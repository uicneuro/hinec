"""Known-answer tests for scripts/tract_assessment (run: pytest tests/python).

Every test builds tractograms whose correct comparison is known analytically,
including the two counterexamples the assessment exists to expose: a
translated copy (same geometry, different place) and a manufactured step-size
error (known convergence order).
"""
import sys
from pathlib import Path

import h5py
import nibabel as nib
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))

from tract_assessment import studies                                  # noqa: E402
from tract_assessment.compare import Features, Grid, compare           # noqa: E402
from tract_assessment.geometry import PRIMARY, connection, track_summary  # noqa: E402
from tract_assessment.io import Tractogram, load_tractogram, to_world  # noqa: E402
from tract_assessment.report import reading                           # noqa: E402

AFFINE = np.diag([2.0, 2.0, 2.0, 1.0])


def helix(r=6.0, c=2.0, turns=1.5, n=600, centre=(40, 40, 20), phase=0.0):
    t = np.linspace(0, 2 * np.pi * turns, n) + phase
    return np.column_stack([centre[0] + r * np.cos(t), centre[1] + r * np.sin(t), centre[2] + c * t])


def make(label, world_curves, seeds_world=None, affine=AFFINE):
    """Tractogram from world-mm curves (stored as one-based voxel coordinates)."""
    inv = np.linalg.inv(affine)
    vox = lambda w: (np.c_[w, np.ones(len(w))] @ inv.T)[:, :3] + 1.0
    tracks = [vox(c) for c in world_curves]
    seeds = None if seeds_world is None else np.array([vox(s[None])[0] for s in seeds_world])
    return Tractogram(label=label, track_file=Path(label), reference=Path('ref'), affine=affine,
                      shape=(64, 64, 64), tracks=tracks, seeds=seeds, n_seeds=len(tracks))


def bundle(n=12, shift=(0, 0, 0), seed_index=300, **kw):
    curves = [helix(centre=(40 + (i % 4) * 1.5, 40 + (i // 4) * 1.5, 20), phase=0.1 * i, **kw) + np.array(shift)
              for i in range(n)]
    return curves, [c[seed_index] for c in curves]


def card(a, b):
    grid = Grid([a, b])
    return compare(Features(a, grid), Features(b, grid))[0]


# ----------------------------------------------------------------- geometry

def test_helix_curvature_and_torsion_match_closed_form():
    r, c = 6.0, 2.0
    kappa, tau = r / (r * r + c * c), c / (r * r + c * c)
    m = connection(helix(r, c, turns=3, n=3000), ds=.25, window_mm=4)
    core = slice(len(m['s']) // 4, -len(m['s']) // 4)
    assert np.nanmedian(m['kappa_w'][core]) == pytest.approx(kappa, rel=.02)
    assert np.nanmedian(abs(m['tau'][core])) == pytest.approx(tau, rel=.02)
    # Frenet gauge on a helix: w13 ~ 0, w23 ~ tau
    assert np.nanmedian(abs(m['w13'][core])) < 1e-3
    assert np.nanmedian(abs(m['w23'][core])) == pytest.approx(tau, rel=.02)


def test_short_track_is_an_exclusion_not_a_zero():
    assert track_summary(np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0.]]), PRIMARY) is None


# ----------------------------------------------------------------- engine

def test_identical_tractograms_agree_perfectly():
    curves, seeds = bundle()
    c = card(make('a', curves, seeds), make('b', curves, seeds))
    assert c['spatial']['occupied_voxel_jaccard'] == 1.0
    assert c['spatial']['endpoint_voxel_jaccard'] == 1.0
    assert c['paths']['correspondence'] == 'seed' and c['paths']['pairs'] == 12
    assert c['paths']['separation_median_mm']['p95'] == pytest.approx(0, abs=1e-9)
    g = c['geometry'][PRIMARY.tag]['kappa_w']
    assert g['relative_median_change'] == pytest.approx(0, abs=1e-12)
    assert g['paired']['spearman'] == pytest.approx(1.0)


def test_translated_copy_has_same_geometry_but_different_place():
    """The study's counterexample: shape statistics cannot certify location."""
    curves, seeds = bundle()
    moved, moved_seeds = bundle(shift=(0, 0, 25))
    c = card(make('a', curves, seeds), make('b', moved, moved_seeds))
    g = c['geometry'][PRIMARY.tag]['kappa_w']
    assert g['relative_median_change'] < 1e-9
    assert c['spatial']['occupied_voxel_jaccard'] < 0.05
    text = ' '.join(reading(studies.headline(c)))
    assert 'geometry agrees but location does not' in text


def test_reversed_storage_order_is_matched_to_the_right_half():
    curves, seeds = bundle()
    c = card(make('a', curves, seeds), make('b', [x[::-1] for x in curves], seeds))
    assert c['paths']['separation_median_mm']['p95'] == pytest.approx(0, abs=1e-9)


def test_nearest_pairing_without_seeds_recovers_shuffled_correspondence():
    curves, _ = bundle()
    order = np.random.default_rng(0).permutation(len(curves))
    a = make('a', curves)
    b = make('b', [curves[k][::-1] for k in order])
    grid = Grid([a, b])
    c, rows = compare(Features(a, grid), Features(b, grid))
    assert c['paths']['correspondence'] == 'nearest'
    inverse = np.argsort(order)
    assert all(r['track_b'] - 1 == inverse[r['track_a'] - 1] for r in rows)
    assert c['paths']['separation_median_mm']['p95'] == pytest.approx(0, abs=1e-9)


def test_convergence_recovers_manufactured_order():
    """Level h perturbs every curve by C*h^2: observed order must be ~2."""
    base, seeds = bundle()
    levels = [0.4, 0.2, 0.1, 0.05]
    ladder = []
    for h in levels + [0.0]:
        curves = []
        for c in base:
            s = np.linspace(0, 1, len(c))
            curves.append(c + (h ** 2) * np.column_stack([np.sin(3 * s), s, 0 * s]))
        ladder.append(make(f'h{h}', curves, seeds))
    result, _ = studies.convergence(ladder, levels + [0.0], 'integrator.step')
    orders = [o['order'] for o in result['observed_order']]
    # The coarsest pair is pre-asymptotic (the perturbation also changes arc
    # length, which arc-aligned comparison sees); the refined pairs are exact.
    assert np.allclose(orders[1:], 2.0, atol=.01), orders
    assert result['reference']['independent_of_ladder'] is False


def test_world_grid_when_affines_differ():
    curves, seeds = bundle()
    a = make('a', curves, seeds)
    b = make('b', curves, seeds, affine=np.diag([1.5, 1.5, 1.5, 1.0]))
    grid = Grid([a, b])
    assert grid.kind == 'world'
    c = compare(Features(a, grid), Features(b, grid))[0]
    assert c['spatial']['occupied_voxel_jaccard'] == 1.0     # same world curves


# ----------------------------------------------------------------- loader

def test_load_matlab_v73_layout_and_affine(tmp_path):
    affine = np.array([[-1.5, 0, 0, 90], [0, 1.5, 0, -120], [0, 0, 2.0, -60], [0, 0, 0, 1]])
    nib.save(nib.Nifti1Image(np.zeros((20, 20, 20, 2), np.float32), affine), tmp_path / 'dwi.nii.gz')
    tracks = [np.column_stack([np.linspace(2, 9, n), np.full(n, 5.0), np.full(n, 6.0)]) for n in (3, 7)]
    seeds = np.array([t[1] for t in tracks])
    with h5py.File(tmp_path / 'tracks_hinec_x.mat', 'w') as f:
        refs = f.create_group('#refs#')
        ref_ids = []
        for k, t in enumerate(tracks):
            ref_ids.append(refs.create_dataset(str(k), data=t.T).ref)   # MATLAB N x 3 -> HDF5 3 x N
        f.create_dataset('tracks', data=np.array([ref_ids], dtype=h5py.ref_dtype))
        f.create_dataset('track_meta/seed_points', data=seeds.T)
        f.create_dataset('track_meta/n_seeds', data=np.array([[4.0]]))
        mask = np.zeros((20, 20, 20), bool)
        mask[4, 4, 5] = True
        f.create_dataset('options/seed_mask', data=mask.T.astype(np.uint8))
    t = load_tractogram(tmp_path / 'tracks_hinec_x.mat', reference=tmp_path / 'dwi.nii.gz')
    assert [len(x) for x in t.tracks] == [3, 7]                 # a 3-point track is not transposed
    assert np.allclose(t.tracks[1], tracks[1])
    assert np.allclose(t.seeds, seeds) and t.n_seeds == 4
    assert t.seed_mask_voxels == {(5, 5, 6)}
    assert np.allclose(t.world[0][0], affine[:3, :3] @ (tracks[0][0] - 1) + affine[:3, 3])
    assert np.allclose(to_world(np.array([[1, 1, 1.]]), affine), affine[:3, 3])
