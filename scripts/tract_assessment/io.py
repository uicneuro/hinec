"""Load a saved tractogram from any HINEC tracker into one common form.

Accepts a run directory (hinec_runs/run_*/) or a tracks .mat file. Tracks are
stored as one-based voxel coordinates; world millimetres are
    world = affine @ (p - 1)
using the full NIfTI affine of a reference image on the DWI grid (the same
convention as scripts/hinec_to_trk.py).

Seed identity is optional. Trackers that record meta.seed_points (hinec,
template, stitching) allow exact per-seed pairing; trackers that do not
(standard, mmf) are paired by nearest streamline instead (see compare.py).
"""
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import nibabel as nib


@dataclass
class Tractogram:
    label: str
    track_file: Path
    reference: Path
    affine: np.ndarray
    shape: tuple
    tracks: list                      # one-based voxel coordinates, N x 3 each
    seeds: np.ndarray = None          # T x 3 one-based voxel seed of each track, or None
    n_seeds: int = None               # seeds attempted (denominator for yield)
    seed_mask_voxels: set = None      # effective seed voxels, one-based
    run_dir: Path = None
    config: dict = field(default_factory=dict)
    sha256: str = ''
    _world: list = None

    @property
    def world(self):
        if self._world is None:
            self._world = [to_world(t, self.affine) for t in self.tracks]
        return self._world

    @property
    def voxel_mm(self):
        return np.linalg.norm(self.affine[:3, :3], axis=0)


def to_world(points, affine):
    points = np.asarray(points, float)
    return (np.c_[points - 1.0, np.ones(len(points))] @ np.asarray(affine).T)[:, :3]


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _read_mat(path):
    """Return dict with tracks, seeds (per track or per seed), n_seeds, seed_mask."""
    import h5py
    out = {}
    try:
        mat = h5py.File(path, 'r')
    except OSError:
        return _read_mat_v7(path)
    with mat:
        refs = np.asarray(mat['tracks']).ravel()
        tracks = []
        for ref in refs:
            arr = np.asarray(mat[ref])
            if arr.ndim != 2 or arr.shape[0] != 3:       # empty MATLAB cell
                continue
            tracks.append(arr.T.astype(float))           # HDF5 stores MATLAB arrays transposed
        out['tracks'] = tracks
        for group in ('track_meta', ''):
            prefix = group + '/' if group else ''
            if prefix + 'seed_points' in mat:
                seeds = np.asarray(mat[prefix + 'seed_points'], float)
                if seeds.ndim == 2 and seeds.shape[0] == 3:
                    out['seeds'] = seeds.T
                if prefix + 'seed_index' in mat:
                    out['seed_index'] = np.asarray(mat[prefix + 'seed_index']).ravel().astype(int)
                break
        if 'track_meta/n_seeds' in mat:
            out['n_seeds'] = int(np.asarray(mat['track_meta/n_seeds']).ravel()[0])
        if 'options/seed_mask' in mat:
            mask = np.asarray(mat['options/seed_mask']).T.astype(bool)
            out['seed_mask'] = {tuple(v + 1) for v in np.argwhere(mask)}
    return out


def _read_mat_v7(path):
    from scipy.io import loadmat
    mat = loadmat(path, simplify_cells=True)
    tracks = [np.asarray(t, float) for t in np.atleast_1d(mat['tracks']) if np.size(t) >= 3]
    out = {'tracks': [t.reshape(-1, 3) for t in tracks]}
    meta = mat.get('track_meta') or {}
    if isinstance(meta, dict) and 'seed_points' in meta:
        out['seeds'] = np.asarray(meta['seed_points'], float).reshape(-1, 3)
        if 'n_seeds' in meta:
            out['n_seeds'] = int(meta['n_seeds'])
    opts = mat.get('options') or {}
    if isinstance(opts, dict) and 'seed_mask' in opts:
        out['seed_mask'] = {tuple(v + 1) for v in np.argwhere(np.asarray(opts['seed_mask'], bool))}
    return out


def _seeds_per_track(tracks, seeds, seed_index=None, tol=1e-6):
    """Return a T x 3 seed per kept track, or None if seed identity is not
    recoverable. Handles per-track seed lists and per-seed lists."""
    if seeds is None or not len(tracks):
        return None
    if len(seeds) == len(tracks):
        per_track = seeds
    elif seed_index is not None and len(seed_index) == len(tracks) and seed_index.max() <= len(seeds):
        per_track = seeds[seed_index - 1]
    else:
        # Per-seed list with dropped seeds: assign each track to the seed it passes through.
        from scipy.spatial import cKDTree
        tree = cKDTree(seeds)
        per_track = np.empty((len(tracks), 3))
        for i, t in enumerate(tracks):
            d, j = tree.query(t)
            per_track[i] = seeds[j[int(d.argmin())]]
    for t, s in zip(tracks, per_track):
        if np.linalg.norm(t - s, axis=1).min() > tol:
            return None      # recorded seed is not on the streamline: do not trust it
    return per_track


def _find_reference(run_dir, track_file):
    inter = run_dir / 'intermediate'
    source = inter / 'SOURCE.txt'
    if source.exists():
        m = re.search(r'dwi_reference:\s*(\S+)', source.read_text())
        if m and (inter / Path(m.group(1)).name).exists():
            return inter / Path(m.group(1)).name
    candidates = sorted(inter.glob('*.nii.gz')) + sorted(inter.glob('*.nii'))
    output_names = {p.stem for p in (run_dir / 'output').glob('*.mat')}
    for c in candidates:
        if c.name.split('.nii')[0] in output_names:
            return c
    four_d = [c for c in candidates if len(nib.load(c).shape) == 4]
    if four_d:
        return four_d[0]
    raise FileNotFoundError(f'No reference NIfTI under {inter}; pass --reference')


def _read_config(run_dir):
    info = {}
    cfg = run_dir / 'config.yml'
    if cfg.exists():
        info['config_file'] = str(cfg)
    overrides = run_dir / 'overrides.txt'
    if overrides.exists():
        info['overrides'] = [l.strip() for l in overrides.read_text().splitlines() if l.strip()]
    return info


def load_tractogram(path, reference=None, label=None, tracks_file=None):
    """Load a run directory or tracks .mat.

    path        run dir, or a .mat with a `tracks` cell array
    reference   NIfTI on the tracking grid (required for a bare .mat unless it
                sits in a run dir's tractography/ folder)
    tracks_file which tracks_*.mat in a run dir, if more than one
    """
    path = Path(path).resolve()
    run_dir = None
    if path.is_dir():
        run_dir = path
        if tracks_file:
            track_file = run_dir / 'tractography' / tracks_file
        else:
            found = sorted((run_dir / 'tractography').glob('tracks_*.mat'))
            if not found:
                raise FileNotFoundError(f'No tractography/tracks_*.mat in {run_dir}')
            if len(found) > 1:
                print(f'[assess] {run_dir.name}: {len(found)} track files, using newest {found[-1].name}')
            track_file = found[-1]
    else:
        track_file = path
        if path.parent.name == 'tractography':
            run_dir = path.parent.parent
    if reference is None:
        if run_dir is None:
            raise ValueError(f'{path}: a bare .mat needs --reference <nifti on the tracking grid>')
        reference = _find_reference(run_dir, track_file)
    reference = Path(reference).resolve()
    image = nib.load(reference)
    raw = _read_mat(track_file)
    tracks = raw['tracks']
    bad = [i for i, t in enumerate(tracks) if not np.all(np.isfinite(t)) or np.iscomplexobj(t)]
    if bad:
        raise ValueError(f'{track_file}: {len(bad)} tracks contain non-finite coordinates')
    seeds = _seeds_per_track(tracks, raw.get('seeds'), raw.get('seed_index'))
    n_seeds = raw.get('n_seeds')
    if n_seeds is None and raw.get('seeds') is not None and len(raw['seeds']) >= len(tracks):
        n_seeds = len(raw['seeds'])
    if label is None:
        label = run_dir.name if run_dir else track_file.stem
    return Tractogram(label=label, track_file=track_file, reference=reference,
                      affine=np.asarray(image.affine, float), shape=tuple(image.shape[:3]),
                      tracks=tracks, seeds=seeds, n_seeds=n_seeds,
                      seed_mask_voxels=raw.get('seed_mask'), run_dir=run_dir,
                      config=_read_config(run_dir) if run_dir else {},
                      sha256=_sha256(track_file))
