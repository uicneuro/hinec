"""Compare two tractograms: the single engine behind every study type.

A comparison produces one "agreement card" with four tiers, reported
separately because they disagree in practice (docs/TRACT_ASSESSMENT.md):

  yield        did both runs produce streamlines from the seeds they tried?
  spatial      do they occupy the same voxels and end in the same places?
  paths        do corresponding streamlines follow the same route?
  geometry     do they have the same MMF (moving-frame) shape statistics?

Stable geometry with unstable space is the characteristic failure the
tract-geometry study found, so no tier may stand in for another.

Streamline correspondence:
  'seed'     both runs record the seed of every streamline (hinec, template,
             stitching): pair by identical seed coordinate. Exact.
  'nearest'  otherwise: pair each streamline with its nearest neighbour by
             mean direct-flip distance (MDF). A heuristic - a pair is "the most
             similar streamline", not "the same seed" - and is labelled so.
"""
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wasserstein_distance

from .geometry import PRIMARY, GEOMETRY_METRICS, track_summary

SEPARATION_THRESHOLD_MM = 1.0
NEAREST_POINTS = 20
# Deliberately small: this machine is shared. Raise with --workers when it is free.
WORKERS = max(1, min(4, (os.cpu_count() or 2) // 2))


def _summaries(job):
    curves, protocol = job
    return [track_summary(c, protocol) for c in curves]


# ----------------------------------------------------------------- helpers

def _length(curve):
    return float(np.linalg.norm(np.diff(curve, axis=0), axis=1).sum()) if len(curve) > 1 else 0.0


def _quantiles(values, qs=(0.5, 0.95)):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {('median' if q == .5 else f'p{int(q*100)}'): None for q in qs}
    return {('median' if q == .5 else f'p{int(q*100)}'): float(np.quantile(values, q)) for q in qs}


def _jaccard(a, b):
    union = a | b
    return len(a & b) / len(union) if union else 1.0


def _dice(a, b):
    return 2 * len(a & b) / (len(a) + len(b)) if (a or b) else 1.0


def _densify(curve, step):
    """Insert points so no segment exceeds `step` (so voxel occupancy has no gaps)."""
    if len(curve) < 2:
        return curve
    seg = np.linalg.norm(np.diff(curve, axis=0), axis=1)
    n = np.maximum(1, np.ceil(seg / step).astype(int))
    if n.max() == 1:
        return curve
    parts = [curve[i] + np.outer(np.arange(k) / k, curve[i + 1] - curve[i]) for i, k in enumerate(n)]
    return np.vstack(parts + [curve[-1:]])


def _arc_resample(curve, n):
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(curve, axis=0), axis=1))]
    if s[-1] <= 0:
        return np.repeat(curve[:1], n, axis=0)
    t = np.linspace(0, s[-1], n)
    return np.column_stack([np.interp(t, s, curve[:, k]) for k in range(3)])


def subsample_indices(n, cap, salt=0):
    """Deterministic, evenly spread subset (reproducible across runs)."""
    if cap is None or n <= cap:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, cap).round().astype(int))


# ----------------------------------------------------------------- grids

class Grid:
    """Where spatial overlap is counted. 'native' when both runs share one
    voxel grid (exactly the tracking voxels); otherwise an isotropic world-mm
    grid, which is the only fair common ground for different acquisitions."""

    def __init__(self, tractograms, spacing_mm=None):
        first = tractograms[0]
        same = all(t.shape == first.shape and np.array_equal(t.affine, first.affine) for t in tractograms)
        if same and spacing_mm is None:
            self.kind, self.spacing = 'native', float(first.voxel_mm.min())
            self._inv = np.linalg.inv(first.affine)
        else:
            self.kind = 'world'
            self.spacing = float(spacing_mm or min(t.voxel_mm.min() for t in tractograms))
            self._inv = None

    def describe(self):
        return {'kind': self.kind, 'spacing_mm': self.spacing}

    def voxels(self, world_curve):
        dense = _densify(world_curve, self.spacing / 2)
        if self.kind == 'native':
            ijk = (np.c_[dense, np.ones(len(dense))] @ self._inv.T)[:, :3]
        else:
            ijk = dense / self.spacing
        return np.rint(ijk).astype(np.int64)


# ----------------------------------------------------------------- per-run features (cached)

class Features:
    """Everything about one tractogram that does not depend on its partner."""

    def __init__(self, tractogram, grid, protocols=(PRIMARY,), max_tracks=None):
        self.t = tractogram
        world = tractogram.world
        self.lengths = np.array([_length(c) for c in world])
        occupied, endpoints, density = set(), set(), {}
        for c in world:
            vox = grid.voxels(c)
            keys = set(map(tuple, vox))
            occupied |= keys
            for k in keys:
                density[k] = density.get(k, 0) + 1
            endpoints.add(tuple(vox[0]))
            endpoints.add(tuple(vox[-1]))
        self.occupied, self.endpoints, self.density = occupied, endpoints, density
        # Geometry: the distribution uses a deterministic subset (all tracks
        # unless capped); paired comparisons measure whatever tracks they need.
        self.protocols = tuple(protocols)
        self.measured = subsample_indices(len(world), max_tracks)
        self._geometry = {p.tag: {} for p in self.protocols}
        self.ensure_geometry(self.measured)
        self.seed_keys = None
        if tractogram.seeds is not None:
            self.seed_keys = [tuple(np.round(s, 6)) for s in tractogram.seeds]
        self._nearest = None

    def ensure_geometry(self, indices):
        world = self.t.world
        for p in self.protocols:
            cache = self._geometry[p.tag]
            missing = [int(i) for i in indices if int(i) not in cache]
            if len(missing) < 400 or WORKERS == 1:
                for i in missing:
                    cache[i] = track_summary(world[i], p)
                continue
            chunks = [missing[k:k + 200] for k in range(0, len(missing), 200)]
            with ProcessPoolExecutor(WORKERS) as pool:
                for chunk, result in zip(chunks, pool.map(_summaries, [([world[i] for i in c], p) for c in chunks])):
                    cache.update(zip(chunk, result))

    def geometry(self, tag, indices):
        """(values [len(indices) x metrics], excluded mask) for one protocol."""
        self.ensure_geometry(indices)
        cache = self._geometry[tag]
        values = np.full((len(indices), len(GEOMETRY_METRICS)), np.nan)
        excluded = np.zeros(len(indices), bool)
        for r, i in enumerate(indices):
            summary = cache[int(i)]
            if summary is None:
                excluded[r] = True
            else:
                values[r] = [summary[m] for m in GEOMETRY_METRICS]
        return values, excluded

    def nearest_signature(self):
        if self._nearest is None:
            self._nearest = np.stack([_arc_resample(c, NEAREST_POINTS) for c in self.t.world])
        return self._nearest


# ----------------------------------------------------------------- correspondence

def _pairs_by_seed(fa, fb):
    def unique(keys):
        seen, dup = {}, set()
        for i, k in enumerate(keys):
            if k in seen:
                dup.add(k)
            seen[k] = i
        return {k: i for k, i in seen.items() if k not in dup}, len(dup)
    ua, da = unique(fa.seed_keys)
    ub, db = unique(fb.seed_keys)
    common = sorted(ua.keys() & ub.keys())
    return [(ua[k], ub[k]) for k in common], {'ambiguous_seed_keys': da + db}


def _mdf(a, b):
    direct = np.linalg.norm(a - b, axis=1).mean()
    flipped = np.linalg.norm(a - b[::-1], axis=1).mean()
    return min(direct, flipped)


def _pairs_by_nearest(fa, fb, cap):
    sa, sb = fa.nearest_signature(), fb.nearest_signature()
    ia = subsample_indices(len(sa), cap)
    if not len(sb):
        return [], {}
    flat_b = sb.reshape(len(sb), -1)
    flat_b_flip = sb[:, ::-1].reshape(len(sb), -1)
    tree = cKDTree(np.vstack([flat_b, flat_b_flip]))
    k = min(8, 2 * len(sb))
    pairs, distances = [], []
    for i in ia:
        _, cand = tree.query(sa[i].reshape(-1), k=k)
        cand = np.unique(np.atleast_1d(cand) % len(sb))
        d = [_mdf(sa[i], sb[j]) for j in cand]
        j = int(cand[int(np.argmin(d))])
        pairs.append((int(i), j))
        distances.append(float(min(d)))
    return pairs, {'nearest_mdf_mm': distances}


# ----------------------------------------------------------------- path agreement

def _split_at_seed(curve_vox, curve_world, seed):
    i = int(np.linalg.norm(curve_vox - seed, axis=1).argmin())
    return curve_world[i::-1], curve_world[i:]


def _arc_param(side):
    return np.r_[0, np.cumsum(np.linalg.norm(np.diff(side, axis=0), axis=1))]


def _separation(side_a, side_b, step):
    """Positions at equal arc distance from the seed; returns (s, distance)."""
    sa, sb = _arc_param(side_a), _arc_param(side_b)
    support = min(sa[-1], sb[-1])
    if support <= 0:
        return np.zeros(1), np.zeros(1)
    s = np.arange(0, support + 1e-9, step)
    pa = np.column_stack([np.interp(s, sa, side_a[:, k]) for k in range(3)])
    pb = np.column_stack([np.interp(s, sb, side_b[:, k]) for k in range(3)])
    return s, np.linalg.norm(pa - pb, axis=1)


def seed_path_agreement(ta, tb, i, j, step=0.25, threshold=SEPARATION_THRESHOLD_MM):
    """Arc-length-aligned comparison of two streamlines grown from one seed.

    Each half is compared from the seed outward at equal arc length (tracker
    agnostic: no assumption about step size or integration variable). Half
    assignment (which end is 'backward') is chosen by early agreement, because
    two runs can legitimately pick opposite initial directions.
    """
    halves_a = _split_at_seed(ta.tracks[i], ta.world[i], ta.seeds[i])
    halves_b = _split_at_seed(tb.tracks[j], tb.world[j], tb.seeds[j])
    options = []
    for order in ((0, 1), (1, 0)):
        comps = [_separation(halves_a[k], halves_b[order[k]], step) for k in range(2)]
        early = np.mean([d[min(len(d) - 1, int(2.0 / step))] for _, d in comps])
        options.append((early, order, comps))
    _, order, comps = min(options, key=lambda x: x[0])
    all_d = np.concatenate([d for _, d in comps])
    divergence = []
    for s, d in comps:
        over = np.flatnonzero(d > threshold)
        divergence.append(float(s[over[0]]) if len(over) else np.nan)   # NaN = never diverged
    endpoint = [float(np.linalg.norm(halves_a[k][-1] - halves_b[order[k]][-1])) for k in range(2)]
    return {'separation_median_mm': float(np.median(all_d)),
            'separation_p95_mm': float(np.quantile(all_d, .95)),
            'separation_max_mm': float(all_d.max()),
            'diverged': bool(np.any(np.isfinite(divergence))),
            'first_divergence_arc_mm': float(np.nanmin(divergence)) if np.any(np.isfinite(divergence)) else None,
            'endpoint_distance_mm': float(max(endpoint)),
            'length_difference_mm': abs(_length(ta.world[i]) - _length(tb.world[j]))}


def nearest_path_agreement(fa, fb, i, j, mdf):
    a, b = fa.t.world[i], fb.t.world[j]
    ends_direct = max(np.linalg.norm(a[0] - b[0]), np.linalg.norm(a[-1] - b[-1]))
    ends_flip = max(np.linalg.norm(a[0] - b[-1]), np.linalg.norm(a[-1] - b[0]))
    return {'separation_median_mm': mdf, 'separation_p95_mm': None, 'separation_max_mm': None,
            'diverged': bool(mdf > SEPARATION_THRESHOLD_MM), 'first_divergence_arc_mm': None,
            'endpoint_distance_mm': float(min(ends_direct, ends_flip)),
            'length_difference_mm': abs(fa.lengths[i] - fb.lengths[j])}


# ----------------------------------------------------------------- the card

def _distribution_change(one, two):
    one, two = one[np.isfinite(one)], two[np.isfinite(two)]
    if len(one) < 2 or len(two) < 2:
        return {'medians': None, 'relative_median_change': None, 'wasserstein_over_iqr': None}
    m1, m2 = float(np.median(one)), float(np.median(two))
    pooled = np.r_[one, two]
    iqr = float(np.quantile(pooled, .75) - np.quantile(pooled, .25))
    return {'medians': [m1, m2],
            'relative_median_change': abs(m1 - m2) / max((abs(m1) + abs(m2)) / 2, 1e-12),
            'wasserstein_over_iqr': float(wasserstein_distance(one, two) / max(iqr, 1e-12))}


def _paired(one, two):
    ok = np.isfinite(one) & np.isfinite(two)
    if ok.sum() < 3:
        return {'n': int(ok.sum()), 'spearman': None, 'median_abs_difference': None}
    rho = spearmanr(one[ok], two[ok]).statistic
    return {'n': int(ok.sum()), 'spearman': None if not np.isfinite(rho) else float(rho),
            'median_abs_difference': float(np.median(abs(one[ok] - two[ok])))}


def compare(fa, fb, max_pairs=None, keep_rows=True):
    """Agreement card for tractograms A and B (Features objects on one Grid)."""
    ta, tb = fa.t, fb.t
    card = {'a': ta.label, 'b': tb.label}

    # yield
    card['yield'] = {
        'tracks': [len(ta.tracks), len(tb.tracks)],
        'seeds_attempted': [ta.n_seeds, tb.n_seeds],
        'fraction': [len(t.tracks) / t.n_seeds if t.n_seeds else None for t in (ta, tb)],
    }
    if ta.seed_mask_voxels is not None and tb.seed_mask_voxels is not None:
        card['yield']['seed_region_dice'] = _dice(ta.seed_mask_voxels, tb.seed_mask_voxels)

    # spatial
    common = fa.density.keys() & fb.density.keys()
    na, nb = sum(fa.density.values()), sum(fb.density.values())
    weighted = (sum(fa.density[k] / na + fb.density[k] / nb for k in common) /
                max(sum(v / na for v in fa.density.values()) + sum(v / nb for v in fb.density.values()), 1e-12)
                if na and nb else None)
    card['spatial'] = {
        'occupied_voxel_jaccard': _jaccard(fa.occupied, fb.occupied),
        'endpoint_voxel_jaccard': _jaccard(fa.endpoints, fb.endpoints),
        'density_weighted_dice': weighted,
        'coverage_of_a_by_b': len(fa.occupied & fb.occupied) / max(len(fa.occupied), 1),
        'coverage_of_b_by_a': len(fa.occupied & fb.occupied) / max(len(fb.occupied), 1),
    }

    # correspondence + paths
    if fa.seed_keys is not None and fb.seed_keys is not None:
        mode = 'seed'
        pairs, extra = _pairs_by_seed(fa, fb)
        if max_pairs and len(pairs) > max_pairs:
            pairs = [pairs[k] for k in subsample_indices(len(pairs), max_pairs)]
        paths = [seed_path_agreement(ta, tb, i, j) for i, j in pairs]
    else:
        mode = 'nearest'
        pairs, extra = _pairs_by_nearest(fa, fb, max_pairs)
        paths = [nearest_path_agreement(fa, fb, i, j, d) for (i, j), d in zip(pairs, extra.pop('nearest_mdf_mm', []))]
    col = lambda k: np.array([p[k] if p[k] is not None else np.nan for p in paths], float)
    card['paths'] = {
        'correspondence': mode,
        'pairs': len(pairs),
        'separation_median_mm': _quantiles(col('separation_median_mm')),
        'separation_p95_mm': _quantiles(col('separation_p95_mm')),
        'endpoint_distance_mm': _quantiles(col('endpoint_distance_mm')),
        'length_difference_mm': _quantiles(col('length_difference_mm')),
        'fraction_diverged': float(np.mean([p['diverged'] for p in paths])) if paths else None,
        'fraction_p95_over_1mm': (float(np.nanmean(col('separation_p95_mm') > SEPARATION_THRESHOLD_MM))
                                  if mode == 'seed' and paths else None),
        'first_divergence_arc_mm': _quantiles(col('first_divergence_arc_mm')),
        **extra,
    }

    # geometry (distribution + paired)
    card['geometry'] = {}
    ia = np.array([i for i, _ in pairs], int)
    ib = np.array([j for _, j in pairs], int)
    for p in fa.protocols:
        tag = p.tag
        ga, xa = fa.geometry(tag, fa.measured)
        gb, xb = fb.geometry(tag, fb.measured)
        pa, _ = fa.geometry(tag, ia)
        pb, _ = fb.geometry(tag, ib)
        entry = {'excluded_short_tracks': [int(xa.sum()), int(xb.sum())],
                 'measured_tracks': [len(fa.measured), len(fb.measured)]}
        entry['length_mm'] = {**_distribution_change(fa.lengths, fb.lengths),
                              **({'paired': _paired(fa.lengths[ia], fb.lengths[ib])} if len(ia) else {})}
        for m, name in enumerate(GEOMETRY_METRICS):
            entry[name] = _distribution_change(ga[:, m], gb[:, m])
            if len(ia):
                entry[name]['paired'] = _paired(pa[:, m], pb[:, m])
        card['geometry'][tag] = entry

    rows = None
    if keep_rows:
        rows = []
        ga, _ = fa.geometry(fa.protocols[0].tag, ia)
        gb, _ = fb.geometry(fb.protocols[0].tag, ib)
        for r, ((i, j), p) in enumerate(zip(pairs, paths)):
            row = {'a': ta.label, 'b': tb.label, 'track_a': i + 1, 'track_b': j + 1, **p,
                   'length_a_mm': fa.lengths[i], 'length_b_mm': fb.lengths[j]}
            if ta.seeds is not None:
                row.update(seed_x=ta.seeds[i][0], seed_y=ta.seeds[i][1], seed_z=ta.seeds[i][2])
            for m, name in enumerate(GEOMETRY_METRICS):
                row[f'{name}_a'], row[f'{name}_b'] = ga[r, m], gb[r, m]
            rows.append(row)
    return card, rows
