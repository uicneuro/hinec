"""The three study designs, all built from compare.compare().

  convergence    ONE input, ONE pipeline, a numerical knob refined (step size,
                 seed density, ...). Every level is compared with the finest
                 level. Answers: does the output settle, how fast, and does it
                 settle for every streamline or only for the typical one?
  repeatability  the SAME pipeline on independent repeats of the same thing
                 (scan-rescan acquisitions, or random-seed reruns). Answers: how
                 much of the output survives re-acquisition?
  consistency    DIFFERENT pipelines/settings on the SAME input (trilinear vs
                 cubic, hinec vs mmf, DTI vs CSD). Answers: how much does the
                 output depend on a methodological choice?

The engine is identical; only which pairs are compared, and how the result is
read, differ. See docs/TRACT_ASSESSMENT.md for interpretation.
"""
import itertools

import numpy as np

from .compare import Features, Grid, compare
from .geometry import PRIMARY, SENSITIVITY


def _features(tractograms, grid_mm, max_tracks, sensitivity):
    grid = Grid(tractograms, grid_mm)
    protocols = (PRIMARY,) + (SENSITIVITY if sensitivity else ())
    feats = []
    for t in tractograms:
        print(f'[assess] measuring {t.label}: {len(t.tracks)} tracks', flush=True)
        feats.append(Features(t, grid, protocols, max_tracks))
    return grid, feats


def headline(card):
    """The handful of numbers that summarise one comparison, one per tier."""
    g = card['geometry'][PRIMARY.tag]
    p = card['paths']
    kappa = g['kappa_w']
    return {
        'a': card['a'], 'b': card['b'],
        'tracks_a': card['yield']['tracks'][0], 'tracks_b': card['yield']['tracks'][1],
        'occupied_jaccard': card['spatial']['occupied_voxel_jaccard'],
        'endpoint_jaccard': card['spatial']['endpoint_voxel_jaccard'],
        'weighted_dice': card['spatial']['density_weighted_dice'],
        'correspondence': p['correspondence'], 'pairs': p['pairs'],
        'path_separation_median_mm': p['separation_median_mm']['median'],
        'path_separation_tail_mm': p['separation_median_mm']['p95'],
        'fraction_pairs_p95_over_1mm': p['fraction_p95_over_1mm'],
        'endpoint_distance_median_mm': p['endpoint_distance_mm']['median'],
        'fraction_pairs_diverged': p['fraction_diverged'],
        'length_rel_median_change': g['length_mm']['relative_median_change'],
        'kappa_w_rel_median_change': kappa['relative_median_change'],
        'kappa_w_paired_spearman': kappa.get('paired', {}).get('spearman'),
        'tau_rel_median_change': g['tau']['relative_median_change'],
        'tau_paired_spearman': g['tau'].get('paired', {}).get('spearman'),
        'w23_rel_median_change': g['w23']['relative_median_change'],
    }


def convergence(tractograms, values, parameter, reference=None, grid_mm=None, max_tracks=None,
                max_pairs=None, sensitivity=False):
    """Compare every level with the reference (default: the finest = last).

    values are the knob settings, ordered coarse -> fine. The reference is part
    of the ladder unless given separately; say so when reporting, because
    self-convergence is not accuracy."""
    ladder = list(tractograms)
    independent = reference is not None
    if reference is None:
        reference = ladder[-1]
        levels, level_values = ladder[:-1], list(values[:-1])
    else:
        levels, level_values = ladder, list(values)
    grid, feats = _features(levels + [reference], grid_mm, max_tracks, sensitivity)
    ref = feats[-1]
    cards, rows, heads = [], [], []
    for value, f in zip(level_values, feats[:-1]):
        card, r = compare(f, ref, max_pairs)
        card['level'] = value
        cards.append(card)
        rows.extend(dict(level=value, **x) for x in r)
        heads.append(dict(level=value, **headline(card)))
    # Cauchy check between consecutive levels: does each refinement change less?
    successive = []
    for k in range(len(level_values) - 1):
        card, _ = compare(feats[k], feats[k + 1], max_pairs, keep_rows=False)
        successive.append({'from': level_values[k], 'to': level_values[k + 1], **headline(card)})
    orders = _observed_orders(level_values, [h['path_separation_median_mm'] for h in heads])
    return {'study': 'convergence', 'parameter': parameter, 'values': list(values),
            'reference': {'label': reference.label, 'value': None if independent else values[-1],
                          'independent_of_ladder': independent},
            'grid': grid.describe(), 'cards': cards, 'headline': heads,
            'successive': successive, 'observed_order': orders}, rows


def _observed_orders(values, errors):
    try:
        h = np.asarray(values, float)
    except (TypeError, ValueError):
        return None                      # categorical knob: no order
    e = np.asarray([np.nan if x is None else x for x in errors], float)
    out = []
    for k in range(len(h) - 1):
        if e[k] > 0 and e[k + 1] > 0 and h[k] != h[k + 1]:
            out.append({'between': [float(h[k]), float(h[k + 1])],
                        'order': float(np.log(e[k] / e[k + 1]) / np.log(h[k] / h[k + 1]))})
    return out


def repeatability(groups, grid_mm=None, max_tracks=None, max_pairs=None, sensitivity=False):
    """groups: {name: [tractogram, ...]} — every pair within a group is compared.
    The GROUP (subject/site/session) is the independent unit, not the streamline."""
    cards, rows, heads = [], [], []
    for name, members in groups.items():
        grid, feats = _features(members, grid_mm, max_tracks, sensitivity)
        for fa, fb in itertools.combinations(feats, 2):
            card, r = compare(fa, fb, max_pairs)
            card['group'] = name
            card['grid'] = grid.describe()
            cards.append(card)
            rows.extend(dict(group=name, **x) for x in r)
            heads.append(dict(group=name, **headline(card)))
    return {'study': 'repeatability', 'groups': {k: [t.label for t in v] for k, v in groups.items()},
            'cards': cards, 'headline': heads, 'across_groups': _ranges(heads)}, rows


def consistency(tractograms, grid_mm=None, max_tracks=None, max_pairs=None, sensitivity=False):
    """All pairwise comparisons between methods run on one input."""
    grid, feats = _features(tractograms, grid_mm, max_tracks, sensitivity)
    cards, rows, heads = [], [], []
    for fa, fb in itertools.combinations(feats, 2):
        card, r = compare(fa, fb, max_pairs)
        cards.append(card)
        rows.extend(r)
        heads.append(headline(card))
    return {'study': 'consistency', 'methods': [t.label for t in tractograms], 'grid': grid.describe(),
            'cards': cards, 'headline': heads, 'across_pairs': _ranges(heads)}, rows


def _ranges(heads):
    out = {}
    for key in heads[0] if heads else []:
        vals = [h[key] for h in heads if isinstance(h[key], (int, float)) and h[key] is not None]
        if vals and key not in ('tracks_a', 'tracks_b', 'pairs'):
            out[key] = {'min': float(min(vals)), 'median': float(np.median(vals)), 'max': float(max(vals))}
    return out
