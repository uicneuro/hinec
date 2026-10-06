#!/usr/bin/env python3
"""Assess the convergence, repeatability or consistency of tractography outputs.

Works for any HINEC tracker (standard, hinec, mmf, stitching, template): inputs
are run directories (hinec_runs/run_*) or tracks .mat files. Normally launched
through bin/run_assessment.sh, which can also generate the runs. See
docs/TRACT_ASSESSMENT.md.

  compare        A B                           one agreement card
  convergence    --param KEY --runs R1 R2 ...  each level vs the finest (or --reference)
  repeatability  --group NAME=R1,R2 ...        all pairs within each group
  consistency    --runs R1 R2 ...              all pairs between methods on one input

A bare .mat needs the NIfTI of its tracking grid: pass --image, or write the
input as tracks.mat::grid.nii.gz.
"""
import argparse
import datetime
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tract_assessment import studies, report                     # noqa: E402
from tract_assessment import compare as compare_module           # noqa: E402
from tract_assessment.geometry import PRIMARY, SENSITIVITY       # noqa: E402
from tract_assessment.io import load_tractogram                  # noqa: E402

REPO = Path(__file__).resolve().parent.parent


def _load(spec, image, label=None):
    path, _, ref = spec.partition('::')
    return load_tractogram(path, reference=ref or image, label=label)


def _override_value(tractogram, key):
    """Read KEY's value from a run's overrides.txt (written by run_tractography.sh --set)."""
    for line in tractogram.config.get('overrides', []):
        k, _, v = line.partition('=')
        if k.strip() == key or k.strip().split('.')[-1] == key.split('.')[-1]:
            try:
                return float(v)
            except ValueError:
                return v.strip()
    return None


def _inputs(tractograms):
    return [{'label': t.label, 'track_file': str(t.track_file), 'sha256': t.sha256,
             'reference': str(t.reference), 'tracks': len(t.tracks),
             'seed_identity': t.seeds is not None, **t.config} for t in tractograms]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='study', required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--out', type=Path, help='output directory (default hinec_runs/assess_<time>_<study>)')
    common.add_argument('--name', default='', help='suffix for the default output directory name')
    common.add_argument('--image', help='reference NIfTI for bare .mat inputs')
    common.add_argument('--grid-mm', type=float, help='force a world grid of this spacing for spatial overlap')
    common.add_argument('--max-tracks', type=int, default=20000,
                        help='cap on streamlines measured per tractogram for geometry distributions (deterministic subset)')
    common.add_argument('--max-pairs', type=int, default=5000, help='cap on matched pairs per comparison')
    common.add_argument('--workers', type=int, help='processes for geometry measurement (default 4)')
    common.add_argument('--sensitivity', action='store_true',
                        help='also measure geometry at the two sensitivity protocols (0.25/4 mm, 0.5/2 mm)')
    p = sub.add_parser('compare', parents=[common]); p.add_argument('a'); p.add_argument('b')
    p = sub.add_parser('convergence', parents=[common])
    p.add_argument('--param', required=True, help='the refined knob, e.g. integrator.step')
    p.add_argument('--runs', nargs='+', required=True)
    p.add_argument('--values', nargs='+', help='knob value per run (default: read from each overrides.txt)')
    p.add_argument('--finer', choices=['smaller', 'larger'],
                   help='which direction is finer (default: larger for *density*, else smaller)')
    p.add_argument('--reference', help='independent reference run (default: the finest level)')
    p = sub.add_parser('repeatability', parents=[common])
    p.add_argument('--group', action='append', required=True, metavar='NAME=R1,R2[,R3]')
    p = sub.add_parser('consistency', parents=[common])
    p.add_argument('--runs', nargs='+', required=True)
    p.add_argument('--labels', nargs='+')
    args = ap.parse_args(argv)

    if args.workers:
        compare_module.WORKERS = max(1, args.workers)
    kw = dict(grid_mm=args.grid_mm, max_tracks=args.max_tracks, max_pairs=args.max_pairs,
              sensitivity=args.sensitivity)
    if args.study == 'compare':
        ts = [_load(args.a, args.image), _load(args.b, args.image)]
        result, rows = studies.consistency(ts, **kw)
        result['study'] = 'consistency'
    elif args.study == 'convergence':
        ts = [_load(r, args.image) for r in args.runs]
        values = args.values or [_override_value(t, args.param) for t in ts]
        if any(v is None for v in values):
            ap.error('could not read --param values from overrides.txt; pass --values')
        values = [float(v) if _numeric(v) else v for v in values]
        if all(isinstance(v, float) for v in values):
            finer = args.finer or ('larger' if 'density' in args.param else 'smaller')
            order = sorted(range(len(ts)), key=lambda k: values[k], reverse=(finer == 'smaller'))
            ts, values = [ts[k] for k in order], [values[k] for k in order]
        ref = _load(args.reference, args.image) if args.reference else None
        result, rows = studies.convergence(ts, values, args.param, reference=ref, **kw)
        if ref:
            ts = ts + [ref]
    elif args.study == 'repeatability':
        groups = {}
        for g in args.group:
            name, _, members = g.partition('=')
            groups[name] = [_load(m, args.image) for m in members.split(',') if m]
            if len(groups[name]) < 2:
                ap.error(f'group {name} needs at least two runs')
        result, rows = studies.repeatability(groups, **kw)
        ts = [t for v in groups.values() for t in v]
    else:
        labels = args.labels or [None] * len(args.runs)
        ts = [_load(r, args.image, l) for r, l in zip(args.runs, labels)]
        result, rows = studies.consistency(ts, **kw)

    out = args.out or REPO / 'hinec_runs' / (
        f"assess_{datetime.datetime.now():%Y%m%d_%H%M%S}_{result['study']}" + (f'_{args.name}' if args.name else ''))
    protocols = [PRIMARY.as_dict()] + ([s.as_dict() for s in SENSITIVITY] if args.sensitivity else [])
    report.write(out, result, rows, _inputs(ts), protocols)
    print(f'[assess] wrote {out}/report.html')
    for h in result['headline']:
        print('  ' + ' | '.join(f'{k}={_fmt(v)}' for k, v in h.items()
                                if k in ('level', 'group', 'occupied_jaccard', 'path_separation_median_mm',
                                         'kappa_w_rel_median_change', 'kappa_w_paired_spearman')))
    return out


def _numeric(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def _fmt(v):
    return f'{v:.4g}' if isinstance(v, float) else str(v)


if __name__ == '__main__':
    main()
