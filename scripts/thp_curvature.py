#!/usr/bin/env python3
"""Traveling Human Phantom (ds000206): is tract curvature a reproducible property of a person?

Protocol (frozen before outcomes): report/KSIAM/20261026/thp_protocol.md. Inputs are the run directories of
    ./bin/run_tractography.sh ds000206_systems --source data/ds000206/nim/<scan>.mat
Each streamline is assigned to a unit by the JHU label at its seed point. Curvature uses the frozen protocol
(tract_assessment.geometry.connection: 0.5 mm resampling, cubic Savitzky-Golay over 4 mm, ends trimmed), in world mm.

Per scan and unit:
  Y      median over streamlines (>= 20 mm) of the per-streamline median kappa   (primary)
  Y90    median over streamlines of the per-streamline 90th percentile of kappa
  Phi    median net turn (degrees; first vs last 10 mm)
  FA     mean FA in the unit's seed label (positive control)
  count  number of qualifying streamlines (yield control)
Per unit and outcome: ICC(1) (one-way random, unbalanced k0, F-based 95% CI) and a nested variance decomposition
person / session within person / run within session.

Usage: python scripts/thp_curvature.py --out report/KSIAM/20261026/thp            (spline, primary)
       python scripts/thp_curvature.py --out report/KSIAM/20261026/thp_cubic --runs 'hinec_runs/run_*_ds000206_systems_methodcubic'
"""
import argparse
import csv
import glob
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from tract_assessment.io import load_tractogram        # noqa: E402
from tract_assessment.geometry import connection        # noqa: E402

UNITS = {'CC': (3, 4, 5), 'fornix': (6,), 'CST_L': (8,), 'CST_R': (7,), 'SS_L': (32,), 'SS_R': (31,),
         'UF_L': (48,), 'UF_R': (47,)}
OUTCOMES = ('Y', 'Y90', 'Phi', 'FA', 'count', 'L')
CENTRAL_MM = None   # --central-mm: measure only the central stretch of each streamline (length-matched check)
MIN_LEN_MM, MIN_STREAMLINES, END_MM = 20.0, 20, 10.0
SCAN_RE = re.compile(r'(thp\d{4})_([a-z]+)_r(\d{2})')


def source_nim(run_dir):
    txt = (Path(run_dir) / 'intermediate' / 'SOURCE.txt').read_text()
    m = re.search(r'(\S*data/ds000206/nim/\S+\.mat)', txt)
    if not m: raise ValueError(f'no ds000206 nim named in {run_dir}/intermediate/SOURCE.txt')
    return ROOT / m.group(1) if not m.group(1).startswith('/') else Path(m.group(1))


def net_turn(p, ds=0.5):
    n = int(END_MM / ds)
    if len(p) < 3 * n: return np.nan
    a, b = p[n] - p[0], p[-1] - p[-1 - n]
    return float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))


def scan_rows(run_dir):
    nim_path = source_nim(run_dir); tag = SCAN_RE.search(nim_path.name)
    person, site, run = tag.group(1), tag.group(2), tag.group(3)
    m = sio.loadmat(str(nim_path), variable_names=['nim'], squeeze_me=True, struct_as_record=False)['nim']
    lab, fa = np.asarray(m.parcellation_mask), np.asarray(m.FA, float)
    tg = load_tractogram(str(run_dir))
    if Path(tg.reference).name != nim_path.stem + '.nii.gz':   # a borrowed reference = wrong world coordinates
        raise ValueError(f'{run_dir}: DWI reference {Path(tg.reference).name} is not this scan\'s own')
    if tg.seeds is None: raise ValueError(f'{run_dir}: seed identity not recoverable')
    sv = np.clip(np.rint(tg.seeds).astype(int) - 1, 0, np.array(lab.shape) - 1)
    seed_label = lab[tuple(sv.T)]
    rows = []
    for unit, ids in UNITS.items():
        inlab = np.isin(lab, ids)
        sel = np.flatnonzero(np.isin(seed_label, ids))
        med, p90, phi, lens = [], [], [], []
        for i in sel:
            w = np.asarray(tg.world[i], float)
            seg = np.linalg.norm(np.diff(w, axis=0), axis=1); L = seg.sum()
            if L < MIN_LEN_MM: continue
            if CENTRAL_MM:
                if L < CENTRAL_MM: continue
                s_ = np.concatenate([[0], np.cumsum(seg)]); a0 = (L - CENTRAL_MM) / 2
                keep = (s_ >= a0) & (s_ <= a0 + CENTRAL_MM)
                if keep.sum() < 4: continue
                w = w[keep]
            lens.append(L)
            g = connection(np.asarray(w, float), diagnostics=False)
            if g is None: continue
            k = g['kappa'][np.isfinite(g['kappa'])]
            if not len(k): continue
            med.append(np.median(k)); p90.append(np.percentile(k, 90)); phi.append(net_turn(g['p']))
        rows.append(dict(person=person, site=site, run=run, unit=unit, run_dir=str(run_dir), label_voxels=int(inlab.sum()),
                         count=len(med), Y=float(np.median(med)) if med else np.nan, L=float(np.median(lens)) if lens else np.nan,
                         Y90=float(np.median(p90)) if p90 else np.nan, Phi=float(np.nanmedian(phi)) if phi else np.nan,
                         FA=float(fa[inlab].mean()) if inlab.any() else np.nan))
    return rows


def exclusions(rows):
    """Pre-set rules: < 20 streamlines; label < 50% of that person's median label size for the unit."""
    med = defaultdict(list)
    for r in rows: med[(r['person'], r['unit'])].append(r['label_voxels'])
    for r in rows:
        r['excluded'] = ('few_streamlines' if r['count'] < MIN_STREAMLINES else
                         'truncated_label' if r['label_voxels'] < 0.5 * np.median(med[(r['person'], r['unit'])]) else '')
    return rows


def icc1(groups, alpha=0.05):
    """One-way random ICC(1) for unbalanced groups; F-based CI with k0 (Searle 1971)."""
    groups = [np.asarray(g, float) for g in groups if len(g) >= 2]
    a = len(groups); n = np.array([len(g) for g in groups]); N = n.sum()
    if a < 3: return dict(icc=np.nan, lo=np.nan, hi=np.nan, a=a, N=int(N))
    grand = np.concatenate(groups).mean()
    msb = sum(len(g) * (g.mean() - grand) ** 2 for g in groups) / (a - 1)
    msw = sum(((g - g.mean()) ** 2).sum() for g in groups) / (N - a)
    k0 = (N - (n ** 2).sum() / N) / (a - 1)
    F = msb / msw if msw > 0 else np.inf
    fl, fu = F / stats.f.ppf(1 - alpha / 2, a - 1, N - a), F / stats.f.ppf(alpha / 2, a - 1, N - a)
    return dict(icc=float((msb - msw) / (msb + (k0 - 1) * msw)), lo=float((fl - 1) / (fl + k0 - 1)),
                hi=float((fu - 1) / (fu + k0 - 1)), a=a, N=int(N), k0=float(k0), F=float(F))


def nested_sd(rows, key):
    """SDs from a nested ANOVA (method of moments, unbalanced handled by per-level means):
    run within session (repeatability), session within person (cross-site), person."""
    by_sess = defaultdict(list)
    for r in rows: by_sess[(r['person'], r['site'])].append(r[key])
    within = [np.var(v, ddof=1) for v in by_sess.values() if len(v) == 2]
    s_run = np.sqrt(np.mean(within)) if within else np.nan
    sess_means = defaultdict(list)
    for (p, s), v in by_sess.items(): sess_means[p].append(np.mean(v))
    v_sess = np.mean([np.var(v, ddof=1) for v in sess_means.values() if len(v) >= 2]) - (s_run ** 2) / 2
    person_means = [np.mean(v) for v in sess_means.values()]
    nbar = np.mean([len(v) for v in sess_means.values()])
    v_person = np.var(person_means, ddof=1) - max(v_sess, 0) / nbar - (s_run ** 2) / (2 * nbar)
    sd = lambda v: float(np.sqrt(v)) if v > 0 else 0.0
    return dict(sd_run=float(s_run), sd_session=sd(v_sess), sd_person=sd(v_person), grand_mean=float(np.mean([r[key] for r in rows])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True); ap.add_argument('--runs', default='hinec_runs/run_*_ds000206_systems',
                    help="glob of run dirs; default = spline runs only. Cubic check: 'hinec_runs/run_*_ds000206_systems_methodcubic'")
    ap.add_argument('--central-mm', type=float, default=None,
                    help='length-matched check: keep streamlines of at least this length and measure only their central stretch of it')
    a = ap.parse_args(); out = ROOT / a.out; out.mkdir(parents=True, exist_ok=True)
    global CENTRAL_MM; CENTRAL_MM = a.central_mm
    runs = sorted(glob.glob(str(ROOT / a.runs)))
    latest = {}                                   # one run per scan: the newest
    for d in runs:
        try: latest[source_nim(d).name] = d
        except (OSError, ValueError) as e: print('skip', d, e)
    rows = []
    for name, d in sorted(latest.items()):
        if not list(Path(d, 'tractography').glob('tracks*.mat')): print('skip (no tracks)', d); continue
        try: rows += scan_rows(d); print('done', name)
        except ValueError as e: print('skip', e)
    rows = exclusions(rows)
    with open(out / 'per_scan.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary = {}
    for unit in UNITS:
        R = [r for r in rows if r['unit'] == unit and not r['excluded'] and np.isfinite(r['Y'])]
        summary[unit] = dict(scans=len(R), excluded={k: sum(1 for r in rows if r['unit'] == unit and r['excluded'] == k)
                                                      for k in ('few_streamlines', 'truncated_label')})
        for key in OUTCOMES:
            g = defaultdict(list)
            for r in R: g[r['person']].append(r[key])
            summary[unit][key] = dict(icc1=icc1(list(g.values())), **(nested_sd(R, key) if R else {}))
    (out / 'summary.json').write_text(json.dumps(dict(protocol='report/KSIAM/20261026/thp_protocol.md', n_scans=len(latest),
                                                      units=summary), indent=1))
    print(f"{'unit':7s} {'scans':>5s}  " + '  '.join(f'{k:>20s}' for k in OUTCOMES))
    for unit, s in summary.items():
        print(f"{unit:7s} {s['scans']:5d}  " + '  '.join(
            f"{s[k]['icc1']['icc']:5.2f} [{s[k]['icc1']['lo']:5.2f},{s[k]['icc1']['hi']:5.2f}]" for k in OUTCOMES))


if __name__ == '__main__':
    main()
