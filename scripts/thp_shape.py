#!/usr/bin/env python3
"""Traveling Human Phantom (ds000206): is the 3-D SHAPE of a bundle reproducible across sites, per tracker?

Protocol (frozen before outcomes): report/KSIAM/20261026/thp_protocol.md, addendum (b). Inputs are the run
directories of  ./bin/run_tractography.sh ds000206_systems --source data/ds000206/nim/<scan>.mat [--set algorithm=mmf].
Streamlines are assigned to a unit by the JHU label at the seed, exactly as in scripts/thp_curvature.py.

Bundle B(scan, unit): streamlines >= 20 mm, world mm, resampled to K = 20 points, at most M = 200 (fixed stride).
  MDF(a, b) = min(mean_k |a_k - b_k|, mean_k |a_k - b_{K+1-k}|)                        (mm)
  BMD(A, B) = 1/2 [mean_a min_b MDF + mean_b min_a MDF]                                  (mm)
  d(A, B)   = BMD after rigid streamline-based linear registration (SLR, dipy), averaged over both directions.
Per unit and arm: W = median d over within-person cross-site pairs (run 01 of each session); R = median
between-person d / W. Decision rule in the protocol.

Usage: python scripts/thp_shape.py --out report/KSIAM/20261026/thp_shape [--workers 16]
       addendum (c): --primary rel --streamline 'hinec_runs/run_*_ds000206_systems_angle_max45'
                     --mmf 'hinec_runs/run_*_ds000206_systems_algorithmmmf_angle_max45' --out report/KSIAM/20261026/thp_shape_angle45
"""
import argparse
import json
import sys
from collections import defaultdict
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import glob
import numpy as np
import scipy.io as sio
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from tract_assessment.io import load_tractogram        # noqa: E402
from thp_curvature import UNITS, MIN_LEN_MM, MIN_STREAMLINES, SCAN_RE, source_nim   # noqa: E402

K, M = 20, 200
ARMS = {'streamline': 'hinec_runs/run_*_ds000206_systems', 'mmf': 'hinec_runs/run_*_ds000206_systems_algorithmmmf'}
PRIMARY = 'd'   # 'd' (addendum b) or 'rel' = d / mean streamline length of the two bundles (addendum c)


def resample(p, k=K):
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s[-1], k)
    return np.stack([np.interp(t, s, p[:, j]) for j in range(3)], axis=1)


def bundles_of(run_dir):
    """{unit: (n, K, 3) array} for one scan, or None for a unit below the streamline floor."""
    nim_path = source_nim(run_dir)
    lab = np.asarray(sio.loadmat(str(nim_path), variable_names=['nim'], squeeze_me=True,
                                 struct_as_record=False)['nim'].parcellation_mask)
    tg = load_tractogram(str(run_dir))
    if Path(tg.reference).name != nim_path.stem + '.nii.gz':
        raise ValueError(f'{run_dir}: DWI reference is not this scan\'s own')
    sv = np.clip(np.rint(tg.seeds).astype(int) - 1, 0, np.array(lab.shape) - 1)
    seed_label = lab[tuple(sv.T)]
    out = {}
    for unit, ids in UNITS.items():
        keep = []
        for i in np.flatnonzero(np.isin(seed_label, ids)):
            w = np.asarray(tg.world[i], float)
            if len(w) >= 2 and np.linalg.norm(np.diff(w, axis=0), axis=1).sum() >= MIN_LEN_MM:
                keep.append(resample(w))
        if len(keep) < MIN_STREAMLINES:
            out[unit] = None; continue
        stride = max(1, int(np.ceil(len(keep) / M)))
        sub = keep[::stride][:M]
        out[unit] = (np.stack(sub), float(np.mean([np.linalg.norm(np.diff(b, axis=0), axis=1).sum() for b in sub])))
    return out


def bmd(A, B):
    direct = np.linalg.norm(A[:, None] - B[None], axis=3).mean(axis=2)
    flip = np.linalg.norm(A[:, None] - B[None, :, ::-1], axis=3).mean(axis=2)
    D = np.minimum(direct, flip)
    return 0.5 * (D.min(axis=1).mean() + D.min(axis=0).mean())


def shape_distance(pair):
    from dipy.align.streamlinear import StreamlineLinearRegistration
    A, B = pair
    ds = []
    for S, Mv in ((A, B), (B, A)):
        srm = StreamlineLinearRegistration(x0='rigid').optimize(static=list(S), moving=list(Mv))
        ds.append(bmd(S, np.stack(srm.transform(list(Mv)))))
    return float(np.mean(ds)), float(bmd(A, B))


def sessions(arm_glob):
    """{(person, site): run_dir} for run 01 of each session; the newest run dir per scan."""
    latest = {}
    for d in sorted(glob.glob(str(ROOT / arm_glob))):
        try: latest[source_nim(d).name] = d
        except (ValueError, FileNotFoundError): continue
    out = {}
    for name, d in latest.items():
        m = SCAN_RE.search(name)
        if m.group(3) == '01' and list(Path(d, 'tractography').glob('tracks*.mat')):
            out[(m.group(1), m.group(2))] = d
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--streamline', default=ARMS['streamline'], help='run-dir glob of the streamline arm')
    ap.add_argument('--mmf', default=ARMS['mmf'], help='run-dir glob of the MMF arm')
    ap.add_argument('--primary', choices=['d', 'rel'], default='d', help='rel = length-normalised d (addendum c)')
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    ARMS['streamline'], ARMS['mmf'] = a.streamline, a.mmf
    sess = {arm: sessions(g) for arm, g in ARMS.items()}
    keys = sorted(set(sess['streamline']) & set(sess['mmf']))
    print(f'sessions: streamline {len(sess["streamline"])}, mmf {len(sess["mmf"])}, both {len(keys)}', flush=True)
    with Pool(min(a.workers, 8)) as pool:
        B = {arm: dict(zip(keys, pool.map(bundles_of, [sess[arm][k] for k in keys]))) for arm in ARMS}

    rows, summary = [], {}
    with Pool(a.workers) as pool:
        for unit in UNITS:
            ok = [k for k in keys if B['streamline'][k][unit] is not None and B['mmf'][k][unit] is not None]
            pairs = [(p, q) for p, q in combinations(ok, 2) if p[1] != q[1]]          # different sites only
            res = {}
            for arm in ARMS:
                r = pool.map(shape_distance, [(B[arm][p][unit][0], B[arm][q][unit][0]) for p, q in pairs])
                L = [0.5 * (B[arm][p][unit][1] + B[arm][q][unit][1]) for p, q in pairs]
                res[arm] = [(x[0] / l if a.primary == 'rel' else x[0], x[1], x[0], l) for x, l in zip(r, L)]
            for i, (p, q) in enumerate(pairs):
                rows.append(dict(unit=unit, a=f'{p[0]}_{p[1]}', b=f'{q[0]}_{q[1]}', within=p[0] == q[0],
                                 **{f'd_{arm}': res[arm][i][0] for arm in ARMS},
                                 **{f'raw_{arm}': res[arm][i][1] for arm in ARMS},
                                 **{f'dmm_{arm}': res[arm][i][2] for arm in ARMS},
                                 **{f'len_{arm}': res[arm][i][3] for arm in ARMS}))
            u = {}
            within = [i for i, (p, q) in enumerate(pairs) if p[0] == q[0]]
            between = [i for i, (p, q) in enumerate(pairs) if p[0] != q[0]]
            for arm in ARMS:
                d = np.array([r[0] for r in res[arm]])
                W, Bt = float(np.median(d[within])), float(np.median(d[between]))
                u[arm] = dict(W=W, between=Bt, R=Bt / W, control_ok=Bt > W, n_within=len(within), n_between=len(between))
            dm = np.array([res['mmf'][i][0] for i in within]); dh = np.array([res['streamline'][i][0] for i in within])
            w = stats.wilcoxon(dm, dh)
            u['paired'] = dict(median_diff=float(np.median(dm - dh)), mmf_lower=int((dm < dh).sum()),
                               mmf_higher=int((dm > dh).sum()), p=float(w.pvalue))
            u['better'] = bool(u['mmf']['W'] < u['streamline']['W'] and w.pvalue < 0.05 and np.median(dm - dh) < 0)
            u['worse'] = bool(u['mmf']['W'] > u['streamline']['W'])
            summary[unit] = u
            print(f"{unit:7} W stream {u['streamline']['W']:.2f}  W mmf {u['mmf']['W']:.2f} mm | between "
                  f"{u['streamline']['between']:.2f}/{u['mmf']['between']:.2f} | paired diff {u['paired']['median_diff']:+.3f} "
                  f"({u['paired']['mmf_lower']}/{u['paired']['mmf_higher']}) p={w.pvalue:.2g} | better={u['better']}", flush=True)

    nb, nw = sum(u['better'] for u in summary.values()), sum(u['worse'] for u in summary.values())
    verdict = 'supported' if nb >= 6 else 'falsified' if nw >= 4 else 'no difference'
    (out / 'summary.json').write_text(json.dumps(dict(protocol='report/KSIAM/20261026/thp_protocol.md (b)/(c)', primary=a.primary,
                                                      arms=dict(ARMS),
                                                      sessions=len(keys), units=summary, n_better=nb, n_worse=nw,
                                                      verdict=verdict), indent=1))
    import csv
    with open(out / 'pairs.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
    print(f'better in {nb}/8, worse in {nw}/8 -> {verdict}')


if __name__ == '__main__':
    main()
