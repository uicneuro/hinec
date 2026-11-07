#!/usr/bin/env python3
"""Check the lifts of field: lift (nim_lift_field) against the ISMRM 2015 ground truth, voxel by voxel.

Ground truth: the Renauld-2023 bundles (data/ismrm2015/scoring_data_Renauld2023/bundles/*.trk), mapped to
1-based voxel coordinates by the inverse of the DWI affine. In every voxel x and for every bundle b with
streamline points in x, d_b(x) = the principal eigenvector of sum t t' over b's unit tangents t in x; the
d_b(x) are clustered greedily into directions more than 30 deg apart, n_dir(x) = their number (capped at 3).

For every voxel with n_dir >= 1 and every ground-truth direction d there:
    v1 error   = arccos |v1(x) . d|              (nim.evec(:,:,:,:,1), the DTI principal eigenvector)
    lift error = min over the lifts u of x of arccos |u . d|  (90 when x has no lift)
and the number of lifts against n_dir. Both fields are in voxel-axis coordinates.

Usage: python scripts/lift_field_check.py <lift cache .mat> [nim .mat] [--z z0 z1]
       (--z: only voxels with 1-based slice index z0..z1, for lifts built on a slab)
       (the cache is data/ismrm2015/ismrm2015_lift_<hash>.mat, written by nim_field; the nim defaults to
        data/ismrm2015/ismrm2015.mat)
"""
import sys, glob
import numpy as np, pandas as pd, nibabel as nib, h5py

ROOT = '/'.join(__file__.split('/')[:-2]) if '/' in __file__ else '..'
GT = ROOT + '/data/ismrm2015/scoring_data_Renauld2023/bundles/'


def load_h5(path, names):
    with h5py.File(path, 'r') as h:
        return {n: np.array(h[n]).T for n in names if n in h}


def gt_directions(aff, shape):
    inv = np.linalg.inv(aff); rows = []
    for f in sorted(glob.glob(GT + '*.trk')):
        name = f.split('/')[-1][:-4]
        for sl in nib.streamlines.load(f).streamlines:
            v = (inv[:3, :3] @ np.asarray(sl).T).T + inv[:3, 3]           # 0-based voxel coordinates
            if len(v) < 3: continue
            t = np.gradient(v, axis=0); t /= np.linalg.norm(t, axis=1, keepdims=True) + 1e-12
            ijk = np.rint(v).astype(int); ok = np.all((ijk >= 0) & (ijk < shape), 1)
            rows.append(pd.DataFrame({'v': np.ravel_multi_index(ijk[ok].T, shape), 'b': name,
                                      'xx': t[ok, 0] ** 2, 'yy': t[ok, 1] ** 2, 'zz': t[ok, 2] ** 2,
                                      'xy': t[ok, 0] * t[ok, 1], 'xz': t[ok, 0] * t[ok, 2], 'yz': t[ok, 1] * t[ok, 2]}))
    g = pd.concat(rows).groupby(['v', 'b']).sum().reset_index()
    M = np.stack([np.stack([g.xx, g.xy, g.xz], -1), np.stack([g.xy, g.yy, g.yz], -1), np.stack([g.xz, g.yz, g.zz], -1)], -1)
    d = np.linalg.eigh(M)[1][:, :, -1]
    out = {}
    for v, idx in g.groupby('v').indices.items():
        reps = []
        for u in d[idx]:
            if all(abs(u @ r) < np.cos(np.radians(30)) for r in reps): reps.append(u)
        out[v] = reps
    return out


def main(lift_mat, nim_mat, zr=None):
    L = load_h5(lift_mat, ['U', 'ok', 'mass'])
    U, ok = L['U'], L['ok'].astype(bool)
    try:                                                    # a -v7.3 nim
        with h5py.File(nim_mat, 'r') as h:
            nim = h['nim'] if 'nim' in h else h
            evec = np.array(nim['evec']).T
    except OSError:                                         # a -v7 nim
        import scipy.io as sio
        evec = sio.loadmat(nim_mat, variable_names=['nim'], struct_as_record=False, squeeze_me=True)['nim'].evec
    shape = ok.shape[:3]
    v1 = evec[..., :, 0]
    aff = nib.load(ROOT + '/data/ismrm2015/ismrm2015.nii.gz').affine
    dirs = gt_directions(aff, shape)
    rows = []
    for v, reps in dirs.items():
        i = np.unravel_index(v, shape)
        if zr is not None and not (zr[0] - 1 <= i[2] <= zr[1] - 1): continue
        if not ok[i].any() and not np.any(v1[i]): continue
        nd = min(len(reps), 3); nl = int(ok[i].sum())
        for r in reps:
            e_v1 = np.degrees(np.arccos(min(1.0, abs(float(v1[i] @ r)))))
            e_l = np.degrees(np.arccos(min(1.0, float(np.abs(U[i][ok[i]] @ r).max())))) if nl else 90.0
            rows.append(dict(ndir=nd, nlift=nl, v1=e_v1, lift=e_l))
    R = pd.DataFrame(rows)
    print('(voxel, ground-truth direction) pairs by n_dir: median error and share > 20 deg; lifts per voxel')
    for nd, q in R.groupby('ndir'):
        print(f'  n_dir {nd}: {len(q):6d} pairs | v1 {q.v1.median():5.1f} deg ({100 * (q.v1 > 20).mean():4.1f}% > 20) | '
              f'lifts {q.lift.median():5.1f} deg ({100 * (q.lift > 20).mean():4.1f}% > 20) | '
              f'lifts per voxel: ' + ' '.join(f'{k}:{100 * (q.nlift == k).mean():.0f}%' for k in range(4)))


if __name__ == '__main__':
    a = sys.argv[1:]; zr = None
    if '--z' in a:
        k = a.index('--z'); zr = (int(a[k + 1]), int(a[k + 2])); a = a[:k] + a[k + 3:]
    main(a[0], a[1] if len(a) > 1 else ROOT + '/data/ismrm2015/ismrm2015.mat', zr)
