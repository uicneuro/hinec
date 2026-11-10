"""DTI stage of HINEC, ported from src/nim_utils/nim_read.m, src/nim_calculation/{nim_dt_spd,nim_eig,nim_fa}.m
and src/nim_preprocessing/{preproc_mask_improvement,preproc_tissue_segmentation}.m.

Conventions kept identical to the MATLAB code:
  * b0 volumes are bval < B0Threshold (default 5); img_b0 is their mean.
  * Y_k = log(S0/S_k)/b_k, H = [gx^2 gy^2 gz^2 2gxgy 2gygz 2gzgx], D = H \\ Y (unweighted LSF).
  * D is stored as the 6-vector [Dxx Dyy Dzz Dxy Dyz Dxz] (nim_reshape_d ordering).
  * eigenvalues sorted descending; evec[..., :, 0] is the primary eigenvector.
  * Tensors with NaN/Inf, negative eigenvalues (after the SPD refit) or eigenvalues > 0.01 are zeroed.
  * FA computed only inside the brain mask, zero elsewhere.
  * b-vectors are used exactly as given in the .bvec file (voxel-index frame), as HINEC does.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import minimize


# --------------------------------------------------------------------------------------
# nim_read
# --------------------------------------------------------------------------------------
@dataclass
class Nim:
    img: np.ndarray                  # (X,Y,Z,N) raw DWI
    affine: np.ndarray
    zooms: tuple
    bval: np.ndarray                 # (N,)
    bvec: np.ndarray                 # (N,3)
    thrsh_b0: float
    img_b0: np.ndarray               # (X,Y,Z) mean b0
    img_bi: np.ndarray               # (X,Y,Z,K) diffusion-weighted volumes
    mask: np.ndarray                 # (X,Y,Z) uint8
    DT: np.ndarray | None = None     # (X,Y,Z,6)
    evec: np.ndarray | None = None   # (X,Y,Z,3,3)
    eval: np.ndarray | None = None   # (X,Y,Z,3)
    FA: np.ndarray | None = None     # (X,Y,Z)
    wm_mask: np.ndarray | None = None
    gm_mask: np.ndarray | None = None
    csf_mask: np.ndarray | None = None
    info: dict = field(default_factory=dict)

    @property
    def dims(self):
        return self.img.shape[:3]


def _read_bval(path):
    lines = [ln for ln in open(path).read().splitlines() if ln.strip()]
    if len(lines) == 1:
        return np.array([float(v) for v in lines[0].split()])
    return np.array([float(ln.strip()) for ln in lines])


def _read_bvec(path):
    lines = [ln for ln in open(path).read().splitlines() if ln.strip()]
    first = [float(v) for v in lines[0].split()]
    if len(first) == 3 and len(lines) > 3:            # IronTract format: one vector per line
        return np.array([[float(v) for v in ln.split()] for ln in lines])
    gx, gy, gz = (np.array([float(v) for v in lines[i].split()]) for i in range(3))
    return np.stack([gx, gy, gz], axis=1)


def nim_read(prefix: str, mask_path: str | None = None, b0_threshold: float = 5.0) -> Nim:
    """Mirror of nim_read.m. `prefix` is the path without extension (expects prefix.nii.gz, .bval, .bvec).
    If mask_path is None and prefix_M.nii.gz is missing, mask is left empty (caller must supply one)."""
    img = nib.load(prefix + ".nii.gz")
    data = np.asanyarray(img.dataobj)
    bval = _read_bval(prefix + ".bval")
    bvec = _read_bvec(prefix + ".bvec")
    is_b0 = bval < b0_threshold
    img_b0 = data[..., is_b0].astype(np.float64).mean(axis=3)
    img_bi = data[..., ~is_b0].astype(np.float64)
    mask_file = mask_path or (prefix + "_M.nii.gz")
    try:
        mask = np.asanyarray(nib.load(mask_file).dataobj).astype(np.uint8)
    except FileNotFoundError:
        mask = np.zeros(data.shape[:3], np.uint8)
    return Nim(img=data, affine=img.affine, zooms=img.header.get_zooms()[:3], bval=bval, bvec=bvec,
               thrsh_b0=b0_threshold, img_b0=img_b0, img_bi=img_bi, mask=mask,
               info=dict(prefix=prefix, n_b0=int(is_b0.sum()), n_dwi=int((~is_b0).sum())))


# --------------------------------------------------------------------------------------
# Brain extraction (stand-in for FSL BET -f 0.3 used by preproc_brain_extraction.m)
# --------------------------------------------------------------------------------------
def brain_mask_from_b0(img_b0: np.ndarray, median_radius: int = 2, numpass: int = 4) -> np.ndarray:
    """dipy median_otsu on the mean b0, keep the largest component, fill holes.
    FSL BET is not available in this environment; this is the closest drop-in."""
    from dipy.segment.mask import median_otsu
    _, mask = median_otsu(img_b0, median_radius=median_radius, numpass=numpass)
    lab, n = ndi.label(mask)
    if n > 1:
        sizes = ndi.sum(mask, lab, index=np.arange(1, n + 1))
        mask = lab == (1 + int(np.argmax(sizes)))
    mask = ndi.binary_fill_holes(mask)
    return mask.astype(np.uint8)


# --------------------------------------------------------------------------------------
# nim_dt_spd
# --------------------------------------------------------------------------------------
def _tensor_from_vec(d: np.ndarray) -> np.ndarray:
    """nim_reshape_d for an (N,6) array -> (N,3,3)."""
    D = np.empty(d.shape[:-1] + (3, 3), d.dtype)
    D[..., 0, 0] = d[..., 0]; D[..., 1, 1] = d[..., 1]; D[..., 2, 2] = d[..., 2]
    D[..., 0, 1] = D[..., 1, 0] = d[..., 3]
    D[..., 1, 2] = D[..., 2, 1] = d[..., 4]
    D[..., 0, 2] = D[..., 2, 0] = d[..., 5]
    return D


def _vec_from_tensor(D: np.ndarray) -> np.ndarray:
    return np.stack([D[..., 0, 0], D[..., 1, 1], D[..., 2, 2], D[..., 0, 1], D[..., 1, 2], D[..., 0, 2]], axis=-1)


def _eig_desc(D: np.ndarray):
    """eigh + descending sort (stable, like MATLAB maxk on ties). Returns (evals (N,3), evecs (N,3,3) columns)."""
    w, v = np.linalg.eigh(D)
    order = np.argsort(-w, axis=-1, kind="stable")
    w_s = np.take_along_axis(w, order, axis=-1)
    v_s = np.take_along_axis(v, order[..., None, :], axis=-1)
    return w_s, v_s


def _spd_refit(S: np.ndarray, S0: float, b: np.ndarray, g: np.ndarray, D_init: np.ndarray) -> np.ndarray:
    """SPD-constrained tensor fit for one voxel: minimise sum_k (log(S_k/S0) + b_k g_k^T D g_k)^2 with D = L L^T.

    This is the objective of lib/bfgs/vox_dt_bfgs.m (objf). The MATLAB gradient (gradf) indexes g(1),g(2),g(3)
    — the x-components of the first three gradients — instead of the per-direction components, so its BFGS
    iterations do not descend on objf. We use the correct analytic gradient with L-BFGS-B instead.
    Units follow the MATLAB file: b in 10^3 s/mm^2 during the fit, D returned in mm^2/s.
    """
    b = b / 1000.0
    y = np.log(S / S0)                      # (K,)
    GG = np.einsum("ki,kj->kij", g, g)      # (K,3,3)

    # Cholesky-parameterised start: project the LSF tensor onto SPD.
    w, v = np.linalg.eigh(D_init * 1000.0)
    w = np.clip(w, 1e-4, None)
    L0 = np.linalg.cholesky((v * w) @ v.T)
    x0 = L0[np.tril_indices(3)]

    def unpack(x):
        L = np.zeros((3, 3)); L[np.tril_indices(3)] = x
        return L

    def fg(x):
        L = unpack(x); D = L @ L.T
        q = np.einsum("kij,ij->k", GG, D)   # g^T D g
        r = y + b * q
        f = np.sum(r * r)
        dD = np.einsum("k,kij->ij", 2.0 * r * b, GG)     # df/dD (symmetric)
        dL = (dD + dD.T) @ L                              # df/dL
        return f, dL[np.tril_indices(3)]

    res = minimize(fg, x0, jac=True, method="L-BFGS-B", options=dict(maxiter=500))
    L = unpack(res.x)
    return (L @ L.T) / 1000.0


def nim_dt_spd(nim: Nim, verbose: bool = True) -> Nim:
    t0 = time.time()
    b = nim.bval[nim.bval >= nim.thrsh_b0]
    g = nim.bvec[nim.bval >= nim.thrsh_b0, :]
    gx, gy, gz = g[:, 0], g[:, 1], g[:, 2]
    H = np.stack([gx ** 2, gy ** 2, gz ** 2, 2 * gx * gy, 2 * gy * gz, 2 * gz * gx], axis=1)   # (K,6)

    X, Y, Z = nim.dims
    K = nim.img_bi.shape[3]
    DT = np.zeros((X, Y, Z, 6))
    evec = np.zeros((X, Y, Z, 3, 3))
    eval_ = np.zeros((X, Y, Z, 3))

    fit_vox = (nim.mask > 0) & (nim.img_b0 != 0)
    idx = np.flatnonzero(fit_vox.ravel())
    S0 = nim.img_b0.reshape(-1)[idx]                              # (N,)
    S = nim.img_bi.reshape(-1, K)[idx]                            # (N,K)
    with np.errstate(divide="ignore", invalid="ignore"):
        Yv = np.log(S0[:, None] / S) / b[None, :]                 # (N,K)
    Yv = np.nan_to_num(Yv, nan=np.nan, posinf=np.inf, neginf=-np.inf)

    # Unweighted least squares for all voxels at once (H \ Y in MATLAB)
    finite = np.all(np.isfinite(Yv), axis=1)
    Dv = np.full((len(idx), 6), np.nan)
    Dv[finite] = np.linalg.lstsq(H, Yv[finite].T, rcond=None)[0].T

    # SPD refit where any eigenvalue is negative (vox_dt_bfgs path)
    Dm = _tensor_from_vec(np.nan_to_num(Dv))
    w, _ = _eig_desc(Dm)
    need = finite & np.any(w < 0, axis=1)
    n_bfgs = int(need.sum())
    if verbose:
        print(f"  LSF tensors: {len(idx)} voxels, {n_bfgs} need SPD refit ({100 * n_bfgs / max(1, len(idx)):.1f}%)")
    for n in np.flatnonzero(need):
        try:
            Dv[n] = _vec_from_tensor(_spd_refit(S[n], S0[n], b, g, Dm[n]))
        except Exception:
            Dv[n] = np.nan

    # Eigendecomposition, validation (mirrors the checks in nim_dt_spd.m)
    bad = ~np.all(np.isfinite(Dv), axis=1)
    Dv[bad] = 0.0
    Dm = _tensor_from_vec(Dv)
    w, v = _eig_desc(Dm)
    bad |= np.any(w < 0, axis=1) | np.any(w > 0.01, axis=1)
    Dv[bad] = 0.0; w[bad] = 0.0; v[bad] = 0.0

    DT.reshape(-1, 6)[idx] = Dv
    evec.reshape(-1, 3, 3)[idx] = v
    eval_.reshape(-1, 3)[idx] = w
    nim.DT, nim.evec, nim.eval = DT, evec, eval_
    nim.info.update(n_fit_voxels=int(len(idx)), n_spd_refit=n_bfgs, n_invalid_tensors=int(bad.sum()),
                    dt_seconds=time.time() - t0)
    if verbose:
        print(f"  invalid tensors zeroed: {int(bad.sum())};  elapsed {time.time() - t0:.1f}s")
    return nim


# --------------------------------------------------------------------------------------
# nim_eig / nim_fa
# --------------------------------------------------------------------------------------
def nim_eig(nim: Nim) -> Nim:
    """Recompute evec/eval from DT for mask==1 voxels (as nim_eig.m does; identical result to nim_dt_spd's own)."""
    X, Y, Z = nim.dims
    sel = nim.mask == 1
    Dm = _tensor_from_vec(nim.DT[sel])
    w, v = _eig_desc(Dm)
    evec = np.zeros((X, Y, Z, 3, 3)); eval_ = np.zeros((X, Y, Z, 3))
    evec[sel] = v; eval_[sel] = w
    nim.evec, nim.eval = evec, eval_
    return nim


def nim_fa(nim: Nim) -> Nim:
    l = nim.eval
    denom = np.sqrt((l ** 2).sum(-1))
    num = np.sqrt(0.5) * np.sqrt((l[..., 0] - l[..., 1]) ** 2 + (l[..., 1] - l[..., 2]) ** 2 + (l[..., 2] - l[..., 0]) ** 2)
    FA = np.zeros(nim.dims)
    ok = (nim.mask == 1) & (denom > 0)
    FA[ok] = num[ok] / denom[ok]
    nim.FA = FA
    return nim


# --------------------------------------------------------------------------------------
# preproc_mask_improvement (Strategy 2: fslmaths -fillh -s 0.5 -thr 0.5 -bin)
# --------------------------------------------------------------------------------------
def mask_improvement(mask: np.ndarray, FA: np.ndarray, zooms) -> tuple[np.ndarray, dict]:
    m = mask > 0.5
    total = m.size
    pct = 100 * m.sum() / total
    fa_out = FA[~m]
    mean_fa_out = float(np.nanmean(fa_out)) if fa_out.size else 0.0
    info = dict(orig_pct=pct, mean_fa_outside=mean_fa_out, strategy1_triggered=bool(pct > 50 or mean_fa_out > 0.1))
    # Strategy 1 (re-run BET -f 0.3 -R) cannot be reproduced without FSL; it would only trigger for a poor mask.
    filled = ndi.binary_fill_holes(m)
    sigma_vox = 0.5 / np.asarray(zooms, float)               # -s 0.5 is in mm
    smooth = ndi.gaussian_filter(filled.astype(np.float32), sigma=sigma_vox)
    filled = smooth > 0.5
    fpct = 100 * filled.sum() / total
    if 15 < fpct < 60:
        m = filled
        info["morph_applied"] = True
    else:
        info["morph_applied"] = False
    info["final_pct"] = 100 * m.sum() / total
    info["final_mean_fa_outside"] = float(np.nanmean(FA[~m]))
    return m.astype(np.uint8), info


# --------------------------------------------------------------------------------------
# preproc_tissue_segmentation (adaptive FA thresholds, WM eroded with fslmaths default 3x3x3 kernel)
# --------------------------------------------------------------------------------------
def tissue_segmentation(FA: np.ndarray, brain_mask: np.ndarray):
    bm = brain_mask > 0.5
    fa_in = FA[bm]
    wm_thr = max(np.percentile(fa_in, 60), 0.1)
    csf_thr = min(np.percentile(fa_in, 20), 0.1)
    if wm_thr <= csf_thr:
        wm_thr = csf_thr + 0.05
    wm_prelim = (FA > wm_thr) & bm
    wm = ndi.binary_erosion(wm_prelim, structure=np.ones((3, 3, 3), bool))
    gm = (FA > csf_thr) & (FA <= wm_thr) & bm
    csf = (FA <= csf_thr) & bm
    info = dict(wm_thr=float(wm_thr), csf_thr=float(csf_thr),
                wm_pct=100 * wm.sum() / bm.sum(), gm_pct=100 * gm.sum() / bm.sum(), csf_pct=100 * csf.sum() / bm.sum(),
                wm_prelim_pct=100 * wm_prelim.sum() / bm.sum())
    return wm.astype(np.uint8), gm.astype(np.uint8), csf.astype(np.uint8), info
