"""Corpus-callosum ROI and streamline selection (atlas-free).

HINEC's own CC handling is either the JHU atlas labels (needs FSL registration, unavailable here) or the crude
5x5x5 box in tests/nim_tests/nim_test_corpus_callosum.m. Here the CC is defined the standard atlas-free way:
a midsagittal slab (±halfwidth voxels around the brain midline along the left-right axis) restricted to
high-FA voxels whose primary eigenvector is left-right oriented, keeping the largest 26-connected component.
Streamlines are then selected if they pass through that ROI and cross the midsagittal plane.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def lr_axis_from_affine(affine: np.ndarray) -> int:
    """Index of the voxel axis that is most left-right (world x)."""
    import nibabel as nib
    codes = nib.aff2axcodes(affine)
    for i, c in enumerate(codes):
        if c in ("L", "R"):
            return i
    return 0


def midline_index(brain_mask: np.ndarray, axis: int = 0) -> float:
    """Left-right centre of mass of the brain mask (continuous voxel index)."""
    com = ndi.center_of_mass(brain_mask > 0)
    return float(com[axis])


def build_cc_roi(FA: np.ndarray, v1: np.ndarray, brain_mask: np.ndarray, x_mid: float, lr_axis: int = 0,
                 halfwidth: int = 1, fa_thr: float = 0.4, lr_thr: float = 0.8):
    dims = FA.shape
    x = np.arange(dims[lr_axis])
    slab_1d = np.abs(x - x_mid) <= halfwidth + 0.5
    shape = [1, 1, 1]; shape[lr_axis] = dims[lr_axis]
    slab = np.broadcast_to(slab_1d.reshape(shape), dims)
    lr = np.abs(v1[..., lr_axis])
    cand = slab & (brain_mask > 0) & (FA > fa_thr) & (lr > lr_thr)
    lab, n = ndi.label(cand, structure=np.ones((3, 3, 3), bool))
    if n == 0:
        raise RuntimeError("No candidate CC voxels found; relax fa_thr / lr_thr")
    sizes = ndi.sum(cand, lab, index=np.arange(1, n + 1))
    order = np.argsort(sizes)[::-1]
    roi = lab == (1 + int(order[0]))
    info = dict(x_mid=x_mid, halfwidth=halfwidth, fa_thr=fa_thr, lr_thr=lr_thr, n_candidates=int(cand.sum()),
                n_components=int(n), component_sizes=sorted(sizes.astype(int).tolist(), reverse=True)[:5],
                roi_voxels=int(roi.sum()))
    return roi, info


def crosses_midline(tracks: list[np.ndarray], x_mid: float, lr_axis: int = 0) -> np.ndarray:
    out = np.zeros(len(tracks), bool)
    for i, t in enumerate(tracks):
        xs = t[:, lr_axis]
        out[i] = (xs.min() < x_mid) and (xs.max() > x_mid)
    return out


def track_lengths_mm(tracks: list[np.ndarray], voxel_size) -> np.ndarray:
    vs = np.asarray(voxel_size, float)
    return np.array([np.sqrt(((np.diff(t.astype(float), axis=0) * vs) ** 2).sum(1)).sum() for t in tracks])
