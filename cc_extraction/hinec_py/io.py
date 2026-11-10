"""Writers: NIfTI masks/maps, TrackVis .trk, and MATLAB .mat files in HINEC's own layout.

HINEC track convention: `tracks` is a cell array, each cell an N x 3 matrix of **1-based voxel coordinates**
(MATLAB indexing). Our internal streamlines are 0-based, so +1 is applied on export and the same nim-compatible
struct (FA, evec, eval, mask, tissue masks) is written so src/nim_visualization/visualizeTractography.m can load it.
"""
from __future__ import annotations

import numpy as np
import nibabel as nib
from nibabel.streamlines import Tractogram, TrkFile
from scipy.io import savemat


def save_nifti(data: np.ndarray, affine: np.ndarray, path: str, dtype=None):
    if dtype is not None:
        data = data.astype(dtype)
    img = nib.Nifti1Image(data, affine)
    img.header.set_xyzt_units("mm")
    nib.save(img, path)


def save_trk(tracks: list[np.ndarray], affine: np.ndarray, ref_img_shape, voxel_size, path: str):
    """Streamlines in 0-based voxel coordinates (voxel centre at integer index) -> TRK in RAS mm."""
    tg = Tractogram([np.asarray(t, np.float32) for t in tracks], affine_to_rasmm=affine)
    hdr = {"voxel_sizes": np.asarray(voxel_size, np.float32), "dimensions": np.asarray(ref_img_shape[:3], np.int16),
           "voxel_to_rasmm": affine.astype(np.float32), "voxel_order": "".join(nib.aff2axcodes(affine))}
    TrkFile(tg, header=hdr).save(path)


def save_tracks_mat(tracks: list[np.ndarray], path: str, options: dict, elapsed_time: float, algorithm: str = "hinec"):
    cell = np.empty((len(tracks), 1), dtype=object)
    for i, t in enumerate(tracks):
        cell[i, 0] = np.asarray(t, np.float64) + 1.0          # 1-based voxel coords as in MATLAB
    opts = {k: (v if not isinstance(v, bool) else float(v)) for k, v in options.items() if v is not None}
    savemat(path, {"tracks": cell, "options": opts, "elapsed_time": float(elapsed_time), "algorithm": algorithm},
            do_compression=True, long_field_names=True)


def save_nim_mat(nim, path: str, extra: dict | None = None, cc_roi: np.ndarray | None = None):
    """nim struct for visualizeTractography.m / nim_plot.m. If cc_roi is given it becomes parcellation_mask label 1
    (with nim.labels) so `visualizeTractography(tracks, nim, 'region', 1)` works. DT is omitted to keep the file small."""
    X, Y, Z = nim.dims
    d = {
        "FA": nim.FA, "evec": nim.evec, "eval": nim.eval, "mask": nim.mask.astype(np.uint8),
        "xdim": float(X), "ydim": float(Y), "zdim": float(Z), "size3": float(X * Y * Z),
        "bval": nim.bval.reshape(-1, 1), "bvec": nim.bvec, "thrsh_b0": float(nim.thrsh_b0),
        "hdr": {"ImageSize": np.array([X, Y, Z, nim.img.shape[3]], float), "PixelDimensions": np.asarray(nim.zooms, float)},
    }
    for k in ("wm_mask", "gm_mask", "csf_mask"):
        v = getattr(nim, k, None)
        if v is not None:
            d[k] = v.astype(np.uint8)
    if cc_roi is not None:
        d["parcellation_mask"] = cc_roi.astype(np.int32)
        d["labels"] = np.array(["Corpus callosum (midsagittal ROI)"], dtype=object).reshape(1, 1)
    if extra:
        d.update(extra)
    savemat(path, {"nim": d}, do_compression=True, long_field_names=True)
