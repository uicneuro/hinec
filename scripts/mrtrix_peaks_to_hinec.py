#!/usr/bin/env python3
"""Convert an MRtrix sh2peaks NIfTI image to a HINEC CSD peak MAT file.

MRtrix peak vectors are in scanner RAS coordinates. HINEC tracks in the
reference image's voxel-axis coordinates. Peak magnitudes are preserved while
directions are transformed, and the source image is sampled on the reference
grid by nearest neighbor. This script reads no tractogram or scoring data.
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.io import savemat


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("peaks", type=Path, help="MRtrix sh2peaks NIfTI (3 channels per peak)")
    parser.add_argument("reference", type=Path, help="HINEC DWI or scalar NIfTI defining the voxel grid")
    parser.add_argument("output", type=Path, help="Output MAT file for tractography.csd.peaks_file")
    parser.add_argument("--min-amplitude", type=float, default=0.01,
                        help="Count peaks above this absolute amplitude (default: 0.01)")
    args = parser.parse_args()
    if not np.isfinite(args.min_amplitude) or args.min_amplitude < 0:
        parser.error("--min-amplitude must be finite and nonnegative")

    source = nib.load(str(args.peaks))
    reference = nib.load(str(args.reference))
    if len(source.shape) != 4 or source.shape[3] < 3 or source.shape[3] % 3:
        parser.error("peaks image must have 3, 6, 9, ... channels")
    dims = reference.shape[:3]
    if len(dims) != 3 or any(d < 2 for d in dims):
        parser.error("reference must have a three-dimensional voxel grid")
    nslots = source.shape[3] // 3
    source_vectors = source.get_fdata(dtype=np.float32).reshape(*source.shape[:3], nslots, 3)

    ijk = np.indices(dims, dtype=np.int32).reshape(3, -1).T
    world = nib.affines.apply_affine(reference.affine, ijk)
    source_ijk = np.rint(nib.affines.apply_affine(
        np.linalg.inv(source.affine), world)).astype(np.int32)
    inside = np.all((source_ijk >= 0) &
                    (source_ijk < np.asarray(source.shape[:3])), axis=1)
    vectors = np.zeros((len(ijk), nslots, 3), dtype=np.float32)
    valid_ijk = source_ijk[inside]
    vectors[inside] = source_vectors[valid_ijk[:, 0], valid_ijk[:, 1], valid_ijk[:, 2]]

    magnitude = np.linalg.norm(vectors, axis=-1)
    inverse_linear = np.linalg.inv(reference.affine[:3, :3]).astype(np.float32)
    voxel_vectors = np.einsum("ij,nkj->nki", inverse_linear, vectors)
    voxel_length = np.linalg.norm(voxel_vectors, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        voxel_vectors *= (magnitude / np.maximum(voxel_length, 1e-30))[..., None]
    voxel_vectors[voxel_length == 0] = 0

    peaks = voxel_vectors.reshape(*dims, nslots, 3)
    peak_w = magnitude.reshape(*dims, nslots)
    npeaks = np.count_nonzero(np.isfinite(peak_w) &
                              (peak_w > args.min_amplitude), axis=-1).astype(np.uint8)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    savemat(args.output, {"peaks": peaks, "npeaks": npeaks, "peak_w": peak_w},
            do_compression=True)
    print(f"Saved {args.output}: {np.count_nonzero(npeaks)} active voxels, "
          f"{nslots} peak slots")


if __name__ == "__main__":
    main()
