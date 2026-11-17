#!/usr/bin/env python3
"""
trackeryoga_prepare.py - lay out the TrackerYoga phantoms for HINEC.

TrackerYoga (D. B. Aydogan, https://github.com/baranaydogan/TrackerYoga, BSD-3-Clause) holds the
synthetic phantoms of Aydogan & Shi, Parallel Transport Tractography, IEEE TMI 40(2), 2021,
Sec. III-A, simulated with MITK Diffusion: one 168 x 32 x 4 volume of 1 mm voxels with six arcs,
straight bundles crossing at 15, 30, 45, 60, 75 and 90 degrees and two fanning bundles side by
side, single shells at b = 1000, 2000, 3000 with and without noise, a multi-shell and an HCP-like
set, a mask, and the ground-truth tracts as VTK files (world mm of the image affine).

This writes, for each data set chosen with --sets (default b1000_no_noise b1000_with_noise):

    <out>/trackeryoga_M.nii.gz               the phantom mask           (preprocessing.mask_file)
    <out>/hinec_pipeline/ty_<tag>_raw.nii.gz DWI                        (run_hinec.sh input)
    <out>/hinec_pipeline/ty_<tag>.bval       one line
    <out>/hinec_pipeline/ty_<tag>.bvec       3 rows, as released
    <out>/trackeryoga.json                   grid, padding, tags, sha256 of the sources

with <tag> = b1000 for b1000_no_noise, b1000n for b1000_with_noise, and so on.

--pad N (default 1) adds N empty voxels on every side of the grid, in the DWI and the mask, and
moves the origin so that the world coordinates of the released voxels (and so of the ground
truth) do not change.  HINEC reads its fields at the voxel centres 1..N and accepts a point only
in [m, N + 1 - m] along each axis (m = 1.5 with the spline kernel, 1.1 trilinear, 1 in the
moving-frames tracker); the released slab is 4 voxels thick and several bundles touch the
y = 0 face, so without padding tracts would be cut half a voxel inside the outermost voxel
centres.

The gradient table is written as released.  Checked on b1000_no_noise with a log-linear tensor
fit: the principal direction is within 10 degrees of the local ground-truth direction in all
mask voxels with a single ground-truth direction (median 1.0 degree); flipping x or y or
swapping them gives a median of about 45 degrees.

Requires numpy and nibabel.  Usage (from the hinec repo root):
    python3 scripts/trackeryoga_prepare.py --src data/trackeryoga/src --out data/trackeryoga
"""
import argparse, hashlib, json, os, sys
import numpy as np
try:
    import nibabel as nib
except ImportError:
    sys.exit("nibabel is required:  pip install nibabel")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save(data, affine, path, dtype):
    img = nib.Nifti1Image(np.asarray(data, dtype=dtype), affine)
    img.set_sform(affine, code=1)
    img.set_qform(affine, code=1)
    nib.save(img, path)


def tag_of(name):
    b, noise = name.split("_", 1)
    return b + ("n" if noise == "with_noise" else "") if noise in ("no_noise", "with_noise") else name


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--src", required=True, help="the TrackerYoga folder (DWI/, Mask/, GroundTruthTracts/)")
    ap.add_argument("--out", required=True, help="output folder, e.g. data/trackeryoga")
    ap.add_argument("--sets", nargs="+", default=["b1000_no_noise", "b1000_with_noise"],
                    help="data sets of DWI/phantom_<set>.nii.gz (default: b1000_no_noise b1000_with_noise)")
    ap.add_argument("--pad", type=int, default=1, help="empty voxels added on every side (default 1)")
    a = ap.parse_args()
    if a.pad < 0:
        sys.exit("--pad must be >= 0")
    pipe = os.path.join(a.out, "hinec_pipeline")
    os.makedirs(pipe, exist_ok=True)
    mpath = os.path.join(a.src, "Mask", "phantom_mask.nii.gz")
    mimg = nib.load(mpath)
    mask = np.asanyarray(mimg.dataobj) > 0
    p = a.pad
    affine = mimg.affine.copy()
    affine[:3, 3] = affine[:3, 3] - p * affine[:3, :3].sum(1)      # voxel (i+p, j+p, k+p) keeps the world position of (i, j, k)
    pad3 = lambda v: np.pad(v, [(p, p)] * 3 + [(0, 0)] * (v.ndim - 3), mode="constant", constant_values=0)
    out_mask = os.path.join(a.out, "trackeryoga_M.nii.gz")
    save(pad3(mask.astype(np.int16)), affine, out_mask, np.int16)
    info = {"source": os.path.abspath(a.src), "pad": p, "mask": out_mask, "mask_voxels": int(mask.sum()),
            "grid": {"released_shape": list(mask.shape), "shape": [s + 2 * p for s in mask.shape],
                     "affine": [[float(v) for v in row] for row in affine]},
            "sha256": {"Mask/phantom_mask.nii.gz": sha256(mpath)}, "sets": {}}
    for name in a.sets:
        base = os.path.join(a.src, "DWI", "phantom_%s" % name)
        img = nib.load(base + ".nii.gz")
        if img.shape[:3] != mask.shape or not np.allclose(img.affine, mimg.affine, atol=1e-5):
            sys.exit("%s: grid differs from the mask's" % name)
        dwi = np.asanyarray(img.dataobj)
        bval = np.loadtxt(base + ".bvals").ravel()
        bvec = np.loadtxt(base + ".bvecs")
        if bvec.shape[0] != 3:
            bvec = bvec.T
        if bvec.shape != (3, dwi.shape[3]) or bval.size != dwi.shape[3]:
            sys.exit("%s: %d volumes, %d b-values, bvec %s" % (name, dwi.shape[3], bval.size, bvec.shape))
        tag = tag_of(name)
        stem = os.path.join(pipe, "ty_%s" % tag)
        save(pad3(dwi), affine, stem + "_raw.nii.gz", img.get_data_dtype())
        with open(stem + ".bval", "w") as fh:
            fh.write(" ".join("%g" % b for b in bval) + "\n")
        with open(stem + ".bvec", "w") as fh:
            for row in bvec:
                fh.write(" ".join("%.8g" % v for v in row) + "\n")
        info["sets"][tag] = {"set": name, "dwi": stem + "_raw.nii.gz", "volumes": int(dwi.shape[3]),
                             "b0": int((bval < 50).sum()), "bvalues": sorted(set(int(round(b)) for b in bval if b >= 50))}
        for ext in (".nii.gz", ".bvals", ".bvecs"):
            info["sha256"]["DWI/phantom_%s%s" % (name, ext)] = sha256(base + ext)
        print("ty_%s: %s, %d volumes (%d b0), b %s" % (tag, name, dwi.shape[3], info["sets"][tag]["b0"], info["sets"][tag]["bvalues"]))
    with open(os.path.join(a.out, "trackeryoga.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    print("grid %s -> %s (padding %d), mask %d voxels -> %s" % (mask.shape, info["grid"]["shape"], p, mask.sum(), out_mask))


if __name__ == "__main__":
    main()
