#!/usr/bin/env python3
"""
fibercup_prepare.py - lay out the FiberCup phantom for HINEC and build the scoring ROIs.

Takes the files of the Tractometer's FiberCup package (all_data_3x3x3_no_rois_mask.zip,
extracted: dwi/acq-averaged_b-1500.nii.gz, dwi/grad64.txt, dwi/wm.nii.gz) and writes what
config/fibercup.yml and bin/run_hinec.sh expect, plus the endpoint regions for scoring:

    <out>/fibercup_M.nii.gz                    white-matter mask   (preprocessing.mask_file)
    <out>/hinec_pipeline/fibercup_raw.nii.gz   DWI, b0 + 64 directions
    <out>/hinec_pipeline/fibercup.bval         0 1500 1500 ...     (one line)
    <out>/hinec_pipeline/fibercup.bvec         3 rows x 65         (voxel axes, as released)
    <out>/hinec_pipeline/fibercup_mrtrix.b     x y z b per line    (for tckgen -grad, baselines)
    <out>/scoring/fibercup_rois.nii.gz         labels 1..12, the two ends of the 7 bundles
    <out>/scoring/fibercup_rois.json           names, boxes, the 7 valid pairs, provenance
    <out>/scoring/fibercup_rois.png            the regions over the b0 image (needs matplotlib)

Usage:
    python3 scripts/fibercup_prepare.py --src ~/Downloads/fibercup --out data/fibercup

--pad-z N (default 1) adds N empty slices below and above the three acquired ones, in the DWI,
the mask and the region image, and moves the origin so that the world coordinates of the acquired
voxels do not change.  The padding is outside the mask: no tensor is fitted there and nothing is
seeded there.  It is there because HINEC reads its fields at the voxel centres 1..N and, along
each axis, accepts a point only in [m, N + 1 - m]: m = 1.1 with trilinear interpolation, 1.5 with
cubic or spline, 1 in the moving-frames tracker.  With N = 3 slices that is a slab 1 to 2 voxels
thick around the middle slice, and every seed outside it is discarded; one padding slice puts
the whole acquired slab inside the interval.  Use --pad-z 0 for the grid as released.

The gradient table is written as released, in the voxel axes.  With HINEC's tensor fit the
principal direction then follows the bundles; flipping or swapping an in-plane axis does not.

THE SCORING REGIONS ARE A RE-CREATION.  The Tractometer's own region mask for the FiberCup is
not in its download package (the data page offers the white-matter mask only; its source has
the region mask commented out as unavailable), so the regions here are boxes around the 12
bundle ends of the released white-matter mask, 3 to 5 voxels deep and reaching 2 to 4 voxels
beyond the mask.  Scores obtained with them are comparable between runs scored
with this file; they are not the numbers of the original Tractometer regions.

Requires numpy and nibabel.
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np

try:
    import nibabel as nib
except ImportError:
    sys.exit("nibabel is required:  pip install nibabel")

EXPECTED_SHAPE = (64, 64, 3)

# Bundle-end regions: boxes of inclusive 0-based voxel ranges (i0, i1, j0, j1) in the released grid
# (i = first array index = x, j = second = y; x to the right and y up shows the phantom as in the
# Tractometer's figure: the small U on the left, the fan on the right).  Every slice.
# The right end of bundle 3 and the upper-right end of bundle 6 touch; the border between their
# regions follows the mask (row 38 belongs to bundle 3 up to i = 48, to bundle 6 from i = 49).
ROI_BOXES = {
    1:  ("b1_left",    [(20, 26, 47, 53)]),                      # bundle 1, the large U: left arm
    2:  ("b1_right",   [(34, 41, 47, 53)]),                      # bundle 1: right arm
    3:  ("b2_top",     [(12, 19, 43, 48)]),                      # bundle 2: upper-left end
    4:  ("b3_left",    [(7, 13, 36, 42)]),                       # bundle 3, the straight one: left end
    5:  ("b3_right",   [(47, 53, 39, 42), (47, 48, 38, 38)]),    # bundle 3: right end
    6:  ("b6_right",   [(46, 54, 34, 37), (49, 54, 38, 38)]),    # bundle 6: upper-right end
    7:  ("b7_top",     [(42, 47, 43, 49)]),                      # bundle 7: top-right end
    8:  ("b5_right",   [(50, 56, 20, 25)]),                      # bundle 5: right end
    9:  ("b4_upper",   [(6, 11, 21, 26)]),                       # bundle 4, the small U: upper arm
    10: ("b4_lower",   [(8, 13, 15, 19)]),                       # bundle 4: lower arm
    11: ("fan_bottom", [(19, 27, 3, 10)]),                       # common lower-left end of bundles 5, 6, 7
    12: ("b2_bottom",  [(35, 43, 3, 10)]),                       # bundle 2: lower-right end
}
BUNDLES = {            # the 7 ground-truth bundles, numbered as in the Tractometer's FiberCup figure: (label, label)
    "bundle_1": (1, 2),
    "bundle_2": (3, 12),
    "bundle_3": (4, 5),
    "bundle_4": (9, 10),
    "bundle_5": (8, 11),
    "bundle_6": (6, 11),
    "bundle_7": (7, 11),
}


def find_file(src, names):
    hits = []
    for root, _, files in os.walk(src):
        for f in files:
            if f in names:
                hits.append(os.path.join(root, f))
    if not hits:
        sys.exit("Not found under %s: %s" % (src, " or ".join(names)))
    hits.sort(key=lambda p: (p.count(os.sep), p))
    return hits[0]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pad_z(vol, n):
    """n empty slices below and above (axis 2)."""
    if n == 0:
        return vol
    width = [(0, 0)] * vol.ndim
    width[2] = (n, n)
    return np.pad(vol, width, mode="constant", constant_values=0)


def save(data, affine, path, dtype):
    img = nib.Nifti1Image(np.asarray(data, dtype=dtype), affine)
    img.set_data_dtype(dtype)
    img.set_sform(affine, code=1)
    img.set_qform(affine, code=1)
    img.header.set_xyzt_units("mm", "sec")       # the same units in every file: HINEC compares them
    nib.save(img, path)


def build_rois(shape):
    lab = np.zeros(shape, dtype=np.int16)
    for label, (_, boxes) in ROI_BOXES.items():
        for i0, i1, j0, j1 in boxes:
            if (lab[i0:i1 + 1, j0:j1 + 1, :] != 0).any():
                sys.exit("internal error: ROI boxes overlap at label %d" % label)
            lab[i0:i1 + 1, j0:j1 + 1, :] = label
    return lab


def bundle_lengths(wm2d, lab2d, voxel_mm):
    """Longest shortest path inside the mask between the two end regions of every bundle, in mm:
    the length of the bundle from tip to tip (8-connected Dijkstra on the in-plane mask)."""
    import heapq
    nb = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
          (-1, -1, 2 ** 0.5), (-1, 1, 2 ** 0.5), (1, -1, 2 ** 0.5), (1, 1, 2 ** 0.5)]
    out = {}
    for name, (la, lb) in BUNDLES.items():
        best = 0.0
        targets = wm2d & (lab2d == lb)
        for src in np.argwhere(wm2d & (lab2d == la)):
            dist = np.full(wm2d.shape, np.inf)
            dist[tuple(src)] = 0.0
            heap = [(0.0, int(src[0]), int(src[1]))]
            while heap:
                d, i, j = heapq.heappop(heap)
                if d > dist[i, j]:
                    continue
                for di, dj, w in nb:
                    a, b = i + di, j + dj
                    if 0 <= a < wm2d.shape[0] and 0 <= b < wm2d.shape[1] and wm2d[a, b] and d + w < dist[a, b]:
                        dist[a, b] = d + w
                        heapq.heappush(heap, (d + w, a, b))
            reach = dist[targets]
            if not np.isfinite(reach).all():
                sys.exit("internal error: the ends of %s are not connected inside the mask" % name)
            best = max(best, float(reach.max()))
        out[name] = round(best * voxel_mm, 1)
    return out


def picture(b0, wm, lab, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        print("  (matplotlib not installed: no picture of the regions)")
        return
    fig, ax = plt.subplots(figsize=(8.2, 8))
    ax.imshow(b0.T, origin="lower", cmap="gray", interpolation="nearest")
    ax.contour(wm.T.astype(float), levels=[0.5], colors="#3aa0ff", linewidths=0.8)
    colors = plt.get_cmap("tab20")
    cx, cy = 31.5, 28.0
    for label, (name, boxes) in ROI_BOXES.items():
        c = colors((label - 1) % 20)
        for i0, i1, j0, j1 in boxes:
            ax.add_patch(Rectangle((i0 - 0.5, j0 - 0.5), i1 - i0 + 1, j1 - j0 + 1, fill=True, fc=c, ec=c, alpha=0.45, lw=1.2))
        i0, i1, j0, j1 = boxes[0]
        mx, my = (i0 + i1) / 2, (j0 + j1) / 2                    # label on the side away from the centre
        dx, dy = mx - cx, my - cy
        if abs(dx) > abs(dy):
            ax.text(i1 + 1.0 if dx > 0 else i0 - 1.0, my, "%d %s" % (label, name), color="white", fontsize=8,
                    ha="left" if dx > 0 else "right", va="center")
        else:
            ax.text(mx, j1 + 1.0 if dy > 0 else j0 - 1.0, "%d %s" % (label, name), color="white", fontsize=8,
                    ha="center", va="bottom" if dy > 0 else "top")
    ax.set_xlabel("i  (first array index, 0-based)")
    ax.set_ylabel("j  (second array index, 0-based)")
    ax.set_title("FiberCup: scoring regions (boxes), white-matter mask (blue), b0")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="folder holding the extracted Tractometer FiberCup files (searched recursively)")
    ap.add_argument("--out", default="data/fibercup", help="output folder (default: data/fibercup)")
    ap.add_argument("--acq", default="averaged", choices=["averaged", "1", "2"], help="acquisition (default: averaged, as in Aydogan & Shi 2021)")
    ap.add_argument("--bvalue", type=int, default=1500, choices=[650, 1500, 2000], help="b-value of the file to use (default 1500)")
    ap.add_argument("--pad-z", type=int, default=1, help="empty slices added below and above (default 1; 0 = grid as released)")
    ap.add_argument("--force", action="store_true", help="overwrite existing output files")
    a = ap.parse_args()
    if a.pad_z < 0:
        sys.exit("--pad-z must be >= 0")

    dwi_path = find_file(a.src, ["acq-%s_b-%d.nii.gz" % (a.acq, a.bvalue)])
    grad_path = find_file(a.src, ["grad64.txt"])
    wm_path = find_file(a.src, ["wm.nii.gz"])
    print("DWI      :", dwi_path)
    print("gradients:", grad_path)
    print("WM mask  :", wm_path)

    dwi_img = nib.load(dwi_path)
    wm_img = nib.load(wm_path)
    dwi = np.asanyarray(dwi_img.dataobj)
    wm = np.asanyarray(wm_img.dataobj) > 0
    grad = np.loadtxt(grad_path)

    # ---- checks ------------------------------------------------------------------------------
    if dwi.ndim != 4 or dwi.shape[:3] != EXPECTED_SHAPE:
        sys.exit("unexpected DWI shape %s (expected 64 x 64 x 3 x 65)" % (dwi.shape,))
    if wm.shape != EXPECTED_SHAPE:
        sys.exit("unexpected mask shape %s" % (wm.shape,))
    if not np.allclose(dwi_img.affine, wm_img.affine, atol=1e-5):
        sys.exit("DWI and mask have different affines")
    if grad.shape != (dwi.shape[3], 3):
        sys.exit("gradient table is %s, the DWI has %d volumes" % (grad.shape, dwi.shape[3]))
    norms = np.linalg.norm(grad, axis=1)
    if norms[0] > 1e-6 or np.abs(norms[1:] - 1).max() > 1e-3:
        sys.exit("gradient table: expected (0,0,0) first, then unit vectors")
    zooms = tuple(float(z) for z in dwi_img.header.get_zooms()[:3])
    print("grid     : %s, voxel %s mm, %d volumes; mask %d voxels (%s per slice)" % (
        dwi.shape[:3], zooms, dwi.shape[3], int(wm.sum()), [int(wm[:, :, k].sum()) for k in range(wm.shape[2])]))

    # ---- padding along z ---------------------------------------------------------------------
    p = a.pad_z
    affine = dwi_img.affine.copy()
    affine[:3, 3] = affine[:3, 3] - p * affine[:3, 2]          # voxel (i, j, k + p) keeps the world position of (i, j, k)
    dwi_p = pad_z(dwi, p)
    wm_p = pad_z(wm.astype(np.int16), p)
    lab = build_rois(EXPECTED_SHAPE)
    lab_p = pad_z(lab, p)
    if p:
        lab_p[:, :, :p] = lab[:, :, [0]]                        # the regions extend through the padding:
        lab_p[:, :, -p:] = lab[:, :, [-1]]                      # an end point there is still at that bundle end

    # every bundle end of the mask must lie inside its region
    covered = {label: int((wm & (lab == label)).sum()) for label in ROI_BOXES}
    if min(covered.values()) < 20:
        sys.exit("a scoring region holds almost no white matter (%s): is this the released wm.nii.gz?" % covered)

    # ---- write -------------------------------------------------------------------------------
    pipe = os.path.join(a.out, "hinec_pipeline")
    scor = os.path.join(a.out, "scoring")
    os.makedirs(pipe, exist_ok=True)
    os.makedirs(scor, exist_ok=True)
    out = {
        "mask": os.path.join(a.out, "fibercup_M.nii.gz"),
        "dwi": os.path.join(pipe, "fibercup_raw.nii.gz"),
        "bval": os.path.join(pipe, "fibercup.bval"),
        "bvec": os.path.join(pipe, "fibercup.bvec"),
        "mrtrix": os.path.join(pipe, "fibercup_mrtrix.b"),
        "rois": os.path.join(scor, "fibercup_rois.nii.gz"),
        "json": os.path.join(scor, "fibercup_rois.json"),
        "png": os.path.join(scor, "fibercup_rois.png"),
    }
    existing = [v for v in out.values() if os.path.exists(v)]
    if existing and not a.force:
        sys.exit("refusing to overwrite (use --force):\n  " + "\n  ".join(existing))
    mine = {os.path.basename(out[k]) for k in ("dwi", "bval", "bvec", "mrtrix")}
    stale = sorted(os.path.join(pipe, f) for f in os.listdir(pipe) if f not in mine)

    save(dwi_p, affine, out["dwi"], dwi_img.get_data_dtype())
    save(wm_p, affine, out["mask"], np.int16)
    save(lab_p, affine, out["rois"], np.int16)
    bvals = np.r_[0, np.full(dwi.shape[3] - 1, a.bvalue)]
    with open(out["bval"], "w") as fh:
        fh.write(" ".join("%d" % b for b in bvals) + "\n")
    with open(out["bvec"], "w") as fh:
        for c in range(3):
            fh.write(" ".join("%.6f" % v for v in grad[:, c]) + "\n")
    with open(out["mrtrix"], "w") as fh:
        for g, b in zip(grad, bvals):
            fh.write("%.6f %.6f %.6f %d\n" % (g[0], g[1], g[2], b))

    meta = {
        "description": "FiberCup endpoint regions for Tractometer-style scoring (re-created from the released WM mask; "
                       "not the Tractometer's original region mask)",
        "grid": {"shape": list(lab_p.shape), "voxel_mm": list(zooms), "pad_z": p,
                 "affine": [[float(v) for v in row] for row in affine]},
        "rois": {str(k): {"name": n, "boxes_ij_inclusive_0based": [list(b) for b in bs], "wm_voxels_inside": covered[k]}
                 for k, (n, bs) in ROI_BOXES.items()},
        "bundles": {k: list(v) for k, v in BUNDLES.items()},
        "bundle_length_mm": bundle_lengths(wm.any(2), lab[:, :, 0], zooms[0]),
        "source": {"dwi": os.path.basename(dwi_path), "dwi_sha256": sha256(dwi_path),
                   "gradients_sha256": sha256(grad_path), "wm_mask_sha256": sha256(wm_path),
                   "acquisition": a.acq, "bvalue": a.bvalue},
    }
    with open(out["json"], "w") as fh:
        json.dump(meta, fh, indent=1)
    picture(np.asarray(dwi[:, :, 1, 0], dtype=float), wm.any(2), lab, out["png"])

    print("\nWritten:")
    for k in ("mask", "dwi", "bval", "bvec", "mrtrix", "rois", "json"):
        print("  " + out[k])
    if os.path.exists(out["png"]):
        print("  " + out["png"])
    print("\ngrid written: %s (%d padding slice(s) below and above; acquired slices are k = %d..%d, 1-based)" % (
        lab_p.shape, p, p + 1, p + 3))
    if stale:
        print("\nWARNING: %s holds files of an earlier run. HINEC skips preprocessing when fibercup.nii.gz is\n"
              "there, and would then use the old grid. Move these away before running run_hinec.sh:\n  %s" % (
                  pipe, "\n  ".join(stale)))
    print("\nNext (FA floor under the phantom's FA, whose median in the mask is 0.12; every streamline kept;\n"
          "49 seeds per voxel = about 100 000 seeds):\n"
          "  ./bin/run_hinec.sh %s/fibercup fibercup.mat config/fibercup.yml \\\n"
          "      --set termination.fa_min=0.03 --set termination.min_arc=0 --set seeding.density=49" % pipe)


if __name__ == "__main__":
    main()
