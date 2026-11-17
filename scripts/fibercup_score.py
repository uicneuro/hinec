#!/usr/bin/env python3
"""
fibercup_score.py - Tractometer connectivity scores of a tractogram on the FiberCup phantom.

The protocol of Cote et al. 2013 as used by Aydogan & Shi 2021 (Table I): every streamline is
classified by where its two END POINTS are.

    VC  valid connection    the ends are in the two regions of one of the 7 bundles
    IC  invalid connection  the ends are in two different regions that no bundle joins
    NC  no connection       anything else (an end outside every region, or both in the same one)
    VB  valid bundles       bundles with at least one VC streamline (0..7)
    IB  invalid bundles     distinct pairs of regions joined by IC streamlines
    overall = sqrt(IC^2 + NC^2 + (1 - VC)^2), with the three rates in [0, 1]; 0 is perfect

Usage:
    python3 scripts/fibercup_score.py hinec_runs/run_20261110_101500_fibercup
    python3 scripts/fibercup_score.py tracks.mat  --ref data/fibercup/hinec_pipeline/fibercup.nii.gz
    python3 scripts/fibercup_score.py sd_stream.tck                      # MRtrix / Trekker baselines
    python3 scripts/fibercup_score.py <run> <run> <run>                  # one table for several

A run directory is scored on its newest tractography/tracks*.mat, and the result is written to
<run>/scoring/fibercup/scores.json.  HINEC tracks are voxel coordinates (1-based) on the grid of
the preprocessed DWI; they are taken to world millimetres with the affine of --ref (default: the
DWI copied into <run>/intermediate, else the grid of the region image), and from there to the
grid of the region image.  .trk and .tck files are read in world millimetres.

The classification looks at end points only, as the protocol does, so a streamline that joins the
right two regions by a wrong route still counts as VC.  Those are reported separately ("VC by a
detour": VC streamlines longer than 1.3 times their bundle, whose tip-to-tip length inside the
mask is stored with the regions) and are NOT taken out of VC.

--export writes vc.trk, ic.trk and nc.trk beside the scores, for a picture like the paper's Fig. 4.
--min-length-mm L scores only streamlines of at least L mm (default 0: all of them; the paper
removed none).  --tolerance T also accepts an end point whose nearest region voxel is within T
voxels (default 0: the end point must be inside the region).

The regions come from scripts/fibercup_prepare.py and are a re-creation, not the Tractometer's
original mask, so compare runs scored here with each other rather than digit for digit with
published tables.

Requires numpy, nibabel and h5py (scipy only for MATLAB files older than v7.3).
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

try:
    import nibabel as nib
except ImportError:
    sys.exit("nibabel is required:  pip install nibabel")

DETOUR = 1.3                      # a VC streamline longer than this many bundle lengths took a wrong route

PAPER_TABLE_I = [                 # Aydogan & Shi, IEEE TMI 40(2), 2021, Table I: IB, NC%, IC%, VC%, overall
    ("SD_STREAM", 10, 80.4, 1.0, 18.6, 1.14),
    ("iFOD1", 23, 73.4, 7.7, 19.0, 1.10),
    ("iFOD1 (-rk4)", 18, 71.0, 4.3, 24.7, 1.04),
    ("iFOD2", 21, 77.8, 7.5, 14.7, 1.16),
    ("PTT (proposed)", 17, 65.6, 4.9, 29.4, 0.96),
]


def load_mat_tracks(path):
    """HINEC tracks .mat -> list of (N, 3) arrays, 1-based voxel coordinates as stored."""
    try:
        import h5py
    except ImportError:
        sys.exit("h5py is required to read HINEC's -v7.3 tracks:  pip install h5py")
    if h5py.is_hdf5(path):
        out = []
        with h5py.File(path, "r") as f:
            name = "tracks" if "tracks" in f else "tracts" if "tracts" in f else None
            if name is None:
                sys.exit("no 'tracks' in %s" % path)
            refs = f[name][()]
            if refs.dtype != object:                       # a single numeric track
                arr = np.asarray(refs, dtype=float)
                return [arr.T] if arr.ndim == 2 else []
            for ref in refs.ravel():
                arr = f[ref][()]
                # MATLAB writes N x 3 column-major, so h5py always sees 3 x N: transpose, also when N == 3
                if arr.ndim == 2 and arr.shape[0] == 3 and arr.shape[1] >= 2:
                    out.append(np.asarray(arr, dtype=float).T)
        return out
    from scipy.io import loadmat
    m = loadmat(path)
    cell = m.get("tracks", m.get("tracts"))
    if cell is None:
        sys.exit("no 'tracks' in %s" % path)
    return [np.asarray(c, dtype=float) for c in cell.ravel() if c.ndim == 2 and c.shape[1] == 3 and c.shape[0] >= 2]


def resolve(target):
    """-> (tracks_file, run_dir or None)"""
    if os.path.isdir(target):
        cands = sorted(glob.glob(os.path.join(target, "tractography", "tracks*.mat")), key=os.path.getmtime)
        if not cands:
            print("skipped %s: no tractography/tracks*.mat (run not finished, or failed)" % target)
            return None, target
        return cands[-1], target
    if not os.path.isfile(target):
        sys.exit("not found: %s" % target)
    return target, None


def default_ref(run_dir):
    if run_dir is None:
        return None
    for f in sorted(glob.glob(os.path.join(run_dir, "intermediate", "*.nii.gz"))):
        base = os.path.basename(f)
        if base.endswith("_M.nii.gz") or "mask" in base or "parcellation" in base or base == "b0.nii.gz":
            continue
        return f
    return None


def load_world(tracks_file, ref_file, roi_img):
    """-> list of (N, 3) arrays in world mm, and a note on how they were read."""
    low = tracks_file.lower()
    if low.endswith((".trk", ".tck")):
        t = nib.streamlines.load(tracks_file)
        return [np.asarray(s, dtype=float) for s in t.streamlines], "world mm from the file"
    tracks = load_mat_tracks(tracks_file)
    if ref_file:
        ref = nib.load(ref_file)
        affine, shape, note = ref.affine, ref.shape[:3], "voxel -> world with %s" % ref_file
    else:
        affine, shape, note = roi_img.affine, roi_img.shape[:3], "voxel -> world with the grid of the region image"
    if tracks:
        allp = np.vstack(tracks)
        lo, hi = allp.min(0), allp.max(0)
        if (lo < 0.49).any() or (hi > np.array(shape) + 0.51).any():
            sys.exit("track coordinates %s .. %s do not fit the %s grid of the reference: wrong --ref?" % (
                np.round(lo, 2), np.round(hi, 2), shape))
    A = affine[:3, :3].T
    t0 = affine[:3, 3]
    return [(t - 1.0) @ A + t0 for t in tracks], note


def end_labels(world_tracks, roi_img, tol):
    lab = np.asanyarray(roi_img.dataobj).astype(np.int32)
    inv = np.linalg.inv(roi_img.affine)
    ends = np.array([[t[0], t[-1]] for t in world_tracks]).reshape(-1, 3)
    vox = ends @ inv[:3, :3].T + inv[:3, 3]
    idx = np.floor(vox + 0.5).astype(int)
    # The phantom is planar and the regions are prisms through the slab: an end point just beyond
    # the first or last slice is still at that bundle end, whatever the padding of the region image.
    idx[:, 2] = np.clip(idx[:, 2], 0, lab.shape[2] - 1)
    inside = ((idx >= 0) & (idx < np.array(lab.shape))).all(1)
    out = np.zeros(len(ends), dtype=np.int32)
    out[inside] = lab[idx[inside, 0], idx[inside, 1], idx[inside, 2]]
    if tol > 0:
        cen = np.argwhere(lab > 0)
        cl = lab[lab > 0]
        for k in np.nonzero(out == 0)[0]:
            v = vox[k].copy()
            v[2] = min(max(v[2], 0.0), lab.shape[2] - 1.0)
            d2 = ((cen - v) ** 2).sum(1)
            j = int(np.argmin(d2))
            if d2[j] <= tol * tol:
                out[k] = cl[j]
    return out.reshape(-1, 2)


def score(world_tracks, roi_img, meta, tol=0.0, min_length_mm=0.0):
    names = {int(k): v["name"] for k, v in meta["rois"].items()}
    bundles = {b: tuple(sorted(p)) for b, p in meta["bundles"].items()}
    pair2bundle = {p: b for b, p in bundles.items()}
    lengths = np.array([np.sqrt((np.diff(t, axis=0) ** 2).sum(1)).sum() for t in world_tracks])
    n_all = len(world_tracks)
    keep = lengths >= min_length_mm
    tracks = [t for t, k in zip(world_tracks, keep) if k]
    n = len(tracks)
    res = {"n_streamlines": n, "n_removed_short": int(n_all - n), "min_length_mm": min_length_mm, "tolerance_vox": tol}
    cls = np.zeros(n, dtype="U2")
    if n == 0:
        res.update(VC=0.0, IC=0.0, NC=0.0, VB=0, IB=0, overall=float(np.sqrt(2.0)))
        return res, cls, keep
    ab = end_labels(tracks, roi_img, tol)
    lo, hi = ab.min(1), ab.max(1)
    per_bundle = {b: 0 for b in bundles}
    bundle_len = {b: [] for b in bundles}
    kept_len = lengths[keep]
    invalid = {}
    nc = {"no_end_in_a_region": 0, "one_end_in_a_region": 0, "both_ends_in_the_same_region": 0}
    for i in range(n):
        a, b = int(lo[i]), int(hi[i])
        if a > 0 and a != b:
            if (a, b) in pair2bundle:
                per_bundle[pair2bundle[(a, b)]] += 1
                bundle_len[pair2bundle[(a, b)]].append(kept_len[i])
                cls[i] = "VC"
            else:
                invalid[(a, b)] = invalid.get((a, b), 0) + 1
                cls[i] = "IC"
        else:
            cls[i] = "NC"
            nc["both_ends_in_the_same_region" if a > 0 else "one_end_in_a_region" if b > 0 else "no_end_in_a_region"] += 1
    # VC streamlines much longer than their bundle: the right two regions by a wrong route
    gt_len = meta.get("bundle_length_mm", {})
    detour = {b: int(sum(v > DETOUR * gt_len[b] for v in bundle_len[b])) for b in bundles if b in gt_len}
    vc = float((cls == "VC").sum()) / n
    ic = float((cls == "IC").sum()) / n
    ncr = float((cls == "NC").sum()) / n
    L = kept_len
    res.update(
        VC=vc, IC=ic, NC=ncr,
        VB=int(sum(c > 0 for c in per_bundle.values())), IB=len(invalid),
        overall=float(np.sqrt(ic ** 2 + ncr ** 2 + (1 - vc) ** 2)),
        VC_over_VC_plus_IC=(vc / (vc + ic) if vc + ic > 0 else None),
        VC_by_detour=(float(sum(detour.values())) / n if detour else None),
        VB_without_detours=(int(sum(per_bundle[b] - detour.get(b, 0) > 0 for b in bundles)) if detour else None),
        VC_by_detour_per_bundle=detour, detour_factor=DETOUR, bundle_length_mm=gt_len,
        valid_bundles={b: per_bundle[b] for b in sorted(per_bundle)},
        valid_bundle_length_mm={b: ({"median": float(np.median(v)), "min": float(np.min(v)), "max": float(np.max(v))} if v else None)
                                for b, v in sorted(bundle_len.items())},
        invalid_bundles={"%s--%s" % (names[a], names[b]): c for (a, b), c in sorted(invalid.items(), key=lambda kv: -kv[1])},
        no_connection=nc,
        length_mm={"mean": float(L.mean()), "median": float(np.median(L)), "max": float(L.max())},
    )
    return res, cls, keep


def export_trk(world_tracks, keep, cls, roi_img, out_dir):
    hdr = nib.streamlines.TrkFile.create_empty_header()
    hdr["voxel_sizes"] = roi_img.header.get_zooms()[:3]
    hdr["dimensions"] = roi_img.shape[:3]
    hdr["voxel_to_rasmm"] = roi_img.affine
    hdr["voxel_order"] = "".join(nib.aff2axcodes(roi_img.affine))
    kept = [t.astype(np.float32) for t, k in zip(world_tracks, keep) if k]
    for c in ("VC", "IC", "NC"):
        sl = [t for t, cc in zip(kept, cls) if cc == c]
        path = os.path.join(out_dir, c.lower() + ".trk")
        nib.streamlines.save(nib.streamlines.Tractogram(sl, affine_to_rasmm=np.eye(4)), path, header=hdr)
        print("  wrote %s (%d streamlines)" % (path, len(sl)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets", nargs="+", help="HINEC run directory, tracks .mat, or .trk / .tck file")
    ap.add_argument("--rois", default="data/fibercup/scoring/fibercup_rois.nii.gz", help="region image from fibercup_prepare.py")
    ap.add_argument("--ref", default=None, help="NIfTI on the grid the .mat tracks live on (default: see above)")
    ap.add_argument("--tolerance", type=float, default=0.0, help="voxels; accept an end point this close to a region (default 0)")
    ap.add_argument("--min-length-mm", type=float, default=0.0, help="score only streamlines at least this long (default 0)")
    ap.add_argument("--export", action="store_true", help="write vc.trk, ic.trk, nc.trk")
    ap.add_argument("--out", default=None, help="output folder (default: <run>/scoring/fibercup, or beside the tracks file)")
    ap.add_argument("--no-paper", action="store_true", help="do not print the paper's Table I under the results")
    a = ap.parse_args()

    if not os.path.isfile(a.rois):
        sys.exit("region image not found: %s (run scripts/fibercup_prepare.py, or pass --rois)" % a.rois)
    meta_path = a.rois.replace(".nii.gz", ".json")
    if not os.path.isfile(meta_path):
        sys.exit("region description not found: %s" % meta_path)
    roi_img = nib.load(a.rois)
    with open(meta_path) as fh:
        meta = json.load(fh)

    if a.out and len(a.targets) > 1:
        sys.exit("--out takes one target; several targets each get their own folder")
    rows = []
    for target in a.targets:
        tracks_file, run_dir = resolve(target)
        if tracks_file is None:
            continue
        ref = a.ref or default_ref(run_dir)
        world, note = load_world(tracks_file, ref, roi_img)
        res, cls, keep = score(world, roi_img, meta, a.tolerance, a.min_length_mm)
        res["tracks_file"] = os.path.abspath(tracks_file)
        res["rois"] = os.path.abspath(a.rois)
        res["coordinates"] = note
        out_dir = a.out or (os.path.join(run_dir, "scoring", "fibercup") if run_dir
                            else os.path.splitext(tracks_file)[0] + "_fibercup")
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "scores.json"), "w") as fh:
            json.dump(res, fh, indent=1)
        label = os.path.basename(os.path.normpath(run_dir)) if run_dir else os.path.basename(tracks_file)
        rows.append((label, res))
        print("\n%s\n  tracks: %s (%s)" % (label, tracks_file, note))
        print("  streamlines %d%s | VB %d/7  IB %d | NC %.1f%%  IC %.1f%%  VC %.1f%% | overall %.3f" % (
            res["n_streamlines"], " (%d shorter than %g mm removed)" % (res["n_removed_short"], a.min_length_mm) if res["n_removed_short"] else "",
            res["VB"], res["IB"], 100 * res["NC"], 100 * res["IC"], 100 * res["VC"], res["overall"]))
        if res["n_streamlines"]:
            print("  VC streamlines per bundle : " + "  ".join("%s:%d" % (b[-1], c) for b, c in res["valid_bundles"].items()))
            print("  their length, mm (min-median-max): " + "  ".join(
                "%s:%s" % (b[-1], "%.0f-%.0f-%.0f" % (v["min"], v["median"], v["max"]) if v else "-")
                for b, v in res["valid_bundle_length_mm"].items()))
            if res["VC_by_detour"] is not None:
                print("  VC by a detour (longer than %.1f x the bundle, %s mm): %d streamlines = %.1f%% of all; "
                      "bundles found without them: %d/7" % (
                          DETOUR, "/".join("%.0f" % res["bundle_length_mm"][b] for b in sorted(res["bundle_length_mm"])),
                          sum(res["VC_by_detour_per_bundle"].values()), 100 * res["VC_by_detour"], res["VB_without_detours"]))
            nc = res["no_connection"]
            print("  NC: no end in a region %d, one end %d, both ends in the same region %d" % (
                nc["no_end_in_a_region"], nc["one_end_in_a_region"], nc["both_ends_in_the_same_region"]))
            top = list(res["invalid_bundles"].items())[:5]
            if top:
                print("  largest invalid bundles: " + ", ".join("%s (%d)" % kv for kv in top))
        print("  -> %s" % os.path.join(out_dir, "scores.json"))
        if a.export:
            export_trk(world, keep, cls, roi_img, out_dir)

    if not rows:
        sys.exit("nothing to score")
    print("\n%-44s %9s %4s %4s %7s %7s %7s %8s" % ("", "n", "VB", "IB", "NC %", "IC %", "VC %", "overall"))
    for label, r in rows:
        print("%-44s %9d %4d %4d %7.1f %7.1f %7.1f %8.2f" % (label[:44], r["n_streamlines"], r["VB"], r["IB"],
                                                           100 * r["NC"], 100 * r["IC"], 100 * r["VC"], r["overall"]))
    if not a.no_paper:
        print("\nAydogan & Shi 2021, Table I (their regions and FODs, 100 000 streamlines; all found 7 valid bundles):")
        for name, ib, ncp, icp, vcp, ov in PAPER_TABLE_I:
            print("%-44s %9s %4d %4d %7.1f %7.1f %7.1f %8.2f" % (name, "100000", 7, ib, ncp, icp, vcp, ov))


if __name__ == "__main__":
    main()
