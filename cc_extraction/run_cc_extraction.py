#!/usr/bin/env python3
"""End-to-end: HINEC pipeline (Python port) on a dMRI sample -> whole-brain tractography -> corpus callosum extraction.

Usage:
    python run_cc_extraction.py --prefix data/original_sample/sample --config config/hinec_default.yml --out cc_extraction
"""
from __future__ import annotations

import argparse, json, os, pickle, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hinec_py import dti, tracking, cc, io as hio, config as hcfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--interp", default=None, help="override interp_method (trilinear|cubic)")
    ap.add_argument("--align-signs", action="store_true", help="non-HINEC: sign-align eigenvectors before interpolation")
    ap.add_argument("--seed-fa-threshold", type=float, default=0.05)
    ap.add_argument("--max-seeds", type=int, default=0, help="debug: limit number of seeds")
    ap.add_argument("--chunk", type=int, default=12000)
    ap.add_argument("--subsample", type=float, default=0.02)
    ap.add_argument("--stage-cache", default=None, help="pickle of a nim with DTI stage done (skips recomputation)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    fig_dir = os.path.join(args.out, "figures"); os.makedirs(fig_dir, exist_ok=True)
    T0 = time.time()
    log = {}

    cfg = hcfg.load_config_yaml(args.config)
    tcfg = cfg["tractography"]
    print("Tractography config:", tcfg, flush=True)

    # ---------------- DTI stage ----------------
    if args.stage_cache and os.path.exists(args.stage_cache):
        nim = pickle.load(open(args.stage_cache, "rb"))
        print("Loaded cached DTI stage", flush=True)
    else:
        nim = dti.nim_read(args.prefix)
        if nim.mask.sum() == 0:
            print("No brain mask next to the data: generating one from the mean b0 (median_otsu; FSL BET unavailable)", flush=True)
            nim.mask = dti.brain_mask_from_b0(nim.img_b0)
        nim = dti.nim_dt_spd(nim)
        nim = dti.nim_eig(nim)
        nim = dti.nim_fa(nim)
        nim.mask, minfo = dti.mask_improvement(nim.mask, nim.FA, nim.zooms)
        nim.wm_mask, nim.gm_mask, nim.csf_mask, tinfo = dti.tissue_segmentation(nim.FA, nim.mask)
        nim.info.update(mask_improvement=minfo, tissue=tinfo)
        if args.stage_cache:
            pickle.dump(nim, open(args.stage_cache, "wb"), protocol=4)
    log["dti"] = {k: (v if not isinstance(v, dict) else {kk: float(vv) if isinstance(vv, (np.floating, float)) else vv for kk, vv in v.items()})
                  for k, v in nim.info.items() if k != "prefix"}
    base = os.path.basename(args.prefix)
    hio.save_nifti(nim.mask, nim.affine, os.path.join(args.out, f"{base}_M.nii.gz"), np.uint8)
    for k, tag in (("wm_mask", "WM"), ("gm_mask", "GM"), ("csf_mask", "CSF")):
        hio.save_nifti(getattr(nim, k), nim.affine, os.path.join(args.out, f"{base}_{tag}_mask.nii.gz"), np.uint8)
    hio.save_nifti(nim.FA.astype(np.float32), nim.affine, os.path.join(args.out, f"{base}_FA_hinecpy.nii.gz"))
    hio.save_nifti(nim.evec[..., :, 0].astype(np.float32), nim.affine, os.path.join(args.out, f"{base}_V1_hinecpy.nii.gz"))

    # ---------------- CC ROI ----------------
    lr = cc.lr_axis_from_affine(nim.affine)
    x_mid = cc.midline_index(nim.mask, lr)
    roi, rinfo = cc.build_cc_roi(nim.FA, nim.evec[..., :, 0], nim.mask, x_mid, lr)
    log["cc_roi"] = rinfo
    hio.save_nifti(roi, nim.affine, os.path.join(args.out, "cc_roi.nii.gz"), np.uint8)
    print("CC ROI:", rinfo, flush=True)

    # ---------------- Tractography ----------------
    seed_mask = tracking.make_seed_mask(nim.mask, nim.FA, args.seed_fa_threshold)
    seeds, sinfo = tracking.generate_seed_points(seed_mask, tcfg["seed_density"], tcfg["seed_strategy"])
    if args.max_seeds:
        seeds = seeds[np.random.default_rng(0).choice(len(seeds), args.max_seeds, replace=False)]
    over = {}
    if args.interp:
        over["interp_method"] = args.interp
    if args.align_signs:
        over["align_signs"] = True
    opts = tracking.TrackOptions.from_config(tcfg, **over)
    print(f"Seeding: {int(seed_mask.sum())} voxels -> {len(seeds)} seeds ({sinfo['description']})", flush=True)
    print(f"Tracker options: {opts}", flush=True)
    tr = tracking.HinecTracker(nim.FA, nim.evec, seed_mask, opts, nim.wm_mask, nim.gm_mask, nim.csf_mask,
                               roi_masks={"cc": roi}, voxel_size=nim.zooms)
    res = tr.run(seeds, chunk_seeds=args.chunk, keep_roi="cc", subsample_prob=args.subsample)
    print(f"Tracking done in {res['elapsed_s'] / 60:.1f} min: {res['n_valid']} valid tracks from {len(seeds)} seeds", flush=True)

    # ---------------- CC selection ----------------
    kept = res["kept"]
    cross = cc.crosses_midline(kept, x_mid, lr)
    cc_tracks = [t for t, c in zip(kept, cross) if c]
    L = cc.track_lengths_mm(cc_tracks, nim.zooms)
    min_len = float(tcfg.get("min_length", 0))
    cc_long = [t for t, l in zip(cc_tracks, L) if l >= min_len]

    log["tracking"] = dict(n_seeds=int(len(seeds)), seeds_per_voxel=sinfo["seeds_per_voxel"], n_valid_tracks=int(res["n_valid"]),
                           half_tracks_discarded_csf=int(res["half_tracks_discarded_csf"]),
                           term_counts=res["term_counts"], elapsed_min=res["elapsed_s"] / 60,
                           wholebrain_len_mm=dict(mean=float(res["lengths_mm"].mean()), median=float(np.median(res["lengths_mm"])),
                                                  p90=float(np.percentile(res["lengths_mm"], 90)), max=float(res["lengths_mm"].max())),
                           options=opts.__dict__)
    log["cc"] = dict(roi_hit=int(len(kept)), roi_hit_and_cross_midline=int(len(cc_tracks)),
                     min_length_mm=min_len, cc_ge_min_length=int(len(cc_long)),
                     cc_len_mm=dict(mean=float(L.mean()) if len(L) else 0, median=float(np.median(L)) if len(L) else 0,
                                    p10=float(np.percentile(L, 10)) if len(L) else 0, p90=float(np.percentile(L, 90)) if len(L) else 0,
                                    max=float(L.max()) if len(L) else 0))
    print("CC:", log["cc"], flush=True)

    # ---------------- Exports ----------------
    shape = nim.dims
    hio.save_trk(cc_tracks, nim.affine, shape, nim.zooms, os.path.join(args.out, "tracks_cc.trk"))
    hio.save_trk(cc_long, nim.affine, shape, nim.zooms, os.path.join(args.out, f"tracks_cc_min{int(min_len)}mm.trk"))
    hio.save_trk(res["subsample"], nim.affine, shape, nim.zooms, os.path.join(args.out, "tracks_wholebrain_subsample.trk"))
    hio.save_tracks_mat(cc_tracks, os.path.join(args.out, "tracks_cc.mat"), opts.__dict__, res["elapsed_s"])
    hio.save_tracks_mat(res["subsample"], os.path.join(args.out, "tracks_wholebrain_subsample.mat"), opts.__dict__, res["elapsed_s"])
    hio.save_nim_mat(nim, os.path.join(args.out, "nim_hinecpy.mat"), extra={"cc_roi": roi.astype(np.uint8), "x_mid_1based": float(x_mid + 1)})
    pickle.dump(dict(cc_tracks=cc_tracks, cc_lengths=L, subsample=res["subsample"], lengths_mm=res["lengths_mm"],
                     lengths_pts=res["lengths_pts"], x_mid=x_mid, lr=lr), open(os.path.join(args.out, "tracks_cache.pkl"), "wb"), protocol=4)
    log["elapsed_total_min"] = (time.time() - T0) / 60
    json.dump(log, open(os.path.join(args.out, "summary.json"), "w"), indent=2, default=float)
    print("All done in %.1f min" % log["elapsed_total_min"], flush=True)


if __name__ == "__main__":
    main()
