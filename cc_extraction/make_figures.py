#!/usr/bin/env python3
"""Figures for the corpus-callosum extraction: projections over FA, 3-D view, length histogram, whole-brain context."""
from __future__ import annotations

import argparse, os, pickle, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def seg_colors(t):
    d = np.diff(t, axis=0)
    n = np.linalg.norm(d, axis=1, keepdims=True); n[n == 0] = 1
    return np.abs(d / n)


def add_projection(ax, tracks, ax_pair, bg, extent, alpha=0.35, lw=0.6, max_tracks=6000, rng=None):
    """ax_pair = (horizontal axis index, vertical axis index) in voxel space."""
    ax.imshow(bg.T, cmap="gray", origin="lower", extent=extent, vmin=0, vmax=0.9, interpolation="nearest")
    rng = rng or np.random.default_rng(0)
    if len(tracks) > max_tracks:
        tracks = [tracks[i] for i in rng.choice(len(tracks), max_tracks, replace=False)]
    segs, cols = [], []
    for t in tracks:
        pts = t[:, list(ax_pair)]
        segs.append(np.stack([pts[:-1], pts[1:]], axis=1))
        cols.append(seg_colors(t))
    if segs:
        lc = LineCollection(np.concatenate(segs), colors=np.concatenate(cols), linewidths=lw, alpha=alpha)
        ax.add_collection(lc)
    ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3]); ax.set_aspect("equal"); ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--nim", required=True)
    ap.add_argument("--title", default="HINEC (Python port), hinec_default.yml")
    args = ap.parse_args()
    fig_dir = os.path.join(args.out, "figures"); os.makedirs(fig_dir, exist_ok=True)
    nim = pickle.load(open(args.nim, "rb"))
    cache = pickle.load(open(os.path.join(args.out, "tracks_cache.pkl"), "rb"))
    summ = json.load(open(os.path.join(args.out, "summary.json")))
    cc_tracks, L, sub, x_mid = cache["cc_tracks"], cache["cc_lengths"], cache["subsample"], cache["x_mid"]
    FA = nim.FA; X, Y, Z = FA.shape
    nim.img = None; nim.img_bi = None; nim.evec = None; nim.DT = None
    import gc; gc.collect()
    roi = np.asanyarray(__import__("nibabel").load(os.path.join(args.out, "cc_roi.nii.gz")).dataobj) > 0
    rz = np.argwhere(roi); yc = int(round(rz[:, 1].mean())); zc = int(round(rz[:, 2].mean())); xm = int(round(x_mid))
    rng = np.random.default_rng(0)

    # 1. three projections of the CC bundle over FA
    fig, ax = plt.subplots(1, 3, figsize=(18, 6.5))
    add_projection(ax[0], cc_tracks, (1, 2), FA[xm], (-0.5, Y - 0.5, -0.5, Z - 0.5), rng=rng)
    ax[0].contour(roi[xm].T, levels=[0.5], colors="w", linewidths=0.8, extent=(-0.5, Y - 0.5, -0.5, Z - 0.5))
    ax[0].set_title(f"Sagittal projection (FA at midline x={xm}, ROI in white)")
    add_projection(ax[1], cc_tracks, (0, 2), FA[:, yc, :], (-0.5, X - 0.5, -0.5, Z - 0.5), rng=rng)
    ax[1].set_title(f"Coronal projection (FA at y={yc})")
    add_projection(ax[2], cc_tracks, (0, 1), FA[:, :, zc], (-0.5, X - 0.5, -0.5, Y - 0.5), rng=rng)
    ax[2].set_title(f"Axial projection (FA at z={zc})")
    fig.suptitle(f"Corpus callosum streamlines — {len(cc_tracks)} tracks — {args.title}\n"
                 f"colour = local orientation (red L–R, green A–P, blue S–I); voxel coordinates, 2 mm isotropic", fontsize=12)
    plt.tight_layout(); plt.savefig(os.path.join(fig_dir, "cc_tracks_projections.png"), dpi=130); plt.close(fig)

    # 2. pseudo-3-D views: orthographic projection after rotation (cheap; matplotlib's 3-D collections are memory-hungry)
    def rotate(pts, elev_deg, azim_deg):
        e, a = np.deg2rad(elev_deg), np.deg2rad(azim_deg)
        Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        Rx = np.array([[1, 0, 0], [0, np.cos(e), -np.sin(e)], [0, np.sin(e), np.cos(e)]])
        return (pts - np.array([X, Y, Z]) / 2.0) @ Rz.T @ Rx.T
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))
    views = [(25, -35, "oblique (anterior-left-superior)"), (0, 0, "coronal, from the front"), (90, 0, "axial, from above")]
    n3d = min(4000, len(cc_tracks))
    pick = [cc_tracks[i] for i in rng.choice(len(cc_tracks), n3d, replace=False)] if len(cc_tracks) else []
    for a, (el, az, name) in zip(axes, views):
        segs, cols, depth = [], [], []
        for t in pick:
            r = rotate(t.astype(float), el, az)
            segs.append(np.stack([r[:-1, :2], r[1:, :2]], axis=1)); cols.append(seg_colors(t)); depth.append(r[:-1, 2])
        if segs:
            segs = np.concatenate(segs); cols = np.concatenate(cols); depth = np.concatenate(depth)
            order = np.argsort(depth)
            a.add_collection(LineCollection(segs[order], colors=cols[order], linewidths=0.5, alpha=0.45))
        a.set_aspect("equal"); a.autoscale(); a.axis("off"); a.set_title(name)
    fig.suptitle(f"Corpus callosum — {n3d} of {len(cc_tracks)} streamlines (colour = local orientation)")
    plt.tight_layout(); plt.savefig(os.path.join(fig_dir, "cc_tracks_3d.png"), dpi=120); plt.close(fig)

    # 3. whole-brain context vs CC (axial + coronal)
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.5))
    add_projection(ax[0], sub, (0, 1), FA[:, :, zc], (-0.5, X - 0.5, -0.5, Y - 0.5), alpha=0.25, lw=0.4, max_tracks=8000, rng=rng)
    ax[0].set_title(f"Whole-brain tractography (random {min(8000, len(sub))} of {summ['tracking']['n_valid_tracks']} tracks), axial")
    add_projection(ax[1], cc_tracks, (0, 1), FA[:, :, zc], (-0.5, X - 0.5, -0.5, Y - 0.5), rng=rng)
    ax[1].set_title(f"Corpus callosum subset ({len(cc_tracks)} tracks), axial")
    plt.tight_layout(); plt.savefig(os.path.join(fig_dir, "wholebrain_vs_cc_axial.png"), dpi=120); plt.close(fig)

    # 4. length histogram + termination reasons
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    if len(L):
        ax[0].hist(L, bins=60, color="#4c72b0"); ax[0].axvline(summ["cc"]["min_length_mm"], color="r", ls="--", label=f"config min_length = {summ['cc']['min_length_mm']:.0f} mm (not enforced by HINEC)")
        ax[0].set_xlabel("CC streamline length (mm)"); ax[0].set_ylabel("count"); ax[0].legend(); ax[0].set_title("Corpus callosum track lengths")
    tc = summ["tracking"]["term_counts"]; keys = [k for k, v in tc.items() if v]
    ax[1].barh(keys, [tc[k] for k in keys], color="#55a868"); ax[1].set_title("Half-track termination reasons (whole brain)"); ax[1].set_xlabel("count")
    plt.tight_layout(); plt.savefig(os.path.join(fig_dir, "cc_lengths_and_terminations.png"), dpi=120); plt.close(fig)
    print("figures written to", fig_dir)


if __name__ == "__main__":
    main()
