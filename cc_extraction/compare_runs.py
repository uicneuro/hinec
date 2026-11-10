#!/usr/bin/env python3
"""Compare two extraction runs (e.g. faithful HINEC interpolation vs sign-aligned eigenvector interpolation)."""
import argparse, json, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import nibabel as nib


def seg_colors(t):
    d = np.diff(t, axis=0); n = np.linalg.norm(d, axis=1, keepdims=True); n[n == 0] = 1
    return np.abs(d / n)


def proj(ax, tracks, pair, bg, extent, rng, max_tracks=5000):
    ax.imshow(bg.T, cmap="gray", origin="lower", extent=extent, vmin=0, vmax=0.9)
    if len(tracks) > max_tracks:
        tracks = [tracks[i] for i in rng.choice(len(tracks), max_tracks, replace=False)]
    segs = [np.stack([t[:-1, list(pair)], t[1:, list(pair)]], 1) for t in tracks]
    cols = [seg_colors(t) for t in tracks]
    if segs:
        ax.add_collection(LineCollection(np.concatenate(segs), colors=np.concatenate(cols), linewidths=0.5, alpha=0.35))
    ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3]); ax.set_aspect("equal"); ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True); ap.add_argument("--b", required=True)
    ap.add_argument("--label-a", default="HINEC as written (raw component interpolation)")
    ap.add_argument("--label-b", default="sign-aligned interpolation (proposed fix)")
    ap.add_argument("--fa", required=True, help="FA NIfTI"); ap.add_argument("--out", required=True)
    args = ap.parse_args()
    FA = np.asanyarray(nib.load(args.fa).dataobj); X, Y, Z = FA.shape
    runs = []
    for d, lab in ((args.a, args.label_a), (args.b, args.label_b)):
        c = pickle.load(open(os.path.join(d, "tracks_cache.pkl"), "rb"))
        s = json.load(open(os.path.join(d, "summary.json")))
        runs.append((lab, c, s))
    roi = np.asanyarray(nib.load(os.path.join(args.a, "cc_roi.nii.gz")).dataobj) > 0
    rz = np.argwhere(roi); yc = int(round(rz[:, 1].mean())); xm = int(round(runs[0][1]["x_mid"]))
    rng = np.random.default_rng(0)

    fig, ax = plt.subplots(2, 3, figsize=(18, 11))
    for r, (lab, c, s) in enumerate(runs):
        proj(ax[r, 0], c["cc_tracks"], (1, 2), FA[xm], (-0.5, Y - 0.5, -0.5, Z - 0.5), rng)
        ax[r, 0].set_title(f"{lab}\nsagittal — {len(c['cc_tracks'])} CC tracks", fontsize=10)
        proj(ax[r, 1], c["cc_tracks"], (0, 2), FA[:, yc, :], (-0.5, X - 0.5, -0.5, Z - 0.5), rng)
        ax[r, 1].set_title("coronal", fontsize=10)
    # histograms + stats
    a = ax[0, 2]
    for (lab, c, s), col in zip(runs, ("#c44e52", "#4c72b0")):
        a.hist(c["cc_lengths"], bins=60, range=(0, 220), histtype="step", lw=1.8, color=col, label=f"{lab}: n={len(c['cc_lengths'])}, median {np.median(c['cc_lengths']):.0f} mm")
    a.set_xlabel("CC streamline length (mm)"); a.set_ylabel("count"); a.legend(fontsize=8); a.set_title("Corpus callosum track lengths")
    a = ax[1, 2]
    for (lab, c, s), col in zip(runs, ("#c44e52", "#4c72b0")):
        a.hist(c["lengths_mm"], bins=80, range=(0, 200), histtype="step", lw=1.8, color=col, log=True,
               label=f"{lab}: n={s['tracking']['n_valid_tracks']}, median {np.median(c['lengths_mm']):.1f} mm")
    a.set_xlabel("whole-brain streamline length (mm)"); a.set_ylabel("count (log)"); a.legend(fontsize=8); a.set_title("Whole-brain track lengths")
    plt.tight_layout(); plt.savefig(args.out, dpi=120)

    rows = []
    for lab, c, s in runs:
        t = s["tracking"]; k = s["cc"]
        rows.append(dict(run=lab, valid_tracks=t["n_valid_tracks"], wb_median_mm=float(np.median(c["lengths_mm"])),
                         wb_mean_mm=t["wholebrain_len_mm"]["mean"], no_direction=t["term_counts"]["no_direction"],
                         angle=t["term_counts"]["angle"], gm=t["term_counts"]["gm"], csf=t["term_counts"]["csf"],
                         cc_tracks=k["roi_hit_and_cross_midline"], cc_ge45=k["cc_ge_min_length"],
                         cc_median_mm=k["cc_len_mm"]["median"], cc_p90_mm=k["cc_len_mm"]["p90"], minutes=t["elapsed_min"]))
    json.dump(rows, open(os.path.splitext(args.out)[0] + ".json", "w"), indent=2)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
