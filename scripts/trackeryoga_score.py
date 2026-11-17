#!/usr/bin/env python3
"""
trackeryoga_score.py - score tractography of the TrackerYoga crossings against the ground truth.

TrackerYoga's crossings (GroundTruthTracts/phantom_cross15 ... cross90.vtk) are pairs of straight
bundles, 25.6 mm long and 2.1 mm wide, 100 fibres each, crossing at their centres at 15, 30, 45,
60, 75 and 90 degrees. Each bundle's TUBE is its own volume: the points within R + 0.5 mm of its
axis (R = 1.06 mm, the largest distance of its fibres from the axis) and within 1 mm of its length.

PASSAGE TEST (the headline numbers). The tracts seeded on an ARM of a bundle - inside its tube and
outside the other bundle's tube, so beside the crossing rather than in it - run through the
crossing on one side of the seed. Each is classified by where it goes:
  complete  both ends at the caps of its own bundle (within 1 mm of them, on opposite sides of
            the crossing): through the crossing on its bundle
  switched  an end lies in the other bundle's tube only: it turned into the other bundle
  exited    an end lies outside both tubes: it left the bundles
  stopped   the rest: it ended inside its own tube short of the far cap (a tract that leaves
            the phantom's mask stops at once, so a tract pushed out of its bundle's side near
            the crossing ends here too, within 0.5 mm of the bundle)
For each tract the arc position of its far end (the end of the half that went towards the
crossing) along its bundle, from the crossing centre, is recorded as far_t (> 0: past the centre).
ALL TRACTS whose midpoint (by arc length) lies in the crossing's box (ground truth +- 1 mm) are
also scored, wherever they were seeded (complete / short = inside one tube all along, with
coverage >= 90 % / below; switched; off = the rest), with
  H         symmetric Hausdorff distance to the nearest ground-truth fibre, the measure of
            Aydogan & Shi (IEEE TMI 2021, Sec. III-A): "acceptable" when H < 4 mm
  d_fib     largest distance of the tract from the single ground-truth fibre it stays closest
            to (directed Hausdorff, tract -> fibre)
  coverage  the part of its bundle's length the tract spans (projection on the bundle axis)

Tracks are read from a HINEC run directory (tractography/tracks*.mat, 1-based voxel
coordinates, mapped to world mm with the DWI in <run>/intermediate, as fibercup_score.py does);
the ground truth is in world mm of the released image.  Writes <run>/scoring/trackeryoga/scores.json.

Usage (from the hinec repo root):
    python3 scripts/trackeryoga_score.py hinec_runs/run_*_trackeryoga_*
Requires numpy, nibabel and h5py.
"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fibercup_score as FS                     # load_mat_tracks, resolve, default_ref, load_world

ANGLES = (15, 30, 45, 60, 75, 90)


def load_tracks_seeds(tracks_file, ref_file):
    """tracks (list of (n, 3)) and seed points ((T, 3) or None) in world mm, aligned by index"""
    import h5py, nibabel as nib
    tracks, seeds = [], None
    with h5py.File(tracks_file, "r") as f:
        for ref in f["tracks"][()].ravel():
            a = np.asarray(f[ref][()], float)
            tracks.append(a.T if (a.ndim == 2 and a.shape[0] == 3) else np.zeros((0, 3)))
        if "track_meta" in f and "seed_points" in f["track_meta"]:
            sp = np.asarray(f["track_meta"]["seed_points"][()], float)
            seeds = sp.T if sp.shape[0] == 3 else sp
            if len(seeds) != len(tracks):
                seeds = None
    aff = nib.load(ref_file).affine
    A, t0 = aff[:3, :3].T, aff[:3, 3]
    world = [(t - 1.0) @ A + t0 for t in tracks]
    return world, (None if seeds is None else (seeds - 1.0) @ A + t0)


def read_vtk_lines(path):
    """ASCII legacy VTK polydata (POINTS + LINES) -> list of (n, 3) arrays."""
    tok = open(path).read().split()
    i = tok.index("POINTS"); n = int(tok[i + 1])
    pts = np.array(tok[i + 3:i + 3 + 3 * n], float).reshape(n, 3)
    j = tok.index("LINES"); nl = int(tok[j + 1]); k = j + 3; lines = []
    for _ in range(nl):
        c = int(tok[k]); lines.append(pts[np.array(tok[k + 1:k + 1 + c], int)]); k += 1 + c
    return lines


def crossing(gt_dir, ang):
    L = read_vtk_lines(os.path.join(gt_dir, "phantom_cross%d.vtk" % ang))
    A = np.array([l[0] for l in L]); B = np.array([l[-1] for l in L])
    d = (B - A) / np.linalg.norm(B - A, axis=1)[:, None]
    lab = (np.abs(d @ d[0]) > np.cos(np.radians(ang / 2.0))).astype(int)   # 1: bundle of fibre 0
    lab = 1 - lab                                                        # 0 = bundle A (fibre 0), 1 = bundle B
    u, cb = [], []
    for b in (0, 1):
        db = d[lab == b]; db = db * np.sign(db @ db[0])[:, None]; v = db.mean(0); u.append(v / np.linalg.norm(v))
        cb.append(((A[lab == b] + B[lab == b]) / 2).mean(0))
    R = []
    for b in (0, 1):                                    # bundle radius: fibre distance from the axis
        q = np.concatenate([L[i] for i in np.nonzero(lab == b)[0]]) - cb[b]
        R.append(float(np.linalg.norm(q - (q @ u[b])[:, None] * u[b], axis=1).max()))
    P = np.concatenate(L)
    return {"angle": ang, "A": A, "B": B, "bundle": lab, "u": np.array(u), "cb": np.array(cb), "R": np.array(R), "centre": P.mean(0),
            "half": float(np.median(np.linalg.norm(B - A, axis=1)) / 2), "lo": P.min(0) - 1.0, "hi": P.max(0) + 1.0,
            "measured_angle": float(np.degrees(np.arccos(abs(u[0] @ u[1]))))}


def seg_dist(P, A, B):
    """distances [len(P), len(A)] from the points P to the segments A-B"""
    AB = B - A; L2 = (AB ** 2).sum(1)
    t = np.clip(((P[:, None, :] - A[None]) * AB[None]).sum(2) / L2[None], 0.0, 1.0)
    return np.linalg.norm(P[:, None, :] - (A[None] + t[..., None] * AB[None]), axis=2)


def arclen(P):
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]


def in_tube(P, c, b, rtol=0.5, stol=1.0):
    """bool per point: inside bundle b's tube (radius R_b + rtol, length 2h + 2 stol)"""
    q = P - c["cb"][b]; t = q @ c["u"][b]
    d = np.linalg.norm(q - t[:, None] * c["u"][b], axis=1)
    return (d <= c["R"][b] + rtol) & (np.abs(t) <= c["half"] + stol)


def axial(p, c, b):
    q = p - c["cb"][b]; t = q @ c["u"][b]
    return t, np.linalg.norm(q - np.outer(t, c["u"][b]) if q.ndim == 2 else q - t * c["u"][b], axis=-1)


def arm_seed(x, c):
    """the bundle (0 / 1) whose arm the seed x lies on, or -1"""
    for b in (0, 1):
        t, d = axial(x, c, b)
        to, do = axial(x, c, 1 - b)
        if d <= c["R"][b] and abs(t) <= c["half"] - 1.0 and not (do <= c["R"][1 - b] + 0.5 and abs(to) <= c["half"] + 1.0):
            return b
    return -1


def passage(P, c, b, seed):
    """class of a tract seeded on the arm of bundle b, and the arc position of its far end"""
    inb, ino = in_tube(P, c, b), in_tube(P, c, 1 - b)
    t, d = axial(P[[0, -1]], c, b)
    ts, _ = axial(seed, c, b)
    far_t = float(np.max(-np.sign(ts) * t))
    e_b, e_o = inb[[0, -1]], ino[[0, -1]]
    cap = (np.abs(t) >= c["half"] - 1.0) & e_b
    if cap.all() and t[0] * t[1] < 0:
        cls = "complete"
    elif (e_o & ~e_b).any():
        cls = "switched"
    elif (~e_b & ~e_o).any():
        cls = "exited"
    else:
        cls = "stopped"
    return cls, far_t, bool(inb.all())


def score_tract(P, c):
    s = arclen(P)
    stride = max(1, int(round(0.5 / max(np.median(np.diff(s)), 1e-6))))
    Ps = P[::stride]
    if len(Ps) < 2 or (Ps[-1] != P[-1]).any():
        Ps = np.vstack([Ps, P[-1]])
    D = seg_dist(Ps, c["A"], c["B"])                              # [points, fibres]
    t2g = D.max(0)                                                # directed Hausdorff tract -> fibre
    best = np.argsort(t2g)[:5]
    H = np.inf
    for g in best:
        n = max(2, int(np.ceil(2 * c["half"] / 0.25)) + 1)
        Q = c["A"][g] + np.linspace(0, 1, n)[:, None] * (c["B"][g] - c["A"][g])
        g2t = np.sqrt(((Q[:, None, :] - P[None]) ** 2).sum(2)).min(1).max()
        H = min(H, max(t2g[g], g2t))
    d_fib = float(t2g.min())
    inA, inB = in_tube(P, c, 0), in_tube(P, c, 1)
    near = int(c["bundle"][int(np.argmin(t2g))])                 # bundle of the nearest fibre
    if inA.all() and inB.all():
        b = near
    elif inA.all():
        b = 0
    elif inB.all():
        b = 1
    else:
        b = -1
    ends = tuple(0 if (a and not bb) else 1 if (bb and not a) else -1 for a, bb in ((inA[0], inB[0]), (inA[-1], inB[-1])))
    bc = b if b >= 0 else near
    proj = (P - c["cb"][bc]) @ c["u"][bc]
    h = c["half"]
    cov = max(0.0, min(proj.max(), h) - max(proj.min(), -h)) / (2 * h)
    if b >= 0:
        cls = "complete" if cov >= 0.9 else "short"
    elif min(ends) >= 0 and ends[0] != ends[1]:
        cls = "switched"
    else:
        cls = "off"
    return {"cls": cls, "H": float(H), "d_fib": d_fib, "length": float(s[-1]), "coverage": float(cov),
            "bundle": int(bc), "ends": ends}


def label_of(run_dir):
    src = ""
    f = os.path.join(run_dir, "intermediate", "SOURCE.txt")
    if os.path.exists(f):
        for line in open(f):
            if line.startswith("source_nim:"):
                src = os.path.basename(line.split(":", 1)[1].strip()).replace(".mat", "")
    ov = os.path.join(run_dir, "overrides.txt")
    over = " ".join(l.strip() for l in open(ov)) if os.path.exists(ov) else ""
    alg = "RK4"
    if "algorithm=mmf" in over:
        alg = {"none": "MMF-4a", "hold": "MMF-4h", "march": "MMF-4s"}.get(
            next((w.split("=")[1] for w in over.split() if w.startswith("mmf.hold=")), "march"), "MMF")
    return src, alg, over


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--gt", default="data/trackeryoga/src/GroundTruthTracts")
    a = ap.parse_args()
    C = [crossing(a.gt, ang) for ang in ANGLES]
    for run in a.runs:
        tf, rd = FS.resolve(run)
        if tf is None:
            continue
        src, alg, over = label_of(rd)
        ref = FS.default_ref(rd)
        world, seeds = load_tracks_seeds(tf, ref)
        res = {"run": rd, "source": src, "tracker": alg, "overrides": over, "tracks_file": tf, "reference": ref,
               "n_tracts": len(world), "crossings": {}}
        mids = np.array([t[np.searchsorted(arclen(t), arclen(t)[-1] / 2)] if len(t) > 1 else (t[0] if len(t) else np.full(3, np.nan))
                         for t in world])
        for c in C:
            inbox = np.nonzero(((mids >= c["lo"]) & (mids <= c["hi"])).all(1))[0]
            inbox = [i for i in inbox if len(world[i]) > 1]
            rows = [score_tract(world[i], c) for i in inbox]
            n = len(rows)
            cnt = {k: sum(r["cls"] == k for r in rows) for k in ("complete", "short", "switched", "off")}
            Hs = np.array([r["H"] for r in rows]) if rows else np.zeros(0)
            pas = []
            if seeds is not None:
                for i in range(len(world)):
                    if len(world[i]) < 2 or not ((seeds[i] >= c["lo"]) & (seeds[i] <= c["hi"])).all():
                        continue
                    b = arm_seed(seeds[i], c)
                    if b >= 0:
                        cls, far_t, intube = passage(world[i], c, b, seeds[i])
                        pas.append({"i": i, "bundle": b, "cls": cls, "far_t": far_t, "in_tube": intube,
                                    "length": float(arclen(world[i])[-1])})
            pc = {k: sum(r["cls"] == k for r in pas) for k in ("complete", "switched", "exited", "stopped")}
            pc["complete_strayed"] = sum(r["cls"] == "complete" and not r["in_tube"] for r in pas)
            st = [r["far_t"] for r in pas if r["cls"] == "stopped"]
            pc["stopped_median_far_t"] = float(np.median(st)) if st else None
            res["crossings"][str(c["angle"])] = {
                "n": n, **cnt,
                "acceptable": int((Hs < 4.0).sum()), "mean_H": float(Hs.mean()) if n else None,
                "median_length": float(np.median([r["length"] for r in rows])) if n else None,
                "median_coverage": float(np.median([r["coverage"] for r in rows])) if n else None,
                "passage": {"n": len(pas), **pc}, "passage_tracts": pas,
                "tracts": [{"i": int(i), **r} for i, r in zip(inbox, rows)]}
        out = os.path.join(rd, "scoring", "trackeryoga")
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "scores.json"), "w") as fh:
            json.dump(res, fh)
        print("\n%s  %s  (%s)  %d tracts" % (src, alg, os.path.basename(rd), len(world)))
        print("   angle | passage test: n  complete switched  exited stopped | all tracts in the box: n  complete   short switched     off  H<4  mean H")
        for c in C:
            r = res["crossings"][str(c["angle"])]; q = r["passage"]
            n = max(r["n"], 1); m = max(q["n"], 1)
            print("   %4d  |  %5d   %5.1f%%  %5.1f%%  %5.1f%%  %5.1f%%  |  %5d   %5.1f%%  %5.1f%%  %5.1f%%  %5.1f%%  %5.1f%%  %5.2f" % (
                c["angle"], q["n"], 100 * q["complete"] / m, 100 * q["switched"] / m, 100 * q["exited"] / m, 100 * q["stopped"] / m,
                r["n"], 100 * r["complete"] / n, 100 * r["short"] / n, 100 * r["switched"] / n, 100 * r["off"] / n,
                100 * r["acceptable"] / n, r["mean_H"] or 0))


if __name__ == "__main__":
    main()
