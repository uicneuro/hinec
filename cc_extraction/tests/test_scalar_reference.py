"""Reference check: a literal, per-streamline transcription of track_fiber_hinec() (RK4 path) vs the vectorised tracker."""
import pickle, sys
import numpy as np
sys.path.insert(0, "/home/claude/hinec_py")
from hinec_py import tracking, config

TERM = tracking.TERM_CODES


def interp1(itp, pos, ref=None):
    d, f, ok = itp(np.asarray(pos)[None, :], None if ref is None else np.asarray(ref)[None, :])
    return (d[0] if ok[0] else None), float(f[0])


def rk4(itp, pos, dirv, h):
    k1 = dirv
    k2, _ = interp1(itp, pos + 0.5 * h * k1); k2 = k1 if k2 is None else (-k2 if np.dot(k2, dirv) < 0 else k2)
    k3, _ = interp1(itp, pos + 0.5 * h * k2); k3 = k2 if k3 is None else (-k3 if np.dot(k3, dirv) < 0 else k3)
    k4, _ = interp1(itp, pos + h * k3);       k4 = k3 if k4 is None else (-k4 if np.dot(k4, dirv) < 0 else k4)
    return pos + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def track_scalar(tr: tracking.HinecTracker, seed, direction):
    o = tr.opts; dims = tr.dims
    pos = np.array(seed, float); track = [pos.copy()]
    dirv, fa = interp1(tr.interp, pos)
    if dirv is None or fa < o.termination_fa:
        return np.array(track), "no_direction"
    dirv = dirv * direction
    cos_thr = tr.cos_thr
    while True:
        if fa < o.termination_fa:
            return np.array(track), "fa"
        if len(track) > 1:
            step = track[-1] - track[-2]
            if np.linalg.norm(step) > 1e-6:
                if np.dot(dirv, step / np.linalg.norm(step)) < cos_thr:
                    return np.array(track), "angle"
        if len(track) >= o.max_steps:
            return np.array(track), "max_steps"
        nxt = rk4(tr.interp, pos, dirv, o.step_size)
        if np.any(nxt < 0) or np.any(nxt > dims - 1):
            return np.array(track), "outside"
        v = np.floor(nxt + 0.5).astype(int)
        if not tr.dilated[v[0], v[1], v[2]]:
            return np.array(track), "outside"
        tcode = tr.tissue[v[0], v[1], v[2]] if tr.tissue is not None else 0
        if tcode == 3:
            return np.zeros((0, 3)), "csf"
        track.append(nxt.copy())
        if tcode == 2:
            return np.array(track), "gm"
        pos = nxt
        nd, fa = interp1(tr.interp, pos, dirv)
        if nd is None:
            return np.array(track), "no_direction"
        if np.dot(dirv, nd) < 0:
            nd = -nd
        dirv = nd


if __name__ == "__main__":
    nim = pickle.load(open("/home/claude/cc_work/nim_stage2.pkl", "rb"))
    cfg = config.load_config_yaml("/mnt/user-data/uploads/hinec/config/hinec_default.yml")["tractography"]
    seed_mask = tracking.make_seed_mask(nim.mask, nim.FA, 0.05)
    seeds, _ = tracking.generate_seed_points(seed_mask, cfg["seed_density"])
    rng = np.random.default_rng(7)
    # bias the sample towards white matter so we test long tracks too
    fa_at = nim.FA[tuple(np.floor(seeds + 0.5).astype(int).T)]
    cand = np.flatnonzero(fa_at > 0.3)
    sel = np.concatenate([rng.choice(cand, 150, replace=False), rng.choice(len(seeds), 50, replace=False)])
    sub = seeds[sel]
    for method in ("trilinear", "cubic"):
        opts = tracking.TrackOptions.from_config(cfg, interp_method=method)
        tr = tracking.HinecTracker(nim.FA, nim.evec, seed_mask, opts, nim.wm_mask, nim.gm_mask, nim.csf_mask, voxel_size=nim.zooms)
        both = np.concatenate([sub, sub]); signs = np.concatenate([np.ones(len(sub)), -np.ones(len(sub))])
        res = tr.track_batch(both, signs)
        n_mismatch = 0; max_dev = 0.0; n_pts_total = 0
        for i in range(len(both)):
            ref, reason = track_scalar(tr, both[i], signs[i])
            n = res["npts"][i]
            vec = res["buf"][i, :n]
            n_pts_total += len(ref)
            if len(ref) != n or tracking.TERM_NAMES[int(res["term"][i])] != reason:
                n_mismatch += 1
                print(f"  MISMATCH seed {i}: scalar {len(ref)} pts/{reason} vs vector {n} pts/{tracking.TERM_NAMES[int(res['term'][i])]}")
            elif n > 0:
                max_dev = max(max_dev, float(np.abs(vec.astype(float) - ref).max()))
        print(f"[{method}] {len(both)} half-tracks, {n_pts_total} points: mismatches={n_mismatch}, max coordinate deviation={max_dev:.2e}")
