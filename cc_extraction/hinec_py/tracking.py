"""Vectorised port of src/nim_tractography/nim_tractography_hinec.m (+ the seeding logic of runTractography.m).

All streamlines of a chunk are advanced together; each streamline follows exactly the per-fibre control flow of
track_fiber_hinec():

    (dir, fa) = interp(seed);  no dir or fa < termination_fa  -> 'no_direction'  (1-point track)
    loop:
        fa < termination_fa                         -> 'fa'
        angle(dir, last step) > angle_thresh        -> 'angle'
        n_points >= max_steps                       -> 'max_steps'
        next = advance(pos, dir)                    (Euler / RK2 / RK4 / RKF45, per options.integration_order)
        next outside volume or outside dilated seed mask -> 'outside'
        tissue(round(next)): CSF -> discard half ('csf'); GM -> append & stop ('gm'); WM/UNKNOWN -> append
        (dir, fa) = interp(next); none -> 'no_direction'; sign-align to previous dir

Coordinates are 0-based voxel indices internally (MATLAB's 1-based positions minus 1); HINEC's margin/bounds tests
are translated accordingly. Interpolation reproduces griddedInterpolant: 'linear' (trilinear) or 'cubic'
(Keys cubic convolution, a = -1/2). Eigenvector components are interpolated independently and then normalised,
exactly as HINEC does (`align_signs=True` is an optional, non-HINEC improvement that flips the neighbouring
eigenvectors to agree with the incoming direction before interpolation).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage as ndi

TERM_CODES = {"running": 0, "no_direction": 1, "fa": 2, "angle": 3, "max_steps": 4, "outside": 5,
              "csf": 6, "gm": 7, "rkf_failure": 8, "unknown": 9}
TERM_NAMES = {v: k for k, v in TERM_CODES.items()}


# --------------------------------------------------------------------------------------
# Interpolation of (FA, v1x, v1y, v1z)
# --------------------------------------------------------------------------------------
class FieldInterpolator:
    def __init__(self, FA: np.ndarray, v1: np.ndarray, method: str = "linear", termination_fa: float = 0.2,
                 align_signs: bool = False):
        self.dims = np.asarray(FA.shape, dtype=np.int64)
        self.method = "cubic" if str(method).lower() == "cubic" else "linear"
        self.margin = 1.5 if self.method == "cubic" else 1.1
        self.termination_fa = float(termination_fa)
        self.align_signs = bool(align_signs)
        fields = np.concatenate([FA[..., None], v1], axis=-1).astype(np.float64)
        self.flat = np.ascontiguousarray(fields.reshape(-1, 4))
        self.sx = int(self.dims[1] * self.dims[2]); self.sy = int(self.dims[2])
        if self.method == "linear":
            self.offsets = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        else:
            self.offsets = [(a, b, c) for a in range(-1, 3) for b in range(-1, 3) for c in range(-1, 3)]

    @staticmethod
    def _keys_weights(t: np.ndarray) -> np.ndarray:
        """Cubic-convolution weights (Keys 1981, a=-1/2) for taps at -1,0,1,2 relative to floor(pos)."""
        t2, t3 = t * t, t * t * t
        return np.stack([-0.5 * t3 + t2 - 0.5 * t,
                         1.5 * t3 - 2.5 * t2 + 1.0,
                         -1.5 * t3 + 2.0 * t2 + 0.5 * t,
                         0.5 * t3 - 0.5 * t2], axis=1)

    def __call__(self, pos: np.ndarray, ref_dir: np.ndarray | None = None):
        """pos: (N,3) 0-based continuous voxel coords. Returns dir (N,3), fa (N,), valid (N,)."""
        N = pos.shape[0]
        direction = np.zeros((N, 3)); fa = np.zeros(N); valid = np.zeros(N, bool)
        if N == 0:
            return direction, fa, valid
        m = self.margin
        inside = np.all(pos >= m - 1.0, axis=1) & np.all(pos <= self.dims[None, :] - m, axis=1)
        if not inside.any():
            return direction, fa, valid
        p = pos[inside]
        i0 = np.floor(p).astype(np.int64)
        f = p - i0
        acc = np.zeros((p.shape[0], 4))
        if self.method == "linear":
            w1 = f; w0 = 1.0 - f
            for (a, b, c) in self.offsets:
                w = (w1[:, 0] if a else w0[:, 0]) * (w1[:, 1] if b else w0[:, 1]) * (w1[:, 2] if c else w0[:, 2])
                idx = (i0[:, 0] + a) * self.sx + (i0[:, 1] + b) * self.sy + (i0[:, 2] + c)
                vals = self.flat[idx]
                if self.align_signs and ref_dir is not None:
                    r = ref_dir[inside]
                    s = np.where(np.einsum("ij,ij->i", vals[:, 1:], r) < 0, -1.0, 1.0)
                    vals = vals.copy(); vals[:, 1:] *= s[:, None]
                acc += vals * w[:, None]
        else:
            wx, wy, wz = (self._keys_weights(f[:, k]) for k in range(3))
            hi = (self.dims - 1)
            for (a, b, c) in self.offsets:
                w = wx[:, a + 1] * wy[:, b + 1] * wz[:, c + 1]
                ix = np.clip(i0[:, 0] + a, 0, hi[0]); iy = np.clip(i0[:, 1] + b, 0, hi[1]); iz = np.clip(i0[:, 2] + c, 0, hi[2])
                vals = self.flat[ix * self.sx + iy * self.sy + iz]
                if self.align_signs and ref_dir is not None:
                    r = ref_dir[inside]
                    s = np.where(np.einsum("ij,ij->i", vals[:, 1:], r) < 0, -1.0, 1.0)
                    vals = vals.copy(); vals[:, 1:] *= s[:, None]
                acc += vals * w[:, None]
        fa_i = acc[:, 0]
        ok_fa = np.isfinite(fa_i) & (fa_i >= self.termination_fa)
        v = acc[:, 1:]
        nrm = np.linalg.norm(v, axis=1)
        ok = ok_fa & np.all(np.isfinite(v), axis=1) & (nrm > 1e-6)
        v_out = np.zeros_like(v)
        v_out[ok] = v[ok] / nrm[ok, None]
        fa_out = np.where(ok_fa, fa_i, 0.0)
        direction[inside] = v_out; fa[inside] = fa_out; valid[inside] = ok
        return direction, fa, valid


# --------------------------------------------------------------------------------------
# Seeding (generate_seed_points_hinec, 'uniform' lattice) and seed-mask logic (runTractography.m)
# --------------------------------------------------------------------------------------
def make_seed_mask(nim_mask: np.ndarray, FA: np.ndarray, seed_fa_threshold: float = 0.05) -> np.ndarray:
    """Strategy 1 of runTractography.m: preprocessed brain mask refined with FA > seed_fa_threshold."""
    return (nim_mask > 0.5) & (FA > seed_fa_threshold)


def generate_seed_points(seed_mask: np.ndarray, seed_density: float, strategy: str = "uniform", rng=None):
    idx = np.argwhere(seed_mask)                       # (V,3) 0-based voxel indices
    density = max(1.0, float(seed_density))
    if strategy.lower() == "random":
        rng = rng or np.random.default_rng(0)
        if density <= 1:
            return idx + (rng.random((len(idx), 3)) - 0.5) * 0.4, dict(description="random jitter", seeds_per_voxel=1)
        d = max(1, int(round(density)))
        seeds = np.repeat(idx, d, axis=0) + (rng.random((len(idx) * d, 3)) - 0.5) * 0.8
        return seeds, dict(description="random jitter", seeds_per_voxel=d)
    per_axis = max(1, int(np.ceil(density ** (1.0 / 3.0))))
    edges = np.linspace(-0.5, 0.5, per_axis + 1)
    off1 = (edges[:-1] + edges[1:]) / 2.0
    ox, oy, oz = np.meshgrid(off1, off1, off1, indexing="ij")
    offsets = np.stack([ox.ravel(), oy.ravel(), oz.ravel()], axis=1)    # (P,3)
    seeds = (idx[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
    info = dict(description=f"uniform grid ({per_axis}x{per_axis}x{per_axis} sub-voxel lattice)",
                seeds_per_voxel=len(offsets), voxel_spacing=1.0 / per_axis)
    return seeds, info


# --------------------------------------------------------------------------------------
# Tracker
# --------------------------------------------------------------------------------------
@dataclass
class TrackOptions:
    step_size: float = 0.2
    integration_order: int = 4
    interp_method: str = "trilinear"
    termination_fa: float = 0.2
    angle_thresh: float = 35.0
    max_steps: int = 800
    seed_density: float = 4
    seed_strategy: str = "uniform"
    adaptive_step: bool = False
    rkf_tolerance: float = 0.01
    rkf_safety: float = 0.9
    step_min: float = 0.01
    step_max: float = 1.0
    act_enabled: bool = True
    align_signs: bool = False              # non-HINEC option, off by default

    @classmethod
    def from_config(cls, tcfg: dict, **overrides):
        keys = cls.__dataclass_fields__.keys()
        kw = {k: tcfg[k] for k in keys if k in tcfg}
        kw.update(overrides)
        return cls(**kw)


def _round_half_away(x):
    """MATLAB round(): halves away from zero (positions are positive here)."""
    return np.floor(x + 0.5).astype(np.int64)


# Dormand–Prince RK5(4)7M tableau (as in rkf45_integration_step)
_DP_A = np.zeros((7, 7))
_DP_A[1, 0] = 1 / 5
_DP_A[2, :2] = [3 / 40, 9 / 40]
_DP_A[3, :3] = [44 / 45, -56 / 15, 32 / 9]
_DP_A[4, :4] = [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729]
_DP_A[5, :5] = [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656]
_DP_A[6, :6] = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]
_DP_BHAT = np.array([35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0])
_DP_B = np.array([5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40])


class HinecTracker:
    def __init__(self, FA, evec, seed_mask, opts: TrackOptions, wm_mask=None, gm_mask=None, csf_mask=None,
                 roi_masks: dict | None = None, voxel_size=(1.0, 1.0, 1.0)):
        self.dims = np.asarray(FA.shape)
        self.opts = opts
        v1 = evec[..., :, 0]
        self.interp = FieldInterpolator(FA, v1, opts.interp_method, opts.termination_fa, opts.align_signs)
        self.seed_mask = seed_mask.astype(bool)
        self.dilated = ndi.binary_dilation(self.seed_mask, structure=np.ones((3, 3, 3), bool))
        self.act = opts.act_enabled and all(m is not None for m in (wm_mask, gm_mask, csf_mask))
        if self.act:
            # tissue code volume: 0 unknown, 1 WM, 2 GM, 3 CSF (priority CSF > GM > WM as in check_tissue_type)
            tissue = np.zeros(FA.shape, np.int8)
            tissue[wm_mask > 0.5] = 1
            tissue[gm_mask > 0.5] = 2
            tissue[csf_mask > 0.5] = 3
            self.tissue = tissue
        else:
            self.tissue = None
        self.roi_masks = roi_masks or {}
        self.roi_names = list(self.roi_masks.keys())
        self.roi_stack = np.stack([self.roi_masks[k].astype(bool) for k in self.roi_names], axis=-1) if self.roi_names else None
        self.voxel_size = np.asarray(voxel_size, float)
        self.cos_thr = np.cos(np.deg2rad(opts.angle_thresh))

    # ---- integration -------------------------------------------------------------------------------------
    def _stage(self, pos, k_prev, ref):
        d, _, ok = self.interp(pos, ref)
        d = np.where(ok[:, None], d, k_prev)                        # fallback to previous k if interpolation fails
        flip = np.einsum("ij,ij->i", d, ref) < 0
        d[flip] *= -1.0
        return d

    def advance(self, pos, dirv, h):
        """Returns next_pos, accepted(bool), new_h — mirrors advance_position()."""
        o = self.opts.integration_order
        hh = h[:, None]
        if o == 1:
            return pos + hh * dirv, np.ones(len(pos), bool), h
        if o == 2:
            k2 = self._stage(pos + 0.5 * hh * dirv, dirv, dirv)
            return pos + hh * k2, np.ones(len(pos), bool), h
        if o == 4:
            k1 = dirv
            k2 = self._stage(pos + 0.5 * hh * k1, k1, dirv)
            k3 = self._stage(pos + 0.5 * hh * k2, k2, dirv)
            k4 = self._stage(pos + hh * k3, k3, dirv)
            return pos + (hh / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4), np.ones(len(pos), bool), h
        if o == 5:
            ks = [dirv]
            for s in range(1, 7):
                inc = sum(_DP_A[s, j] * ks[j] for j in range(s))
                ks.append(self._stage(pos + hh * inc, ks[-1], dirv))
            y_hat = pos + hh * sum(_DP_BHAT[j] * ks[j] for j in range(7))
            y_til = pos + hh * sum(_DP_B[j] * ks[j] for j in range(7))
            if not self.opts.adaptive_step:
                return y_hat, np.ones(len(pos), bool), h
            err = np.linalg.norm(y_hat - y_til, axis=1)
            tol, saf = self.opts.rkf_tolerance, self.opts.rkf_safety
            new_h = np.where(err > 1e-10, saf * h * (tol / np.maximum(err, 1e-300)) ** 0.2, saf * h * 2.0)
            new_h = np.clip(new_h, self.opts.step_min, self.opts.step_max)
            new_h = np.minimum(new_h, 2.0 * h)
            acc = err <= tol
            return np.where(acc[:, None], y_hat, pos), acc, new_h
        raise ValueError(f"unsupported integration_order {o}")

    # ---- main loop ---------------------------------------------------------------------------------------
    def track_batch(self, seeds: np.ndarray, signs: np.ndarray):
        """Track N half-streamlines (seed_i, sign_i). Returns dict with point buffer, counts, term codes, flags."""
        N = len(seeds)
        maxp = self.opts.max_steps + 1
        buf = np.zeros((N, maxp, 3), np.float32)
        buf[:, 0, :] = seeds
        npts = np.ones(N, np.int64)
        term = np.zeros(N, np.int8)
        discard = np.zeros(N, bool)
        roi_hit = np.zeros((N, len(self.roi_names)), bool)
        if self.roi_stack is not None:
            r = _round_half_away(seeds)
            roi_hit |= self.roi_stack[r[:, 0], r[:, 1], r[:, 2]]

        pos = seeds.astype(np.float64).copy()
        dirv, fa, ok = self.interp(pos, None)
        alive = ok & (fa >= self.opts.termination_fa)
        term[~alive] = TERM_CODES["no_direction"]
        dirv = dirv * signs[:, None]
        h = np.full(N, float(self.opts.step_size))
        dims = self.dims.astype(np.float64)
        A = np.flatnonzero(alive)

        for _ in range(maxp + 5):
            if A.size == 0:
                break
            p, d, f = pos[A], dirv[A], fa[A]
            # --- termination checks before advancing
            stop = np.zeros(A.size, np.int8)
            stop[f < self.opts.termination_fa] = TERM_CODES["fa"]
            has_prev = npts[A] > 1
            if has_prev.any():
                ia = A[has_prev]
                last = buf[ia, npts[ia] - 1, :].astype(np.float64) - buf[ia, npts[ia] - 2, :].astype(np.float64)
                ln = np.linalg.norm(last, axis=1)
                m = ln > 1e-6
                cosang = np.einsum("ij,ij->i", d[has_prev], last / np.maximum(ln, 1e-12)[:, None])
                bad = m & (cosang < self.cos_thr) & (stop[has_prev] == 0)
                tmp = stop[has_prev]; tmp[bad] = TERM_CODES["angle"]; stop[has_prev] = tmp
            stop[(npts[A] >= self.opts.max_steps) & (stop == 0)] = TERM_CODES["max_steps"]
            cont = stop == 0
            term[A[~cont]] = stop[~cont]
            A = A[cont]; p, d, f = p[cont], d[cont], f[cont]
            if A.size == 0:
                break
            # --- advance (with RKF retry loop, max 5 retries, as in track_fiber_hinec)
            hA = h[A]
            nxt, acc, newh = self.advance(p, d, hA)
            retries = 0
            while (~acc).any() and retries < 5:
                retries += 1
                r = np.flatnonzero(~acc)
                hA[r] = newh[r]
                n2, a2, h2 = self.advance(p[r], d[r], hA[r])
                nxt[r] = n2; acc[r] = a2; newh[r] = h2
            h[A] = newh
            fail = ~acc
            term[A[fail]] = TERM_CODES["rkf_failure"]
            keep = ~fail
            A, p, d, nxt = A[keep], p[keep], d[keep], nxt[keep]
            if A.size == 0:
                break
            # --- volume bounds and dilated brain mask
            out = np.any(nxt < 0, axis=1) | np.any(nxt > dims[None, :] - 1, axis=1)
            rv = np.zeros((A.size, 3), np.int64)
            rv[~out] = _round_half_away(nxt[~out])
            inb = ~out
            inb_idx = np.flatnonzero(inb)
            dm = self.dilated[rv[inb_idx, 0], rv[inb_idx, 1], rv[inb_idx, 2]]
            out2 = np.zeros(A.size, bool); out2[inb_idx[~dm]] = True
            outside = out | out2
            term[A[outside]] = TERM_CODES["outside"]
            keep = ~outside
            A, p, d, nxt, rv = A[keep], p[keep], d[keep], nxt[keep], rv[keep]
            if A.size == 0:
                break
            # --- ACT tissue check at next position
            if self.tissue is not None:
                tcode = self.tissue[rv[:, 0], rv[:, 1], rv[:, 2]]
            else:
                tcode = np.zeros(A.size, np.int8)
            is_csf = tcode == 3
            is_gm = tcode == 2
            # CSF: discard this half-track (track_length = 0)
            term[A[is_csf]] = TERM_CODES["csf"]; discard[A[is_csf]] = True; npts[A[is_csf]] = 0
            # GM and WM/UNKNOWN both append the point; GM then stops
            app = ~is_csf
            ia = A[app]
            buf[ia, npts[ia], :] = nxt[app]
            npts[ia] += 1
            if self.roi_stack is not None:
                rr = rv[app]
                roi_hit[ia] |= self.roi_stack[rr[:, 0], rr[:, 1], rr[:, 2]]
            term[A[is_gm]] = TERM_CODES["gm"]
            cont = app & ~is_gm
            A, p, d, nxt = A[cont], p[cont], d[cont], nxt[cont]
            if A.size == 0:
                break
            # --- move and re-interpolate direction
            pos[A] = nxt
            nd, nf, nok = self.interp(nxt, d)
            term[A[~nok]] = TERM_CODES["no_direction"]
            A2 = A[nok]
            nd = nd[nok]; nf = nf[nok]; dprev = d[nok]
            flip = np.einsum("ij,ij->i", dprev, nd) < 0
            nd[flip] *= -1.0
            dirv[A2] = nd; fa[A2] = nf
            A = A2
        return dict(buf=buf, npts=npts, term=term, discard=discard, roi_hit=roi_hit)

    def run(self, seeds: np.ndarray, chunk_seeds: int = 12000, keep_roi: str | None = None,
            subsample_prob: float = 0.0, rng=None, verbose: bool = True, min_length_mm: float = 0.0):
        """Whole run: bidirectional tracking of all seeds in chunks.

        Returns a dict with: lengths_pts, lengths_mm, term_counts, n_valid, kept streamlines (list of (n,3) 0-based
        voxel coords) for those hitting ROI `keep_roi`, and an optional random subsample of all valid streamlines.
        """
        rng = rng or np.random.default_rng(0)
        n_seeds = len(seeds)
        t0 = time.time()
        lengths_pts, lengths_mm = [], []
        kept, kept_idx, sub = [], [], []
        roi_counts = np.zeros(len(self.roi_names), np.int64)
        term_counts = {k: 0 for k in TERM_CODES}
        n_valid = 0
        half_discarded = 0
        ki = self.roi_names.index(keep_roi) if keep_roi else None
        for c0 in range(0, n_seeds, chunk_seeds):
            s = seeds[c0:c0 + chunk_seeds]
            n = len(s)
            both = np.concatenate([s, s], axis=0)
            signs = np.concatenate([np.ones(n), -np.ones(n)])
            res = self.track_batch(both, signs)
            buf, npts, term, disc, hit = res["buf"], res["npts"], res["term"], res["discard"], res["roi_hit"]
            for code in term[:n]:
                term_counts[TERM_NAMES[int(code)]] += 1
            for code in term[n:]:
                term_counts[TERM_NAMES[int(code)]] += 1
            nf, nb = npts[:n], npts[n:]
            half_discarded += int(disc.sum())
            for i in range(n):
                f_pts = nf[i]; b_pts = nb[i]
                fwd = buf[i, 1:f_pts] if f_pts > 1 else None
                bwd = buf[n + i, 1:b_pts][::-1] if b_pts > 1 else None
                if fwd is None and bwd is None:
                    continue
                parts = [x for x in (bwd, s[i][None, :].astype(np.float32), fwd) if x is not None]
                tr = np.concatenate(parts, axis=0)
                if tr.shape[0] <= 1:
                    continue
                n_valid += 1
                seg = np.diff(tr.astype(np.float64), axis=0) * self.voxel_size[None, :]
                L = float(np.sqrt((seg ** 2).sum(1)).sum())
                lengths_pts.append(tr.shape[0]); lengths_mm.append(L)
                if L < min_length_mm:
                    continue
                h = hit[i] | hit[n + i]
                roi_counts += h
                if ki is not None and h[ki]:
                    kept.append(tr); kept_idx.append(c0 + i)
                if subsample_prob > 0 and rng.random() < subsample_prob:
                    sub.append(tr)
            if verbose:
                done = min(c0 + chunk_seeds, n_seeds)
                el = time.time() - t0
                print(f"  seeds {done}/{n_seeds} ({100 * done / n_seeds:.1f}%)  valid tracks so far {n_valid}"
                      f"  kept({keep_roi}) {len(kept)}  elapsed {el / 60:.1f} min  ETA {el / done * (n_seeds - done) / 60:.1f} min", flush=True)
        return dict(lengths_pts=np.array(lengths_pts), lengths_mm=np.array(lengths_mm), term_counts=term_counts,
                    n_valid=n_valid, n_seeds=n_seeds, half_tracks_discarded_csf=half_discarded,
                    roi_counts=dict(zip(self.roi_names, roi_counts.tolist())), kept=kept, kept_seed_idx=np.array(kept_idx),
                    subsample=sub, elapsed_s=time.time() - t0)
