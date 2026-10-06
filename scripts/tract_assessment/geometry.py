"""Moving-frame (MMF / Peng) geometry of a single streamline.

Curve-level connection coefficients following the Frenet/Bishop prescription in
Chun & Peng's discussion note. Convention: frame vectors are ROWS and
    w_ij = <d e_i / ds, e_j>
so an ideal Frenet frame has w12 = kappa, w13 = 0, w23 = tau.

Which quantities mean what (see docs/TRACT_ASSESSMENT.md):
  kappa_w = sqrt(w12^2 + w13^2)  curvature magnitude; invariant to how the normal
                                 plane is oriented -> the PRIMARY descriptor.
  tau                            Frenet torsion; defined only where curvature is
                                 above threshold -> secondary.
  w23                            frame twist in the declared hybrid gauge; changes
                                 when the gauge changes -> diagnostic only.

Input curves must already be in WORLD MILLIMETRES (full NIfTI affine applied).
Measuring voxel-index curves on an anisotropic, rotated or reflected grid gives
different numbers and can flip the sign of torsion.
"""
import math
from dataclasses import dataclass, asdict
from functools import lru_cache

import numpy as np
from scipy.ndimage import convolve1d
from scipy.signal import savgol_coeffs


@dataclass(frozen=True)
class GeometryProtocol:
    """A frozen measurement definition. Values from different protocols are
    different measurements and must not be pooled."""
    spacing_mm: float = 0.5        # arc-length resampling before differentiation
    window_mm: float = 4.0         # Savitzky-Golay derivative window (physical)
    kappa_threshold: float = 1e-3  # Frenet above, Bishop transport below (mm^-1)
    summary_quantile: float = 0.9  # per-track summary = this quantile of |value|

    def as_dict(self):
        return asdict(self)

    @property
    def tag(self):
        return f"ds{self.spacing_mm:g}_window{self.window_mm:g}"


PRIMARY = GeometryProtocol()
SENSITIVITY = (GeometryProtocol(spacing_mm=0.25), GeometryProtocol(window_mm=2.0))
GEOMETRY_METRICS = ("kappa_w", "tau", "w23")


@lru_cache(maxsize=None)
def _savgol_operators(window, order, deriv, delta):
    """Cached linear operators equal to SciPy's Savitzky-Golay filter in 'interp' mode:
    the interior convolution kernel, plus the matrices that evaluate the
    least-squares edge polynomials (fitted to the first/last `window` samples)."""
    kernel = savgol_coeffs(window, order, deriv=deriv, delta=delta)
    fit = np.linalg.pinv(np.vander(np.arange(window, dtype=float), order + 1))

    def evaluate(positions):
        rows = np.zeros((len(positions), order + 1))
        for r, p in enumerate(positions):
            for k in range(order + 1):
                power = order - k
                if power >= deriv:
                    rows[r, k] = math.perm(power, deriv) * p ** (power - deriv)
        return rows @ fit / delta ** deriv

    half = window // 2
    return kernel, evaluate(range(half)), evaluate(range(window - half, window))


def _savgol(x, window, deriv, delta, order=3):
    """Savitzky-Golay derivative along axis 0 ('interp' edges), with the kernel and
    edge fits precomputed; equal to scipy.signal's filter up to rounding."""
    kernel, left, right = _savgol_operators(window, order, deriv, float(delta))
    y = convolve1d(x, kernel, axis=0, mode='constant')
    half = window // 2
    flat, out = x.reshape(len(x), -1), y.reshape(len(y), -1)
    out[:half] = left @ flat[:window]
    out[-half:] = right @ flat[-window:]
    return y


def _unit(x):
    return x / np.maximum(np.linalg.norm(x, axis=-1, keepdims=True), 1e-14)


def resample(points, ds):
    p = np.asarray(points, float)
    p = p[np.r_[True, np.linalg.norm(np.diff(p, axis=0), axis=1) > 1e-10]]
    if len(p) < 2:
        return None
    s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    targets = np.arange(0, s[-1] + 1e-10, ds)
    return np.column_stack([np.interp(targets, s, p[:, i]) for i in range(3)])


def _reference_normal(tangent):
    axis = np.zeros(3)
    axis[int(np.argmin(abs(tangent)))] = 1
    return _unit(axis - tangent * np.dot(axis, tangent))


def _transport(normal, t0, t1, eps=1e-10):
    cross = np.cross(t0, t1)
    sn = np.linalg.norm(cross)
    cs = float(np.clip(np.dot(t0, t1), -1, 1))
    if sn < eps:
        if cs >= 0:
            return _unit(normal - t1 * np.dot(normal, t1))
        axis = _reference_normal(t0)
        return _unit(2 * axis * np.dot(axis, normal) - normal)
    axis = cross / sn
    moved = normal * cs + np.cross(axis, normal) * sn + axis * np.dot(axis, normal) * (1 - cs)
    return _unit(moved - t1 * np.dot(moved, t1))


def connection(points, ds=.5, window_mm=4., kappa_threshold=1e-3, diagnostics=True):
    """Pointwise frame and connection coefficients of a world-mm polyline.

    Returns None when the curve is too short for the derivative window.
    Samples that cannot be measured (window edges, frame-switch windows,
    torsion where curvature is below threshold) are NaN, never zero.
    diagnostics=False skips the frame-ODE residual (not needed for summaries).
    """
    q = resample(points, ds)
    if q is None:
        return None
    w = max(5, int(round(window_mm / ds)) + 1)
    w += 1 - w % 2
    if len(q) < w + 4:
        return None
    d1 = _savgol(q, w, 1, ds)
    d2 = _savgol(q, w, 2, ds)
    d3 = _savgol(q, w, 3, ds)
    speed = np.linalg.norm(d1, axis=1)
    T = _unit(d1)
    K = (d2 - T * np.sum(d2 * T, axis=1)[:, None]) / np.maximum(speed[:, None] ** 2, 1e-14)
    kappa = np.linalg.norm(K, axis=1)
    N = _unit(K)
    high = kappa >= kappa_threshold
    # Frenet samples take the principal normal directly; only below-threshold
    # samples need the Bishop transport, which is sequential, so loop over those
    # alone (each reads the already-final e2 of its predecessor).
    e2 = N.copy()
    low = np.flatnonzero(~high)
    if low.size and low[0] == 0:
        e2[0] = _reference_normal(T[0])
    for i in low[low > 0]:
        e2[i] = _transport(e2[i - 1], T[i - 1], T[i])
    e2 = _unit(e2 - T * np.sum(e2 * T, axis=1)[:, None])
    e3 = _unit(np.cross(T, e2))
    frames = np.stack([T, e2, e3], axis=1)
    derivative = _savgol(frames, w, 1, ds)
    raw = np.einsum('nic,njc->nij', derivative, frames)
    omega = .5 * (raw - raw.transpose(0, 2, 1))
    cross = np.cross(d1, d2)
    tau = np.sum(cross * d3, axis=1) / np.maximum(np.sum(cross * cross, axis=1), 1e-20)
    margin = w // 2 + 1
    valid = np.ones(len(q), bool)
    valid[:margin] = False
    valid[-margin:] = False
    valid &= speed > .1
    # The gauge changes where curvature crosses the threshold. A derivative
    # window spanning such a switch represents neither smooth frame.
    switch_valid = np.ones(len(q), bool)
    for i in np.flatnonzero(high[1:] != high[:-1]) + 1:
        switch_valid[max(0, i - margin):min(len(q), i + margin + 1)] = False
    coefficient_valid = valid & switch_valid
    frenet_valid = valid & high
    tau[~frenet_valid] = np.nan
    kappa[~valid] = np.nan
    w12, w13, w23 = omega[:, 0, 1], omega[:, 0, 2], omega[:, 1, 2]
    for value in (w12, w13, w23):
        value[~coefficient_valid] = np.nan
    out = dict(p=q, s=np.arange(len(q)) * ds, frames=frames,
               w12=w12, w13=w13, w23=w23, kappa_w=np.hypot(w12, w13),
               kappa=kappa, tau=tau, coefficient_valid=coefficient_valid,
               frenet_valid=frenet_valid,
               convention='row frames; w_ij=<d e_i/ds,e_j>')
    if diagnostics:
        reconstructed = np.einsum('nij,njc->nic', omega, frames)
        out['frame_ode_residual'] = np.linalg.norm(derivative - reconstructed, axis=(1, 2))
    return out


def track_summary(curve_world, protocol=PRIMARY):
    """Per-track summary: the protocol quantile of |value| for each metric.

    Returns None if the track is too short to measure (counted as an exclusion,
    never silently dropped)."""
    measured = connection(curve_world, ds=protocol.spacing_mm, window_mm=protocol.window_mm,
                          kappa_threshold=protocol.kappa_threshold, diagnostics=False)
    if measured is None:
        return None
    out = {}
    for name in GEOMETRY_METRICS:
        values = measured[name]
        values = abs(values[np.isfinite(values)])
        out[name] = float(np.quantile(values, protocol.summary_quantile)) if len(values) else np.nan
    return out
