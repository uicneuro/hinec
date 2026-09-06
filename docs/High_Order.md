# High-Order Tractography

`nim_tractography_hinec` extends the FACT baseline (`nim_tractography_standard`) along three
axes: **sub-voxel interpolation** of the direction field, **high-order integration** of the
streamline ODE, and **anatomically constrained tractography (ACT)**. The first two are what
turn voxel-hopping into the numerical solution of an initial value problem; the third
constrains where that solution is allowed to start and stop.

The three are configured under `tractography.interpolation`, `tractography.integrator` and
`tractography.act`. Every parameter named here is declared in
`src/nim_utils/nim_config_schema.m`, and the generated reference is
[YAML Config](YAML_CONFIG.md).

---

## 1. Interpolation of the direction field

FACT reads the principal eigenvector of the voxel a point falls in, so the direction field
is piecewise constant and jumps at every voxel face. HINEC interpolates it, which both
smooths the trajectory and — more importantly — makes the right-hand side of the ODE
differentiable enough for a high-order integrator to be worth using.

### What is interpolated: the dyadic, not the eigenvector

\( v_1 \) is a **line field**: the eigensolver's sign for each voxel is arbitrary, and two
adjacent voxels can describe the same orientation with opposite signs. Interpolating the
signed components between them computes a *difference* — the result shrinks toward zero and
its direction is meaningless. That is an artefact of the storage, not of the anatomy.

HINEC therefore interpolates the six unique components of the **dyadic**
\( v_1 v_1^{\mathsf T} \), which is invariant under \( v_1 \to -v_1 \), and extracts the
principal eigenvector of the interpolated dyadic (`nim_principal_dir`). This averages
*orientations* rather than vectors. A sign for the result is then chosen to point along the
incoming tangent; that is a choice of representative, not a modification of the field, so
the returned direction stays a pure function of position.

Under `field: csd` the field is multi-valued, and the peak nearest the incoming tangent is
taken before the spatial blend. That reduction is structural — a multi-valued field must be
collapsed to one direction somehow — and is not a tunable steering term.

### Kernels

`interpolation.method` selects the kernel. They differ in **smoothness**, and smoothness is
what caps the order an integrator can actually reach: a Runge–Kutta method of order \(p\)
needs a right-hand side with \(p\) continuous derivatives, no matter how many stages the
tableau has.

| `interpolation.method` | MATLAB kernel | Smoothness | Observed RK4 order |
|---|---|---|---|
| `trilinear` (default) | `linear` | \(C^0\) — kinked at every voxel face | 2.00 |
| `cubic` | `cubic` | \(C^1\) — Keys cubic **convolution**, not a spline | 3.06 |
| `spline` | `spline` | \(C^2\) — a genuine cubic spline | 4.00 |

The observed orders are measured on ISMRM 2015 data; the ladders, error metric and threats
to validity are in [Solution Verification](CONVERGENCE.md).

Interpolants are built once per run as `griddedInterpolant` objects with `'none'`
extrapolation, so a query outside the domain returns `NaN` and the streamline terminates
there. `cubic` and `spline` need a wider stencil, so their in-bounds test is
correspondingly tighter.

### Spatial sampling (`interpolation.upsample`)

`interpolation.upsample` resamples the field onto a grid of spacing `1/upsample` voxels
before the interpolants are built: above 1 refines, below 1 coarsens, 1 is the acquisition
grid. The coordinate frame is untouched, so positions, step sizes and reported lengths stay
in native voxel units and runs at different factors are directly comparable. The
\( u \to \infty \) limit is the continuous field implied by the measured samples — not
ground-truth anatomy.

---

## 2. High-order integration

FACT advances to the next voxel boundary along a constant direction, which is Euler stepping
with a geometrically chosen step. HINEC integrates

\[
\frac{dx}{ds} = v(x), \qquad |v| = 1
\]

with a fixed step `integrator.step` (in voxels) and a scheme chosen by
`integrator.method`:

| `integrator.method` | Scheme | Stages | Step |
|---|---|---|---|
| `euler` | forward Euler | 1 | fixed |
| `rk2` | midpoint | 2 | fixed |
| `rk4` (default) | classical Runge–Kutta | 4 | fixed |
| `rkf45` | Dormand–Prince embedded pair | 7 | adaptive (`integrator.adaptive`) |

RK4, the default:

\[
\begin{aligned}
k_1 &= v(r_n), \\
k_2 &= v(r_n + \tfrac{1}{2}h k_1), \\
k_3 &= v(r_n + \tfrac{1}{2}h k_2), \\
k_4 &= v(r_n + h k_3), \\
r_{n+1} &= r_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).
\end{aligned}
\]

Each stage is an interpolation of the direction field, sign-aligned to the incoming tangent;
a stage whose interpolation fails falls back to the previous stage's vector. The embedded
pair and its step-size control are covered in [RKF45](RKF.md).

Because the direction is a pure function of position — HINEC has no direction-dependent
steering term — classical Runge–Kutta order theory applies, and the observed orders in the
table above mean what they say.

### Termination

| Criterion | Config key | Note |
|---|---|---|
| FA floor | `termination.fa_min` | also defines the propagation mask, together with the brain mask |
| Turn budget | `termination.angle_max` | degrees per **voxel of arc**; the budget for one step is `angle_max × step` |
| Arc length | `termination.max_arc` | in voxels; `max_steps` is derived as `ceil(max_arc/step)` |
| Minimum length | `termination.min_arc` | tracks shorter than this arc length are discarded |
| Domain | — | leaving the propagation mask, or the interpolation domain, ends the track |

The turn budget is a **rate**, which is what makes it step-invariant: refining the step
tightens the per-step allowance in proportion. Consecutive tangents are sign-aligned, so a
measured turn never exceeds 90° and any budget above that is inert rather than merely loose.
`angle_max: 0` disables the criterion outright, which is what a control run should use.

Note that the propagation mask is derived from FA and the brain mask **independently of the
seed mask** — where a track may go is not where it may start. See
[Seeding Strategy](SEEDING_STRATEGY.md).

---

## 3. Anatomically constrained tractography (ACT)

ACT adds tissue-based rules so that streamlines propagate in white matter and terminate in
grey matter rather than in CSF or outside the brain.

**It is off by default** (`tractography.act: false`). It also requires all three tissue masks
— `nim.wm_mask`, `nim.gm_mask`, `nim.csf_mask` — which `main.m` produces via
`preproc_tissue_segmentation`. The tracker decides whether ACT is active purely from whether
it received the masks; with any of them missing, tracking falls back to FA-based and
mask-based termination only.

| Tissue at the new position | Action |
|---|---|
| WM | continue |
| GM | append the point and terminate — a valid endpoint |
| CSF | terminate — an invalid endpoint |
| outside the brain | terminate |

Termination reasons are counted per run and reported in the run's termination analysis
(`gm`, `csf`, `outside`, `fa`, `angle`, `max_steps`, `no_direction`), which is the quickest
way to see which criterion is actually binding.

---

## Summary

| Component | FACT (`standard`) | HINEC (`hinec`) |
|---|---|---|
| Direction lookup | discrete voxel eigenvector | interpolated dyadic, or nearest FOD peak |
| Position update | jump to the next voxel boundary | fixed step, Euler / RK2 / RK4 / RKF45 |
| Field smoothness | piecewise constant | \(C^0\), \(C^1\) or \(C^2\) by kernel |
| Termination | FA, angle, length | FA, angle, length, and optional ACT tissue rules |

Configs: `config/standard_dti.yml`, `config/hinec_dti.yml`, `config/hinec_csd.yml`, plus the
`hinec_dti_*` variants. See [Tractography Methods](TRACTOGRAPHY_METHODS.md) for how the
trackers are dispatched and [YAML Config](YAML_CONFIG.md) for every parameter.
