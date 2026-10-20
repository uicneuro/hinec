# Moving frames tractography (`algorithm: mmf`)

`nim_tractography_mmf` follows a fibre as a curve that carries its own
orthonormal frame. Where the streamline trackers integrate the direction field
`dx/ds = v(x)` and re-read `v` from the data at every stage, the moving-frames
tracker advances the whole frame `E = (e1; e2; e3)` with the **connection
1-forms** of the data's frame field, evaluated on the direction of travel, and
the position with the first frame vector:

```
dx/ds = e1,        dE/ds = Omega E,        Omega_ij = w_ij(e1) = sum_k wij_k e1_k.
```

Omega is skew-symmetric with the two curvature components `w_12(e1)`,
`w_13(e1)` and the twist `w_23(e1)`. One step is the closed-form exponential of
this motion (a rotation of the frame together with a circular arc of the
position, by the Rodrigues formulas), so the frame stays exactly orthonormal
without re-orthogonalisation and the step is **exact for a planar arc of
constant curvature, at any step size**. Higher order comes from Lie group
integrators built around that step: Lie–Euler, the exponential midpoint rule and
the commutator-free fourth-order method CF4 (`integrator.method: euler | rk2 |
rk4`). On a fibre with torsion the tracker drops the twist (it is not a datum of
DTI, see below), so there Omega is not constant along the fibre in the carried
frame and the step is exact only to the order of the integrator. The scheme
itself is `advanceFrame.m`, shared with the moving-frames test suite (the
`Moving Frame Tractography/matlab` folder), whose notes
`doc/MMFtractography_notes.tex` give the mathematics, the gauge, the variants and
the verification on analytic fields in two and three dimensions. This page
covers only what connects the scheme to the voxel data of a nim.

## The frame field and its connection, from the voxels

The DTI principal eigenvector is a **line** field: the sign the eigensolver
writes for each voxel is arbitrary. The direction between voxels is therefore
read from the interpolated dyadic `v1 v1'`, which is sign-invariant — the same
six interpolants, kernel (`interpolation.method`) and `upsample` as the hinec
tracker — and the data's frame at a point `x` is

```
e1(x)   = principal axis of D(x) = interp(v1 v1')
e2, e3  = a completion of e1 (reference-axis projection)
w_1j,k  = e_j . d e1 / dx_k      j = 2, 3     (the curvature components)
w_23,k  = 0                                    (the twist is not a datum of DTI)
```

`d e1/dx_k` is a sign-aligned central difference of `e1` of the *interpolated*
field (`mmf.gradient_delta` voxels), so the connection is that of the field the
streamline trackers integrate. With the spline kernel that field is C2 and the
connection continuous; with the cubic and trilinear kernels the derivative of
the field jumps at the voxel faces, the connection is discontinuous there, and
the step — exact only for a connection that is constant along the arc — is
then as accurate as any integrator on a kinked right-hand side (see the
numbers below). `config/mmf_dti.yml` therefore uses the spline kernel. The
scheme runs with `'Twist' 'zero'` (a rotation-minimizing normal plane), so the
arbitrary completion `e2, e3` never enters the position: a step uses only the
curvature vector `(e1 . grad) e1` of the field, expressed in the carried frame
through the gauge `G = E Et'` that `advanceFrame` computes at every evaluation.
That gauge is also what absorbs the sign flips of the axial data — the "flip if
the dot product is negative" of streamline tracking, extended to the whole
frame.

## Modes (`mmf.mode`)

| mode | what the carried frame is | when |
|---|---|---|
| `transport` | carried by the connection only; the data enter through the curvature and the gauge | the frame scheme proper; the heading, once in error, is never corrected |
| `reanchor` | replaced by the data frame after every step (signs chosen by the transported frame) | the closest relative of streamline tracking, with an arc instead of a chord; the direction cannot drift from the data |
| `adaptive` (default) | transported, re-anchored only after a step at whose end the carried `e1` has drifted from the data's `e1` by `mmf.reanchor_tol` radians or more | trusts the transported heading while it agrees with the data and snaps it back where it does not |

The drift at every point and the re-anchoring decisions are recorded
(`meta.mmf.reanchored_fraction`; per step in `meta.trace` with `debug.trace`).
The default tolerance `1e-3` rad is an order of magnitude below the heading
error of a one-voxel CF4 step across a sharp bend and above the step-to-step
wobble of interpolated data at normal steps (see the suite's notes, section on
variants).

## Termination

The connection is returned non-finite wherever a point is outside the volume,
outside the propagation mask (or the brain mask), below the FA floor
`termination.fa_min` or without a well-conditioned principal axis;
`advanceFrame` stops a tract **before** a step at whose stages the connection
is not finite. The stages are the start of a step, its midpoint (orders 2 and
4) and its end (order 4), so with orders 1 and 2 the scheme never examines
the end of a step: the tracker drops a last node that lies outside tissue.
Every node of a track is therefore a valid integration node inside tissue.
Leaving the volume stops the tract as well, and the cut boundary point is
dropped, so tracks consist of integration nodes only, like those of the other
trackers.

The angle budget `termination.angle_max` (degrees per voxel of arc, through
`nim_angle_limit`: `min(angle_max * h, 90)` degrees per step, none for
`angle_max <= 0`) becomes `advanceFrame`'s `MaxTurn`. It bounds two turns: the
rotation `h |kappa|` of a step, evaluated at every stage, and the turn of the
heading as the data see it — the data's `e1` at the end of a step against the
heading at its start (modulo sign), in every mode. In the `reanchor` and
`adaptive` modes that is the turn a re-anchoring would realise, the
moving-frames counterpart of the angle between consecutive directions that the
streamline trackers test, which could otherwise put a kink of up to 90 degrees
into a track where the data's direction changes abruptly (a crossing seen in
DTI, an interface between bundles); in `transport` mode, which never
re-anchors, it stops a tract whose heading the data contradict beyond the
budget instead of letting it run through them. Where the interpolated
principal axis is piecewise constant with a swap plane, as across such an
interface, the connection is zero except within `gradient_delta` of the plane
and a step whose stages straddle it sees no curvature at all; this second test
is what stops the tract there. Unlike the streamline trackers, whose sign-aligned turn cannot
exceed 90 degrees, the 90 degree ceiling of `nim_angle_limit` is a real cap
here. Note also that the curvature of the interpolated principal-axis field is
large wherever the dyadic is nearly degenerate (crossings, planar regions:
the derivative of the principal axis grows like `1/(lambda_1 - lambda_2)`), so
the budget stops more half-tracks in such regions than the node-to-node angle
test of the streamline trackers does. `termination.max_arc` is `advanceFrame`'s
`smax`.

The reasons are counted over the half-tracks of the seeds that were tracked in
`meta.termination` (`max_arc`, `angle`, `fa`, `outside_mask`, `outside`,
`outside_domain`, `no_direction`); an `advanceFrame` stop on a non-finite
connection is resolved into one of these by probing the stages of the next
step and their difference stencils. Seeds outside tissue or without a
direction are not tracked and are counted by reason in
`meta.mmf.seeds_skipped`. ACT masks are not used by this tracker.

With `debug.trace` the per-step records of the traced seeds have the layout of
the hinec tracker's (`pos`, `dir` = the carried `e1`, `fa`, `h`, `chord`,
`turn` in degrees, `flip` always false, `n_steps`, `termination`), so
`nim_track_diagnose` reads them, plus `drift` (angle to the data's `e1`),
`reanchored`, `kappa` and `s`.

## What is not supported

- `field: csd`. The connection handles of `advanceFrame` are functions of
  position alone, and choosing an FOD peak needs the incoming tangent. Use the
  hinec tracker for CSD.
- `integrator.method: rkf45`. The exponential step has no embedded error
  estimate; the tracker rejects it (`mmf:integrator`).

## Verification

`tests/unit/TestMmfTracker.m` runs the tracker on the rigid-rotation phantom of
`TestIntegratorOrder` (circles of every radius, FA 1 in an annulus) from a seed
on the circle of radius 15.5, for 60 voxels of arc in each direction with
`h = 0.25`, and on the helix of `tests/fixtures/make_synthetic_nim.m` (radius
15, pitch 6, constant curvature and torsion) for 40 voxels each way:

| field | kernel | mode | integrator | max distance from the true curve (voxels) |
|---|---|---|---|---|
| circle | spline | any | euler / rk2 / rk4 | 1.1e-5 / 6e-6 / 5e-6 |
| circle | cubic | transport | euler / rk2 / rk4 | 3.7e-3 / 7e-4 / 7e-4 |
| circle | trilinear | transport | euler / rk2 / rk4 | 2.0e-2 / 1.5e-2 / 1.6e-2 |
| circle | cubic / trilinear | reanchor | rk4 | 4e-5 / 1.4e-4 |
| circle | trilinear | (Euler template) | | 0.48 |
| helix | spline | adaptive | euler / rk2 / rk4 | 1.0e-3 / 3.5e-4 / 5e-6 |

On the circle the exponential step is exact, so with the spline kernel the
residual is the interpolation error of the dyadic field and every mode and
order gives the same 1e-5; the curvature returned with the tract is `1/15.508`
to seven digits and nothing re-anchors at the default tolerance (the drift
stays at `4e-8` rad). With the cubic and trilinear kernels the residual of the
transport mode is dominated by the integration of the connection that jumps at
the voxel faces, which is why the orders separate there and why re-anchoring to
the data after every step (which follows the interpolated circle instead) is a
hundred times better on trilinear data; the adaptive mode keeps the transported
heading on these fields (no step drifts by `1e-3` rad) and coincides with
transport. On the helix the twist is dropped and the residual is the integration
error, with the orders 1, 2 and 4 cleanly separated; the curvature is recovered
to `4e-8`.

The test also pins the tracker contract (`nim_check_tracker_output`, seed
lattice and order, `seed_density`, random jitter, `min_length`, a run across
the chunk boundary of 4000 seeds, a seed that cannot be tracked), the straight
field (a line to roundoff), the sign alignment of the connection across the
flip locus of the principal axis, the precedence of `propagation_mask` over
the brain mask, the termination reasons (`max_arc`, `angle`, `fa`,
`outside_mask`, `outside`) with every integrator, the angle budget on both
the step and the data's heading (a 90-degree interface between bundles with
the stages straddling the swap plane, and the crossing fixture at eight
sub-voxel seed phases: every mode stops with `angle` and no track has a kink
above the budget), the safety ceiling when the budget is disabled, the
one-sided stencil at a volume face, a non-finite voxel in a spline field, the
per-step trace and the rejection of `rkf45` and `csd`.

On the sample dataset (`data/parcellation_sample/sample.mat`, 128 x 128 x 50,
`config/mmf_dti.yml` at `seeding.density: 1`, 147k seeds of which 100k lie in
tissue) the tracker produces 69.7k tracks of at least 15 voxels in about 11
minutes, every node inside the mask and no node-to-node turn above the budget;
`hinec_dti.yml` on the same seeds gives 74k in about 10 minutes, and for the
68.4k seeds both keep the two tracks agree to a median nearest-node distance of
0.25 voxel. The half-tracks end on the FA floor (77%), the volume faces (10%,
the slab is cut at both ends), the angle budget (9%, against 0.5% for hinec,
for the reason given above), the mask (3%) and a missing axis (1%).

## Files

- `src/nim_tractography/nim_tractography_mmf.m` — the tracker (seeds, field,
  connection handles, termination, metadata)
- `src/nim_tractography/advanceFrame.m` — the scheme (identical to the suite's copy)
- `config/mmf_dti.yml` — tracker config: `algorithm: mmf`, CF4, spline kernel,
  adaptive mode
- `tests/unit/TestMmfTracker.m`
