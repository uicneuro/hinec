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

By default (`mmf.connection: nodal`) the connection is computed at the voxel nodes
from the nodal `v1`, by fourth-order differences of the dyadic, and interpolated
between them with the kernel. With `mmf.connection: derivative`, `d e1/dx_k` is a
sign-aligned central difference of `e1` of the *interpolated*
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

The angle limit `termination.angle_max` (degrees per voxel of arc, through
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
re-anchors, it stops a tract whose heading the data contradict beyond what
`angle_max` allows instead of letting it run through them. Where the interpolated
principal axis is piecewise constant with a swap plane, as across such an
interface, the connection is zero except within `gradient_delta` of the plane
and a step whose stages straddle it sees no curvature at all; this second test
is what stops the tract there. Unlike the streamline trackers, whose sign-aligned turn cannot
exceed 90 degrees, the 90 degree ceiling of `nim_angle_limit` is a real cap
here. Note also that the curvature of the interpolated principal-axis field is
large wherever the dyadic is nearly degenerate (crossings, planar regions:
the derivative of the principal axis grows like `1/(lambda_1 - lambda_2)`), so
`angle_max` stops more half-tracks in such regions than the node-to-node angle
test of the streamline trackers does. `termination.max_arc` is `advanceFrame`'s
`smax`.

The reasons are counted over the half-tracks of the seeds that were tracked in
`meta.termination` (`max_arc`, `angle`, `fa`, `outside_mask`, `outside`,
`outside_domain`, `no_direction`, and with `mmf.hold` `coast`, `exit` and `focus`,
see below); an `advanceFrame` stop on a non-finite
connection is resolved into one of these by probing the stages of the next
step and their difference stencils. Seeds outside tissue or without a
direction are not tracked and are counted by reason in
`meta.mmf.seeds_skipped`. ACT masks are not used by this tracker.

With `debug.trace` the per-step records of the traced seeds have the layout of
the hinec tracker's (`pos`, `dir` = the carried `e1`, `fa`, `h`, `chord`,
`turn` in degrees, `flip` always false, `n_steps`, `termination`), so
`nim_track_diagnose` reads them, plus `drift` (angle to the data's `e1`),
`reanchored`, `kappa`, `s`, `coasted` and `confidence`.

## Degenerate data: MMF-4h and MMF-4s (`mmf.hold`)

In a voxel where two bundles cross, the single tensor is a pancake and its
principal axis means nothing. A tract that reads the connection there turns into
the other bundle or, more often, stops on `angle_max`. The hold scheme of
`advanceFrame` lets it cross instead: inside a zone of degenerate data it does not
read the data at all but **coasts**.

| `mmf.hold` | scheme | what the tract follows in the zone |
|---|---|---|
| `none` (default) | MMF-4a, the adaptive mode alone | the data, as everywhere |
| `hold` | MMF-4h | the mean connection of its last `hold_k` trusted nodes, in the carried frame, constant along the coast: a circle |
| `march` | MMF-4s | that curvature and fan-out, marched by the flatness of the connection: exact on circles, logarithmic spirals and rays, and turning the curvature by the torsion on a helix |

The presets `config/mmf_dti_hold.yml` and `config/mmf_dti_march.yml` are `mmf_dti`
with the scheme on, `reanchor_tol: 0.2` and `angle_max: 0` (see
[below](#angle-max) why). See [the scores](#on-data) below: on the
ISMRM 2015 phantom the scheme halves the wrong-path streamlines at the price of
shorter tracks, and most of what the presets gain in F1 comes from `reanchor_tol`
alone.

**The zone.** The confidence is the linear shape `CL = (l1 - l2)/l1` of the voxel
tensors (`nim.eval`), interpolated linearly and eroded by `hold_guard` voxels (the
minimum over the point and 14 points at that distance); a point is in the zone
where it is below `hold_min`. For a single bundle of eigenvalues 1.7, 0.3, 0.3
(CL 0.82) the default 0.2 is an equal crossing of about 80 degrees, 0.4 one of
about 70. But CL also falls with the anisotropy of the data, so **the floor belongs
to the dataset**: 0.15 and 0.2 were the best of the floors scored on the phantom,
and on a dataset of lower anisotropy 0.2 can be the tissue itself rather than its
crossings. The tracker reports the share of the tissue voxels below the floor, before
the erosion (`meta.mmf.zone_fraction`, and on the console), and warns above 40%
(`mmf:holdZone`). A voxel without a fit, outside the mask or below the FA floor is
no data: CL is interpolated over the voxels of tissue alone and the erosion reads
nothing else, so neither the edge of the mask nor the FA floor is a zone.

**Leaving it.** In the `adaptive` and `reanchor` modes, a tract reaching a trusted
node is re-anchored if its heading is within `exit_tol` of the data; otherwise the
coast has left the bundle and it stops (`exit`). The heading test of
`angle_max` comes first: data at the end of the step further from the heading
at its start than `angle_max` allows stop the tract with `angle`. In `transport` mode nothing is re-anchored or tested on leaving: the tract
keeps its heading and only `angle_max` applies. `exit_tol` is apart from
`reanchor_tol` because they ask different things of the data: a test on leaving has
to sit above the disagreement of neighbouring voxels. In the tissue of the sample
dataset (FA at least 0.2) a quarter of the neighbouring voxels differ by more than
0.2 rad, and with the 0.2 rad of the suite's noise-free fields as the test 23% of its
half-tracks ended on it. `exit_delay` carries the coast on for that
arc beyond the edge before the test; those steps neither count towards `coast_max`
nor reset it. A coast longer than `coast_max` voxels stops (`coast`), and so does a
marched coast whose fan closes to a tenth of its cross-section, the density ten
times the ring's (`focus`). `angle_max`, when it is active, applies to the coast's
own rotation.
**The volume, the mask and the FA floor stay stops**: the confidence is NaN outside
tissue, which `advanceFrame` takes for no data, and the connection is not finite
there. A seed in the zone has no trusted node and coasts straight, along the axis
of the filled field described next; where that axis ties, the seed is skipped as
`no_direction`.

**The trusted nodes must not read the zone.** The held connection comes from the
nodes just before it, and these are exactly where the voxel pipeline would mix the
degenerate voxels in: the spline kernel is global and the jump of the axis at the
edge rings over several voxels, the nodal stencil reaches two voxels, and the
derivative of a jump is of the order of 1/voxel against a curvature of 0.05. With a
scheme on, therefore, the principal axis of the voxels of the zone is first replaced
by the one of the mean dyadic of the trusted voxels around them, from the outside
in; a voxel outside the zone takes its differences from voxels outside it (central
or one-sided); and the nodal tensor is continued into the zone in the same way. The
last two apply to the nodal connection, the default; with `mmf.connection:
derivative` only the axis is filled. Without these the phantoms below read 0.3 to 5
voxels instead of 0.02 to 0.35, and the march gains nothing over the hold.

**The same bundle never crosses.** In trusted data two tracts of one bundle cannot
cross, being integral curves of one field; two coasts can, since each follows its own
held or marched connection. In the suite's crossing runs MMF-4h and MMF-4s crossed
373 and 370 pairs of same-family tracts, nearly all in the core, against 13 for the
adaptive mode and none for RK4. `advanceFrame` therefore takes `'Bundle'`, the label
of each seed's bundle, and stops both tracts of a bundle on a step on which they
would cross while one of them coasts (reason `cross`). The tracker does not pass it
yet: its seeds carry no bundle labels, and a bundle would have to be defined at the
entrance of each zone (the tracts entering from the same side, close together, with
headings within `exit_tol`).

**The march in 3-D.** The planar march of the suite's notes, `dsigma/ds = -sigma^2 +
D`, `dkappa/ds = -kappa sigma`, has a matrix form. In the carried frame, whose
normal plane is rotation-minimizing along the coast, the curvature vector
`a = (Omega_12, Omega_13)` and the transverse block `B_pq = e_p . d e1/d e_q`
(p, q = 2, 3) obey

```
dB/ds = -B^2 + D          the flatness of the connection (the trace is Raychaudhuri's equation)
da/ds = -B' a             the curvature 1-form is closed, as div(kappa e1) = 0 says in the plane
```

and the tracker takes `B_0` and `a_0` as the means over the ring and `D = 0`, so
that with `U = I + (s - s_0) B_0`, `B = B_0 U^-1`, `a = U^-T a_0` and
`u = det U`, which is `1 + sigma_0 (s - s_0)` on a planar field; the coast's density
is `1/u` of the ring's and `u < 0.1` is the stop `focus`. The matrix is not a
luxury: a closure on the trace of `B` would let a fan-out *across* the curvature
change the curvature (circles in the half-planes through an axis have `B = diag(0,
cot(theta)/r)` and constant curvature). On a helix `B_0` turns `a` by the torsion to
first order, which is all that a principal axis says about the torsion. The ring is
kept in the carried frame, so a re-anchoring, which replaces the normal plane by the
data's, turns the ring with it. The derivation and the checks on analytic fields
are in the suite (`doc/MMFtractography_notes.tex`, section on the hold scheme;
`testadvanceFrame` section 11).

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
`outside_mask`, `outside`) with every integrator, `angle_max` on both
the step and the data's heading (a 90-degree interface between bundles with
the stages straddling the swap plane, and the crossing fixture at eight
sub-voxel seed phases: every mode stops with `angle` and no track has a kink
above `angle_max`), the safety ceiling when `angle_max` is disabled, the
one-sided stencil at a volume face, a non-finite voxel in a spline field, the
per-step trace and the rejection of `rkf45` and `csd`.

For `mmf.hold` it runs the phantoms with the voxels of a zone turned into pancakes,
CL 0 and the axis 90 degrees off, with the presets' settings. Without the scheme
every one of these tracts stops at the zone (`angle`):

| field | zone | seeds | max distance from the true curve (voxels): hold / march |
|---|---|---|---|
| circles | sector of 20 degrees, 4 to 8 voxels of arc | 32 phases at 4 radii | 0.04 median, 0.10 worst / the same |
| log spirals, 60 degrees, outwards | annulus of 3 voxels, 6 of arc | 24 phases at 3 radii | 0.21, 0.25 / 0.02, 0.05 |
| the same, inwards | | 24 phases at 3 radii | 0.09, 0.17 / 0.02, 0.07 |
| helix of the fixture | sector of 40 degrees, 11 of arc | 1, at the end of the coast | 0.35 / 0.04 |
| crossing of the fixture | the 7 voxels of band B | 8 phases | 0, straight to the far face / the same |

Without the filling of the zone and the one-sided stencils of the section above the
circles read 1.0 and 1.5 voxels, the spirals outwards 1.2 and 1.8 for the hold and
the march alike, and the helix 2.6. The test also pins that the scheme changes
nothing without a zone; `exit` (band A turned by 30 degrees beyond the crossing, with
`exit_tol` 0.2, and no stop with `exit_tol` above that); `coast`; `focus` (rays through
a point in 3-D, at r/r_0 = 0.32 where the cross-section is a tenth; the hold has no
fan and runs on to the FA floor); that a tract beside low-CL voxels below the FA floor
or outside the mask does not coast; the ring across a re-anchoring (circles in a tilted
plane, two sectors fewer trusted nodes apart than the ring holds: the hold stays in
the plane to 1e-6 voxel and the march to 0.04, against 0.4 to 0.6 with the ring left
behind); the FA floor and the mask inside the zone, with `exit_delay` as well; the
share of the zone and its warning at 40%; that `hold_min` and `hold_k` reach the
scheme; the signs of the axial data; the trace and the rejected inputs.

On the sample dataset (`data/parcellation_sample/sample.mat`, 128 x 128 x 50,
`config/mmf_dti.yml` at `seeding.density: 1`, 147k seeds of which 100k lie in
tissue) the tracker produces 69.7k tracks of at least 15 voxels in about 11
minutes, every node inside the mask and no node-to-node turn above `angle_max`;
`hinec_dti.yml` on the same seeds gives 74k in about 10 minutes, and for the
68.4k seeds both keep the two tracks agree to a median nearest-node distance of
0.25 voxel. The half-tracks end on the FA floor (77%), the volume faces (10%,
the slab is cut at both ends), `angle_max` (9%, against 0.5% for hinec,
for the reason given above), the mask (3%) and a missing axis (1%).

<a id="angle-max"></a>
### angle_max is not needed with the scheme

With the scheme on, the tract transports its frame wherever the data are degenerate,
and these are the places where the data's axis turns sharply. `angle_max` then
has little left to do, and the scheme's own exit test does that. Of the half-tracks
of `mmf_dti_hold` on the ISMRM 2015 phantom with the `angle_max` of `mmf_dti` (45 degrees
per step), 14.7% stopped on it. In a trace of 8000 seeds, 86% of those stops were steps
leaving the zone whose drift (1.1 rad median) the exit test would have stopped anyway,
since the heading test of `angle_max` runs first; 3.5% had the zone within a voxel
ahead, the axis swinging before the eroded confidence falls; 10% were in trusted data
with trusted data ahead. The rotation of a step never reached the limit from `angle_max` (0.17 rad at
most). With `angle_max` off those stops become exits, and the scores do not move:

| `angle_max` | streamlines | valid | wrong path | overlap | overreach | mean F1 | bundles | mean length |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| `mmf_dti_hold`, 225 | 48487 | 53.9% | 6.1% | 0.483 | 0.337 | 0.497 | 23 | 45.6 |
| `mmf_dti_hold`, 0 (the preset) | 48588 | 53.5% | 6.8% | 0.483 | 0.339 | 0.497 | 23 | 46.8 |
| `mmf_dti_march`, 225 | 47874 | 53.5% | 5.8% | 0.476 | 0.325 | 0.495 | 23 | 45.3 |
| `mmf_dti_march`, 0 (the preset) | 47979 | 53.2% | 6.4% | 0.477 | 0.329 | 0.495 | 23 | 46.5 |

So the presets run with `angle_max: 0`. Without the scheme it is a different matter: the
adaptive mode does not transport through degenerate data, it re-anchors onto them,
and `angle_max` is the only thing that stops that snap. `mmf_dti` with `reanchor_tol`
0.2 and `angle_max: 0` gives 57634 streamlines of mean length 58.5, 53.2% valid and
13.5% wrong path against 11.2%, F1 0.483 against 0.492. `mmf_dti` keeps its `angle_max`.

<a id="on-data"></a>
### MMF-4h and MMF-4s on data

The ISMRM 2015 phantom has a ground truth. These rows come from the scorer of
`bin/run_ismrm_scoring.sh` (scilpy 2.3, `scil_tractogram_segment_with_ROI_and_score`;
point `ISMRM_VENV` at an environment that has it), with `--save_wpc_separately` for the
wrong-path column. A streamline can count for several bundles there, so valid and wrong path are
unique streamlines as a share of all of them, a wrong-path one being valid for no
bundle. All seeds with FA at least 0.1 in the mask, `termination.fa_min: 0.08`,
`min_arc: 15`, `max_arc: 400`:

| config | streamlines | valid | wrong path | overlap | overreach | mean F1 | bundles | mean length |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| `mmf_dti` (`reanchor_tol` 1e-3) | 56823 | 58.2% | 10.3% | 0.391 | 0.226 | 0.436 | 22 | 55.9 |
| `mmf_dti` with `reanchor_tol` 0.2 | 56703 | 54.5% | 11.2% | 0.509 | 0.405 | 0.492 | 23 | 54.0 |
| `mmf_dti_hold` (with `angle_max` 225) | 48487 | 53.9% | 6.1% | 0.483 | 0.337 | 0.497 | 23 | 45.6 |
| `mmf_dti_march` (with `angle_max` 225) | 47874 | 53.5% | 5.8% | 0.476 | 0.325 | 0.495 | 23 | 45.3 |

`reanchor_tol: 0.2` alone (MMF-4a in the notes) lifts F1 by far more than the scheme
does, by trusting the smoother transported heading between the voxels, at the cost of
more overreach (these two rows and the next were run with the `angle_max` of 225 that the
presets have since dropped, see above). On top of it the hold halves the wrong-path streamlines (6367 to 2973,
a share down by 45%) and cuts the overreach by a sixth, for 14% fewer streamlines, 16%
shorter ones and 5% less overlap, so F1 hardly moves. It is a trade of recall for
precision rather than a gain. The march adds nothing to the hold here: the fan-out of
the last trusted nodes is a derivative of noisy data, and its gain on the phantoms
above is for clean curved bundles. The settings were chosen on half the seeds, by a
local reimplementation of the scorer's rules that read F1 about 0.01 lower and kept the
order: `hold_min` 0.15 to 0.2, `hold_guard` 0.5, `exit_tol` 0.5 to 0.8 and `coast_max`
12 to 20 were the best of those tried, 0.15 with the higher F1 and 0.2 with fewer wrong
paths; `hold_min` 0.25 and `hold_guard` 1 cost overlap. `hold_min` 0.2 is 39% of this
phantom's tissue voxels.

On the sample dataset there is no ground truth, and its anisotropy is low: the median
CL over the tissue is 0.2, so the default floor marks 51% of the tissue voxels and the
tracker warns. With a quarter of the seeds, `mmf_dti` with `reanchor_tol` 0.2 gives
17745 tracks of at least 15 voxels and a mean of 44.1; `mmf_dti_hold` gives 11900 and
32.5 with 49% of the steps coasted, and with `hold_min: 0.1`, a zone of 12.5%, 14937
and 40.0 with 10%. The fewer tracks are not a verdict either way, as the phantom shows,
but a floor that coasts half of the arc is not detecting crossings. Lower it there.

## Files

- `src/nim_tractography/nim_tractography_mmf.m` — the tracker (seeds, field,
  connection handles, termination, metadata)
- `src/nim_tractography/advanceFrame.m` — the scheme (identical to the suite's copy)
- `config/mmf_dti.yml` — tracker config: `algorithm: mmf`, CF4, spline kernel,
  adaptive mode
- `config/mmf_dti_hold.yml`, `config/mmf_dti_march.yml` — the same with MMF-4h and
  MMF-4s
- `tests/unit/TestMmfTracker.m`
