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
numbers below). `config/mmf_dti_march.yml` therefore uses the spline kernel. The
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
The default tolerance is `0.2` rad, the one of MMF-4s (see [the
scores](#on-data) for what it does on data).

## Termination

The connection is returned non-finite wherever a point is outside the volume,
outside the propagation mask (or the brain mask), below the FA floor
`termination.fa_min` or without a well-conditioned principal axis.

**Since 2026-11-16 a non-finite connection does not stop a tract.** `advanceFrame`
takes it as zero at that stage of the step and takes the step (line 945 of
`advanceFrame.m`, `maxRot(bad) = NaN;`, is switched off, so that the file is the
moving-frames test suite's). A tract that leaves the mask, drops below the FA floor
or meets a point without a principal axis runs on **straight** along its carried
`e1`, through non-tissue, until it re-enters tissue (where the data's axis is read
again, and the adaptive mode re-anchors onto whatever bundle it meets), leaves the
volume (`outside`) or reaches `max_arc`. The tracker drops only a last node that
lies outside tissue, so tracks contain these straight runs outside the brain.
`fa`, `outside_mask` and `outside_domain` no longer occur as stop reasons. With
`angle_max: 0` (the preset of MMF-4s) nothing limits the heading at
re-entry. Before, `advanceFrame` stopped a tract **before** a step at whose stages
the connection was not finite, and every node of a track was an integration node
inside tissue; reverting the commit that switched the line off (the line, the
tests and these notes) restores it. `advanceFrame.m`'s own help still describes the
`nan` stop, as the suite's copy does. Leaving the volume stops the tract, and the
cut boundary point is dropped.

**Since 2026-11-27 the mask stops a tract again** (`mmf.mask_stop`, default `true`):
each half-track ends at its last node inside the propagation mask (or the brain
mask), termination `outside_mask`. The tracking is causal, so this is the start of
the track that runs on without the key, and the straight runs (or coasts) outside
the brain are gone; inside the mask the FA floor and points without an axis are
still no stop. `mask_stop: false` restores the tracker of 2026-11-16 to 2026-11-26,
whose tracts ended at the face of the volume (ISMRM 2015 on 2026-11-26: MMF-4a F1 0,
99.9 % of its half-tracks leaving the volume).

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
see below; since 2026-11-16 `fa` and `outside_domain` do not occur, and
`outside_mask` only with `mmf.mask_stop`, the default since 2026-11-27, see above). Seeds outside tissue or without a
direction are not tracked and are counted by reason in
`meta.mmf.seeds_skipped`, and so, with `mmf.hold` and `mmf.seed_trusted: true` (the
default; `false` on 2026-11-25 only), are the seeds in the zone of the scheme
(`zone`, see below). ACT masks are
not used by this tracker.

With `debug.trace` the per-step records of the traced seeds have the layout of
the hinec tracker's (`pos`, `dir` = the carried `e1`, `fa`, `h`, `chord`,
`turn` in degrees, `flip` always false, `n_steps`, `termination`), so
`nim_track_diagnose` reads them, plus `drift` (angle to the data's `e1`),
`reanchored`, `kappa`, `s`, `coasted` and `confidence`.

## Degenerate data: MMF-4s (`mmf.hold`)

In a voxel where two bundles cross, the single tensor is a pancake and its
principal axis means nothing. A tract that reads the connection there turns into
the other bundle or, more often, stops on `angle_max`. The hold scheme of
`advanceFrame` lets it cross instead: inside a zone of degenerate data it does not
read the data at all but **coasts**.

With `mmf.hold: march`, the default, this is **MMF-4s**: in the zone the tract follows
the curvature and the fan-out of its last `hold_k` trusted nodes, marched by the
flatness of the connection. That is exact on circles, logarithmic spirals and rays, and
turns the curvature by the torsion on a helix. With it comes `reanchor_tol: 0.2`: a
config that only selects `algorithm: mmf` runs MMF-4s, with whatever `angle_max` it has
(the schema default is 225). The preset `config/mmf_dti_march.yml` sets
`angle_max: 0` (see [below](#angle-max) why); see [the scores](#on-data) below.
(`mmf.hold` also accepts `none`, no coasting, and `hold`, the connection held constant
along the coast. They are not presented here.)

The confidence is the linear shape `CL = (l1 - l2)/l1` of the voxel
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

**The confidence from the number of CSD peaks** (`mmf.confidence: csd_peaks`,
2026-11-27; `cl` stays the default). On FiberCup (3 mm, b 1500, median FA about
0.15) CL does not single out the X crossing: its median is 0.18 inside the X and
0.15 in the rest of the phantom, so no `hold_min` makes the zone the crossing. With
`csd_peaks` the confidence is 1 in a voxel with one CSD peak, 0 in a voxel with two
or more and NaN in a voxel with none, interpolated and eroded as CL is; a peak counts
only when its weight is at least `mmf.peak_ratio` times the largest of the voxel (0,
the default, counts every peak CSD kept, `tractography.csd.peak_thresh` and
`max_peaks`). `nim_field` provisions the peaks (the same cached sidecar as
`field: csd`) while the tracked field stays dti. The confidence takes the values 0
and 1 at the voxels, so any `hold_min` in (0, 1] gives the same zone at the voxels
(the interpolation between them still depends on it).

**Leaving it.** In the `adaptive` and `reanchor` modes, a tract reaching a trusted
node is always re-anchored on the data, whatever its heading. The exit test of
`exit_tol`, which stopped a tract whose heading was more than `exit_tol` off the data
there (`exit`, the coast having left the bundle), is obsolete since 2026-11-23: the
key is still read, so that existing configs load, but ignored whatever its value
(the tracker warns, `mmf:obsoleteExitTol`, when it is finite), and `exit` no longer
occurs. The rest of this paragraph keeps the history of the test. Since 2026-11-24
non-tissue (outside the mask, below the FA floor) counts as the zone: a coast that
leaves the zone into non-tissue goes on coasting through it (its steps count
towards `coast_max`), and the tract is re-anchored where it reaches trusted
tissue; only `angle_max` limits that snap, nothing with the preset's
`angle_max: 0`. (From 2026-11-16 until then such a coast was not re-anchored in
non-tissue, ran on straight there, and its `coast_max` count restarted.)
Since 2026-11-17 `exit_tol` and
`coast_max` are `Inf` by default: no exit test, a tract leaving the zone is always
re-anchored, and no cap of the coast. The scores below were measured with
`exit_tol: 0.5` and `coast_max: 12`, the defaults until then, and were not measured
again; `coast_max: 12` can still be set, but since 2026-11-23 no setting gives the
exit test of `exit_tol: 0.5` back. The heading test of
`angle_max` comes first: data at the end of the step further from the heading
at its start than `angle_max` allows stop the tract with `angle`. In `transport` mode nothing is re-anchored or tested on leaving: the tract
keeps its heading and only `angle_max` applies. `exit_tol` is apart from
`reanchor_tol` because they ask different things of the data: a test on leaving has
to sit above the disagreement of neighbouring voxels. In the tissue of the sample
dataset (FA at least 0.2) a quarter of the neighbouring voxels differ by more than
0.2 rad, and with the 0.2 rad of the suite's noise-free fields as the test 23% of its
half-tracks ended on it. `exit_delay`, which carried the coast on for that
arc beyond the edge before the test (steps that neither counted towards `coast_max`
nor reset it), is obsolete since 2026-11-23 too: accepted and ignored, with the
warning `mmf:obsoleteExitDelay` when it is not 0; the coast ends at the first trusted
node. A coast longer than `coast_max` voxels stops (`coast`), and so does a
marched coast whose fan closes to a tenth of its cross-section, the density ten
times the ring's (`focus`). `angle_max`, when it is active, applies to the coast's
own rotation.
**The mask and the FA floor are crossed like the zone** (leaving the volume still
stops the tract, `outside`): the confidence is NaN outside tissue, which
`advanceFrame` takes for no data, and since 2026-11-24 no data counts as the zone.
The tract crosses a slab of non-tissue on its buffer, held (`hold`) or marched
(`march`), its steps there count towards `coast_max`, they never enter the buffer,
and it is re-anchored where it reaches trusted tissue. From 2026-11-16 until then
it crossed such a slab straight (the connection there, not finite, being taken as
zero, see Termination), and that zero connection entered the buffer, so that a
coast starting within `hold_k` steps after such a run held a reduced curvature, and
a straight one after a run of at least `hold_k` steps; on the crossing-bundles test
fields of the moving-frames suite this put the tracts whose fibre passes through a
box without data just before the crossing onto the other family (up to 2.7 boxes
off their fibre, against 0.65 at most with the change). A seed in the zone has no trusted node
and would coast straight, along the
axis of the filled field described next. **The grace at the seed**
(`mmf.seed_grace`, voxels of arc, an option since 2026-11-25, default `0` = none
since 2026-11-26; advanceFrame `'SeedGrace'`; it matters only with
`mmf.seed_trusted: false`): each half-track reads the data for this arc from its seed whatever
the confidence. The connection is read there, the frame is re-anchored as in trusted
data, and those nodes fill the buffer (a grace node whose own confidence is below
`hold_min` gives the buffer its connection but not its transverse sigma: the data
there are the blended crossing, whose sigma made MMF-4s stop at a false focus);
beyond the grace the scheme applies by the
confidence, so a seed in the zone coasts or marches through it on its own data
rather than along a straight line. Non-tissue stays untrusted in the grace, and the
grace applies to every seed (one within `seed_grace` of the zone reads it there). On
the crossing-bundles test fields of the moving-frames suite, with the seeds on the
entry walls (also those in the zone), a grace of 0.25 box gave a mean largest MMF-4s
error of 1.55 boxes over the 15 cases and 8 of 150 tracts more than a box off their
fibre, against 3.84 and 29 without it; 0.5 gave 1.85 and 9, 1 gave 3.41 and 48: a
longer grace reads the blend of the two bundles longer and sends more seeds onto the
other one. The tracts whose fibre meets the crossing right after the wall stay off
their fibre with any grace. On 2026-11-25 only, `mmf.seed_trusted` was `false` and
`seed_grace` `0.25` by default. No variant tried in the suite (the grace, following
the data until the first trusted node, a pointwise confidence along the seed walls)
fixed the wrong curvature and heading that the blended data give at a seed in the
zone, so since 2026-11-26 the seeds stay in trusted data by default: seed where the
fibres are coherent.
**Seeds only in trusted data**
(`mmf.seed_trusted: true`, the default): such a seed is not tracked
but skipped as `zone` in `meta.mmf.seeds_skipped`, its trace record (with
`debug.trace`) a zero-step one that keeps its confidence. The test is the scheme's
own at a first node, the confidence read as at every node (CL interpolated over the
tissue and eroded by `hold_guard`) finite and below `hold_min`, so that the erosion
also takes the seeds within about `hold_guard` of the zone; seeds at the edge of
the tissue are not taken, as the edge is no zone. The key shares its aim, tracts that
start only in trusted data, and its test with the `TrustedSeeds` of the moving frames
test suite (`runtractographycrossingtest.m`), but where the suite replaces all the
seeds of such a family by as many on the trusted part of its entry wall,
`seed_trusted` drops each seed in the zone and adds none, so fewer tracks result (the
dropped seeds stay counted in `meta.n_seeds`). It differs from
`seeding.cl_min`, which drops whole voxels by their own CL before any tracker runs
(no interpolation, no erosion, absent from `meta.n_seeds`). Neither key drops the
voxels the tensor refit fills (`preprocessing.refit_zeroed_tensors`): at the noisy
edge of a mask their floored tensors are close to linear, CL near 1, and read as
trusted. With
`seed_trusted: false` and `seed_grace: 0` the seed coasts straight, as before
2026-11-16; where the axis of the filled field ties, the seed is skipped as
`no_direction` either way. `seed_trusted: true, seed_grace: 0`, the defaults, give the
tracker of 2026-11-24.

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
voxels instead of 0.02 to 0.1.

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
the brain mask, the termination reasons (`max_arc`, `angle`, `outside`) with every
integrator, and, since 2026-11-16, the straight runs through non-tissue (beyond the
mask, through an FA cliff, with no `fa` or `outside_mask` stop), `angle_max` on both
the step and the data's heading (a 90-degree interface between bundles with
the stages straddling the swap plane, and the crossing fixture at eight
sub-voxel seed phases: every mode stops with `angle` and no track has a kink
above `angle_max`), the safety ceiling when `angle_max` is disabled (a stage within
`gradient_delta` of a swap plane stops with `no_direction`; a stage exactly on it,
where the connection is not finite, is now crossed), the
one-sided stencil at a volume face, a non-finite voxel in a spline field, the
per-step trace and the rejection of `rkf45` and `csd`.

For `mmf.hold` it runs the phantoms with the voxels of a zone turned into pancakes,
CL 0 and the axis 90 degrees off, with the preset's settings. Without the scheme
every one of these tracts stops at the zone (`angle`):

| field | zone | seeds | max distance from the true curve (voxels) |
|---|---|---|---|
| circles | sector of 20 degrees, 4 to 8 voxels of arc | 32 phases at 4 radii | 0.04 median, 0.10 worst |
| log spirals, 60 degrees, outwards | annulus of 3 voxels, 6 of arc | 24 phases at 3 radii | 0.02, 0.05 |
| the same, inwards | | 24 phases at 3 radii | 0.02, 0.07 |
| helix of the fixture | sector of 40 degrees, 11 of arc | 1, at the end of the coast | 0.04 |
| crossing of the fixture | the 7 voxels of band B | 8 phases | 0, straight to the far face |

Without the filling of the zone and the one-sided stencils of the section above the
circles read 1.0 and 1.5 voxels, the spirals outwards 1.2 and 1.8, and the helix 2.6. The test also pins that the scheme changes
nothing without a zone; that `exit_tol` is ignored (band A turned by 30 degrees beyond
the crossing: the same tracts, re-anchored and going on, with `exit_tol` 0, 0.2 and 1 as
without it; until 2026-11-23 `exit_tol` 0.2 stopped them, `exit`); `coast`; `focus` (rays through
a point in 3-D, at r/r_0 = 0.32 where the cross-section is a tenth); that a tract beside
low-CL voxels below the FA floor
or outside the mask does not coast; the ring across a re-anchoring (circles in a tilted
plane, two sectors fewer trusted nodes apart than the ring holds: the march stays in
the plane to 0.04 voxel, against 0.4 to 0.6 with the ring left behind); that the FA floor and the mask inside the zone are crossed; that `exit_delay` is
ignored; the
share of the zone and its warning at 40%; that `hold_min` and `hold_k` reach the
scheme; the signs of the axial data; the trace and the rejected inputs.

<a id="angle-max"></a>
### angle_max is not needed with the scheme

With the scheme on, the tract transports its frame wherever the data are degenerate,
and these are the places where the data's axis turns sharply. `angle_max` then
has little left to do, and the scheme's own exit test did that until it was made
obsolete on 2026-11-23; since then nothing tests the direction on leaving the zone
with `angle_max: 0`. (The numbers of this
section and the next were measured before 2026-11-16, when tracts stopped at the mask
and the FA floor; since then they run on through non-tissue, see Termination, and the
scores would have to be measured again.) With `angle_max` off the stops it made become
exits, and the scores do not move:

| `angle_max` | streamlines | valid | wrong path | overlap | overreach | mean F1 | bundles | mean length |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| `mmf_dti_march`, 225 | 47874 | 53.5% | 5.8% | 0.476 | 0.325 | 0.495 | 23 | 45.3 |
| `mmf_dti_march`, 0 (the preset) | 47979 | 53.2% | 6.4% | 0.477 | 0.329 | 0.495 | 23 | 46.5 |

So the preset runs with `angle_max: 0`.

<a id="on-data"></a>
### MMF-4s on data

The ISMRM 2015 phantom has a ground truth. These rows come from the scorer of
`bin/run_ismrm_scoring.sh` (scilpy 2.3, `scil_tractogram_segment_with_ROI_and_score`;
point `ISMRM_VENV` at an environment that has it), with `--save_wpc_separately` for the
wrong-path column. The phantom has 25 ground-truth bundles, each defined by two end
regions (head and tail), a containment mask and, for some, an inclusion mask or length
bounds; the scorer sorts every streamline against those definitions and compares the
accepted ones with the true bundles voxel by voxel. The columns, here and in the
`angle_max` table above:

- **streamlines** — the tracks the run produced and kept.
- **valid** — streamlines assigned to a bundle: both ends in its two end regions,
  every point inside its containment mask, its extra conditions met. A streamline
  can qualify for several bundles, so this is unique streamlines as a share of all.
- **wrong path** — streamlines that connect the right two end regions of a bundle but
  leave its mask on the way or fail its length or inclusion rule: right endpoints,
  wrong route. Unique, and valid for no bundle. The rest, not shown, have an end in
  no end region at all (stopped in white matter, or in the wrong region).
- **overlap** (OL) — per bundle, the fraction of the true bundle's voxels that its
  valid streamlines pass through; averaged over the bundles, a bundle without a
  valid streamline counting 0. Recall of the bundle's volume.
- **overreach** (OR) — per bundle, the voxels those streamlines visit outside the
  true bundle as a fraction of the true bundle's voxel count (so it can exceed 1);
  averaged the same way. The volume added that is not there.
- **mean F1** — per bundle the Dice score `2 TP / (2 TP + FP + FN)` on voxels (TP the
  true voxels covered, FP the voxels covered outside the bundle, FN the true voxels
  missed), averaged over the bundles; a missing bundle counts 0.
- **bundles** (VB) — how many of the 25 bundles have at least one valid streamline.
- **mean length** — the mean arc length of the streamlines, voxels.

Valid and overlap are recall-like, wrong path and overreach precision-like, and F1
balances the two per bundle. All seeds with FA at least 0.1 in the mask,
`termination.fa_min: 0.08`, `min_arc: 15`, `max_arc: 400`:

| config | streamlines | valid | wrong path | overlap | overreach | mean F1 | bundles | mean length |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| `mmf_dti_march` (with `angle_max` 225) | 47874 | 53.5% | 5.8% | 0.476 | 0.325 | 0.495 | 23 | 45.3 |

This row was run with the `angle_max` of 225 that the preset has since dropped (see
above). The settings were chosen on half the seeds, by a
local reimplementation of the scorer's rules that read F1 about 0.01 lower and kept the
order: `hold_min` 0.15 to 0.2, `hold_guard` 0.5, `exit_tol` 0.5 to 0.8 (obsolete since 2026-11-23) and `coast_max`
12 to 20 were the best of those tried, 0.15 with the higher F1 and 0.2 with fewer wrong
paths; `hold_min` 0.25 and `hold_guard` 1 cost overlap. `hold_min` 0.2 is 39% of this
phantom's tissue voxels.

On the sample dataset there is no ground truth, and its anisotropy is low: the median
CL over the tissue is 0.2, so the default floor marks 51% of the tissue voxels and the
tracker warns. A floor that coasts half of the arc is not detecting crossings: lower
it there.

## Files

- `src/nim_tractography/nim_tractography_mmf.m` — the tracker (seeds, field,
  connection handles, termination, metadata)
- `src/nim_tractography/advanceFrame.m` — the scheme
- `config/mmf_dti_march.yml` — tracker config: `algorithm: mmf`, CF4, spline kernel,
  MMF-4s
- `tests/unit/TestMmfTracker.m`
