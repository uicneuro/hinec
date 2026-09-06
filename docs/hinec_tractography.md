# HINEC Tracker Implementation Notes

`src/nim_tractography/nim_tractography_hinec.m` is the interpolated streamline tracker. It
began as a copy of the FACT baseline `nim_tractography_standard.m` with the discrete
voxel-boundary machinery replaced by a numerical ODE solve, and the two still share their
seeding, propagation-mask and length-filter conventions. This page is a map of the file for
anyone modifying it; for the method, see
[High-Order Tractography](High_Order.md), and for how it is selected and configured,
[Tractography Methods](TRACTOGRAPHY_METHODS.md).

---

## What changed relative to FACT

| | `standard` (FACT) | `hinec` |
|---|---|---|
| Direction lookup | principal eigenvector of the containing voxel | interpolated dyadic \( v_1 v_1^{\mathsf T} \), or nearest FOD peak |
| Position update | intersect the current voxel's boundary and jump there | fixed step of `integrator.step` voxels, Euler / RK2 / RK4 / RKF45 |
| Direction source | DTI only | `field: dti` or `field: csd` |
| Termination | FA, angle, length | FA, angle, arc length, interpolation domain, optional ACT tissue rules |

The voxel-boundary intersection has no counterpart in `hinec`: nothing in the file computes
an exit face, because the step length is set by the integrator rather than by the geometry.

---

## Structure of a run

1. **Parse options.** Every parameter has a default set at the top of the file, so the
   tracker is callable with a bare `nim`. In the sanctioned path the values come from the
   config instead — see [Options plumbing](#options-plumbing).
2. **Build the interpolants.** The six unique components of the dyadic
   \( v_1 v_1^{\mathsf T} \), plus FA, are wrapped in `griddedInterpolant` objects with
   `'none'` extrapolation, using the kernel named by `interpolation.method`. When
   `interpolation.upsample ≠ 1` the fields are resampled onto a finer or coarser grid first,
   leaving the coordinate frame untouched. The raw signed `v1_*` interpolants are kept only
   for the CSD path's initial reference direction.
3. **CSD peaks**, if `field: csd`: `nim.peaks` / `nim.npeaks` must already be on the `nim`.
   `runTractography` provisions them before dispatch; calling the tracker directly without
   them is an error.
4. **Propagation mask.** `FA > termination.fa_min` intersected with the brain mask, dilated
   by one voxel — deliberately **not** derived from the seed mask.
5. **Seeds.** `generate_seed_points_hinec` places `seeding.density` deterministic sub-voxel
   offsets (`nim_seed_offsets`) in every seeded voxel.
6. **Track.** Each seed is tracked bidirectionally, over `parfor` when Parallel Computing
   Toolbox is present and a plain loop otherwise. The two halves are joined, and the combined
   track is kept only if its arc length reaches `termination.min_arc`.
7. **Report.** Timing breakdown, seeding layout, and the termination analysis — the counts of
   ends that fell to FA, to the angle budget, to ACT tissue rules, or to having no initial
   direction.

`runTractography` then applies the ROI filters (`nim_filter_tracks_roi`) and, when
`output.arc_step > 0`, resamples the saved streamlines to a fixed arc-length spacing
(`nim_resample_track_arc`) before writing the `.mat`.

---

## The main functions

| Function | Role |
|---|---|
| `track_fiber_hinec` | one direction of one streamline: the step loop, all termination tests, and the per-step timing |
| `advance_position` | dispatch on the integration method; returns the new position, whether the step was accepted, and the next step size |
| `rk4_integration_step` | the four-stage step |
| `rkf45_integration_step` | the seven-stage Dormand–Prince pair, returning both solutions' difference as the error estimate |
| `interpolate_direction_trilinear` | direction and FA at a sub-voxel position; the single point where the field is read |
| `check_tissue_type` | ACT classification of a position: `WM`, `GM`, `CSF`, `OUTSIDE`, `UNKNOWN` |
| `generate_seed_points_hinec` | seeded voxels × sub-voxel offsets |
| `interp_peak_traj_h`, `dominant_peak_h` | the CSD path: nearest-peak selection, then the spatial blend |

Despite its name, `interpolate_direction_trilinear` honours whichever kernel was selected;
the name predates `cubic` and `spline`.

### Direction lookup

For `field: dti` the six dyadic components are interpolated and `nim_principal_dir` returns
the principal eigenvector of the result; the sign is then chosen to point along the incoming
tangent. Interpolating the signed eigenvector components instead would average across the
arbitrary per-voxel sign of a line field.

For `field: csd` the peak nearest the incoming tangent is chosen at each contributing voxel
and the selected peaks are blended spatially, over a 2×2×2 stencil for `trilinear` or 4×4×4
for `cubic`.

Outside the interpolation domain the interpolants return `NaN`, which the caller turns into
an empty direction and the step loop treats as a termination. `cubic` and `spline` need a
wider stencil, so their in-bounds test is correspondingly tighter.

### Integration and its fallbacks

Every Runge–Kutta stage is a field lookup that can fail near the domain edge. A failed stage
reuses the previous stage's vector rather than aborting the step, and each stage is
sign-aligned to the incoming tangent before being used. For `rkf45` a rejected step returns
the original position and the loop retries with the reduced step, up to five times, after
which that end of the track terminates.

### Termination

Checked per step, in the loop: FA below `termination.fa_min`; the turn exceeding the budget
`termination.angle_max × step` (`nim_angle_limit`); leaving the propagation mask or the
interpolation domain; the derived `max_steps`; and, with ACT active, entering CSF or GM. Each
end carries a reason string, which is what the termination analysis counts.

---

## Options plumbing

The config surface is nested and validated (`nim_config_schema`), but the tracker still reads
flat option names. `nim_config_to_options` is the one place that translates between them, and
it is where the historical names survive:

| Config path | Option the tracker reads |
|---|---|
| `integrator.method` | `integration_order` (1, 2, 4, 5 — a method selector despite the name) |
| `integrator.step` | `step_size` |
| `integrator.adaptive` | `adaptive_step` (forced false for the fixed-step methods) |
| `interpolation.method` | `interp_method` |
| `termination.fa_min` | `termination_fa` |
| `termination.angle_max` | `angle_thresh` |
| `termination.min_arc` | `min_length` |
| `termination.max_arc` | `max_steps`, derived as `ceil(max_arc/step)` |
| `act` | `act_enabled` |

When the trackers are updated to read canonical names, this function shrinks and then
disappears. Until then, adding a knob means adding it to `nim_config_schema` **and** mapping
it here — never hardcoding a default in the tracker alone.

---

## Invariants worth preserving

- **The direction is a pure function of position.** No steering term may depend on the
  incoming tangent beyond choosing a sign (DTI) or selecting a peak (CSD); otherwise the ODE
  becomes \( dx/ds = v(x, dx/ds) \) and classical order theory no longer applies to it. This
  is why `sel_power` was removed.
- **Propagation domain ≠ seed mask.** Deriving one from the other silently caps track length
  by the extent of the seed region, which is invisible under whole-brain seeding and fatal
  under ROI seeding.
- **Step-invariant termination.** The angle limit is a rate per voxel of arc and the step
  ceiling is an arc length; both must stay independent of the step size, or refining the step
  changes which tracks survive.
- **Deterministic seeding.** `uniform` seeding uses no RNG, which is what makes per-streamline
  comparison across a parameter ladder possible.

---

## Tests

| Test | Covers |
|---|---|
| `tests/unit/TestDyadicInterpolation.m` | sign-invariance of the dyadic direction lookup |
| `tests/unit/TestIntegratorOrder.m` | observed convergence order per integrator on a smooth field |
| `tests/unit/TestAngleLimit.m`, `TestAngleCriterionBehaviour.m` | the turn budget, its step-invariance and its 90° ceiling |
| `tests/unit/TestRoiSeeding.m` | ROI seed masks and the propagation-mask separation |
| `tests/unit/TestTrackResample.m` | arc-length resampling of saved streamlines |
| `tests/unit/TestConfigSchema.m`, `TestConfigYaml.m` | the schema, its defaults, migration and validation |

```matlab
results = runtests('tests/unit/TestIntegratorOrder.m');
```
