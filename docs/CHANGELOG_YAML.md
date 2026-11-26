# YAML Configuration Changelog

## 2026-11-26 - `mmf.seed_grace`: an optional grace at the seed; seeds stay in trusted data

- New key `tractography.mmf.seed_grace` (numeric, voxels of arc, default `0` = none;
  flat option name `mmf_seed_grace`), read by `algorithm: mmf` with `mmf.hold` `hold`
  or `march`: each half-track reads the data for this arc from its seed whatever the
  confidence (advanceFrame's new option `'SeedGrace'`), so that the nodes there fill
  the buffer, and the scheme applies by the confidence beyond it. A grace node whose
  own confidence is below `hold_min` gives the buffer its connection but no sigma.
  Non-tissue stays untrusted. `meta.mmf` gains `seed_grace`. It matters mainly with
  `mmf.seed_trusted: false`, where it lets a seed in the zone start from its own
  data instead of coasting straight.
- `tractography.mmf.seed_trusted` stays `true` by default: a seed in the zone is not
  tracked. On the moving-frames test suite's crossing cases no remedy for seeds in
  the zone (a grace of 0.25 to 1.5 box, following the data until the first trusted
  node, a pointwise confidence along the seed walls) matched seeding in trusted data,
  since the blended data at such a seed give a wrong curvature and heading. The
  suite's crossing drivers re-seed in trusted data by default again (2026-11-26).
- With the defaults the tracker is that of 2026-11-24: no behaviour change.
- The ISMRM 2015 scores of `docs/MMF_TRACTOGRAPHY.md` were measured before the
  changes of 2026-11-16 to 2026-11-24 and have not been re-scored.

## 2026-11-24 - `algorithm: mmf` with `mmf.hold`: non-tissue counts as the zone

- No new key. In `src/nim_tractography/advanceFrame.m` a node or stage whose
  confidence is not finite (no data: outside the volume or the mask, below the FA
  floor) counts as the zone of the hold scheme, as one below `mmf.hold_min` does.
- **Behaviour change** for runs with `mmf.hold: hold` or `march` (MMF-4h, MMF-4s):
  a tract crosses non-tissue on its buffer (held or marched) instead of running
  straight, non-tissue nodes no longer enter the buffer with a zero connection, its
  steps there count towards `mmf.coast_max` and are reported as coasted, and it is
  re-anchored where it reaches trusted tissue. Runs without `mmf.hold` are
  unchanged. The same change is made in the moving-frames test suite's copy of
  `advanceFrame.m`, where it removed the deviations of up to 2.7 boxes of the
  tracts whose fibre crosses a box without data just before the crossing.

## 2026-11-23 - `mmf.exit_delay` is obsolete too: no coast beyond the zone

- `tractography.mmf.exit_delay` (flat `mmf_exit_delay`) is still a key, so that
  existing configs load, but it has no effect: whatever its value the coast ends at
  the first trusted node beyond the zone, where the tract is re-anchored, the
  behaviour of the default `0`. `nim_tractography_mmf` warns
  (`mmf:obsoleteExitDelay`) when it is not 0; `meta.mmf` no longer has `exit_delay`.
  A config that sets it changes meaning.
- `advanceFrame.m`: the option `'ExitDelay'` is accepted and ignored, with the warning
  `advanceFrame:obsoleteExitDelay`; `info.hold` no longer has `exitDelay`. The code
  of the delay is commented out, not deleted. The same change is made in the
  moving-frames test suite's copy of `advanceFrame.m`.

## 2026-11-23 - `mmf.exit_tol` is obsolete: no exit test, whatever its value

- `tractography.mmf.exit_tol` (flat `mmf_exit_tol`) is still a key, so that existing
  configs load, but it has no effect: whatever its value, a tract leaving the zone is
  always re-anchored on the data, the behaviour of the default `Inf` since
  2026-11-17, and the stop reason `exit` no longer occurs. `nim_tractography_mmf`
  warns (`mmf:obsoleteExitTol`) when it is finite. `meta.mmf` no longer has
  `exit_tol`.
- **A config that sets a finite `exit_tol` changes meaning**: until now it stopped
  tracts leaving the zone more than `exit_tol` off the data; now they are re-anchored
  and go on. With `termination.angle_max: 0`, as in `config/mmf_dti_march.yml`,
  nothing tests the direction on leaving the zone; set `angle_max` for a test. The
  ISMRM 2015 scores of [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md), measured with
  `exit_tol: 0.5`, can no longer be reproduced by a setting.
- `advanceFrame.m`: the option `'ExitTol'` is accepted and ignored, with the warning
  `advanceFrame:obsoleteExitTol`; `info.hold` no longer has `exitTol`. The code of the
  test is commented out, not deleted. (`'ExitDelay'` was kept by this change; see the
  entry above.) The same change is made
  in the moving-frames test suite's copy of `advanceFrame.m`.

## 2026-11-17 - MMF-4s is the moving-frames scheme presented

- `config/mmf_dti.yml` (MMF-4a) and `config/mmf_dti_hold.yml` (MMF-4h) are no longer
  shipped configs: they are test fixtures now, `tests/fixtures/mmf_dti_4a.yml` and
  `tests/fixtures/mmf_dti_4h.yml`. `./bin/run_tractography.sh mmf_dti` and
  `mmf_dti_hold` therefore no longer find a config; `config/mmf_dti_march.yml`
  (MMF-4s) is the moving-frames preset.
- No key changed: `mmf.hold` still accepts `none` and `hold`, and a config that sets
  them runs as before. The docs present MMF-4s only: the descriptions, the score rows
  and the comparisons of MMF-4a and MMF-4h are removed from
  [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md) and the config tables. The entries below
  that name the two schemes or their configs describe the code as it was; the last
  commit with them in the docs is `db45c54`.

## 2026-11-17 - advanceFrame: 'BundleRule' and 'BundleLaw' removed

- No key changed: the tracker never passed them. `src/nim_tractography/advanceFrame.m`
  no longer has the options `'BundleRule'` and `'BundleLaw'` (2026-11-03 and
  2026-11-04), the stop reason `cross`, or `info.hold.bundleRule` / `bundleLaw`; a call
  with either option is rejected as an unknown option. Everything else in the file is
  as before, the `'ExitTol'` default of `Inf` and the non-finite connection included.
- `advanceFrame.m` is therefore no longer the moving-frames test suite's version byte
  for byte: the suite keeps the two options. The last commit with them is `dd5dee3`.

## 2026-11-17 - MMF-ml (multilayer moving frames) removed

- Removed: the MMF-ml scheme of 2026-11-07 and everything that existed only for it.
  `field: lift` is no longer a value of `tractography.field` (it is `dti | csd`
  again), for `algorithm: mmf` and for `algorithm: hinec`, where the lifts were a
  peak field. The whole group `tractography.lift.*` is gone, and so are
  `tractography.mmf.veto_arc` and `tractography.mmf.veto_tol` (the history veto,
  which only the MMF-ml copy of the scheme had). A config that sets any of them is
  rejected as an unknown key or value.
- Files removed: `src/nim_tractography/advanceFrame_ml.m`,
  `nim_tractography_mmf_ml.m`, `nim_lift_connection.m`,
  `src/nim_calculation/nim_lift_field.m`, `config/mmf_ml.yml`, `docs/MMF_ML.md`,
  `scripts/lift_field_check.py`, `tests/fixtures/make_lift_crossing_nim.m` and the
  tests `TestAdvanceFrameMl`, `TestLiftField`, `TestMmfHistoryVeto`.
  `runTractography` sends every `algorithm: mmf` run to `nim_tractography_mmf`.
- Kept: `seeding.cl_min` and `seeding.coherence_min`, which came with the same commit
  but are not part of the scheme. Cached lifts beside a dataset
  (`<nim>_lift_<hash>.mat`) are not deleted and are no longer read.
- The entries below that mention MMF-ml, `field: lift`, the veto or the `_ml` files
  describe the code as it was. The last commit with MMF-ml is `9c81caf`.

## 2026-11-17 - `mmf.exit_tol` and `mmf.coast_max` are `Inf` by default

- The defaults of `tractography.mmf.exit_tol` and `tractography.mmf.coast_max` are now
  `Inf` (they were `0.5` rad and `12` voxels of arc): no exit test and no cap of the
  coast. A tract leaving the degenerate zone is always re-anchored on the data, and a
  tract still in the zone coasts on until it leaves the volume or reaches `max_arc`. The tracker's own defaults (`nim_tractography_mmf`) follow. The ranges of the two keys are `[0, Inf]`.
- **Behaviour change** for every config with `mmf.hold: hold` or `march` (the default)
  that sets neither key: the stops `exit` and `coast` no longer occur. With
  `termination.angle_max: 0`, as in `config/mmf_dti_hold.yml` and
  `config/mmf_dti_march.yml`, nothing then tests the direction of a tract leaving the
  zone. Set `exit_tol: 0.5` and `coast_max: 12` for the former behaviour.
- The scores of [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md)
  were measured with the former defaults and were not measured again.
- A numeric key whose range has `Inf` as a bound takes `inf`, `Inf` or `.inf` in a
  config file and with `--set`; any other numeric key still rejects it.
- `advanceFrame.m` itself now defaults `'ExitTol'` to `Inf`
  (it was `'ReanchorTol'`), as the moving-frames test suite's version does. No run
  through the tracker changes: it always passes `mmf.exit_tol`. Only a direct call
  of `advanceFrame` with `'Confidence'` and no `'ExitTol'` loses its exit test.

## 2026-11-16 - `algorithm: mmf`: a non-finite connection no longer stops a tract

- No new key. `src/nim_tractography/advanceFrame.m` takes the moving-frames test
  suite's version, which differs from the one before in one line: `maxRot(bad) = NaN;`
  in `lieStep` is switched off. A non-finite connection at a stage of a step (outside
  the volume or the mask, below the FA floor, no principal axis) is taken as zero and
  the step is taken.
- **Behaviour change** for every run of `algorithm: mmf` (`nim_tractography_mmf`; not
  `nim_tractography_mmf_ml`, whose `advanceFrame_ml.m` keeps the line): a tract that
  leaves the tissue is no longer stopped there but runs on straight through non-tissue
  until it re-enters tissue, leaves the volume (`outside`) or reaches `max_arc`, and
  the track keeps those nodes. `fa`, `outside_mask` and `outside_domain` no longer
  occur in `meta.termination`. With `mmf.hold`, a coast that leaves the zone into
  non-tissue is not exit-tested and its `coast_max` count restarts, and the zero
  connection of such a run enters the hold ring. `fa`, `mask` and
  `termination.fa_min` no longer stop this tracker (their descriptions say so).
  Reverting the commit (the line, the tests and the notes) restores the tracker; the
  numbers of MMF_TRACTOGRAPHY.md that were measured with the stop are marked so.

## 2026-11-16 - `preprocessing.refit_zeroed_tensors`: the voxels the fit left at zero get a tensor

- New key `preprocessing.refit_zeroed_tensors` (logical, default `true`). The voxels
  of the mask that `nim_dt_spd` left at zero (no tensor, FA 0: no tissue to every
  tracker) get their log-linear least-squares tensor with the eigenvalues raised to
  1e-7 mm²/s (`src/nim_calculation/nim_dt_refit_zeroed.m`). Nearly all of them are the
  voxels whose least-squares tensor has a negative eigenvalue: `nim_dt_spd` sends them
  to BFGS (`lib/bfgs/vox_dt_bfgs.m`), whose tensor comes back 1e6 or more too large
  (it multiplies by 1000 where it should divide, and its gradient reads the wrong
  gradient directions) and is rejected by the bound on the eigenvalues. A voxel whose
  signals do not determine a tensor stays at zero. `vox_dt_bfgs.m` is not changed.
- Applied by `main` after `nim_fa` and again before the save (the voxels the improved
  mask adds), by `runTractography` to the nim it loads, so that a dataset made before
  reaches the trackers without holes (in memory: the file is not rewritten), and by
  `cc_extraction/extract_cc_tracks.m` after its `nim_fa`.
- **Behaviour change** for every pipeline and tractography run: the voxels that were
  holes now carry a tensor where their signals allow it. On the TrackerYoga phantom all
  352 holes (12% of the mask, b = 1000; 1697 of 1698 in the noisy set), with the
  tensors of the suite's refit to the last bit. On the ISMRM 2015 phantom 2681 of its
  8562 holes; the other 5881 stay at zero (5878 have no eigenvalue above 1e-7, 3 have
  b0 <= 0). 2196 of the 2681 lie at the edge of the mask, and they come out close to
  linear (FA median 0.98, 90th percentile 0.999, against a median of 0.07 for the
  voxels of the mask with a tensor): in least squares their third eigenvalue is
  negative in all of them and their second in 1354, and both are floored. These
  voxels now pass every FA threshold: seeding (`seeding.fa_min`), the FA stop, the
  mask improvement and the segmentation of `main`, and with `mmf.hold` their CL near 1
  makes them trusted (and seeded under `mmf.seed_trusted`); `seeding.cl_min` does not
  drop them either. The segmentation runs before the second call, so the voxels the
  improved mask adds are segmented with FA 0.
- `refit_zeroed_tensors: false` leaves them all at zero, as before. `main` reads it
  from the pipeline config and `runTractography` from the tracking config (a call
  without a config always refits): set it in both.

## 2026-11-16 - `mmf.seed_trusted`: seeds only in trusted data

- New key `tractography.mmf.seed_trusted` (logical, default `true`; flat option name
  `mmf_seed_trusted`), read by `algorithm: mmf` with `mmf.hold` `hold` or `march`: a
  seed whose confidence, read as the scheme reads it at a node (CL interpolated over
  the tissue and eroded by `mmf.hold_guard`), is below `mmf.hold_min` is not tracked.
  Such a seed has no trusted node and would coast straight along the filled axis.
  It is counted in `meta.mmf.seeds_skipped` under the new reason `zone`, its
  `debug.trace` record keeps its confidence, and `meta.mmf` gains `seed_trusted`.
  See [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md).
- **Behaviour change** for every run of `nim_tractography_mmf` with the scheme on, the
  schema default (MMF-4s) and `config/mmf_dti_hold.yml`, `config/mmf_dti_march.yml`
  included: the seeds in the zone, and within about `hold_guard` of it, are no longer
  tracked. `seed_trusted: false` restores the tracker as it was. `hold: none`
  (`config/mmf_dti.yml`) is unaffected, and so are the runs of
  `nim_tractography_mmf_ml` (`field: lift`, as `config/mmf_ml.yml`, or
  `mmf.veto_arc > 0`), which do not read the key.

## 2026-11-10 - `algorithm: mmf` is MMF-4s by default

- The defaults of `tractography.mmf.hold` and `tractography.mmf.reanchor_tol` are now
  `march` and `0.2` rad (they were `none` and `1e-3`): a config that only selects
  `algorithm: mmf`, or `--set algorithm=mmf` on the config of another tracker, runs
  MMF-4s. The trackers' own defaults (`nim_tractography_mmf`,
  `nim_tractography_mmf_ml`, which `runTractography(data, 'mmf')` reads) follow. The
  march reads `nim.eval`.
- **Behaviour change** for a config that selects `algorithm: mmf` and sets neither
  key. `config/mmf_dti.yml` now sets `hold: none` beside its `reanchor_tol: 0.001`
  and is the tracker it was (MMF-4a); `config/mmf_dti_hold.yml`,
  `config/mmf_dti_march.yml` and `config/mmf_ml.yml` set both keys and are unchanged.
  `termination.angle_max` is not part of this: the presets of MMF-4h and MMF-4s
  switch it off, a config that only selects the algorithm keeps its own.

## 2026-11-04 - advanceFrame: 'BundleRule' and 'BundleLaw'

- `src/nim_tractography/advanceFrame.m` takes the suite's version with the
  options `'BundleRule'` (labels of the seeds' bundles; two tracts of one
  bundle stop with reason `cross` on a step on which they would cross while
  one of them coasts; the option added as `'Bundle'` the day before) and
  `'BundleLaw'` (the divergence term of the march from the neighbouring
  tracts of the bundle; planar connection only). The tracker passes neither:
  no config key, no change in behaviour. See
  [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md).

## 2026-11-03 - the word "budget" is gone

- Every mention of the angle criterion as a "budget" in the help texts, docs,
  configs and tests now names the parameter: `termination.angle_max` (the limit
  per step being `min(angle_max * step, 90)` degrees). No key changed. Test
  methods were renamed accordingly (for example `angleMaxStopsASharpTurn`), the
  error identifier `stitching:budget` is `stitching:maxFragments`, and the
  non-angle uses say "limit" (`stitching.max_fragments`, `stitching.neighbors`).

## 2026-11-02 - `mmf.hold`: MMF-4h and MMF-4s

- New keys of `tractography.mmf`, for the tracker `algorithm: mmf`: `hold`
  (`none | hold | march`, default `none`: nothing changes for existing configs),
  `hold_min` (0.2), `hold_k` (3), `hold_guard` (0.5 voxel), `coast_max` (12 voxels
  of arc), `exit_tol` (0.5 rad) and `exit_delay` (0 voxels of arc); flat option names
  `mmf_hold`, `mmf_hold_min`, `mmf_hold_k`, `mmf_hold_guard`, `mmf_coast_max`,
  `mmf_exit_tol`, `mmf_exit_delay`. With `hold` or `march` a tract coasts through
  degenerate data (a crossing fitted as a pancake, CL below `hold_min`), on the held
  connection of its last trusted nodes (MMF-4h) or on that curvature and fan-out
  marched by the flatness of the connection (MMF-4s). See
  [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md), with scores on the ISMRM 2015 phantom.
- New termination reasons of the tracker: `coast`, `exit`, `focus`. `meta.mmf`
  gains the settings above, `zone_fraction` and `coasted_fraction`; the `debug.trace`
  records gain `coasted` and `confidence`.
- New warning `mmf:holdZone`: more than 40% of the tissue voxels below `hold_min`.
- New tracker configs `config/mmf_dti_hold.yml` and `config/mmf_dti_march.yml`
  (`mmf_dti` with `hold` set, `reanchor_tol: 0.2` and `termination.angle_max: 0`:
  with the scheme the stops of `angle_max` are the exit test's, and the scores do not
  change without it, see [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md#angle-max)).
- `src/nim_tractography/advanceFrame.m` takes the suite's version with the whole hold
  scheme ('Confidence', the 'Hold*' options, 'CoastMax', 'ExitTol', 'ExitDelay', and the
  `info` fields `coasted`, `confidence`, `held`, `spread`), and 'HoldCurvature' 'march'
  works with three frame vectors there (the matrix law of the flat connection).

## 2026-10-20 - `algorithm: mmf` selects the moving-frames tracker

- `algorithm: mmf` is a permitted value again and now selects the new moving
  frames tracker `nim_tractography_mmf` (scheme `advanceFrame.m`: the carried
  orthonormal frame is transported by the connection 1-forms of the data's frame
  field with a closed-form exponential step, exact for a planar arc of
  constant curvature at any step size, to the order of the integrator on a
  fibre with torsion; see [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md)).
  It is unrelated to the tracker removed below; only the name is reused.
  `algorithm` takes `hinec | standard | mmf | stitching | template`.
- New group `tractography.mmf`: `mode` (`transport | reanchor | adaptive`,
  default `adaptive`), `reanchor_tol` (default `1e-3` rad) and `gradient_delta`
  (default `1e-3` voxels); flat option names `mmf_mode`, `mmf_reanchor_tol`,
  `mmf_gradient_delta`.
- For `algorithm: mmf`, `integrator.method: euler | rk2 | rk4` selects the Lie
  group integrator of order 1, 2, 4 (Lie-Euler, exponential midpoint, CF4);
  `rkf45` is rejected. `field` must be `dti` (`csd` is rejected). ACT masks are
  not used by this tracker.
- The "used by" column of the reference changes accordingly: `debug.trace` and
  `debug.trace_max` now list `mmf`; the `csd.*` keys and the rkf45 keys
  (`integrator.tolerance`, `step_min`, `step_max`) no longer do.
- For `algorithm: mmf`, `termination.angle_max` (through `nim_angle_limit`,
  `min(angle_max * step, 90)` degrees per step) bounds both the rotation of a
  frame step and the turn of the heading as the data see it (the data's `e1`
  at the end of a step against the heading at its start, in every mode), so a
  re-anchoring cannot exceed `angle_max`; the 90 degree ceiling is a real cap
  for this tracker (warning `nim_config:angleMaxCapped` instead of
  `angleMaxInert`), and `angle_max: 0` keeps a safety ceiling of 180 degrees
  per step at singular data (reported as `no_direction`).
- With `algorithm: mmf`, the `debug.trace` records have the hinec layout
  (`pos`, `dir`, `fa`, `h`, `chord`, `turn`, `flip`, `n_steps`, `termination`)
  plus `drift`, `reanchored`, `kappa`, `s`; `nim_track_diagnose` reads both.
- New tracker config `config/mmf_dti.yml` (CF4, spline kernel, adaptive mode).
  Output file `tracks_mmf_<timestamp>.mat`.

## 2026-10-20 - mmf tracker removed

- The `mmf` tracking algorithm (`algorithm: mmf`, the connection-form Method of
  Moving Frames tracer) has been removed, together with the per-run moving-frame
  geometry step it required. `algorithm: mmf` is no longer a permitted value;
  `algorithm` takes `standard | hinec | stitching | template`. A config that still
  says `algorithm: mmf` must be changed to the tracker it meant (before July 2026
  that label denoted today's `hinec`).
- `field: dwi` (frame and curvature fitted directly to the DW signal) has been
  removed. `field` takes `dti | csd`.
- The `tractography.mmf` group (`mmf.anchor`, option `mmf_anchor`) is gone.
- The stitching tracker keeps its direction mode only: `stitching.geometry` and
  `stitching.max_curvature` have been removed. `stitching.curvature_tolerance`
  and `stitching.torsion_tolerance` remain in the schema; the endpoint graph
  applies them only when a per-endpoint curvature/torsion table is handed to
  `nim_stitch_fragments`, which the pipeline no longer supplies.
- The configs `mmf_dti.yml`, `mmf_csd.yml`, `ismrm_wholebrain_mmf.yml`,
  `ismrm_wholebrain_mmf_csd.yml`, `stitching_mmf_dti.yml` and
  `stitching_mmf_csd.yml` have been deleted.
- The "MMF geometry" tier of the tract-assessment pipeline
  (`scripts/tract_assessment/geometry.py`) is unaffected: it is a per-streamline
  curvature statistic of any tractogram, not the tracker.

## Schema rework - canonical nested configuration (current)

The configuration surface is now declared once, in
`src/nim_utils/nim_config_schema.m`. Defaults, validation, unknown-key
rejection, legacy migration, the flat option names the trackers read, and the
[configuration reference](YAML_CONFIG.md) are all derived from it. Old configs
still load; each superseded key produces a deprecation warning naming its
replacement.

### Structure

- Config gained **two levels of nesting** (`section.group.key`) and inline lists.
  A third level is a parse error. The previous parser was indentation-blind and
  silently misparsed nested YAML rather than rejecting it.
- Unknown keys are now an **error**, not a silent no-op - which is how
  `rkf_tolerance` (hinec) and `rkf_tol` (mmf) came to mean the same thing under two
  names.
- Every parameter, `preprocessing.*` included, is reachable from the CLI via
  `--set <path>=<value>` on both launchers.
- [YAML_CONFIG.md](YAML_CONFIG.md) is now **generated** from
  `src/nim_utils/nim_config_schema.m` by `src/nim_utils/nim_config_docs.m`, with
  a test that fails if code and docs drift. Do not edit it by hand.

### Renamed and redefined keys

- `integration_order` -> `integrator.method`. The old key **spelled a method
  selector as a number**: `5` selected RKF45. That value was not wrong
  numerically - the implementation uses Dormand-Prince coefficients and advances
  on the 5th-order solution, keeping the embedded 4th-order one for error control
  - but a method name belongs in a method key. The values are now named: `euler | rk2 | rk4 | rkf45`.
- `max_steps` -> `termination.max_arc`. The old key counted integration steps, so
  halving the step size silently halved how far a track could travel. `max_arc` is an
  arc length in voxels; `max_steps` is derived as `ceil(max_arc / step)`.
- `min_length` -> `termination.min_arc`, an arc length in **voxels**, and it is
  now actually enforced. It had been defaulted and then never read by the hinec
  and standard trackers - a schema-validated key that every config set was
  honoured by exactly one of the three.
- `termination.angle_max` is documented, and implemented, as degrees of turning
  per **voxel of arc**, not per step: the limit for one step is
  `angle_max x step`, so refining the step does not loosen the constraint.
  Because tangents are sign-aligned (`v1` is a line field), a measured turn
  never exceeds 90 degrees, so any `angle_max` above `90/step` is *inert* rather than
  merely loose; the loader warns when a config lands there. `angle_max: 0`
  disables the criterion outright. See `src/nim_utils/nim_angle_limit.m`.
- `interpolation.method` gained `spline` alongside `trilinear` and `cubic`. The
  three differ in smoothness - C0, C1 (Keys cubic convolution, *not* a spline)
  and C2 - which caps the order a Runge-Kutta method can reach. See
  [Solution Verification](CONVERGENCE.md) for the measured orders.
- `interpolation.upsample` is implemented: the direction field is sampled on a
  grid of spacing `1/upsample` voxels before the interpolants are built. Above 1
  refines, below 1 coarsens. The coordinate frame is unchanged, so step sizes and
  lengths stay comparable across factors.
- `tractography.act` now defaults to **false**. It is read only by the hinec
  tracker.
- New filter predicates `filter.endpoints_in` (one end in each of two regions)
  and `filter.contained_in` (every point inside a corridor). Together with
  `filter.include_roi` and `filter.exclude_roi` these express the ISMRM 2015
  bundle definition.

### Retired keys

Retired keys are accepted with a warning and dropped, never silently migrated -
mapping a key that had no effect onto one that does would change behaviour.
`src/nim_utils/nim_config_retired.m` holds the list and the reason for each.

- Eight keys no tracker ever read: `gate_power`, `crossing_cp`, `curv_beta`,
  `crossing_detect`, `swing_ratio_max`, `transport_gate`, `transport_strength`,
  `bishop_eps`.
- `fa_threshold`, which was printed but never used to terminate. Use
  `termination.fa_min` to stop tracking, or `seeding.fa_min` to restrict seeding.
- `order`, a legacy compatibility key read by no tracker as an option.
- `sel_power`, **removed from the hinec tracker entirely**. It re-weighted each
  stencil voxel by `alignment^sel_power`, biasing interpolation toward the
  incoming direction. For DTI there is one principal eigenvector per voxel, so
  there was nothing to disambiguate; it simply bent a single-valued field toward
  the current heading, which is self-reinforcing, and it made the ODE
  direction-dependent (`dx/ds = v(x, dx/ds)`), so classical Runge-Kutta order
  theory no longer applied. HINEC tracking is now interpolation and integration
  only. CSD still selects the FOD peak nearest the incoming tangent, because
  reducing a multi-valued field to one direction is structurally required, but
  the alignment exponent on top of it is gone. The old MMF tracker's separate
  `mmf.frame_sel_power`, the same mechanism applied when building its moving
  frame, has since been retired as well.

---

## History

Everything below predates the schema rework and is kept as a record of how the
configuration system was introduced. It describes the **flat key names and
validation rules that the rework replaced**; for current behaviour read the
section above, or the generated [configuration reference](YAML_CONFIG.md).

### Overview

Implemented a comprehensive YAML-based parameter configuration system for the HINEC pipeline, making parameter management easy and reproducible.

### What Changed

#### New Files Created

**Configuration System**:

- `nim_utils/load_config_yaml.m` - YAML parser with validation
- `config/hinec_default.yml` - Default balanced configuration
- `config/hinec_dti_cubic.yml` - Publication-quality preset
- `config/hinec_dti_fast.yml` - Quick testing preset
- `config/irontract.yml` - IronTract challenge preset

**Documentation**:

- `docs/YAML_CONFIG.md` - Complete configuration reference (124 KB)
- `config/README.md` - Quick reference for config directory
- `test_yaml_config.m` - Validation test script

#### Modified Files

**Core Pipeline**:

- `run_hinec.sh` - Added optional config file argument
- `main.m` - Added YAML config support (backward compatible)
- `runTractography.m` - Added YAML config support (backward compatible)
- `README.md` - Added YAML quick start section

> Paths in the two lists above predate the source reorganisation: the MATLAB
> sources moved under `src/`, the launchers under `bin/`, and the test script to
> `tests/test_yaml_config.m`. `docs/YAML_CONFIG.md` is now generated rather than
> written by hand.

### Key Features

#### 1. Easy Parameter Management

**Before** (editing source code):
```matlab
% Edit runTractography.m:88-93
options.seed_density = 4;
options.step_size = 0.5;
options.termination_fa = 0.15;
options.angle_thresh = 35;
options.max_steps = 1000;
options.min_length = 35;
```

**After** (editing YAML file):
```yaml
# config/my_experiment.yml
tractography:
  seed_density: 4
  step_size: 0.5
  termination_fa: 0.15
  angle_thresh: 35
  max_steps: 1000
  min_length: 35
```

> **Superseded.** These flat key names no longer load without a deprecation
> warning. The canonical form is `tractography.seeding.density`,
> `tractography.integrator.step`, `tractography.termination.fa_min`,
> `termination.angle_max`, `termination.max_arc` and `termination.min_arc`.

#### 2. Preset Configurations

Four ready-to-use presets for common scenarios:

| Preset | Integration | Speed | Use Case |
|--------|-------------|-------|----------|
| `hinec_default.yml` | RK4 | Baseline | Standard analysis |
| `hinec_dti_cubic.yml` | RKF45 adaptive | 2-3x slower | Publications |
| `hinec_dti_fast.yml` | RK2 | 3-5x faster | Parameter testing |
| `irontract.yml` | RK4 | Baseline | IronTract challenge |

> **Superseded.** The config set has since been reorganised around one naming
> rule per family - tracker configs `<algorithm>_<field>[_<variant>].yml`,
> dataset configs `<dataset>[_variant].yml`. See `config/README.md` for the
> current list and what each preset is for.

#### 3. Automatic Validation

Parameters are validated when loaded:
```matlab
✓ Parameter validation passed
```

**Checks include**:

- Step size must be positive
- Angle threshold in (0, 180]
- Integration order in {1, 2, 4, 5}
- RKF bounds: step_min < step_max
- RKF safety factor in (0, 1]
- RKF tolerance > 0

> **Superseded.** Ranges and permitted values now come from the `range` and
> `allowed` fields of `src/nim_utils/nim_config_schema.m`, so this list cannot
> drift from the code. Two entries are no longer true as written: integration
> order is gone (`integrator.method` takes `euler | rk2 | rk4 | rkf45`), and
> `termination.angle_max` is a rate in degrees per voxel of arc with a range of
> `[0, 3600]`, where `0` means "criterion disabled".

#### 4. Backward Compatibility

All existing code continues to work:

```matlab
% Legacy usage - still works
runTractography('data.mat', 'hinec');

% New YAML usage
config = load_config_yaml('config/hinec_dti_cubic.yml');
runTractography('data.mat', config);
```

### Usage Examples

> **Superseded.** The launchers now live under `bin/`, so every `./run_hinec.sh`
> below is `./bin/run_hinec.sh` today, and the programmatic example sets flat
> field names that have since been nested. Tractography-only iteration on an
> already-preprocessed dataset goes through `bin/run_tractography.sh`, which did
> not exist yet.

#### Basic Usage

```bash
# Use default config
./run_hinec.sh data/subject processed.mat

# Use custom config
./run_hinec.sh data/subject processed.mat config/hinec_dti_cubic.yml
```

#### Create Custom Config

```bash
# Copy template
cp config/hinec_default.yml config/my_experiment.yml

# Edit parameters
nano config/my_experiment.yml

# Run with custom config
./run_hinec.sh data/subject processed.mat config/my_experiment.yml
```

#### Programmatic Usage

```matlab
% Load and modify config
config = load_config_yaml('config/hinec_default.yml');
config.tractography.seed_density = 8;
config.tractography.termination_fa = 0.05;

% Run with modified config
runTractography('data.mat', config);
```

### Parameter Organization

#### Preprocessing Section

- Denoising settings
- Motion/eddy correction
- Brain mask improvement
- Atlas selection
- T1 registration (optional)

#### Tractography Section
- Algorithm selection (standard/hinec)
- Integration method (Euler/RK2/RK4/RKF45)
- Step size control
- RKF45 adaptive parameters
- Seeding strategy
- Termination criteria
- ACT configuration

#### Output Section
- Output directory
- Filename timestamps
- Auto-visualization

> **Superseded.** There is no top-level `output:` section: output directories
> and timestamps are handled by the run-directory system, and `runTractography`
> no longer draws figures. The only output key is
> `tractography.output.arc_step`. The tractography section has also gained a
> third algorithm (`mmf`), an `interpolation` group, ROI seeding and the track
> `filter` group.

### Testing

Run validation test:

```matlab
test_yaml_config
```

> The script now lives at `tests/test_yaml_config.m`, and the schema itself is
> covered by `tests/unit/TestConfigSchema.m` and `tests/unit/TestConfigYaml.m`.
> The transcript below is the original output and no longer matches.

Expected output:
```
=== HINEC YAML Configuration System Test ===

Test 1: Loading default configuration...
  ✓ Default config loaded successfully

Test 2: Loading high precision configuration...
  ✓ High precision config loaded

Test 3: Parameter validation...
  ✓ All parameter ranges validated

Test 4: Checking preset configurations...
  ✓ All preset configs found

Test 5: RKF45 parameter validation...
  ✓ RKF45 parameters validated

Test 6: Comparing configurations...
  [comparison table]

========================================
ALL TESTS PASSED ✓
========================================
```

### Benefits

#### For Users
1. **Easy parameter changes** - Edit YAML files instead of source code
2. **Reproducibility** - Version control configurations in git
3. **Presets** - Quick access to validated parameter sets
4. **Safety** - Automatic validation prevents invalid parameters

#### For Developers
1. **Centralized configuration** - All parameters in one place
2. **Type safety** - Validation catches errors early
3. **Backward compatibility** - Legacy code still works
4. **Extensibility** - Easy to add new parameters

### Documentation

- **Complete reference**: [YAML_CONFIG.md](YAML_CONFIG.md)
- **Quick reference**: `config/README.md` in the repository
- **RKF45 guide**: [RKF_Usage.md](RKF_Usage.md)
- **Main README**: `README.md` in the repository

### Migration Guide

#### For Existing Users

No changes required. Existing usage patterns continue to work:

```matlab
% This still works exactly as before
runTractography('data.mat', 'hinec');
```

To use YAML configs:

```matlab
% New usage (recommended)
config = load_config_yaml('config/hinec_default.yml');
runTractography('data.mat', config);
```

#### For Script Users

Update shell scripts to use config files:

**Before**:
```bash
./run_hinec.sh data/subject processed.mat
```

**After** (optional, but recommended):
```bash
./run_hinec.sh data/subject processed.mat config/my_experiment.yml
```

### Future Enhancements

Potential improvements:

1. JSON support (in addition to YAML)
2. GUI config editor
3. Config validation CLI tool
4. Parameter optimization suggestions
5. Automatic config generation from logs

### Technical Details

#### YAML Parser

Simple MATLAB-compatible YAML parser that handles:

- Key-value pairs
- Nested sections
- Comments
- Booleans (`true`/`false`)
- Numbers (integers and floats)
- Strings (quoted and unquoted)

**Limitations**:

- No support for arrays/lists
- No support for complex YAML features (anchors, aliases, etc.)
- For HINEC config needs only

For complex YAML needs, consider visiting: [yamlmatlab](https://github.com/ewiger/yamlmatlab)

> **Superseded.** Parsing now lives in `src/nim_utils/nim_yaml_parse.m`. It
> supports inline lists (`roi: [41, 42]`), which ROI seeding and the track
> filters need, and it is indentation-aware: two levels below a section are
> accepted and a third is a parse error rather than a silent misparse. Anchors
> and aliases are still unsupported.

#### Validation Strategy

Two-phase validation:

1. **Syntax validation** - During YAML parsing
2. **Semantic validation** - After loading complete config

**Validation errors halt execution early** - prevents invalid tractography runs.

### Credits

Implementation Date: 2025-01-18

Based on code review findings and parameter management requirements.

### See Also

- [RKF Implementation](RKF.md) - RKF45 technical details
- [Tractography Guide](TRACTOGRAPHY.md) - the FACT tracker
- [YAML Configuration Reference](YAML_CONFIG.md) - the current, generated
  parameter reference
