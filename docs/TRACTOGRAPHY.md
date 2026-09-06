# Standard Tractography: the FACT Tracker

`nim_tractography_standard` implements FACT (Fiber Assignment by Continuous
Tracking, Mori et al. 1999): deterministic streamline tracking on the discrete
voxel tensor field, with **no interpolation**. Direction is read from the
principal eigenvector of the voxel the streamline is currently in, and the
streamline advances by intersecting the voxel boundary rather than by taking a
fixed step. It is the baseline tracker in HINEC — the reference against which
the interpolated trackers (`nim_tractography_hinec`,
`nim_tractography_mmf_connframe`) are compared.

Selected with `algorithm: standard` in the config. FACT is DTI-only by
construction: `field: csd` has no meaning here, because there is no
interpolation stage in which to resolve multiple peaks.

| | |
|---|---|
| Source | `src/nim_tractography/nim_tractography_standard.m` |
| Driver | `runTractography.m` (seeding, masks, filtering, saving) |
| Config | `config/standard_dti.yml` |
| Interpolated alternatives | [High-Order Methods](High_Order.md), [MMF Connection-Form](MMF_TRACTOGRAPHY.md) |

## Algorithm

FACT does not step by a fixed increment. From the current position it finds
where the ray along the current direction leaves the current voxel, jumps to
that exit point, and reads a new direction from the voxel it has entered. Step
length therefore varies: it depends on where in the voxel the streamline
entered and how obliquely it is travelling.

```
place seed anywhere inside a seeded voxel
direction <- principal eigenvector of the seed voxel, signed by +1 / -1
loop:
    stop if FA at the current voxel < termination.fa_min
    stop if the turn from the previous segment exceeds the angle budget
    stop if the crossing count reaches max_steps
    exit_point <- intersection of the ray with the current voxel boundary
    stop if there is no forward intersection, or the exit point leaves the volume
    stop if the exit voxel is outside the propagation mask
    append exit_point; move there
    direction <- principal eigenvector of the new voxel, sign-aligned to the old
```

Two passes are run from every seed, one with the initial direction negated. The
backward pass is reversed and concatenated with the forward pass, so **one seed
yields at most one track**, spanning the seed in both directions.

Implementation notes:

- The exit point is nudged `1e-4` along the direction before being stored.
  Without it, `round()` maps a boundary coordinate of the form `n.5` back to
  the voxel just left (MATLAB rounds half away from zero) and the streamline
  stalls.
- Voxel membership is `round(position)` throughout, so a voxel is the unit cube
  centred on its integer index.
- Direction is sign-aligned to the previous direction at every crossing
  (`dot < 0` flips it), because `v1` is a line field defined only up to sign.

## Inputs

Required fields on the `nim` struct:

| Field | Shape | Use |
|---|---|---|
| `nim.evec` | `[x y z 3 3]` | principal eigenvector `evec(:,:,:,:,1)` is the direction field |
| `nim.FA` | `[x y z]` | seeding, propagation and termination |
| `nim.eval` | `[x y z 3]` | optional; used only for a sanity warning at the centre voxel |
| `nim.mask` | `[x y z]` | optional; intersected into the propagation mask |

The eigenvector array is validated to be `3 x 3` in its last two dimensions,
and the eigenvalue ordering is checked at the centre voxel: if `eval(1)` is not
the largest, a warning is issued that the principal eigenvector may not
correspond to the largest eigenvalue.

## Seeding

**Seeding is decided in `runTractography.m`, not in the tracker.** The tracker
errors out if `options.seed_mask` is empty — it executes a seeding decision, it
does not make one. `runTractography` picks a mask by priority:

1. `seeding.roi` — explicit atlas regions, optionally dilated, then intersected
   with the brain mask and `FA > seeding.fa_min`.
2. `nim.mask` refined by `FA > seeding.fa_min`.
3. `nim.parcellation_mask` dilated by 3 voxels, then `FA > 0.05`.
4. `FA > 0.10` as a last resort, with a warning that low-anisotropy structures
   (fornix, cingulum) will be missed.

Within each seeded voxel, `seeding.density` seeds are placed by
`nim_seed_offsets`:

- `strategy: uniform` (default) — a deterministic sub-voxel lattice. The count
  is honoured **exactly**: a perfect cube gives the full `n x n x n` lattice,
  and any other count gives a deterministic farthest-point subset of the next
  larger lattice, so seeds stay spread rather than clustering in one corner. No
  RNG is involved, so runs are reproducible and streamlines can be matched
  one-to-one across a parameter sweep.
- `strategy: random` — uniform jitter within the voxel. Not reproducible.

See [Seeding Strategy](SEEDING_STRATEGY.md) for the full account.

## The propagation mask

Where a track may **go** is not where it may **start**. The propagation domain
is derived from the data:

```matlab
prop_mask = nim.FA > options.termination_fa;
if isfield(nim, 'mask') && ~isempty(nim.mask)
    prop_mask = prop_mask & (nim.mask > 0.5);
end
nim.dilated_brain_mask = imdilate(prop_mask, ones(3,3,3));
```

This was previously a dilation of the *seed* mask. With whole-brain seeding the
two coincide and nothing looks wrong; with ROI seeding the propagation domain
collapsed to a one-voxel skin around the ROI and streamlines died almost
immediately. The same bug had been fixed in `nim_tractography_hinec` earlier and
the fix was not carried across. A caller may still override the domain by
passing `options.propagation_mask` directly.

## Termination

| Criterion | Config key | Reason string |
|---|---|---|
| FA at the seed voxel below floor | `termination.fa_min` | `initial_fa` |
| FA at the current voxel below floor | `termination.fa_min` | `fa` |
| Turn exceeds the angle budget | `termination.angle_max` | `angle` |
| Crossing count reached | derived from `termination.max_arc` | `max_steps` |
| No forward boundary intersection | — | `no_exit` |
| Exit voxel outside the propagation mask | — | `brain_mask` |
| No valid eigenvector in the new voxel | — | `no_direction` |

Counts per reason are printed at the end of the run when
`tractography.diagnostics` is true.

After tracking, a track is kept only if its arc length — the sum of segment
lengths of the combined forward+backward polyline, in **voxels** — is at least
`termination.min_arc`. This filter was previously defaulted at the top of the
file and then never read; the tracker carried a comment stating it applied "no
filters at all", which was literally true. A schema-validated key that every
config sets was honoured by exactly one of the three trackers.

### The angle budget

`termination.angle_max` is in **degrees of turning per voxel of arc**, not per
step. All three trackers convert it through `nim_utils/nim_angle_limit.m`. It is
a minimum radius of curvature in disguise: `R = 57.3 / angle_max` voxels.

**FACT is budgeted per voxel transition.** Direction is piecewise constant
within a voxel and changes only at boundaries, so one transition is charged as
one voxel of arc — `nim_angle_limit(angle_thresh, 1)`. This is deliberately
*not* the length of the segment just travelled: a streamline that clips a voxel
corner covers a tiny arc, but the turn waiting for it at the boundary is a
property of the data, not of how much of the voxel it happened to cross.
Budgeting by the clipped length would terminate a streamline for the geometry
of the intersection.

**The 90 degree ceiling.** `v1` is a line field, defined only up to sign, so
every tangent is sign-aligned to its predecessor before the turn between them is
measured. That confines the measured turn to `[0, 90]` degrees. A budget above
90 degrees can never be exceeded, so the criterion there is not loose — it is
**absent**. `nim_config_to_options` warns when a config lands in that regime.

`angle_max: 0` disables the criterion outright. Use that for a control run
rather than setting a very large value: a large value is a silently inert
criterion, which looks exactly like a satisfied one.

For the interpolated trackers the budget for one step is `angle_max x step`, so
refining the step tightens the budget in proportion and the same physical
curvature bound is enforced at every step size. The budget is always taken from
the *nominal* step, never the realised chord — chord-based budgeting is
method-dependent, since Euler and RK2 advance by `h x (unit vector)` and have a
chord of exactly `h`, while RK4 averages four stage vectors and its chord is
shorter. `config/standard_dti.yml` sets `angle_max: 45`, the classic FACT
turning angle, evaluated per voxel transition.

## Configuration

Only the keys FACT actually reads are listed. The full surface, with defaults,
is in the generated [YAML Configuration Reference](YAML_CONFIG.md); the mapping
from canonical config path to the flat option names the tracker reads lives in
`src/nim_utils/nim_config_to_options.m`.

| Config path | Tracker option | Effect on FACT |
|---|---|---|
| `tractography.algorithm: standard` | — | selects this tracker |
| `tractography.seeding.density` | `seed_density` | seeds per seeded voxel |
| `tractography.seeding.strategy` | `seed_strategy` | `uniform` lattice or `random` jitter |
| `tractography.seeding.fa_min` | `seed_fa_threshold` | seed mask FA floor (applied in `runTractography`) |
| `tractography.seeding.roi` / `roi_dilate` | `seed_roi` / `seed_roi_dilate` | restrict seeding to atlas regions |
| `tractography.termination.fa_min` | `termination_fa` | termination and propagation-mask FA floor |
| `tractography.termination.angle_max` | `angle_thresh` | turn budget per voxel transition |
| `tractography.termination.max_arc` | `max_steps` (derived) | see the caveat below |
| `tractography.termination.min_arc` | `min_length` | discard shorter tracks, in voxels |
| `tractography.filter.*` | — | applied by `runTractography` after tracking |
| `tractography.output.arc_step` | `output_arc_step` | resample saved polylines |
| `tractography.diagnostics` | `enable_diagnostics` | timing and failure reports |

Keys FACT ignores by design: `integrator.method` and the RKF45 group (there is
no numerical integrator here — steps are exact boundary intersections),
`interpolation.*` (no interpolation), `csd.*` and `field`, `mmf.*`, and
`tractography.act` (ACT is read only by `nim_tractography_hinec`; the schema
marks it `hinec`).

**`integrator.step` caveat.** FACT's geometry does not use a step size, but the
step is not inert: `nim_config_to_options` derives `max_steps =
ceil(max_arc / step)` for every tracker, and FACT compares that number against
its **crossing count**. So under FACT `max_arc` and `step` together set a cap on
the number of voxel boundaries crossed, not on arc length in voxels; the two
coincide only if each crossing happens to advance exactly `step` of arc. With
`config/standard_dti.yml` (`max_arc: 2500`, `step: 0.5`) the cap is 5000
crossings.

## Output

`nim_tractography_standard` returns `tracks` and writes nothing itself. Saving,
ROI filtering and decimation are `runTractography`'s job, in that order:

1. `nim_filter_tracks_roi` applies `filter.include_roi`, `filter.exclude_roi`,
   `filter.endpoints_in` (one end in each of two regions) and
   `filter.contained_in` (every point inside a corridor). Filtering runs on the
   full-resolution polyline, because it tests which voxels a track visits.
2. `nim_resample_track_arc` decimates to `output.arc_step` voxels of spacing, if
   set. This decouples file size from step size and does not affect the tracking
   itself.
3. The result is saved to `<run_dir>/tractography/tracks_standard_<timestamp>.mat`
   (or `tractography_results/` without a run directory), containing `tracks`,
   `options`, `elapsed_time`, `algorithm` and `track_meta`.

### Track data structure

```matlab
tracks = cell(N, 1);
tracks{i} = [x1, y1, z1;    % complete trajectory, not just endpoints
             x2, y2, z2;
             ...
             xN, yN, zN];   % voxel-space coordinates
```

Every track is the full polyline; row counts vary per track. Coordinates stay in
voxel space throughout — conversion to world/RAS coordinates happens only at
export, in `scripts/hinec_to_trk.py`, using the DWI affine.

## Diagnostics

With `tractography.diagnostics: true` the tracker prints, after the run:

- a timing report — total, eigenvector verification, seed generation, tracking,
  boundary-check time, total voxel steps, average voxels per track, voxels per
  second;
- a failure analysis — seeds with no initial direction, successful tracks,
  failures during generation;
- the termination-reason breakdown tabulated above, sorted by frequency;
- a warning when the success rate falls below 10%, which usually means the FA
  floor or the seed mask is wrong.

`runTractography` additionally writes `track_statistics.txt` into the run
directory's diagnostics folder.

## Limitations

- **Single tensor model.** One principal eigenvector per voxel cannot represent
  crossing fibres, and streamlines terminate or take the dominant branch at a
  crossing. Use `algorithm: hinec` with `field: csd` for multi-peak direction
  resolution.
- **No interpolation.** Direction is piecewise constant, so tracks are visibly
  faceted and turn only at voxel boundaries. This is FACT's definition, not a
  gap; it is also why FACT has no meaningful order of accuracy to verify — see
  [Solution Verification](CONVERGENCE.md), which covers the interpolated
  trackers.
- **Deterministic only.** There is no probabilistic variant; low-probability
  pathways are not sampled.
- **Voxel-space output.** Anything requiring physical units (millimetre lengths,
  world-space overlays) must apply the affine downstream.

## See also

- [Methods Overview](TRACTOGRAPHY_METHODS.md) — how the three trackers relate
- [High-Order Methods](High_Order.md) — the interpolated tracker
- [MMF Connection-Form](MMF_TRACTOGRAPHY.md) — moving-frame tracking
- [Seeding Strategy](SEEDING_STRATEGY.md)
- [YAML Configuration Reference](YAML_CONFIG.md) and its
  [changelog](CHANGELOG_YAML.md)
- [Visualization Guide](VISUALIZATION_GUIDE.md) — `visualizeTractography`,
  `visualizeTractographySlices` and the slice-cache viewer

## Reference

Mori, S., Crain, B. J., Chacko, V. P., & van Zijl, P. C. (1999).
Three-dimensional tracking of axonal projections in the brain by magnetic
resonance imaging. *Annals of Neurology*, 45(2), 265-269.
