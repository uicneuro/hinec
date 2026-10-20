# Writing your own tracker

This page is the contract between `runTractography.m` and a tracker: exactly
what a tracker is handed, exactly what it must return, and how to plug a new
one in. Everything here is also **logged by every run** — the block headed
`=== Tracker input (step 4 hand-off) ===` in the run log, saved as
`<run_dir>/tractography/tracker_input.txt` — so if this page and a log ever
disagree, the log is right and this page needs fixing.

The shortest route: copy
[`src/nim_tractography/nim_tractography_template.m`](../src/nim_tractography/nim_tractography_template.m)
(≈100 lines, runnable as `algorithm: template`), replace its two marked
sections, and register the name. The rest of this page is what that file
assumes.

## Where a tracker sits in the pipeline

`runTractography` is a fixed sequence; a tracker is step 4 and nothing else.

| step | function | what it does | tracker sees the result as |
|---|---|---|---|
| 1 | `nim_load_nim` | load the processed dataset from disk | `nim` (dataset fields) |
| 2 | `nim_field` | build the direction field the config asks for (`field: csd` → FOD peaks; `dti` needs nothing) | `nim.peaks…` |
| 3 | *(in `runTractography`)* | decide **where** seeds may go: brain mask ∩ FA ≥ `seeding.fa_min`, or the ROI in `seeding.roi` | `options.seed_mask` |
| 3b | *(in `runTractography`)* | copy the WM/GM/CSF masks over when `act: true` | `options.wm_mask`, `gm_mask`, `csf_mask` (or `[]`) |
| **4** | **`nim_tractography_<algorithm>`** | **integrate streamlines** | — |
| 5 | `nim_filter_tracks_roi`, `nim_resample_track_arc` | ROI **selection** (keep/drop flags only — nothing is removed or cut), then decimation to `output.arc_step` | — |
| 6 | `save` | `tracks`, `options`, `elapsed_time`, `algorithm`, `track_meta` → `tractography/tracks_<algorithm>_<timestamp>.mat` — always the **full** tractogram; when a `filter.*` key is set, the selection goes beside it as `tractography/roi_selection.mat` (`keep`, the kept `tracks`, subset `track_meta`, per-criterion drop counts in `stats`) | — |

Two consequences worth stating:

* A tracker **generates its own seed points** from `options.seed_mask` (a voxel
  mask) and `options.seed_density` (seeds per voxel). It does not decide where
  seeding is allowed — that is step 3's job, so that every tracker seeds the
  same voxels for the same config. Use `nim_seed_offsets(density)` for the
  sub-voxel positions; it is deterministic, which is what makes per-streamline
  comparison across runs possible.
* A tracker does **not** filter by region or resample its output. It returns
  every integration point of every track that passes the length rule; steps 5
  and 6 do the rest.

## The call

```matlab
[tracks, meta] = nim_tractography_<algorithm>(nim, options)
```

**Coordinates.** Everything is voxel space, 1-based, continuous: voxel
`(i,j,k)` of the arrays is centred at position `[i j k]`, its faces at ±0.5. A
step of `h` moves `h` voxels. Angles are in degrees. Conversion to scanner
millimetres happens once, at scoring time (`scripts/hinec_to_trk.py`, using
`nim.hdr`), never inside a tracker.

## Input 1 — `nim`

`nim` is a struct. The fields below are grouped by **who wrote them**, because
that determines whether a field is always there.

### Dataset fields — written by `main.m`, on disk, config-independent

Always present for any processed dataset.

| field | layout | meaning |
|---|---|---|
| `hdr` | struct | NIfTI header (`niftiinfo`): `ImageSize`, `PixelDimensions`, `Transform` (voxel → world affine) |
| `xdim`, `ydim`, `zdim` | scalars | volume size `X`, `Y`, `Z` |
| `img` | `[X Y Z V]` single | raw DWI signal, `V` volumes |
| `bval` | `[V 1]` | b-value per volume |
| `bvec` | `[V 3]` | unit gradient direction per volume |
| `img_b0`, `img_bi` | `[X Y Z]`, `[X Y Z V-nb0]` | mean b0 image; diffusion-weighted volumes only |
| `mask` | `[X Y Z]` | brain mask, 1 inside. Tracking stops on leaving it |
| `DT` | `[X Y Z 6]` | diffusion tensor as `[Dxx Dyy Dzz Dxy Dyz Dxz]`; `nim_reshape_d` rebuilds the 3×3 |
| `evec` | `[X Y Z 3 3]` | eigenvectors: `evec(x,y,z,:,k)` is the k-th, **k = 1 is the principal**. Sign is arbitrary voxel to voxel — see below |
| `eval` | `[X Y Z 3]` | eigenvalues, descending. `[0 0 0]` where the fit was skipped (mask edge) |
| `FA` | `[X Y Z]` | fractional anisotropy; the termination criterion |
| `parcellation_mask` | `[X Y Z]` uint16 | region label per voxel, 0 = none |
| `atlas_labels` | struct | label index → name |
| `roi_masks` | `containers.Map` | bundle name → `[X Y Z]` logical (ISMRM scoring ROIs); read by `nim_roi_mask` |
| `wm_mask`, `gm_mask`, `csf_mask` | `[X Y Z]` | tissue masks for ACT (present only if `main.m` ran tissue segmentation) |

Plus bookkeeping you will not need: `size3`, `size_b0`, `size_bi`,
`thrsh_b0`, `parcellation_mask_file`, `*_mask_file`, `parcellation_mask_jhu`,
`atlas_labels_jhu`, `atlas_type`, `roi_source`.

**The sign of `evec`.** The principal eigenvector is a *line* field: the
eigensolver returns `v` or `−v` at random from voxel to voxel. Interpolating
the components directly averages a vector with its negative and collapses
toward zero. Either pick the sign per voxel to agree with your incoming
direction (what FACT does — `nim_tractography_standard`), or interpolate the
dyadic `v·vᵀ` and take its principal axis afterwards (what `hinec` and the
template do — `nim_principal_dir` gives you the closed form).

### Direction-field fields — built per run by `nim_field` (step 2)

Present only when the config asks for that field. Native fitted fields may be
cached beside the dataset because they cost minutes; they are not saved inside
the nim. For `field: csd`, `tractography.csd.peaks_file` can supply a peak MAT
on the same voxel grid instead of fitting native CSD.

| `field:` | field | layout | meaning |
|---|---|---|---|
| `csd` | `peaks` | `[X Y Z P 3]` | FOD peak directions, unit, largest first. `P = csd.max_peaks` |
| `csd` | `npeaks` | `[X Y Z]` | how many of the `P` slots are valid at each voxel |
| `csd` | `peak_w` | `[X Y Z P]` | FOD amplitude of each peak |
| `csd` | `fod_sh` | `[X Y Z C]` | FOD spherical-harmonic coefficients (optional) |

A tracker that needs one of these should `error` if it is absent, with the
config key that provides it, as `nim_tractography_hinec` does:
`HINEC field=csd needs nim.peaks/npeaks (run nim_csd).`

### Per-run derived fields

If your algorithm needs a derived field of its own — one that depends on the
tractography config — build it in a function called from `runTractography`
under your algorithm's name, after step 2 and before seeding; do not compute it
inside `main.m` (the dataset must stay config-independent) and do not save it
into the nim.

## Input 2 — `options`

`options` is a flat struct. Every field comes from one YAML key via
`nim_config_to_options` (the schema in `nim_config_schema.m` is the single
source of truth for defaults, ranges and help text), except the four that
`runTractography` adds in steps 3 and 3b. The run log prints each option
**with its value, its config key and its meaning**; the table here gives the
same for the fields a tracker is expected to honour.

### Read these

| option | from | meaning |
|---|---|---|
| `seed_mask` | step 3 | `[X Y Z]` logical. **The** seed input: place `seed_density` seeds in each true voxel |
| `seed_density` | `seeding.density` | seeds per seeded voxel — honour it exactly |
| `seed_strategy` | `seeding.strategy` | `uniform` (deterministic lattice, `nim_seed_offsets`) or `random` |
| `step_size` | `integrator.step` | `h` in voxels (initial step for adaptive schemes) |
| `termination_fa` | `termination.fa_min` | stop when the interpolated FA drops below this |
| `angle_thresh` | `termination.angle_max` | maximum turn in **degrees per voxel of arc**, i.e. minimum radius of curvature `57.3 / angle_max` voxels. For a fixed step, the test is `turn_this_step / h > angle_thresh` — so the limit means the same thing at every `h` |
| `max_arc` | `termination.max_arc` | stop a half-track after this arc length in voxels (`max_steps = ceil(max_arc/step)` is the same limit as a count) |
| `min_length` | `termination.min_arc` | discard a finished track whose chord length `sum(‖diff(track)‖)` is below this |
| `wm_mask`, `gm_mask`, `csf_mask` | step 3b | `[X Y Z]` when ACT is on, `[]` when off. Trackers decide whether ACT is active from these, not from `act_enabled` |
| `field` | `tractography.field` | which direction source to read: `dti` → `nim.evec`, `csd` → `nim.peaks` |

### Read these if the scheme applies

| option | from | meaning |
|---|---|---|
| `integrator` / `integration_order` | `integrator.method` | scheme name (`euler`, `rk2`, `rk4`, `rkf45`) and the same choice as a number (1, 2, 4, 5 — a selector, not an order claim) |
| `adaptive_step`, `rkf_tolerance` (= `rkf_tol`), `step_min`, `step_max`, `rkf_safety` | `integrator.*` | adaptive-step control, `rkf45` only |
| `interp_method`, `upsample` | `interpolation.*` | kernel (`trilinear` C⁰, `cubic` C¹, `spline` C²) and pre-sampling factor for the direction field |
| `trace`, `trace_max` | `debug.*` | if you support tracing, record a per-step trace for `trace_max` evenly spaced seeds into `meta.trace` (see `nim_tractography_hinec` for the record layout and `docs/CONVERGENCE.md` for what it is used for) |
| `mmf_mode` | `mmf.mode` | how the carried frame relates to the data frame (`transport`, `reanchor`, `adaptive`); `mmf` only |
| `mmf_reanchor_tol` | `mmf.reanchor_tol` | drift (radians) between the carried and the data's `e1` at which the adaptive mode re-anchors; `mmf` only |
| `mmf_gradient_delta` | `mmf.gradient_delta` | central-difference step (voxels) for the derivative of the interpolated principal axis that forms the connection; `mmf` only |

### Present but not yours

These are in the struct because one struct is handed to everything, but they
are consumed elsewhere: `csd_*`, including `csd_peaks_file` (step 2), `include_roi`, `exclude_roi`,
`roi_filter_*`, `endpoints_in`, `contained_in`, `any_in`, `length*` (step 5),
`output_arc_step` (step 6), `seed_roi`, `seed_roi_dilate`, `seed_roi_info`,
`seed_fa_threshold` (already applied to `seed_mask`), `algorithm`,
`enable_diagnostics`. Ignore them.

### Adding a knob

A new algorithm usually wants one or two parameters of its own. Add them to
`nim_config_schema.m` under a group named after the algorithm (as the
`stitching.*` and `mmf.*` keys are), list your algorithm in the entry's `algos`, and map them in
`nim_config_to_options`. Then they are validated, defaulted, documented in
`docs/YAML_CONFIG.md` (generated — run `nim_config_docs`) and overridable from
the shell with `run_tractography.sh --set <group>.<key>=…`. Do not read
`config` directly inside a tracker and do not invent defaults inside the
tracker — the schema is where defaults live.

## Output — `tracks` and `meta`

```
tracks : cell {T x 1}
         tracks{t} is N_t x 3 double, N_t >= 2, voxel coordinates, finite,
         one row per integration point, ordered
             [ backward half, reversed ; seed ; forward half ]
         so the polyline reads end to end and the seed is on it exactly once.
         Include ONLY tracks whose chord length >= options.min_length.

meta   : struct, may be struct() with no fields. Recognised fields:
         seed_index   1 x T   index into your seed list, one per KEPT track
         seed_points  T x 3   the seed position of each kept track
         n_seeds      scalar  how many seeds you attempted
         trace        optional per-step record (debug.trace)
```

`runTractography` runs `nim_check_tracker_output` on the return values before
anything downstream touches them; a wrong shape fails there with the rule it
broke (`tracker:outputContract`). Note that `seed_index` and `seed_points` are
indexed by **track**, not by seed — `roi_selection.mat` subsets them
alongside its `tracks`, which only works if they are already aligned.

`meta` is saved verbatim as `track_meta` in the output `.mat`, so anything
else you put in it (termination-reason counts, timings) is preserved.

## Registering the algorithm

Three edits, all one line:

1. `src/nim_utils/nim_config_schema.m` — add the name to the
   `tractography.algorithm` enum.
2. `runTractography.m`, step 4 — add an `elseif strcmpi(algorithm, '<name>')`
   branch that calls your function and sets `output_filename`.
3. A config: copy `config/hinec_dti.yml` to `config/<name>_dti.yml` and set
   `algorithm: <name>` (naming rule in `config/README.md`).

Then:

```bash
./bin/run_tractography.sh <name>_dti --score
```

runs it on the preprocessed ISMRM 2015 dataset and scores it. The run dir
contains `tractography/tracker_input.txt` (what your function was given),
`tractography/tracks_<name>_*.mat` (what it returned, after steps 5–6) and
`scoring/renauld2023/results.json` (`mean_f1`, `mean_OL`, `mean_OR_gt`, `VB`).
For a fast first check, seed one bundle:

```bash
./bin/run_tractography.sh <name>_dti --set seeding.roi='[Cingulum_right]'
```

## The template, annotated

`nim_tractography_template.m` is the minimal honest tracker: Euler steps along
the trilinearly interpolated principal-eigenvector dyadic, with FA / mask /
edge / turn / arc termination and the length rule. Its structure is the one
every shipped tracker follows:

```
seeds        <- options.seed_mask, options.seed_density, nim_seed_offsets
direction    <- a function  v = direction_at(pos, previous_direction)   [replace]
step         <- a function  pos' = step(pos, v, h)                       [replace]
per seed     <- forward half, backward half (start with -v), join through the seed
length rule  <- keep if chord >= options.min_length
meta         <- seed_index, seed_points, n_seeds
```

`tests/unit/TestTrackerInterface.m` runs it on a synthetic rotation field and
checks the contract; running that test class after copying the template is the
quickest way to know a new tracker is at least well-formed.

## Checklist

- [ ] reads seeds from `options.seed_mask`, honours `seed_density` exactly
- [ ] handles the eigenvector sign (per-voxel flip or dyadic interpolation)
- [ ] stops on FA, mask, volume edge, turn per voxel of arc, `max_arc`
- [ ] applies `min_length` to the chord of the joined track
- [ ] returns `tracks` as `[bwd reversed; seed; fwd]`, N×3, no NaN
- [ ] returns `meta.seed_index` / `seed_points` indexed by kept track
- [ ] `error`s with the config key when a field it needs is missing
- [ ] registered in the schema enum and the step-4 dispatch; has a config
- [ ] `runtests('tests/unit/TestTrackerInterface.m')` passes

### Short tractlet stitching

`algorithm: stitching` dispatches to `nim_tractography_stitching(nim, options)`.
`options.stitching` contains the schema-defined fragment and graph settings;
the tracker runs on `field: dti` or `field: csd`.
See [STITCHING.md](STITCHING.md) for supported kernels/integration and units.

For assembled chains, `meta.seed_index` and `meta.seed_points` identify the first
source fragment. `meta.fragment_seed_indices` preserves all contributing seeds,
and `meta.graph.fragment_ids` preserves the assembly provenance. With tracing,
`meta.fragments` stores the original fragments; graph endpoint IDs `2*i-1` and
`2*i` refer to the first and last point of fragment `i`. These chains remain
ordinary finite N-by-3 voxel-coordinate polylines at the output boundary.

### Moving frames tractography

`algorithm: mmf` dispatches to `nim_tractography_mmf(nim, options)`, which
advances a carried orthonormal frame with the connection of the data's frame
field through the exponential step of `advanceFrame.m`. It reads
`options.mmf_mode`, `options.mmf_reanchor_tol` and `options.mmf_gradient_delta`
(the `mmf.*` group), runs on `field: dti` only (`csd` is rejected), maps
`integrator` `euler` / `rk2` / `rk4` to Lie group order 1 / 2 / 4 (`rkf45` is
rejected), and ignores the ACT masks. The angle budget goes through
`nim_angle_limit` like everywhere else and bounds the rotation of a step and
the turn of the heading as the data see it, in every mode. Besides the contract fields, `meta` carries
`termination` (reason / count over the half-tracks tracked), `mmf` (the scheme
settings, the fraction of re-anchored steps, `seeds_alive` and the
`seeds_skipped` reasons) and, with `debug.trace`, a per-step `trace` whose
records have the hinec layout above (`pos`, `dir`, `fa`, `h`, `chord`, `turn`,
`flip`, `n_steps`, `termination`; `flip` is always false) plus `drift`,
`reanchored`, `kappa` and `s`, so `nim_track_diagnose` reads them. See
[MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md).
