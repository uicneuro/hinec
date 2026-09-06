# HINEC Configuration Reference

<!-- GENERATED FILE - do not edit by hand.
     Regenerate with:  nim_config_docs('docs/YAML_CONFIG.md')
     Source of truth:  src/nim_utils/nim_config_schema.m -->

Configuration is YAML with **two levels of nesting** below a section
(`section` -> `group` -> `key`). A third level is a parse error, deliberately:
deeper nesting was judged more confusing than helpful.

Every key is optional. Anything you omit takes the default below, so a working
config can be very short.

## Minimal example

```yaml
tractography:
  algorithm: hinec
  seeding:
    roi: [SLF_L, SLF_R]
```

## Overriding from the command line

Every parameter in this reference can be overridden with `--set`, on both
`bin/run_hinec.sh` and `bin/run_tractography.sh`:

```bash
./bin/run_tractography.sh hinec_dti --set tractography.integrator.step=0.05
./bin/run_tractography.sh hinec_dti --set integrator.step=0.05    # section assumed
./bin/run_tractography.sh hinec_dti --set upsample=2              # bare leaf, if unique
./bin/run_hinec.sh data/x x.mat config/y.yml --set preprocessing.run_eddy=false
./bin/run_tractography.sh hinec_dti --set seeding.roi=[41,42]     # lists
```

Overrides are checked against the schema. An unknown or mistyped key is an
error, never a silent no-op. A bare leaf name that is ambiguous (`method`,
`fa_min`) is rejected with the candidate full paths.

## Parameters

The **Applies to** column says which tracking algorithms actually read the
parameter. A key marked `hinec` is ignored by `standard` and `mmf`.

### `tractography`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `algorithm` | string | `hinec` | all | Tracking algorithm. |
| `field` | string | `dti` | hinec, mmf | Direction source: DTI principal eigenvector, or CSD FOD peaks. |
| `act` | logical | `false` | hinec | Anatomically constrained tracking using WM/GM/CSF masks. |
| `diagnostics` | logical | `true` | all | Write per-run diagnostic reports. |

#### `tractography.integrator`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `method` | string | `rk4` | hinec, mmf | Numerical stepping scheme. NOTE: this is a method NAME, not an order claim. rkf45 here is Dormand-Prince: it advances on the 5th-order solution and uses the embedded 4th-order one only to size the next step. |
| `step` | numeric | `0.2` | all | Integration step in voxels. Fixed step, or the initial step for rkf45. |
| `tolerance` | numeric | `0.01` | hinec, mmf | Adaptive error tolerance in voxels. rkf45 ONLY. |
| `step_min` | numeric | `0.01` | hinec, mmf | Minimum adaptive step. rkf45 ONLY. |
| `step_max` | numeric | `1` | hinec, mmf | Maximum adaptive step. rkf45 ONLY. |
| `safety` | numeric | `0.9` | hinec | Adaptive step safety factor. rkf45 ONLY. |
| `adaptive` | logical | `true` | hinec | Adaptive step-size control. rkf45 ONLY. rkf45 with adaptive:false is fixed-step RKF45 (a real third mode, kept for completeness). |

#### `tractography.interpolation`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `method` | string | `trilinear` | hinec, mmf | Spatial interpolation kernel for the direction field. These differ in SMOOTHNESS, which caps how much of an integrator formal order is reachable: a Runge-Kutta method of order p needs a right-hand side with p continuous derivatives. trilinear is C0 (kinked at every voxel face), cubic is C1 (Keys cubic convolution, NOT a spline - its second derivative jumps), spline is C2 (a genuine cubic spline). Measured here: RK4 reaches observed order 2.00 on trilinear, 3.06 on cubic and 4.00 on spline - one order per continuous derivative, on an unchanged tableau. |
| `upsample` | numeric | `1` | hinec, mmf | Spatial sampling factor for the direction field: the field is sampled on a grid of spacing 1/upsample voxels before the interpolants are built. 1 = the acquisition grid. Above 1 refines toward the continuous field the samples imply; BELOW 1 coarsens, discarding spatial information, which is how the space axis of a convergence study is swept. The coordinate frame is unchanged, so positions, step sizes and lengths stay in native voxel units and runs at different factors are directly comparable. Note the u -> infinity limit is the native-resolution interpolant, not ground-truth anatomy. |

#### `tractography.seeding`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `density` | numeric | `8` | all | Seeds per seeded voxel - honoured exactly. Placed on a deterministic sub-voxel lattice (perfect cubes) or a spread subset of one (other counts), so runs are reproducible. |
| `strategy` | string | `uniform` | all | Seed placement. uniform = deterministic lattice (reproducible); random = jittered. |
| `fa_min` | numeric | `0.05` | all | Minimum FA for a voxel to be seeded. Default 0.05 excludes CSF only. |
| `roi` | list | `[]` | all | Restrict seeding to these atlas regions. JHU indices and/or names, freely mixed. Empty = whole brain mask. |
| `roi_dilate` | numeric | `0` | all | Dilate the seed ROI by this many voxels before seeding. |

#### `tractography.termination`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `fa_min` | numeric | `0.1` | all | Stop tracking below this FA. (NOT the legacy fa_threshold - see nim_config_retired.) |
| `angle_max` | numeric | `225` | all | Maximum turn in DEGREES PER VOXEL OF ARC, i.e. a minimum radius of curvature R = 57.3/angle_max voxels. The default 225 is the classic 45 deg/step evaluated at the default step of 0.2. Step-invariant: the budget for one step is angle_max x the NOMINAL step arc (never the realised chord - that would make the budget method-dependent), so refining the step does not loosen the constraint. CEILING: tangents are sign-aligned because v1 is a line field, so a measured turn never exceeds 90 deg - any angle_max above 90/step is INERT, not merely loose, and 225 goes inert for any step >= 0.4. 0 disables the criterion outright, which is what a control run should use. |
| `max_arc` | numeric | `400` | all | Maximum track arc length in voxels. Step-invariant: max_steps is derived as ceil(max_arc/step). |
| `min_arc` | numeric | `15` | all | Discard tracks shorter than this arc length in voxels. |

#### `tractography.filter`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `include_roi` | list | `[]` | all | Keep only tracks touching these regions (waypoint). |
| `exclude_roi` | list | `[]` | all | Discard tracks touching these regions. |
| `mode` | string | `all` | all | Whether a track must touch ALL include regions or ANY of them. |
| `roi_dilate` | numeric | `0` | all | Dilate the include/exclude masks by this many voxels before testing. |
| `endpoints_in` | list | `[]` | all | Two regions; keep a track only if one END lands in the first and the other END in the second, either way round. This is an ENDPOINT test, not a waypoint test - include_roi asks whether a track passes through a region, this asks where it stops. It is half of how the ISMRM 2015 scorer defines a bundle (head + tail). |
| `contained_in` | list | `[]` | all | Keep a track only if EVERY point lies inside these regions. The other half of the ISMRM bundle definition (all_mask): a streamline that wanders outside the corridor is not that bundle, however it ends. |

#### `tractography.output`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `arc_step` | numeric | `0` | all | Resample saved streamlines to this arc-length spacing in voxels. 0 = store every integration step. Decouples file size from step size; integration accuracy is unaffected. |

#### `tractography.csd`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `lmax` | numeric | `6` | hinec, mmf | Spherical harmonic order for CSD. |
| `max_peaks` | numeric | `3` | hinec, mmf | Maximum FOD peaks retained per voxel. |
| `peak_thresh` | numeric | `0.5` | hinec, mmf | Relative amplitude threshold for accepting an FOD peak. |
| `peak_min_sep` | numeric | `45` | hinec, mmf | Minimum angular separation between FOD peaks (degrees). |
| `n_iter` | numeric | `50` | hinec, mmf | CSD deconvolution iterations. |

#### `tractography.mmf`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `anchor` | numeric | `0` | mmf | Re-anchor strength of e1 toward the field tangent. 0 = pure connection-form evolution. |
| `frame_sel_power` | numeric | `16` | mmf | Directional selectivity used when building the moving frame. |

### `preprocessing`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `run_denoising` | logical | `true` | - | Run denoising before tensor fitting. |
| `denoise_method` | string | `dwidenoise` | - | Denoising backend. Falls back to gaussian when MRtrix3 is absent. |
| `run_motion_correction` | logical | `true` | - | Run motion correction (FSL). |
| `run_eddy` | logical | `true` | - | Run eddy-current correction (FSL). |
| `improve_mask` | logical | `true` | - | Refine the brain mask using FA. |
| `atlas_type` | string | `jhu` | - | Atlas used for parcellation. |
| `bundle_roi_dir` | string | `` | - | Optional directory of bundle masks to use as the parcellation INSTEAD of the atlas, e.g. data/ismrm2015/scoring_data_Renauld2023/ROI. Must contain all_masks/ (containment corridors, which become the labels) and may contain endpoints/ and any_masks/ (gates, addressable by name but not labels). The atlas parcellation is preserved as parcellation_mask_<atlas_type> rather than discarded. Empty = atlas only. Set this when ROIs should be addressed by the names a scorer uses: an atlas label and a bundle of the same name are not the same region. |
| `t1_available` | logical | `false` | - | A T1 volume is present alongside the DWI. |
| `use_t1_registration` | logical | `false` | - | Register T1 to DWI space. |
| `register_to_mni` | logical | `false` | - | Register to MNI space. |

## Migrating from the old flat config

Old configs still load; each superseded key produces a deprecation warning
naming its replacement.

| Old key | Replacement |
|---|---|
| `act_enabled` | `tractography.act` |
| `integrator` | `tractography.integrator.method` |
| `step_size` | `tractography.integrator.step` |
| `rkf_tolerance` | `tractography.integrator.tolerance` |
| `rkf_tol` | `tractography.integrator.tolerance` |
| `step_min` | `tractography.integrator.step_min` |
| `step_max` | `tractography.integrator.step_max` |
| `rkf_safety` | `tractography.integrator.safety` |
| `adaptive_step` | `tractography.integrator.adaptive` |
| `interp_method` | `tractography.interpolation.method` |
| `seed_density` | `tractography.seeding.density` |
| `seed_strategy` | `tractography.seeding.strategy` |
| `seed_fa_threshold` | `tractography.seeding.fa_min` |
| `termination_fa` | `tractography.termination.fa_min` |
| `min_length` | `tractography.termination.min_arc` |
| `csd_lmax` | `tractography.csd.lmax` |
| `csd_max_peaks` | `tractography.csd.max_peaks` |
| `csd_peak_thresh` | `tractography.csd.peak_thresh` |
| `csd_peak_min_sep` | `tractography.csd.peak_min_sep` |
| `mmf_anchor` | `tractography.mmf.anchor` |
| `frame_sel_power` | `tractography.mmf.frame_sel_power` |
| `enable_diagnostics` | `tractography.diagnostics` |
| `integration_order: 1\|2\|4\|5` | `integrator.method: euler\|rk2\|rk4\|rkf45` |
| `max_steps` | `termination.max_arc` (converted as `max_steps x step`) |

`integration_order` was a **method selector wearing a numeric-order name**:
the value `5` selected RKF45. That value was not wrong numerically - the
implementation uses Dormand-Prince coefficients and advances on the 5th-order
solution, using the embedded 4th-order one only for error control - but a
method selector should not be spelled as a number. `integrator.method` names
the method instead.

`max_steps` counted integration steps, so halving the step size silently
halved how far a track could travel. `termination.max_arc` is an arc length
in voxels and is step-invariant; `max_steps` is now derived as
`ceil(max_arc / step)`.

## Retired keys

These were accepted by the old loader but read by no tracker. They are
ignored with a warning rather than silently accepted.

| Key | Why it is gone |
|---|---|
| `gate_power` | Never read by any tracker. Defaulted and validated by the old loader but had no effect. |
| `crossing_cp` | Never read by any tracker. |
| `curv_beta` | Never read by any tracker. |
| `crossing_detect` | Never read by any tracker. |
| `swing_ratio_max` | Never read by any tracker. |
| `transport_gate` | Never read by any tracker. |
| `transport_strength` | Never read by any tracker. |
| `bishop_eps` | Referenced only in a comment; no tracker reads it as an option. |
| `fa_threshold` | Functionally dead: printed by hinec/standard, and only a seed-mask fallback in mmf/highorder that never fires because runTractography always supplies a seed mask. Use termination.fa_min to stop tracking, or seeding.fa_min to restrict seeding. NOTE: in 5 of the 12 original configs fa_threshold and termination_fa held DIFFERENT values, so they were never interchangeable. |
| `order` | Legacy backward-compatibility key, read by no tracker as an option. |
| `sel_power` | Removed. HINEC is now pure interpolation + integration. It biased interpolation toward the incoming direction (weight = alignment^sel_power), which has no justification for DTI (one eigenvector per voxel, nothing to disambiguate) and made the ODE direction-dependent. For CSD, nearest-peak selection is retained because it is structural, but the alignment exponent is gone. |
