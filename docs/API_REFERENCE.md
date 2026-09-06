# Function Reference

Signatures, arguments, return values and defaults for the public functions in
HINEC, organized by module. Parameter *semantics* for the YAML configuration
surface live in [YAML_CONFIG.md](YAML_CONFIG.md), which is generated from
`src/nim_utils/nim_config_schema.m`; this page covers the MATLAB interfaces.

Two conventions run through the codebase and are worth stating once:

- **Voxel coordinates.** Streamlines, seeds and masks are all in 1-based voxel
  index space of the DWI grid. Only `nim_nifti_affine` and `nim_read_trk` deal in
  world (RAS) millimetres.
- **Options structs, not name-value pairs.** Most tractography and plotting
  functions take a single `options` struct in `varargin{1}`. The visualization
  entry points (`visualizeTractography`, `visualizeTractographySlices`,
  `nim_plot`) are the exception and use `inputParser` name-value pairs.

---

## 1. Entry Points

### `main(imgpath, nimpath, varargin)`

**File**: `main.m`

Builds the processed `nim`: data loading, preprocessing, tensor fitting,
eigendecomposition, FA, registration, parcellation and tissue segmentation.

| Argument | Type | Required | Description |
|---|---|---|---|
| `imgpath` | char/string | Yes | Path to the NIfTI image, without extension |
| `nimpath` | char/string | Yes | Path to save the processed `.mat` (must end in `.mat`) |

The trailing arguments are **positional and type-dispatched**, not name-value:

| Position | Accepted | Effect |
|---|---|---|
| 3rd | struct with a `preprocessing` field | YAML config; `t1_file` is derived from `preprocessing.t1_available` |
| 3rd | struct with a `run_dir` field | `run_info`; redirects all writes into the run directory |
| 3rd | struct (other) | Legacy options struct; `.t1_file` honored |
| 3rd | char/string | Legacy T1 file path |
| 4th | struct with a `run_dir` field | `run_info`, when the 3rd argument was a config |

**Returns**: nothing; writes the `nim` struct to `nimpath`.

**Behavior**:

- Short-circuits if the output `.mat` already exists, and skips preprocessing if
  the preprocessed `.nii.gz` exists. Delete the target to force reprocessing.
- Detects raw data by `_raw` in the filename and triggers preprocessing.
- Runs `nim_read` → `nim_dt_spd` → `nim_eig` → `nim_fa` → optional
  `nim_registration` → parcellation → mask improvement → tissue segmentation →
  `nim_save`.
- Adds every `src/nim_*` directory and `lib/spm12` to the MATLAB path, and
  `char()`-converts path inputs.

```matlab
main('data/sample', 'output/sample.mat');
main('data/sample_raw', 'output/sample.mat', 'data/sample_T1.nii.gz');

config   = load_config_yaml('config/ismrm2015.yml');
run_info = create_run_directory('config/ismrm2015.yml');
main('data/ismrm2015/ismrm2015', 'ismrm2015.mat', config, run_info);
```

---

### `runTractography(data_path, varargin)`

**File**: `runTractography.m`

Loads a processed `nim`, resolves the seed mask and ACT masks, dispatches to one
of the three trackers, applies ROI filtering and output decimation, and saves.

| Argument | Type | Required | Description |
|---|---|---|---|
| `data_path` | char | Yes | Path to the `.mat` holding the `nim` struct |

Again positional and type-dispatched:

| Position | Accepted | Effect |
|---|---|---|
| 2nd | struct with a `tractography` field | YAML config; the tracker comes from `tractography.algorithm` |
| 2nd | struct with a `run_dir` field | `run_info` |
| 2nd | `'standard'` / `'hinec'` / `'mmf'` | Legacy algorithm selector |
| 2nd | `'IronTract'` | IronTract submission mode; requires an injection file and output directory as the 3rd and 4th arguments |
| 3rd | struct with a `run_dir` field | `run_info`, when the 2nd argument was a config |

**Returns**: nothing.

**Output**: `tracks_{algorithm}_{yyyy-mm-dd_HH_MM_SS}.mat`, written to
`run_info.tractography_dir` when a `run_info` is supplied and to
`tractography_results/` otherwise. The file contains:

| Variable | Description |
|---|---|
| `tracks` | Cell array of $N \times 3$ voxel-coordinate polylines |
| `options` | The flat option struct actually handed to the tracker |
| `elapsed_time` | Tracking wall time, seconds |
| `algorithm` | `'standard'`, `'hinec'` or `'mmf'` |
| `track_meta` | `.seed_index`, `.seed_points`, `.n_seeds` — populated by `hinec`, empty otherwise |

**Seed mask priority** (first match wins, then intersected with the brain mask
where available and with `seeding.fa_min`):

1. `seeding.roi` — explicit atlas regions, via `nim_roi_mask`
2. `nim.mask` — the preprocessed brain mask
3. `nim.parcellation_mask`, dilated by 3 voxels
4. `nim.FA > 0.10` — last-resort fallback

**CSD provisioning**: with `field: csd`, FOD peaks are computed by `nim_csd`
before dispatch and cached alongside the source nim as `<source>_csd.mat`, so
they are computed once per dataset and reused by every subsequent config.

```matlab
runTractography('output/sample.mat');
runTractography('output/sample.mat', 'hinec');
runTractography('output/sample.mat', config, run_info);
```

---

### `runhinec`

**File**: `runhinec.m`

Demo **script** (not a function). Loads `sample_parcellated.mat` from the working
directory and calls `nim_plot(nim, 'mode', 'parcels')`.

---

## 2. Configuration (`src/nim_utils/`)

`nim_config_schema` is the single source of truth for the configuration surface.
Every other function in this section derives from it, which is what keeps the
loader, the CLI overrides and the generated reference from drifting apart.

### `nim_config_schema()`

**Returns**: a struct array, one entry per configurable parameter, with fields
`path` (canonical dotted path), `default`, `type`
(`numeric`/`string`/`logical`/`list`), `allowed` (enum values), `range`
(`[lo hi]` for numerics), `legacy` (superseded flat key names that migrate here),
`algos` (which trackers read it) and `help`.

Add or change a knob here, not in the trackers.

### `load_config_yaml(config_file)`

Parses a YAML config, applies schema defaults, validates types, enums and ranges,
migrates legacy flat keys, and rejects unknown keys. Returns a nested `config`
struct with `.preprocessing` and `.tractography` sections.

### `nim_yaml_parse(filename)`

Indentation-aware reader for the YAML subset HINEC uses: `section:` at column 0,
`key: value` one level in, and `group:` / `key: value` two levels in. A third
level of nesting is a hard error. Handles numbers, booleans, `null`/`~`, quoted
and unquoted strings, inline lists `[a, b, c]`, and `#` comments outside quotes.
Returns the raw nested struct with no defaults and no validation.

### `nim_config_to_options(config)`

Flattens the canonical nested config into the flat option names the trackers
read. This is the one place the legacy names live. Notable derivations:

- `integrator.method` → `options.integration_order` (`euler`→1, `rk2`→2,
  `rk4`→4, `rkf45`→5) — a method selector wearing a numeric name.
- `options.adaptive_step` is true only for `rkf45` with `integrator.adaptive`.
- `options.max_steps = ceil(termination.max_arc / integrator.step)`, so refining
  the step cannot truncate tracks.
- Warns when `angle_max × step ≥ 90`, where the angle criterion is inert.
- `sel_power` and `fa_threshold` are deliberately **not** emitted (see
  `nim_config_retired`).

### `nim_config_apply_overrides(config, sets)`

Applies CLI `--set key=value` overrides. Returns `[config, applied]`. Keys resolve
in order: full canonical path (`tractography.integrator.step`), path with
`tractography.` assumed (`integrator.step`), unique bare leaf name (`upsample`),
then legacy alias (`step_size`, with a warning). An ambiguous bare leaf or an
unknown key is an error, never a silent no-op. Values parse as YAML scalars plus
lists (`[41,42]`, `41,42`, `SLF_L,SLF_R`).

### `nim_config_write(config, out_file, opts)`

Writes a config back out as nested YAML. By default only values differing from
the schema default are emitted; `opts.all = true` writes every key (for run
manifests), and `opts.header` is a cell array of comment lines for the top.

### `nim_config_docs(out_file)`

Generates the configuration reference from the schema. `nim_config_docs()`
returns the markdown; `nim_config_docs('docs/YAML_CONFIG.md')` writes it.

!!! warning "`docs/YAML_CONFIG.md` is generated"
    Do not edit it by hand — regenerate it from the schema instead.

### `nim_config_retired()`

Returns the keys that were once accepted and now do nothing, each with the reason.
The loader warns and drops them rather than silently remapping them, since they
had no functional effect and mapping them anywhere would change behavior.

### `create_run_directory(config_file, varargin)`

Creates a timestamped, self-contained run directory.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `config_file` | char | — | YAML config used for this run (copied in) |
| `'base_dir'` | char | `'hinec_runs'` | Parent directory |
| `'run_name'` | char | auto | Custom run name |
| `'description'` | char | `''` | Recorded in the run metadata |

**Returns**: `run_info` with `.run_dir`, `.config_file`, `.logs_dir`,
`.intermediate_dir`, `.output_dir`, `.tractography_dir`, `.diagnostics_dir`,
`.run_id`, `.timestamp`. Creates the subdirectory tree, copies the config,
records git information and updates a `latest` symlink.

Code that writes outputs branches on `use_run_dir = ~isempty(fieldnames(run_info))`;
preserve that branch when extending it.

---

## 3. Data I/O (`src/nim_utils/`)

### `nim_read(datapath, opts)`

Reads NIfTI-1 diffusion data together with its acquisition parameters. `opts` are
name-value arguments.

| Option | Type | Default | Description |
|---|---|---|---|
| `BvalPath`, `BvecPath`, `MaskPath` | string | derived from `datapath` | Explicit overrides |
| `Mask`, `Bval`, `Bvec` | `"on"`/`"off"` | `"on"` | Whether to read each sidecar |
| `B0Threshold` | double | `5` | Volumes with `bval < B0Threshold` are b0 |

**Returns** `nim` with:

| Field | Description |
|---|---|
| `.img` | 4D image data, `[X Y Z N]` |
| `.img_b0` | Mean of the b0 volumes |
| `.img_bi` | The diffusion-weighted volumes, as double |
| `.bval`, `.bvec` | B-values and gradient directions |
| `.hdr` | `niftiinfo` header |
| `.mask` | Brain mask, if `{name}_M.nii.gz` exists |
| `.xdim`, `.ydim`, `.zdim`, `.size3` | Volume dimensions and their product |
| `.size_b0`, `.size_bi` | Number of b0 and diffusion-weighted volumes |
| `.thrsh_b0` | The b0 threshold actually used |

Both the standard single-line `.bval`/`.bvec` layout and the IronTract
one-value-per-line layout are detected automatically.

### `nim_save(nim, nimpath)`

Saves the `nim` struct. Uses MATLAB v7.3 (HDF5) when the struct exceeds 1900 MB,
which is the practical limit of the default format, and the default format below
that.

### `nim_load_nim(file_prefix)`

Loads a previously saved `nim` together with its mask, parcellation mask and atlas
labels.

### `nim_load_labels(nim)`

Populates `nim.atlas_labels` from the atlas referenced by
`nim.parcellation_mask_file`, searching for FSL XML and `.mat` label files.

### `nim_load_atlas_labels(atlas_type)`

Reads an FSL atlas XML directly. `atlas_type` is `'HarvardOxford'`, `'JHU'` or
`'JHU-tract'`. Returns a struct with `.map` (a `containers.Map`) and `.atlas_type`.

### `nim_atlas_label_map(nim)`

Returns the atlas index → region name map for a `nim`, resolving in order: a map
embedded in the nim, the sidecar `<prefix>_atlas_labels.mat` if recorded, then
`nim_load_atlas_labels`. Errors with an actionable message rather than returning
an empty map.

### `nim_nifti_affine(f)`

Returns the $4 \times 4$ voxel-to-world affine of a NIfTI-1 file, mapping
**0-based** voxel indices to world RAS millimetres. Prefers the sform and falls
back to the qform, since roughly a third of the ISMRM 2015 ground-truth ROI masks
set `sform_code = 0` and describe their geometry only through the qform
quaternion; reading `srow` blindly returns zeros for those and silently places
them at the origin. Read from the raw header rather than through a toolbox,
because toolboxes disagree on whether their affine is 0- or 1-based.

### `nim_read_trk(trk_file, target_nii)`

Reads a TrackVis `.trk` file into a cell array of $N \times 3$ tracks, plus the
header. With `target_nii`, streamlines are resampled into that volume's voxel grid
via world RAS, which is the only frame the two files agree on:

$$
\text{ras} = \mathbf{M}_{\text{vox}\to\text{ras}}\,[\,p / \text{voxel size}\,],
\qquad
\text{vox}_{\text{target}} = \mathbf{S}^{-1}\,\text{ras} + 1 .
$$

Without it, points come back in the TRK's own voxel units. The transform matters:
the ISMRM 2015 ground truth is 1 mm $180\times216\times180$ while this pipeline's
DWI is 2 mm $90\times108\times90$, and plotting the two together untransformed
draws the ground truth at half scale in a corner of the volume.

---

## 4. Calculation (`src/nim_calculation/`)

### `nim_dt_spd(nim, opts)`

Tensor fitting with an SPD constraint.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nim` | struct | — | Needs `.img`, `.bval`, `.bvec` |
| `opts.Mask` | string | `"on"` | Restrict the fit to the brain mask |
| `opts.Steps` | int | 2000 | BFGS iteration cap |
| `opts.FontSize` | int | 16 | Diagnostic plot font size |

**Returns** `nim` with `.DT` `[X Y Z 6]`, `.evec` `[X Y Z 3 3]` and `.eval`
`[X Y Z 3]` sorted descending.

Per voxel: compute the ADC $Y = \ln(S_0/S)/b$, build the design matrix $H$, solve
$D = H \backslash Y$, and if any eigenvalue is negative re-solve by BFGS over a
Cholesky parameterization (`lib/bfgs/vox_dt_bfgs.m`), which cannot produce a
non-SPD tensor. NaN, Inf and unrealistic eigenvalues (> 0.01) are rejected;
eigenvalues are then sorted with `maxk`.

### `nim_dt(nim, opts)`

Plain least-squares tensor fitting with no SPD enforcement. Returns `nim` with
`.DT` and reports the count of non-positive-definite tensors without correcting
them.

### `nim_eig(nim, opts)`

Eigendecomposition of `.DT`. Returns `.evec` (columns ordered by eigenvalue) and
`.eval` (descending).

### `nim_fa(nim, opts)`

Fractional anisotropy from `.eval`. Returns `.FA` in $[0, 1]$, with FA set to 0
where all eigenvalues vanish.

### `nim_csd(nim, options)`

Single-shell, single-tissue constrained spherical deconvolution (Tournier 2007),
implemented natively in MATLAB: real spherical-harmonic fit of the DW signal, a
single-fiber response estimated from high-FA voxels, iterative non-negativity
deconvolution, then peak extraction.

| Option | Default | Description |
|---|---|---|
| `lmax` | 6 | SH order, capped by the number of directions |
| `sf_fa` | 0.7 | FA threshold for response-function voxels |
| `fod_lambda` | 1.0 | Non-negativity regularisation weight |
| `fod_thresh` | 0.1 × mean | FOD amplitude treated as negative |
| `n_iter` | 50 | Deconvolution iterations |
| `peak_thresh` | 0.1 | Minimum peak amplitude, as a fraction of the max FOD |
| `peak_min_sep` | 25 | Minimum angular peak separation, degrees |
| `max_peaks` | 3 | Peaks kept per voxel |
| `response` | — | Precomputed zonal response, to skip estimation |

**Returns** `nim` with `.peaks`, `.npeaks`, `.peak_w`, `.fod_sh` and `.response`
— the same interface the trackers consume for DTI, so the two direction sources
are interchangeable upstream of tractography.

!!! note "Defaults differ between `nim_csd` and the config"
    `runTractography` passes the schema's `csd.*` values, which are stricter than
    `nim_csd`'s own defaults: `peak_thresh` 0.5 and `peak_min_sep` 45.

### `nim_mmf_geometry(nim, options)`

Builds the moving-frame geometry over the whole domain and stores it in the `nim`
(Chun & Peng, in preparation). Constructs the denoised tangent field $e_1$, the
Frenet normal $e_2 = (de_1/ds)/\lVert de_1/ds \rVert$ with a reference-axis
fallback where the curvature vanishes, $e_3 = e_1 \times e_2$, and the connection
1-form of the frame field.

| Option | Default | Description |
|---|---|---|
| `frame_sel_power` | 16 | Trajectory-dependent denoising selectivity |
| `field` | `'dti'` | `'csd'` builds the per-peak connection from FOD peaks instead |

**Returns** `nim` with `.mmf_frames` `[X Y Z 3 3]`, `.mmf_kappa` `[X Y Z 3]`,
`.mmf_tau` `[X Y Z]`, `.mmf_built` and `.mmf_field` (stamped from what was
actually built, so the tracer can tell DTI geometry from CSD geometry).

### `nim_build_frames(nim, options)` and `nim_connection_form(frames, mask, options)`

The two components `nim_mmf_geometry` composes: the frame field, and the
connection 1-form derived from it.

---

## 5. Preprocessing (`src/nim_preprocessing/`)

### `nim_preprocessing(file_prefix, varargin)`

FSL-driven preprocessing orchestrator. Takes an options struct (or the legacy
positional `run_eddy`, `atlas_type` form) and writes its outputs to disk; it
returns nothing.

| Option | Default | Description |
|---|---|---|
| `.run_denoising` | true | Denoise before tensor fitting |
| `.denoise_method` | `'dwidenoise'` | `'dwidenoise'`, `'nlmeans'` or `'gaussian'` |
| `.run_motion_correction` | true | MCFLIRT motion correction |
| `.run_eddy` | true | FSL `eddy` |
| `.improve_mask` | true | FA-based mask refinement |
| `.atlas_type` | `'HarvardOxford'` | Atlas for parcellation |
| `.phase_encoding_direction` | `""` | Phase-encoding axis for `eddy`, e.g. `'j-'` |
| `.total_readout_time` | `[]` | Readout time, seconds |
| `.eddy_index_vector` | `[]` | Per-volume acquisition indices |

### Individual preprocessing functions

| Function | Signature | FSL tool | Purpose |
|---|---|---|---|
| `preproc_denoising` | `(dwi_file, file_prefix, method)` | dwidenoise / SUSAN | Remove thermal noise |
| `preproc_extract_b0` | `(dwi_file, output_dir)` | fslroi | Extract the b0 reference |
| `preproc_motion_correction` | `(dwi_file, bvec_file, bval_file, file_prefix)` | mcflirt | Correct head motion |
| `preproc_eddy_correction` | `(dwi_file, brain_mask_file, bvec_file, bval_file, acqp_file, index_file, file_prefix)` | eddy | Correct eddy currents |
| `preproc_brain_extraction` | `(dwi_or_b0_file, output_dir, brain_mask_file)` | BET (`-f 0.3`) | Generate the brain mask |
| `preproc_t1_brain_extraction` | `(t1_file, output_dir)` | BET (`-f 0.4 -R`) | Extract the T1 brain |
| `preproc_fieldmap_correction` | `(dwi_file, fieldmap_file, output_prefix, varargin)` | FUGUE | Correct field inhomogeneity |
| `preproc_mask_improvement` | `(brain_mask_file, fa_data, file_prefix)` | BET + fslmaths | FA-based mask refinement |
| `preproc_tissue_segmentation` | `(fa_data, brain_mask_file, file_prefix, t1_file)` | fslmaths | WM/GM/CSF classification |
| `preproc_create_dwi_reference` | `(dwi_file, brain_mask_file, output_dir, file_prefix)` | fslroi + fslmaths | Clean DWI reference volume |
| `preproc_t1_dwi_registration` | `(dwi_ref_file, t1_file, t1_brain_file, output_dir, file_prefix)` | epi_reg | T1 ↔ DWI transform |
| `preproc_t1_mni_registration` | `(t1_file, t1_brain_file, output_dir, file_prefix)` | FLIRT + FNIRT | T1 → MNI warp |
| `preproc_atlas_resampling` | `(reference_file, output_dir, file_prefix, atlas_type, options)` | applywarp / flirt | Atlas into DWI space |
| `preproc_cleanup` | `(output_dir, file_prefix, keep_files)` | — | Remove intermediate files |

Returns worth noting: `preproc_tissue_segmentation` returns
`[wm_mask_file, gm_mask_file, csf_mask_file]`, and `preproc_atlas_resampling`
returns `[parcellation_mask_output, atlas_labels_file]`.

**Tissue thresholds**: WM `FA > 0.2` (then eroded), GM `0.05 < FA ≤ 0.2`, CSF
`FA ≤ 0.05`.

### `preproc_fieldmap_correction` optional parameters

| Name | Type | Default | Description |
|---|---|---|---|
| `'TE'` | double | auto-detect | Echo time difference |
| `'dwell_time'` | double | 0.00058 | Effective echo spacing, seconds |
| `'phase_dir'` | char | `'y'` | Phase-encoding direction |
| `'units'` | char | auto-detect | `'Hz'` or `'rad/s'` |
| `'smooth_sigma'` | double | 2.0 | Field map smoothing, mm |
| `'mask_file'` | char | auto-generate | Brain mask path |

---

## 6. Registration (`src/nim_registration/`)

### `nim_registration(dwi_file, t1_file, varargin)`

Multi-modal registration orchestrator. Returns a `registration_data` struct
carrying `.transforms`, `.registered_images` and `.quality_metrics`.

### `register_dti_to_t1(registration_data, options)`

DTI → T1. FSL backend uses FLIRT with 6 DOF and a correlation-ratio cost; an SPM
backend is also available. Returns the updated `registration_data`.

### `register_t1_to_mni(registration_data, options)`

T1 → MNI152 in two stages: linear FLIRT (12 DOF), then optional non-linear FNIRT.
Returns the updated `registration_data`, including inverse warps.

### `nim_apply_transforms(input_file, registration_data, transform_chain, varargin)`

Applies a chain of transforms between coordinate spaces.

| Parameter | Type | Description |
|---|---|---|
| `input_file` | char | Image to transform |
| `registration_data` | struct | Transforms from the registration pipeline |
| `transform_chain` | cell | Sequence, e.g. `{'mni_to_t1', 't1_to_dti'}` |
| `options` | struct | `.output_file`, `.interpolation` (`'linear'`/`'nearest'`/`'spline'`), `.method` |

**Returns**: `output_file`, the path to the transformed image.

### Supporting functions

| Function | Signature | Purpose |
|---|---|---|
| `extract_reference_volumes` | `(registration_data, options)` | Extract the b0 volume from the DWI |
| `compute_registration_quality` | `(registration_data, options)` | Score alignment |
| `compute_normalized_mutual_information` | `(img1, img2)` | NMI between two volumes |
| `save_registration_data` | `(registration_data, options)` | Persist the results |
| `generate_registration_report` | `(registration_data, options)` | HTML quality report |

---

## 7. Parcellation and ROIs

### `nim_parcellation(nim, parcellation_file, varargin)`

**File**: `src/nim_parcellation/nim_parcellation.m`

Loads a parcellation NIfTI, or generates one from an MNI atlas if the file is
absent. `varargin{1}` is an options struct with `.dwi_file` (required for
generation), `.atlas_type` (`'jhu'` default, plus `'HarvardOxford'` and `'aal'`)
and `.force_regenerate`. Returns `nim` with `.parcellation_mask`,
`.parcellation_mask_file` and `.labels`, handling dimension mismatches
automatically.

### `nim_parcellation_registered(nim, registration_data, parcellation_mask_file)`

**File**: `src/nim_parcellation/nim_parcellation_registered.m`

Parcellation through the full registration chain (MNI → T1 → DWI). Returns `nim`
with `.parcellation_mask` and `.parcellation_info` (quality metrics, atlas type,
coverage).

### `nim_parcellation_from_masks(mask_dirs, ref_nii, target_dims)`

**File**: `src/nim_utils/nim_parcellation_from_masks.m`

Builds a parcellation from a directory of binary mask volumes, resampled
nearest-neighbor through the two sforms into the target grid.

**Returns** a struct with `.labels` (integer label volume), `.map` (label →
name), `.masks` (name → logical volume, overlaps intact) and `.overlap` (a report
on how much the source masks overlap).

Both `.labels` and `.masks` exist because anatomical bundles genuinely overlap —
a temporal-stem voxel can belong to the uncinate, the inferior longitudinal
fasciculus and the corpus callosum at once — and an integer label volume cannot
represent that. Overlap in `.labels` is resolved smallest-region-wins, so a large
structure cannot swallow every small bundle crossing it.

### `nim_roi_mask(nim, roi_spec, dilate)`

**File**: `src/nim_utils/nim_roi_mask.m`

Builds a binary mask from a list of atlas regions. `roi_spec` freely mixes numeric
indices and names; names match case-insensitively by exact match, then unique
substring, then short alias (`SLF_L`, `CST_R`, `CC_genu`). `dilate` (default 0)
grows the result with a spherical structuring element — often necessary rather
than cosmetic, since the JHU parcellation is thin (9,641 labelled voxels in the
ISMRM data, of which the uncinate is 23–24).

**Returns** `[mask, info]`, where `info` carries per-label voxel counts at each
stage so an empty result can be diagnosed rather than guessed at.

---

## 8. Tractography (`src/nim_tractography/`)

All three trackers take a `nim` struct or a path to one, plus a single flat
`options` struct. In normal use `runTractography` produces that struct from the
config via `nim_config_to_options`; the tables below give both the config key you
set and the option field the tracker reads.

### `nim_tractography_standard(data_path, varargin)`

FACT — Fiber Assignment by Continuous Tracking (Mori et al., 1999). Discrete
per-voxel tensors with no interpolation; each step is an exact ray-box
intersection with the current voxel boundary, so the step length is determined by
the geometry rather than configured.

| Config key | Option field | Default | Description |
|---|---|---|---|
| `seeding.density` | `seed_density` | 1 | Seeds per voxel |
| `seeding.strategy` | `seed_strategy` | `"uniform"` | Seed placement |
| `termination.fa_min` | `termination_fa` | 0.05 | FA floor for termination |
| `termination.angle_max` | `angle_thresh` | 60 | Turn budget, **degrees per voxel of arc** |
| `termination.max_arc` | `max_steps` | 5000 | Derived step cap (boundary crossings) |
| `termination.min_arc` | `min_length` | 10 | Minimum track length |
| — | `seed_mask` | `[]` | Seed mask, supplied by `runTractography` |
| — | `step_size` | 0.5 | Unused: exits are computed exactly |
| — | `use_fa_seed_filter` | false | Filter seeds by FA |
| `diagnostics` | `enable_diagnostics` | true | Timing diagnostics |

**Returns**: `tracks`, a cell array of $N \times 3$ matrices.

### `nim_tractography_hinec(data_path, varargin)`

Interpolated streamline tractography: a continuous direction field, Runge-Kutta
integration and optional ACT. **Returns** `[tracks, meta]`, where `meta` carries
`.seed_index`, `.seed_points` and `.n_seeds` — the seed correspondence that makes
per-streamline comparison across a refinement ladder possible, since the tracker
compacts out seeds that produce nothing.

Everything in the standard table applies, plus:

| Config key | Option field | Default | Description |
|---|---|---|---|
| `field` | `field` | `'dti'` | `'dti'` (principal eigenvector) or `'csd'` (FOD peaks) |
| `interpolation.method` | `interp_method` | `'trilinear'` | `'trilinear'` ($C^0$), `'cubic'` (Keys convolution, $C^1$), `'spline'` ($C^2$) |
| `interpolation.upsample` | `upsample` | 1 | Sample the field at $1/u$ voxel spacing before building interpolants |
| `integrator.method` | `integration_order` | `'rk4'` → 4 | `euler`→1, `rk2`→2, `rk4`→4, `rkf45`→5 |
| `integrator.step` | `step_size` | 0.2 | Step in voxels; the initial step for `rkf45` |
| `integrator.adaptive` | `adaptive_step` | true (rkf45 only) | Adaptive step-size control |
| `integrator.tolerance` | `rkf_tolerance` | 0.01 | RKF45 error tolerance, voxels |
| `integrator.safety` | `rkf_safety` | 0.9 | RKF45 step safety factor |
| `integrator.step_min` | `step_min` | 0.01 | RKF45 minimum step |
| `integrator.step_max` | `step_max` | 1.0 | RKF45 maximum step |
| `act` | `act_enabled` | **false** | Anatomically constrained tracking |
| — | `wm_mask`, `gm_mask`, `csf_mask` | `[]` | Tissue masks; supplied by `runTractography` when ACT is on |

!!! note "Two different ACT defaults"
    `tractography.act` defaults to **false** in the config schema. The tracker's
    own internal default for `act_enabled`, used only when it is called directly
    with a hand-built options struct, is `true`. ACT is additionally inert
    without WM/GM/CSF masks in the `nim`.

**Precomputation**: `griddedInterpolant` objects for FA and for the six unique
components of the dyadic $\mathbf{v}_1\mathbf{v}_1^\top$ are built once per run,
with `'none'` extrapolation so out-of-volume queries return `NaN` and terminate
the streamline. This is 2–5× faster than repeated `interp3` calls. The dyadic,
not the signed eigenvector components, is what DTI tracking interpolates — see
[MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md).

### `nim_tractography_mmf_connframe(nim, options)`

Method of Moving Frames tracking (Chun & Peng, in preparation). Reads the geometry
built by `nim_mmf_geometry`, evolves the carried frame by the structure equation
and advances $d\mathbf{x}/ds = e_1$. The direction comes purely from the
connection form; the integrator is only the stepping scheme. With `field: 'csd'`,
a per-peak connection yields multiple pathways through a crossing.

| Config key | Option field | Default | Description |
|---|---|---|---|
| `field` | `field` | `'dti'` | `'dti'` or `'csd'` (per-peak, multi-frame) |
| `integrator.method` | `integrator` | `'rk4'` | `'rk4'` fixed step, or `'rkf45'` adaptive Dormand-Prince |
| `integrator.tolerance` | `rkf_tol` | 0.02 | RKF45 error tolerance |
| `integrator.step` | `step_size` | 0.2 | Step in voxels |
| `integrator.step_min` / `step_max` | `step_min` / `step_max` | 0.02 / 0.5 | Adaptive step bounds |
| `interpolation.method` | `interp_method` | `'trilinear'` | How $\kappa$, $\tau$ and $e_1$ are sampled; peak *directions* always use nearest neighbor, since sign ambiguity breaks under averaging |
| `mmf.anchor` | `mmf_anchor` | 0 | Re-anchor strength of $e_1$ toward the field tangent; 0 is the pure connection-form evolution |
| `mmf.frame_sel_power` | `frame_sel_power` | 16 | Selectivity used when building the frames |

Requires `nim.evec`, `.FA` and `.mask`; `field: 'csd'` additionally requires
`nim.peaks` and `nim.npeaks`. The geometry is rebuilt in place if it is absent or
was baked for a different field. **Returns** `[tracks, info]`.

### `nim_tractography_highorder(data_path, varargin)`

Superseded predecessor of the `hinec` tracker, kept for reference. Options
`order` (3), `step_size` (0.2), `fa_threshold` (0.1), `angle_thresh` (60),
`max_steps` (5000), `min_length` (10), `termination_fa` (0.05), `interp_method`
(`"spline"`), `seed_density` (1). Not reachable from `algorithm:` — new work
should use `nim_tractography_hinec`.

### Tractography support functions

#### `nim_seed_offsets(density)`

Returns `[offsets, info]`: exactly `density` sub-voxel seed offsets in
$[-0.5, +0.5]^3$. A perfect cube uses the full $n^3$ lattice in `ndgrid` order;
any other count takes a deterministic farthest-point subset of the next-larger
lattice, so seeds stay spread through the voxel rather than clustering in a
corner, with no RNG and therefore reproducible. `info` reports `.requested`,
`.actual`, `.per_axis`, `.lattice_points`, `.exact_lattice` and `.spacing`.

#### `nim_angle_limit(angle_rate, arc)`

Returns `[allowed_deg, cos_allowed]`, the turn budget for one step.
`angle_rate` is degrees of turning per voxel of arc — a minimum radius of
curvature $R = (180/\pi)/\texttt{angle\_rate}$ voxels in disguise — and `arc` is
the arc the step advances.

Pass the **nominal** step size, never the realised chord. Under
$d\mathbf{x}/ds = \mathbf{v}$ with $|\mathbf{v}| = 1$ a step advances `arc` of arc
by construction; the chord is shorter mostly because the RK stage vectors
disagree, not because less arc was covered, so scaling by it would make the
criterion method-dependent — Euler and RK2 have a chord of exactly $h$ while
RK4's falls as low as $0.25h$ on this data.

Because $\mathbf{v}_1$ is a line field and tangents are sign-aligned before the
turn is measured, the measured turn never exceeds 90°: a budget at or above that
is inert, not merely loose.

#### `nim_principal_dir(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz)`

Unit principal eigenvector of a symmetric $3\times3$, in closed form; `[]` for a
degenerate or non-finite matrix. The sign is arbitrary — what it returns is a
line, and the caller picks the representative. Cardano's trigonometric solution
via the deviatoric part, then the eigenvector from the null space of
$\mathbf{D} - \lambda\mathbf{I}$, taken as the largest cross product of its rows
so the result stays conditioned when a row is near zero. Measured at 1.7 µs
against 2.5 µs for `eig` on the same matrix, with a worst disagreement of
1.1e-13 over 20,000 random PSD matrices.

#### `nim_resample_track_arc(tracks, arc_step)`

Resamples streamlines to a fixed arc-length spacing, decoupling the stored
representation from the integration step — every step is otherwise kept, so
storage grows as $1/h$. Endpoints are preserved exactly, since they carry the
connectivity that scoring is built on. `arc_step ≤ 0` returns the input
untouched. Accepts a single $N \times 3$ track or a cell array.

#### `nim_filter_tracks_roi(tracks, nim, options)`

Keeps or discards streamlines by the regions they touch. Returns
`[tracks, stats]`, with `stats.keep` a logical index into the input so
`track_meta` can be kept aligned.

| Option field | Config key | Description |
|---|---|---|
| `include_roi` | `filter.include_roi` | Keep tracks touching these regions (waypoint) |
| `exclude_roi` | `filter.exclude_roi` | Discard tracks touching these regions |
| `roi_filter_mode` | `filter.mode` | `'all'` (every include region) or `'any'` |
| `roi_filter_dilate` | `filter.roi_dilate` | Dilate the filter masks first |
| `endpoints_in` | `filter.endpoints_in` | Two regions; one end in each, either way round |
| `contained_in` | `filter.contained_in` | Every point must lie inside these regions |

With none of these set, the tracks are returned untouched. `endpoints_in` and
`contained_in` together are how the ISMRM 2015 scorer defines a bundle: an
endpoint pair plus a containment corridor. Neither is expressible as a waypoint
include, which is why they are separate keys.

#### `nim_sift(tracks, nim, options)`

SIFT-style filtering with per-lobe attribution (an analog of Smith 2013). Each
streamline segment is attributed to the FOD peak it best aligns with, and
streamlines running through over-represented lobes are removed, so that
streamline density on each lobe matches that lobe's FOD amplitude. This separates
a legitimately dense bundle from spray overloading a weak lobe, which a
per-voxel density match cannot. Requires `nim.peaks` and `nim.peak_w` from CSD.
Options: `sift_batch_frac` (0.03), `sift_n_iter` (60), `sift_min_keep` (0.20),
`sift_align_pow` (1). Returns `[tracks_out, info]`.

---

## 9. Solution Verification (`src/nim_utils/`)

### `nim_convergence_error(test_file, ref_file, opts)`

Per-streamline geometric error of one tractography run against a reference run.
Streamlines are matched by `track_meta.seed_index`, not by position in the cell
array, since the tracker compacts out unproductive seeds and "streamline 7" is a
different fiber in the two runs.

Two error families, and the distinction governs which one an order may be fitted
to:

- **`prefix`** — compared over $\pm$ `opts.prefix_arc` voxels of arc *either side
  of the seed* (default 20). Anchoring at the seed matters: a track is stored as
  `[reversed backward half; seed; forward half]`, so its first point is a
  termination point. This is the family to fit an observed order to.
- **`whole`** — the full common arc, plus endpoint distance and arc-length
  difference. Dominated by termination, which is discontinuous in the step size.
  Characterises termination sensitivity; not an order.

Statistics are reported as median and p95, never a mean: crossing regions make
the distribution heavy-tailed. Returned fields include `.n_test`, `.n_ref`,
`.n_matched`, `.matched_frac`, `.seed_index`, `.prefix`, `.prefix_max`,
`.whole`, `.whole_max`, `.endpoint` and `.length_diff`.

See [CONVERGENCE.md](CONVERGENCE.md) for the measured rates.

---

## 10. Tractography Plotting (`src/nim_tractography/`)

These take an options **struct** in `varargin{1}`.

### `nim_plot_tractography(tracks, nim, varargin)`

| Field | Default | Description |
|---|---|---|
| `color_mode` | `"direction"` | `'direction'`, `'fa'`, `'parcellation'` |
| `show_anatomy` | true | FA slice background |
| `track_subset` | all | Indices of tracks to draw |
| `line_width` | 1 | Line width |
| `alpha` | 0.8 | Transparency |

Caps the display at 2000 tracks.

### `nim_plot_tractography_region(tracks, nim, region_id, varargin)`

Draws the tracks associated with one parcellation region. `nim` must carry
`.parcellation_mask`.

| Field | Default | Description |
|---|---|---|
| `filter_mode` | `'pass_through'` | `'pass_through'`, `'start_in'`, `'end_in'`, `'connect_to'` |
| `min_overlap` | 0.3 | Minimum fraction of the track inside the region |
| `show_region` | true | Region isosurface overlay |
| `region_alpha` | 0.3 | Overlay transparency |
| `track_color` | `'direction'` | `'direction'`, `'fa'`, `'uniform'`, `'region'` |
| `max_tracks` | `inf` | Display cap |
| `line_width` | 1.2 | Line width |
| `min_track_length` | 10 | Minimum track points |
| `show_stats` | true | Print statistics |

### `nim_plot_connectivity_matrix(tracks, nim, varargin)`

| Field | Default | Description |
|---|---|---|
| `min_track_length` | 10 | Minimum track points |
| `normalize` | true | Scale to $[0, 1]$ |
| `symmetric` | true | Symmetrise |

**Returns** the $N \times N$ connectivity matrix, and draws a four-panel figure:
heatmap, histogram, node strengths, summary.

### `nim_plot_vector_field(nim, varargin)`

| Field | Default | Description |
|---|---|---|
| `slice` | middle | Slice index |
| `downsample` | 2 | Downsample factor |
| `scale` | 1 | Vector scale |
| `axis_view` | `"axial"` | `'axial'`, `'coronal'`, `'sagittal'` |

Masks vectors by `FA > 0.2` and overlays them on an FA background.

### `nim_excitation_time_map(nim, seed_points, varargin)`

Excitation propagation time from `seed_points` ($N \times 3$), with velocity
shaped by DTI anisotropy.

| Field | Default | Description |
|---|---|---|
| `conduction_speed` | 1.0 | Base conduction speed |
| `fa_scaling` | true | Scale velocity by FA |
| `max_time` | 100 | Time cap |
| `method` | `"fast_marching"` | `'fast_marching'` or `'dijkstra'` |

**Returns** `[time_map, velocity_field]`.

---

## 11. Visualization (`src/nim_visualization/`)

### `visualizeTractography(varargin)`

Unified viewer with four modes and figure export. The first arguments are either
a run directory (files auto-detected, wildcards accepted) or an explicit
`tracks_file` and `nim_file` pair, followed by name-value options.

| Name | Default | Description |
|---|---|---|
| `'mode'` | `'whole'` | `'whole'`, `'region'`, `'grid'`, `'sequential'` |
| `'region'` | `[]` | Region ID or array of IDs |
| `'export_dir'` | `''` | Directory for image export |
| `'export_format'` | `'png'` | `'png'`, `'pdf'`, `'eps'`, `'fig'` |
| `'export_dpi'` | 300 | Export resolution |
| `'filter_mode'` | `'pass_through'` | `'pass_through'`, `'start_in'`, `'end_in'`, `'connect_to'` |
| `'region_filter'` | `'pass_through'` | Slice-viewer-compatible spelling of the above |
| `'min_overlap'` | 0.05 | Minimum region overlap ratio |
| `'min_region_points'` | 3 | Minimum in-region points for pass-through |
| `'max_tracts'` | `[]` | Display cap (2000 in whole mode, unlimited otherwise) |
| `'color_mode'` | `'direction'` | `'direction'`, `'fa'`, `'uniform'`, `'region'` |
| `'show_anatomy'` | true | FA background |
| `'show_region'` | true | Region overlay |
| `'line_width'` | 1.2 | Line width |
| `'region_alpha'` | 0.3 | Overlay transparency |
| `'grid_cols'` | auto | Grid columns |
| `'subplot_size'` | `[400, 350]` | Panel size |
| `'show_all_tracks'` | true | Force all tracks in grid mode |
| `'silent_export'` | true | Do not display figures while exporting |

### `visualizeTractographySlices(tracks_file, nim_file, x, y, z, varargin)`

Command-line viewer for three orthogonal slices. The slice positions `x` (sagittal),
`y` (coronal) and `z` (axial) are **required positional** arguments.

| Name | Default | Description |
|---|---|---|
| `'save'` | `''` | Output PNG filename |
| `'tolerance'` | 2 | Slice thickness in voxels for track selection |
| `'show_crosshairs'` | true | Draw slice intersections |
| `'show_anatomy'` | true | FA background |
| `'color_mode'` | `'direction'` | `'direction'` or `'uniform'` |
| `'alpha'` | 0.6 | FA background transparency |
| `'regions'` | `[]` | Region ID(s) to filter by |
| `'region_filter'` | `'pass_through'` | `'pass_through'`, `'start_in'`, `'end_in'` |
| `'min_overlap'` | 0.1 | Minimum region overlap |

### `visualizeTractographyAngles(varargin)`

Generates and exports 3D tractography figures from the standard anatomical
viewing angles (axial, sagittal, coronal, oblique). Takes the same
run-directory-or-file-pair first arguments as `visualizeTractography`, plus
`'output_dir'`.

### `nim_plot_bundles(items, opts)`

Renders one or more streamline sets in a shared voxel frame, and returns the
figure handle. `items` is a struct array with `.tracks`, `.name`, `.color`
(`[]` for per-segment direction coloring), `.alpha` (0.35), `.width` (1.2) and
`.clip` — an optional logical volume that keeps only the longest contiguous run of
each streamline inside it and reports the discarded fraction in
`items(i).frac_outside`. (Clipping used to be a separate helper; it is now this
field.)

`opts` takes `.views` (default sagittal / coronal / axial / oblique), `.title`,
`.save` (a path stem; writes `<stem>.png`), `.dpi` (150) and `.maxper` (900
streamlines drawn per item).

Per-segment direction coloring is what makes a whole-brain tractogram render as
anatomy rather than a monotone hairball.

### `nim_plot_vs_groundtruth(tracks_src, gt_trk, ref_nii, opts)`

Draws a run's streamlines (orange) against an ISMRM 2015 ground-truth bundle
(grey) in one shared voxel frame, mapping the ground truth through world RAS with
`nim_read_trk`. `tracks_src` is a run directory or a tracks `.mat`.

`opts`: `.max_gt` (700), `.max_ours` (700), `.title`, `.save` (writes
`<stem>.png` and `<stem>.fig`), `.roi` (JHU label(s) to outline as the seed
region) and `.nim` (needed only with `.roi`).

### Slice cache pipeline

| Function | Signature | Purpose |
|---|---|---|
| `generateSlices` | `(tracks_file, nim_file, output_dir, options)` | Server-side slice generation; returns `cache_dir` |
| `generateTractographySliceCache` | `(tracks_file, nim_file, cache_dir, options)` | The full cache build, with parallel rendering and validation |
| `TractographyCacheManager` | class, constructed from `output_dir` | Cache directory management with JSON metadata |
| `buildOptimizedTrackSliceLookup` | `(tracks, dims)` | Track↔slice intersection lookup, two-pass count-then-fill |
| `optimizedSliceRenderer` | `(ax, orientation, slice_value, slices, nim, tracks, slice_lookup, opts, dims)` | Fast single-slice render with precomputed color tables |
| `launchFastViewer` | `(tracks_file, nim_file, options)` | Generates the cache if needed and launches the Python viewer |

`generateSlices` options (a struct, not name-value pairs):

| Field | Default | Description |
|---|---|---|
| `.tolerance` | 2 | Slice thickness, voxels |
| `.color_mode` | `'direction'` | `'direction'` or `'uniform'` |
| `.show_anatomy` | true | FA background |
| `.alpha` | 0.6 | FA transparency |
| `.image_resolution` | `[1024, 1024]` | `[width, height]` in pixels |
| `.image_format` | `'png'` | `'png'` or `'jpg'` |
| `.compression_level` | 6 (PNG) / 90 (JPG) | Compression |
| `.parallel_workers` | auto | Worker count |
| `.regions` | `[]` (all) | Region IDs to include |

`output_dir` defaults to `./tractography_cache`.

`TractographyCacheManager` methods: `initializeCache(tracks_file, nim_file, options)`,
`saveMetadata()`, `loadMetadata()`, `validateCache()`, `getCacheStats()`.

`launchFastViewer` options: `.cache_dir`, `.force_regenerate` (false),
`.python_path` (`'python'`), `.viewer_path` (auto), `.auto_launch` (true),
`.background_mode` (false).

See [DISTRIBUTED_WORKFLOW.md](DISTRIBUTED_WORKFLOW.md) for the server → transfer →
local workflow these support.

---

## 12. General Plotting and Utilities

### `nim_plot(nim, varargin)`

**File**: `src/nim_plots/nim_plot.m`

The consolidated eigenvector plotter, replacing the former `nim_plotall`,
`nim_plotparcelall` and `nim_plotparcellation`. Name-value options:

| Name | Default | Description |
|---|---|---|
| `'mode'` | `'single'` | `'single'`, `'parcel'` (one parcel), `'parcels'` (all, one figure each) |
| `'indx'`, `'indy'`, `'indz'` | `[]` | Voxel index ranges, `'single'` mode |
| `'parcel_id'` | — | Parcel to draw, `'parcel'` mode |
| `'figindex'` | 1 | Figure number |
| `'downsample_factor'` | 1 | Downsample factor |
| `'show_figure'` | true | Create a new figure |
| `'show_progress'` | true | Pause and report progress |

A first argument that is a filename is loaded before plotting. Arrow color
follows the DTI convention: red left–right, green anterior–posterior, blue
superior–inferior.

### `src/nim_utils/` helpers

| Function | Signature | Description |
|---|---|---|
| `nim_reshape_d` | `(d)` | `[Dxx Dyy Dzz Dxy Dxz Dyz]` → the symmetric $3\times3$ matrix |
| `vector_to_color` | `(Vx, Vy, Vz)` | Unit-normalized direction vectors → $N \times 3$ RGB, taking the absolute value of each component |
| `nim_interp` | `(nim, p)` | Element-wise sub-voxel interpolation of the eigenvector field at $(p+1)^3$ uniform nodes per voxel (`p` default 3, nodes from `zwuni`); writes `nim_interp.mat` |
| `zwgll` | `(p)` | `[z, w]`: the $p+1$ Gauss-Lobatto-Legendre nodes on $[-1, 1]$ and their weights. Provided, but not currently called by `nim_interp` |
| `zwuni` | `(N)` | `[z, w]`: $N+1$ uniform nodes on $[-1, 1]$ and their weights |
| `nim_vis_eig` | `(nim, plane, slice, opts)` | Eigenvector figure for one slice; returns the figure handle |
| `nim_vis_fa` | `(nim, plane, slice, opts)` | FA figure for one slice; returns the figure handle |
| `plot_nim_interp` | `(grid, indx, indy, indz)` | Plots an interpolated grid produced by `nim_interp` |

`hdr.m`, `gen_vis_eig.m` and `runnimplot.m` are **scripts**, not functions:
plotting defaults, a batch of `nim_vis_eig` exports, and a `nim_plot` example
respectively.

---

## 13. Challenge Submissions (`src/nim_challenges/`)

### `nim_irontract_submit(nim_file, injection_file, output_dir, varargin)`

Generates IronTract Challenge submission volumes: binary NIfTI voxel masks from
tractography results, filtered by injection site, swept over a parameter range for
ROC analysis. Headers are taken from the injection mask so the volumes align.

| Option | Default | Description |
|---|---|---|
| `angle_thresholds` | `[30, 45, 60, 75, 90]` | Parameter sweep |
| `tracks_file` | `''` | Precomputed tracks; empty regenerates with injection-site seeding |
| `base_options` | — | Base tractography options |

**Output**: `submission_001.nii.gz`, `submission_002.nii.gz`, … in `output_dir`.

Reference: Maffei et al. (2022), *Insights from the IronTract challenge*,
NeuroImage.

---

## 14. Python Scripts (`scripts/`)

| Script | Purpose |
|---|---|
| `FastTractographyViewer.py` | Tkinter GUI for instant slice navigation over a pre-computed cache. `python3 scripts/FastTractographyViewer.py <cache_dir>`. Needs Python 3.7+, `tkinter`, `pillow`, `numpy` |
| `tractography_slice_gui.py` | Interactive slice GUI |
| `hinec_to_trk.py` | Convert HINEC `.mat` tracks to TrackVis `.trk` |
| `validate_ismrm_tractography.py` | Validate against the ISMRM 2015 ground truth |
| `compare_ismrm_results.py` | Compare scoring results across runs |
| `build_ismrm_scoring_config.py` | Build the scilpy scorer configuration |
| `run_pft_dipy.py` | DIPY particle-filtering tractography, for comparison |

`requirements.txt` covers the Python dependencies.

---

## Cross-References

- Architecture and data flow: [ARCHITECTURE.md](ARCHITECTURE.md)
- Mathematical foundations: [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md)
- Configuration parameter reference: [YAML_CONFIG.md](YAML_CONFIG.md)
- Run directory system: [RUN_DIRECTORY_SYSTEM.md](RUN_DIRECTORY_SYSTEM.md)
- Visualization guide: [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)
- User guide: [USER_GUIDE.md](USER_GUIDE.md)
