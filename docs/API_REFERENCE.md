# Complete Function Reference

Comprehensive reference for every public function in HINEC. Organized by module with function signatures, parameters, return values, and usage examples.

---

## 1. Entry Points

### `main(imgpath, nimpath, varargin)`

**File**: `main.m`

Core pipeline: data loading, preprocessing, DTI calculation, registration, and parcellation.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `imgpath` | char | Yes | Path to NIfTI image file (without extension) |
| `nimpath` | char | Yes | Path to save processed .mat file (must end in `.mat`) |
| `'t1_file'` | char | No | Path to T1 anatomical file for registration |
| `'config'` | struct | No | YAML-loaded configuration structure |
| `'run_info'` | struct | No | Run directory organization info |

**Returns**: None (saves nim structure to `nimpath`)

**Behavior**:

- Skips processing if output file already exists (caching)
- Auto-detects raw data (filename contains `_raw`) and triggers preprocessing
- Sequentially runs: nim_read → nim_preprocessing → nim_dt_spd → nim_eig → nim_fa → nim_registration → nim_parcellation → tissue segmentation
- Adds all nim_* subdirectories and SPM12 to MATLAB path

**Example**:
```matlab
main('data/sample', 'output/sample.mat');
main('data/sample_raw', 'output/sample.mat', 't1_file', 'data/t1.nii.gz');
```

---

### `runTractography(data_path, varargin)`

**File**: `runTractography.m`

Fiber tractography execution with algorithm selection and seeding strategy.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `data_path` | char | Yes | Path to .mat file with processed nim structure |
| `'algorithm'` | char | No | `'standard'` (FACT) or `'hinec'` (high-order). Default: `'standard'` |
| `'config'` | struct | No | YAML configuration structure |
| `'run_info'` | struct | No | Run directory organization |
| `'IronTract'` | flag | No | Enable IronTract Challenge submission mode |
| `'injection_file'` | char | No | Injection site mask for IronTract |
| `'irontract_output_dir'` | char | No | Output directory for IronTract submissions |

**Returns**: None (saves tracks to `tractography_results/` directory)

**Output file**: `tracks_{algorithm}_{YYYY-MM-DD_HH_MM_SS}.mat` containing:

- `tracks` — cell array of fiber tracks
- `options` — parameters used
- `elapsed_time` — computation time
- `algorithm` — algorithm name

**Example**:
```matlab
runTractography('output/sample.mat');
runTractography('output/sample.mat', 'algorithm', 'hinec');
```

---

### `runhinec()`

**File**: `runhinec.m`

Quick-start demo script. Loads `sample_parcellated.mat` and visualizes parcellations.

---

## 2. Data I/O (`src/nim_utils/`)

### `nim_read(imgpath)`

**File**: `src/nim_utils/nim_read.m`

Read NIfTI-1 diffusion MRI data with acquisition parameters.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `imgpath` | char | Yes | Path to NIfTI file (without extension) |

**Returns**: `nim` struct with fields:

- `.img` — 4D image data [X, Y, Z, N]
- `.img_b0` — B0 (non-diffusion) volumes
- `.img_bi` — Diffusion-weighted volumes
- `.bval` — B-values [1, N]
- `.bvec` — B-vectors [3, N] or [N, 3]
- `.hdr` — NIfTI header
- `.mask` — Brain mask (if `{name}_M.nii.gz` exists)
- `.xdim`, `.ydim`, `.zdim` — Volume dimensions
- `.size3` — Product of 3D dimensions
- `.size_bi` — Number of diffusion-weighted volumes
- `.thrsh_b0` — B-value threshold for b0 classification (50)

**Behavior**:
- Automatically separates b0 (bval ≤ 50) from diffusion-weighted volumes
- Supports both standard format (`{name}.bval/.bvec`) and IronTract format
- Loads brain mask from `{name}_M.nii.gz` if present

---

### `nim_save(nim, filepath)`

**File**: `src/nim_utils/nim_save.m`

Save processed nim structure to MAT file.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM data structure |
| `filepath` | char | Yes | Output .mat file path |

**Behavior**: Uses v7.3 format for files >1900MB, standard format otherwise.

---

### `nim_load_nim(file_prefix)`

**File**: `src/nim_utils/nim_load_nim.m`

Load a previously saved nim structure with associated masks and labels.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `file_prefix` | char | Yes | File path prefix for data files |

**Returns**: `nim` struct with image data, mask, parcellation mask, and atlas labels.

---

### `nim_load_labels(nim)`

**File**: `src/nim_utils/nim_load_labels.m`

Load parcellation labels into nim structure. Searches for XML (FSL atlas) and MAT file formats.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure with `parcellation_mask_file` field |

**Returns**: `nim` updated with `.atlas_labels` (containers.Map) field.

---

### `nim_load_atlas_labels(atlas_type)`

**File**: `src/nim_utils/nim_load_atlas_labels.m`

Load atlas labels from FSL XML file.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `atlas_type` | char | Yes | `'HarvardOxford'`, `'JHU'`, or `'JHU-tract'` |

**Returns**: `labels` struct with `.map` (containers.Map) and `.atlas_type`.

---

### `load_config_yaml(config_file)`

**File**: `src/nim_utils/load_config_yaml.m`

Parse YAML configuration file into MATLAB struct.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `config_file` | char | Yes | Path to YAML config file |

**Returns**: `config` struct with `.preprocessing` and `.tractography` sections. Validates tractography parameters and sets defaults for missing fields.

---

### `create_run_directory(config_file, varargin)`

**File**: `src/nim_utils/create_run_directory.m`

Create organized, timestamped run directory for reproducibility.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `config_file` | char | Yes | Path to YAML config file |
| `'base_dir'` | char | No | Base directory. Default: `'hinec_runs'` |
| `'run_name'` | char | No | Custom run name |
| `'description'` | char | No | Run description |

**Returns**: `run_info` struct with fields:

- `.run_dir` — Full path to created directory
- `.config_file` — Path to copied config
- `.logs_dir`, `.intermediate_dir`, `.output_dir`, `.tractography_dir`, `.diagnostics_dir`
- `.run_id` — Unique identifier
- `.timestamp` — ISO timestamp

**Behavior**: Creates subdirectory structure, copies config, records git info, creates `latest` symlink.

---

## 3. Calculation (`src/nim_calculation/`)

### `nim_dt_spd(nim, opts)`

**File**: `src/nim_calculation/nim_dt_spd.m`

Calculate diffusion tensors with SPD (symmetric positive-definite) constraint.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure with `.img`, `.bval`, `.bvec` |
| `opts.Mask` | string | No | `"on"` (default) or `"off"` — use brain mask |
| `opts.Steps` | int | No | BFGS max steps. Default: 2000 |
| `opts.FontSize` | int | No | Plot font size. Default: 16 |

**Returns**: `nim` updated with:

- `.DT` — [X, Y, Z, 6] diffusion tensor components
- `.evec` — [X, Y, Z, 3, 3] sorted eigenvectors
- `.eval` — [X, Y, Z, 3] sorted eigenvalues (descending)

**Behavior**:

1. Compute ADC: Y = ln(S₀/S) / b
2. Build design matrix H from gradient directions
3. Least squares fit: D = H \ Y for each voxel
4. Check eigenvalues; if any negative → BFGS optimization
5. Validate: reject NaN/Inf and unrealistic eigenvalues (> 0.01)
6. Sort eigenvalues descending via `maxk()`

---

### `nim_dt(nim, opts)`

**File**: `src/nim_calculation/nim_dt.m`

Basic diffusion tensor calculation without SPD enforcement.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure with `.img`, `.bval`, `.bvec` |
| `opts.Mask` | string | No | `"on"` (default) or `"off"` |

**Returns**: `nim` updated with `.DT` — [X, Y, Z, 6] diffusion tensor. Reports count of non-PSD tensors.

---

### `nim_eig(nim, opts)`

**File**: `src/nim_calculation/nim_eig.m`

Eigendecomposition of diffusion tensors.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure with `.DT` |
| `opts.Mask` | string | No | `"on"` (default) or `"off"` |

**Returns**: `nim` updated with:

- `.evec` — [X, Y, Z, 3, 3] eigenvectors (columns sorted by eigenvalue)
- `.eval` — [X, Y, Z, 3] eigenvalues (sorted descending via `maxk()`)

---

### `nim_fa(nim, opts)`

**File**: `src/nim_calculation/nim_fa.m`

Compute fractional anisotropy from eigenvalues.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure with `.eval` |
| `opts.Mask` | string | No | `"on"` (default) or `"off"` |

**Returns**: `nim` updated with `.FA` — [X, Y, Z] fractional anisotropy values in [0, 1].

---

## 4. Preprocessing (`src/nim_preprocessing/`)

### `nim_preprocessing(file_prefix, varargin)`

**File**: `src/nim_preprocessing/nim_preprocessing.m`

Orchestrator for the 10-step FSL-based preprocessing pipeline (815 lines).

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `file_prefix` | char | Yes | Path prefix for input files |
| `options` | struct | No | Preprocessing options from YAML config |
| `run_info` | struct | No | Run directory information |

**Returns**: `result` struct with `.preprocessed_file`, `.brain_mask_file`, `.parcellation_file`.

### Individual Preprocessing Functions

| Function | Signature | FSL Tool | Purpose |
|---|---|---|---|
| `preproc_denoising` | `(dwi_file, file_prefix, method)` | dwidenoise/SUSAN | Remove thermal noise |
| `preproc_extract_b0` | `(dwi_file, output_dir)` | fslroi | Extract b0 reference |
| `preproc_motion_correction` | `(dwi_file, bvec_file, bval_file, file_prefix)` | mcflirt | Correct head motion |
| `preproc_eddy_correction` | `(dwi_file, mask_file, bvec_file, bval_file, acqp_file, index_file, prefix)` | eddy | Correct eddy currents |
| `preproc_brain_extraction` | `(dwi_or_b0_file, output_dir, brain_mask_file)` | BET (-f 0.3) | Generate brain mask |
| `preproc_t1_brain_extraction` | `(t1_file, output_dir)` | BET (-f 0.4 -R) | Extract T1 brain |
| `preproc_fieldmap_correction` | `(dwi_file, fieldmap_file, output_prefix, varargin)` | FUGUE | Correct field inhomogeneity |
| `preproc_mask_improvement` | `(brain_mask_file, fa_data, file_prefix)` | BET + fslmaths | FA-based mask refinement |
| `preproc_tissue_segmentation` | `(fa_data, brain_mask_file, file_prefix)` | fslmaths -ero | WM/GM/CSF classification |
| `preproc_create_dwi_reference` | `(dwi_file, brain_mask_file, output_dir, file_prefix)` | fslroi + fslmaths | Create clean reference |
| `preproc_t1_dwi_registration` | `(dwi_ref, t1, t1_brain, output_dir, prefix)` | epi_reg | Register T1 to DWI |
| `preproc_t1_mni_registration` | `(t1_file, t1_brain, output_dir, prefix)` | FLIRT + FNIRT | Register T1 to MNI |
| `preproc_atlas_resampling` | `(ref_file, output_dir, prefix, atlas_type, options)` | applywarp/flirt | Atlas to DWI space |
| `preproc_cleanup` | `(output_dir, file_prefix, keep_files)` | — | Remove intermediate files |

### `preproc_fieldmap_correction` Optional Parameters (varargin)

| Name | Type | Default | Description |
|---|---|---|---|
| `'TE'` | double | auto-detect | Echo time difference |
| `'dwell_time'` | double | 0.00058 | Effective echo spacing (seconds) |
| `'phase_dir'` | char | `'y'` | Phase encoding direction |
| `'units'` | char | auto-detect | `'Hz'` or `'rad/s'` |
| `'smooth_sigma'` | double | 2.0 | Field map smoothing (mm) |
| `'mask_file'` | char | auto-generate | Brain mask path |

### `preproc_tissue_segmentation` Return Values

Returns three file paths: `[wm_mask_file, gm_mask_file, csf_mask_file]`

Thresholds: WM: FA > 0.2 (eroded), GM: 0.05 < FA ≤ 0.2, CSF: FA ≤ 0.05.

---

## 5. Registration (`src/nim_registration/`)

### `nim_registration(nim, options)`

**File**: `src/nim_registration/nim_registration.m`

Multi-modal registration pipeline orchestrator.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure with `.FA` |
| `options` | struct | Yes | Registration options with `.t1_file`, `.method` |

**Returns**: `registration_data` struct with `.transforms`, `.registered_images`, `.quality_metrics`.

---

### `register_dti_to_t1(registration_data, options)`

**File**: `src/nim_registration/register_dti_to_t1.m`

Register DTI to T1 space. Supports FSL (FLIRT, 6 DOF, correlation ratio) and SPM methods.

**Returns**: Updated `registration_data` with DTI↔T1 transforms.

---

### `register_t1_to_mni(registration_data, options)`

**File**: `src/nim_registration/register_t1_to_mni.m`

Register T1 to MNI152 template. Two-stage: linear (FLIRT 12 DOF) + optional nonlinear (FNIRT).

**Returns**: Updated `registration_data` with T1↔MNI transforms and inverse warps.

---

### `nim_apply_transforms(input_file, registration_data, transform_chain, varargin)`

**File**: `src/nim_registration/nim_apply_transforms.m`

Apply a chain of transforms between coordinate spaces.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `input_file` | char | Yes | Input image file |
| `registration_data` | struct | Yes | Registration data with transforms |
| `transform_chain` | cell | Yes | Transform sequence, e.g. `{'mni_to_t1', 't1_to_dti'}` |
| `options` | struct | No | `.output_file`, `.interpolation` (`'linear'`/`'nearest'`/`'spline'`), `.method` |

**Returns**: `output_file` — path to transformed image.

---

### `registration_utils.m` Functions

| Function | Signature | Purpose |
|---|---|---|
| `extract_reference_volumes` | `(registration_data, options)` | Extract b0 volume from DWI |
| `compute_registration_quality` | `(registration_data, options)` | Compute NMI quality scores |
| `save_registration_data` | `(registration_data, options)` | Save registration results |
| `generate_registration_report` | `(registration_data, options)` | Create HTML quality report |

---

## 6. Parcellation (`src/nim_parcellation/`)

### `nim_parcellation(nim, mask_file)`

**File**: `src/nim_parcellation/nim_parcellation.m`

Atlas-based brain parcellation.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM structure |
| `mask_file` | char | Yes | Path to parcellation mask NIfTI file |

**Returns**: `nim` updated with `.parcellation_mask`, `.parcellation_mask_file`, `.labels`.

**Supported Atlases**: JHU, HarvardOxford, AAL. Handles dimension mismatches with automatic fixes.

---

### `nim_parcellation_registered(nim, registration_data, parcellation_mask_file)`

**File**: `src/nim_parcellation/nim_parcellation_registered.m`

Enhanced parcellation using multi-modal registration chain (MNI → T1 → DWI).

**Returns**: `nim` updated with `.parcellation_mask`, `.parcellation_info` (quality metrics, atlas type, coverage).

---

## 7. Tractography (`src/nim_tractography/`)

### `nim_tractography_standard(data_path, varargin)`

**File**: `src/nim_tractography/nim_tractography_standard.m` (674 lines)

FACT (Fiber Assignment by Continuous Tracking) deterministic tractography.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `data_path` | char or struct | Yes | Path to .mat file or nim struct |
| `options` | struct | No | Tractography parameters (see below) |

**Options**:

| Field | Type | Default | Description |
|---|---|---|---|
| `seed_density` | int | 1 | Seeds per voxel dimension |
| `seed_strategy` | string | `"uniform"` | Seeding strategy |
| `step_size` | double | 0.5 | Step size in voxel units |
| `fa_threshold` | double | 0.1 | FA threshold for seeding |
| `termination_fa` | double | 0.05 | FA threshold for termination |
| `angle_thresh` | double | 60 | Maximum angle (degrees) |
| `max_steps` | int | 5000 | Maximum boundary crossings |
| `min_length` | int | 10 | Minimum track points |
| `seed_mask` | 3D array | [] | Custom seed mask |
| `use_fa_seed_filter` | bool | false | Filter seeds by FA |
| `enable_diagnostics` | bool | true | Enable timing diagnostics |

**Returns**: `tracks` — cell array where each element is an Nx3 matrix of [x, y, z] coordinates.

---

### `nim_tractography_hinec(data_path, varargin)`

**File**: `src/nim_tractography/nim_tractography_hinec.m` (1347 lines)

High-order deterministic tractography with interpolation, numerical integration, and ACT.

**Parameters**: Same as standard, plus:

| Field | Type | Default | Description |
|---|---|---|---|
| `interp_method` | string | `"trilinear"` | `"trilinear"`, `"cubic"`, or `"none"` |
| `integration_order` | int | 4 | 1=Euler, 2=RK2, 4=RK4, 5=RKF45 |
| `step_size` | double | 0.2 | Fixed step size (voxel units) |
| `adaptive_step` | bool | false | Enable RKF45 adaptive stepping |
| `rkf_tolerance` | double | 0.01 | RKF error tolerance (voxels) |
| `rkf_safety` | double | 0.9 | RKF safety factor |
| `step_min` | double | 0.01 | RKF minimum step size |
| `step_max` | double | 1.0 | RKF maximum step size |
| `act_enabled` | bool | true | Enable ACT |
| `wm_mask` | 3D array | [] | White matter mask |
| `gm_mask` | 3D array | [] | Gray matter mask |
| `csf_mask` | 3D array | [] | CSF mask |

**Returns**: `tracks` — cell array of Nx3 track matrices.

**Pre-computation**: Creates `griddedInterpolant` objects for FA and eigenvector components before tracking for 2-5× performance improvement.

---

### `nim_tractography_highorder(nim, options)`

**File**: `src/nim_tractography/nim_tractography_highorder.m` (439 lines)

Legacy high-order tractography with spline/linear interpolation and RK4/RK2 integration.

**Parameters**:

| Field | Type | Default | Description |
|---|---|---|---|
| `order` | int | 3 | Interpolation order (high-order: spline, 1: linear) |
| `step_size` | double | 0.2 | Integration step size |
| `fa_threshold` | double | 0.2 | FA threshold |
| `angle_thresh` | double | 45 | Angle threshold (degrees) |
| `seed_density` | int | 2 | Seeds per voxel dimension |

**Returns**: `tracks` — cell array of Nx3 track matrices.

---

## 8. Tractography Visualization (`src/nim_tractography/`)

### `nim_plot_tractography(tracks, nim, varargin)`

**File**: `src/nim_tractography/nim_plot_tractography.m`

3D visualization of fiber tracks.

**Options**:

| Field | Type | Default | Description |
|---|---|---|---|
| `color_mode` | char | `'direction'` | `'direction'`, `'fa'`, `'parcellation'` |
| `show_anatomy` | bool | true | Show FA slice background |
| `track_subset` | array | all | Indices of tracks to display |
| `line_width` | double | 1 | Line width |
| `alpha` | double | 0.8 | Transparency |

Limits display to 2000 tracks for performance.

---

### `nim_plot_tractography_region(tracks, nim, region_id, varargin)`

**File**: `src/nim_tractography/nim_plot_tractography_region.m`

Visualize tracks for a specific brain region with filtering.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `tracks` | cell | Yes | Fiber tracks |
| `nim` | struct | Yes | NIM with `.parcellation_mask` |
| `region_id` | int | Yes | Parcellation region index |

**Options**:

| Field | Type | Default | Description |
|---|---|---|---|
| `filter_mode` | char | `'pass_through'` | `'pass_through'`, `'start_in'`, `'end_in'`, `'connect_to'` |
| `min_overlap` | double | 0.3 | Minimum fraction in region |
| `show_region` | bool | true | Show region isosurface overlay |
| `region_alpha` | double | 0.3 | Overlay transparency |
| `track_color` | char | `'direction'` | `'direction'`, `'fa'`, `'uniform'`, `'region'` |
| `max_tracks` | double | inf | Maximum tracks to display |
| `min_track_length` | int | 10 | Minimum track points |
| `show_stats` | bool | true | Display statistics |

---

### `nim_plot_connectivity_matrix(tracks, nim, varargin)`

**File**: `src/nim_tractography/nim_plot_connectivity_matrix.m`

Compute and visualize region-to-region connectivity matrix.

**Options**:

| Field | Type | Default | Description |
|---|---|---|---|
| `min_track_length` | int | 10 | Minimum track points |
| `normalize` | bool | true | Normalize to [0, 1] |
| `symmetric` | bool | true | Symmetrize matrix |

**Returns**: `connectivity_matrix` — NxN connectivity matrix.

Creates 4-subplot figure: heatmap, histogram, node strengths, and summary.

---

### `nim_plot_vector_field(nim, varargin)`

**File**: `src/nim_tractography/nim_plot_vector_field.m`

Visualize eigenvector field on a 2D slice.

**Options**:

| Field | Type | Default | Description |
|---|---|---|---|
| `slice` | int | middle | Slice index |
| `downsample` | int | 2 | Downsample factor |
| `scale` | double | 1 | Vector scale factor |
| `axis_view` | char | `'axial'` | `'axial'`, `'coronal'`, `'sagittal'` |

Masks vectors by FA > 0.2 and overlays on FA background.

---

### `nim_excitation_time_map(nim, seed_points, varargin)`

**File**: `src/nim_tractography/nim_excitation_time_map.m`

Compute excitation propagation time map based on DTI anisotropy.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `nim` | struct | Yes | NIM with `.FA` and `.evec` |
| `seed_points` | Nx3 | Yes | Seed coordinates |

**Options**:

| Field | Type | Default | Description |
|---|---|---|---|
| `conduction_speed` | double | 1.0 | Base conduction speed |
| `fa_scaling` | bool | true | Scale velocity by FA |
| `max_time` | double | 100 | Maximum time cap |
| `method` | char | `'fast_marching'` | `'fast_marching'` or `'dijkstra'` |

**Returns**: `[time_map, velocity_field]` — 3D time map and 3D×3 velocity field.

---

## 9. Advanced Visualization (`src/nim_visualization/`)

### `visualizeTractography(tracks_file, nim_file, varargin)`

**File**: `src/nim_visualization/visualizeTractography.m` (1609 lines)

Unified tractography visualization with 4 modes. Supports run directory auto-detection.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `tracks_file` | char | Yes | Path to tracks .mat file (supports wildcards and run dirs) |
| `nim_file` | char | Yes | Path to nim .mat file |

**Name-Value Options**:

| Name | Type | Default | Description |
|---|---|---|---|
| `'mode'` | char | `'whole'` | `'whole'`, `'region'`, `'grid'`, `'sequential'` |
| `'regions'` | array | all | Region IDs to visualize |
| `'filter_mode'` | char | `'pass_through'` | `'pass_through'`, `'start_in'`, `'end_in'`, `'connect_to'` |
| `'color_mode'` | char | `'direction'` | `'direction'`, `'fa'`, `'uniform'`, `'region'` |
| `'max_tracks'` | int | 5000 | Maximum tracks to display |
| `'export'` | char | `''` | Export path (`.png`, `.pdf`, `.eps`, `.fig`) |
| `'show_stats'` | bool | true | Display track statistics |

---

### `visualizeTractographySlices(tracks_file, nim_file, varargin)`

**File**: `src/nim_visualization/visualizeTractographySlices.m` (536 lines)

Command-line 2D slice viewer with three orthogonal views.

**Name-Value Options**:

| Name | Type | Default | Description |
|---|---|---|---|
| `'slice_axial'` | int | middle | Axial slice index |
| `'slice_sagittal'` | int | middle | Sagittal slice index |
| `'slice_coronal'` | int | middle | Coronal slice index |
| `'tolerance'` | double | 1.0 | Slice thickness for track filtering |
| `'regions'` | array | all | Region filter |
| `'export_dir'` | char | `''` | Directory for batch PNG export |

---

### `generateSlices(tracks_file, nim_file, output_dir, varargin)`

**File**: `src/nim_visualization/generateSlices.m` (246 lines)

Server-side pre-generation of slice images for fast distributed viewing.

**Parameters**:

| Parameter | Type | Required | Description |
|---|---|---|---|
| `tracks_file` | char | Yes | Path to tracks .mat file |
| `nim_file` | char | Yes | Path to nim .mat file |
| `output_dir` | char | Yes | Output directory for cache |

**Name-Value Options**:

| Name | Type | Default | Description |
|---|---|---|---|
| `'resolution'` | [1x2] | [800, 600] | Image resolution [width, height] |
| `'format'` | char | `'png'` | `'png'` or `'jpg'` |
| `'quality'` | int | 85 | JPEG quality (1-100) |
| `'parallel'` | bool | auto | Use parallel workers |

---

### `generateTractographySliceCache(tracks, nim, cache_dir, varargin)`

**File**: `src/nim_visualization/generateTractographySliceCache.m` (569 lines)

Complete cache generation pipeline with parallel processing and quality validation.

---

### `TractographyCacheManager`

**File**: `src/nim_visualization/TractographyCacheManager.m` (230 lines)

Class for managing slice cache directories with JSON metadata.

**Methods**:

| Method | Description |
|---|---|
| `TractographyCacheManager(cache_dir)` | Constructor |
| `initialize(tracks_info, nim_info, params)` | Create cache structure |
| `validate()` | Check cache integrity |
| `getStatistics()` | Image counts and total size |

---

### `buildOptimizedTrackSliceLookup(tracks, dims)`

**File**: `src/nim_visualization/buildOptimizedTrackSliceLookup.m` (146 lines)

Pre-compute vectorized track-slice intersection lookup tables. Two-pass algorithm: count → pre-allocate → fill. Returns lookup struct with per-slice track indices.

---

### `optimizedSliceRenderer(tracks, nim, orientation, slice_idx, varargin)`

**File**: `src/nim_visualization/optimizedSliceRenderer.m` (331 lines)

High-performance slice rendering with pre-computed color lookup tables and FA background.

---

### `launchFastViewer(tracks_file, nim_file, varargin)`

**File**: `src/nim_visualization/launchFastViewer.m` (257 lines)

MATLAB bridge to Python FastTractographyViewer. Auto-generates cache, checks Python dependencies, launches GUI.

---

## 10. General Plotting (`src/nim_plots/`)

### `nim_plot(nim, varargin)`

**File**: `src/nim_plots/nim_plot.m` (348 lines)

DTI eigenvector visualization.

**Name-Value Options**:

| Name | Type | Default | Description |
|---|---|---|---|
| `'mode'` | char | `'single'` | `'single'`, `'parcel'`, `'parcels'` |
| `'region_id'` | int | 1 | Region to display (parcel mode) |
| `'downsample'` | int | 1 | Downsample factor |

Color coding follows DTI convention: Red = Left/Right (X), Green = Anterior/Posterior (Y), Blue = Superior/Inferior (Z).

---

## 11. Utilities (`src/nim_utils/`)

### `nim_reshape_d(d)`

Reshape 6-element diffusion vector to 3×3 symmetric matrix.

**Input**: `d` — [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]

**Output**: `D` — 3×3 symmetric tensor matrix.

---

### `vector_to_color(Vx, Vy, Vz)`

Convert direction vectors to RGB colors using DTI convention.

**Returns**: `colors` — Nx3 RGB matrix where R=|Vx|, G=|Vy|, B=|Vz|.

---

### `nim_interp(nim, p)`

High-order spectral interpolation of eigenvectors using GLL quadrature nodes.

**Parameters**: `nim` — NIM structure, `p` — interpolation order (default: 3).

**Output**: Saves interpolated grids to `nim_interp.mat`.

---

### `zwgll(p)`

Compute Gauss-Lobatto-Legendre quadrature nodes and weights.

**Input**: `p` — polynomial degree.

**Returns**: `[z, w]` — (p+1)×1 nodes on [-1, 1] and weights.

---

### `zwuni(N)`

Compute uniform quadrature nodes and weights.

**Input**: `N` — number of intervals.
**Returns**: `[z, w]` — (N+1)×1 uniform nodes on [-1, 1] and weights.

---

### `hdr(nim)`

NIfTI header utility functions.

---

### `nim_vis_eig(nim)`, `nim_vis_fa(nim)`, `gen_vis_eig(nim)`

Eigenvector and FA visualization helper functions.

---

### `plot_nim_interp(nim)`, `runnimplot(nim)`

Interpolation result plotting and plot workflow utilities.

---

## 12. Python Scripts (`scripts/`)

### `FastTractographyViewer.py`

Python GUI for instant slice navigation using pre-computed images. Requires: Python 3.7+, tkinter, Pillow, numpy.

**Launch**: `python FastTractographyViewer.py <cache_directory>`

### `hinec_to_trk.py`

Convert HINEC .mat track output to TrackVis .trk format.

### `validate_ismrm_tractography.py`

Validate tractography results against ISMRM 2015 ground truth.

---

## 13. Challenge Submissions (`src/nim_challenges/`)

### `nim_irontract_submit(tracks, nim, output_dir, options)`

**File**: `src/nim_challenges/nim_irontract_submit.m`

Format and export tractography results for IronTract Challenge submission.

---

## Cross-References

- Architecture and data flow: [ARCHITECTURE.md](ARCHITECTURE.md)
- Mathematical foundations: [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md)
- Configuration system: [YAML_CONFIG.md](YAML_CONFIG.md)
- Visualization guide: [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)
- User guide: [USER_GUIDE.md](USER_GUIDE.md)
