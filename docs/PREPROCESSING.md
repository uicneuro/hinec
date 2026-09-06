# HINEC Preprocessing Pipeline

This document covers the HINEC preprocessing pipeline: its workflow steps, the fallback path
taken when a step's preferred tool is unavailable, and the error conditions each step can
raise.

## Overview

`nim_preprocessing.m` prepares raw diffusion-weighted MRI for tensor estimation and
tractography. It runs ten sequential steps, from raw NIfTI through to an atlas resampled into
DWI space, using FSL throughout. `main.m` decides whether it needs to run at all, based on
which files are present.

Where a T1 anatomical image is available, the pipeline uses it for brain extraction and for
atlas registration through an MNI→T1→DWI transformation chain rather than a direct MNI→DWI
registration. The pipeline itself produces no white matter mask; tissue segmentation for ACT
happens later, in `main.m` (see [Tissue Segmentation for ACT](#tissue-segmentation-for-act)).

## Architecture

### Preprocessing Detection Logic (main.m)

The main pipeline uses a hierarchical file detection system:

```matlab
% Detection order:
1. Check for processed file: {name}.nii.gz
2. If not found, check for raw file: {name}_raw.nii.gz
3. If raw file exists, trigger preprocessing
4. If neither exists, throw error
```

**File naming conventions:**

- Raw data: `{name}_raw.nii.gz`, `{name}.bval`, `{name}.bvec`
- T1 structural data: `{name}_T1.nii.gz` (optional for enhanced processing)
- Field map data: `{name}_fmap_Hz.nii.gz` (optional for distortion correction)
- Processed data: `{name}.nii.gz` (output of preprocessing)
- Brain mask: `{name}_M.nii.gz`
- Acquisition parameters: `{name}_acqp.txt`, `{name}_index.txt` (optional for advanced eddy correction)

### Pipeline Entry Points

1. **Direct preprocessing call:** `nim_preprocessing(file_prefix, options)`
2. **Automatic from main:** Called when raw data detected in `main.m`

## Preprocessing Workflow

### Phase 1: Environment Setup and Validation

#### Step 1: B0 Volume Extraction
**Function:** `preproc_extract_b0.m`

**Purpose:** Extract b=0 volume for brain extraction and registration

```bash
fslroi {input_dwi} {b0_output} 0 1
```

**Fallback:** Use first volume if b=0 not identifiable

**Error condition:** DWI file unreadable or empty

#### Step 2: Brain Extraction

**Function:** `preproc_t1_brain_extraction.m` (when T1 available) or `preproc_brain_extraction.m`

**Purpose:** Create the brain mask, using T1 structural data when available

**T1-Based Processing (Preferred):**

1. T1 brain extraction using FSL BET with T1-optimized parameters (`-f 0.4`)
2. Boundary-based registration of T1 brain mask to DWI space using FSL epi_reg
3. Transfer the T1-derived brain mask into DWI space

```bash
# T1 brain extraction
bet {T1_file} {T1_brain} -R -m -f 0.4

# T1-DWI registration
epi_reg --epi={b0_volume} --t1={T1_file} --t1brain={T1_brain} --out={T1_to_dwi}

# Transform T1 mask to DWI space
flirt -in {T1_brain_mask} -ref {b0_volume} -applyxfm -init {T1_to_dwi.mat} -interp nearestneighbour -out {dwi_brain_mask}
```

**DWI-Only Fallback:**
```bash
bet {b0_volume} {brain_extracted} -m -f 0.3
```

**T1 Detection:** Automatic detection of `{name}_T1.nii.gz` in the same directory

**Rationale:** T1 gives sharper tissue contrast at the brain boundary than a b0 volume, so the
mask follows the actual boundary more closely.

**Fallback:** DWI-based BET if T1 unavailable or T1 processing fails

**Error condition:** All brain extraction methods fail

### Phase 2: Distortion and Artifact Correction

#### Step 3: Denoising (Optional)

**Function:** `preproc_denoising.m`

**Methods:**

- `dwidenoise` (default) - MRtrix3 denoising
- `nlmeans` - Non-local means
- `gaussian` - Gaussian filtering

**Configuration:**
```matlab
options.run_denoising = true;  % Enable/disable
options.denoise_method = 'dwidenoise';  % Method selection
```

**Fallback:** Skip if denoising tools unavailable

**Error condition:** All denoising methods fail

#### Step 4: Field Map Distortion Correction (Optional)

**Function:** `preproc_fieldmap_correction.m`

**Purpose:** Correct susceptibility-induced distortions

**Configuration:**
```matlab
options.run_fieldmap_correction = true;
options.fieldmap_file = 'auto';  % Auto-detect or specify path
options.phase_encoding_dir = 'y';  % Phase encoding direction
options.dwell_time = 0.00058;  % Effective echo spacing (seconds)
```

**Process:**

1. Auto-detect field map units (Hz vs rad/s)
2. Convert rad/s to Hz if necessary
3. Apply brain mask and smooth field map
4. Apply FUGUE distortion correction

**Fallback:** Skip field map correction if field map missing/invalid

**Error conditions:**

- Field map file corrupt or wrong format
- FUGUE command fails
- Incompatible field map dimensions

#### Step 5: Motion Correction (Optional)

**Function:** `preproc_motion_correction.m`

**Purpose:** Correct head motion between volumes

**Process:**

1. Extract b=0 volumes as motion reference
2. Run FSL mcflirt for volume-to-volume alignment
3. Rotate b-vectors according to motion parameters
4. Generate motion quality metrics

```bash
mcflirt -in {dwi} -out {corrected} -refvol {b0_index} -plots -mats -rmsrel -rmsabs
```

**Configuration:**
```matlab
options.run_motion_correction = true;
```

**Fallback:** Copy original data if motion correction fails

**Error conditions:**

- No b=0 volumes found
- mcflirt fails
- Motion parameters corrupt

**Quality thresholds:**

- Translation warning: >3mm
- Rotation warning: >3 degrees
- RMS displacement warning: >1mm

#### Step 6: Eddy Current Correction (Optional)

**Function:** `preproc_eddy_correction.m`

**Purpose:** Correct eddy current distortions and residual motion

**Automatic Method Selection:**
```matlab
if exist(acqp_file) && exist(index_file)
    method = 'eddy';  % Full FSL eddy
else
    method = 'eddy_correct';  % Basic correction
end
```

**Full `eddy` (Preferred):**
```bash
eddy --imain={dwi} --mask={mask} --bvecs={bvec} --bvals={bval}
     --out={output} --acqp={acqp} --index={index}
     --repol --cnr_maps --residuals
```

**Basic Fallback:**
```bash
eddy_correct {dwi} {output} 0
```

**Parameter File Generation:**
If acquisition parameter files missing, pipeline attempts to create them from options:
```matlab
options.phase_encoding_direction = 'y-';  % BIDS format
options.total_readout_time = 0.05;  % seconds
options.eddy_index_vector = [];  % Auto-generate if empty
```

**Fallback hierarchy:**

1. Use full `eddy` with existing parameter files
2. Generate parameter files from options, then use full `eddy`
3. Use basic eddy_correct
4. Skip eddy correction (with warning)

**Error conditions:**

- Both eddy methods fail
- Parameter files corrupt
- Insufficient system resources for eddy

### Phase 3: Finalization

#### Step 7: Data Finalization

**Purpose:** Copy processed data to final locations with proper file naming

**Process:**

1. Copy brain mask to final location as `{name}_M.nii.gz`
2. Copy processed DWI data to final output file
3. Copy final b-vectors (motion/eddy corrected)
4. Update preprocessing report with final file locations

**No Configuration Required** - Always executed

#### Step 8: T1 Preprocessing and Registration (When Available)

**Functions:** `preproc_t1_dwi_registration.m`, `preproc_t1_mni_registration.m`, `preproc_create_dwi_reference.m`

**Purpose:** Complete T1-based registration workflow for enhanced atlas processing

**T1 Registration Chain (When T1 Available):**

1. **T1-DWI Registration:** Refine boundary-based registration with final processed DWI
2. **T1-MNI Registration:** Create nonlinear T1-to-MNI transformation using FSL FNIRT
3. **DWI Reference Creation:** Generate distortion-corrected DWI reference volume

**Process Flow:**
```bash
# T1-DWI registration refinement
epi_reg --epi={dwi_ref} --t1={T1_file} --t1brain={T1_brain} --out={T1_to_dwi_refined}

# T1-MNI nonlinear registration
flirt -in {T1_brain} -ref {MNI_template} -omat {T1_to_MNI_linear.mat} -out {T1_to_MNI_linear}
fnirt --in={T1_file} --aff={T1_to_MNI_linear.mat} --cout={T1_to_MNI_warp} --config=T1_2_MNI152_2mm

# Inverse warp for atlas transformation
invwarp -w {T1_to_MNI_warp} -r {T1_file} -o {MNI_to_T1_warp}
```

**Rationale:** T1-to-MNI registration is a same-contrast, undistorted problem and converges
far better than a direct MNI-to-DWI registration, which must cope with EPI distortion and
diffusion contrast at once. Routing the atlas through T1 inherits that better alignment.

**Configuration.** Both keys are derived, not chosen: when `main.m` calls preprocessing it
sets `t1_available` and `use_t1_registration` from whether `{name}_T1.nii.gz` actually exists,
overriding whatever the config carried. To take the direct path, call `nim_preprocessing`
yourself with the keys left unset (absent means off).

```matlab
options.use_t1_registration = true;
options.t1_available = true;
options.t1_file = '/path/to/subject_T1.nii.gz';
```

**Fallback:** Skip T1 registration if T1 data unavailable or T1 processing fails

**Error condition:** T1 registration chain fails (falls back to direct atlas registration)

#### Step 9: Atlas Processing

**Function:** `preproc_atlas_resampling.m`

**Purpose:** Resample atlas to DWI space for parcellation using optimal registration method

**T1-Based Atlas Registration (When T1 Available):**
Uses the complete MNI→T1→DWI transformation chain created in Step 8:
```bash
# Apply composite transformation
applywarp -i {atlas_MNI} -r {dwi_ref} --warp={MNI_to_T1_warp} --postmat={T1_to_dwi.mat} --interp=nn -o {atlas_dwi}
```

**Direct Registration Fallback (DWI-Only):**
```bash
# Direct atlas resampling (fallback method)
flirt -in {atlas_MNI} -ref {dwi_ref} -out {atlas_dwi} -interp nearestneighbour
```

**Supported atlases:**

- `HarvardOxford` (HarvardOxford-cort-maxprob-thr0-1mm) — the default in `nim_preprocessing`,
  and the fallback for an unrecognised name
- `JHU` (JHU-ICBM-labels-1mm)
- `JHU-tract` (JHU-ICBM-tracts-maxprob-thr0-1mm)

All three are read from `$FSLDIR/data/atlases/`.

**Configuration:**
```matlab
options.atlas_type = 'HarvardOxford';
options.use_t1_registration = true;  % Enables T1-guided atlas registration
```

Note the two layers disagree on the default: `nim_preprocessing`'s own default is
`HarvardOxford`, while the YAML schema's `preprocessing.atlas_type` defaults to `jhu`. A run
launched from a config therefore gets JHU unless it says otherwise.

**Atlas Quality Validation:**

- Label value range checking (ensures integer labels preserved)
- Voxel coverage assessment (validates successful registration)
- Spatial consistency verification

**Why the T1 route:** it avoids registering an MNI template directly onto distorted,
diffusion-contrast data, and it preserves label integrity by resampling with nearest-neighbour
interpolation throughout.

**Fallback:** Use direct FLIRT registration if T1-based method unavailable

**Error condition:** All atlas registration methods fail

#### Step 10: Finalization and Cleanup

**Function:** `preproc_cleanup.m`

**Purpose:** Organize final outputs and remove temporary files

**Final outputs:**

- `{name}.nii.gz` - Processed DWI data
- `{name}.bvec` - Final b-vectors (motion/eddy corrected)
- `{name}.bval` - B-values (unchanged)
- `{name}_M.nii.gz` - Brain mask
- `parcellation_mask.nii.gz` - Atlas mask

## Configuration Options

### Default Configuration

These are `nim_preprocessing.m`'s own defaults, applied to any key the caller omits:

```matlab
default_options = struct(...
    'run_denoising', true, ...
    'denoise_method', 'dwidenoise', ...
    'run_motion_correction', true, ...
    'run_fieldmap_correction', true, ...
    'fieldmap_file', '', ...
    'fieldmap_units', 'auto', ...
    'phase_encoding_dir', 'y', ...
    'dwell_time', 0.00058, ...
    'eddy_method', 'auto', ...
    'run_eddy', true, ...
    'improve_mask', true, ...
    'atlas_type', 'HarvardOxford', ...
    'phase_encoding_direction', "", ...
    'total_readout_time', [], ...
    'eddy_index_vector', [] ...
);
```

The T1 keys (`t1_available`, `use_t1_registration`, `t1_file`) are deliberately absent from
this struct: an unset key means the T1 path is off, and `main.m` sets all three from whether
`{name}_T1.nii.gz` exists.

The YAML surface is narrower than this struct and is documented separately. The `preprocessing`
section of a config exposes `run_denoising`, `denoise_method`, `run_motion_correction`,
`run_eddy`, `improve_mask`, `atlas_type`, `t1_available`, `use_t1_registration` and
`register_to_mni`; `src/nim_utils/nim_config_schema.m` is the single source of truth for it,
and [YAML_CONFIG.md](YAML_CONFIG.md) is generated from that file.

### T1 Integration

The pipeline detects T1 structural data by the naming convention `{name}_T1.nii.gz` in the
same directory as the DWI, and turns on the T1 path when it is found. A T1 changes three
things: brain extraction runs on structural contrast rather than a b0 volume; the atlas
reaches DWI space through MNI→T1→DWI instead of directly; and tissue segmentation for ACT can
use FSL FAST instead of the FA-tertile fallback.

### Legacy Support

Two positional forms are still accepted for backward compatibility:

```matlab
nim_preprocessing(file_prefix, options);                   % current
nim_preprocessing(file_prefix, run_eddy);                  % legacy
nim_preprocessing(file_prefix, run_eddy, atlas_type);      % legacy
```

Both legacy forms force `run_denoising`, `run_motion_correction` and `improve_mask` on, and
`denoise_method` to `dwidenoise`.

## Error Handling and Recovery

### Comprehensive Error Reporting
All errors and warnings are logged in a preprocessing report:
```matlab
preprocessing_report = struct();
preprocessing_report.errors = {};      % Critical errors
preprocessing_report.warnings = {};    % Non-fatal issues
preprocessing_report.steps_completed = {};  % Successful steps
```

### Common Error Scenarios

#### 1. FSL Environment Issues

**Symptoms:** FSLDIR not set, FSL commands not found

**Recovery:** Install FSL, set environment variables

**Pipeline behavior:** Fatal error, cannot proceed

#### 2. Insufficient Disk Space

**Symptoms:** Preprocessing steps fail with I/O errors

**Recovery:** Free disk space, resume from last successful step

**Pipeline behavior:** Stop and save partial results

#### 3. Memory Limitations

**Symptoms:** Large dataset processing fails

**Recovery:** Reduce processing options, use cluster computing

**Pipeline behavior:** Fall back to less memory-intensive methods

#### 4. Corrupted Input Data

**Symptoms:** NIfTI reading errors, dimension mismatches

**Recovery:** Re-export data from scanner, verify file integrity

**Pipeline behavior:** Fatal error with detailed diagnostics

#### 5. Parameter File Issues

**Symptoms:** Eddy correction fails due to missing acqp/index files

**Recovery:** Automatic parameter generation from options

**Pipeline behavior:**

1. Try to generate parameters from options
2. Fall back to basic eddy_correct
3. Skip eddy correction if all methods fail

### Quality Control Metrics

#### Motion Assessment

- Maximum translation (mm)
- Maximum rotation (degrees)
- Mean relative RMS displacement (mm)
- Frame-to-frame displacement spikes

#### Eddy Current Assessment

- Outlier slices detected and corrected
- CNR (Contrast-to-Noise Ratio) maps
- Residual maps for quality assessment

#### Field Map Assessment

- Field map coverage (% of brain)
- Distortion magnitude statistics
- Correction effectiveness metrics

## Performance Considerations

### Processing Time Estimates

- Basic preprocessing (no eddy): 3-10 minutes
- With advanced eddy correction: 2-8 hours
- With field map correction: Additional 30-60 minutes

### System Requirements

- RAM: Minimum 8GB, recommended 16GB+
- Disk space: 3-5x input data size for temporary files
- CPU: Multi-core beneficial for FSL operations

### Optimization Tips

1. Use cluster computing for large datasets
2. Enable parallel processing in FSL
3. Use SSDs for temporary file storage
4. Pre-compute parameter files for batch processing

## Integration with Main Pipeline

The preprocessing module integrates seamlessly with the main HINEC pipeline:

1. **Automatic Detection:** Main pipeline detects raw vs processed data
2. **Transparent Processing:** Preprocessing runs automatically when needed
3. **Quality Validation:** Results validated before proceeding to DTI calculation
4. **Error Propagation:** Preprocessing errors properly handled by main pipeline

## Troubleshooting Guide

### Common Issues and Solutions

#### "FSLDIR not set" error
```bash
export FSLDIR=/usr/local/fsl
export PATH=${FSLDIR}/bin:${PATH}
source ${FSLDIR}/etc/fslconf/fsl.sh
```

#### Field map correction fails

1. Check field map units (Hz vs rad/s)
2. Verify phase encoding direction
3. Ensure field map covers brain region
4. Check dwell time parameter accuracy

#### Eddy correction takes too long

1. Use basic eddy_correct instead of full eddy
2. Reduce number of iterations
3. Use cluster computing with multiple cores
4. Consider GPU-accelerated eddy

#### Motion correction produces artifacts

1. Check motion parameters for excessive movement
2. Verify b-vector rotation is applied correctly
3. Consider excluding high-motion volumes
4. Use more robust reference volume selection

## Tissue Segmentation for ACT

Tissue masks are produced by `main.m`, not by `nim_preprocessing`, and are stored on the `nim`
as `.wm_mask`, `.gm_mask` and `.csf_mask` (with the corresponding NIfTI paths in
`.wm_mask_file` and friends).

**Primary method** (`preproc_tissue_segmentation.m`): FSL FAST on the anatomical T1, resampled
into DWI space through the two images' **world affines** (`flirt -usesqform`). Deliberately no
registration step is involved, so the masks cannot be corrupted by a T1→DWI registration that
is disabled or has converged poorly. This is the right approach whenever the T1 is already
world-aligned to the diffusion data.

**Fallback** (no usable T1): FA-tertile binning. This bins **anisotropy, not tissue**, so ACT
driven by it terminates streamlines mid-crossing. It exists for DWI-only datasets and should
be read as a degraded mode, not an equivalent one.

Output files are `{name}_WM_mask.nii.gz`, `{name}_GM_mask.nii.gz` and
`{name}_CSF_mask.nii.gz`. Note that these masks are generated on every run, but ACT itself is
off unless `tractography.act` is set true.

## Field Map Correction Reference

### Setup Example

```matlab
options = struct();
options.run_fieldmap_correction = true;
options.fieldmap_file = 'subject_fmap_Hz.nii.gz';
options.phase_encoding_dir = 'y';      % acquisition-dependent
options.dwell_time = 0.00058;          % scanner-dependent, seconds
options.eddy_method = 'eddy_correct';  % force the basic method

nim_preprocessing('path/to/subject', options);
main('path/to/subject', 'subject.mat');
runTractography('subject.mat');
```

### Typical Dwell Times by Vendor

| Vendor | Effective echo spacing |
|---|---|
| Siemens | ~0.00058 s |
| GE | ~0.000476 s |
| Philips | ~0.000694 s |

These are representative values only; read the actual figure from your sequence protocol.

## T1-Based Registration Implementation Details

The atlas reaches DWI space through an MNI→T1→DWI chain rather than a direct MNI→DWI
registration, because the direct route must absorb EPI distortion and diffusion contrast at
once. The chain is assembled from four steps:

1. **T1 Brain Extraction** (`preproc_t1_brain_extraction.m`):

   - Use FSL BET with robust skull stripping parameters
   - Generate T1 brain mask and brain-extracted T1 image

2. **T1-DWI Registration** (`preproc_t1_dwi_registration.m`):

   - Use FSL epi_reg for boundary-based registration
   - Handle EPI distortions and generate transformation matrix

3. **T1-MNI Registration** (`preproc_t1_mni_registration.m`):

   - Perform linear pre-alignment using FLIRT
   - Execute nonlinear registration using FNIRT
   - Generate and invert warp fields for MNI→T1 transformation

4. **Composite Transform Application**:

   - Use applywarp to combine MNI→T1 warp with T1→DWI affine
   - Apply nearest neighbor interpolation to preserve label values

### Coordinate System Considerations

- **MATLAB vs FSL**: Handle coordinate system differences
- **Voxel vs World**: All internal processing in voxel coordinates
- **Label Preservation**: Use nearest neighbor interpolation for atlases

## Known Limitations

- Motion correction is volume-to-volume; there is no slice-to-volume correction, so
  within-volume motion is not recovered.
- The FA-tertile tissue fallback is not a tissue segmentation and should not be used to drive
  ACT on data where a T1 is available.
- Eddy correction is not GPU-accelerated here; on a large dataset the full `eddy` path
  dominates the wall-clock time.
- There is no automatic BIDS ingestion; inputs follow the `{name}_*` prefix convention
  described above.

---

For behaviour not covered here, the function files under `src/nim_preprocessing/` are the
authority.