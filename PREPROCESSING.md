# HINEC Preprocessing Pipeline

This document provides a comprehensive overview of the HINEC preprocessing pipeline, including the workflow steps, fallback mechanisms, and error conditions.

## Overview

The HINEC preprocessing pipeline is a modular system that prepares raw diffusion-weighted MRI (dMRI) data for tensor calculation and tractography analysis. The pipeline automatically detects whether data requires preprocessing and applies appropriate corrections using FSL tools.

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
- Processed data: `{name}.nii.gz` (output of preprocessing)
- Brain mask: `{name}_M.nii.gz`
- Acquisition parameters: `{name}_acqp.txt`, `{name}_index.txt` (optional for advanced eddy correction)

### Pipeline Entry Points

1. **Direct preprocessing call:** `nim_preprocessing(file_prefix, options)`
2. **Automatic from main:** Called when raw data detected in `main.m`

## Preprocessing Workflow

### Phase 1: Environment Setup and Validation

#### Step 1: FSL Environment Initialization
```matlab
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set...');
end
system('source $FSLDIR/etc/fslconf/fsl.sh');
```

**Fallback:** None - FSL is mandatory for preprocessing
**Error condition:** Missing FSLDIR environment variable

#### Step 2: Input File Verification
Required files:
- `{name}_raw.nii.gz` - Raw DWI data
- `{name}.bvec` - Gradient directions (3×N or N×3 format)
- `{name}.bval` - B-values

**Fallback:** None - all files mandatory
**Error condition:** Any required file missing

#### Step 3: Field Map Auto-Detection
```matlab
% Search patterns:
- {name}_fmap_Hz.nii.gz
- {name}_fmap_RadPerSec.nii.gz
```

**Fallback:** Field map correction disabled if no field map found
**Warning condition:** Field map files not found

### Phase 2: Initial Processing

#### Step 4: B0 Volume Extraction
**Function:** `preproc_extract_b0.m`
**Purpose:** Extract b=0 volume for brain extraction and registration

```bash
fslroi {input_dwi} {b0_output} 0 1
```

**Fallback:** Use first volume if b=0 not identifiable
**Error condition:** DWI file unreadable or empty

#### Step 5: Initial Brain Extraction
**Function:** `preproc_brain_extraction.m`
**Purpose:** Create initial brain mask using FSL BET

```bash
bet {b0_volume} {brain_extracted} -m -f 0.3
```

**Fallback:** Adjust BET threshold (-f parameter) if initial extraction fails
**Error condition:** BET fails completely

### Phase 3: Distortion and Artifact Correction

#### Step 6: Denoising (Optional)
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

#### Step 7: Field Map Distortion Correction (Optional)
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

#### Step 8: Motion Correction (Optional)
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

#### Step 9: Enhanced Eddy Current Correction (Optional)
**Function:** `preproc_eddy_correction.m`
**Purpose:** Correct eddy current distortions and residual motion

**Automatic Method Selection:**
```matlab
if exist(acqp_file) && exist(index_file)
    method = 'eddy';  % Advanced FSL eddy
else
    method = 'eddy_correct';  % Basic correction
end
```

**Advanced Eddy (Preferred):**
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
1. Use advanced eddy with existing parameter files
2. Generate parameter files from options → use advanced eddy
3. Use basic eddy_correct
4. Skip eddy correction (with warning)

**Error conditions:**
- Both eddy methods fail
- Parameter files corrupt
- Insufficient system resources for eddy

### Phase 4: Segmentation and Finalization

#### Step 10: White Matter Segmentation (Optional)
**Purpose:** Create white matter mask for enhanced tractography

**Process:**
1. Run preliminary DTI calculation with dtifit
2. Create FA map
3. Threshold FA > 0.2 for white matter
4. Apply morphological erosion to avoid boundaries

```bash
dtifit -k {dwi} -o {temp_dti} -m {mask} -r {bvec} -b {bval}
fslmaths {temp_dti_FA} -thr 0.2 -bin {wm_mask}
fslmaths {wm_mask} -ero -ero {wm_mask}
```

**Configuration:**
```matlab
options.create_wm_mask = true;
```

**Fallback:** Use brain mask if WM segmentation fails
**Error condition:** DTI calculation fails

#### Step 11: Atlas Processing
**Function:** `preproc_atlas_resampling.m`
**Purpose:** Resample atlas to DWI space for parcellation

**Supported atlases:**
- HarvardOxford (default)
- AAL
- JHU-tract
- Destrieux

**Configuration:**
```matlab
options.atlas_type = 'HarvardOxford';
```

**Fallback:** Use HarvardOxford if specified atlas unavailable
**Error condition:** No atlas files found

#### Step 12: Finalization and Cleanup
**Function:** `preproc_cleanup.m`
**Purpose:** Organize final outputs and remove temporary files

**Final outputs:**
- `{name}.nii.gz` - Processed DWI data
- `{name}.bvec` - Final b-vectors (motion/eddy corrected)
- `{name}.bval` - B-values (unchanged)
- `{name}_M.nii.gz` - Brain mask
- `parcellation_mask.nii.gz` - Atlas mask
- `{name}_WM_mask.nii.gz` - White matter mask (if created)

## Configuration Options

### Default Configuration
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
    'create_wm_mask', true, ...
    'improve_mask', true, ...
    'atlas_type', 'HarvardOxford', ...
    'phase_encoding_direction', "", ...
    'total_readout_time', [], ...
    'eddy_index_vector', [] ...
);
```

### Legacy Support
The pipeline maintains backward compatibility:
```matlab
% New usage
nim_preprocessing(file_prefix, options);

% Legacy usage
nim_preprocessing(file_prefix, run_eddy, atlas_type);
```

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
- Basic preprocessing (no eddy): 5-15 minutes
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

## Future Enhancements

### Planned Improvements
1. GPU acceleration for eddy correction
2. Deep learning denoising methods
3. Automatic quality assessment scoring
4. Integration with BIDS data organization
5. Real-time processing monitoring
6. Automatic parameter optimization

### Research Features
1. Advanced distortion correction methods
2. Slice-to-volume motion correction
3. Multi-shell optimization
4. Cardiac/respiratory noise removal

---

*This documentation reflects the current state of the HINEC preprocessing pipeline. For the latest updates and detailed API documentation, consult the individual function files in the `nim_preprocessing/` directory.*