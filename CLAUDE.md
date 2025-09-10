# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HINEC (HIgh-order NEural Connectivity) is a MATLAB-based pipeline for processing and analyzing diffusion-weighted MRI (dMRI) data. The pipeline processes raw NIfTI diffusion MRI data through preprocessing, diffusion tensor calculation, fractional anisotropy computation, parcellation, and tractography analysis.

## Key Commands

### Main Execution
- `main('path/to/data', 'output.mat')` - Core pipeline execution
- `runhinec` - Complete processing workflow with visualization
- `runTractography('data.mat')` - Run tractography analysis on processed data

### Data Processing Functions
- `nim_read(imgpath)` - Load NIfTI diffusion data
- `nim_dt_spd(nim)` - Calculate diffusion tensors (SPD-constrained)
- `nim_fa(nim)` - Compute fractional anisotropy
- `nim_parcellation(nim, mask_file)` - Apply brain parcellation
- `nim_tractography_standard(nim, options)` - Standard tractography

### Visualization
- `nim_plotall(nim)` - Comprehensive visualization
- `nim_plotparcelall(nim)` - Parcellation-specific plots
- `nim_plotparcellation(nim)` - Parcellation mask visualization

## Architecture

### Core Modules
- **nim_calculation/**: Diffusion tensor and FA calculations
- **nim_parcellation/**: Brain region segmentation
- **nim_preprocessing/**: Raw data preprocessing (FSL integration)
- **nim_tractography/**: Fiber tractography algorithms
- **nim_plots/**: Visualization functions
- **nim_utils/**: Utility functions for data I/O and manipulation

### Data Flow
1. **Input**: Raw NIfTI files with bval/bvec (e.g., `sample_raw.nii.gz`)
2. **Preprocessing**: Motion/eddy correction via FSL (if raw data provided)
3. **Processing**: Diffusion tensor → Eigenvalues/vectors → FA calculation
4. **Parcellation**: Atlas-based brain region segmentation
5. **Output**: Processed `.mat` file with all computed metrics
6. **Analysis**: Tractography and visualization

### File Naming Conventions
- Raw data: `{name}_raw.nii.gz`, `{name}.bval`, `{name}.bvec`
- Acquisition params: `{name}_acqp.txt`, `{name}_index.txt`
- Processed: `{name}.nii.gz` (after preprocessing)
- Output: User-specified `.mat` file

## Dependencies

### MATLAB Toolboxes (Required)
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Tools for NIfTI and ANALYZE image

### External Software
- **SPM12**: Must be in `spm12/` directory (included in repo)
- **FSL**: Required for preprocessing, must be initialized before use

## Development Notes

### Testing
- SPM12 includes comprehensive test suite in `spm12/tests/`
- No dedicated HINEC test framework present
- Manual testing through sample data in `nifti_sample/`

### Sample Data
- `nifti_sample/original_sample/`: Basic diffusion data
- `nifti_sample/parcellation_sample/`: Data with parcellation masks
- Pre-computed results: `sample_parcellated.mat`

### Tractography
- Standard and high-order methods available
- Configurable parameters: step size, FA threshold, angle threshold
- Output includes track coordinates and statistics

## Important Implementation Details

### Data Structure (nim)
The main data structure contains:
- `.img`: 4D diffusion image data
- `.evec`: Eigenvectors from diffusion tensors
- `.eval`: Eigenvalues
- `.FA`: Fractional anisotropy maps
- `.parcellation_mask`: Brain region labels
- `.labels`: Anatomical region names

### Path Management
The `main.m` function automatically adds all necessary paths:
- All nim_* subdirectories
- SPM12 with subpaths
- Utility directories

### Preprocessing Integration
- Automatic detection of raw vs. preprocessed data
- FSL preprocessing called when raw data detected
- Preprocessing includes motion correction, eddy current correction, brain extraction