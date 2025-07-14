# HINEC: HIgh-order NEural Connectivity

## What HINEC Does

HINEC is a brain white matter tractography analysis pipeline for processing and analyzing diffusion MRI data. It enables researchers to:

- Process raw diffusion MRI data into a format suitable for analysis
- Calculate diffusion tensors and their properties (like fractional anisotropy)
- Apply brain parcellation to analyze connectivity between specific brain regions
- Visualize diffusion tensor data and brain connectivity metrics
- Extract quantitative measurements from diffusion data across brain regions

The program facilitates neuroimaging research by providing tools to investigate structural brain connectivity and white matter properties.

## Tech Stack & Requirements

### Core Technology
- **MATLAB** - Primary programming environment

### Required MATLAB Toolboxes
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

### External Software Dependencies
- **SPM12** (Statistical Parametric Mapping)
  - Must be added to the root directory
- **FSL** (FMRIB Software Library)
  - Must be initialized before use

### Input Data Format
- NiFTI format (.nii.gz) diffusion MRI data
- Associated metadata files (.bval, .bvec, .json)
- Parcellation masks for region-based analysis

## Pipeline Overview

### 1. Data Preparation
- Organize raw diffusion data with consistent naming convention
- Ensure required files exist: raw NiFTI, bvec, bval, index, and acquisition parameters

### 2. Preprocessing
- `nim_preprocessing.m` handles:
  - Motion and eddy current correction
  - Brain extraction
  - Generation of brain mask
  - Preparation for diffusion tensor calculation

### 3. Diffusion Tensor Calculation
- `nim_dt_spd.m` calculates the diffusion tensor for each voxel
- `nim_fa.m` computes fractional anisotropy maps

### 4. Parcellation Analysis
- `nim_parcellation.m` applies brain region masks
- `nim_extract_parcellation.m` extracts metrics by brain region
- Region labels are loaded with `nim_load_labels.m`

### 5. Visualization
- `nim_plotparcelall.m` displays metrics by parcellation
- `nim_plotparcellation.m` visualizes the parcellation system
- Various plot functions for diffusion data visualization

### 6. Data Export
- Results saved as MATLAB .mat files
- Contains processed data, calculated metrics, and parcellation information

The complete pipeline can be executed through `runhinec.m` or by using the `main` function with appropriate parameters.