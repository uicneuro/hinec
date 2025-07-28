# HINEC Pipeline Overview

This document provides an overview of the HINEC pipeline, a workflow for processing and analyzing diffusion-weighted MRI (dMRI) data.

## Main Entry Points

The pipeline has three main entry points in the root directory:

1.  **`runhinec.m`**: Main entry point for DTI processing. Sets file paths and calls `main.m` for core processing, then generates visualizations.

2.  **`main.m`**: Core DTI processing pipeline. Handles data loading, tensor calculation, and parcellation.

3.  **`runTractography.m`**: Entry point for fiber tractography. Loads processed DTI data and runs tractography with visualization.

4.  **`visualizeTractography.m`**: Standalone visualization of saved tractography results without re-running tracking.

## Pipeline Stages

The HINEC pipeline can be broken down into three main stages:

### 1. Preprocessing

-   The `main.m` script first checks for the existence of a preprocessed NIfTI file.
-   If a preprocessed file is not found, it looks for a "raw" NIfTI file (e.g., `sample_raw.nii.gz`).
-   If the raw file is found, the `nim_preprocessing.m` function is automatically executed. This function leverages external neuroimaging tools (like FSL) to perform essential preprocessing steps, including:
    -   Motion correction
    -   Eddy current correction
    -   Brain extraction

### 2. Core Data Processing

Once the data is preprocessed, `main.m` executes the following sequence of operations:

-   **`nim_read`**: Loads the preprocessed NIfTI data, including the dMRI image, b-values, and b-vectors, into a MATLAB structure (`nim`).
-   **`nim_dt_spd`**: Calculates the diffusion tensor for each voxel. The diffusion tensor is a 3x3 matrix that models the local diffusion of water molecules in the brain tissue.
-   **`nim_fa`**: Computes the Fractional Anisotropy (FA) from the diffusion tensors. FA is a scalar value between 0 and 1 that measures the directionality of diffusion, which is particularly useful for characterizing white matter tracts.
-   **`nim_parcellation`**: Segments the brain into distinct anatomical regions based on a provided parcellation mask (e.g., an atlas). This enables the analysis of diffusion metrics within specific brain structures.
-   **`nim_load_labels`**: Loads the corresponding anatomical labels for the brain regions defined by the parcellation mask.
-   **`nim_save`**: Saves all the computed data—including the diffusion tensors, FA maps, parcellation information, and labels—into a single output `.mat` file for future use.

### 3. Tractography (Optional)

After DTI processing, fiber tractography can be performed:

-   **`runTractography.m`**: Loads processed DTI data and runs fiber tracking using the improved standard algorithm
-   **Automatic Saving**: Tractography results are automatically saved to `tractography_results/` with timestamp
-   **Integrated Visualization**: Creates comprehensive plots including tracks, FA background, and statistics

### 4. Visualization

**DTI Visualization** (via `runhinec.m`):
-   `nim_plotall`: Creates comprehensive DTI plots for the entire brain
-   `nim_plotparcelall`: Generates plots for each parcellated region  
-   `nim_plotparcellation`: Visualizes the parcellation mask

**Tractography Visualization**:
-   **`visualizeTractography.m`**: Standalone visualization of saved tractography results
-   **Automatic**: Can load most recent tractography results without arguments
-   **Comprehensive**: Shows 3D tracks, length distribution, and seed points

## Testing and Diagnostics

The `nim_tests/` directory contains diagnostic tools:

-   **`test_functions.m`**: Quick test script to verify tractography functions work
-   **`nim_diagnostic_check.m`**: Validates tensor quality and eigenvector orientation
-   **`nim_test_corpus_callosum.m`**: Tests tracking in known high-anisotropy region

## Data Flow

```
Raw DWI → main.m → DTI Data (.mat) → runTractography.m → Tracks (.mat) → visualizeTractography.m
                                  ↓
                              nim plots
```

---

In summary, the HINEC pipeline provides an end-to-end solution for dMRI analysis, from raw data to fiber tractography, with robust diagnostic tools and comprehensive visualization capabilities.
