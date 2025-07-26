# HINEC Pipeline Overview

This document provides an overview of the HINEC pipeline, a workflow for processing and analyzing diffusion-weighted MRI (dMRI) data.

## Execution Flow

The pipeline is primarily executed through two main scripts: `runhinec.m` and `main.m`.

1.  **`runhinec.m`**: This is the main entry point for a user. It sets the necessary file paths for the input NIfTI data and the desired output `.mat` file. It then calls the `main.m` function to perform the core data processing. After the processing is complete, `runhinec.m` loads the resulting `.mat` file and generates various visualizations.

2.  **`main.m`**: This script contains the core logic of the pipeline. It orchestrates the entire workflow from data loading to final analysis, calling various sub-functions to perform specific tasks.

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

### 3. Visualization

After the core processing is finished and the `.mat` file is saved, the `runhinec.m` script generates a series of plots to visualize the results. This is accomplished by calling several plotting functions:

-   `nim_plotall`: Creates a comprehensive set of plots for the entire brain.
-   `nim_plotparcelall`: Generates plots specifically for each parcellated region.
-   `nim_plotparcellation`: Visualizes the parcellation mask itself overlaid on the brain.

---

In summary, the HINEC pipeline provides an end-to-end solution for dMRI analysis, transforming raw data into detailed, region-specific diffusion metrics and insightful visualizations.
