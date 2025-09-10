# HINEC: HIgh-order NEural Connectivity

Human brain white matter tractography pipeline

## Instructions

1.  Launch MATLAB and use `cd` in MATLAB to navigate to the HINEC directory.
2.  Invoke the `main` function.

### Usage

#### Initialization

To run any command without running the main function, the paths must be added with the command `addpaths`.

#### Data Preparation

To run HINEC from scratch, you must provide a raw NIfTI file (`.nii.gz`), a b-vector file (`.bvec`), and a b-value file (`.bval`). These files should share the same prefix and be organized as follows:

*   `{prefix}_raw.nii.gz`
*   `{prefix}.bvec`
*   `{prefix}.bval`

For example, if your prefix is `my_data`, you should have:

*   `my_data_raw.nii.gz`
*   `my_data.bvec`
*   `my_data.bval`

#### Running HINEC

To process your data and save the output to a `.mat` file, run the following command in MATLAB, replacing `{data_location}` with the path to your files and `{prefix}` with your chosen prefix:

```matlab
main('{data_location}/{prefix}', 'output.mat')
```

For instance, if your data is in a folder named `input_data`, the command would be:

```matlab
main('input_data/my_data', 'output.mat')
```

The pipeline will automatically perform all necessary preprocessing steps, including brain extraction, parcellation, and diffusion tensor calculation, and save the final results to `output.mat`.

### Viewing Parcellation Results

To view the data by each parcellation, use the following commands:

```matlab
load('output.mat');
nim_plotparcelall(nim);
```

## Requirements

Addons:

-   `Image Processing Toolbox`
-   `Statistics and Machine Learning Toolbox`
-   `Tools for NIfTI and ANALYZE image`

External Softwares:

-   `Statistical Parametric Mapping`
    -   Must add folder SPM12 in the root directory
-   `FSL`
    -   Must be initialized before use
