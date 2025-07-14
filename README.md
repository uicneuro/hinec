# HINEC: HIgh-order NEural Connectivity

Human brain white matter tractography pipeline

## Instructions

1. Launch MATLAB and use `cd` in MATLAB to navigate to the HINEC directory.
2. Invoke the `main` function.

### Usage

#### Initialization

To run any command without running the main function, the paths must be added with the command `addpaths`.

#### Running HINEC

To run the program on the provided sample NiFTI file, and save it to "sample.mat", run the following in MATLAB.

```matlab
main nifti_sample/sample sample.mat
```

### Usage (Parcellation)

#### Data Preparation

To run hinec with parcellation, prepare the following files:

- raw NIFTI file
- bvec and bval file
- index file (index.txt)
- acquisition parameters (acqp.txt)

Name them all with same prefix. The naming scheme is as below:

- {name}\_raw.nii.gz
- {name}.bvec
- {name}.bval
- {name}\_index.txt
- {name}\_acqp.txt

#### Running HINEC

To run data processing only, use the following command:

```matlab
main {data location} {nim data file}
```

- `data location` refers to the location of the data along with the prefix
  - e.g. nifti_sample/parcellation_sample/sample
- `nim data file` refers to the location of save file generated (must be .mat)
  - e.g. sample.mat

To view the data by each parcellation, use the following commands:

```matlab
load {name}.mat
nim_plotparcelall(nim);
```

The total process can be executed by editing and running `runhinec`.

## Requirements

Addons:

- `Image Processing Toolbox`
- `Statistics and Machine Learning Toolbox`

External Softwares:

- `Statistical Parametric Mapping`
  - Must add folder SPM12 in the root directory
- `FSL`
  - Must be initialized before use
