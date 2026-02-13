# System Architecture and Data Structures

This document provides a technical overview of the HINEC (HIgh-order NEural Connectivity) pipeline architecture: how modules connect, the complete data flow, the central `nim` data structure, and dependencies between components.

---

## 1. System Overview

HINEC is a modular MATLAB pipeline that transforms raw diffusion-weighted MRI (dMRI) data into fiber tractography results and connectivity analysis. The pipeline follows a sequential processing model with optional stages.

```
Input (NIfTI)
    |
    v
+------------------+     +------------------+     +------------------+
|  Preprocessing   | --> | DTI Calculation  | --> |   Registration   |
|  (FSL-based)     |     | (Tensor + FA)    |     |  (DWI/T1/MNI)   |
+------------------+     +------------------+     +------------------+
                                                         |
    +----------------------------------------------------+
    |
    v
+------------------+     +------------------+     +------------------+
|  Parcellation    | --> |  Tractography    | --> |  Visualization   |
|  (Atlas-based)   |     | (FACT / HINEC)   |     | (3D/Slice/Fast)  |
+------------------+     +------------------+     +------------------+
```

### Entry Points

| Entry Point | Purpose | Typical Use |
|---|---|---|
| `main.m` | Core pipeline: preprocessing through parcellation | Process raw or preprocessed dMRI data |
| `runTractography.m` | Tractography execution with algorithm selection | Run fiber tracking on processed data |
| `runhinec.m` | Quick demo script | Load sample data and visualize |
| `bin/run_hinec.sh` | Shell wrapper for full pipeline | Batch processing with YAML config |
| `bin/run_tractography.sh` | Shell wrapper for tractography | Launch tractography jobs |

---

## 2. Directory Structure

```
hinec/
├── main.m                          # Core pipeline entry point
├── runTractography.m               # Tractography entry point
├── runhinec.m                      # Demo/quick-start script
├── CLAUDE.md                       # Project context for Claude Code
├── README.md                       # Project overview
├── requirements.txt                # Python dependencies (fast viewer)
│
├── config/                         # YAML configuration presets
│   ├── hinec_default.yml           #   Default balanced parameters
│   ├── high_precision.yml          #   RKF45 adaptive, publication quality
│   ├── fast_exploration.yml        #   RK2, quick testing
│   ├── standard_fact.yml           #   Baseline FACT algorithm
│   └── irontract.yml              #   IronTract challenge preset
│
├── src/                            # Source code modules
│   ├── nim_calculation/            #   DTI tensor, eigendecomposition, FA
│   │   ├── nim_dt_spd.m            #     Diffusion tensor (SPD-constrained)
│   │   ├── nim_dt.m                #     Diffusion tensor (basic LSF)
│   │   ├── nim_eig.m               #     Eigenvalue/eigenvector decomposition
│   │   └── nim_fa.m                #     Fractional anisotropy
│   │
│   ├── nim_preprocessing/          #   FSL-based preprocessing (16 files)
│   │   ├── nim_preprocessing.m     #     Orchestrator (10-step pipeline)
│   │   ├── preproc_denoising.m     #     MP-PCA / SUSAN / Gaussian
│   │   ├── preproc_extract_b0.m    #     B0 reference extraction
│   │   ├── preproc_motion_correction.m   #  MCFLIRT motion correction
│   │   ├── preproc_eddy_correction.m     #  Eddy current correction
│   │   ├── preproc_brain_extraction.m    #  BET brain extraction (DWI)
│   │   ├── preproc_t1_brain_extraction.m #  BET brain extraction (T1)
│   │   ├── preproc_fieldmap_correction.m #  FUGUE distortion correction
│   │   ├── preproc_mask_improvement.m    #  FA-based mask refinement
│   │   ├── preproc_tissue_segmentation.m #  WM/GM/CSF segmentation
│   │   ├── preproc_create_dwi_reference.m #  DWI reference volume
│   │   ├── preproc_t1_dwi_registration.m #  EPI_REG T1-DWI alignment
│   │   ├── preproc_t1_mni_registration.m #  FLIRT+FNIRT T1-MNI alignment
│   │   ├── preproc_atlas_resampling.m    #  Atlas to DWI space
│   │   ├── preproc_cleanup.m       #     Intermediate file removal
│   │   └── preprocessing_config_example.m # Config template
│   │
│   ├── nim_registration/           #   Multi-modal image registration
│   │   ├── nim_registration.m      #     Registration orchestrator
│   │   ├── register_dti_to_t1.m    #     DTI-to-T1 (FLIRT/SPM)
│   │   ├── register_t1_to_mni.m    #     T1-to-MNI (FLIRT+FNIRT/SPM)
│   │   ├── nim_apply_transforms.m  #     Transform chain application
│   │   └── registration_utils.m    #     Quality metrics, reporting
│   │
│   ├── nim_parcellation/           #   Brain region segmentation
│   │   ├── nim_parcellation.m      #     Basic parcellation (atlas/mask file)
│   │   └── nim_parcellation_registered.m  # Registration-enhanced parcellation
│   │
│   ├── nim_tractography/           #   Fiber tractography algorithms
│   │   ├── nim_tractography_standard.m    # FACT algorithm
│   │   ├── nim_tractography_hinec.m       # HINEC high-order tractography
│   │   ├── nim_tractography_highorder.m   # Legacy high-order methods
│   │   ├── nim_plot_tractography.m        # 3D track visualization
│   │   ├── nim_plot_tractography_region.m # Region-filtered visualization
│   │   ├── nim_plot_connectivity_matrix.m # Connectivity matrix
│   │   ├── nim_plot_vector_field.m        # Eigenvector field display
│   │   └── nim_excitation_time_map.m      # Excitation propagation
│   │
│   ├── nim_visualization/          #   Advanced tractography viewers
│   │   ├── visualizeTractography.m        # Unified 3D viewer (4 modes)
│   │   ├── visualizeTractographySlices.m  # Interactive 2D slice viewer
│   │   ├── visualizeTractographyAngles.m  # Angle-based visualization
│   │   ├── generateSlices.m               # Server-side slice generation
│   │   ├── generateTractographySliceCache.m # Cache generation pipeline
│   │   ├── TractographyCacheManager.m     # Cache metadata management
│   │   ├── buildOptimizedTrackSliceLookup.m # Vectorized slice lookup
│   │   ├── optimizedSliceRenderer.m       # High-performance renderer
│   │   ├── launchFastViewer.m             # MATLAB-to-Python bridge
│   │   └── testFastViewer.m              # Fast viewer test suite
│   │
│   ├── nim_plots/                  #   General plotting
│   │   └── nim_plot.m              #     DTI eigenvector visualization
│   │
│   ├── nim_presentation/           #   Research figure generation
│   │   ├── visualize_tractography_example.m
│   │   ├── visualize_integration_methods.m
│   │   ├── visualize_interpolation_methods.m
│   │   └── visualize_tractography_slice.m
│   │
│   ├── nim_utils/                  #   Utility functions (18 files)
│   │   ├── nim_read.m              #     NIfTI file reader
│   │   ├── nim_save.m              #     MAT file writer
│   │   ├── nim_load_nim.m          #     NIM structure loader
│   │   ├── nim_load_labels.m       #     Parcellation label loader
│   │   ├── nim_load_atlas_labels.m #     FSL atlas XML parser
│   │   ├── load_config_yaml.m      #     YAML configuration parser
│   │   ├── create_run_directory.m  #     Timestamped run organization
│   │   ├── nim_interp.m            #     High-order interpolation
│   │   ├── nim_reshape_d.m         #     6-vector to 3x3 tensor
│   │   ├── vector_to_color.m       #     Direction-to-RGB mapping
│   │   ├── nim_vis_eig.m           #     Eigenvector visualization helper
│   │   ├── nim_vis_fa.m            #     FA visualization helper
│   │   ├── gen_vis_eig.m           #     Eigenvector figure generation
│   │   ├── hdr.m                   #     NIfTI header utilities
│   │   ├── zwgll.m                 #     Gauss-Lobatto-Legendre nodes
│   │   ├── zwuni.m                 #     Uniform quadrature nodes
│   │   ├── plot_nim_interp.m       #     Interpolation result plotting
│   │   └── runnimplot.m            #     Plot workflow runner
│   │
│   └── nim_challenges/             #   Challenge/benchmark submissions
│       └── nim_irontract_submit.m  #     IronTract challenge formatter
│
├── bin/                            # Executable shell scripts
│   ├── run_hinec.sh                #   Full pipeline launcher
│   ├── run_tractography.sh         #   Tractography launcher
│   ├── run_generateSlices.sh       #   Slice cache generator
│   ├── viewSlices.sh               #   Fast Python viewer launcher
│   └── download.sh                 #   Sample data download
│
├── scripts/                        # Auxiliary scripts
│   ├── FastTractographyViewer.py   #   Python fast slice viewer GUI
│   ├── hinec_to_trk.py             #   HINEC-to-TRK format converter
│   ├── tractography_slice_gui.py   #   Slice viewer GUI (Python)
│   ├── validate_ismrm_tractography.py  # ISMRM validation
│   ├── diagnose_irontract.m        #   IronTract diagnostics
│   └── generate_presentation_figures.m # Presentation figures
│
├── tests/                          # Test suites
│   ├── test_yaml_config.m          #   YAML config validation
│   ├── test_irontract_submissions.m #  IronTract format tests
│   └── nim_tests/
│       ├── test_functions.m        #   Core function tests
│       ├── nim_diagnostic_check.m  #   Diagnostic checks
│       └── nim_test_corpus_callosum.m # Corpus callosum validation
│
├── data/                           # Sample datasets
│   ├── original_sample/            #   Basic diffusion data
│   └── parcellation_sample/        #   Data with parcellation masks
│
├── docs/                           # Documentation (19 existing files)
│
├── lib/                            # External libraries
│   ├── spm12/                      #   SPM12 toolbox (large)
│   └── bfgs/                       #   BFGS optimization
│       ├── bfgs_loop.m             #     Optimization loop
│       ├── bfgs_test.m             #     Optimization tests
│       └── vox_dt_bfgs.m           #     Voxel-wise tensor optimization
│
└── experiments/                    # Research outputs
    ├── hinec_runs/                 #   Pipeline run results
    ├── hinec_figures/              #   Generated figures
    ├── tractography_results/       #   Tractography outputs
    ├── tractography_figures/       #   Algorithm comparisons
    ├── validation_results/         #   Validation test results
    ├── ISMRM_data/                 #   ISMRM benchmark data
    ├── ironTract_data/             #   IronTract challenge data
    └── CovidNeuro_data/            #   COVID neuroimaging data
```

---

## 3. The `nim` Data Structure

The `nim` struct is the central data container that flows through the entire pipeline. Fields are progressively populated as each stage completes.

### Complete Field Reference

#### Metadata Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.hdr` | struct | — | NIfTI header (ImageSize, dime, etc.) | `nim_read` |
| `.xdim` | int | scalar | X dimension (columns) | `nim_read` |
| `.ydim` | int | scalar | Y dimension (rows) | `nim_read` |
| `.zdim` | int | scalar | Z dimension (slices) | `nim_read` |
| `.size3` | [1x3] | [xdim, ydim, zdim] | 3D volume dimensions | `nim_read` |
| `.run_info` | struct | — | Run directory metadata | `create_run_directory` |

#### Image Data Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.img` | double | [X, Y, Z, N] | 4D diffusion-weighted images (all volumes) | `nim_read` |
| `.img_b0` | double | [X, Y, Z, B0] | B0 (non-diffusion-weighted) volumes | `nim_read` |
| `.img_bi` | double | [X, Y, Z, Bi] | Diffusion-weighted volumes only | `nim_read` |

#### Acquisition Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.bval` | double | [1, N] | B-values for each volume | `nim_read` |
| `.bvec` | double | [3, N] | Gradient directions for each volume | `nim_read` |

#### Mask Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.mask` | logical | [X, Y, Z] | Brain mask (binary) | `nim_read` / preprocessing |
| `.wm_mask` | logical | [X, Y, Z] | White matter mask (FA > 0.2) | `preproc_tissue_segmentation` |
| `.gm_mask` | logical | [X, Y, Z] | Gray matter mask (0.05 < FA <= 0.2) | `preproc_tissue_segmentation` |
| `.csf_mask` | logical | [X, Y, Z] | CSF mask (FA <= 0.05) | `preproc_tissue_segmentation` |

#### DTI Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.DT` | double | [X, Y, Z, 6] | Diffusion tensor (6-element: Dxx,Dyy,Dzz,Dxy,Dxz,Dyz) | `nim_dt_spd` / `nim_dt` |
| `.evec` | double | [X, Y, Z, 3, 3] | Eigenvectors (sorted by eigenvalue) | `nim_eig` |
| `.eval` | double | [X, Y, Z, 3] | Eigenvalues (descending: lambda1 >= lambda2 >= lambda3) | `nim_eig` |
| `.FA` | double | [X, Y, Z] | Fractional anisotropy (0 to 1) | `nim_fa` |

#### Parcellation Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.parcellation_mask` | double | [X, Y, Z] | Integer labels per voxel (0 = background) | `nim_parcellation` |
| `.parcellation_mask_file` | char | — | Path to parcellation mask file | `nim_parcellation` |
| `.labels` | cell array | {N x 1} | Anatomical region names | `nim_load_labels` |
| `.atlas_labels` | containers.Map | — | Index-to-name label mapping | `nim_load_labels` |
| `.parcellation_info` | struct | — | Parcellation metadata and quality metrics | `nim_parcellation_registered` |

#### Registration Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.registration` | struct | — | Multi-modal registration data | `nim_registration` |
| `.registration.transforms` | struct | — | Transformation matrices (DTI/T1/MNI) | `register_dti_to_t1`, `register_t1_to_mni` |
| `.registration.quality_metrics` | struct | — | NMI scores and quality assessment | `registration_utils` |

### How the Structure Grows

```
nim_read           → .img, .img_b0, .img_bi, .bval, .bvec, .hdr, .mask, dimensions
    |
nim_dt_spd         → .DT (6-element tensor per voxel)
    |
nim_eig            → .evec (eigenvectors), .eval (eigenvalues)
    |
nim_fa             → .FA (fractional anisotropy map)
    |
nim_registration   → .registration (transforms, quality metrics)
    |
nim_parcellation   → .parcellation_mask, .labels, .atlas_labels
    |
tissue_segmentation → .wm_mask, .gm_mask, .csf_mask
    |
[Saved as .mat]    → Input to tractography
```

---

## 4. Module Dependency Graph

Arrows indicate "depends on" (must execute after).

```
nim_read
    |
    +---> nim_preprocessing (optional, if raw data)
    |         |
    |         +---> preproc_denoising
    |         +---> preproc_extract_b0
    |         +---> preproc_motion_correction
    |         +---> preproc_brain_extraction
    |         +---> preproc_fieldmap_correction (optional, if fieldmap available)
    |         +---> preproc_t1_brain_extraction (optional, if T1 available)
    |         +---> preproc_t1_dwi_registration (optional, if T1 available)
    |         +---> preproc_t1_mni_registration (optional, if T1 available)
    |         +---> preproc_atlas_resampling
    |         +---> preproc_cleanup
    |
    +---> nim_dt_spd (requires: .img, .bval, .bvec, .mask)
              |
              +---> nim_eig (requires: .DT)
                        |
                        +---> nim_fa (requires: .eval)
                                  |
    +-----------------------------+
    |
    +---> nim_registration (optional, requires: .FA, T1 file)
    |         |
    |         +---> register_dti_to_t1 (requires: b0 volume, T1)
    |         +---> register_t1_to_mni (requires: T1)
    |         +---> nim_apply_transforms
    |
    +---> nim_parcellation (requires: .FA or mask file)
    |         |
    |         +---> nim_parcellation_registered (if registration available)
    |
    +---> preproc_mask_improvement (requires: .FA, brain mask)
    |
    +---> preproc_tissue_segmentation (requires: .FA, brain mask)
    |
    +---> [nim saved to .mat]
              |
              +---> nim_tractography_standard (requires: .evec, .FA, .mask)
              |         OR
              +---> nim_tractography_hinec (requires: .evec, .FA, .mask, optional: tissue masks)
              |
              +---> visualizeTractography (requires: tracks, nim)
              +---> nim_plot_connectivity_matrix (requires: tracks, nim with parcellation)
              +---> generateSlices (requires: tracks, nim)
```

---

## 5. Pipeline Stages in Detail

### Stage 1: Data Loading (`nim_read.m`)

Reads NIfTI-1 format diffusion MRI data with associated acquisition parameters.

**Inputs**: NIfTI file path (without extension), optional brain mask and b-value/b-vector files.

**Processing**:
1. Load NIfTI image using SPM12's `load_nii()`
2. Extract header metadata (dimensions, voxel sizes)
3. Load b-values from `.bval` file and b-vectors from `.bvec` file
4. Separate b0 volumes (bval <= 50) from diffusion-weighted volumes
5. Load brain mask if provided (`.nii.gz` with `_M` suffix)

**Output**: Initialized `nim` struct with image data, acquisition parameters, and metadata.

### Stage 2: Preprocessing (`nim_preprocessing.m`)

Orchestrates a configurable 10-step FSL-based preprocessing pipeline. Automatically detects raw data (filename contains `_raw`) and triggers preprocessing.

**Steps** (each independently toggleable via config):

| Step | Function | FSL Tool | Purpose |
|---|---|---|---|
| 1 | `preproc_denoising` | dwidenoise/SUSAN | Remove thermal noise |
| 2 | `preproc_extract_b0` | fslroi | Extract b0 reference volume |
| 3 | `preproc_motion_correction` | mcflirt | Correct head motion artifacts |
| 4 | `preproc_eddy_correction` | eddy | Correct eddy current distortions |
| 5 | `preproc_brain_extraction` | BET | Generate brain mask |
| 6 | `preproc_fieldmap_correction` | FUGUE | Correct B0 field inhomogeneity |
| 7 | `preproc_t1_brain_extraction` | BET | Extract T1 brain (if T1 provided) |
| 8 | `preproc_t1_dwi_registration` | epi_reg | Register T1 to DWI space |
| 9 | `preproc_t1_mni_registration` | FLIRT+FNIRT | Register T1 to MNI space |
| 10 | `preproc_atlas_resampling` | applywarp/flirt | Resample atlas to DWI space |

**Output**: Preprocessed NIfTI file, brain mask, and parcellation mask in DWI space.

### Stage 3: DTI Calculation

Three sequential functions compute the diffusion tensor model:

1. **`nim_dt_spd.m`**: Fits diffusion tensors using least squares, then enforces positive-definiteness via BFGS optimization with Cholesky parameterization. Stores 6-element tensor representation in `nim.DT`.

2. **`nim_eig.m`**: Eigendecomposition of each voxel's 3x3 tensor. Sorts eigenvalues in descending order (lambda1 >= lambda2 >= lambda3). Stores eigenvectors and eigenvalues.

3. **`nim_fa.m`**: Computes fractional anisotropy from eigenvalues. Values range 0 (isotropic) to 1 (perfectly anisotropic).

### Stage 4: Registration (`nim_registration.m`, optional)

Multi-modal registration chain when a T1 anatomical image is available:

```
DWI native space  <--FLIRT 6DOF-->  T1 native space  <--FLIRT+FNIRT-->  MNI standard space
```

- DTI-to-T1: Uses FLIRT (correlation ratio) or SPM coregistration
- T1-to-MNI: Linear (FLIRT 12 DOF) + optional nonlinear (FNIRT)
- Quality assessed via normalized mutual information (NMI)

### Stage 5: Parcellation (`nim_parcellation.m`)

Atlas-based brain region segmentation. Two pathways:

1. **Basic**: Loads a provided parcellation mask file or generates one using direct atlas registration
2. **Registration-enhanced** (`nim_parcellation_registered.m`): Transforms MNI atlas through T1 to DWI space using the registration chain, producing more accurate parcellation

**Supported Atlases**: JHU white matter labels, JHU white matter tracts, HarvardOxford cortical, AAL.

### Stage 6: Tissue Segmentation (`preproc_tissue_segmentation.m`)

FA-based tissue classification for Anatomically Constrained Tractography (ACT):

| Tissue | FA Criterion | Post-processing |
|---|---|---|
| White matter (WM) | FA > 0.2 | Eroded (remove boundary voxels) |
| Gray matter (GM) | 0.05 < FA <= 0.2 | None |
| CSF | FA <= 0.05 | None |

Validates mutual exclusivity and expected tissue proportions (WM ~40-45%, GM ~40-45%, CSF ~10-20%).

### Stage 7: Tractography

Two primary algorithms are available:

**Standard FACT** (`nim_tractography_standard.m`):
- Discrete voxel-by-voxel tracking without interpolation
- Variable step size determined by voxel boundary intersection
- Bidirectional tracking from seed points
- Seeding: uniform lattice with optional random jitter

**HINEC High-Order** (`nim_tractography_hinec.m`):
- Trilinear interpolation via pre-computed `griddedInterpolant` objects
- Multiple integration methods: Euler, RK2, RK4, RKF45 (adaptive)
- Anatomically Constrained Tractography (ACT) with tissue-type termination
- Parallel processing support (`parfor` with `DataQueue` progress)

Both algorithms share these stopping criteria:
- FA below threshold (typically 0.1-0.2)
- Angle between consecutive steps exceeds threshold (typically 35-60 degrees)
- Track exits brain mask
- Maximum step count reached

**Seeding strategy** (3-tier priority):
1. Preprocessed brain mask (best)
2. Expanded parcellation mask (dilated regions)
3. FA-threshold fallback (FA > 0.10)

**Output**: Cell array of tracks (each track is Nx3 matrix of 3D coordinates), options struct, timing, and statistics.

### Stage 8: Visualization

Three tiers of visualization capability:

1. **MATLAB 3D** (`visualizeTractography.m`): Unified viewer with whole-brain, region, grid, and sequential modes. Supports direction/FA/region coloring and export.

2. **MATLAB 2D Slices** (`visualizeTractographySlices.m`): Interactive orthogonal slice viewer with crosshair synchronization. 5-30 second latency per slice update.

3. **Fast Distributed Viewer**: Pre-generate slice cache on server (`generateSlices.m`), transfer, view locally with Python GUI (`FastTractographyViewer.py`). Sub-100ms slice transitions.

See [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md) for complete details.

---

## 6. Configuration System

HINEC uses YAML configuration files for reproducible parameter management.

### Configuration Structure

```yaml
preprocessing:
  denoising: true/false
  motion_correction: true/false
  eddy_correction: true/false
  brain_extraction: true/false
  fieldmap_correction: true/false

tractography:
  algorithm: standard | hinec
  integration_order: 1 (Euler) | 2 (RK2) | 4 (RK4) | 5 (RKF45)
  interpolation: none | trilinear | cubic | spline
  step_size: 0.1 - 1.0 (voxel units)
  angle_threshold: 20 - 90 (degrees)
  fa_threshold: 0.05 - 0.3
  seed_density: 1 - 10 (seeds per voxel dimension)
  adaptive_step: true/false (RKF45 only)
  adaptive_tolerance: 0.001 - 0.1
  act_enabled: true/false
```

### Available Presets

| Preset | Algorithm | Integration | Step Size | Seed Density | Use Case |
|---|---|---|---|---|---|
| `hinec_default.yml` | hinec | RK4 | 0.2 | 4 | Balanced default |
| `high_precision.yml` | hinec | RKF45 | 0.2 | 4 | Publication quality |
| `fast_exploration.yml` | hinec | RK2 | 0.3 | 2 | Quick testing |
| `standard_fact.yml` | standard | Euler | 0.5 | 4 | Baseline FACT |
| `irontract.yml` | hinec | RK4 | 0.2 | 4 | IronTract challenge |

See [YAML_CONFIG.md](YAML_CONFIG.md) for complete configuration reference.

---

## 7. Run Directory System

Each pipeline execution creates a timestamped directory for reproducibility:

```
hinec_runs/
└── run_YYYY-MM-DD_HH-MM-SS/
    ├── run_info.txt              # Metadata (git hash, timestamp, description)
    ├── config.yml                # Copied configuration file
    ├── logs/                     # Execution logs
    ├── intermediate/             # Intermediate processing files
    ├── output/                   # Final processed nim.mat
    ├── tractography/             # Track results
    └── diagnostics/              # Quality reports
```

A `latest` symlink points to the most recent run. See [RUN_DIRECTORY_SYSTEM.md](RUN_DIRECTORY_SYSTEM.md) for details.

---

## 8. External Dependencies

### FSL (FMRIB Software Library)

Required for all preprocessing operations. Must be installed and initialized (`FSLDIR` set) before running.

| FSL Tool | Used By | Purpose |
|---|---|---|
| `mcflirt` | `preproc_motion_correction` | Motion correction |
| `eddy` | `preproc_eddy_correction` | Eddy current correction |
| `bet` | `preproc_brain_extraction` | Brain extraction |
| `fugue` | `preproc_fieldmap_correction` | Distortion correction |
| `flirt` | Registration functions | Linear registration (6/12 DOF) |
| `fnirt` | `register_t1_to_mni` | Nonlinear registration |
| `epi_reg` | `preproc_t1_dwi_registration` | EPI-to-structural registration |
| `applywarp` | `preproc_atlas_resampling` | Apply warp fields |
| `invwarp` | `preproc_t1_mni_registration` | Invert warp fields |
| `fslroi` | `preproc_extract_b0` | Volume extraction |
| `fslmaths` | Multiple functions | Arithmetic, morphology |
| `fslcpgeom` | `preproc_atlas_resampling` | Copy geometry headers |

### SPM12

Included in `lib/spm12/`. Used for:
- NIfTI file I/O (`load_nii`, `save_nii`)
- Optional registration method (alternative to FSL)
- Image coregistration and normalization

### BFGS Optimization Library

Included in `lib/bfgs/`. Used by `nim_dt_spd.m` for enforcing positive-definiteness of diffusion tensors through constrained optimization.

### MATLAB Toolboxes

| Toolbox | Required By | Purpose |
|---|---|---|
| Image Processing Toolbox | `preproc_denoising` (nlmeans), visualization | Image filtering and display |
| Statistics and Machine Learning | DTI calculations | Statistical computations |
| Tools for NIfTI and ANALYZE image | `nim_read`, `nim_save` | NIfTI file format support |

### Python (Optional)

Required only for the fast distributed slice viewer:
- Python 3.7+
- tkinter (GUI framework)
- Pillow (image loading)
- numpy (array operations)

---

## 9. File Naming Conventions

### Input Files

| Pattern | Description | Example |
|---|---|---|
| `{name}_raw.nii.gz` | Raw diffusion data (triggers preprocessing) | `sample_raw.nii.gz` |
| `{name}.nii.gz` | Preprocessed diffusion data | `sample.nii.gz` |
| `{name}.bval` | B-values (space-separated, one line) | `sample.bval` |
| `{name}.bvec` | B-vectors (3 rows x N columns) | `sample.bvec` |
| `{name}_acqp.txt` | Acquisition parameters (for eddy) | `sample_acqp.txt` |
| `{name}_index.txt` | Volume index file (for eddy) | `sample_index.txt` |
| `{name}_M.nii.gz` | Brain mask | `sample_M.nii.gz` |

### Output Files

| Pattern | Description | Example |
|---|---|---|
| `{name}.mat` | Processed nim structure | `sample_parcellated.mat` |
| `tracks_{algo}_{timestamp}.mat` | Tractography results | `tracks_hinec_2025-01-15_14_30_00.mat` |

### Intermediate Files (cleaned up after preprocessing)

| Pattern | Description |
|---|---|
| `b0.nii.gz` | Extracted b0 volume |
| `nodif_brain.nii.gz` | Brain-extracted b0 |
| `nodif_brain_mask.nii.gz` | Brain mask from BET |
| `parcellation_mask.nii.gz` | Atlas in DWI space |

---

## Cross-References

- Configuration details: [YAML_CONFIG.md](YAML_CONFIG.md)
- Run directory system: [RUN_DIRECTORY_SYSTEM.md](RUN_DIRECTORY_SYSTEM.md)
- Preprocessing details: [PREPROCESSING.md](PREPROCESSING.md)
- Tractography algorithms: [TRACTOGRAPHY.md](TRACTOGRAPHY.md)
- Scientific foundations: [SCIENCE.md](SCIENCE.md)
- Mathematical methods: [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md)
- Complete function reference: [API_REFERENCE.md](API_REFERENCE.md)
- Visualization guide: [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)
