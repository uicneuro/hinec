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
|  (FSL-based)     |     | (Tensor + FA)    |     |  (DWI/T1/MNI)    |
+------------------+     +------------------+     +------------------+
                                                         |
    +----------------------------------------------------+
    |
    v
+------------------+     +---------------------+  +------------------+
|  Parcellation    | --> |    Tractography     |->|  Visualization   |
|  (Atlas-based)   |     | (FACT/HINEC/MMF)    |  | (3D/Slice/Fast)  |
+------------------+     +---------------------+  +------------------+
```

### Entry Points

| Entry Point | Purpose | Typical Use |
|---|---|---|
| `main.m` | Core pipeline: preprocessing through tissue segmentation | Process raw or preprocessed dMRI data |
| `runTractography.m` | Seeding, tracker dispatch, ROI filtering | Run fiber tracking on a processed `nim` |
| `runhinec.m` | Scratch script | Load a processed `.mat` and plot it |
| `bin/run_hinec.sh` | Full pipeline launcher | Preprocess and track in one background job |
| `bin/run_tractography.sh` | Tractography-only launcher | Iterate tracking configs on a built `nim` |
| `bin/run_ismrm_scoring.sh` | Scoring | Convert tracks to TRK and run the scilpy scorer |
| `bin/run_visualization.sh` | Headless figure export | Render figures into `<run_dir>/figures/` |

---

## 2. Directory Structure

```
hinec/
├── main.m                          # Core pipeline entry point
├── runTractography.m               # Tractography entry point
├── runhinec.m                      # Scratch plotting script
├── README.md                       # Project overview
├── mkdocs.yml                      # Documentation site configuration
├── requirements.txt                # Python dependencies (viewer, scoring, TRK export)
│
├── config/                         # YAML configuration presets
│   ├── README.md                   #   Naming conventions
│   ├── hinec_default.yml           #   Generic full-pipeline fallback
│   ├── standard_dti.yml            #   FACT baseline
│   ├── hinec_dti.yml               #   Interpolated streamlines, DTI field
│   ├── hinec_csd.yml               #   Interpolated streamlines, CSD field
│   ├── hinec_dti_cubic.yml         #   Cubic interpolation + RKF45
│   ├── hinec_dti_cubic_recall.yml  #   Cubic, recall-tuned termination
│   ├── hinec_dti_fast.yml          #   RK2, quick screening
│   ├── hinec_dti_euler.yml         #   Euler, lowest-order comparison point
│   ├── mmf_dti.yml                 #   Connection-form MMF, DTI field
│   ├── mmf_csd.yml                 #   Connection-form MMF, CSD field
│   ├── ismrm2015.yml               #   ISMRM 2015 challenge dataset
│   ├── irontract.yml               #   IronTract challenge dataset
│   └── reference.yml               #   Convergence-ladder template
│
├── src/                            # Source code modules
│   ├── nim_calculation/            #   Tensor, eigendecomposition, FA, CSD, MMF geometry
│   │   ├── nim_dt_spd.m            #     Diffusion tensor (SPD-constrained)
│   │   ├── nim_dt.m                #     Diffusion tensor (basic LSF)
│   │   ├── nim_eig.m               #     Eigenvalue/eigenvector decomposition
│   │   ├── nim_fa.m                #     Fractional anisotropy
│   │   ├── nim_csd.m               #     Constrained spherical deconvolution + FOD peaks
│   │   ├── nim_build_frames.m      #     Moving-frame construction
│   │   ├── nim_connection_form.m   #     Connection 1-form (curvature, torsion)
│   │   └── nim_mmf_geometry.m      #     Frame field + connection baked into the nim
│   │
│   ├── nim_preprocessing/          #   FSL-based preprocessing
│   │   ├── nim_preprocessing.m     #     Orchestrator (10-step pipeline)
│   │   ├── preproc_denoising.m     #     dwidenoise / nlmeans / Gaussian
│   │   ├── preproc_extract_b0.m    #     B0 reference extraction
│   │   ├── preproc_motion_correction.m   #  MCFLIRT motion correction
│   │   ├── preproc_eddy_correction.m     #  Eddy current correction
│   │   ├── preproc_brain_extraction.m    #  BET brain extraction (DWI)
│   │   ├── preproc_t1_brain_extraction.m #  BET brain extraction (T1)
│   │   ├── preproc_fieldmap_correction.m #  FUGUE distortion correction
│   │   ├── preproc_mask_improvement.m    #  FA-based mask refinement
│   │   ├── preproc_tissue_segmentation.m #  WM/GM/CSF via FAST-on-T1, FA fallback
│   │   ├── preproc_create_dwi_reference.m #  DWI reference volume
│   │   ├── preproc_t1_dwi_registration.m #  epi_reg T1-DWI alignment
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
│   │   ├── extract_reference_volumes.m    # Reference volume extraction
│   │   ├── compute_registration_quality.m # Quality assessment
│   │   ├── compute_normalized_mutual_information.m # NMI metric
│   │   ├── generate_registration_report.m # Reporting
│   │   └── save_registration_data.m       # Persistence
│   │
│   ├── nim_parcellation/           #   Brain region segmentation
│   │   ├── nim_parcellation.m      #     Basic parcellation (atlas/mask file)
│   │   └── nim_parcellation_registered.m  # Registration-enhanced parcellation
│   │
│   ├── nim_tractography/           #   Fiber tractography
│   │   ├── nim_tractography_standard.m    # FACT
│   │   ├── nim_tractography_hinec.m       # Interpolated streamlines
│   │   ├── nim_tractography_mmf_connframe.m # Connection-form MMF
│   │   ├── nim_tractography_highorder.m   # Legacy high-order methods
│   │   ├── mmf_bishop_update.m            # Bishop frame update
│   │   ├── mmf_gram_schmidt.m             # Frame reorthonormalization
│   │   ├── mmf_reference_axis_frame.m     # Reference-axis frame
│   │   ├── nim_seed_offsets.m             # Sub-voxel seed lattice
│   │   ├── nim_filter_tracks_roi.m        # Waypoint / endpoint / containment gates
│   │   ├── nim_resample_track_arc.m       # Arc-length resampling
│   │   ├── nim_sift.m                     # Track filtering
│   │   ├── nim_plot_tractography.m        # 3D track visualization
│   │   ├── nim_plot_tractography_region.m # Region-filtered visualization
│   │   ├── nim_plot_connectivity_matrix.m # Connectivity matrix
│   │   ├── nim_plot_vector_field.m        # Eigenvector field display
│   │   └── nim_excitation_time_map.m      # Excitation propagation
│   │
│   ├── nim_visualization/          #   Tractography viewers and slice cache
│   │   ├── visualizeTractography.m        # Unified 3D viewer (4 modes)
│   │   ├── visualizeTractographySlices.m  # Interactive 2D slice viewer
│   │   ├── visualizeTractographyAngles.m  # Angle-based visualization
│   │   ├── nim_plot_bundles.m             # Per-bundle rendering
│   │   ├── nim_plot_vs_groundtruth.m      # Tracks against ground truth
│   │   ├── generateSlices.m               # Server-side slice generation
│   │   ├── generateTractographySliceCache.m # Cache generation pipeline
│   │   ├── TractographyCacheManager.m     # Cache metadata management
│   │   ├── buildOptimizedTrackSliceLookup.m # Vectorized slice lookup
│   │   ├── optimizedSliceRenderer.m       # Slice renderer
│   │   ├── launchFastViewer.m             # MATLAB-to-Python bridge
│   │   └── testFastViewer.m               # Fast viewer test suite
│   │
│   ├── nim_plots/                  #   General plotting
│   │   └── nim_plot.m              #     Consolidated eigenvector/parcel plotter
│   │
│   ├── nim_presentation/           #   Research figure generation
│   │   ├── visualize_tractography_example.m
│   │   ├── visualize_integration_methods.m
│   │   ├── visualize_interpolation_methods.m
│   │   └── visualize_tractography_slice.m
│   │
│   ├── nim_utils/                  #   Utilities
│   │   ├── nim_read.m              #     NIfTI file reader
│   │   ├── nim_save.m              #     MAT file writer (v7.3 for large structs)
│   │   ├── nim_load_nim.m          #     nim structure loader
│   │   ├── nim_load_labels.m       #     Parcellation label loader
│   │   ├── nim_load_atlas_labels.m #     FSL atlas XML parser
│   │   ├── nim_atlas_label_map.m   #     Label index <-> name mapping
│   │   ├── nim_parcellation_from_masks.m # Parcellation from binary mask volumes
│   │   ├── nim_roi_mask.m          #     Resolve ROI names/indices to a mask
│   │   ├── nim_nifti_affine.m      #     sform/qform affine extraction
│   │   ├── nim_read_trk.m          #     TRK reader
│   │   ├── nim_config_schema.m     #     SINGLE SOURCE OF TRUTH for all parameters
│   │   ├── nim_yaml_parse.m        #     YAML parser (two-level nesting)
│   │   ├── load_config_yaml.m      #     Defaults, validation, legacy migration
│   │   ├── nim_config_to_options.m #     Nested config -> flat tracker options
│   │   ├── nim_config_apply_overrides.m # CLI --set resolution
│   │   ├── nim_config_write.m      #     Config serialization
│   │   ├── nim_config_docs.m       #     Generates docs/YAML_CONFIG.md from the schema
│   │   ├── nim_config_retired.m    #     Keys that are accepted and dropped
│   │   ├── nim_angle_limit.m       #     Curvature budget from angle_max and step
│   │   ├── nim_principal_dir.m     #     Principal direction extraction
│   │   ├── nim_convergence_error.m #     Track-to-reference error metrics
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
│   ├── run_tractography.sh         #   Tractography-only launcher
│   ├── run_ismrm_scoring.sh        #   TRK conversion + scilpy scorer
│   ├── run_visualization.sh        #   Headless figure export
│   ├── run_generateSlices.sh       #   Slice cache generator
│   ├── run_sift.sh                 #   Track filtering
│   ├── viewSlices.sh               #   Fast Python viewer launcher
│   └── download.sh                 #   Sample data download
│
├── scripts/                        # Auxiliary scripts
│   ├── FastTractographyViewer.py   #   Python fast slice viewer GUI
│   ├── tractography_slice_gui.py   #   Slice viewer GUI (Python)
│   ├── hinec_to_trk.py             #   HINEC-to-TRK format converter
│   ├── validate_ismrm_tractography.py  # ISMRM validation
│   ├── build_ismrm_scoring_config.py   # Scorer configuration
│   ├── compare_ismrm_results.py    #   Cross-run score comparison
│   ├── run_pft_dipy.py             #   DIPY PFT reference tracking
│   ├── diagnose_irontract.m        #   IronTract diagnostics
│   └── generate_presentation_figures.m # Presentation figures
│
├── tests/                          # Test suites
│   ├── unit/                       #   Unit tests (config schema, angle limit,
│   │                               #   frames, interpolation, seeding, resampling...)
│   ├── integration/                #   TestFullPipeline.m
│   ├── fsl/                        #   TestFSLPreprocessing.m
│   ├── nim_tests/                  #   Quick diagnostics
│   ├── fixtures/                   #   Synthetic phantom and nim builders
│   ├── test_yaml_config.m          #   Exercises every YAML preset
│   └── test_irontract_submissions.m #  IronTract format tests
│
├── data/                           # Datasets (large binaries are gitignored)
│   ├── original_sample/            #   Preprocessed diffusion sample
│   ├── parcellation_sample/        #   Raw diffusion sample
│   ├── ismrm2015/                  #   ISMRM 2015 challenge data + canonical nim
│   ├── irontract/                  #   IronTract challenge data
│   └── covid/                      #   COVID neuroimaging data
│
├── docs/                           # MkDocs documentation source
│
├── hinec_runs/                     # Timestamped run outputs (gitignored)
│
└── lib/                            # External libraries
    ├── spm12/                      #   SPM12 — install here; not vendored
    └── bfgs/                       #   BFGS optimization (vendored)
        ├── bfgs_loop.m             #     Optimization loop
        ├── bfgs_test.m             #     Optimization tests
        └── vox_dt_bfgs.m           #     Voxel-wise tensor optimization
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
| `.img_b0` | double | [X, Y, Z] | Mean of the b0 (non-diffusion-weighted) volumes | `nim_read` |
| `.img_bi` | double | [X, Y, Z, Bi] | Diffusion-weighted volumes only | `nim_read` |
| `.size_b0`, `.size_bi` | int | scalar | Volume counts either side of the b0 threshold | `nim_read` |
| `.thrsh_b0` | double | scalar | b-value below which a volume counts as b0 (default 5) | `nim_read` |

#### Acquisition Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.bval` | double | [1, N] | B-values for each volume | `nim_read` |
| `.bvec` | double | [3, N] | Gradient directions for each volume | `nim_read` |

#### Mask Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.mask` | logical | [X, Y, Z] | Brain mask (binary) | `nim_read` / preprocessing |
| `.wm_mask` | double | [X, Y, Z] | White matter mask for ACT | `preproc_tissue_segmentation` |
| `.gm_mask` | double | [X, Y, Z] | Gray matter mask for ACT | `preproc_tissue_segmentation` |
| `.csf_mask` | double | [X, Y, Z] | CSF mask for ACT | `preproc_tissue_segmentation` |
| `.wm_mask_file`, `.gm_mask_file`, `.csf_mask_file` | char | — | Paths to the mask NIfTIs | `main.m` |

The tissue masks are always generated. They come from FSL FAST on the anatomical T1 where one
is available, and from FA-tertile binning otherwise; the fallback bins anisotropy rather than
tissue and should be treated as degraded. ACT is nonetheless off unless `tractography.act` is
set true.

#### DTI Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.DT` | double | [X, Y, Z, 6] | Diffusion tensor (6-element: Dxx,Dyy,Dzz,Dxy,Dxz,Dyz) | `nim_dt_spd` / `nim_dt` |
| `.evec` | double | [X, Y, Z, 3, 3] | Eigenvectors (sorted by eigenvalue) | `nim_eig` |
| `.eval` | double | [X, Y, Z, 3] | Eigenvalues (descending: lambda1 >= lambda2 >= lambda3) | `nim_eig` |
| `.FA` | double | [X, Y, Z] | Fractional anisotropy (0 to 1) | `nim_fa` |

#### CSD Fields (`field: csd` only)

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.peaks` | double | [X, Y, Z, P, 3] | Unit FOD peak directions | `nim_csd` |
| `.npeaks` | double | [X, Y, Z] | Number of accepted peaks per voxel | `nim_csd` |
| `.peak_w` | double | [X, Y, Z, P] | Peak amplitudes | `nim_csd` |
| `.fod_sh` | double | — | Spherical harmonic FOD coefficients | `nim_csd` |
| `.response` | — | — | Estimated response function | `nim_csd` |

These are computed by `runTractography` on demand and cached as `<source>_csd.mat`.

#### Moving-Frame Fields (`algorithm: mmf`)

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.mmf_frames` | double | [X, Y, Z, 3, 3] | Frame field {e1, e2, e3} | `nim_mmf_geometry` |
| `.mmf_kappa` | double | [X, Y, Z, 3] | Curvature vector de1/ds | `nim_mmf_geometry` |
| `.mmf_tau` | double | [X, Y, Z] | Torsion, omega_23(e1) | `nim_mmf_geometry` |
| `.mmf_field` | char | — | Which field the geometry was built from (`dti` or `csd`) | `nim_mmf_geometry` |
| `.mmf_built` | logical | — | Geometry present | `nim_mmf_geometry` |
| `.mmf_peakdirs`, `.mmf_kappa_p`, `.mmf_npeaks`, `.mmf_multi` | — | — | Per-peak variants, CSD field only | `nim_mmf_geometry` |

#### Parcellation Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.parcellation_mask` | double | [X, Y, Z] | Integer labels per voxel (0 = background) | `nim_parcellation` |
| `.parcellation_mask_file` | char | — | Path to parcellation mask file | `nim_parcellation` |
| `.labels` | cell array | {N x 1} | Anatomical region names | `nim_load_labels` |
| `.atlas_labels` | containers.Map | — | Index-to-name label mapping | `nim_load_labels` |
| `.parcellation_info` | struct | — | Parcellation metadata and quality metrics | `nim_parcellation_registered` |
| `.roi_masks` | containers.Map | — | Region name -> logical volume, overlaps intact | `nim_parcellation_from_masks` |

`.parcellation_mask` is an integer label volume, so every voxel has exactly one owner. That
representation cannot express genuinely overlapping structures — a temporal-stem voxel belongs
to the uncinate, the inferior longitudinal fasciculus and, near midline, the corpus callosum
at once. `.roi_masks`, where present, holds the regions as separate binary volumes with their
overlaps intact, and `nim_roi_mask` prefers it when resolving ROI names for seeding and
filtering.

#### Registration Fields

| Field | Type | Dimensions | Description | Populated By |
|---|---|---|---|---|
| `.registration` | struct | — | Multi-modal registration data | `nim_registration` |
| `.registration.transforms` | struct | — | Transformation matrices (DTI/T1/MNI) | `register_dti_to_t1`, `register_t1_to_mni` |
| `.registration.quality_metrics` | struct | — | NMI scores and quality assessment | `compute_registration_quality` |

### How the Structure Grows

```
nim_read            → .img, .img_b0, .img_bi, .bval, .bvec, .hdr, .mask, dimensions
    |
nim_dt_spd          → .DT (6-element tensor per voxel)
    |
nim_eig             → .evec (eigenvectors), .eval (eigenvalues)
    |
nim_fa              → .FA (fractional anisotropy map)
    |
nim_mmf_geometry    → .mmf_frames, .mmf_kappa, .mmf_tau, .mmf_built
    |
nim_registration    → .registration (transforms, quality metrics)   [optional]
    |
nim_parcellation    → .parcellation_mask, .labels, .atlas_labels
    |
mask improvement    → .mask (refined using FA)
    |
tissue segmentation → .wm_mask, .gm_mask, .csf_mask
    |
[Saved as .mat]     → Input to tractography
    |
nim_csd             → .peaks, .npeaks, .peak_w   [on demand, field: csd]
```

---

## 4. Module Dependency Graph

Arrows indicate "depends on" (must execute after).

```
nim_read
    |
    +---> nim_preprocessing (optional, if raw data)
    |         |
    |         +---> preproc_extract_b0
    |         +---> preproc_brain_extraction  |  preproc_t1_brain_extraction (if T1)
    |         +---> preproc_denoising (optional)
    |         +---> preproc_fieldmap_correction (optional, if fieldmap available)
    |         +---> preproc_motion_correction (optional)
    |         +---> preproc_eddy_correction (optional)
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
                                  +---> nim_mmf_geometry (requires: .evec or .peaks)
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
              +---> runTractography (resolves the seed mask, then dispatches)
              |         |
              |         +---> nim_csd (if field: csd; cached as <source>_csd.mat)
              |         |
              |         +---> nim_tractography_standard (requires: .evec, .FA, .mask)
              |         |         OR
              |         +---> nim_tractography_hinec  (requires: .evec or .peaks, .FA,
              |         |                              .mask; optional: tissue masks)
              |         |         OR
              |         +---> nim_tractography_mmf_connframe (requires: mmf geometry)
              |         |
              |         +---> nim_filter_tracks_roi (requires: .parcellation_mask or .roi_masks)
              |         +---> nim_resample_track_arc (if output.arc_step > 0)
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

1. Load the image and header with MATLAB's `niftiread` / `niftiinfo`
2. Extract header metadata (dimensions, voxel sizes)
3. Load b-values from the `.bval` file and b-vectors from the `.bvec` file, accepting both the
   one-line and one-value-per-line layouts
4. Split volumes at the b0 threshold (`B0Threshold`, default 5), averaging the b0 volumes into
   `.img_b0`. `.bval` and `.bvec` are kept unfiltered; downstream code filters using
   `.thrsh_b0`
5. Load the brain mask if present (`.nii.gz` with the `_M` suffix)

**Output**: Initialized `nim` struct with image data, acquisition parameters, and metadata.

### Stage 2: Preprocessing (`nim_preprocessing.m`)

Orchestrates a configurable 10-step FSL-based preprocessing pipeline. Automatically detects raw data (filename contains `_raw`) and triggers preprocessing.

**Steps**, in execution order. The optional ones are toggled by config:

| Step | Function | FSL Tool | Purpose |
|---|---|---|---|
| 1 | `preproc_extract_b0` | fslroi | Extract the b0 reference volume |
| 2 | `preproc_brain_extraction` / `preproc_t1_brain_extraction` | BET (+ epi_reg) | Generate the brain mask, from T1 when available |
| 3 | `preproc_denoising` (optional) | dwidenoise / nlmeans / Gaussian | Remove thermal noise |
| 4 | `preproc_fieldmap_correction` (optional) | FUGUE | Correct B0 field inhomogeneity |
| 5 | `preproc_motion_correction` (optional) | mcflirt | Correct head motion, rotating b-vectors to match |
| 6 | `preproc_eddy_correction` (optional) | eddy / eddy_correct | Correct eddy current distortions |
| 7 | — | — | Copy processed data, mask and b-vectors to their final paths |
| 8 | `preproc_t1_dwi_registration`, `preproc_t1_mni_registration` (if T1) | epi_reg, FLIRT+FNIRT, invwarp | Build the MNI→T1→DWI chain |
| 9 | `preproc_atlas_resampling` | applywarp / flirt | Resample the atlas into DWI space |
| 10 | `preproc_cleanup` | — | Finalize outputs, remove temporaries |

**Output**: Preprocessed NIfTI file, brain mask, and parcellation mask in DWI space.

Tissue segmentation is *not* part of this pipeline; it runs from `main.m` (Stage 6).

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

**Supported Atlases**: `jhu` (JHU-ICBM-labels-1mm), `jhu-tract`
(JHU-ICBM-tracts-maxprob-thr0-1mm) and `harvardoxford` (HarvardOxford-cort-maxprob-thr0-1mm),
all read from `$FSLDIR/data/atlases/`. An unrecognised name warns and falls back to
HarvardOxford.

A third pathway exists for datasets that ship their own region definitions:
`nim_parcellation_from_masks` builds a parcellation from a directory of binary mask volumes,
resampling nearest-neighbour through the two sforms so the masks need not be on the DWI grid.
It returns both an integer label volume (overlaps resolved smallest-region-wins, so a large
structure cannot swallow the small specific bundle crossing it) and the untouched per-region
masks, which is what ROI selection should use.

### Stage 6: Tissue Segmentation (`preproc_tissue_segmentation.m`)

Produces the WM/GM/CSF masks that Anatomically Constrained Tractography needs. Called from
`main.m`, not from the preprocessing pipeline.

**Primary**: FSL FAST on the anatomical T1, resampled into DWI space through the two images'
world affines (`flirt -usesqform`). No registration step is involved, so the masks cannot be
corrupted by a T1→DWI registration that is disabled or poorly converged — correct whenever the
T1 is already world-aligned to the diffusion data.

**Fallback** (no usable T1): FA-tertile binning. This bins anisotropy, not tissue, so ACT
driven by it terminates streamlines mid-crossing. It exists for DWI-only data.

The masks are written as `{name}_WM_mask.nii.gz`, `{name}_GM_mask.nii.gz` and
`{name}_CSF_mask.nii.gz`, and loaded onto the `nim`. Failure is non-fatal: tractography
proceeds without ACT.

### Stage 7: Tractography

Two primary algorithms are available:

**Standard FACT** (`nim_tractography_standard.m`):

- Discrete voxel-by-voxel tracking without interpolation
- Variable step size determined by voxel boundary intersection
- Bidirectional tracking from seed points
- Seeding: uniform lattice with optional random jitter

**HINEC interpolated streamlines** (`nim_tractography_hinec.m`):

- Interpolated direction field via pre-computed `griddedInterpolant` objects; the kernel is
  `trilinear` (C0), `cubic` (C1 Keys convolution) or `spline` (C2)
- `interpolation.upsample` refines or coarsens the sampling grid without changing the
  coordinate frame, which is how the spatial axis of a convergence study is swept
- Integration methods: Euler, RK2, RK4, RKF45 (adaptive Dormand-Prince)
- Direction from the DTI principal eigenvector or, with `field: csd`, the FOD peak nearest the
  incoming tangent
- Anatomically Constrained Tractography (ACT) with tissue-type termination, off by default
- Parallel processing support (`parfor` with `DataQueue` progress)

**MMF connection-form** (`nim_tractography_mmf_connframe.m`):

- Advances `dx/ds = e1` while evolving the carried frame by the structure equation, driven by
  the interpolated connection field built in `nim_mmf_geometry`
- `integrator` (`rk4` fixed or `rkf45` adaptive) chooses only the numerical scheme; the
  direction comes entirely from the connection form
- With `field: csd`, the geometry is built per FOD peak, giving multiple pathways

All three trackers share these stopping criteria:

- FA below `termination.fa_min`
- Curvature exceeds `termination.angle_max`, measured in degrees per voxel of arc rather than
  per step. Because the principal direction is a line field, tangents are sign-aligned and a
  measured turn never exceeds 90°, so a budget above `90 / step` is inert rather than loose;
  `0` disables the criterion
- Arc length exceeds `termination.max_arc` voxels (`max_steps` is derived from it)
- Track exits the brain mask
- With ACT on, the track reaches GM (accept) or CSF (reject)

Tracks shorter than `termination.min_arc` voxels of arc are discarded.

**Seeding strategy**, resolved in `runTractography` (4-tier priority):

1. Explicit ROI (`seeding.roi`) — named atlas regions, optionally dilated
2. Preprocessed brain mask
3. Expanded parcellation mask (dilated 3 voxels)
4. FA-threshold fallback (FA > 0.10)

Every tier then intersects with `seeding.fa_min`. Within a seeded voxel, `seeding.density`
seeds are placed on a deterministic sub-voxel lattice (`uniform`) or jittered (`random`).

**Output**: cell array of tracks (each an Nx3 matrix of voxel-space coordinates spanning the
complete trajectory), plus the flat options struct, the algorithm name, elapsed time and
per-track metadata.

### Stage 8: Visualization

Three tiers of visualization capability:

1. **MATLAB 3D** (`visualizeTractography.m`): Unified viewer with whole-brain, region, grid, and sequential modes. Supports direction/FA/region coloring and export.

2. **MATLAB 2D Slices** (`visualizeTractographySlices.m`): interactive orthogonal slice viewer with crosshair synchronization. Each slice update re-renders from the track data, so navigation is slow on large tractograms.

3. **Fast Distributed Viewer**: pre-generate a slice cache on the machine holding the data (`generateSlices.m`), transfer it, and browse locally with the Python GUI (`FastTractographyViewer.py`). Navigation reads pre-rendered images, so no MATLAB is needed locally.

See [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md) for complete details.

---

## 6. Configuration System

HINEC uses YAML configuration files for reproducible parameter management.

### Configuration Structure

Configuration nests exactly two levels below a section (`section` → `group` → `key`); a third
level is a parse error. Every key is optional and falls back to the schema default, so a
working config can be three lines.

```yaml
preprocessing:
  run_denoising: true|false
  denoise_method: dwidenoise | nlmeans | gaussian
  run_motion_correction: true|false
  run_eddy: true|false
  improve_mask: true|false
  atlas_type: jhu | harvardoxford | jhu-tract
  t1_available: true|false
  use_t1_registration: true|false
  register_to_mni: true|false

tractography:
  algorithm: hinec | standard | mmf       # which tracker
  field: dti | csd                        # where the local direction comes from
  act: true|false                         # default FALSE
  integrator:
    method: euler | rk2 | rk4 | rkf45     # a METHOD name, not an order
    step: 0.2                             # voxel units; initial step for rkf45
    tolerance, step_min, step_max, safety, adaptive   # rkf45 only
  interpolation:
    method: trilinear | cubic | spline    # C0 | C1 | C2
    upsample: 1                           # >1 refines, <1 coarsens
  seeding:
    density, strategy, fa_min, roi, roi_dilate
  termination:
    fa_min: 0.10
    angle_max: 225                        # DEGREES PER VOXEL OF ARC; 0 disables
    max_arc: 400                          # voxels; max_steps derived from it
    min_arc: 15                           # voxels
  filter:
    include_roi, exclude_roi, mode, roi_dilate, endpoints_in, contained_in
  output:
    arc_step: 0                           # resample saved streamlines
  csd:
    lmax, max_peaks, peak_thresh, peak_min_sep, n_iter
  mmf:
    anchor, frame_sel_power
  diagnostics: true|false
```

`src/nim_utils/nim_config_schema.m` declares every parameter exactly once — its canonical
path, default, type, permitted values, legacy aliases, the algorithms that actually read it,
and a one-line description. Everything else derives from it: `load_config_yaml` takes its
defaults, validation, unknown-key rejection and legacy migration from the schema;
`nim_config_to_options` flattens to the legacy option names the trackers still read; and
[YAML_CONFIG.md](YAML_CONFIG.md) is generated by `nim_config_docs`, so the reference cannot
drift from the code.

Two categories of legacy key are handled differently. **Migrated** keys are transformed:
`integration_order` (1/2/4/5) becomes `integrator.method`, and `max_steps` becomes
`max_arc = max_steps × step`. **Retired** keys are warned about and dropped rather than
remapped, because remapping them would change behaviour — `fa_threshold` was functionally dead
and is not equivalent to either `termination.fa_min` or `seeding.fa_min`, and `sel_power` has
been removed from HINEC entirely.

Any parameter can be overridden on the command line with `--set key=value`, resolved against
the schema by canonical path, by `group.key`, or by a bare leaf name when unique.

### Available Presets

Configs come in two families, and the family determines which launcher consumes them.
**Tracker** configs are named `<algorithm>_<field>[_<variant>].yml` and drive
`bin/run_tractography.sh`; **dataset** configs are named `<dataset>[_variant].yml` and drive
`bin/run_hinec.sh`.

| Preset | Family | Algorithm | Field | Integrator | Interpolation | Step |
|---|---|---|---|---|---|---|
| `standard_dti.yml` | tracker | standard | dti | euler | — | 0.5 |
| `hinec_dti.yml` | tracker | hinec | dti | rkf45 | trilinear | 0.2 |
| `hinec_csd.yml` | tracker | hinec | csd | rkf45 | trilinear | 0.2 |
| `hinec_dti_cubic.yml` | tracker | hinec | dti | rkf45 | cubic | 0.2 |
| `hinec_dti_cubic_recall.yml` | tracker | hinec | dti | rkf45 | cubic | 0.2 |
| `hinec_dti_fast.yml` | tracker | hinec | dti | rk2 | trilinear | 0.3 |
| `hinec_dti_euler.yml` | tracker | hinec | dti | euler | trilinear | 0.5 |
| `mmf_dti.yml` | tracker | mmf | dti | rk4 | trilinear | 0.2 |
| `mmf_csd.yml` | tracker | mmf | csd | rk4 | trilinear | 0.2 |
| `ismrm2015.yml` | dataset | hinec | dti | rkf45 | trilinear | 0.2 |
| `irontract.yml` | dataset | hinec | dti | rkf45 | trilinear | 0.2 |
| `hinec_default.yml` | fallback | hinec | dti | rk4 | cubic | 0.2 |
| `reference.yml` | template | — | — | fixed step, set on the CLI | — | — |

Presets carry only their non-default values; everything else comes from the schema. There are
three orthogonal axes — `algorithm` (which tracker), `field` (`dti` or `csd`) and the
integrator — and, for hinec, direction is further shaped by two composable factors:
`interpolation.method` (spatial kernel) and `interpolation.upsample` (sampling density).
Settings irrelevant to a chosen combination are ignored by design.

See [YAML_CONFIG.md](YAML_CONFIG.md) for the complete parameter reference and
`config/README.md` for the naming rules.

---

## 7. Run Directory System

Each pipeline execution creates a timestamped directory for reproducibility:

```
hinec_runs/
└── run_YYYYMMDD_HHMMSS_<config>/
    ├── run_info.txt              # Metadata (git hash, timestamp, description)
    ├── config.yml                # Copied configuration file
    ├── overrides.txt             # --set overrides, when any were given
    ├── logs/pipeline.log         # Live execution log
    ├── intermediate/             # Preprocessing artifacts
    ├── output/                   # Copy of the processed nim
    └── tractography/             # Track results
        └── diagnostics/          # Quality reports
```

A `latest` symlink points to the most recent run (a text file containing the path, on
Windows).

Note the division of layers: the canonical processed `nim`, the preprocessed references and
the CSD cache live in the **data directory** next to the input; run directories hold
tractography outputs. See [RUN_DIRECTORY_SYSTEM.md](RUN_DIRECTORY_SYSTEM.md) for details.

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

Expected at `lib/spm12/`, which `main.m` adds with `addpath(genpath('lib/spm12'))`. It is
**not** vendored in this repository — install it there. Used for:

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

Not needed for the MATLAB pipeline itself. Required for the fast distributed slice viewer
(`scripts/FastTractographyViewer.py`, `scripts/tractography_slice_gui.py`) and for the TRK
export and validation tooling (`scripts/hinec_to_trk.py`,
`scripts/validate_ismrm_tractography.py`). `requirements.txt` covers all of them; the viewer
additionally needs tkinter, which normally ships with Python.

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
| `{name}_T1.nii.gz` | T1 anatomical (optional) | `sample_T1.nii.gz` |
| `{name}_fmap_Hz.nii.gz` | Field map in Hz (optional) | `sample_fmap_Hz.nii.gz` |

### Output Files

| Pattern | Description | Example |
|---|---|---|
| `{name}.mat` | Processed nim structure | `ismrm2015.mat` |
| `{name}_csd.mat` | Cached CSD FOD peaks | `ismrm2015_csd.mat` |
| `{name}_WM_mask.nii.gz` and friends | Tissue masks for ACT | `sample_GM_mask.nii.gz` |
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
- Tractography methods overview: [TRACTOGRAPHY_METHODS.md](TRACTOGRAPHY_METHODS.md)
- Standard FACT: [TRACTOGRAPHY.md](TRACTOGRAPHY.md)
- MMF connection-form tractography: [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md)
- Scientific foundations: [SCIENCE.md](SCIENCE.md)
- Mathematical methods: [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md)
- Complete function reference: [API_REFERENCE.md](API_REFERENCE.md)
- Visualization guide: [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)
