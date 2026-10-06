# HINEC: HIgh-order NEural Connectivity

Human brain white matter tractography pipeline with **YAML-based parameter configuration**.

## Quick Start

### Using Shell Script (Recommended)

```bash
# Usage: ./bin/run_hinec.sh <data_prefix> <output_mat> [config_file]
#
# <data_prefix> is the shared path of your input files, without extensions.
# Given files: my_data/subject_raw.nii.gz, my_data/subject.bval, my_data/subject.bvec
# The prefix is: my_data/subject

# Default configuration
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat

# High precision (publication quality)
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/hinec_dti_cubic.yml

# Fast exploration (parameter testing)
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/hinec_dti_fast.yml

# Export visualization figures (after pipeline completes)
./bin/run_visualization.sh hinec_runs/run_20260330_*/ figures/my_output
```

**NEW**: All runs are now automatically organized in `hinec_runs/run_YYYYMMDD_HHMMSS_<config>/` directories! No more scattered files cluttering your workspace.

### Using MATLAB

```matlab
% Load configuration
config = load_config_yaml('config/hinec_default.yml');

% Create organized run directory
run_info = create_run_directory('config/hinec_default.yml');

% Run pipeline with organized output
main('data/parcellation_sample/sample', 'sample.mat', config, run_info);
output_mat = fullfile(run_info.output_dir, 'sample.mat');
runTractography(output_mat, config, run_info);
```

## NEW Features

### 🗂️ Organized Run Directories

Every pipeline run creates a timestamped directory with all outputs:

```
hinec_runs/run_20250118_143045_hinec_dti_cubic/
├── config.yml          # Copy of config used
├── run_info.txt        # Run metadata
├── logs/               # Pipeline logs
├── intermediate/       # Preprocessing artifacts
├── output/             # Final .mat file
└── tractography/       # Tracks + diagnostics
```

**Benefits**:
- ✅ Clean workspace - no scattered files
- ✅ Full reproducibility - config and metadata saved
- ✅ Easy comparison between runs
- ✅ Simple cleanup - delete entire run directory

**Documentation**: See [docs/RUN_DIRECTORY_SYSTEM.md](docs/RUN_DIRECTORY_SYSTEM.md)

### ⚙️ YAML Configuration System

**Easy parameter management** - no more editing source code!

**Benefits**:
- Human-readable configuration files
- Version control your experiment parameters
- Switch between presets instantly
- Automatic parameter validation

**Available Presets**:
- `hinec_default.yml` - Balanced performance (recommended)
- `hinec_dti_cubic.yml` - Publication quality (RKF45 adaptive)
- `hinec_dti_fast.yml` - Quick parameter testing (3-5x faster)
- `hinec_dti_euler.yml` - Euler with trilinear interpolation
- `irontract.yml` - IronTract challenge optimized

**Documentation**: See [docs/YAML_CONFIG.md](docs/YAML_CONFIG.md)

---

## Requirements

Addons:

-   `Image Processing Toolbox`
-   `Statistics and Machine Learning Toolbox`
-   `Tools for NIfTI and ANALYZE image`

External Softwares:

-   `Statistical Parametric Mapping`
    -   Must add folder SPM12 in the `lib` directory
-   `FSL`
    -   Must be initialized before use
