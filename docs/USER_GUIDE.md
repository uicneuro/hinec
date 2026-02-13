# Getting Started and Usage Guide

Practical guide for installing, configuring, running the HINEC pipeline, and interpreting results.

---

## 1. Prerequisites

### MATLAB

- **Version**: R2020b or later (for `arguments` block syntax and `griddedInterpolant`)
- **Required Toolboxes**:
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Tools for NIfTI and ANALYZE image (included in `lib/spm12/`)

Verify toolboxes:
```matlab
ver('images')    % Image Processing Toolbox
ver('stats')     % Statistics and Machine Learning Toolbox
```

### FSL (FMRIB Software Library)

Required for preprocessing. Not needed if you provide already-preprocessed data.

1. Install FSL from [fsl.fmrib.ox.ac.uk](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)
2. Initialize FSL before running HINEC:
```bash
# Add to your shell profile (~/.bashrc or ~/.zshrc)
export FSLDIR=/usr/local/fsl
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH
```
3. Verify: `flirt -version` should print the FSL version

### SPM12

Included in the repository at `lib/spm12/`. No separate installation needed — `main.m` automatically adds it to the MATLAB path.

### Python (Optional)

Required only for the fast distributed slice viewer:
- Python 3.7+
- `pip install Pillow numpy`
- tkinter (usually included with Python)

---

## 2. Installation

```bash
# Clone the repository
git clone <repository-url> hinec
cd hinec

# Verify sample data exists
ls data/original_sample/
ls data/parcellation_sample/
```

No build step is required. MATLAB source files are used directly.

---

## 3. Quick Start

Process the sample data end-to-end in under 5 minutes:

```matlab
% In MATLAB, navigate to the hinec directory
cd /path/to/hinec

% Process sample data (preprocessed, no FSL needed)
main('data/original_sample/sample', 'output/sample.mat');

% Run tractography
runTractography('output/sample.mat');

% Load and view results
load('output/sample.mat');
nim_plot(nim, 'mode', 'single');
```

Or use the pre-computed sample:

```matlab
% Load pre-computed results
load('data/parcellation_sample/sample_motion_corrected.mat');

% Visualize parcellations
nim_plot(nim, 'mode', 'parcels');
```

### Expected Output

After `main()` completes, the output .mat file contains a `nim` struct with:
- `.img` — Original image data
- `.DT` — Diffusion tensors (6 components per voxel)
- `.evec` — Eigenvectors (fiber directions)
- `.eval` — Eigenvalues (diffusion magnitudes)
- `.FA` — Fractional anisotropy map
- `.parcellation_mask` — Brain region labels
- `.labels` — Region names

After `runTractography()`, a tracks file appears in `tractography_results/` containing:
- `tracks` — Cell array of fiber pathways (each track is an Nx3 matrix)
- `options` — Parameters used
- `elapsed_time` — Processing time

---

## 4. Running the Full Pipeline

### Option A: MATLAB Direct

```matlab
% Basic usage (preprocessed data)
main('path/to/data/subject', 'output/subject.mat');

% With T1 for registration-enhanced parcellation
main('path/to/data/subject_raw', 'output/subject.mat', ...
     't1_file', 'path/to/data/subject_t1.nii.gz');

% With YAML configuration
config = load_config_yaml('config/high_precision.yml');
main('path/to/data/subject', 'output/subject.mat', 'config', config);

% With run directory organization
config = load_config_yaml('config/hinec_default.yml');
run_info = create_run_directory('config/hinec_default.yml', ...
    'description', 'First processing run');
main('path/to/data/subject', 'output/subject.mat', ...
     'config', config, 'run_info', run_info);
```

### Option B: Shell Script

```bash
# Full pipeline with configuration
./bin/run_hinec.sh path/to/data/subject.nii.gz config/hinec_default.yml

# The script:
# 1. Validates input files exist
# 2. Creates a timestamped run directory
# 3. Runs preprocessing (Phase 1)
# 4. Runs tractography (Phase 2)
# 5. Logs everything to the run directory
```

### Input Data Requirements

| File | Required | Format | Description |
|---|---|---|---|
| `{name}.nii.gz` | Yes | NIfTI-1 | 4D diffusion-weighted images |
| `{name}.bval` | Yes | Text | B-values (space-separated, one line) |
| `{name}.bvec` | Yes | Text | Gradient directions (3 rows × N columns) |
| `{name}_M.nii.gz` | No | NIfTI-1 | Brain mask (binary) |
| `{name}_raw.nii.gz` | Alt | NIfTI-1 | Raw data (triggers preprocessing) |
| `{name}_acqp.txt` | No | Text | Acquisition parameters (for eddy correction) |
| `{name}_index.txt` | No | Text | Volume index (for eddy correction) |

**Raw vs preprocessed**: If the filename contains `_raw`, HINEC automatically runs FSL preprocessing. Otherwise, it assumes data is already preprocessed.

---

## 5. Configuration

### Using YAML Presets

HINEC provides 5 configuration presets in `config/`:

| Preset | File | Algorithm | Integration | Best For |
|---|---|---|---|---|
| Default | `hinec_default.yml` | hinec | RK4 | General use |
| High Precision | `high_precision.yml` | hinec | RKF45 adaptive | Publication figures |
| Fast Exploration | `fast_exploration.yml` | hinec | RK2 | Quick parameter testing |
| Standard FACT | `standard_fact.yml` | standard | Euler | Baseline comparison |
| IronTract | `irontract.yml` | hinec | RK4 | IronTract challenge |

```matlab
% Load a preset
config = load_config_yaml('config/high_precision.yml');

% Use it
main('data/subject', 'output/subject.mat', 'config', config);
runTractography('output/subject.mat', 'config', config);
```

### Key Tractography Parameters

| Parameter | Range | Default | Impact |
|---|---|---|---|
| `step_size` | 0.1 – 1.0 | 0.2 (HINEC), 0.5 (FACT) | Smaller = smoother, slower |
| `fa_threshold` | 0.05 – 0.3 | 0.1 | Lower = more tracks, more noise |
| `termination_fa` | 0.01 – 0.2 | 0.05 | Lower = longer tracks |
| `angle_threshold` | 20 – 90 | 60 | Higher = allows sharper turns |
| `seed_density` | 1 – 10 | 4 | Higher = denser coverage, slower |
| `integration_order` | 1, 2, 4, 5 | 4 | Higher = more accurate |

### Creating Custom Configurations

Copy an existing preset and modify:

```bash
cp config/hinec_default.yml config/my_config.yml
# Edit my_config.yml with your parameters
```

See [YAML_CONFIG.md](YAML_CONFIG.md) for complete parameter documentation.

---

## 6. Running Tractography

### Standard FACT

```matlab
% Basic FACT tractography
runTractography('output/subject.mat');

% Or explicitly
runTractography('output/subject.mat', 'algorithm', 'standard');
```

### HINEC High-Order

```matlab
% RK4 with trilinear interpolation (default HINEC)
runTractography('output/subject.mat', 'algorithm', 'hinec');

% RKF45 adaptive stepping (highest accuracy)
config = load_config_yaml('config/high_precision.yml');
runTractography('output/subject.mat', 'algorithm', 'hinec', 'config', config);
```

### Direct Function Calls (Advanced)

For more control, call tractography functions directly:

```matlab
% Load processed data
data = load('output/subject.mat');
nim = data.nim;

% Configure options
options = struct();
options.step_size = 0.2;
options.integration_order = 4;
options.fa_threshold = 0.15;
options.angle_thresh = 45;
options.seed_density = 3;
options.act_enabled = true;
options.wm_mask = nim.wm_mask;
options.gm_mask = nim.gm_mask;
options.csf_mask = nim.csf_mask;

% Run
tracks = nim_tractography_hinec(nim, options);
```

### Understanding Track Output

Each track is an Nx3 matrix representing a complete fiber pathway:

```matlab
% Load tracks
data = load('tractography_results/tracks_hinec_2025-01-15_14_30_00.mat');
tracks = data.tracks;

fprintf('Total tracks: %d\n', length(tracks));

% Examine a single track
track = tracks{1};
fprintf('Points in track: %d\n', size(track, 1));
fprintf('Start: [%.1f, %.1f, %.1f]\n', track(1,:));
fprintf('End:   [%.1f, %.1f, %.1f]\n', track(end,:));

% Track statistics
lengths = cellfun(@(t) size(t, 1), tracks);
fprintf('Mean length: %.1f points\n', mean(lengths));
fprintf('Max length:  %d points\n', max(lengths));
```

---

## 7. Visualization

### Quick 3D View

```matlab
% Load data
load('output/subject.mat');

% Eigenvector visualization
nim_plot(nim, 'mode', 'single');        % Single region
nim_plot(nim, 'mode', 'parcel', 'region_id', 5);  % Specific region
nim_plot(nim, 'mode', 'parcels');       % All regions (separate figures)
```

### Tractography Visualization

```matlab
% Whole brain
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'whole');

% Single region
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'region', 'regions', [5]);

% All regions in grid
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'grid');

% Export for publication
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'whole', ...
    'export', 'figures/whole_brain.png');
```

### Interactive Slice Viewer

```matlab
visualizeTractographySlices('tracks.mat', 'nim.mat');
```

### Fast Distributed Viewer

For instant navigation (sub-100ms), pre-generate slices:

```matlab
% On server (MATLAB): generate cache
generateSlices('tracks.mat', 'nim.mat', '/path/to/cache');
```

```bash
# Transfer cache to local machine
rsync -avz server:/path/to/cache/ ~/local/cache/

# View locally (Python only, no MATLAB needed)
./bin/viewSlices.sh ~/local/cache/
```

See [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md) for complete visualization reference and [DISTRIBUTED_WORKFLOW.md](DISTRIBUTED_WORKFLOW.md) for distributed setup.

---

## 8. Working with Results

### Loading Saved Data

```matlab
% Load processed nim structure
data = load('output/subject.mat');
nim = data.nim;

% Access fields
fa_map = nim.FA;              % 3D fractional anisotropy
evecs = nim.evec;             % Eigenvectors
evals = nim.eval;             % Eigenvalues
mask = nim.parcellation_mask; % Brain region labels
labels = nim.labels;          % Region names
```

### Extracting Connectivity Matrices

```matlab
% Load tracks and nim
tracks_data = load('tractography_results/tracks_hinec.mat');
nim_data = load('output/subject.mat');

% Compute connectivity
C = nim_plot_connectivity_matrix(tracks_data.tracks, nim_data.nim);

% Access the matrix
fprintf('Regions: %d x %d\n', size(C));
fprintf('Connections: %d\n', nnz(C));
```

### Comparing Runs

```matlab
% Load two runs
run1 = load('hinec_runs/run_2025-01-15/output/subject.mat');
run2 = load('hinec_runs/run_2025-01-16/output/subject.mat');

% Compare FA maps
fa_diff = run2.nim.FA - run1.nim.FA;
fprintf('Mean FA difference: %.4f\n', mean(fa_diff(:), 'omitnan'));

% Compare track counts
t1 = load('hinec_runs/run_2025-01-15/tractography/tracks.mat');
t2 = load('hinec_runs/run_2025-01-16/tractography/tracks.mat');
fprintf('Run 1 tracks: %d, Run 2 tracks: %d\n', length(t1.tracks), length(t2.tracks));
```

---

## 9. Troubleshooting

### FSL Not Found

**Symptom**: `'flirt' is not recognized` or `command not found`

**Solution**:
```bash
# Check FSL is installed
echo $FSLDIR
which flirt

# If not set, add to shell profile:
export FSLDIR=/usr/local/fsl
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH
```

In MATLAB, ensure FSL is accessible:
```matlab
[status, result] = system('flirt -version');
if status ~= 0
    error('FSL not found. Set FSLDIR and source fsl.sh before launching MATLAB.');
end
```

### Out of Memory

**Symptom**: MATLAB crashes or reports `Out of memory`

**Solutions**:
1. Process smaller regions: use parcellation to focus on specific areas
2. Reduce seed density: `options.seed_density = 2` instead of 4
3. Use v7.3 MAT files: `nim_save` automatically uses v7.3 for large files
4. Close other applications to free RAM
5. Use the `fast_exploration.yml` preset for initial testing

### String vs Char Array Errors

**Symptom**: `Error using fullfile` or path concatenation errors

**Cause**: MATLAB string (`"text"`) vs char array (`'text'`) incompatibility.

**Solution**: Convert strings to char arrays when passing to file operations:
```matlab
path = char(string_path);  % Convert if needed
```

See [STRING_CHAR_FIXES.md](STRING_CHAR_FIXES.md) for known issues and fixes.

### No Tracks Generated

**Symptom**: `runTractography` produces 0 tracks

**Possible causes**:
1. **FA threshold too high**: Try lowering `fa_threshold` to 0.1 or 0.05
2. **No valid seed mask**: Check that brain mask or parcellation mask exists
3. **Data quality**: Verify FA map has reasonable values (max > 0.3)
4. **Preprocessing issues**: Check that eigenvectors are computed correctly

**Diagnostic**:
```matlab
data = load('output/subject.mat');
fprintf('Max FA: %.3f\n', max(data.nim.FA(:)));
fprintf('Mask voxels: %d\n', sum(data.nim.mask(:)));
fprintf('Has eigenvectors: %d\n', isfield(data.nim, 'evec'));
```

### Preprocessing Fails

**Symptom**: Error during `nim_preprocessing`

**Common issues**:
1. **Missing FSL tools**: Ensure all required FSL tools are installed (mcflirt, bet, flirt, fnirt)
2. **Missing input files**: Verify `.bval` and `.bvec` files exist alongside the NIfTI
3. **Permission errors**: Check write permissions in the output directory
4. **Incompatible data format**: Ensure NIfTI files are valid (try `fslinfo input.nii.gz`)

### Registration Quality Issues

**Symptom**: Parcellation regions appear misaligned

**Solutions**:
1. Check registration quality metrics (NMI scores)
2. Provide a T1 image for better registration accuracy
3. Verify T1 and DWI are from the same subject/session
4. Try different registration method (FSL vs SPM)

---

## Cross-References

- Architecture overview: [ARCHITECTURE.md](ARCHITECTURE.md)
- Complete API reference: [API_REFERENCE.md](API_REFERENCE.md)
- Configuration system: [YAML_CONFIG.md](YAML_CONFIG.md)
- Visualization guide: [VISUALIZATION_GUIDE.md](VISUALIZATION_GUIDE.md)
- Preprocessing details: [PREPROCESSING.md](PREPROCESSING.md)
- IronTract workflow: [IRONTRACT_WORKFLOW.md](IRONTRACT_WORKFLOW.md)
- Run directory system: [RUN_DIRECTORY_SYSTEM.md](RUN_DIRECTORY_SYSTEM.md)
