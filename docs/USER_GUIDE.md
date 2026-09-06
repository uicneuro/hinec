# Getting Started and Usage Guide

Practical guide for installing, configuring, running the HINEC pipeline, and interpreting results.

---

## 1. Prerequisites

### MATLAB

- **Version**: R2020b or later (for `arguments` block syntax and `griddedInterpolant`)
- **Required Toolboxes**:
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Tools for NIfTI and ANALYZE image (MATLAB File Exchange)

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

Used for NIfTI I/O and as an optional registration backend. SPM12 is **not** vendored in this
repository: install it at `lib/spm12/`, which is the path `main.m` adds with
`addpath(genpath('lib/spm12'))`.

### Python (Optional)

Not needed for the MATLAB pipeline. Required for the fast distributed slice viewer, for TRK
export, and for the ISMRM scoring tooling:

- Python 3.7+
- `pip install -r requirements.txt`
- tkinter for the viewer's GUI (usually included with Python)

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

### Using Shell Scripts (Recommended)

The simplest way to run HINEC end-to-end — preprocessing, DTI, and tractography in one command:

```bash
cd /path/to/hinec

# Process sample data with default configuration
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat

# Export visualization figures (after pipeline completes)
./bin/run_visualization.sh hinec_runs/latest/
```

The first argument is a **data prefix** — this is the shared path and name that all your input files have in common, without any suffix or extension. For example, given these files:

```
data/parcellation_sample/
├── sample_raw.nii.gz    ← raw diffusion data
├── sample.bval          ← b-values
├── sample.bvec          ← b-vectors
└── sample_T1.nii.gz     ← T1 anatomical (optional; not shipped with this sample)
```

The data prefix is `data/parcellation_sample/sample` — HINEC appends `_raw.nii.gz`, `.bval`, `.bvec`, etc. automatically.

More examples:

| Your files | Data prefix |
|---|---|
| `data/ISMRM/ISMRM_raw.nii.gz`, `ISMRM.bval`, `ISMRM.bvec` | `data/ISMRM/ISMRM` |
| `subjects/sub01/dwi_raw.nii.gz`, `dwi.bval`, `dwi.bvec` | `subjects/sub01/dwi` |
| `my_data.nii.gz`, `my_data.bval`, `my_data.bvec` | `my_data` |

If your files end in `_raw.nii.gz`, HINEC runs FSL preprocessing automatically. If they end in just `.nii.gz` (no `_raw`), it assumes data is already preprocessed.

### Using MATLAB

```matlab
cd /path/to/hinec

% Process sample data (preprocessed, no FSL needed)
main('data/original_sample/sample', 'output/sample.mat');

% Run tractography
runTractography('output/sample.mat');

% Load and view results
load('output/sample.mat');
nim_plot(nim, 'mode', 'single');
```

### Expected Output

After `main()` completes, the output `.mat` file contains a `nim` struct with:

- `.img` — Original image data
- `.DT` — Diffusion tensors (6 components per voxel)
- `.evec` — Eigenvectors (fiber directions)
- `.eval` — Eigenvalues (diffusion magnitudes)
- `.FA` — Fractional anisotropy map
- `.wm_mask`, `.gm_mask`, `.csf_mask` — Tissue masks for ACT
- `.parcellation_mask` — Brain region labels
- `.labels` — Region names

After `runTractography()`, a tracks file is written to the run directory's `tractography/`
(or to `tractography_results/` when no run directory is in use), containing:

- `tracks` — Cell array of fiber pathways (each track is an Nx3 matrix)
- `options` — Flat option struct actually passed to the tracker
- `algorithm` — Tracker that produced the tracks
- `elapsed_time` — Processing time
- `track_meta` — Per-track metadata

---

## 4. Running the Full Pipeline

### Option A: MATLAB Direct

`main` and `runTractography` take **positional** optional arguments, not name-value pairs;
the argument's type selects its meaning (see the parsing block at the top of each file).

```matlab
% Basic usage (preprocessed data)
main('path/to/data/subject', 'output/subject.mat');

% With T1 for registration-enhanced parcellation (3rd argument = T1 path)
main('path/to/data/subject_raw', 'output/subject.mat', ...
     'path/to/data/subject_T1.nii.gz');

% With YAML configuration (3rd argument = config struct)
config = load_config_yaml('config/hinec_dti_cubic.yml');
main('path/to/data/subject', 'output/subject.mat', config);

% With run directory organization (4th argument = run_info struct)
config   = load_config_yaml('config/hinec_default.yml');
run_info = create_run_directory('config/hinec_default.yml', ...
    'description', 'First processing run');
main('path/to/data/subject', 'output/subject.mat', config, run_info);
runTractography(fullfile(run_info.output_dir, 'subject.mat'), config, run_info);
```

### Option B: Shell Script

```bash
# Usage: ./bin/run_hinec.sh <data_prefix> <output_mat> [config_file] [--set key=value ...]

# Default config (cubic interpolation + RK4 + ACT)
./bin/run_hinec.sh path/to/data/subject subject.mat

# Standard FACT for comparison
./bin/run_hinec.sh path/to/data/subject subject.mat config/standard_dti.yml

# Cubic interpolation with adaptive RKF45
./bin/run_hinec.sh path/to/data/subject subject.mat config/hinec_dti_cubic.yml

# Override any config parameter on the command line
./bin/run_hinec.sh path/to/data/subject subject.mat config/ismrm2015.yml \
    --set seeding.density=8 --set integrator.method=rk4
```

The script validates the input files (`<prefix>_raw.nii.gz`, `.bval`, `.bvec`), creates a
timestamped run directory under `hinec_runs/`, then runs preprocessing plus DTI calculation
followed by tractography in a single background `matlab -batch` process, logging to
`<run_dir>/logs/pipeline.log`.

`<data_prefix>` is the shared path without `_raw.nii.gz` — see [Quick Start](#3-quick-start)
for examples.

### Option C: Preprocess Once, Iterate Tractography

Preprocessing is the expensive stage, and it does not change when you vary tracking
parameters. `bin/run_tractography.sh` runs tractography **only**, against an
already-preprocessed `nim`:

```bash
# 1. Preprocess once. main.m writes the canonical nim next to the input data.
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml

# 2. Run any number of tractography configs against it — no re-preprocessing.
./bin/run_tractography.sh hinec_csd --score

# 3. Sweep a parameter without writing new YAML files.
for d in 2 4 8; do
    ./bin/run_tractography.sh hinec_dti --set seeding.density=$d
done
```

`<config>` accepts `config/<name>.yml` or the bare `<name>`. Each invocation produces its own
run directory, tagged with any `--set` overrides and recording them in `overrides.txt`.
`--source <nim|dir>` selects a different preprocessed source; `--score` runs
`bin/run_ismrm_scoring.sh` on the finished run directory.

Keep the processed `nim`, the preprocessed references and the CSD cache in the data directory;
run directories hold tractography **outputs** only.

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
| `{name}_T1.nii.gz` | No | NIfTI-1 | T1 anatomical; enables T1-based brain extraction, registration and FAST tissue segmentation |
| `{name}_fmap_Hz.nii.gz` | No | NIfTI-1 | Field map in Hz, for susceptibility distortion correction |

**Raw vs preprocessed**: if the filename contains `_raw`, HINEC runs FSL preprocessing.
Otherwise it assumes the data is already preprocessed. `main.m` short-circuits aggressively:
an existing output `.mat` skips everything, and an existing preprocessed `.nii.gz` skips
preprocessing. To force reprocessing, delete the target file.

---

## 5. Configuration

### Using YAML Presets

Configs in `config/` fall into two families, distinguished by which launcher consumes them.

**Tracker configs** — `<algorithm>_<field>[_<variant>].yml`, for `bin/run_tractography.sh`:

| File | Algorithm | Field | Integrator | Interpolation |
|---|---|---|---|---|
| `standard_dti.yml` | standard (FACT) | dti | euler | none (discrete voxel) |
| `hinec_dti.yml` | hinec | dti | rkf45 | trilinear |
| `hinec_csd.yml` | hinec | csd | rkf45 | trilinear |
| `hinec_dti_cubic.yml` | hinec | dti | rkf45 | cubic |
| `hinec_dti_cubic_recall.yml` | hinec | dti | rkf45 | cubic (recall-tuned termination) |
| `hinec_dti_fast.yml` | hinec | dti | rk2 | trilinear |
| `hinec_dti_euler.yml` | hinec | dti | euler | trilinear |
| `mmf_dti.yml` | mmf | dti | rk4 | trilinear |
| `mmf_csd.yml` | mmf | csd | rk4 | trilinear |

**Dataset configs** — `<dataset>[_variant].yml`, for `bin/run_hinec.sh`: `ismrm2015.yml` and
`irontract.yml`. `hinec_default.yml` is the generic full-pipeline fallback and is the only
config pairing a full preprocessing block with a general-purpose tractography block.
`reference.yml` is a template for step-size convergence ladders, driven entirely from
`--set`.

```matlab
% Load a preset and pass it positionally
config = load_config_yaml('config/hinec_dti_cubic.yml');
main('data/subject', 'output/subject.mat', config);
runTractography('output/subject.mat', config);
```

### Key Tractography Parameters

Configuration nests exactly two levels below a section (`section` → `group` → `key`); a third
level is a parse error. Defaults come from the schema, so a working config is usually short.

| Parameter | Range | Default | Effect |
|---|---|---|---|
| `integrator.method` | euler, rk2, rk4, rkf45 | rk4 | Stepping scheme — a method name, not an order |
| `integrator.step` | 1e-6 – 10 | 0.2 | Step in voxels; the initial step for rkf45 |
| `interpolation.method` | trilinear, cubic, spline | trilinear | Kernel smoothness: C0, C1 (Keys convolution), C2 (true spline) |
| `interpolation.upsample` | 0.05 – 8 | 1 | Field sampling factor: >1 refines, <1 coarsens |
| `seeding.density` | 1 – 1000 | 8 | Seeds per seeded voxel, honoured exactly |
| `seeding.fa_min` | 0 – 1 | 0.05 | Minimum FA for a voxel to be seeded |
| `termination.fa_min` | 0 – 1 | 0.10 | Stop tracking below this FA |
| `termination.angle_max` | 0 – 3600 | 225 | Maximum turn in **degrees per voxel of arc**; `0` disables |
| `termination.max_arc` | 1 – 100000 | 400 | Maximum arc length in voxels |
| `termination.min_arc` | 0 – 100000 | 15 | Discard tracks shorter than this arc length |
| `act` | true/false | **false** | Anatomically constrained tracking |

Two of these behave in ways worth stating explicitly.

`termination.angle_max` is a curvature budget, not a per-step limit: it fixes a minimum radius
of curvature `R = 57.3 / angle_max` voxels, so refining the step does not loosen the
constraint. Because the principal direction is a line field, tangents are sign-aligned and a
measured turn never exceeds 90 degrees — any budget above `90 / step` is **inert**, not merely
loose. Set `angle_max: 0` when you want the criterion off deliberately.

`act` defaults to false. `main.m` always generates the WM/GM/CSF tissue masks and stores them
on the `nim`, so ACT being off is a configuration choice rather than a missing input;
`runTractography` prints which of the two reasons applies for a given run.

### Creating Custom Configurations

Copy an existing preset and modify:

```bash
cp config/hinec_default.yml config/my_config.yml
# Edit my_config.yml with your parameters
```

For a one-off sweep, prefer `--set key=value` over a new file. Any parameter is reachable by
its canonical path (`tractography.integrator.step`), by its group and key
(`integrator.step` — the `tractography.` prefix is assumed), or by a bare leaf name when that
name is unique across the schema.

`src/nim_utils/nim_config_schema.m` is the single source of truth for every parameter, its
default, and its permitted values. [YAML_CONFIG.md](YAML_CONFIG.md) is generated from it, so
the reference cannot drift from the code.

---

## 6. Running Tractography

`runTractography`'s optional third argument is either a tracker name or a config struct.
When a config is given, `tractography.algorithm` selects the tracker.

### By Tracker Name

```matlab
% Standard FACT (the default when no argument is given)
runTractography('output/subject.mat');
runTractography('output/subject.mat', 'standard');

% HINEC interpolated streamlines
runTractography('output/subject.mat', 'hinec');
```

### By Configuration

```matlab
% Cubic interpolation with adaptive RKF45
config = load_config_yaml('config/hinec_dti_cubic.yml');
runTractography('output/subject.mat', config);

% Connection-form Method of Moving Frames
config = load_config_yaml('config/mmf_dti.yml');
runTractography('output/subject.mat', config);
```

### Direct Function Calls (Advanced)

The trackers themselves still read a **flat** option struct with legacy names.
`nim_config_to_options` is the single translation point from the nested config surface to
those names; call it rather than hand-assembling a struct, so a direct call stays consistent
with what the launchers do.

```matlab
data = load('output/subject.mat');
nim  = data.nim;

config  = load_config_yaml('config/hinec_dti_cubic.yml');
options = nim_config_to_options(config);

% A seed mask is required — the trackers do no seeding decisions of their own
options.seed_mask = nim.mask > 0.5 & nim.FA >= options.seed_fa_threshold;

% ACT is opt-in: the tracker decides purely from whether it received tissue masks
if options.act_enabled
    options.wm_mask  = nim.wm_mask;
    options.gm_mask  = nim.gm_mask;
    options.csf_mask = nim.csf_mask;
end

tracks = nim_tractography_hinec(nim, options);
```

Note that `options.integration_order` in the flat struct is a **method selector**
(1/2/4/5 → euler/rk2/rk4/rkf45) despite its name, and `options.max_steps` is derived as
`ceil(max_arc / step)`.

### Understanding Track Output

`tracks` is a cell array. Each cell is an Nx3 matrix of **voxel-space** coordinates holding
the complete trajectory, not just the endpoints; N varies per fiber.

```matlab
data = load('tractography_results/tracks_hinec_2025-01-15_14_30_00.mat');
tracks = data.tracks;

fprintf('Total tracks: %d\n', numel(tracks));

% Examine a single track
track = tracks{1};
fprintf('Points in track: %d\n', size(track, 1));
fprintf('Start: [%.1f, %.1f, %.1f]\n', track(1,:));
fprintf('End:   [%.1f, %.1f, %.1f]\n', track(end,:));

% Arc lengths in voxels. Point counts are not comparable across runs when
% output.arc_step or an adaptive integrator changes the point spacing.
arclen = cellfun(@(t) sum(vecnorm(diff(t, 1, 1), 2, 2)), tracks);
fprintf('Mean arc length: %.1f voxels\n', mean(arclen));
fprintf('Max arc length:  %.1f voxels\n', max(arclen));
```

Scoring tools expect RAS world coordinates in millimetres. `bin/run_ismrm_scoring.sh` converts
a run's tracks to TRK using the DWI affine from the run's `intermediate/` directory before
handing them to the scorer.

---

## 7. Visualization

### Quick 3D View

```matlab
% Load data
load('output/subject.mat');

% Eigenvector visualization
nim_plot(nim, 'mode', 'single');                  % Voxel-indexed view (default)
nim_plot(nim, 'mode', 'parcel', 'parcel_id', 5);  % One specific parcel
nim_plot(nim, 'mode', 'parcels');                 % All parcels (separate figures)
```

### Tractography Visualization

```matlab
% Whole brain
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'whole');

% Single region, or several
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'region', 'region', 5);
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'region', 'region', [5 10 15]);

% All regions in grid
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'grid');

% Auto-detect tracks and nim from a run directory
visualizeTractography('hinec_runs/latest');

% Export
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'whole', ...
    'export_dir', 'figures/', 'export_format', 'png', 'export_dpi', 600);
```

### CLI Visualization Export (Shell Script)

Export figures headlessly without an interactive MATLAB session:

```bash
# From a run directory (auto-detects tracks + nim files)
./bin/run_visualization.sh hinec_runs/latest/

# Restrict to specific regions
./bin/run_visualization.sh hinec_runs/latest/ png '23,24,25,26'

# PDF at high DPI
./bin/run_visualization.sh hinec_runs/latest/ pdf '' 600

# Explicit file paths and output directory
./bin/run_visualization.sh tracks.mat nim.mat figures/euler png
```

Two argument forms:

```
./bin/run_visualization.sh <run_dir>                     [format] [region] [dpi]
./bin/run_visualization.sh <tracks_file> <nim_file> <output_dir> [format] [region] [dpi]
```

| Argument | Default | Options |
|----------|---------|---------|
| format | `png` | `png`, `pdf`, `eps` |
| region | all tracks | Region ID or comma-separated list, e.g. `'5,10,15'` |
| dpi | `300` | Any positive integer |

The script renders all eight anatomical viewing angles (superior, inferior, anterior,
posterior, left, right, oblique left, oblique right). In run-directory mode, output goes to
`<run_dir>/figures/`.

### Interactive Slice Viewer

```matlab
visualizeTractographySlices('tracks.mat', 'nim.mat');
```

### Fast Distributed Viewer

Slice images can be pre-rendered on the machine that holds the data, then browsed locally
without MATLAB.

```matlab
% On the server, in MATLAB: generate the cache
generateSlices('tracks.mat', 'nim.mat', '/path/to/cache');
```

```bash
# Or from the shell, as a background batch job
./bin/run_generateSlices.sh tracks.mat nim.mat /path/to/cache

# Transfer the cache to the local machine
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

Because tractography-only runs share one preprocessed `nim`, comparing two tracking
configurations means comparing their track files:

```matlab
runs = dir('hinec_runs/run_*');
a = load(fullfile('hinec_runs', runs(1).name, 'tractography', 'tracks.mat'));
b = load(fullfile('hinec_runs', runs(2).name, 'tractography', 'tracks.mat'));
fprintf('Tracks: %d vs %d\n', numel(a.tracks), numel(b.tracks));
```

Their exact parameters are recorded in each run directory's `config.yml` and, for `--set`
sweeps, `overrides.txt`:

```bash
diff hinec_runs/run_A/config.yml hinec_runs/run_B/config.yml
cat hinec_runs/run_*/overrides.txt
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

1. Restrict seeding to specific regions with `seeding.roi`
2. Reduce seeding density: `--set seeding.density=2`
3. Increase `output.arc_step` to resample saved streamlines — this shrinks the track file
   without affecting integration accuracy
4. Close other applications to free RAM
5. Use the `hinec_dti_fast.yml` preset for initial screening

`nim_save` switches to the v7.3 MAT format automatically for large `nim` structs, so the
2 GB per-variable limit of the older format is not what you are hitting.

### String vs Char Array Errors

**Symptom**: `Error using fullfile` or path concatenation errors

**Cause**: MATLAB string (`"text"`) vs char array (`'text'`) incompatibility.

**Solution**: Convert strings to char arrays when passing to file operations:
```matlab
path = char(string_path);  % Convert if needed
```

See `docs/STRING_CHAR_FIXES.md` in the repository for the string/char
post-mortem; it is a development note and is not part of this site.

### No Tracks Generated

**Symptom**: `runTractography` produces 0 tracks

**Possible causes**:

1. **Termination FA too high**: lower `termination.fa_min` (the seeding threshold,
   `seeding.fa_min`, is a separate knob — they are not interchangeable)
2. **`min_arc` discarding everything**: `termination.min_arc` is an arc length in voxels and
   is enforced; a value above the tracks your parameters actually produce silently empties
   the output
3. **Empty seed mask**: check that the brain mask, or the `seeding.roi` regions you named,
   resolve to a non-empty set of voxels
4. **Data quality**: verify the FA map has reasonable values (max > 0.3)
5. **Preprocessing issues**: check that eigenvectors were computed correctly

The retired key `fa_threshold` no longer does anything. If a config still carries it, the
loader warns and drops it rather than migrating it, because it was never equivalent to either
of the two thresholds that replaced it.

**Diagnostic**:
```matlab
data = load('output/subject.mat');
fprintf('Max FA: %.3f\n', max(data.nim.FA(:)));
fprintf('Mask voxels: %d\n', sum(data.nim.mask(:)));
fprintf('Has eigenvectors: %d\n', isfield(data.nim, 'evec'));
```

### ACT Reported as Disabled

**Symptom**: `runTractography` prints `ACT disabled`

There are two distinct reasons, and the log line says which one applies:

- `ACT disabled by config (tractography.act: false)` — `act` defaults to false. The tissue
  masks are present on the `nim`, and the line reports their voxel counts. Set `act: true`
  under `tractography:` to enable it.
- `ACT disabled: tissue masks not present in the nim` — the `nim` predates tissue
  segmentation, or segmentation failed. Re-run `main()` to generate the masks.

Regenerating masks you already have will not fix the first case.

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
