# IronTract Challenge Workflow

Complete workflow for processing IronTract challenge data and generating submission files.

## Overview

HINEC reads IronTract challenge data directly: `nim_read` detects its gradient-table layout automatically, and `nim_irontract_submit` packages tractography results as the binary volumes the challenge evaluates. The only difference from a standard DTI dataset is the `.bval`/`.bvec` layout.

## File Format Differences

### Standard DTI Format
```
# .bval (single line)
0 1000 1000 1000 1000 ...

# .bvec (3 lines: X, Y, Z components)
0 0.2 0.4 0.6 ...
0 0 0.2 0.4 ...
0 0 0 0.2 ...
```

### IronTract Format
```
# .bval (one value per line)
0
1601.25
1601.25
1601.25
...

# .bvec (one 3D vector per line)
0.2 0 0
0 0.4 0
0 0 0.6
...
```

**Format Detection**: `nim_read.m` automatically detects both formats - no user action required.

## Complete Workflow

### Step 1: Preprocess IronTract Data

IronTract data arrives already preprocessed, but still needs a brain mask,
parcellation and tensor fit:

```bash
./bin/run_hinec.sh ironTract/sub-MR243 ironTract.mat config/irontract.yml
```

The MATLAB entry point is equivalent:

```matlab
config   = load_config_yaml('config/irontract.yml');
run_info = create_run_directory('config/irontract.yml');
main('ironTract/sub-MR243', 'ironTract.mat', config, run_info);
```

**What happens**:

- Detects IronTract file format automatically
- Generates brain mask using FSL BET
- Applies MNI atlas parcellation
- Computes diffusion tensors and FA maps
- Saves complete `nim` structure to `ironTract.mat`

**Important Files Generated**:

- `ironTract.mat` - Processed nim structure with all metrics
- `ironTract_mask_improved.nii.gz` - Brain mask
- `parcellation_mask.nii.gz` - Atlas parcellation
- `atlas_labels.xml` - Region labels

### Step 2: Run Tractography

**Option A: Using Shell Script (Recommended)**

```bash
# Tractography only, config-driven
./bin/run_tractography.sh irontract --source ironTract.mat

# IronTract submission mode (legacy positional form)
./bin/run_tractography.sh ironTract.mat IronTract ironTract/injection.nii.gz ironTract_submissions/
```

**Option B: Using MATLAB Directly**

```matlab
% Standard tractography only
runTractography('ironTract.mat');

% IronTract submission mode
runTractography('ironTract.mat', 'IronTract', 'ironTract/injection.nii.gz', 'ironTract_submissions/');
```

**What happens**:

- Loads the processed `nim` structure
- Tracks according to the config (or, in IronTract mode, with `nim_tractography_standard` — FACT)
- Writes tracks under the run directory's `tractography/`
- **In IronTract mode**: hands off to `nim_irontract_submit` for the parameter sweep

### Step 3: IronTract Submission Generation

**Automatic in IronTract mode.** For each angle threshold in the sweep,
`nim_irontract_submit`:

1. **Seeds including the injection site** — the brain mask (or, failing that, the
   parcellation mask, or FA > 0.1) is unioned with the injection mask, so the
   site is guaranteed to be seeded.
2. **Tracks** with `nim_tractography_standard` at that angle threshold.
3. **Filters** to streamlines with at least one point inside the injection site.
4. **Rasterises** the survivors into a binary voxel mask.
5. **Writes** `submission_001.nii.gz`, `submission_002.nii.gz`, …

**Sweep defaults**:

- Angle thresholds: `[30, 45, 60, 75, 90]` degrees, giving 5 volumes for ROC
  analysis.
- Each volume is a **separate tracking run** at its own angle threshold. Passing
  `options.tracks_file` short-circuits that: all volumes then reuse the one
  supplied track set, and the sweep degenerates to a single result repeated.

**Manual Submission Generation** (if needed):
```matlab
% Basic usage
nim_irontract_submit('ironTract.mat', 'ironTract/injection.nii.gz', 'submissions/');

% Custom angle sweep
options.angle_thresholds = [20, 40, 60, 80];
nim_irontract_submit('ironTract.mat', 'ironTract/injection.nii.gz', 'submissions/', options);

% Use pre-computed tracks
options.tracks_file = 'tractography_results/tracks_standard.mat';
nim_irontract_submit('ironTract.mat', 'ironTract/injection.nii.gz', 'submissions/', options);
```

## Complete Example

```bash
# 1. Process IronTract data
./bin/run_hinec.sh ironTract/sub-MR243 ironTract.mat config/irontract.yml

# 2. Run tractography and generate submissions
./bin/run_tractography.sh ironTract.mat IronTract ironTract/injection.nii.gz ironTract_submissions/

# 3. Monitor progress
tail -f tractography_*.log
```

## Output Files

### Standard Tractography
```
<run_dir>/tractography/
  tracks_standard.mat         # Complete tractography results
    - tracks: cell array of fiber pathways
    - options: tractography parameters
    - elapsed_time: processing time
```

### IronTract Submissions
```
ironTract_submissions/
  submission_001.nii.gz       # Angle threshold = 30°
  submission_002.nii.gz       # Angle threshold = 45°
  submission_003.nii.gz       # Angle threshold = 60°
  submission_004.nii.gz       # Angle threshold = 75°
  submission_005.nii.gz       # Angle threshold = 90°
```

Each submission file is:

- Binary mask (0 or 1)
- 1 = voxel visited by streamlines from injection site
- Matches injection mask header for proper alignment
- Ready for IronTract challenge evaluation

## Tractography Parameters

Parameters come from the YAML config, in the canonical nested schema. Only
non-default keys need to appear; the rest are filled in from
`nim_config_schema`. `config/irontract.yml` sets:

```yaml
tractography:
  integrator:
    method: rkf45
    tolerance: 0.02
    step_min: 0.02
    step_max: 0.5
  termination:
    angle_max: 225      # degrees per VOXEL OF ARC; 0 disables the criterion
    max_arc: 600        # voxels of arc, not steps
```

Override any of them per run without a new file:

```bash
./bin/run_tractography.sh irontract --source ironTract.mat \
    --set termination.angle_max=175 --set seeding.density=4
```

Submission mode is different: `nim_irontract_submit` calls
`nim_tractography_standard` (FACT), which walks voxel-to-voxel with no
interpolation, so the integrator and interpolation keys do not apply. It sets
the tracker's internal `angle_thresh` from each entry of the sweep and takes the
rest of its options from `base_options`.

**Seeding**: the brain mask unioned with the injection site — comprehensive
rather than FA-thresholded.

## Implementation Details

### Format Detection (nim_read.m)
```matlab
% .bval format detection
if length(bval_lines) == 1
    % Standard format: single line
    nim.bval = transpose(str2num(string(bval_lines(1))));
elseif length(bval_lines) > 1
    % IronTract format: one value per line
    nim.bval = zeros(length(bval_lines), 1);
    for i = 1:length(bval_lines)
        nim.bval(i) = str2double(string(bval_lines(i)));
    end
end

% .bvec format detection
first_line_vals = str2num(string(lines(1)));
if length(first_line_vals) == 3 && length(lines) > 3
    % IronTract format: one vector per line
    bvec_all = zeros(length(lines), 3);
    for i = 1:length(lines)
        bvec_all(i, :) = str2num(string(lines(i)));
    end
elseif length(first_line_vals) > 3 || (length(first_line_vals) >= 1 && length(lines) == 3)
    % Standard format: 3 lines of components
    gx = str2num(string(lines(1)));
    gy = str2num(string(lines(2)));
    gz = str2num(string(lines(3)));
    bvec_all = transpose([gx; gy; gz]);
end
```

### Injection Site Filtering (nim_irontract_submit.m)
```matlab
% For each track
for i = 1:length(tracks)
    track = tracks{i};
    track_voxels = round(track);  % Convert to voxel indices

    % Check if any point intersects injection site
    for j = 1:size(track_voxels, 1)
        pos = track_voxels(j, :);
        if all(pos >= 1) && all(pos <= dims)
            if injection_mask(pos(1), pos(2), pos(3)) > 0
                filtered_tracks{end+1} = track;
                break;
            end
        end
    end
end
```

### Binary Mask Generation (nim_irontract_submit.m)
```matlab
% Mark all voxels visited by filtered tracks
voxel_mask = zeros(dims);
for i = 1:length(filtered_tracks)
    track = filtered_tracks{i};
    track_voxels = round(track);

    for j = 1:size(track_voxels, 1)
        pos = track_voxels(j, :);
        if all(pos >= 1) && all(pos <= dims)
            voxel_mask(pos(1), pos(2), pos(3)) = 1;
        end
    end
end
```

## Troubleshooting

| symptom | cause and resolution |
|---|---|
| "Dimension mismatch" | b-value/b-vector arrays filtered inconsistently. `nim` stores the unfiltered arrays; downstream functions filter as needed, so do not pre-filter them yourself. |
| "Compressed NIfTI files are not supported" | SPM cannot read `.nii.gz` directly. The preprocessing functions gunzip to a temporary file automatically; the error means the wrapper was bypassed. |
| "Variable not saved" on large volumes | The default MAT format caps at 2 GB. `nim_save` switches to v7.3 above ~1900 MB. |
| "Label file not found" | Parcellation emits XML labels, not `.mat`. `nim_load_labels` reads both. |

## References

- Maffei et al. (2022). Insights from the IronTract challenge: Optimal methods for mapping the human corticospinal tract. NeuroImage.
- IronTract Challenge: https://tractometer.org/irontract/

## Where the code lives

| file | role |
|---|---|
| `src/nim_utils/nim_read.m` | gradient-table format detection |
| `src/nim_utils/nim_load_labels.m` | XML atlas label loading |
| `src/nim_utils/nim_save.m` | v7.3 MAT format for large nims |
| `src/nim_preprocessing/preproc_mask_improvement.m` | compressed-NIfTI handling |
| `src/nim_challenges/nim_irontract_submit.m` | submission packaging |
| `runTractography.m` | dispatch into IronTract mode |
| `bin/run_tractography.sh` | shell wrapper |
