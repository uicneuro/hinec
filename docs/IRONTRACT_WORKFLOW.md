# IronTract Challenge Workflow

Complete workflow for processing IronTract challenge data and generating submission files.

## Overview

The HINEC pipeline now supports IronTract challenge data with automatic format detection and specialized submission generation. IronTract data uses a different file format than standard DTI data.

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

IronTract data is already preprocessed, but needs brain mask and parcellation:

```matlab
% Process the data (generates brain mask, parcellation, DTI metrics)
main('ironTract/sub-MR243', 'ironTract.mat');
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
# Standard tractography
./run_tractography.sh ironTract.mat

# IronTract submission mode
./run_tractography.sh ironTract.mat IronTract ironTract/injection.nii.gz ironTract_submissions/
```

**Option B: Using MATLAB Directly**

```matlab
% Standard tractography only
runTractography('ironTract.mat');

% IronTract submission mode
runTractography('ironTract.mat', 'IronTract', 'ironTract/injection.nii.gz', 'ironTract_submissions/');
```

**What happens**:

- Loads processed `nim` structure
- Runs FACT tractography with optimized parameters
- Saves tracks to `tractography_results/tracks_standard.mat`
- **If IronTract mode**: Generates submission files with parameter sweep

### Step 3: IronTract Submission Generation

**Automatic in IronTract Mode**:
When you run tractography with IronTract arguments, it automatically:

1. **Filters Tracks**: Only keeps streamlines passing through injection site
2. **Parameter Sweep**: Generates multiple volumes with varying angle thresholds
3. **Binary Masks**: Converts tracks to binary voxel masks
4. **Saves Submissions**: Creates numbered files `submission_001.nii.gz`, `submission_002.nii.gz`, etc.

**Parameter Sweep Defaults**:

- Angle thresholds: `[30, 45, 60, 75, 90]` degrees
- Generates 5 submission volumes for ROC analysis
- Each volume uses same tracks but different filtering

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
matlab -batch "addpath(genpath('.')); main('ironTract/sub-MR243', 'ironTract.mat');"

# 2. Run tractography and generate submissions
./run_tractography.sh ironTract.mat IronTract ironTract/injection.nii.gz ironTract_submissions/

# 3. Monitor progress
tail -f tractography_*.log
```

## Output Files

### Standard Tractography
```
tractography_results/
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

Current FACT parameters (optimized for comprehensive coverage):

```matlab
options.seed_density = 4;           % 4 seeds per voxel
options.step_size = 0.5;            % 0.5mm Euler steps
options.fa_threshold = 0.15;        % FA threshold for seeding
options.termination_fa = 0.15;      % FA stopping criterion
options.angle_thresh = 35;          % 35° maximum turning angle
options.max_steps = 1000;           % Maximum steps per track
options.min_length = 35;            % 35mm minimum track length
options.order = 1;                  % First-order Euler integration
options.interp_method = 'none';     % FACT uses nearest neighbor
```

**Seeding Strategy**: Brain mask-based comprehensive seeding (not FA-threshold seeding)

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
elseif length(first_line_vals) > 3 || length(lines) == 3
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

### Issue: "Dimension mismatch" errors

**Cause**: Incorrect filtering of bval/bvec arrays

**Solution**: Updated code maintains unfiltered arrays in `nim` structure, downstream functions filter when needed

### Issue: "Compressed NIfTI files are not supported"

**Cause**: SPM doesn't handle .nii.gz directly

**Solution**: Automatic gunzip wrapper added to preprocessing functions

### Issue: "Variable not saved" for large files

**Cause**: Default MAT-file format limited to 2GB

**Solution**: Automatic v7.3 format for files >1900MB

### Issue: Label file not found

**Cause**: Parcellation generates XML labels, not MAT files

**Solution**: Added XML label loading support to `nim_load_labels.m`

## References

- Maffei et al. (2022). Insights from the IronTract challenge: Optimal methods for mapping the human corticospinal tract. NeuroImage.
- IronTract Challenge: https://tractometer.org/irontract/

## File Changes Summary

Modified files for IronTract support:

1. `nim_utils/nim_read.m` - Multi-format detection
2. `nim_utils/nim_load_labels.m` - XML label support
3. `nim_utils/nim_save.m` - Large file support (v7.3)
4. `nim_preprocessing/preproc_mask_improvement.m` - Compressed NIfTI handling
5. `runTractography.m` - IronTract mode integration, visualization removed
6. `run_tractography.sh` - IronTract mode shell wrapper
7. `nim_challenges/nim_irontract_submit.m` - Submission generation (already existed, now integrated)
