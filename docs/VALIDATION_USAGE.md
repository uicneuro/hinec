# ISMRM Tractography Validation - Usage Guide

## Overview

The `validate_ismrm_tractography.py` script provides comprehensive validation of HINEC tractography against ISMRM 2015 ground truth data.

**What it does**:
1. ✅ Converts HINEC tracks from MATLAB format to world coordinates
2. ✅ Segments bundles using ROI-based anatomical constraints
3. ✅ Compares with ISMRM ground truth tractography
4. ✅ Computes accuracy metrics (coverage, overreach, Dice)
5. ✅ Generates HTML report with visualizations
6. ✅ Validates coordinate system alignment

## Quick Start

### Prerequisites

```bash
# Install required Python packages
pip install nibabel scipy dipy
```

### Basic Usage

```bash
python validate_ismrm_tractography.py \
    --nim-file sample_parcellated.mat \
    --tracks-file tractography_results/tracks_standard.mat \
    --scoring-dir ISMRM_data/scoring_data_Renauld2023 \
    --output-dir validation_results/
```

### Expected Output

```
validation_results/
├── validation_report.html          # Interactive HTML report
├── validation_results.json         # Machine-readable results
├── CA_segmented.trk                # Segmented bundles (one per bundle)
├── CC_temporal_segmented.trk
├── Cingulum_left_segmented.trk
└── ...
```

## Detailed Usage

### Command-Line Arguments

| Argument | Required | Description | Example |
|----------|----------|-------------|---------|
| `--nim-file` | Yes | HINEC nim structure (.mat) | `sample_parcellated.mat` |
| `--tracks-file` | Yes | HINEC tractography results (.mat) | `tractography_results/tracks_standard.mat` |
| `--scoring-dir` | Yes | ISMRM scoring data directory | `ISMRM_data/scoring_data_Renauld2023` |
| `--output-dir` | No | Output directory (default: `validation_results`) | `my_validation/` |

### Input Requirements

**NIM File** (`--nim-file`):
- MATLAB `.mat` file with nim structure
- Must contain `hdr` field with Transform information
- Used for voxel-to-world coordinate conversion

**Tracks File** (`--tracks-file`):
- MATLAB `.mat` file with tractography results
- Must contain `tracks` field (cell array of Nx3 matrices)
- Tracks stored as 1-based MATLAB voxel indices

**Scoring Directory** (`--scoring-dir`):
- Must contain ISMRM scoring data structure:
  ```
  scoring_data_Renauld2023/
  ├── bundles/              # Ground truth .trk files
  ├── ROI/                  # Anatomical masks
  │   ├── all_masks/
  │   ├── any_masks/
  │   └── endpoints/
  ├── config_file_segmentation.json
  ├── config_file_tractometry.json
  └── t1.nii.gz             # Reference T1 image
  ```

## Output Files

### 1. HTML Report (`validation_report.html`)

Interactive report with:
- **Summary statistics**: Valid bundles, coverage, overreach
- **Per-bundle table**: Detailed metrics for each bundle
- **Coordinate diagnostics**: Voxel and world coordinate ranges
- **Interpretation guide**: How to read the metrics

**Open in browser**:
```bash
open validation_results/validation_report.html
```

### 2. JSON Results (`validation_results.json`)

Machine-readable results for further analysis:

```json
{
  "bundles": {
    "CA": {
      "name": "CA",
      "user_count": 45,
      "gt_count": 540,
      "coverage": 78.5,
      "overreach": 12.3,
      "dice": 0.82,
      "avg_length_user": 45.2,
      "avg_length_gt": 48.1,
      "valid": true
    },
    ...
  },
  "summary": {
    "total_bundles": 25,
    "valid_bundles": 18,
    "total_user_tracks": 45289,
    "total_gt_tracks": 180345,
    "avg_coverage": 65.3,
    "avg_overreach": 18.7
  },
  "diagnostics": {
    "input_coords": {...},
    "world_coords": {...}
  }
}
```

### 3. Segmented Bundles (`*_segmented.trk`)

TrackVis format files for each bundle that passed ROI filtering:
- Can be visualized in TrackVis, MI-Brain, or DSI Studio
- Used for visual validation of segmentation accuracy
- Compatible with further analysis tools

## Metrics Explained

### Coverage (%)
**Definition**: Percentage of ground truth streamlines that are covered by user tractography

**Formula**: `(# GT streamlines within 10mm of user tracks) / (total GT streamlines) × 100`

**Interpretation**:
- **>70%**: Excellent coverage - capturing most of ground truth
- **50-70%**: Good coverage - major pathways captured
- **30-50%**: Moderate coverage - missing some pathways
- **<30%**: Poor coverage - significant anatomical gaps

### Overreach (%)
**Definition**: Percentage of user streamlines that don't match ground truth

**Formula**: `(# user streamlines >10mm from any GT) / (total user streamlines) × 100`

**Interpretation**:
- **<20%**: Excellent specificity - few spurious tracks
- **20-35%**: Good specificity - moderate false positives
- **35-50%**: Moderate specificity - significant false positives
- **>50%**: Poor specificity - many invalid tracks

### Dice Coefficient
**Definition**: Spatial overlap measure (0 = no overlap, 1 = perfect overlap)

**Formula**: `2 × (# matched streamlines) / (# user + # GT streamlines)`

**Interpretation**:
- **>0.8**: Excellent agreement
- **0.6-0.8**: Good agreement
- **0.4-0.6**: Moderate agreement
- **<0.4**: Poor agreement

### Valid Bundles
**Definition**: Bundles with at least one streamline passing all ROI constraints

**ROI Constraints** (per bundle):
- **All mask**: All points must pass through
- **Any mask**: At least one point must intersect
- **Endpoints**: Start/end in correct anatomical regions
- **Length**: Geometric constraints (total length, axis-specific extents)

## Bundle Segmentation Details

### ROI-Based Filtering Process

For each bundle (e.g., "Cingulum_left"):

1. **Load Anatomical Masks**:
   - `all_masks/Cingulum_left.nii.gz` - Required pathway
   - `endpoints/Cingulum_left_head.nii.gz` - Start region
   - `endpoints/Cingulum_left_tail.nii.gz` - End region

2. **Apply Constraints**:
   ```
   For each streamline:
     ✓ All points pass through all_mask?
     ✓ First/last points in head/tail regions?
     ✓ Length within [min, max] range?
     ✓ X-axis extent within constraints?
     → Keep if ALL checks pass
   ```

3. **Compare with Ground Truth**:
   - Load `bundles/Cingulum_left.trk`
   - Compute distance between user and GT streamlines
   - Calculate coverage and overreach metrics

### Example: Corpus Callosum U-shaped

**Config** (`config_file_segmentation.json`):
```json
{
  "CC_u_shaped": {
    "all_mask": "ROI/all_masks/CC_u_shaped.nii.gz",
    "any_mask": "ROI/any_masks/CC_u_shaped_inclusion_mask.nii.gz",
    "head": "ROI/endpoints/CC_striatal_left.nii.gz",
    "tail": "ROI/endpoints/CC_striatal_right.nii.gz",
    "length": [70, 1000],       // Total length 70-1000mm
    "length_y": [0, 32],         // Anterior extent <32mm
    "length_x_abs": [35, 1000]   // Lateral span >35mm
  }
}
```

**Filtering Logic**:
1. Keep streamlines passing through CC u-shaped mask (entire pathway)
2. Keep streamlines with at least one point in inclusion mask
3. Keep streamlines with endpoints in left and right striatal regions
4. Keep streamlines with 70-1000mm total length
5. Keep streamlines with <32mm anterior extent
6. Keep streamlines spanning >35mm laterally

## Coordinate System Validation

### Automatic Diagnostics

The script automatically checks coordinate system alignment:

**Input Coordinates** (0-based voxels):
```
X: [0.5, 95.8]
Y: [1.2, 110.5]
Z: [0.3, 94.7]
```

**World Coordinates** (after affine transform):
```
X: [-63.5, 68.2] mm
Y: [-5.6, 72.1] mm
Z: [6.3, 72.8] mm
```

**Warning Signs** (check diagnostics section in report):
- ⚠️ World coordinates outside expected brain range (-100 to +100 mm)
- ⚠️ Voxel coordinates exceeding image dimensions
- ⚠️ All bundles showing 0% coverage (complete misalignment)

### Manual Validation

If you get suspicious results (all 0% coverage):

1. **Visual Check**:
   ```bash
   # Convert one bundle to TRK for visualization
   python hinec_to_trk.py \
       tractography_results/tracks_standard.mat \
       ISMRM_data/scoring_data_Renauld2023/t1.nii.gz \
       test_tracks.trk --space rasmm

   # Load in TrackVis with T1 overlay
   trackvis test_tracks.trk -v ISMRM_data/scoring_data_Renauld2023/t1.nii.gz
   ```

2. **Check Alignment**:
   - Do tracks follow white matter pathways?
   - Are tracks inside the brain?
   - Do major bundles (corpus callosum, cingulum) look anatomically correct?

## Example Workflows

### Workflow 1: Basic Validation

```bash
# 1. Run HINEC tractography (if not done)
matlab -batch "runTractography('sample_parcellated.mat')"

# 2. Validate against ISMRM
python validate_ismrm_tractography.py \
    --nim-file sample_parcellated.mat \
    --tracks-file tractography_results/tracks_standard.mat \
    --scoring-dir ISMRM_data/scoring_data_Renauld2023 \
    --output-dir validation_results/

# 3. View report
open validation_results/validation_report.html
```

### Workflow 2: Parameter Optimization

Test different tractography parameters:

```bash
# Run tractography with different parameters
for fa_thresh in 0.10 0.15 0.20; do
    for angle in 35 45 60; do
        output_dir="validation_fa${fa_thresh}_angle${angle}"

        # Run tractography with specific parameters
        matlab -batch "options.fa_threshold=$fa_thresh; options.angle_thresh=$angle; \
                       tracks=nim_tractography_standard('data.mat', options); \
                       save('tracks_fa${fa_thresh}_a${angle}.mat', 'tracks');"

        # Validate
        python validate_ismrm_tractography.py \
            --nim-file data.mat \
            --tracks-file "tracks_fa${fa_thresh}_a${angle}.mat" \
            --scoring-dir ISMRM_data/scoring_data_Renauld2023 \
            --output-dir "$output_dir"
    done
done

# Compare results
python compare_validation_results.py validation_*/validation_results.json
```

### Workflow 3: Bundle-Specific Analysis

Focus on specific bundles:

```bash
# 1. Run full validation
python validate_ismrm_tractography.py ...

# 2. Extract specific bundle for detailed analysis
python -c "
import nibabel as nib
from nibabel.streamlines import load
import sys

# Load segmented bundle
sft = load('validation_results/Cingulum_left_segmented.trk', 'same')
print(f'Cingulum Left:')
print(f'  Streamlines: {len(sft.streamlines)}')
print(f'  Length range: {[len(s) for s in sft.streamlines]}')

# Save statistics
with open('cingulum_stats.txt', 'w') as f:
    f.write(f'Bundle: Cingulum Left\n')
    f.write(f'Count: {len(sft.streamlines)}\n')
"
```

## Troubleshooting

### Issue: All bundles show 0% coverage

**Cause**: Coordinate system misalignment

**Solutions**:
1. Check coordinate diagnostics in report
2. Visually validate using TrackVis
3. Verify T1 reference image matches tractography
4. Check if nim.hdr.Transform is correctly set

### Issue: "No tracks found in .mat file"

**Cause**: Field name mismatch

**Solutions**:
1. Check .mat file structure: `scipy.io.whosmat('tracks.mat')`
2. Ensure field is named `tracks` or `tracts`
3. Verify tracks are stored as cell array, not plain array

### Issue: High overreach (>50%) but good coverage

**Interpretation**: Algorithm is capturing ground truth but also generating many extra streamlines

**Solutions**:
1. Increase FA threshold to reduce spurious tracks
2. Reduce step size for more anatomically constrained tracking
3. Apply stricter angle constraints
4. Use more conservative termination criteria

### Issue: Low coverage (<30%) but low overreach

**Interpretation**: Algorithm is conservative, missing valid pathways

**Solutions**:
1. Decrease FA threshold to track in lower anisotropy regions
2. Increase seeding density
3. Relax angle constraints
4. Check if seeding covers all brain regions

## Advanced Analysis

### Export Results to CSV

```python
import json
import csv

# Load results
with open('validation_results/validation_results.json') as f:
    results = json.load(f)

# Export per-bundle metrics
with open('bundle_metrics.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Bundle', 'User Tracks', 'GT Tracks', 'Coverage %', 'Overreach %', 'Dice', 'Avg Length'])

    for bundle_name, metrics in results['bundles'].items():
        writer.writerow([
            bundle_name,
            metrics['user_count'],
            metrics['gt_count'],
            f"{metrics['coverage']:.1f}",
            f"{metrics['overreach']:.1f}",
            f"{metrics['dice']:.3f}",
            f"{metrics['avg_length_user']:.1f}"
        ])
```

### Compare Multiple Validation Runs

```python
import json
import numpy as np

runs = ['run1', 'run2', 'run3']
coverages = []

for run in runs:
    with open(f'{run}/validation_results.json') as f:
        data = json.load(f)
        coverages.append(data['summary']['avg_coverage'])

print(f"Average coverage across runs: {np.mean(coverages):.1f}% ± {np.std(coverages):.1f}%")
```

## References

- ISMRM 2015 Tractography Challenge: http://tractometer.dinf.usherbrooke.ca/ismrm_2015_challenge/
- Validation methodology: See `ISMRM_SCORING_ANALYSIS.md`
- HINEC tractography: See `CLAUDE.md` for pipeline details
