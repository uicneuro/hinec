# `validate_ismrm_tractography.py` — Usage

`scripts/validate_ismrm_tractography.py` is a standalone diagnostic: it segments
a HINEC tractogram with the ISMRM 2015 ROI gates, compares each segmented bundle
with the ground truth by streamline distance, and writes an HTML report and a
JSON summary.

!!! warning "This is a diagnostic, not the scorer"
    The authoritative scoring path is `bin/run_ismrm_scoring.sh`, which runs
    scilpy's Renauld 2023 scorer — see the
    [Validation Guide](VALIDATION_GUIDE.md). This script reimplements the ROI
    gating in Python with its own metrics; its coverage, overreach and Dice are
    **streamline-count** quantities computed from pairwise streamline distances,
    whereas the challenge metrics are **voxel-mask** quantities. The two do not
    produce the same numbers and neither this script's output nor its Dice may
    be reported as a challenge score. Known limitations are listed at the bottom
    of this page.

    The same whole-brain caveat applies here: the ROI gates define 26 bundles,
    and a tractogram seeded from a single ROI can populate at most one of them.

What the script does:

1. Loads HINEC tracks from a MATLAB `.mat` (v7 or v7.3).
2. Segments bundles using the ROI constraints in
   `config_file_segmentation.json`.
3. Converts each segmented bundle to world coordinates.
4. Compares with the ground-truth `.trk` bundles by streamline distance.
5. Writes an HTML report, a JSON summary, and one `.trk` per non-empty bundle.
6. Reports coordinate ranges before and after the world transform, as a
   space-alignment check.

## Quick Start

### Prerequisites

```bash
pip install nibabel scipy dipy
```

`dipy` is optional but required for the distance metrics: without it the script
still segments bundles and writes TRK files, but coverage, overreach and Dice
all stay at 0.

### Basic Usage

```bash
python scripts/validate_ismrm_tractography.py \
    --nim-file sample_parcellated.mat \
    --tracks-file tractography_results/tracks_standard.mat \
    --scoring-dir data/ismrm2015/scoring_data_Renauld2023 \
    --output-dir validation_results/
```

### Expected Output

```
validation_results/
├── validation_report.html          # HTML report
├── validation_results.json         # machine-readable results
├── CA_segmented.trk                # one per bundle that kept any streamline
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
| `--scoring-dir` | Yes | ISMRM scoring data directory | `data/ismrm2015/scoring_data_Renauld2023` |
| `--output-dir` | No | Output directory (default: `validation_results`) | `my_validation/` |

### Input Requirements

**Tracks file** (`--tracks-file`)

- MATLAB `.mat` with a `tracks` field: a cell array of N×3 matrices.
- Coordinates are 1-based MATLAB voxel indices in DWI space.

**Nim file** (`--nim-file`)

- Required by the argument parser, but used only as a provenance label in the
  HTML report. The voxel→world transform is taken from the scoring
  `t1.nii.gz`, not from the nim header.

**Scoring directory** (`--scoring-dir`)

```
scoring_data_Renauld2023/
├── bundles/              # ground-truth .trk files
├── ROI/
│   ├── all_masks/
│   ├── any_masks/
│   └── endpoints/
├── config_file_segmentation.json
├── config_file_tractometry.json
└── t1.nii.gz             # reference T1, 1 mm
```

## Output Files

### 1. HTML report (`validation_report.html`)

Summary statistics, a per-bundle metric table, the coordinate diagnostics, and a
short guide to reading the numbers.

```bash
open validation_results/validation_report.html
```

### 2. JSON results (`validation_results.json`)

The schema, with illustrative values — these are not measured results:

```json
{
  "bundles": {
    "CA": {
      "name": "CA",
      "user_count": 45,
      "gt_count": 430,
      "coverage": 78.5,
      "overreach": 12.3,
      "dice": 0.82,
      "avg_length_user": 45.2,
      "avg_length_gt": 48.1,
      "valid": true
    }
  },
  "summary": {
    "total_bundles": 26,
    "valid_bundles": 18,
    "total_user_tracks": 45289,
    "total_gt_tracks": 180345,
    "avg_coverage": 65.3,
    "avg_overreach": 18.7
  },
  "diagnostics": {
    "input_coords": {},
    "world_coords": {}
  }
}
```

`gt_count` is fixed by the data: `CA` has 430 ground-truth streamlines, `CP`
365, `MCP` 21008.

### 3. Segmented bundles (`*_segmented.trk`)

One TrackVis file per bundle that retained at least one streamline, in world
coordinates. Loadable in TrackVis, MI-Brain or DSI Studio for visual inspection
of the segmentation.

## Metrics

All three metrics are computed from a pairwise streamline-distance matrix
(Dipy's `bundles_distances_mam`, after resampling every streamline to 20
points), with a fixed 10 mm match threshold.

| metric | definition |
|---|---|
| coverage (%) | ground-truth streamlines with a produced streamline within 10 mm, over all ground-truth streamlines |
| overreach (%) | produced streamlines with no ground-truth streamline within 10 mm, over all produced streamlines |
| dice | `2 × (matched ground-truth streamlines) / (produced + ground-truth streamline counts)` |
| valid | the bundle kept at least one streamline through the ROI gates |

Two properties are worth keeping in mind when reading them. The Dice numerator
counts matched *ground-truth* streamlines only, so the measure is asymmetric and
is not the voxel Dice the challenge reports. And all three depend on the
streamline counts on both sides, so a sparsely seeded run and a densely seeded
one are not directly comparable.

There are no calibrated thresholds for these numbers on this phantom. Use them
to compare runs of this script against each other, and use
`bin/run_ismrm_scoring.sh` for anything that needs to be comparable to published
work.

## Bundle segmentation

For each entry in `config_file_segmentation.json` — e.g. `Cingulum_left`, whose
gates are `all_masks/Cingulum_left.nii.gz`,
`endpoints/Cingulum_left_head.nii.gz`, `endpoints/Cingulum_left_tail.nii.gz` and
`length_x: [0, 40]` — a streamline is kept only if:

```
✓ every point lies inside all_mask
✓ at least one point lies inside any_mask   (when the bundle declares one)
✓ one endpoint in head and the other in tail (either way round)
✓ length and per-axis extents within the declared ranges
```

Surviving streamlines are then converted to world coordinates and compared with
`bundles/<name>.trk`.

## Coordinate diagnostics

The script prints voxel and world coordinate ranges, both of which land in the
JSON under `diagnostics`. Expect world coordinates within roughly ±100 mm of the
origin, and voxel coordinates inside the image dimensions. All bundles reporting
0% coverage is the signature of a space mismatch rather than a tracking failure.

To check by eye:

```bash
python scripts/hinec_to_trk.py \
    tractography_results/tracks_standard.mat \
    data/ismrm2015/ismrm2015.nii.gz \
    test_tracks.trk --space rasmm

trackvis test_tracks.trk -v data/ismrm2015/scoring_data_Renauld2023/t1.nii.gz
```

Then confirm the tracks sit inside the brain and that the corpus callosum and
cingulum look anatomically plausible.

## Example workflows

### Basic validation

```bash
# 1. Tractography (config-driven; writes into a timestamped run dir)
./bin/run_tractography.sh hinec_dti_cubic

# 2. Validate against ISMRM
python scripts/validate_ismrm_tractography.py \
    --nim-file data/ismrm2015/ismrm2015.mat \
    --tracks-file hinec_runs/run_<TIMESTAMP>_hinec_dti_cubic/tractography/tracks_hinec_<TS>.mat \
    --scoring-dir data/ismrm2015/scoring_data_Renauld2023 \
    --output-dir validation_results/

# 3. View report
open validation_results/validation_report.html
```

### Parameter sweep

Override any config key on the command line with `--set`; each run lands in its
own tagged run directory.

```bash
for fa in 0.10 0.15 0.20; do
  for angle in 175 225 300; do
    ./bin/run_tractography.sh hinec_dti_cubic \
        --set termination.fa_min=$fa --set termination.angle_max=$angle
  done
done
```

Then validate each run directory, or score them properly with
`./bin/run_ismrm_scoring.sh <run_dir>` and tabulate with
`scripts/compare_ismrm_results.py`.

### Bundle-specific inspection

```bash
python -c "
from nibabel.streamlines import load
sft = load('validation_results/Cingulum_left_segmented.trk')
print('streamlines:', len(sft.streamlines))
print('points per streamline:', min(len(s) for s in sft.streamlines),
      '-', max(len(s) for s in sft.streamlines))
"
```

## Troubleshooting

### All bundles show 0% coverage

Almost always a coordinate mismatch. Check the coordinate diagnostics in the
report, inspect the tracks in TrackVis against the scoring T1, and confirm the
tracks were written from the DWI the nim was built on.

### "No tracks found in .mat file"

The field name did not match. Inspect with `scipy.io.whosmat('tracks.mat')`; the
script expects `tracks` or `tracts`, stored as a cell array rather than a plain
numeric array.

### High overreach with good coverage

The tracker is finding the bundles and producing a great deal besides. Tighten
termination: raise `termination.fa_min`, lower `termination.angle_max`.

### Low coverage with low overreach

The tracker is too conservative. Lower `termination.fa_min`, raise
`seeding.density`, loosen `termination.angle_max`, and check that seeding covers
the regions of interest (`seeding.roi`, `seeding.roi_dilate`).

## Coordinate handling

Streamlines are mapped into the ROI masks through **world (RAS) coordinates**,
using each file's own affine — the DWI reference's to leave voxel space, and the
mask's to enter it. Pass `--dwi-ref <dwi.nii.gz>` to name the grid the tracks are
stored on; if omitted the script looks for a sibling `*_dwi_ref.nii.gz` beside
`--nim-file` and fails loudly if it finds none.

Length constraints (`length`, `length_x`, `length_y`, `length_x_abs`) are
specified in millimetres in `config_file_segmentation.json` and are evaluated on
those world coordinates.

!!! note "Both of these were wrong until recently"
    The script used to multiply track coordinates by a hardcoded 2.0 to index the
    1 mm masks. On the ISMRM data that happens to agree with the correct mapping
    — the DWI is 2 mm at zero offset, the masks are 1 mm at a half-voxel offset,
    and rounding absorbs the difference — so it produced no visible error. It was
    still wrong for any other pair of grids.

    The length bug was not harmless. Bounds in millimetres were compared against
    lengths measured in 2 mm voxel units, so every geometric gate fired at half
    its threshold: a 70 mm minimum rejected any streamline shorter than 140 mm.
    This affected the four bundles that carry length constraints — `CC_u_shaped`,
    `Cingulum_left`, `Cingulum_right` and `MCP`.

This script remains a **diagnostic**. `bin/run_ismrm_scoring.sh` is authoritative
for scores.


## References

- [Validation Guide](VALIDATION_GUIDE.md) — the authoritative scoring path
- [ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md) — how a bundle is defined
- [ISMRM 2015 Tractography Challenge](http://tractometer.dinf.usherbrooke.ca/ismrm_2015_challenge/)
