# Complete Visualization Reference

Unified guide to all visualization capabilities in HINEC. Covers DTI visualization, 3D tractography, interactive slice viewing, fast distributed viewing, connectivity matrices, and publication-quality exports.

---

## 1. Visualization Overview

HINEC provides multiple visualization tools optimized for different tasks:

| Tool | Type | Use case |
|---|---|---|
| `nim_plot` | DTI eigenvectors | DTI quality check |
| `visualizeTractography` | 3D tractography | whole-brain and per-region exploration |
| `visualizeTractographyAngles` | 3D, batch export | headless figure export from 8 anatomical views |
| `visualizeTractographySlices` | 2D slices | slice-by-slice inspection |
| `generateSlices` + `FastTractographyViewer.py` | pre-computed 2D | fast slice navigation without MATLAB |
| `nim_plot_bundles` | 3D streamline sets | comparing bundles in one shared frame |
| `nim_plot_vs_groundtruth` | 3D overlay | ours against an ISMRM ground-truth bundle |
| `nim_plot_tractography` | basic 3D | quick look at a track set |
| `nim_plot_connectivity_matrix` | 2D heatmap | connectivity analysis |
| `nim_plot_vector_field` | 2D quiver | direction-field inspection |

### Decision Tree: Which Tool to Use

```
What do you need?
├── Quick DTI quality check → nim_plot
├── 3D whole-brain tractography → visualizeTractography('mode', 'whole')
├── Specific brain region analysis → visualizeTractography('mode', 'region')
├── All regions at once → visualizeTractography('mode', 'grid')
├── Step through regions one by one → visualizeTractography('mode', 'sequential')
├── Slice-by-slice inspection → visualizeTractographySlices
├── Fast slice navigation → generateSlices + FastTractographyViewer.py
├── One or more bundles side by side → nim_plot_bundles
├── Ours vs ISMRM ground truth → nim_plot_vs_groundtruth
├── Connectivity between regions → nim_plot_connectivity_matrix
├── Direction field on a slice → nim_plot_vector_field
└── Publication figures → visualizeTractography with 'export_dir'
```

---

## 2. DTI Eigenvector Visualization (`nim_plot`)

**File**: `src/nim_plots/nim_plot.m`

Displays principal eigenvectors color-coded by direction. This is the primary tool for verifying DTI processing quality.

### DTI Color Convention

The standard neuroimaging color coding:

- **Red** = Left/Right (X axis)
- **Green** = Anterior/Posterior (Y axis)
- **Blue** = Superior/Inferior (Z axis)

Colors are derived from the absolute value of the primary eigenvector components using `vector_to_color()`.

### Modes

**Single region** (default):
```matlab
nim_plot(nim, 'mode', 'single');
```
Shows eigenvectors for a single region or the whole volume.

**Single parcel**:
```matlab
nim_plot(nim, 'mode', 'parcel', 'parcel_id', 5);
```
Shows eigenvectors within a specific parcellation region.

**All parcels** (separate figures):
```matlab
nim_plot(nim, 'mode', 'parcels');
```
Creates one figure per parcellation region. Useful for surveying all regions.

### Performance

For large datasets, downsample:
```matlab
nim_plot(nim, 'mode', 'single', 'downsample_factor', 3);  % every 3rd voxel
```

`indx`/`indy`/`indz` restrict `mode: single` to an explicit voxel index range;
`figindex`, `show_figure` and `show_progress` control figure handling.

---

## 3. 3D Tractography Visualization (`visualizeTractography`)

**File**: `src/nim_visualization/visualizeTractography.m`

The primary tool for 3D tractography visualization. It supersedes the older
`visualizeTractographyRegion`, `visualizeTractographyAllRegions` and
`nim_plot_tractography_region` entry points, and supports four modes, several
colour schemes, track filtering, run-directory auto-detection and export.

### Input Flexibility

```matlab
% Direct file paths
visualizeTractography('path/to/tracks.mat', 'path/to/nim.mat');

% From a run directory (auto-detects nim + latest tracks; wildcards allowed)
visualizeTractography('hinec_runs/run_<TIMESTAMP>_<config>/');

% Wildcard support
visualizeTractography('tractography_results/tracks_hinec_*.mat', 'output/*.mat');
```

### Mode: Whole Brain

Complete 3D view of all fiber tracks with summary statistics.

```matlab
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'whole');

% With custom settings
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'whole', ...
    'max_tracts', 10000, ...
    'color_mode', 'direction');
```

### Mode: Single/Multi Region

Focused view of tracks associated with specific brain regions.

```matlab
% Single region
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', ...
    'region', 5);

% Multiple regions
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', ...
    'region', [5, 10, 15]);
```

The parameter is `region` (singular). `visualizeTractographySlices` uses
`regions` for the same idea, so check which function you are calling.

### Mode: Grid

All parcellation regions displayed simultaneously in a grid layout. Useful for comparing regions at a glance.

```matlab
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'grid');
```

### Mode: Sequential

Step through regions one at a time with pause between each.

```matlab
visualizeTractography('tracks.mat', 'nim.mat', 'mode', 'sequential');
```

### Track Coloring

| Color Mode | Description | Best For |
|---|---|---|
| `'direction'` | RGB from eigenvector (R=L/R, G=A/P, B=S/I) | Anatomical orientation |
| `'fa'` | Hot colormap based on mean FA along track | White matter integrity |
| `'uniform'` | Single color for all tracks | Simplicity |
| `'region'` | Random color per parcellation region | Region comparison |

```matlab
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'whole', 'color_mode', 'fa');
```

### Track Filtering

Control which tracks are displayed based on their relationship to brain regions:

| Filter Mode | Description |
|---|---|
| `'pass_through'` | Track passes through the region (default) |
| `'start_in'` | Track starts in the region |
| `'end_in'` | Track ends in the region |
| `'connect_to'` | Track connects region to another specified region |

```matlab
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', ...
    'region', 5, ...
    'filter_mode', 'start_in');
```

See [REGION_VISUALIZATION_EXAMPLES.md](REGION_VISUALIZATION_EXAMPLES.md) for detailed filtering examples.

### Export

Export takes a *directory* plus a format, not a filename: `export_dir`,
`export_format` (`png`, `pdf`, `eps`, `fig`) and `export_dpi`. Subdirectories are
created automatically, and figures are not displayed while exporting unless
`silent_export` is set false.

```matlab
% PNG (raster, for presentations)
visualizeTractography('tracks.mat', 'nim.mat', 'export_dir', 'figures/');

% PDF (vector, for publications)
visualizeTractography('tracks.mat', 'nim.mat', ...
    'export_dir', 'figures/', 'export_format', 'pdf');

% MATLAB figure, for later editing
visualizeTractography('tracks.mat', 'nim.mat', ...
    'export_dir', 'figures/', 'export_format', 'fig');
```

---

## 4. Orthogonal Slice Viewer (`visualizeTractographySlices`)

**File**: `src/nim_visualization/visualizeTractographySlices.m`

Displays tractography on three orthogonal anatomical slices: sagittal, coronal
and axial. The three slice positions are **required positional arguments**, in
`x, y, z` order.

### Basic Usage

```matlab
% x = sagittal, y = coronal, z = axial
visualizeTractographySlices('tracks.mat', 'nim.mat', 60, 50, 45);
```

### Parameters

```matlab
visualizeTractographySlices('tracks.mat', 'nim.mat', 60, 50, 45, ...
    'tolerance', 2.0, ...          % slice thickness in voxels for track inclusion
    'regions', [5, 10], ...        % filter by parcellation region
    'region_filter', 'start_in', ...% 'pass_through' | 'start_in' | 'end_in'
    'min_overlap', 0.1, ...        % minimum track-region overlap
    'color_mode', 'direction', ... % 'direction' | 'uniform'
    'show_anatomy', true, ...      % FA background
    'alpha', 0.6, ...              % FA background transparency
    'show_crosshairs', true, ...   % draw the slice intersections
    'save', 'figures/slices.png'); % write a PNG instead of only displaying
```

### Features

- **Three views**: sagittal, coronal and axial planes at the requested positions
- **Crosshairs**: the slice intersections drawn on each panel
- **FA background**: anatomy overlay showing tissue structure
- **Region filtering**: display only tracks related to given parcellation regions
- **PNG export**: via `save`

### Performance

Each render recomputes track–slice intersections, so stepping through many slice
positions interactively is slow. For that, pre-compute a cache and use the
distributed workflow below.

---

## 5. Fast Distributed Viewer Workflow

Pre-compute every slice image once on a server, then page through them locally with a lightweight Python GUI. Navigation cost drops to an image load, and the local machine needs no MATLAB.

### Step 1: Generate Cache (Server, MATLAB)

Options are passed as a **struct**, not name-value pairs.

```matlab
% Basic generation
addpath('src/nim_visualization');
generateSlices('tracks.mat', 'nim.mat', '/export/slices');

% With custom settings
opts.image_resolution = [1024, 768];
opts.image_format     = 'png';       % 'png' | 'jpg'
opts.tolerance        = 2;           % slice thickness in voxels
opts.regions          = [];          % [] = all regions
generateSlices('tracks.mat', 'nim.mat', '/export/slices', opts);
```

The cache is keyed by dataset and by rendering parameters, so several parameter
sets can coexist under one output directory:

```
/export/slices/
└── datasets/
    └── <dataset_hash>/
        ├── metadata.json
        └── parameters/
            └── <param_hash>/
                ├── config.json
                ├── axial/*.png
                ├── sagittal/*.png
                └── coronal/*.png
```

Generation time depends on brain dimensions and track count. Workers are used
when the Parallel Computing Toolbox is available (`opts.parallel_workers`).

### Step 2: Transfer (rsync/scp)

```bash
# Transfer cache to local machine
rsync -avz server:/export/slices/ ~/local/slices/

# Or use scp
scp -r server:/export/slices/ ~/local/slices/
```

### Step 3: View (Local, Python only)

```bash
# Using the launcher script
./bin/viewSlices.sh ~/local/slices/

# Or directly
python scripts/FastTractographyViewer.py ~/local/slices/
```

### Python Viewer Features

- **Three synchronized panels**: Axial, sagittal, coronal views
- **Keyboard navigation**: Arrow keys to scroll through slices
- **Fast transitions**: slices are pre-rendered, so paging is an image load
- **Image caching**: In-memory cache for smooth scrolling
- **Performance monitoring**: Real-time display of load times
- **Export**: Save current view or batch export

### MATLAB Bridge

Alternatively, launch the Python viewer directly from MATLAB:

```matlab
launchFastViewer('tracks.mat', 'nim.mat');
% Automatically generates cache if needed, then launches Python GUI
```

See [DISTRIBUTED_WORKFLOW.md](DISTRIBUTED_WORKFLOW.md) for complete setup instructions.

---

## 6. Connectivity Matrix Visualization (`nim_plot_connectivity_matrix`)

**File**: `src/nim_tractography/nim_plot_connectivity_matrix.m`

Computes and displays the structural connectivity matrix between parcellation regions.

```matlab
% Compute and visualize
data_t = load('tracks.mat');
data_n = load('nim.mat');
C = nim_plot_connectivity_matrix(data_t.tracks, data_n.nim);
```

### Output

Creates a 4-panel figure:

1. **Heatmap**: Region-by-region connectivity strength
2. **Histogram**: Distribution of connection strengths
3. **Node strengths**: Bar chart showing total connectivity per region
4. **Summary statistics**: Track counts, region counts

### Options

```matlab
options.min_track_length = 15;   % only count tracks with 15+ points (default 10)
options.normalize = true;        % normalise to [0, 1]
options.symmetric = true;        % C(i,j) = C(j,i)
C = nim_plot_connectivity_matrix(tracks, nim, options);
```

---

## 7. Vector Field Visualization (`nim_plot_vector_field`)

**File**: `src/nim_tractography/nim_plot_vector_field.m`

Displays the primary eigenvector direction field on a 2D anatomical slice. Useful for verifying DTI direction estimates and understanding local fiber architecture.

Options are passed as a struct.

```matlab
nim_plot_vector_field(nim);                    % default: axial, middle slice

opts.axis_view = 'sagittal';
nim_plot_vector_field(nim, opts);

opts = struct('slice', 30, 'downsample', 3);   % specific slice, sparser quiver
nim_plot_vector_field(nim, opts);
```

Vectors are masked at FA > 0.2, so only white matter is drawn, and overlaid on an
FA background image.

---

## 8. Bundle Rendering (`nim_plot_bundles`, `nim_plot_vs_groundtruth`)

**Files**: `src/nim_visualization/nim_plot_bundles.m`,
`src/nim_visualization/nim_plot_vs_groundtruth.m`

These two render *streamline sets* rather than a whole tractogram, in a shared
voxel frame with fixed limits across panels, and they are what produced the
bundle figures in [ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md).

### `nim_plot_bundles`

Takes a struct array of items, one per streamline set, and draws each across four
views (sagittal, coronal, axial, oblique) by default.

```matlab
items(1).tracks = gt_tracks;              % cell array of Nx3 voxel polylines
items(1).name   = 'ground truth';
items(1).color  = [0.62 0.63 0.66];       % flat grey
items(1).alpha  = 0.2;

items(2).tracks = our_tracks;
items(2).name   = 'ours';
items(2).color  = [];                     % [] = colour each segment by direction
items(2).clip   = bundle_mask;            % optional: trim to a corridor

opts = struct('title', 'UF right', 'save', 'docs/img/bundle_uf_right', 'dpi', 150);
fig = nim_plot_bundles(items, opts);
```

Per item: `tracks`, `name`, `color` (`[]` for direction encoding — |dx| red, |dy|
green, |dz| blue, applied per segment), `alpha`, `width`, and `clip`. Per figure:
`views`, `title`, `subtitle`, `save` (a path stem; `.png` is appended), `dpi`,
`maxper` (streamlines drawn per item, default 900).

!!! note "`clip` trims, it does not reject"
    With `clip` set, each streamline is reduced to its longest contiguous run
    inside the mask, and the fraction of length that fell outside is returned in
    `items(i).frac_outside`. This is deliberately *not* the scorer's rule, which
    discards a streamline entirely if any point leaves the corridor — that rule
    makes an overshooting streamline vanish, so the figure shows a bundle
    stopping short when it in fact reached too far. Rejection lives in
    `nim_filter_tracks_roi` and is what scoring uses.

### `nim_plot_vs_groundtruth`

Overlays our streamlines (orange) on an ISMRM ground-truth bundle (grey) across
the same four views.

```matlab
opts = struct('max_gt', 700, 'max_ours', 700, ...
              'title', 'UF right', 'save', 'figures/uf_right');
nim_plot_vs_groundtruth(run_dir_or_tracks_mat, gt_trk, ref_nii, opts);
```

The ground truth arrives on the 1 mm 180×216×180 scoring grid and is mapped
through world RAS into the DWI's 2 mm 90×108×90 voxel space by `nim_read_trk`.
Skipping that step draws the ground truth at half scale in a corner, which looks
like a tracking failure and is a units bug.

### Figures in the docs

`docs/img/` holds the rendered outputs used across this documentation:
`tractogram_whole_brain.png`, `tractogram_bundles.png`, and one `bundle_*.png`
per bundle (`bundle_bps_right`, `bundle_cc_u_shaped`, `bundle_cingulum_right`,
`bundle_ilf_right`, `bundle_slf_right`, `bundle_uf_right`).

---

## 9. Presentation and Publication Figures

### Method Comparison Visualizations

The `src/nim_presentation/` module provides scripts for generating publication-quality comparison figures:

```matlab
% Compare integration methods (Euler vs RK2 vs RK4 vs RKF45)
visualize_integration_methods;

% Compare interpolation methods (none vs trilinear vs cubic)
visualize_interpolation_methods;

% Example tractography visualization
visualize_tractography_example;

% Slice-based comparison
visualize_tractography_slice;
```

### Export Settings for Publications

```matlab
% High-resolution PNG
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'whole', ...
    'max_tracts', 20000, ...
    'export_dir', 'figures/', 'export_dpi', 600);

% Vector format for LaTeX
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', 'region', 5, ...
    'export_dir', 'figures/', 'export_format', 'eps');

% MATLAB figure for further customisation
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'whole', ...
    'export_dir', 'figures/', 'export_format', 'fig');
% Then open and adjust:
openfig('figures/tractography_whole.fig');
set(gcf, 'PaperSize', [8 6]);
print('-dpdf', '-r300', 'figures/fig1_final.pdf');
```

### Recommended Figure Sizes

| Journal Format | Width | Height | DPI |
|---|---|---|---|
| Single column | 3.5 in | 2.5 in | 300 |
| Double column | 7.0 in | 5.0 in | 300 |
| Full page | 7.0 in | 9.0 in | 300 |
| Poster | 12.0 in | 9.0 in | 150 |

---

## CLI Export with Shell Script

Export tractography figures from the command line without an interactive MATLAB
session. The script wraps `visualizeTractographyAngles()`, which renders the same
track set from eight standard anatomical viewpoints, and runs MATLAB in the
background.

### Usage

```bash
# From a run directory (auto-detects tracks + nim; figures go to <run_dir>/figures/)
./bin/run_visualization.sh <run_dir> [format] [region] [dpi]

# With explicit file paths (output_dir is required in this form)
./bin/run_visualization.sh <tracks_file> <nim_file> <output_dir> [format] [region] [dpi]
```

There is no mode argument: the script always generates all eight views
(superior, inferior, anterior, posterior, left, right, oblique_left,
oblique_right).

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `output_dir` | `<run_dir>/figures/` in run-dir form; required in explicit form | where to save exported images |
| `format` | `png` | `png`, `pdf`, `eps` |
| `region` | — | region ID(s) to filter by, e.g. `5` or `5,10,15`; omitted means all tracks |
| `dpi` | `300` | export resolution |

### Examples

```bash
# Whole brain, all eight views
./bin/run_visualization.sh hinec_runs/run_<TIMESTAMP>_standard_dti/

# PDF at high DPI
./bin/run_visualization.sh hinec_runs/run_<TIMESTAMP>_standard_dti/ pdf '' 600

# Specific regions
./bin/run_visualization.sh hinec_runs/run_<TIMESTAMP>_standard_dti/ png '5,10,15'

# Explicit files and output directory
./bin/run_visualization.sh tracks.mat nim.mat publication_figures/ pdf '' 600
```

### Output

Files land flat in the output directory, one per view:

```
<output_dir>/
├── tractography_superior.png
├── tractography_inferior.png
├── tractography_anterior.png
├── ...
└── tractography_oblique_right.png
```

With a region filter the names carry a region suffix. In run-directory form logs
go to `<run_dir>/logs/visualization_<timestamp>.log`; otherwise
`visualization_<timestamp>.log` in the working directory. Monitor with
`tail -f`.

---

## Cross-References

- Distributed workflow setup: [DISTRIBUTED_WORKFLOW.md](DISTRIBUTED_WORKFLOW.md)
- Region visualization examples: [REGION_VISUALIZATION_EXAMPLES.md](REGION_VISUALIZATION_EXAMPLES.md)
- Complete API reference: [API_REFERENCE.md](API_REFERENCE.md)
- User guide: [USER_GUIDE.md](USER_GUIDE.md)
- Architecture overview: [ARCHITECTURE.md](ARCHITECTURE.md)
