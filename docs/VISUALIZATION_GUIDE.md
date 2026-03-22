# Complete Visualization Reference

Unified guide to all visualization capabilities in HINEC. Covers DTI visualization, 3D tractography, interactive slice viewing, fast distributed viewing, connectivity matrices, and publication-quality exports.

---

## 1. Visualization Overview

HINEC provides multiple visualization tools optimized for different tasks:

| Tool | Type | Speed | Use Case |
|---|---|---|---|
| `nim_plot` | DTI eigenvectors | Instant | Quick DTI quality check |
| `visualizeTractography` | 3D tractography | 5-30s | Comprehensive 3D exploration |
| `visualizeTractographySlices` | 2D slices | 5-30s/slice | Detailed slice-by-slice inspection |
| `generateSlices` + `FastTractographyViewer.py` | Pre-computed 2D | <100ms/slice | Fast daily navigation |
| `nim_plot_tractography` | Basic 3D | 2-10s | Quick track visualization |
| `nim_plot_connectivity_matrix` | 2D heatmap | 1-5s | Connectivity analysis |
| `nim_plot_vector_field` | 2D quiver | 1-2s | Direction field inspection |

### Decision Tree: Which Tool to Use

```
What do you need?
├── Quick DTI quality check → nim_plot
├── 3D whole-brain tractography → visualizeTractography('mode', 'whole')
├── Specific brain region analysis → visualizeTractography('mode', 'region')
├── All regions at once → visualizeTractography('mode', 'grid')
├── Step through regions one by one → visualizeTractography('mode', 'sequential')
├── Slice-by-slice inspection (interactive) → visualizeTractographySlices
├── Fast daily slice navigation → generateSlices + FastTractographyViewer.py
├── Connectivity between regions → nim_plot_connectivity_matrix
├── Direction field on a slice → nim_plot_vector_field
└── Publication figures → visualizeTractography with 'export' option
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
nim_plot(nim, 'mode', 'parcel', 'region_id', 5);
```
Shows eigenvectors within a specific parcellation region.

**All parcels** (separate figures):
```matlab
nim_plot(nim, 'mode', 'parcels');
```
Creates one figure per parcellation region. Useful for surveying all regions.

### Performance

For large datasets, use downsampling:
```matlab
nim_plot(nim, 'mode', 'single', 'downsample', 3);  % Show every 3rd voxel
```

---

## 3. 3D Tractography Visualization (`visualizeTractography`)

**File**: `src/nim_visualization/visualizeTractography.m` (1609 lines)

The primary tool for comprehensive 3D tractography visualization. Supports 4 modes, multiple color schemes, track filtering, run directory auto-detection, and export.

### Input Flexibility

```matlab
% Direct file paths
visualizeTractography('path/to/tracks.mat', 'path/to/nim.mat');

% From run directory (auto-detects latest tracks)
visualizeTractography('hinec_runs/run_2025-01-15/', 'path/to/nim.mat');

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
    'max_tracks', 10000, ...
    'color_mode', 'direction');
```

### Mode: Single/Multi Region

Focused view of tracks associated with specific brain regions.

```matlab
% Single region
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', ...
    'regions', [5]);

% Multiple regions
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', ...
    'regions', [5, 10, 15]);
```

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
    'regions', [5], ...
    'filter_mode', 'start_in');
```

See [REGION_VISUALIZATION_EXAMPLES.md](REGION_VISUALIZATION_EXAMPLES.md) for detailed filtering examples.

### Export

```matlab
% PNG (raster, best for presentations)
visualizeTractography('tracks.mat', 'nim.mat', 'export', 'figures/brain.png');

% PDF (vector, best for publications)
visualizeTractography('tracks.mat', 'nim.mat', 'export', 'figures/brain.pdf');

% EPS (vector, for LaTeX)
visualizeTractography('tracks.mat', 'nim.mat', 'export', 'figures/brain.eps');

% MATLAB figure (for later editing)
visualizeTractography('tracks.mat', 'nim.mat', 'export', 'figures/brain.fig');
```

---

## 4. Interactive Slice Viewer (`visualizeTractographySlices`)

**File**: `src/nim_visualization/visualizeTractographySlices.m` (536 lines)

Displays tractography on three orthogonal anatomical slices: axial (top-down), sagittal (side), and coronal (front).

### Basic Usage

```matlab
visualizeTractographySlices('tracks.mat', 'nim.mat');
```

### Parameters

```matlab
visualizeTractographySlices('tracks.mat', 'nim.mat', ...
    'slice_axial', 45, ...       % Specific axial slice
    'slice_sagittal', 60, ...    % Specific sagittal slice
    'slice_coronal', 50, ...     % Specific coronal slice
    'tolerance', 2.0, ...        % Slice thickness (mm) for track inclusion
    'regions', [5, 10]);         % Filter by regions
```

### Features

- **Three synchronized views**: Axial, sagittal, and coronal planes
- **Crosshair synchronization**: Click in one view to update the others
- **FA background**: Anatomy overlay showing tissue structure
- **Region filtering**: Display only tracks through specific brain regions
- **Batch export**: Save all slices as PNG images

### Performance

Each slice update takes 5-30 seconds due to real-time computation of track-slice intersections. For faster navigation, use the distributed workflow below.

### Batch Export

```matlab
visualizeTractographySlices('tracks.mat', 'nim.mat', ...
    'export_dir', 'figures/slices/');
```

---

## 5. Fast Distributed Viewer Workflow

For sub-100ms slice navigation, pre-compute all slice images on a server and view them locally with a lightweight Python GUI.

### Step 1: Generate Cache (Server, MATLAB)

```matlab
% Basic generation
addpath('src/nim_visualization');
generateSlices('tracks.mat', 'nim.mat', '/export/slices');

% With custom settings
generateSlices('tracks.mat', 'nim.mat', '/export/slices', ...
    'resolution', [1024, 768], ...
    'format', 'png', ...
    'parallel', true);
```

This creates a cache directory structure:
```
/export/slices/
├── metadata.json          # Dataset info, parameters
├── axial/                 # Axial slice images
│   ├── slice_001.png
│   ├── slice_002.png
│   └── ...
├── sagittal/              # Sagittal slice images
│   └── ...
└── coronal/               # Coronal slice images
    └── ...
```

Generation time depends on brain dimensions and number of tracks. Parallel processing is used when available.

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
- **Sub-100ms transitions**: Pre-computed images load instantly
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
options.min_track_length = 15;   % Only count tracks with 15+ points
options.normalize = true;        % Normalize to [0, 1]
options.symmetric = true;        % C(i,j) = C(j,i)
C = nim_plot_connectivity_matrix(tracks, nim, options);
```

---

## 7. Vector Field Visualization (`nim_plot_vector_field`)

**File**: `src/nim_tractography/nim_plot_vector_field.m`

Displays the primary eigenvector direction field on a 2D anatomical slice. Useful for verifying DTI direction estimates and understanding local fiber architecture.

```matlab
nim_plot_vector_field(nim);                                  % Default: axial, middle slice
nim_plot_vector_field(nim, 'axis_view', 'sagittal');          % Sagittal view
nim_plot_vector_field(nim, 'slice', 30, 'downsample', 3);    % Specific slice, sparser
```

Vectors are masked by FA > 0.2 (only shown in white matter) and overlaid on an FA background image.

---

## 8. Presentation and Publication Figures

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

For publication-quality figures:

```matlab
% High-resolution PNG
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'whole', ...
    'max_tracks', 20000, ...
    'export', 'figures/fig1_whole_brain.png');

% Vector format for LaTeX
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'region', ...
    'regions', [5], ...
    'export', 'figures/fig2_region.eps');

% MATLAB figure for further customization
visualizeTractography('tracks.mat', 'nim.mat', ...
    'mode', 'whole', ...
    'export', 'figures/fig1.fig');
% Then open and customize:
openfig('figures/fig1.fig');
set(gcf, 'PaperSize', [8 6]);  % Adjust paper size
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

## Cross-References

- Distributed workflow setup: [DISTRIBUTED_WORKFLOW.md](DISTRIBUTED_WORKFLOW.md)
- Region visualization examples: [REGION_VISUALIZATION_EXAMPLES.md](REGION_VISUALIZATION_EXAMPLES.md)
- Complete API reference: [API_REFERENCE.md](API_REFERENCE.md)
- User guide: [USER_GUIDE.md](USER_GUIDE.md)
- Architecture overview: [ARCHITECTURE.md](ARCHITECTURE.md)
