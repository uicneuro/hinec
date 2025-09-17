# Tractography Visualization Documentation

## Overview

The **`visualizeTractography.m`** function is a comprehensive, unified solution for tractography visualization in the HINEC pipeline. This single function consolidates all previous visualization functionality into one powerful interface that supports multiple visualization modes, advanced filtering, and automatic image export.

## Key Features

- ✅ **4 Visualization Modes**: Whole brain, single region, grid layout, sequential
- ✅ **Advanced Track Filtering**: Multiple filtering strategies with overlap control
- ✅ **Multiple Color Schemes**: Direction-based, FA-based, uniform, and region-based coloring
- ✅ **Automatic Image Export**: PNG, PDF, EPS, FIG formats with organized directory structure
- ✅ **Performance Optimization**: Configurable track limits and smart rendering
- ✅ **Legacy Compatibility**: Supports all parameters from previous functions

## Architecture

The function is organized into **12 major sections**:

1. **Input Validation & Parsing** - Comprehensive parameter validation
2. **Data Loading & Validation** - Robust data loading with error handling
3. **Visualization Modes** - Four distinct visualization implementations
4. **Core Plotting Functions** - Shared plotting utilities
5. **Track Filtering** - Advanced filtering algorithms
6. **Utility Functions** - Helper functions for names, grids, layouts
7. **Statistical Analysis** - Comprehensive statistics and reporting
8. **Image Export** - Professional image export with metadata logging

## Visualization Modes

### 1. Whole Brain Mode (`mode='whole'`)

Comprehensive overview of all tractography tracks with detailed analytics.

**Features:**
- Main 3D visualization with anatomical background
- Track length histogram
- Seed point distribution analysis
- Performance optimization for large datasets
- Comprehensive statistics output

**Usage:**
```matlab
% Basic whole brain visualization
visualizeTractography('tracks.mat', 'sample_parcellated.mat')

% With custom settings
visualizeTractography('tracks.mat', 'sample_parcellated.mat', ...
    'color_mode', 'fa', 'max_tracks', 5000, 'export_dir', 'figures/')
```

**Output:**
- 2x3 subplot layout with main 3D view
- Performance statistics
- Track summary information

### 2. Single Region Mode (`mode='region'` or `'region', ID`)

Detailed analysis of tractography for a specific brain region.

**Features:**
- Region-specific track filtering
- Region overlay visualization
- Detailed region statistics
- Multiple filtering strategies
- Comprehensive track analysis

**Usage:**
```matlab
% Single region visualization
visualizeTractography('tracks.mat', 'sample_parcellated.mat', 'region', 5)

% With advanced filtering
visualizeTractography('tracks.mat', 'sample_parcellated.mat', ...
    'region', 5, 'filter_mode', 'start_in', 'min_overlap', 0.2)
```

**Output:**
- Full-screen 3D visualization
- Region information display
- Track filtering statistics
- Regional connectivity analysis

### 3. Grid Layout Mode (`mode='grid'`)

Overview of all brain regions in an organized grid layout.

**Features:**
- Automatic grid layout optimization
- Individual region subplots
- Comprehensive region statistics
- ALL tracks display capability
- Publication-ready layouts

**Usage:**
```matlab
% Grid layout for all regions
visualizeTractography('tracks.mat', 'sample_parcellated.mat', 'mode', 'grid')

% Custom grid configuration
visualizeTractography('tracks.mat', 'sample_parcellated.mat', ...
    'mode', 'grid', 'grid_cols', 4, 'subplot_size', [300, 250])
```

**Output:**
- Optimized grid layout (2x2 to 6x5 depending on region count)
- Individual region titles and track counts
- Comprehensive statistics table
- Overall summary information

### 4. Sequential Mode (`mode='sequential'`)

Interactive sequential visualization of all regions.

**Features:**
- Step-through individual regions
- User-controlled progression
- Individual exports for each region
- Detailed per-region analysis
- Interactive pause functionality

**Usage:**
```matlab
% Sequential visualization
visualizeTractography('tracks.mat', 'sample_parcellated.mat', 'mode', 'sequential')

% With automatic export
visualizeTractography('tracks.mat', 'sample_parcellated.mat', ...
    'mode', 'sequential', 'export_dir', 'sequential_output/')
```

**Output:**
- Individual full-screen visualizations
- User interaction prompts
- Optional individual image exports
- Progress tracking

## Parameter Reference

### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | string | `'whole'` | Visualization mode: `'whole'`, `'region'`, `'grid'`, `'sequential'` |
| `region` | integer | `[]` | Region ID for single region mode (alternative to `mode='region'`) |

### Export Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `export_dir` | string | `''` | Directory for image export (auto-creates subdirectories) |
| `export_format` | string | `'png'` | Image format: `'png'`, `'pdf'`, `'eps'`, `'fig'` |
| `export_dpi` | integer | `300` | Export resolution in DPI |
| `silent_export` | logical | `true` | Hide figures and skip pauses when exporting |

### Track Filtering Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `filter_mode` | string | `'pass_through'` | Track filtering: `'pass_through'`, `'start_in'`, `'end_in'`, `'connect_to'` |
| `min_overlap` | double | `0.1` | Minimum region overlap ratio (0-1) |
| `max_tracks` | integer | varies | Maximum tracks to display (2000 for whole, inf for others) |

### Visualization Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `color_mode` | string | `'direction'` | Track coloring: `'direction'`, `'fa'`, `'uniform'`, `'region'` |
| `show_anatomy` | logical | `true` | Show FA background slices |
| `show_region` | logical | `true` | Show region boundary overlays |
| `line_width` | double | `1.2` | Track line width |
| `region_alpha` | double | `0.3` | Region overlay transparency (0-1) |

### Grid Layout Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `grid_cols` | integer | `0` | Grid columns (0 = auto-calculate) |
| `subplot_size` | array | `[400, 350]` | Subplot dimensions [width, height] |
| `show_all_tracks` | logical | `true` | Force display of ALL tracks in grid mode |

## Filtering Strategies

### `'pass_through'` (Default)
Includes tracks that pass through the region at any point.
- **Use case**: General connectivity analysis
- **Performance**: Fastest filtering
- **Typical tracks**: 50-200% of other modes

### `'start_in'`
Includes only tracks that start within the region.
- **Use case**: Efferent connectivity analysis
- **Performance**: Fast filtering
- **Typical tracks**: 20-50% of pass_through

### `'end_in'`
Includes only tracks that end within the region.
- **Use case**: Afferent connectivity analysis
- **Performance**: Fast filtering
- **Typical tracks**: 20-50% of pass_through

### `'connect_to'`
Includes tracks that connect the region to other regions.
- **Use case**: Inter-regional connectivity
- **Performance**: Moderate filtering
- **Typical tracks**: 30-70% of pass_through

## Color Schemes

### `'direction'` (Default)
RGB colors based on track direction vectors.
- **Red**: Left-Right movement
- **Green**: Anterior-Posterior movement
- **Blue**: Superior-Inferior movement
- **Best for**: Understanding fiber orientation

### `'fa'`
Hot colormap based on Fractional Anisotropy values.
- **Dark red**: Low FA values
- **Yellow/White**: High FA values
- **Best for**: Tissue integrity assessment

### `'uniform'`
Single blue color for all tracks.
- **Color**: Light blue `[0.2, 0.6, 0.8]`
- **Best for**: Clean presentations, focus on structure

### `'region'`
Distinct colors based on track's primary region.
- **Colors**: Consistent per region (seeded random)
- **Best for**: Regional connectivity analysis

## Export System

### Directory Structure
```
export_dir/
├── whole/          # Whole brain visualizations
├── region/         # Single region visualizations
├── grid/           # Grid layout visualizations
├── sequential/     # Sequential mode visualizations
└── metadata/       # Export logs and metadata
    └── export_log.txt
```

### Naming Convention
- **Whole brain**: `tractography_whole.png`
- **Single region**: `tractography_region-05.png` (where 05 is the region ID)
- **Grid layout**: `tractography_grid.png`
- **Sequential**: `sequential_region-05.png` (individual regions)

### Export Formats
- **PNG**: High-quality raster (default, 300 DPI)
- **PDF**: Vector format for publications
- **EPS**: Vector format for journals
- **FIG**: MATLAB figure format (editable)

## Migration from Legacy Functions

### From `visualizeTractography.m` (original)
```matlab
% OLD
visualizeTractography('tracks.mat', 'data.mat')

% NEW (same)
visualizeTractography('tracks.mat', 'data.mat')
```

### From `visualizeTractographyRegion.m`
```matlab
% OLD
visualizeTractographyRegion(5, 'tracks.mat', 'data.mat', 'start_in')

% NEW
visualizeTractography('tracks.mat', 'data.mat', 'region', 5, 'filter_mode', 'start_in')
```

### From `visualizeTractographyAllRegions.m`
```matlab
% OLD
visualizeTractographyAllRegions('tracks.mat', 'data.mat', ...
    'filter_mode', 'pass_through', 'show_all_tracks', true)

% NEW
visualizeTractography('tracks.mat', 'data.mat', 'mode', 'grid', ...
    'filter_mode', 'pass_through', 'show_all_tracks', true)
```

## Example Workflows

### Basic Analysis Workflow
```matlab
% 1. Overview
visualizeTractography('tracks.mat', 'data.mat', 'export_dir', 'analysis_2024/')

% 2. Grid comparison
visualizeTractography('tracks.mat', 'data.mat', 'mode', 'grid', 'export_dir', 'analysis_2024/')

% 3. Detailed region analysis
visualizeTractography('tracks.mat', 'data.mat', 'region', 5, 'export_dir', 'analysis_2024/')
```

### Publication Workflow
```matlab
% High-quality exports for publication
visualizeTractography('tracks.mat', 'data.mat', ...
    'mode', 'grid', 'export_dir', 'publication/', ...
    'export_format', 'pdf', 'export_dpi', 600, ...
    'color_mode', 'direction', 'show_anatomy', true)
```

### Connectivity Analysis Workflow
```matlab
% Compare different filtering strategies
for filter_type = {'pass_through', 'start_in', 'end_in', 'connect_to'}
    visualizeTractography('tracks.mat', 'data.mat', ...
        'region', 5, 'filter_mode', filter_type{1}, ...
        'export_dir', ['connectivity_' filter_type{1}])
end
```

### Batch Processing Workflow
```matlab
% Silent batch processing for multiple datasets
datasets = {'subj001', 'subj002', 'subj003'};
for i = 1:length(datasets)
    subj = datasets{i};

    % Process all modes without user interaction
    visualizeTractography([subj '_tracks.mat'], [subj '_data.mat'], ...
        'mode', 'whole', 'export_dir', ['batch_output/' subj '/']);

    visualizeTractography([subj '_tracks.mat'], [subj '_data.mat'], ...
        'mode', 'grid', 'export_dir', ['batch_output/' subj '/']);

    fprintf('Completed processing for %s\n', subj);
end
% All figures exported silently - no user interaction required!
```

## Performance Tips

- Use `max_tracks` for large datasets in whole brain mode
- Enable export only when needed
- Use `'uniform'` color for fastest rendering
- Consider `filter_mode` impact on track count
- **Silent Export**: When `export_dir` is specified, figures are hidden and pauses are skipped for faster batch processing

## Troubleshooting

### Common Issues

**Slow Performance**
- Reduce `max_tracks` parameter
- Use `'uniform'` color mode
- Disable `show_anatomy` for faster rendering

**Memory Issues**
- Process regions individually instead of grid mode
- Reduce track dataset size before visualization
- Use sequential mode for very large datasets

**Export Issues**
- Check disk space in export directory
- Verify write permissions
- Use PNG format if PDF export fails

---

*This documentation covers the consolidated `visualizeTractography.m` function that replaces all previous visualization functions in the HINEC pipeline.*