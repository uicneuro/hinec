# Region-Specific Tractography Visualization

This guide shows how to visualize tractography tracks for specific brain regions using the new region visualization functions.

## Quick Start

First, make sure you have processed data with tractography and parcellation:

```matlab
% Run the full pipeline if you haven't already
main('nifti_sample/sample', 'sample_parcellated.mat');
runTractography('sample_parcellated.mat');
```

## Basic Usage

### Step 1: List Available Regions
```matlab
% Use your nim file from main()
listBrainRegions('sample_parcellated.mat');
```

### Step 2: Find Your Tracks File
After running `runTractography()`, your tracks will be in `tractography_results/`. 

Use `ls tractography_results/` to find your specific tracks file name, e.g.:
- `tracks_2024-01-01_12-00-00.mat`

### Step 3: Visualize a Specific Region
```matlab
% Basic usage - ALL arguments required
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat');

% With different filter modes
visualizeTractographyRegion(12, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'start_in');
visualizeTractographyRegion(8, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'end_in');
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'connect_to');
```

## Advanced Options

### Custom Visualization Settings
```matlab
% Advanced example with custom settings
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'pass_through', ...
    'color_mode', 'fa', ...           % Color by FA values
    'max_tracks', 500, ...            % Limit to 500 tracks for performance
    'show_region', true, ...          % Show region overlay
    'region_alpha', 0.5, ...          % Region transparency
    'min_overlap', 0.2);              % 20% minimum overlap
```

### Different Color Modes
```matlab
% Direction-based coloring (default)
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'pass_through', ...
    'color_mode', 'direction');

% FA-based coloring (hot colormap)
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'pass_through', ...
    'color_mode', 'fa');

% Single color for all tracks
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'pass_through', ...
    'color_mode', 'uniform');

% Region-based coloring
visualizeTractographyRegion(5, 'tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat', 'pass_through', ...
    'color_mode', 'region');
```

## Filter Modes Explained

1. **'pass_through'** (default): Shows tracks that pass through the region at any point
2. **'start_in'**: Shows only tracks that start within the region
3. **'end_in'**: Shows only tracks that end within the region  
4. **'connect_to'**: Shows tracks that connect the region to other regions

## Function Reference

### visualizeTractographyRegion()
Main function for region-specific visualization.

**Syntax:**
```matlab
visualizeTractographyRegion(region_id, tracks_file, nim_file)
visualizeTractographyRegion(region_id, tracks_file, nim_file, filter_mode)
visualizeTractographyRegion(region_id, tracks_file, nim_file, filter_mode, Name, Value, ...)
```

**Parameters:**
- `region_id`: Integer index of brain region (from parcellation mask) - **REQUIRED**
- `tracks_file`: Path to tracks .mat file - **REQUIRED**
- `nim_file`: Path to .mat file containing nim structure with parcellation - **REQUIRED**
- `filter_mode`: How to filter tracks ('pass_through', 'start_in', 'end_in', 'connect_to')
- `min_overlap`: Minimum fraction of track in region (0-1)
- `show_region`: Show region as 3D overlay (true/false)
- `color_mode`: Track coloring ('direction', 'fa', 'uniform', 'region')
- `max_tracks`: Maximum tracks to display (default: unlimited)
- `region_alpha`: Region overlay transparency (0-1)

### listBrainRegions()
Lists all available brain regions with statistics.

**Syntax:**
```matlab
listBrainRegions(nim_file)
region_info = listBrainRegions(nim_file)
```

**Arguments:**
- `nim_file`: Path to .mat file containing nim structure with parcellation_mask - **REQUIRED**

**Returns:**
- `region_info`: Struct array with region ID, name, voxel count, percentage

### visualizeTractography()
Original visualization function (also requires all arguments).

**Syntax:**
```matlab
visualizeTractography(tracks_file, nim_file)
```

**Arguments:**
- `tracks_file`: Path to tracks .mat file - **REQUIRED**
- `nim_file`: Path to nim .mat file - **REQUIRED**

## Tips and Best Practices

### Finding Your Files
```matlab
% List available tracks files
ls tractography_results/

% List .mat files in current directory
ls *.mat
```

### Performance Tips
- Use `max_tracks` parameter to limit display for better performance  
- Start with smaller regions (higher region IDs) for faster visualization
- Use `'uniform'` color mode for fastest rendering

### Interpretation Tips
- **Direction coloring**: Red=Left-Right, Green=Anterior-Posterior, Blue=Superior-Inferior
- **FA coloring**: Bright colors indicate high anisotropy (well-organized fibers)
- **Region overlay**: Red transparent surface shows the selected region
- Track statistics are printed to console for quantitative analysis

### Common Use Cases

1. **Explore region connectivity**:
   ```matlab
   visualizeTractographyRegion(region_id, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat', 'connect_to');
   ```

2. **Find outgoing pathways**:
   ```matlab
   visualizeTractographyRegion(region_id, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat', 'start_in');
   ```

3. **Find incoming pathways**:
   ```matlab
   visualizeTractographyRegion(region_id, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat', 'end_in');
   ```

4. **Study tract organization**:
   ```matlab
   visualizeTractographyRegion(region_id, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat', 'pass_through', ...
       'color_mode', 'fa');
   ```

## Troubleshooting

### Common Errors

**"visualizeTractographyRegion requires 3 arguments"**
- You must specify ALL three arguments: region_id, tracks_file, nim_file
- Example: `visualizeTractographyRegion(5, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat')`

**"Region X does not exist"**
- Run `listBrainRegions('sample_parcellated.mat')` to see available regions
- Check that parcellation was successful

**"Tracks file not found"**
- Check the exact filename in `tractography_results/` directory
- Run `runTractography('sample_parcellated.mat')` first to generate tracks

**"Nim file not found"**
- Check that your .mat file exists and contains nim structure
- Run `main()` first to generate nim data with parcellation

**"No tracks found for region X"**
- Try lowering `min_overlap` parameter
- Use different filter mode (e.g., 'pass_through' instead of 'start_in')
- Check if region is too small or in non-white-matter area

**Poor performance/visualization**
- Reduce `max_tracks` parameter
- Use `'uniform'` color mode
- Set `show_region` to false

## Example Workflow

```matlab
% 1. List available regions
listBrainRegions('sample_parcellated.mat');

% 2. Find your tracks file
ls tractography_results/

% 3. Pick an interesting region (e.g., region 15) and your actual tracks file
region_id = 15;
tracks_file = 'tractography_results/tracks_2024-01-01_12-00-00.mat';  % Use your actual file
nim_file = 'sample_parcellated.mat';

% 4. Explore different aspects of this region
visualizeTractographyRegion(region_id, tracks_file, nim_file, 'pass_through');  % All related tracks
visualizeTractographyRegion(region_id, tracks_file, nim_file, 'start_in');      % Outgoing tracks  
visualizeTractographyRegion(region_id, tracks_file, nim_file, 'end_in');        % Incoming tracks
visualizeTractographyRegion(region_id, tracks_file, nim_file, 'connect_to');    % Inter-region connections

% 5. Study with different color modes
visualizeTractographyRegion(region_id, tracks_file, nim_file, 'pass_through', ...
    'color_mode', 'fa');
visualizeTractographyRegion(region_id, tracks_file, nim_file, 'pass_through', ...
    'color_mode', 'direction');
```

This gives you a comprehensive view of how the selected brain region connects to the rest of the brain through white matter pathways.

## Key Differences from Original Functions

- **No auto-detection**: ALL file arguments must be specified explicitly
- **Clear error messages**: Functions will tell you exactly what arguments are missing
- **Consistent with visualizeTractography**: Both functions now require all arguments
- **No guessing**: You always know exactly which files are being used