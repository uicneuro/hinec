# Tractography Slice Viewer

A simple two-component system for viewing tractography slices:

1. **MATLAB Command-Line Function**: `visualizeTractographySlices.m`
2. **Python GUI**: `tractography_slice_gui.py`

## Quick Start

### Option 1: MATLAB Command Line (Direct)

```matlab
% Basic usage - display slice at position x=64, y=64, z=32
visualizeTractographySlices('tractography_results/tracks_standard.mat', 'sample_parcellated.mat', 64, 64, 32)

% Save to PNG file
visualizeTractographySlices('tractography_results/tracks_standard.mat', 'sample_parcellated.mat', 64, 64, 32, 'save', 'my_slice.png')

% With additional options
visualizeTractographySlices('tracks.mat', 'nim.mat', 80, 90, 45, ...
    'tolerance', 3, ...
    'show_crosshairs', false, ...
    'show_anatomy', true, ...
    'color_mode', 'direction', ...
    'alpha', 0.7, ...
    'save', 'output.png')
```

### Option 2: Python GUI (Interactive)

```bash
# Run the Python GUI
python tractography_slice_gui.py
```

The GUI provides:
- File selection for tracks and NIM files
- Interactive sliders for X, Y, Z slice positions
- Display options (crosshairs, anatomy background, color mode, transparency)
- Preview button (display only)
- Render & Save button (save to PNG)

## Parameters

### Required Arguments
- `tracks_file`: Path to tracks .mat file (e.g., `'tractography_results/tracks_standard.mat'`)
- `nim_file`: Path to nim structure .mat file (e.g., `'sample_parcellated.mat'`)
- `x`: Sagittal slice position (X coordinate)
- `y`: Coronal slice position (Y coordinate)
- `z`: Axial slice position (Z coordinate)

### Optional Parameters
- `'save'`: Output PNG filename (if not specified, just displays figure)
- `'tolerance'`: Slice thickness in voxels (default: 2)
- `'show_crosshairs'`: Show slice intersections (default: true)
- `'show_anatomy'`: Show FA background (default: true)
- `'color_mode'`: `'direction'` or `'uniform'` (default: `'direction'`)
- `'alpha'`: FA background transparency 0-1 (default: 0.6)

## Output

The viewer shows three orthogonal slice views:
1. **Axial (Z)**: Top-down view (X-Y plane)
2. **Sagittal (X)**: Side view (Y-Z plane)
3. **Coronal (Y)**: Front view (X-Z plane)
4. **Info Panel**: Track statistics and settings

## Performance

- Fast command-line execution
- Pre-computed track lookup tables for efficiency
- No laggy UI interactions
- Clear terminal feedback during processing

## Requirements

- MATLAB with required toolboxes
- Python 3.x with tkinter (for GUI)
- Processed tractography data from HINEC pipeline

## Architecture Benefits

✅ **Simple**: No complex interactive UI code
✅ **Fast**: Direct rendering without UI overhead
✅ **Flexible**: Command-line allows scripting and automation
✅ **User-Friendly**: Python GUI for interactive exploration
✅ **Reliable**: Clear separation between UI and processing logic
✅ **Extensible**: Easy to add new features to either component

## Example Workflow

1. Run tractography: `runTractography('sample_parcellated.mat')`
2. Open Python GUI: `python tractography_slice_gui.py`
3. Set file paths to tracks and NIM files
4. Adjust slice positions with sliders
5. Click "Preview" to display, or "Render & Save" to export PNG

## Troubleshooting

- Ensure MATLAB is in your system PATH for Python integration
- Check that tracks and NIM files exist and have correct structure
- Slice coordinates are automatically clamped to valid ranges
- High resolution PNG output available with `-r300` flag