# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HINEC (HIgh-order NEural Connectivity) is a MATLAB-based pipeline for processing and analyzing diffusion-weighted MRI (dMRI) data. The pipeline processes raw NIfTI diffusion MRI data through preprocessing, diffusion tensor calculation, fractional anisotropy computation, parcellation, and tractography analysis.

## Key Commands

### Main Execution
- `main('path/to/data', 'output.mat')` - Core pipeline execution
- `runhinec` - Complete processing workflow with visualization
- `runTractography('data.mat')` - Run tractography analysis on processed data

### Data Processing Functions
- `nim_read(imgpath)` - Load NIfTI diffusion data
- `nim_dt_spd(nim)` - Calculate diffusion tensors (SPD-constrained)
- `nim_fa(nim)` - Compute fractional anisotropy
- `nim_parcellation(nim, mask_file)` - Apply brain parcellation
- `nim_tractography_standard(nim, options)` - Standard tractography

### Visualization
- `nim_plotall(nim)` - Comprehensive visualization
- `nim_plotparcelall(nim)` - Parcellation-specific plots
- `nim_plotparcellation(nim)` - Parcellation mask visualization

### Tractography Visualization (nim_visualization/)
- `visualizeTractography(tracks_file, nim_file)` - Unified tractography viewer with multiple modes
- `visualizeTractographySlices(tracks_file, nim_file)` - Interactive 2D slice viewer
- `generateSlices(tracks_file, nim_file, output_dir)` - Server-side slice generation for distributed viewing
- `launchFastViewer(tracks_file, nim_file)` - MATLAB bridge to Python fast viewer
- `FastTractographyViewer.py` - Python GUI for instant slice navigation (no MATLAB required)

## Architecture

### Core Modules
- **nim_calculation/**: Diffusion tensor and FA calculations
- **nim_parcellation/**: Brain region segmentation
- **nim_preprocessing/**: Raw data preprocessing (FSL integration)
- **nim_tractography/**: Fiber tractography algorithms
- **nim_visualization/**: Tractography visualization and fast slice viewer
- **nim_plots/**: Visualization functions
- **nim_utils/**: Utility functions for data I/O and manipulation

### Data Flow
1. **Input**: Raw NIfTI files with bval/bvec (e.g., `sample_raw.nii.gz`)
2. **Preprocessing**: Motion/eddy correction via FSL (if raw data provided)
3. **Processing**: Diffusion tensor → Eigenvalues/vectors → FA calculation
4. **Parcellation**: Atlas-based brain region segmentation
5. **Output**: Processed `.mat` file with all computed metrics
6. **Analysis**: Tractography and visualization

### File Naming Conventions
- Raw data: `{name}_raw.nii.gz`, `{name}.bval`, `{name}.bvec`
- Acquisition params: `{name}_acqp.txt`, `{name}_index.txt`
- Processed: `{name}.nii.gz` (after preprocessing)
- Output: User-specified `.mat` file

## Dependencies

### MATLAB Toolboxes (Required)
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Tools for NIfTI and ANALYZE image

### External Software
- **SPM12**: Must be in `spm12/` directory (included in repo)
- **FSL**: Required for preprocessing, must be initialized before use

## Development Notes

### Testing
- SPM12 includes comprehensive test suite in `spm12/tests/`
- No dedicated HINEC test framework present
- Manual testing through sample data in `nifti_sample/`

### Sample Data
- `nifti_sample/original_sample/`: Basic diffusion data
- `nifti_sample/parcellation_sample/`: Data with parcellation masks
- Pre-computed results: `sample_parcellated.mat`

### Tractography
- **FACT Algorithm**: Fiber Assignment by Continuous Tracking with voxel boundary intersection
- **High-order Methods**: Available in `nim_tractography_highorder.m`
- **Seeding Strategies**: FA-based (current) and grid-based (recommended for uniform coverage)
- **Output**: Tracks as cell arrays, each track is N×3 matrix of 3D coordinates
- **Parameters**: FA thresholds, angle constraints, step size, seed density
- **File Structure**: Saved as `.mat` files with tracks, options, and statistics

## Important Implementation Details

### Data Structure (nim)
The main data structure contains:
- `.img`: 4D diffusion image data
- `.evec`: Eigenvectors from diffusion tensors
- `.eval`: Eigenvalues
- `.FA`: Fractional anisotropy maps
- `.parcellation_mask`: Brain region labels
- `.labels`: Anatomical region names

### Path Management
The `main.m` function automatically adds all necessary paths:
- All nim_* subdirectories
- SPM12 with subpaths
- Utility directories

### Preprocessing Integration
- Automatic detection of raw vs. preprocessed data
- FSL preprocessing called when raw data detected
- Preprocessing includes motion correction, eddy current correction, brain extraction

## Visualization Components (nim_visualization/)

### Traditional Visualization (MATLAB)
- **`visualizeTractography.m`**: Unified visualization system with multiple modes:
  - `'whole'` - Complete 3D brain view with track statistics
  - `'region'` - Single region analysis with filtering options
  - `'grid'` - All regions in grid layout
  - `'sequential'` - Region-by-region sequential viewing
  - Export capabilities in PNG, PDF, EPS, FIG formats

- **`visualizeTractographySlices.m`**: Interactive 2D slice viewer:
  - Three orthogonal views: Axial (top), Sagittal (side), Coronal (front)
  - Interactive sliders for slice navigation
  - Crosshair synchronization across views
  - Adjustable slice thickness and track filtering
  - Real-time updates with FA background overlay
  - **Performance**: 5-30 second delays per slice update (real-time computation)

### Fast Slice Viewer (Distributed Workflow)
For high-performance viewing with instant navigation, use the distributed workflow:

**Server (MATLAB)**: Pre-generate all slice images
```matlab
addpath('nim_visualization');
generateSlices('tracks.mat', 'nim.mat', '/export/slices');
```

**Transfer**: Copy cache directory to local computer
```bash
rsync -avz server:/export/slices/ ~/local/slices/
```

**Local (Python only - no MATLAB)**: View with instant navigation
```bash
./viewSlices.sh ~/local/slices/
```

**Performance**: Sub-100ms slice transitions (pre-computed images)

**Key Files**:
- `generateSlices.m` - Server-side batch image generation
- `generateTractographySliceCache.m` - Optimized rendering pipeline
- `TractographyCacheManager.m` - Cache directory management
- `FastTractographyViewer.py` - Python GUI for instant viewing
- `viewSlices.sh` - Local viewer launcher

See `DISTRIBUTED_WORKFLOW.md` for complete setup instructions.

### Track Data Structure
Each track in the `tracks` cell array represents a complete fiber pathway:
```matlab
track = [x1, y1, z1;    % First point (start)
         x2, y2, z2;    % Second point
         ...             % Intermediate points along fiber
         xN, yN, zN];   % Final point (end)
```
- **Not just endpoints**: Each track contains the complete 3D trajectory
- **Variable length**: Tracks contain 10-1000+ points depending on fiber length
- **Coordinate system**: Voxel space coordinates (X, Y, Z)

### Current Seeding Strategy (Needs Improvement)
- FA-based seeding creates uneven track distribution
- Only places seeds where FA exceeds threshold (typically 0.15-0.2)
- Results in sparse coverage in low-anisotropy regions
- **Recommendation**: Switch to uniform grid-based seeding for complete brain coverage

## Development Recommendations

### Tractography Improvements Needed
1. **Grid-based seeding**: Replace FA-threshold seeding with uniform 3D grid
2. **Dense coverage**: Use very dense seed spacing (0.5-1mm intervals)
3. **Brain mask constraint**: Limit seeds to brain tissue only (not CSF/background)
4. **Performance optimization**: Pre-compute seed grids and track intersections