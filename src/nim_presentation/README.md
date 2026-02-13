# HINEC Presentation Figure Generation

Automated figure generation for HINEC academic poster presentation. Creates high-quality methodological diagrams demonstrating interpolation, integration, and tractography techniques.

## Quick Start

### Generate All Figures

From MATLAB command line in the HINEC root directory:

```matlab
generate_presentation_figures()
```

This creates 7 figures in the `presentation_figures/` directory.

### Custom Output Directory

```matlab
generate_presentation_figures('my_poster_figures')
```

## Generated Figures

### 0. Grid-Based Tracking Concept
**File**: `0_grid_tracking_simple.png`

**Content**:
- Panel A: Discrete voxel grid with direction vectors at voxel centers
- Panel B: Smooth interpolated track path through direction field
- Panel C: Combined view showing both discrete and continuous aspects

**Use in presentation**: Introduction to basic tractography concept - how discrete voxel directions create smooth continuous paths

---

### 1. Interpolation Comparison (3D)
**File**: `1_interpolation_comparison.png`

**Content**:
- Panel A: FACT (no interpolation) - discrete voxel jumps
- Panel B: Trilinear interpolation - smooth transitions
- Panel C: Cubic interpolation - smoothest trajectories

**Use in presentation**: Main comparison of interpolation methods

---

### 2. Interpolation Field (2D Slice)
**File**: `2_interpolation_field_comparison.png`

**Content**:
- Direction field visualization at Z=3 slice
- Shows how interpolation creates smooth vector fields
- Demonstrates sub-voxel precision

**Use in presentation**: Detailed view of interpolation mechanics

---

### 3. Integration Methods (Algorithms)
**File**: `3_integration_method_comparison.png`

**Content**:
- Panel A: RK4 algorithm pseudocode and properties
- Panel B: RKF45 algorithm with adaptive step control
- Panel C: Visual comparison at high-curvature region

**Use in presentation**: Explain RK4 vs RKF45 differences

---

### 4. Step Size Adaptation
**File**: `4_step_size_adaptation.png`

**Content**:
- Panel A: RK4 fixed step sizes (bar chart)
- Panel B: RKF45 adaptive step sizes with curvature annotation

**Use in presentation**: Demonstrate adaptive stepping advantage

---

### 5. Tractography Step-by-Step
**File**: `5_tractography_step_by_step.png`

**Content** (6 panels):
1. Seed placement on FA map
2. Direction extraction from tensor field
3. Integration (complete track)
4. ACT tissue constraints (WM/GM/CSF)
5. Termination criteria (FA threshold)
6. Final result with statistics

**Use in presentation**: Complete workflow explanation

---

### 6. 3D Tractography Visualization
**File**: `6_tractography_3d_visualization.png`

**Content**:
- Panel A: Multiple fiber tracks in 3D space
- Panel B: Track density map (connectivity matrix)

**Use in presentation**: Show realistic multi-track reconstruction

---

## Individual Function Usage

### Generate Only Grid Tracking Concept

```matlab
addpath('nim_presentation');
visualize_grid_tracking('output_dir');
```

Generates figure 0.

### Generate Only Interpolation Figures

```matlab
addpath('nim_presentation');
visualize_interpolation_methods('output_dir');
```

Generates figures 1-2.

### Generate Only Integration Figures

```matlab
addpath('nim_presentation');
visualize_integration_methods('output_dir');
```

Generates figures 3-4.

### Generate Only Tractography Example

```matlab
addpath('nim_presentation');
visualize_tractography_example('output_dir');
```

Generates figures 5-6.

---

## Poster Layout Recommendations

### Option 1: Horizontal Layout with Introduction

```
┌──────────────────────────────────────────┐
│          Figure 0 (Introduction)         │
│         Grid Tracking Concept            │
├─────────────┬─────────────┬──────────────┤
│   Figure 1  │   Figure 3  │   Figure 5   │
│Interpolation│ Integration │ Tractography │
│  Comparison │    Methods  │   Workflow   │
├─────────────┼─────────────┼──────────────┤
│   Figure 2  │   Figure 4  │   Figure 6   │
│Field Detail │  Step Size  │  3D Tracks   │
└─────────────┴─────────────┴──────────────┘
```

### Option 2: Vertical Layout (Methods Focus)

```
┌───────────────────────────┐
│      Introduction         │
│   Figure 0: Grid Concept  │
├───────────────────────────┤
│   Figure 1 + Figure 2     │
│   Interpolation Methods   │
├───────────────────────────┤
│   Figure 3 + Figure 4     │
│   Integration Methods     │
├───────────────────────────┤
│   Figure 5 + Figure 6     │
│   Complete Tractography   │
├───────────────────────────┤
│   Results & Conclusions   │
└───────────────────────────┘
```

---

## Customization

### Modify Figure Parameters

Edit the individual visualization scripts in `nim_presentation/`:

- **Grid size**: Change `grid_size` variable (default: 5)
- **Track length**: Modify `max_steps` (default: 20-25)
- **Colors**: Adjust `'r'`, `'b'`, `'g'` color specifications
- **Font sizes**: Change `FontSize` properties

### Example: Larger Grid

```matlab
% In visualize_interpolation_methods.m, line 23
grid_size = 10;  % Change from 5 to 10
```

### Example: More Tracks

```matlab
% In visualize_tractography_example.m, line 174
n_seeds = 20;  % Change from 10 to 20
```

---

## Figure Quality

### Resolution
- **Default**: MATLAB default DPI (96-150 DPI)
- **For printing**: Export as vector format

### Vector Export

```matlab
% After running generate_presentation_figures()
cd presentation_figures

% Convert to vector PDF
fig = openfig('1_interpolation_comparison.fig');
print(fig, '1_interpolation_comparison', '-dpdf', '-vector');
```

### High-DPI PNG

```matlab
% Modify saveas calls in visualization scripts
print(fig1, fullfile(output_dir, '1_interpolation_comparison.png'), ...
    '-dpng', '-r300');  % 300 DPI
```

---

## Troubleshooting

### Figures Don't Appear

**Issue**: MATLAB running in headless mode
**Solution**: Enable display

```matlab
set(0, 'DefaultFigureVisible', 'on');
generate_presentation_figures();
```

### Out of Memory

**Issue**: Large grid sizes or too many tracks
**Solution**: Reduce parameters

```matlab
% In visualization scripts
grid_size = 3;      % Smaller grid
n_seeds = 5;        % Fewer tracks
max_steps = 15;     % Shorter tracks
```

### Path Errors

**Issue**: Functions not found
**Solution**: Ensure correct working directory

```matlab
cd /path/to/hinec
addpath(genpath('.'));
generate_presentation_figures();
```

---

## Dependencies

### Required MATLAB Toolboxes
- **Core MATLAB**: No additional toolboxes required
- **Built-in functions used**:
  - `interp3` (interpolation)
  - `quiver`, `quiver3` (vector field plotting)
  - `imagesc` (heatmap visualization)
  - `plot3` (3D plotting)

### No External Dependencies
All visualization code is self-contained and uses only MATLAB built-in functions.

---

## File Structure

```
hinec/
├── generate_presentation_figures.m    # Main executor
├── nim_presentation/
│   ├── README.md                      # This file
│   ├── visualize_grid_tracking.m
│   ├── visualize_interpolation_methods.m
│   ├── visualize_integration_methods.m
│   └── visualize_tractography_example.m
└── presentation_figures/              # Output directory
    ├── 0_grid_tracking_simple.png
    ├── 1_interpolation_comparison.png
    ├── 2_interpolation_field_comparison.png
    ├── 3_integration_method_comparison.png
    ├── 4_step_size_adaptation.png
    ├── 5_tractography_step_by_step.png
    └── 6_tractography_3d_visualization.png
```

---

## Academic Use

### Citation
When using these figures in presentations or publications:

```
HINEC (High-order Neural Connectivity) Tractography Pipeline
[Your Institution], [Year]
Figures generated using automated visualization tools
```

### Figure Captions (Suggested)

**Figure 1**: *Comparison of interpolation methods in fiber tractography. (A) FACT algorithm with no interpolation shows discrete voxel jumps. (B) Trilinear interpolation provides smooth direction transitions. (C) Cubic interpolation achieves highest smoothness.*

**Figure 3**: *Integration method comparison. (A) RK4 uses fixed 4th-order integration. (B) RKF45 employs 5th-order integration with adaptive step size control. (C) Behavioral difference at high-curvature regions shows RKF45 adapts step size to local geometry.*

**Figure 5**: *Complete tractography workflow. Step-by-step demonstration of (1) seed placement, (2) direction extraction, (3) integration, (4) anatomical constraints, (5) termination criteria, and (6) final track reconstruction.*

---

## Future Enhancements

Potential additions for future versions:

- [ ] Real DTI data integration (currently synthetic)
- [ ] Crossing fiber region examples
- [ ] Quantitative comparison metrics
- [ ] Animation generation (GIF/MP4)
- [ ] Interactive 3D viewers
- [ ] Batch processing for multiple datasets

---

## Support

For questions or issues:
1. Check HINEC main documentation: `../README.md`
2. Review presentation outline: `../docs/presentation_outline.md`
3. Check MATLAB version compatibility (R2018b or later recommended)

---

**Last Updated**: 2025-01-19
**HINEC Version**: 1.0
**Author**: HINEC Development Team
