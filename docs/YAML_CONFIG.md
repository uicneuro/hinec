# HINEC YAML Configuration System

## Overview

The HINEC pipeline now supports **YAML-based parameter configuration**, making it easy to:

- Manage complex parameter sets in human-readable files
- Version control your experimental configurations
- Share reproducible analysis pipelines
- Switch between different parameter presets instantly

## Quick Start

### Basic Usage

```bash
# Use default configuration
./run_hinec.sh data/parcellation_sample/sample sample.mat

# Use custom configuration
./run_hinec.sh data/parcellation_sample/sample sample.mat config/hinec_dti_cubic.yml
```

### Available Presets

| Config File | Purpose | Speed | Quality |
|-------------|---------|-------|---------|
| `hinec_default.yml` | Balanced performance | Medium | Good |
| `hinec_dti_cubic.yml` | Publication quality | Slow | Excellent |
| `hinec_dti_fast.yml` | Quick testing | Fast | Moderate |
| `irontract.yml` | IronTract challenge | Medium | High |

## Configuration File Structure

### Complete Example

```yaml
# my_experiment.yml
preprocessing:
  run_denoising: true
  denoise_method: 'dwidenoise'
  run_motion_correction: true
  run_eddy: true
  improve_mask: true
  atlas_type: 'jhu'

tractography:
  algorithm: 'hinec'
  integration_order: 4
  interp_method: 'trilinear'
  step_size: 0.2
  seed_density: 4
  termination_fa: 0.15
  angle_thresh: 35
  max_steps: 1000
  min_length: 35
  act_enabled: true

output:
  tractography_dir: 'tractography_results'
  include_timestamp: true
```

## Parameter Reference

### Preprocessing Section

```yaml
preprocessing:
  # Denoising
  run_denoising: true|false
  denoise_method: 'dwidenoise'|'susan'|'nlmeans'|'gaussian'

  # Correction
  run_motion_correction: true|false
  run_eddy: true|false

  # Masks and atlases
  improve_mask: true|false
  atlas_type: 'jhu'|'aal'|'desikan'

  # T1 registration (optional)
  t1_available: true|false
  use_t1_registration: true|false
  register_to_mni: true|false
```

### Tractography Section

#### Algorithm Selection

```yaml
tractography:
  algorithm: 'standard'|'hinec'
```

**Options**:

- `standard`: FACT algorithm (Euler integration, no interpolation)
- `hinec`: High-order (RK2/RK4/RKF45, trilinear interpolation, ACT)

#### Integration Method

```yaml
  integration_order: 1|2|4|5
  interp_method: 'trilinear'|'none'
```

**Integration Orders**:

- `1`: Euler (FACT default)
- `2`: Runge-Kutta 2nd order
- `4`: Runge-Kutta 4th order (HINEC default)
- `5`: RKF45 adaptive (highest accuracy)

#### Step Size Control

```yaml
  step_size: 0.2          # Fixed step (voxel units)

  # RKF45 adaptive parameters (integration_order=5 only)
  adaptive_step: true|false
  rkf_tolerance: 0.01     # Error tolerance (voxels)
  rkf_safety: 0.9         # Safety factor (0-1)
  step_min: 0.01          # Minimum step (voxels)
  step_max: 1.0           # Maximum step (voxels)
```

**Guidelines**:

- `step_size`: 0.1-0.2 (precise), 0.2-0.3 (balanced), 0.3-0.5 (fast)
- `rkf_tolerance`: 0.005 (high precision), 0.01 (balanced), 0.02 (fast)

#### Seeding Parameters

```yaml
  seed_density: 4         # Seeds per voxel (1-8 recommended)
  seed_strategy: 'uniform'
```

**Seed Density Guidelines**:

- `1-2`: Fast exploration
- `4`: Balanced coverage (default)
- `6-8`: Dense coverage for critical regions

#### Termination Criteria

```yaml
  termination_fa: 0.15    # FA threshold for termination
  angle_thresh: 35        # Max curvature (degrees)
  max_steps: 1000         # Max integration steps
  min_length: 35          # Min track length (mm)
```

**Guidelines**:

- `termination_fa`: 0.05 (long tracks), 0.15 (balanced), 0.20 (short/fast)
- `angle_thresh`: 30° (strict), 35-45° (balanced), 60° (permissive)

#### ACT Configuration

```yaml
  act_enabled: true|false
```

ACT (Anatomically Constrained Tractography) automatically uses tissue masks generated during preprocessing (WM/GM/CSF).

### Output Section

```yaml
output:
  tractography_dir: 'tractography_results'
  include_timestamp: true|false
  auto_visualize: true|false
```

## Creating Custom Configurations

### Step 1: Copy Template

```bash
cp config/hinec_default.yml config/my_experiment.yml
```

### Step 2: Edit Parameters

```bash
nano config/my_experiment.yml
```

### Step 3: Run with Custom Config

```bash
./run_hinec.sh data/mydata data_processed.mat config/my_experiment.yml
```

## Common Configuration Patterns

### Publication-Quality Results

```yaml
tractography:
  integration_order: 5       # RKF45 adaptive
  adaptive_step: true
  rkf_tolerance: 0.005       # Tight tolerance
  step_size: 0.1             # Small initial step
  seed_density: 6            # Dense seeding
  termination_fa: 0.1        # Continue in low FA
  angle_thresh: 30           # Strict curvature
```

**Expected runtime**: 2-3x slower than default

**Benefits**: Maximum accuracy, smooth tracks, fewer spurious connections

---

### Fast Exploration

```yaml
tractography:
  integration_order: 2       # RK2
  adaptive_step: false
  step_size: 0.3             # Larger steps
  seed_density: 2            # Sparse seeding
  termination_fa: 0.2        # Stop earlier
  angle_thresh: 45           # More permissive
```

**Expected runtime**: 3-5x faster than default

**Benefits**: Quick parameter testing, exploratory analysis

---

### IronTract Challenge

```yaml
tractography:
  integration_order: 4       # RK4
  step_size: 0.2
  seed_density: 4
  termination_fa: 0.1        # Complete tracking
  angle_thresh: 35
  min_length: 35
```

**Expected runtime**: Similar to default

**Benefits**: Optimized for phantom validation, proven effective

---

### Low-Anisotropy Regions (Fornix, Cingulum)

```yaml
tractography:
  integration_order: 4
  step_size: 0.15            # Smaller steps
  seed_density: 6            # Dense seeding
  termination_fa: 0.05       # Very low threshold
  angle_thresh: 45           # More permissive
  min_length: 20             # Shorter tracks OK
```

**Expected runtime**: Slower (more seeds, longer tracks)

**Benefits**: Better coverage in challenging white matter regions

---

## Parameter Validation

The configuration system includes automatic validation:

```matlab
✓ Parameter validation passed
```

**Validation checks**:

- `step_size > 0`
- `angle_thresh ∈ (0, 180]`
- `integration_order ∈ {1, 2, 4, 5}`
- `rkf_tolerance > 0` (if RKF45 enabled)
- `step_min < step_max` (if RKF45 enabled)
- `rkf_safety ∈ (0, 1]` (if RKF45 enabled)

**Warnings**:

- Unusual `seed_density` outside 1-8 range

## Legacy Compatibility

The YAML config system is **fully backward compatible**:

```matlab
% Legacy usage still works
runTractography('data.mat', 'hinec');

% New YAML usage
config = load_config_yaml('config/my_experiment.yml');
runTractography('data.mat', config);
```

## Troubleshooting

### Config File Not Found

```bash
Error: Configuration file not found: config/missing.yml
Available configs:
  config/hinec_default.yml
  config/hinec_dti_cubic.yml
  config/hinec_dti_fast.yml
  config/irontract.yml
```

**Solution**: Check filename and path

---

### YAML Parsing Error

```
Error parsing YAML at line 15: Invalid value
```

**Common issues**:

- Missing quotes around strings with special characters
- Incorrect indentation (use 2 spaces, not tabs)
- Missing colons after keys

**Example fix**:

```yaml
# WRONG
tractography:
algorithm: hinec

# CORRECT
tractography:
  algorithm: 'hinec'
```

---

### Parameter Validation Failed

```
Error: step_size must be positive, got: -0.5
```

**Solution**: Check parameter values against validation rules in this document

---

## Advanced Usage

### Programmatic Configuration

```matlab
% Load and modify config in MATLAB
config = load_config_yaml('config/hinec_default.yml');

% Override specific parameters
config.tractography.seed_density = 8;
config.tractography.termination_fa = 0.05;

% Run with modified config
runTractography('data.mat', config);
```

### Batch Processing

```bash
#!/bin/bash
# Process multiple subjects with same config

for subject in subject_001 subject_002 subject_003; do
    ./run_hinec.sh data/${subject} processed/${subject}.mat config/hinec_dti_cubic.yml
done
```

### Config Versioning

```bash
# Track configs in git for reproducibility
git add config/experiment_v1.yml
git commit -m "Add config for Experiment 1 (high precision, dense seeding)"
```

## Best Practices

1. **Start with presets**: Use `hinec_default.yml` or `hinec_dti_cubic.yml` as starting points
2. **Document changes**: Add comments explaining why parameters were changed
3. **Version control**: Keep configs in git for reproducibility
4. **Test incrementally**: Change one parameter at a time when optimizing
5. **Monitor logs**: Check `hinec_*.log` files for parameter confirmation and diagnostics

## See Also

- [RKF_Usage.md](RKF_Usage.md) - Detailed RKF45 adaptive integration guide
- [TRACTOGRAPHY.md](TRACTOGRAPHY.md) - Tractography algorithm documentation
- [PREPROCESSING.md](../PREPROCESSING.md) - Preprocessing pipeline details
