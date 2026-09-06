# YAML Configuration Changelog

## Schema rework - canonical nested configuration (current)

Old configs still load; each superseded key produces a deprecation warning naming
its replacement.

- `integration_order` -> `integrator.method`. The old key was a **method selector
  wearing a numeric-order name**: `5` meant RKF45, a 4(5) embedded pair whose order
  is 4, not 5.
- `max_steps` -> `termination.max_arc`. The old key counted integration steps, so
  halving the step size silently halved how far a track could travel. `max_arc` is an
  arc length in voxels; `max_steps` is derived as `ceil(max_arc / step)`.
- Config gained **two levels of nesting** (`section.group.key`) and inline lists.
  A third level is a parse error. The previous parser was indentation-blind and
  silently misparsed nested YAML rather than rejecting it.
- Unknown keys are now an **error**, not a silent no-op - which is how
  `rkf_tolerance` (hinec) and `rkf_tol` (mmf) came to mean the same thing under two
  names.
- Retired eight keys no tracker read (`gate_power`, `crossing_cp`, `curv_beta`,
  `crossing_detect`, `swing_ratio_max`, `transport_gate`, `transport_strength`,
  `bishop_eps`) plus `fa_threshold`, which was printed but never used to terminate.
- Every parameter, `preprocessing.*` included, is reachable from the CLI via
  `--set <path>=<value>` on both launchers.
- `docs/YAML_CONFIG.md` is now **generated** from `src/nim_utils/nim_config_schema.m`,
  with a test that fails if code and docs drift.

Everything below predates the rework and is kept as history.

---

# YAML Configuration System - Implementation Summary

## Overview

Implemented a comprehensive YAML-based parameter configuration system for the HINEC pipeline, making parameter management easy and reproducible.

## What Changed

### New Files Created

**Configuration System**:

- `nim_utils/load_config_yaml.m` - YAML parser with validation
- `config/hinec_default.yml` - Default balanced configuration
- `config/hinec_dti_cubic.yml` - Publication-quality preset
- `config/hinec_dti_fast.yml` - Quick testing preset
- `config/irontract.yml` - IronTract challenge preset

**Documentation**:

- `docs/YAML_CONFIG.md` - Complete configuration reference (124 KB)
- `config/README.md` - Quick reference for config directory
- `test_yaml_config.m` - Validation test script

### Modified Files

**Core Pipeline**:

- `run_hinec.sh` - Added optional config file argument
- `main.m` - Added YAML config support (backward compatible)
- `runTractography.m` - Added YAML config support (backward compatible)
- `README.md` - Added YAML quick start section

## Key Features

### 1. Easy Parameter Management

**Before** (editing source code):
```matlab
% Edit runTractography.m:88-93
options.seed_density = 4;
options.step_size = 0.5;
options.termination_fa = 0.15;
options.angle_thresh = 35;
options.max_steps = 1000;
options.min_length = 35;
```

**After** (editing YAML file):
```yaml
# config/my_experiment.yml
tractography:
  seed_density: 4
  step_size: 0.5
  termination_fa: 0.15
  angle_thresh: 35
  max_steps: 1000
  min_length: 35
```

### 2. Preset Configurations

Four ready-to-use presets for common scenarios:

| Preset | Integration | Speed | Use Case |
|--------|-------------|-------|----------|
| `hinec_default.yml` | RK4 | Baseline | Standard analysis |
| `hinec_dti_cubic.yml` | RKF45 adaptive | 2-3x slower | Publications |
| `hinec_dti_fast.yml` | RK2 | 3-5x faster | Parameter testing |
| `irontract.yml` | RK4 | Baseline | IronTract challenge |

### 3. Automatic Validation

Parameters are validated when loaded:
```matlab
✓ Parameter validation passed
```

**Checks include**:

- Step size must be positive
- Angle threshold in (0, 180]
- Integration order in {1, 2, 4, 5}
- RKF bounds: step_min < step_max
- RKF safety factor in (0, 1]
- RKF tolerance > 0

### 4. Backward Compatibility

All existing code continues to work:

```matlab
% Legacy usage - still works
runTractography('data.mat', 'hinec');

% New YAML usage
config = load_config_yaml('config/hinec_dti_cubic.yml');
runTractography('data.mat', config);
```

## Usage Examples

### Basic Usage

```bash
# Use default config
./run_hinec.sh data/subject processed.mat

# Use custom config
./run_hinec.sh data/subject processed.mat config/hinec_dti_cubic.yml
```

### Create Custom Config

```bash
# Copy template
cp config/hinec_default.yml config/my_experiment.yml

# Edit parameters
nano config/my_experiment.yml

# Run with custom config
./run_hinec.sh data/subject processed.mat config/my_experiment.yml
```

### Programmatic Usage

```matlab
% Load and modify config
config = load_config_yaml('config/hinec_default.yml');
config.tractography.seed_density = 8;
config.tractography.termination_fa = 0.05;

% Run with modified config
runTractography('data.mat', config);
```

## Parameter Organization

### Preprocessing Section

- Denoising settings
- Motion/eddy correction
- Brain mask improvement
- Atlas selection
- T1 registration (optional)

### Tractography Section
- Algorithm selection (standard/hinec)
- Integration method (Euler/RK2/RK4/RKF45)
- Step size control
- RKF45 adaptive parameters
- Seeding strategy
- Termination criteria
- ACT configuration

### Output Section
- Output directory
- Filename timestamps
- Auto-visualization

## Testing

Run validation test:

```matlab
test_yaml_config
```

Expected output:
```
=== HINEC YAML Configuration System Test ===

Test 1: Loading default configuration...
  ✓ Default config loaded successfully

Test 2: Loading high precision configuration...
  ✓ High precision config loaded

Test 3: Parameter validation...
  ✓ All parameter ranges validated

Test 4: Checking preset configurations...
  ✓ All preset configs found

Test 5: RKF45 parameter validation...
  ✓ RKF45 parameters validated

Test 6: Comparing configurations...
  [comparison table]

========================================
ALL TESTS PASSED ✓
========================================
```

## Benefits

### For Users
1. **Easy parameter changes** - Edit YAML files instead of source code
2. **Reproducibility** - Version control configurations in git
3. **Presets** - Quick access to validated parameter sets
4. **Safety** - Automatic validation prevents invalid parameters

### For Developers
1. **Centralized configuration** - All parameters in one place
2. **Type safety** - Validation catches errors early
3. **Backward compatibility** - Legacy code still works
4. **Extensibility** - Easy to add new parameters

## Documentation

- **Complete reference**: [docs/YAML_CONFIG.md](docs/YAML_CONFIG.md)
- **Quick reference**: [config/README.md](config/README.md)
- **RKF45 guide**: [docs/RKF_Usage.md](docs/RKF_Usage.md)
- **Main README**: [README.md](README.md)

## Migration Guide

### For Existing Users

No changes required! Existing usage patterns continue to work:

```matlab
% This still works exactly as before
runTractography('data.mat', 'hinec');
```

To use YAML configs:

```matlab
% New usage (recommended)
config = load_config_yaml('config/hinec_default.yml');
runTractography('data.mat', config);
```

### For Script Users

Update shell scripts to use config files:

**Before**:
```bash
./run_hinec.sh data/subject processed.mat
```

**After** (optional, but recommended):
```bash
./run_hinec.sh data/subject processed.mat config/my_experiment.yml
```

## Future Enhancements

Potential improvements:

1. JSON support (in addition to YAML)
2. GUI config editor
3. Config validation CLI tool
4. Parameter optimization suggestions
5. Automatic config generation from logs

## Technical Details

### YAML Parser

Simple MATLAB-compatible YAML parser that handles:

- Key-value pairs
- Nested sections
- Comments
- Booleans (`true`/`false`)
- Numbers (integers and floats)
- Strings (quoted and unquoted)

**Limitations**:

- No support for arrays/lists
- No support for complex YAML features (anchors, aliases, etc.)
- For HINEC config needs only

For complex YAML needs, consider visiting: [yamlmatlab](https://github.com/ewiger/yamlmatlab)

### Validation Strategy

Two-phase validation:

1. **Syntax validation** - During YAML parsing
2. **Semantic validation** - After loading complete config

**Validation errors halt execution early** - prevents invalid tractography runs.

## Credits

Implementation Date: 2025-01-18

Based on code review findings and parameter management requirements.

## See Also

- [Code Review Report](docs/CODE_REVIEW_2025-01-18.md) - Full analysis
- [RKF Implementation](docs/RKF.md) - RKF45 technical details
- [Tractography Guide](docs/TRACTOGRAPHY.md) - Algorithm documentation
