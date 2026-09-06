# YAML Configuration Changelog

## Schema rework - canonical nested configuration (current)

The configuration surface is now declared once, in
`src/nim_utils/nim_config_schema.m`. Defaults, validation, unknown-key
rejection, legacy migration, the flat option names the trackers read, and the
[configuration reference](YAML_CONFIG.md) are all derived from it. Old configs
still load; each superseded key produces a deprecation warning naming its
replacement.

### Structure

- Config gained **two levels of nesting** (`section.group.key`) and inline lists.
  A third level is a parse error. The previous parser was indentation-blind and
  silently misparsed nested YAML rather than rejecting it.
- Unknown keys are now an **error**, not a silent no-op - which is how
  `rkf_tolerance` (hinec) and `rkf_tol` (mmf) came to mean the same thing under two
  names.
- Every parameter, `preprocessing.*` included, is reachable from the CLI via
  `--set <path>=<value>` on both launchers.
- [YAML_CONFIG.md](YAML_CONFIG.md) is now **generated** from
  `src/nim_utils/nim_config_schema.m` by `src/nim_utils/nim_config_docs.m`, with
  a test that fails if code and docs drift. Do not edit it by hand.

### Renamed and redefined keys

- `integration_order` -> `integrator.method`. The old key **spelled a method
  selector as a number**: `5` selected RKF45. That value was not wrong
  numerically - the implementation uses Dormand-Prince coefficients and advances
  on the 5th-order solution, keeping the embedded 4th-order one for error control
  - but a method name belongs in a method key. The values are now named: `euler | rk2 | rk4 | rkf45`.
- `max_steps` -> `termination.max_arc`. The old key counted integration steps, so
  halving the step size silently halved how far a track could travel. `max_arc` is an
  arc length in voxels; `max_steps` is derived as `ceil(max_arc / step)`.
- `min_length` -> `termination.min_arc`, an arc length in **voxels**, and it is
  now actually enforced. It had been defaulted and then never read by the hinec
  and standard trackers - a schema-validated key that every config set was
  honoured by exactly one of the three.
- `termination.angle_max` is documented, and implemented, as degrees of turning
  per **voxel of arc**, not per step: the budget for one step is
  `angle_max x step`, so refining the step does not loosen the constraint.
  Because tangents are sign-aligned (`v1` is a line field), a measured turn
  never exceeds 90 degrees, so any budget above `90/step` is *inert* rather than
  merely loose; the loader warns when a config lands there. `angle_max: 0`
  disables the criterion outright. See `src/nim_utils/nim_angle_limit.m`.
- `interpolation.method` gained `spline` alongside `trilinear` and `cubic`. The
  three differ in smoothness - C0, C1 (Keys cubic convolution, *not* a spline)
  and C2 - which caps the order a Runge-Kutta method can reach. See
  [Solution Verification](CONVERGENCE.md) for the measured orders.
- `interpolation.upsample` is implemented: the direction field is sampled on a
  grid of spacing `1/upsample` voxels before the interpolants are built. Above 1
  refines, below 1 coarsens. The coordinate frame is unchanged, so step sizes and
  lengths stay comparable across factors.
- `tractography.act` now defaults to **false**. It is read only by the hinec
  tracker.
- New filter predicates `filter.endpoints_in` (one end in each of two regions)
  and `filter.contained_in` (every point inside a corridor). Together with
  `filter.include_roi` and `filter.exclude_roi` these express the ISMRM 2015
  bundle definition.

### Retired keys

Retired keys are accepted with a warning and dropped, never silently migrated -
mapping a key that had no effect onto one that does would change behaviour.
`src/nim_utils/nim_config_retired.m` holds the list and the reason for each.

- Eight keys no tracker ever read: `gate_power`, `crossing_cp`, `curv_beta`,
  `crossing_detect`, `swing_ratio_max`, `transport_gate`, `transport_strength`,
  `bishop_eps`.
- `fa_threshold`, which was printed but never used to terminate. Use
  `termination.fa_min` to stop tracking, or `seeding.fa_min` to restrict seeding.
- `order`, a legacy compatibility key read by no tracker as an option.
- `sel_power`, **removed from the hinec tracker entirely**. It re-weighted each
  stencil voxel by `alignment^sel_power`, biasing interpolation toward the
  incoming direction. For DTI there is one principal eigenvector per voxel, so
  there was nothing to disambiguate; it simply bent a single-valued field toward
  the current heading, which is self-reinforcing, and it made the ODE
  direction-dependent (`dx/ds = v(x, dx/ds)`), so classical Runge-Kutta order
  theory no longer applied. HINEC tracking is now interpolation and integration
  only. CSD still selects the FOD peak nearest the incoming tangent, because
  reducing a multi-valued field to one direction is structurally required, but
  the alignment exponent on top of it is gone. The MMF tracker keeps a separate
  `mmf.frame_sel_power`, used when building the moving frame.

---

## History

Everything below predates the schema rework and is kept as a record of how the
configuration system was introduced. It describes the **flat key names and
validation rules that the rework replaced**; for current behaviour read the
section above, or the generated [configuration reference](YAML_CONFIG.md).

### Overview

Implemented a comprehensive YAML-based parameter configuration system for the HINEC pipeline, making parameter management easy and reproducible.

### What Changed

#### New Files Created

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

#### Modified Files

**Core Pipeline**:

- `run_hinec.sh` - Added optional config file argument
- `main.m` - Added YAML config support (backward compatible)
- `runTractography.m` - Added YAML config support (backward compatible)
- `README.md` - Added YAML quick start section

> Paths in the two lists above predate the source reorganisation: the MATLAB
> sources moved under `src/`, the launchers under `bin/`, and the test script to
> `tests/test_yaml_config.m`. `docs/YAML_CONFIG.md` is now generated rather than
> written by hand.

### Key Features

#### 1. Easy Parameter Management

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

> **Superseded.** These flat key names no longer load without a deprecation
> warning. The canonical form is `tractography.seeding.density`,
> `tractography.integrator.step`, `tractography.termination.fa_min`,
> `termination.angle_max`, `termination.max_arc` and `termination.min_arc`.

#### 2. Preset Configurations

Four ready-to-use presets for common scenarios:

| Preset | Integration | Speed | Use Case |
|--------|-------------|-------|----------|
| `hinec_default.yml` | RK4 | Baseline | Standard analysis |
| `hinec_dti_cubic.yml` | RKF45 adaptive | 2-3x slower | Publications |
| `hinec_dti_fast.yml` | RK2 | 3-5x faster | Parameter testing |
| `irontract.yml` | RK4 | Baseline | IronTract challenge |

> **Superseded.** The config set has since been reorganised around one naming
> rule per family - tracker configs `<algorithm>_<field>[_<variant>].yml`,
> dataset configs `<dataset>[_variant].yml`. See `config/README.md` for the
> current list and what each preset is for.

#### 3. Automatic Validation

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

> **Superseded.** Ranges and permitted values now come from the `range` and
> `allowed` fields of `src/nim_utils/nim_config_schema.m`, so this list cannot
> drift from the code. Two entries are no longer true as written: integration
> order is gone (`integrator.method` takes `euler | rk2 | rk4 | rkf45`), and
> `termination.angle_max` is a rate in degrees per voxel of arc with a range of
> `[0, 3600]`, where `0` means "criterion disabled".

#### 4. Backward Compatibility

All existing code continues to work:

```matlab
% Legacy usage - still works
runTractography('data.mat', 'hinec');

% New YAML usage
config = load_config_yaml('config/hinec_dti_cubic.yml');
runTractography('data.mat', config);
```

### Usage Examples

> **Superseded.** The launchers now live under `bin/`, so every `./run_hinec.sh`
> below is `./bin/run_hinec.sh` today, and the programmatic example sets flat
> field names that have since been nested. Tractography-only iteration on an
> already-preprocessed dataset goes through `bin/run_tractography.sh`, which did
> not exist yet.

#### Basic Usage

```bash
# Use default config
./run_hinec.sh data/subject processed.mat

# Use custom config
./run_hinec.sh data/subject processed.mat config/hinec_dti_cubic.yml
```

#### Create Custom Config

```bash
# Copy template
cp config/hinec_default.yml config/my_experiment.yml

# Edit parameters
nano config/my_experiment.yml

# Run with custom config
./run_hinec.sh data/subject processed.mat config/my_experiment.yml
```

#### Programmatic Usage

```matlab
% Load and modify config
config = load_config_yaml('config/hinec_default.yml');
config.tractography.seed_density = 8;
config.tractography.termination_fa = 0.05;

% Run with modified config
runTractography('data.mat', config);
```

### Parameter Organization

#### Preprocessing Section

- Denoising settings
- Motion/eddy correction
- Brain mask improvement
- Atlas selection
- T1 registration (optional)

#### Tractography Section
- Algorithm selection (standard/hinec)
- Integration method (Euler/RK2/RK4/RKF45)
- Step size control
- RKF45 adaptive parameters
- Seeding strategy
- Termination criteria
- ACT configuration

#### Output Section
- Output directory
- Filename timestamps
- Auto-visualization

> **Superseded.** There is no top-level `output:` section: output directories
> and timestamps are handled by the run-directory system, and `runTractography`
> no longer draws figures. The only output key is
> `tractography.output.arc_step`. The tractography section has also gained a
> third algorithm (`mmf`), an `interpolation` group, ROI seeding and the track
> `filter` group.

### Testing

Run validation test:

```matlab
test_yaml_config
```

> The script now lives at `tests/test_yaml_config.m`, and the schema itself is
> covered by `tests/unit/TestConfigSchema.m` and `tests/unit/TestConfigYaml.m`.
> The transcript below is the original output and no longer matches.

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

### Benefits

#### For Users
1. **Easy parameter changes** - Edit YAML files instead of source code
2. **Reproducibility** - Version control configurations in git
3. **Presets** - Quick access to validated parameter sets
4. **Safety** - Automatic validation prevents invalid parameters

#### For Developers
1. **Centralized configuration** - All parameters in one place
2. **Type safety** - Validation catches errors early
3. **Backward compatibility** - Legacy code still works
4. **Extensibility** - Easy to add new parameters

### Documentation

- **Complete reference**: [YAML_CONFIG.md](YAML_CONFIG.md)
- **Quick reference**: `config/README.md` in the repository
- **RKF45 guide**: [RKF_Usage.md](RKF_Usage.md)
- **Main README**: `README.md` in the repository

### Migration Guide

#### For Existing Users

No changes required. Existing usage patterns continue to work:

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

#### For Script Users

Update shell scripts to use config files:

**Before**:
```bash
./run_hinec.sh data/subject processed.mat
```

**After** (optional, but recommended):
```bash
./run_hinec.sh data/subject processed.mat config/my_experiment.yml
```

### Future Enhancements

Potential improvements:

1. JSON support (in addition to YAML)
2. GUI config editor
3. Config validation CLI tool
4. Parameter optimization suggestions
5. Automatic config generation from logs

### Technical Details

#### YAML Parser

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

> **Superseded.** Parsing now lives in `src/nim_utils/nim_yaml_parse.m`. It
> supports inline lists (`roi: [41, 42]`), which ROI seeding and the track
> filters need, and it is indentation-aware: two levels below a section are
> accepted and a third is a parse error rather than a silent misparse. Anchors
> and aliases are still unsupported.

#### Validation Strategy

Two-phase validation:

1. **Syntax validation** - During YAML parsing
2. **Semantic validation** - After loading complete config

**Validation errors halt execution early** - prevents invalid tractography runs.

### Credits

Implementation Date: 2025-01-18

Based on code review findings and parameter management requirements.

### See Also

- [RKF Implementation](RKF.md) - RKF45 technical details
- [Tractography Guide](TRACTOGRAPHY.md) - the FACT tracker
- [YAML Configuration Reference](YAML_CONFIG.md) - the current, generated
  parameter reference
