# HINEC Configuration Files

This directory contains YAML configuration files for the HINEC tractography pipeline.

## Available Configurations

### `hinec_default.yml` (Recommended)
**Balanced performance and quality**
- Integration: RK4
- Step size: 0.2 voxels
- Seed density: 4 seeds/voxel
- Runtime: Baseline

**Use for**: Standard analysis, most datasets

```bash
./run_hinec.sh data/subject processed.mat
# Or explicitly:
./run_hinec.sh data/subject processed.mat config/hinec_default.yml
```

---

### `high_precision.yml`
**Maximum accuracy for publication**
- Integration: RKF45 adaptive (tol=0.005)
- Step size: 0.1 voxels (initial)
- Seed density: 6 seeds/voxel
- Runtime: 2-3x slower

**Use for**: Publication figures, critical connectivity analysis

```bash
./run_hinec.sh data/subject processed.mat config/high_precision.yml
```

---

### `fast_exploration.yml`
**Quick parameter testing**
- Integration: RK2
- Step size: 0.3 voxels
- Seed density: 2 seeds/voxel
- Runtime: 3-5x faster

**Use for**: Parameter screening, exploratory analysis

```bash
./run_hinec.sh data/subject processed.mat config/fast_exploration.yml
```

---

### `irontract.yml`
**IronTract Challenge optimized**
- Integration: RK4
- Step size: 0.2 voxels
- Seed density: 4 seeds/voxel
- Runtime: Baseline

**Use for**: IronTract phantom validation, ISMRM challenges

```bash
./run_hinec.sh ironTract/subject processed.mat config/irontract.yml
```

---

## Creating Custom Configs

1. **Copy a template**:
   ```bash
   cp hinec_default.yml my_experiment.yml
   ```

2. **Edit parameters**:
   ```bash
   nano my_experiment.yml
   ```

3. **Run with custom config**:
   ```bash
   ./run_hinec.sh data/subject processed.mat config/my_experiment.yml
   ```

## Key Parameters Quick Reference

### Integration Quality
```yaml
integration_order: 1  # Euler (FACT)
integration_order: 2  # RK2 (fast)
integration_order: 4  # RK4 (balanced) ← DEFAULT
integration_order: 5  # RKF45 (adaptive, highest quality)
```

### Seed Density
```yaml
seed_density: 2  # Sparse (fast exploration)
seed_density: 4  # Balanced ← DEFAULT
seed_density: 6  # Dense (complete coverage)
seed_density: 8  # Very dense (critical regions)
```

### Termination FA
```yaml
termination_fa: 0.05  # Continue in very low FA (long tracks)
termination_fa: 0.15  # Balanced ← DEFAULT
termination_fa: 0.20  # Stop earlier (short tracks, faster)
```

### Step Size
```yaml
step_size: 0.1   # Precise (slow)
step_size: 0.2   # Balanced ← DEFAULT
step_size: 0.3   # Fast
```

## Parameter Validation

Configs are automatically validated when loaded:
- ✓ Step size must be positive
- ✓ Angle threshold must be in (0, 180]
- ✓ Integration order must be {1, 2, 4, 5}
- ✓ RKF parameters validated if adaptive_step enabled

## Documentation

See [docs/YAML_CONFIG.md](../docs/YAML_CONFIG.md) for complete reference.

## Examples

### High-precision, dense seeding
```yaml
tractography:
  integration_order: 5
  adaptive_step: true
  rkf_tolerance: 0.005
  seed_density: 6
  termination_fa: 0.1
```

### Fast exploration, sparse seeding
```yaml
tractography:
  integration_order: 2
  step_size: 0.3
  seed_density: 2
  termination_fa: 0.2
```

### Low-anisotropy regions (fornix, cingulum)
```yaml
tractography:
  integration_order: 4
  step_size: 0.15
  seed_density: 6
  termination_fa: 0.05
  angle_thresh: 45
```
