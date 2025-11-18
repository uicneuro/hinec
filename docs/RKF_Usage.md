# RKF45 Adaptive Tractography - Usage Guide

## Overview

The HINEC tractography pipeline now supports **RKF45 (Runge-Kutta-Fehlberg)** adaptive step size integration using the Dormand-Prince embedded RK pair. This provides superior accuracy with automatic error control compared to fixed-step methods (Euler, RK2, RK4).

## Key Features

### 1. **Embedded RK Pair**
- **5th-order solution**: Used for position advancement (high accuracy)
- **4th-order solution**: Used for error estimation
- **7-stage integration**: Single set of k_i evaluations for both solutions
- **Dormand-Prince coefficients**: Industry-standard RK5(4)7M method

### 2. **Adaptive Step Size Control**
- **Automatic error estimation**: `error = ||y_hat - y_tilde||`
- **Dynamic step adjustment**: `h_new = safety * h * (tol/err)^(1/5)`
- **Step bounds**: `[step_min, step_max]` = `[0.01, 1.0]` voxels
- **Safety factor**: 0.9 (conservative)
- **Growth limits**: Maximum 2× increase per step

### 3. **Integration Hierarchy**
```
Order 1 (Euler):    1st-order, fixed step
Order 2 (RK2):      2nd-order, fixed step
Order 4 (RK4):      4th-order, fixed step
Order 5 (RKF45):    5th-order, adaptive step  ← NEW
```

## Usage Examples

### Basic RKF45 with Default Settings

```matlab
% Load preprocessed data
main('data/subject_raw', 'output/subject.mat');

% Run RKF45 tractography with adaptive stepping
options = struct();
options.integration_order = 5;        % Enable RKF45
options.adaptive_step = true;         % Enable adaptive stepping
options.rkf_tolerance = 0.01;         % Error tolerance (voxels)

tracks = nim_tractography_hinec('output/subject.mat', options);
```

### High-Precision Tracking (Recommended)

```matlab
% For publication-quality results
options = struct();
options.integration_order = 5;
options.adaptive_step = true;
options.rkf_tolerance = 0.005;        % Tighter tolerance (higher accuracy)
options.step_size = 0.1;              % Initial step size (smaller = safer start)
options.step_min = 0.005;             % Minimum allowed step
options.step_max = 0.5;               % Maximum allowed step
options.rkf_safety = 0.85;            % More conservative safety factor

tracks = nim_tractography_hinec('data.mat', options);
```

### Fast Exploration Mode

```matlab
% For quick exploratory analysis
options = struct();
options.integration_order = 5;
options.adaptive_step = true;
options.rkf_tolerance = 0.02;         % Looser tolerance (faster)
options.step_size = 0.2;              % Standard initial step
options.seed_density = 2;             % Moderate density

tracks = nim_tractography_hinec('data.mat', options);
```

### Fixed-Step RKF45 (No Adaptivity)

```matlab
% Use RKF45 accuracy without adaptive stepping
options = struct();
options.integration_order = 5;
options.adaptive_step = false;        % Disable adaptivity
options.step_size = 0.2;              % Fixed step size

tracks = nim_tractography_hinec('data.mat', options);
```

## Parameter Guidelines

### Error Tolerance (`rkf_tolerance`)

| Tolerance | Use Case | Expected Behavior |
|-----------|----------|-------------------|
| 0.001 | Ultra-high precision | Very small steps, slow, highest accuracy |
| 0.005 | High precision | Small steps, moderate speed, excellent accuracy |
| **0.01** | **Recommended** | **Balanced performance and accuracy** |
| 0.02 | Fast exploration | Larger steps, faster, good accuracy |
| 0.05 | Quick screening | Large steps, fastest, acceptable accuracy |

### Initial Step Size (`step_size`)

| Step Size | Recommendation | Notes |
|-----------|---------------|-------|
| 0.05 | Conservative | Very safe, may be slow |
| **0.1-0.2** | **Recommended** | Good balance |
| 0.3-0.5 | Aggressive | Faster but may reject more steps |

### Safety Factor (`rkf_safety`)

| Safety | Behavior | Trade-off |
|--------|----------|-----------|
| 0.7-0.8 | Aggressive | Faster adaptation, more rejections |
| **0.9** | **Recommended** | Conservative, fewer rejections |
| 0.95 | Very safe | Slowest adaptation, safest |

## Performance Characteristics

### Computational Cost

**RKF45 vs RK4**:
- **Stages**: 7 vs 4 (1.75× more interpolations)
- **Error estimation**: Included (no extra cost)
- **Adaptive stepping**: May reduce total steps
- **Overall**: ~1.5-2× slower than RK4, but higher accuracy

**Expected Timing** (typical brain, 10K seeds):
- RK4 fixed: ~5-10 minutes
- RKF45 adaptive (tol=0.01): ~8-15 minutes
- RKF45 adaptive (tol=0.005): ~12-20 minutes

### Accuracy Gains

**Local Error Scaling**:
- Euler (order 1): Error ~ h²
- RK2 (order 2): Error ~ h³
- RK4 (order 4): Error ~ h⁵
- **RKF45 (order 5)**: Error ~ h⁶ ← Best

**Global Error** (accumulation over N steps):
- RK4: Error ~ N·h⁵
- RKF45: Error ~ N·h⁶ (superior)

### Step Rejection Statistics

**Normal Operation**:
- Rejection rate: 5-15% of steps
- Retry rate: 1-3% of steps
- Most rejected steps succeed on first retry

**High rejection rates** (>20%) indicate:
- Tolerance too tight for step size
- Complex direction field
- Consider: reduce tolerance or increase initial step_size

## Diagnostic Output

### Console Output Example

```
Starting HINEC High-Order Tractography...
Parameters: step=0.20, FA_thresh=0.10, angle_thresh=60.0°, integration_order=5
Integration: RKF45 (Dormand-Prince) (adaptive, tol=0.0100)
Interpolation: trilinear, ACT: true

=== HINEC TIMING REPORT ===
Integration method: Order 5 (RKF45 (Dormand-Prince))
Adaptive stepping: ENABLED (tolerance=0.0100 voxels)
Total time: 842.31 seconds
Eigenvector extraction: 12.45 seconds (1.5%)
Seed generation: 3.21 seconds (0.4%)
HINEC tracking: 826.65 seconds (98.1%)
  - Interpolation + integration overhead included
  - Brain boundary checks: 45.23 seconds (5.5% of tracking)
  - RKF step rejections: 1247 (8.3% of steps)
  - RKF retry attempts: 213
Total integration steps: 15023
Average steps per track: 75.1
Integration steps per second: 18.2
========================
```

### Interpretation

- **RKF step rejections**: 8.3% is normal (5-15% expected)
- **RKF retry attempts**: 213 retries for 1247 rejections = ~17% needed retry (most succeed immediately)
- **Steps per second**: 18.2 is reasonable for RKF45 with ACT

## Advanced Usage

### Combining with ACT

```matlab
% RKF45 + ACT for maximum biological plausibility
main('data/subject_raw', 'subject.mat');  % Generates tissue masks

options = struct();
options.integration_order = 5;
options.adaptive_step = true;
options.rkf_tolerance = 0.01;
% ACT automatically enabled if tissue masks available

tracks = nim_tractography_hinec('subject.mat', options);
```

### Varying Tolerance for Different Regions

```matlab
% Not directly supported - use multiple runs
% High precision for critical regions:
options_precise = struct();
options_precise.integration_order = 5;
options_precise.adaptive_step = true;
options_precise.rkf_tolerance = 0.005;
% ... restrict seed_mask to critical regions ...

% Standard precision for exploratory regions:
options_standard = struct();
options_standard.integration_order = 5;
options_standard.adaptive_step = true;
options_standard.rkf_tolerance = 0.02;
% ... different seed_mask ...
```

## Validation & Quality Assurance

### Recommended Validation

1. **Compare with Fixed RK4**:
```matlab
% Run both and compare track counts, lengths, coverage
tracks_rk4 = nim_tractography_hinec('data.mat', struct('integration_order', 4));
tracks_rkf = nim_tractography_hinec('data.mat', struct('integration_order', 5, 'adaptive_step', true));
```

2. **Check Statistics**:
- Similar track counts (±10%)
- RKF tracks may be slightly longer (more accurate integration)
- Lower rejection rate is better

3. **Visualize Differences**:
```matlab
visualizeTractography('tracks_rk4_*.mat', 'data.mat');
visualizeTractography('tracks_rkf_*.mat', 'data.mat');
```

### Known Limitations

1. **Adaptive stepping overhead**: 1.5-2× slower than fixed RK4
2. **Memory**: Same as RK4 (no additional memory required)
3. **Not recommended** for very coarse initial steps (>0.5 voxels)
4. **Best for**: High-quality final analysis, not rapid screening

## Troubleshooting

### Problem: High Rejection Rate (>25%)

**Diagnosis**: Tolerance too tight or initial step too large

**Solutions**:
```matlab
% Option 1: Relax tolerance
options.rkf_tolerance = 0.02;  % instead of 0.01

% Option 2: Reduce initial step
options.step_size = 0.1;  % instead of 0.2

% Option 3: Adjust safety factor
options.rkf_safety = 0.95;  % more conservative
```

### Problem: Very Slow Execution

**Diagnosis**: Tolerance too tight or step bounds too restrictive

**Solutions**:
```matlab
% Increase tolerance
options.rkf_tolerance = 0.015;

% Increase minimum step
options.step_min = 0.02;  % instead of 0.01

% Allow larger maximum step
options.step_max = 1.5;  % instead of 1.0
```

### Problem: Tracks Too Short

**Diagnosis**: Overly conservative stepping or termination

**Solutions**:
```matlab
% Increase step bounds
options.step_max = 2.0;

% Relax termination criteria
options.termination_fa = 0.03;  % instead of 0.05

% Increase max steps
options.max_steps = 10000;  % instead of 5000
```

## References

1. **Dormand, J. R., & Prince, P. J.** (1980). A family of embedded Runge-Kutta formulae. *Journal of Computational and Applied Mathematics*, 6(1), 19-26.

2. **Hairer, E., Norsett, S. P., & Wanner, G.** (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer.

3. **Basser, P. J., Pajevic, S., Pierpaoli, C., Duda, J., & Aldroubi, A.** (2000). In vivo fiber tractography using DT-MRI data. *Magnetic Resonance in Medicine*, 44(4), 625-632.

## Summary

RKF45 adaptive integration provides **superior numerical accuracy** with **minimal user intervention** through automatic error control. Recommended for:

✅ **Publication-quality results**
✅ **High-precision connectivity analysis**
✅ **Challenging fiber configurations**
✅ **Regions with complex geometry**

Use fixed-step RK4 for:
✅ **Rapid exploration**
✅ **Parameter screening**
✅ **Large-scale batch processing**

**Default recommendation**: `integration_order=5, adaptive_step=true, rkf_tolerance=0.01`
