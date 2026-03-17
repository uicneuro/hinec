# HINEC: HIgh-order NEural Connectivity

Human brain white matter tractography pipeline with YAML-based parameter configuration.

This file is for testing mkdocs stuff.
![External Image Testing](https://pythonfix.com/pkg/m/mkdocs/mkdocs-banner.webp)
![In-Repository Image Testing](img/hinec_test_image.png)

## What is HINEC?

HINEC processes raw diffusion-weighted MRI (dMRI) data through a complete pipeline:

1. **Preprocessing** — motion/eddy correction via FSL
2. **Diffusion tensor estimation** — SPD-constrained fitting
3. **FA computation** — fractional anisotropy maps
4. **Parcellation** — atlas-based brain region segmentation
5. **Tractography** — high-order fiber tracking (FACT, RK4, RKF45)
6. **Visualization** — 3D and fast slice-based viewing

## Quick Start

```bash
# Default run — auto-creates organized run directory
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat

# High precision (publication quality)
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/high_precision.yml

# Fast exploration
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/fast_exploration.yml
```

Or from MATLAB:

```matlab
main('path/to/data', 'output.mat')
runhinec
```

## Navigation

| Section | Description |
|---|---|
| [Pipeline](PIPELINE.md) | End-to-end data flow and module overview |
| [Tractography](TRACTOGRAPHY.md) | Standard and high-order fiber tracking |
| [Configuration](YAML_CONFIG.md) | YAML parameter system |
| [Visualization](VISUALIZATION_GUIDE.md) | 3D and distributed slice viewing |
| [API Reference](API_REFERENCE.md) | Complete function reference |
| [Math Foundations](MATHEMATICAL_FOUNDATIONS.md) | Formulas and numerical methods |
