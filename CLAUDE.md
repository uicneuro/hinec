# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HINEC (HIgh-order NEural Connectivity) is a MATLAB pipeline for processing diffusion-weighted MRI (dMRI) and running fiber tractography. Inputs are raw NIfTI dMRI + bval/bvec; outputs are a processed `nim` struct (`.mat`) plus tractography tracks.

## Key Commands

### Primary entry point (shell)

The supported workflow is the shell launcher, which loads a YAML config, creates a timestamped run directory, and chains preprocessing + tractography in one `matlab -batch` background process:

```bash
./bin/run_hinec.sh <data_prefix> <output_mat> [config_file]
# e.g.
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/high_precision.yml
```

`<data_prefix>` is the shared path without extension. Expected files: `<prefix>_raw.nii.gz`, `<prefix>.bval`, `<prefix>.bvec`, optionally `<prefix>_T1.nii.gz`.

Other shell wrappers:
- `bin/run_tractography.sh <input.mat> [standard|hinec]` — tractography only on an already-processed `.mat`. Also supports `IronTract <injection.nii.gz> <output_dir>`.
- `bin/run_visualization.sh <run_dir> [format] [region] [dpi]` — headless figure export to `<run_dir>/figures/`.
- `bin/run_generateSlices.sh`, `bin/viewSlices.sh` — slice cache generation and Python viewer.

### MATLAB entry points

```matlab
config   = load_config_yaml('config/hinec_default.yml');
run_info = create_run_directory('config/hinec_default.yml');
main(data_prefix, 'output.mat', config, run_info);          % preprocessing + DTI + FA + parcellation
runTractography(fullfile(run_info.output_dir,'output.mat'), config, run_info);
```

Both `main` and `runTractography` accept legacy positional arg forms (no config / no `run_info`) — see the arg-parsing blocks at the top of each file.

### Tests

MATLAB unit/integration tests live under `tests/` (`unit/`, `integration/`, `fsl/`, `nim_tests/`). Run a single test class with MATLAB's runner, e.g.:

```matlab
results = runtests('tests/unit/TestNimFa.m');
```

`tests/test_yaml_config.m` is a standalone script that exercises every YAML preset in `config/`. There is no top-level `make test` target — the `Makefile` only drives `mkdocs`.

### Docs site

`make serve` (mkdocs) at `127.0.0.1:8000` — see `mkdocs.yml` and `docs/`.

## Architecture

### Source layout (everything is under `src/`)

- `src/nim_preprocessing/` — FSL-driven preprocessing (denoise, motion/eddy, brain extraction, T1↔DWI↔MNI registration, tissue segmentation for ACT).
- `src/nim_calculation/` — `nim_dt_spd` (SPD-constrained tensor fit), `nim_eig`, `nim_fa`.
- `src/nim_parcellation/` — atlas-based parcellation; `_registered` variant uses the registration transforms.
- `src/nim_registration/` — DTI↔T1↔MNI registration (FSL or SPM backend).
- `src/nim_tractography/` — `nim_tractography_standard` (FACT) and `nim_tractography_hinec` (high-order: interpolation + RK4/RKF45 + ACT). `nim_tractography_highorder` is the shared kernel.
- `src/nim_challenges/nim_irontract_submit.m` — IronTract Challenge submission packager.
- `src/nim_visualization/` — 3D viewer (`visualizeTractography`), interactive slice viewer (`visualizeTractographySlices`), and the server-side cache pipeline (`generateSlices` → `generateTractographySliceCache` → `TractographyCacheManager`) feeding the Python `FastTractographyViewer`.
- `src/nim_plots/nim_plot.m` — consolidated plotter. Use `nim_plot(nim, 'mode', 'parcels'|'grid'|...)`. (The old `nim_plotall` / `nim_plotparcelall` / `nim_plotparcellation` names are gone.)
- `src/nim_utils/` — I/O (`nim_read`, `nim_save`, `nim_load_nim`, `nim_load_labels`), `load_config_yaml`, `create_run_directory`, interpolation, GLL/uniform quadrature.

### Pipeline data flow

1. **`main.m`** orchestrates: detects raw vs. preprocessed → preprocessing (if raw) → `nim_read` → `nim_dt_spd` → `nim_eig` → `nim_fa` → optional registration → parcellation → mask improvement → tissue segmentation (WM/GM/CSF) for ACT → `nim_save`.
2. **`runTractography.m`** loads the saved `nim`, dispatches to `nim_tractography_standard` or `nim_tractography_hinec` per config, and writes tracks under `<run_dir>/tractography/`.
3. `main.m` short-circuits aggressively: if the output `.mat` exists it skips everything; if the preprocessed `.nii.gz` exists it skips preprocessing. To force reprocessing, delete the target file. This behavior is intentional — do not "fix" it by always reprocessing.

### Run directory system (`hinec_runs/run_YYYYMMDD_HHMMSS_<config>/`)

`create_run_directory` (and the shell launcher) build a self-contained run with `config.yml`, `logs/`, `intermediate/`, `output/`, `tractography/`, and (after visualization) `figures/`. When a `run_info` struct is passed through, both `main` and `runTractography` redirect their writes into it rather than scattering files alongside the input. Code paths that touch outputs almost always branch on `use_run_dir = ~isempty(fieldnames(run_info))` — preserve that branch when modifying.

### YAML configuration

`config/*.yml` drives everything. Two top-level keys: `preprocessing:` and `tractography:`. Presets in repo: `hinec_default`, `high_precision` (RKF45 adaptive), `fast_exploration`, `euler_interpolated`, `standard_fact`, `irontract`. `load_config_yaml` validates parameter ranges; prefer adding/changing knobs there over hardcoding in MATLAB.

### `nim` struct (key fields)

`.img` (4D DWI), `.evec`/`.eval` (per-voxel eigendecomp), `.FA`, `.mask`, `.wm_mask`/`.gm_mask`/`.csf_mask` (for ACT), `.parcellation_mask`, `.labels`, `.parcellation_mask_file`, optional `.registration`, `.run_info`.

### Track data structure

`tracks` is a cell array. Each cell is an `N×3` matrix of voxel-space coordinates representing the complete fiber trajectory (not just endpoints); N varies per fiber.

## Dependencies

- MATLAB toolboxes: Image Processing, Statistics and Machine Learning, "Tools for NIfTI and ANALYZE image".
- **SPM12**: expected at `lib/spm12/` (added via `addpath(genpath('lib/spm12'))` in `main.m`). Not vendored in this repo — install it there.
- **FSL**: required for preprocessing; must be initialized in the shell before `matlab -batch` runs.
- **`lib/bfgs/`** is vendored (BFGS solver used by tensor fitting).
- Python: `scripts/FastTractographyViewer.py`, `scripts/tractography_slice_gui.py`, `scripts/hinec_to_trk.py`, `scripts/validate_ismrm_tractography.py`. `requirements.txt` covers these.

## Things that bite

- **Path setup**: `main.m` does its own `addpath('src/nim_*')` and `addpath(genpath('lib/spm12'))`. Scripts run outside `main` need to call `addpath(genpath('.'))` or the equivalent — the shell launchers already do this.
- **Old vs. new function names**: the unified `nim_plot` replaced `nim_plotall`/`nim_plotparcelall`/`nim_plotparcellation`. Old call sites in notes/docs are stale.
- **Tractography seeding** is uniform grid-based (`seed_strategy: uniform`, `seed_density` per voxel) — the historical "FA-threshold seeding needs replacing" recommendation has already been done.
- **ACT** requires WM/GM/CSF masks produced by `preproc_tissue_segmentation` during `main.m`. If those masks are missing, tractography falls back to unconstrained.
- **String vs. char**: `main.m` aggressively `char(...)`-converts inputs because mixing MATLAB `string` and `char` types broke earlier file-path code (`docs/STRING_CHAR_FIXES.md`). Keep this when editing.

## Reference docs

`docs/` contains topical guides — consult these before re-deriving behavior:
`ARCHITECTURE.md`, `PIPELINE.md`, `PREPROCESSING.md`, `TRACTOGRAPHY.md`, `High_Order.md`, `RKF.md`, `SEEDING_STRATEGY.md`, `YAML_CONFIG.md`, `RUN_DIRECTORY_SYSTEM.md`, `DISTRIBUTED_WORKFLOW.md`, `VISUALIZATION_GUIDE.md`, `IRONTRACT_WORKFLOW.md`, `USER_GUIDE.md`.
