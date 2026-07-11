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
./bin/run_hinec.sh data/parcellation_sample/sample sample.mat config/hinec_default.yml
```

`<data_prefix>` is the shared path without extension. Expected files: `<prefix>_raw.nii.gz`, `<prefix>.bval`, `<prefix>.bvec`, optionally `<prefix>_T1.nii.gz`.

Other shell wrappers:
- `bin/run_tractography.sh <config> [--score] [--set k=v ...] [--source <nim|dir>]` — **config-driven, tractography-ONLY** iteration on already-preprocessed data. Companion to `run_hinec.sh`: preprocess ONCE, then run any number of tractography configs on the processed nim WITHOUT re-preprocessing. `<config>` = `config/<name>.yml` or bare `<name>`. **Source defaults to the canonical nim in the data dir (`data/ismrm2015/ismrm2015.mat`)** — `main.m` writes it there — so you never point at a run dir; override with `--source`. It COPIES (freezes) the DWI reference into the new `run_<ts>_<config>/intermediate/` + writes `intermediate/SOURCE.txt` (so a later reprocess can't corrupt an old run's scoring reference), tracks, and with `--score` runs `run_ismrm_scoring.sh`. **`--set key=value` (DSE)** overrides any `config.tractography` param on the CLI without a new YAML (e.g. `--set sel_power=32 --set seed_density=8`); the run dir name is tagged with the overrides and `overrides.txt` records them. **This is the sanctioned way to run tractography experiments — write/copy a config and run this, or sweep with `--set`. Do NOT hand-write ad-hoc MATLAB batch scripts or hand-create run dirs.** Legacy `<input.mat> standard|hinec` and `IronTract <inj> <out>` forms still work.
- `bin/run_ismrm_scoring.sh <run_dir>` — the single scoring script: converts `tractography/tracks*.mat` → TRK (RAS world via the `intermediate/` DWI affine, with a `data/ismrm2015/ismrm2015.nii.gz` fallback) and runs the scilpy Renauld-2023 scorer → `<run_dir>/scoring/renauld2023/results.json` (headline keys `mean_f1`, `mean_OL`, `mean_OR_gt`, `VB`).
- `bin/run_visualization.sh <run_dir> [format] [region] [dpi]` — headless figure export to `<run_dir>/figures/`.
- `bin/run_generateSlices.sh`, `bin/viewSlices.sh` — slice cache generation and Python viewer.

**Config-driven experiment workflow (READ THIS before running experiments):**
1. Preprocess once: `./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml`. `main.m` writes the canonical processed nim to the **data layer** (`data/ismrm2015/ismrm2015.mat`) plus a copy in the run dir. (It reuses `data/ismrm2015/`'s preprocessed DWI if present, so repeats are cheap.)
2. Iterate tractography (no re-preprocess, auto-sourced from the data-layer nim):
   - single config: `./bin/run_tractography.sh hinec_csd --score`
   - DSE sweep without new YAMLs: `for s in 8 16 32; do ./bin/run_tractography.sh hinec_csd --score --set sel_power=$s; done`
   Each run → its own tagged, scorable run dir; compare by `scoring/renauld2023/results.json`.
- **Data layer vs run dirs:** the nim + preprocessed refs + CSD cache (`data/ismrm2015/ismrm2015_csd.mat`) live in `data/ismrm2015/`; `hinec_runs/run_*/` hold tractography OUTPUTS only. Never store nims inside run dirs.
- **Only committed `bin/*.sh` are the real pipeline.** Ignore/never create ad-hoc launchers like `run_one_config.sh`/`run_all_ismrm.sh` (prior untracked pollution that re-preprocessed every run — the exact thing the run_hinec/run_tractography split exists to avoid).

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
- `src/nim_tractography/` — three trackers, one shared nim:
  - `nim_tractography_standard` — FACT (discrete voxel tensors, no interpolation).
  - `nim_tractography_hinec` — **interpolated streamline tractography**: trajectory-dependent `sel_power` directional interpolation, `field: dti` (principal eigenvector) OR `field: csd` (FOD peaks, crossing resolution by approach direction), RK2/RK4/RKF45 (`integration_order`), ACT. This is the interpolated tracker (what was previously mislabelled `mmf` + `integrator: rkf45`). Configs: `config/hinec_dti.yml`, `config/hinec_csd.yml`.
  - `nim_tractography_mmf_connframe` — **GENUINE Method of Moving Frames** (Chun-Peng, `DiscussionWithPeng_Winter2026.pdf`): `nim_mmf_geometry` builds the moving-frame field + connection 1-form (curvature+torsion; per-peak for CSD → multiple pathways) into the nim (`main.m` Step 2b, Eq 6-9); the tracer evolves the carried frame by the structure equation and advances `dx/ds=e1` (Eq 10-11). Selected by `algorithm: 'mmf'`; the direction comes purely from the connection form, while `integrator` (`rk4` fixed | `rkf45` adaptive Dormand-Prince) is the numerical stepping scheme. Options: `field` (`dti`|`csd`), `integrator`, `mmf_anchor` (0 = pure Eq.10-11), `frame_sel_power`. Configs: `config/mmf_dti.yml`, `config/mmf_csd.yml`.
  - Dispatch by `algorithm:` (`standard`|`hinec`|`mmf`). `field: csd` requires FOD peaks — `runTractography` computes them via `nim_csd`, cached as `<source>_csd.mat` (`csd_lmax`/`csd_max_peaks`/`csd_peak_thresh`/`csd_peak_min_sep`).
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

`config/*.yml` drives everything. Two top-level keys: `preprocessing:` and `tractography:`. Configs follow one naming rule per family (see `config/README.md`): **tracker** configs are `<algorithm>_<field>[_<variant>].yml` — `standard_dti`, `hinec_dti`, `hinec_csd`, `mmf_dti`, `mmf_csd`, plus reusable hinec DTI variants `hinec_dti_cubic` (cubic+RKF45), `hinec_dti_cubic_recall` (cubic, recall-tuned), `hinec_dti_euler` (order 1), `hinec_dti_fast` (RK2) — all for `run_tractography.sh`; **dataset** configs are `<dataset>[_variant].yml` — `ismrm2015`, `ismrm2015_t1reg`, `irontract` (for `run_hinec.sh`); plus `hinec_default` (the generic full-pipeline fallback). The three orthogonal axes are `algorithm` (which tracker) / `field` (`dti`|`csd`, direction source) / integrator (stepping scheme: `integration_order` 1/2/4/5 for hinec/standard, `integrator` rk4/rkf45 for mmf). For hinec, direction is further shaped by two *composable* factors: `interp_method` (spatial kernel: `trilinear`|`cubic`) × `sel_power` (directional steering; 0 = plain interp) — cubic and sel_power combine, they are not mutually exclusive. **This is a research tool: keep useful experimental configs.** For a reusable knob-set make a `_<variant>` config; for a one-off sweep use `run_tractography.sh --set key=value` (don't delete existing configs just because a knob is expressible via `--set`). Settings irrelevant to the chosen combination are ignored by design. `load_config_yaml` validates parameter ranges; prefer adding/changing knobs there over hardcoding in MATLAB.

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
