# ISMRM 2015 Tractography Challenge — data layout (live)

This directory holds the ISMRM 2015 dataset and both flavors of scoring data.
On this machine the layout is already populated — see "Current state" below.

## Current state

```
data/ismrm2015/
├── challenge_data/
│   └── ISMRM_2015_Tracto_challenge_data/
│       ├── Diffusion.nii.gz          ← raw DWI
│       ├── Diffusion.bvals
│       ├── Diffusion.bvecs
│       ├── T1.nii.gz                 ← T1 anatomical
│       ├── fmap_Hz.nii.gz
│       ├── fmap_RadPerSec.nii.gz
│       └── README.txt
├── scoring_data_Renauld2023/         ← Renauld 2023 scoring (scilpy)
│   ├── bundles/                      ← 25 ground-truth TRK bundles
│   ├── ROI/                          ← all_masks / any_masks / endpoints
│   ├── config_file_segmentation.json
│   ├── config_file_tractometry.json
│   └── t1.nii.gz                     ← 1 mm reference
├── scoring_data_2015/                ← Original 2015 scoring (dedicated scorer)
│   ├── bundles/
│   ├── masks/{rois,bundles}/
│   └── gt_bundles_attributes.json
├── ISMRM_2015_Tracto_challenge_data_v1_1.zip
├── scoring_data_Renauld2023.zip
└── scoring_data_tractography_challenge.tar.gz
```

## Rename for HINEC (one-time)

HINEC's `main.m` and `bin/run_hinec.sh` expect files with a shared prefix:
`<prefix>_raw.nii.gz`, `<prefix>.bval`, `<prefix>.bvec`, optionally `<prefix>_T1.nii.gz`.
The challenge uses `.bvals` / `.bvecs` and lacks the `_raw` suffix, so the
canonical layout for `./bin/run_hinec.sh data/ismrm2015/ismrm2015 …` is:

```bash
cd data/ismrm2015
src=challenge_data/ISMRM_2015_Tracto_challenge_data
ln -sf "$src/Diffusion.nii.gz" ismrm2015_raw.nii.gz
ln -sf "$src/Diffusion.bvals"  ismrm2015.bval
ln -sf "$src/Diffusion.bvecs"  ismrm2015.bvec
ln -sf "$src/T1.nii.gz"        ismrm2015_T1.nii.gz
```

(Symlinks — saves disk and stays bit-identical to the originals.)

## Run HINEC (when MATLAB + FSL are available)

```bash
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml
# or for comparison runs:
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/high_precision.yml
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/hinec_default.yml
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/fast_exploration.yml
```

## Score against ISMRM ground truth

After a run finishes (output in `hinec_runs/run_<TIMESTAMP>_<config>/`):

```bash
./bin/run_ismrm_scoring.sh hinec_runs/run_*_ismrm2015/
```

Writes Renauld 2023 metrics to `<run_dir>/scoring/renauld2023/results.json`
and dedicated-2015 metrics to `<run_dir>/scoring/dedicated2015/scores/`.

## Acquisition notes

- DWI: 2 mm isotropic, 32 directions, b = 1000 s/mm².
- Scoring ROI masks live in **T1 (1 mm) space**, not DWI space. The
  `config/ismrm2015.yml` preset turns on T1 registration so the affine
  saved with the NIfTI maps voxels → RAS world correctly. Without it,
  the scorer rejects most tracks because they land in the wrong space.
- The challenge data page warns about a **MRtrix gradient-vector orientation**
  issue. If you preprocess with `dwidenoise`, sanity-check bvec signs.

## What's installed on this machine

- `~/venvs/ismrm/` — Python 3.11 venv with scilpy, nibabel, dipy, scipy, h5py,
  and the dedicated `ismrm_2015_tractography_challenge_scoring` package.
- `~/fsl/` — FSL 6.0.7.16 local install (when complete; check
  `~/installers/fslinstall.log`). Activate with `export FSLDIR=~/fsl` +
  `source $FSLDIR/etc/fslconf/fsl.sh`.
- `~/ismrm_2015_tractography_challenge_scoring/` — cloned source for the
  dedicated 2015 scorer.
- **Not installed**: MATLAB (impossible without license + interactive installer),
  SPM12 (needs MATLAB), MRtrix3 (`dwidenoise` — pending; can install via the
  micromamba that FSL brought in once FSL finishes).
