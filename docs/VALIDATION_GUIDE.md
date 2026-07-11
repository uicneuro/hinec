# Validating Tractography Against the ISMRM 2015 Ground Truth

This guide is a step-by-step runbook for **quantitatively testing** a HINEC
tractography result against a known ground truth. It uses the ISMRM 2015
Tractography Challenge phantom — a synthetic brain where the true white-matter
bundles are known — and the Renauld et al. 2023 scoring system (via `scilpy`).

The end result is a `results.json` with per-bundle **Dice / overlap / overreach**
scores and the number of **valid bundles** recovered.

---

## 0. What you get

For every run you score, you get:

| metric | meaning |
|---|---|
| **Dice / F1** | overall agreement with the true bundle (0–1); `2·OL/(OL+OR+1)` |
| **OL** (overlap) | fraction of the true bundle you covered (recall) |
| **OR** (overreach) | fraction of your streamlines outside the true bundle |
| **VB** | how many of the 25 ground-truth bundles you recovered |

Reference points on this exact scoring (for calibration): best-ever submission
≈ 0.64 Dice, best 2015 original ≈ 0.58, mean of 96 submissions ≈ 0.41.

---

## 1. One-time environment setup

You need these available (install locations are suggestions; scripts default to
`$HOME/...` and can be overridden):

| tool | why | notes |
|---|---|---|
| **MATLAB** + Image Processing & Statistics toolboxes | runs HINEC | R2024a+; Parallel Computing Toolbox optional (speeds up tracking) |
| **FSL** | preprocessing (bet, eddy, flirt) | e.g. `~/fsl`; `source $FSLDIR/etc/fslconf/fsl.sh` |
| **SPM12** | I/O in the MATLAB pipeline | at `lib/spm12/` |
| **Python venv** | scoring | `~/venvs/ismrm` with `scilpy`, `dipy`, `nibabel` |

Create the Python scoring environment once:

```bash
python3.11 -m venv ~/venvs/ismrm        # scilpy needs Python 3.10/3.11
~/venvs/ismrm/bin/pip install scilpy dipy nibabel numpy scipy h5py
```

---

## 2. Get the data (one-time)

Full details in [`data/ismrm2015/README.md`](../data/ismrm2015/README.md).
In short, download and unpack into `data/ismrm2015/`:

- **Challenge DWI** — `ISMRM_2015_Tracto_challenge_data_v1_1.zip`
- **Scoring data (Renauld 2023)** — `scoring_data_Renauld2023.zip`

Then symlink the DWI into HINEC's `<prefix>_raw.nii.gz` naming:

```bash
cd data/ismrm2015
src=challenge_data/ISMRM_2015_Tracto_challenge_data
ln -sf "$src/Diffusion.nii.gz" ismrm2015_raw.nii.gz
ln -sf "$src/Diffusion.bvals"  ismrm2015.bval
ln -sf "$src/Diffusion.bvecs"  ismrm2015.bvec
ln -sf "$src/T1.nii.gz"        ismrm2015_T1.nii.gz
```

### Build the scoring config (one-time)

The scorer needs a config that lists each bundle's `gt_mask`. Merge the two
provided configs into one:

```bash
~/venvs/ismrm/bin/python scripts/build_ismrm_scoring_config.py
# -> writes data/ismrm2015/scoring_data_Renauld2023/config_file_merged.json
```

Without the merged config you get bundle counts but not Dice/overlap.

---

## 3. Run HINEC

```bash
# <data_prefix> <output.mat> <config>
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml
```

Output lands in a timestamped `hinec_runs/run_<TIMESTAMP>_ismrm2015/`:
- `output/ismrm2015.mat` — the processed `nim` struct
- `tractography/tracks_hinec_*.mat` — the streamlines (voxel coords)

`config/ismrm2015.yml` and `config/hinec_dti_cubic_recall.yml` are the two presets tuned
for this phantom (hinec_dti_cubic_recall scored highest in our tests).

> **Tip:** on a shared machine, cap the parallel pool with
> `HINEC_MAX_WORKERS=8` so tracking doesn't oversubscribe cores.

---

## 4. Score against ground truth

```bash
./bin/run_ismrm_scoring.sh hinec_runs/run_<TIMESTAMP>_ismrm2015/
```

This does three things:
1. Converts `tracks_*.mat` → `tracks.trk` in RAS world coordinates
   (`scripts/hinec_to_trk.py`, using the DWI affine, scoring-T1 reference).
2. Runs `scil_tractogram_segment_with_ROI_and_score` against the Renauld 2023
   ground truth.
3. Writes results under `<run_dir>/scoring/renauld2023/`.

Read the headline numbers:

```bash
~/venvs/ismrm/bin/python -c "import json; \
r=json.load(open('hinec_runs/run_<TIMESTAMP>_ismrm2015/scoring/renauld2023/results.json')); \
print('Dice=%.3f  OL=%.3f  OR=%.3f  VB=%.0f'%(r['mean_f1'],r['mean_OL'],r['mean_OR_gt'],r['VB']))"
```

---

## 5. Compare several runs/configs

```bash
~/venvs/ismrm/bin/python scripts/compare_ismrm_results.py \
    hinec_runs/run_*_ismrm2015 hinec_runs/run_*_hinec_dti_cubic_recall
```

Prints a side-by-side table of headline metrics and per-bundle Dice.

---

## 6. Interpreting the score

- **Higher Dice is better.** It balances *finding* the bundle (OL) against
  *staying inside* it (OR) — you can't game it by flooding streamlines.
- **Low OL, low OR** → under-coverage (too conservative); loosen termination
  (lower `termination_fa`, raise `angle_thresh`) or seed more densely.
- **High OL, high OR** → over-production; tighten the same knobs.
- **VB counts** which bundles you found at all; the small/deep ones
  (CA, CP, commissures) are the hardest.

---

## 7. Reproducing on your own data

The same flow validates any DWI that has a matching ground truth. Swap the
`data/ismrm2015/` inputs for your dataset (keeping the
`<prefix>_raw.nii.gz` / `.bval` / `.bvec` naming), point the scorer at your
ground-truth bundles/config, and the rest is identical.

---

## Troubleshooting

| symptom | cause / fix |
|---|---|
| scorer errors "incompatible headers" | tracks and scoring masks in different spaces — the scoring script attaches the scoring T1 reference; make sure you pass a run dir with an `intermediate/*.nii.gz` DWI |
| all bundles score 0 | usually a coordinate/space mismatch, or (for non-HINEC tools) a wrong gradient table — verify with `mrinfo`/`dwigradcheck` before trusting a low score |
| `results.json` has "No gt_mask" per bundle | you scored with `config_file_segmentation.json`; use the merged config (Step 2) for Dice |
| very short streamlines / low coverage | check `termination_fa`, `angle_thresh`, `min_length` in the config |
