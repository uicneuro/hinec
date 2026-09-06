# Validating Tractography Against the ISMRM 2015 Ground Truth

This is a runbook for scoring a HINEC tractography result against a known ground
truth. It uses the ISMRM 2015 Tractography Challenge phantom — a synthetic brain
whose true white-matter bundles are known — and the Renauld et al. 2023 scoring
system (via `scilpy`). The output is a `results.json` carrying per-bundle Dice,
overlap and overreach, plus the number of valid bundles recovered.

!!! warning "Score a whole-brain tractogram, or the numbers mean nothing"
    The scorer takes **one whole-brain tractogram** and segments it into 26
    bundles using ROI gates (head + tail endpoint ROIs and an `all_mask`
    containment corridor). Its headline metrics are averages over all 26. A
    tractogram seeded from a single ROI can populate at most one of them, so
    `VB` caps at 1/26 and `mean_f1` is dominated by the 25 bundles that were
    never attempted. Such a run is a useful bundle-level diagnostic, but it is
    **not** a challenge score and must not be compared with published ones. See
    [ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md) for how a bundle is defined.

---

## 0. What you get

| metric | meaning |
|---|---|
| `mean_f1` | per-bundle Dice, averaged over the 26 bundles |
| `mean_OL` | overlap: fraction of the true bundle's voxels covered (recall) |
| `mean_OR_gt` | overreach: produced-but-wrong volume, as a fraction of the true bundle's volume |
| `VB` | valid bundles: how many of the 26 were recovered at all |

All three are computed on **voxel masks**, not streamline counts. Writing
$A$ for the ground-truth bundle mask and $B$ for the segmented result:

$$
\mathrm{OL} = \frac{|A \cap B|}{|A|}, \qquad
\mathrm{OR_{gt}} = \frac{|B \setminus A|}{|A|}, \qquad
\mathrm{OR_{vs}} = \frac{|B \setminus A|}{|B|}
$$

and `f1` is the harmonic mean of recall $\mathrm{OL}$ and precision
$1 - \mathrm{OR_{vs}}$. Note that `mean_OR_gt`, the reported overreach, is
normalised by the *ground-truth* volume and can therefore exceed 1, while the
overreach entering `f1` is normalised by the produced volume and cannot.
`results.json` also carries `VS`, `VS_ratio` and `IS_ratio` (valid and invalid
streamline counts and their ratios).

Published reference points for this scorer, for calibration: the best submission
recorded is ≈ 0.64 Dice, the best of the original 2015 entries ≈ 0.58, and the
mean over 96 submissions ≈ 0.41. These are whole-brain submissions; only a
whole-brain run of your own is comparable to them.

---

## 1. One-time environment setup

Install locations below are suggestions; the scripts default to `$HOME/...` and
can be overridden (`ISMRM_VENV` for the scoring venv).

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

The scorer needs a config that pairs each bundle's ROI gates with its `gt_mask`.
Merge the two provided configs into one:

```bash
~/venvs/ismrm/bin/python scripts/build_ismrm_scoring_config.py
# -> writes data/ismrm2015/scoring_data_Renauld2023/config_file_merged.json
```

Without the merged config you get bundle counts but no Dice or overlap.

---

## 3. Run HINEC

Preprocess once:

```bash
# <data_prefix> <output.mat> <config>
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml
```

This writes the canonical processed `nim` to `data/ismrm2015/ismrm2015.mat` and
a copy into a timestamped `hinec_runs/run_<TIMESTAMP>_ismrm2015/`, alongside
`tractography/tracks_hinec_*.mat` (streamlines in voxel coordinates).

Then iterate tractography on the processed `nim` without re-preprocessing, each
config landing in its own scorable run directory:

```bash
./bin/run_tractography.sh hinec_dti_cubic --score
./bin/run_tractography.sh hinec_dti_cubic_recall --score
```

`config/ismrm2015.yml` is the dataset preset for this phantom;
`config/hinec_dti_cubic_recall.yml` is a variant that deliberately trades
precision for coverage (denser seeding, looser termination). Which config scores
best has not been established on a whole-brain submission — see the warning at
the top of this page.

> **Tip:** on a shared machine, cap the parallel pool with
> `HINEC_MAX_WORKERS=8` so tracking doesn't oversubscribe cores.

---

## 4. Score against ground truth

```bash
./bin/run_ismrm_scoring.sh hinec_runs/run_<TIMESTAMP>_ismrm2015/
```

This does three things:

1. Converts the newest `tractography/tracks*.mat` to `scoring/tracks.trk` in RAS
   world millimetres (`scripts/hinec_to_trk.py`, using the run's preprocessed DWI
   affine for placement and the scoring T1 as the saved TRK reference).
2. Runs `scil_tractogram_segment_with_ROI_and_score` against the merged Renauld
   2023 config, writing `<run_dir>/scoring/renauld2023/results.json`.
3. If `data/ismrm2015/scoring_data_2015/` and the original challenge scorer are
   also present, runs the dedicated 2015 scorer into
   `<run_dir>/scoring/dedicated2015/` as a cross-check. This step is skipped
   silently when either is missing.

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

Prints a side-by-side table of headline metrics and per-bundle F1, labelled by
the config name embedded in each run directory.

---

## 6. Interpreting the score

- **Higher Dice is better.** It balances *finding* the bundle (OL) against
  *staying inside* it (OR), so flooding the brain with streamlines does not
  raise it.
- **Low OL, low OR** → under-coverage. Loosen termination (lower
  `termination.fa_min`, raise `termination.angle_max`) or raise
  `seeding.density`.
- **High OL, high OR** → over-production. Tighten the same knobs.
- **VB counts** which bundles were recovered at all; the small and deep ones
  (`CA`, `CP`, the commissures) are hardest — `CP` has 365 ground-truth
  streamlines against `MCP`'s 21008.
- **Read `VB` first.** A `VB` far below 26 means most bundles contributed a zero
  to `mean_f1`, and the average says more about coverage than about quality.

---

## 7. Reproducing on your own data

The same flow validates any DWI that has a matching ground truth. Swap the
`data/ismrm2015/` inputs for your dataset (keeping the
`<prefix>_raw.nii.gz` / `.bval` / `.bvec` naming), point the scorer at your
ground-truth bundles and config, and the rest is identical.

---

## Troubleshooting

| symptom | cause / fix |
|---|---|
| scorer errors "incompatible headers" | tracks and scoring masks are in different spaces — the scoring script attaches the scoring T1 reference, so make sure you pass a run dir with an `intermediate/*.nii.gz` DWI |
| all bundles score 0 | usually a coordinate/space mismatch, or (for non-HINEC tools) a wrong gradient table — verify with `mrinfo`/`dwigradcheck` before trusting a low score |
| `results.json` has "No gt_mask" per bundle | you scored with `config_file_segmentation.json`; use the merged config (Step 2) for Dice |
| `VB` is 1 and `mean_f1` is near zero | the tractogram was seeded from a single ROI — see the warning at the top |
| very short streamlines / low coverage | check `termination.fa_min`, `termination.angle_max` and `termination.min_arc` in the config |
