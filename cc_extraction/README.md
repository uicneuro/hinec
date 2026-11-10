# Corpus-callosum tractography from HINEC (`data/original_sample`)

This folder contains a corpus-callosum (CC) extraction produced by running the HINEC pipeline on
`data/original_sample/sample` with `config/hinec_default.yml`, plus everything needed to reproduce it.

Because neither MATLAB nor FSL/SPM was available in the environment where this was run, the pipeline was executed
through a **faithful Python port of the HINEC MATLAB code** (`hinec_py/`), and a **MATLAB script that runs the real
HINEC kernels without FSL/SPM** (`extract_cc_tracks.m`) is provided for reproduction in MATLAB R2026a.
The run date is 2026-08-25.

## Results at a glance (HINEC as written)

| stage | result |
|---|---|
| data | 128 × 128 × 70 voxels, 2 mm isotropic, 1 b0 + 30 DWI at b = 1000 s/mm² (LAS orientation) |
| brain mask | 169,606 voxels (median-Otsu on the b0; FSL BET was unavailable) |
| tensor fit | 169,603 voxels, 369 needed an SPD refit (0.2 %), 471 zeroed as invalid |
| V1 check | median angle to the shipped FSL `sample_V1.nii.gz`: **0.00°**, 99.9 % of FA > 0.3 voxels within 10° |
| tissue masks (ACT) | WM threshold FA > 0.231, CSF threshold FA ≤ 0.098 (adaptive percentiles, as in the MATLAB code) |
| seeding | 163,467 seed voxels × 8 seeds (2×2×2 lattice; `seed_density: 4` → `ceil(4^(1/3)) = 2` per axis) = 1,307,736 seeds |
| tracking | RK4, step 0.2 voxel, cubic interpolation, angle 35°, termination FA 0.2, ACT on — 635,255 valid tracks in 9.0 min |
| CC ROI | midsagittal slab (±1 voxel), FA > 0.4, \|v1·LR\| > 0.8, largest component: 343 voxels |
| **corpus callosum** | **23,796 tracks** pass through the ROI and cross the midline (24,870 hit the ROI); 17,956 are ≥ 45 mm |
| CC lengths | median 54.5 mm, 10–90 % range 33–83 mm, max 204 mm |

Figures: `results/figures/cc_tracks_projections.png` (sagittal / coronal / axial over FA), `cc_tracks_3d.png`,
`wholebrain_vs_cc_axial.png`, `cc_lengths_and_terminations.png`, and QC images for the brain mask, tissue masks and
the ROI on the colour-FA map.

## Files

```
results/                         faithful run (HINEC as written)
  tracks_cc.trk                  CC streamlines, TrackVis format, RAS mm  (MRtrix / DSI Studio / TrackVis / dipy)
  tracks_cc.mat                  same tracks in HINEC's own layout: cell array of N×3 1-based voxel coordinates
  tracks_cc_min45mm.trk          subset ≥ 45 mm (the config's min_length — see note 5 below)
  tracks_wholebrain_subsample.*  random 2 % of all whole-brain tracks, for context
  nim_hinecpy.mat                nim struct (FA, evec, eval, mask, WM/GM/CSF masks, parcellation_mask = CC ROI, labels)
  cc_roi.nii.gz                  the CC ROI
  sample_M.nii.gz                brain mask  — use it as <prefix>_M.nii.gz / 'MaskPath' for nim_read
  sample_{WM,GM,CSF}_mask.nii.gz tissue masks used for ACT
  sample_FA_hinecpy.nii.gz, sample_V1_hinecpy.nii.gz   FA and primary eigenvector from the port
  summary.json                   all counts, thresholds, termination statistics, timings
results_sign_aligned/            comparison run with sign-aligned interpolation (see "Findings")
comparison_raw_vs_sign_aligned.{png,json}
extract_cc_tracks.m              MATLAB reproduction with the real HINEC functions (no FSL/SPM needed)
hinec_py/, run_cc_extraction.py, make_figures.py, compare_runs.py, tests/   the Python port
```

### Using the results in MATLAB (HINEC's own viewer)

```matlab
addpath(genpath('.'));                                   % from the repo root
visualizeTractography('cc_extraction/results/tracks_cc.mat', 'cc_extraction/results/nim_hinecpy.mat', ...
                      'export_dir', 'cc_extraction/results/figures_matlab');
visualizeTractography('cc_extraction/results/tracks_cc.mat', 'cc_extraction/results/nim_hinecpy.mat', ...
                      'region', 1, 'min_region_points', 1);   % region 1 = the CC ROI
load('cc_extraction/results/tracks_cc.mat', 'tracks');   % tracks{i} is N×3, 1-based voxel coordinates
```

### Reproducing with the real HINEC code in MATLAB

```matlab
addpath(genpath('.'));
result = extract_cc_tracks('data/original_sample/sample', 'config/hinec_default.yml', ...
                           'cc_extraction/matlab', 'mask_file', 'cc_extraction/results/sample_M.nii.gz');
```

The script calls `load_config_yaml → nim_read → nim_dt_spd → nim_eig → nim_fa → nim_tractography_hinec` exactly as
`main.m`/`runTractography.m` do, replaces the three FSL/SPM steps (BET, `fslmaths -ero`, atlas registration) with the
supplied mask, MATLAB `imerode`, and the midsagittal ROI, and then applies the same CC selection rule.
**It has not been executed** (no MATLAB here) — please treat it as untested until the first run. Expect a longer runtime
than the Python port: `nim_dt_spd` loops over voxels and `nim_tractography_hinec` runs 1.3 M seeds through `parfor`
(set `HINEC_MAX_WORKERS` if you want more than 8 workers).

### Re-running the Python port

```bash
pip install numpy scipy nibabel dipy matplotlib
python cc_extraction/run_cc_extraction.py --prefix data/original_sample/sample \
       --config config/hinec_default.yml --out cc_extraction/results          # ~10 min on 2 cores
python cc_extraction/run_cc_extraction.py ... --align-signs                    # the comparison run
python cc_extraction/make_figures.py --out cc_extraction/results --nim <nim pickle written with --stage-cache>
python cc_extraction/tests/test_scalar_reference.py                            # vectorised tracker vs literal transcription
```

## What the port does, and where it deviates from the MATLAB code

The port mirrors the MATLAB control flow function by function (`dti.py` ↔ `nim_read/nim_dt_spd/nim_eig/nim_fa/
preproc_mask_improvement/preproc_tissue_segmentation`, `tracking.py` ↔ `nim_tractography_hinec` + the seeding block
of `runTractography.m`). The vectorised tracker was checked against a literal per-streamline transcription of
`track_fiber_hinec` on 400 half-tracks: identical point counts and termination reasons, coordinates equal to float32
precision (`tests/test_scalar_reference.py`). Known deviations:

1. **Brain extraction** — `dipy.segment.mask.median_otsu` on the b0 instead of `bet -f 0.3`. The mask is 14.8 % of the
   image; because that is just under the 15 % floor in `preproc_mask_improvement.m`, its morphological clean-up step
   is skipped, exactly as the MATLAB code would skip it.
2. **SPD refit** — `lib/bfgs/vox_dt_bfgs.m` minimises the right objective but its gradient (`gradf`) uses `g(1), g(2),
   g(3)`, i.e. the x-components of the first three gradient directions, for every direction, so its BFGS iterations do not
   descend on the objective. The port minimises the same objective with the correct analytic gradient (L-BFGS-B on a
   Cholesky factor). It affects 369 of 169,603 voxels (0.2 %), essentially noise voxels.
3. **Cubic interpolation** — Keys cubic convolution (a = −½), which is what `griddedInterpolant(...,'cubic')` implements,
   with replicate padding at the volume edge (far from the brain).
4. **Eigenvector signs** — NumPy's `eigh` and MATLAB's `eig` pick eigenvector signs independently, so a voxel-by-voxel
   numerical match of MATLAB tracks is impossible in principle (see Finding 1); statistically the runs should agree.
5. **Nothing else was "fixed"**: `min_length` is not enforced (HINEC parses it but never uses it), a CSF hit discards
   only that half of a track (as in the MATLAB code), FA terminations are reported as `no_direction` (as in the MATLAB
   code, because `interpolate_direction_trilinear` returns an empty direction when FA < `termination_fa`).

## Findings about the HINEC code that came out of this run

1. **Eigenvector sign ambiguity dominates the tracking.** `nim_tractography_hinec` interpolates the three components
   of the primary eigenvector independently (`griddedInterpolant` on `v1_x, v1_y, v1_z`). Eigenvectors have arbitrary
   sign, and in this dataset 26 % of neighbouring white-matter voxel pairs (FA > 0.3) have opposite signs (18 % are pure
   flips of nearly parallel vectors). Interpolating across such a pair partially cancels the vector, which deflects the
   step and trips the 35° angle test. Aligning the eight (or 64) neighbouring eigenvectors to the incoming direction
   before weighting (`--align-signs`, `FieldInterpolator(align_signs=True)`) changes the outcome dramatically:

   | | HINEC as written | sign-aligned |
   |---|---|---|
   | whole-brain `angle` terminations | 445,982 | 116 |
   | whole-brain median track length | 6.4 mm | 11.2 mm |
   | CC tracks (ROI ∧ midline) | 23,796 | 42,382 |
   | CC median / 90th-pct length | 54.5 / 83 mm | 95.6 / 142 mm |

   With the fix, callosal fibres reach the cortex (`comparison_raw_vs_sign_aligned.png`); without it they stop about
   half-way. The same effect will be present in the MATLAB run because it is inherent to `eig`, not to NumPy.
   A three-line fix in `interpolate_direction_trilinear` (flip corner vectors whose dot product with `dir_vec` is
   negative before interpolating) — or interpolating the dyadic tensor v vᵀ — removes it.
2. **`min_length` is never applied** in either kernel (`nim_tractography_hinec.m`, `nim_tractography_standard.m`); the
   45 mm in `hinec_default.yml` has no effect. Two-point tracks are kept as "valid".
3. **`main.m` cannot run without FSL**, even for preprocessed data: `preproc_mask_improvement` and
   `preproc_tissue_segmentation` `error()` when `FSLDIR` is unset (the first is not inside a try/catch), and
   `nim_parcellation` shells out to `flirt`. On a machine without FSL the practical path is `extract_cc_tracks.m`.
4. **WM erosion for ACT** uses `fslmaths -ero` with its default 3×3×3 kernel; at 2 mm this leaves only 5 % of the brain
   labelled WM (the cores of the major tracts). Eroded-away voxels fall to `UNKNOWN`, which continues tracking, so the
   effective ACT rule is simply "stop in GM (FA ≤ 0.231), discard the half-track in CSF (FA ≤ 0.098)".
5. **`vox_dt_bfgs.m` gradient bug** (see deviation 2).

## CC definition

Atlas-free (the JHU labels need FSL registration): a slab of ±1 voxel around the brain-mask centre of mass along the
left–right axis (x = 63.55 in 0-based voxel units, i.e. 64.55 in MATLAB indexing), restricted to FA > 0.4 and
|v1 · LR| > 0.8, largest 26-connected component (343 voxels — the genu-body-splenium arch; the transverse pontine fibres
form a separate, smaller component and are excluded). A streamline is callosal if any of its points rounds into the
ROI **and** it has points on both sides of the midline. Thresholds are arguments of `cc.build_cc_roi` /
`extract_cc_tracks(..., 'roi_fa_thr', 'roi_lr_thr', 'roi_halfwidth')`.
