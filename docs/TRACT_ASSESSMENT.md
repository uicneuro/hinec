# Tractogram assessment: convergence, repeatability, consistency

One pipeline answers three different questions about **any** HINEC tracker
(`standard`, `hinec`, `mmf`, `stitching`, `template`):

| Study | What is compared | Question it answers |
|---|---|---|
| **Convergence** | One input and one pipeline, run at several settings of a numerical knob (step size, seed density, …) from coarse to fine. Each level is compared with the finest level. | Does the output settle as the setting is refined? Does it settle for every streamline, or only for the typical one? |
| **Repeatability** | The same pipeline run on independent repeats of the same subject (scan and rescan). | How much of the output survives re-acquisition? |
| **Consistency** | Different methods or settings run on the same input (trilinear vs cubic, hinec vs mmf, DTI vs CSD). | How much does the output depend on the method you chose? |

All three use one comparison engine. Only the choice of which pairs to compare
changes, along with how the result is read.

## The four tiers of agreement

Every comparison of two tractograms produces one agreement card with four
tiers. The tiers are reported separately because in real data they disagree.
The tract-geometry study showed this directly: geometry statistics repeated to
within 2–5% while the streamlines themselves occupied only 44–55% of the same
voxels.

| Tier | Measures | Headline numbers |
|---|---|---|
| **Yield** | Did both runs produce streamlines from the seeds they tried? | tracks, seeds attempted, seed-region Dice |
| **Space** | Do they pass through, and end in, the same voxels? | occupied-voxel Jaccard, endpoint-voxel Jaccard, density-weighted Dice |
| **Paths** | Do corresponding streamlines follow the same route? | median separation (mm), worst-5% separation, fraction that diverge by >1 mm, arc distance to first divergence |
| **MMF geometry** | Do they have the same shape? | κ<sub>w</sub>, τ, w<sub>23</sub>: shift in the bundle median (distribution) and paired Spearman ρ (per streamline) |

**Only the Space tier tells you the tractograms are in the same place.** A
translated copy of a tractogram has identical geometry and zero spatial
overlap; `tests/python/test_tract_assessment.py` checks this case.

### Streamline correspondence

- **seed**: used when both runs record the seed of every streamline
  (`meta.seed_points`: `hinec`, `template`, `stitching`, and custom `.mat` files
  with `seed_points`). Streamlines with the same seed coordinate are paired.
  The pairing is exact. Each streamline is then split at its seed, and the two
  halves are compared at equal arc length outward from the seed. Nothing is
  assumed about step size or the integration variable, so runs with different
  `integrator.step` can be compared.
- **nearest**: used for `standard` and `mmf`, which do not record seeds. Each
  streamline is paired with the streamline closest to it by mean direct-flip
  distance. This is a heuristic ("the most similar streamline", not "the same
  seed"). The report always states which mode was used.

## The MMF geometry measurement

Each streamline is transformed to **world millimetres with the full NIfTI
affine**, resampled every 0.5 mm, and differentiated over a 4 mm window. A
moving frame is built along it: the Frenet normal where curvature is at least
0.001 mm⁻¹, Bishop transport elsewhere. Its connection coefficients are then
measured, `w_ij = <de_i/ds, e_j>`, with frames stored as rows. Each streamline
contributes the 90th percentile of |value|. The code is in
`scripts/tract_assessment/geometry.py`.

| Quantity | Meaning | Use |
|---|---|---|
| κ<sub>w</sub> = √(w₁₂² + w₁₃²) | curvature magnitude, independent of how the normal plane is oriented | **primary** |
| τ | Frenet torsion, only where curvature is valid | secondary |
| w₂₃ | frame twist in the declared hybrid gauge; it changes when the gauge changes | diagnostic only |

These settings are the measurement definition. Values measured with other
settings are different measurements. `--sensitivity` adds the two nearby
definitions: 0.25 mm spacing with a 4 mm window, and 0.5 mm spacing with a
2 mm window. Streamlines too short for the window are counted as exclusions,
never as zeros.

## Running it

```bash
# Convergence: generate a step-size ladder with run_tractography.sh, then assess it
./bin/run_assessment.sh convergence --config hinec_dti \
    --sweep integrator.step=0.4,0.2,0.1,0.05,0.025 --set 'seeding.roi=[Fornix]' --set seeding.density=1

# Consistency: one knob, or whole configs, on the same data
./bin/run_assessment.sh consistency --config hinec_dti --sweep interpolation.method=trilinear,cubic,spline
./bin/run_assessment.sh consistency --configs hinec_dti mmf_dti standard_dti --set 'seeding.roi=[Fornix]' --set seeding.density=1

# Repeatability: independent acquisitions, each preprocessed with run_hinec.sh and tracked
./bin/run_assessment.sh repeatability \
    --group subj1=hinec_runs/run_A_scan1,hinec_runs/run_A_scan2 \
    --group subj2=hinec_runs/run_B_scan1,hinec_runs/run_B_scan2

# Any existing runs or .mat files (a bare .mat needs the NIfTI of its grid)
./bin/run_assessment.sh convergence --param integrator.step --runs <run> <run> ... [--reference <run>]
./bin/run_assessment.sh compare runA runB
./bin/run_assessment.sh consistency --runs a.mat::grid.nii.gz b.mat::grid.nii.gz
```

The output goes to `hinec_runs/assess_<timestamp>_<study>[_<name>]/`:

- `report.html`: a plain-language reading of each tier, one figure, and a table
  of every comparison.
- `summary.json`: every tier of every comparison, the inputs with SHA-256
  hashes, and the measurement protocol.
- `pairs.csv`: one row per matched streamline pair (separation, divergence,
  endpoints, and geometry for both streamlines).

Useful flags: `--name`, `--sensitivity`, `--grid-mm` (force a world grid),
`--max-tracks` / `--max-pairs` (deterministic subsets for whole-brain runs),
`--workers`, `--reference` (an independent convergence reference), and
`--finer larger` (for knobs such as seed density, where a larger value is finer).

**Resources.** Generated trackings run with `HINEC_MAX_WORKERS=4` and the
analysis uses 4 Python processes. Raise them (`HINEC_MAX_WORKERS=16`,
`--workers 16`) only when the machine is free. Before a sweep, check the seed
count in a single run's log: a sweep multiplies it by the number of levels.

## Reading the result

The report uses the words *high*, *moderate* and *low*:

| Kind of number | high | moderate | low |
|---|---|---|---|
| overlap (Jaccard, Dice) | ≥ 0.8 | 0.5 to 0.8 | < 0.5 |
| rank agreement (Spearman) | ≥ 0.8 | 0.5 to 0.8 | < 0.5 |
| median shift | ≤ 2% | 2% to 10% | > 10% |

These bands are **reading conventions, not validated pass/fail thresholds**.
No universal cutoff is supported by the evidence so far.

Rules that come from the study:

1. **Geometry agreement is not agreement.** When geometry agrees but space does
   not, the report prints a warning. The spatial tier must always be reported
   next to the geometry tier.
2. **Distribution agreement is weaker evidence than paired agreement.** A small
   shift in the median can coexist with weak per-streamline Spearman ρ.
3. **Convergence of the typical streamline is not convergence of every
   streamline.** Read the worst-5% separation and the fraction of pairs over
   1 mm, not only the median.
4. **Self-convergence is not accuracy.** If no independent `--reference` is
   given, the finest level is the reference. The report says so.
5. **Repeatability has one unit per subject/session group.** It is not
   measured per streamline. Two subjects are two data points.

## Caveats

- `seeding.strategy: random` calls MATLAB `rand` without seeding the generator.
  Every `matlab -batch` session starts from the same default seed, so repeated
  "random" runs are expected to place the same seeds. Do not use it as a source
  of independent repeats.
- Spatial overlap is counted on the native voxel grid when both runs share it.
  Otherwise (for example across sites) it is counted on an isotropic world-mm
  grid. `summary.json` records which grid was used.
- The pipeline measures agreement. It does not measure anatomical correctness.
  For correctness, use `bin/run_ismrm_scoring.sh` or another source of ground
  truth.

## Validation

- `tests/python/test_tract_assessment.py` (pytest) contains known-answer tests:
  helix curvature and torsion, identical and translated tractograms, reversed
  storage order, nearest pairing with seeds shuffled, a manufactured
  O(h²) convergence ladder, and a MATLAB v7.3 layout with a reflected,
  anisotropic affine.
- Rerunning the tract-geometry study's inputs through this pipeline reproduces
  its separately computed results exactly:
  - four traveling-human pairs: occupied Jaccard 0.552/0.440/0.515/0.459,
    κ<sub>w</sub> shift 2.28/2.15/4.72/5.03%, paired ρ 0.811/0.369/0.563/0.490;
  - MADI MMF step ladders: 2/66 cubic and 7/66 trilinear seeds with p95
    separation above 1 mm at h = 0.0125.
