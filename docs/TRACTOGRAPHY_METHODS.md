# Tractography Methods

HINEC ships **three trackers** over one shared `nim` struct, selected by the
`algorithm:` key in the config's `tractography:` section. A second axis, `field:`,
chooses whether direction information comes from the **DTI tensor** or from **CSD FOD
peaks**. This page is the map; each tracker links to its own deep-dive.

```
runTractography(nim, config)
        │  dispatch on config.tractography.algorithm   (runTractography.m:265)
        ├── 'standard' → nim_tractography_standard      FACT, discrete voxel tensors
        ├── 'hinec'    → nim_tractography_hinec          interpolated streamlines (RK4/RKF45)
        └── 'mmf'      → nim_tractography_mmf_connframe   connection-form Method of Moving Frames
```

---

## What it produces

Seeding from a region rather than the whole brain gives bundle-specific
tractography. Below, six ISMRM 2015 bundles spanning the major fibre classes,
each seeded from its own mask and filtered by the endpoint pair and containment
corridor that define it — see [`seeding.roi`](YAML_CONFIG.md) and
[`filter.endpoints_in`](YAML_CONFIG.md).

![Six bundles](img/tractogram_bundles.png)

Each is reconstructed independently; they are drawn together here only to show
their spatial relationship. Per-bundle comparisons against ground truth are in
[ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md#reconstruction-against-ground-truth).

---

## The three trackers

| Algorithm | File | Model | Integration | Interpolation | `field` |
|---|---|---|---|---|---|
| `standard` | `nim_tractography_standard.m` | single tensor | Euler (FACT) | none (discrete voxel tensor) | dti |
| `hinec` | `nim_tractography_hinec.m` | tensor **or** CSD FOD | RK2 / RK4 / RKF45 | trilinear + trajectory-dependent (`sel_power`) | dti \| csd |
| `mmf` | `nim_tractography_mmf_connframe.m` | tensor **or** CSD FOD | RK4 on the frame system | connection 1-form field | dti \| csd |

### `standard` — FACT

Classic **F**iber **A**ssignment by **C**ontinuous **T**racking. Steps voxel-to-voxel
along the discrete principal eigenvector with no sub-voxel interpolation. Fast, simple,
the reference baseline. Full internals: **[Standard Tractography](TRACTOGRAPHY.md)**.

### `hinec` — interpolated streamlines

The workhorse interpolated tracker:

- **Trilinear + trajectory-dependent interpolation.** `sel_power` (default `0` = plain
  trilinear) weights each contributing voxel by how well its direction aligns with the
  incoming trajectory, \( w = |n\cdot v|^{\text{sel\_power}} \). Higher `sel_power`
  sharpens crossings — the streamline follows the population it is already on instead of
  averaging across the crossing. `config/hinec_dti.yml` and `config/hinec_csd.yml` use 16.
- **High-order integration.** `integration_order` = 1 (Euler), 2 (RK2), 4 (RK4), or 5
  (**RKF45** Dormand–Prince adaptive step; set `adaptive_step: true`). See
  [RKF45](RKF.md).
- **`field: dti`** tracks the interpolated principal eigenvector.
  **`field: csd`** tracks **FOD peaks**: at each step it selects the peak best aligned with
  the approach direction, resolving crossings by trajectory. Requires the CSD peak field
  (see [below](#csd-fod-reconstruction)).
- **ACT** (Anatomically Constrained Tractography) when WM/GM/CSF masks are present.

!!! note "Renaming (July 2026)"
    `hinec` is the tracker previously **mislabelled** `algorithm: mmf` with
    `integrator: rkf45`. It has nothing to do with the Method of Moving Frames. The name
    `mmf` now belongs to the genuine connection-form tracker below.

### `mmf` — connection-form Method of Moving Frames

Chun & Peng (in preparation). Builds an orthonormal **moving-frame
field** {e1, e2, e3} plus its **connection 1-form** (curvature + torsion) into the `nim`
(`main.m` Step 2b), then traces by evolving a carried frame with the connection structure
equation — not by re-sampling a direction field. `field: csd` builds a *per-peak*
connection, giving multiple pathways through crossings. Full internals:
**[MMF Connection-Form Tractography](MMF_TRACTOGRAPHY.md)**.

---

## The `field:` axis — DTI vs CSD

Both `hinec` and `mmf` accept `field: dti` (default) or `field: csd`.

- **`dti`** — direction(s) from the diffusion-tensor principal eigenvector. One direction
  per voxel; cannot represent crossings.
- **`csd`** — direction(s) from **Constrained Spherical Deconvolution** FOD peaks. Up to
  `csd_max_peaks` fiber populations per voxel; resolves crossings.

### CSD FOD reconstruction

`src/nim_calculation/nim_csd.m` — a native single-shell, single-tissue CSD (Tournier 2007):

1. real spherical-harmonic fit of the DW signal (design matrix from bvecs, `lmax`),
2. single-fiber **response function** estimated from high-FA voxels,
3. constrained spherical deconvolution (iterative non-negativity),
4. **peak extraction** from the FOD.

It writes peaks into the **same interface** the multi-direction trackers consume, so CSD and
the DTI field are interchangeable upstream of tracking:

| `nim` field | Meaning |
|---|---|
| `nim.peaks`   `[X Y Z maxK 3]` | up to `maxK` unit peak directions per voxel |
| `nim.npeaks`  `[X Y Z]`        | number of valid peaks per voxel |
| `nim.peak_w`  `[X Y Z maxK]`   | peak amplitudes (used by [SIFT](#sift-tractogram-filtering)) |
| `nim.fod_sh`                   | SH FOD coefficients |

**Provisioning & caching.** When `algorithm: mmf` and `field: csd`, `runTractography`
(mmf branch, `runTractography.m:265`) computes the FOD once with `nim_csd` and caches it next
to the source nim as `<source>_csd.mat` (e.g. `data/ismrm2015/ismrm2015_csd.mat`), reusing
the cache on subsequent runs.

!!! warning "Known limitation — `hinec` + `field: csd` peak provisioning"
    The CSD compute/load block currently lives **only** in `runTractography`'s `mmf` branch.
    The `hinec` branch calls `nim_tractography_hinec` directly, which **errors**
    (`HINEC field=csd needs nim.peaks/npeaks (run nim_csd)`) if `nim.peaks` is not already on
    the loaded `nim` — which the standard flow does not populate. So `config/hinec_csd.yml`
    only runs if peaks are provisioned first (e.g. by an earlier `mmf`+`csd` run, or by
    hoisting the CSD block to cover both branches). `mmf_csd.yml` is unaffected.

**CSD parameters** (config `tractography:` keys; defaults as applied by `runTractography`):

| Key | Default | Meaning |
|---|---|---|
| `csd_lmax` | 6 | SH order (capped by #directions) |
| `csd_max_peaks` | 3 | max peaks kept per voxel |
| `csd_peak_thresh` | 0.5 | min peak amplitude (fraction of max FOD) |
| `csd_peak_min_sep` | 45 | min angular separation between peaks (deg) |

---

## SIFT — tractogram filtering

`src/nim_tractography/nim_sift.m` + `bin/run_sift.sh`. A **post-processing** filter
(Smith 2013 analog, per-lobe): it prunes streamlines so that streamline density on each FOD
**lobe** matches that lobe's FOD amplitude. Each streamline segment is attributed to the FOD
peak it best aligns with; lobes carrying more density than their amplitude warrants are
over-represented, and streamlines running through them are removed. This cuts **overreach**
(false-positive spray) for precision, and — unlike a per-voxel density match — distinguishes
a legitimately dense bundle from invalid spray loading a weak/wrong lobe.

Requires a CSD FOD field (`nim.peaks` + `nim.peak_w`). No re-tracking — it filters an
existing run's tractogram into a **new scorable run dir**:

```bash
./bin/run_sift.sh <run_dir> [--score] [--set sift_key=val ...] [--csd <csd.mat>]
```

| Key | Default | Meaning |
|---|---|---|
| `sift_batch_frac` | 0.03 | fraction of streamlines removed per iteration |
| `sift_n_iter` | 60 | max iterations |
| `sift_min_keep` | 0.20 | floor on fraction of streamlines kept |
| `sift_align_pow` | 1 | weight excess by mis-alignment (0 = off) |

Output → `hinec_runs/<run>_sift/`. With `--score`, runs the
[ISMRM scorer](ISMRM_SCORING_ANALYSIS.md) on the filtered result.

---

## PFT — external baseline (DIPY)

`scripts/run_pft_dipy.py` — **Particle Filtering Tractography** (Girard 2014) via DIPY, the
winning ISMRM-2015 approach per Renauld 2023: CSD FOD + probabilistic tracking + CMC
anatomical priors + particle-filter backtracking. It exists as a **comparison baseline** for
the recall/overlap deficit of deterministic tracking, and emits a `.trk` (RAS world) for the
same scilpy scorer. Tissue masks are the DWI-space masks from HINEC preprocessing (no
T1→DWI registration to misalign).

```bash
python scripts/run_pft_dipy.py --out tracks.trk [--density 1] [--seed wm|interface] [--rng K]
```

---

## Track data structure

Every tracker returns the same structure: `tracks` is a cell array; each cell is an `N×3`
matrix of **voxel-space** coordinates giving the *complete* trajectory (not just endpoints);
`N` varies per fiber. Saved under `<run_dir>/tractography/tracks_<algorithm>_<timestamp>.mat`
together with `options`, `elapsed_time`, and `algorithm`.

---

## Config-driven experiment workflow

**Preprocess once, then iterate tractography** without re-preprocessing. This is the
sanctioned way to run experiments — do **not** hand-write MATLAB batch scripts or hand-create
run dirs.

**1. Preprocess once** (`main.m` writes the canonical nim to the *data layer*, plus a copy in
the run dir, and reuses an existing preprocessed DWI if present):

```bash
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml
# → data/ismrm2015/ismrm2015.mat  (canonical processed nim)
```

**2. Iterate tractography** (auto-sourced from the data-layer nim; each run gets its own
tagged, scorable run dir):

```bash
# single config
./bin/run_tractography.sh hinec_csd --score

# design-space sweep with --set (DSE), no new YAML per point
for s in 8 16 32; do
  ./bin/run_tractography.sh hinec_csd --score --set sel_power=$s
done
```

- `<config>` accepts `config/<name>.yml`, `<name>.yml`, or bare `<name>`.
- **Source** defaults to `data/ismrm2015/ismrm2015.mat`; override with
  `--source <nim.mat | run_dir>`. The DWI reference is **copied (frozen)** into the new run's
  `intermediate/` (+ `SOURCE.txt`) so a later reprocess can't corrupt an old run's scoring
  reference.
- **`--set key=value`** overrides any `config.tractography` param on the CLI; the run-dir
  name is tagged with the overrides and `overrides.txt` records them.
- **`--score`** runs `bin/run_ismrm_scoring.sh` on the run dir →
  `scoring/renauld2023/results.json` (headline keys `mean_f1`, `mean_OL`, `mean_OR_gt`,
  `VB`). Compare runs by that JSON.

**Data layer vs run dirs.** The nim + preprocessed refs + CSD cache
(`data/ismrm2015/ismrm2015_csd.mat`) live in `data/ismrm2015/`; `hinec_runs/run_*/` hold
tractography **outputs only**. Never store a nim inside a run dir.

See also [Run Directory System](RUN_DIRECTORY_SYSTEM.md) and
[ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md).

---

## Choosing a tracker

| Goal | Start with |
|---|---|
| Fast baseline / sanity check | `standard` (`config/standard_dti.yml`) |
| Best deterministic single-tensor result | `hinec` + `field: dti` (`config/hinec_dti.yml`; add `--set angle_thresh=60 --set termination_fa=0.05` to push recall) |
| Crossings resolved by FOD | `hinec` + `field: csd` (`config/hinec_csd.yml`) |
| Intrinsic curvature/torsion geometry | `mmf` (`config/mmf_dti.yml` / `mmf_csd.yml`) |
| Cut overreach after tracking | any of the above, then [`run_sift.sh`](#sift-tractogram-filtering) |
| Probabilistic recall reference | PFT baseline (`scripts/run_pft_dipy.py`) |

See the [config catalog](YAML_CONFIG.md#preset-catalog) for the full preset list and every
parameter.
