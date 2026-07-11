# HINEC Configuration Files

YAML configs for the HINEC pipeline. Every file has two top-level sections —
`preprocessing:` and `tractography:` — parsed by `src/nim_utils/load_config_yaml.m`
(flat key/value only; the parser has no lists or nesting). Full parameter reference:
[docs/YAML_CONFIG.md](../docs/YAML_CONFIG.md).

## The naming convention (one rule per family)

There are exactly **two kinds** of config, and each follows one naming rule:

| Family | Drives | Name rule | The `preprocessing:` block | The `tractography:` block |
|---|---|---|---|---|
| **Dataset** | `run_hinec.sh` (build the nim, once) | `<dataset>[_variant].yml` | **meaningful** (denoise, eddy, registration…) | a sane default |
| **Tracker** | `run_tractography.sh` (iterate on a built nim) | `<algorithm>_<field>.yml` | minimal (just `denoise_method`) | **meaningful** |

Plus one special file, **`hinec_default.yml`** — the fallback used when no config is
given. It is the only config that pairs a full preprocessing block with a general-purpose
hinec DTI tractography block (cubic interpolation + RK4 + ACT).

## The three orthogonal axes

A tracker config's identity is three independent knobs. Its **filename encodes the
first two**; the third is a knob with a per-algorithm default.

1. **`algorithm:`** — *which tracker*. `standard` (FACT) | `hinec` (interpolated
   streamlines) | `mmf` (connection-form Method of Moving Frames). This is the dispatch key.
2. **`field:`** — *where the local direction comes from*. `dti` (tensor principal
   eigenvector) | `csd` (CSD FOD peaks; needs the `csd_*` params, cached as
   `<source>_csd.mat`). `hinec` and `mmf` accept either; `standard` is DTI-only.
3. **integrator** — *the numerical stepping scheme (how, and how far, to advance)*.
   This is **not** the direction. `hinec`/`standard` spell it `integration_order:`
   (`1`=Euler, `2`=RK2, `4`=RK4, `5`=RKF45 adaptive); `mmf` spells it
   `integrator:` (`rk4` fixed | `rkf45` adaptive). Same concept, historically
   different key names.

> Direction (`field`, `sel_power`, `interp_method`) and integrator are genuinely
> independent: you can hold one fixed and sweep the other. Settings that don't apply
> to the chosen combination (e.g. `csd_lmax` when `field: dti`, `sel_power` when
> `algorithm: mmf`) are simply ignored — that is expected, not a bug.

### Direction shaping: `interp_method` and `sel_power` compose (not either/or)

For `hinec`, the local direction is a weighted blend over neighbor voxels where each
neighbor's weight is a **product of two independent factors**:

```
weight = spatial_kernel(distance)  ×  alignment^sel_power
         └── interp_method ──┘         └──── sel_power ────┘
```

- **`interp_method`** is the *spatial interpolation kernel* — `trilinear` (2×2×2) or
  `cubic` (4×4×4 tricubic). It decides each neighbor's contribution by distance.
- **`sel_power`** is *directional steering* — it reweights neighbors by how well they
  align with the incoming tangent (`0` = no steering = plain interpolation; higher =
  sharper crossing resolution). It is **not** a spatial scheme.

So `sel_power: 16` **and** `interp_method: cubic` combine — tricubic spatial weights, each
sharpened by alignment. `interp_method: cubic` drives the direction at any `sel_power`
(and at `sel_power: 0` it uses the fast cubic `griddedInterpolant`). This is why cubic and
sel_power are two knobs, not a choice between them.

## The configs

### Tracker configs — `<algorithm>_<field>[_<variant>].yml`  (for `run_tractography.sh`)

|            | `field: dti`      | `field: csd`      |
|------------|-------------------|-------------------|
| **standard** (FACT) | `standard_dti.yml` | — *(FACT is DTI-only)* |
| **hinec** (interpolated) | `hinec_dti.yml` (+ variants below) | `hinec_csd.yml` |
| **mmf** (connection-form) | `mmf_dti.yml` | `mmf_csd.yml` |

- `standard_dti` — FACT, discrete voxel tensors, no interpolation. Deterministic baseline.
- `hinec_dti` / `hinec_csd` — interpolated streamlines, trajectory-dependent
  (`sel_power`) direction, RKF45 + trilinear by default, optional ACT. `hinec_csd` computes
  its own FOD peaks (`nim_csd`) at track time — no separate step needed.
- `mmf_dti` / `mmf_csd` — genuine connection-form tracer (Chun & Peng): the moving-frame
  field + connection 1-form are built into the nim and evolved by the structure equation;
  stepping via `integrator: rk4|rkf45`. `mmf_csd` builds a per-peak connection → multiple
  pathways through crossings. See [docs/MMF_TRACTOGRAPHY.md](../docs/MMF_TRACTOGRAPHY.md).

**hinec DTI variants** — same tracker, distinct *reusable* knob-sets (the `_<variant>` suffix):

| Config | What makes it distinct |
|---|---|
| `hinec_dti_cubic` | `interp_method: cubic` + RKF45 — the cubic high-precision config |
| `hinec_dti_cubic_recall` | cubic + aggressive termination (`termination_fa 0.05`, `angle 60`) to push coverage/recall |
| `hinec_dti_euler` | Euler (`integration_order 1`) + trilinear — didactic / FACT-vs-high-order comparison |
| `hinec_dti_fast` | RK2 (`integration_order 2`), coarse steps & seeding — quick parameter screening |

### Dataset configs — `<dataset>[_variant].yml`  (for `run_hinec.sh`)

| Config | What it builds |
|---|---|
| `ismrm2015.yml` | ISMRM-2015 phantom, no T1 registration |
| `ismrm2015_t1reg.yml` | ISMRM-2015 phantom **with** T1→MNI registration |
| `irontract.yml` | IronTract Challenge dataset |

### Default

| Config | Role |
|---|---|
| `hinec_default.yml` | Generic full-pipeline fallback: full preprocessing + a general hinec DTI tractography block (cubic + RK4 + ACT). Used when `run_hinec.sh` / `load_config_yaml` get no config. |

## How to run

```bash
# 1. Preprocess ONCE (build the nim) with a dataset config
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml

# 2. Iterate tractography on the built nim (NO re-preprocess), and score
./bin/run_tractography.sh hinec_csd --score
./bin/run_tractography.sh mmf_dti   --score
```

`run_tractography.sh` accepts `config/<name>.yml`, `<name>.yml`, or bare `<name>`, and
defaults its source to the canonical nim in `data/ismrm2015/ismrm2015.mat`.

## Named variant configs vs `--set` (both are fine — this is a research tool)

Two ways to change knobs; use whichever fits. Neither is "more correct":

- **Named variant config** (`<algorithm>_<field>_<variant>.yml`) — for a setup you'll
  *reuse* and want in version control, like the `hinec_dti_*` variants above. Keep it.
- **`--set key=value`** — for a *one-off* sweep, no new file. The run dir is tagged with
  the override so the experiment stays identifiable.

```bash
# one-off sweep — no file needed
for s in 8 16 32; do ./bin/run_tractography.sh hinec_csd --score --set sel_power=$s; done

# cubic on any hinec config, at any sel_power (composes — see "Direction shaping" above)
./bin/run_tractography.sh hinec_dti --score --set interp_method=cubic
./bin/run_tractography.sh hinec_csd --score --set interp_method=cubic --set sel_power=16
```

Rule of thumb: **reuse it → make a `_<variant>` config; one-shot → `--set`.** Nothing is
off-limits — this is a research tool, so keep whatever experimental configs are useful.

## Creating a new config

Copy the nearest one and keep the naming rule (`<algorithm>_<field>[_<variant>]` for
trackers, `<dataset>[_variant]` for datasets):

```bash
cp config/hinec_csd.yml config/hinec_csd_myvariant.yml
# edit; keep 2-space indent, flat key/value
./bin/run_tractography.sh hinec_csd_myvariant --score
```

## Key parameters (quick reference)

```yaml
tractography:
  algorithm: 'hinec'         # standard | hinec | mmf          (axis 1: which tracker)
  field: 'csd'               # dti | csd                       (axis 2: direction source)
  # --- axis 3: integrator (numerical stepping scheme) ---
  integration_order: 5       # hinec/standard: 1 Euler | 2 RK2 | 4 RK4 | 5 RKF45
  # integrator: 'rk4'        # mmf only:       rk4 (fixed) | rkf45 (adaptive)
  adaptive_step: true        # RKF45 only
  rkf_tolerance: 0.02        # RKF45 error tol (hinec); mmf uses rkf_tol
  # --- direction shaping ---
  sel_power: 16              # hinec: 0 = plain interp; >0 = trajectory-dependent selection
  interp_method: 'cubic'     # hinec: trilinear | cubic | none
  frame_sel_power: 16        # mmf: selectivity for building the moving frames
  mmf_anchor: 0              # mmf: 0 = pure structure equation; >0 re-anchors e1 to field
  # --- seeding / termination ---
  seed_density: 8            # seeds/voxel (1-8)
  step_size: 0.2
  angle_thresh: 45           # max turn per step (deg)
  termination_fa: 0.10       # FA floor
  min_length: 15             # mm
  act_enabled: false         # WM/GM/CSF-constrained tracking (needs tissue masks)
  # --- field: csd only ---
  csd_lmax: 6
  csd_max_peaks: 3
  csd_peak_thresh: 0.5
  csd_peak_min_sep: 45
```

Configs are validated on load (`step_size > 0`, `angle_thresh ∈ (0,180]`,
`integration_order ∈ {1,2,4,5}`, RKF/MMF ranges). See
[docs/YAML_CONFIG.md](../docs/YAML_CONFIG.md).
