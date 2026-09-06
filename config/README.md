# HINEC Configuration Files

YAML configs for the HINEC pipeline. Every file has two top-level sections —
`preprocessing:` and `tractography:` — parsed by `src/nim_utils/load_config_yaml.m`.

**The parameter reference is [docs/YAML_CONFIG.md](../docs/YAML_CONFIG.md)**, which is
*generated* from `src/nim_utils/nim_config_schema.m` and is the only place parameters
are documented. This file covers naming conventions only; do not duplicate parameter
tables here, or they will drift.

Configs use **two levels of nesting** below a section (`section` → `group` → `key`);
a third level is a parse error. Inline lists (`roi: [41, 42]`) are supported. Every
key is optional, so a working config can be three lines:

```yaml
tractography:
  algorithm: hinec
  seeding:
    roi: [SLF_L, SLF_R]
```

Any parameter can be overridden on the command line, including `preprocessing.*`:

```bash
./bin/run_tractography.sh hinec_dti --set integrator.step=0.05
./bin/run_hinec.sh data/x x.mat config/y.yml --set preprocessing.run_eddy=false
```

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
3. **`integrator.method:`** — *the numerical stepping scheme (how, and how far,
   to advance)*. This is **not** the direction. `euler` | `rk2` | `rk4` |
   `rkf45`, spelled the same way for every tracker. The old
   `integration_order: 1|2|4|5` is a legacy alias and migrates automatically —
   it spelled a method selector as a number. The value `5` selected RKF45, and
   that was not wrong numerically - the implementation uses Dormand-Prince
   coefficients and advances on the 5th-order solution, keeping the embedded
   4th-order one for error control - but a method name belongs in a method key.

> Direction (`field`, `interpolation.method`) and integrator are genuinely
> independent: you can hold one fixed and sweep the other. Settings that don't apply
> to the chosen combination (e.g. `csd.lmax` when `field: dti`) are simply ignored
> — that is expected, not a bug.

### Direction shaping: the interpolation kernel

For `hinec`, the direction at a sub-voxel point is interpolated from the voxel
grid. Because `v1` is a **line** field — the sign the eigensolver writes for each
voxel is arbitrary — the tracker interpolates the dyadic `v1*v1'`, which is
invariant under `v -> -v`, and recovers the principal eigenvector from it.
Interpolating the signed components directly would average across disagreeing
signs and collapse the result toward zero.

`interpolation.method` selects the spatial kernel, and the three differ in
**smoothness**, which is what caps the reachable order of the integrator:

| | class | what it is |
|---|:--:|---|
| `trilinear` | C⁰ | straight lines between samples; kinked at every voxel face |
| `cubic` | C¹ | Keys cubic *convolution*, a fixed 4×4×4 stencil — **not** a spline |
| `spline` | C² | a genuine cubic spline; one global solve, then local evaluation |

A Runge–Kutta method needs enough continuous derivatives in the right-hand side
to express its formal order. Measured on ISMRM 2015, RK4 reaches order **2.00**
on trilinear, **3.06** on cubic and **4.00** on spline — one order per
derivative. See [Convergence & Verification](../docs/CONVERGENCE.md).

`interpolation.upsample` is the **space axis**: the field is sampled on a grid of
spacing `1/upsample` voxels before the interpolants are built. Above 1 refines,
below 1 coarsens; the coordinate frame is unchanged, so runs at different factors
stay directly comparable.

> **`sel_power` has been removed.** It reweighted neighbours by alignment with
> the incoming tangent, which made the ODE's right-hand side depend on the
> trajectory that reached it — there is no justification for that with DTI, where
> each voxel holds one eigenvector and there is nothing to disambiguate. HINEC is
> now interpolation + integration only. For CSD, nearest-peak selection is kept
> because it is structural, but the alignment exponent is gone.

## The configs

### Tracker configs — `<algorithm>_<field>[_<variant>].yml`  (for `run_tractography.sh`)

|            | `field: dti`      | `field: csd`      |
|------------|-------------------|-------------------|
| **standard** (FACT) | `standard_dti.yml` | — *(FACT is DTI-only)* |
| **hinec** (interpolated) | `hinec_dti.yml` (+ variants below) | `hinec_csd.yml` |
| **mmf** (connection-form) | `mmf_dti.yml` | `mmf_csd.yml` |

- `standard_dti` — FACT, discrete voxel tensors, no interpolation. Deterministic baseline.
- `hinec_dti` / `hinec_csd` — interpolated streamlines, trajectory-dependent
  direction, RKF45 + trilinear by default, optional ACT. `hinec_csd` computes
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
| `hinec_dti_euler` | Euler (`integrator.method: euler`) + trilinear — didactic / FACT-vs-high-order comparison |
| `hinec_dti_fast` | RK2 (`integrator.method: rk2`), coarse steps & seeding — quick parameter screening |

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
for h in 0.5 0.25 0.125; do ./bin/run_tractography.sh hinec_csd --score --set integrator.step=$h; done

# swap the interpolation kernel on any hinec config
./bin/run_tractography.sh hinec_dti --score --set interp_method=cubic
./bin/run_tractography.sh hinec_csd --score --set interpolation.method=spline
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
  algorithm: hinec           # hinec | standard | mmf     (which tracker)
  field: csd                 # dti | csd                  (direction source)
  act: false                 # WM/GM/CSF-constrained tracking (hinec; needs tissue masks)

  # --- integrator: ONE scheme for every tracker (mmf no longer differs) ---
  integrator:
    method: rkf45            # euler | rk2 | rk4 | rkf45
    step: 0.2                # voxels; fixed step, or initial step for rkf45
    tolerance: 0.02          # rkf45 ONLY - setting it with a fixed-step method is an error
    step_min: 0.02
    step_max: 0.5

  # --- direction shaping ---
  interpolation:
    method: spline           # trilinear (C0) | cubic (C1) | spline (C2)
    upsample: 1              # spatial sampling factor; >1 refines, <1 coarsens

  # --- seeding ---
  seeding:
    density: 8               # seeds/voxel
    fa_min: 0.05
    roi: []                  # JHU indices and/or names, e.g. [41, 42] or [SLF_L]

  # --- termination: arc length, so refining the step cannot truncate tracks ---
  termination:
    fa_min: 0.10
    angle_max: 225           # deg per VOXEL OF ARC (not per step), so the budget
                             # for one step is angle_max * step. Tangents are
                             # sign-aligned, so a measured turn never exceeds 90
                             # deg and any budget above that is INERT. 0 disables.
    max_arc: 400             # voxels; max_steps is DERIVED as ceil(max_arc/step)
    min_arc: 15              # voxels

  # --- field: csd only ---
  csd:
    lmax: 6
    max_peaks: 3
    peak_thresh: 0.5
    peak_min_sep: 45

  # --- algorithm: mmf only ---
  mmf:
    frame_sel_power: 16
    anchor: 0
```

Every key, its type, range and default is generated from the schema
(`src/nim_utils/nim_config_schema.m`) into
[docs/YAML_CONFIG.md](../docs/YAML_CONFIG.md) — that file is the reference, and
this one is the orientation. Configs are validated on load; unknown keys and
retired keys are reported rather than silently accepted.
