# MMF Connection-Form Tractography

**Method of Moving Frames (MMF) tractography** — Chun & Peng, in preparation.

This is the *genuine* moving-frames tracker: it builds an orthonormal **moving-frame
field** {e1, e2, e3} and its **connection 1-form** (curvature + torsion) over the whole
brain, then traces streamlines by *evolving a carried frame with the connection structure
equation* rather than by re-sampling a direction field at each step.

!!! note "Not to be confused with the interpolated tracker"
    The interpolated streamline tracker (spatial interpolation, RK4/RKF45, CSD-peak
    resampling) is `algorithm: hinec` — see
    [Tractography Methods](TRACTOGRAPHY_METHODS.md). *That* tracker was previously
    mislabelled `mmf` + `integrator: rkf45`. The tracker documented here — `algorithm: mmf`
    — is the real connection-form Method of Moving Frames and shares none of that code.

| | |
|---|---|
| **Dispatch** | `algorithm: mmf` → `runTractography` → `nim_tractography_mmf_connframe` |
| **Geometry build** | `nim_mmf_geometry`, called by `runTractography` **step 3**, every `mmf` run |
| **Configs** | `config/mmf_dti.yml` (DTI field), `config/mmf_csd.yml` (CSD field) |
| **Reference** | Chun & Peng, in preparation. Equation numbers below follow that formulation; the equations themselves are written out on this page. |

---

## Why moving frames?

Classical streamline tractography integrates `dx/ds = v(x)`, where `v` is a diffusion
direction re-interpolated at every point. It carries no memory of *how* the fiber is
turning. The Method of Moving Frames instead attaches an orthonormal frame to each point of
the trajectory and describes the fiber intrinsically through the **connection 1-form** of
that frame field — the same object that, in the Frenet–Serret picture, packages curvature
and torsion:

$$
\frac{d}{ds}\begin{bmatrix} e_1 \\ e_2 \\ e_3 \end{bmatrix}
=\begin{bmatrix} 0 & \omega_{12} & \omega_{13} \\ -\omega_{12} & 0 & \omega_{23} \\ -\omega_{13} & -\omega_{23} & 0 \end{bmatrix}
\begin{bmatrix} e_1 \\ e_2 \\ e_3 \end{bmatrix}
\qquad(\text{Eq. 10})
$$

$$
\frac{dx}{ds} = e_1 \qquad(\text{Eq. 11})
$$

The connection coefficients are the **curvature** \( \kappa = \nabla_{e_1} e_1 \)
(with \( \omega_{12} = \kappa\cdot e_2 \), \( \omega_{13} = \kappa\cdot e_3 \)) and the
**torsion** \( \tau = \omega_{23}(e_1) \). Those are properties of the *space*: they
depend only on the direction field, not on any particular streamline, so HINEC builds the
whole field **once per run**, before any streamline starts, and the tracer only samples it.
They are **not** part of the `nim` on disk — they depend on `tractography.field`, and the
nim is the dataset. The build costs 4.9 s (DTI) / 8.1 s (CSD), which is why it is redone
every run instead of being cached and version-stamped.

---

## Stage 1 — Building the geometry (`nim_mmf_geometry`)

Called from `runTractography` **step 3**, after the direction field (step 2) and only when
`algorithm: mmf` — `hinec` and `standard` never build it:

```matlab
%% Step 3 - geometry: MMF moving frames + connection 1-form (Eq 6-9), mmf ONLY
if strcmpi(algorithm, 'mmf')
    nim = nim_mmf_geometry(nim, options);   % honours field; no denoising parameter
end
```

With `field: dwi` the frame and curvature come from `nim_mmf_from_dwi`, which step 2
(`nim_field`) runs; `nim_mmf_geometry` only **consumes** `nim.mmf_e1_dwi` /
`nim.mmf_kappa_dwi` and errors if they are absent.

The build follows the Frenet construction of the formulation (steps 1–3, Eq. 6–9):

**1. Alignment-selective tangent field `e1`.**
The raw tangent is the tensor principal eigenvector (DTI) or, when a CSD FOD peak field is
present *and* requested, the **dominant FOD peak** (so `field: csd` genuinely builds the
connection from CSD data — see [CSD field](#csd-field-multiple-pathways) below). It is then
used directly, with **no denoising step**. An alignment-selective filter used to sit here
(`mmf_traj_denoise`, each 3×3×3 neighbour weighted \( |n\cdot e_1|^{\text{sel}} \) with
selectivity `mmf.frame_sel_power`), but that was the same mechanism as the `sel_power` term
removed from `hinec` — a free exponent with no principled value — and it has been removed
with it. Measured against the curvature of the ISMRM ground-truth curves on 5699
`Cingulum_right` voxels, the exponent changed nothing: correlation with the true curvature
was 0.229 / 0.239 / 0.224 / 0.217 at sel 0 / 2 / 16 / 64.

**2. Curvature vector \( \kappa = \nabla_{e_1} e_1 \) (Eq. 7 source).**
Computed basis-free through the connection form itself: with a reference-axis frame
`(e2a, e3a)`, \( \kappa = \omega_{12}(e_1)\,e_2^a + \omega_{13}(e_1)\,e_3^a \).

**3. Frenet normal `e2` (Eq. 7) and binormal `e3` (Eq. 8).**
\( e_2 = \kappa/\lVert\kappa\rVert \) where curvature is well-defined; where
\( \lVert\kappa\rVert \approx 0 \) (straight fiber) it falls back to the **reference-axis
projection** (Eq. 6, `mmf_reference_axis_frame`), which is robust to the
\( \lambda_2 \approx \lambda_3 \) degeneracy that makes the DTI e2/e3 arbitrary. Then
\( e_3 = e_1 \times e_2 \).

**4. Torsion \( \tau = \omega_{23}(e_1) \) (Eq. 9).**
A second pass of `nim_connection_form` on the completed Frenet frame gives the torsion.

**Returned on the `nim` for this run** (in memory only — never saved to disk):

| Field | Shape | Meaning |
|---|---|---|
| `nim.mmf_frames` | `[X Y Z 3 3]` | frame field, `(:,:,:,c,i)` = component *c* of \(e_i\) |
| `nim.mmf_kappa`  | `[X Y Z 3]`   | curvature vector \( \nabla_{e_1} e_1 \) |
| `nim.mmf_tau`    | `[X Y Z]`     | torsion \( \omega_{23}(e_1) \) |

There are no build/version stamps (`mmf_built`, `mmf_field`, `mmf_geom_version` are gone):
the geometry is rebuilt unconditionally on every `mmf` run, so there is nothing to
invalidate. The log line says which field it was actually built from.

### The connection 1-form (`nim_connection_form`)

The connection 1-form \( [\omega] = dA\,A^{T} \) (with \( A=[e_1\ e_2\ e_3] \)) records how
the orthonormal frame rotates through space. It is antisymmetric, so it reduces to three
scalar 1-forms \( (\omega_{12},\omega_{13},\omega_{23}) \). Evaluated along a frame
direction \(e_k\):

$$
\omega_{ij}(e_k) = \sum_c e_j^{\,c}\,\big(\nabla e_i^{\,c}\cdot e_k\big)
\qquad(\text{component form of the connection 1-form})
$$

i.e. take the spatial Jacobian of each component field of \(e_i\), contract with \(e_k\) to
get the directional derivative of \(e_i\) along \(e_k\), then dot with \(e_j\).
`nim_connection_form` returns `wijk[X,Y,Z,i,j,k] = ω_ij(e_k)` with antisymmetry enforced.

### Frame helpers

| Function | Role |
|---|---|
| `nim_build_frames` | sign-consistent `e1` field + reference-axis `e2/e3` |
| `mmf_reference_axis_frame` | complete `e2, e3` from `e1` alone (robust to λ₂≈λ₃) |
| `mmf_gram_schmidt` | re-orthonormalize a drifted frame during integration |
| `mmf_bishop_update` | rotation-minimizing (Bishop) parallel transport of a frame vector |

---

## Stage 2 — Tracing (`nim_tractography_mmf_connframe`)

The geometry arrives already built by step 3 — the tracer never builds it and asserts that
`nim.mmf_frames` is present. It wraps each field in a `griddedInterpolant` and traces:

1. **Seeds** are placed on the seed mask (`seeding.density` sub-voxel offsets, on a
   deterministic lattice). For DTI, each
   seed's initial tangent is the principal direction; for CSD, **one seed per FOD peak**, so
   crossing populations are all launched.
2. **Bidirectional** tracking (`track_bi`) from each seed.
3. Each streamline **carries a moving frame**: `e1 = tangent`, `e2/e3` initialized by the
   reference-axis projection at the seed (Eq. 6), then **evolved by the structure equation
   (Eq. 10)** while advancing `dx/ds = e1` (Eq. 11).
4. The coupled \((x, e_1, e_2, e_3)\) system is integrated by `mmf_step` with **RK4**
   (fixed step) or **RKF45** (adaptive Dormand–Prince), the frame **re-orthonormalized**
   (`mmf_gram_schmidt`) each step. The curvature and torsion are interpolated from the
   stored connection field at each substep.
5. **Termination**: the turn budget (`termination.angle_max`), the propagation mask, which
   is derived from FA (`termination.fa_min`) and the brain mask **independently of the seed
   mask**, and ACT tissue rules — stop on entering CSF or leaving the brain, terminate on
   reaching GM. Tracks shorter than `termination.min_arc` (an arc length in voxels) are
   discarded.

!!! info "The path equation"
    Because \( dT/ds = \kappa N \) in the Frenet picture, the streamline **path** is a
    curvature-vector-field streamline: \( de_1/ds = \kappa(x) \). The full frame *and*
    torsion are evolved and available downstream, but do not feed back into `dx/ds`.
    `mmf_anchor` optionally re-anchors `e1` toward the field tangent (see below).

### `mmf_anchor` — faithful vs. anchored

```yaml
tractography:
  mmf:
    anchor: 0     # 0 = pure Eq.10-11 (faithful); >0 re-anchors e1 to the field
```

- `anchor: 0` — the **pure** connection-form formulation: `e1` is driven only by the
  integrated curvature. Most faithful to the formulation.
- `anchor` in `(0, 1]` — after each step, blend `e1` toward the interpolated field
  tangent: \( e_1 \leftarrow (1-a)\,e_1 + a\,e_1^{\text{field}} \), then re-orthonormalize.
  A stabilizer that trades faithfulness for robustness against accumulated curvature error.

---

## CSD field: multiple pathways

With `field: csd`, `nim_mmf_geometry` builds a **per-peak connection**: for every FOD peak
\(p\) (a distinct fiber population) it computes a curvature vector
\( \nabla_{e_{1p}} e_{1p} \) by **trajectory-aligned differencing** — each neighbour
contributes the peak best aligned with the centre's peak \(p\) (peak matching across the
crossing). Stored as:

| Field | Shape | Meaning |
|---|---|---|
| `nim.mmf_peakdirs` | `[X Y Z maxK 3]` | unit peak directions |
| `nim.mmf_kappa_p`  | `[X Y Z maxK 3]` | per-peak curvature vectors |
| `nim.mmf_multi`    | `true`           | multi-frame flag |

At trace time the tracer selects the peak aligned with the **incoming** tangent, so two
streamlines entering one voxel from different approach directions follow *different*
curvatures → different continuations. That is how a single moving-frame formulation resolves
crossings into multiple pathways.

`runTractography` computes/loads the FOD peaks (via `nim_csd`, cached as `<source>_csd.mat`)
before tracking — see [CSD](TRACTOGRAPHY_METHODS.md#csd-fod-reconstruction).

---

## Configuration

Two sanctioned configs; run them with the [config-driven
workflow](TRACTOGRAPHY_METHODS.md#config-driven-experiment-workflow):

```bash
./bin/run_tractography.sh mmf_dti --score        # DTI-field connection
./bin/run_tractography.sh mmf_csd --score        # CSD-field connection (multiple pathways)
```

### Key parameters

Paths below are canonical config paths under `tractography:`; defaults come from
`src/nim_utils/nim_config_schema.m`.

| Parameter | Default | Meaning |
|---|---|---|
| `algorithm` | `hinec` | must be `mmf` to select this tracker (the dispatch key) |
| `field` | `dti` | `dti` = tensor principal direction; `csd` = per-peak FOD connection |
| `integrator.method` | `rk4` | numerical stepping scheme; `rk4` (fixed) or `rkf45` (adaptive) |
| `integrator.step` | `0.2` | step in voxels (initial step for `rkf45`) |
| `interpolation.method` | `trilinear` | kernel used to sample the stored connection field |
| `mmf.anchor` | `0` | `0` = pure Eq.10-11; `>0` re-anchors `e1` to the field |
| `termination.angle_max` | `225` | turn budget in degrees **per voxel of arc** |
| `termination.fa_min` | `0.10` | FA floor for propagation |
| `termination.min_arc` | `15` | minimum track arc length, in **voxels** |
| `termination.max_arc` | `400` | maximum track arc length in voxels; `max_steps` is derived as `ceil(max_arc/step)` |
| `csd.lmax`, `csd.max_peaks`, `csd.peak_thresh`, `csd.peak_min_sep` | see [CSD](TRACTOGRAPHY_METHODS.md#csd-fod-reconstruction) | FOD peak extraction (only `field: csd`) |

!!! warning "`angle_max` is a rate, and it has a ceiling"
    `termination.angle_max` is degrees of turning per **voxel of arc**, so the budget for
    one step is `angle_max × step` and the criterion is step-invariant. Consecutive tangents
    are sign-aligned — \(e_1\) names a line, not a ray — so a measured turn never exceeds
    90°. Any budget above that is **inert**, not merely loose: the shipped `mmf_dti.yml`
    setting of 225°/voxel goes inert for any step ≥ 0.4. `angle_max: 0` disables the
    criterion outright.

!!! note "`algorithm` is the only dispatch key"
    `runTractography` dispatches purely on `algorithm`. `integrator.method` selects the
    stepping scheme *within* the chosen tracker; changing it never changes which tracker
    runs. Earlier configs carried `integrator: mmf` as a readability marker — that is no
    longer an accepted value.

---

## Tests

- `tests/unit/TestConnectionForm.m` — `nim_connection_form` on analytic frame fields.
- `tests/unit/TestBishopFrame.m` — `mmf_bishop_update` transport / orthonormality.
- `tests/fixtures/make_synthetic_nim.m` — synthetic `nim` (known curvature/torsion) for the above.

## See also

- [Tractography Methods](TRACTOGRAPHY_METHODS.md) — the three trackers and how they dispatch.
- [YAML Config](YAML_CONFIG.md) — full parameter reference.
- [Mathematical Foundations](MATHEMATICAL_FOUNDATIONS.md) — frames, connections, curvature.
