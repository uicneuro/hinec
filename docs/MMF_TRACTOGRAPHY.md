# MMF Connection-Form Tractography

**Method of Moving Frames (MMF) tractography** — Chun & Peng, `DiscussionWithPeng_Winter2026.pdf`.

This is the *genuine* moving-frames tracker: it builds an orthonormal **moving-frame
field** {e1, e2, e3} and its **connection 1-form** (curvature + torsion) over the whole
brain, then traces streamlines by *evolving a carried frame with the connection structure
equation* rather than by re-sampling a direction field at each step.

!!! note "Not to be confused with the interpolated tracker"
    The interpolated streamline tracker (`sel_power` directional interpolation, RK4/RKF45,
    CSD-peak resampling) is `algorithm: hinec` — see
    [Tractography Methods](TRACTOGRAPHY_METHODS.md). *That* tracker was previously
    mislabelled `mmf` + `integrator: rkf45`. The tracker documented here — `algorithm: mmf`
    — is the real connection-form Method of Moving Frames and shares none of that code.

| | |
|---|---|
| **Dispatch** | `algorithm: mmf` → `runTractography.m:265` → `nim_tractography_mmf_connframe` |
| **Geometry build** | `nim_mmf_geometry` (`main.m` Step 2b) → stored into the `nim` |
| **Configs** | `config/mmf_dti.yml` (DTI field), `config/mmf_csd.yml` (CSD field) |
| **Reference** | `DiscussionWithPeng_Winter2026.pdf` (equation numbers below refer to it) |

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
**torsion** \( \tau = \omega_{23}(e_1) \). Because those are properties of the *space*
(they depend only on the tensor/FOD field, not on any particular streamline), HINEC builds
them **once** and stores them in the `nim`, exactly like FA or the eigenvectors.

---

## Stage 1 — Building the geometry (`nim_mmf_geometry`)

Called from `main.m` **Step 2b**, right after `nim_dt_spd`/`nim_eig`/`nim_fa`:

```matlab
%% Step 2b: MMF moving-frame geometry (frame field + connection 1-form) baked into the nim
nim = nim_mmf_geometry(nim, config.tractography);   % honours frame_sel_power / field
```

The build is faithful to the spec's Frenet construction (steps 1–3 / Eq. 6–9):

**1. Trajectory-dependent tangent field `e1`.**
The raw tangent is the tensor principal eigenvector (DTI) or, when a CSD FOD peak field is
present *and* requested, the **dominant FOD peak** (so `field: csd` genuinely builds the
connection from CSD data — see [CSD field](#csd-field-multiple-pathways) below). It is then
denoised by an **alignment-selective** filter — *not* an isotropic Gaussian
(`mmf_traj_denoise`): each 3×3×3 neighbour is weighted by \( |n\cdot e_1|^{\text{sel}} \),
so aligned neighbours dominate and misaligned (crossing) neighbours get ≈0 weight. This
denoises without blurring across crossings. The selectivity is `frame_sel_power`
(default 16).

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

**Stored into the `nim`:**

| Field | Shape | Meaning |
|---|---|---|
| `nim.mmf_frames` | `[X Y Z 3 3]` | frame field, `(:,:,:,c,i)` = component *c* of \(e_i\) |
| `nim.mmf_kappa`  | `[X Y Z 3]`   | curvature vector \( \nabla_{e_1} e_1 \) |
| `nim.mmf_tau`    | `[X Y Z]`     | torsion \( \omega_{23}(e_1) \) |
| `nim.mmf_field`  | `'dti'`\|`'csd'` | which field the geometry was **actually** built from |
| `nim.mmf_built`  | `true`        | build flag |

The build is wrapped in a `try/catch` in `main.m`; if it fails, the tracker rebuilds it
on demand.

### The connection 1-form (`nim_connection_form`)

The connection 1-form \( [\omega] = dA\,A^{T} \) (with \( A=[e_1\ e_2\ e_3] \)) records how
the orthonormal frame rotates through space. It is antisymmetric, so it reduces to three
scalar 1-forms \( (\omega_{12},\omega_{13},\omega_{23}) \). Evaluated along a frame
direction \(e_k\):

$$
\omega_{ij}(e_k) = \sum_c e_j^{\,c}\,\big(\nabla e_i^{\,c}\cdot e_k\big)
\qquad(\text{Book B-3.96})
$$

i.e. take the spatial Jacobian of each component field of \(e_i\), contract with \(e_k\) to
get the directional derivative of \(e_i\) along \(e_k\), then dot with \(e_j\).
`nim_connection_form` returns `wijk[X,Y,Z,i,j,k] = ω_ij(e_k)` with antisymmetry enforced.

### Frame helpers

| Function | Role |
|---|---|
| `nim_build_frames` | sign-consistent `e1` field + reference-axis `e2/e3`; `frame_smooth_sigma` |
| `mmf_reference_axis_frame` | complete `e2, e3` from `e1` alone (robust to λ₂≈λ₃) |
| `mmf_gram_schmidt` | re-orthonormalize a drifted frame during integration |
| `mmf_bishop_update` | rotation-minimizing (Bishop) parallel transport of a frame vector |

---

## Stage 2 — Tracing (`nim_tractography_mmf_connframe`)

The tracer reads the precomputed geometry from the `nim` (rebuilding it if absent, or if it
was baked for a different `field` — e.g. DTI geometry baked by `main.m` but this is a CSD
run). It wraps each stored field in a `griddedInterpolant` and traces:

1. **Seeds** are placed on the seed mask (`seed_density` sub-voxel offsets). For DTI, each
   seed's initial tangent is the principal direction; for CSD, **one seed per FOD peak**, so
   crossing populations are all launched.
2. **Bidirectional** tracking (`track_bi`) from each seed.
3. Each streamline **carries a moving frame**: `e1 = tangent`, `e2/e3` initialized by the
   reference-axis projection at the seed (Eq. 6), then **evolved by the structure equation
   (Eq. 10)** while advancing `dx/ds = e1` (Eq. 11).
4. The coupled \((x, e_1, e_2, e_3)\) system is integrated with **RK4** (`mmf_step`), the
   frame **re-orthonormalized** (`mmf_gram_schmidt`) each step. The curvature and torsion
   are interpolated from the stored connection field at each substep.
5. **Termination**: angle threshold (`angle_thresh`), propagation-mask / FA floor
   (`termination_fa`), and ACT tissue rules — stop on entering CSF/OUTSIDE, terminate on
   reaching GM.

!!! info "The path equation"
    Because \( dT/ds = \kappa N \) in the Frenet picture, the streamline **path** is a
    curvature-vector-field streamline: \( de_1/ds = \kappa(x) \). The full frame *and*
    torsion are evolved and available downstream, but do not feed back into `dx/ds`.
    `mmf_anchor` optionally re-anchors `e1` toward the field tangent (see below).

### `mmf_anchor` — faithful vs. anchored

```yaml
mmf_anchor: 0     # 0 = pure Eq.10-11 (faithful); >0 re-anchors e1 to the field
```

- `mmf_anchor: 0` — the **pure** connection-form formulation: `e1` is driven only by the
  integrated curvature. Most faithful to the spec.
- `mmf_anchor` in `(0, 1]` — after each RK4 step, blend `e1` toward the interpolated field
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
./bin/run_tractography.sh mmf_csd --score    # CSD-field connection (multiple pathways)
```

### Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `algorithm` | — | must be `mmf` (the dispatch key) |
| `field` | `dti` | `dti` = tensor principal direction; `csd` = per-peak FOD connection |
| `integrator` | — | descriptive marker (`mmf`) in the sanctioned configs; **dispatch is by `algorithm`**, not this key |
| `mmf_anchor` | `0` | `0` = pure Eq.10-11; `>0` re-anchors `e1` to the field |
| `frame_sel_power` | `16` | alignment selectivity of the tangent denoise when building the frames |
| `step_size` | `0.2` | RK4 step (voxels) |
| `angle_thresh` | `45` | max turn per step (deg) |
| `termination_fa` | `0.05` | FA floor for propagation |
| `min_length` | `15` | min streamline length (mm) |
| `curv_beta` | `2.0` | curvature-gate strength (shared crossing-aware knob) |
| `csd_lmax`, `csd_max_peaks`, `csd_peak_thresh`, `csd_peak_min_sep` | see [CSD](TRACTOGRAPHY_METHODS.md#csd-fod-reconstruction) | FOD peak extraction (only `field: csd`) |

!!! warning "`integrator` is not the dispatch key"
    `runTractography` dispatches purely on `algorithm` (`mmf` → connection-form tracer). The
    sanctioned configs *also* set `integrator: mmf` as a readability marker, and the tracer
    stamps `'mmf'` into its returned `info` struct — but changing `integrator` alone does
    **not** change which tracker runs. Set `algorithm: mmf`.

---

## Tests

- `tests/unit/TestConnectionForm.m` — `nim_connection_form` on analytic frame fields.
- `tests/unit/TestBishopFrame.m` — `mmf_bishop_update` transport / orthonormality.
- `tests/fixtures/make_synthetic_nim.m` — synthetic `nim` (known curvature/torsion) for the above.

## See also

- [Tractography Methods](TRACTOGRAPHY_METHODS.md) — the three trackers and how they dispatch.
- [YAML Config](YAML_CONFIG.md) — full parameter reference.
- [Mathematical Foundations](MATHEMATICAL_FOUNDATIONS.md) — frames, connections, curvature.
