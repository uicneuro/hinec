# HINEC Pipeline Overview

HINEC (HIgh-order NEural Connectivity) is a MATLAB pipeline that takes raw NIfTI diffusion
MRI through preprocessing, diffusion tensor estimation, fractional anisotropy, parcellation
and tractography. This document describes each stage and the mathematics behind it.

## Main Entry Points

Two MATLAB functions in the repository root do the work, and shell launchers in `bin/` drive
them:

| Entry point | Role |
|---|---|
| `main.m` | Preprocessing → tensor fit → eigendecomposition → FA → MMF geometry → optional registration → parcellation → mask improvement → tissue segmentation → save |
| `runTractography.m` | Loads the saved `nim`, resolves the seeding strategy, dispatches to a tracker, filters and writes the tracks |
| `bin/run_hinec.sh` | Full pipeline: builds a run directory, then chains `main` and `runTractography` in one `matlab -batch` process |
| `bin/run_tractography.sh` | Tractography only, against an already-preprocessed `nim` |

`runhinec.m` is a short scratch script that loads a processed `.mat` and plots it; it is not
part of the pipeline. Visualization entry points live in `src/nim_visualization/`.

## Pipeline Stages

The pipeline consists of four main stages.

### 1. Preprocessing Pipeline

The preprocessing pipeline (`nim_preprocessing.m`) runs ten sequential steps, using T1
structural data where it is available to improve brain extraction and atlas registration.

#### **Step 1: B0 Extraction**
Extract the first volume ($b=0$) as reference:

$$
B_0(x,y,z) = \text{DWI}(x,y,z,0)
$$

#### **Step 2: Brain Extraction**
Brain mask creation, using T1 structural data when available:

**T1-Based Method (Preferred):**

$$
M_{\text{T1}}(x,y,z) = \text{BET}(T1(x,y,z),\; f=0.4)
$$

$$
M_0(x,y,z) = \text{Transform}(M_{\text{T1}},\; T1 \to \text{DWI})
$$

**DWI-Only Fallback:**

$$
M_0(x,y,z) = \text{BET}(B_0(x,y,z),\; f=0.3)
$$

**Registration Chain:**
```
T_T1→DWI  = epi_reg(B0, T1, T1_brain)
M_DWI     = apply_transform(M_T1, T_T1→DWI)
```

#### **Step 3: Denoising (Optional)**
MP-PCA denoising or Gaussian smoothing:

$$
S_{\text{denoised}}(x,y,z,b) = S_{\text{raw}}(x,y,z,b) \otimes G(\sigma)
$$

where $G(\sigma)$ is a Gaussian kernel with standard deviation $\sigma$.

#### **Step 4: Field Map Distortion Correction**
Susceptibility distortion correction using field maps:

**Field Map Processing:**

$$
\Delta B_0(x,y,z) = \text{fieldmap}_{\text{Hz}}(x,y,z) \quad [\text{Hz}]
$$

**FUGUE Distortion Correction:**

$$
S_{\text{corrected}}(x',y',z',b) = S_{\text{raw}}(x,y,z,b)
$$

where the spatial transformation is:

$$
x' = x + \Delta B_0(x,y,z) \times t_{\text{dwell}} \times \hat{\mathbf{n}}_{\text{PE}}
$$

**Parameters:**

- $t_{\text{dwell}}$: Effective echo spacing (typically 0.58 ms)
- $\hat{\mathbf{n}}_{\text{PE}}$: Phase encoding direction vector
- $\Delta B_0$: Field inhomogeneity in Hz

#### **Step 5: Motion Correction**
Rigid body motion correction with b-vector rotation:

```
R_b = mcflirt(DWI_volumes)
```

The b-vectors must be rotated to match the corrected orientations:

$$
\mathbf{g}_i' = \mathbf{R}_b \, \mathbf{g}_i \quad \forall \, i \in [1, N_{\text{directions}}]
$$

#### **Step 6: Eddy Current Correction**
Eddy current correction, with a fallback when acquisition parameter files are absent:

**Method 1 (Preferred):** FSL eddy with acquisition parameters
**Method 2 (Fallback):** FSL eddy_correct for datasets without acqp/index files

Mathematical model for eddy currents:

$$
S_{\text{corrected}}(i) = S_{\text{raw}}(i) \circ T_{\text{eddy}}(\mathbf{g}_i)
$$

where $T_{\text{eddy}}$ represents the eddy-induced geometric distortion.

#### **Step 7: Copy Processed Data to Final Location**
Write the corrected DWI, the brain mask (`{name}_M.nii.gz`) and the motion/eddy-corrected
b-vectors to their final paths, and record them in the preprocessing report.

Tissue segmentation for ACT is **not** part of this pipeline; it runs later, from `main.m`
(see [Tissue Segmentation](#tissue-segmentation-for-act) below).

#### **Step 8: T1 Preprocessing and Registration (When Available)**
Complete T1-based registration workflow for enhanced atlas processing:

**T1-MNI Registration Chain:**
```
T_linear    = FLIRT(T1_brain → MNI152_1mm)
W_nonlinear = FNIRT(T1 → MNI152, init=T_linear)
W_inverse   = invwarp(W_nonlinear)              [MNI → T1]
```

**DWI Reference Creation:**
```
DWI_ref        = fslroi(DWI_processed, 0, 1)    [Extract first volume]
DWI_ref_masked = DWI_ref × M0                   [Apply brain mask]
```

#### **Step 9: Atlas Registration**
T1-guided atlas transformation through a composite registration chain:

**T1-Based Atlas Registration (Preferred):**
```
Atlas_DWI = applywarp(Atlas_MNI,
                     warp=W_inverse,
                     postmat=T_T1_to_DWI,
                     ref=DWI_ref,
                     interp=nn)
```

**Mathematical Transformation:**

$$
\text{Atlas}_{\text{DWI}}(\mathbf{r}) = \text{Atlas}_{\text{MNI}}\!\left(\mathbf{W}_{\text{inverse}}\!\left(\mathbf{T}_{T1 \to \text{DWI}}^{-1}(\mathbf{r})\right)\right)
$$

**Direct Registration Fallback:**

$$
\text{Atlas}_{\text{DWI}} = \text{FLIRT}(\text{Atlas}_{\text{MNI}} \to \text{DWI}_{\text{ref}},\; \text{interp} = \text{nn})
$$

#### **Step 10: Finalization**
- Copy processed data to standard locations
- Quality validation and reporting
- Cleanup and report generation

### 2. Core Data Processing

DTI processing with robust tensor estimation:

#### **Diffusion Tensor Calculation (`nim_dt_spd`)**
Symmetric positive definite (SPD) constrained tensor estimation:

$$
\mathbf{D} = \begin{bmatrix} D_{xx} & D_{xy} & D_{xz} \\ D_{xy} & D_{yy} & D_{yz} \\ D_{xz} & D_{yz} & D_{zz} \end{bmatrix} \in \mathcal{S}_+^3
$$

**Log-linear fitting:**

$$
\ln(S_i / S_0) = -b_i \, \mathbf{g}_i^\top \mathbf{D} \, \mathbf{g}_i
$$

**SPD Constraint:** Ensure $\mathbf{D}$ has positive eigenvalues: $\lambda_1 \geq \lambda_2 \geq \lambda_3 > 0$

#### **Fractional Anisotropy (`nim_fa`)**

$$
\text{FA} = \sqrt{\frac{3}{2}} \cdot \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}
$$

where $\bar{\lambda} = (\lambda_1 + \lambda_2 + \lambda_3)/3$ is the mean diffusivity.

**Tensor Eigendecomposition:**

$$
\mathbf{D} = \mathbf{V} \, \mathbf{\Lambda} \, \mathbf{V}^\top
$$

where:

- $\mathbf{V} = [\mathbf{v}_1 \; \mathbf{v}_2 \; \mathbf{v}_3]$ are eigenvectors (fiber directions)
- $\mathbf{\Lambda} = \mathrm{diag}(\lambda_1, \lambda_2, \lambda_3)$ are eigenvalues

#### **Moving-Frame Geometry (`nim_mmf_geometry`)**

Immediately after the eigendecomposition, `main.m` builds the moving-frame field and its
connection 1-form into the `nim`. This is a property of the direction field, not of any
individual streamline, so it is computed once here rather than per track: $\mathbf{e}_1$ is
the denoised tangent, $\mathbf{e}_2 = d\mathbf{e}_1/ds$ the Frenet normal,
$\mathbf{e}_3 = \mathbf{e}_1 \times \mathbf{e}_2$, and the connection
$[\omega] = d[A]\,A^\top$ carries curvature and torsion. The `mmf` tracker traces through this
field; the other trackers ignore it. Failure here is non-fatal — the tracker rebuilds the
geometry on demand.

For the derivation see Chun & Peng, in preparation, and
[MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md).

#### **Tissue Segmentation for ACT**

`main.m` generates the WM/GM/CSF masks that Anatomically Constrained Tractography needs, and
stores them on the `nim` as `.wm_mask`, `.gm_mask` and `.csf_mask`.

The primary method is FSL FAST on the anatomical T1, resampled into DWI space through the two
images' **world affines** (`flirt -usesqform`). Going through the world affines rather than
the T1→DWI registration means the masks do not depend on — and cannot be corrupted by — a
registration step that may be disabled or poorly converged. This is correct whenever the T1 is
already world-aligned to the diffusion data.

The fallback, used only when no T1 is present, bins FA into tertiles. It bins **anisotropy,
not tissue**, so ACT driven by it will terminate streamlines mid-crossing. Treat it as a
degraded mode rather than an equivalent one.

Note that these masks are always generated; ACT is nonetheless **off by default**
(`tractography.act` defaults to false) and must be enabled in the config.

### 3. Tractography

`runTractography.m` makes every seeding decision, then dispatches to one of three trackers
selected by `tractography.algorithm`. The trackers themselves do no seeding.

#### **Seeding**

Which voxels are seeded is decided by the first strategy that applies:

1. **Explicit ROI** — `seeding.roi` names atlas regions; seeds go only there. This is the
   basis of single-bundle experiments and of the convergence ladder's fixed seed set.
2. **Preprocessed brain mask** (`nim.mask`) — the usual path.
3. **Expanded parcellation mask** — the parcellation dilated by 3 voxels.
4. **FA threshold** (FA > 0.10) — a fallback when no mask exists at all. It misses
   low-anisotropy structures such as the fornix and cingulum.

Every strategy then intersects with `seeding.fa_min` (default 0.05, which excludes CSF
only). Within each seeded voxel, `seeding.density` seeds are placed on a deterministic
sub-voxel lattice (`seeding.strategy: uniform`) or jittered (`random`). The count is honoured
exactly, so seeding is reproducible run to run.

#### **Trackers**

| `algorithm` | Function | Direction at a point |
|---|---|---|
| `standard` | `nim_tractography_standard` | FACT: the discrete voxel's principal eigenvector, no interpolation |
| `hinec` | `nim_tractography_hinec` | Interpolated direction field, integrated with a Runge-Kutta method |
| `mmf` | `nim_tractography_mmf_connframe` | The carried frame, evolved by the connection 1-form |

**FACT integration** advances along the voxel's own direction until the streamline crosses a
voxel face:

$$
\mathbf{r}_{i+1} = \mathbf{r}_i + \Delta s \cdot \mathbf{e}_1(\mathbf{r}_i)
$$

**HINEC integration** replaces $\mathbf{e}_1(\mathbf{r}_i)$ with an interpolated field
$\mathbf{v}(\mathbf{r})$ and integrates $d\mathbf{x}/ds = \mathbf{v}(\mathbf{x})$ with
`integrator.method` ∈ {`euler`, `rk2`, `rk4`, `rkf45`}. Two things shape the right-hand side:

- `interpolation.method` sets the kernel's smoothness — `trilinear` is C0 (kinked at every
  voxel face), `cubic` is C1 (Keys cubic convolution, whose second derivative jumps),
  `spline` is C2. This caps the attainable order: a Runge-Kutta method of order $p$ needs a
  right-hand side with $p$ continuous derivatives, which is why RK4 measures the same observed
  order on trilinear and on cubic.
- `interpolation.upsample` sets the spacing at which the field is sampled before the
  interpolants are built ($1/u$ voxels). The coordinate frame does not change, so positions,
  step sizes and lengths stay in native voxel units and runs at different factors are directly
  comparable. Note the $u \to \infty$ limit is the native-resolution interpolant, not
  ground-truth anatomy.

With `field: csd`, the direction comes from CSD FOD peaks rather than the tensor, and the peak
nearest the incoming tangent is selected. That selection is structural — it is what reduces a
multi-valued field to one direction — not a tunable steering term.

**MMF integration** advances $d\mathbf{x}/ds = \mathbf{e}_1$ while evolving the full carried
frame by the structure equation, driven by the interpolated connection field. Here
`integrator` chooses only the numerical scheme (`rk4` fixed or `rkf45` adaptive
Dormand-Prince); the direction comes entirely from the connection form.

#### **Termination**

All three trackers share the same criteria, expressed in step-invariant units:

- $\text{FA}(\mathbf{r}) < \texttt{termination.fa\_min}$
- Curvature exceeds `termination.angle_max`, in **degrees per voxel of arc**. This fixes a
  minimum radius of curvature $R = 57.3 / \texttt{angle\_max}$ voxels, so refining the step
  does not loosen the constraint. Because the principal direction is a line field, tangents
  are sign-aligned and a measured turn never exceeds 90°; a budget above $90/\Delta s$ is
  therefore **inert**, not merely loose. `angle_max: 0` disables the criterion.
- Arc length exceeds `termination.max_arc` voxels (`max_steps` is derived as
  $\lceil \texttt{max\_arc} / \Delta s \rceil$)
- The streamline leaves the brain mask
- With ACT enabled, the streamline reaches GM (accept) or CSF (reject)

Completed tracks shorter than `termination.min_arc` voxels of arc are discarded.

#### **Filtering and Output**

After tracking, `nim_filter_tracks_roi` applies the configured region gates: `include_roi` and
`exclude_roi` are waypoint tests (does the track pass through), while `endpoints_in` and
`contained_in` are the two halves of a bundle definition — where a track *stops*, and whether
it stays inside a corridor. `output.arc_step` optionally resamples the saved streamlines to a
fixed arc spacing, which decouples file size from step size without affecting integration
accuracy.

### 4. Mathematical Framework

#### **Diffusion Signal Model**
The Stejskal-Tanner equation:

$$
S(b, \mathbf{g}) = S_0 \exp\!\left(-b \, \mathbf{g}^\top \mathbf{D} \, \mathbf{g}\right)
$$

#### **Tensor Metrics**

- **Mean Diffusivity**: $\text{MD} = (\lambda_1 + \lambda_2 + \lambda_3) / 3$
- **Radial Diffusivity**: $\text{RD} = (\lambda_2 + \lambda_3) / 2$
- **Axial Diffusivity**: $\text{AD} = \lambda_1$

#### **Field Map Correction Theory**
Susceptibility-induced distortion model:

$$
k_{\text{actual}} = k_{\text{nominal}} + \gamma \, \Delta B_0 \, T_E
$$

where:

- $\gamma$: Gyromagnetic ratio ($42.58$ MHz/T)
- $T_E$: Echo time
- $\Delta B_0$: $B_0$ field inhomogeneity

## Data Flow

```
                        main.m                          runTractography.m
  ┌─────────────────────────────────────────────┐   ┌────────────────────────┐
  │ nim_read                                    │   │ seed mask resolution   │
  │   → nim_dt_spd → nim_eig → nim_fa           │   │   → tracker dispatch   │
  │   → nim_mmf_geometry  (frame + connection)  │──▶│   → ROI filtering      │
  │   → nim_registration  (optional, needs T1)  │   │   → arc resampling     │
  │   → nim_parcellation                        │   │   → tracks_*.mat       │
  │   → mask improvement                        │   └────────────────────────┘
  │   → tissue segmentation (WM/GM/CSF for ACT) │                │
  │   → nim_save                                │                ▼
  └─────────────────────────────────────────────┘        visualization /
             ▲                                            scoring / export
             │
   nim_preprocessing (FSL): b0 → brain extraction → denoise → fieldmap
                          → motion → eddy → T1 registration → atlas
```

## File Naming Convention

**Input Files:**

- `{name}_raw.nii.gz` - Raw DWI data
- `{name}.bvec` - B-vectors
- `{name}.bval` - B-values
- `{name}_T1.nii.gz` - T1 structural data (optional, enables enhanced processing)
- `{name}_fmap_Hz.nii.gz` - Field map in Hz (optional)
- `{name}_acqp.txt` - Acquisition parameters (optional)
- `{name}_index.txt` - Volume indices (optional)

**Output Files:**

- `{name}.nii.gz` - Preprocessed DWI data
- `{name}_M.nii.gz` - Brain mask
- `{name}_WM_mask.nii.gz`, `{name}_GM_mask.nii.gz`, `{name}_CSF_mask.nii.gz` - Tissue masks for ACT
- `parcellation_mask.nii.gz` - Atlas resampled into DWI space
- `{name}_preprocessing_report.mat` - Processing report

## Running the Pipeline

### Requirements

**MATLAB Toolboxes:**

- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Tools for NIfTI and ANALYZE image

**External Software:**

- **SPM12** — NIfTI I/O and an optional registration backend. Not vendored here; install it
  at `lib/spm12/`, the path `main.m` adds.
- **FSL** — required for preprocessing, and must be initialized in the shell before
  `matlab -batch` runs.
- **`lib/bfgs/`** — vendored; used by the SPD-constrained tensor fit.

### Data Preparation

To run from raw data you must provide:

- `{prefix}_raw.nii.gz` - Raw NIfTI file
- `{prefix}.bvec` - B-vector file
- `{prefix}.bval` - B-value file
- `{prefix}_T1.nii.gz` - T1 structural data (optional; enables T1-based brain extraction,
  registration and FAST tissue segmentation)

### Invocation

```bash
# Full pipeline, via the shell launcher
./bin/run_hinec.sh data/ismrm2015/ismrm2015 ismrm2015.mat config/ismrm2015.yml

# Tractography only, against the nim that run produced
./bin/run_tractography.sh hinec_dti --score
```

```matlab
% Or from MATLAB
config   = load_config_yaml('config/hinec_default.yml');
run_info = create_run_directory('config/hinec_default.yml');
main('input_data/my_data', 'output.mat', config, run_info);
runTractography(fullfile(run_info.output_dir, 'output.mat'), config, run_info);
```

`main.m` short-circuits aggressively: if the output `.mat` exists it skips everything, and if
the preprocessed `.nii.gz` exists it skips preprocessing. Delete the target file to force
reprocessing.

### Viewing Results

```matlab
load('output.mat');
nim_plot(nim, 'mode', 'parcels');                    % Parcellation
visualizeTractography('tracks.mat', 'output.mat');   % Tractography
```

## Development Notes

### Project Structure

- `main.m` orchestrates preprocessing, registration, parcellation and tissue segmentation;
  `runTractography.m` handles seeding, tracking and filtering
- Source modules live under `src/`: `nim_preprocessing/`, `nim_registration/`,
  `nim_calculation/`, `nim_parcellation/`, `nim_tractography/`, `nim_visualization/`,
  `nim_plots/`, `nim_utils/`, `nim_challenges/`
- Sample data: `data/original_sample/` (preprocessed) and `data/parcellation_sample/` (raw)

### Coding Style

- Four-space indentation
- One MATLAB function per file, named identically to the function
- lowerCamelCase for pipeline entry points (`runTractography`, `visualizeTractography`)
- snake_case for utilities (`nim_diagnostic_check`)
- Descriptive variable names, uppercase for constants
- Structure arguments preferred over long positional lists

### Testing

Tests live under `tests/` in four groups — `unit/`, `integration/`, `fsl/` and `nim_tests/`.
Run one class with MATLAB's runner:

```matlab
results = runtests('tests/unit/TestNimFa.m');
```

`tests/test_yaml_config.m` is a standalone script that exercises every YAML preset in
`config/`. There is no top-level `make test` target; the `Makefile` only drives `mkdocs`.