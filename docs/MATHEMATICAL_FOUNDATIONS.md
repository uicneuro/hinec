# Mathematical Foundations in HINEC

This document consolidates all mathematical formulas, derivations, and numerical methods used throughout the HINEC pipeline. It serves as a self-contained reference for researchers who want to understand the mathematics without reading source code.

---

## 1. DTI Signal Model

### Stejskal-Tanner Equation

The diffusion-weighted MRI signal measured with gradient direction **g** and b-value *b* is:

```
S(g, b) = S₀ · exp(-b · gᵀDg)
```

where:
- **S₀** is the non-diffusion-weighted signal (b = 0 image)
- **g** = [gx, gy, gz]ᵀ is the unit gradient direction vector
- **D** is the 3×3 symmetric diffusion tensor
- *b* is the b-value (s/mm²), encoding gradient strength and timing

### Log-Linear Formulation

Taking the natural logarithm and rearranging:

```
ln(S₀ / S) / b = gᵀDg
```

Expanding the quadratic form gᵀDg with the symmetric tensor:

```
gᵀDg = Dxx·gx² + Dyy·gy² + Dzz·gz² + 2Dxy·gx·gy + 2Dyz·gy·gz + 2Dxz·gz·gx
```

This is linear in the 6 unique tensor elements. Defining:
- **Y** = ln(S₀ / S) / b (the apparent diffusion coefficient along direction **g**)
- **d** = [Dxx, Dyy, Dzz, Dxy, Dyz, Dxz]ᵀ (6-element tensor vector)

We obtain the linear system:

```
Y = Hd
```

### Design Matrix Construction

The design matrix **H** (N × 6, where N = number of gradient directions) is constructed from gradient directions:

```
H = | gx₁²   gy₁²   gz₁²   2·gx₁·gy₁   2·gy₁·gz₁   2·gz₁·gx₁ |
    | gx₂²   gy₂²   gz₂²   2·gx₂·gy₂   2·gy₂·gz₂   2·gz₂·gx₂ |
    |  ...     ...     ...      ...          ...          ...       |
    | gxN²   gyN²   gzN²   2·gxN·gyN   2·gyN·gzN   2·gzN·gxN |
```

**Implementation** (`nim_dt.m:27-31`, `nim_dt_spd.m:27`): Both tensor calculation functions construct H identically:
```matlab
H = [gx.^2  gy.^2  gz.^2  2.*gx.*gy  2.*gy.*gz  2.*gz.*gx];
```

The signal Y is computed from log-ratio of b0 and diffusion-weighted images (`nim_dt_spd.m:22-24`):
```matlab
Y(:,:,:,t) = log(nim.img_b0 ./ nim.img_bi(:,:,:,t)) ./ b(t);
```

---

## 2. Tensor Estimation

### Least Squares Solution

Given the overdetermined system Y = Hd (typically N = 30-60 equations, 6 unknowns), the least squares solution minimizes ||Y - Hd||²:

```
d = (HᵀH)⁻¹Hᵀy
```

In MATLAB, this is computed via the backslash operator (`nim_dt.m:48`, `nim_dt_spd.m:73`):
```matlab
D = H \ Y_i;
```

which uses QR decomposition internally for numerical stability.

### The Non-Positive-Definite Problem

The diffusion tensor must be symmetric positive-definite (SPD) — all eigenvalues must be positive — because diffusion coefficients are physical quantities that cannot be negative. However, noise in the measured signal can produce tensors with negative eigenvalues.

In `nim_dt.m`, non-SPD tensors are noted but not corrected. In `nim_dt_spd.m`, they are corrected using BFGS optimization.

### SPD Constraint via BFGS Optimization

When least squares produces a non-SPD tensor (any eigenvalue < 0), `nim_dt_spd.m` invokes BFGS optimization (`vox_dt_bfgs.m`) to find the nearest SPD tensor that fits the data.

**BFGS** (Broyden-Fletcher-Goldfarb-Shanno) is a quasi-Newton optimization algorithm that:

1. Starts from the least squares solution
2. Parameterizes the tensor using Cholesky decomposition (D = LLᵀ), which guarantees positive-definiteness for any L
3. Minimizes the fitting error ||Y - Hd||² subject to the SPD constraint
4. Updates an approximate Hessian matrix at each iteration for efficient convergence

**Implementation flow** (`nim_dt_spd.m:76-98`):
```
1. Compute LSF solution: D = H \ Y
2. Check eigenvalues: [Q, l] = eig(reshape(D))
3. If any(l < 0):
   a. Run BFGS: [D, nsteps] = vox_dt_bfgs(nim, x, y, z)
   b. Re-compute eigenvalues from optimized tensor
4. Validate: reject if NaN, Inf, or eigenvalues outside [0, 0.01]
5. Sort eigenvalues descending: [lM, ilM] = maxk(l, 3)
6. Store sorted eigenvectors and eigenvalues
```

### 6-Element Tensor Representation

The symmetric 3×3 tensor is stored as a 6-element vector:

```
d = [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
```

The reconstruction to a full matrix is (`nim_reshape_d.m`):

```
D = | d(1)  d(4)  d(5) |     | Dxx  Dxy  Dxz |
    | d(4)  d(2)  d(6) |  =  | Dxy  Dyy  Dyz |
    | d(5)  d(6)  d(3) |     | Dxz  Dyz  Dzz |
```

---

## 3. Eigendecomposition

### Symmetric Matrix Eigenvalue Problem

For a real symmetric 3×3 matrix D, the eigenvalue decomposition is:

```
Dv = λv
```

yielding three real eigenvalue-eigenvector pairs: (λ₁, **e₁**), (λ₂, **e₂**), (λ₃, **e₃**).

**Sorting convention**: Eigenvalues are sorted in descending order (`nim_eig.m`, `nim_dt_spd.m:104`):

```
λ₁ ≥ λ₂ ≥ λ₃ ≥ 0
```

using MATLAB's `maxk()`:
```matlab
[lM, ilM] = maxk(l, 3);
QM = Q(:, ilM);
```

### Physical Interpretation

| Component | Symbol | Physical Meaning |
|---|---|---|
| Primary eigenvalue | λ₁ | Axial diffusivity (along fiber) |
| Primary eigenvector | **e₁** | Estimated fiber direction |
| Secondary eigenvalue | λ₂ | Intermediate diffusivity |
| Secondary eigenvector | **e₂** | Secondary diffusion axis |
| Tertiary eigenvalue | λ₃ | Minimum diffusivity (perpendicular to fiber) |
| Tertiary eigenvector | **e₃** | Direction of maximum restriction |

### Tensor Reconstruction from Eigendecomposition

The tensor can be reconstructed from its eigendecomposition:

```
D = λ₁·e₁e₁ᵀ + λ₂·e₂e₂ᵀ + λ₃·e₃e₃ᵀ = QΛQᵀ
```

where Q = [**e₁** **e₂** **e₃**] and Λ = diag(λ₁, λ₂, λ₃).

### Tensor Geometry

The eigenvalues define the shape of the diffusion ellipsoid:

- **Prolate** (λ₁ >> λ₂ ≈ λ₃): Cigar-shaped, coherent single-fiber bundle
- **Oblate** (λ₁ ≈ λ₂ >> λ₃): Pancake-shaped, crossing fibers or planar structure
- **Spherical** (λ₁ ≈ λ₂ ≈ λ₃): Isotropic diffusion, CSF or gray matter

---

## 4. Diffusion Scalar Metrics

### Fractional Anisotropy (FA)

FA quantifies the degree of anisotropy on a normalized scale [0, 1]:

```
FA = √(1/2) · √[(λ₁ - λ₂)² + (λ₂ - λ₃)² + (λ₃ - λ₁)²] / √(λ₁² + λ₂² + λ₃²)
```

**Derivation**: FA is the normalized variance of the eigenvalues. Defining the mean eigenvalue λ̄ = (λ₁ + λ₂ + λ₃)/3:

```
FA = √(3/2) · √[(λ₁ - λ̄)² + (λ₂ - λ̄)² + (λ₃ - λ̄)²] / √(λ₁² + λ₂² + λ₃²)
```

Both formulations are mathematically equivalent. The first is used in the implementation.

**Implementation** (`nim_fa.m:23-24`):
```matlab
denom = sqrt(l(1)^2 + l(2)^2 + l(3)^2);
if denom > 0
    nim.FA(x,y,z) = sqrt(1/2) * sqrt((l(1)-l(2))^2 + (l(2)-l(3))^2 + (l(3)-l(1))^2) / denom;
end
```

Division by zero is protected: FA = 0 when all eigenvalues are zero.

### Mean Diffusivity (MD)

```
MD = (λ₁ + λ₂ + λ₃) / 3 = trace(D) / 3
```

Typical values: ~0.7-1.0 × 10⁻³ mm²/s in healthy white matter, ~3.0 × 10⁻³ mm²/s in CSF.

### Axial Diffusivity (AD) and Radial Diffusivity (RD)

```
AD = λ₁
RD = (λ₂ + λ₃) / 2
```

Note: FA can be expressed in terms of AD and RD:
```
FA = √(1/2) · √[(AD - RD)² · 3] / √(AD² + 2·RD²)    (when λ₂ = λ₃ = RD)
```

---

## 5. FACT Tractography Algorithm

The FACT (Fiber Assignment by Continuous Tracking) algorithm, implemented in `nim_tractography_standard.m`, uses discrete voxel directions with exact boundary intersection.

### Algorithm Overview

```
Input: Seed point p₀, eigenvector field e₁(x,y,z), FA field
Output: Track T = [p₀, p₁, p₂, ...]

1. Initialize: p = p₀, previous_direction = e₁(voxel(p₀))
2. While termination criteria not met:
   a. Get current voxel indices: (i, j, k) = floor(p)
   b. Get direction: d = e₁(i, j, k)  [discrete, no interpolation]
   c. Ensure consistency: if dot(d, previous_direction) < 0, flip d = -d
   d. Compute boundary exit: p_exit = find_voxel_boundary_exit(p, d)
   e. Step to boundary: p = p_exit + ε·d  [enter next voxel]
   f. Check termination: FA(new voxel) < threshold? angle(d, d_prev) > limit?
   g. Store point: append p to track T
   h. Update: previous_direction = d
3. Repeat steps 1-2 tracking in opposite direction from p₀
4. Concatenate: T = [reverse(T_backward), T_forward]
```

### Voxel Boundary Intersection

Given position **p** = (px, py, pz) and direction **d** = (dx, dy, dz) within a voxel with boundaries [xmin, xmax] × [ymin, ymax] × [zmin, zmax]:

For each axis, compute the parametric distance to the near and far boundaries:

```
tx_near = (xmin - px) / dx    tx_far = (xmax - px) / dx    (swap if dx < 0)
ty_near = (ymin - py) / dy    ty_far = (ymax - py) / dy    (swap if dy < 0)
tz_near = (zmin - pz) / dz    tz_far = (zmax - pz) / dz    (swap if dz < 0)
```

The exit parameter is the minimum positive far-boundary parameter:

```
t_exit = min(tx_far, ty_far, tz_far)    where t > 0
```

The exit point is:

```
p_exit = p + t_exit · d
```

### Direction Consistency

When entering a new voxel, the eigenvector e₁ may point in the opposite direction (eigenvectors are bidirectional). The dot product check ensures consistent propagation:

```
if dot(d_new, d_prev) < 0:
    d_new = -d_new
```

### Bidirectional Tracking

Each seed point generates two half-tracks: one forward, one backward. The backward track starts by negating the initial direction. The final track is the concatenation:

```
T = [reverse(T_backward); T_forward]
```

---

## 6. Numerical Integration Methods

HINEC's high-order tractography (`nim_tractography_hinec.m`) models streamline tracking as an initial value ODE:

```
dy/ds = f(y)    where y = [x, y, z] position, s = arc length, f = interpolated direction field
```

### Euler Method (Order 1, `integration_order = 1`)

The simplest explicit method:

```
y_{n+1} = y_n + h · f(y_n)
```

- Error per step: O(h²) — local truncation error
- Global error: O(h) — first-order accuracy
- One function evaluation per step
- Fast but least accurate; equivalent to standard step-based tracking

### RK2 Midpoint Method (Order 2, `integration_order = 2`)

Uses a trial step to the midpoint for better accuracy:

```
k₁ = f(y_n)
k₂ = f(y_n + h/2 · k₁)
y_{n+1} = y_n + h · k₂
```

- Error per step: O(h³)
- Global error: O(h²) — second-order accuracy
- Two function evaluations per step
- Significantly more accurate than Euler for curved trajectories

### RK4 Classical Method (Order 4, `integration_order = 4`)

The standard fourth-order Runge-Kutta method:

```
k₁ = f(y_n)
k₂ = f(y_n + h/2 · k₁)
k₃ = f(y_n + h/2 · k₂)
k₄ = f(y_n + h · k₃)
y_{n+1} = y_n + h/6 · (k₁ + 2k₂ + 2k₃ + k₄)
```

- Error per step: O(h⁵)
- Global error: O(h⁴) — fourth-order accuracy
- Four function evaluations per step
- Excellent accuracy-to-cost ratio; the default HINEC method

### RKF45 Adaptive Method (Order 4-5, `integration_order = 5`)

The Dormand-Prince embedded pair provides automatic step size control by computing both a 4th-order and 5th-order solution using shared function evaluations.

**Stage evaluations** (7 stages):

```
k₁ = f(y_n)
k₂ = f(y_n + h·a₂₁·k₁)
k₃ = f(y_n + h·(a₃₁·k₁ + a₃₂·k₂))
k₄ = f(y_n + h·(a₄₁·k₁ + a₄₂·k₂ + a₄₃·k₃))
k₅ = f(y_n + h·(a₅₁·k₁ + a₅₂·k₂ + a₅₃·k₃ + a₅₄·k₄))
k₆ = f(y_n + h·(a₆₁·k₁ + a₆₂·k₂ + a₆₃·k₃ + a₆₄·k₄ + a₆₅·k₅))
k₇ = f(y_n + h·(a₇₁·k₁ + ... + a₇₆·k₆))
```

**Solutions**:

```
y₅ = y_n + h·(b₁·k₁ + b₃·k₃ + b₄·k₄ + b₅·k₅ + b₆·k₆ + b₇·k₇)     (5th order)
y₄ = y_n + h·(b̂₁·k₁ + b̂₃·k₃ + b̂₄·k₄ + b̂₅·k₅ + b̂₆·k₆)              (4th order)
```

**Error estimate**:

```
err = ||y₅ - y₄|| = h·||Σᵢ (bᵢ - b̂ᵢ)·kᵢ||
```

**Step size control**:

```
h_new = safety · h · (tolerance / err)^(1/5)
h_new = max(step_min, min(step_max, h_new))
```

If err ≤ tolerance: accept step, use y₅, adjust h for next step.
If err > tolerance: reject step, reduce h, retry.

**HINEC defaults**: tolerance = 0.01 voxels, safety = 0.9, step_min = 0.01, step_max = 1.0.

For the complete Dormand-Prince coefficients and derivation, see [RKF.md](RKF.md). For usage examples and parameter guidelines, see [RKF_Usage.md](RKF_Usage.md).

### Comparison of Methods

| Method | Order | Evaluations/Step | Relative Cost | Relative Accuracy | Best For |
|---|---|---|---|---|---|
| Euler | 1 | 1 | 1× | Low | Quick exploration |
| RK2 | 2 | 2 | 2× | Moderate | Fast, reasonable accuracy |
| RK4 | 4 | 4 | 4× | High | Default, good cost/accuracy |
| RKF45 | 4-5 | 6-7 | 6-7× | Highest (adaptive) | Publication quality |

The higher-order methods are especially valuable in regions of high fiber curvature, where lower-order methods accumulate significant tracking error.

---

## 7. Interpolation Methods

### Trilinear Interpolation

The default interpolation method in HINEC. For a query point **p** = (x, y, z) inside a voxel, the interpolated value is a weighted average of the 8 surrounding voxel centers.

Given voxel centers at integer coordinates with values V(i,j,k), and fractional coordinates:

```
xd = x - floor(x),  yd = y - floor(y),  zd = z - floor(z)
```

The trilinear interpolation is:

```
V(x,y,z) = V₀₀₀(1-xd)(1-yd)(1-zd) + V₁₀₀(xd)(1-yd)(1-zd) +
            V₀₁₀(1-xd)(yd)(1-zd)   + V₁₁₀(xd)(yd)(1-zd)   +
            V₀₀₁(1-xd)(1-yd)(zd)   + V₁₀₁(xd)(1-yd)(zd)   +
            V₀₁₁(1-xd)(yd)(zd)     + V₁₁₁(xd)(yd)(zd)
```

**Implementation**: HINEC uses MATLAB's `griddedInterpolant` with `'linear'` method for 2-5× faster interpolation compared to `interp3` (`nim_tractography_hinec.m:239-242`):

```matlab
nim.FA_interp = griddedInterpolant(grid_vectors, nim.FA, 'linear', 'none');
nim.v1_x_interp = griddedInterpolant(grid_vectors, nim.v1_x, 'linear', 'none');
nim.v1_y_interp = griddedInterpolant(grid_vectors, nim.v1_y, 'linear', 'none');
nim.v1_z_interp = griddedInterpolant(grid_vectors, nim.v1_z, 'linear', 'none');
```

Each eigenvector component is interpolated independently, and the result is renormalized to unit length.

### Cubic Spline Interpolation

Available via `interp_method = 'cubic'`. Uses cubic B-splines for C²-continuous interpolation, providing smoother direction fields at the cost of additional computation. Uses the same `griddedInterpolant` framework with `'cubic'` method.

### High-Order Spectral Interpolation

The `nim_interp.m` utility implements spectral element interpolation using Gauss-Lobatto-Legendre (GLL) quadrature nodes:

1. Expand eigenvectors from cell centers to vertices
2. Place (p+1)³ interpolation nodes per voxel using GLL quadrature
3. Perform 3D spline interpolation within each voxel

GLL nodes are computed by `zwgll.m`: roots of (1-x²)P'_N(x) = 0, where P'_N is the derivative of the Nth Legendre polynomial. These nodes cluster near element boundaries, providing optimal interpolation accuracy with minimal Runge phenomenon.

Uniform nodes from `zwuni.m` provide an alternative with equally-spaced points.

### Direction Field Normalization

After interpolation, the direction vector must be renormalized:

```
d_interp = interpolate(v1_x, v1_y, v1_z, position)
d_normalized = d_interp / ||d_interp||
```

This is necessary because linear interpolation of unit vectors does not produce unit vectors.

---

## 8. Registration Mathematics

### Affine Transformation

An affine transformation maps coordinates between spaces using a 4×4 matrix (12 degrees of freedom):

```
| x' |   | a₁₁  a₁₂  a₁₃  t₁ |   | x |
| y' | = | a₂₁  a₂₂  a₂₃  t₂ | · | y |
| z' |   | a₃₁  a₃₂  a₃₃  t₃ |   | z |
| 1  |   |  0    0    0    1  |   | 1 |
```

Decomposed into rotation (3 DOF), translation (3 DOF), scaling (3 DOF), and shear (3 DOF).

### Registration Cost Functions

HINEC supports two registration backends:

**FSL FLIRT** uses correlation ratio for same-modality and mutual information for cross-modality registration:

```
Correlation Ratio (CR) = 1 - Var(Y|X) / Var(Y)
```

where Y is the target image intensity and X is the source image intensity, partitioned into bins.

**Normalized Mutual Information** (used for quality assessment in `registration_utils.m`):

```
NMI = (H(A) + H(B)) / H(A, B)
```

where H(A) is the entropy of image A, H(B) is the entropy of image B, and H(A,B) is their joint entropy. Computed using 64-bin histograms.

### Multi-Modal Registration Chain

HINEC constructs a chain of transforms between three coordinate spaces:

```
DWI native  <---(6 DOF FLIRT)---->  T1 native  <---(12 DOF FLIRT + FNIRT)---->  MNI standard
```

**Forward chain** (DWI → MNI): Compose transforms
```
T_total = T_t1_to_mni · T_dti_to_t1
```

**Inverse chain** (MNI → DWI): Apply inverse transforms in reverse order
```
T⁻¹_total = T⁻¹_dti_to_t1 · T⁻¹_t1_to_mni
```

This chain is used by `nim_apply_transforms.m` to bring atlas labels from MNI space to DWI space.

### Nearest-Neighbor Resampling for Labels

When transforming parcellation labels, nearest-neighbor interpolation is required to preserve discrete integer values:

```
label(x', y', z') = label(round(x), round(y), round(z))
```

Linear or spline interpolation would produce non-integer values that corrupt label identity.

---

## 9. Tissue Segmentation Thresholds

FA-based tissue classification (`preproc_tissue_segmentation.m`):

```
WM:  FA(x,y,z) > 0.2         (eroded by 1-voxel sphere)
GM:  0.05 < FA(x,y,z) ≤ 0.2
CSF: FA(x,y,z) ≤ 0.05
```

The WM mask is eroded (morphological erosion with spherical kernel, radius=1) to remove boundary voxels where partial volume effects may place gray matter voxels inside the white matter mask.

Masks are mutually exclusive and constrained by the brain mask:

```
WM ∩ GM = ∅,  WM ∩ CSF = ∅,  GM ∩ CSF = ∅
WM ∪ GM ∪ CSF ⊆ brain_mask
```

Expected tissue proportions: WM ~40-45%, GM ~40-45%, CSF ~10-20% of brain volume.

---

## 10. Connectivity Matrix Construction

Given a set of tracks T = {t₁, t₂, ..., tM} and a parcellation with N regions (`nim_plot_connectivity_matrix.m`):

```
For each track tₖ:
    start_region = parcellation_mask(tₖ[1])      # First point
    end_region   = parcellation_mask(tₖ[end])     # Last point

    if start_region ≠ 0 and end_region ≠ 0:
        C(start_region, end_region) += 1
        if symmetric:
            C(end_region, start_region) += 1

if normalize:
    C = C / max(C(:))
```

The resulting matrix C ∈ ℝ^(N×N) represents the structural connectivity between brain regions, where C(i,j) is the number of streamlines connecting region i to region j.

---

## Cross-References

- RKF45 complete derivation and coefficients: [RKF.md](RKF.md)
- RKF45 usage guide and parameter selection: [RKF_Usage.md](RKF_Usage.md)
- High-order methods implementation details: [High_Order.md](High_Order.md)
- Scientific context for these methods: [SCIENCE.md](SCIENCE.md)
- Architecture and data flow: [ARCHITECTURE.md](ARCHITECTURE.md)
- Complete function reference: [API_REFERENCE.md](API_REFERENCE.md)
