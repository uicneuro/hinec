# Mathematical Foundations in HINEC

This document consolidates all mathematical formulas, derivations, and numerical methods used throughout the HINEC pipeline. It serves as a self-contained reference for researchers who want to understand the mathematics without reading source code.

---

## 1. Diffusion Physics and Signal Model

### The Bloch-Torrey Equation

The physical foundation of diffusion MRI is the Bloch-Torrey equation, which extends the Bloch equation to include the effects of molecular diffusion on the transverse magnetization $M_{xy}$:

$$
\frac{\partial M_{xy}}{\partial t} = -i \gamma \mathbf{G}(t) \cdot \mathbf{r} \, M_{xy} + \nabla \cdot (\mathbf{D} \nabla M_{xy})
$$

where:

- $\gamma$ is the gyromagnetic ratio ($\gamma / 2\pi \approx 42.58$ MHz/T for hydrogen)
- $\mathbf{G}(t)$ is the time-varying magnetic field gradient
- $\mathbf{r}$ is the spatial position vector
- $\mathbf{D}$ is the local diffusion tensor

The first term encodes spatial position via phase accumulation from the applied gradient. The second term models signal attenuation from molecular diffusion. In regions where diffusion is anisotropic (e.g., white matter), $\mathbf{D}$ is a full $3 \times 3$ tensor rather than a scalar.

### Stejskal-Tanner Equation

For a pulsed gradient spin-echo (PGSE) sequence, integrating the Bloch-Torrey equation under the narrow-pulse approximation yields the Stejskal-Tanner equation:

$$
S(\mathbf{g}, b) = S_0 \exp(-b \, \mathbf{g}^\top \mathbf{D} \mathbf{g})
$$

where:

- $S_0$ is the non-diffusion-weighted signal ($b = 0$ image)
- $\mathbf{g} = [g_x, g_y, g_z]^\top$ is the unit gradient direction vector
- $\mathbf{D}$ is the $3 \times 3$ symmetric positive-definite diffusion tensor
- $b$ is the b-value (s/mm$^2$), encoding gradient strength and timing

### The B-value

The b-value is not a single parameter but a composite that encodes the gradient waveform:

$$
b = \gamma^2 G^2 \delta^2 \left( \Delta - \frac{\delta}{3} \right)
$$

where:

- $G$ is the gradient amplitude (T/m)
- $\delta$ is the gradient pulse duration (s)
- $\Delta$ is the time between the leading edges of the two gradient pulses (s)
- The factor $(\Delta - \delta/3)$ is the effective diffusion time, accounting for the finite pulse width

Typical clinical DTI uses $b \approx 1000$ s/mm$^2$, corresponding to $G \approx 40$ mT/m, $\delta \approx 30$ ms, $\Delta \approx 40$ ms.

### Apparent Diffusion Coefficient

The quadratic form $\mathbf{g}^\top \mathbf{D} \mathbf{g}$ represents the **apparent diffusion coefficient** (ADC) along direction $\mathbf{g}$:

$$
\text{ADC}(\mathbf{g}) = \mathbf{g}^\top \mathbf{D} \mathbf{g} = \sum_{i=1}^{3} \sum_{j=1}^{3} g_i \, D_{ij} \, g_j
$$

The ADC is the directional projection of the diffusion tensor --- it measures how fast water diffuses along a specific direction. In white matter, the ADC is large along the fiber direction and small perpendicular to it.

### Log-Linear Formulation

Taking the natural logarithm of the Stejskal-Tanner equation and rearranging:

$$
\frac{\ln(S_0 / S)}{b} = \mathbf{g}^\top \mathbf{D} \mathbf{g}
$$

Expanding the quadratic form with the symmetric tensor $\mathbf{D}$:

$$
\mathbf{g}^\top \mathbf{D} \mathbf{g} = D_{xx} g_x^2 + D_{yy} g_y^2 + D_{zz} g_z^2 + 2D_{xy} g_x g_y + 2D_{yz} g_y g_z + 2D_{xz} g_z g_x
$$

This expression is linear in the 6 unique tensor elements. Defining:

- $\mathbf{Y}$: the vector of observed ADC values, $Y_i = \ln(S_0 / S_i) / b_i$
- $\mathbf{d} = [D_{xx}, D_{yy}, D_{zz}, D_{xy}, D_{yz}, D_{xz}]^\top$: the 6 independent tensor elements

We obtain the linear system:

$$
\mathbf{Y} = \mathbf{H} \mathbf{d}
$$

### Design Matrix Construction

The design matrix $\mathbf{H} \in \mathbb{R}^{N \times 6}$ (where $N$ is the number of gradient directions) is:

$$
\mathbf{H} = \begin{bmatrix}
g_{x_1}^2 & g_{y_1}^2 & g_{z_1}^2 & 2g_{x_1}g_{y_1} & 2g_{y_1}g_{z_1} & 2g_{z_1}g_{x_1} \\
g_{x_2}^2 & g_{y_2}^2 & g_{z_2}^2 & 2g_{x_2}g_{y_2} & 2g_{y_2}g_{z_2} & 2g_{z_2}g_{x_2} \\
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots \\
g_{x_N}^2 & g_{y_N}^2 & g_{z_N}^2 & 2g_{x_N}g_{y_N} & 2g_{y_N}g_{z_N} & 2g_{z_N}g_{x_N}
\end{bmatrix}
$$

**Implementation** (`nim_dt.m:27-31`, `nim_dt_spd.m:27`): Both tensor calculation functions construct $\mathbf{H}$ identically:
```matlab
H = [gx.^2  gy.^2  gz.^2  2.*gx.*gy  2.*gy.*gz  2.*gz.*gx];
```

The signal $\mathbf{Y}$ is computed from the log-ratio (`nim_dt_spd.m:22-24`):
```matlab
Y(:,:,:,t) = log(nim.img_b0 ./ nim.img_bi(:,:,:,t)) ./ b(t);
```

---

## 2. Tensor Estimation

### Ordinary Least Squares Solution

Given the overdetermined system $\mathbf{Y} = \mathbf{H}\mathbf{d}$ (typically $N = 30$--$60$ equations, 6 unknowns), the ordinary least squares (OLS) solution minimizes the sum of squared residuals:

$$
\hat{\mathbf{d}} = \arg\min_{\mathbf{d}} \| \mathbf{Y} - \mathbf{H}\mathbf{d} \|^2
$$

Setting the gradient to zero yields the normal equations:

$$
\mathbf{H}^\top \mathbf{H} \hat{\mathbf{d}} = \mathbf{H}^\top \mathbf{Y}
$$

with solution:

$$
\hat{\mathbf{d}} = (\mathbf{H}^\top \mathbf{H})^{-1} \mathbf{H}^\top \mathbf{Y}
$$

In MATLAB, this is computed via the backslash operator (`nim_dt.m:48`, `nim_dt_spd.m:73`):
```matlab
D = H \ Y_i;
```
which internally uses QR decomposition for numerical stability rather than explicitly forming $(\mathbf{H}^\top \mathbf{H})^{-1}$.

### Weighted Least Squares

In the log-linear model, the noise in $\mathbf{Y}$ is heteroscedastic --- the variance depends on the signal magnitude. A weighted least squares (WLS) formulation accounts for this:

$$
\hat{\mathbf{d}}_{\text{WLS}} = (\mathbf{H}^\top \mathbf{W} \mathbf{H})^{-1} \mathbf{H}^\top \mathbf{W} \mathbf{Y}
$$

where $\mathbf{W} = \text{diag}(S_1^2, S_2^2, \ldots, S_N^2)$ weights each observation by its signal intensity squared. HINEC currently uses OLS, which is a reasonable approximation at typical clinical SNR levels.

### The Non-Positive-Definite Problem

The diffusion tensor must be symmetric positive-definite (SPD) --- all eigenvalues must be strictly positive --- because diffusion coefficients are physical quantities that cannot be negative. However, noise in the measured signal can produce tensors with negative eigenvalues via least squares fitting.

A tensor $\mathbf{D}$ is SPD if and only if:

$$
\mathbf{x}^\top \mathbf{D} \mathbf{x} > 0 \quad \forall \, \mathbf{x} \neq \mathbf{0}
$$

Equivalently, $\mathbf{D}$ is SPD if all its eigenvalues $\lambda_1, \lambda_2, \lambda_3 > 0$.

In `nim_dt.m`, non-SPD tensors are noted but not corrected. In `nim_dt_spd.m`, they are corrected using BFGS optimization.

### SPD Constraint via Cholesky Parameterization

The key insight for enforcing the SPD constraint is the **Cholesky decomposition**: any SPD matrix can be uniquely written as:

$$
\mathbf{D} = \mathbf{L} \mathbf{L}^\top
$$

where $\mathbf{L}$ is a lower triangular matrix with positive diagonal entries:

$$
\mathbf{L} = \begin{bmatrix}
l_{11} & 0 & 0 \\
l_{21} & l_{22} & 0 \\
l_{31} & l_{32} & l_{33}
\end{bmatrix}, \quad l_{11}, l_{22}, l_{33} > 0
$$

This gives the tensor elements:

$$
\begin{aligned}
D_{xx} &= l_{11}^2 \\
D_{yy} &= l_{21}^2 + l_{22}^2 \\
D_{zz} &= l_{31}^2 + l_{32}^2 + l_{33}^2 \\
D_{xy} &= l_{11} l_{21} \\
D_{xz} &= l_{11} l_{31} \\
D_{yz} &= l_{21} l_{31} + l_{22} l_{32}
\end{aligned}
$$

By optimizing over the 6 free parameters of $\mathbf{L}$ rather than directly over $\mathbf{D}$, positive-definiteness is guaranteed for any value of the parameters.

### BFGS Optimization

When least squares produces a non-SPD tensor, `nim_dt_spd.m` invokes BFGS optimization (`vox_dt_bfgs.m`) to find the nearest SPD tensor that fits the data. The objective function minimizes the fitting error:

$$
f(\mathbf{l}) = \| \mathbf{Y} - \mathbf{H} \, \mathbf{d}(\mathbf{l}) \|^2
$$

where $\mathbf{d}(\mathbf{l})$ maps the Cholesky parameters to the 6-element tensor vector.

**BFGS** (Broyden-Fletcher-Goldfarb-Shanno) is a quasi-Newton method that maintains an approximate inverse Hessian $\mathbf{B}_k$, updated at each iteration by:

$$
\mathbf{B}_{k+1} = \left(\mathbf{I} - \frac{\mathbf{s}_k \mathbf{y}_k^\top}{\mathbf{y}_k^\top \mathbf{s}_k}\right) \mathbf{B}_k \left(\mathbf{I} - \frac{\mathbf{y}_k \mathbf{s}_k^\top}{\mathbf{y}_k^\top \mathbf{s}_k}\right) + \frac{\mathbf{s}_k \mathbf{s}_k^\top}{\mathbf{y}_k^\top \mathbf{s}_k}
$$

where:

- $\mathbf{s}_k = \mathbf{l}_{k+1} - \mathbf{l}_k$ is the parameter step
- $\mathbf{y}_k = \nabla f(\mathbf{l}_{k+1}) - \nabla f(\mathbf{l}_k)$ is the gradient change

The search direction is $\mathbf{p}_k = -\mathbf{B}_k \nabla f(\mathbf{l}_k)$, followed by a line search to determine the step length.

**Implementation flow** (`nim_dt_spd.m:76-98`):

1. Compute OLS solution: $\mathbf{D} = \mathbf{H} \backslash \mathbf{Y}$
2. Check eigenvalues: $[\mathbf{Q}, \mathbf{\lambda}] = \text{eig}(\text{reshape}(\mathbf{D}))$
3. If $\exists \, \lambda_i < 0$: run BFGS with Cholesky parameterization
4. Validate: reject if NaN, Inf, or eigenvalues outside $[0, \, 0.01]$
5. Sort eigenvalues descending: $\lambda_1 \geq \lambda_2 \geq \lambda_3$
6. Store sorted eigenvectors and eigenvalues

### 6-Element Tensor Representation

The symmetric $3 \times 3$ tensor is stored as a 6-element vector:

$$
\mathbf{d} = [D_{xx}, D_{yy}, D_{zz}, D_{xy}, D_{xz}, D_{yz}]
$$

The reconstruction to the full matrix (`nim_reshape_d.m`) is:

$$
\mathbf{D} = \begin{bmatrix}
D_{xx} & D_{xy} & D_{xz} \\
D_{xy} & D_{yy} & D_{yz} \\
D_{xz} & D_{yz} & D_{zz}
\end{bmatrix}
$$

---

## 3. Eigendecomposition

### Symmetric Matrix Eigenvalue Problem

For a real symmetric $3 \times 3$ matrix $\mathbf{D}$, the eigenvalue decomposition solves:

$$
\mathbf{D} \mathbf{v} = \lambda \mathbf{v}
$$

The eigenvalues are roots of the **characteristic polynomial**:

$$
\det(\mathbf{D} - \lambda \mathbf{I}) = -\lambda^3 + I_1 \lambda^2 - I_2 \lambda + I_3 = 0
$$

where the tensor invariants are:

$$
\begin{aligned}
I_1 &= \text{tr}(\mathbf{D}) = D_{xx} + D_{yy} + D_{zz} \\
I_2 &= \frac{1}{2}\left[(\text{tr}\,\mathbf{D})^2 - \text{tr}(\mathbf{D}^2)\right] = D_{xx}D_{yy} + D_{yy}D_{zz} + D_{zz}D_{xx} - D_{xy}^2 - D_{yz}^2 - D_{xz}^2 \\
I_3 &= \det(\mathbf{D})
\end{aligned}
$$

These invariants are rotationally invariant --- they do not change when the coordinate system is rotated. For a real symmetric matrix, all three eigenvalues are guaranteed to be real.

This yields three real eigenvalue-eigenvector pairs: $(\lambda_1, \mathbf{e}_1)$, $(\lambda_2, \mathbf{e}_2)$, $(\lambda_3, \mathbf{e}_3)$.

**Sorting convention**: Eigenvalues are sorted in descending order (`nim_eig.m`, `nim_dt_spd.m:104`):

$$
\lambda_1 \geq \lambda_2 \geq \lambda_3 \geq 0
$$

using MATLAB's `maxk()`:
```matlab
[lM, ilM] = maxk(l, 3);
QM = Q(:, ilM);
```

### Physical Interpretation

| Component | Symbol | Physical Meaning |
|---|---|---|
| Primary eigenvalue | $\lambda_1$ | Axial diffusivity (along fiber) |
| Primary eigenvector | $\mathbf{e}_1$ | Estimated fiber direction |
| Secondary eigenvalue | $\lambda_2$ | Intermediate diffusivity |
| Secondary eigenvector | $\mathbf{e}_2$ | Secondary diffusion axis |
| Tertiary eigenvalue | $\lambda_3$ | Minimum diffusivity (perpendicular to fiber) |
| Tertiary eigenvector | $\mathbf{e}_3$ | Direction of maximum restriction |

### Tensor Reconstruction from Eigendecomposition

The tensor can be reconstructed from its spectral decomposition:

$$
\mathbf{D} = \lambda_1 \mathbf{e}_1 \mathbf{e}_1^\top + \lambda_2 \mathbf{e}_2 \mathbf{e}_2^\top + \lambda_3 \mathbf{e}_3 \mathbf{e}_3^\top = \mathbf{Q} \mathbf{\Lambda} \mathbf{Q}^\top
$$

where $\mathbf{Q} = [\mathbf{e}_1 \; \mathbf{e}_2 \; \mathbf{e}_3]$ and $\mathbf{\Lambda} = \text{diag}(\lambda_1, \lambda_2, \lambda_3)$.

### Tensor Geometry

The eigenvalues define the shape of the **diffusion ellipsoid** --- the surface of constant mean squared displacement:

$$
\mathbf{x}^\top \mathbf{D}^{-1} \mathbf{x} = 6t
$$

where $t$ is the diffusion time. The semi-axes of this ellipsoid have lengths $\sqrt{6t\lambda_i}$ along directions $\mathbf{e}_i$.

- **Prolate** ($\lambda_1 \gg \lambda_2 \approx \lambda_3$): Cigar-shaped, coherent single-fiber bundle
- **Oblate** ($\lambda_1 \approx \lambda_2 \gg \lambda_3$): Pancake-shaped, crossing fibers or planar structure
- **Spherical** ($\lambda_1 \approx \lambda_2 \approx \lambda_3$): Isotropic diffusion, CSF or gray matter

### Westin Shape Measures

Westin et al. (1997) defined three shape measures that sum to 1, providing a geometric decomposition of the tensor:

$$
\begin{aligned}
C_L &= \frac{\lambda_1 - \lambda_2}{\lambda_1} \quad &\text{(linear anisotropy --- fiber-like)} \\
C_P &= \frac{\lambda_2 - \lambda_3}{\lambda_1} \quad &\text{(planar anisotropy --- crossing/sheet-like)} \\
C_S &= \frac{\lambda_3}{\lambda_1} \quad &\text{(spherical isotropy --- ball-like)}
\end{aligned}
$$

with $C_L + C_P + C_S = 1$. In a coherent white matter voxel, $C_L \gg C_P, C_S$. At a fiber crossing, $C_P$ becomes significant.

---

## 4. Diffusion Scalar Metrics

### Fractional Anisotropy (FA)

FA quantifies the degree of anisotropy on a normalized scale $[0, 1]$:

$$
\text{FA} = \sqrt{\frac{1}{2}} \cdot \frac{\sqrt{(\lambda_1 - \lambda_2)^2 + (\lambda_2 - \lambda_3)^2 + (\lambda_3 - \lambda_1)^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}
$$

**Equivalence with the variance formulation**: Defining the mean eigenvalue $\bar{\lambda} = (\lambda_1 + \lambda_2 + \lambda_3)/3$, an equivalent expression is:

$$
\text{FA} = \sqrt{\frac{3}{2}} \cdot \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}
$$

**Proof of equivalence**: Expanding both forms shows they differ only by algebraic rearrangement. Note that:

$$
(\lambda_1 - \lambda_2)^2 + (\lambda_2 - \lambda_3)^2 + (\lambda_3 - \lambda_1)^2 = 3(\lambda_1^2 + \lambda_2^2 + \lambda_3^2) - (\lambda_1 + \lambda_2 + \lambda_3)^2
$$

and:

$$
3\left[(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2\right] = 3(\lambda_1^2 + \lambda_2^2 + \lambda_3^2) - (\lambda_1 + \lambda_2 + \lambda_3)^2
$$

which confirms $\frac{1}{2} \times (\text{pairwise form}) = \frac{3}{2} \times (\text{variance form})$.

**In terms of tensor invariants**: FA can also be written using the tensor invariants, making it manifestly rotation-invariant:

$$
\text{FA} = \sqrt{\frac{3}{2}} \cdot \sqrt{1 - \frac{I_2}{I_1^2 / 3}}
$$

where $I_1 = \text{tr}(\mathbf{D})$ and $I_2$ is as defined in Section 3.

**Implementation** (`nim_fa.m:23-24`):
```matlab
denom = sqrt(l(1)^2 + l(2)^2 + l(3)^2);
if denom > 0
    nim.FA(x,y,z) = sqrt(1/2) * sqrt((l(1)-l(2))^2 + (l(2)-l(3))^2 + (l(3)-l(1))^2) / denom;
end
```

Division by zero is protected: $\text{FA} = 0$ when all eigenvalues are zero.

### Mean Diffusivity (MD)

$$
\text{MD} = \frac{\lambda_1 + \lambda_2 + \lambda_3}{3} = \frac{\text{tr}(\mathbf{D})}{3}
$$

MD equals one-third of the trace, which is rotationally invariant. Typical values: $\sim 0.7$--$1.0 \times 10^{-3}$ mm$^2$/s in healthy white matter, $\sim 3.0 \times 10^{-3}$ mm$^2$/s in CSF.

### Axial Diffusivity (AD) and Radial Diffusivity (RD)

$$
\text{AD} = \lambda_1, \qquad \text{RD} = \frac{\lambda_2 + \lambda_3}{2}
$$

FA can be expressed in terms of AD and RD in the special case $\lambda_2 = \lambda_3 = \text{RD}$:

$$
\text{FA} = \sqrt{\frac{1}{2}} \cdot \frac{\sqrt{3} \, |\text{AD} - \text{RD}|}{\sqrt{\text{AD}^2 + 2 \, \text{RD}^2}}
$$

---

## 5. FACT Tractography Algorithm

The FACT (Fiber Assignment by Continuous Tracking) algorithm, implemented in `nim_tractography_standard.m`, uses discrete voxel directions with exact boundary intersection.

### Streamline ODE Formulation

Deterministic tractography formulates streamline tracing as an initial value problem (IVP) for an ordinary differential equation:

$$
\frac{d\mathbf{r}}{ds} = \mathbf{e}_1(\mathbf{r}(s)), \qquad \mathbf{r}(0) = \mathbf{r}_0
$$

where:

- $\mathbf{r}(s) = [x(s), y(s), z(s)]$ is the streamline position parameterized by arc length $s$
- $\mathbf{e}_1(\mathbf{r})$ is the primary eigenvector field (a unit vector field defined at every voxel)
- $\mathbf{r}_0$ is the seed point

The FACT algorithm solves this ODE using a piecewise-constant approximation of $\mathbf{e}_1$ within each voxel.

### Algorithm Overview

```
Input: Seed point r_0, eigenvector field e_1(x,y,z), FA field
Output: Track T = [r_0, r_1, r_2, ...]

1. Initialize: r = r_0, d_prev = e_1(voxel(r_0))
2. While termination criteria not met:
   a. Get current voxel indices: (i, j, k) = floor(r)
   b. Get direction: d = e_1(i, j, k)  [discrete, no interpolation]
   c. Ensure consistency: if dot(d, d_prev) < 0, flip d = -d
   d. Compute boundary exit: r_exit = find_boundary_exit(r, d)
   e. Step to boundary: r = r_exit + epsilon * d  [enter next voxel]
   f. Check termination: FA(new voxel) < threshold? angle > limit?
   g. Store point: append r to T
   h. Update: d_prev = d
3. Repeat tracking in opposite direction from r_0
4. Concatenate: T = [reverse(T_backward), T_forward]
```

### Voxel Boundary Intersection

Given position $\mathbf{p} = (p_x, p_y, p_z)$ and direction $\mathbf{d} = (d_x, d_y, d_z)$ within a voxel with boundaries $[x_\min, x_\max] \times [y_\min, y_\max] \times [z_\min, z_\max]$, the ray-box intersection is computed parametrically.

For each axis $\alpha \in \{x, y, z\}$, compute the parametric distances to the boundaries:

$$
t_\alpha^{\text{near}} = \frac{\alpha_\min - p_\alpha}{d_\alpha}, \qquad t_\alpha^{\text{far}} = \frac{\alpha_\max - p_\alpha}{d_\alpha}
$$

(swapping near/far if $d_\alpha < 0$). The exit parameter is:

$$
t_{\text{exit}} = \min\!\left(t_x^{\text{far}}, \, t_y^{\text{far}}, \, t_z^{\text{far}}\right), \quad t > 0
$$

The exit point is:

$$
\mathbf{p}_{\text{exit}} = \mathbf{p} + t_{\text{exit}} \cdot \mathbf{d}
$$

### Direction Consistency

Eigenvectors are bidirectional ($\mathbf{e}_1$ and $-\mathbf{e}_1$ are both valid). The dot product check ensures consistent propagation:

$$
\text{if } \mathbf{d}_{\text{new}} \cdot \mathbf{d}_{\text{prev}} < 0 \implies \mathbf{d}_{\text{new}} \leftarrow -\mathbf{d}_{\text{new}}
$$

### Angular Stopping Criterion

The angle $\theta$ between consecutive directions is checked against a threshold $\theta_\max$:

$$
\cos\theta = \frac{|\mathbf{d}_{\text{new}} \cdot \mathbf{d}_{\text{prev}}|}{|\mathbf{d}_{\text{new}}| \, |\mathbf{d}_{\text{prev}}|}
$$

If $\theta > \theta_\max$ (typically 25--60 degrees), the streamline is terminated. Pre-computing $\cos\theta_\max$ avoids the expensive $\arccos$ operation at every step.

### Bidirectional Tracking

Each seed point generates two half-tracks. The backward track starts by negating the initial direction. The final track is:

$$
T = [\text{reverse}(T_{\text{backward}}); \; T_{\text{forward}}]
$$

### Track Length

The total arc length of a track with $N$ points is:

$$
L = \sum_{i=1}^{N-1} \|\mathbf{r}_{i+1} - \mathbf{r}_i\|_2
$$

Tracks shorter than a minimum length (typically 10--30 mm) are filtered out as likely spurious.

---

## 6. Numerical Integration Methods

HINEC's high-order tractography (`nim_tractography_hinec.m`) replaces the piecewise-constant FACT approach with numerical ODE integration on a continuous interpolated direction field:

$$
\frac{d\mathbf{y}}{ds} = \mathbf{f}(\mathbf{y}), \qquad \mathbf{y} = [x, y, z], \quad s = \text{arc length}, \quad \mathbf{f} = \text{interpolated direction field}
$$

### Euler Method (Order 1, `integration_order = 1`)

The simplest explicit method:

$$
\mathbf{y}_{n+1} = \mathbf{y}_n + h \, \mathbf{f}(\mathbf{y}_n)
$$

- Local truncation error: $O(h^2)$
- Global error: $O(h)$ --- first-order accuracy
- One function evaluation per step
- Fast but least accurate; equivalent to standard step-based tracking

### RK2 Midpoint Method (Order 2, `integration_order = 2`)

Uses a trial step to the midpoint for better accuracy:

$$
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(\mathbf{y}_n + \frac{h}{2} \mathbf{k}_1\right) \\
\mathbf{y}_{n+1} &= \mathbf{y}_n + h \, \mathbf{k}_2
\end{aligned}
$$

- Local truncation error: $O(h^3)$
- Global error: $O(h^2)$ --- second-order accuracy
- Two function evaluations per step

### RK4 Classical Method (Order 4, `integration_order = 4`)

The standard fourth-order Runge-Kutta method:

$$
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(\mathbf{y}_n) \\
\mathbf{k}_2 &= \mathbf{f}\!\left(\mathbf{y}_n + \frac{h}{2} \mathbf{k}_1\right) \\
\mathbf{k}_3 &= \mathbf{f}\!\left(\mathbf{y}_n + \frac{h}{2} \mathbf{k}_2\right) \\
\mathbf{k}_4 &= \mathbf{f}(\mathbf{y}_n + h \, \mathbf{k}_3) \\
\mathbf{y}_{n+1} &= \mathbf{y}_n + \frac{h}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)
\end{aligned}
$$

- Local truncation error: $O(h^5)$
- Global error: $O(h^4)$ --- fourth-order accuracy
- Four function evaluations per step
- Excellent accuracy-to-cost ratio; the default HINEC method

### RKF45 Adaptive Method (Order 4-5, `integration_order = 5`)

The Dormand-Prince embedded pair provides automatic step size control by computing both a 4th-order and 5th-order solution using shared function evaluations.

**Stage evaluations** (7 stages, $i = 1, \ldots, 7$):

$$
\mathbf{k}_i = \mathbf{f}\!\left(\mathbf{y}_n + h \sum_{j=1}^{i-1} a_{ij} \mathbf{k}_j\right)
$$

where the coefficients $a_{ij}$ are given by the Dormand-Prince Butcher tableau.

**Two solutions from shared stages**:

$$
\begin{aligned}
\hat{\mathbf{y}}_{n+1} &= \mathbf{y}_n + h \sum_{i=1}^{7} \hat{b}_i \mathbf{k}_i \quad &\text{(5th order)} \\
\tilde{\mathbf{y}}_{n+1} &= \mathbf{y}_n + h \sum_{i=1}^{7} b_i \mathbf{k}_i \quad &\text{(4th order)}
\end{aligned}
$$

**Error estimate**:

$$
\text{err} = \|\hat{\mathbf{y}}_{n+1} - \tilde{\mathbf{y}}_{n+1}\| = h \left\| \sum_{i=1}^{7} (\hat{b}_i - b_i) \mathbf{k}_i \right\|
$$

**Step size control**:

$$
h_{\text{new}} = \text{safety} \cdot h \cdot \left(\frac{\text{tol}}{\text{err}}\right)^{1/5}
$$

$$
h_{\text{new}} = \max(h_\min, \min(h_\max, h_{\text{new}}))
$$

If $\text{err} \leq \text{tol}$: accept step, use $\hat{\mathbf{y}}_{n+1}$, adjust $h$ for next step.
If $\text{err} > \text{tol}$: reject step, reduce $h$, retry.

**HINEC defaults**: $\text{tol} = 0.01$ voxels, $\text{safety} = 0.9$, $h_\min = 0.01$, $h_\max = 1.0$.

For the complete Dormand-Prince coefficients and derivation, see [RKF.md](RKF.md). For usage examples and parameter guidelines, see [RKF_Usage.md](RKF_Usage.md).

### Error Analysis Summary

The **local truncation error** (LTE) for an order-$p$ method scales as:

$$
\text{LTE} = C \cdot h^{p+1}
$$

where $C$ depends on the $(p+1)$-th derivative of the solution. The **global error** after $N = L/h$ steps over total arc length $L$ accumulates to:

$$
\text{Global error} \sim N \cdot \text{LTE} = \frac{L}{h} \cdot C \cdot h^{p+1} = C L \cdot h^p
$$

This explains why the global accuracy order equals $p$, one less than the local order.

### Comparison of Methods

| Method | Order | Evaluations/Step | Relative Cost | Relative Accuracy | Best For |
|---|---|---|---|---|---|
| Euler | 1 | 1 | $1\times$ | Low | Quick exploration |
| RK2 | 2 | 2 | $2\times$ | Moderate | Fast, reasonable accuracy |
| RK4 | 4 | 4 | $4\times$ | High | Default, good cost/accuracy |
| RKF45 | 4-5 | 6-7 | $6$--$7\times$ | Highest (adaptive) | Publication quality |

Higher-order methods are especially valuable in regions of high fiber curvature, where lower-order methods accumulate significant tracking error.

---

## 7. Interpolation Methods

### Trilinear Interpolation

The default interpolation method in HINEC. For a query point $\mathbf{p} = (x, y, z)$ inside a voxel, the interpolated value is a weighted average of the 8 surrounding voxel centers.

Given voxel centers at integer coordinates with values $V(i,j,k)$, and fractional coordinates:

$$
x_d = x - \lfloor x \rfloor, \quad y_d = y - \lfloor y \rfloor, \quad z_d = z - \lfloor z \rfloor
$$

The trilinear interpolation can be written as a tensor product of 1D linear interpolations:

$$
V(x,y,z) = \sum_{\alpha=0}^{1} \sum_{\beta=0}^{1} \sum_{\gamma=0}^{1} V_{\alpha\beta\gamma} \, w_\alpha(x_d) \, w_\beta(y_d) \, w_\gamma(z_d)
$$

where $w_0(t) = 1 - t$ and $w_1(t) = t$. Expanding:

$$
\begin{aligned}
V(x,y,z) =\; & V_{000}(1-x_d)(1-y_d)(1-z_d) + V_{100} \, x_d(1-y_d)(1-z_d) \\
          +\; & V_{010}(1-x_d) \, y_d(1-z_d) + V_{110} \, x_d \, y_d(1-z_d) \\
          +\; & V_{001}(1-x_d)(1-y_d) \, z_d + V_{101} \, x_d(1-y_d) \, z_d \\
          +\; & V_{011}(1-x_d) \, y_d \, z_d + V_{111} \, x_d \, y_d \, z_d
\end{aligned}
$$

**Implementation**: HINEC uses MATLAB's `griddedInterpolant` with `'linear'` method for 2--5$\times$ faster interpolation compared to `interp3` (`nim_tractography_hinec.m:239-242`):

```matlab
nim.FA_interp = griddedInterpolant(grid_vectors, nim.FA, 'linear', 'none');
nim.v1_x_interp = griddedInterpolant(grid_vectors, nim.v1_x, 'linear', 'none');
nim.v1_y_interp = griddedInterpolant(grid_vectors, nim.v1_y, 'linear', 'none');
nim.v1_z_interp = griddedInterpolant(grid_vectors, nim.v1_z, 'linear', 'none');
```

Each eigenvector component is interpolated independently, and the result is renormalized to unit length.

### Cubic Spline Interpolation

Available via `interp_method = 'cubic'`. Uses cubic B-splines for $C^2$-continuous interpolation, providing smoother direction fields at the cost of additional computation. The cubic B-spline basis function is:

$$
B_3(t) = \frac{1}{6} \begin{cases}
(2-|t|)^3 & 1 \leq |t| < 2 \\
4 - 6t^2 + 3|t|^3 & |t| < 1 \\
0 & |t| \geq 2
\end{cases}
$$

Uses the same `griddedInterpolant` framework with `'cubic'` method.

### High-Order Spectral Interpolation

The `nim_interp.m` utility implements spectral element interpolation using **Gauss-Lobatto-Legendre** (GLL) quadrature nodes:

1. Expand eigenvectors from cell centers to vertices
2. Place $(p+1)^3$ interpolation nodes per voxel using GLL quadrature
3. Perform 3D spline interpolation within each voxel

GLL nodes are the roots of:

$$
(1 - \xi^2) P_N'(\xi) = 0
$$

where $P_N'$ is the derivative of the $N$-th Legendre polynomial. These nodes include the endpoints $\xi = \pm 1$ and cluster near element boundaries, providing optimal interpolation accuracy with minimal Runge phenomenon.

The Lagrange interpolant through GLL nodes $\xi_0, \ldots, \xi_N$ is:

$$
f(\xi) \approx \sum_{j=0}^{N} f(\xi_j) \, \ell_j(\xi), \qquad \ell_j(\xi) = \prod_{\substack{k=0 \\ k \neq j}}^{N} \frac{\xi - \xi_k}{\xi_j - \xi_k}
$$

Uniform nodes from `zwuni.m` provide an alternative with equally-spaced points.

### Direction Field Normalization

After interpolation, the direction vector must be renormalized because linear interpolation of unit vectors does not produce unit vectors:

$$
\mathbf{d}_{\text{interp}} = \text{interpolate}(\mathbf{e}_1, \mathbf{r}), \qquad \hat{\mathbf{d}} = \frac{\mathbf{d}_{\text{interp}}}{\|\mathbf{d}_{\text{interp}}\|}
$$

The interpolation error in the direction (angular deviation) is bounded by:

$$
\theta_{\text{error}} \leq \arcsin\!\left(\frac{\|\nabla \mathbf{e}_1\| \, \Delta x}{2}\right)
$$

where $\Delta x$ is the voxel size. This is why trilinear interpolation works well: for smoothly varying fiber directions, the angular error per voxel is small.

---

## 8. Registration Mathematics

### Affine Transformation

An affine transformation maps coordinates between spaces using a $4 \times 4$ matrix (12 degrees of freedom):

$$
\begin{bmatrix} x' \\ y' \\ z' \\ 1 \end{bmatrix} = \begin{bmatrix} a_{11} & a_{12} & a_{13} & t_1 \\ a_{21} & a_{22} & a_{23} & t_2 \\ a_{31} & a_{32} & a_{33} & t_3 \\ 0 & 0 & 0 & 1 \end{bmatrix} \begin{bmatrix} x \\ y \\ z \\ 1 \end{bmatrix}
$$

The $3 \times 3$ submatrix decomposes into rotation (3 DOF), scaling (3 DOF), and shear (3 DOF). The translation vector $\mathbf{t} = [t_1, t_2, t_3]^\top$ adds 3 DOF.

A rigid-body transformation (6 DOF) restricts to rotation + translation only:

$$
\mathbf{r}' = \mathbf{R} \mathbf{r} + \mathbf{t}, \qquad \mathbf{R} \in SO(3), \quad \det(\mathbf{R}) = 1
$$

### Registration Cost Functions

HINEC supports two registration backends:

**Correlation Ratio** (FSL FLIRT, same-modality):

$$
\text{CR} = 1 - \frac{\text{Var}(Y | X)}{\text{Var}(Y)} = 1 - \frac{\sum_k n_k \, \sigma_k^2}{\sigma_Y^2}
$$

where $Y$ is the target image, $X$ is the source partitioned into $k$ bins, $n_k$ is the count in bin $k$, $\sigma_k^2$ is the variance of $Y$ within bin $k$, and $\sigma_Y^2$ is the total variance.

**Normalized Mutual Information** (quality assessment in `registration_utils.m`):

$$
\text{NMI} = \frac{H(A) + H(B)}{H(A, B)}
$$

where the entropies are:

$$
H(A) = -\sum_a p(a) \log p(a), \quad H(A, B) = -\sum_{a,b} p(a,b) \log p(a,b)
$$

$p(a)$ and $p(a,b)$ are estimated from 64-bin histograms. NMI $= 1$ indicates no shared information; NMI $> 1$ indicates increasing alignment.

### Multi-Modal Registration Chain

HINEC constructs a chain of transforms between three coordinate spaces:

$$
\text{DWI native} \xleftrightarrow{\text{6 DOF FLIRT}} \text{T1 native} \xleftrightarrow{\text{12 DOF FLIRT + FNIRT}} \text{MNI standard}
$$

**Forward chain** (DWI $\to$ MNI): compose transforms

$$
\mathbf{T}_{\text{total}} = \mathbf{T}_{\text{T1} \to \text{MNI}} \cdot \mathbf{T}_{\text{DWI} \to \text{T1}}
$$

**Inverse chain** (MNI $\to$ DWI): apply inverse transforms in reverse order

$$
\mathbf{T}_{\text{total}}^{-1} = \mathbf{T}_{\text{DWI} \to \text{T1}}^{-1} \cdot \mathbf{T}_{\text{T1} \to \text{MNI}}^{-1}
$$

This chain is used by `nim_apply_transforms.m` to bring atlas labels from MNI space to DWI space.

### Nearest-Neighbor Resampling for Labels

When transforming parcellation labels, nearest-neighbor interpolation is required to preserve discrete integer values:

$$
\text{label}(\mathbf{r}') = \text{label}(\text{round}(\mathbf{T}^{-1} \mathbf{r}'))
$$

Linear or spline interpolation would produce non-integer values that corrupt label identity.

---

## 9. Tissue Segmentation Thresholds

FA-based tissue classification (`preproc_tissue_segmentation.m`) uses empirical thresholds:

$$
\text{WM}: \text{FA}(\mathbf{r}) > 0.2, \qquad \text{GM}: 0.05 < \text{FA}(\mathbf{r}) \leq 0.2, \qquad \text{CSF}: \text{FA}(\mathbf{r}) \leq 0.05
$$

The WM mask is eroded (morphological erosion with spherical kernel, radius=1) to remove boundary voxels where partial volume effects contaminate the mask:

$$
\text{WM}_{\text{eroded}} = \text{WM} \ominus \mathcal{B}_1
$$

where $\mathcal{B}_1$ is a unit ball structuring element and $\ominus$ denotes morphological erosion.

Masks are mutually exclusive and constrained by the brain mask:

$$
\text{WM} \cap \text{GM} = \emptyset, \quad \text{WM} \cap \text{CSF} = \emptyset, \quad \text{GM} \cap \text{CSF} = \emptyset
$$

$$
\text{WM} \cup \text{GM} \cup \text{CSF} \subseteq \mathcal{M}_{\text{brain}}
$$

Expected tissue proportions: WM $\sim$40--45%, GM $\sim$40--45%, CSF $\sim$10--20% of brain volume.

---

## 10. Connectivity Matrix Construction

Given a set of tracks $\mathcal{T} = \{t_1, t_2, \ldots, t_M\}$ and a parcellation with $N$ regions (`nim_plot_connectivity_matrix.m`):

$$
C_{ij} = \#\{t_k \in \mathcal{T} : \text{start}(t_k) \in R_i \text{ and } \text{end}(t_k) \in R_j\}
$$

where $R_i$ denotes parcellation region $i$, and $\text{start}(t_k)$, $\text{end}(t_k)$ are the first and last points of track $t_k$.

For a symmetric (undirected) connectivity matrix:

$$
C_{ij}^{\text{sym}} = C_{ij} + C_{ji}
$$

Optional normalization:

$$
\hat{C}_{ij} = \frac{C_{ij}}{\max_{p,q} C_{pq}}
$$

The resulting matrix $\mathbf{C} \in \mathbb{R}^{N \times N}$ represents the structural connectivity between brain regions, where $C_{ij}$ is the number of streamlines connecting region $i$ to region $j$.

---

## Cross-References

- RKF45 complete derivation and coefficients: [RKF.md](RKF.md)
- RKF45 usage guide and parameter selection: [RKF_Usage.md](RKF_Usage.md)
- High-order methods implementation details: [High_Order.md](High_Order.md)
- Scientific context for these methods: [SCIENCE.md](SCIENCE.md)
- Architecture and data flow: [ARCHITECTURE.md](ARCHITECTURE.md)
- Complete function reference: [API_REFERENCE.md](API_REFERENCE.md)
