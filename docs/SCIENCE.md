# The Neuroscience of Diffusion MRI and Tractography

This document explains the scientific foundations behind HINEC for researchers and students who are new to diffusion MRI or need a refresher on the biology and physics that motivate every computation in the pipeline.

---

## 1. Introduction to Brain White Matter

### Neurons and Axons

The brain contains approximately 86 billion neurons. Each neuron communicates with others through electrical signals transmitted along axons --- long, thin projections that can extend centimeters across the brain. While neuronal cell bodies (gray matter) perform computation, axons bundled together form white matter tracts that serve as the brain's communication highways.

White matter appears white in dissected brain tissue because axons are wrapped in myelin, a fatty insulating sheath produced by oligodendrocytes. Myelin dramatically increases the speed of electrical signal conduction (from $\sim$1 m/s unmyelinated to $\sim$100 m/s myelinated).

### White Matter Tracts and Connectivity

Axons traveling between the same brain regions tend to bundle together into coherent fiber tracts. Major examples include:

- **Corpus callosum**: Connects left and right hemispheres ($\sim$200 million axons)
- **Corticospinal tract**: Motor commands from cortex to spinal cord
- **Arcuate fasciculus**: Links language areas (Broca's and Wernicke's)
- **Cingulum**: Connects limbic structures for emotion and memory
- **Uncinate fasciculus**: Connects temporal and frontal lobes

The pattern of these connections --- the brain's structural connectome --- fundamentally shapes how information flows and is processed.

### Why Mapping Connectivity Matters

**Clinical applications**:

- Pre-surgical planning: Identifying critical tracts to avoid during tumor resection
- Traumatic brain injury: Detecting diffuse axonal injury invisible on standard MRI
- Neurodegenerative disease: Tracking white matter degeneration in multiple sclerosis, ALS, and Alzheimer's
- Stroke recovery: Understanding which connections survive and how they reorganize

**Research applications**:

- Mapping the human structural connectome
- Understanding brain development and aging
- Relating structural connectivity to cognitive function
- Investigating psychiatric disorders (schizophrenia, depression, autism)

---

## 2. Physics of Diffusion MRI

### Brownian Motion and the Einstein Relation

Water molecules in tissue undergo random thermal motion (Brownian motion). Einstein (1905) showed that for free (unrestricted) diffusion in three dimensions, the mean squared displacement grows linearly with time:

$$
\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle = 6Dt
$$

where $D$ is the scalar diffusion coefficient (units: mm$^2$/s) and $t$ is the diffusion time. For pure water at body temperature ($37°$C), $D \approx 3.0 \times 10^{-3}$ mm$^2$/s.

At a typical diffusion time of $\Delta \approx 40$ ms, the root-mean-square displacement is:

$$
\sqrt{\langle r^2 \rangle} = \sqrt{6 D \Delta} \approx \sqrt{6 \times 3.0 \times 10^{-3} \times 0.04} \approx 0.027 \text{ mm} = 27 \; \mu\text{m}
$$

This length scale is comparable to the diameter of axons ($1$--$20 \; \mu\text{m}$) and the spacing between them, which is exactly why diffusion MRI is sensitive to microstructure.

### Anisotropic Diffusion in Tissue

In biological structures, molecular motion is constrained. In white matter, the cylindrical geometry of axons restricts diffusion perpendicular to the fiber direction while allowing relatively free diffusion along the fiber axis. For anisotropic diffusion, the scalar $D$ is replaced by the diffusion tensor $\mathbf{D}$, and the mean squared displacement along any direction $\hat{\mathbf{n}}$ is:

$$
\langle (\hat{\mathbf{n}} \cdot \Delta\mathbf{r})^2 \rangle = 2 \, \hat{\mathbf{n}}^\top \mathbf{D} \hat{\mathbf{n}} \, t
$$

Along the fiber direction ($\hat{\mathbf{n}} = \mathbf{e}_1$), diffusion is fast ($\lambda_1$ large). Perpendicular to fibers, diffusion is restricted ($\lambda_2, \lambda_3$ small). This directional dependence is the fundamental measurement that enables fiber tracking.

### How MRI Measures Diffusion

Standard MRI creates images based on water proton density and relaxation times. Diffusion MRI adds magnetic field gradients that make the signal sensitive to molecular displacement.

The Stejskal-Tanner pulsed gradient spin-echo sequence works as follows:

1. A 90-degree RF pulse excites protons, aligning their magnetic moments
2. A magnetic field gradient pulse is applied in direction $\mathbf{g}$ with strength $G$, causing protons at different positions to precess at different frequencies --- encoding their positions in their phase
3. After a diffusion time $\Delta$, a second gradient pulse of equal magnitude is applied
4. Stationary molecules experience equal and opposite phase shifts that cancel perfectly, preserving signal
5. Molecules that moved along the gradient direction accumulate a net phase $\phi$:

$$
\phi = \gamma \, \mathbf{g} \cdot [\mathbf{r}(t_2) - \mathbf{r}(t_1)] \cdot G \delta
$$

The ensemble dephasing from many molecules causes signal attenuation. The greater the molecular motion along the gradient direction, the greater the signal loss.

### B-values and Gradient Directions

Two key acquisition parameters control diffusion sensitivity:

**B-value** (units: s/mm$^2$): Controls the strength of diffusion weighting. It combines gradient amplitude, duration, and timing:

$$
b = \gamma^2 G^2 \delta^2 \left(\Delta - \frac{\delta}{3}\right)
$$

- $b = 0$: No diffusion weighting (reference image, "b0")
- $b = 1000$: Standard clinical DTI
- $b = 2000$--$3000$: Higher angular resolution, better for crossing fibers

**Signal-to-noise trade-off**: Higher $b$-values provide stronger diffusion contrast but lower SNR because the signal decays exponentially:

$$
\text{SNR}(b) \propto \frac{S_0 \exp(-b \cdot \text{ADC})}{\sigma_{\text{noise}}}
$$

**Gradient directions** (b-vectors): Each diffusion-weighted volume applies the gradient in a specific 3D direction. A minimum of 6 non-collinear directions is needed to fit the diffusion tensor (6 unknowns), but 30--60+ directions are typical for quality DTI, providing better angular coverage and noise averaging.

In the HINEC pipeline, b-values are stored in `.bval` files and gradient directions in `.bvec` files, read by `nim_read.m`. Volumes with $b \leq 50$ s/mm$^2$ are classified as b0 references.

---

## 3. Diffusion Tensor Imaging (DTI)

### The Diffusion Tensor Model

DTI models the diffusion process at each voxel using a $3 \times 3$ symmetric positive-definite matrix called the diffusion tensor $\mathbf{D}$:

$$
\mathbf{D} = \begin{bmatrix}
D_{xx} & D_{xy} & D_{xz} \\
D_{xy} & D_{yy} & D_{yz} \\
D_{xz} & D_{yz} & D_{zz}
\end{bmatrix}
$$

Being symmetric, $\mathbf{D}$ has 6 unique parameters. This is why a minimum of 6 gradient directions plus one b0 image are required for DTI.

### The Signal Equation

The relationship between the measured diffusion-weighted signal and the diffusion tensor is the Stejskal-Tanner equation:

$$
S(\mathbf{g}, b) = S_0 \exp\!\left(-b \, \mathbf{g}^\top \mathbf{D} \mathbf{g}\right)
$$

where:

- $S(\mathbf{g}, b)$ = signal with gradient direction $\mathbf{g}$ and b-value $b$
- $S_0$ = signal without diffusion weighting (b0 image)
- $\mathbf{g} = [g_x, g_y, g_z]^\top$ = unit gradient direction vector
- $\mathbf{g}^\top \mathbf{D} \mathbf{g}$ = apparent diffusion coefficient along direction $\mathbf{g}$

Taking the logarithm yields a linear equation solvable by least squares:

$$
\ln\!\left(\frac{S}{S_0}\right) = -b \, \mathbf{g}^\top \mathbf{D} \mathbf{g}
$$

In HINEC, `nim_dt.m` implements direct least squares fitting, while `nim_dt_spd.m` adds BFGS optimization to enforce positive-definiteness. See [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md) for the complete derivation.

### Eigendecomposition

The diffusion tensor $\mathbf{D}$ is a real symmetric matrix, so it has three real eigenvalues ($\lambda_1 \geq \lambda_2 \geq \lambda_3 \geq 0$) and three orthogonal eigenvectors ($\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3$):

$$
\mathbf{D} = \lambda_1 \mathbf{e}_1 \mathbf{e}_1^\top + \lambda_2 \mathbf{e}_2 \mathbf{e}_2^\top + \lambda_3 \mathbf{e}_3 \mathbf{e}_3^\top
$$

**Physical interpretation**:

- $\lambda_1$ (largest eigenvalue): Magnitude of diffusion along the primary direction. In coherent white matter, this is the **axial diffusivity** --- diffusion parallel to fibers.
- $\mathbf{e}_1$ (primary eigenvector): Points along the direction of maximum diffusion. In white matter, this estimates the **local fiber direction**.
- $\lambda_2, \lambda_3$: Diffusion magnitudes in the two perpendicular directions. Their average $(\lambda_2 + \lambda_3)/2$ is the **radial diffusivity**.

The tensor can be visualized as an ellipsoid with semi-axes of length $\sqrt{\lambda_i}$ along directions $\mathbf{e}_i$. In white matter (anisotropic diffusion), the ellipsoid is elongated along the fiber direction. In CSF (isotropic diffusion), it approaches a sphere.

In HINEC, eigendecomposition is computed by `nim_eig.m`, which sorts eigenvalues in descending order.

### Limitations of the Single-Tensor Model

DTI assumes a single fiber population per voxel. This model breaks down in regions with:

- **Crossing fibers**: Two or more tracts intersecting (estimated in 60--90% of white matter voxels)
- **Kissing fibers**: Tracts running parallel then diverging
- **Fanning fibers**: A single tract splaying into multiple branches

In these regions, the tensor becomes oblate (pancake-shaped) rather than prolate (cigar-shaped). The primary eigenvector no longer reliably indicates any single fiber direction --- instead it points between the true fiber directions. Mathematically, for two crossing fibers with directions $\mathbf{v}_1$ and $\mathbf{v}_2$, the observed tensor is approximately:

$$
\mathbf{D}_{\text{observed}} \approx f_1 \mathbf{D}_1 + f_2 \mathbf{D}_2
$$

where $f_1 + f_2 = 1$ are volume fractions. This is a fundamental limitation that affects tractography accuracy.

---

## 4. Scalar Metrics from DTI

### Fractional Anisotropy (FA)

FA is the most widely used DTI metric. It quantifies how directionally dependent the diffusion is, on a scale from 0 (perfectly isotropic) to 1 (perfectly anisotropic):

$$
\text{FA} = \sqrt{\frac{1}{2}} \cdot \frac{\sqrt{(\lambda_1 - \lambda_2)^2 + (\lambda_2 - \lambda_3)^2 + (\lambda_3 - \lambda_1)^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}
$$

**Interpretation by tissue type**:

| FA Range | Tissue Type | Interpretation |
|---|---|---|
| 0.0 -- 0.1 | CSF, background | Free isotropic diffusion or noise |
| 0.1 -- 0.25 | Gray matter | Low anisotropy, complex cellular architecture |
| 0.25 -- 0.5 | Mixed/crossing WM | Moderate anisotropy, possible fiber crossings |
| 0.5 -- 0.8 | Coherent white matter | High anisotropy, well-organized fiber bundles |
| 0.8 -- 1.0 | Highly coherent WM | Very high anisotropy (e.g., corpus callosum) |

FA is sensitive to white matter pathology: demyelination, axonal loss, edema, and inflammation all tend to reduce FA. However, FA is non-specific --- it cannot distinguish between different types of pathology.

In HINEC, FA is computed by `nim_fa.m` and serves as a key input for tractography (seed placement, stopping criteria) and tissue segmentation.

### Mean Diffusivity (MD)

The average diffusion across all directions:

$$
\text{MD} = \frac{\lambda_1 + \lambda_2 + \lambda_3}{3} = \frac{\text{tr}(\mathbf{D})}{3}
$$

MD is elevated in CSF (free water) and edema, and reduced in densely packed cellular tissue. Unlike FA, MD is rotationally invariant and provides a measure of overall diffusion magnitude.

### Axial Diffusivity (AD) and Radial Diffusivity (RD)

$$
\text{AD} = \lambda_1, \qquad \text{RD} = \frac{\lambda_2 + \lambda_3}{2}
$$

These directional components can help distinguish types of white matter damage:

- **Reduced AD**: May indicate axonal damage or degeneration
- **Elevated RD**: May indicate demyelination

However, these interpretations should be applied cautiously and validated in context.

---

## 5. Fiber Tractography

Tractography reconstructs 3D fiber pathways from the voxel-wise direction estimates provided by DTI. It is the core computational task in HINEC.

### Deterministic vs Probabilistic Approaches

**Deterministic tractography** (used in HINEC): Given a seed point, follows the local fiber direction step-by-step to trace a single streamline. The result is a specific, reproducible 3D curve for each seed point. Simple and fast, but does not capture uncertainty.

**Probabilistic tractography** (not implemented in HINEC): Samples from a distribution of possible fiber directions at each step, generating many streamlines per seed point. Produces a probability map of connectivity rather than specific tracts. Better handles uncertainty and crossing fibers, but computationally expensive.

### Streamline Tracking as an ODE

The fundamental mathematical formulation: streamline tracking solves an initial value problem for an ordinary differential equation:

$$
\frac{d\mathbf{r}}{ds} = \mathbf{e}_1(\mathbf{r}), \qquad \mathbf{r}(0) = \mathbf{r}_0
$$

Starting from seed point $\mathbf{r}_0$, follow the primary eigenvector $\mathbf{e}_1$ step by step. At each point along the streamline:

1. Determine the local fiber direction (primary eigenvector $\mathbf{e}_1$)
2. Advance a small step along that direction
3. At the new position, determine the updated fiber direction
4. Repeat until a stopping criterion is met

The result is a polyline (sequence of 3D points) that approximates the trajectory of a fiber tract through the brain.

### FACT Algorithm (Standard Tractography)

HINEC's standard tractography uses FACT (Fiber Assignment by Continuous Tracking), implemented in `nim_tractography_standard.m`.

FACT is a discrete voxel-by-voxel tracking method that uses a piecewise-constant direction field:

1. Start at a seed point within a voxel
2. Assign the voxel's primary eigenvector as the tracking direction
3. Compute the intersection of this direction with the voxel boundary (ray-box intersection)
4. Step to the boundary intersection point, entering the next voxel
5. In the new voxel, use that voxel's eigenvector (checking direction consistency via dot product)
6. Repeat until stopping criteria are met

The effective integration step for FACT is the Euler method with variable step size (determined by voxel boundary distance):

$$
\mathbf{r}_{n+1} = \mathbf{r}_n + t_{\text{exit}} \cdot \mathbf{e}_1(\text{voxel}_n)
$$

where $t_{\text{exit}}$ is the parametric distance to the voxel boundary. This makes FACT fast and straightforward but can produce angular trajectories at voxel boundaries.

### High-Order Methods (HINEC Tractography)

HINEC's high-order tractography (`nim_tractography_hinec.m`) improves upon FACT with two enhancements:

**Interpolation**: Instead of using discrete voxel directions, trilinear interpolation creates a continuous direction field. At any sub-voxel position, the local direction is a weighted average of the 8 surrounding voxel centers. This produces smoother, more accurate trajectories.

**Numerical integration**: Instead of simple stepping, numerical ODE integration methods propagate the streamline with controlled accuracy:

| Method | Order | Description | Global Error |
|---|---|---|---|
| Euler | 1 | Simple forward step | $O(h)$ |
| RK2 | 2 | Midpoint method | $O(h^2)$ |
| RK4 | 4 | Classical Runge-Kutta | $O(h^4)$ |
| RKF45 | 4-5 | Adaptive Dormand-Prince | Adaptive error control |

RKF45 uses an embedded pair of 4th and 5th order methods to estimate local error and automatically adjust the step size:

$$
h_{\text{new}} = 0.9 \cdot h \left(\frac{\text{tol}}{\text{err}}\right)^{1/5}
$$

using smaller steps in regions of high curvature and larger steps in straight regions. See [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md) for the complete numerical methods, and [RKF.md](RKF.md) for the full RKF45 derivation.

### Seeding Strategies

Where tractography begins dramatically affects the results:

**FA-based seeding**: Place seeds only where FA exceeds a threshold. This concentrates seeds in well-organized white matter but misses regions with crossing fibers or partial volume effects.

**Uniform grid seeding** (preferred): Place seeds on a regular 3D grid throughout the brain mask, regardless of local FA. This ensures uniform spatial coverage. The seed density parameter controls the number of seeds per voxel dimension (e.g., density=4 means $4 \times 4 \times 4 = 64$ seeds per voxel).

In HINEC, `runTractography.m` implements a 3-tier seeding priority system that prefers the brain mask, falls back to expanded parcellation, and uses FA-threshold as last resort.

### Stopping Criteria

Tractography terminates a streamline when any condition is met:

1. **FA threshold**: $\text{FA}(\mathbf{r}) < \text{FA}_\min$ (typically 0.1--0.2), indicating the streamline has left organized white matter
2. **Angle threshold**: The angle between consecutive steps exceeds $\theta_\max$ (typically 25--60 degrees): $\arccos(|\mathbf{d}_n \cdot \mathbf{d}_{n-1}|) > \theta_\max$
3. **Brain mask boundary**: The streamline exits the brain
4. **Maximum steps**: Safety limit to prevent infinite loops

---

## 6. Anatomically Constrained Tractography (ACT)

Standard stopping criteria based solely on FA and angle can produce anatomically implausible results --- streamlines that terminate in white matter (biologically, fibers should end in gray matter) or enter CSF (where no fibers exist).

ACT, implemented in `nim_tractography_hinec.m`, uses tissue segmentation to enforce anatomical plausibility:

### Tissue Classification

HINEC segments brain tissue using FA values (computed in `preproc_tissue_segmentation.m`):

| Tissue | FA Criterion | Tractography Role |
|---|---|---|
| White matter (WM) | $\text{FA} > 0.2$ | Continue tracking |
| Gray matter (GM) | $0.05 < \text{FA} \leq 0.2$ | Valid termination |
| CSF | $\text{FA} \leq 0.05$ | Invalid --- reject track |

### ACT Rules

1. **Continue tracking** while the streamline remains in white matter
2. **Accept termination** if the streamline enters gray matter (biologically valid --- axons synapse in gray matter)
3. **Reject the streamline** if it enters CSF (biologically invalid --- no fibers in fluid)
4. Standard FA and angle criteria still apply as additional constraints

ACT significantly improves the anatomical plausibility of tractography results by reducing false positive connections through non-neural tissue.

---

## 7. Brain Parcellation and Connectivity

### Atlas-Based Parcellation

Brain parcellation divides the brain into anatomically defined regions using a standardized atlas. HINEC supports multiple atlases:

| Atlas | Type | Regions | Resolution |
|---|---|---|---|
| JHU ICBM Labels | White matter | 48 regions | 1mm |
| JHU ICBM Tracts | White matter tracts | 20 tracts | 1mm |
| HarvardOxford Cortical | Cortical gray matter | 48 regions | 1mm |
| AAL (Automated Anatomical Labeling) | Whole brain | 116 regions | 1mm |

The atlas, defined in MNI standard space, is transformed to each subject's native DWI space via the registration chain (MNI $\to$ T1 $\to$ DWI) using `preproc_atlas_resampling.m`. Nearest-neighbor interpolation preserves discrete integer labels.

### Structural Connectivity Matrices

Once both tractography and parcellation are complete, a structural connectivity matrix quantifies the strength of connections between brain regions:

1. For each streamline, identify the parcellation label at its start and end points
2. Increment the corresponding cell in an $N \times N$ matrix ($N$ = number of regions)
3. Optionally normalize and symmetrize

$$
C_{ij} = \#\{t \in \mathcal{T} : \text{start}(t) \in R_i, \; \text{end}(t) \in R_j\}
$$

The resulting matrix represents the **structural connectome** --- a complete map of white matter connections between brain regions. This is computed by `nim_plot_connectivity_matrix.m`.

### Network Neuroscience Applications

Connectivity matrices enable graph-theoretic analysis of brain networks:

- **Hub identification**: Regions with many strong connections
- **Modularity**: Communities of densely interconnected regions
- **Path length**: How many connections separate any two regions
- **Network efficiency**: How efficiently information can be transmitted

---

## 8. Preprocessing Rationale

Each preprocessing step in HINEC addresses a specific source of artifact. Uncorrected artifacts propagate through the pipeline and degrade DTI and tractography quality.

### Motion Correction

**Problem**: Head motion during the $\sim$20-minute scan causes misalignment between diffusion volumes. Even 1--2 mm of motion can corrupt tensor estimation.

**Solution**: `preproc_motion_correction.m` uses FSL's MCFLIRT to rigidly align each volume to a reference b0 image. The b-vectors are rotated to match the corrected orientations:

$$
\mathbf{g}_i' = \mathbf{R}_i \, \mathbf{g}_i
$$

where $\mathbf{R}_i$ is the rotation matrix from motion correction of volume $i$. This is a critical step, since the gradient directions must correspond to the actual image orientations.

### Eddy Current Correction

**Problem**: Rapid switching of strong diffusion gradients induces eddy currents in the scanner's conductive components. These currents create transient magnetic fields that distort the images, with different distortion patterns for each gradient direction.

**Solution**: `preproc_eddy_correction.m` wraps FSL's `eddy` tool, which models and corrects these distortions while simultaneously correcting for motion.

### Brain Extraction

**Problem**: Non-brain tissue (skull, scalp, eyes) can interfere with registration and tensor estimation.

**Solution**: `preproc_brain_extraction.m` uses FSL's BET (Brain Extraction Tool) to create a binary brain mask. A conservative threshold ($f=0.3$, more inclusive than default 0.5) is used for DWI data, which has lower contrast than structural images.

### Field Map Correction

**Problem**: Magnetic field inhomogeneities near tissue-air interfaces (e.g., sinuses, ear canals) cause geometric distortions in EPI-based diffusion images.

**Solution**: `preproc_fieldmap_correction.m` uses FSL's FUGUE with a measured field map $\Delta B_0(\mathbf{r})$ to unwarp these distortions:

$$
x' = x + \Delta B_0(x, y, z) \times t_{\text{dwell}} \times \hat{\mathbf{n}}_{\text{PE}}
$$

where $t_{\text{dwell}}$ is the effective echo spacing and $\hat{\mathbf{n}}_{\text{PE}}$ is the phase-encoding direction.

### Denoising

**Problem**: Diffusion MRI has inherently low signal-to-noise ratio, especially at high b-values. Noise in individual voxels propagates to tensor estimates and FA values.

**Solution**: `preproc_denoising.m` supports multiple methods with a fallback chain: MRtrix3 MP-PCA (preferred, exploits redundancy across diffusion directions) $\to$ FSL SUSAN (edge-preserving smoothing) $\to$ Gaussian smoothing (conservative).

### Tissue Segmentation

**Problem**: Tractography needs to know tissue boundaries to avoid producing biologically implausible results.

**Solution**: `preproc_tissue_segmentation.m` classifies voxels as white matter, gray matter, or CSF based on FA thresholds:

$$
\text{WM}: \text{FA} > 0.2, \qquad \text{GM}: 0.05 < \text{FA} \leq 0.2, \qquad \text{CSF}: \text{FA} \leq 0.05
$$

This enables Anatomically Constrained Tractography (Section 6).

---

## 9. References

### Foundational DTI and Tractography

- Basser PJ, Mattiello J, LeBihan D. (1994). "MR diffusion tensor spectroscopy and imaging." *Biophysical Journal*, 66(1):259-267. --- Original diffusion tensor formulation.

- Basser PJ, Pierpaoli C. (1996). "Microstructural and physiological features of tissues elucidated by quantitative-diffusion-tensor MRI." *Journal of Magnetic Resonance, Series B*, 111(3):209-219. --- FA and other scalar metrics.

- Mori S, Crain BJ, Chacko VP, van Zijl PCM. (1999). "Three-dimensional tracking of axonal projections in the brain by magnetic resonance imaging." *Annals of Neurology*, 45(2):265-269. --- FACT algorithm.

- Westin CF, Maier SE, Mamata H, Nabavi A, Jolesz FA, Kikinis R. (2002). "Processing and visualization for diffusion tensor MRI." *Medical Image Analysis*, 6(2):93-108. --- Tensor shape measures ($C_L$, $C_P$, $C_S$).

### Advanced Tractography

- Behrens TEJ, et al. (2003). "Characterization and propagation of uncertainty in diffusion-weighted MR imaging." *Magnetic Resonance in Medicine*, 50(5):1077-1088. --- Probabilistic tractography.

- Smith RE, Tournier J-D, Calamante F, Connelly A. (2012). "Anatomically-constrained tractography: Improved diffusion MRI streamlines tractography through effective use of anatomical information." *NeuroImage*, 62(3):1924-1938. --- ACT framework.

### Preprocessing

- Andersson JLR, Sotiropoulos SN. (2016). "An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging." *NeuroImage*, 125:1063-1078. --- FSL eddy.

- Jenkinson M, Beckmann CF, Behrens TEJ, Woolrich MW, Smith SM. (2012). "FSL." *NeuroImage*, 62(2):782-790. --- FSL software suite.

### Brain Atlases

- Mori S, et al. (2008). "Stereotaxic white matter atlas based on diffusion tensor imaging in an ICBM template." *NeuroImage*, 40(2):570-582. --- JHU atlas.

- Desikan RS, et al. (2006). "An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest." *NeuroImage*, 31(3):968-980. --- Atlas-based parcellation framework.

---

## Cross-References

- Mathematical formulas and derivations: [MATHEMATICAL_FOUNDATIONS.md](MATHEMATICAL_FOUNDATIONS.md)
- System architecture and data flow: [ARCHITECTURE.md](ARCHITECTURE.md)
- Pipeline implementation details: [PIPELINE.md](PIPELINE.md)
- Tractography algorithm details: [TRACTOGRAPHY.md](TRACTOGRAPHY.md)
- Preprocessing step details: [PREPROCESSING.md](PREPROCESSING.md)
- RKF45 adaptive integration: [RKF.md](RKF.md)
