# Solution Verification

Streamline tractography is the numerical solution of an ordinary differential
equation, $d\mathbf{x}/ds = \mathbf{v}(\mathbf{x})$. This page establishes that
HINEC's solver converges under refinement of both discretisations — the
integration step and the spatial sampling of the direction field — and measures
the rate of each.

!!! abstract "Result"
    Both axes converge monotonically at stable power-law rates. Under step
    refinement the integrators reach their formal orders: Euler $0.99$, RK2
    $2.01$, RK4 $\mathbf{4.00}$ with a $C^2$ spline interpolant. Under spatial
    refinement the order is $\approx 2.25$. **Spatial resolution is the binding
    constraint**: past a moderately coarse step there is little left to gain from
    refining time.

!!! note "Scope — verification, not validation"
    Verification asks whether the equations are being solved correctly.
    Validation asks whether the answer corresponds to anatomy. This page is the
    former only, and no claim about anatomical accuracy follows from anything
    below. The ISMRM ground-truth bundles appear here solely to define a seed
    region; for how the streamlines actually measure against them, see
    [ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md), where between 28% and 62% of
    produced streamline length leaves the bundle it was seeded in.

---

## Method

### The problem

Streamline tractography integrates $d\mathbf{x}/ds = \mathbf{v}(\mathbf{x})$,
where $\mathbf{v}$ is the principal eigenvector field and $|\mathbf{v}| = 1$, so
the parameter $s$ is arc length. Seeds are the 8,384 voxels of the ISMRM
ground-truth `UF_right` bundle mask with FA > 0.15.

The angle termination criterion is **disabled** (`angle_max: 0`). It is a
stopping rule whose budget scales with the step size, and a rule that changes
with the refinement parameter is not a fixed problem to converge to.

![Right uncinate fasciculus: ours against the ISMRM ground truth](img/bundle_uf_right.png)

*The seed region. Grey is the ISMRM ground-truth `UF_right` bundle, red is ours,
clipped to the bundle corridor for display. Shown only to establish what the
convergence ladders are tracking through — it is a validation figure and the
ladders do not depend on it.*

!!! note "Do not read this as an accuracy claim"
    Applying the scorer's full definition (endpoint pair plus containment) leaves
    streamlines that are 98.5% inside the ground-truth envelope, but that number
    describes the survivors of a filter that discards most of what was produced.
    Measured on the unfiltered output, **62% of `UF_right` streamline length
    leaves the bundle corridor**. See
    [ISMRM Scoring](ISMRM_SCORING_ANALYSIS.md#reconstruction-against-ground-truth).

### The error metric

Each run is compared against a reference run. Both are streamline *sets*, so the
comparison must decide which streamline corresponds to which, and then which
point corresponds to which. `nim_convergence_error` does this in six steps:

1. **Pair by seed.** Matched through `track_meta.seed_index`, not by position in
   the array — the tracker compacts out seeds that produce nothing, so
   "streamline 7" is a different fibre in the two runs.
2. **Anchor at the seed.** A track is stored as
   `[reversed backward half; seed; forward half]`, so index 1 is a *termination*
   point. Arc is measured outward from the seed.
3. **Parameterise by nominal arc.** The point $n$ steps from the seed is at arc
   $nh$ *by construction*, since each step advances exactly $h$ of arc.
   Cumulative chord length is not used: it underestimates true arc by
   $h^2/(24R^2)$ per unit length, putting the two runs on parameters that
   disagree at $O(h^2)$.
4. **Window to ±10 voxels** about the seed, using only streamlines that reach
   that in both directions in both runs. This excludes termination, which is a
   discontinuous function of $h$.
5. **Compare the test track's own stored points** against a splined reference.
   The test points are exact integrator output and are not interpolated;
   resampling them linearly imposes a chord-error floor of $\approx h^2/(8R)$
   that has nothing to do with the integrator.
6. **Reduce.** Median over the streamlines measurable at *every* rung — a fixed
   population. The 95th percentile is tracked separately.

!!! warning "Steps 3 and 5 both matter at $O(h^2)$"
    Either one, done naively, imposes a floor no integrator can beat and reports
    order 2 for every method regardless of the tableau — which is
    indistinguishable from a genuine physical cap **unless methods with known
    formal orders are run alongside as a control**.

### Reference solutions

One shared reference per experiment, always the most accurate configuration
available, **never a refined copy of the method under test**. Referencing Euler
to Euler leaves the reference carrying ~25% of the error being measured at the
finest rung, and inflates Euler's observed order from 0.99 to 1.24.

---

## Refinement in time

Nine step sizes, refined in ratios of $\sqrt{2}$ from 0.5 to 0.03125 voxels; reference RK4
at $h = 0.0078125$; fixed population $n \approx 1300$.

| method | field | formal | $h=0.5$ | $h=0.031$ | **observed** |
|---|:--:|:--:|--:|--:|--:|
| Euler | — | 1 | 1.33e-1 | 8.41e-3 | **0.99** |
| RK2 | — | 2 | 3.43e-3 | 1.30e-5 | **2.01** |
| RK4 trilinear | $C^0$ | 4 | 5.05e-4 | 1.91e-6 | **2.00** |
| RK4 cubic | $C^1$ | 4 | 6.58e-5 | 1.35e-8 | **3.06** |
| RK4 spline | $C^2$ | 4 | 2.12e-5 | 3.05e-10 | **4.00** |

![Convergence under step refinement](img/convergence_time.png)

Euler and RK2 land on their formal orders. Because those two answers are known
in advance they serve as a **control on the measurement itself**: a metric that
reports 0.99 and 2.01 for methods whose orders are 1 and 2 is measuring the
integrator and not an artefact of its own.

**The observed order of RK4 is set by the smoothness of the interpolant** — one
order per continuous derivative, until the tableau's formal order is reached.
The three RK4 rows differ *only* in the interpolation kernel.

Local slopes are stable across every triplet (the asymptotic-range condition),
and the 95th percentile converges at the same rate as the median, so the whole
distribution converges rather than a well-behaved majority.

---

## Refinement in space

The direction field is sampled on a grid of spacing $1/u$ voxels before the
interpolants are built — see [`interpolation.upsample`](YAML_CONFIG.md). The
coordinate frame is unchanged, so positions, step sizes and lengths stay in
native voxel units and runs at different factors are directly comparable.

The integration step is pinned at $h = 0.125$ with RK4 and spline interpolation,
where the temporal error is $\approx 8\times10^{-8}$ voxels — four to six orders
below the spatial errors being measured.

| spacing (vox) | median | p95 | local $p$ |
|--:|--:|--:|--:|
| 8.000 | 1.0713 | 6.709 | — |
| 5.657 | 1.0167 | 5.148 | 0.15 |
| 4.000 | 0.5303 | 3.899 | 1.88 |
| 2.828 | 0.2924 | 2.520 | 1.72 |
| 2.000 | 0.1188 | 1.839 | 2.60 |
| 1.414 | 0.0476 | 0.193 | 2.64 |

![Convergence under spatial refinement](img/convergence_space.png)

Fitted order **2.25** over the four asymptotic rungs ($n = 761$), rising to
2.43–2.62 over the finest three. Error falls by a factor of 22 as the grid
refines 5.7×, and the p95 falls with it (6.71 → 0.19). The two coarsest rungs
are pre-asymptotic: at 8-voxel spacing the field is barely resolved and the
error saturates near one voxel.

!!! note "The limit is the native-resolution interpolant, not the anatomy"
    Sampling *above* the acquisition grid returns error that is exactly zero — a
    spline through samples of itself reproduces itself — which says the 2 mm
    samples already carry all the information there is. This ladder demonstrates
    grid convergence *of the algorithm*; whether 2 mm resolves the anatomy is a
    validation question requiring different data.

---

## Which axis binds

RK4 with spline interpolation converges at order 4 in time and $\approx 2.3$ in
space. **Resolution is the binding constraint.** Beyond a moderately coarse step
there is little to gain from refining time; effort is better spent on the
spatial representation of the direction field.

Practically: pair `interpolation.method: spline` with `integrator.method: rk4`.
At $h = 0.5$ that combination is 162× more accurate than RK2 and 6300× more
accurate than Euler, so the larger step fourth-order accuracy affords can be
taken without giving the accuracy back.

---

## What does not limit the observed order

Three candidate explanations were tested directly and **rejected**. The negative
results constrain what the numbers above can mean.

| tested | result |
|---|---|
| **Field noise** — a field whose every voxel direction is rotated through a random angle $\sigma$, swept 0–20° against the real field's measured median of 4.84° between adjacent voxels, with random per-voxel sign flips | RK4 stayed above order 2.5 and 2–756× ahead of RK2 throughout |
| **Coherent discontinuities** — a surface across which the principal direction jumps by up to 60°, crossed mid-integration (the shape a fibre crossing takes) | RK4 still finished near order 4 and 9× ahead of RK2 |
| **Stage fallbacks** — `rk4_integration_step` substitutes the previous stage when interpolation refuses a probe below the FA floor | 0.80% of probes, independent of $h$ — present, but far too small to govern the observed order |

---

## Threats to validity

- **One bundle, one subject.** All results are the right uncinate fasciculus of
  the ISMRM 2015 phantom. Nothing here establishes generality.
- **Self-convergence, not accuracy.** The reference is our own refined solution.
  The pipeline converges to *something*; that it is the right thing is a
  separate claim.
- **The spatial ladder is not fully asymptotic.** Local order is still rising at
  the finest rung, so 2.25 is a lower bound.
- **Termination is quantised to one step.** Measured arc-length difference is
  0.87–0.98 $h$ across a 16× refinement, so the endpoint converges at $O(h)$
  while the interior converges at the integrator's order. The windowed metric
  excludes this deliberately; it remains a real defect.
- **MMF is untested here.** The angle criterion in
  `nim_tractography_mmf_connframe` is covered only by a source-text check, not a
  behavioural test.
- **The whole-brain benchmark has not been run.** Every ISMRM score collected so
  far comes from a single-ROI submission, which is not a valid use of a scorer
  that segments a whole-brain tractogram into 26 bundles.

---

## Reproduction

```bash
# time ladder, one rung
./bin/run_tractography.sh reference \
    --set seeding.roi=UF_right --set seeding.fa_min=0.15 \
    --set interpolation.method=spline --set integrator.method=rk4 \
    --set integrator.step=0.125 --set output.arc_step=0

# space ladder, one rung — add
    --set upsample=0.5
```

Analysis is `nim_convergence_error(test_run, reference_run, struct('prefix_arc', 10))`,
which returns per-seed errors so a fixed population can be intersected across
rungs.

Verification suite: **82/82** passing, including behavioural cases that
integrate real streamlines through synthetic fields of known analytic curvature —
order verification for each integrator, the angle criterion firing exactly at
the field's true turning rate, sign-invariance of the dyadic interpolation, and
the endpoint/containment bundle gates.
