# Multilayer moving frames (MMF-ml: `algorithm: mmf`, `field: lift`)

MMF-ml is the moving-frames tracker of [MMF_TRACTOGRAPHY.md](MMF_TRACTOGRAPHY.md)
(MMF-4h: CF4 steps of `advanceFrame` with the hold scheme) run on a **multilayer
field** instead of the principal eigenvector. The DTI field has one direction per
voxel, and where two bundles cross that direction belongs to neither of them.
MMF-ml gives every bundle through a voxel its own **layer**: a direction, a frame
and a connection 1-form. A tract reads only the layer that continues its own
heading. Everything else is unchanged: the step, the integrator, the
termination, the coast through missing data and the exit test are those of
MMF-4h.

**MMF-ml leaves the professor's scheme untouched.** Sehun Chun's
`advanceFrame.m` and `nim_tractography_mmf.m` are not edited. MMF-ml runs on
copies of them, `advanceFrame_ml.m` and `nim_tractography_mmf_ml.m`, which add
one thing: the connection, the frames and the confidence may be functions
f(X, T) of the position and the heading. `runTractography` sends
`algorithm: mmf` to the copy only for `field: lift` (or for the history veto,
`mmf.veto_arc > 0`). Every other `mmf` run uses the originals.
`tests/unit/TestAdvanceFrameMl.m` checks that, with neither addition in use,
the copies take exactly the same steps as the originals.

```
./bin/run_tractography.sh mmf_ml --score          # config/mmf_ml.yml
```

The layers are the **lifts** of the bundles. A bundle's direction field, lifted
to the space of positions and directions, is a sheet over the brain. Code and
config keys keep that name: `nim_lift_field`, `nim_lift_connection`,
`field: lift`, `tractography.lift.*`.

## Notation

| symbol | meaning | domain / units |
|---|---|---|
| Omega | the voxel grid of the scan, nodes x in Z^3 | voxel coordinates |
| RP^2 | lines through the origin, S^2 / {u ~ -u} (a fibre direction has no sign) | |
| M = Omega x RP^2 | positions and fibre directions; pi(x, u) = x | 5-dimensional |
| g_i, b_i | unit gradient direction and b-value of measurement i | s/mm^2 |
| y_i(x) = S_i(x) / S_0(x) | the attenuation at node x | dimensionless |
| n_j, j = 1..J | dictionary directions, near-uniform on the half sphere (Fibonacci), J = `lift.ndirs` | unit vectors |
| a_j(x) >= 0 | amount of fibre along n_j | |
| u_k(x) | direction of layer k at node x, k = 1..P (P = `lift.max_lifts`) | unit vector, sign free |
| e1, e2, e3 | the moving frame; e1 = tangent | unit vectors |
| w_ij | connection 1-forms, w_ij,k = e_j . d e_i / dx_k | rad / voxel |
| c(X, T) | read weight of the layer of heading T at the point X | in [0, 1] |
| CL = (l1 - l2) / l1 | linear shape of the diffusion tensor (eigenvalues l1 >= l2 >= l3) | in [0, 1] |

## 1. From the signal to the layers (`nim_lift_field`)

**Forward model.** A voxel holds sticks of fibre along the n_j, plus isotropic
compartments that carry no direction:

```
y_i(x) = sum_j a_j(x) K(g_i, n_j) + sum_d c_d(x) exp(-b_i d)
K(g, n) = exp(-b (l_perp + (l_par - l_perp) (g . n)^2))
```

`l_par`, `l_perp` are the single-fibre diffusivities: the median DTI l1 and
(l2 + l3)/2 of the tissue voxels whose FA is above the `lift.ref_fa_pct`-th
percentile. The isotropic compartments are d in {l_perp, 3e-3 mm^2/s}
(`lift.iso: auto`), the hindered water and the free water. Without them, the
stick sum would build isotropic signal out of many sticks, and those sticks
would become false layers.

**Deconvolution.** At every tissue node, non-negative least squares
(`lsqnonneg`) with a prior of few layers, mu = `lift.sparsity`:

```
min over a, c >= 0 of  || K a + K_iso c - y ||^2 + mu (sum_j a_j)^2
```

The prior is one extra row of the system, sqrt(mu) (1 ... 1 0 0) against 0.
It is needed because a sum of many sticks can reproduce broad signal as well as
the isotropic compartments can. mu = 0.03 was the best value tried on ISMRM 2015
(see [Settings](#settings-of-configmmf_mlyml)).

**Orientation density.** p(x, u) = sum_j a_j(x) W(u . n_j), with the Watson
kernel W(t) = exp(kappa (t^2 - 1)), kappa = `lift.kappa` (20).

**The layers at a node** are the local maxima of p(x, .):

- **Search.** Maxima are found over the dictionary (12 degree neighbourhood),
  then refined to the continuous maximum by the fixed point of the gradient
  condition u <- sum_j a_j W(u . n_j)(u . n_j) n_j / | . |.
- **Kept when all of these hold:**
  - p >= `lift.rel` x max p;
  - lobe mass >= `lift.mass_min`. The lobe mass m = sum_{|u.n_j| >= cos 20 deg} a_j / sum_j a_j estimates the volume fraction.
  - at least `lift.sep` degrees from a stronger kept layer;
  - strongest first, at most `lift.max_lifts`.

**Continuity through lost lobes** (`lift.bridge`, 2 nodes). A layer's lobe can
vanish for m <= 2 consecutive nodes along u while the layer continues beyond
them within 20 degrees. The gap nodes then receive the mean direction of the
two ends (`L.filled`), because fibres do not end inside white matter. On the
ISMRM 2015 coronal slice y = 54, 113 of the 154 layers that are not maxima of
their own node's density are such fills.

The output is `nim.lift_U` [X Y Z P 3] (unit directions), `nim.lift_ok`
[X Y Z P] (which exist) and `nim.lift_mass`. `nim_field` builds them once and
caches them next to the nim. The cache is `<nim>_lift_<hash>.mat`, where the
hash covers every setting that changes the layers, the version of
`nim_lift_field` and the nim file.

## 2. The connection of each layer (`nim_lift_connection`)

The connection is computed **at the nodes, before tracking** ("field first,
then the ODE"), exactly as the nodal connection of `field: dti`. The difference
is that it is computed per layer.

**Branch.** Each layer must first be continued to its neighbours. A layer u at
node x continues to x +- e_k (k = 1, 2, 3) as the layer v of that neighbour
closest to u, with |u . v| >= cos 30 deg; if no layer qualifies, the branch
ends. It continues to x +- 2 e_k as the layer closest to that match. This
matching separates crossing bundles: their layers are more than 30 degrees
apart, so a branch never jumps from one bundle to the other.

**Derivative** (`lift.stencil: fd4`, the stencil of `field: dti`). The
dyadic M = u u' is sign free, and on the matched branch:

```
dM/dx_k = (M_-2 - 8 M_-1 + 8 M_+1 - M_+2) / 12     branch reaches +-2 on both sides
        = (M_+1 - M_-1) / 2                        it reaches +-1
        = least-squares slope over the reached t   otherwise
dU/dx_k = (dM/dx_k) u
T(a, b, k) = u_a dU_b/dx_k                         (the nodal connection tensor, sign free)
```

`lift.stencil: ls` (the code default) uses the least-squares slope everywhere.
It is smoother and reads less curvature on bends.

## 3. Reading a layer at a point

At the point X with incoming heading T, each of the 8 corner nodes of X
contributes the layer u closest to T if its angle to T is within `lift.read_tol`
(30 deg); a corner without one contributes nothing. With trilinear weights
omega_q, summed over the contributing corners, the read weight is:

```
c(X, T) = sum_{corners q that contribute} omega_q(X)                 in [0, 1]
```

If c < `lift.read_min` (0.125) nothing is read. Otherwise the dyadic and the
tensor T are interpolated with the weights renormalised over the contributing
corners:

```
e1      = principal axis of sum_q omega_q u_q u_q' / c,  signed along T
e2      = the Cartesian axis least aligned with e1, projected on its normal plane;  e3 = e1 x e2
w_1j,k  = e_j . sum_a e1_a T(a, :, k),  j = 2, 3;     w_23 = 0  ('Twist' 'zero')
```

With one layer per voxel and u_1 = v1, this is the scheme of `field: dti`:
MMF-ml with P = 1 is MMF-4h on DTI. `lift.interp: sheet_linear | sheet_cubic`
reads the whole 4 x 4 x 4 neighbourhood as one sheet instead. It scored no
better (0.539 against 0.540) and is off.

## 4. Tracking

The tract solves the moving-frames ODE of `advanceFrame` with the field read on
its own heading:

```
dx/ds = e1,   dE/ds = Omega(x, e1) E,   Omega_1j = w_1j(x, e1) . e1   (CF4, h = integrator.step)
```

**The hold** (`mmf.hold: hold`). In MMF-4h the confidence is the tensor's CL.
In MMF-ml it is the read weight c(X, T): a point where c < `mmf.hold_min`
(0.15) has lost the tract's own layer. There the tract coasts on the mean
connection of its last `mmf.hold_k` trusted nodes, at most `mmf.coast_max`
voxels. On reaching trusted data again it is re-anchored if its heading is
within `mmf.exit_tol` (1.2 rad) of the layer read there. Otherwise it stops
(`exit`).

**Seeds** start on the strongest layer of their voxel.

**Stops.** The tract stops outside the volume or the propagation mask, below
`termination.fa_min`, after `termination.max_arc`, and with `no_lift` where no
corner carries a layer within 30 degrees of the heading.

## Settings of `config/mmf_ml.yml`

| key | value | why |
|---|---|---|
| everything outside `lift` | as the scored MMF-4h runs | MMF-ml differs from MMF-4h only in the field |
| `lift.sparsity` | 0.03 | F1 0.537 against 0.522 at 0.1 |
| `lift.smooth` | 0 | smoothing along straight lines pulls the layers of curved bundles upstream: 0.537 at 0 against 0.521 at 1 voxel |
| `lift.stencil` | fd4 | the stencil of the DTI field; 0.540 against 0.537 with `ls` |
| `lift.read_min` | 0.125 | a crossing bundle's layer is missing at scattered nodes; its continuation is read from the nodes that carry it |
| `lift.read_tol` | 30 | 0.509 against 0.490 at 20 degrees, which cuts real layers on bends |
| `mmf.hold_min` | 0.15 | floor of the read weight below which the tract coasts; 0.541 against 0.540 at 0.25 |

The F1 values in this table are single ISMRM 2015 runs that change one key at
a time. They use the setup of [Results](#results) with `mmf_dti` and
`mmf.hold: hold`, `mmf.hold_min` 0.25 (except in its own row), and the other
lift keys at their preset values when they were known. The `read_tol` pair was
run with `lift.stencil: ls` and the code-default smoothing.

The other `lift` options are experimental and off: `smooth_kernel: curved`,
`continuation`, `fit`, `debias` and `interp: sheet_*`. Their help in
`docs/YAML_CONFIG.md` states what each did on the data.

## Results

ISMRM 2015 phantom, scored by `bin/run_ismrm_scoring.sh` (Renauld et al. 2023
rules). Per-bundle F1 = 2 TP / (2 TP + FP + FN) on the voxels of the bundle, and
mean F1 over the 26 bundles of the ground truth. OL = overlap, OR = overreach,
in the scorer's definition. Every row uses seed density 1, `seeding.fa_min` 0.1,
`termination.fa_min` 0.10 and `max_arc` 400.

| tracker | field | mean F1 | OL | OR | valid bundles |
|---|---|---|---|---|---|
| hinec (RK4 streamline) | DTI | 0.443 | 0.387 | 0.201 | 22 |
| hinec (RK4 streamline) | layers | 0.475 | 0.444 | 0.246 | 22 |
| MMF-4h | DTI | 0.503 | 0.499 | 0.348 | 23 |
| **MMF-ml** | layers | **0.541** | 0.559 | 0.459 | 24 |

```
./bin/run_tractography.sh mmf_dti_hold --score --set seeding.density=1 --set seeding.fa_min=0.1 \
    --set termination.max_arc=400 --set mmf.hold_min=0.15 --set mmf.exit_tol=1.2      # MMF-4h
./bin/run_tractography.sh mmf_ml --score                                               # MMF-ml
./bin/run_tractography.sh mmf_dti --score --set algorithm=hinec --set seeding.density=1 \
    --set seeding.fa_min=0.1 --set termination.max_arc=400 \
    [--set field=lift --set lift.sparsity=0.03 --set lift.smooth=0]                   # hinec
```

**Where it gains and loses, MMF-ml against MMF-4h, per bundle:**
- **Gains, where bundles cross:**

  | bundle | F1 change |
  |---|---|
  | CA | 0 → 0.48 |
  | OR_right | +0.15 |
  | BPS left / right | +0.14 / +0.17 |
  | CC_ventro_striatal2 | +0.13 |
  | MCP | +0.09 |

- **Losses, curved left association bundles:**

  | bundle | F1 change |
  |---|---|
  | OR_left | −0.15 |
  | Cingulum_left | −0.06 |
  | UF_left | −0.06 |
  | ILF_right | −0.04 |

- **Overall:** 14 bundles up and 11 down. A paired Wilcoxon signed-rank test
  over the 26 bundles gives p = 0.17, so on one phantom the gain in mean F1 is
  not yet shown to be more than bundle-to-bundle scatter.
- **Where the losses come from:** the cingulum loss stays when `hold_min`,
  `angle_max`, the stencil, the read kernel and the dictionary size are varied.
  It is consistent with a property of the layers themselves. The cause is not
  yet located on individual streamlines.

**With the streamline tracker (hinec):** the layers also help: +0.033 mean F1,
17 bundles up and 5 down, Wilcoxon p = 0.013. For MMF the gain is +0.037. The
interaction (the MMF gain minus the hinec gain, per bundle) is +0.004, with
p = 0.24. So the field's gain and the tracker's gain add up; neither replaces
the other. The hinec tracker reads the layers as a peak field: the layer
nearest the tangent, interpolated trilinearly.

**Seed density.** At density 2, MMF-ml scores 0.546 (OL 0.601, OR 0.571) and
MMF-4h 0.516.

## Files

- `src/nim_calculation/nim_lift_field.m`: signal → layers (section 1)
- `src/nim_tractography/nim_lift_connection.m`: branches, nodal connection, the read at a point (sections 2–3)
- `src/nim_tractography/nim_field.m`: builds and caches the layers for `field: lift`
- `src/nim_tractography/nim_tractography_mmf_ml.m`: the tracker; the layers replace the DTI frame field, and the read weight replaces CL as the hold's confidence
- `src/nim_tractography/advanceFrame_ml.m`: the scheme, with handles that may read the heading, f(X, T)
- `src/nim_tractography/nim_tractography_hinec.m`: the streamline tracker on the layers (nearest layer to the tangent), the control above
- `config/mmf_ml.yml`: the preset
- `tests/unit/TestLiftField.m`, `tests/fixtures/make_lift_crossing_nim.m`: synthetic crossings
- `scripts/lift_field_check.py`: compares the layers with the ISMRM 2015 ground-truth directions, for checking only; the tracker never reads ground truth
