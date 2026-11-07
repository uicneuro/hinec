# HINEC Configuration Reference

<!-- GENERATED FILE - do not edit by hand.
     Regenerate with:  nim_config_docs('docs/YAML_CONFIG.md')
     Source of truth:  src/nim_utils/nim_config_schema.m -->

Configuration is YAML with **two levels of nesting** below a section
(`section` -> `group` -> `key`). A third level is a parse error, deliberately:
deeper nesting was judged more confusing than helpful.

Every key is optional. Anything you omit takes the default below, so a working
config can be very short.

## Minimal example

```yaml
tractography:
  algorithm: hinec
  seeding:
    roi: [SLF_L, SLF_R]
```

## Overriding from the command line

Every parameter in this reference can be overridden with `--set`, on both
`bin/run_hinec.sh` and `bin/run_tractography.sh`:

```bash
./bin/run_tractography.sh hinec_dti --set tractography.integrator.step=0.05
./bin/run_tractography.sh hinec_dti --set integrator.step=0.05    # section assumed
./bin/run_tractography.sh hinec_dti --set upsample=2              # bare leaf, if unique
./bin/run_hinec.sh data/x x.mat config/y.yml --set preprocessing.run_eddy=false
./bin/run_tractography.sh hinec_dti --set seeding.roi=[41,42]     # lists
```

Overrides are checked against the schema. An unknown or mistyped key is an
error, never a silent no-op. A bare leaf name that is ambiguous (`method`,
`fa_min`) is rejected with the candidate full paths.

## Parameters

The **Applies to** column says which tracking algorithms actually read the
parameter. A key marked `hinec` is ignored by `standard`, `mmf` and `stitching`.

### `tractography`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `algorithm` | string | `hinec` | all | Tracking algorithm. hinec = interpolated streamlines (RK); mmf = moving frames tractography (advanceFrame: the frame is transported by the connection of the data, exponential step; docs/MMF_TRACTOGRAPHY.md); standard = FACT; stitching = short tractlets; template = the minimal worked example of the tracker interface (docs/TRACKER_INTERFACE.md). |
| `field` | string | `dti` | hinec, mmf, stitching | Direction source. dti = DTI principal eigenvector; csd = CSD FOD peaks (hinec and stitching); lift = the lifts of the bundles on M = Omega x RP^2 built from the raw signal (nim_lift_field), one frame and one connection 1-form per bundle where bundles cross (with mmf this is MMF-ml, multilayer moving frames, run by nim_tractography_mmf_ml; see docs/MMF_ML.md and tractography.lift); with hinec the lifts are a peak field, the lift nearest the incoming tangent at every stencil voxel as for csd. |
| `act` | logical | `false` | hinec | Anatomically constrained tracking using WM/GM/CSF masks. |
| `diagnostics` | logical | `true` | all | Write per-run diagnostic reports. |

#### `tractography.integrator`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `method` | string | `rk4` | hinec, mmf, stitching | Numerical stepping scheme. NOTE: this is a method NAME, not an order claim. rkf45 here is Dormand-Prince: it advances on the 5th-order solution and uses the embedded 4th-order one only to size the next step. For mmf: euler, rk2, rk4 select the Lie group integrator of order 1, 2, 4 (Lie-Euler, exponential midpoint, CF4) around the exponential step; rkf45 is rejected (no embedded error estimate exists for it). |
| `step` | numeric | `0.2` | all | Integration step in voxels. Fixed step, or the initial step for rkf45. |
| `tolerance` | numeric | `0.01` | hinec | Adaptive error tolerance in voxels. rkf45 ONLY. |
| `step_min` | numeric | `0.01` | hinec | Minimum adaptive step. rkf45 ONLY. |
| `step_max` | numeric | `1` | hinec | Maximum adaptive step. rkf45 ONLY. |
| `safety` | numeric | `0.9` | hinec | Adaptive step safety factor. rkf45 ONLY. |
| `adaptive` | logical | `true` | hinec | Adaptive step-size control. rkf45 ONLY. rkf45 with adaptive:false is fixed-step RKF45 (a real third mode, kept for completeness). |

#### `tractography.interpolation`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `method` | string | `trilinear` | hinec, mmf, stitching | Spatial interpolation kernel for the direction field. These differ in SMOOTHNESS, which caps how much of an integrator formal order is reachable: a Runge-Kutta method of order p needs a right-hand side with p continuous derivatives. trilinear is C0 (kinked at every voxel face), cubic is C1 (Keys cubic convolution, NOT a spline - its second derivative jumps), spline is C2 (a genuine cubic spline). Measured here: RK4 reaches observed order 2.00 on trilinear, 3.06 on cubic and 4.00 on spline - one order per continuous derivative, on an unchanged tableau. |
| `upsample` | numeric | `1` | hinec, mmf | Spatial sampling factor for the direction field: the field is sampled on a grid of spacing 1/upsample voxels before the interpolants are built. 1 = the acquisition grid. Above 1 refines toward the continuous field the samples imply; BELOW 1 coarsens, discarding spatial information, which is how the space axis of a convergence study is swept. The coordinate frame is unchanged, so positions, step sizes and lengths stay in native voxel units and runs at different factors are directly comparable. Note the u -> infinity limit is the native-resolution interpolant, not ground-truth anatomy. |

#### `tractography.seeding`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `density` | numeric | `8` | all | Seeds per seeded voxel - honoured exactly. Placed on a deterministic sub-voxel lattice (perfect cubes) or a spread subset of one (other counts), so runs are reproducible. |
| `strategy` | string | `uniform` | all | Seed placement. uniform = deterministic lattice (reproducible); random = jittered. |
| `fa_min` | numeric | `0.05` | all | Minimum FA for a voxel to be seeded. Default 0.05 excludes CSF only. |
| `cl_min` | numeric | `0` | all | Minimum linear shape CL = (lambda_1 - lambda_2)/lambda_1 of the voxel tensor for a voxel to be seeded (nim_seed_cl_mask). 0 = off. A crossing fitted as one tensor has CL near 0, so e.g. 0.15 (the mmf.hold_min of the hold presets) keeps seeds out of crossings. |
| `coherence_min` | numeric | `0` | all | Minimum coherence of the principal direction for a voxel to be seeded (nim_seed_coherence_mask): the largest eigenvalue of the mean dyadic v1*v1' over the 3 x 3 x 3 block, in [1/3, 1] (1 = all neighbours point the same way). 0 = off. |
| `roi` | list | `[]` | all | Restrict seeding to these atlas regions. JHU indices and/or names, freely mixed. Empty = whole brain mask. |
| `roi_dilate` | numeric | `0` | all | Dilate the seed ROI by this many voxels before seeding. |

#### `tractography.termination`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `fa_min` | numeric | `0.1` | all | Stop tracking below this FA. (NOT the legacy fa_threshold - see nim_config_retired.) |
| `fa_mode` | string | `absolute` | all | How the FA stop threshold is set. absolute = stop below termination.fa_min (unchanged behaviour). percentile = stop below the termination.fa_percentile-th percentile of THIS scan's FA inside its brain mask, resolved per scan in runTractography before any tracker runs; fa_min is then ignored. Removes scanner-level FA offsets (the same absolute FA cuts different fractions of the brain on different scanners). The resolved value is logged and kept in options.termination_fa, with the configured one in options.termination_fa_config. |
| `fa_percentile` | numeric | `26.3` | all | Percentile of the scan's in-mask FA used as the stop threshold when fa_mode = percentile. Default 26.3 = the median percentile rank of FA 0.10 over the 76 ds000206 scans (the reference protocol), so it keeps the meaning of fa_min = 0.10. |
| `angle_max` | numeric | `225` | all | Maximum turn in DEGREES PER VOXEL OF ARC, i.e. a minimum radius of curvature R = 57.3/angle_max voxels. The default 225 is the classic 45 deg/step evaluated at the default step of 0.2. Step-invariant: the limit for one step is angle_max x the NOMINAL step arc (never the realised chord - that would make the limit method-dependent), so refining the step does not loosen the constraint. CEILING: tangents are sign-aligned because v1 is a line field, so a measured turn never exceeds 90 deg - any angle_max above 90/step is INERT, not merely loose, and 225 goes inert for any step >= 0.4. 0 disables the criterion outright, which is what a control run should use. |
| `max_arc` | numeric | `400` | all | Maximum track arc length in voxels. Step-invariant: max_steps is derived as ceil(max_arc/step). |
| `min_arc` | numeric | `15` | all | Discard tracks shorter than this arc length in voxels. |

#### `tractography.filter`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `include_roi` | list | `[]` | all | Keep only tracks touching these regions (waypoint). |
| `exclude_roi` | list | `[]` | all | Discard tracks touching these regions. |
| `mode` | string | `all` | all | Whether a track must touch ALL include regions or ANY of them. |
| `roi_dilate` | numeric | `0` | all | Dilate the include/exclude masks by this many voxels before testing. |
| `endpoints_in` | list | `[]` | all | Two regions; keep a track only if one END lands in the first and the other END in the second, either way round. This is an ENDPOINT test, not a waypoint test - include_roi asks whether a track passes through a region, this asks where it stops. It is half of how the ISMRM 2015 scorer defines a bundle (head + tail). |
| `contained_in` | list | `[]` | all | Keep a track only if EVERY point lies inside these regions. The other half of the ISMRM bundle definition (all_mask): a streamline that wanders outside the corridor is not that bundle, however it ends. Tracks are kept or dropped whole, never cut. Like every filter.* key this is a SELECTION: tracks.mat still holds the full tractogram (what gets scored) and the selection is saved beside it as roi_selection.mat with the per-criterion drop counts. |
| `any_in` | list | `[]` | all | Keep a track only if it touches these regions in AT LEAST ONE point (the scorer's any_mask). Weaker than include_roi with mode all: one region, one point. |
| `length` | list | `[]` | all | Keep tracks whose total length in MILLIMETRES is within [min max]. Millimetres, not voxels - the scorer measures in RAS mm. |
| `length_x` | list | `[]` | all | Keep tracks whose NET displacement along x is within [min max] mm. Net: a track that doubles back cancels itself out. |
| `length_y` | list | `[]` | all | Keep tracks whose NET displacement along y is within [min max] mm. |
| `length_z` | list | `[]` | all | Keep tracks whose NET displacement along z is within [min max] mm. |
| `length_x_abs` | list | `[]` | all | Keep tracks whose TOTAL travel along x is within [min max] mm. Total: a track that doubles back adds to it. |
| `length_y_abs` | list | `[]` | all | Keep tracks whose TOTAL travel along y is within [min max] mm. |
| `length_z_abs` | list | `[]` | all | Keep tracks whose TOTAL travel along z is within [min max] mm. |

#### `tractography.debug`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `trace` | logical | `false` | hinec, mmf | Record a PER-STEP trace of the tracker: position, direction used, interpolated FA, turn angle, nominal vs realised step arc, sign flips, and termination reason. None of this survives to the saved polyline, which output.arc_step decimates, so it cannot be reconstructed afterwards. Off by default; see nim_track_diagnose. The mmf tracker writes the same record fields (pos, dir, fa, h, chord, turn, flip, n_steps, termination; flip is always false there, the gauge absorbs sign flips) plus drift, reanchored, kappa, s, coasted and confidence, so nim_track_diagnose reads both. |
| `trace_max` | numeric | `200` | hinec, mmf | How many seeds to trace, sampled EVENLY across the seed list (a prefix would be one corner of the volume). 0 traces every seed, which is large: allow roughly max_steps x 12 doubles per seed per direction. |

#### `tractography.output`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `arc_step` | numeric | `0` | all | Resample saved streamlines to this arc-length spacing in voxels. 0 = store every integration step. Decouples file size from step size; integration accuracy is unaffected. |

#### `tractography.csd`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `peaks_file` | string | `` | hinec, stitching | Optional MAT file with peaks [X Y Z P 3], npeaks [X Y Z], and peak_w [X Y Z P] in dataset voxel-axis coordinates. When set, use this deterministic field instead of fitting native CSD. The file must match the dataset grid; no GT or scoring data is needed. |
| `lmax` | numeric | `4` | hinec, stitching | Spherical harmonic order for CSD. Default 4, set from the data rather than convention: with the 32 directions of this acquisition lmax=6 needs 28 coefficients from 32 measurements and the response falls to r_4 = -0.035 against r_0 = 2.73, so deconvolution divides by a near-zero and amplifies noise - its primary FOD peak misses the tensor v1 by 26.6 deg in UNAMBIGUOUS single-fibre voxels, against 6.6 deg at lmax=4. lmax=2 spans exactly the tensor and cannot represent a crossing at all. |
| `max_peaks` | numeric | `3` | hinec, stitching | Maximum FOD peaks retained per voxel. |
| `peak_thresh` | numeric | `0.2` | hinec, stitching | Minimum FOD peak amplitude as a fraction of the voxel maximum. The old default of 0.5 required a second fibre population to be at least half as strong as the first, which discards ordinary unequal crossings; MRtrix uses 0.1 for the equivalent. Peaks must now also be local maxima of the FOD, so a lower threshold admits genuine fibres rather than points on the shoulder of the first peak. |
| `peak_min_sep` | numeric | `45` | hinec, stitching | Minimum angular separation between FOD peaks (degrees). |
| `n_iter` | numeric | `50` | hinec, stitching | CSD deconvolution iterations. |

#### `tractography.lift`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `smooth` | numeric | `0` | hinec, mmf | Half-width T (voxels) of the regularisation of the orientation density along the fibres: each direction n_j is smoothed along x + t n_j, t = -T..T, weights T + 1 - \|t\|; 0 = off (nim_lift_field). Default 0: the kernel runs along straight lines, so on a bend of curvature k it mixes directions up to k T apart and pulls the lobe toward the upstream direction (ISMRM cingulum, k ~ 0.5 per voxel: the lifted curvature fell to 0.02 at T = 2). ISMRM mmf lift+hold F1 0.509 / 0.521 / 0.537 at T = 2 / 1 / 0. |
| `smooth_kernel` | string | `straight` | hinec, mmf | Path of the regularisation along the fibres (lift.smooth > 0): straight = the line x + t n_j; curved = the best of the circular arcs leaving (x, n_j) with curvature lift.smooth_kappa in lift.smooth_normals normal directions (and the straight line), the density read at each arc point in the dictionary direction nearest the arc tangent (Duits and Franken, IJCV 92:231-264, 2011: smoothing on position x orientation follows curves, not lines). |
| `smooth_kappa` | numeric | `` | hinec, mmf | Curvatures (per voxel) of the arcs of lift.smooth_kernel: curved, besides the straight line; the ISMRM cingulum bends at ~0.5 per voxel. |
| `smooth_normals` | numeric | `8` | hinec, mmf | Number of normal directions m = cos(q) e2 + sin(q) e3, q = 2 pi i / smooth_normals, in which the arcs of lift.smooth_kernel: curved bend. |
| `bridge` | numeric | `2` | hinec, mmf | Nodes a lift is continued through where its lobe is lost, when it continues beyond within 20 degrees (fibres do not end inside white matter); 0 = off. |
| `rel` | numeric | `0.1` | hinec, mmf | A lift is a local maximum of p(x,.) of at least rel times the voxel maximum. |
| `mass_min` | numeric | `0.1` | hinec, mmf | Minimum lobe mass of a lift: the share of the deconvolved fibre within 20 degrees of it, an estimate of its volume fraction. |
| `max_lifts` | numeric | `3` | hinec, mmf | Maximum lifts per voxel. |
| `sep` | numeric | `15` | hinec, mmf | Minimum angle (degrees) between two lifts of a voxel. |
| `ndirs` | numeric | `300` | hinec, mmf | Dictionary directions on the half sphere for the non-negative deconvolution. |
| `kappa` | numeric | `20` | hinec, mmf | Concentration of the Watson kernel W(t) = exp(kappa (t^2 - 1)) of the orientation density p = sum_j a_j W(u.n_j). |
| `ref_fa_pct` | numeric | `90` | hinec, mmf | The single-fibre kernel is the median DTI l1 and (l2 + l3)/2 of the voxels above this FA percentile of the tissue. |
| `iso` | string | `auto` | hinec, mmf | Isotropic compartments of the forward model, which carry no direction: auto = the hindered (l_perp) and the free water (3e-3 mm^2/s); none = sticks only (false lifts where the tissue is isotropic). |
| `continuation` | numeric | `0` | hinec, mmf | Continuation prior beta (nim_lift_field): a second deconvolution whose prior of few lifts weighs direction n_j by c_j = 1 - beta max W(n_j . u^), u^ the moving-frame continuation of each neighbouring lift to the voxel (its nodal connection, curvature included); the directions a neighbouring fibre continues into cost less. 0 = off; needs lift.sparsity > 0. |
| `sparsity` | numeric | `0.1` | hinec, mmf | Prior of few lifts: the deconvolution minimises \|\|K a + K_iso c - y\|\|^2 + sparsity (sum_j a_j)^2, so broad or isotropic signal goes to the isotropic compartments instead of becoming false lifts. ISMRM 2015 slab (32 directions): 0.1 keeps one lift in 85% of the single-bundle voxels at the v1 error (7.3 vs 7.5 deg) and halves the error in two-bundle crossings (12.8 vs 24.9 deg); 0 = none. |
| `read_min` | numeric | `0.25` | mmf | The tracker reads a lift at a point when the corners carrying a lift within 30 degrees of the heading hold at least this share of the trilinear weight (renormalised over them). Inside crossings the weaker lift is missing at scattered nodes; at 0.5 31.9% of ISMRM 2015 half-tracks stopped there with 0.25-0.5 of the weight still on their lift. |
| `interp` | string | `corner` | mmf | How the lifted field is read at a point (nim_lift_connection): corner = at each of the 8 corners the lift nearest the heading, trilinear (C0, jumps where a corner switches lift); sheet_linear / sheet_cubic = the local sheet of the heading on the 4 x 4 x 4 nodes (branch-matched from the nearest corner), interpolated trilinearly / with the Keys cubic kernel (C1, interpolation.method: cubic). |
| `stencil` | string | `ls` | mmf | Derivative of each lift along its matched branch for the lifted nodal connection (nim_lift_connection): ls = least-squares slope over t = -2..2; fd4 = the stencil of the nodal connection of field: dti, (M_-2 - 8 M_-1 + 8 M_+1 - M_+2) / 12. On the ISMRM ground truth, bends with kappa > 0.15: ls reads 0.37-0.61 of the curvature, fd4 0.45-0.71. |
| `read_tol` | numeric | `30` | mmf | Degrees: a node contributes its lift to the point the tracker reads only when that lift is within read_tol of the heading. With mmf.hold the share of the weight so read is the confidence: the tract coasts where it falls below mmf.hold_min. |
| `fit` | logical | `false` | hinec, mmf | Refit the stick model of each voxel from its lifts (Levenberg-Marquardt). Off: on synthetic crossings it undoes the regularisation along the fibres. |

#### `tractography.mmf`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `mode` | string | `adaptive` | mmf | How the carried frame relates to the data frame along the tract (advanceFrame 'Mode'). transport = carried by the connection only, the data enter through the curvature and the gauge; reanchor = replaced by the data frame after every step (the direction from the data, the curvature of the step from the transport: the closest relative of streamline tracking, with an arc instead of a chord); adaptive = transport, re-anchored only after a step at whose end the carried e1 has drifted from the data's e1 by mmf.reanchor_tol or more. |
| `reanchor_tol` | numeric | `0.001` | mmf | Drift tolerance of the adaptive mode, radians: the angle between the carried e1 and the data's e1 (modulo sign) at the end of a step at or above which the frame is re-read from the data. 0 re-anchors after every step like mode reanchor; a value above pi never, like transport. |
| `gradient_delta` | numeric | `0.001` | mmf | Central-difference step, in voxels, for the derivative of the interpolated principal axis from which the connection (curvature) components are formed. |
| `connection` | string | `nodal` | mmf | How the connection 1-forms are obtained (MMFtractography_notes, Sec. 7). nodal = computed at the voxel nodes from the nodal frames and INTERPOLATED between them with interpolation.method (the scheme of the notes: Computeconnection + interpGLL); derivative = central difference (step gradient_delta) of the principal axis of the interpolated dyadic at every query point (the adapter of 2026-10-20, kept to reproduce its runs). |
| `max_turn` | numeric | `0` | mmf | advanceFrame 'MaxTurn' in radians per step (the notes' default is 0.5). 0 = derive it from termination.angle_max as before (min(angle_max*h, 90) degrees per step; pi when angle_max <= 0). |
| `hold` | string | `none` | mmf | How a tract crosses degenerate data, where the tensor is a crossing fitted as a pancake and its principal axis means nothing (see mmf.hold_min). none = it does not: the tract reads the data, follows the axis of the pancake or stops on angle_max. hold = MMF-4h: it coasts through the zone on the mean connection of its last trusted nodes (advanceFrame 'Confidence'), a circle when the twist is zero. march = MMF-4s: it coasts on the curvature and the fan-out of those nodes marched by the flatness law of the connection ('HoldCurvature' 'march'), exact on circles, logarithmic spirals and rays, and turning the curvature by the torsion on a helix. The volume, the mask and the FA floor stay stops inside the zone. config/mmf_dti_hold.yml and config/mmf_dti_march.yml turn them on; docs/MMF_TRACTOGRAPHY.md has the scores against the ground truth of the ISMRM 2015 phantom. |
| `hold_min` | numeric | `0.2` | mmf | Confidence floor of mmf.hold (advanceFrame 'HoldMin'): a point is degenerate where the linear shape CL = (lambda_1 - lambda_2)/lambda_1 of the voxel tensors, interpolated linearly and eroded by mmf.hold_guard, is below this. For a single bundle of eigenvalues 1.7, 0.3, 0.3 (CL 0.82) 0.2 is an equal crossing of about 80 degrees and 0.4 of about 70, but CL also falls with the anisotropy of the data, so the floor belongs to the dataset: 0.15 and 0.2 were the best of the floors scored on the ISMRM 2015 phantom, and on tissue of lower anisotropy 0.2 can mark half of the voxels. The tracker reports the share of the tissue voxels below the floor, before the erosion (meta.mmf.zone_fraction), and warns above 40 percent (mmf:holdZone): lower the floor there. |
| `hold_k` | numeric | `3` | mmf | Trusted nodes in the running mean the coast starts from (advanceFrame 'HoldK'), a positive integer. |
| `hold_guard` | numeric | `0.5` | mmf | Erosion radius of the confidence, voxels: its minimum over the point and 14 points at this distance, which keeps the last trusted nodes, whose connection the coast starts from, off the voxels next to the zone. CL is interpolated over the voxels of tissue alone and the minimum reads nothing else, so that neither the edge of the mask nor the FA floor is a zone. A larger guard is a longer coast. |
| `coast_max` | numeric | `12` | mmf | Longest coast, voxels of arc (advanceFrame 'CoastMax' = coast_max/step steps): a tract still in the zone after this stops with reason coast. |
| `exit_tol` | numeric | `0.5` | mmf | Largest angle, radians, between the carried e1 and the data's e1 with which a tract may leave the zone (advanceFrame 'ExitTol'): it is re-anchored then, and above this the coast has left the bundle and the tract stops with reason exit (with reason angle if the data at the end of the step are further from the heading at its start than angle_max allows). It must sit above the disagreement of neighbouring voxels, a quarter of which differ by more than 0.2 rad on real data, so it is not mmf.reanchor_tol; 0.5 rad is 29 degrees. Adaptive and reanchor modes only: in transport mode nothing is re-anchored or tested on leaving. |
| `exit_delay` | numeric | `0` | mmf | Arc of coast beyond the edge of the zone before the exit test, voxels (advanceFrame 'ExitDelay' = exit_delay/step steps). 0 = test at the first trusted node. These steps neither count towards coast_max nor reset it, and the volume, the mask and the FA floor stop them. |
| `veto_arc` | numeric | `0` | mmf | History veto of mmf.hold, look-back in voxels of arc (advanceFrame_ml 'VetoK' = round(veto_arc/step); a run with it goes to nim_tractography_mmf_ml); 0 = off. A trusted step whose end the data see more than mmf.veto_tol from the heading veto_arc back is not taken: the tract returns there and coasts on its ring (a crossing whose tensor looks linear along the other bundle). |
| `veto_tol` | numeric | `0.5236` | mmf | History veto: largest angle, radians, between the data e1 at a step end and the heading mmf.veto_arc back (advanceFrame_ml 'VetoTol'; 0.5236 = 30 deg). |

#### `tractography.stitching`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `fragment_arc` | numeric | `3` | stitching | Total bidirectional fragment arc in voxels; 3 is 6 mm at 2 mm resolution. |
| `min_fragment_arc` | numeric | `1` | stitching | Discard local fragments shorter than this voxel arc before graph assembly. |
| `join` | logical | `true` | stitching | Enable endpoint stitching. False is the fragment-only ablation. |
| `join_radius` | numeric | `1` | stitching | Maximum endpoint gap in voxels. |
| `join_angle` | numeric | `30` | stitching | Maximum endpoint tangent and forward-gap alignment angle, degrees. |
| `bridge_angle` | numeric | `45` | stitching | Maximum bridge tangent disagreement with the local field, degrees. |
| `search` | string | `knn` | stitching | Endpoint search: knn caps before compatibility; forward searches the full forward cone and caps compatible candidates across distance/direction sectors. |
| `neighbors` | numeric | `16` | stitching | Endpoint candidate limit: raw nearest endpoints for knn; compatible sector-balanced endpoints for forward. |
| `curvature_tolerance` | numeric | `0.25` | stitching | Maximum endpoint curvature difference, inverse voxels, applied by the endpoint graph when a per-endpoint curvature/torsion table is supplied to nim_stitch_fragments (the pipeline currently supplies none). |
| `torsion_tolerance` | numeric | `0.5` | stitching | Maximum endpoint torsion difference when curvature is nonzero, inverse voxels; same condition as curvature_tolerance. |
| `max_fragments` | numeric | `200000` | stitching | Fail before generation if seed count times possible peaks exceeds this explicit memory limit. |

### `preprocessing`

| Key | Type | Default | Applies to | Description |
|---|---|---|---|---|
| `run_denoising` | logical | `true` | - | Run denoising before tensor fitting. |
| `denoise_method` | string | `dwidenoise` | - | Denoising backend. Falls back to gaussian when MRtrix3 is absent. |
| `run_motion_correction` | logical | `true` | - | Run motion correction (FSL). |
| `run_eddy` | logical | `true` | - | Run eddy-current correction (FSL). |
| `improve_mask` | logical | `true` | - | Refine the brain mask using FA. |
| `mask_file` | string | `` | - | Optional dataset-supplied brain or tracking mask. Empty runs the normal T1/DWI brain extraction stage. |
| `atlas_file` | string | `` | - | Optional external 3-D integer label atlas (.nii/.nii.gz), already aligned to the DWI grid. Overrides atlas_type and bypasses built-in atlas registration. Background is 0; no automatic resampling or registration. |
| `atlas_labels_file` | string | `` | - | Optional UTF-8 TSV with index and name columns for atlas_file. Indices are actual voxel values, without offsets. Empty generates Region_<index> names. Requires atlas_file. |
| `atlas_type` | string | `jhu` | - | Atlas used for parcellation. |
| `bundle_roi_dir` | string | `` | - | Optional directory of bundle masks to use as the parcellation INSTEAD of the atlas, e.g. data/ismrm2015/scoring_data_Renauld2023/ROI. Must contain all_masks/ (containment corridors, which become the labels) and may contain endpoints/ and any_masks/ (gates, addressable by name but not labels). The atlas parcellation is preserved as parcellation_mask_<atlas_type> rather than discarded. Empty = atlas only. Set this when ROIs should be addressed by the names a scorer uses: an atlas label and a bundle of the same name are not the same region. |
| `t1_available` | logical | `false` | - | A T1 volume is present alongside the DWI. |
| `use_t1_registration` | logical | `false` | - | Register T1 to DWI space. |
| `register_to_mni` | logical | `false` | - | Register to MNI space. |
| `dwell_time` | numeric | `0.00058` | - | Field-map (FUGUE) correction: effective echo spacing in SECONDS (time between phase-encode lines). Set it from the acquisition; the default is a generic guess. Displacement in voxels = field[Hz] x dwell_time x phase-encode lines, so a wrong value scales the correction. ISMRM 2015: 0.001. |
| `phase_encoding_dir` | string | `y` | - | Field-map (FUGUE) correction: phase-encoding axis and polarity (FUGUE --unwarpdir). ISMRM 2015: y (checked against GT: y- makes the distortion worse). |
| `fieldmap_units` | string | `auto` | - | Units of the field map file. auto = from the file name (fmap_Hz / RadPerSec). Hz maps are converted to rad/s, which is what FUGUE expects. |

## Migrating from the old flat config

Old configs still load; each superseded key produces a deprecation warning
naming its replacement.

| Old key | Replacement |
|---|---|
| `act_enabled` | `tractography.act` |
| `integrator` | `tractography.integrator.method` |
| `step_size` | `tractography.integrator.step` |
| `rkf_tolerance` | `tractography.integrator.tolerance` |
| `rkf_tol` | `tractography.integrator.tolerance` |
| `step_min` | `tractography.integrator.step_min` |
| `step_max` | `tractography.integrator.step_max` |
| `rkf_safety` | `tractography.integrator.safety` |
| `adaptive_step` | `tractography.integrator.adaptive` |
| `interp_method` | `tractography.interpolation.method` |
| `seed_density` | `tractography.seeding.density` |
| `seed_strategy` | `tractography.seeding.strategy` |
| `seed_fa_threshold` | `tractography.seeding.fa_min` |
| `termination_fa` | `tractography.termination.fa_min` |
| `min_length` | `tractography.termination.min_arc` |
| `csd_peaks_file` | `tractography.csd.peaks_file` |
| `csd_lmax` | `tractography.csd.lmax` |
| `csd_max_peaks` | `tractography.csd.max_peaks` |
| `csd_peak_thresh` | `tractography.csd.peak_thresh` |
| `csd_peak_min_sep` | `tractography.csd.peak_min_sep` |
| `enable_diagnostics` | `tractography.diagnostics` |
| `integration_order: 1\|2\|4\|5` | `integrator.method: euler\|rk2\|rk4\|rkf45` |
| `max_steps` | `termination.max_arc` (converted as `max_steps x step`) |

`integration_order` was a **method selector wearing a numeric-order name**:
the value `5` selected RKF45. That value was not wrong numerically - the
implementation uses Dormand-Prince coefficients and advances on the 5th-order
solution, using the embedded 4th-order one only for error control - but a
method selector should not be spelled as a number. `integrator.method` names
the method instead.

`max_steps` counted integration steps, so halving the step size silently
halved how far a track could travel. `termination.max_arc` is an arc length
in voxels and is step-invariant; `max_steps` is now derived as
`ceil(max_arc / step)`.

## Retired keys

These were accepted by the old loader but read by no tracker. They are
ignored with a warning rather than silently accepted.

| Key | Why it is gone |
|---|---|
| `gate_power` | Never read by any tracker. Defaulted and validated by the old loader but had no effect. |
| `crossing_cp` | Never read by any tracker. |
| `curv_beta` | Never read by any tracker. |
| `crossing_detect` | Never read by any tracker. |
| `swing_ratio_max` | Never read by any tracker. |
| `transport_gate` | Never read by any tracker. |
| `transport_strength` | Never read by any tracker. |
| `bishop_eps` | Referenced only in a comment; no tracker reads it as an option. |
| `frame_sel_power` | Removed with the alignment-weighted denoise it controlled. Same mechanism as the retired sel_power: a free exponent with no principled value. Measured against ground-truth curvature it also did nothing - correlation 0.217 to 0.239 across sel 0 to 64, while every setting suppressed curvature ~17x below the true 0.1538/vox. |
| `fa_threshold` | Functionally dead: printed by hinec/standard, and only a seed-mask fallback in the removed trackers that never fired because runTractography always supplies a seed mask. Use termination.fa_min to stop tracking, or seeding.fa_min to restrict seeding. NOTE: in 5 of the 12 original configs fa_threshold and termination_fa held DIFFERENT values, so they were never interchangeable. |
| `order` | Legacy backward-compatibility key, read by no tracker as an option. |
| `sel_power` | Removed. HINEC is now pure interpolation + integration. It biased interpolation toward the incoming direction (weight = alignment^sel_power), which has no justification for DTI (one eigenvector per voxel, nothing to disambiguate) and made the ODE direction-dependent. For CSD, nearest-peak selection is retained because it is structural, but the alignment exponent is gone. |
