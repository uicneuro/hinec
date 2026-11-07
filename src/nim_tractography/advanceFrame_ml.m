function [tracts, frames, info] = advanceFrame_ml(w, mf, dom, seeds, dir0, h, smax, varargin)
%ADVANCEFRAME_ML  advanceFrame with two additions, for MMF-ml and the history veto.
%
%   A copy of Sehun Chun's advanceFrame.m (commit 7fb4c75), which stays unchanged. The scheme, its
%   options and their defaults are the same; the copy adds two things, both off unless asked for:
%     - HEADING-AWARE FIELDS. w, mf and 'Confidence' may be handles of two arguments f(X, T), T the
%       e1 of the frame being evaluated, so that a multilayer field (MMF-ml, docs/MMF_ML.md) is read
%       on the layer of the tract's own heading. Handles of one argument are called as before.
%     - THE HISTORY VETO ('VetoK', 'VetoTol'), part of the hold scheme; see below.
%   With handles of one argument and 'VetoK' 0 it takes the same steps as advanceFrame
%   (TestAdvanceFrameMl). nim_tractography_mmf_ml calls it; nim_tractography_mmf calls advanceFrame.
%
%   The original description follows.
%
%ADVANCEFRAME  Moving frames tractography: transport of the full moving frame
%              with the connection 1-forms (de = [w] e along the fibre).
%
%   [tracts, frames, info] = advanceFrame(w, mf, dom, seeds, dir0, h, smax)
%
%   The fibre through a seed is followed as a curve x(s) (s = arc length)
%   that carries its moving frame E(s) = (e_1; e_2[; e_3]) (one frame
%   vector per row):
%
%       dx/ds = e_1,     dE/ds = Omega E,     Omega_ij = w_ij(e_1) = sum_k wij_k e_1k.
%
%   The frame is advanced by the connection 1-forms w_ij = < d e_i, e_j >
%   evaluated on the direction of travel (the change of the frame is
%   de_i = sum_j w_ij e_j), and the position by the first frame vector.
%   Omega is skew-symmetric with the three entries
%       kappa_1 = w_12(e_1),  kappa_2 = w_13(e_1)     curvature components,
%       tau     = w_23(e_1)                           twist of the normal plane,
%   |kappa| = sqrt(kappa_1^2 + kappa_2^2) is the curvature of the fibre.
%
%   Step: with Omega frozen over a step the pose (x, E) has the closed-form
%   solution (exponential of the pose, a rotation plus an arc/helix)
%       E_{n+1} = exp(h Omega) E_n,
%       x_{n+1} = x_n + sum_j [ int_0^h exp(s Omega) ds ]_1j e_j,
%   evaluated with the Rodrigues formulas. The frame therefore stays
%   exactly orthonormal and right-handed without re-orthogonalization, and
%   the step is EXACT wherever Omega is constant along the fibre (circular
%   arcs in 2-D, helices in 3-D), for any step size. 'Order' selects
%       1  Omega frozen at the start of the step (Lie-Euler),
%       2  exponential midpoint rule (Omega at the half step), 2nd order,
%       4  commutator-free 4th-order Lie group method CF4 (Celledoni,
%          Marthinsen & Owren 2003): four connection evaluations and two
%          exponentials per step, like RK4 in cost.
%   All of them keep the exactness for constant Omega.
%
%   Gauge: the connection components wij_k are those of the frame field of
%   the data, Et(x). The transported frame E differs from Et(x) by a
%   rotation G = E Et' (sign flips of the axial data, drift, a different
%   normal-plane orientation), and the frame field G Et has the connection
%   G w G' + dG G'. The step uses Omega = G Omegat G' with Omegat the
%   contraction of the data's wij_k with the direction of travel, so the
%   components are always applied in the transported frame's own
%   coordinates (for axial data this is the "flip if the dot product is
%   negative" of streamline tractography, extended to the whole frame).
%
%   Inputs
%       w      connection components wij_k = e_j . d(e_i)/dx_k of the data
%              frames: the nodal array [N+1 x N+1 x Nx x Ny x m x m x expdim]
%              (2-D domain) or [N+1 x N+1 x N+1 x Nx x Ny x Nz x m x m x
%              expdim] (3-D domain) of Computeconnection (evaluated along
%              the tract with interpGLL), or a function handle wq = w(X)
%              returning [nq x m x m x expdim] for query points X [nq x expdim]
%       mf     the data frames: the nodal array [N+1 x N+1 x Nx x Ny x m x
%              spacedim] or [N+1 x N+1 x N+1 x Nx x Ny x Nz x m x spacedim]
%              of diffvectomf (interpolated with interpGLL and
%              re-orthonormalized), or a function handle Eq = mf(X)
%              returning [nq x m x spacedim]. They give the frame at the
%              seeds and the gauge G along the tract. May be [] when
%              'Frame0' is given (then G = I is assumed).
%       dom    structure of TractDomain2D or TractDomain3D (needed for nodal
%              w, mf), or for function handles any structure with the
%              domain extent xlim, ylim [, zlim]; a tract stops when it
%              leaves this extent
%       seeds  [nseed x expdim] seed points
%       dir0   [1 x spacedim] or [nseed x spacedim] initial direction: the
%              sign of e_1 at the seed is chosen to continue it (e_2 is
%              flipped with e_1, so the frame stays right-handed)
%       h      step in arc length (the frame scheme has no stability limit;
%              h is set by accuracy and by the resolution of the data: keep
%              h*theta below 'MaxTurn', see below)
%       smax   maximum arc length followed from the seed; a tract stops
%              earlier when it leaves the domain (the last segment is cut
%              at the boundary), when the next step would turn the frame by
%              more than 'MaxTurn' (or when the data's e_1 at its end is
%              more than 'MaxTurn' from the heading at its start), or when
%              the connection is not finite
%
%   Options (name/value pairs)
%       'Order'   1, 2 or 4 (default 4), see above
%       'Mode'    'transport' (default): the frame is carried along by the
%                 connection only; the data enter through wij_k and the
%                 gauge, and the direction of travel is the transported e_1.
%                 'reanchor': after every step the frame is replaced by the
%                 data frame at the new point with the signs of its vectors
%                 chosen by the transported frame, i.e. the transport fixes
%                 the orientation and provides the curvature of the step,
%                 the data fix the direction (closest to streamline
%                 tractography; limits drift when the connection is noisy).
%                 'adaptive': transport, but the frame is re-anchored (as
%                 in 'reanchor') after a step at whose end the transported
%                 e_1 has drifted from the data's e_1 by 'ReanchorTol' or
%                 more; elsewhere the transported frame is kept. The
%                 transported heading is trusted while it agrees with the
%                 data and snapped back to them when it does not.
%       'ReanchorTol'
%                 drift angle in radians (default 1e-3) at or above which
%                 the 'adaptive' mode re-anchors: the angle between the
%                 transported e_1 and the data's e_1 (modulo sign) at the
%                 end of the step. 0 re-anchors after every step like
%                 'reanchor', Inf never, like 'transport'. Ignored by the
%                 other modes, and the default of 'ExitTol' of the hold
%                 scheme. The default is an order of magnitude below
%                 the heading error of a one-voxel CF4 step across a bend
%                 (~1e-2 rad) and above the step-to-step wobble of the
%                 interpolated test data at normal steps, so it fires on
%                 every one-voxel step across a sharp bend (80% of them on
%                 the interpolated data), on about 30% of the steps at
%                 h = 0.16 there, where it still lowers the tract error,
%                 and on at most a fifth of the steps elsewhere, without
%                 effect on the error there.
%       'Twist'   'data' (default): tau = w_23(e_1) from the data;
%                 'zero': tau is dropped, the normal plane (e_2, e_3) is
%                 transported without rotation about e_1 (rotation-
%                 minimizing / Bishop transport), which is the right choice
%                 where the data's normal-plane frame is arbitrary
%                 (lambda_2 ~ lambda_3). The fibre followed is the same in
%                 the limit h -> 0 (the curvature vector kappa_1 e_2 +
%                 kappa_2 e_3 is the same vector in space), but with tau
%                 dropped Omega is no longer constant along a helix in the
%                 transported frame, so the step is exact there only to the
%                 order of the integrator (exact helix of TractDomain3D,
%                 h = 0.1: 1e-2 for 'Order' 1, 8e-7 for 'Order' 4, against
%                 1e-16 with 'data'). Note
%                 that a re-anchoring ('reanchor', and 'adaptive' where it
%                 fires) replaces the whole frame, so the rotation-
%                 minimizing normal plane is reset to the data's normal
%                 plane at every re-anchored point and kept only between
%                 re-anchorings.
%       'Frame0'  explicit seed frames, [m x spacedim] (all seeds) or
%                 [nseed x m x spacedim], used instead of the data frame
%                 at the seed (dir0 is then ignored)
%       (w, mf and 'Confidence' may also be handles of two arguments, f(X, T), T [n x spacedim] the
%       heading of the frame being evaluated: a LIFTED field on positions and directions, read on the
%       sheet of the incoming heading where bundles cross; see nim_lift_connection.)
%       'Confidence'
%                 the MMF-4h hold scheme (default [], off): a handle
%                 c = cl(X) returning the confidence of the data at the rows
%                 of X [nq x expdim] (for DTI the linear shape measure
%                 CL = (lambda_1 - lambda_2)/lambda_1, or the anisotropy of
%                 a crossing's averaged tensor), or a nodal field of the
%                 layout [N+1 N+1 Nx Ny] / [N+1 N+1 N+1 Nx Ny Nz] of the
%                 domain. Where c < 'HoldMin' the data are degenerate (a
%                 crossing fitted as a pancake): at such a stage point the
%                 connection is NOT read from the data but HELD, Omega_hold
%                 = the running mean of the connection at the last 'HoldK'
%                 TRUSTED nodes (outside the zone, with rotation rate theta
%                 <= 'HoldCap'; a node of non-finite confidence and finite
%                 connection is one), taken in the carried frame - the only frame in
%                 which "the same curvature as before" is defined - with
%                 its twist Omega_23 set to zero ('HoldTwist' 'zero',
%                 default) or kept ('hold'), or with its curvature set to
%                 zero as well ('HoldCurvature' 'zero': a straight coast;
%                 default 'hold'). Omega being constant along the coast the
%                 exponential step is exact there: the pose follows the
%                 osculating circle of its last trusted pose (or the
%                 osculating helix with the twist held), through the
%                 degenerate zone, and nothing in the integrator changes.
%                 Without a trusted node yet the coast is straight. In the
%                 re-anchoring modes the frame is never re-anchored at a
%                 node with c < 'HoldMin', and the heading test of
%                 'MaxTurn' against the data's e_1 is not applied there (the
%                 data's axis means nothing in the zone). On leaving the
%                 zone (a node with c >= 'HoldMin' after one below) the
%                 frame is re-anchored only if its drift is at most
%                 'ExitTol'; a larger drift means the coast has left the
%                 bundle and the tract stops ('exit'; a data direction more
%                 than 'MaxTurn' from the heading at the start of the step
%                 is already a 'turn', and in 'transport' mode there is
%                 neither the re-anchoring nor the test); with 'ExitDelay' n > 0
%                 the coast goes on for n more steps beyond the zone's edge,
%                 on the held connection, and the drift is tested at the end
%                 of them, where the data's direction is read clear of the
%                 swing ahead of a degenerate zone (these nodes count as
%                 coasted in info.coasted, do not enter the ring and do not
%                 count towards 'CoastMax'; 2026-11-02). A coast of more than
%                 'CoastMax' consecutive steps started in the zone stops the
%                 tract ('coast'; default Inf). A non-finite confidence is
%                 no data: the connection is used as it is (NaN stops the
%                 tract as usual), during the steps of 'ExitDelay' as well.
%                 The ring is kept in the carried frame, so a re-anchoring
%                 with three frame vectors, which changes its normal plane,
%                 turns the ring with it. The nodes in the zone are reported in
%                 info.coasted, the confidence in info.confidence and the
%                 Omega used at every node in info.held.
%                 MMF-4s ('HoldCurvature' 'march'): instead of being held
%                 constant, the curvature kappa =
%                 Omega_12 and the transverse component sigma = w_12(e_2) =
%                 div e_1 of the data's connection (the fan-out of the
%                 fibres, which transports the axon density lambda by
%                 d ln lambda/ds = -sigma) are MARCHED through the zone,
%                   d sigma/ds = -sigma^2 + div(kappa e_2),  d kappa/ds = -kappa sigma,
%                 the first the flatness law of the planar connection, the
%                 second the closure div(kappa e_1) = 0, which scales the
%                 curvature with the fan-out, from their means over the
%                 ring of trusted nodes, anchored
%                 at the mean arc length of those nodes. The divergence
%                 term D = div(kappa e_2) = d kappa/d e_2 - kappa^2, which
%                 the zone cannot give, is zero ('HoldDivergence' 'zero',
%                 default) or held at its mean over the ring ('hold': at
%                 every trusted node D is measured by central differences
%                 of the data's curvature along e_2 over the step h, and
%                 left out where a probe falls into the zone). With
%                 D = 0 the ratio kappa/sigma is constant along the coast,
%                 kappa(s) = kappa_0/u and sigma(s) = sigma_0/u with u = 1
%                 + sigma_0 (s - s_0): the coast is a logarithmic spiral (a
%                 circle for sigma_0 = 0, a straight ray for kappa_0 = 0),
%                 exact for the fields of concentric circles, of
%                 logarithmic spirals and of rays through a point; with D
%                 held, u'' = D u (u a cosh or a cos). The connection so
%                 prescribed along the coast is evaluated at the stage
%                 points of the step, so the CF4 integrates it to fourth
%                 order. A coast whose u falls below 0.1 (converging fibres
%                 approaching the focus of their fan, where the marched
%                 curvature would exceed ten times the trusted one) stops
%                 the tract ('focus'). The transverse component used at
%                 every node - the data's at a trusted node, the held mean
%                 or the marched value in the zone - is returned in
%                 info.spread (for every hold mode).
%                 With three frame vectors (m = 3; 2026-11-02) the law is
%                 the matrix one. In the carried frame, rotation-minimizing
%                 along the coast, the curvature vector a = (Omega_12,
%                 Omega_13) and the transverse block B_pq = e_p . d e_1/d e_q
%                 (p, q = 2, 3) of the data's connection obey
%                   dB/ds = -B^2 + D,   da/ds = -B' a,
%                 the first the flatness of the connection (the Riccati
%                 equation of the congruence; its trace is Raychaudhuri's),
%                 the second the closure d(a^flat) = 0 that div(kappa e_1) =
%                 0 is for a planar field. D is taken zero, so that with
%                 U = I + (s - s_0) B_0, B = B_0 U^-1 and a = U^-T a_0: the
%                 fan-out across the curvature (a fibre bundle turned about
%                 an axis) leaves the curvature alone, and a helix's B_0
%                 turns a by its torsion to first order. u = det U, the
%                 inverse of the density the coast carries, is 1 + sigma_0
%                 (s - s_0) on a planar field and its u < 0.1 stops the
%                 tract ('focus'); info.spread is the trace of B. The law
%                 holds in a rotation-minimizing frame, so m = 3 needs 'Twist'
%                 'zero'; 'HoldTwist' and 'HoldDivergence' 'hold' are not
%                 defined for it.
%       'HoldMin' confidence below which the data are degenerate (default 0.2)
%       'HoldK'   length of the running mean of trusted connections (default 3)
%       'HoldCap' largest rotation rate theta of a trusted connection that
%                 enters the mean (default Inf)
%       'CoastMax'
%                 largest number of consecutive steps started in the zone
%                 (default Inf)
%       'ExitTol' largest drift, in radians, of a tract leaving the zone
%                 (default 'ReanchorTol'). Trusted tracking and this test
%                 want different things of noisy data: a tolerance small
%                 enough to keep to the data, and one above the disagreement
%                 of neighbouring voxels.
%       'ExitDelay'
%                 steps of coast beyond the zone's edge before the exit
%                 test (default 0: the test at the first trusted node; in
%                 'transport' mode the coast is extended the same way and
%                 no test is applied)
%       'VetoK'   the history veto (default 0: off; needs the hold scheme,
%                 mf and a re-anchoring mode). A crossing in which the
%                 tract's bundle is the minor one is fitted as a linear
%                 tensor along the other bundle: its confidence is high and
%                 the data's e_1 is wrong, and re-anchoring follows it a
%                 few degrees per step. The veto compares the data's e_1 at
%                 the end of a trusted step with the heading e_1 at the node
%                 VetoK steps back, d = arccos|e_1(data, x_n) . e_1(n - VetoK)|
%                 (not before the end of the tract's last vetoed coast).
%                 Where d > 'VetoTol' the step is not taken: the tract
%                 returns to that node, with its frame and its ring as they
%                 were there, and coasts from it as in the zone (the steps
%                 count towards 'CoastMax'). The vetoed coast goes on at
%                 least to the node where the veto fired (before it the data
%                 agree with the heading they are about to leave, and an end
%                 there would re-anchor onto the turn), and then while the
%                 data's e_1 at the end of a step is more than 'VetoTol'
%                 from the coast's heading; it leaves the zone through
%                 'ExitTol'. info.vetoes counts the returns per tract.
%                 (2026-11-06)
%       'VetoTol' largest disagreement d of the veto, radians (default pi/6)
%       'BundleRule'
%                 labels of the seeds' bundles, [nseed x 1] (default []):
%                 the rule that the fibres of one bundle never cross. Two
%                 coasts do not know of each other, so a step on which two
%                 tracts of the same bundle, at least one of them coasting,
%                 would cross stops both ('cross'). They cross when the
%                 offset between them, in the plane normal to their mean
%                 heading, reverses or vanishes over the step (in the
%                 plane, when their order across the bundle swaps). Tracts
%                 of different bundles may cross. (2026-11-03)
%       'BundleLaw'
%                 labels of the seeds' bundles, [nseed x 1] (default []):
%                 the law behind that rule, for the march of MMF-4s with
%                 the planar connection. The fibres of one bundle are the
%                 integral curves of one field, whose connection is flat:
%                 the divergence term D = d kappa/d e_2 - kappa^2 of the
%                 march, which the zone cannot give and which 'HoldDiver-
%                 gence' takes as zero or from the data at the ring, is
%                 taken at every step from the bundle itself, d kappa/d e_2
%                 from the curvature of the neighbouring tracts of the
%                 bundle (the nearest on either side of the tract along its
%                 e_2, the coast's curvature where they coast and the
%                 data's where they are trusted, the part of their
%                 difference due to the offset along the fibre removed by
%                 the march's own d kappa/ds = -kappa sigma; one-sided
%                 with one neighbour, zero with none). The coasts of a
%                 bundle then march as one sheet, each step of the march
%                 taking the law u'' = D u from the state (kappa, sigma)
%                 at its start with D frozen over the step, and 'focus'
%                 counting u from the entry of the zone. Needs
%                 'HoldCurvature' 'march', m = 2 and 'HoldDivergence'
%                 'zero'; may be combined with 'BundleRule'. (2026-11-04)
%       'HoldTwist', 'HoldCurvature'
%                 'zero' | 'hold' for the twist (default 'zero') and
%                 'hold' | 'zero' | 'march' for the curvature (default
%                 'hold') of the held connection; 'march' is the MMF-4s
%                 scheme described under 'Confidence'
%       'HoldDivergence'
%                 'zero' (default) | 'hold': the divergence term of the
%                 march of MMF-4s, zero or held at its mean over the
%                 trusted nodes of the ring (needs the data frames mf;
%                 m = 2 only)
%       'MaxTurn' largest rotation of the frame allowed in one step, in
%                 radians (default 0.5): a tract stops, at its current
%                 point, before a step in which h*theta exceeds it at any
%                 stage, theta = sqrt(kappa_1^2 + kappa_2^2 + tau^2) being
%                 the rotation rate of the frame (|Omega| in the sense
%                 theta^2 = -tr(Omega^2)/2). The exponential step is exact
%                 for any h*theta, but beyond ~0.5 radian the step no longer
%                 resolves the turning of the fibre (the data are sampled
%                 once or four times over a large fraction of a turn, at
%                 2*pi the arc closes on itself), and a large theta usually
%                 means the data are singular there (crossing, isotropic
%                 voxel). This is the frame counterpart of the angle
%                 threshold of streamline tractography. Inf disables it.
%                 'MaxTurn' also applies to the heading as the data see
%                 it: a step at whose end the data's e_1 is more than
%                 MaxTurn (modulo sign, so at most pi/2) from the heading
%                 at its start is not taken either, with the same reason
%                 'turn'. In the re-anchoring modes that is the turn the
%                 re-anchoring would realise, the counterpart of the angle
%                 test between consecutive directions of streamline
%                 tracking, which a snap to the data could otherwise
%                 bypass; in 'transport' mode it stops a tract whose
%                 transported heading the data contradict beyond
%                 'MaxTurn' (not evaluated when mf is []). A tract also stops
%                 before a step whose connection is not finite at some
%                 stage.
%
%   Outputs
%       tracts  cell {nseed x 1} of [npts x spacedim] tract points, the
%               seed first (the position has the components of the frame
%               vectors; the first expdim are the domain coordinates)
%       frames  cell {nseed x 1} of [npts x m x spacedim] transported
%               frames at the tract points (at a clipped boundary point:
%               the frame at the uncut end of the last step, beyond the
%               boundary, where the data are extrapolated)
%       info    structure:
%               stop     cell {nseed x 1} of char, why each tract ended:
%                        'boundary' (left the domain, last point is the cut
%                        on the boundary), 'length' (reached smax), 'turn'
%                        (the next step would turn by more than MaxTurn,
%                        or the data's e_1 at its end is more than MaxTurn
%                        from the heading at its start), 'coast', 'exit',
%                        'focus' and 'cross' (the hold scheme: coast cap, the zone
%                        left with a drift above ExitTol, the march
%                        of MMF-4s reaching the focus of its fan, or a
%                        step crossing a tract of the same 'BundleRule'),
%                        'nan' (the connection is not finite at the next
%                        step); for 'turn' and 'nan' the last point is the
%                        last admissible integration node
%               clipped  nseed x 1 logical, true when stop is 'boundary'
%               s        cell, arc length at the points
%               kappa, tau, turn
%                        cell, curvature, twist and rotation per step
%                        h*theta at the points (NaN at a clipped point)
%               drift    cell, angle between the transported e_1 arriving
%                        at a point and the data's e_1 there (modulo sign),
%                        before any re-anchoring; 0 at the seed; NaN where
%                        there is no data frame (mf = [])
%               reanchored
%                        cell, logical per point: the frame stored at the
%                        point is the data frame (re-anchored) rather than
%                        the transported one; false at the seed. Like the
%                        frame, drift and reanchored at a clipped boundary
%                        point refer to the uncut end of the last step.
%               coasted, confidence, held
%                        the hold scheme (see 'Confidence'): per point,
%                        whether the node lies in the zone, the confidence
%                        there (NaN without the scheme) and the Omega used
%                        at the node [npts x m x m] (NaN at a clipped point,
%                        empty without the scheme)
%               spread   cell, the transverse component sigma = w_12(e_2)
%                        [+ w_13(e_3)] = div e_1 of the data's connection
%                        at the points along the direction of travel (the
%                        data's frame aligned with the carried one; the
%                        carried e_2 when mf is []): the data's value at a
%                        trusted node, the held mean of the ring or the
%                        marched value (MMF-4s) at a node in the zone; NaN
%                        at a clipped point and without the hold scheme
%               hold     the hold options (on = false without 'Confidence')
%               order, mode, twist, maxTurn, reanchorTol (the options given)
%
%   Examples
%       dom = TractDomain2D('arc');
%       mf  = diffvectomf(dom.gll.Vx, dom.gll.Vy);        % 2-D frames
%       w   = Computeconnection(mf, dom);
%       [tr, fr, info] = advanceFrame(w, mf, dom, [0.25 0; 0.5 0], [0 1], 0.01, 4);
%       plot(tr{1}(:,1), tr{1}(:,2))
%       [tr1] = advanceFrame(w, mf, dom, [0.25 0], [0 1], 0.1, 4, 'Order', 1);
%       [tr4] = advanceFrame(w, mf, dom, [0.25 0], [0 1], 0.1, 4, 'Mode', 'reanchor');
%       [trA, ~, infA] = advanceFrame(w, mf, dom, [0.25 0], [0 1], 0.1, 4, 'Mode', 'adaptive', 'ReanchorTol', 1e-3);
%       nnz(infA.reanchored{1})                           % steps at which the frame was re-anchored
%       [tr, fr, info] = advanceFrame(w, mf, dom, [0.75 0], [0 1], 0.25, 4, 'MaxTurn', Inf);
%       info.stop{1}                                      % 'boundary', 'length', 'turn' or 'nan'
%
%   Requires interpGLL.m for nodal w and mf (and TractDomain2D.m or
%   TractDomain3D.m, diffvectomf.m, Computeconnection.m to produce them).
%
%   HINEC / moving-frames tractography test domain. Created 2026-10-09;
%   'adaptive' mode and 'ReanchorTol' added 2026-10-13. The linear-curvature
%   step ('Model', 'linear', 2026-10-14) was removed on 2026-10-19; the version
%   that had it is archived as MMF4vca_old.m. Since 2026-10-20 'MaxTurn' also
%   bounds the turn of the heading as the data see it (the data's e_1 at the
%   end of a step against the heading at its start), in every mode, so that a
%   re-anchoring cannot exceed 'MaxTurn'; drift is evaluated in every mode.
%   The MMF-4h hold scheme ('Confidence' and the 'Hold*' options, 2026-10-26)
%   coasts through degenerate data on the held connection of the last
%   trusted nodes; MMF-4s ('HoldCurvature' 'march', 2026-10-30) marches the
%   curvature and the transverse component through the zone by the flatness
%   law of the planar connection instead, and info.spread was added; with
%   three frame vectors by the matrix law of the flat connection since
%   2026-11-02.

%% --- options --------------------------------------------------------------
order   = 4;
mode    = 'transport';
twist   = 'data';
frame0  = [];
maxTurn = 0.5;
reanchorTol = 1e-3;
exitTol    = [];             % hold scheme: the test on leaving the zone ([]: reanchorTol)
confidence = [];             % the hold scheme: confidence of the data (handle or nodal field)
holdMin    = 0.2;
holdK      = 3;
holdCap    = Inf;
coastMax   = Inf;
exitDelay  = 0;              % hold scheme: steps of coast beyond the zone's edge before the exit test
bundleRule = [];             % labels of the seeds' bundles: the tracts of one bundle never cross (the rule)
bundleLaw  = [];             % labels of the seeds' bundles: the march's divergence term from the bundle (the law)
holdTwist  = 'zero';
holdCurv   = 'hold';
holdDiv    = 'zero';         % MMF-4s: the divergence term of the march
vetoK      = 0;              % hold scheme: the history veto's look-back in steps (0: off)
vetoTol    = pi / 6;         % hold scheme: the history veto's largest disagreement, radians
if mod(numel(varargin), 2) ~= 0
    error('advanceFrame:badOptions', 'Options must come in name/value pairs.');
end
for n = 1:2:numel(varargin)
    name = varargin{n};
    val  = varargin{n+1};
    if ~ischar(name) && ~(isstring(name) && isscalar(name))
        error('advanceFrame:badOptions', 'Option names must be character vectors.');
    end
    switch lower(char(name))
        case 'order'
            if ~isnumeric(val) || ~isscalar(val) || ~any(val == [1 2 4])
                error('advanceFrame:badOptions', '''Order'' must be 1, 2 or 4.');
            end
            order = double(val);
        case 'mode'
            mode = validatestring(val, {'transport', 'reanchor', 'adaptive'}, 'advanceFrame', '''Mode''');
        case 'reanchortol'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnegative', 'nonnan'}, 'advanceFrame', '''ReanchorTol''');
            reanchorTol = double(val);
        case 'twist'
            twist = validatestring(val, {'data', 'zero'}, 'advanceFrame', '''Twist''');
        case 'frame0'
            frame0 = double(val);
        case 'maxturn'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive'}, 'advanceFrame', '''MaxTurn''');
            maxTurn = double(val);
        case 'confidence'
            if ~(isempty(val) || isa(val, 'function_handle') || isnumeric(val))
                error('advanceFrame:badOptions', '''Confidence'' must be a function handle, a nodal field or [].');
            end
            confidence = val;
        case 'holdmin'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnan'}, 'advanceFrame', '''HoldMin''');
            holdMin = double(val);
        case 'holdk'
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'positive'}, 'advanceFrame', '''HoldK''');
            holdK = double(val);
        case 'holdcap'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive', 'nonnan'}, 'advanceFrame', '''HoldCap''');
            holdCap = double(val);
        case 'coastmax'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnegative', 'nonnan'}, 'advanceFrame', '''CoastMax''');
            coastMax = double(val);
        case 'exittol'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'nonnegative', 'nonnan'}, 'advanceFrame', '''ExitTol''');
            exitTol = double(val);
        case 'bundlerule'
            bundleRule = double(val(:));
        case 'bundlelaw'
            bundleLaw = double(val(:));
        case 'exitdelay'
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, 'advanceFrame', '''ExitDelay''');
            exitDelay = double(val);
        case 'holdtwist'
            holdTwist = validatestring(val, {'zero', 'hold'}, 'advanceFrame', '''HoldTwist''');
        case 'holdcurvature'
            holdCurv = validatestring(val, {'hold', 'zero', 'march'}, 'advanceFrame', '''HoldCurvature''');
        case 'holddivergence'
            holdDiv = validatestring(val, {'zero', 'hold'}, 'advanceFrame', '''HoldDivergence''');
        case 'vetok'
            validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, 'advanceFrame', '''VetoK''');
            vetoK = double(val);
        case 'vetotol'
            validateattributes(val, {'numeric'}, {'scalar', 'real', 'positive', 'nonnan'}, 'advanceFrame', '''VetoTol''');
            vetoTol = double(val);
        otherwise
            error('advanceFrame:badOptions', 'Unknown option ''%s''.', char(name));
    end
end
validateattributes(seeds, {'numeric'}, {'2d', 'real', 'finite', 'nonempty'}, 'advanceFrame', 'seeds');
validateattributes(h,     {'numeric'}, {'scalar', 'real', 'positive', 'finite'}, 'advanceFrame', 'h');
validateattributes(smax,  {'numeric'}, {'scalar', 'real', 'positive', 'finite'}, 'advanceFrame', 'smax');
seeds  = double(seeds);
nseed  = size(seeds, 1);
expdim = size(seeds, 2);

%% --- the data: connection and frames as functions of the position ----------
if isa(w, 'function_handle') && nargin(w) == 2        % a lifted field w(X, T): read on the sheet of the heading T
    evalW = w;
elseif isa(w, 'function_handle')
    evalW = @(X, T) w(X);
elseif isnumeric(w) && any(ndims(w) == [7 9])
    gdim = (ndims(w) - 3) / 2;                   % 2: [N+1 N+1 Nx Ny m m 2], 3: [N+1 N+1 N+1 Nx Ny Nz m m 3]
    if size(w, ndims(w)) ~= expdim || gdim ~= expdim
        error('advanceFrame:badConnection', ...
            'w is the nodal connection of a %d-D domain with %d differentiation directions, but the seeds have %d coordinates.', ...
            gdim, size(w, ndims(w)), expdim);
    end
    if gdim == 2
        evalW = @(X, T) interpGLL(w, dom, X(:,1), X(:,2));
    else
        evalW = @(X, T) interpGLL(w, dom, X(:,1), X(:,2), X(:,3));
    end
else
    error('advanceFrame:badConnection', ...
        ['w must be the nodal connection [N+1 N+1 Nx Ny m m 2] or [N+1 N+1 N+1 Nx Ny Nz m m 3] ' ...
         'or a function handle.']);
end
if isa(mf, 'function_handle') && nargin(mf) == 2      % lifted frames mf(X, T)
    evalMF = @(X, T) orthonormalize(mf(X, T));
elseif isa(mf, 'function_handle')
    evalMF = @(X, T) orthonormalize(mf(X));
elseif isnumeric(mf) && any(ndims(mf) == [6 8])
    gdim = (ndims(mf) - 2) / 2;                  % 2: [N+1 N+1 Nx Ny m spacedim], 3: [N+1 N+1 N+1 Nx Ny Nz m spacedim]
    if gdim ~= expdim
        error('advanceFrame:badFrames', ...
            'mf is the nodal frame field of a %d-D domain, but the seeds have %d coordinates.', gdim, expdim);
    end
    if gdim == 2
        evalMF = @(X, T) orthonormalize(interpGLL(mf, dom, X(:,1), X(:,2)));
    else
        evalMF = @(X, T) orthonormalize(interpGLL(mf, dom, X(:,1), X(:,2), X(:,3)));
    end
elseif isempty(mf)
    evalMF = [];
else
    error('advanceFrame:badFrames', ...
        ['mf must be the nodal frames [N+1 N+1 Nx Ny m spacedim] or [N+1 N+1 N+1 Nx Ny Nz m spacedim], ' ...
         'a function handle, or [].']);
end
if isempty(evalMF) && isempty(frame0)
    error('advanceFrame:noFrames', 'Give the data frames mf or the seed frames ''Frame0''.');
end
if isa(confidence, 'function_handle') && nargin(confidence) == 2   % a heading-aware confidence c(X, T)
    evalCL = confidence;
elseif isa(confidence, 'function_handle')
    evalCL = @(X, T) confidence(X);
elseif isnumeric(confidence) && ~isempty(confidence)
    gdim = ndims(confidence) / 2;                % 2: [N+1 N+1 Nx Ny], 3: [N+1 N+1 N+1 Nx Ny Nz]
    if gdim ~= expdim
        error('advanceFrame:badOptions', ...
            '''Confidence'' is a nodal field of a %d-D domain, but the seeds have %d coordinates.', gdim, expdim);
    end
    if gdim == 2
        evalCL = @(X, T) interpGLL(confidence, dom, X(:,1), X(:,2));
    else
        evalCL = @(X, T) interpGLL(confidence, dom, X(:,1), X(:,2), X(:,3));
    end
else
    evalCL = [];
end
useHold = ~isempty(evalCL);
if any(strcmp(mode, {'reanchor', 'adaptive'})) && isempty(evalMF)
    error('advanceFrame:noFrames', '''Mode'', ''%s'' needs the data frames mf.', mode);
end
lims = [dom.xlim(:).'; dom.ylim(:).'];
if isfield(dom, 'zlim') && ~isempty(dom.zlim)
    lims = [lims; dom.zlim(:).'];
end
if size(lims, 1) < expdim
    error('advanceFrame:badDomain', 'dom must give the extent of all %d coordinates.', expdim);
end
lims = lims(1:expdim, :);

%% --- seed frames -----------------------------------------------------------
if isempty(frame0)
    E = evalMF(seeds, repmat(double(dir0), nseed / size(dir0, 1), 1));   % [nseed x m x spacedim]
    m        = size(E, 2);
    spacedim = size(E, 3);
    validateattributes(dir0, {'numeric'}, {'2d', 'real', 'finite'}, 'advanceFrame', 'dir0');
    if size(dir0, 2) ~= spacedim || ~any(size(dir0, 1) == [1 nseed])
        error('advanceFrame:badDirection', ...
            'dir0 must be [1 x %d] or [%d x %d].', spacedim, nseed, spacedim);
    end
    if size(dir0, 1) == 1
        dir0 = dir0(ones(nseed, 1), :);
    end
    sgn = sign(sum(reshape(E(:,1,:), nseed, spacedim) .* double(dir0), 2));
    sgn(sgn == 0) = 1;
    E(:,1,:) = sgn .* E(:,1,:);                  % flip e_1 and e_2 together:
    E(:,2,:) = sgn .* E(:,2,:);                  % the frame stays right-handed
else
    if ndims(frame0) == 2
        frame0 = reshape(frame0, [1 size(frame0)]);
        frame0 = frame0(ones(nseed, 1), :, :);
    end
    if size(frame0, 1) ~= nseed
        error('advanceFrame:badFrame0', '''Frame0'' must be [m x spacedim] or [nseed x m x spacedim].');
    end
    E        = orthonormalize(frame0);
    m        = size(E, 2);
    spacedim = size(E, 3);
end
if spacedim < expdim
    error('advanceFrame:badFrames', ...
        'The frame vectors have %d components but the seeds %d coordinates.', spacedim, expdim);
end
if m < 2 || m > 3 || m > spacedim
    error('advanceFrame:badFrames', 'There must be 2 or 3 frame vectors (m = %d, spacedim = %d).', m, spacedim);
end
dropTwist = strcmp(twist, 'zero') && m == 3;
useGauge  = ~isempty(evalMF);
omegaFcn  = @(X, E) connectionMatrix(X, E, evalW, evalMF, useGauge, expdim, m, dropTwist);
marchHold = useHold && strcmp(holdCurv, 'march');    % MMF-4s
if marchHold && m == 3 && (~dropTwist || strcmp(holdTwist, 'hold') || strcmp(holdDiv, 'hold'))
    error('advanceFrame:badOptions', ['''HoldCurvature'' ''march'' (MMF-4s) with 3 frame vectors needs a rotation-minimizing ' ...
        'frame and the divergence term zero: use ''Twist'' ''zero''; ''HoldTwist'' and ''HoldDivergence'' ''hold'' are not defined for it.']);
end
if isempty(exitTol)
    exitTol = reanchorTol;
end
holdDivOn = marchHold && strcmp(holdDiv, 'hold');
vetoOn = vetoK > 0;
if vetoOn && (~useHold || isempty(evalMF) || strcmp(mode, 'transport'))
    error('advanceFrame:badOptions', '''VetoK'' needs the hold scheme (''Confidence''), the data frames mf and a re-anchoring mode.');
end
if holdDivOn && isempty(evalMF)
    error('advanceFrame:badOptions', '''HoldDivergence'' ''hold'' needs the data frames mf.');
end
uMin     = 0.1;                                       % MMF-4s: the coast stops ('focus') where u falls below this (kappa_0/kappa for m = 2, det U for m = 3)
sigmaFcn = @(X, E) transverseOf(X, E, evalW, evalMF, expdim, m, h, holdDivOn, evalCL, holdMin);

%% --- march all seeds together ------------------------------------------------
X = [seeds, zeros(nseed, spacedim - expdim)];
tracts  = cell(nseed, 1);
frames  = cell(nseed, 1);
kappa   = cell(nseed, 1);
tau     = cell(nseed, 1);
turn    = cell(nseed, 1);
drift   = cell(nseed, 1);                   % heading drift of the transported frame per point
reanch  = cell(nseed, 1);                   % whether the stored frame was re-anchored to the data
stop    = repmat({'length'}, nseed, 1);
coasted = cell(nseed, 1);                   % hold scheme: the node is in the degenerate zone
confid  = cell(nseed, 1);                   % hold scheme: confidence at the points
held    = cell(nseed, 1);                   % hold scheme: the Omega used at the node [npts x m x m]
hbuf    = NaN(nseed, holdK, m, m);          % hold scheme: the last holdK trusted connections (ring)
sbuf    = NaN(nseed, holdK, m - 1, m - 1);  % hold scheme: the transverse block of the connection at the ring's nodes (m = 2: sigma)
dbuf    = NaN(nseed, holdK);                % MMF-4s: the divergence term D = div(kappa e_2) at the ring's nodes
abuf    = NaN(nseed, holdK);                % hold scheme: the arc length of the ring's nodes
spread  = cell(nseed, 1);                   % the transverse component used at the nodes
hcount  = zeros(nseed, 1);
coast   = zeros(nseed, 1);                  % consecutive steps started in the zone
grace   = zeros(nseed, 1);                  % steps of coast left beyond the zone's edge ('ExitDelay')
if ~isempty(bundleRule)
    validateattributes(bundleRule, {'numeric'}, {'vector', 'numel', nseed, 'real', 'nonnan'}, 'advanceFrame', '''BundleRule''');
end
useLaw = ~isempty(bundleLaw);
nn      = ones(nseed, 1);                   % nodes of each tract so far (the history veto rolls tracts back)
vco     = false(nseed, 1);                  % the history veto: the tract coasts since a rollback
vfloor  = ones(nseed, 1);                   % the history veto: the earliest node a rollback may return to
vcount  = zeros(nseed, 1);                  % the history veto: rollbacks per tract
vneed   = zeros(nseed, 1);                  % the history veto: the node a vetoed coast must reach before it may end
if vetoOn                                   % the ring at each of the last vetoK + 1 nodes, to restore on a rollback
    rhH = NaN(nseed * (vetoK + 1), holdK * m * m);  rhS = NaN(nseed * (vetoK + 1), holdK * (m - 1) ^ 2);
    rhD = NaN(nseed * (vetoK + 1), holdK);          rhA = rhD;  rhC = zeros(nseed * (vetoK + 1), 1);
end
if useLaw && vetoOn
    error('advanceFrame:badOptions', '''VetoK'' is not defined with ''BundleLaw''.');
end
if useLaw
    validateattributes(bundleLaw, {'numeric'}, {'vector', 'numel', nseed, 'real', 'nonnan'}, 'advanceFrame', '''BundleLaw''');
    if ~marchHold || m ~= 2 || holdDivOn
        error('advanceFrame:badOptions', ...
            '''BundleLaw'' needs ''HoldCurvature'' ''march'' with the planar connection (m = 2) and ''HoldDivergence'' ''zero''.');
    end
end
lawOm = zeros(nseed, m, m);                 % the law: the coast's connection, transverse component and u at the current node
lawSg = zeros(nseed, 1);
lawU  = ones(nseed, 1);
lawOn = false(nseed, 1);                    % the tract coasts under the law since its last trusted node
if useHold
    c0 = evalCL(X(:, 1:expdim), reshape(E(:, 1, :), nseed, []));
else
    c0 = NaN(nseed, 1);
end
for k = 1:nseed
    tracts{k} = X(k, :);
    frames{k} = E(k, :, :);
    drift{k}  = 0;
    reanch{k} = false;
    confid{k} = c0(k);
    coasted{k} = useHold && isfinite(c0(k)) && c0(k) < holdMin;
end
active = true(nseed, 1);
nstep  = ceil(smax / h);
for n = 1:nstep
    idx = find(active);
    if isempty(idx)
        break
    end
    Xa = X(idx, :);
    Ea = E(idx, :, :);
    fo = false(numel(idx), 1);                  % MMF-4s: rows whose coast reaches the focus of its fan
    cx = false(numel(idx), 1);                  % rows whose step crosses a tract of the same bundle ('BundleRule')
    if useHold
        sA = h * (nn(idx) - 1);                 % arc length of the step's start node
        [sgM, dvM, blM] = sigmaFcn(Xa, Ea);     % the data's transverse component (and block, divergence term) at the start nodes
        [Om0, bl0, dv0, sb0, sg0] = ringMeans(hbuf(idx, :, :, :), sbuf(idx, :, :, :), dbuf(idx, :), abuf(idx, :), hcount(idx), m, holdTwist, holdCurv);
        forcedS = grace(idx) > 0 | vco(idx);    % rows still coasting beyond the zone's edge ('ExitDelay') or vetoed
        if marchHold && useLaw                  % the law: the march from the state at the start node, D from the bundle
            [OmS, sgS, uS] = marchedOmega(Om0, bl0, dv0, sA - sb0, uMin);    % where a coast starts: the plain march's state
            sgS = reshape(sgS, [], 1);  uS = reshape(uS, [], 1);
            on = lawOn(idx);
            OmS(on, :, :) = lawOm(idx(on), :, :);  sgS(on) = lawSg(idx(on));  uS(on) = lawU(idx(on));
            cS  = evalCL(Xa(:, 1:expdim), reshape(Ea(:, 1, :), numel(idx), []));
            zS0 = isfinite(cS(:)) & (cS(:) < holdMin | forcedS);
            dvB = bundleDivergence(Xa, Ea, OmS, sgS, zS0, bundleLaw(idx), omegaFcn, expdim);
            holdOm = @(c) OmS ./ max(marchLaw(sgS, dvB, c * h * ones(numel(idx), 1)), uMin);
            [uh, uph] = marchLaw(sgS, dvB, h * ones(numel(idx), 1));
            uE = uS .* uh;                      % u counted from the entry of the zone
            lawNext = struct('Om', OmS ./ max(uh, uMin), 'sg', uph ./ uh, 'u', uE);
        elseif marchHold                        % the connection prescribed along the coast, stage by stage
            holdOm = @(c) marchedOmega(Om0, bl0, dv0, sA + c * h - sb0, uMin);
            [~, sgS] = marchedOmega(Om0, bl0, dv0, sA - sb0, uMin);           % sigma marched to the start node
            [~, ~, uE] = marchedOmega(Om0, bl0, dv0, sA + h - sb0, uMin);     % u at the end node
        else
            holdOm = Om0;  sgS = sg0;  uE = ones(numel(idx), 1);
        end
        [Xn, En, Om1, maxRot, zoneS] = lieStep(Xa, Ea, h, order, omegaFcn, evalCL, holdMin, holdOm, expdim, forcedS);
        zoneD = zoneS & (~forcedS | vco(idx));  % the start node is degenerate by the confidence or vetoed (counts towards 'CoastMax')
        cE = evalCL(Xn(:, 1:expdim), reshape(En(:, 1, :), numel(idx), []));   % confidence at the end of the step
        zoneC = isfinite(cE) & cE < holdMin;    % the end node is degenerate by the confidence
        g = grace(idx);                         % 'ExitDelay': the coast goes on beyond the zone's edge
        g(zoneS & ~zoneC & g == 0 & ~vco(idx)) = exitDelay + 1;      % the edge is reached: exitDelay more steps of coast
        g = max(g - 1, 0);                      % one of them spent
        g(zoneC) = 0;                           % back in the degenerate zone
        grace(idx) = g;
        zoneE = zoneC | g > 0;                  % the end node still coasts
        sgS(~zoneS) = sgM(~zoneS);              % the transverse component used at the start node: the data's where trusted
        fo = marchHold & zoneE & uE < uMin;
        if ~isempty(bundleRule)
            cx = bundleCross(Xa, Xn, Ea, zoneS, bundleRule(idx), expdim);
        end
    else
        [Xn, En, Om1, maxRot] = lieStep(Xa, Ea, h, order, omegaFcn);
        zoneS = false(numel(idx), 1);  zoneE = zoneS;  zoneD = zoneS;  cE = NaN(numel(idx), 1);
        sgS = NaN(numel(idx), 1);
    end
    Eb = En;                                    % the transported frames, before any re-anchoring
    dr = NaN(numel(idx), 1);                    % drift of the transported heading at Xn
    dt = NaN(numel(idx), 1);                    % turn of the heading over the step as the data see it
    re = false(numel(idx), 1);                  % rows whose frame is re-anchored at Xn
    ex = false(numel(idx), 1);                  % rows leaving the zone with a drift above the tolerance
    vfire = false(numel(idx), 1);  vend = false(numel(idx), 1);  jv = nn(idx);   % the history veto (see below)
    if ~isempty(evalMF)
        Et = evalMF(Xn(:, 1:expdim), reshape(En(:, 1, :), numel(idx), []));   % the data frame at the end of the step
        dr = headingAngle(En, Et);              % transported e_1 at Xn against the data's e_1 there
        dt = headingAngle(Ea, Et);              % e_1 at the start of the step against the data's e_1 at its end
        if vetoOn
            % The history veto. A vetoed coast goes on while the data's e_1 at its end is more than vetoTol
            % from the coast's heading (or not finite); then it leaves the zone through the exit test.
            vrow = vco(idx);
            vcont = vrow & (nn(idx) + 1 < vneed(idx) | ~(dr <= vetoTol));   % not before the node where the veto fired
            vend  = vrow & ~vcont;                  % vetoed coasts that end here
            zoneE = zoneE | vcont;
            fo = marchHold & zoneE & uE < uMin;
            % A trusted step whose end the data see more than vetoTol from the heading vetoK steps back
            % (not before the end of the last vetoed coast) is not taken: the tract returns to that node
            % and coasts from it on the ring it had there.
            cand = ~zoneS & ~forcedS & ~zoneE;
            jv = max(nn(idx) + 1 - vetoK, vfloor(idx));
            vfire = false(numel(idx), 1);
            if any(cand)
                rows = find(cand);
                Ej = zeros(numel(rows), m, spacedim);
                for q = 1:numel(rows)
                    Ej(q, :, :) = frames{idx(rows(q))}(jv(rows(q)), :, :);
                end
                vfire(rows) = headingAngle(Ej, Et(rows, :, :)) > vetoTol;
            end
        end
        switch mode
            case 'reanchor',  re = isfinite(dr);            % after every step (where there is a data frame)
            case 'adaptive',  re = dr >= reanchorTol;       % only where the heading has drifted
        end
        if useHold && ~strcmp(mode, 'transport')
            leaving = (zoneS | forcedS) & ~zoneE;   % back in trusted data: re-anchor only if the drift is small
            re(leaving) = isfinite(dr(leaving)) & dr(leaving) <= exitTol;
            ex(leaving) = dr(leaving) > exitTol;
            re(zoneE) = false;                  % never in the zone
            if vetoOn, re(vfire) = false; end
        end
        En = reanchorRows(En, Et, re);
    end
    [ka, ta, th] = curvatureTwist(Om1, m);
    for q = 1:numel(idx)
        k = idx(q);
        kappa{k}(end+1, 1) = ka(q);                                    %#ok<AGROW>
        tau{k}(end+1, 1)   = ta(q);                                    %#ok<AGROW>
        turn{k}(end+1, 1)  = h * th(q);                                %#ok<AGROW>
        if useHold
            held{k}(end+1, :, :) = Om1(q, :, :);                       %#ok<AGROW>
            spread{k}(end+1, 1) = sgS(q);                              %#ok<AGROW>
            if zoneD(q)                                                % consecutive steps started in the degenerate zone
                coast(k) = coast(k) + 1;
            elseif ~forcedS(q)                                         % the steps of 'ExitDelay' neither count nor reset
                coast(k) = 0;
            end
            if ~zoneS(q) && ~forcedS(q) && isfinite(maxRot(q)) && th(q) <= holdCap   % a trusted node enters the ring
                slot = mod(hcount(k), holdK) + 1;
                hbuf(k, slot, :, :) = Om1(q, :, :);
                sbuf(k, slot, :, :) = blM(q, :, :);  dbuf(k, slot) = dvM(q);  abuf(k, slot) = sA(q);
                hcount(k) = hcount(k) + 1;
            end
        end
        if ~isfinite(maxRot(q))                                        % the step is not taken:
            stop{k}   = 'nan';                                         % connection not finite
            active(k) = false;
        elseif vfire(q)                                                % the step is not taken: the history
            j = jv(q);                                                 % veto returns the tract to node j,
            tracts{k} = tracts{k}(1:j, :);  frames{k} = frames{k}(1:j, :, :);      % which coasts on the ring it had there
            drift{k} = drift{k}(1:j);  reanch{k} = reanch{k}(1:j);
            coasted{k} = coasted{k}(1:j);  confid{k} = confid{k}(1:j);    % node j stays trusted: the coast starts after it
            kappa{k} = kappa{k}(1:j-1);  tau{k} = tau{k}(1:j-1);  turn{k} = turn{k}(1:j-1);
            held{k} = held{k}(1:j-1, :, :);  spread{k} = spread{k}(1:j-1);
            r = k + (mod(j - 1, vetoK + 1)) * nseed;
            hbuf(k, :, :, :) = reshape(rhH(r, :), [1 holdK m m]);  sbuf(k, :, :, :) = reshape(rhS(r, :), [1 holdK m-1 m-1]);
            dbuf(k, :) = rhD(r, :);  abuf(k, :) = rhA(r, :);  hcount(k) = rhC(r);
            X(k, :) = tracts{k}(j, :);  E(k, :, :) = frames{k}(j, :, :);
            vneed(k) = nn(k) + 1;                                      % the end of the vetoed step
            nn(k) = j;  vco(k) = true;  coast(k) = 0;  grace(k) = 0;  vcount(k) = vcount(k) + 1;
        elseif maxRot(q) > maxTurn                                     % the step is not taken:
            stop{k}   = 'turn';                                        % too large a rotation
            active(k) = false;
        elseif ~zoneE(q) && dt(q) > maxTurn                            % the step is not taken: the data's
            stop{k}   = 'turn';                                        % e_1 at its end is more than MaxTurn
            active(k) = false;                                         % from the heading at its start
        elseif useHold && coast(k) > coastMax                          % the step is not taken: coasted too long
            stop{k}   = 'coast';
            active(k) = false;
        elseif ex(q)                                                   % the step is not taken: the coast has
            stop{k}   = 'exit';                                        % left the bundle
            active(k) = false;
        elseif cx(q)                                                   % the step is not taken: it would cross
            stop{k}   = 'cross';                                       % a tract of the same bundle
            active(k) = false;
        elseif fo(q)                                                   % the step is not taken: the march of
            stop{k}   = 'focus';                                       % MMF-4s reaches the focus of its fan
            active(k) = false;
        elseif all(Xn(q, 1:expdim) >= lims(:, 1).' & Xn(q, 1:expdim) <= lims(:, 2).')
            if useLaw                                                  % the law's state at the new node
                lawOn(k) = zoneE(q);
                if zoneE(q)
                    lawOm(k, :, :) = lawNext.Om(q, :, :);  lawSg(k) = lawNext.sg(q);  lawU(k) = lawNext.u(q);
                end
            end
            tracts{k}(end+1, :)    = Xn(q, :);                         %#ok<AGROW>
            frames{k}(end+1, :, :) = En(q, :, :);                      %#ok<AGROW>
            drift{k}(end+1, 1)  = dr(q);                               %#ok<AGROW>
            reanch{k}(end+1, 1) = re(q);                               %#ok<AGROW>
            if useHold
                coasted{k}(end+1, 1) = zoneE(q);                       %#ok<AGROW>
                confid{k}(end+1, 1)  = cE(q);                          %#ok<AGROW>
            end
            X(k, :)    = Xn(q, :);
            E(k, :, :) = En(q, :, :);
            nn(k) = nn(k) + 1;
            if vend(q), vco(k) = false; vfloor(k) = nn(k); end         % the vetoed coast has ended
        else                                                           % leaves the domain
            tracts{k}(end+1, :)    = clipToBoundary(X(k, :), Xn(q, :), lims, expdim);   %#ok<AGROW>
            frames{k}(end+1, :, :) = En(q, :, :);                      %#ok<AGROW>
            drift{k}(end+1, 1)  = dr(q);                               %#ok<AGROW>
            reanch{k}(end+1, 1) = re(q);                               %#ok<AGROW>
            if useHold
                coasted{k}(end+1, 1) = zoneE(q);                       %#ok<AGROW>
                confid{k}(end+1, 1)  = cE(q);                          %#ok<AGROW>
            end
            kappa{k}(end+1, 1) = NaN;                                  %#ok<AGROW>
            tau{k}(end+1, 1)   = NaN;                                  %#ok<AGROW>
            turn{k}(end+1, 1)  = NaN;                                  %#ok<AGROW>
            if useHold
                held{k}(end+1, :, :) = NaN;                            %#ok<AGROW>
                spread{k}(end+1, 1) = NaN;                             %#ok<AGROW>
            end
            stop{k}   = 'boundary';
            active(k) = false;
        end
    end
    if useHold && m == 3 && any(re)             % the ring follows the frame onto the re-anchored normal plane
        rows = idx(re);
        [hbuf(rows, :, :, :), sbuf(rows, :, :, :)] = turnRing(hbuf(rows, :, :, :), sbuf(rows, :, :, :), En(re, :, :), Eb(re, :, :));
    end
    if vetoOn                                   % the ring at the new node of every tract that took its step
        rows = idx(active(idx) & nn(idx) == cellfun(@(P) size(P, 1), tracts(idx)) & ~vfire);
        r = rows + mod(nn(rows) - 1, vetoK + 1) * nseed;
        rhH(r, :) = reshape(hbuf(rows, :, :, :), numel(rows), []);  rhS(r, :) = reshape(sbuf(rows, :, :, :), numel(rows), []);
        rhD(r, :) = dbuf(rows, :);  rhA(r, :) = abuf(rows, :);  rhC(r) = hcount(rows);
    end
end
% curvature, twist and turn at the last node of the tracts that reached smax
idx = find(active);
if ~isempty(idx)
    OmL = omegaFcn(X(idx, :), E(idx, :, :));
    sgL = NaN(numel(idx), 1);
    if useHold                                  % the held (or marched) connection where the last node is in the zone
        cL = evalCL(X(idx, 1:expdim), reshape(E(idx, 1, :), numel(idx), []));
        zL = isfinite(cL) & (cL < holdMin | grace(idx) > 0 | vco(idx));
        sgL = sigmaFcn(X(idx, :), E(idx, :, :));
        [Om0, bl0, dv0, sb0, sg0] = ringMeans(hbuf(idx, :, :, :), sbuf(idx, :, :, :), dbuf(idx, :), abuf(idx, :), hcount(idx), m, holdTwist, holdCurv);
        if marchHold
            [holdOm, sgH] = marchedOmega(Om0, bl0, dv0, h * (nn(idx) - 1) - sb0, uMin);
            if useLaw
                on = lawOn(idx);
                holdOm(on, :, :) = lawOm(idx(on), :, :);  sgH(on) = lawSg(idx(on));
            end
        else
            holdOm = Om0;  sgH = sg0;
        end
        OmL(zL, :, :) = holdOm(zL, :, :);
        sgL(zL) = sgH(zL);
    end
    [ka, ta, th] = curvatureTwist(OmL, m);
    for q = 1:numel(idx)
        k = idx(q);
        kappa{k}(end+1, 1) = ka(q);                                    %#ok<AGROW>
        tau{k}(end+1, 1)   = ta(q);                                    %#ok<AGROW>
        turn{k}(end+1, 1)  = h * th(q);                                %#ok<AGROW>
        if useHold
            held{k}(end+1, :, :) = OmL(q, :, :);                       %#ok<AGROW>
            spread{k}(end+1, 1) = sgL(q);                              %#ok<AGROW>
        end
    end
end
if ~useHold
    spread  = cellfun(@(P) NaN(size(P, 1), 1), tracts, 'UniformOutput', false);
    confid  = spread;
    coasted = cellfun(@(P) false(size(P, 1), 1), tracts, 'UniformOutput', false);
end

info = struct();
info.stop    = stop;
info.clipped = strcmp(stop, 'boundary');
info.s       = cell(nseed, 1);
for k = 1:nseed
    npts = size(tracts{k}, 1);
    s = h * (0:npts-1).';
    if info.clipped(k) && npts > 1              % the cut point: chord length of the last piece
        s(end) = s(end-1) + norm(tracts{k}(end, :) - tracts{k}(end-1, :));
    end
    info.s{k} = s;
end
info.kappa   = kappa;
info.tau     = tau;
info.turn    = turn;
info.drift      = drift;
info.reanchored = reanch;
info.coasted    = coasted;
info.confidence = confid;
info.held       = held;
info.spread     = spread;
info.hold       = struct('on', useHold, 'holdMin', holdMin, 'holdK', holdK, 'holdCap', holdCap, ...
                         'coastMax', coastMax, 'exitTol', exitTol, 'exitDelay', exitDelay, 'holdTwist', holdTwist, 'holdCurvature', holdCurv, 'holdDivergence', holdDiv, ...
                         'bundleRule', ~isempty(bundleRule), 'bundleLaw', useLaw, 'vetoK', vetoK, 'vetoTol', vetoTol);
info.vetoes  = vcount;
info.order   = order;
info.mode    = mode;
info.twist   = twist;
info.maxTurn = maxTurn;
info.reanchorTol = reanchorTol;
end

% =========================================================================
function [Xn, En, Om1, maxRot, zone1] = lieStep(X, E, h, order, omegaFcn, evalCL, holdMin, holdOm, expdim, forced)
% One step of the Lie group integrator for the pose (X, E) with the
% exponential poseExp; Om1 is Omega at the start of the step and maxRot
% the largest rotation h*theta over the stages of the step (NaN when the
% connection is not finite at some stage; such an Omega is replaced by
% zero so that the remaining stages stay finite, and the caller rejects
% the step). With the hold scheme (evalCL and holdOm given) a stage whose
% confidence is below holdMin uses the held connection holdOm of its tract
% instead of the data - an array [n x m x m], or (MMF-4s) a function of the
% stage's abscissa c in [0, 1] along the step returning that array, the
% connection prescribed along the coast; the rows forced (logical, 'ExitDelay')
% coast at every stage whatever the finite confidence; zone1 flags the tracts
% whose start node is such.
if nargin < 6
    evalCL = [];  holdMin = NaN;  holdOm = [];  expdim = size(X, 2);
end
if nargin < 10, forced = false(size(X, 1), 1); end
if isa(holdOm, 'function_handle'), heldAt = holdOm; else, heldAt = @(c) holdOm; end
stage = @(Xs, Es, c) stageOmega(Xs, Es, omegaFcn, evalCL, holdMin, heldAt(c), expdim, forced);
[Om1, bad, th, zone1] = stage(X, E, 0);
rates = th;
switch order
    case 1                                      % Lie-Euler: Omega frozen
        [Xn, En] = poseExp(X, E, Om1, 1, h);
    case 2                                      % exponential midpoint rule
        [X2, E2] = poseExp(X, E, Om1, 1, h / 2);
        [Om2, b2, r2] = stage(X2, E2, 1/2);
        bad = bad | b2;  rates = [rates, r2];
        [Xn, En] = poseExp(X, E, Om2, 1, h);
    case 4                                      % CF4 (Celledoni, Marthinsen & Owren 2003)
        [X2, E2] = poseExp(X, E, Om1, 1, h / 2);
        [Om2, b2, r2] = stage(X2, E2, 1/2);
        [X3, E3] = poseExp(X, E, Om2, 1, h / 2);
        [Om3, b3, r3] = stage(X3, E3, 1/2);
        [X4, E4] = poseExp(X2, E2, Om3 - Om1 / 2, 1/2, h);
        [Om4, b4, r4] = stage(X4, E4, 1);
        bad = bad | b2 | b3 | b4;  rates = [rates, r2, r3, r4];
        [Xn, En] = poseExp(X,  E,  (3*Om1 + 2*Om2 + 2*Om3 -   Om4) / 12, 1/2, h);
        [Xn, En] = poseExp(Xn, En, ( -Om1 + 2*Om2 + 2*Om3 + 3*Om4) / 12, 1/2, h);
end
maxRot = h * max(rates, [], 2);
maxRot(bad) = NaN;
end

function [Om, bad, th, zone] = stageOmega(X, E, omegaFcn, evalCL, holdMin, holdOm, expdim, forced)
% The connection at a stage: from the data, or the held one of the tract
% where the confidence of the data is below holdMin (hold scheme) or the
% tract still coasts beyond the zone's edge (forced, 'ExitDelay'), the data
% not being read there. A non-finite confidence is no data at all: the
% connection is read, and stops the tract if it is not finite.
zone = false(size(X, 1), 1);
if ~isempty(evalCL)
    c = evalCL(X(:, 1:expdim), reshape(E(:, 1, :), size(X, 1), []));
    zone = isfinite(c(:)) & (c(:) < holdMin | forced(:));
end
if any(zone)
    Om = holdOm;
    if ~all(zone)
        Om(~zone, :, :) = omegaFcn(X(~zone, :), E(~zone, :, :));
    end
else
    Om = omegaFcn(X, E);
end
[Om, bad, th] = finiteOmega(Om);
end

function Om = heldOmega(hbuf, hcount, m, holdTwist, holdCurv)
% The held connection of every tract: the mean of its ring of trusted
% connections (zero before the first one), with the twist zeroed or kept
% and the curvature kept or zeroed.
n  = size(hbuf, 1);
Om = zeros(n, m, m);
has = hcount > 0;
if any(has)
    Om(has, :, :) = reshape(mean(hbuf(has, :, :, :), 2, 'omitnan'), nnz(has), m, m);
end
if strcmp(holdCurv, 'zero')
    Om(:, 1, :) = 0;  Om(:, :, 1) = 0;
end
if m == 3 && strcmp(holdTwist, 'zero')
    Om(:, 2, 3) = 0;  Om(:, 3, 2) = 0;
end
Om(~isfinite(Om)) = 0;
end

function [Om, bl, dv, sb, sg] = ringMeans(hbuf, sbuf, dbuf, abuf, hcount, m, holdTwist, holdCurv)
% The held connection of every tract (heldOmega) with the means over its
% ring of the transverse block of the connection (m = 2: the transverse
% component; sg is the trace), of the divergence term of MMF-4s and of the
% arc length of the ring's nodes (all zero before the first trusted node,
% so that such a tract coasts straight with sigma = 0).
Om = heldOmega(hbuf, hcount, m, holdTwist, holdCurv);
n  = size(hbuf, 1);
bl = zeros(n, m - 1, m - 1);  dv = zeros(n, 1);  sb = zeros(n, 1);
has = hcount > 0;
bl(has, :, :) = reshape(mean(sbuf(has, :, :, :), 2, 'omitnan'), nnz(has), m - 1, m - 1);
dv(has) = mean(dbuf(has, :), 2, 'omitnan');
sb(has) = mean(abuf(has, :), 2, 'omitnan');
bl(~isfinite(bl)) = 0;  dv(~isfinite(dv)) = 0;  sb(~isfinite(sb)) = 0;
sg = bl(:, 1);
if m == 3
    sg = bl(:, 1, 1) + bl(:, 2, 2);
end
end

function [hb, sb] = turnRing(hb, sb, En, Eo)
% The rings of tracts whose frames Eo were re-anchored to En [n x 3 x 3], in
% the coordinates of En: the curvature vectors hb(:,:,1,2:3) and the blocks
% sb [n x K x 2 x 2] turned by the rotation nearest to the projection between
% the two normal planes, as in transverseOf (empty slots stay NaN; the twist
% of a held connection is unchanged by a turn of the plane).
n = size(En, 1);
G = zeros(n, 2, 2);
for p = 1:2
    for q = 1:2
        G(:, p, q) = sum(En(:, p + 1, :) .* Eo(:, q + 1, :), 3);
    end
end
c = G(:, 1, 1) + G(:, 2, 2);
d = G(:, 1, 2) - G(:, 2, 1);
r = hypot(c, d);
c = c ./ r;  d = d ./ r;
a2 = hb(:, :, 1, 2);  a3 = hb(:, :, 1, 3);
k2 = c .* a2 + d .* a3;
k3 = -d .* a2 + c .* a3;
hb(:, :, 1, 2) = k2;  hb(:, :, 2, 1) = -k2;
hb(:, :, 1, 3) = k3;  hb(:, :, 3, 1) = -k3;
b11 = sb(:, :, 1, 1);  b12 = sb(:, :, 1, 2);  b21 = sb(:, :, 2, 1);  b22 = sb(:, :, 2, 2);
sb(:, :, 1, 1) = c .^ 2 .* b11 + c .* d .* (b12 + b21) + d .^ 2 .* b22;
sb(:, :, 1, 2) = c .^ 2 .* b12 + c .* d .* (b22 - b11) - d .^ 2 .* b21;
sb(:, :, 2, 1) = c .^ 2 .* b21 + c .* d .* (b22 - b11) - d .^ 2 .* b12;
sb(:, :, 2, 2) = d .^ 2 .* b11 - c .* d .* (b12 + b21) + c .^ 2 .* b22;
end

function [Om, sg, u] = marchedOmega(Om0, sg0, dv0, tau, uMin)
% MMF-4s: the connection and the transverse component marched over the arc
% length tau from the ring's centroid by the flatness law of the planar
% connection with the divergence term held, d sigma/ds = -sigma^2 + D, and
% the closure div(kappa e_1) = 0, d kappa/ds = -kappa sigma, i.e. with
% sigma = u'/u, u'' = D u, u(0) = 1,
% u'(0) = sigma_0: Omega = Omega_0/u (every component, so the ratio
% kappa/sigma stays constant for D = 0, where u = 1 + sigma_0 tau) and sigma
% = u'/u. The connection is floored at u = uMin (the caller stops the tract
% before the focus of a converging fan); a row without a trusted node
% (Omega_0 = 0, sigma_0 = 0) coasts straight. With three frame vectors
% sg0 is the transverse block [n x 2 x 2], see marchedOmega3.
if size(Om0, 2) == 3
    [Om, sg, u] = marchedOmega3(Om0, sg0, tau, uMin);
    return
end
[u, up] = marchLaw(sg0, dv0, tau);
sg = up ./ u;
sg(~isfinite(sg)) = 0;
Om = Om0 ./ max(u, uMin);
end

function [Om, sg, u] = marchedOmega3(Om0, B0, tau, uMin)
% MMF-4s with three frame vectors, the divergence term zero: dB/ds = -B^2
% and da/ds = -B' a for the transverse block B and the curvature vector
% a = (Omega_12, Omega_13), so U = I + tau B_0, B = B_0 U^-1 and a = U^-T a_0
% (the adjugate of U' over u = det U = 1 + tau tr B_0 + tau^2 det B_0, floored
% at uMin as in marchedOmega), and sigma = tr B = u'/u. The rest of the
% held connection is left as it is.
b11 = B0(:, 1, 1);  b12 = B0(:, 1, 2);  b21 = B0(:, 2, 1);  b22 = B0(:, 2, 2);
tr  = b11 + b22;
dt  = b11 .* b22 - b12 .* b21;
u   = 1 + tau .* tr + tau .^ 2 .* dt;
sg  = (tr + 2 * tau .* dt) ./ u;
sg(~isfinite(sg)) = 0;
uf  = max(u, uMin);
a2  = Om0(:, 1, 2);  a3 = Om0(:, 1, 3);
k2  = ((1 + tau .* b22) .* a2 - tau .* b21 .* a3) ./ uf;
k3  = ((1 + tau .* b11) .* a3 - tau .* b12 .* a2) ./ uf;
Om  = Om0;
Om(:, 1, 2) = k2;  Om(:, 2, 1) = -k2;
Om(:, 1, 3) = k3;  Om(:, 3, 1) = -k3;
end

function [u, up] = marchLaw(sg0, dv0, tau)
% u and u' of u'' = D u, u(0) = 1, u'(0) = sigma_0, at tau (columns):
% 1 + sigma_0 tau for D = 0, cosh/sinh for D > 0, cos/sin for D < 0.
u  = 1 + sg0 .* tau;
up = sg0 + 0 * tau;
pos = dv0 > 0;  neg = dv0 < 0;
if any(pos)
    a = sqrt(dv0(pos));  t = tau(pos);  s0 = sg0(pos);
    u(pos)  = cosh(a .* t) + s0 ./ a .* sinh(a .* t);
    up(pos) = a .* sinh(a .* t) + s0 .* cosh(a .* t);
end
if any(neg)
    a = sqrt(-dv0(neg));  t = tau(neg);  s0 = sg0(neg);
    u(neg)  = cos(a .* t) + s0 ./ a .* sin(a .* t);
    up(neg) = -a .* sin(a .* t) + s0 .* cos(a .* t);
end
end

function [sig, dv, bl] = transverseOf(X, E, evalW, evalMF, expdim, m, delta, wantDiv, evalCL, holdMin)
% The transverse component sigma = sum_{j>1} w_1j(e_j) = div e_1 of the
% data's connection at the rows of X along the direction of travel: the
% data's frame with its e_1 pointing the way the carried frame's does (the
% carried frame itself when there is no data frame), and, with wantDiv,
% the divergence term D = div(kappa e_2) = d kappa/d e_2 - kappa^2 of
% MMF-4s, kappa = w_12(e_1) the data's curvature, its derivative along e_2
% by central differences over delta (planar connection; NaN where a probe
% leaves the data or falls into the zone, c < holdMin, so that a node next
% to the zone contributes no divergence to the ring). bl is what the ring
% keeps of it: sig for m = 2, for m = 3 the transverse block B_pq =
% e_p . d e_1/d e_q, p, q = 2, 3, of the data's frame along the travel [n x
% 2 x 2], turned onto the carried normal plane by the rotation nearest to
% the projection between the two planes (the carried frame being the only
% one in which it means the same from node to node).
n  = size(X, 1);
wq = evalW(X(:, 1:expdim), reshape(E(:, 1, :), n, []));      % [n x m x m x expdim]
if isempty(evalMF)
    Et = E;  s = ones(n, 1);
else
    Et = evalMF(X(:, 1:expdim), reshape(E(:, 1, :), n, []));  % [n x m x spacedim]
    s  = sign(sum(Et(:, 1, 1:expdim) .* E(:, 1, 1:expdim), 3));  s(s == 0) = 1;
end
sig = zeros(n, 1);
for j = 2:m
    for k = 1:expdim
        sig = sig + wq(:, 1, j, k) .* Et(:, j, k);
    end
end
sig = s .* sig;
bl  = sig;
if m == 3
    Et(:, 3, :) = s .* Et(:, 3, :);                           % with e_1 along the travel, right-handed
    bl = zeros(n, 2, 2);
    for p = 1:2
        for q = 1:2
            for k = 1:expdim
                bl(:, p, q) = bl(:, p, q) + wq(:, 1, p + 1, k) .* Et(:, q + 1, k);
            end
        end
    end
    bl(:, 1, :) = s .* bl(:, 1, :);                           % the data's block on its own normal plane
    G  = zeros(n, 2, 2);                                      % that plane turned onto the carried one
    for p = 1:2
        for q = 1:2
            G(:, p, q) = sum(E(:, p + 1, 1:expdim) .* Et(:, q + 1, 1:expdim), 3);
        end
    end
    c  = G(:, 1, 1) + G(:, 2, 2);
    d  = G(:, 1, 2) - G(:, 2, 1);
    r  = hypot(c, d);
    R  = zeros(n, 2, 2);
    R(:, 1, 1) = c ./ r;   R(:, 1, 2) = d ./ r;
    R(:, 2, 1) = -d ./ r;  R(:, 2, 2) = c ./ r;
    bl = pagemul(pagemul(R, bl), permute(R, [1 3 2]));
end
dv  = zeros(n, 1);
if wantDiv
    e1 = s .* reshape(Et(:, 1, 1:expdim), n, expdim);         % the data's frame along the travel
    e2 = s .* reshape(Et(:, 2, 1:expdim), n, expdim);
    Y  = X(:, 1:expdim);
    kp = dataCurvature(Y + delta * e2, e1, evalW, evalMF, expdim);
    km = dataCurvature(Y - delta * e2, e1, evalW, evalMF, expdim);
    k0 = dataCurvature(Y, e1, evalW, evalMF, expdim);
    dv = (kp - km) / (2 * delta) - k0 .^ 2;
    cp = evalCL(Y + delta * e2, e1);  cm = evalCL(Y - delta * e2, e1);
    dv((isfinite(cp) & cp < holdMin) | (isfinite(cm) & cm < holdMin)) = NaN;
end
end

function kap = dataCurvature(Y, e1ref, evalW, evalMF, expdim)
% The data's curvature kappa = w_12(e_1) at the rows of Y, for the data's
% e_1 pointing the way e1ref does.
n  = size(Y, 1);
wq = evalW(Y, e1ref);
Et = evalMF(Y, e1ref);
e1 = reshape(Et(:, 1, 1:expdim), n, expdim);
s  = sign(sum(e1 .* e1ref, 2));  s(s == 0) = 1;
kap = zeros(n, 1);
for k = 1:expdim
    kap = kap + wq(:, 1, 2, k) .* e1(:, k);
end
kap = s .* kap;
end

function [Om, bad, th] = finiteOmega(Om)
% Rotation rate th = sqrt(-tr(Om^2)/2) of the skew matrices Om [n x m x m];
% rows with a non-finite Om are flagged in bad and replaced by zero.
[n, m, ~] = size(Om);
bad = any(any(~isfinite(Om), 2), 3);
Om(bad, :, :) = 0;
th2 = zeros(n, 1);
for i = 1:m-1
    for j = i+1:m
        th2 = th2 + Om(:, i, j) .^ 2;
    end
end
th = sqrt(th2);
end

% =========================================================================
function [Xn, En] = poseExp(X, E, Om, cvel, h)
% Exponential of the pose with the constant angular velocity Om (skew,
% [n x m x m]) and the velocity cvel*e_1 over the arc length h:
%     En = exp(h Om) E,    Xn = X + cvel * sum_j [int_0^h exp(s Om) ds]_1j e_j,
% with exp(h Om) = I + A Om + B Om^2 and int_0^h exp(s Om) ds = h I + B Om
% + D Om^2 (Rodrigues; theta^2 = -tr(Om^2)/2 is the squared rotation rate).
[n, m, ~] = size(E);
th2 = zeros(n, 1);
for i = 1:m-1
    for j = i+1:m
        th2 = th2 + Om(:, i, j) .^ 2;
    end
end
th = sqrt(th2);
a  = h * th;
A  = zeros(n, 1);  B = zeros(n, 1);  D = zeros(n, 1);
big = a > 1e-3;
A(big) = sin(a(big)) ./ th(big);
B(big) = (1 - cos(a(big))) ./ th2(big);
D(big) = (a(big) - sin(a(big))) ./ (th2(big) .* th(big));
sm = ~big;
a2 = a(sm) .^ 2;
A(sm) = h   * (1   - a2/6   + a2.^2/120);
B(sm) = h^2 * (1/2 - a2/24  + a2.^2/720);
D(sm) = h^3 * (1/6 - a2/120 + a2.^2/5040);

Om2 = zeros(n, m, m);
for i = 1:m
    for j = 1:m
        for l = 1:m
            Om2(:, i, j) = Om2(:, i, j) + Om(:, i, l) .* Om(:, l, j);
        end
    end
end
R = A .* Om + B .* Om2;
for i = 1:m
    R(:, i, i) = R(:, i, i) + 1;
end
V1 = B .* Om(:, 1, :) + D .* Om2(:, 1, :);      % first row of the integral, [n x 1 x m]
V1(:, 1, 1) = V1(:, 1, 1) + h;

En = zeros(size(E));
Xn = X;
for j = 1:m
    Ej = E(:, j, :);                            % [n x 1 x spacedim]
    for i = 1:m
        En(:, i, :) = En(:, i, :) + R(:, i, j) .* Ej;
    end
    Xn = Xn + cvel * V1(:, 1, j) .* reshape(Ej, n, []);
end
end

% =========================================================================
function Om = connectionMatrix(X, E, evalW, evalMF, useGauge, expdim, m, dropTwist)
% Omega_ij = w_ij(e_1) at the positions X for the frames E, expressed in
% the coordinates of E: Omega = G Omegat G' with G = E Et' the rotation
% from the data frame Et(x) to E, and Omegat = sum_k wij_k(x) e_1k.
n  = size(X, 1);
wq = evalW(X(:, 1:expdim), reshape(E(:, 1, :), n, []));   % [n x m x m x expdim]
if size(wq, 2) ~= m || size(wq, 3) ~= m
    error('advanceFrame:badConnection', 'The connection has %d frames but the frames array %d.', size(wq, 2), m);
end
t  = reshape(E(:, 1, 1:expdim), n, expdim);     % direction of travel
Om = zeros(n, m, m);
for k = 1:expdim
    Om = Om + wq(:, :, :, k) .* t(:, k);
end
if useGauge
    Et = evalMF(X(:, 1:expdim), reshape(E(:, 1, :), n, []));   % [n x m x spacedim]
    G  = zeros(n, m, m);
    for i = 1:m
        for j = 1:m
            G(:, i, j) = sum(E(:, i, :) .* Et(:, j, :), 3);
        end
    end
    Om = pagemul(pagemul(G, Om), permute(G, [1 3 2]));
end
Om = (Om - permute(Om, [1 3 2])) / 2;           % exactly skew (removes the differentiation error)
if dropTwist
    Om(:, 2, 3) = 0;
    Om(:, 3, 2) = 0;
end
end

function C = pagemul(A, B)
% Page-wise product C(:,i,j) = sum_l A(:,i,l) B(:,l,j) for [n x m x m] arrays.
[n, m, ~] = size(A);
C = zeros(n, m, m);
for i = 1:m
    for j = 1:m
        for l = 1:m
            C(:, i, j) = C(:, i, j) + A(:, i, l) .* B(:, l, j);
        end
    end
end
end

% =========================================================================
function E = orthonormalize(E)
% Gram-Schmidt of the frames E [n x m x spacedim] keeping the direction of
% e_1; for three frames in 3-D, e_3 = e_1 x e_2 (right-handed).
[n, m, d] = size(E);
nonfinite_in = any(any(~isfinite(E), 3), 2);
e1 = reshape(E(:, 1, :), n, d);
e1 = e1 ./ sqrt(sum(e1 .^ 2, 2));
e2 = reshape(E(:, 2, :), n, d);
e2 = e2 - sum(e2 .* e1, 2) .* e1;
e2 = e2 ./ sqrt(sum(e2 .^ 2, 2));
E(:, 1, :) = reshape(e1, n, 1, d);
E(:, 2, :) = reshape(e2, n, 1, d);
if m == 3
    if d == 3
        e3 = cross(e1, e2, 2);
    else
        e3 = reshape(E(:, 3, :), n, d);
        e3 = e3 - sum(e3 .* e1, 2) .* e1 - sum(e3 .* e2, 2) .* e2;
        e3 = e3 ./ sqrt(sum(e3 .^ 2, 2));
    end
    E(:, 3, :) = reshape(e3, n, 1, d);
end
% A frame that was not finite on input (a data frame outside the data, e.g.
% beyond a mask) is returned as NaN: the connection there is not finite
% either, and the march stops before the step. A finite frame that cannot be
% orthonormalized (zero or parallel vectors) is an error.
bad = any(any(~isfinite(E), 3), 2);
if any(bad & ~nonfinite_in)
    error('advanceFrame:degenerateFrame', 'A frame could not be orthonormalized (zero or parallel vectors).');
end
E(bad, :, :) = NaN;
end

% =========================================================================
function ang = headingAngle(A, B)
% Angle between the first vectors of the frames A and B, row by row, modulo
% sign (the data's e_1 is a line): accurate for small angles, at most pi/2,
% NaN where a frame is not finite.
g   = sum(A(:, 1, :) .* B(:, 1, :), 3);
r   = A(:, 1, :) - g .* B(:, 1, :);                 % part of e_1 normal to the other
ang = atan2(sqrt(sum(r .^ 2, 3)), abs(g));
end

function D = bundleDivergence(Xa, Ea, OmS, sgS, zone, label, omegaFcn, expdim)
% The law ('BundleLaw'): the divergence term D = d kappa/d e_2 - kappa^2 of
% the march at the coasting rows (zone), with d kappa/d e_2 from the
% neighbouring tracts of the bundle. A neighbour j of the row i, at the
% offset r, gives kappa_j - kappa_i = (d kappa/d e_1)(r . e_1) + (d kappa/d
% e_2)(r . e_2) to first order, with d kappa/d e_1 = -kappa sigma the march's
% own law along the fibre (the neighbours of a row are at its arc length,
% not on its normal); the nearest neighbour on either side of the row along
% its e_2 (the coast's curvature where that tract coasts, the data's where it
% is trusted) gives an estimate, combined as the central difference of
% unequal spacings (exact for a quadratic), one-sided with one neighbour,
% and D = 0 without any. Zero at the trusted rows.
n   = size(Xa, 1);
D   = zeros(n, 1);
kap = reshape(OmS(:, 1, 2), n, 1);
if any(~zone)
    OmT = omegaFcn(Xa(~zone, :), Ea(~zone, :, :));
    kap(~zone) = OmT(:, 1, 2);
end
e1 = reshape(Ea(:, 1, 1:expdim), n, expdim);
e2 = reshape(Ea(:, 2, 1:expdim), n, expdim);
for g = unique(label(:)).'
    rows = find(label == g);
    if numel(rows) < 2, continue; end
    for i = rows(zone(rows)).'
        r  = Xa(rows, 1:expdim) - Xa(i, 1:expdim);                % offsets of the bundle's tracts from the row
        d  = r * e2(i, :).';  a = r * e1(i, :).';
        k0 = kap(i);
        est = (kap(rows) - k0 + k0 * sgS(i) * a) ./ d;            % d kappa/d e_2 implied by each neighbour
        ok = rows ~= i & isfinite(est) & d ~= 0;
        dp = d;  dp(~ok | d <= 0) = Inf;  [dP, jp] = min(dp);
        dm = -d; dm(~ok | d >= 0) = Inf;  [dM, jm] = min(dm);
        if isfinite(dP) && isfinite(dM)
            dk = (est(jp) * dM + est(jm) * dP) / (dP + dM);
        elseif isfinite(dP)
            dk = est(jp);
        elseif isfinite(dM)
            dk = est(jm);
        else
            continue
        end
        D(i) = dk - k0 ^ 2;
    end
end
end

function cx = bundleCross(Xa, Xn, Ea, zone, label, expdim)
% The rows whose step from Xa to Xn crosses the step of another tract of
% the same bundle (equal label), at least one of the two coasting (zone):
% the offset between the two, in the plane normal to their mean heading at
% the start of the step, reverses or vanishes over the step. In the plane
% this is the swap of their order across the bundle; in space two fibres
% that pass each other keep their offset, so they are not stopped.
n  = size(Xa, 1);
cx = false(n, 1);
e1 = reshape(Ea(:, 1, 1:expdim), n, expdim);
for g = unique(label(:)).'
    rows = find(label == g);
    for a = 1:numel(rows)
        for b = a + 1:numel(rows)
            i = rows(a);  j = rows(b);
            if ~(zone(i) || zone(j)), continue; end
            m = e1(i, :) + e1(j, :);
            if norm(m) < 1e-12, continue; end
            m  = m / norm(m);
            wo = Xa(i, 1:expdim) - Xa(j, 1:expdim);  wo = wo - (wo * m.') * m;
            wn = Xn(i, 1:expdim) - Xn(j, 1:expdim);  wn = wn - (wn * m.') * m;
            if wo * wn.' <= 0
                cx(i) = true;  cx(j) = true;
            end
        end
    end
end
end

function E = reanchorRows(E, Et, sel)
% Re-anchoring: the rows sel of E are replaced by the data frames Et with
% the sign of each vector chosen by the transported frame (an even number
% of flips, so right-handed).
if ~any(sel)
    return
end
[n, m, ~] = size(E);
g = zeros(n, m);
for i = 1:m
    g(:, i) = sum(E(:, i, :) .* Et(:, i, :), 3);
end
s = sign(g);
s(s == 0) = 1;
if m == 3
    odd = prod(s, 2) < 0;                        % make det = +1: flip the least aligned normal vector
    if any(odd)
        [~, imin] = min(abs(g(odd, 2:3)), [], 2);
        rows = find(odd);
        for q = 1:numel(rows)
            s(rows(q), imin(q) + 1) = -s(rows(q), imin(q) + 1);
        end
    end
end
for i = 1:m
    E(sel, i, :) = s(sel, i) .* Et(sel, i, :);
end
end

% =========================================================================
function [kappa, tau, theta] = curvatureTwist(Om, m)
% Curvature |kappa| = sqrt(w_12(e_1)^2 + w_13(e_1)^2), twist tau = w_23(e_1)
% and rotation rate theta = sqrt(kappa^2 + tau^2) of the frame.
if m == 3
    kappa = sqrt(Om(:, 1, 2) .^ 2 + Om(:, 1, 3) .^ 2);
    tau   = Om(:, 2, 3);
else
    kappa = abs(Om(:, 1, 2));
    tau   = zeros(size(kappa));
end
theta = sqrt(kappa .^ 2 + tau .^ 2);
end

% =========================================================================
function xb = clipToBoundary(xa, xn, lims, expdim)
% Point where the chord xa -> xn (xa inside, xn outside) crosses the boundary.
t = 1;
for c = 1:expdim
    for e = 1:2
        dx = xn(c) - xa(c);
        if dx ~= 0
            tc = (lims(c, e) - xa(c)) / dx;
            if tc >= 0 && tc < t
                t = tc;
            end
        end
    end
end
xb = xa + t * (xn - xa);
end
