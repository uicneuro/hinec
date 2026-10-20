function [tracts, frames, info] = advanceFrame(w, mf, dom, seeds, dir0, h, smax, varargin)
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
%              more than 'MaxTurn' (in a re-anchoring mode: also when the
%              re-anchoring at its end would), or when the connection is
%              not finite
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
%                 other modes. The default is an order of magnitude below
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
%                 In the re-anchoring modes the re-anchoring at the end of
%                 a step is a rotation of the frame too: a step after which
%                 the data's e_1 differs from the transported one by more
%                 than MaxTurn (the drift, at most pi/2) is not taken
%                 either, with the same reason 'turn' - the counterpart of
%                 the angle test between consecutive directions of
%                 streamline tracking, which a snap to the data would
%                 otherwise bypass. A tract also stops before a step whose
%                 connection is not finite at some stage.
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
%                        (the next step, or the re-anchoring at its end,
%                        would turn by more than MaxTurn),
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
%                        before any re-anchoring; 0 at the seed; in
%                        'transport' mode NaN after the seed (not evaluated
%                        there)
%               reanchored
%                        cell, logical per point: the frame stored at the
%                        point is the data frame (re-anchored) rather than
%                        the transported one; false at the seed. Like the
%                        frame, drift and reanchored at a clipped boundary
%                        point refer to the uncut end of the last step.
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
%   that had it is archived as MMF4vca_old.m. Since 2026-10-20 a re-anchoring
%   larger than 'MaxTurn' stops the tract ('turn') instead of being applied.

%% --- options --------------------------------------------------------------
order   = 4;
mode    = 'transport';
twist   = 'data';
frame0  = [];
maxTurn = 0.5;
reanchorTol = 1e-3;
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
if isa(w, 'function_handle')
    evalW = w;
elseif isnumeric(w) && any(ndims(w) == [7 9])
    gdim = (ndims(w) - 3) / 2;                   % 2: [N+1 N+1 Nx Ny m m 2], 3: [N+1 N+1 N+1 Nx Ny Nz m m 3]
    if size(w, ndims(w)) ~= expdim || gdim ~= expdim
        error('advanceFrame:badConnection', ...
            'w is the nodal connection of a %d-D domain with %d differentiation directions, but the seeds have %d coordinates.', ...
            gdim, size(w, ndims(w)), expdim);
    end
    if gdim == 2
        evalW = @(X) interpGLL(w, dom, X(:,1), X(:,2));
    else
        evalW = @(X) interpGLL(w, dom, X(:,1), X(:,2), X(:,3));
    end
else
    error('advanceFrame:badConnection', ...
        ['w must be the nodal connection [N+1 N+1 Nx Ny m m 2] or [N+1 N+1 N+1 Nx Ny Nz m m 3] ' ...
         'or a function handle.']);
end
if isa(mf, 'function_handle')
    evalMF = @(X) orthonormalize(mf(X));
elseif isnumeric(mf) && any(ndims(mf) == [6 8])
    gdim = (ndims(mf) - 2) / 2;                  % 2: [N+1 N+1 Nx Ny m spacedim], 3: [N+1 N+1 N+1 Nx Ny Nz m spacedim]
    if gdim ~= expdim
        error('advanceFrame:badFrames', ...
            'mf is the nodal frame field of a %d-D domain, but the seeds have %d coordinates.', gdim, expdim);
    end
    if gdim == 2
        evalMF = @(X) orthonormalize(interpGLL(mf, dom, X(:,1), X(:,2)));
    else
        evalMF = @(X) orthonormalize(interpGLL(mf, dom, X(:,1), X(:,2), X(:,3)));
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
    E = evalMF(seeds);                           % [nseed x m x spacedim]
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
for k = 1:nseed
    tracts{k} = X(k, :);
    frames{k} = E(k, :, :);
    drift{k}  = 0;
    reanch{k} = false;
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
    [Xn, En, Om1, maxRot] = lieStep(Xa, Ea, h, order, omegaFcn);
    dr = NaN(numel(idx), 1);                    % drift of the transported heading at Xn
    re = false(numel(idx), 1);                  % rows whose frame is re-anchored at Xn
    switch mode
        case 'reanchor'                         % after every step
            [En, dr, re] = reanchorFrames(Xn, En, evalMF, expdim, -Inf);
        case 'adaptive'                         % only where the heading has drifted
            [En, dr, re] = reanchorFrames(Xn, En, evalMF, expdim, reanchorTol);
    end
    [ka, ta, th] = curvatureTwist(Om1, m);
    for q = 1:numel(idx)
        k = idx(q);
        kappa{k}(end+1, 1) = ka(q);                                    %#ok<AGROW>
        tau{k}(end+1, 1)   = ta(q);                                    %#ok<AGROW>
        turn{k}(end+1, 1)  = h * th(q);                                %#ok<AGROW>
        if ~isfinite(maxRot(q))                                        % the step is not taken:
            stop{k}   = 'nan';                                         % connection not finite
            active(k) = false;
        elseif maxRot(q) > maxTurn                                     % the step is not taken:
            stop{k}   = 'turn';                                        % too large a rotation
            active(k) = false;
        elseif re(q) && dr(q) > maxTurn                                % the step is not taken: the
            stop{k}   = 'turn';                                        % re-anchoring at its end would
            active(k) = false;                                         % turn e_1 by more than MaxTurn
        elseif all(Xn(q, 1:expdim) >= lims(:, 1).' & Xn(q, 1:expdim) <= lims(:, 2).')
            tracts{k}(end+1, :)    = Xn(q, :);                         %#ok<AGROW>
            frames{k}(end+1, :, :) = En(q, :, :);                      %#ok<AGROW>
            drift{k}(end+1, 1)  = dr(q);                               %#ok<AGROW>
            reanch{k}(end+1, 1) = re(q);                               %#ok<AGROW>
            X(k, :)    = Xn(q, :);
            E(k, :, :) = En(q, :, :);
        else                                                           % leaves the domain
            tracts{k}(end+1, :)    = clipToBoundary(X(k, :), Xn(q, :), lims, expdim);   %#ok<AGROW>
            frames{k}(end+1, :, :) = En(q, :, :);                      %#ok<AGROW>
            drift{k}(end+1, 1)  = dr(q);                               %#ok<AGROW>
            reanch{k}(end+1, 1) = re(q);                               %#ok<AGROW>
            kappa{k}(end+1, 1) = NaN;                                  %#ok<AGROW>
            tau{k}(end+1, 1)   = NaN;                                  %#ok<AGROW>
            turn{k}(end+1, 1)  = NaN;                                  %#ok<AGROW>
            stop{k}   = 'boundary';
            active(k) = false;
        end
    end
end
% curvature, twist and turn at the last node of the tracts that reached smax
idx = find(active);
if ~isempty(idx)
    [ka, ta, th] = curvatureTwist(omegaFcn(X(idx, :), E(idx, :, :)), m);
    for q = 1:numel(idx)
        k = idx(q);
        kappa{k}(end+1, 1) = ka(q);                                    %#ok<AGROW>
        tau{k}(end+1, 1)   = ta(q);                                    %#ok<AGROW>
        turn{k}(end+1, 1)  = h * th(q);                                %#ok<AGROW>
    end
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
info.order   = order;
info.mode    = mode;
info.twist   = twist;
info.maxTurn = maxTurn;
info.reanchorTol = reanchorTol;
end

% =========================================================================
function [Xn, En, Om1, maxRot] = lieStep(X, E, h, order, omegaFcn)
% One step of the Lie group integrator for the pose (X, E) with the
% exponential poseExp; Om1 is Omega at the start of the step and maxRot
% the largest rotation h*theta over the stages of the step (NaN when the
% connection is not finite at some stage; such an Omega is replaced by
% zero so that the remaining stages stay finite, and the caller rejects
% the step).
[Om1, bad, th] = finiteOmega(omegaFcn(X, E));
rates = th;
switch order
    case 1                                      % Lie-Euler: Omega frozen
        [Xn, En] = poseExp(X, E, Om1, 1, h);
    case 2                                      % exponential midpoint rule
        [X2, E2] = poseExp(X, E, Om1, 1, h / 2);
        [Om2, b2, r2] = finiteOmega(omegaFcn(X2, E2));
        bad = bad | b2;  rates = [rates, r2];
        [Xn, En] = poseExp(X, E, Om2, 1, h);
    case 4                                      % CF4 (Celledoni, Marthinsen & Owren 2003)
        [X2, E2] = poseExp(X, E, Om1, 1, h / 2);
        [Om2, b2, r2] = finiteOmega(omegaFcn(X2, E2));
        [X3, E3] = poseExp(X, E, Om2, 1, h / 2);
        [Om3, b3, r3] = finiteOmega(omegaFcn(X3, E3));
        [X4, E4] = poseExp(X2, E2, Om3 - Om1 / 2, 1/2, h);
        [Om4, b4, r4] = finiteOmega(omegaFcn(X4, E4));
        bad = bad | b2 | b3 | b4;  rates = [rates, r2, r3, r4];
        [Xn, En] = poseExp(X,  E,  (3*Om1 + 2*Om2 + 2*Om3 -   Om4) / 12, 1/2, h);
        [Xn, En] = poseExp(Xn, En, ( -Om1 + 2*Om2 + 2*Om3 + 3*Om4) / 12, 1/2, h);
end
maxRot = h * max(rates, [], 2);
maxRot(bad) = NaN;
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
wq = evalW(X(:, 1:expdim));                     % [n x m x m x expdim]
if size(wq, 2) ~= m || size(wq, 3) ~= m
    error('advanceFrame:badConnection', 'The connection has %d frames but the frames array %d.', size(wq, 2), m);
end
t  = reshape(E(:, 1, 1:expdim), n, expdim);     % direction of travel
Om = zeros(n, m, m);
for k = 1:expdim
    Om = Om + wq(:, :, :, k) .* t(:, k);
end
if useGauge
    Et = evalMF(X(:, 1:expdim));                % [n x m x spacedim]
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
function [E, drift, sel] = reanchorFrames(X, E, evalMF, expdim, tol)
% Re-anchoring: the rows of E whose heading drift is tol or more (sel) are
% replaced by the data frame at X with the sign of each vector chosen by
% the transported frame (an even number of flips, so right-handed); tol =
% -Inf replaces every row, tol = Inf none (the drift is at most pi/2).
% drift is the angle between the transported e_1 and the data's e_1
% (modulo sign) at X, for every row, before the replacement.
[n, m, ~] = size(E);
Et = evalMF(X(:, 1:expdim));
g  = zeros(n, m);
for i = 1:m
    g(:, i) = sum(E(:, i, :) .* Et(:, i, :), 3);
end
r     = E(:, 1, :) - g(:, 1) .* Et(:, 1, :);       % part of e_1 normal to the data's e_1
drift = atan2(sqrt(sum(r .^ 2, 3)), abs(g(:, 1))); % accurate for small angles
sel   = drift >= tol;
if ~any(sel)
    return
end
s = sign(g);
s(s == 0) = 1;
if m == 3
    odd = prod(s, 2) < 0;                        % make det = +1: flip the least aligned vector
    if any(odd)
        [~, imin] = min(abs(g(odd, :)), [], 2);
        rows = find(odd);
        for q = 1:numel(rows)
            s(rows(q), imin(q)) = -s(rows(q), imin(q));
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
