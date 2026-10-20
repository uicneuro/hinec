function [tracks, meta] = nim_tractography_mmf(nim, options)
% nim_tractography_mmf: Moving frames tractography on the voxel data (advanceFrame).
%
%   [tracks, meta] = nim_tractography_mmf(nim, options)
%
% The fibre is followed as a curve that carries its own orthonormal frame
% E = (e1; e2; e3), advanced by the connection 1-forms of the data's frame field
% on the direction of travel,
%
%     dx/ds = e1,     dE/ds = Omega E,     Omega_ij = w_ij(e1) = sum_k wij_k e1_k,
%
% with the closed-form exponential step of advanceFrame.m (a rotation of the
% frame together with a circular arc of the position: exact for a planar arc of
% constant curvature at any step size; on a fibre with torsion exact to the
% order of the integrator, since the twist is dropped, see below) and the Lie
% group integrators of order 1, 2 and 4 (Lie-Euler, exponential midpoint, CF4).
% The mathematics and the verification are in the moving-frames test suite
% (Moving Frame Tractography/matlab, doc/advanceFrame_notes.tex); this file only
% connects the scheme to the voxel data of a nim.
%
% THE DATA'S FRAME FIELD AND CONNECTION, FROM THE VOXELS. The principal
% eigenvector is a LINE field (the sign the eigensolver writes per voxel is
% arbitrary), so the direction between voxels is read from the interpolated
% dyadic v1*v1' (sign-invariant, the same six interpolants the hinec tracker
% uses, with the same kernel 'interpolation.method' and the same 'upsample'):
%
%     e1(x)  = principal axis of D(x) = interp(v1 v1'),
%     e2, e3 = a completion of e1 (reference-axis projection; the data's normal
%              plane is arbitrary in DTI, so nothing is read from it),
%     w_1j,k = e_j . d e1/dx_k,  j = 2, 3   (the curvature components)
%     w_23,k = 0                            (the twist is not a datum here)
%
% d e1/dx_k is a central difference of e1 of the interpolated field (sign
% aligned, step 'mmf.gradient_delta' voxels), so the connection is that of
% the field the streamline trackers integrate: with the spline kernel a C2
% field, with cubic or trilinear a field whose derivative jumps at the voxel
% faces (where the exponential step is then only as accurate as any
% integrator on a kinked right-hand side; prefer spline). The scheme runs
% with 'Twist' 'zero' (rotation-minimizing normal plane), so the arbitrary
% completion never enters the position: the step uses only the curvature
% vector (e1 . grad) e1 of the field, expressed in the carried frame through
% the gauge G = E Et' that advanceFrame computes at every evaluation. That
% gauge also absorbs the sign flips of the axial data.
%
% TERMINATION. The connection is returned non-finite wherever a point is
% outside the volume, outside the propagation mask (or the brain mask), below
% the FA floor or without a well-conditioned principal axis; advanceFrame
% stops a tract BEFORE a step at whose stages the connection is not finite.
% The stages are the start of the step, its midpoint (orders 2 and 4) and its
% end (order 4), so with orders 1 and 2 the scheme never examines the end of a
% step: a last node that lies outside tissue is dropped here, and every node
% of a track is a valid integration node inside tissue. Leaving the volume
% stops a tract as well; the cut boundary point is dropped, so tracks consist
% of integration nodes only, like those of the other trackers. The angle
% budget termination.angle_max (degrees per voxel of arc; through
% nim_angle_limit it is min(angle_max*h, 90) degrees per step, none for
% angle_max <= 0) is advanceFrame's 'MaxTurn', applied to the rotation
% h*|kappa| of a step and to the re-anchoring at its end (the snap of the
% carried heading to the data's e1, which could otherwise bypass the budget).
% Unlike the streamline trackers, whose sign-aligned turn cannot exceed 90
% degrees, the 90 degree ceiling is a real cap here. The curvature of the
% interpolated principal-axis field is large wherever the dyadic is nearly
% degenerate (crossings, planar regions), so the budget stops more half-tracks
% there than the node-to-node angle test of the streamline trackers does.
%
% OPTIONS READ (all defaulted and validated by nim_config_schema; see
% docs/YAML_CONFIG.md)
%   field              'dti' only. CSD peaks need the incoming tangent to
%                      choose a peak, and the connection handles of
%                      advanceFrame are functions of position alone.
%   integrator         'euler' | 'rk2' | 'rk4' -> Lie group order 1 | 2 | 4.
%                      'rkf45' has no counterpart (the exponential step has no
%                      embedded error estimate) and is rejected.
%   step_size          h, the arc length of a step in voxels (fixed)
%   interp_method, upsample
%                      kernel and spatial sampling of the dyadic field, as hinec
%   mmf_mode           'transport' | 'reanchor' | 'adaptive' (advanceFrame 'Mode')
%   mmf_reanchor_tol   drift tolerance of the adaptive mode, radians
%   mmf_gradient_delta central-difference step for d e1/dx_k, voxels
%   angle_thresh       degrees of turning per VOXEL OF ARC (nim_angle_limit),
%                      applied as advanceFrame's 'MaxTurn', see above
%   max_arc            arc length followed from the seed in each direction
%                      (max_steps * step_size when absent, legacy calls)
%   termination_fa, min_length, seed_mask, seed_density, seed_strategy,
%   propagation_mask, trace, trace_max                      as for hinec
%   ACT masks (wm/gm/csf) are not used by this tracker.
%
% OUTPUT (the tracker contract, docs/TRACKER_INTERFACE.md)
%   tracks  cell {T x 1}; tracks{t} is N x 3, [backward half, reversed; seed;
%           forward half], integration nodes only, arc length >= min_length
%   meta    seed_index (1 x T), seed_points (T x 3), n_seeds;
%           termination: reason/count over the half-tracks of the live seeds
%           (max_arc, angle, fa, outside_mask, outside, outside_domain,
%           no_direction; empty when no seed is alive);
%           mmf: the scheme settings, reanchored_fraction of the steps,
%           seeds_alive, and seeds_skipped (reason/count of the seeds outside
%           the volume or the mask, below the FA floor or without a direction);
%           trace (debug.trace): per traced seed, forward/backward records with
%           the fields the hinec tracker writes - pos, dir (the carried e1),
%           fa (interpolated), h, chord, turn (degrees: the rotation of the
%           step starting at the node), flip (always false: the gauge absorbs
%           the sign flips), n_steps, termination - plus drift (angle to the
%           data's e1), reanchored, kappa and s (arc length); nim_track_diagnose
%           reads them. A skipped seed has a zero-step record with its reason.
%
% Requires advanceFrame.m, nim_seed_offsets.m, nim_angle_limit.m.

dims = size(nim.FA);
h    = options.step_size;
o    = @(f, d) getdef(options, f, d);

% ---- options ---------------------------------------------------------------
fld_name = lower(char(string(o('field', 'dti'))));
if ~strcmp(fld_name, 'dti')
    error('mmf:field', ['nim_tractography_mmf: field ''%s'' is not supported. The connection handles of ' ...
        'advanceFrame are functions of position alone, and selecting a CSD peak needs the incoming ' ...
        'tangent; run field: dti, or use the hinec tracker for CSD.'], fld_name);
end
integ = lower(char(string(o('integrator', 'rk4'))));
switch integ
    case 'euler', order = 1;
    case 'rk2',   order = 2;
    case 'rk4',   order = 4;
    case 'rkf45'
        error('mmf:integrator', ['nim_tractography_mmf: integrator rkf45 has no moving-frames ' ...
            'counterpart (the exponential step carries no embedded error estimate). Use euler, rk2 or rk4.']);
    otherwise
        error('mmf:integrator', 'nim_tractography_mmf: unknown integrator ''%s''.', integ);
end
mode  = lower(char(string(o('mmf_mode', 'adaptive'))));
rtol  = o('mmf_reanchor_tol', 1e-3);
delta = o('mmf_gradient_delta', 1e-3);
maxTurn = deg2rad(nim_angle_limit(options.angle_thresh, h));   % min(angle_max*h, 90) deg per step; Inf disables
smax  = o('max_arc', o('max_steps', 1000) * h);
min_length     = o('min_length', 0);
termination_fa = options.termination_fa;
act_given = ~isempty(o('wm_mask', [])) || ~isempty(o('gm_mask', [])) || ~isempty(o('csf_mask', []));
if act_given
    fprintf('MMF: ACT masks are present but not used by this tracker (termination by FA, mask and volume).\n');
end

% ---- the field: interpolated dyadic, FA, masks (as the hinec tracker) ---------
switch lower(char(string(o('interp_method', 'trilinear'))))
    case 'cubic',  kernel = 'cubic';
    case 'spline', kernel = 'spline';
    otherwise,     kernel = 'linear';
end
upsample = o('upsample', 1);
if isempty(upsample) || upsample <= 0, upsample = 1; end
v1  = reshape(nim.evec(:, :, :, :, 1), [dims 3]);
Dc  = {v1(:,:,:,1).*v1(:,:,:,1), v1(:,:,:,2).*v1(:,:,:,2), v1(:,:,:,3).*v1(:,:,:,3), ...
       v1(:,:,:,1).*v1(:,:,:,2), v1(:,:,:,1).*v1(:,:,:,3), v1(:,:,:,2).*v1(:,:,:,3)};
FAg = double(nim.FA);
gv  = {1:dims(1), 1:dims(2), 1:dims(3)};
if abs(upsample - 1) > 1e-9
    ng   = max(2, round((dims - 1) * upsample) + 1);
    gnew = {linspace(1, dims(1), ng(1)), linspace(1, dims(2), ng(2)), linspace(1, dims(3), ng(3))};
    [GX, GY, GZ] = ndgrid(gnew{1}, gnew{2}, gnew{3});
    for q = 1:6
        gi = griddedInterpolant(gv, Dc{q}, kernel, 'none'); Dc{q} = gi(GX, GY, GZ);
    end
    gi = griddedInterpolant(gv, FAg, kernel, 'none'); FAg = gi(GX, GY, GZ);
    gv = gnew;
    fprintf('MMF: upsample %.4g -> field sampled at %.4g voxel spacing (%d x %d x %d)\n', ...
        upsample, (dims(1)-1)/(ng(1)-1), ng(1), ng(2), ng(3));
end
DI = cell(1, 6);
for q = 1:6, DI{q} = griddedInterpolant(gv, Dc{q}, kernel, 'none'); end
FAI = griddedInterpolant(gv, FAg, 'linear', 'none');      % FA linear: a stop test, not a field (see hinec)
prop = o('propagation_mask', []);
if ~isempty(prop)
    prop = logical(prop > 0);
elseif isfield(nim, 'mask') && ~isempty(nim.mask)
    prop = logical(nim.mask > 0.5);
else
    prop = true(dims);
end
fprintf('MMF: moving frames tractography (advanceFrame), order %d, mode %s%s, %s dyadic field, h = %.3g, MaxTurn %.3g rad/step\n', ...
    order, mode, ifelse(strcmp(mode, 'adaptive'), sprintf(' (tol %.2g rad)', rtol), ''), kernel, h, maxTurn);

% The field, as one value the local functions below receive explicitly (local
% functions, not nested ones: advanceFrame calls the handles from inside its own
% loops, and nested functions would share this workspace).
fld = struct('DI', {DI}, 'FAI', FAI, 'prop', prop, 'dims', dims, 'fa_min', termination_fa, 'delta', delta);
names = reason_names();

% ---- seeds -------------------------------------------------------------------
seed_points = seeds_from_mask(options.seed_mask, options, dims);
n_seeds = size(seed_points, 1);
[e1s, alive_seed, code0] = e1_at(seed_points, fld);
skipped_cnt = arrayfun(@(c) nnz(code0 == c), 1:numel(names))';   % the seeds not tracked, by reason
has = find(skipped_cnt > 0);
fprintf('MMF: %d seeds, %d inside tissue with a direction\n', n_seeds, nnz(alive_seed));

% ---- which seeds are traced per step (debug.trace) ---------------------------
do_trace = false(n_seeds, 1);
if o('trace', false)
    tmax = o('trace_max', 200);
    if tmax <= 0 || tmax >= n_seeds, do_trace(:) = true; else, do_trace(round(linspace(1, n_seeds, tmax))) = true; end
    fprintf('MMF: per-step tracing enabled for %d of %d seeds\n', nnz(do_trace), n_seeds);
end
dirnames = {'forward', 'backward'};

% ---- march all live seeds together, in chunks, both directions -----------------
dom   = struct('xlim', [1 dims(1)], 'ylim', [1 dims(2)], 'zlim', [1 dims(3)]);
wFcn  = @(X) connection_at(X, fld);
mfFcn = @(X) frames_at(X, fld);
live  = find(alive_seed);
half  = cell(n_seeds, 2);                       % forward, backward halves (seed first)
term  = cell(n_seeds, 2);                       % termination reason per half
trace_fwd = cell(n_seeds, 1); trace_bwd = cell(n_seeds, 1);
n_re = 0; n_st = 0;                             % re-anchored steps / steps
chunk = 4000;
progress = nim_tracking_progress('MMF', numel(live));
t0 = tic;
for c0 = 1:chunk:numel(live)
    idx = live(c0:min(c0 + chunk - 1, numel(live)));
    S   = seed_points(idx, :);
    d0  = e1s(idx, :);
    for dirn = 1:2
        sgn = 3 - 2 * dirn;                     % +1 forward, -1 backward
        [tr, fr, info] = advanceFrame(wFcn, mfFcn, dom, S, sgn * d0, h, smax, ...
            'Order', order, 'Mode', mode, 'ReanchorTol', rtol, 'Twist', 'zero', 'MaxTurn', maxTurn);
        % Integration nodes only: drop a clipped boundary point, and a step end
        % the scheme never examined (orders 1 and 2) that lies outside tissue.
        Pq = tr(:);
        for q = 1:numel(idx)
            if info.clipped(q), Pq{q} = Pq{q}(1:end-1, :); end
        end
        ends = cell2mat(cellfun(@(P) P(end, :), Pq, 'UniformOutput', false));
        [~, okend] = e1_at(ends, fld);
        for q = 1:numel(idx)
            k = idx(q);
            P = Pq{q}; xdrop = zeros(0, 3);
            if ~okend(q) && size(P, 1) > 1, xdrop = P(end, :); P = P(1:end-1, :); end
            half{k, dirn} = P;
            term{k, dirn} = why_stopped(info.stop{q}, P, fr{q}, h, fld, xdrop);
            re = info.reanchored{q}(2:size(P, 1));
            n_re = n_re + nnz(re); n_st = n_st + numel(re);
            if do_trace(k)
                rec = trace_record(P, fr{q}, info, q, h, term{k, dirn}, dirnames{dirn}, fld);
                if dirn == 1, trace_fwd{k} = rec; else, trace_bwd{k} = rec; end
            end
        end
    end
    progress(numel(idx), false);
end
progress(numel(live), true);
fprintf('MMF: tracking %.1f s\n', toc(t0));

% ---- combine the halves, apply the length rule ---------------------------------
tracks = cell(n_seeds, 1);
kept   = false(n_seeds, 1);
for k = live(:)'
    fwd = half{k, 1}; bwd = half{k, 2};
    if size(fwd, 1) > 1, fwd = fwd(2:end, :); else, fwd = zeros(0, 3); end
    if size(bwd, 1) > 1, bwd = flipud(bwd(2:end, :)); else, bwd = zeros(0, 3); end
    track = [bwd; seed_points(k, :); fwd];
    if size(track, 1) >= 2 && sum(sqrt(sum(diff(track, 1, 1).^2, 2))) >= min_length
        tracks{k} = track; kept(k) = true;
    end
end
tracks = reshape(tracks(kept), [], 1);

meta = struct();
meta.seed_index  = reshape(find(kept), 1, []);
meta.seed_points = seed_points(kept, :);
meta.n_seeds     = n_seeds;
meta.mmf = struct('order', order, 'mode', mode, 'reanchor_tol', rtol, 'twist', 'zero', ...
                  'max_turn_rad', maxTurn, 'gradient_delta', delta, 'kernel', kernel, ...
                  'reanchored_fraction', n_re / max(n_st, 1), 'seeds_alive', nnz(alive_seed), ...
                  'seeds_skipped', struct('reason', {reshape(names(has), [], 1)}, 'count', {skipped_cnt(has)}));
tt = term(alive_seed, :); tt = tt(:); tt = tt(~cellfun(@isempty, tt));
if ~isempty(tt)
    [u, ~, ii] = unique(tt); cnt = accumarray(ii, 1); [cnt, ord] = sort(cnt, 'descend');
    meta.termination = struct('reason', {u(ord)}, 'count', {cnt});
    fprintf('MMF termination (%d half-tracks):', numel(tt));
    for j = 1:numel(cnt), fprintf(' %s=%.1f%%', u{ord(j)}, 100 * cnt(j) / numel(tt)); end
    fprintf('\n');
else
    meta.termination = struct('reason', {cell(0, 1)}, 'count', {zeros(0, 1)});
end
if any(do_trace)
    ti = find(do_trace(:))';
    for k = ti(~alive_seed(ti))                 % a skipped seed: a zero-step record with its reason
        trace_fwd{k} = trace_skipped(seed_points(k, :), names{code0(k)}, dirnames{1}, fld);
        trace_bwd{k} = trace_skipped(seed_points(k, :), names{code0(k)}, dirnames{2}, fld);
    end
    meta.trace = struct('seed_index', num2cell(ti), ...
                        'seed_point', num2cell(seed_points(ti, :), 2)', ...
                        'forward',    trace_fwd(ti)', ...
                        'backward',   trace_bwd(ti)');
end
fprintf('MMF: %d of %d seeds gave a track >= %.3g voxels; %.1f%% of the steps re-anchored\n', ...
    numel(tracks), n_seeds, min_length, 100 * meta.mmf.reanchored_fraction);

end

% =============================================================================
% The data's frame field and connection, vectorized over query points. Local
% functions with the field passed in explicitly (see the note at fld above).
% =============================================================================
function names = reason_names()
% The failure codes of e1_at, in the order they are tested.
names = {'outside', 'outside_mask', 'outside_domain', 'fa', 'no_direction'};
end

function [e1, ok, code] = e1_at(X, fld)
% Principal axis of the interpolated dyadic at the rows of X. code is 0 where
% the point is usable and otherwise the first failing test: 1 outside the
% volume, 2 outside the propagation mask, 3 FA not finite (outside the
% interpolation domain), 4 below the FA floor, 5 no well-conditioned
% principal axis (reason_names); ok = code == 0.
n  = size(X, 1);
e1 = NaN(n, 3); code = ones(n, 1);
inside = all(X >= 1, 2) & all(X <= fld.dims, 2);
if any(inside)
    xi  = X(inside, :);
    vox = min(max(round(xi), 1), fld.dims);
    inm = fld.prop(sub2ind(fld.dims, vox(:, 1), vox(:, 2), vox(:, 3)));
    fa  = fld.FAI(xi(:, 1), xi(:, 2), xi(:, 3));
    c   = zeros(nnz(inside), 6);
    for q = 1:6, c(:, q) = fld.DI{q}(xi(:, 1), xi(:, 2), xi(:, 3)); end
    v = principal_axis(c(:, 1), c(:, 2), c(:, 3), c(:, 4), c(:, 5), c(:, 6));
    ci = zeros(nnz(inside), 1);
    ci(~inm) = 2;
    ci(ci == 0 & ~isfinite(fa)) = 3;
    ci(ci == 0 & fa < fld.fa_min) = 4;
    ci(ci == 0 & any(~isfinite(v), 2)) = 5;
    v(ci ~= 0, :) = NaN;
    e1(inside, :) = v;
    code(inside)  = ci;
end
ok = code == 0;
end

function E = frames_at(X, fld)
% Data frames [n x 3 x 3]: e1 the principal axis, e2 by reference-axis
% projection (the Cartesian axis least aligned with e1), e3 = e1 x e2.
n  = size(X, 1);
e1 = e1_at(X, fld);
[~, ax] = min(abs(e1), [], 2);
r  = zeros(n, 3); r(sub2ind([n 3], (1:n)', ax)) = 1;
e2 = r - sum(r .* e1, 2) .* e1;
e2 = e2 ./ sqrt(sum(e2 .^ 2, 2));
e3 = cross(e1, e2, 2);
E = zeros(n, 3, 3);
E(:, 1, :) = reshape(e1, n, 1, 3);
E(:, 2, :) = reshape(e2, n, 1, 3);
E(:, 3, :) = reshape(e3, n, 1, 3);
end

function w = connection_at(X, fld)
% Connection components w(n,i,j,k) = e_j . d e_i/dx_k of the data frames:
% the curvature rows (i = 1) from central differences of the interpolated
% principal axis, sign aligned; the twist row w_23 is set to zero (the
% scheme runs with 'Twist' 'zero'). Non-finite where the point is invalid.
n  = size(X, 1);
E  = frames_at(X, fld);
e1 = reshape(E(:, 1, :), n, 3); e2 = reshape(E(:, 2, :), n, 3); e3 = reshape(E(:, 3, :), n, 3);
w  = zeros(n, 3, 3, 3);
for k = 1:3
    dX = zeros(1, 3); dX(k) = fld.delta;
    ep = e1_at(X + dX, fld); em = e1_at(X - dX, fld);
    ep = ep .* sign(sum(ep .* e1, 2) + 1e-300);       % align the signs of the line field
    em = em .* sign(sum(em .* e1, 2) + 1e-300);
    de1 = (ep - em) / (2 * fld.delta);
    w12 = sum(e2 .* de1, 2);
    w13 = sum(e3 .* de1, 2);
    w(:, 1, 2, k) =  w12;  w(:, 2, 1, k) = -w12;
    w(:, 1, 3, k) =  w13;  w(:, 3, 1, k) = -w13;
end
bad = any(~isfinite(reshape(w, n, [])), 2) | any(~isfinite(e1), 2);
w(bad, :, :, :) = NaN;
end

function reason = why_stopped(stop, P, F, hh, fld, xdrop)
% advanceFrame's reasons, with 'nan' resolved into what made the connection
% non-finite: the node dropped for lying outside tissue (xdrop, if any), else
% the first failing point among the stages of the next step (its midpoint and
% end along the carried e1, then the node itself), each with the six probes
% of the central-difference stencil, in the order the tests are applied.
switch stop
    case 'boundary', reason = 'outside';
    case 'length',   reason = 'max_arc';
    case 'turn',     reason = 'angle';
    otherwise
        if isempty(P), reason = 'outside'; return; end
        p   = P(end, :); e = reshape(F(size(P, 1), 1, :), 1, 3);
        pts = [xdrop; p + 0.5 * hh * e; p + hh * e; p];
        d   = fld.delta;
        st  = [0 0 0; d 0 0; -d 0 0; 0 d 0; 0 -d 0; 0 0 d; 0 0 -d];
        Q   = reshape(permute(reshape(pts, [], 1, 3) + reshape(st, 1, [], 3), [2 1 3]), [], 3);
        [~, ~, code] = e1_at(Q, fld);
        c = code(find(code > 0, 1));
        if isempty(c), reason = 'no_direction'; else, names = reason_names(); reason = names{c}; end
end
end

% =============================================================================
function rec = trace_record(P, F, info, q, h, reason, dirname, fld)
% Per-step record of one half-track for debug.trace: the fields of the hinec
% tracker's record (nim_track_diagnose) plus the moving-frames quantities.
n = size(P, 1);
rec = struct('direction', dirname, 'pos', P, 'dir', reshape(F(1:n, 1, :), n, 3), ...
             'fa', fld.FAI(P(:, 1), P(:, 2), P(:, 3)), ...
             'h', [h * ones(n - 1, 1); NaN], 'chord', [sqrt(sum(diff(P, 1, 1) .^ 2, 2)); NaN], ...
             'turn', [rad2deg(info.turn{q}(1:n-1)); NaN], 'flip', false(n, 1), ...
             'drift', info.drift{q}(1:n), 'reanchored', info.reanchored{q}(1:n), ...
             'kappa', info.kappa{q}(1:n), 's', info.s{q}(1:n), 'n_steps', n - 1, 'termination', reason);
end

function rec = trace_skipped(x, reason, dirname, fld)
% The zero-step record of a seed that was not tracked.
rec = struct('direction', dirname, 'pos', x, 'dir', NaN(1, 3), 'fa', fld.FAI(x(1), x(2), x(3)), ...
             'h', NaN, 'chord', NaN, 'turn', NaN, 'flip', false, 'drift', 0, 'reanchored', false, ...
             'kappa', NaN, 's', 0, 'n_steps', 0, 'termination', reason);
end

function seed_points = seeds_from_mask(seed_mask, options, dims)
% Seeds on the sub-voxel lattice of nim_seed_offsets (uniform) or jittered
% (random: +-0.2 voxel at density 1, +-0.4 above), as the hinec tracker.
idx = find(seed_mask);
[x, y, z] = ind2sub(dims, idx);
centres = [x y z];
density = max(1, round(double(getdef(options, 'seed_density', 1))));
if strcmpi(string(getdef(options, 'seed_strategy', 'uniform')), 'random')
    jit = ifelse(density <= 1, 0.4, 0.8);
    seed_points = repelem(centres, density, 1) + (rand(size(centres, 1) * density, 3) - 0.5) * jit;
    return
end
off = nim_seed_offsets(density);
seed_points = zeros(size(centres, 1) * size(off, 1), 3);
for k = 1:size(off, 1)
    seed_points((k-1) * size(centres, 1) + (1:size(centres, 1)), :) = centres + off(k, :);
end
end

function v = principal_axis(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz)
% Unit principal eigenvector of symmetric 3x3 matrices, vectorized over rows
% (the closed form of nim_principal_dir: Cardano on the deviatoric part, then
% the null space of D - lambda I from the largest cross product of its rows).
% NaN rows where the matrix is non-finite or degenerate.
n = numel(Dxx);
v = NaN(n, 3);
fin = isfinite(Dxx) & isfinite(Dyy) & isfinite(Dzz) & isfinite(Dxy) & isfinite(Dxz) & isfinite(Dyz);
p1 = Dxy.^2 + Dxz.^2 + Dyz.^2;
qm = (Dxx + Dyy + Dzz) / 3;
diag_case = fin & p1 < 1e-24;
if any(diag_case)
    [~, kmax] = max([Dxx(diag_case) Dyy(diag_case) Dzz(diag_case)], [], 2);
    vd = zeros(nnz(diag_case), 3); vd(sub2ind(size(vd), (1:nnz(diag_case))', kmax)) = 1;
    v(diag_case, :) = vd;
end
gen = fin & ~diag_case;
if ~any(gen), return; end
Dxx = Dxx(gen); Dyy = Dyy(gen); Dzz = Dzz(gen); Dxy = Dxy(gen); Dxz = Dxz(gen); Dyz = Dyz(gen);
p1 = p1(gen); qm = qm(gen);
p2 = (Dxx-qm).^2 + (Dyy-qm).^2 + (Dzz-qm).^2 + 2*p1;
p  = sqrt(p2 / 6);
Bxx = (Dxx-qm)./p; Byy = (Dyy-qm)./p; Bzz = (Dzz-qm)./p; Bxy = Dxy./p; Bxz = Dxz./p; Byz = Dyz./p;
r = (Bxx.*(Byy.*Bzz - Byz.*Byz) - Bxy.*(Bxy.*Bzz - Byz.*Bxz) + Bxz.*(Bxy.*Byz - Byy.*Bxz)) / 2;
r = min(1, max(-1, r));
lam = qm + 2*p.*cos(acos(r)/3);
r1 = [Dxx-lam, Dxy, Dxz]; r2 = [Dxy, Dyy-lam, Dyz]; r3 = [Dxz, Dyz, Dzz-lam];
c1 = cross(r1, r2, 2); c2 = cross(r2, r3, 2); c3 = cross(r3, r1, 2);
n1 = sum(c1.^2, 2); n2 = sum(c2.^2, 2); n3 = sum(c3.^2, 2);
[nb, kb] = max([n1 n2 n3], [], 2);
vg = c1; vg(kb == 2, :) = c2(kb == 2, :); vg(kb == 3, :) = c3(kb == 3, :);
vg = vg ./ sqrt(nb);
vg(p < 1e-12 | nb < 1e-24, :) = NaN;
v(gen, :) = vg;
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
