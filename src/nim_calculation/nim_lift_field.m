function L = nim_lift_field(nim, opts)
%NIM_LIFT_FIELD  The lifts of the fibre bundles on M = Omega x RP^2, from the raw diffusion signal.
%
%   L = nim_lift_field(nim, opts)
%
%   SPACES. Omega is the voxel domain, RP^2 = S^2/{u ~ -u} the lines through the origin (a fibre
%   direction has no sign), M = Omega x RP^2 the 5-dimensional space of positions and directions,
%   pi(x, u) = x. A bundle k with direction field u_k on its region Omega_k is the LIFT
%       Sigma_k = { (x, u_k(x)) : x in Omega_k },
%   a 3-dimensional submanifold of M meeting every fibre pi^-1(x) at most once. Bundles that cross
%   at x in Omega have different u there, so their lifts are disjoint in M; the connection 1-forms
%   omega_ij = <d e_i, e_j> of each lift's frame are defined where the single-valued principal
%   eigenvector field is a mixture of the bundles (nim_tractography_mmf, field: lift).
%
%   THE DATA. The signal at x for the unit gradient g and b-value b is the projection of the lifts
%   through x along pi, blurred by the single-fibre kernel:
%       S(x, g) / S0(x) = sum_j a_j(x) K(g, n_j),   K(g, n) = exp(-b (l_perp + (l_par - l_perp)(g.n)^2)),
%   + sum_d c_d exp(-b d)  (isotropic compartments, d in opts.iso: 'auto' = [l_perp, 3e-3], the hindered
%   and the free water, which carry no direction: without them the stick sum has to build isotropic signal
%   out of many sticks, which become false lifts), a_j >= 0 the amount of fibre along n_j, {n_j} a near-uniform set of ndirs directions on the half
%   sphere (Fibonacci). l_par, l_perp: opts.lambda, or the median l1 and (l2 + l3)/2 of the DTI
%   eigenvalues (nim.eval) of the voxels whose FA is above the ref_fa_pct-th percentile of the
%   tissue (single-fibre voxels). a(x) is the non-negative least-squares solution (lsqnonneg) of
%       min ||K a + K_iso c - y||^2 + mu (sum_j a_j)^2,   a, c >= 0   (mu = opts.sparsity),
%   the second term the prior of few lifts: without it a sum of many sticks reproduces broad or isotropic
%   signal as exactly as the isotropic compartments do, and its sticks become false lifts. The prior only
%   SELECTS (opts.debias): the atoms within 25 deg of a significant one (a_j >= 0.02 max a) are refitted
%   without it, because the penalty alone pulls the lobes of a crossing toward each other. OFF by default:
%   on ISMRM 2015 (slab z 38-52, scripts/lift_field_check.py) the refit returns the noise atoms beside a
%   lobe as false lifts (one lift in 25% of single-bundle voxels, 11.6 deg, against 85% and 7.3 deg).
%
%   THE DENSITY ON M.  p(x, u) = sum_j a_j(x) W(u.n_j),  W(t) = exp(kappa (t^2 - 1))  (Watson).
%
%   REGULARISATION ALONG THE FIBRES (opts.smooth, half-width in voxels, 0 = off). Each orientation
%   slice is smoothed along its own direction only, never across it:
%       a~_j(x) = sum_{t = -T..T} w_t a_j(x + t n_j),  w_t proportional to T + 1 - |t|,
%   so that a lift is reinforced by its own continuation and an isolated noise lobe is not.
%   OFF by default: the kernel runs along the straight line x + t n_j, while a curved fibre leaves it, its
%   direction at x + t n_j turned by about k |t| (curvature k). Where a straight run enters a bend, the
%   slice of the upstream direction keeps its support from behind and the slice of the local tangent loses
%   it, so the lobe is pulled upstream and a second lobe can form: on the ISMRM cingulum (k ~ 0.5 per voxel)
%   the lifted curvature fell from 0.15 to 0.02 at T = 2 (0.17-0.21 at T = 0, v1 0.22-0.30). ISMRM mmf
%   lift+hold F1 0.509 / 0.521 / 0.537 at T = 2 / 1 / 0.
%   THE CURVED KERNEL (opts.smooth_kernel 'curved'). A fibre through (x, n_j) is a curve in M whose direction
%   turns as it moves, so the support is gathered along circular arcs instead of the line: for curvature k in
%   opts.smooth_kappa and normal m = cos(q) e2 + sin(q) e3, q = 2 pi i / opts.smooth_normals,
%       gamma(t) = x + (sin(k t)/k) n_j + ((1 - cos(k t))/k) m,   d(t) = cos(k t) n_j + sin(k t) m,
%   S_(k,m)(x) = sum_t w_t a_(j(t))(gamma(t)) / sum_t w_t chi(gamma(t)),  j(t) the dictionary direction
%   nearest d(t), chi the tissue mask (trilinear), and a~_j(x) = max over (k, m), the line (k = 0) included.
%   A slice is supported by the continuation of its own fibre whichever way that fibre bends.
%   The directions are in the frame of nim.bvec, which is the frame of nim.evec (nim_dt_spd fits
%   the tensor with nim.bvec), the voxel-axis frame every tracker reads nim.evec in.
%
%   THE CONTINUATION PRIOR (opts.continuation = beta in [0, 1), 0 = off; needs opts.sparsity > 0). The prior
%   of few lifts weighs every dictionary direction alike. With beta > 0 a second deconvolution weighs them by
%   the geometry of the neighbouring lifts: the lifts of a first pass (u_k at the voxel nodes) and their nodal
%   connection (nim_lift_connection, opts.continuation_stencil) give, from each of the 6 neighbours y of x and
%   each lift u_k(y), the moving-frame continuation to x,
%       u^_k(x) = (u_k(y) + (du_k/dx)(y) (x - y)) / | . |      (curvature included, to first order),
%   and the penalty becomes mu (sum_j c_j(x) a_j)^2 with
%       c_j(x) = 1 - beta max over (y, k) of W(n_j . u^_k(x)),  W(t) = exp(kappa (t^2 - 1)),
%   so the directions into which a neighbouring lift continues cost less and the rest as before. Every lift
%   of every neighbour predicts, so the two bundles of a crossing are favoured each along its own continuation;
%   the misfit term still decides where the data point elsewhere. Deterministic: a weight of the penalty, not
%   a probability.
%
%   THE LIFTS. At each voxel, the local maxima of p(x, .) over the dictionary (neighbourhood
%   opts.nb degrees), refined to the continuous maximum by the fixed point of the gradient condition
%       u <- (sum_j a_j W(u.n_j)(u.n_j) n_j) / | . |,
%   kept when p >= opts.rel * max p and the lobe mass m = sum_{|u.n_j| >= cos cap} a_j / sum_j a_j
%   (an estimate of the volume fraction) is at least opts.mass_min, at least opts.sep degrees from a
%   stronger kept lift, at most opts.max_lifts, strongest first.
%
%   THE FORWARD MODEL, FITTED (opts.fit). Where two lobes of p overlap, the maximum of the weaker is
%   pulled toward the stronger. The lifts start a nonlinear least-squares fit of the model the data
%   obey, y_i = sum_k f_k E_k(g_i), E_k(g) = exp(-b (l_perp + dl (g.u_k)^2)), f_k >= 0 (Levenberg-
%   Marquardt, each u_k moved in its own normal plane); the fitted u_k replace the lifts and
%   f_k / sum f their mass. A stick whose mass falls below mass_min or that merges with another is dropped.
%   OFF by default: on the synthetic crossings (make_lift_crossing_nim, 32 directions, SNR 20) the fit of
%   each voxel's own noisy signal undoes the regularisation along the fibres and is worse than the
%   smoothed maxima (bundle A at 45 deg, f_B 0.7: 9.7 deg against 7.4 deg, 10.3 with the exact kernel).
%
%   CONTINUITY THROUGH LOST LOBES (opts.bridge, at most that many nodes). Where a lift's lobe is lost over
%   m <= bridge consecutive nodes along u and the lift continues beyond them (within 20 deg), the gap nodes
%   receive the mean direction of the two ends (L.filled marks them). Fibres do not end inside white matter.
%
%   OUTPUT  L.U [X Y Z P 3] unit lift directions, L.ok [X Y Z P], L.mass [X Y Z P], L.n [X Y Z],
%           L.filled [X Y Z P] (bridged), L.lambda [l_par l_perp], L.opts.
%
%   Options (defaults): ndirs 300, kappa 20, smooth 0, rel 0.1, mass_min 0.1, max_lifts 3, sep 15,
%   nb 12, cap 20, lambda [] (estimated), ref_fa_pct 90, mask [] (nim.mask & FA > 0, the voxels
%   whose lifts are wanted), workers [] (HINEC_MAX_WORKERS or 4; 0 = serial), fit false, fit_iter 25,
%   bridge 2, iso 'auto' (or a list of isotropic diffusivities, [] = none), sparsity 0.1, debias false,
%   smooth_kernel 'straight' ('curved'), smooth_kappa [0.15 0.3 0.5], smooth_normals 8,
%   continuation 0 (beta, off), continuation_stencil 'fd4'.
%
%   See also nim_field, nim_tractography_mmf, lsqnonneg.

if nargin < 2, opts = struct(); end
d = struct('ndirs', 300, 'kappa', 20, 'smooth', 0, 'rel', 0.1, 'mass_min', 0.1, 'max_lifts', 3, ...
           'sep', 15, 'nb', 12, 'cap', 20, 'lambda', [], 'ref_fa_pct', 90, 'mask', [], 'workers', [], ...
           'fit', false, 'fit_iter', 25, 'bridge', 2, 'iso', 'auto', 'sparsity', 0.1, 'debias', false, ...
           'smooth_kernel', 'straight', 'smooth_kappa', [0.15 0.3 0.5], 'smooth_normals', 8, ...
           'continuation', 0, 'continuation_stencil', 'fd4');
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = d.(fn{i}); end
end
dims = size(nim.FA);
if ~isfield(nim, 'img_bi') || ~isfield(nim, 'img_b0') || ~isfield(nim, 'bval') || ~isfield(nim, 'bvec')
    error('nim_lift_field:noSignal', 'field: lift needs the raw signal nim.img_b0, nim.img_bi, nim.bval, nim.bvec.');
end
dw = nim.bval(:) >= nim.thrsh_b0;
g  = double(nim.bvec(dw, :)); g = g ./ max(vecnorm(g, 2, 2), eps);
bb = double(nim.bval(dw)); bb = bb(:);
mask = opts.mask;
if isempty(mask)
    mask = logical(nim.FA > 0);
    if isfield(nim, 'mask') && ~isempty(nim.mask), mask = mask & nim.mask > 0.5; end
end
mask = logical(mask) & double(nim.img_b0) > 0;

% the kernel: the DTI eigenvalues of the single-fibre voxels
lam = opts.lambda;
if isempty(lam)
    ev  = reshape(double(nim.eval), [], 3);
    fa  = double(nim.FA(:)); tis = mask(:) & ev(:, 1) > 0;
    ref = tis & fa >= prctile(fa(tis), opts.ref_fa_pct);
    lam = [median(ev(ref, 1)), median(mean(ev(ref, 2:3), 2))];
end
dirs = halfsphere(opts.ndirs);
K = exp(-bb .* (lam(2) + (lam(1) - lam(2)) * (g * dirs') .^ 2));      % [N x ND]
iso = opts.iso;
if ischar(iso) || isstring(iso)
    if strcmpi(char(iso), 'auto'), iso = [lam(2) 3e-3]; else, iso = []; end   % 'none' = sticks only
end
Ki = exp(-bb .* iso(:)');                                             % [N x ni] isotropic compartments
Kall = [K Ki];
if opts.sparsity > 0   % the prior of few lifts: mu (sum_j a_j)^2, an extra row, isotropic compartments unpenalised
    Kall = [Kall; sqrt(opts.sparsity) * [ones(1, size(K, 2)) zeros(1, size(Ki, 2))]];
end
fprintf('lift field: %d directions at b = %s, kernel l_par %.3g l_perp %.3g mm^2/s, %d dictionary directions, %d voxels\n', ...
    nnz(dw), mat2str(unique(round(bb'))), lam(1), lam(2), size(dirs, 1), nnz(mask));

% the deconvolution, voxel by voxel (lsqnonneg), over a capped pool when one can be had
idx = find(mask);
S0  = double(nim.img_b0(idx));
Si  = reshape(double(nim.img_bi), [], size(nim.img_bi, 4));
Y   = Si(idx, :) ./ max(S0, eps);
Ya  = Y; if opts.sparsity > 0, Ya = [Y zeros(size(Y, 1), 1)]; end
ND  = size(dirs, 1); nv = numel(idx);
Acoef = zeros(nv, ND); Aiso = zeros(nv, numel(iso));
nw = opts.workers;
if isempty(nw)
    nw = str2double(getenv('HINEC_MAX_WORKERS')); if isnan(nw), nw = 4; end
end
lo = optimset('TolX', 1e-10);
t0 = tic;
Afull = deconvolve(K, Ki, Ya, [], opts.sparsity, nw, lo);
Acoef = Afull(:, 1:ND); Aiso = Afull(:, ND + 1:end); clear Afull
if opts.continuation > 0 && opts.sparsity <= 0
    warning('nim_lift_field:continuation', 'continuation needs sparsity > 0; it has no effect.');
elseif opts.continuation > 0                   % the continuation prior: a second pass, weighed by the neighbours
    Cw = continuation_weights(Acoef, dirs, idx, dims, opts);
    fprintf('lift field: continuation prior beta %.2f: penalty weight median %.2f, 5th pct %.2f\n', opts.continuation, ...
        median(Cw(:)), prctile(Cw(:), 5));
    Afull = deconvolve(K, Ki, Ya, Cw, opts.sparsity, nw, lo);
    Acoef = Afull(:, 1:ND); Aiso = Afull(:, ND + 1:end); clear Afull Cw
end
% the prior chooses the support, an unpenalised fit on it estimates (relaxed fit): the penalty would
% otherwise pull the lobes of a crossing toward each other (shrinkage), the support keeps out false lobes
made = false;
if opts.sparsity > 0 && opts.debias
    C = abs(dirs * dirs') >= cosd(25);
    sig = Acoef >= 0.02 * max(Acoef, [], 2) & Acoef > 0;
    sup = (double(sig) * double(C)) > 0;                             % atoms within 25 deg of a significant one
    Kd = [K Ki]; Ad = zeros(nv, ND + numel(iso)); ni = numel(iso);
    if nw > 0 && ~isempty(ver('parallel')) && nv > 5000
        if isempty(gcp('nocreate')), parpool('local', nw); made = true; end
        parfor (v = 1:nv, nw)
            cols = [find(sup(v, :)), ND + (1:ni)]; a = zeros(1, ND + ni);
            a(cols) = lsqnonneg(Kd(:, cols), Y(v, :)', lo)'; Ad(v, :) = a;
        end
        if made, delete(gcp('nocreate')); end
    else
        for v = 1:nv
            cols = [find(sup(v, :)), ND + (1:ni)]; a = zeros(1, ND + ni);
            a(cols) = lsqnonneg(Kd(:, cols), Y(v, :)', lo)'; Ad(v, :) = a;
        end
    end
    Acoef = Ad(:, 1:ND); Aiso = Ad(:, ND + 1:end); clear Ad sup sig
end
fprintf('lift field: deconvolution %.0f s; isotropic share of the signal weight: median %.2f\n', toc(t0), ...
    median(sum(Aiso, 2) ./ max(sum(Aiso, 2) + sum(Acoef, 2), eps)));

% regularisation along the fibres (in the directions' own frame = the voxel axes, see the header)
if opts.smooth > 0 && strcmpi(opts.smooth_kernel, 'curved')
    Acoef = smooth_curved(Acoef, dirs, idx, dims, mask, opts);
elseif opts.smooth > 0
    T  = opts.smooth; ts = -T:T; ws = T + 1 - abs(ts); ws = ws / sum(ws);
    [GX, GY, GZ] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
    Gp = [GX(idx) GY(idx) GZ(idx)];
    vol = zeros(dims);
    Fm = griddedInterpolant({1:dims(1), 1:dims(2), 1:dims(3)}, double(mask), 'linear', 'none');
    for j = 1:ND
        vol(:) = 0; vol(idx) = Acoef(:, j);
        F = griddedInterpolant({1:dims(1), 1:dims(2), 1:dims(3)}, vol, 'linear', 'none');
        acc = zeros(nv, 1); wacc = zeros(nv, 1);
        for q = 1:numel(ts)
            Xq = Gp + ts(q) * dirs(j, :); vq = F(Xq); tq = Fm(Xq);      % tq: share of the sample inside the tissue
            vq(isnan(vq)) = 0; tq(isnan(tq)) = 0;
            acc = acc + ws(q) * vq; wacc = wacc + ws(q) * tq;
        end
        Acoef(:, j) = acc ./ max(wacc, 1e-3);                          % renormalised: no bias at the tissue edge
    end
end

% the lifts: maxima of p over the dictionary, refined continuously, lobe mass, separation
[Uv, okv, Mv] = lifts_of(Acoef, dirs, opts);
P = opts.max_lifts;
% the forward model itself, fitted from the lifts (opts.fit): removes the pull of a strong lobe on a weak one
if opts.fit
    [Uv, Mv, okv] = fit_sticks(Y, g, bb, lam, Uv, Mv, okv, opts);
end
L.U = zeros([dims P 3], 'single'); L.ok = false([dims P]); L.mass = zeros([dims P], 'single');
for s = 1:P
    for c = 1:3
        tmp = zeros(dims, 'single'); tmp(idx) = Uv(:, s, c); L.U(:, :, :, s, c) = tmp;
    end
    tmp = false(dims); tmp(idx) = okv(:, s); L.ok(:, :, :, s) = tmp;
    tmp = zeros(dims, 'single'); tmp(idx) = Mv(:, s); L.mass(:, :, :, s) = tmp;
end
L.filled = false([dims P]);
if opts.bridge > 0
    [L.U, L.ok, L.mass, L.filled] = bridge_gaps(L.U, L.ok, L.mass, mask, opts.bridge, opts);
end
L.n = uint8(sum(L.ok, 4));
L.lambda = lam; L.opts = opts; L.opts.mask = [];
fprintf('lift field: %.1f%% of the voxels with one lift, %.1f%% with two, %.1f%% with three or more\n', ...
    100 * mean(L.n(idx) == 1), 100 * mean(L.n(idx) == 2), 100 * mean(L.n(idx) >= 3));
end

function Afull = deconvolve(K, Ki, Ya, Cw, mu, nw, lo)
% lsqnonneg voxel by voxel of [K Ki] a = y with the prior row sqrt(mu) [c 0] (c = 1, or the rows of Cw), over a
% capped pool when one can be had.
nv = size(Ya, 1); ND = size(K, 2); ni = size(Ki, 2); Kb = [K Ki];
Afull = zeros(nv, ND + ni); made = false;
if isempty(Cw), Kall = Kb; if mu > 0, Kall = [Kb; sqrt(mu) * [ones(1, ND) zeros(1, ni)]]; end, end
use_pool = nw > 0 && ~isempty(ver('parallel')) && nv > 5000;
if use_pool && isempty(gcp('nocreate')), parpool('local', nw); made = true; end
if ~use_pool && nw > 0, fprintf('lift field: no parallel pool, the deconvolution runs serially\n'); end
if isempty(Cw)
    if use_pool
        parfor (v = 1:nv, nw), Afull(v, :) = lsqnonneg(Kall, Ya(v, :)', lo)'; end      % at most nw workers
    else
        for v = 1:nv, Afull(v, :) = lsqnonneg(Kall, Ya(v, :)', lo)'; end
    end
else
    sq = sqrt(mu); z = zeros(1, ni);
    if use_pool
        parfor (v = 1:nv, nw), Afull(v, :) = lsqnonneg([Kb; sq * [Cw(v, :) z]], Ya(v, :)', lo)'; end
    else
        for v = 1:nv, Afull(v, :) = lsqnonneg([Kb; sq * [Cw(v, :) z]], Ya(v, :)', lo)'; end
    end
end
if made, delete(gcp('nocreate')); end
end

function Cw = continuation_weights(A, dirs, idx, dims, o)
% c_j(x) = 1 - beta max over the 6 neighbours y and their lifts k of W(n_j . u^_k(x)) (see THE CONTINUATION PRIOR).
[Uv, okv] = lifts_of(A, dirs, o); P = size(okv, 2); N = prod(dims); nv = numel(idx);
U = zeros([dims P 3]); ok = false([dims P]);
for s = 1:P
    for c = 1:3, tmp = zeros(dims); tmp(idx) = Uv(:, s, c); U(:, :, :, s, c) = tmp; end
    tmp = false(dims); tmp(idx) = okv(:, s); ok(:, :, :, s) = tmp;
end
[~, ~, ~, F] = nim_lift_connection(U, ok, 0.125, 30, o.continuation_stencil, 'corner');
[I, J, K] = ind2sub(dims, idx); sub = [I J K];
best = zeros(nv, size(dirs, 1), 'single'); chunk = 20000;
for m = 1:3
    for sg = [-1 1]
        y = sub; y(:, m) = y(:, m) + sg; in = y(:, m) >= 1 & y(:, m) <= dims(m);
        ly = ones(nv, 1); ly(in) = sub2ind(dims, y(in, 1), y(in, 2), y(in, 3));
        for k = 1:P
            u = reshape(F.U(ly, k, :), nv, 3); has = in & F.ok(ly, k);
            du = zeros(nv, 3);                                       % du_b/dx_m = sum_a u_a T(a, b, m)
            for b = 1:3, for a = 1:3
                du(:, b) = du(:, b) + u(:, a) .* F.Tn(ly, k, sub2ind([3 3 3], a, b, m));
            end, end
            uh = u - sg * du; uh = uh ./ max(vecnorm(uh, 2, 2), eps);   % x - y = -sg e_m
            for c0 = 1:chunk:nv
                r = c0:min(c0 + chunk - 1, nv); h = has(r);
                Wq = single(exp(o.kappa * ((uh(r, :) * dirs') .^ 2 - 1))); Wq(~h, :) = 0;
                best(r, :) = max(best(r, :), Wq);
            end
        end
    end
end
Cw = double(1 - o.continuation * best);
end

function [U, ok, Mass] = lifts_of(A, dirs, o)
% Lifts of every row of A (vectorised over voxels; see the header for the definitions).
nv = size(A, 1); ND = size(dirs, 1); P = o.max_lifts;
C  = abs(dirs * dirs');
W  = exp(o.kappa * (C .^ 2 - 1));
p  = A * W;                                                   % [nv x ND]
nbrs = arrayfun(@(m) find(C(m, :) >= cosd(o.nb) & (1:ND) ~= m), 1:ND, 'UniformOutput', false);
nmax = max(cellfun(@numel, nbrs));
NB = ones(ND, nmax); for m = 1:ND, NB(m, 1:numel(nbrs{m})) = nbrs{m}; NB(m, numel(nbrs{m})+1:end) = m; end
pn = zeros(nv, ND);
for m = 1:ND, pn(:, m) = max(p(:, NB(m, :)), [], 2); end      % largest neighbour value
mx = max(p, [], 2); tot = sum(A, 2);
cand = p >= pn & p >= o.rel * mx & mx > 0;                    % local maxima above rel
pc = p; pc(~cand) = -Inf;
[~, order] = sort(pc, 2, 'descend');
U = zeros(nv, P, 3); ok = false(nv, P); Mass = zeros(nv, P); nk = zeros(nv, 1);
for r = 1:min(ND, 8)                                          % the strongest candidates, in order
    m = order(:, r); live = isfinite(pc(sub2ind([nv ND], (1:nv)', m))) & nk < P;
    if ~any(live), break; end
    u = refine(A(live, :), dirs, dirs(m(live), :), o.kappa);
    vi = find(live);
    sepok = true(numel(vi), 1);
    for s = 1:P
        prev = reshape(U(vi, s, :), [], 3);
        sepok = sepok & ~(ok(vi, s) & abs(sum(prev .* u, 2)) >= cosd(o.sep));
    end
    lobe = (abs(u * dirs') >= cosd(o.cap));
    mass = sum(A(vi, :) .* lobe, 2) ./ max(tot(vi), eps);
    keep = sepok & mass >= o.mass_min;
    for q = find(keep)'
        v = vi(q); s = nk(v) + 1;
        U(v, s, :) = reshape(u(q, :), 1, 1, 3); ok(v, s) = true; Mass(v, s) = mass(q); nk(v) = s;
    end
end
end

function [U, Mass, ok] = fit_sticks(Y, g, b, lam, U, Mass, ok, o)
% Nonlinear least squares of the forward model from the lifts, voxels with the same number of lifts
% together (Levenberg-Marquardt, analytic Jacobian):
%     y_i = sum_k q_k^2 E_k(g_i),  E_k(g) = exp(-b (l_perp + dl (g.u_k)^2)),
% each u_k moved in its own normal plane, u_k <- (u_k + d2 e2_k + d3 e3_k)/|.|, q_k^2 >= 0 its weight.
% The fitted directions replace the lifts and q_k^2 / sum q^2 their mass; a fit that does not lower the
% residual leaves the voxel as it was.
dl = lam(1) - lam(2); nk = sum(ok, 2); N = size(g, 1);
for K = 1:size(U, 2)
    v = find(nk == K); if isempty(v), continue; end
    n = numel(v); y = Y(v, :)';                                         % [N x n]
    u = U(v, 1:K, :); q = sqrt(max(Mass(v, 1:K), 1e-3));                % [n K 3], [n K]
    [~, Ek] = resid(u, q); sc = reshape(sum(Ek .* reshape(q', [1 K n]) .^ 2, 2), N, n);   % scale q to the signal
    q = q .* sqrt(max((sum(y .* sc, 1) ./ max(sum(sc .^ 2, 1), eps))', eps));
    [r0, Ek, c] = resid(u, q); cost = sum(r0 .^ 2, 1)'; cost0 = cost; mu = 1e-2 * ones(n, 1);
    for it = 1:o.fit_iter
        [e2, e3] = normal_plane(u);
        J = zeros(N, 3 * K, n);
        for k = 1:K
            Ekk = reshape(Ek(:, k, :), N, n); ck = reshape(c(:, k, :), N, n); qk = q(:, k)';
            J(:, k, :)         = reshape(2 * qk .* Ekk, N, 1, n);
            dE = -2 * b .* dl .* ck .* Ekk .* qk .^ 2;                   % dy/d(g.u_k)
            J(:, K + 2*k - 1, :) = reshape(dE .* (g * reshape(e2(:, k, :), n, 3)'), N, 1, n);
            J(:, K + 2*k, :)     = reshape(dE .* (g * reshape(e3(:, k, :), n, 3)'), N, 1, n);
        end
        H  = pagemtimes(J, 'transpose', J, 'none'); gr = pagemtimes(J, 'transpose', reshape(r0, N, 1, n), 'none');
        Hd = H + reshape(mu, 1, 1, n) .* (H .* eye(3 * K)) + 1e-12 * eye(3 * K);
        dlt = -reshape(pagemldivide(Hd, gr), 3 * K, n);                 % [3K x n]
        qn = q + dlt(1:K, :)';
        un = u;
        for k = 1:K
            w = reshape(u(:, k, :), n, 3) + dlt(K + 2*k - 1, :)' .* reshape(e2(:, k, :), n, 3) + dlt(K + 2*k, :)' .* reshape(e3(:, k, :), n, 3);
            un(:, k, :) = reshape(w ./ max(vecnorm(w, 2, 2), eps), n, 1, 3);
        end
        [rn, En, cn] = resid(un, qn); cn_ = sum(rn .^ 2, 1)';
        acc = cn_ < cost;
        u(acc, :, :) = un(acc, :, :); q(acc, :) = qn(acc, :); cost(acc) = cn_(acc);
        r0(:, acc) = rn(:, acc); Ek(:, :, acc) = En(:, :, acc); c(:, :, acc) = cn(:, :, acc);
        mu(acc) = mu(acc) / 3; mu(~acc) = mu(~acc) * 4;
    end
    better = cost < cost0;
    f = q .^ 2; f = f ./ max(sum(f, 2), eps);
    U(v(better), 1:K, :) = u(better, :, :); Mass(v(better), 1:K) = f(better, :);
    % a fitted stick whose weight fell below mass_min, or that merged into another (sep), is dropped
    for k = 1:K
        drop = better & Mass(v, k) < o.mass_min;
        for j = 1:k - 1
            drop = drop | (better & abs(sum(reshape(U(v, k, :), [], 3) .* reshape(U(v, j, :), [], 3), 2)) >= cosd(o.sep));
        end
        ok(v(drop), k) = false;
    end
end
    function [r, E, c] = resid(u, q)
        nn = size(u, 1); KK = size(u, 2);
        c = zeros(N, KK, nn); for kk = 1:KK, c(:, kk, :) = reshape(g * reshape(u(:, kk, :), nn, 3)', N, 1, nn); end
        E = exp(-b .* (lam(2) + dl * c .^ 2));
        r = reshape(sum(E .* reshape(q', [1 KK nn]) .^ 2, 2), N, nn) - y;
    end
end

function A = smooth_curved(A, dirs, idx, dims, mask, o)
% The regularisation along the fibres on curved arcs (opts.smooth_kernel 'curved'; see the header).
% For slice j (direction n), an arc of curvature k in the normal direction m = cos(q) e2 + sin(q) e3 is
%   gamma(t) = x + (sin(k t) / k) n + ((1 - cos(k t)) / k) m,   its tangent d(t) = cos(k t) n + sin(k t) m,
% (k = 0: gamma(t) = x + t n, d = n), t = -T..T. Its support is the weighted mean of the density of the
% dictionary direction nearest d(t), trilinear at gamma(t), renormalised by the tissue share as the straight
% kernel; the slice keeps the support of its best arc: a~_j(x) = max over (k, q) of that mean.
T = o.smooth; ts = -T:T; ws = T + 1 - abs(ts); ws = ws / sum(ws);
kap = o.smooth_kappa(:)'; nq = o.smooth_normals; ND = size(dirs, 1); nv = numel(idx);
pad = ceil(T) + 2; pd = dims + 2 * pad;                       % padded grid: every arc point stays inside
[I, J, K] = ind2sub(dims, idx); lp = sub2ind(pd, I + pad, J + pad, K + pad);
V = zeros(prod(pd), ND, 'single'); V(lp, :) = single(A);
Mk = zeros(pd, 'single'); Mk(pad + 1:pad + dims(1), pad + 1:pad + dims(2), pad + 1:pad + dims(3)) = single(mask);
Mk = Mk(:); st = [1, pd(1), pd(1) * pd(2)];
[e2, e3] = normal_plane(reshape(dirs, ND, 1, 3)); e2 = reshape(e2, ND, 3); e3 = reshape(e3, ND, 3);
arcs = [0 0; [repelem(kap(:), nq), repmat((0:nq - 1)' * 2 * pi / nq, numel(kap), 1)]];
out = zeros(nv, ND);
for j = 1:ND
    n = dirs(j, :); best = -inf(nv, 1);
    for c = 1:size(arcs, 1)
        k = arcs(c, 1); m = cos(arcs(c, 2)) * e2(j, :) + sin(arcs(c, 2)) * e3(j, :);
        acc = zeros(nv, 1); wacc = zeros(nv, 1);
        for q = 1:numel(ts)
            t = ts(q);
            if t == 0, acc = acc + ws(q) * double(V(lp, j)); wacc = wacc + ws(q); continue; end
            if k == 0, dp = t * n; dt = n;
            else, dp = (sin(k * t) / k) * n + ((1 - cos(k * t)) / k) * m; dt = cos(k * t) * n + sin(k * t) * m; end
            [~, jp] = max(abs(dirs * dt'));
            f = floor(dp); r = dp - f; v = zeros(nv, 1); tq = zeros(nv, 1);
            for cx = 0:1, for cy = 0:1, for cz = 0:1
                wc = (cx * r(1) + (1 - cx) * (1 - r(1))) * (cy * r(2) + (1 - cy) * (1 - r(2))) * (cz * r(3) + (1 - cz) * (1 - r(3)));
                if wc == 0, continue; end
                li = lp + (f + [cx cy cz]) * st';
                v = v + wc * double(V(li, jp)); tq = tq + wc * double(Mk(li));
            end, end, end
            acc = acc + ws(q) * v; wacc = wacc + ws(q) * tq;
        end
        best = max(best, acc ./ max(wacc, 1e-3));
    end
    out(:, j) = best;
end
A = out;
end

function [e2, e3] = normal_plane(u)
% An orthonormal basis of the plane normal to each u(:, k, :) (reference-axis projection).
[n, K, ~] = size(u); e2 = zeros(n, K, 3); e3 = zeros(n, K, 3);
for k = 1:K
    a = reshape(u(:, k, :), n, 3);
    [~, ax] = min(abs(a), [], 2); r = zeros(n, 3); r(sub2ind([n 3], (1:n)', ax)) = 1;
    b2 = r - sum(r .* a, 2) .* a; b2 = b2 ./ max(vecnorm(b2, 2, 2), eps);
    e2(:, k, :) = reshape(b2, n, 1, 3); e3(:, k, :) = reshape(cross(a, b2, 2), n, 1, 3);
end
end

function [U, ok, M, filled] = bridge_gaps(U, ok, M, mask, gap, o)
% Continuity of a lift through nodes where its lobe is lost (fibres do not end inside white matter).
% From every lift u at node x, walk t = 1..gap+1 nodes along +-u (nearest nodes): when the nodes at t = 1..m
% (m <= gap) are in the mask and carry no lift within o.sep of u, and the node at t = m + 1 carries one
% within 20 deg of u, the m gap nodes receive the mean direction of the two ends (filled = true) and the
% smaller of their masses, unless a lift within o.sep is already there (no duplicates).
dims = size(ok, 1:3); P = size(ok, 4); nv = prod(dims);
filled = false(size(ok));
comp = @(l, q, c) l + (q - 1) * nv + (c - 1) * nv * P;
    function v = liftAt(l, q), v = double([U(comp(l, q, 1)), U(comp(l, q, 2)), U(comp(l, q, 3))]); end
    function q = nearLift(l, u, cosmin)
        q = 0; best = cosmin;
        for qq = 1:P
            if ok(l + (qq - 1) * nv)
                d = abs(liftAt(l, qq) * u'); if d >= best, best = d; q = qq; end
            end
        end
    end
for s = 1:P
    src = find(ok(:, :, :, s) & ~filled(:, :, :, s));
    for ii = 1:numel(src)
        l0 = src(ii); u = liftAt(l0, s); [x, y, z] = ind2sub(dims, l0); x0 = [x y z];
        for sg = [-1 1]
            path = zeros(gap, 1); m = 0; done = false;
            for t = 1:gap + 1
                xt = round(x0 + sg * t * u);
                if any(xt < 1 | xt > dims), break; end
                lt = sub2ind(dims, xt(1), xt(2), xt(3));
                if ~mask(lt), break; end
                if nearLift(lt, u, cosd(o.sep)) > 0
                    if m > 0
                        qe = nearLift(lt, u, cosd(20));
                        if qe > 0
                            v = liftAt(lt, qe); v = v * sign(v * u'); w = (u + v) / norm(u + v);
                            for k = 1:m
                                lk = path(k);
                                if nearLift(lk, w, cosd(o.sep)) > 0, continue; end
                                slot = find(~ok(lk + (0:P - 1) * nv), 1);
                                if isempty(slot), continue; end
                                U(comp(lk, slot, 1)) = w(1); U(comp(lk, slot, 2)) = w(2); U(comp(lk, slot, 3)) = w(3);
                                ok(lk + (slot - 1) * nv) = true; filled(lk + (slot - 1) * nv) = true;
                                M(lk + (slot - 1) * nv) = min(M(l0 + (s - 1) * nv), M(lt + (qe - 1) * nv));
                            end
                        end
                    end
                    done = true; break
                end
                if t > gap, break; end
                m = m + 1; path(m) = lt;
            end
            if ~done, continue; end
        end
    end
end
end

function u = refine(A, dirs, u, kappa)
% Continuous maximum of p(u) = sum_j a_j W(u.n_j) from u: the fixed point of the gradient condition.
for it = 1:30
    t = u * dirs';                                            % [n x ND]
    v = (A .* exp(kappa * (t .^ 2 - 1)) .* t) * dirs;
    v = v ./ max(vecnorm(v, 2, 2), eps);
    done = abs(sum(v .* u, 2)) > 1 - 1e-10;
    u = v;
    if all(done), break; end
end
end

function v = halfsphere(n)
% Fibonacci points of the sphere, the half z > 0: about n near-uniform directions of RP^2.
k = (0:2 * n - 1)' + 0.5;
z = 1 - k / n; r = sqrt(max(1 - z .^ 2, 0)); t = pi * (1 + sqrt(5)) * k;
v = [r .* cos(t), r .* sin(t), z];
v = v(v(:, 3) > 0, :);
end
