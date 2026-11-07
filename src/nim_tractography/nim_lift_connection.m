function [w, mf, conf, F] = nim_lift_connection(U, ok, read_min, read_tol, stencil, interp)
%NIM_LIFT_CONNECTION  The nodal frames and connection 1-forms of the lifts, as heading-aware handles
%                     w(X, T), mf(X, T) for advanceFrame (field: lift; see nim_lift_field).
%
%   U [X Y Z P 3] the unit direction of up to P lifts at the voxel nodes, ok [X Y Z P] where each exists.
%
%   THE BRANCH. The lift through node x with direction u continues to the neighbour nodes along
%   axis k as the branch of the implicit function u(x): at x +- e_k the neighbour's lift closest to u,
%   |u.v| >= cos 30 deg (a branch ends where there is none), and at x +- 2 e_k the lift closest to that
%   match. Its derivative is the least-squares slope of the dyadic M = u u' over the offsets reached,
%       dM/dx_k = sum_t t (M_t - M_0) / sum_t t^2,   t in {-2, -1, 1, 2},
%   and dU/dx_k = (dM/dx_k) u (u . du = 0). With stencil 'fd4' the derivative is instead that of the
%   nodal connection of field: dti (nim_tractography_mmf): (M_-2 - 8 M_-1 + 8 M_+1 - M_+2) / 12 where the
%   branch reaches +-2 on both sides, (M_+1 - M_-1) / 2 where it reaches +-1, the slope elsewhere. The slope
%   ('ls', the default) is a smoother derivative and reads less curvature on bends. The nodal connection tensor T(a,b,k) = u_a dU_b/dx_k is
%   free of the sign of u, as in nim_tractography_mmf's nodal connection (MMFtractography_notes Sec. 7).
%
%   THE QUERY. At a point X with heading T, at each of the 8 nodes around X the lift closest to T,
%   |u.T| >= cos 30 deg (a node without one contributes nothing), is read; its dyadic and T are
%   interpolated with the trilinear weights renormalised over the corners read (no value where they carry
%   less than read_min of the weight, default 0.25: inside a crossing the weaker bundle's lift is missing at
%   scattered nodes, and its continuation through them is read from the nodes that carry it), e1
%   is the principal axis of the interpolated dyadic, signed along T, e2 the projection of the
%   Cartesian axis least aligned with e1 onto the normal plane, e3 = e1 x e2, and
%       w_1j,k = e_j . sum_a e1_a T(a, b, k),   j = 2, 3   (the twist w_23 is zero: 'Twist' 'zero').
%   This is the lift of the R^3 scheme to M = Omega x RP^2: P = 1 with U = v1 is the single-valued field.
%
%   conf(X, T) is the read weight itself, in [0, 1]: the share of the trilinear weight on the corners that
%   carry a lift within read_tol (default 30) degrees of T. It is the confidence of the hold (advanceFrame
%   'Confidence'): below mmf.hold_min the tract's own lift is absent and the tract coasts on its held
%   connection instead of reading another bundle's lift or stopping.
%
%   THE LOCAL SHEET (interp 'sheet_linear' | 'sheet_cubic'; 'corner', the default, is the query above).
%   Over a small neighbourhood V of X the lifts split into sheets, each the graph of one direction field over
%   V. The query reads the sheet of the heading on the 4 x 4 x 4 nodes x_q = floor(X) + o, o in {-1..2}^3, and
%   interpolates that single-valued nodal field as the scheme of field: dti does (nim_tractography_mmf):
%     - the start node is the inner corner nearest X, its lift the one closest to T (|u.T| >= cos read_tol);
%     - every other node, in breadth-first order over the 6-neighbourhood from the start, takes the lift
%       closest to its parent's, |u.u_parent| >= cos 30 deg, signed along it; a node without one is a
%       gap of the sheet and continues its parent's direction and nodal tensor (as fill_zone continues the
%       field into the zone from the known nodes around it);
%     - the dyadic and the nodal tensor T of the sheet are interpolated with the separable weights
%       w(X) = prod_k phi(r_k - o_k), r = X - floor(X): 'sheet_cubic' phi = the Keys cubic convolution
%       kernel (a = -1/2; interpolation.method: cubic, C1, exact for quadratics), 'sheet_linear' phi = the
%       hat function on the inner corners (trilinear);
%     - conf(X, T) = the trilinear weight on the inner corners where the sheet is present (not a gap).
%   Reading one sheet removes the corner-wise choice that made the field jump between the lifts of a crossing.
%
%   F (fourth output) holds the nodal quantities: F.U [N x P x 3] (unit, zero where absent), F.ok [N x P],
%   F.Tn [N x P x 27] the nodal tensor T(a, b, k) at linear index sub2ind([3 3 3], a, b, k), N = prod(dims), so that
%   dU_b/dx_k = sum_a U_a T(a, b, k) (nim_lift_field's continuation prior transports the lifts with it).
%
%   See also nim_lift_field, nim_tractography_mmf, advanceFrame.
if nargin < 5 || isempty(stencil), stencil = 'ls'; end
if nargin < 6 || isempty(interp), interp = 'corner'; end
stencil = lower(char(string(stencil))); interp = lower(char(string(interp)));
if ~any(strcmp(interp, {'corner', 'sheet_linear', 'sheet_cubic'}))
    error('nim_lift_connection:interp', 'interp must be corner, sheet_linear or sheet_cubic, not ''%s''.', interp);
end
dims = size(ok, 1:3); P = size(ok, 4); N = prod(dims);
Uf = reshape(U, N, P, 3); okf = reshape(ok, N, P);
Uf = Uf ./ max(vecnorm(Uf, 2, 3), eps); Uf(repmat(~okf, [1 1 3])) = 0;
[I, J, K] = ind2sub(dims, (1:N)'); sub = [I J K];
Tn = zeros(N, P, 27);
% The derivative of each lift along x_k: least-squares slope of the dyadic over the matched branch at
% t = -2, -1, 1, 2 (the +-2 node matched to the +-1 match), dM/dx_k = sum_t t (M_t - M_0) / sum_t t^2 over the
% offsets reached; a branch is broken where no neighbour lift is within 30 deg.
for k = 1:3
    for s = 1:P
        u0 = reshape(Uf(:, s, :), N, 3); M0 = reshape(u0, N, 3, 1) .* reshape(u0, N, 1, 3);
        num = zeros(N, 3, 3); den = zeros(N, 1); Mb = cell(2, 2); ab = cell(2, 2);
        for sg = [-1 1]
            cur = u0; curli = (1:N)'; alive = okf(:, s);
            for t = 1:2
                nb = sub(curli, :); nb(:, k) = nb(:, k) + sg; in = nb(:, k) >= 1 & nb(:, k) <= dims(k);
                li = ones(N, 1); li(in) = sub2ind(dims, nb(in, 1), nb(in, 2), nb(in, 3));
                best = zeros(N, 1); bu = zeros(N, 3);
                for q = 1:P
                    v = reshape(Uf(li, q, :), N, 3); d = abs(sum(cur .* v, 2)); d(~in | ~okf(li, q)) = 0;
                    pick = d > best; best(pick) = d(pick); bu(pick, :) = v(pick, :) .* sign(sum(cur(pick, :) .* v(pick, :), 2));
                end
                alive = alive & best >= cosd(30);
                Mt = reshape(bu, N, 3, 1) .* reshape(bu, N, 1, 3);
                num(alive, :, :) = num(alive, :, :) + (sg * t) * (Mt(alive, :, :) - M0(alive, :, :));
                den(alive) = den(alive) + t ^ 2;
                Mb{(sg + 3) / 2, t} = Mt; ab{(sg + 3) / 2, t} = alive;
                cur = bu; curli = li;
            end
        end
        G = num ./ max(den, 1);
        if strcmp(stencil, 'fd4')
            % the stencil of the nodal connection of field: dti (nim_tractography_mmf, nodal_jacobian):
            % fourth-order central differences (M_-2 - 8 M_-1 + 8 M_+1 - M_+2) / 12 along the matched branch,
            % second order (M_+1 - M_-1) / 2 where the branch reaches only +-1, the slope above elsewhere
            a4 = ab{1, 1} & ab{1, 2} & ab{2, 1} & ab{2, 2}; a2 = ab{1, 1} & ab{2, 1} & ~a4;
            G(a4, :, :) = (Mb{1, 2}(a4, :, :) - 8 * Mb{1, 1}(a4, :, :) + 8 * Mb{2, 1}(a4, :, :) - Mb{2, 2}(a4, :, :)) / 12;
            G(a2, :, :) = (Mb{2, 1}(a2, :, :) - Mb{1, 1}(a2, :, :)) / 2;
        end
        Jk = squeeze(sum(G .* reshape(u0, N, 1, 3), 3));          % dU/dx_k = (dM/dx_k) u
        for a = 1:3, for b = 1:3
            Tn(:, s, sub2ind([3 3 3], a, b, k)) = u0(:, a) .* Jk(:, b);
        end, end
    end
end
if nargin < 3 || isempty(read_min), read_min = 0.25; end
if nargin < 4 || isempty(read_tol), read_tol = 30; end
F = struct('read_min', read_min, 'cos_tol', cosd(read_tol), 'dims', dims, 'P', P, 'U', Uf, 'ok', okf, 'Tn', Tn);
if strcmp(interp, 'corner')
    w  = @(X, T) query(X, T, 'w', F);
    mf = @(X, T) query(X, T, 'mf', F);
    conf = @(X, T) query(X, T, 'conf', F);
else
    F.cubic = strcmp(interp, 'sheet_cubic'); F.trees = sheet_trees(); F.Tf = reshape(Tn, N * P, 27);
    w  = @(X, T) query_sheet(X, T, 'w', F);
    mf = @(X, T) query_sheet(X, T, 'mf', F);
    conf = @(X, T) query_sheet(X, T, 'conf', F);
end
end

function trees = sheet_trees()
% For each of the 8 inner corners of the 4 x 4 x 4 window (offsets {-1..2}^3), the breadth-first order of the
% 64 nodes over the 6-neighbourhood from that corner and each node's parent.
[a, b, c] = ndgrid(-1:2, -1:2, -1:2); O = [a(:) b(:) c(:)];
trees = struct('O', O, 'order', {cell(8, 1)}, 'parent', {cell(8, 1)}, 'start', zeros(8, 1));
for s = 1:8
    c0 = [bitand(s - 1, 1), bitand(bitshift(s - 1, -1), 1), bitand(bitshift(s - 1, -2), 1)];
    q0 = find(all(O == c0, 2)); par = zeros(64, 1); seen = false(64, 1); seen(q0) = true; order = q0; head = 1;
    while head <= numel(order)
        q = order(head); head = head + 1;
        for d = [eye(3); -eye(3)]'
            m = find(all(O == O(q, :) + d', 2));
            if ~isempty(m) && ~seen(m), seen(m) = true; par(m) = q; order(end + 1) = m; end %#ok<AGROW>
        end
    end
    trees.order{s} = order; trees.parent{s} = par; trees.start(s) = q0;
end
end

function out = query_sheet(X, T, what, F)
% The local sheet of the heading on the 4 x 4 x 4 nodes around X (see THE LOCAL SHEET in the header).
dims = F.dims; P = F.P; n = size(X, 1); N = prod(dims);
T = T(:, 1:3) ./ max(vecnorm(T(:, 1:3), 2, 2), eps);
f = floor(X); r = X - f; O = F.trees.O;
c0 = round(r); sid = 1 + c0(:, 1) + 2 * c0(:, 2) + 4 * c0(:, 3);      % the inner corner nearest X
% separable weights of the 64 nodes
if F.cubic
    phi = @(s) (abs(s) <= 1) .* (1.5 * abs(s) .^ 3 - 2.5 * abs(s) .^ 2 + 1) + ...
               (abs(s) > 1 & abs(s) < 2) .* (-0.5 * abs(s) .^ 3 + 2.5 * abs(s) .^ 2 - 4 * abs(s) + 2);
else
    phi = @(s) max(1 - abs(s), 0);
end
Wt = ones(n, 64); Wl = ones(n, 64);
for k = 1:3
    Wt = Wt .* phi(r(:, k) - O(:, k)');
    Wl = Wl .* max(1 - abs(r(:, k) - O(:, k)'), 0);                 % trilinear: the read weight
end
Usel = zeros(n, 64, 3); src = zeros(n, 64); present = false(n, 64); valid = true(n, 1);
for s = 1:8
    qi = find(sid == s); if isempty(qi), continue; end
    m = numel(qi); order = F.trees.order{s}; par = F.trees.parent{s};
    for j = 1:64
        q = order(j); c = f(qi, :) + O(q, :); in = all(c >= 1 & c <= dims, 2);
        li = ones(m, 1); li(in) = sub2ind(dims, c(in, 1), c(in, 2), c(in, 3));
        if j == 1, ref = T(qi, :); tol = F.cos_tol; else, ref = reshape(Usel(qi, par(q), :), m, 3); tol = cosd(30); end
        best = zeros(m, 1); bs = zeros(m, 1); bu = zeros(m, 3);
        for p = 1:P
            v = reshape(F.U(li, p, :), m, 3); d = sum(v .* ref, 2); al = abs(d); al(~in | ~F.ok(li, p) | al < tol) = -1;
            pick = al > best; best(pick) = al(pick); bs(pick) = p; bu(pick, :) = v(pick, :) .* sign(d(pick, :) + 1e-300);
        end
        has = bs > 0;
        if j == 1, valid(qi) = has; end
        present(qi(has), q) = true; src(qi(has), q) = (bs(has) - 1) * N + li(has);
        Usel(qi(has), q, :) = reshape(bu(has, :), [], 1, 3);
        if j > 1                                                    % a gap continues its parent
            g = ~has; Usel(qi(g), q, :) = Usel(qi(g), par(q), :); src(qi(g), q) = src(qi(g), par(q));
        end
    end
end
wsum = sum(Wl .* present, 2); wsum(~valid) = 0;
if strcmp(what, 'conf'), out = wsum; return; end
bad = ~valid | wsum < F.read_min;
M = zeros(n, 3, 3); TT = zeros(n, 27); ok = valid;
for q = 1:64
    u = reshape(Usel(:, q, :), n, 3); wq = Wt(:, q);
    M = M + wq .* (reshape(u, n, 3, 1) .* reshape(u, n, 1, 3));
    sq = src(:, q); h = ok & sq > 0;
    TT(h, :) = TT(h, :) + wq(h, :) .* F.Tf(sq(h), :);
end
e1 = T;
for it = 1:25
    e1 = reshape(sum(M .* reshape(e1, n, 1, 3), 3), n, 3); e1 = e1 ./ max(vecnorm(e1, 2, 2), eps);
end
e1 = e1 .* sign(sum(e1 .* T, 2) + 1e-300);
[~, ax] = min(abs(e1), [], 2); ref = zeros(n, 3); ref(sub2ind([n 3], (1:n)', ax)) = 1;
e2 = ref - sum(ref .* e1, 2) .* e1; e2 = e2 ./ max(vecnorm(e2, 2, 2), eps); e3 = cross(e1, e2, 2);
if strcmp(what, 'mf')
    out = cat(2, reshape(e1, n, 1, 3), reshape(e2, n, 1, 3), reshape(e3, n, 1, 3)); out(bad, :, :) = NaN; return
end
TT = reshape(TT, n, 3, 3, 3);
out = zeros(n, 3, 3, 3);
for k = 1:3
    de1 = zeros(n, 3);
    for b = 1:3, for a = 1:3, de1(:, b) = de1(:, b) + e1(:, a) .* TT(:, a, b, k); end, end
    w12 = sum(e2 .* de1, 2); w13 = sum(e3 .* de1, 2);
    out(:, 1, 2, k) = w12; out(:, 2, 1, k) = -w12; out(:, 1, 3, k) = w13; out(:, 3, 1, k) = -w13;
end
out(bad, :, :, :) = NaN;
end

function out = query(X, T, what, F)
dims = F.dims; P = F.P; n = size(X, 1);
T = T(:, 1:3) ./ max(vecnorm(T(:, 1:3), 2, 2), eps);
f = floor(X); r = X - f;
M = zeros(n, 3, 3); TT = zeros(n, 27); wsum = zeros(n, 1);
for dx = 0:1, for dy = 0:1, for dz = 0:1
    c = f + [dx dy dz]; in = all(c >= 1 & c <= dims, 2);
    li = ones(n, 1); li(in) = sub2ind(dims, c(in, 1), c(in, 2), c(in, 3));
    wt = (dx * r(:, 1) + (1 - dx) * (1 - r(:, 1))) .* (dy * r(:, 2) + (1 - dy) * (1 - r(:, 2))) .* (dz * r(:, 3) + (1 - dz) * (1 - r(:, 3)));
    best = zeros(n, 1); bs = zeros(n, 1);
    for s = 1:P
        al = abs(sum(reshape(F.U(li, s, :), n, 3) .* T, 2)); al(~in | ~F.ok(li, s)) = -1;
        al(al < F.cos_tol) = -1;                              % a corner without a lift near the heading gives none
        pick = al > best; best(pick) = al(pick); bs(pick) = s;
    end
    has = bs > 0; wt(~has) = 0;
    for s = 1:P
        q = has & bs == s; if ~any(q), continue; end
        us = reshape(F.U(li(q), s, :), [], 3); us = us .* sign(sum(us .* T(q, :), 2) + 1e-300);
        M(q, :, :) = M(q, :, :) + wt(q) .* (reshape(us, [], 3, 1) .* reshape(us, [], 1, 3));
        TT(q, :) = TT(q, :) + wt(q) .* reshape(F.Tn(li(q), s, :), [], 27);
    end
    wsum = wsum + wt;
end, end, end
if strcmp(what, 'conf'), out = wsum; return; end         % the read weight: the confidence of field: lift
bad = wsum < F.read_min;
Mn = M ./ max(wsum, eps); e1 = T;                       % principal axis: power iteration from the heading,
for it = 1:25                                           % which lies within 30 deg of it at every corner read
    e1 = reshape(sum(Mn .* reshape(e1, n, 1, 3), 3), n, 3); e1 = e1 ./ max(vecnorm(e1, 2, 2), eps);
end
e1 = e1 .* sign(sum(e1 .* T, 2) + 1e-300);
[~, ax] = min(abs(e1), [], 2); ref = zeros(n, 3); ref(sub2ind([n 3], (1:n)', ax)) = 1;   % as frames_at of nim_tractography_mmf
e2 = ref - sum(ref .* e1, 2) .* e1; e2 = e2 ./ max(vecnorm(e2, 2, 2), eps); e3 = cross(e1, e2, 2);
if strcmp(what, 'mf')
    out = cat(2, reshape(e1, n, 1, 3), reshape(e2, n, 1, 3), reshape(e3, n, 1, 3)); out(bad, :, :) = NaN; return
end
TT = reshape(TT ./ max(wsum, eps), n, 3, 3, 3);
out = zeros(n, 3, 3, 3);
for k = 1:3
    de1 = zeros(n, 3);
    for b = 1:3, for a = 1:3, de1(:, b) = de1(:, b) + e1(:, a) .* TT(:, a, b, k); end, end
    w12 = sum(e2 .* de1, 2); w13 = sum(e3 .* de1, 2);
    out(:, 1, 2, k) = w12; out(:, 2, 1, k) = -w12; out(:, 1, 3, k) = w13; out(:, 3, 1, k) = -w13;
end
out(bad, :, :, :) = NaN;
end
