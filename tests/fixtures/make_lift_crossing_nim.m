function [nim, truth] = make_lift_crossing_nim(varargin)
%MAKE_LIFT_CROSSING_NIM  A nim carrying the raw diffusion signal of two crossing analytic bundles.
%
%   [nim, truth] = make_lift_crossing_nim('alpha', 60, 'fmax', 0.5, ...)
%
%   Bundle A is the arc field of the MMF notes (Sec. 9): circles about the axis through c, direction
%   phi-hat. Bundle B is the logarithmic spiral of pitch alpha about the same axis, crossing every circle
%   at alpha: u_B = -sin(alpha) rho-hat + cos(alpha) phi-hat. B is present in the sector |phi - 90| < 25
%   deg with f_B = fmax cos^2(90 (phi - 90)/25 deg), f_A = 1 - f_B. The signal is the two-stick model
%       S(g)/S0 = f_A exp(-b (l_perp + dl (g.u_A)^2)) + f_B exp(-b (l_perp + dl (g.u_B)^2)),
%   with Rician noise of standard deviation 1/snr (S0 = 1) on every volume, the b = 0 volume included.
%   The DTI fields (evec, eval, FA) are the log-linear tensor fit of the noisy signal.
%
%   Options: 'alpha' (60), 'fmax' (0.5), 'dims' ([80 56 11]), 'center' ([40.5 8]), 'lambda'
%   ([1.7e-3 3e-4]), 'b' (1000), 'bvec' ([] = 32 Fibonacci directions), 'snr' (20, Inf = none),
%   'seed' (0, the noise stream), 'fiso' (0, a free-water fraction with D = 3e-3 mm^2/s in every voxel:
%   S/S0 = (1 - fiso) (the two sticks) + fiso exp(-3e-3 b)), 'rho' ([8 44], the radii of the arcs kept in the
%   mask: the curvature of bundle A is 1 / rho per voxel).
%   truth: .c, .uA, .uB [X Y Z 3], .fB [X Y Z], .alpha, .fmax.

p = inputParser;
p.addParameter('alpha', 60); p.addParameter('fmax', 0.5); p.addParameter('dims', [80 56 11]);
p.addParameter('center', [40.5 8]); p.addParameter('lambda', [1.7e-3 3e-4]); p.addParameter('b', 1000);
p.addParameter('bvec', []); p.addParameter('snr', 20); p.addParameter('seed', 0); p.addParameter('fiso', 0);
p.addParameter('rho', [8 44]);
p.parse(varargin{:}); o = p.Results;
dims = o.dims; c = o.center; lp = o.lambda(2); dl = o.lambda(1) - o.lambda(2);
g = o.bvec;
if isempty(g)
    k = (0:31)' + 0.5; z = 1 - k / 32; r = sqrt(1 - z .^ 2); t = pi * (1 + sqrt(5)) * k;
    g = [r .* cos(t), r .* sin(t), z];                         % 32 directions over the sphere
end
g = g ./ vecnorm(g, 2, 2); N = size(g, 1);
[X, Y, ~] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
rx = X - c(1); ry = Y - c(2); rho = hypot(rx, ry); phi = atan2d(ry, rx);
rh = cat(4, rx ./ rho, ry ./ rho, zeros(dims)); ph = cat(4, -ry ./ rho, rx ./ rho, zeros(dims));
uA = ph; uB = -sind(o.alpha) * rh + cosd(o.alpha) * ph;
mask = rho >= o.rho(1) & rho <= o.rho(2) & ry > 1;
fB = o.fmax * cosd(90 * min(abs(phi - 90) / 25, 1)) .^ 2; fB(abs(phi - 90) >= 25 | ~mask) = 0; fA = 1 - fB;
cA = reshape(reshape(uA, [], 3) * g', [dims N]); cB = reshape(reshape(uB, [], 3) * g', [dims N]);
S = (1 - o.fiso) * (fA .* exp(-o.b * (lp + dl * cA .^ 2)) + fB .* exp(-o.b * (lp + dl * cB .^ 2))) + o.fiso * exp(-o.b * 3e-3);
S0 = ones(dims);
if isfinite(o.snr)
    rs = RandStream('mt19937ar', 'Seed', o.seed); sg = 1 / o.snr;
    S  = abs(S + sg * randn(rs, size(S)) + 1i * sg * randn(rs, size(S)));
    S0 = abs(S0 + sg * randn(rs, dims) + 1i * sg * randn(rs, dims));
end
nim = struct('xdim', dims(1), 'ydim', dims(2), 'zdim', dims(3), 'mask', double(mask));
nim.bval = [0; o.b * ones(N, 1)]; nim.bvec = [0 0 0; g]; nim.thrsh_b0 = 50;
nim.img_b0 = S0 .* mask; nim.img_bi = S .* mask;
% the log-linear DTI fit of the noisy signal
Bm = -o.b * [g(:, 1) .^ 2, g(:, 2) .^ 2, g(:, 3) .^ 2, 2 * g(:, 1) .* g(:, 2), 2 * g(:, 1) .* g(:, 3), 2 * g(:, 2) .* g(:, 3)];
yv = log(max(reshape(S, [], N) ./ max(S0(:), eps), 1e-6));
dv = (Bm \ yv')';
nim.evec = zeros([dims 3 3]); nim.eval = zeros([dims 3]); nim.FA = zeros(dims);
for i = find(mask)'
    Dm = [dv(i, 1) dv(i, 4) dv(i, 5); dv(i, 4) dv(i, 2) dv(i, 6); dv(i, 5) dv(i, 6) dv(i, 3)];
    [V, Lm] = eig((Dm + Dm') / 2); [l, ord] = sort(diag(Lm), 'descend'); V = V(:, ord);
    [ix, iy, iz] = ind2sub(dims, i);
    nim.evec(ix, iy, iz, :, :) = reshape(V, 1, 1, 1, 3, 3); nim.eval(ix, iy, iz, :) = l;
    md = mean(l); nim.FA(ix, iy, iz) = sqrt(1.5 * sum((l - md) .^ 2) / max(sum(l .^ 2), eps));
end
truth = struct('c', c, 'uA', uA, 'uB', uB, 'fB', fB, 'alpha', o.alpha, 'fmax', o.fmax);
end
