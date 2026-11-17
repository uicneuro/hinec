function [nim, n_refit, n_zero] = nim_dt_refit_zeroed(nim, floor_ev)
% nim_dt_refit_zeroed: a tensor for the voxels of the mask that the tensor fit left at zero.
%
%   [nim, n_refit, n_zero] = nim_dt_refit_zeroed(nim)
%   [nim, n_refit, n_zero] = nim_dt_refit_zeroed(nim, floor_ev)
%
% nim_dt_spd leaves a voxel of the mask at zero (DT and eval zero, FA 0 after nim_fa;
% nim_eig makes its evec the identity) when it cannot keep a tensor there, and every
% tracker then reads the voxel as no tissue. Nearly all of these holes are the voxels
% whose least-squares tensor has a negative eigenvalue: nim_dt_spd sends them to the
% BFGS fit of lib/bfgs/vox_dt_bfgs, and the fit's bound on the eigenvalues
% (0.01 mm^2/s) rejects every result. vox_dt_bfgs scales b by 1/1000 and then
% multiplies D by 1000 where it should divide (a factor 1e6), and its gradient (gradf)
% reads the x-components of the first three gradient directions instead of each
% sample's direction, so it does not minimise its objective either: on the TrackerYoga
% holes its largest eigenvalue comes out 6e6 to 1e7 times the least-squares one, and
% fixing the scale alone would not repair it. The holes were 352 of the 2974 mask
% voxels of the TrackerYoga phantom (b = 1000), 1698 of its noisy version and 8562 of
% the 182857 of the ISMRM 2015 phantom.
%
% Each such voxel gets, where its signals allow it, its log-linear least-squares
% tensor with the eigenvalues raised to floor_ev (default 1e-7 mm^2/s), the nearest
% positive definite tensor with the same axes; DT, evec, eval and, when the nim has
% it, FA (nim_fa's formula) are set from it. A voxel is a hole when it is in the mask
% (mask > 0) and its eigenvalues, and its DT where the nim has one, are all zero (a
% voxel whose DT is not zero keeps it, even where nim_eig skipped it for a mask value
% other than 1). It is refitted when its b0 is positive and left at zero when fewer
% than six of its diffusion-weighted samples are positive or their gradients do not
% determine the tensor (the voxels nim_dt_spd counts as too few usable samples), when
% the least-squares tensor has no eigenvalue above floor_ev, or when its largest
% eigenvalue exceeds the fit's bound of 0.01 mm^2/s. Voxels with a tensor and voxels
% outside the mask are not touched. The units are those of nim_dt_spd: b in s/mm^2,
% D in mm^2/s. Where the least-squares tensor has one or two negative eigenvalues, as
% at the noisy edge of a mask, the floored tensor is close to linear (FA near 1).
%
% A nim without eigenvalues or mask, or without holes, is returned as it is, silently.
% A nim with holes but without the raw signals (img_b0, img_bi, bval, bvec,
% thrsh_b0), or with signals that do not match its grid, is returned as it is, with a
% message. n_refit counts the voxels refitted, n_zero the holes found.
%
% main calls it after the fit, and again before the save for the voxels that the
% improved mask adds (outside the mask of the fit, so never fitted), behind
% preprocessing.refit_zeroed_tensors of the pipeline's config; runTractography calls
% it on the dataset it loads, so that a nim made before reaches the trackers without
% holes (in memory: the file is not rewritten), behind the same key of the tracking
% config (a call without a config always refits). The moving-frames test suite did the
% same on the side for the TrackerYoga data (runtrackeryogatest 'Fit' 'refit').

if nargin < 2 || isempty(floor_ev)
    floor_ev = 1e-7;
end
n_refit = 0;
n_zero = 0;
if ~isfield(nim, 'eval') || ~isfield(nim, 'mask') || isempty(nim.eval) || isempty(nim.mask) ...
        || ndims(nim.eval) ~= 4 || size(nim.eval, 4) ~= 3
    return
end
dims = size(nim.eval, 1:3);
if ~isequal(size(nim.mask, 1:3), dims)
    return
end
hole = nim.mask > 0 & all(nim.eval == 0, 4);
if isfield(nim, 'DT') && ~isempty(nim.DT) && isequal(size(nim.DT), [dims 6])
    hole = hole & all(nim.DT == 0, 4);
end
n_zero = nnz(hole);
if n_zero == 0
    return
end
need = {'img_b0', 'img_bi', 'bval', 'bvec', 'thrsh_b0', 'DT', 'evec'};
missing = need(~isfield(nim, need));
if ~isempty(missing)
    fprintf('Tensor refit: %d voxels of the mask have no tensor, but the nim has no %s: left at zero.\n', ...
        n_zero, strjoin(missing, ', '));
    return
end

% The diffusion-weighted volumes, paired with their b-values and gradients as nim_dt_spd pairs them.
sel = nim.bval >= nim.thrsh_b0;
if ~isequal(size(nim.img_b0, 1:3), dims) || ~isequal(size(nim.img_bi, 1:3), dims) || size(nim.img_bi, 4) ~= nnz(sel)
    fprintf('Tensor refit: %d voxels of the mask have no tensor, but the signals of the nim do not match its grid or its b-values: left at zero.\n', n_zero);
    return
end
b = reshape(double(nim.bval(sel)), [], 1);
g = double(nim.bvec(sel, :));
H = [g(:, 1).^2, g(:, 2).^2, g(:, 3).^2, 2 * g(:, 1) .* g(:, 2), 2 * g(:, 2) .* g(:, 3), 2 * g(:, 3) .* g(:, 1)];
has_fa = isfield(nim, 'FA') && ~isempty(nim.FA) && isequal(size(nim.FA, 1:3), dims);
idx = find(hole & nim.img_b0 > 0);
for n = 1:numel(idx)
    [x, y, z] = ind2sub(dims, idx(n));
    s  = reshape(double(nim.img_bi(x, y, z, :)), [], 1);
    s0 = double(nim.img_b0(x, y, z));
    use = s > 0;
    if nnz(use) < 6 || rank(H(use, :)) < 6
        continue
    end
    d = H(use, :) \ (log(s0 ./ s(use)) ./ b(use));    % [Dxx Dyy Dzz Dxy Dyz Dxz], as nim_dt_spd
    [Q, L] = eig(nim_reshape_d(d));
    [l, order] = sort(diag(L).', 'descend');
    Q = Q(:, order);
    if any(~isfinite(l)) || l(1) <= floor_ev || l(1) > 0.01
        continue
    end
    l = max(l, floor_ev);
    Dn = Q * diag(l) * Q.';
    nim.DT(x, y, z, :) = [Dn(1, 1) Dn(2, 2) Dn(3, 3) Dn(1, 2) Dn(2, 3) Dn(1, 3)];
    nim.evec(x, y, z, :, :) = Q;
    nim.eval(x, y, z, :) = l;
    if has_fa
        nim.FA(x, y, z) = sqrt(1/2) * sqrt((l(1) - l(2))^2 + (l(2) - l(3))^2 + (l(3) - l(1))^2) / sqrt(sum(l .^ 2));
    end
    n_refit = n_refit + 1;
end
fprintf('Tensor refit: %d of the %d voxels of the mask without a tensor refitted (least squares, eigenvalues floored at %.3g mm^2/s).\n', ...
    n_refit, n_zero, floor_ev);
end
