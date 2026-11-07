function [keep, coh] = nim_seed_coherence_mask(nim, coh_min)
%NIM_SEED_COHERENCE_MASK Voxels whose principal direction agrees with its neighbours.
%   [keep, coh] = nim_seed_coherence_mask(nim, coh_min) returns keep = coh > coh_min,
%   where the coherence of voxel x is the largest eigenvalue of the mean dyadic
%
%       M(x) = (1/|N(x)|) sum_{y in N(x)} v1(y) v1(y)',
%
%   N(x) the voxels of the 3 x 3 x 3 block around x (x included) that carry a
%   tensor (lambda_1 > 0, inside nim.mask when there is one) and v1 the unit
%   principal eigenvector (nim.evec(:,:,:,:,1)). The dyadic makes it sign
%   invariant (v1 is a line field). coh lies in [1/3, 1]: 1 when every
%   neighbour points the same way, near 1/3 when they point every way. Unlike
%   CL (nim_seed_cl_mask), which reads one voxel's tensor shape, this asks
%   whether the axis is consistent with the voxels around it: a crossing whose
%   tensor looks linear but points along neither bundle has a high CL and a low
%   coherence. A voxel without a tensor has coh = 0 and is dropped.
%
%   See also nim_seed_cl_mask, runTractography (seeding.coherence_min), nim_config_schema.

if ~isfield(nim, 'evec') || isempty(nim.evec) || ~isfield(nim, 'eval') || isempty(nim.eval)
    error('nim_seed_coherence_mask:noEig', 'seeding.coherence_min needs nim.evec and nim.eval.');
end
valid = nim.eval(:, :, :, 1) > 0;
if isfield(nim, 'mask') && ~isempty(nim.mask)
    valid = valid & logical(nim.mask);
end
v = double(nim.evec(:, :, :, :, 1));
v = v ./ max(sqrt(sum(v .^ 2, 4)), eps);
v(repmat(~valid, [1 1 1 3])) = 0;

k   = ones(3, 3, 3);
cnt = convn(double(valid), k, 'same');
pr  = [1 1; 2 2; 3 3; 1 2; 1 3; 2 3];
m   = cell(1, 6);
for j = 1:6
    m{j} = convn(v(:, :, :, pr(j, 1)) .* v(:, :, :, pr(j, 2)), k, 'same') ./ max(cnt, 1);
end
coh = zeros(size(cnt));
coh(valid) = largestEig(m{1}(valid), m{2}(valid), m{3}(valid), m{4}(valid), m{5}(valid), m{6}(valid));
keep = coh > coh_min;
end

function l = largestEig(a, b, c, d, e, f)
% Largest eigenvalue of the symmetric [a d e; d b f; e f c], elementwise
% (the trigonometric closed form for a symmetric 3 x 3 matrix).
q  = (a + b + c) / 3;
p1 = d .^ 2 + e .^ 2 + f .^ 2;
p2 = (a - q) .^ 2 + (b - q) .^ 2 + (c - q) .^ 2 + 2 * p1;
p  = sqrt(p2 / 6);
l  = q;                                  % isotropic: all three equal
nz = p > 1e-12;
ip = 1 ./ p(nz);
B11 = (a(nz) - q(nz)) .* ip;  B22 = (b(nz) - q(nz)) .* ip;  B33 = (c(nz) - q(nz)) .* ip;
B12 = d(nz) .* ip;  B13 = e(nz) .* ip;  B23 = f(nz) .* ip;
r = (B11 .* (B22 .* B33 - B23 .^ 2) - B12 .* (B12 .* B33 - B23 .* B13) + B13 .* (B12 .* B23 - B22 .* B13)) / 2;
phi = acos(min(max(r, -1), 1)) / 3;
l(nz) = q(nz) + 2 * p(nz) .* cos(phi);
end
