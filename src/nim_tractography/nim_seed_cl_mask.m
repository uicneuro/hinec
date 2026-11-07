function [keep, cl] = nim_seed_cl_mask(nim, cl_min)
%NIM_SEED_CL_MASK Voxels whose tensor is linear enough to seed from.
%   [keep, cl] = nim_seed_cl_mask(nim, cl_min) returns keep = CL > cl_min, where
%   CL = (lambda_1 - lambda_2) / lambda_1 is the linear shape of the voxel tensor
%   (nim.eval, eigenvalues in descending order) - the same confidence the hold
%   scheme of the mmf tracker reads (mmf.hold_min), here per voxel, without its
%   interpolation or erosion. A voxel with lambda_1 <= 0 has CL = 0 and is
%   dropped. A crossing fitted as one tensor is a pancake (lambda_1 ~ lambda_2,
%   CL ~ 0), so cl_min > 0 keeps seeds out of crossings, where a seed has no
%   trusted direction to start from. cl_min = 0 keeps every voxel with
%   lambda_1 > lambda_2 (tracking is unchanged when the option is off: the
%   caller skips this function).
%
%   See also runTractography (seeding.cl_min), nim_config_schema.

if ~isfield(nim, 'eval') || isempty(nim.eval)
    error('nim_seed_cl_mask:noEval', 'seeding.cl_min needs nim.eval (the tensor eigenvalues).');
end
l1 = nim.eval(:, :, :, 1);
l2 = nim.eval(:, :, :, 2);
cl = zeros(size(l1));
ok = l1 > 0;
cl(ok) = (l1(ok) - l2(ok)) ./ l1(ok);
keep = cl > cl_min;
end
