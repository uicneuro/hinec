function frames = nim_build_frames(nim, options)
% nim_build_frames: Build a smooth, sign-consistent orthonormal frame FIELD.
%
% This is the Phase-2 precompute for MMF (moving-frames) tractography. It turns
% the per-voxel principal eigenvector e1 into a globally sign-consistent, gently
% smoothed unit vector field, then constructs e2, e3 from e1 alone using the
% reference-axis projection (robust to the lambda2 ~ lambda3 degeneracy that
% makes the DTI e2/e3 arbitrary). See MMF_Tractography_ImplementationPlan §6.2.
%
% Arguments:
%   nim     - NIM structure with nim.evec [D1 D2 D3 3 3] (evec(:,:,:,c,k) =
%             component c of the k-th eigenvector) and, ideally, nim.mask/nim.FA.
%   options - struct; uses:
%               .frame_smooth_sigma  Gaussian sigma in voxels (default 1.0; 0 = off)
%
% Returns:
%   frames  - [D1 D2 D3 3 3], frames(:,:,:,c,i) = component c of frame vector e_i.
%             e1 = sign-consistent, smoothed principal eigenvector.
%             e2, e3 = reference-axis projection of e1 (orthonormal).
%             Frames outside the brain mask are zeroed.

if nargin < 2 || isempty(options)
    options = struct();
end
if ~isfield(options, 'frame_smooth_sigma') || isempty(options.frame_smooth_sigma)
    options.frame_smooth_sigma = 1.0;
end

if ~isfield(nim, 'evec')
    error('nim_build_frames: nim.evec not found. Run nim_eig() first.');
end

dims = size(nim.FA);

% --- Brain mask (propagation/analysis domain) --------------------------------
if isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
    mask = nim.mask > 0.5;
elseif isfield(nim, 'FA')
    mask = nim.FA > 0;
else
    mask = true(dims);
end

% --- Principal eigenvector components ----------------------------------------
ex = squeeze(nim.evec(:,:,:,1,1));
ey = squeeze(nim.evec(:,:,:,2,1));
ez = squeeze(nim.evec(:,:,:,3,1));

% Normalize (guard against any un-normalized / zero voxels)
nrm = sqrt(ex.^2 + ey.^2 + ez.^2);
valid = nrm > 1e-6 & mask;
ex(~valid) = 0; ey(~valid) = 0; ez(~valid) = 0;
nrm(nrm < 1e-6) = 1;
ex = ex ./ nrm; ey = ey ./ nrm; ez = ez ./ nrm;

% --- (a) Global sign consistency (BFS flood over the mask) --------------------
% Eigenvectors are sign-ambiguous; a derivative of e1 (the connection form) is
% meaningless unless neighbouring vectors point the same way. Flood-fill from the
% highest-FA voxel, flipping any neighbour whose dot with the visited voxel < 0.
% A globally consistent field may not exist (topological obstruction) -- residual
% seams are tolerated (the tracer also flips signs locally).
[ex, ey, ez] = flood_sign_consistency(ex, ey, ez, valid, nim, dims);

% --- (b) Optional Gaussian smoothing of e1, then renormalize -----------------
sigma = options.frame_smooth_sigma;
if sigma > 0
    % Zero outside the mask so smoothing does not bleed garbage in; smooth the
    % mask too and divide to approximate a normalized (mask-aware) convolution.
    m = double(valid);
    sm = imgaussfilt3(m, sigma);
    sm(sm < 1e-6) = 1;
    ex = imgaussfilt3(ex .* m, sigma) ./ sm;
    ey = imgaussfilt3(ey .* m, sigma) ./ sm;
    ez = imgaussfilt3(ez .* m, sigma) ./ sm;

    nrm = sqrt(ex.^2 + ey.^2 + ez.^2);
    bad = nrm < 1e-6;
    nrm(bad) = 1;
    ex = ex ./ nrm; ey = ey ./ nrm; ez = ez ./ nrm;
    ex(bad) = 0; ey(bad) = 0; ez(bad) = 0;
end

% --- (c) Build e2, e3 by reference-axis projection (vectorized §2.2) ----------
% r = Cartesian axis least aligned with e1 (smallest |component|).
absx = abs(ex); absy = abs(ey); absz = abs(ez);
pick_x = (absx <= absy) & (absx <= absz);
pick_y = (~pick_x) & (absy <= absz);
pick_z = ~(pick_x | pick_y);

% r as a unit axis vector per voxel; r.e1 = the e1 component along that axis.
rx = double(pick_x); ry = double(pick_y); rz = double(pick_z);
r_dot_e1 = rx .* ex + ry .* ey + rz .* ez;

e2x = rx - r_dot_e1 .* ex;
e2y = ry - r_dot_e1 .* ey;
e2z = rz - r_dot_e1 .* ez;
n2 = sqrt(e2x.^2 + e2y.^2 + e2z.^2);
n2(n2 < 1e-9) = 1;
e2x = e2x ./ n2; e2y = e2y ./ n2; e2z = e2z ./ n2;

% e3 = e1 x e2
e3x = ey .* e2z - ez .* e2y;
e3y = ez .* e2x - ex .* e2z;
e3z = ex .* e2y - ey .* e2x;

% --- (d) Zero frames outside the mask, assemble output -----------------------
z = ~valid;
ex(z)=0; ey(z)=0; ez(z)=0;
e2x(z)=0; e2y(z)=0; e2z(z)=0;
e3x(z)=0; e3y(z)=0; e3z(z)=0;

frames = zeros([dims, 3, 3]);
frames(:,:,:,1,1) = ex;  frames(:,:,:,2,1) = ey;  frames(:,:,:,3,1) = ez;
frames(:,:,:,1,2) = e2x; frames(:,:,:,2,2) = e2y; frames(:,:,:,3,2) = e2z;
frames(:,:,:,1,3) = e3x; frames(:,:,:,2,3) = e3y; frames(:,:,:,3,3) = e3z;

end


function [ex, ey, ez] = flood_sign_consistency(ex, ey, ez, valid, nim, dims)
% BFS flood-fill sign alignment over the valid mask (6-connectivity).

% Seed at the highest-FA valid voxel (fallback: first valid voxel).
if isfield(nim, 'FA') && ~isempty(nim.FA)
    fa = nim.FA;
    fa(~valid) = -Inf;
    [~, seed_lin] = max(fa(:));
else
    seed_lin = find(valid, 1, 'first');
end
if isempty(seed_lin) || ~valid(seed_lin)
    return;  % nothing to do
end

visited = false(dims);
% Manual stack of linear indices (preallocated to number of valid voxels).
nvalid = nnz(valid);
stack = zeros(nvalid, 1);
sp = 1;
stack(sp) = seed_lin;
visited(seed_lin) = true;

% Neighbour linear-index offsets for 6-connectivity.
D1 = dims(1); D2 = dims(2);
offs = [1, -1, D1, -D1, D1*D2, -D1*D2];

while sp > 0
    cur = stack(sp);
    sp = sp - 1;
    exc = ex(cur); eyc = ey(cur); ezc = ez(cur);

    [ci, cj, ck] = ind2sub(dims, cur);
    for o = 1:6
        nb = cur + offs(o);
        if nb < 1 || nb > numel(visited), continue; end
        % Reject wrap-around across the dim-1 face for the +/-1 offsets.
        [ni, nj, nk] = ind2sub(dims, nb);
        if abs(ni-ci) + abs(nj-cj) + abs(nk-ck) ~= 1, continue; end
        if visited(nb) || ~valid(nb), continue; end

        if exc*ex(nb) + eyc*ey(nb) + ezc*ez(nb) < 0
            ex(nb) = -ex(nb); ey(nb) = -ey(nb); ez(nb) = -ez(nb);
        end
        visited(nb) = true;
        sp = sp + 1;
        stack(sp) = nb;
    end
end
end
