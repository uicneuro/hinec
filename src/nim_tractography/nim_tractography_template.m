function [tracks, meta] = nim_tractography_template(nim, options)
% nim_tractography_template: The smallest tracker that satisfies the contract.
%
%   [tracks, meta] = nim_tractography_template(nim, options)
%
% This is a WORKED EXAMPLE of the interface in docs/TRACKER_INTERFACE.md, kept
% deliberately short so the whole hand-off fits on one screen. To write a new
% algorithm, copy this file, rename it nim_tractography_<name>.m, replace the
% two marked sections (the direction rule and the stepping rule), and add
% <name> to the algorithm enum in nim_config_schema and the dispatch in
% runTractography step 5. Select it with  algorithm: <name>  in the YAML.
%
% It is runnable as shipped:  algorithm: template  in any config.
%
% What it does: fixed-step Euler along the DTI principal eigenvector. The
% eigenvector is interpolated as the dyadic v1*v1' (sign-invariant, because the
% eigensolver's sign is arbitrary from voxel to voxel) and the principal axis is
% recovered afterwards. Stops on FA, the brain mask, the volume edge, the turn
% limit and the arc limit. Nothing else - no ACT, no CSD, no adaptive stepping.
%
% INPUT - everything is documented, per field, in the run log ("Tracker input")
% and in docs/TRACKER_INTERFACE.md. The fields this tracker reads:
%   nim.FA         [X Y Z]      termination
%   nim.mask       [X Y Z]      termination (optional)
%   nim.evec       [X Y Z 3 3]  evec(x,y,z,:,1) is the principal eigenvector
%   options.seed_mask       [X Y Z] logical - seed inside these voxels
%   options.seed_density    seeds per voxel
%   options.step_size       h, voxels
%   options.termination_fa  stop below this FA
%   options.angle_thresh    max turn, degrees per voxel of arc
%   options.max_arc         stop a half-track after this arc, voxels
%   options.min_length      discard a track shorter than this, voxels
% Coordinates are voxel indices, 1-based, continuous; voxel (i,j,k) is centred
% at [i j k]. A step of h moves h voxels.
%
% OUTPUT
%   tracks  cell {T x 1}; tracks{t} is N x 3, ordered
%           [backward half, reversed; seed; forward half]
%   meta    .seed_index (1 x T), .seed_points (T x 3), .n_seeds

dims = size(nim.FA);
h    = options.step_size;

% ---- seeds: from the mask runTractography built, nothing else ---------------
[offsets, ~] = nim_seed_offsets(options.seed_density);
[ix, iy, iz] = ind2sub(dims, find(options.seed_mask));
centres = [ix iy iz];
seed_points = zeros(size(centres, 1) * size(offsets, 1), 3);
for k = 1:size(offsets, 1)
    seed_points((k-1)*size(centres,1) + (1:size(centres,1)), :) = centres + offsets(k, :);
end
n_seeds = size(seed_points, 1);
fprintf('template: %d seed voxels x %d = %d seeds, h = %.3g\n', ...
    size(centres, 1), size(offsets, 1), n_seeds, h);

% ---- direction rule (REPLACE THIS for a new algorithm) ----------------------
% Six interpolants for the unique dyadic components; 'none' extrapolation makes
% any point outside the grid read as NaN, which is the leaving-the-volume test.
v1 = reshape(nim.evec(:, :, :, :, 1), [dims 3]);
gv = {1:dims(1), 1:dims(2), 1:dims(3)};
comp = @(a, b) griddedInterpolant(gv, v1(:,:,:,a) .* v1(:,:,:,b), 'linear', 'none');
D  = {comp(1,1), comp(2,2), comp(3,3), comp(1,2), comp(1,3), comp(2,3)};
FA = griddedInterpolant(gv, nim.FA, 'linear', 'none');
if isfield(nim, 'mask') && ~isempty(nim.mask)
    inmask = @(p) nim.mask(round(p(1)), round(p(2)), round(p(3))) > 0.5;
else
    inmask = @(p) true;
end

    function v = direction_at(p, ref)
    % Principal axis of the interpolated dyadic, signed to follow ref.
        c = cellfun(@(f) f(p(1), p(2), p(3)), D);
        if any(~isfinite(c)), v = []; return; end
        v = nim_principal_dir(c(1), c(2), c(3), c(4), c(5), c(6));
        if isempty(v), return; end
        v = v(:)';
        if ~isempty(ref) && dot(v, ref) < 0, v = -v; end
    end

    function ok = alive(p)
        ok = all(p >= 1) && all(p <= dims) && inmask(p) && FA(p(1), p(2), p(3)) >= options.termination_fa;
    end

% ---- stepping rule (REPLACE THIS for a new algorithm) -----------------------
    function pts = half_track(seed, d0)
    % Euler: p <- p + h * v(p). Returns the points AFTER the seed, in order.
        pts = zeros(ceil(options.max_arc / h) + 1, 3);
        n = 0; p = seed; d = d0; arc = 0;
        while arc < options.max_arc
            v = direction_at(p, d);
            if isempty(v), break; end
            turn = acosd(min(1, dot(v, d)));           % degrees, this step
            if turn / h > options.angle_thresh, break; end   % per voxel of arc
            p = p + h * v;
            if ~alive(p), break; end
            n = n + 1; pts(n, :) = p; d = v; arc = arc + h;
        end
        pts = pts(1:n, :);
    end

% ---- one track per seed, both directions, then the length rule -------------
tracks = cell(n_seeds, 1);
kept   = false(n_seeds, 1);
for s = 1:n_seeds
    seed = seed_points(s, :);
    if ~alive(seed), continue; end
    d0 = direction_at(seed, []);
    if isempty(d0), continue; end
    fwd = half_track(seed,  d0);
    bwd = half_track(seed, -d0);
    track = [flipud(bwd); seed; fwd];
    if size(track, 1) >= 2 && sum(sqrt(sum(diff(track, 1, 1).^2, 2))) >= options.min_length
        tracks{s} = track;
        kept(s) = true;
    end
end

tracks = tracks(kept);
meta = struct('seed_index', find(kept)', 'seed_points', seed_points(kept, :), 'n_seeds', n_seeds);
fprintf('template: %d of %d seeds gave a track >= %.3g voxels\n', numel(tracks), n_seeds, options.min_length);
end
