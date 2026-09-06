function [offsets, info] = nim_seed_offsets(density)
% nim_seed_offsets: Sub-voxel seed offsets, honouring the requested count exactly.
%
%   offsets = nim_seed_offsets(8)   -> 8x3, offsets in [-0.5, +0.5]
%
% Returns EXACTLY `density` offsets. The previous implementation computed
%   per_axis = ceil(density^(1/3));  offsets = per_axis^3 points
% which silently rounded UP to the next perfect cube: asking for 2..7 gave 8,
% and 9..26 gave 27. `seed_density: 4` therefore placed 8 seeds per voxel, twice
% what the config said.
%
% Placement:
%   - density is a perfect cube -> the full per_axis^3 lattice, in ndgrid order.
%     Bit-identical to the old behaviour for 1, 8, 27, 64 (which is every value
%     the shipped configs actually use), so existing runs are unchanged.
%   - otherwise -> a deterministic farthest-point subset of the next-larger
%     lattice, so the seeds stay spread through the voxel instead of clustering
%     in one corner. No RNG, so it stays reproducible - which is what makes
%     per-streamline comparison across a convergence ladder possible.
%
% info reports the actual count, the lattice used, and the resulting spacing.

    density = max(1, round(double(density)));

    per_axis = max(1, ceil(density^(1/3) - 1e-12));
    lat = build_lattice(per_axis);

    if size(lat, 1) == density
        offsets = lat;
        exact_lattice = true;
    else
        % Need fewer points than the enclosing lattice provides. Take a spread
        % subset rather than the first N (which would clump at one face).
        offsets = farthest_point_subset(lat, density);
        exact_lattice = false;
    end

    info = struct('requested', density, ...
                  'actual', size(offsets, 1), ...
                  'per_axis', per_axis, ...
                  'lattice_points', size(lat, 1), ...
                  'exact_lattice', exact_lattice, ...
                  'spacing', 1 / per_axis);
end

% =========================================================================
function lat = build_lattice(per_axis)
    edges = linspace(-0.5, 0.5, per_axis + 1);
    c = (edges(1:end-1) + edges(2:end)) / 2;
    [ox, oy, oz] = ndgrid(c, c, c);
    lat = [ox(:), oy(:), oz(:)];
end

function sub = farthest_point_subset(lat, k)
% Deterministic greedy farthest-point sampling. Starts from the lattice point
% closest to the voxel centre (ties broken by index, so it is reproducible) and
% repeatedly adds the point maximising the minimum distance to those chosen.
    n = size(lat, 1);
    k = min(k, n);

    d0 = sum(lat.^2, 2);
    [~, first] = min(d0);            % min() returns the lowest index on ties
    chosen = zeros(1, k);
    chosen(1) = first;

    mind = sum((lat - lat(first, :)).^2, 2);
    mind(first) = -inf;

    for i = 2:k
        [~, nxt] = max(mind);        % max() returns the lowest index on ties
        chosen(i) = nxt;
        mind = min(mind, sum((lat - lat(nxt, :)).^2, 2));
        mind(chosen(1:i)) = -inf;
    end

    sub = lat(sort(chosen), :);
end
