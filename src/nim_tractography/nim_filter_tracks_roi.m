function [tracks, stats] = nim_filter_tracks_roi(tracks, nim, options)
% nim_filter_tracks_roi: Keep/discard streamlines by the regions they touch.
%
%   [tracks, stats] = nim_filter_tracks_roi(tracks, nim, options)
%
% Reads from options:
%   include_roi      - keep only tracks touching these regions (waypoint)
%   roi_filter_mode  - 'all' (must touch every include region) | 'any'
%   exclude_roi      - discard tracks touching any of these regions
%   roi_filter_dilate- dilate the filter masks by this many voxels (default 0)
%
% With no include/exclude set, tracks are returned untouched.
%
% Combined with ROI seeding this gives both requested behaviours:
%   seed_roi: X                       -> tracks seeded in X (bidirectional, so
%                                        this already recovers tracks THROUGH X)
%   seed_roi: <brain>, include_roi: X -> full-brain tracking, then keep only the
%                                        tracks that pass through X

    stats = struct('n_in', numel(tracks), 'n_out', numel(tracks), ...
                   'n_dropped_include', 0, 'n_dropped_exclude', 0, 'applied', false, ...
                   'keep', true(1, numel(tracks)));

    inc = get_opt(options, 'include_roi', {});
    exc = get_opt(options, 'exclude_roi', {});
    eps = get_opt(options, 'endpoints_in', {});
    con = get_opt(options, 'contained_in', {});
    anyi = get_opt(options, 'any_in', {});

    % Per-dimension and total length criteria, in MILLIMETRES - the scorer
    % converts to RAS mm (sft.to_rasmm()) before measuring, so voxel units would
    % be wrong by the voxel size (2x on this data).
    Lspec = struct( ...
        'total',  {num_opt(options, 'length')}, ...
        'net',    {{num_opt(options,'length_x'),     num_opt(options,'length_y'),     num_opt(options,'length_z')}}, ...
        'abs',    {{num_opt(options,'length_x_abs'), num_opt(options,'length_y_abs'), num_opt(options,'length_z_abs')}});
    have_len = ~isempty(Lspec.total) || ...
               any(~cellfun(@isempty, Lspec.net)) || any(~cellfun(@isempty, Lspec.abs));

    if isempty(inc) && isempty(exc) && isempty(eps) && isempty(con) && ...
       isempty(anyi) && ~have_len
        return;
    end

    mode = lower(char(string(get_opt(options, 'roi_filter_mode', 'all'))));
    dil  = get_opt(options, 'roi_filter_dilate', 0);
    dims = size(nim.FA);

    fprintf('\n=== ROI track filtering ===\n');

    % Build one mask per include region so 'all' can be tested per region.
    inc_masks = {};
    if ~isempty(inc)
        if ~iscell(inc), inc = {inc}; end
        fprintf('  include (%s):\n', mode);
        for i = 1:numel(inc)
            inc_masks{end+1} = nim_roi_mask(nim, inc(i), dil); %#ok<AGROW>
        end
    end
    exc_mask = [];
    if ~isempty(exc)
        fprintf('  exclude:\n');
        exc_mask = nim_roi_mask(nim, exc, dil);
    end

    % Endpoint pair + containment corridor: how the ISMRM 2015 scorer actually
    % defines a bundle. Both are ordinary set tests, but on different things -
    % one on the two END points only, one on every point - and neither is
    % expressible as a waypoint include, which is why they are separate keys.
    ep_masks = {};
    if ~isempty(eps)
        if ~iscell(eps), eps = {eps}; end
        if numel(eps) ~= 2
            error('nim_filter_tracks_roi:endpoints', ...
                'endpoints_in needs exactly two regions (head and tail); got %d.', numel(eps));
        end
        fprintf('  endpoints (head/tail):\n');
        ep_masks = {nim_roi_mask(nim, eps(1), dil), nim_roi_mask(nim, eps(2), dil)};
    end
    con_mask = [];
    if ~isempty(con)
        fprintf('  contained in:\n');
        con_mask = nim_roi_mask(nim, con, dil);
    end
    any_mask = [];
    if ~isempty(anyi)
        fprintf('  any (touch at least one point):\n');
        any_mask = nim_roi_mask(nim, anyi, dil);
    end
    stats.n_dropped_endpoints = 0;
    stats.n_dropped_contained = 0;
    stats.n_dropped_any       = 0;
    stats.n_dropped_length    = 0;

    % Voxel -> mm. The ISMRM 2015 affine is diagonal 2 mm isotropic; a rotated
    % affine would make a per-axis scale meaningless, so refuse rather than
    % report a length that is quietly wrong.
    mm = [1 1 1];
    if have_len
        if isfield(nim, 'hdr') && isfield(nim.hdr, 'PixelDimensions')
            mm = double(nim.hdr.PixelDimensions(1:3));
        else
            error('nim_filter_tracks_roi:noVoxelSize', ...
                ['Length criteria are specified in millimetres but the nim carries no ' ...
                 'PixelDimensions, so voxel coordinates cannot be converted.']);
        end
    end

    keep = true(1, numel(tracks));
    for t = 1:numel(tracks)
        tr = tracks{t};
        if isempty(tr) || size(tr, 1) < 2, keep(t) = false; continue; end

        li = track_linear_indices(tr, dims);
        if isempty(li), keep(t) = false; continue; end

        if ~isempty(exc_mask) && any(exc_mask(li))
            keep(t) = false;
            stats.n_dropped_exclude = stats.n_dropped_exclude + 1;
            continue;
        end

        if ~isempty(con_mask)
            % EVERY point, so test the full index list rather than the unique
            % set - a track leaving the corridor once is out.
            if ~all(con_mask(li))
                keep(t) = false;
                stats.n_dropped_contained = stats.n_dropped_contained + 1;
                continue;
            end
        end

        if ~isempty(ep_masks)
            e1 = voxel_index(tr(1,:),   dims);
            e2 = voxel_index(tr(end,:), dims);
            if isempty(e1) || isempty(e2)
                ok_ends = false;
            else
                ok_ends = (ep_masks{1}(e1) && ep_masks{2}(e2)) || ...
                          (ep_masks{2}(e1) && ep_masks{1}(e2));
            end
            if ~ok_ends
                keep(t) = false;
                stats.n_dropped_endpoints = stats.n_dropped_endpoints + 1;
                continue;
            end
        end

        if ~isempty(any_mask) && ~any(any_mask(li))
            keep(t) = false;
            stats.n_dropped_any = stats.n_dropped_any + 1;
            continue;
        end

        if have_len && ~length_ok(tr, mm, Lspec)
            keep(t) = false;
            stats.n_dropped_length = stats.n_dropped_length + 1;
            continue;
        end

        if ~isempty(inc_masks)
            hits = false(1, numel(inc_masks));
            for k = 1:numel(inc_masks)
                hits(k) = any(inc_masks{k}(li));
            end
            if strcmp(mode, 'any'), ok = any(hits); else, ok = all(hits); end
            if ~ok
                keep(t) = false;
                stats.n_dropped_include = stats.n_dropped_include + 1;
            end
        end
    end

    tracks = tracks(keep);
    stats.keep = keep;      % so callers can subset per-track metadata in step
    stats.n_out = numel(tracks);
    stats.applied = true;

    fprintf(['  kept %d / %d tracks (dropped: %d include, %d exclude, %d containment, ' ...
             '%d endpoints, %d any, %d length)\n'], ...
        stats.n_out, stats.n_in, stats.n_dropped_include, stats.n_dropped_exclude, ...
        stats.n_dropped_contained, stats.n_dropped_endpoints, ...
        stats.n_dropped_any, stats.n_dropped_length);
    if stats.n_out == 0
        warning('nim_filter_tracks_roi:empty', ...
            ['ROI filtering removed every track. Check the include/exclude regions, ' ...
             'or raise roi_filter_dilate - the JHU parcellation is thin and a track ' ...
             'can pass beside a region without entering a labelled voxel.']);
    end
end

% =========================================================================
function li = track_linear_indices(tr, dims)
% Voxel indices visited by a track. Tracks are in voxel coordinates, so round to
% the nearest voxel centre and drop anything outside the volume.
    v = round(tr);
    inb = v(:,1) >= 1 & v(:,1) <= dims(1) & ...
          v(:,2) >= 1 & v(:,2) <= dims(2) & ...
          v(:,3) >= 1 & v(:,3) <= dims(3);
    v = v(inb, :);
    if isempty(v), li = []; return; end
    li = unique(sub2ind(dims, v(:,1), v(:,2), v(:,3)));
end

function li = voxel_index(p, dims)
% Linear index of the voxel a single point falls in, or [] if outside.
    v = round(p);
    if any(v < 1) || v(1) > dims(1) || v(2) > dims(2) || v(3) > dims(3)
        li = []; return;
    end
    li = sub2ind(dims, v(1), v(2), v(3));
end

function v = get_opt(options, name, dflt)
    if isfield(options, name) && ~isempty(options.(name))
        v = options.(name);
    else
        v = dflt;
    end
end

function ok = length_ok(tr, mm, L)
% The scorer's length criteria, in millimetres, matching scilpy's
% filter_streamlines_by_total_length_per_dim:
%   length_<d>      -> |sum(d)|   net displacement; a loop cancels itself out
%   length_<d>_abs  -> sum(|d|)   total travel; a loop adds to it
    ok = true;
    d = diff(tr, 1, 1) .* mm(:)';           % per-step displacement, mm
    if isempty(d), ok = false; return; end

    if ~isempty(L.total)
        tot = sum(sqrt(sum(d.^2, 2)));
        if tot < L.total(1) || tot > L.total(2), ok = false; return; end
    end
    net = abs(sum(d, 1));
    tvl = sum(abs(d), 1);
    for k = 1:3
        r = L.net{k};
        if ~isempty(r) && (net(k) < r(1) || net(k) > r(2)), ok = false; return; end
        r = L.abs{k};
        if ~isempty(r) && (tvl(k) < r(1) || tvl(k) > r(2)), ok = false; return; end
    end
end

function v = num_opt(options, name)
% A [min max] pair that may arrive as a numeric vector or as a cell of strings
% (which is what a --set override parses to).
    v = [];
    if ~isfield(options, name) || isempty(options.(name)), return; end
    r = options.(name);
    if iscell(r), r = cellfun(@(x) str2double(string(x)), r); end
    r = double(r(:))';
    if numel(r) ~= 2 || any(isnan(r))
        error('nim_filter_tracks_roi:badLength', ...
            '%s must be a [min max] pair in millimetres.', name);
    end
    v = r;
end
