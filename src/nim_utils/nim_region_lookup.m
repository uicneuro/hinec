function R = nim_region_lookup(nim)
% nim_region_lookup: Region membership that allows a voxel to belong to several
% regions at once.
%
%   R = nim_region_lookup(nim);
%   ids = R.at(pos);            % region ids containing voxel pos ([] if none)
%   ids = R.touched(track);     % every region an Nx3 track passes through
%
% WHY NOT `parcellation_mask(pos)`. That returns ONE label, because a label
% volume has one owner per voxel. Anatomy does not: on the ISMRM bundle masks
% 84.8% of labelled voxels lie in two or more bundles and one voxel lies in ten.
% A streamline point in a voxel shared by the corpus callosum and the cingulum is
% in both, and code that assigns it a single label silently attributes it to
% whichever region happened to win the tie - which, for a connectivity matrix, is
% an edge credited to the wrong pair.
%
% Membership is built once from nim.roi_masks (the regions as defined). Where a
% region has no true mask, its label-volume extent is used and R.exact is false.

    dims = size(nim.parcellation_mask);
    lab  = nim_atlas_label_map(nim);
    ids  = cell2mat(keys(lab));
    names = cellfun(@char, values(lab), 'uni', 0);

    have = isfield(nim, 'roi_masks') && isa(nim.roi_masks, 'containers.Map');
    masks = cell(1, numel(ids));
    exact = true(1, numel(ids));
    for i = 1:numel(ids)
        if have && isKey(nim.roi_masks, names{i})
            masks{i} = logical(nim.roi_masks(names{i}));
        else
            masks{i} = (nim.parcellation_mask == ids(i));
            exact(i) = false;
        end
    end

    % Flatten to a per-voxel membership list so lookups are cheap.
    NV = prod(dims);
    cnt = zeros(NV, 1, 'uint8');
    for i = 1:numel(ids), cnt(masks{i}(:)) = cnt(masks{i}(:)) + 1; end

    R = struct();
    R.ids = ids; R.names = {names{:}}; R.dims = dims; R.exact = exact;
    R.n_multi = sum(cnt > 1); R.n_labelled = sum(cnt > 0);
    R.masks = masks;
    R.any   = reshape(cnt > 0, dims);   % voxels belonging to ANY region
    R.count = reshape(cnt, dims);       % how many regions own each voxel
    R.at = @(p) region_at(masks, ids, dims, p);
    R.touched = @(t) region_touched(masks, ids, dims, t);
end

% =========================================================================
function out = region_at(masks, ids, dims, p)
    v = round(p(:))';
    if any(v < 1) || any(v > dims), out = []; return; end
    li = sub2ind(dims, v(1), v(2), v(3));
    keep = false(1, numel(ids));
    for i = 1:numel(ids), keep(i) = masks{i}(li); end
    out = ids(keep);
end

function out = region_touched(masks, ids, dims, track)
    v = round(track);
    ok = all(v >= 1, 2) & v(:,1) <= dims(1) & v(:,2) <= dims(2) & v(:,3) <= dims(3);
    v = v(ok, :);
    if isempty(v), out = []; return; end
    li = unique(sub2ind(dims, v(:,1), v(:,2), v(:,3)));
    keep = false(1, numel(ids));
    for i = 1:numel(ids), keep(i) = any(masks{i}(li)); end
    out = ids(keep);
end
