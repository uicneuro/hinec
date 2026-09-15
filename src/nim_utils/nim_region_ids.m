function ids = nim_region_ids(nim)
% nim_region_ids: The parcellation's region indices, from the label MAP.
%
%   ids = nim_region_ids(nim)
%
% WHY NOT unique(nim.parcellation_mask(:)). Because a region can be entirely
% absent from the label volume. Regions overlap, the label volume keeps one owner
% per voxel, and a region every one of whose voxels is claimed by a smaller
% neighbour has no voxels left in it. On the ISMRM bundle masks the label volume
% retains a median 43% of each region, so listing regions by what survives there
% under-reports the parcellation and can drop a region completely.
%
% The label map is the definitive list of what the parcellation contains.

    if isfield(nim, 'atlas_labels') && isstruct(nim.atlas_labels) ...
            && isfield(nim.atlas_labels, 'map') ...
            && isa(nim.atlas_labels.map, 'containers.Map') ...
            && nim.atlas_labels.map.Count > 0
        ids = sort(cell2mat(keys(nim.atlas_labels.map)));
        ids = ids(ids > 0);
        ids = ids(:);
        return;
    end
    if ~isfield(nim, 'parcellation_mask') || isempty(nim.parcellation_mask)
        ids = zeros(0, 1); return;
    end
    ids = unique(nim.parcellation_mask(:));
    ids = double(ids(ids > 0));
end
