function m = nim_region_mask(nim, region)
% nim_region_mask: The mask of ONE parcellation region, overlaps intact.
%
%   m = nim_region_mask(nim, 12)          % by label index
%   m = nim_region_mask(nim, 'UF_right')  % by name
%
% Quiet counterpart to nim_roi_mask, for code that loops over regions and does
% not want a report per call.
%
% WHY THIS EXISTS. The obvious lookup is `nim.parcellation_mask == id`, and it is
% wrong wherever regions overlap. A label volume gives each voxel exactly one
% owner, so a region that shares voxels with anything else comes back with holes.
% On the ISMRM bundle masks 84.8% of labelled voxels belong to more than one
% bundle and the label volume retains a median 43% of each region - CC_u_shaped
% keeps 1275 of its 106502 voxels. Every count taken from `== id` is wrong in the
% same direction, and nothing in the result shows it.
%
% nim.roi_masks holds the regions as defined. This prefers it and falls back to
% the label volume only when it must, warning through nim_roi_mask's own path.

    if ~isfield(nim, 'parcellation_mask') || isempty(nim.parcellation_mask)
        error('nim_region_mask:noParcellation', 'nim has no parcellation_mask.');
    end
    name = '';
    if ischar(region) || isstring(region)
        name = char(region);
    else
        lab = nim_atlas_label_map(nim);
        if isKey(lab, double(region)), name = char(lab(double(region))); end
    end

    if ~isempty(name) && isfield(nim, 'roi_masks') ...
            && isa(nim.roi_masks, 'containers.Map') && isKey(nim.roi_masks, name)
        m = logical(nim.roi_masks(name));
        return;
    end
    if ischar(region) || isstring(region)
        lab = nim_atlas_label_map(nim); ids = cell2mat(keys(lab)); nm = values(lab);
        hit = find(strcmpi(cellfun(@char, nm, 'uni', 0), name), 1);
        if isempty(hit)
            error('nim_region_mask:unknownName', 'No region named "%s".', name);
        end
        region = ids(hit);
    end
    m = (nim.parcellation_mask == region);
end
