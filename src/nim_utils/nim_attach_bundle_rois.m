function nim = nim_attach_bundle_rois(nim, roi_root, ref_nii)
% nim_attach_bundle_rois: Replace the atlas parcellation with one built from a
% directory of bundle masks, and attach the masks themselves for ROI selection.
%
%   nim = nim_attach_bundle_rois(nim, 'data/ismrm2015/scoring_data_Renauld2023/ROI', ...
%                                'data/ismrm2015/ismrm2015_dwi_ref.nii.gz');
%   nim_save(nim, 'data/ismrm2015/ismrm2015.mat');
%
% roi_root must contain an `all_masks/` directory (the bundle containment
% corridors, which become the parcellation labels) and may contain `endpoints/`
% and `any_masks/` (gates, addressable by name but deliberately NOT labels).
%
% WHY THIS EXISTS. Seeding from an atlas label and scoring against a bundle of
% the same name is not a like-for-like comparison: JHU label 47, "Uncinate
% fasciculus R", is 24 voxels on the 2mm grid against the ISMRM bundle-density
% mask's 1503, a Dice of 0.018. Atlas labels are deep white-matter cores for
% tract statistics; a bundle is an endpoint pair plus a containment corridor.
% This attaches the challenge's own regions so ROIs can be addressed by the names
% the scorer uses.
%
% The previous parcellation is preserved as nim.parcellation_mask_<atlas_type>
% rather than discarded - it remains the right atlas for tract-core statistics.
%
% Both a winner-takes-all label volume and the true masks are stored, because
% bundles genuinely overlap (84.8% of labelled voxels belong to more than one on
% the ISMRM data, up to ten on a single voxel) and a label volume can keep only
% one owner per voxel. nim_roi_mask prefers nim.roi_masks, so a request for a
% bundle returns all of it rather than the part no larger structure claimed.

    if nargin < 3 || isempty(ref_nii)
        error('nim_attach_bundle_rois:ref', ...
            ['A reference NIfTI defining the nim''s voxel grid is required - the ' ...
             'masks are resampled into it through both sforms.']);
    end
    all_dir = fullfile(roi_root, 'all_masks');
    if ~isfolder(all_dir)
        error('nim_attach_bundle_rois:layout', 'No all_masks/ under %s', roi_root);
    end

    dims = size(nim.FA);
    P = nim_parcellation_from_masks({all_dir}, ref_nii, dims);

    % Keep the outgoing parcellation under its atlas name instead of losing it.
    if isfield(nim, 'parcellation_mask') && ~isempty(nim.parcellation_mask)
        prev = 'atlas';
        if isfield(nim, 'atlas_type') && ~isempty(nim.atlas_type)
            prev = lower(char(nim.atlas_type));
        end
        nim.(sprintf('parcellation_mask_%s', prev)) = nim.parcellation_mask;
        if isfield(nim, 'atlas_labels')
            nim.(sprintf('atlas_labels_%s', prev)) = nim.atlas_labels;
        end
        fprintf('  previous parcellation preserved as parcellation_mask_%s\n', prev);
    end

    nim.parcellation_mask = P.labels;
    nim.atlas_labels = struct('map', P.map, 'atlas_type', 'bundle_masks', ...
                              'atlas_variant', all_dir);
    nim.atlas_type   = 'bundle_masks';

    % Gates are addressable by name but are not parcellation labels.
    nim.roi_masks = containers.Map('KeyType','char','ValueType','any');
    for k = keys(P.masks), nim.roi_masks(k{1}) = P.masks(k{1}); end
    for sub = {'endpoints', 'any_masks'}
        d = fullfile(roi_root, sub{1});
        if ~isfolder(d), continue; end
        G = nim_parcellation_from_masks({d}, ref_nii, dims);
        for k = keys(G.masks), nim.roi_masks(k{1}) = G.masks(k{1}); end
    end
    nim.roi_source = roi_root;

    fprintf('  %d bundle labels, %d addressable ROIs\n', P.map.Count, nim.roi_masks.Count);
end
