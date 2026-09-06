function P = nim_parcellation_from_masks(mask_dirs, ref_nii, target_dims)
% nim_parcellation_from_masks: Build a parcellation from a directory of binary
% mask volumes, resampled into a target voxel grid.
%
%   P = nim_parcellation_from_masks({'.../all_masks'}, 'dwi.nii.gz', size(nim.FA))
%
% Returns a struct with
%   .labels    - integer label volume on the target grid (winner-takes-all)
%   .map       - containers.Map, label index -> region name
%   .masks     - containers.Map, region name -> logical volume (the TRUE masks,
%                overlaps intact)
%   .overlap   - report on how much the source masks overlap each other
%
% WHY BOTH .labels AND .masks. Anatomical bundles genuinely overlap: a voxel in
% the temporal stem belongs to the uncinate AND the inferior longitudinal
% fasciculus AND, near midline, the corpus callosum. An integer label volume
% cannot represent that - every voxel gets exactly one owner - so building one
% silently deletes the overlaps. .labels exists because a lot of code expects a
% parcellation to look like FreeSurfer's; .masks is what ROI selection should
% actually use, because it is the data as given.
%
% Overlap in .labels is resolved SMALLEST-REGION-WINS. The alternative,
% first-come or largest-wins, lets a big structure like the corpus callosum
% swallow every small bundle that crosses it, which is the opposite of useful:
% the small specific bundle is the one you wanted to name.
%
% Resampling is nearest-neighbour through the two sforms, so it works when the
% masks are on a different grid from the DWI - which they are here: the ISMRM
% ground truth is 1mm 180x216x180 and the DWI is 2mm 90x108x90.

    if ischar(mask_dirs) || isstring(mask_dirs), mask_dirs = {char(mask_dirs)}; end

    A_ref = nim_nifti_affine(ref_nii);
    files = [];
    for d = 1:numel(mask_dirs)
        files = [files; dir(fullfile(mask_dirs{d}, '*.nii.gz'))]; %#ok<AGROW>
        files = [files; dir(fullfile(mask_dirs{d}, '*.nii'))];    %#ok<AGROW>
    end
    if isempty(files)
        error('nim_parcellation_from_masks:empty', ...
            'No .nii or .nii.gz masks under %s', strjoin(mask_dirs, ', '));
    end

    % Target voxel centres, expressed as indices into each source mask.
    [I, J, K] = ndgrid(1:target_dims(1), 1:target_dims(2), 1:target_dims(3));
    V = [I(:)-1, J(:)-1, K(:)-1, ones(numel(I),1)]';     % 0-based homogeneous
    world = A_ref * V;

    P = struct();
    P.labels = zeros(target_dims, 'uint16');
    P.map    = containers.Map('KeyType','double','ValueType','char');
    P.masks  = containers.Map('KeyType','char','ValueType','any');
    sizes    = zeros(numel(files),1);
    vols     = cell(numel(files),1);
    names    = cell(numel(files),1);

    for f = 1:numel(files)
        fp = fullfile(files(f).folder, files(f).name);
        names{f} = erase(erase(files(f).name, '.nii.gz'), '.nii');
        A_src = nim_nifti_affine(fp);
        M = niftiread(fp) > 0;
        sd = size(M);
        Q = round(A_src \ world)' + 1;                    % target centres -> source voxels
        ok = all(Q >= 1, 2) & Q(:,1) <= sd(1) & Q(:,2) <= sd(2) & Q(:,3) <= sd(3);
        v = false(prod(target_dims), 1);
        v(ok) = M(sub2ind(sd, Q(ok,1), Q(ok,2), Q(ok,3)));
        vols{f}  = reshape(v, target_dims);
        sizes(f) = sum(v);
        P.masks(names{f}) = vols{f};
    end

    % Smallest region wins: paint in descending size order so small regions,
    % painted last, survive.
    [~, order] = sort(sizes, 'descend');
    for oi = 1:numel(order)
        f = order(oi);
        P.labels(vols{f}) = uint16(f);
    end
    for f = 1:numel(files)
        P.map(f) = names{f};
    end

    % How much did the label volume have to throw away?
    stack = cat(4, vols{:});
    n_per = sum(stack, 4);
    P.overlap = struct( ...
        'n_regions',        numel(files), ...
        'voxels_labelled',  sum(n_per(:) > 0), ...
        'voxels_multi',     sum(n_per(:) > 1), ...
        'max_regions_per_voxel', max(n_per(:)), ...
        'region_sizes',     containers.Map(names, num2cell(sizes)));

    fprintf('parcellation from masks: %d regions, %d labelled voxels, %d in >1 region (%.1f%%), max %d regions on one voxel\n', ...
        P.overlap.n_regions, P.overlap.voxels_labelled, P.overlap.voxels_multi, ...
        100*P.overlap.voxels_multi/max(P.overlap.voxels_labelled,1), P.overlap.max_regions_per_voxel);
end
