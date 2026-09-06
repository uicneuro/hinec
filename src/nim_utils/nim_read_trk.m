function [tracks, hdr] = nim_read_trk(trk_file, target_nii)
% nim_read_trk: Read a TrackVis .trk file, optionally into another volume's
% voxel grid.
%
%   [tracks, hdr] = nim_read_trk(f)              % streamlines in the TRK's own voxel units
%   [tracks, hdr] = nim_read_trk(f, 'dwi.nii.gz') % resampled into that volume's voxel grid
%
% tracks is a cell array of Nx3 matrices, matching the convention used by the
% trackers in this repo (1-based voxel coordinates when target_nii is given).
%
% WHY THE SECOND ARGUMENT MATTERS. The ISMRM 2015 ground-truth bundles live on a
% 1mm 180x216x180 grid; this pipeline's DWI is 2mm 90x108x90. Plotting the two
% together without a transform silently draws the ground truth at half scale in
% the corner of the volume, which looks like a catastrophic tracking failure and
% is in fact a units bug. The mapping goes through world (RAS) space, which is
% the only frame both files agree on:
%
%   TRK stores points in "voxmm" - millimetres in ITS voxel frame, so its voxel
%   index is point/voxel_size - and carries a vox_to_ras matrix. The target NIfTI
%   carries an srow (sform). So:
%
%       ras          = trk_vox_to_ras * [point ./ trk_voxel_size; 1]
%       target_vox0  = srow \ [ras; 1]                  (0-based)
%       tracks       = target_vox0 + 1                  (MATLAB is 1-based)

    fid = fopen(trk_file, 'r', 'l');
    if fid < 0, error('nim_read_trk:open', 'Cannot open %s', trk_file); end
    cleanup = onCleanup(@() fclose(fid));

    raw = fread(fid, 1000, '*uint8');
    hdr = struct();
    hdr.dim         = double(typecast(raw(7:12),   'int16'))';
    hdr.voxel_size  = double(typecast(raw(13:24),  'single'))';
    hdr.n_scalars   = double(typecast(raw(37:38),  'int16'));
    hdr.n_properties= double(typecast(raw(239:240),'int16'));
    hdr.vox_to_ras  = reshape(double(typecast(raw(441:504), 'single')), 4, 4)';
    hdr.voxel_order = char(raw(949:951))';
    hdr.n_count     = double(typecast(raw(989:992), 'int32'));

    % Stream the streamlines. n_count can be 0 in files written by tools that did
    % not backfill the header, so read until EOF rather than trusting it.
    tracks = {};
    while true
        n = fread(fid, 1, 'int32');
        if isempty(n) || n <= 0, break; end
        blk = fread(fid, n * (3 + hdr.n_scalars), 'single');
        if numel(blk) < n * (3 + hdr.n_scalars), break; end
        blk = reshape(blk, 3 + hdr.n_scalars, n)';
        tracks{end+1} = blk(:, 1:3); %#ok<AGROW>
        if hdr.n_properties > 0, fread(fid, hdr.n_properties, 'single'); end
    end
    tracks = tracks(:);

    if nargin < 2 || isempty(target_nii), return; end

    if isnumeric(target_nii)
        srow = target_nii;                       % caller supplied the affine directly
    else
        srow = nim_nifti_affine(target_nii);
    end
    vs   = hdr.voxel_size; vs(vs == 0) = 1;

    for i = 1:numel(tracks)
        p    = tracks{i} ./ vs;                       % voxmm -> TRK voxel index
        ras  = (hdr.vox_to_ras * [p, ones(size(p,1),1)]')';
        vox0 = (srow \ ras')';
        tracks{i} = vox0(:, 1:3) + 1;                 % 0-based -> MATLAB 1-based
    end
end
