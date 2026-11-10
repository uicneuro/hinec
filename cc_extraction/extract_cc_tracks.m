function result = extract_cc_tracks(data_prefix, config_file, out_dir, varargin)
% EXTRACT_CC_TRACKS  Corpus-callosum tractography with the real HINEC kernels, without FSL/SPM.
%
%   result = extract_cc_tracks(data_prefix, config_file, out_dir)
%   result = extract_cc_tracks(..., 'mask_file', 'cc_extraction/sample_M.nii.gz')
%
% Pipeline (same call sequence as main.m + runTractography.m, minus the FSL/SPM steps):
%   nim_read -> nim_dt_spd -> nim_eig -> nim_fa
%     -> FA-based WM/GM/CSF masks (same adaptive thresholds as preproc_tissue_segmentation.m,
%        with MATLAB imerode(…, ones(3,3,3)) standing in for `fslmaths -ero`, whose default kernel is that 3x3x3 box)
%     -> seed mask = brain mask & FA > 0.05   (Strategy 1 of runTractography.m)
%     -> nim_tractography_hinec(nim, options)  (options straight from the YAML preset, ACT enabled)
%     -> corpus-callosum ROI: midsagittal slab (+-1 voxel around the brain midline), FA > 0.4,
%        |v1 . left-right| > 0.8, largest 26-connected component
%     -> keep tracks that pass through the ROI AND cross the midsagittal plane.
%
% Why not just main.m?  main.m calls preproc_brain_extraction (bet), preproc_mask_improvement (fslmaths, spm_vol)
% and nim_parcellation (flirt) unconditionally, so it errors without FSL/SPM. This script takes the brain mask from
% a file instead (the one produced by the Python run is next to this script) and never touches the atlas.
%
% Inputs
%   data_prefix : path without extension; expects <prefix>.nii.gz, <prefix>.bval, <prefix>.bvec
%   config_file : YAML preset, e.g. 'config/hinec_default.yml'
%   out_dir     : output folder (created); writes tracks_cc.mat, tracks_hinec_wholebrain.mat, nim_cc.mat,
%                 cc_roi.nii.gz, cc_summary.txt and figures/ (via visualizeTractography)
% Name/value options
%   'mask_file'      brain mask NIfTI (default: <prefix>_M.nii.gz)
%   'roi_fa_thr'     FA threshold for CC ROI voxels            (default 0.4)
%   'roi_lr_thr'     |v1_LR| threshold for CC ROI voxels       (default 0.8)
%   'roi_halfwidth'  slab half-width in voxels around midline  (default 1)
%   'save_wholebrain' also save all whole-brain tracks (large) (default true)
%   'make_figures'   call visualizeTractography for PNG export (default true)
%
% Example (run from the repository root):
%   addpath(genpath('.'));
%   result = extract_cc_tracks('data/original_sample/sample', 'config/hinec_default.yml', ...
%                              'cc_extraction/matlab', 'mask_file', 'cc_extraction/sample_M.nii.gz');
%
% NOTE: This script was written alongside the Python port (cc_extraction/hinec_py) and could not be executed
%       here because MATLAB is not available in the authoring environment. It only calls existing HINEC
%       functions with their documented signatures, but please treat it as untested until the first run.

    %% ---- options -------------------------------------------------------------------------------------------
    p = inputParser;
    addParameter(p, 'mask_file', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'roi_fa_thr', 0.4, @isnumeric);
    addParameter(p, 'roi_lr_thr', 0.8, @isnumeric);
    addParameter(p, 'roi_halfwidth', 1, @isnumeric);
    addParameter(p, 'save_wholebrain', true, @islogical);
    addParameter(p, 'make_figures', true, @islogical);
    parse(p, varargin{:});
    o = p.Results;

    data_prefix = char(data_prefix);
    config_file = char(config_file);
    out_dir     = char(out_dir);
    mask_file   = char(o.mask_file);
    if isempty(mask_file)
        mask_file = [data_prefix '_M.nii.gz'];
    end
    if ~isfile(mask_file)
        error('Brain mask not found: %s  (pass ''mask_file'' or run the Python pipeline first)', mask_file);
    end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    fig_dir = fullfile(out_dir, 'figures');

    t_all = tic;
    fprintf('=== extract_cc_tracks ===\nprefix: %s\nconfig: %s\nmask:   %s\nout:    %s\n', data_prefix, config_file, mask_file, out_dir);

    %% ---- configuration -------------------------------------------------------------------------------------
    config = load_config_yaml(config_file);
    options = config.tractography;

    %% ---- DTI stage (identical calls to main.m) --------------------------------------------------------------
    nim = nim_read(data_prefix, 'MaskPath', mask_file);
    nim = nim_dt_spd(nim);
    nim = nim_eig(nim);
    nim = nim_fa(nim);

    %% ---- tissue masks for ACT (preproc_tissue_segmentation.m logic, no FSL) --------------------------------
    bm = nim.mask > 0.5;
    fa_in = nim.FA(bm);
    wm_thr  = max(prctile(fa_in, 60), 0.1);
    csf_thr = min(prctile(fa_in, 20), 0.1);
    if wm_thr <= csf_thr, wm_thr = csf_thr + 0.05; end
    wm_prelim = (nim.FA > wm_thr) & bm;
    wm  = imerode(wm_prelim, ones(3, 3, 3));                 % fslmaths -ero default kernel = 3x3x3 box
    gm  = (nim.FA > csf_thr) & (nim.FA <= wm_thr) & bm;
    csf = (nim.FA <= csf_thr) & bm;
    nim.wm_mask = double(wm); nim.gm_mask = double(gm); nim.csf_mask = double(csf);
    fprintf('Tissue masks: WM thr %.3f, CSF thr %.3f | WM %.1f%% (eroded), GM %.1f%%, CSF %.1f%% of brain\n', ...
        wm_thr, csf_thr, 100*nnz(wm)/nnz(bm), 100*nnz(gm)/nnz(bm), 100*nnz(csf)/nnz(bm));

    %% ---- seeding + ACT options (runTractography.m, Strategy 1) ---------------------------------------------
    seed_fa_thr = 0.05;
    if isfield(options, 'seed_fa_threshold') && ~isempty(options.seed_fa_threshold)
        seed_fa_thr = options.seed_fa_threshold;
    end
    options.seed_mask = bm & (nim.FA > seed_fa_thr);
    options.wm_mask  = nim.wm_mask;
    options.gm_mask  = nim.gm_mask;
    options.csf_mask = nim.csf_mask;
    fprintf('Seed voxels: %d (%.1f%% of volume), seed_density %d\n', nnz(options.seed_mask), ...
        100*nnz(options.seed_mask)/numel(nim.FA), options.seed_density);

    %% ---- tractography (the real kernel) ---------------------------------------------------------------------
    algorithm = lower(char(options.algorithm));
    tic;
    if strcmp(algorithm, 'hinec')
        tracks = nim_tractography_hinec(nim, options);
    else
        tracks = nim_tractography_standard(nim, options);
    end
    elapsed_time = toc;
    fprintf('Tractography: %d tracks in %.1f s\n', numel(tracks), elapsed_time);

    %% ---- corpus-callosum ROI ---------------------------------------------------------------------------------
    dims = size(nim.FA);
    lr_axis = detect_lr_axis(nim);                            % voxel axis that is left-right (1 for this sample)
    [sx, sy, sz] = ind2sub(dims, find(bm));
    sub = [sx sy sz];
    x_mid = mean(sub(:, lr_axis));                            % midline as brain-mask centre of mass (1-based)
    v1_lr = squeeze(nim.evec(:, :, :, lr_axis, 1));           % LR component of the primary eigenvector
    slab = false(dims);
    idx_lr = find(abs((1:dims(lr_axis)) - x_mid) <= o.roi_halfwidth + 0.5);
    switch lr_axis
        case 1, slab(idx_lr, :, :) = true;
        case 2, slab(:, idx_lr, :) = true;
        case 3, slab(:, :, idx_lr) = true;
    end
    cand = slab & bm & (nim.FA > o.roi_fa_thr) & (abs(v1_lr) > o.roi_lr_thr);
    CC = bwconncomp(cand, 26);
    if CC.NumObjects == 0
        error('No CC candidate voxels; relax roi_fa_thr / roi_lr_thr');
    end
    [~, imax] = max(cellfun(@numel, CC.PixelIdxList));
    roi = false(dims);
    roi(CC.PixelIdxList{imax}) = true;
    fprintf('CC ROI: %d voxels (largest of %d components), midline at %.2f along axis %d\n', nnz(roi), CC.NumObjects, x_mid, lr_axis);

    %% ---- select callosal tracks -------------------------------------------------------------------------------
    keep = false(numel(tracks), 1);
    for i = 1:numel(tracks)
        t = tracks{i};
        r = max(1, min(round(t), dims));
        hits = any(roi(sub2ind(dims, r(:, 1), r(:, 2), r(:, 3))));
        crosses = min(t(:, lr_axis)) < x_mid && max(t(:, lr_axis)) > x_mid;
        keep(i) = hits && crosses;
    end
    tracks_cc = tracks(keep);
    len_mm = cellfun(@(t) sum(sqrt(sum((diff(t, 1, 1) .* nim.hdr.PixelDimensions(1:3)).^2, 2))), tracks_cc);
    fprintf('Corpus callosum: %d tracks (pass through ROI and cross midline); length mean %.1f mm, median %.1f mm\n', ...
        numel(tracks_cc), mean(len_mm), median(len_mm));

    %% ---- save ------------------------------------------------------------------------------------------------
    nim.parcellation_mask = int32(roi);                       % lets visualizeTractography(..., 'region', 1) work
    nim.labels = {'Corpus callosum (midsagittal ROI)'};
    nim.cc_roi_info = struct('x_mid', x_mid, 'lr_axis', lr_axis, 'fa_thr', o.roi_fa_thr, 'lr_thr', o.roi_lr_thr, ...
                             'halfwidth', o.roi_halfwidth, 'roi_voxels', nnz(roi));
    nim_file = fullfile(out_dir, 'nim_cc.mat');
    nim_save(nim, nim_file);

    tracks_all = tracks; %#ok<NASGU>
    tracks = tracks_cc;
    cc_file = fullfile(out_dir, 'tracks_cc.mat');
    save(cc_file, 'tracks', 'options', 'elapsed_time', 'algorithm', '-v7.3');
    if o.save_wholebrain
        tracks = tracks_all;
        save(fullfile(out_dir, 'tracks_hinec_wholebrain.mat'), 'tracks', 'options', 'elapsed_time', 'algorithm', '-v7.3');
        tracks = tracks_cc;
    end
    try
        info = niftiinfo(mask_file);
        info.Datatype = 'uint8'; info.BitsPerPixel = 8; info.ImageSize = dims; info.PixelDimensions = info.PixelDimensions(1:3);
        niftiwrite(uint8(roi), fullfile(out_dir, 'cc_roi'), info, 'Compressed', true);
    catch ME
        warning('Could not write cc_roi.nii.gz: %s', ME.message);
    end

    fid = fopen(fullfile(out_dir, 'cc_summary.txt'), 'w');
    fprintf(fid, 'extract_cc_tracks summary\n=========================\n');
    fprintf(fid, 'prefix: %s\nconfig: %s\nalgorithm: %s\n', data_prefix, config_file, algorithm);
    fprintf(fid, 'seed voxels: %d, seed_density: %d\n', nnz(options.seed_mask), options.seed_density);
    fprintf(fid, 'whole-brain tracks: %d (%.1f s)\n', numel(tracks_all), elapsed_time);
    fprintf(fid, 'CC ROI voxels: %d, midline %.2f (axis %d)\n', nnz(roi), x_mid, lr_axis);
    fprintf(fid, 'CC tracks: %d, length mean %.1f mm, median %.1f mm, max %.1f mm\n', numel(tracks_cc), mean(len_mm), median(len_mm), max(len_mm));
    fclose(fid);

    %% ---- figures (HINEC's own viewer) -----------------------------------------------------------------------
    if o.make_figures
        try
            visualizeTractography(cc_file, nim_file, 'export_dir', fig_dir, 'color_mode', 'direction', 'max_tracts', 5000);
            visualizeTractography(cc_file, nim_file, 'region', 1, 'export_dir', fig_dir, 'min_region_points', 1);
        catch ME
            warning('visualizeTractography failed (%s); tracks are saved regardless.', ME.message);
        end
    end

    result = struct('tracks_cc', {tracks_cc}, 'n_wholebrain', numel(tracks_all), 'roi', roi, 'x_mid', x_mid, ...
                    'lr_axis', lr_axis, 'len_mm', len_mm, 'nim_file', nim_file, 'tracks_file', cc_file);
    fprintf('Done in %.1f min. Outputs in %s\n', toc(t_all)/60, out_dir);
end

function lr_axis = detect_lr_axis(nim)
% Voxel axis whose direction is most left-right (world x), from the NIfTI transform returned by niftiinfo.
% affine3d stores T with [x y z 1] = [i j k 1]*T; affinetform3d (R2022b+) stores A with [x;y;z;1] = A*[i;j;k;1].
% Falls back to axis 1, which is correct for data/original_sample (LAS orientation).
    lr_axis = 1;
    try
        tf = nim.hdr.Transform;
        if isprop(tf, 'A')
            A = tf.A;
            [~, lr_axis] = max(abs(A(1, 1:3)));
        else
            T = tf.T;
            [~, lr_axis] = max(abs(T(1:3, 1)));
        end
    catch
        % keep default
    end
end
