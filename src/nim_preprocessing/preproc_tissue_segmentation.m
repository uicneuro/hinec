function [wm_mask_file, gm_mask_file, csf_mask_file] = preproc_tissue_segmentation(fa_data, brain_mask_file, file_prefix, t1_file)
% preproc_tissue_segmentation: Generate WM/GM/CSF tissue masks for ACT.
%
% PRIMARY method: FSL FAST on the anatomical T1, resampled into DWI space via the
% images' WORLD affines (flirt -usesqform) — no registration step, so it does not
% depend on (and is not corrupted by) the T1->DWI registration. This is the right
% approach when the T1 is already world-aligned to the diffusion (e.g. the ISMRM
% challenge T1). FALLBACK (no T1): FA-tertile binning, which bins anisotropy — NOT
% tissue — so ACT would terminate streamlines mid-crossing. Used only for DWI-only
% data where no anatomical T1 exists.
%
% Arguments:
%   fa_data          - FA array (nim.FA), used only by the FA-tertile fallback
%   brain_mask_file  - brain mask (DWI space); the resample grid + tissue domain
%   file_prefix      - prefix for output files
%   t1_file          - (optional) anatomical T1 (world-aligned to the DWI, any
%                      resolution); enables real FSL FAST tissue segmentation
%
% Returns paths to the WM / GM / CSF masks.
%
% Reference: PIPELINE.md Step 7 (Tractoflow-style anatomical priors for ACT)

% Define output file paths
wm_mask_file = [strrep(file_prefix, '_raw', '') '_WM_mask.nii.gz'];
gm_mask_file = [strrep(file_prefix, '_raw', '') '_GM_mask.nii.gz'];
csf_mask_file = [strrep(file_prefix, '_raw', '') '_CSF_mask.nii.gz'];

% PRIMARY: real anatomical segmentation from the T1 via FSL FAST.
if nargin >= 4 && ~isempty(t1_file) && isfile(t1_file)
    if tissue_from_t1_fast(t1_file, brain_mask_file, wm_mask_file, gm_mask_file, csf_mask_file)
        fprintf('  ✓ Tissue masks from FSL FAST on the anatomical T1 (world-aligned)\n');
        return;
    end
    fprintf('  ⚠ T1 FAST segmentation failed; falling back to FA-tertiles\n');
end

fprintf('Step: FA-based tissue segmentation for ACT (FALLBACK — no usable T1)...\n');

% Verify inputs
if isempty(fa_data)
    error('FA data is empty or not provided');
end

if ~isfile(brain_mask_file)
    error('Brain mask file not found: %s', brain_mask_file);
end

% Ensure FSL is available for morphological operations
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

fprintf('Loading brain mask...\n');

% Load brain mask
temp_mask_file = brain_mask_file;
cleanup_temp_mask = false;
if endsWith(brain_mask_file, '.nii.gz')
    temp_mask_file = gunzip(brain_mask_file, tempdir);
    temp_mask_file = temp_mask_file{1};
    cleanup_temp_mask = true;
end

V_mask = spm_vol(temp_mask_file);
brain_mask_data = spm_read_vols(V_mask);

if cleanup_temp_mask && isfile(temp_mask_file)
    delete(temp_mask_file);
end

% Get FA dimensions
fa_dims = size(fa_data);
mask_dims = size(brain_mask_data);

% Verify dimensions match
if ~isequal(fa_dims, mask_dims)
    error('FA data dimensions %s do not match brain mask dimensions %s', ...
          mat2str(fa_dims), mat2str(mask_dims));
end

% Convert brain mask to binary
brain_mask_binary = brain_mask_data > 0.5;

% Adaptive FA thresholds: fixed FA cutoffs (WM>0.2, CSF<=0.05) are tuned for
% real human-brain DTI. They misclassify low-FA datasets (synthetic phantoms,
% paediatric, low-SNR) as mostly CSF, which makes ACT terminate tracks
% aggressively. Use percentile-based thresholds within the brain mask so the
% partitioning adapts to whatever FA distribution the data actually has,
% with sanity floors so we don't drift too far from biology.
fa_in_brain = fa_data(brain_mask_binary);
adaptive_wm_thresh  = max(prctile(fa_in_brain, 60), 0.1);   % WM ≈ top 40% FA, floor 0.10
adaptive_csf_thresh = min(prctile(fa_in_brain, 20), 0.1);   % CSF ≈ bottom 20% FA, ceil 0.10
% Ensure ordering: wm_thresh > csf_thresh
if adaptive_wm_thresh <= adaptive_csf_thresh
    adaptive_wm_thresh = adaptive_csf_thresh + 0.05;
end

fprintf('Generating tissue masks using adaptive FA thresholds...\n');
fprintf('  FA range in brain: [%.3f, %.3f] (median %.3f)\n', ...
        min(fa_in_brain), max(fa_in_brain), median(fa_in_brain));
fprintf('  WM threshold:  FA > %.3f (60th percentile, floored at 0.10)\n', adaptive_wm_thresh);
fprintf('  GM range:      %.3f < FA <= %.3f\n', adaptive_csf_thresh, adaptive_wm_thresh);
fprintf('  CSF threshold: FA <= %.3f (20th percentile, capped at 0.10)\n', adaptive_csf_thresh);

%% Generate White Matter Mask
% WM: high-FA tissue, constrained to brain, eroded to remove boundary voxels
wm_preliminary = (fa_data > adaptive_wm_thresh) & brain_mask_binary;

% Count preliminary WM voxels
wm_prelim_count = sum(wm_preliminary(:));
fprintf('  Preliminary WM voxels: %d (%.1f%% of brain)\n', ...
        wm_prelim_count, 100*wm_prelim_count/sum(brain_mask_binary(:)));

% Save preliminary WM mask for erosion. SPM can't write .nii.gz directly,
% so write to .nii and let fslmaths (next step) read it.
temp_wm_prelim = [strrep(file_prefix, '_raw', '') '_temp_WM_prelim.nii'];
V_out = V_mask;
V_out.fname = temp_wm_prelim;
V_out.dt = [2 0]; % uint8 binary mask
spm_write_vol(V_out, double(wm_preliminary));

% Apply erosion using FSL to remove boundary voxels
% Use spherical structuring element with radius 1 voxel
fprintf('  Applying erosion to WM mask (removing boundary voxels)...\n');
cmd_erode = sprintf('%s/bin/fslmaths %s -ero %s', ...
                    fsl_path, temp_wm_prelim, wm_mask_file);

[status, cmdout] = system(cmd_erode);

if status ~= 0
    warning('FSL erosion failed, using preliminary WM mask:\n%s', cmdout);
    copyfile(temp_wm_prelim, wm_mask_file);
end

% Clean up temporary file
if isfile(temp_wm_prelim)
    delete(temp_wm_prelim);
end

% Load final WM mask for statistics
temp_wm_final = wm_mask_file;
cleanup_wm_final = false;
if endsWith(wm_mask_file, '.nii.gz')
    temp_wm_final = gunzip(wm_mask_file, tempdir);
    temp_wm_final = temp_wm_final{1};
    cleanup_wm_final = true;
end

V_wm = spm_vol(temp_wm_final);
wm_data = spm_read_vols(V_wm);

if cleanup_wm_final && isfile(temp_wm_final)
    delete(temp_wm_final);
end

wm_count = sum(wm_data(:) > 0.5);
fprintf('  ✓ Final WM voxels (after erosion): %d (%.1f%% of brain)\n', ...
        wm_count, 100*wm_count/sum(brain_mask_binary(:)));

%% Generate Gray Matter Mask
% GM: 0.05 < FA <= 0.2, constrained to brain
gm_mask = (fa_data > adaptive_csf_thresh) & (fa_data <= adaptive_wm_thresh) & brain_mask_binary;

gm_count = sum(gm_mask(:));
fprintf('  ✓ GM voxels: %d (%.1f%% of brain)\n', ...
        gm_count, 100*gm_count/sum(brain_mask_binary(:)));

% Save GM mask. SPM can't write .nii.gz directly — write .nii then gzip.
temp_gm_nii = strrep(gm_mask_file, '.nii.gz', '.nii');
V_out.fname = temp_gm_nii;
spm_write_vol(V_out, double(gm_mask));
gzip(temp_gm_nii);
if isfile(temp_gm_nii); delete(temp_gm_nii); end

%% Generate CSF Mask
% CSF: FA <= 0.05, constrained to brain
csf_mask = (fa_data <= adaptive_csf_thresh) & brain_mask_binary;

csf_count = sum(csf_mask(:));
fprintf('  ✓ CSF voxels: %d (%.1f%% of brain)\n', ...
        csf_count, 100*csf_count/sum(brain_mask_binary(:)));

% Save CSF mask. SPM can't write .nii.gz directly — write .nii then gzip.
temp_csf_nii = strrep(csf_mask_file, '.nii.gz', '.nii');
V_out.fname = temp_csf_nii;
spm_write_vol(V_out, double(csf_mask));
gzip(temp_csf_nii);
if isfile(temp_csf_nii); delete(temp_csf_nii); end

%% Validate tissue masks
fprintf('\nValidating tissue segmentation...\n');

% Check for tissue overlap (should be mutually exclusive)
wm_binary = wm_data > 0.5;
overlap_wm_gm = sum(wm_binary(:) & gm_mask(:));
overlap_wm_csf = sum(wm_binary(:) & csf_mask(:));
overlap_gm_csf = sum(gm_mask(:) & csf_mask(:));

total_overlap = overlap_wm_gm + overlap_wm_csf + overlap_gm_csf;

if total_overlap == 0
    fprintf('  ✅ No tissue overlap detected (mutually exclusive masks)\n');
else
    fprintf('  ⚠ WARNING: Tissue overlap detected:\n');
    fprintf('    WM-GM overlap: %d voxels\n', overlap_wm_gm);
    fprintf('    WM-CSF overlap: %d voxels\n', overlap_wm_csf);
    fprintf('    GM-CSF overlap: %d voxels\n', overlap_gm_csf);
end

% Check coverage (should cover most of brain)
total_tissue = wm_count + gm_count + csf_count;
brain_coverage = 100 * total_tissue / sum(brain_mask_binary(:));

fprintf('  Brain coverage: %.1f%% (%d/%d voxels)\n', ...
        brain_coverage, total_tissue, sum(brain_mask_binary(:)));

if brain_coverage > 90
    fprintf('  ✅ Excellent brain coverage (>90%%)\n');
elseif brain_coverage > 75
    fprintf('  ✓ Good brain coverage (>75%%)\n');
else
    fprintf('  ⚠ WARNING: Low brain coverage (<75%%) - some voxels unclassified\n');
end

% Check tissue proportions (typical brain: ~40-45% WM, ~40-45% GM, ~10-20% CSF)
wm_proportion = 100 * wm_count / sum(brain_mask_binary(:));
gm_proportion = 100 * gm_count / sum(brain_mask_binary(:));
csf_proportion = 100 * csf_count / sum(brain_mask_binary(:));

fprintf('  Tissue proportions:\n');
fprintf('    WM:  %.1f%% (typical: 40-45%%)\n', wm_proportion);
fprintf('    GM:  %.1f%% (typical: 40-45%%)\n', gm_proportion);
fprintf('    CSF: %.1f%% (typical: 10-20%%)\n', csf_proportion);

% Validate proportions are reasonable
if wm_proportion < 20 || wm_proportion > 60
    fprintf('  ⚠ WARNING: WM proportion outside typical range\n');
end
if gm_proportion < 20 || gm_proportion > 60
    fprintf('  ⚠ WARNING: GM proportion outside typical range\n');
end
if csf_proportion > 30
    fprintf('  ⚠ WARNING: CSF proportion unusually high\n');
end

%% Final report
fprintf('\n✅ Tissue segmentation complete:\n');
fprintf('  WM mask:  %s (%.1f MB)\n', wm_mask_file, dir(wm_mask_file).bytes/1024/1024);
fprintf('  GM mask:  %s (%.1f MB)\n', gm_mask_file, dir(gm_mask_file).bytes/1024/1024);
fprintf('  CSF mask: %s (%.1f MB)\n', csf_mask_file, dir(csf_mask_file).bytes/1024/1024);
fprintf('  Ready for Anatomically Constrained Tractography (ACT)\n');

end

% ---------------------------------------------------------------------------
function ok = tissue_from_t1_fast(t1_file, dwi_ref, wm_out, gm_out, csf_out)
% Real WM/GM/CSF via FSL FAST on the anatomical T1, resampled into DWI space via
% the world affines (flirt -usesqform) — NO registration, so a broken/absent
% T1->DWI registration cannot corrupt it (the T1 must be world-aligned to the DWI).
% FAST orders classes by intensity: seg_0=CSF (dark), seg_1=GM, seg_2=WM (bright).
ok = false;
fsldir = getenv('FSLDIR'); if isempty(fsldir), return; end
tmp = tempname;
try
    % brain-extract the T1, then 3-class FAST (full T1 resolution)
    system(sprintf('%s/bin/bet %s %s_brain -f 0.3 -R > /dev/null 2>&1', fsldir, t1_file, tmp));
    bt = [tmp '_brain.nii.gz']; if ~isfile(bt), bt = t1_file; end
    system(sprintf('%s/bin/fast -t 1 -n 3 -g -o %s %s > /dev/null 2>&1', fsldir, tmp, bt));
    if ~isfile([tmp '_seg_2.nii.gz']), return; end
    map = {csf_out, '0'; gm_out, '1'; wm_out, '2'};      % out file <- FAST class idx
    for i = 1:3
        seg = sprintf('%s_seg_%s.nii.gz', tmp, map{i,2});
        % resample T1-space seg -> DWI grid through the world affine, then binarize
        system(sprintf(['%s/bin/flirt -in %s -ref %s -applyxfm -usesqform ' ...
                '-interp nearestneighbour -out %s > /dev/null 2>&1'], fsldir, seg, dwi_ref, map{i,1}));
        system(sprintf('%s/bin/fslmaths %s -mas %s -bin %s > /dev/null 2>&1', ...
                fsldir, map{i,1}, dwi_ref, map{i,1}));
    end
    ok = isfile(wm_out) && isfile(gm_out) && isfile(csf_out);
catch
    ok = false;
end
try, delete([tmp '*']); catch, end
end
