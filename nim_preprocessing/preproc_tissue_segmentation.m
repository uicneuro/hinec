function [wm_mask_file, gm_mask_file, csf_mask_file] = preproc_tissue_segmentation(fa_data, brain_mask_file, file_prefix)
% preproc_tissue_segmentation: Generate tissue-specific masks using FA-based segmentation
%
% This function creates three tissue masks for Anatomically Constrained Tractography (ACT):
%   - White Matter (WM): FA > 0.2, eroded to remove boundary voxels
%   - Gray Matter (GM): 0.05 < FA < 0.2, intermediate anisotropy
%   - CSF: FA < 0.05, low anisotropy
%
% Arguments:
%   fa_data - FA data array (from nim.FA) for tissue classification
%   brain_mask_file - Path to the brain mask (tissue masks will be constrained to brain)
%   file_prefix - Prefix for output files
%
% Returns:
%   wm_mask_file - Path to white matter mask
%   gm_mask_file - Path to gray matter mask
%   csf_mask_file - Path to CSF mask
%
% Reference: PIPELINE.md Step 7 - White Matter Segmentation (extended for ACT)

fprintf('Step: FA-based tissue segmentation for ACT...\n');

% Define output file paths
wm_mask_file = [strrep(file_prefix, '_raw', '') '_WM_mask.nii.gz'];
gm_mask_file = [strrep(file_prefix, '_raw', '') '_GM_mask.nii.gz'];
csf_mask_file = [strrep(file_prefix, '_raw', '') '_CSF_mask.nii.gz'];

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

fprintf('Generating tissue masks using FA thresholds...\n');
fprintf('  WM threshold: FA > 0.2\n');
fprintf('  GM threshold: 0.05 < FA <= 0.2\n');
fprintf('  CSF threshold: FA <= 0.05\n');

%% Generate White Matter Mask
% WM: FA > 0.2, constrained to brain, eroded to remove boundary voxels
wm_preliminary = (fa_data > 0.2) & brain_mask_binary;

% Count preliminary WM voxels
wm_prelim_count = sum(wm_preliminary(:));
fprintf('  Preliminary WM voxels: %d (%.1f%% of brain)\n', ...
        wm_prelim_count, 100*wm_prelim_count/sum(brain_mask_binary(:)));

% Save preliminary WM mask for erosion
temp_wm_prelim = [strrep(file_prefix, '_raw', '') '_temp_WM_prelim.nii.gz'];
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
gm_mask = (fa_data > 0.05) & (fa_data <= 0.2) & brain_mask_binary;

gm_count = sum(gm_mask(:));
fprintf('  ✓ GM voxels: %d (%.1f%% of brain)\n', ...
        gm_count, 100*gm_count/sum(brain_mask_binary(:)));

% Save GM mask
V_out.fname = gm_mask_file;
spm_write_vol(V_out, double(gm_mask));

%% Generate CSF Mask
% CSF: FA <= 0.05, constrained to brain
csf_mask = (fa_data <= 0.05) & brain_mask_binary;

csf_count = sum(csf_mask(:));
fprintf('  ✓ CSF voxels: %d (%.1f%% of brain)\n', ...
        csf_count, 100*csf_count/sum(brain_mask_binary(:)));

% Save CSF mask
V_out.fname = csf_mask_file;
spm_write_vol(V_out, double(csf_mask));

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
