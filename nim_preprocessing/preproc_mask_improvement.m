function improved_mask_file = preproc_mask_improvement(brain_mask_file, fa_file, file_prefix)
% preproc_mask_improvement: Improve brain mask quality using FA data
%
% Arguments:
%   brain_mask_file - Path to the initial brain mask
%   fa_file - Path to the FA map for mask refinement
%   file_prefix - Prefix for output files
%
% Returns:
%   improved_mask_file - Path to the improved brain mask

fprintf('Step: Brain mask improvement and validation...\n');

% Define output file path
improved_mask_file = [strrep(file_prefix, '_raw', '') '_mask_improved.nii.gz'];

% Verify input files exist
if ~isfile(brain_mask_file)
    error('Brain mask file not found: %s', brain_mask_file);
end

if ~isfile(fa_file)
    error('FA file not found: %s', fa_file);
end

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

fprintf('Loading brain mask and FA data...\n');

% Load original mask and FA data for analysis
V_mask = spm_vol(brain_mask_file);
mask_data = spm_read_vols(V_mask);

V_fa = spm_vol(fa_file);
fa_data = spm_read_vols(V_fa);

% Analyze original mask quality
mask_volume = sum(mask_data(:) > 0.5);
total_volume = numel(mask_data);
mask_percentage = 100 * mask_volume / total_volume;

fprintf('Original mask statistics:\n');
fprintf('  Mask volume: %d voxels (%.1f%% of image)\n', mask_volume, mask_percentage);

% Check FA values outside mask (should be close to 0)
fa_outside_mask = fa_data(mask_data < 0.5);
mean_fa_outside = mean(fa_outside_mask(~isnan(fa_outside_mask)));

fprintf('  Mean FA outside mask: %.3f (should be ~0)\n', mean_fa_outside);

% Check for obvious mask issues
if mask_percentage > 60
    fprintf('  ⚠ WARNING: Mask covers >60%% of image - may include non-brain tissue\n');
elseif mask_percentage < 20
    fprintf('  ⚠ WARNING: Mask covers <20%% of image - may be too restrictive\n');
else
    fprintf('  ✓ Mask coverage looks reasonable\n');
end

if mean_fa_outside > 0.1
    fprintf('  ⚠ WARNING: High FA outside mask - mask may be too restrictive\n');
else
    fprintf('  ✓ Low FA outside mask - good mask quality\n');
end

% Strategy 1: Use FSL's bet2 with more aggressive settings if original mask is poor
if mask_percentage > 50 || mean_fa_outside > 0.1
    fprintf('Applying improved brain extraction...\n');
    
    % Extract b0 for better brain extraction
    temp_b0 = [strrep(file_prefix, '_raw', '') '_temp_b0_for_mask.nii.gz'];
    
    % Get the b0 volume from the FA's header info to find the original DWI
    dwi_file = strrep(fa_file, '_FA', '');
    if ~isfile(dwi_file)
        dwi_file = [strrep(file_prefix, '_raw', '') '.nii.gz'];
    end
    
    if isfile(dwi_file)
        % Extract first volume (assumed to be b0 or close to it)
        cmd_extract_b0 = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, temp_b0);
        [status, ~] = system(cmd_extract_b0);
        
        if status == 0
            % Apply bet with more conservative settings
            temp_bet_prefix = [strrep(file_prefix, '_raw', '') '_temp_bet'];
            cmd_bet = sprintf('%s/bin/bet %s %s -m -f 0.3 -R', fsl_path, temp_b0, temp_bet_prefix);
            
            fprintf('Improving brain extraction: %s\n', cmd_bet);
            [status, cmdout] = system(cmd_bet);
            
            if status == 0
                new_mask = [temp_bet_prefix '_mask.nii.gz'];
                if isfile(new_mask)
                    % Test the new mask
                    V_new = spm_vol(new_mask);
                    new_mask_data = spm_read_vols(V_new);
                    new_mask_volume = sum(new_mask_data(:) > 0.5);
                    new_mask_percentage = 100 * new_mask_volume / total_volume;
                    
                    fa_outside_new = fa_data(new_mask_data < 0.5);
                    mean_fa_outside_new = mean(fa_outside_new(~isnan(fa_outside_new)));
                    
                    fprintf('New mask statistics:\n');
                    fprintf('  New mask volume: %d voxels (%.1f%% of image)\n', new_mask_volume, new_mask_percentage);
                    fprintf('  Mean FA outside new mask: %.3f\n', mean_fa_outside_new);
                    
                    % Use new mask if it's better
                    if mean_fa_outside_new < mean_fa_outside && new_mask_percentage > 15 && new_mask_percentage < 55
                        fprintf('✓ Using improved mask\n');
                        copyfile(new_mask, improved_mask_file);
                        mask_data = new_mask_data;
                    else
                        fprintf('Improved mask not better, keeping original\n');
                        copyfile(brain_mask_file, improved_mask_file);
                    end
                else
                    fprintf('New mask file not created, keeping original\n');
                    copyfile(brain_mask_file, improved_mask_file);
                end
            else
                fprintf('Brain extraction failed, keeping original mask\n');
                copyfile(brain_mask_file, improved_mask_file);
            end
            
            % Clean up temporary files
            if isfile(temp_b0)
                delete(temp_b0);
            end
            temp_files = {[temp_bet_prefix '.nii.gz'], [temp_bet_prefix '_mask.nii.gz']};
            for i = 1:length(temp_files)
                if isfile(temp_files{i})
                    delete(temp_files{i});
                end
            end
        else
            fprintf('Failed to extract b0, keeping original mask\n');
            copyfile(brain_mask_file, improved_mask_file);
        end
    else
        fprintf('DWI file not found for mask improvement, keeping original\n');
        copyfile(brain_mask_file, improved_mask_file);
    end
else
    fprintf('Original mask quality is good, keeping it\n');
    copyfile(brain_mask_file, improved_mask_file);
end

% Strategy 2: Fill holes and smooth the mask
fprintf('Applying morphological operations to improve mask...\n');

% Fill holes and smooth using FSL
temp_filled = [strrep(file_prefix, '_raw', '') '_temp_filled.nii.gz'];

% Fill holes and apply morphological operations
cmd_morph = sprintf(['%s/bin/fslmaths %s -fillh -s 0.5 -thr 0.5 -bin %s'], ...
                   fsl_path, improved_mask_file, temp_filled);

[status, ~] = system(cmd_morph);

if status == 0 && isfile(temp_filled)
    % Verify the morphologically processed mask
    V_filled = spm_vol(temp_filled);
    filled_data = spm_read_vols(V_filled);
    
    filled_volume = sum(filled_data(:) > 0.5);
    filled_percentage = 100 * filled_volume / total_volume;
    
    % Use filled mask if it's reasonable
    if filled_percentage > 15 && filled_percentage < 60
        fprintf('✓ Using morphologically improved mask\n');
        copyfile(temp_filled, improved_mask_file);
    else
        fprintf('Morphological operations made mask worse, reverting\n');
    end
    
    % Clean up
    if isfile(temp_filled)
        delete(temp_filled);
    end
end

% Final validation
if isfile(improved_mask_file)
    file_info = dir(improved_mask_file);
    fprintf('✓ Improved brain mask saved: %s (%.1f MB)\n', improved_mask_file, file_info.bytes/1024/1024);
    
    % Final quality check
    V_final = spm_vol(improved_mask_file);
    final_mask_data = spm_read_vols(V_final);
    final_volume = sum(final_mask_data(:) > 0.5);
    final_percentage = 100 * final_volume / total_volume;
    
    fa_outside_final = fa_data(final_mask_data < 0.5);
    mean_fa_outside_final = mean(fa_outside_final(~isnan(fa_outside_final)));
    
    fprintf('Final mask quality:\n');
    fprintf('  Volume: %d voxels (%.1f%% of image)\n', final_volume, final_percentage);
    fprintf('  Mean FA outside mask: %.3f\n', mean_fa_outside_final);
    
    if final_percentage > 15 && final_percentage < 55 && mean_fa_outside_final < 0.15
        fprintf('  ✅ Mask quality assessment: GOOD\n');
    elseif final_percentage > 10 && final_percentage < 65 && mean_fa_outside_final < 0.25
        fprintf('  ⚠ Mask quality assessment: ACCEPTABLE\n');
    else
        fprintf('  ❌ Mask quality assessment: POOR - manual review recommended\n');
    end
else
    error('Failed to create improved mask file');
end

end