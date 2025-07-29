function brain_mask_file = preproc_brain_extraction(b0_file, output_dir, brain_mask_file)
% preproc_brain_extraction: Create brain mask using FSL's BET tool
%
% Arguments:
%   b0_file - Path to the b0 volume
%   output_dir - Directory for intermediate files
%   brain_mask_file - Final path for the brain mask
%
% Returns:
%   brain_mask_file - Path to the created brain mask

fprintf('Step 2: Creating brain mask using BET...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define intermediate file paths
bet_output_prefix = fullfile(output_dir, 'nodif_brain');
bet_mask_file = fullfile(output_dir, 'nodif_brain_mask.nii.gz');

% Run FSL's BET tool to create brain mask
cmd_bet = sprintf('%s/bin/bet %s %s -m', fsl_path, b0_file, bet_output_prefix);

fprintf('Running: %s\n', cmd_bet);
[status, cmdout] = system(cmd_bet);

if status ~= 0
    error('Error in bet (brain extraction): %s', cmdout);
end

% Verify the mask file was created
if ~isfile(bet_mask_file)
    error('BET brain extraction failed: mask file not found at %s', bet_mask_file);
end

% Move the mask file to the final location
try
    movefile(bet_mask_file, brain_mask_file);
    fprintf('✓ Brain mask moved to: %s\n', brain_mask_file);
catch ME
    error('Failed to move brain mask from %s to %s: %s', bet_mask_file, brain_mask_file, ME.message);
end

% Verify the final mask file exists
if ~isfile(brain_mask_file)
    error('Brain mask not found at final location: %s', brain_mask_file);
end

% Get file size for reporting
file_info = dir(brain_mask_file);
fprintf('✓ Brain mask created: %s (%.1f MB)\n', brain_mask_file, file_info.bytes/1024/1024);

% Clean up intermediate brain volume (keep only the mask)
bet_brain_file = fullfile(output_dir, 'nodif_brain.nii.gz');
if isfile(bet_brain_file)
    try
        delete(bet_brain_file);
        fprintf('  Cleaned up intermediate brain volume\n');
    catch ME
        warning('Could not remove intermediate brain volume %s: %s', bet_brain_file, ME.message);
    end
end

end