function brain_mask_file = preproc_brain_extraction(dwi_or_b0_file, output_dir, brain_mask_file)
% preproc_brain_extraction: Create brain mask using FSL's BET
%
% Works with both:
%   - Extracted b0 volume (from preprocessing pipeline)
%   - Full DWI dataset (extracts b0 automatically for preprocessed data)
%
% Arguments:
%   dwi_or_b0_file - Path to b0 volume OR full DWI file
%   output_dir - Directory for intermediate files
%   brain_mask_file - Final path for the brain mask
%
% Returns:
%   brain_mask_file - Path to the created brain mask
%
% Examples:
%   % From preprocessing pipeline (b0 already extracted)
%   mask = preproc_brain_extraction('nodif.nii.gz', 'output/', 'mask.nii.gz');
%
%   % From preprocessed DWI (will extract b0)
%   mask = preproc_brain_extraction('sub-MR256.nii.gz', './', 'sub-MR256_M.nii.gz');

fprintf('Creating brain mask using BET...\n');

%% Validate environment
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR not set. Please ensure FSL is installed and configured.');
end

%% Determine if input is b0 or full DWI
[~, fname, ext] = fileparts(dwi_or_b0_file);
if strcmp(ext, '.gz')
    [~, fname] = fileparts(fname);
end

% Check if we need to extract b0 first
needs_b0_extraction = ~contains(fname, 'nodif') && ~contains(fname, 'b0');

if needs_b0_extraction
    fprintf('  Extracting b0 volume from DWI...\n');
    b0_file = extract_b0_from_dwi(dwi_or_b0_file, output_dir, fsl_path);
else
    b0_file = dwi_or_b0_file;
    fprintf('  Using provided b0 volume\n');
end

% Define intermediate file paths
bet_output_prefix = fullfile(output_dir, 'nodif_brain');
bet_mask_file = fullfile(output_dir, 'nodif_brain_mask.nii.gz');

%% Run FSL BET brain extraction
cmd_bet = sprintf('%s/bin/bet %s %s -m -f 0.3', fsl_path, b0_file, bet_output_prefix);

fprintf('  Running BET...\n');
[status, cmdout] = system(cmd_bet);

if status ~= 0
    error('BET failed: %s', cmdout);
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

%% Cleanup intermediate files
cleanup_files = {[bet_output_prefix '.nii.gz']};
if needs_b0_extraction
    cleanup_files{end+1} = b0_file;
end

for i = 1:length(cleanup_files)
    if isfile(cleanup_files{i})
        delete(cleanup_files{i});
    end
end

fprintf('✓ Brain mask created: %s\n', brain_mask_file);
end

%% Helper Functions

function b0_file = extract_b0_from_dwi(dwi_file, output_dir, fsl_path)
% Extract mean b0 volume from DWI dataset
    b0_file = fullfile(output_dir, 'temp_nodif.nii.gz');

    % Try to find bval file
    bval_file = strrep(dwi_file, '.nii.gz', '.bval');
    if ~isfile(bval_file)
        bval_file = '';
    end

    if ~isempty(bval_file)
        % Extract and average b0 volumes
        bval_data = load(bval_file);
        b0_indices = find(bval_data < 50);
        if isempty(b0_indices)
            b0_indices = 1;
        end

        fprintf('    Found %d b0 volume(s)\n', length(b0_indices));

        % Extract individual b0s
        temp_b0s = cell(length(b0_indices), 1);
        for i = 1:length(b0_indices)
            idx = b0_indices(i) - 1; % FSL: 0-based
            temp_b0s{i} = fullfile(output_dir, sprintf('temp_b0_%03d.nii.gz', i));
            cmd = sprintf('%s/bin/fslroi %s %s %d 1', fsl_path, dwi_file, temp_b0s{i}, idx);
            [status, ~] = system(cmd);
            if status ~= 0
                error('Failed to extract b0 volume %d', i);
            end
        end

        % Average b0 volumes
        if length(temp_b0s) > 1
            cmd = sprintf('%s/bin/fslmaths %s', fsl_path, temp_b0s{1});
            for i = 2:length(temp_b0s)
                cmd = sprintf('%s -add %s', cmd, temp_b0s{i});
            end
            cmd = sprintf('%s -div %d %s', cmd, length(temp_b0s), b0_file);
        else
            cmd = sprintf('%s/bin/fslmaths %s %s', fsl_path, temp_b0s{1}, b0_file);
        end

        [status, cmdout] = system(cmd);
        if status ~= 0
            error('Failed to average b0 volumes: %s', cmdout);
        end

        % Cleanup temp b0 files
        cellfun(@delete, temp_b0s(isfile(string(temp_b0s))));
    else
        % No bval - extract first volume
        warning('No .bval file found. Using first volume as b0.');
        cmd = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, b0_file);
        [status, cmdout] = system(cmd);
        if status ~= 0
            error('Failed to extract first volume: %s', cmdout);
        end
    end
end