function dwi_ref_file = preproc_create_dwi_reference(dwi_file, brain_mask_file, output_dir, file_prefix)
% preproc_create_dwi_reference: Create distortion-corrected DWI reference volume
%
% Arguments:
%   dwi_file - Processed DWI file (after eddy correction)
%   brain_mask_file - Brain mask file
%   output_dir - Directory for output files
%   file_prefix - Prefix for output files
%
% Returns:
%   dwi_ref_file - Path to DWI reference volume

fprintf('Creating DWI reference volume...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output file paths
% Extract just the filename part from file_prefix to avoid duplicate directory paths
[~, filename_only] = fileparts(file_prefix);
b0_temp_file = fullfile(output_dir, [filename_only '_b0_temp.nii.gz']);
dwi_ref_file = fullfile(output_dir, [filename_only '_dwi_ref.nii.gz']);

%% Step 1: Extract B0 volume from processed DWI
fprintf('Extracting B0 volume from processed DWI...\n');
cmd_extract = sprintf('%s/bin/fslroi %s %s 0 1', ...
    fsl_path, dwi_file, b0_temp_file);

fprintf('Running: %s\n', cmd_extract);
[status, cmdout] = system(cmd_extract);

if status ~= 0
    error('B0 extraction failed: %s', cmdout);
end

if ~isfile(b0_temp_file)
    error('B0 extraction failed: output file not found at %s', b0_temp_file);
end

%% Step 2: Apply brain mask to create clean reference
fprintf('Applying brain mask to create DWI reference...\n');
cmd_mask = sprintf('%s/bin/fslmaths %s -mas %s %s', ...
    fsl_path, b0_temp_file, brain_mask_file, dwi_ref_file);

fprintf('Running: %s\n', cmd_mask);
[status, cmdout] = system(cmd_mask);

if status ~= 0
    error('Brain masking failed: %s', cmdout);
end

if ~isfile(dwi_ref_file)
    error('DWI reference creation failed: output file not found at %s', dwi_ref_file);
end

%% Cleanup temporary file
if isfile(b0_temp_file)
    delete(b0_temp_file);
end

%% Validate DWI reference
try
    % Check that the reference has reasonable intensity values
    cmd_stats = sprintf('%s/bin/fslstats %s -M -S', fsl_path, dwi_ref_file);
    [status, stats_output] = system(cmd_stats);
    
    if status == 0
        stats = str2num(stats_output);
        if length(stats) >= 2
            mean_intensity = stats(1);
            std_intensity = stats(2);
            
            fprintf('DWI reference statistics: mean=%.1f, std=%.1f\n', mean_intensity, std_intensity);
            
            if mean_intensity < 10 || std_intensity < 1
                warning('DWI reference may have low intensity values');
            end
        end
    end
catch ME
    warning(ME.identifier, 'Could not validate DWI reference statistics: %s', ME.message);
end

fprintf('✓ DWI reference volume created: %s\n', dwi_ref_file);

end