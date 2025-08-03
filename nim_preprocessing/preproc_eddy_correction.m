function eddy_corrected_file = preproc_eddy_correction(dwi_file, brain_mask_file, bvec_file, bval_file, acqp_file, index_file, file_prefix)
% preproc_eddy_correction: Apply eddy current and motion correction using FSL's eddy tool
%
% Arguments:
%   dwi_file - Path to the raw DWI file
%   brain_mask_file - Path to the brain mask
%   bvec_file - Path to the b-vectors file
%   bval_file - Path to the b-values file
%   acqp_file - Path to acquisition parameters file
%   index_file - Path to the index file
%   file_prefix - Prefix for output files
%
% Returns:
%   eddy_corrected_file - Path to the eddy-corrected DWI file
%
% NOTE: This function is currently commented out as eddy correction is not
%       implemented in the current pipeline. Uncomment and modify as needed.

fprintf('Step 3: Eddy current and motion correction...\n');

% COMMENTED OUT - Eddy current correction is not currently implemented
% Uncomment and modify the following code to enable eddy correction:


% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output file path
eddy_corrected_file = [strrep(file_prefix, '_raw', '') '_eddy_corrected.nii.gz'];

% Verify required input files exist
required_files = {dwi_file, brain_mask_file, bvec_file, bval_file, acqp_file, index_file};
file_descriptions = {'DWI', 'brain mask', 'bvec', 'bval', 'acquisition parameters', 'index'};

for i = 1:length(required_files)
    if ~isfile(required_files{i})
        error('Required file for eddy correction not found: %s (%s)', required_files{i}, file_descriptions{i});
    end
end

% Remove file extension from output for eddy (it adds .nii.gz automatically)
[output_dir, output_name, ~] = fileparts(eddy_corrected_file);
if endsWith(output_name, '.nii')
    output_name = output_name(1:end-4);
end
eddy_output_prefix = fullfile(output_dir, output_name);

% Run FSL's eddy tool
% Note: This is a computationally intensive process that can take hours
cmd_eddy = sprintf(['%s/bin/eddy --imain=%s --mask=%s --bvecs=%s --bvals=%s ' ...
                   '--out=%s --acqp=%s --index=%s --repol --cnr_maps --residuals'], ...
                   fsl_path, dwi_file, brain_mask_file, bvec_file, bval_file, ...
                   eddy_output_prefix, acqp_file, index_file);

fprintf('Running eddy correction (this may take several hours)...\n');
fprintf('Command: %s\n', cmd_eddy);

% Show progress indicator
fprintf('Eddy correction in progress');
tic;
[status, cmdout] = system(cmd_eddy);
elapsed_time = toc;

if status ~= 0
    error('Error in eddy correction: %s', cmdout);
end

% Verify the output file was created
eddy_corrected_file = [eddy_output_prefix '.nii.gz'];
if ~isfile(eddy_corrected_file)
    error('Eddy correction failed: output file not found at %s', eddy_corrected_file);
end

% Get file size for reporting
file_info = dir(eddy_corrected_file);
fprintf('✓ Eddy correction completed: %s (%.1f MB)\n', eddy_corrected_file, file_info.bytes/1024/1024);
fprintf('  Processing time: %.1f minutes\n', elapsed_time/60);

% Additional eddy outputs that may be useful:
eddy_outputs = {
    [eddy_output_prefix '.eddy_rotated_bvecs'],  % Corrected b-vectors
    [eddy_output_prefix '.eddy_movement_rms'],    % Movement parameters
    [eddy_output_prefix '.eddy_outlier_map'],     % Outlier map
    [eddy_output_prefix '.eddy_cnr_maps.nii.gz'] % CNR maps
};

fprintf('Additional eddy outputs:\n');
for i = 1:length(eddy_outputs)
    if isfile(eddy_outputs{i})
        file_info = dir(eddy_outputs{i});
        fprintf('  ✓ %s (%.1f MB)\n', eddy_outputs{i}, file_info.bytes/1024/1024);
    end
end

end