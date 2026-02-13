function [t1_to_dwi_mat] = preproc_t1_dwi_registration(dwi_ref_file, t1_file, t1_brain_file, output_dir, file_prefix)
% preproc_t1_dwi_registration: Register T1 to DWI space using boundary-based registration
%
% Arguments:
%   dwi_ref_file - Reference DWI volume (distortion-corrected B0)
%   t1_file - T1 structural image
%   t1_brain_file - Brain-extracted T1 image
%   output_dir - Directory for output files
%   file_prefix - Prefix for output files
%
% Returns:
%   t1_to_dwi_mat - Path to transformation matrix file

fprintf('T1-DWI registration using FSL epi_reg...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output file paths
% Extract just the filename part from file_prefix to avoid duplicate directory paths
[~, filename_only] = fileparts(file_prefix);
t1_to_dwi_mat = fullfile(output_dir, [filename_only '_T1_to_dwi.mat']);
t1_to_dwi_output = fullfile(output_dir, [filename_only '_T1_to_dwi']);

% Use FSL epi_reg for boundary-based registration
% epi_reg handles EPI distortions and provides robust T1-EPI registration
% --epi: EPI image (DWI reference)
% --t1: T1 image
% --t1brain: Brain-extracted T1
% --out: Output prefix
cmd = sprintf('%s/bin/epi_reg --epi=%s --t1=%s --t1brain=%s --out=%s', ...
    fsl_path, dwi_ref_file, t1_file, t1_brain_file, t1_to_dwi_output);

fprintf('Running: %s\n', cmd);
[status, cmdout] = system(cmd);

if status ~= 0
    error('T1-DWI registration failed: %s', cmdout);
end

% Verify transformation matrix was created
if ~isfile(t1_to_dwi_mat)
    error('T1-DWI registration failed: transformation matrix not found at %s', t1_to_dwi_mat);
end

% Check registration quality by examining the transformation matrix
try
    % Read transformation matrix to validate
    transform_data = dlmread(t1_to_dwi_mat);
    if size(transform_data, 1) ~= 4 || size(transform_data, 2) ~= 4
        warning('T1-DWI transformation matrix has unexpected dimensions');
    end
    
    % Check for reasonable translation values (should be < 50mm typically)
    translation_x = abs(transform_data(1, 4));
    translation_y = abs(transform_data(2, 4));
    translation_z = abs(transform_data(3, 4));
    
    if translation_x > 100 || translation_y > 100 || translation_z > 100
        warning('T1-DWI registration may have failed: large translations detected (%.1f, %.1f, %.1f)', ...
            translation_x, translation_y, translation_z);
    end
    
catch ME
    warning(ME.identifier, 'Could not validate T1-DWI transformation matrix: %s', ME.message);
end

fprintf('✓ T1-DWI registration completed\n');
fprintf('  Transformation matrix: %s\n', t1_to_dwi_mat);

end