function [t1_brain_file, t1_brain_mask_file] = preproc_t1_brain_extraction(t1_file, output_dir)
% preproc_t1_brain_extraction: Extract brain from T1 structural image
%
% Arguments:
%   t1_file - Path to T1 structural NIfTI file
%   output_dir - Directory for output files
%
% Returns:
%   t1_brain_file - Path to brain-extracted T1 image
%   t1_brain_mask_file - Path to brain mask

fprintf('T1 brain extraction using FSL BET...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Get file prefix for naming outputs
[~, base_name, ~] = fileparts(t1_file);
base_name = strrep(base_name, '.nii', ''); % Remove .nii if present

% Define output file paths
t1_brain_file = fullfile(output_dir, [base_name '_brain.nii.gz']);
t1_brain_mask_file = fullfile(output_dir, [base_name '_brain_mask.nii.gz']);

% Use FSL BET with robust parameters for T1 brain extraction
% -R: Use robust brain centre estimation
% -m: Generate brain mask
% -f 0.4: More conservative threshold for T1 (default 0.5)
cmd = sprintf('%s/bin/bet %s %s -R -m -f 0.4', ...
    fsl_path, t1_file, strrep(t1_brain_file, '_brain.nii.gz', '_brain'));

fprintf('Running: %s\n', cmd);
[status, cmdout] = system(cmd);

if status ~= 0
    error('T1 brain extraction failed: %s', cmdout);
end

% Verify output files were created
if ~isfile(t1_brain_file)
    error('T1 brain extraction failed: output file not found at %s', t1_brain_file);
end

if ~isfile(t1_brain_mask_file)
    error('T1 brain extraction failed: mask file not found at %s', t1_brain_mask_file);
end

fprintf('✓ T1 brain extraction completed\n');
fprintf('  Brain image: %s\n', t1_brain_file);
fprintf('  Brain mask: %s\n', t1_brain_mask_file);

end