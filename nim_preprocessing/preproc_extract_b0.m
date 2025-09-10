function b0_file = preproc_extract_b0(dwi_file, output_dir)
% preproc_extract_b0: Extract b0 volume from DWI data
%
% Arguments:
%   dwi_file - Path to the raw DWI file
%   output_dir - Directory to save the b0 volume
%
% Returns:
%   b0_file - Path to the extracted b0 volume

fprintf('Step 1: Extracting b0 volume...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output file path
b0_file = fullfile(output_dir, 'b0.nii.gz');

% Extract b0 image (first volume) from DWI data using fslroi
cmd_extract_b0 = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, b0_file);

fprintf('Running: %s\n', cmd_extract_b0);
[status, cmdout] = system(cmd_extract_b0);

if status ~= 0
    error('Error in fslroi (b0 extraction): %s', cmdout);
end

% Verify the output file was created
if ~isfile(b0_file)
    error('B0 extraction failed: output file not found at %s', b0_file);
end

% Get file size for reporting
file_info = dir(b0_file);
fprintf('✓ B0 volume extracted: %s (%.1f MB)\n', b0_file, file_info.bytes/1024/1024);

end