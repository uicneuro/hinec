function registration_data = extract_reference_volumes(registration_data, options)
% extract_reference_volumes: Extract reference volumes for registration
%
% This function extracts b0 volume from DWI data for registration purposes

fprintf('Extracting reference volumes...\n');

dwi_file = registration_data.input.dwi_file;
output_dir = registration_data.output_dir;

% Extract b0 volume from DWI
fprintf('  Extracting b0 volume from DWI...\n');

[dwi_dir, dwi_name, dwi_ext] = fileparts(dwi_file);
if strcmp(dwi_ext, '.gz')
    [~, dwi_name, dwi_ext2] = fileparts(dwi_name);
    dwi_ext = [dwi_ext2 dwi_ext];
end

b0_file = fullfile(output_dir, [dwi_name '_b0' dwi_ext]);

% Use FSL to extract first volume (assumed to be b0 or low b-value)
fsl_path = getenv('FSLDIR');
if ~isempty(fsl_path)
    cmd_extract = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, b0_file);
    [status, cmdout] = system(cmd_extract);

    if status ~= 0
        error('Failed to extract b0 volume: %s', cmdout);
    end

    fprintf('    ✓ B0 volume extracted: %s\n', b0_file);
else
    % Fallback: copy entire DWI file (not ideal but works)
    warning('FSL not found, using entire DWI file as reference');
    copyfile(dwi_file, b0_file);
end

% Store reference volume information
registration_data.reference_volumes = struct();
registration_data.reference_volumes.b0_file = b0_file;

fprintf('  ✓ Reference volumes extracted\n');

end
