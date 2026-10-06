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

% Pick the b0 from the b-values (same rule as nim_read: b < 5), not volume 0:
% many acquisitions do not start with a b0 (MASiVar: the only b0 is volume 32).
b0_index = 0;
bval_file = fullfile(dwi_dir, [dwi_name '.bval']);
if isfile(bval_file)
    bvals = load(bval_file);
    first_b0 = find(bvals(:) < 5, 1);
    if isempty(first_b0)
        error('extract_reference_volumes:noB0', 'No b0 volume (b < 5) in %s', bval_file);
    end
    b0_index = first_b0 - 1;                    % fslroi is 0-based
else
    warning('No b-value file next to %s; assuming volume 0 is a b0', dwi_file);
end
fprintf('    b0 volume index (0-based): %d\n', b0_index);

fsl_path = getenv('FSLDIR');
if ~isempty(fsl_path)
    cmd_extract = sprintf('%s/bin/fslroi %s %s %d 1', fsl_path, dwi_file, b0_file, b0_index);
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
