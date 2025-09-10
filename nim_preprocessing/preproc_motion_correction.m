function motion_corrected_file = preproc_motion_correction(dwi_file, bvec_file, bval_file, file_prefix)
% preproc_motion_correction: Apply motion correction using FSL's mcflirt
%
% Arguments:
%   dwi_file - Path to the DWI file
%   bvec_file - Path to the b-vectors file
%   bval_file - Path to the b-values file
%   file_prefix - Prefix for output files
%
% Returns:
%   motion_corrected_file - Path to the motion-corrected DWI file

fprintf('Step: Motion correction using FSL mcflirt...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output file path
motion_corrected_file = [strrep(file_prefix, '_raw', '') '_motion_corrected.nii.gz'];

% Verify required input files exist
required_files = {dwi_file, bvec_file, bval_file};
file_descriptions = {'DWI', 'bvec', 'bval'};

for i = 1:length(required_files)
    if ~isfile(required_files{i})
        error('Required file for motion correction not found: %s (%s)', required_files{i}, file_descriptions{i});
    end
end

% Get the output directory and base name
[output_dir, output_name, ~] = fileparts(motion_corrected_file);
if endsWith(output_name, '.nii')
    output_name = output_name(1:end-4);
end
motion_output_prefix = fullfile(output_dir, output_name);

% Extract b=0 volumes for motion reference
fprintf('Extracting b=0 volumes for motion reference...\n');
bvals = load(bval_file);
b0_indices = find(bvals <= 50); % Typically b=0 volumes have values <= 50

if isempty(b0_indices)
    error('No b=0 volumes found in %s', bval_file);
end

% Use the first b=0 volume as reference
ref_vol = b0_indices(1) - 1; % FSL uses 0-based indexing

% Run FSL's mcflirt for motion correction
% mcflirt parameters:
% -in: input 4D image
% -out: output corrected image
% -refvol: reference volume number (0-based)
% -plots: generate plots
% -mats: save transformation matrices
% -rmsrel: save relative RMS displacement
% -rmsabs: save absolute RMS displacement

cmd_mcflirt = sprintf(['%s/bin/mcflirt -in %s -out %s -refvol %d ' ...
                      '-plots -mats -rmsrel -rmsabs'], ...
                      fsl_path, dwi_file, motion_output_prefix, ref_vol);

fprintf('Running motion correction...\n');
fprintf('Command: %s\n', cmd_mcflirt);
fprintf('Reference volume: %d (b=0)\n', ref_vol);

tic;
[status, cmdout] = system(cmd_mcflirt);
elapsed_time = toc;

if status ~= 0
    error('Error in motion correction: %s', cmdout);
end

% Verify the output file was created
if ~isfile(motion_corrected_file)
    error('Motion correction failed: output file not found at %s', motion_corrected_file);
end

% Get file size for reporting
file_info = dir(motion_corrected_file);
fprintf('✓ Motion correction completed: %s (%.1f MB)\n', motion_corrected_file, file_info.bytes/1024/1024);
fprintf('  Processing time: %.1f seconds\n', elapsed_time);

% Report on motion parameters
motion_params_file = [motion_output_prefix '.par'];
motion_rms_rel_file = [motion_output_prefix '_rel_mean.rms'];
motion_rms_abs_file = [motion_output_prefix '_abs_mean.rms'];

if isfile(motion_params_file)
    motion_params = load(motion_params_file);
    max_translation = max(max(abs(motion_params(:, 1:3))));
    max_rotation = max(max(abs(motion_params(:, 4:6)))) * 180/pi; % Convert to degrees
    
    fprintf('Motion analysis:\n');
    fprintf('  Max translation: %.2f mm\n', max_translation);
    fprintf('  Max rotation: %.2f degrees\n', max_rotation);
    
    if max_translation > 3 || max_rotation > 3
        fprintf('  ⚠ WARNING: Significant motion detected (>3mm or >3°)\n');
    else
        fprintf('  ✓ Motion within acceptable limits\n');
    end
end

% Check for relative RMS displacement
if isfile(motion_rms_rel_file)
    rms_rel = load(motion_rms_rel_file);
    mean_rms_rel = mean(rms_rel);
    fprintf('  Mean relative RMS displacement: %.3f mm\n', mean_rms_rel);
    
    if mean_rms_rel > 1
        fprintf('  ⚠ WARNING: High relative motion (>1mm mean RMS)\n');
    end
end

% Rotate b-vectors according to motion correction
fprintf('Rotating b-vectors according to motion correction...\n');
corrected_bvec_file = strrep(bvec_file, '.bvec', '_motion_corrected.bvec');

% Load original b-vectors
bvecs = load(bvec_file);
if size(bvecs, 1) == 3
    bvecs = bvecs'; % Make sure it's Nx3 format
end

% Apply rotations from motion correction
mat_dir = [motion_output_prefix '.mat'];
if isfolder(mat_dir)
    corrected_bvecs = bvecs;
    mat_files = dir(fullfile(mat_dir, 'MAT_*'));
    
    for i = 1:length(mat_files)
        if i <= size(bvecs, 1)
            mat_file = fullfile(mat_files(i).folder, mat_files(i).name);
            if isfile(mat_file)
                % Load transformation matrix as ASCII
transform_mat = load(mat_file, '-ascii');
                if size(transform_mat, 1) == 4 && size(transform_mat, 2) == 4
                    % Extract rotation matrix (top-left 3x3)
                    rotation_mat = transform_mat(1:3, 1:3);
                    % Apply rotation to b-vector
                    if norm(bvecs(i, :)) > 0.1 % Only rotate non-zero b-vectors
                        corrected_bvecs(i, :) = (rotation_mat * bvecs(i, :)')';
                    end
                end
            end
        end
    end
    
    % Save corrected b-vectors
    if size(corrected_bvecs, 1) == 3
        corrected_bvecs = corrected_bvecs'; % Make sure it's 3xN for saving
    else
        corrected_bvecs = corrected_bvecs'; % Transpose to 3xN
    end
    save(corrected_bvec_file, 'corrected_bvecs', '-ascii');
    fprintf('✓ Corrected b-vectors saved to: %s\n', corrected_bvec_file);
else
    fprintf('⚠ Motion transformation matrices not found, copying original b-vectors\n');
    copyfile(bvec_file, corrected_bvec_file);
end

% Additional motion correction outputs
motion_outputs = {
    motion_params_file;              % Motion parameters
    motion_rms_rel_file;             % Relative RMS displacement
    motion_rms_abs_file;             % Absolute RMS displacement
    corrected_bvec_file              % Corrected b-vectors
};

fprintf('Motion correction outputs:\n');
for i = 1:length(motion_outputs)
    if isfile(motion_outputs{i})
        file_info = dir(motion_outputs{i});
        fprintf('  ✓ %s (%.1f KB)\n', motion_outputs{i}, file_info.bytes/1024);
    end
end

end