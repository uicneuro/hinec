function [mni_to_t1_warp] = preproc_t1_mni_registration(t1_file, t1_brain_file, output_dir, file_prefix)
% preproc_t1_mni_registration: Register T1 to MNI space and create inverse warp
%
% Arguments:
%   t1_file - T1 structural image
%   t1_brain_file - Brain-extracted T1 image
%   output_dir - Directory for output files
%   file_prefix - Prefix for output files
%
% Returns:
%   mni_to_t1_warp - Path to inverse warp field (MNI to T1)

fprintf('T1-MNI registration using FSL FLIRT+FNIRT...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define MNI template
mni_template = fullfile(fsl_path, 'data/standard/MNI152_T1_1mm_brain.nii.gz');
if ~isfile(mni_template)
    error('MNI template not found: %s', mni_template);
end

% Define output file paths
% Extract just the filename part from file_prefix to avoid duplicate directory paths
[~, filename_only] = fileparts(file_prefix);
t1_to_mni_lin_mat = fullfile(output_dir, [filename_only '_T1_to_MNI_lin.mat']);
t1_to_mni_lin_img = fullfile(output_dir, [filename_only '_T1_to_MNI_lin.nii.gz']);
t1_to_mni_warp = fullfile(output_dir, [filename_only '_T1_to_MNI_warp']);
mni_to_t1_warp = fullfile(output_dir, [filename_only '_MNI_to_T1_warp.nii.gz']);

%% Step 1: Linear pre-alignment using FLIRT
fprintf('Step 1: Linear T1 to MNI registration...\n');
cmd_linear = sprintf('%s/bin/flirt -in %s -ref %s -omat %s -out %s', ...
    fsl_path, t1_brain_file, mni_template, t1_to_mni_lin_mat, t1_to_mni_lin_img);

fprintf('Running: %s\n', cmd_linear);
[status, cmdout] = system(cmd_linear);

if status ~= 0
    error('Linear T1-MNI registration failed: %s', cmdout);
end

if ~isfile(t1_to_mni_lin_mat)
    error('Linear registration failed: matrix not found at %s', t1_to_mni_lin_mat);
end

%% Step 2: Nonlinear registration using FNIRT
fprintf('Step 2: Nonlinear T1 to MNI registration...\n');
cmd_nonlinear = sprintf('%s/bin/fnirt --in=%s --aff=%s --cout=%s --config=T1_2_MNI152_2mm', ...
    fsl_path, t1_file, t1_to_mni_lin_mat, t1_to_mni_warp);

fprintf('Running: %s\n', cmd_nonlinear);
[status, cmdout] = system(cmd_nonlinear);

if status ~= 0
    error('Nonlinear T1-MNI registration failed: %s', cmdout);
end

% Check if warp field was created (FNIRT creates .nii.gz automatically)
t1_to_mni_warp_file = [t1_to_mni_warp '.nii.gz'];
if ~isfile(t1_to_mni_warp_file)
    error('Nonlinear registration failed: warp field not found at %s', t1_to_mni_warp_file);
end

%% Step 3: Invert warp field to get MNI to T1 transformation
fprintf('Step 3: Inverting warp field...\n');
cmd_invert = sprintf('%s/bin/invwarp -w %s -r %s -o %s', ...
    fsl_path, t1_to_mni_warp_file, t1_file, mni_to_t1_warp);

fprintf('Running: %s\n', cmd_invert);
[status, cmdout] = system(cmd_invert);

if status ~= 0
    error('Warp field inversion failed: %s', cmdout);
end

if ~isfile(mni_to_t1_warp)
    error('Warp field inversion failed: inverse warp not found at %s', mni_to_t1_warp);
end

%% Validate registration quality
try
    % Check if linear transformation is reasonable
    transform_data = dlmread(t1_to_mni_lin_mat);
    if size(transform_data, 1) ~= 4 || size(transform_data, 2) ~= 4
        warning('T1-MNI linear transformation matrix has unexpected dimensions');
    end

    % Check scaling factors (should be close to 1 for good registration)
    scale_x = norm(transform_data(1, 1:3));
    scale_y = norm(transform_data(2, 1:3));
    scale_z = norm(transform_data(3, 1:3));

    if scale_x < 0.8 || scale_x > 1.2 || scale_y < 0.8 || scale_y > 1.2 || scale_z < 0.8 || scale_z > 1.2
        warning('T1-MNI registration may have unusual scaling: %.3f, %.3f, %.3f', ...
            scale_x, scale_y, scale_z);
    end

catch ME
    warning(ME.identifier, 'Could not validate T1-MNI transformation: %s', ME.message);
end

fprintf('✓ T1-MNI registration completed\n');
fprintf('  Linear matrix: %s\n', t1_to_mni_lin_mat);
fprintf('  Nonlinear warp: %s\n', t1_to_mni_warp_file);
fprintf('  Inverse warp: %s\n', mni_to_t1_warp);

end