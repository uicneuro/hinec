function registration_data = register_dti_to_t1(registration_data, options)
% register_dti_to_t1: Register DTI data to T1 anatomical space
%
% This function performs robust registration between DTI and T1 images using
% the b0 volume as the DTI representative image.

fprintf('Registering DTI to T1 anatomical space...\n');

% Define output paths
dti_to_t1_matrix = [registration_data.output_prefix '_dti_to_t1.mat'];
dti_to_t1_transform = [registration_data.output_prefix '_dti_to_t1_transform.txt'];
registered_b0 = [registration_data.output_prefix '_b0_in_t1.nii.gz'];
t1_in_dti = [registration_data.output_prefix '_t1_in_dti.nii.gz'];

% Check if registration already exists and not forcing recompute
if ~options.force_recompute && isfile(dti_to_t1_matrix)
    fprintf('  DTI->T1 registration already exists, loading...\n');
    load(dti_to_t1_matrix, 'dti_to_t1_transform_matrix');
    registration_data.transforms.dti_to_t1_matrix = dti_to_t1_transform_matrix;
    registration_data.transforms.dti_to_t1_file = dti_to_t1_matrix;
    return;
end

% Perform registration based on method
switch lower(options.registration_method)
    case 'fsl'
        registration_data = register_dti_to_t1_fsl(registration_data, options, ...
            dti_to_t1_matrix, dti_to_t1_transform, registered_b0, t1_in_dti);
    case 'spm'
        registration_data = register_dti_to_t1_spm(registration_data, options, ...
            dti_to_t1_matrix, registered_b0, t1_in_dti);
    otherwise
        error('Unknown registration method: %s', options.registration_method);
end

fprintf('  ✓ DTI to T1 registration complete\n');

end

function registration_data = register_dti_to_t1_fsl(registration_data, options, ...
    dti_to_t1_matrix, dti_to_t1_transform, registered_b0, t1_in_dti)
% Register DTI to T1 using FSL tools

fsl_path = getenv('FSLDIR');
b0_file = registration_data.reference_volumes.b0_file;
t1_file = registration_data.input.t1_file;

fprintf('  Using FSL FLIRT for DTI->T1 registration...\n');

% Step 1: Initial linear registration (rigid, 6 DOF)
fprintf('    Running rigid registration (6 DOF)...\n');

% Create brain-extracted versions for better registration
b0_brain = [registration_data.output_prefix '_b0_brain.nii.gz'];
t1_brain = [registration_data.output_prefix '_t1_brain.nii.gz'];

% Brain extract b0
cmd_bet_b0 = sprintf('%s/bin/bet %s %s -f 0.3 -R', fsl_path, b0_file, ...
    strrep(b0_brain, '.nii.gz', ''));
[status, cmdout] = system(cmd_bet_b0);
if status ~= 0
    warning('B0 brain extraction failed: %s', cmdout);
    b0_brain = b0_file; % Use original if brain extraction fails
end

% Brain extract T1
cmd_bet_t1 = sprintf('%s/bin/bet %s %s -f 0.5 -B', fsl_path, t1_file, ...
    strrep(t1_brain, '.nii.gz', ''));
[status, cmdout] = system(cmd_bet_t1);
if status ~= 0
    warning('T1 brain extraction failed: %s', cmdout);
    t1_brain = t1_file; % Use original if brain extraction fails
end

% Run FLIRT registration
% Use correlation ratio as cost function (good for DTI->T1)
cmd_flirt = sprintf(['%s/bin/flirt -in %s -ref %s -out %s ' ...
                    '-omat %s -dof %d -cost corratio -searchrx -90 90 ' ...
                    '-searchry -90 90 -searchrz -90 90 -interp trilinear'], ...
                    fsl_path, b0_brain, t1_brain, registered_b0, ...
                    dti_to_t1_transform, options.dti_reg_dof);

fprintf('    Command: %s\n', cmd_flirt);
[status, cmdout] = system(cmd_flirt);

if status ~= 0
    error('FLIRT registration failed: %s', cmdout);
end

% Step 2: Register T1 to DTI space (inverse transform)
fprintf('    Computing inverse transform (T1->DTI)...\n');
t1_to_dti_transform = [registration_data.output_prefix '_t1_to_dti_transform.txt'];

cmd_convert = sprintf('%s/bin/convert_xfm -omat %s -inverse %s', ...
    fsl_path, t1_to_dti_transform, dti_to_t1_transform);
[status, cmdout] = system(cmd_convert);

if status ~= 0
    error('Transform inversion failed: %s', cmdout);
end

% Apply inverse transform to get T1 in DTI space
cmd_apply_inverse = sprintf('%s/bin/flirt -in %s -ref %s -out %s -init %s -applyxfm', ...
    fsl_path, t1_brain, b0_file, t1_in_dti, t1_to_dti_transform);
[status, cmdout] = system(cmd_apply_inverse);

if status ~= 0
    warning('Inverse transform application failed: %s', cmdout);
end

% Step 3: Load and store transformation matrix
if isfile(dti_to_t1_transform)
    % Read FSL transformation matrix
    fsl_matrix = load(dti_to_t1_transform, '-ascii');
    
    % Convert FSL matrix to standard 4x4 format and save
    dti_to_t1_transform_matrix = fsl_matrix;
    save(dti_to_t1_matrix, 'dti_to_t1_transform_matrix');
    
    % Store in registration data
    registration_data.transforms.dti_to_t1_matrix = dti_to_t1_transform_matrix;
    registration_data.transforms.dti_to_t1_file = dti_to_t1_matrix;
    registration_data.transforms.dti_to_t1_fsl_file = dti_to_t1_transform;
    registration_data.transforms.t1_to_dti_fsl_file = t1_to_dti_transform;
    
    fprintf('    ✓ Registration matrix saved: %s\n', dti_to_t1_matrix);
else
    error('Registration failed: transform file not created');
end

% Store registered image paths
registration_data.registered_images.b0_in_t1 = registered_b0;
registration_data.registered_images.t1_in_dti = t1_in_dti;

% Clean up temporary brain-extracted files
if isfile(b0_brain) && ~strcmp(b0_brain, b0_file)
    delete(b0_brain);
end
if isfile(t1_brain) && ~strcmp(t1_brain, t1_file)
    delete(t1_brain);
end

fprintf('    ✓ DTI->T1 registration successful\n');

end

function registration_data = register_dti_to_t1_spm(registration_data, options, ...
    dti_to_t1_matrix, registered_b0, t1_in_dti)
% Register DTI to T1 using SPM tools

fprintf('  Using SPM for DTI->T1 registration...\n');

b0_file = registration_data.reference_volumes.b0_file;
t1_file = registration_data.input.t1_file;

try
    % Load images into SPM format
    VG = spm_vol(t1_file);  % Reference (T1)
    VF = spm_vol(b0_file);  % Source (b0)
    
    fprintf('    Running SPM coregistration...\n');
    
    % Set up coregistration parameters
    flags = struct();
    flags.cost_fun = 'nmi';  % Normalized mutual information
    flags.sep = [4 2];       % Multi-resolution sampling
    flags.tol = [0.02 0.02 0.02 0.001 0.001 0.001]; % Tolerance
    flags.fwhm = [7 7];      % Smoothing
    
    % Run coregistration
    M = spm_coreg(VG, VF, flags);
    
    % Create full transformation matrix
    dti_to_t1_transform_matrix = M * VF.mat / VG.mat;
    
    % Save transformation matrix
    save(dti_to_t1_matrix, 'dti_to_t1_transform_matrix');
    
    % Apply transformation to create registered b0
    fprintf('    Creating registered b0 image...\n');
    
    % Update header with new transformation
    VF_reg = VF;
    VF_reg.fname = registered_b0;
    VF_reg.mat = dti_to_t1_transform_matrix * VF.mat;
    
    % Write registered image
    VF_reg = spm_create_vol(VF_reg);
    for i = 1:VF_reg.dim(3)
        img = spm_slice_vol(VF, spm_matrix([0 0 i]), VF_reg.dim(1:2), 1);
        spm_write_plane(VF_reg, img, i);
    end
    
    % Store in registration data
    registration_data.transforms.dti_to_t1_matrix = dti_to_t1_transform_matrix;
    registration_data.transforms.dti_to_t1_file = dti_to_t1_matrix;
    registration_data.registered_images.b0_in_t1 = registered_b0;
    
    fprintf('    ✓ SPM registration successful\n');
    
catch ME
    error('SPM registration failed: %s', ME.message);
end

end