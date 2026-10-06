function registration_data = register_t1_to_mni(registration_data, options)
% register_t1_to_mni: Register T1 anatomical image to MNI template space
%
% This function performs linear and/or nonlinear registration between 
% T1 anatomical image and MNI152 template space.

fprintf('Registering T1 to MNI template space...\n');

% Define output paths
t1_to_mni_matrix = [registration_data.output_prefix '_t1_to_mni.mat'];
t1_to_mni_transform = [registration_data.output_prefix '_t1_to_mni_transform.txt'];
t1_to_mni_warp = [registration_data.output_prefix '_t1_to_mni_warp.nii.gz'];
mni_to_t1_warp = [registration_data.output_prefix '_mni_to_t1_warp.nii.gz'];
registered_t1_mni = [registration_data.output_prefix '_t1_in_mni.nii.gz'];
mni_in_t1 = [registration_data.output_prefix '_mni_in_t1.nii.gz'];

% Check if registration already exists and not forcing recompute
if ~options.force_recompute && isfile(t1_to_mni_matrix)
    fprintf('  T1->MNI registration already exists, loading...\n');
    load(t1_to_mni_matrix, 't1_to_mni_data');
    registration_data.transforms.t1_to_mni = t1_to_mni_data;
    return;
end

% Perform registration based on method
switch lower(options.registration_method)
    case 'fsl'
        registration_data = register_t1_to_mni_fsl(registration_data, options, ...
            t1_to_mni_matrix, t1_to_mni_transform, t1_to_mni_warp, mni_to_t1_warp, ...
            registered_t1_mni, mni_in_t1);
    otherwise
        error('Unknown registration method: %s', options.registration_method);
end

fprintf('  ✓ T1 to MNI registration complete\n');

end

function registration_data = register_t1_to_mni_fsl(registration_data, options, ...
    t1_to_mni_matrix, t1_to_mni_transform, t1_to_mni_warp, mni_to_t1_warp, ...
    registered_t1_mni, mni_in_t1)
% Register T1 to MNI using FSL tools (FLIRT + FNIRT)

fsl_path = getenv('FSLDIR');
t1_file = registration_data.input.t1_file;
mni_template = registration_data.input.mni_template;

fprintf('  Using FSL for T1->MNI registration...\n');

% Step 1: Brain extraction for better registration
fprintf('    Extracting brain from T1...\n');
t1_brain = [registration_data.output_prefix '_t1_brain.nii.gz'];
t1_brain_mask = [registration_data.output_prefix '_t1_brain_mask.nii.gz'];

cmd_bet = sprintf('%s/bin/bet %s %s -f 0.5 -B -m', fsl_path, t1_file, ...
    strrep(t1_brain, '.nii.gz', ''));
[status, cmdout] = system(cmd_bet);

if status ~= 0 || ~isfile(t1_brain)
    % No silent fallback: registering the whole-head T1 to the brain template is
    % the same template/input mismatch that corrupted the registration.
    error('register_t1_to_mni:bet', 'T1 brain extraction failed: %s', cmdout);
end

% Step 2: Linear registration (FLIRT)
fprintf('    Running linear registration (FLIRT)...\n');
t1_to_mni_linear = [registration_data.output_prefix '_t1_to_mni_linear.nii.gz'];

% Brain-extracted T1 must be registered to the BRAIN template: against the
% whole-head template FLIRT scales the brain up to fill the skull (MASiVar pilot:
% x1.6 scaling and ~90 deg axis swap). Same recipe as preproc_t1_mni_registration.
mni_brain_template = regexprep(mni_template, '(\.nii)?(\.gz)?$', '_brain.nii.gz');
if ~isfile(mni_brain_template)
    error('register_t1_to_mni:brainTemplate', 'Brain template not found: %s', mni_brain_template);
end
cmd_flirt = sprintf(['%s/bin/flirt -in %s -ref %s -out %s -omat %s ' ...
                    '-cost corratio -dof 12 -searchrx -90 90 ' ...
                    '-searchry -90 90 -searchrz -90 90 -interp trilinear'], ...
                    fsl_path, t1_brain, mni_brain_template, t1_to_mni_linear, t1_to_mni_transform);

[status, cmdout] = system(cmd_flirt);
if status ~= 0
    error('FLIRT linear registration failed: %s', cmdout);
end

% Step 3: Nonlinear registration (FNIRT) if requested
if strcmp(options.t1_mni_reg_type, 'nonlinear')
    fprintf('    Running nonlinear registration (FNIRT)...\n');
    
    % FNIRT configuration for T1->MNI
    fnirt_config = fullfile(fsl_path, 'etc', 'flirtsch', 'T1_2_MNI152_2mm.cnf');
    if ~isfile(fnirt_config)
        fprintf('    Standard FNIRT config not found, using default parameters\n');
        fnirt_config = '';
    end
    
    % FNIRT follows FSL's standard recipe (as preproc_t1_mni_registration): the
    % WHOLE-HEAD T1 as input, and the config's own reference and reference mask
    % (MNI152_T1_2mm + its dilated brain mask). Overriding --ref with the 1 mm
    % template paired it with the config's 2 mm mask. The warp is defined on the
    % 2 mm grid; applywarp resamples it onto any output grid.
    fnirt_ref = fullfile(fsl_path, 'data', 'standard', 'MNI152_T1_2mm.nii.gz');
    fnirt_ref_mask = fullfile(fsl_path, 'data', 'standard', 'MNI152_T1_2mm_brain_mask_dil.nii.gz');
    jacobian_file = regexprep(t1_to_mni_warp, '(\.nii)?(\.gz)?$', '_jacobian.nii.gz');

    % Run FNIRT
    if isempty(fnirt_config)
        cmd_fnirt = sprintf(['%s/bin/fnirt --in=%s --ref=%s --aff=%s ' ...
                            '--iout=%s --fout=%s --jout=%s --refmask=%s ' ...
                            '--warpres=10,10,10 --subsamp=8,4,2,1 --miter=5,5,5,5 --lambda=240,120,90,30 ' ...
                            '--ssqlambda=1 --regmod=bending_energy --estint=1,1,1 --applyrefmask=0,0,1 ' ...
                            '--applyinmask=0,0,1 --verbose'], ...
                            fsl_path, t1_file, fnirt_ref, t1_to_mni_transform, ...
                            registered_t1_mni, t1_to_mni_warp, jacobian_file, fnirt_ref_mask);
    else
        cmd_fnirt = sprintf(['%s/bin/fnirt --in=%s --aff=%s ' ...
                            '--config=%s --iout=%s --fout=%s --jout=%s'], ...
                            fsl_path, t1_file, t1_to_mni_transform, ...
                            fnirt_config, registered_t1_mni, t1_to_mni_warp, jacobian_file);
    end
    
    [status, cmdout] = system(cmd_fnirt);
    if status ~= 0
        warning('FNIRT nonlinear registration failed: %s\nFalling back to linear only', cmdout);
        options.t1_mni_reg_type = 'linear';
        copyfile(t1_to_mni_linear, registered_t1_mni);
    else
        fprintf('    ✓ Nonlinear registration successful\n');
    end
else
    % Use linear registration result
    copyfile(t1_to_mni_linear, registered_t1_mni);
end

% Step 4: Create inverse transform (MNI->T1)
if strcmp(options.t1_mni_reg_type, 'nonlinear') && isfile(t1_to_mni_warp)
    fprintf('    Computing inverse nonlinear transform...\n');
    cmd_invwarp = sprintf('%s/bin/invwarp --ref=%s --warp=%s --out=%s', ...
        fsl_path, t1_file, t1_to_mni_warp, mni_to_t1_warp);
    [status, cmdout] = system(cmd_invwarp);
    
    if status ~= 0
        warning('Inverse warp computation failed: %s', cmdout);
    end
    
    % Apply inverse warp to get MNI in T1 space
    cmd_applywarp = sprintf('%s/bin/applywarp --ref=%s --in=%s --warp=%s --out=%s', ...
        fsl_path, t1_file, mni_template, mni_to_t1_warp, mni_in_t1);
    [status, cmdout] = system(cmd_applywarp);
    
    if status ~= 0
        warning('Inverse warp application failed: %s', cmdout);
    end
else
    % Linear inverse transform
    fprintf('    Computing inverse linear transform...\n');
    mni_to_t1_transform = [registration_data.output_prefix '_mni_to_t1_transform.txt'];
    
    cmd_convert = sprintf('%s/bin/convert_xfm -omat %s -inverse %s', ...
        fsl_path, mni_to_t1_transform, t1_to_mni_transform);
    [status, cmdout] = system(cmd_convert);
    
    if status == 0
        % Apply inverse transform
        cmd_apply_inverse = sprintf('%s/bin/flirt -in %s -ref %s -out %s -init %s -applyxfm', ...
            fsl_path, mni_template, t1_file, mni_in_t1, mni_to_t1_transform);
        [status, cmdout] = system(cmd_apply_inverse);
        
        if status ~= 0
            warning('Inverse transform application failed: %s', cmdout);
        end
    end
end

% Step 5: Store transformation data
t1_to_mni_data = struct();
t1_to_mni_data.type = options.t1_mni_reg_type;
t1_to_mni_data.linear_transform_file = t1_to_mni_transform;

if isfile(t1_to_mni_transform)
    t1_to_mni_data.linear_matrix = load(t1_to_mni_transform, '-ascii');
end

if strcmp(options.t1_mni_reg_type, 'nonlinear')
    t1_to_mni_data.forward_warp = t1_to_mni_warp;
    t1_to_mni_data.inverse_warp = mni_to_t1_warp;
end

% Save transformation data
save(t1_to_mni_matrix, 't1_to_mni_data');

% Store in registration data
registration_data.transforms.t1_to_mni = t1_to_mni_data;
registration_data.transforms.t1_to_mni_file = t1_to_mni_matrix;
registration_data.registered_images.t1_in_mni = registered_t1_mni;
registration_data.registered_images.mni_in_t1 = mni_in_t1;

% Clean up temporary files
if isfile(t1_brain) && ~strcmp(t1_brain, t1_file)
    delete(t1_brain);
end
if isfile(t1_brain_mask)
    delete(t1_brain_mask);
end
if isfile(t1_to_mni_linear)
    delete(t1_to_mni_linear);
end

fprintf('    ✓ T1->MNI registration successful\n');

end
