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
    case 'spm'
        registration_data = register_t1_to_mni_spm(registration_data, options, ...
            t1_to_mni_matrix, registered_t1_mni, mni_in_t1);
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

if status ~= 0
    warning('T1 brain extraction failed: %s', cmdout);
    t1_brain = t1_file; % Use original if brain extraction fails
end

% Step 2: Linear registration (FLIRT)
fprintf('    Running linear registration (FLIRT)...\n');
t1_to_mni_linear = [registration_data.output_prefix '_t1_to_mni_linear.nii.gz'];

cmd_flirt = sprintf(['%s/bin/flirt -in %s -ref %s -out %s -omat %s ' ...
                    '-cost corratio -dof 12 -searchrx -90 90 ' ...
                    '-searchry -90 90 -searchrz -90 90 -interp trilinear'], ...
                    fsl_path, t1_brain, mni_template, t1_to_mni_linear, t1_to_mni_transform);

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
    
    % Run FNIRT
    if isempty(fnirt_config)
        cmd_fnirt = sprintf(['%s/bin/fnirt --in=%s --ref=%s --aff=%s ' ...
                            '--iout=%s --fout=%s --jout=%s --refmask=%s/data/standard/MNI152_T1_2mm_brain_mask_dil ' ...
                            '--warpres=10,10,10 --subsamp=8,4,2,1 --miter=5,5,5,5 --lambda=240,120,90,30 ' ...
                            '--ssqlambda=1 --regmod=bending_energy --estint=1,1,1 --applyrefmask=0,0,1 ' ...
                            '--applyinmask=0,0,1 --verbose'], ...
                            fsl_path, t1_file, mni_template, t1_to_mni_transform, ...
                            registered_t1_mni, t1_to_mni_warp, mni_to_t1_warp, fsl_path);
    else
        cmd_fnirt = sprintf(['%s/bin/fnirt --in=%s --ref=%s --aff=%s ' ...
                            '--config=%s --iout=%s --fout=%s --jout=%s'], ...
                            fsl_path, t1_file, mni_template, t1_to_mni_transform, ...
                            fnirt_config, registered_t1_mni, t1_to_mni_warp, mni_to_t1_warp);
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

function registration_data = register_t1_to_mni_spm(registration_data, options, ...
    t1_to_mni_matrix, registered_t1_mni, mni_in_t1)
% Register T1 to MNI using SPM tools

fprintf('  Using SPM for T1->MNI registration...\n');

t1_file = registration_data.input.t1_file;
mni_template = registration_data.input.mni_template;

try
    % Load images
    VG = spm_vol(mni_template);  % Reference (MNI)
    VF = spm_vol(t1_file);       % Source (T1)
    
    fprintf('    Running SPM normalization...\n');
    
    if strcmp(options.t1_mni_reg_type, 'nonlinear')
        % Use SPM12 unified segmentation/normalization
        fprintf('    Using SPM12 unified segmentation...\n');
        
        % Set up batch for unified segmentation
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[t1_file ',1']};
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
        
        % Tissue probability maps
        for c = 1:6
            matlabbatch{1}.spm.spatial.preproc.tissue(c).tpm = {sprintf('%s,%d', ...
                fullfile(spm('Dir'), 'tpm', 'TPM.nii'), c)};
            matlabbatch{1}.spm.spatial.preproc.tissue(c).ngaus = [1 1 2 3 4 2];
            matlabbatch{1}.spm.spatial.preproc.tissue(c).native = [1 0];  % Native space
            matlabbatch{1}.spm.spatial.preproc.tissue(c).warped = [0 0];  % No warped output
        end
        
        % Warping options
        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];  % Forward and inverse
        
        % Run batch
        spm_jobman('run', matlabbatch);
        
        % Get deformation field paths
        [pth, nam, ext] = fileparts(t1_file);
        forward_deformation = fullfile(pth, ['y_' nam ext]);
        inverse_deformation = fullfile(pth, ['iy_' nam ext]);
        
        % Apply deformation to create registered T1
        if isfile(forward_deformation)
            fprintf('    Applying deformation field...\n');
            
            % Set up batch for applying deformation
            clear matlabbatch;
            matlabbatch{1}.spm.spatial.normalise.write.subj.def = {forward_deformation};
            matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {t1_file};
            matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
            matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
            
            spm_jobman('run', matlabbatch);
            
            % Move result to expected location
            warped_t1 = fullfile(pth, ['w' nam ext]);
            if isfile(warped_t1)
                movefile(warped_t1, registered_t1_mni);
            end
        end
        
        % Store transformation data
        t1_to_mni_data = struct();
        t1_to_mni_data.type = 'nonlinear';
        t1_to_mni_data.forward_deformation = forward_deformation;
        t1_to_mni_data.inverse_deformation = inverse_deformation;
        
    else
        % Linear registration only
        fprintf('    Running linear coregistration...\n');
        
        % Set up coregistration parameters
        flags = struct();
        flags.cost_fun = 'nmi';
        flags.sep = [4 2];
        flags.tol = [0.02 0.02 0.02 0.001 0.001 0.001];
        flags.fwhm = [7 7];
        
        % Run coregistration
        M = spm_coreg(VG, VF, flags);
        
        % Create transformation matrix
        t1_to_mni_transform_matrix = M * VF.mat / VG.mat;
        
        % Apply transformation to create registered T1
        VF_reg = VF;
        VF_reg.fname = registered_t1_mni;
        VF_reg.mat = t1_to_mni_transform_matrix * VF.mat;
        
        % Write registered image
        VF_reg = spm_create_vol(VF_reg);
        for i = 1:VF_reg.dim(3)
            img = spm_slice_vol(VF, spm_matrix([0 0 i]), VF_reg.dim(1:2), 1);
            spm_write_plane(VF_reg, img, i);
        end
        
        % Store transformation data
        t1_to_mni_data = struct();
        t1_to_mni_data.type = 'linear';
        t1_to_mni_data.transform_matrix = t1_to_mni_transform_matrix;
    end
    
    % Save transformation data
    save(t1_to_mni_matrix, 't1_to_mni_data');
    
    % Store in registration data
    registration_data.transforms.t1_to_mni = t1_to_mni_data;
    registration_data.transforms.t1_to_mni_file = t1_to_mni_matrix;
    registration_data.registered_images.t1_in_mni = registered_t1_mni;
    
    fprintf('    ✓ SPM registration successful\n');
    
catch ME
    error('SPM registration failed: %s', ME.message);
end

end