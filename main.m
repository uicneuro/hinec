function main(imgpath, nimpath, varargin)
% main_with_registration: Enhanced HINEC main function with proper registration
%
% Arguments:
%   imgpath - Path to NiFTI image file (without extension)
%   nimpath - Path to save processed .mat file (must end in `.mat`)
%   t1_file - Path to T1 anatomical file (optional)
%   options - Structure with processing options (optional):
%     .enable_registration - Enable multi-modal registration (default: false if no T1, true if T1 provided)
%     .registration_method - 'fsl' or 'spm' (default: 'fsl')
%     .register_to_mni - Register to MNI space (default: true if registration enabled)
%     .force_recompute_registration - Force recomputation of registration (default: false)
%     .preprocessing_options - Options for preprocessing (see nim_preprocessing)

arguments
    % Path to NiFTI image file
    imgpath string
    
    % Path to save processed .mat file (must end in `.mat`)
    nimpath string
end

arguments (Repeating)
    % Variable arguments for advanced options (T1 file, options struct)
    varargin
end

% Convert inputs to char arrays for consistent handling
imgpath = char(imgpath);
nimpath = char(nimpath);

% Parse variable arguments
t1_file = '';
options = struct();
config = struct();
run_info = struct();

if nargin >= 3
    arg1 = varargin{1};

    % Check if argument is a YAML config structure
    if isstruct(arg1) && isfield(arg1, 'preprocessing')
        % YAML config structure
        config = arg1;
        fprintf('Using YAML configuration for preprocessing\n');

        % Extract preprocessing options from config
        if isfield(config, 'preprocessing')
            options.preprocessing_options = config.preprocessing;
        end

        % Check for T1 in preprocessing config
        if isfield(config.preprocessing, 't1_available') && config.preprocessing.t1_available
            % Look for T1 file based on imgpath
            t1_file = [char(imgpath) '_T1.nii.gz'];
            if ~isfile(t1_file)
                warning('T1 configured but file not found: %s', t1_file);
                t1_file = '';
            end
        end

        % Check for run_info in second varargin
        if nargin >= 4 && isstruct(varargin{2}) && isfield(varargin{2}, 'run_dir')
            run_info = varargin{2};
            fprintf('Using run directory: %s\n', run_info.run_dir);
        end
    elseif ischar(arg1) || isstring(arg1)
        % Third argument is T1 file (legacy)
        t1_file = char(arg1);
        if nargin >= 4 && isstruct(varargin{2})
            options = varargin{2};
        end
    elseif isstruct(arg1)
        % Check if this is run_info structure
        if isfield(arg1, 'run_dir')
            run_info = arg1;
            fprintf('Using run directory: %s\n', run_info.run_dir);
        else
            % Third argument is options (legacy)
            options = arg1;
            if isfield(options, 't1_file')
                t1_file = options.t1_file;
            end
        end
    end
end

% Check if using run directory organization
use_run_dir = ~isempty(fieldnames(run_info));

% Set registration defaults
if ~isfield(options, 'enable_registration')
    options.enable_registration = ~isempty(t1_file) && any(isfile(t1_file));
end

if ~isfield(options, 'registration_method')
    options.registration_method = 'fsl';
end

if ~isfield(options, 'register_to_mni')
    options.register_to_mni = options.enable_registration;
end

if ~isfield(options, 'force_recompute_registration')
    options.force_recompute_registration = false;
end

% Include folders to path
addpath('nim_preprocessing/');
addpath('nim_plots');
addpath('nim_utils');
addpath('nim_calculation');
addpath('nim_parcellation');
addpath('nim_tractography');
addpath('nim_registration');  % Add registration module
addpath('nim_tests');
addpath(genpath('spm12'));
addpath('utils');
addpath('nifti_sample');
addpath('bfgs');

fprintf('=== HINEC Enhanced Pipeline with Registration ===\n');
if options.enable_registration
    fprintf('Registration enabled: %s method\n', options.registration_method);
    fprintf('T1 file: %s\n', t1_file);
else
    fprintf('Registration disabled - using original pipeline\n');
end

%% Data Type Detection and Preprocessing
[img_file, raw_file, t1_file, mask_file] = setup_file_paths(imgpath, run_info);
[is_raw_data, is_preprocessed_data] = detect_data_type(img_file, raw_file);

if is_preprocessed_data
    handle_preprocessed_data(img_file, mask_file, imgpath, run_info);
elseif is_raw_data
    handle_raw_data(img_file, raw_file, t1_file, imgpath, options, run_info);
else
    error('No valid data found. Expected either:\n  - Preprocessed: %s\n  - Raw: %s', char(img_file), char(raw_file));
end

start_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC START: %s\n", string(start_time));

%% Step 1: Load NiFTI image
nim = nim_read(imgpath);

%% Step 2: Calculate diffusion tensor and metrics
nim = nim_dt_spd(nim);
nim = nim_eig(nim);
nim = nim_fa(nim);

%% Step 3: Registration (if enabled)
registration_data = [];
if options.enable_registration
    fprintf('\n=== Running Multi-Modal Registration ===\n');
    
    % Check if registration already exists
    [output_dir, ~, ~] = fileparts(imgpath);
    registration_file = fullfile(output_dir, [strrep(char(imgpath), '_raw', '') '_registration.mat']);
    
    if ~options.force_recompute_registration && any(isfile(registration_file))
        fprintf('Loading existing registration data...\n');
        reg_data = load(registration_file);
        registration_data = reg_data.registration_data;
    else
        % Run registration
        reg_options = struct();
        reg_options.output_dir = output_dir;
        reg_options.registration_method = options.registration_method;
        reg_options.register_to_mni = options.register_to_mni;
        reg_options.force_recompute = options.force_recompute_registration;
        
        registration_data = nim_registration(img_file, t1_file, reg_options);
    end
    
    % Store registration data in nim structure
    nim.registration = registration_data;
    fprintf('✓ Registration data integrated into nim structure\n');
end

%% Step 4: Enhanced Parcellation (with proper registration)
if use_run_dir
    parcellation_mask_file = fullfile(run_info.intermediate_dir, 'parcellation_mask.nii.gz');
else
    [output_dir, ~, ~] = fileparts(imgpath);
    parcellation_mask_file = fullfile(output_dir, 'parcellation_mask.nii.gz');
end

%% Parcellation: Load or generate as needed
if options.enable_registration
    fprintf('\n=== Running Enhanced Parcellation with Registration ===\n');
    nim = nim_parcellation_registered(nim, registration_data, parcellation_mask_file);
else
    fprintf('\n=== Parcellation ===\n');

    % Setup options for generation if needed
    parc_opts = struct();
    parc_opts.dwi_file = img_file;
    if isfield(options, 'atlas_type')
        parc_opts.atlas_type = options.atlas_type;
    else
        parc_opts.atlas_type = 'jhu';
    end

    % nim_parcellation will generate if file doesn't exist
    nim = nim_parcellation(nim, parcellation_mask_file, parc_opts);
end

% Store parcellation mask file path for reference
nim.parcellation_mask_file = parcellation_mask_file;

% Load parcellation labels
nim = nim_load_labels(nim);

%% Step 5: Brain mask improvement using FA data (final step)
fprintf("Improving brain mask using FA data...\n");

% Use run directory if specified, otherwise use imgpath location
if use_run_dir
    [~, base_name, ~] = fileparts(imgpath);
    brain_mask_file = fullfile(run_info.intermediate_dir, [base_name '_M.nii.gz']);
    output_prefix = fullfile(run_info.intermediate_dir, base_name);
else
    brain_mask_file = [char(imgpath) '_M.nii.gz'];
    output_prefix = char(imgpath);
end

if isfile(brain_mask_file)
    % Improve the brain mask using FA data directly from nim structure
    improved_mask_file = preproc_mask_improvement(brain_mask_file, nim.FA, output_prefix);

    % Update the mask in the nim structure if improvement was successful
    if isfile(improved_mask_file)
        % Load the improved mask
        V_improved = spm_vol(improved_mask_file);
        improved_mask_data = spm_read_vols(V_improved);

        % Update the brain mask in nim structure if it exists
        if isfield(nim, 'mask')
            nim.mask = improved_mask_data;
        end

        fprintf("✓ Brain mask improved and updated\n");

        % Replace the original mask file with the improved one
        copyfile(improved_mask_file, brain_mask_file);
    else
        fprintf("⚠ Brain mask improvement failed, keeping original\n");
    end
else
    fprintf("⚠ No brain mask found to improve\n");
end

%% Step 5b: Tissue segmentation for Anatomically Constrained Tractography (ACT)
fprintf("\nGenerating tissue masks for ACT (WM, GM, CSF)...\n");
if isfile(brain_mask_file)
    try
        % Generate tissue-specific masks using FA-based segmentation
        [wm_mask_file, gm_mask_file, csf_mask_file] = ...
            preproc_tissue_segmentation(nim.FA, brain_mask_file, output_prefix);

        % Load tissue masks into nim structure for tractography
        if isfile(wm_mask_file)
            V_wm = spm_vol(wm_mask_file);
            nim.wm_mask = spm_read_vols(V_wm);
            fprintf("  ✓ WM mask loaded: %s\n", wm_mask_file);
        end

        if isfile(gm_mask_file)
            V_gm = spm_vol(gm_mask_file);
            nim.gm_mask = spm_read_vols(V_gm);
            fprintf("  ✓ GM mask loaded: %s\n", gm_mask_file);
        end

        if isfile(csf_mask_file)
            V_csf = spm_vol(csf_mask_file);
            nim.csf_mask = spm_read_vols(V_csf);
            fprintf("  ✓ CSF mask loaded: %s\n", csf_mask_file);
        end

        % Store file paths for reference
        nim.wm_mask_file = wm_mask_file;
        nim.gm_mask_file = gm_mask_file;
        nim.csf_mask_file = csf_mask_file;

        fprintf("✅ Tissue segmentation complete and loaded into nim structure\n");
        fprintf("   Ready for high-order tractography with ACT\n");

    catch ME
        fprintf("⚠ WARNING: Tissue segmentation failed: %s\n", ME.message);
        fprintf("   Tractography will proceed without ACT constraints\n");
    end
else
    fprintf("⚠ No brain mask found for tissue segmentation\n");
    fprintf("   Tractography will proceed without ACT constraints\n");
end

%% Step 6: Save enhanced nim structure
% Redirect output to run directory if specified
if use_run_dir
    % Extract just the filename (no directory) from nimpath
    [~, output_name, output_ext] = fileparts(nimpath);
    fprintf('DEBUG: nimpath = %s\n', nimpath);
    fprintf('DEBUG: output_name = %s\n', output_name);
    fprintf('DEBUG: output_ext = %s\n', output_ext);
    fprintf('DEBUG: run_info.output_dir = %s\n', run_info.output_dir);
    final_nimpath = fullfile(run_info.output_dir, [output_name output_ext]);
    fprintf('DEBUG: final_nimpath = %s\n', final_nimpath);
    nim_save(nim, final_nimpath);
    fprintf('✓ Saved to run directory: %s\n', final_nimpath);

    % Also save run_info to nim structure for reference
    nim.run_info = run_info;
else
    nim_save(nim, nimpath);
end

end_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC END: %s\n", string(end_time));

print_pipeline_summary(options, registration_data, imgpath, nimpath, run_info);
end

%% Helper Functions

function [img_file, raw_file, t1_file, mask_file] = setup_file_paths(imgpath, run_info)
% Setup standard file paths based on input prefix
% If run_info provided, use intermediate directory for generated files

    % Convert imgpath to char array first to ensure compatibility
    imgpath = char(imgpath);

    % Input files are always in original location
    raw_file = [imgpath '_raw.nii.gz'];
    t1_file = [imgpath '_T1.nii.gz'];

    % Output files go to run directory if specified
    if ~isempty(fieldnames(run_info))
        % Use run directory for intermediate files
        [~, base_name, ~] = fileparts(imgpath);
        img_file = char(fullfile(run_info.intermediate_dir, [base_name '.nii.gz']));
        mask_file = char(fullfile(run_info.intermediate_dir, [base_name '_M.nii.gz']));
    else
        % Legacy: keep in same directory as input
        img_file = [imgpath '.nii.gz'];
        mask_file = [imgpath '_M.nii.gz'];
    end
end

function [is_raw, is_preprocessed] = detect_data_type(img_file, raw_file)
% Detect whether data is raw or preprocessed
    % Ensure scalar logical values
    is_raw = any(isfile(raw_file));
    is_preprocessed = any(isfile(img_file)) && ~is_raw;
end

function handle_preprocessed_data(img_file, mask_file, imgpath, run_info)
% Handle preprocessed data: generate auxiliary files only
    % Convert imgpath to char array for consistent handling
    imgpath = char(imgpath);

    % Original file location (where the actual preprocessed file is)
    original_img_file = [imgpath '.nii.gz'];
    original_mask_file = [imgpath '_M.nii.gz'];

    fprintf("=== PREPROCESSED DATA DETECTED ===\n");
    fprintf("Found: %s\n", original_img_file);
    fprintf("Strategy: Generate auxiliary files (no NIfTI modification)\n\n");

    % Generate brain mask if needed (in original location)
    if ~isfile(original_mask_file)
        fprintf("--- Generating Brain Mask ---\n");
        [output_dir, ~, ~] = fileparts(imgpath);
        if isempty(output_dir)
            output_dir = pwd;
        end
        preproc_brain_extraction(original_img_file, output_dir, original_mask_file);
    else
        fprintf("✓ Brain mask exists: %s\n", original_mask_file);
    end

    % If using run directory, copy files to run directory
    use_run_dir = ~isempty(fieldnames(run_info));
    if use_run_dir
        fprintf("\n--- Copying preprocessed files to run directory ---\n");

        % Copy main preprocessed file
        if isfile(original_img_file)
            copyfile(original_img_file, img_file);
            fprintf("  Copied: %s -> %s\n", original_img_file, img_file);
        end

        % Copy brain mask
        if isfile(original_mask_file)
            copyfile(original_mask_file, mask_file);
            fprintf("  Copied: %s -> %s\n", original_mask_file, mask_file);
        end

        fprintf("✓ Files copied to run directory\n");
    end

    fprintf("\n✓ Preprocessed data ready for DTI analysis\n");
end

function handle_raw_data(img_file, raw_file, t1_file, imgpath, options, run_info)
% Handle raw data: run full preprocessing pipeline
    % Convert imgpath to char array for consistent handling
    imgpath = char(imgpath);

    fprintf("=== RAW DATA DETECTED ===\n");
    fprintf("Found: %s\n", raw_file);
    fprintf("Strategy: Full preprocessing pipeline\n\n");

    % Check for T1 data
    t1_available = isfile(t1_file);
    if t1_available
        fprintf("✓ T1 structural data: %s\n", t1_file);
    else
        fprintf("ℹ No T1 structural data\n");
    end

    % Setup preprocessing options
    preproc_options = setup_preprocessing_options(options, t1_available, t1_file);

    % Run preprocessing to original location (preprocessing doesn't support output_dir redirection)
    fprintf("\n--- Starting Preprocessing ---\n");
    nim_preprocessing(imgpath, preproc_options);

    % Original preprocessed file location
    original_img_file = [imgpath '.nii.gz'];

    % Verify preprocessing succeeded
    if ~isfile(original_img_file)
        error('Preprocessing failed. Output not found: %s', original_img_file);
    end
    fprintf("✓ Preprocessing complete: %s\n", original_img_file);

    % If using run directory, copy preprocessed files to run directory
    use_run_dir = ~isempty(fieldnames(run_info));
    if use_run_dir
        fprintf("\n--- Copying preprocessed files to run directory ---\n");

        % Copy main preprocessed file
        copyfile(original_img_file, img_file);
        fprintf("  Copied: %s -> %s\n", original_img_file, img_file);

        % Copy brain mask if exists
        original_mask = [imgpath '_M.nii.gz'];
        if isfile(original_mask)
            [~, base_name, ~] = fileparts(imgpath);
            dest_mask = char(fullfile(run_info.intermediate_dir, [base_name '_M.nii.gz']));
            copyfile(original_mask, dest_mask);
            fprintf("  Copied: %s -> %s\n", original_mask, dest_mask);
        end

        fprintf("✓ Files copied to run directory\n");
    end
end

function preproc_options = setup_preprocessing_options(options, t1_available, t1_file)
% Setup preprocessing options with T1 integration if available
    if isfield(options, 'preprocessing_options')
        preproc_options = options.preprocessing_options;
    else
        preproc_options = struct();
        preproc_options.run_denoising = true;
        preproc_options.denoise_method = 'dwidenoise';
        preproc_options.run_motion_correction = true;
        preproc_options.run_eddy = true;
        preproc_options.improve_mask = true;
        preproc_options.atlas_type = 'jhu';
    end

    % Add T1 integration
    preproc_options.t1_available = t1_available;
    preproc_options.t1_file = t1_file;
    preproc_options.use_t1_registration = t1_available;
end

function print_pipeline_summary(options, registration_data, imgpath, nimpath, run_info)
% Print pipeline completion summary
    if options.enable_registration
        fprintf('\n=== Registration-Enhanced Pipeline Complete ===\n');
        fprintf('Enhanced features:\n');
        fprintf('  • DTI ↔ T1 ↔ MNI registration\n');
        fprintf('  • Atlas-based parcellation\n');
        fprintf('  • Cross-modal ROI transformations\n');
        fprintf('  • Improved tractography seeds\n');

        if isfield(registration_data, 'quality_metrics')
            fprintf('\nQuality metrics:\n');
            if isfield(registration_data.quality_metrics, 'dti_t1_nmi')
                fprintf('  DTI→T1 NMI: %.4f\n', registration_data.quality_metrics.dti_t1_nmi);
            end
            if isfield(registration_data.quality_metrics, 't1_mni_nmi')
                fprintf('  T1→MNI NMI: %.4f\n', registration_data.quality_metrics.t1_mni_nmi);
            end
        end

        fprintf('\nReport: %s_registration_report.html\n', strrep(char(imgpath), '_raw', ''));
    else
        fprintf('\n=== Standard Pipeline Complete ===\n');
        fprintf('To enable registration, provide T1:\n');
        fprintf('  main(''%s'', ''%s'', ''path/to/T1.nii.gz'')\n', imgpath, nimpath);
    end

    % Print run directory info if used
    if ~isempty(fieldnames(run_info))
        fprintf('\n=== Run Directory Organization ===\n');
        fprintf('Run ID: %s\n', run_info.run_id);
        fprintf('Location: %s\n', run_info.run_dir);
        fprintf('  Logs: %s\n', run_info.logs_dir);
        fprintf('  Intermediate: %s\n', run_info.intermediate_dir);
        fprintf('  Output: %s\n', run_info.output_dir);
        fprintf('  Tractography: %s\n', run_info.tractography_dir);
    end
end