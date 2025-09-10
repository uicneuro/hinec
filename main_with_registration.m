function main_with_registration(imgpath, nimpath, varargin)
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

% Parse variable arguments
t1_file = '';
options = struct();

if nargin >= 3
    if ischar(varargin{1}) || isstring(varargin{1})
        % Third argument is T1 file
        t1_file = char(varargin{1});
        if nargin >= 4 && isstruct(varargin{2})
            options = varargin{2};
        end
    elseif isstruct(varargin{1})
        % Third argument is options
        options = varargin{1};
        if isfield(options, 't1_file')
            t1_file = options.t1_file;
        end
    end
end

% Set registration defaults
if ~isfield(options, 'enable_registration')
    options.enable_registration = ~isempty(t1_file) && isfile(t1_file);
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

% Define the file extensions and suffixes
img_file = [char(imgpath) '.nii.gz'];
raw_file = [char(imgpath) '_raw.nii.gz'];

% Check if the processed NIfTI file exists
if isfile(img_file)
    fprintf("Found processed NIfTI image: %s\n", img_file);
else
    fprintf("Processed NIfTI image not found: %s\n", img_file);
    
    % Check if the raw data exists
    if isfile(raw_file)
        fprintf("Found raw data: %s\n", raw_file);
        fprintf("Starting preprocessing...\n");
        
        % Set up preprocessing options
        if isfield(options, 'preprocessing_options')
            preproc_options = options.preprocessing_options;
        else
            preproc_options = struct();
            preproc_options.run_denoising = true;
            preproc_options.denoise_method = 'dwidenoise';
            preproc_options.run_motion_correction = true;
            preproc_options.run_eddy = true;
            preproc_options.improve_mask = true;
            preproc_options.atlas_type = 'JHU-tract';
        end
        
        % Preprocess the raw data
        nim_preprocessing(imgpath, preproc_options);
        
        % Check if the preprocessing was successful
        if isfile(img_file)
            fprintf("Preprocessing completed. Processed file: %s\n", img_file);
        else
            error('Preprocessing failed. Processed NIfTI image not found: %s\n', img_file);
        end
    else
        error('Raw data not found: %s\n', raw_file);
    end
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
    
    if ~options.force_recompute_registration && isfile(registration_file)
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
[output_dir, ~, ~] = fileparts(imgpath);
parcellation_mask_file = fullfile(output_dir, 'parcellation_mask.nii.gz');

if options.enable_registration
    fprintf('\n=== Running Enhanced Parcellation with Registration ===\n');
    nim = nim_parcellation_registered(nim, registration_data, parcellation_mask_file);
else
    fprintf('\n=== Running Standard Parcellation ===\n');
    nim = nim_parcellation(nim, parcellation_mask_file);
end

% Store parcellation mask file path for reference
nim.parcellation_mask_file = parcellation_mask_file;

% Load parcellation labels
nim = nim_load_labels(nim);

%% Step 5: Brain mask improvement using FA data (final step)
fprintf("Improving brain mask using FA data...\n");
brain_mask_file = [char(imgpath) '_M.nii.gz'];
if isfile(brain_mask_file)
    % Improve the brain mask using FA data directly from nim structure
    improved_mask_file = preproc_mask_improvement(brain_mask_file, nim.FA, char(imgpath));
    
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

%% Step 6: Save enhanced nim structure
nim_save(nim, nimpath);

end_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC END: %s\n", string(end_time));

if options.enable_registration
    fprintf('\n=== Registration-Enhanced Pipeline Complete ===\n');
    fprintf('Enhanced features available:\n');
    fprintf('  • Proper DTI ↔ T1 ↔ MNI registration\n');
    fprintf('  • Accurate atlas-based parcellation\n');
    fprintf('  • Cross-modal ROI transformations\n');
    fprintf('  • Improved tractography seed masks\n');
    
    if isfield(registration_data, 'quality_metrics')
        fprintf('\nRegistration Quality Summary:\n');
        if isfield(registration_data.quality_metrics, 'dti_t1_nmi')
            fprintf('  DTI→T1 NMI: %.4f\n', registration_data.quality_metrics.dti_t1_nmi);
        end
        if isfield(registration_data.quality_metrics, 't1_mni_nmi')
            fprintf('  T1→MNI NMI: %.4f\n', registration_data.quality_metrics.t1_mni_nmi);
        end
    end
    
    fprintf('\nRegistration report: %s_registration_report.html\n', strrep(char(imgpath), '_raw', ''));
else
    fprintf('\n=== Standard Pipeline Complete ===\n');
    fprintf('To enable registration features, provide T1 file:\n');
    fprintf('  main_with_registration(''%s'', ''%s'', ''path/to/T1.nii.gz'')\n', imgpath, nimpath);
end

end