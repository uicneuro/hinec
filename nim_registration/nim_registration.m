function registration_data = nim_registration(dwi_file, t1_file, varargin)
% nim_registration: Comprehensive registration pipeline for DTI analysis
%
% This function performs multi-modal registration between DWI, T1, and MNI spaces
% enabling proper anatomical alignment for tractography and parcellation.
%
% Arguments:
%   dwi_file - Path to preprocessed DWI file (e.g., 'sample.nii.gz')
%   t1_file - Path to T1 anatomical file (e.g., 'sample_T1.nii.gz')
%   options - Structure with registration options (optional):
%     .output_dir - Output directory (default: same as DWI file)
%     .register_to_mni - Boolean for MNI registration (default: true)
%     .mni_template - Path to MNI template (default: FSL MNI152)
%     .registration_method - 'fsl' or 'spm' (default: 'fsl')
%     .dti_reg_dof - DOF for DTI->T1 registration (default: 6 for rigid)
%     .t1_mni_reg_type - 'linear' or 'nonlinear' (default: 'nonlinear')
%     .force_recompute - Force recomputation of existing transforms (default: false)
%
% Returns:
%   registration_data - Structure containing:
%     .transforms - All transformation matrices and warp fields
%     .registered_images - Paths to registered images
%     .quality_metrics - Registration quality measures
%     .spaces - Information about different coordinate spaces

fprintf('=== HINEC Registration Pipeline ===\n');

% Parse input arguments and set defaults
if nargin < 2
    error('Both DWI and T1 files are required');
end

% Handle options
if nargin > 2 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default options
options = set_registration_defaults(options, dwi_file);

% Validate input files
validate_input_files(dwi_file, t1_file, options);

% Initialize registration data structure
registration_data = initialize_registration_data(dwi_file, t1_file, options);

fprintf('Registration setup:\n');
fprintf('  DWI file: %s\n', dwi_file);
fprintf('  T1 file: %s\n', t1_file);
fprintf('  Output directory: %s\n', options.output_dir);
fprintf('  Method: %s\n', options.registration_method);
fprintf('  Register to MNI: %s\n', string(options.register_to_mni));

%% Step 1: Extract reference volumes
fprintf('\n--- Step 1: Extracting reference volumes ---\n');
registration_data = extract_reference_volumes(registration_data, options);

%% Step 2: DTI to T1 registration
fprintf('\n--- Step 2: DTI to T1 registration ---\n');
registration_data = register_dti_to_t1(registration_data, options);

%% Step 3: T1 to MNI registration (if requested)
if options.register_to_mni
    fprintf('\n--- Step 3: T1 to MNI registration ---\n');
    registration_data = register_t1_to_mni(registration_data, options);
end

%% Step 4: Compute quality metrics
fprintf('\n--- Step 4: Computing registration quality ---\n');
registration_data = compute_registration_quality(registration_data, options);

%% Step 5: Save registration data
fprintf('\n--- Step 5: Saving registration data ---\n');
save_registration_data(registration_data, options);

%% Step 6: Generate registration report
fprintf('\n--- Step 6: Generating registration report ---\n');
generate_registration_report(registration_data, options);

fprintf('\n=== Registration Pipeline Complete ===\n');
fprintf('Registration data saved to: %s\n', registration_data.output_file);

end

function options = set_registration_defaults(options, dwi_file)
% Set default registration options

% Get output directory from DWI file if not specified
if ~isfield(options, 'output_dir') || isempty(options.output_dir)
    [dwi_dir, ~, ~] = fileparts(dwi_file);
    if isempty(dwi_dir)
        options.output_dir = pwd;
    else
        options.output_dir = dwi_dir;
    end
end

% Registration method
if ~isfield(options, 'registration_method')
    options.registration_method = 'fsl';  % FSL is more robust for DTI
end

% MNI registration
if ~isfield(options, 'register_to_mni')
    options.register_to_mni = true;
end

% MNI template
if ~isfield(options, 'mni_template') || isempty(options.mni_template)
    fsl_path = getenv('FSLDIR');
    if ~isempty(fsl_path)
        options.mni_template = fullfile(fsl_path, 'data', 'standard', 'MNI152_T1_1mm.nii.gz');
    else
        warning('FSL not found, MNI registration may fail');
        options.register_to_mni = false;
    end
end

% Registration parameters
if ~isfield(options, 'dti_reg_dof')
    options.dti_reg_dof = 6;  % Rigid registration for DTI->T1
end

if ~isfield(options, 't1_mni_reg_type')
    options.t1_mni_reg_type = 'nonlinear';  % Nonlinear for T1->MNI
end

if ~isfield(options, 'force_recompute')
    options.force_recompute = false;
end

% Quality control
if ~isfield(options, 'generate_qc_images')
    options.generate_qc_images = true;
end

end

function validate_input_files(dwi_file, t1_file, options)
% Validate that input files exist and are readable

if ~isfile(dwi_file)
    error('DWI file not found: %s', dwi_file);
end

if ~isfile(t1_file)
    error('T1 file not found: %s', t1_file);
end

if options.register_to_mni && ~isfile(options.mni_template)
    error('MNI template not found: %s', options.mni_template);
end

% Check FSL availability for FSL-based registration
if strcmp(options.registration_method, 'fsl')
    fsl_path = getenv('FSLDIR');
    if isempty(fsl_path)
        error('FSL not found. Please set FSLDIR environment variable or use registration_method = ''spm''');
    end
end

% Create output directory if it doesn't exist
if ~exist(options.output_dir, 'dir')
    mkdir(options.output_dir);
    fprintf('Created output directory: %s\n', options.output_dir);
end

end

function registration_data = initialize_registration_data(dwi_file, t1_file, options)
% Initialize the registration data structure

registration_data = struct();

% Input files
registration_data.input = struct();
registration_data.input.dwi_file = dwi_file;
registration_data.input.t1_file = t1_file;
registration_data.input.mni_template = options.mni_template;

% Output paths
[~, dwi_name, ~] = fileparts(dwi_file);
dwi_name = strrep(dwi_name, '.nii', ''); % Handle .nii.gz
registration_data.output_dir = options.output_dir;
registration_data.output_prefix = fullfile(options.output_dir, dwi_name);
registration_data.output_file = [registration_data.output_prefix '_registration.mat'];

% Initialize transform storage
registration_data.transforms = struct();
registration_data.registered_images = struct();
registration_data.quality_metrics = struct();
registration_data.spaces = struct();

% Timestamps
registration_data.created = datetime('now');
registration_data.options = options;

end