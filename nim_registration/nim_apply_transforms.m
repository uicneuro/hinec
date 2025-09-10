function output_file = nim_apply_transforms(input_file, registration_data, transform_chain, varargin)
% nim_apply_transforms: Apply transformation chain to image or data
%
% This function applies a series of transformations to move data between
% different coordinate spaces (DTI, T1, MNI).
%
% Arguments:
%   input_file - Path to input image/data file
%   registration_data - Registration data structure from nim_registration
%   transform_chain - Cell array defining transformation chain, e.g.:
%                    {'dti_to_t1', 't1_to_mni'} - DTI -> T1 -> MNI
%                    {'mni_to_t1', 't1_to_dti'} - MNI -> T1 -> DTI
%   options - Structure with application options (optional):
%     .output_file - Output file path (auto-generated if not specified)
%     .interpolation - 'linear', 'nearest', 'spline' (default: 'linear')
%     .reference_space - Reference image for output space
%     .method - 'fsl', 'spm' (default: same as registration)
%
% Returns:
%   output_file - Path to transformed output file

fprintf('=== Applying Transform Chain ===\n');
fprintf('Input: %s\n', input_file);
fprintf('Transform chain: %s\n', strjoin(transform_chain, ' -> '));

% Parse options
if nargin > 3 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set defaults
options = set_transform_defaults(options, input_file, transform_chain, registration_data);

% Validate inputs
validate_transform_inputs(input_file, registration_data, transform_chain);

% Generate output filename if not provided
if ~isfield(options, 'output_file') || isempty(options.output_file)
    options.output_file = generate_output_filename(input_file, transform_chain);
end

% Apply transformation chain
switch lower(options.method)
    case 'fsl'
        output_file = apply_transforms_fsl(input_file, registration_data, transform_chain, options);
    case 'spm'
        output_file = apply_transforms_spm(input_file, registration_data, transform_chain, options);
    otherwise
        error('Unknown transformation method: %s', options.method);
end

fprintf('✓ Transform complete: %s\n', output_file);

end

function options = set_transform_defaults(options, input_file, transform_chain, registration_data)
% Set default options for transformation

if ~isfield(options, 'interpolation')
    options.interpolation = 'linear';
end

if ~isfield(options, 'method')
    % Use same method as registration
    if isfield(registration_data, 'options') && isfield(registration_data.options, 'registration_method')
        options.method = registration_data.options.registration_method;
    else
        options.method = 'fsl';
    end
end

% Determine reference space for output
if ~isfield(options, 'reference_space') || isempty(options.reference_space)
    final_transform = transform_chain{end};
    switch final_transform
        case 'dti_to_t1'
            options.reference_space = registration_data.input.t1_file;
        case 't1_to_mni'
            options.reference_space = registration_data.input.mni_template;
        case 't1_to_dti'
            options.reference_space = registration_data.input.dwi_file;
        case 'mni_to_t1'
            options.reference_space = registration_data.input.t1_file;
        otherwise
            error('Cannot determine reference space for transform: %s', final_transform);
    end
end

end

function validate_transform_inputs(input_file, registration_data, transform_chain)
% Validate transformation inputs

if ~isfile(input_file)
    error('Input file not found: %s', input_file);
end

if ~isstruct(registration_data) || ~isfield(registration_data, 'transforms')
    error('Invalid registration data structure');
end

% Validate each transform in chain exists
for i = 1:length(transform_chain)
    transform = transform_chain{i};
    
    switch transform
        case 'dti_to_t1'
            if ~isfield(registration_data.transforms, 'dti_to_t1_matrix')
                error('DTI to T1 transform not found in registration data');
            end
        case 't1_to_dti'
            if ~isfield(registration_data.transforms, 'dti_to_t1_matrix')
                error('DTI to T1 transform (needed for inverse) not found in registration data');
            end
        case 't1_to_mni'
            if ~isfield(registration_data.transforms, 't1_to_mni')
                error('T1 to MNI transform not found in registration data');
            end
        case 'mni_to_t1'
            if ~isfield(registration_data.transforms, 't1_to_mni')
                error('T1 to MNI transform (needed for inverse) not found in registration data');
            end
        otherwise
            error('Unknown transform: %s', transform);
    end
end

end

function output_file = generate_output_filename(input_file, transform_chain)
% Generate output filename based on input and transform chain

[pth, nam, ext] = fileparts(input_file);
if strcmp(ext, '.gz')
    [~, nam2, ext2] = fileparts(nam);
    nam = nam2;
    ext = [ext2 ext];
end

% Create suffix from transform chain
suffix = '';
for i = 1:length(transform_chain)
    transform = transform_chain{i};
    switch transform
        case 'dti_to_t1'
            suffix = [suffix '_in_t1'];
        case 't1_to_dti'
            suffix = [suffix '_in_dti'];
        case 't1_to_mni'
            suffix = [suffix '_in_mni'];
        case 'mni_to_t1'
            suffix = [suffix '_in_t1'];
    end
end

output_file = fullfile(pth, [nam suffix ext]);

end

function output_file = apply_transforms_fsl(input_file, registration_data, transform_chain, options)
% Apply transformation chain using FSL tools

fprintf('Using FSL for transformation...\n');

fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSL not found. Please set FSLDIR environment variable.');
end

current_file = input_file;
temp_files = {};

% Apply each transform in sequence
for i = 1:length(transform_chain)
    transform = transform_chain{i};
    fprintf('  Applying transform %d/%d: %s\n', i, length(transform_chain), transform);
    
    % Create temporary filename for intermediate results (except last)
    if i == length(transform_chain)
        next_file = options.output_file;
    else
        [pth, nam, ext] = fileparts(current_file);
        next_file = fullfile(pth, sprintf('%s_temp_%d%s', nam, i, ext));
        temp_files{end+1} = next_file;
    end
    
    % Apply specific transform
    switch transform
        case 'dti_to_t1'
            apply_dti_to_t1_fsl(current_file, next_file, registration_data, options, fsl_path);
            
        case 't1_to_dti'
            apply_t1_to_dti_fsl(current_file, next_file, registration_data, options, fsl_path);
            
        case 't1_to_mni'
            apply_t1_to_mni_fsl(current_file, next_file, registration_data, options, fsl_path);
            
        case 'mni_to_t1'
            apply_mni_to_t1_fsl(current_file, next_file, registration_data, options, fsl_path);
            
        otherwise
            error('Unknown transform: %s', transform);
    end
    
    current_file = next_file;
end

% Clean up temporary files
for i = 1:length(temp_files)
    if isfile(temp_files{i})
        delete(temp_files{i});
    end
end

output_file = options.output_file;

end

function apply_dti_to_t1_fsl(input_file, output_file, registration_data, options, fsl_path)
% Apply DTI to T1 transform using FSL

if isfield(registration_data.transforms, 'dti_to_t1_fsl_file')
    transform_file = registration_data.transforms.dti_to_t1_fsl_file;
else
    error('FSL DTI to T1 transform file not found');
end

% Set interpolation method
interp_method = get_fsl_interpolation(options.interpolation);

cmd = sprintf('%s/bin/flirt -in %s -ref %s -out %s -init %s -applyxfm -interp %s', ...
    fsl_path, input_file, options.reference_space, output_file, transform_file, interp_method);

fprintf('    Command: %s\n', cmd);
[status, cmdout] = system(cmd);

if status ~= 0
    error('FSL DTI to T1 transform failed: %s', cmdout);
end

end

function apply_t1_to_dti_fsl(input_file, output_file, registration_data, options, fsl_path)
% Apply T1 to DTI transform using FSL (inverse of DTI to T1)

if isfield(registration_data.transforms, 't1_to_dti_fsl_file')
    transform_file = registration_data.transforms.t1_to_dti_fsl_file;
else
    error('FSL T1 to DTI transform file not found');
end

% Set interpolation method
interp_method = get_fsl_interpolation(options.interpolation);

cmd = sprintf('%s/bin/flirt -in %s -ref %s -out %s -init %s -applyxfm -interp %s', ...
    fsl_path, input_file, options.reference_space, output_file, transform_file, interp_method);

fprintf('    Command: %s\n', cmd);
[status, cmdout] = system(cmd);

if status ~= 0
    error('FSL T1 to DTI transform failed: %s', cmdout);
end

end

function apply_t1_to_mni_fsl(input_file, output_file, registration_data, options, fsl_path)
% Apply T1 to MNI transform using FSL

t1_mni_data = registration_data.transforms.t1_to_mni;

if strcmp(t1_mni_data.type, 'nonlinear') && isfield(t1_mni_data, 'forward_warp')
    % Use nonlinear warp
    cmd = sprintf('%s/bin/applywarp --ref=%s --in=%s --warp=%s --out=%s --interp=%s', ...
        fsl_path, options.reference_space, input_file, t1_mni_data.forward_warp, ...
        output_file, options.interpolation);
else
    % Use linear transform
    interp_method = get_fsl_interpolation(options.interpolation);
    cmd = sprintf('%s/bin/flirt -in %s -ref %s -out %s -init %s -applyxfm -interp %s', ...
        fsl_path, input_file, options.reference_space, output_file, ...
        t1_mni_data.linear_transform_file, interp_method);
end

fprintf('    Command: %s\n', cmd);
[status, cmdout] = system(cmd);

if status ~= 0
    error('FSL T1 to MNI transform failed: %s', cmdout);
end

end

function apply_mni_to_t1_fsl(input_file, output_file, registration_data, options, fsl_path)
% Apply MNI to T1 transform using FSL (inverse of T1 to MNI)

t1_mni_data = registration_data.transforms.t1_to_mni;

if strcmp(t1_mni_data.type, 'nonlinear') && isfield(t1_mni_data, 'inverse_warp')
    % Use inverse nonlinear warp
    cmd = sprintf('%s/bin/applywarp --ref=%s --in=%s --warp=%s --out=%s --interp=%s', ...
        fsl_path, options.reference_space, input_file, t1_mni_data.inverse_warp, ...
        output_file, options.interpolation);
else
    % Use inverse linear transform - need to compute it
    error('Linear MNI to T1 inverse transform not implemented yet');
end

fprintf('    Command: %s\n', cmd);
[status, cmdout] = system(cmd);

if status ~= 0
    error('FSL MNI to T1 transform failed: %s', cmdout);
end

end

function interp_method = get_fsl_interpolation(interpolation)
% Convert interpolation method to FSL format

switch lower(interpolation)
    case 'linear'
        interp_method = 'trilinear';
    case 'nearest'
        interp_method = 'nearestneighbour';
    case 'spline'
        interp_method = 'spline';
    otherwise
        warning('Unknown interpolation method: %s, using trilinear', interpolation);
        interp_method = 'trilinear';
end

end

function output_file = apply_transforms_spm(input_file, registration_data, transform_chain, options)
% Apply transformation chain using SPM tools

fprintf('Using SPM for transformation...\n');
error('SPM transformation chain not implemented yet');

end