function [parcellation_mask_output, atlas_labels_file] = preproc_atlas_resampling_fixed(reference_file, output_dir, file_prefix, atlas_type)
% preproc_atlas_resampling_fixed: Enhanced atlas resampling with better validation
%
% This version addresses common atlas resampling issues:
% 1. Better atlas file validation
% 2. Improved resampling with label preservation
% 3. Enhanced error checking and reporting
% 4. Atlas-specific optimizations

arguments
    reference_file  % Reference DWI file for dimensions
    output_dir     % Directory for output files
    file_prefix    % Prefix for output files
    atlas_type     % Atlas type ('HarvardOxford', 'JHU', 'JHU-tract')
end

fprintf('Step 4: Enhanced atlas resampling and validation...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output paths
resampled_parcellation_mask = fullfile(output_dir, 'resampled_parcellation_mask.nii.gz');
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');
atlas_labels_file = file_prefix + "_atlas_labels.mat";

% Get and validate atlas path
atlas_path = get_atlas_path_validated(fsl_path, atlas_type);
fprintf('Using validated atlas: %s\n', atlas_path);

% Check reference file exists
if ~isfile(reference_file)
    error('Reference DWI file not found: %s', reference_file);
end

% Get reference dimensions for validation
try
    if exist('niftiinfo', 'file')
        ref_info = niftiinfo(reference_file);
        ref_dims = ref_info.ImageSize(1:3);
    else
        ref_V = spm_vol(reference_file);
        ref_dims = ref_V.dim;
    end
    fprintf('Reference DWI dimensions: %d x %d x %d\n', ref_dims);
catch ME
    warning('Could not determine reference dimensions: %s', ME.message);
    ref_dims = [];
end

% Load original atlas to check its properties
fprintf('Analyzing original atlas properties...\n');
try
    if exist('niftiread', 'file')
        orig_atlas_data = niftiread(atlas_path);
        orig_atlas_info = niftiinfo(atlas_path);
    else
        orig_V = spm_vol(atlas_path);
        orig_atlas_data = spm_read_vols(orig_V);
    end
    
    orig_unique_labels = unique(orig_atlas_data(:));
    orig_unique_labels = orig_unique_labels(orig_unique_labels > 0);
    
    fprintf('Original atlas statistics:\n');
    fprintf('  Dimensions: %d x %d x %d\n', size(orig_atlas_data));
    fprintf('  Unique labels: %d\n', length(orig_unique_labels));
    fprintf('  Label range: %d to %d\n', min(orig_unique_labels), max(orig_unique_labels));
    
    % Validate atlas labels are reasonable
    validate_original_atlas_labels(atlas_type, orig_unique_labels);
    
catch ME
    warning('Could not analyze original atlas: %s', ME.message);
    orig_unique_labels = [];
end

% Perform enhanced resampling
fprintf('Performing enhanced atlas resampling...\n');
success = perform_enhanced_resampling(fsl_path, atlas_path, reference_file, resampled_parcellation_mask, atlas_type);

if ~success
    error('Atlas resampling failed');
end

% Validate resampled atlas
fprintf('Validating resampled atlas...\n');
validate_resampled_atlas(resampled_parcellation_mask, orig_unique_labels, ref_dims, atlas_type);

% Copy resampled atlas to final location
try
    copyfile(resampled_parcellation_mask, parcellation_mask_output);
    fprintf('✓ Parcellation mask created: %s\n', parcellation_mask_output);
catch ME
    error('Failed to copy parcellation mask: %s', ME.message);
end

% Load and save atlas labels with enhanced validation
fprintf('Loading and validating atlas labels...\n');
atlas_labels_file = load_and_validate_atlas_labels(atlas_type, file_prefix, atlas_labels_file);

% Final validation report
generate_final_validation_report(parcellation_mask_output, atlas_labels_file, atlas_type);

end

function atlas_path = get_atlas_path_validated(fsl_path, atlas_type)
% Get atlas path with enhanced validation
fprintf('Validating atlas files for type: %s\n', atlas_type);

switch lower(atlas_type)
    case 'jhu'
        potential_paths = {
            fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz'),
            fullfile(fsl_path, 'data/atlases/JHU/JHU-WhiteMatter-labels-1mm.nii.gz')
        };
        
    case 'jhu-tract'
        potential_paths = {
            fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz'),
            fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-prob-1mm.nii.gz')
        };
        
    case 'harvardoxford'
        potential_paths = {
            fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'),
            fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cortl-maxprob-thr0-1mm.nii.gz')
        };
        
    otherwise
        error('Unknown atlas type: %s', atlas_type);
end

% Find the first existing atlas
atlas_path = '';
for i = 1:length(potential_paths)
    if isfile(potential_paths{i})
        atlas_path = potential_paths{i};
        fprintf('✓ Found atlas: %s\n', atlas_path);
        break;
    end
end

if isempty(atlas_path)
    fprintf('✗ Atlas files checked:\n');
    for i = 1:length(potential_paths)
        fprintf('  - %s (not found)\n', potential_paths{i});
    end
    error('No atlas file found for type: %s', atlas_type);
end
end

function validate_original_atlas_labels(atlas_type, unique_labels)
% Validate original atlas has reasonable labels
expected_counts = containers.Map({'jhu', 'jhu-tract', 'harvardoxford'}, {48, 20, 96});

if isempty(unique_labels)
    error('Original atlas contains no non-zero labels!');
end

atlas_key = lower(atlas_type);
if expected_counts.isKey(atlas_key)
    expected_count = expected_counts(atlas_key);
    actual_count = length(unique_labels);
    
    if actual_count < 0.5 * expected_count
        warning('Original atlas has fewer labels than expected (%d vs ~%d for %s)', ...
                actual_count, expected_count, atlas_type);
    else
        fprintf('✓ Original atlas label count looks reasonable (%d labels)\n', actual_count);
    end
end
end

function success = perform_enhanced_resampling(fsl_path, atlas_path, reference_file, output_path, atlas_type)
% Enhanced resampling with better error handling
success = false;

% Method 1: Standard FLIRT with nearest neighbor
fprintf('Attempting FLIRT resampling (nearest neighbor)...\n');
cmd_resample = sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyxfm -usesqform -interp nearestneighbour', ...
    fsl_path, atlas_path, reference_file, output_path);

fprintf('Running: %s\n', cmd_resample);
[status, cmdout] = system(cmd_resample);

if status == 0 && isfile(output_path)
    fprintf('✓ FLIRT resampling successful\n');
    success = true;
    return;
else
    fprintf('✗ FLIRT failed: %s\n', cmdout);
end

% Method 2: Try with different interpolation
fprintf('Attempting FLIRT with different options...\n');
cmd_resample2 = sprintf('%s/bin/flirt -in %s -ref %s -out %s -interp nearestneighbour -omat /tmp/identity.mat', ...
    fsl_path, atlas_path, reference_file, output_path);

[status, cmdout] = system(cmd_resample2);
if status == 0 && isfile(output_path)
    fprintf('✓ Alternative FLIRT successful\n');
    success = true;
    return;
else
    fprintf('✗ Alternative FLIRT failed: %s\n', cmdout);
end

% Method 3: If FLIRT fails, try manual resampling (last resort)
fprintf('FLIRT failed, attempting manual resampling...\n');
try
    success = manual_atlas_resampling(atlas_path, reference_file, output_path);
catch ME
    fprintf('✗ Manual resampling failed: %s\n', ME.message);
end
end

function success = manual_atlas_resampling(atlas_path, reference_file, output_path)
% Manual resampling as fallback
success = false;

try
    % Load atlas and reference
    if exist('niftiread', 'file')
        atlas_data = niftiread(atlas_path);
        atlas_info = niftiinfo(atlas_path);
        ref_info = niftiinfo(reference_file);
        ref_dims = ref_info.ImageSize(1:3);
    else
        atlas_V = spm_vol(atlas_path);
        atlas_data = spm_read_vols(atlas_V);
        ref_V = spm_vol(reference_file);
        ref_dims = ref_V.dim;
    end
    
    % Simple nearest neighbor resampling
    atlas_dims = size(atlas_data);
    
    % Create coordinate grids
    [X, Y, Z] = meshgrid(1:atlas_dims(2), 1:atlas_dims(1), 1:atlas_dims(3));
    [Xq, Yq, Zq] = meshgrid(linspace(1, atlas_dims(2), ref_dims(2)), ...
                           linspace(1, atlas_dims(1), ref_dims(1)), ...
                           linspace(1, atlas_dims(3), ref_dims(3)));
    
    % Resample using nearest neighbor
    resampled_data = interp3(X, Y, Z, double(atlas_data), Xq, Yq, Zq, 'nearest', 0);
    resampled_data = int32(round(resampled_data));
    
    % Save using SPM
    ref_V = spm_vol(reference_file);
    out_V = ref_V;
    out_V.fname = output_path;
    out_V.dt = [16 0]; % int32
    out_V = spm_write_vol(out_V, resampled_data);
    
    fprintf('✓ Manual resampling completed\n');
    success = true;
    
catch ME
    fprintf('✗ Manual resampling error: %s\n', ME.message);
end
end

function validate_resampled_atlas(resampled_file, orig_unique_labels, ref_dims, atlas_type)
% Validate the resampled atlas
try
    if exist('niftiread', 'file')
        resampled_data = niftiread(resampled_file);
    else
        V = spm_vol(resampled_file);
        resampled_data = spm_read_vols(V);
    end
    
    resampled_unique_labels = unique(resampled_data(:));
    resampled_unique_labels = resampled_unique_labels(resampled_unique_labels > 0);
    
    fprintf('Resampled atlas validation:\n');
    fprintf('  Dimensions: %d x %d x %d\n', size(resampled_data));
    fprintf('  Unique labels: %d\n', length(resampled_unique_labels));
    
    if ~isempty(resampled_unique_labels)
        fprintf('  Label range: %d to %d\n', min(resampled_unique_labels), max(resampled_unique_labels));
    end
    
    % Check if we lost too many labels during resampling
    if ~isempty(orig_unique_labels)
        label_retention = length(resampled_unique_labels) / length(orig_unique_labels);
        fprintf('  Label retention: %.1f%% (%d/%d)\n', ...
                100*label_retention, length(resampled_unique_labels), length(orig_unique_labels));
        
        if label_retention < 0.7
            warning('Significant label loss during resampling (%.1f%% retained)', 100*label_retention);
        else
            fprintf('✓ Good label retention during resampling\n');
        end
    end
    
    % Check dimensions match reference
    if ~isempty(ref_dims) && ~isequal(size(resampled_data), ref_dims)
        warning('Resampled dimensions do not match reference: %dx%dx%d vs %dx%dx%d', ...
                size(resampled_data), ref_dims);
    else
        fprintf('✓ Resampled dimensions match reference\n');
    end
    
catch ME
    warning('Could not validate resampled atlas: %s', ME.message);
end
end

function atlas_labels_file = load_and_validate_atlas_labels(atlas_type, file_prefix, atlas_labels_file)
% Load atlas labels with validation
try
    % Add path to nim_utils if needed
    if ~exist('nim_load_atlas_labels', 'file')
        addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'nim_utils'));
    end
    
    % Load the atlas labels
    atlas_labels = nim_load_atlas_labels(atlas_type);
    
    if atlas_labels.map.Count > 0
        fprintf('✓ Loaded %d atlas labels for %s\n', atlas_labels.map.Count, atlas_type);
        
        % Add enhanced metadata
        atlas_labels.atlas_type = atlas_type;
        atlas_labels.creation_time = datetime('now');
        atlas_labels.file_prefix = char(file_prefix);
        
        % For JHU atlases, add specific variant info
        if strcmpi(atlas_type, 'JHU') || strcmpi(atlas_type, 'JHU-tract')
            atlas_labels.atlas_variant = atlas_type;
        end
        
        % Save the atlas labels
        save(atlas_labels_file, 'atlas_labels');
        
        % Validate saved file
        if isfile(atlas_labels_file)
            file_info = dir(atlas_labels_file);
            fprintf('✓ Atlas labels saved: %s (%.1f KB)\n', atlas_labels_file, file_info.bytes/1024);
        else
            error('Failed to save atlas labels file');
        end
        
    else
        warning('No atlas labels loaded for %s', atlas_type);
        atlas_labels_file = ''; % Return empty if failed
    end
    
catch ex
    warning('Error loading atlas labels for %s: %s', atlas_type, ex.message);
    atlas_labels_file = ''; % Return empty if failed
end
end

function generate_final_validation_report(parcellation_mask_output, atlas_labels_file, atlas_type)
% Generate final validation report
fprintf('\n=== ATLAS RESAMPLING VALIDATION REPORT ===\n');

% Check parcellation mask
if isfile(parcellation_mask_output)
    try
        if exist('niftiread', 'file')
            mask_data = niftiread(parcellation_mask_output);
        else
            V = spm_vol(parcellation_mask_output);
            mask_data = spm_read_vols(V);
        end
        
        unique_labels = unique(mask_data(:));
        unique_labels = unique_labels(unique_labels > 0);
        
        file_info = dir(parcellation_mask_output);
        fprintf('✓ Parcellation Mask: %s\n', parcellation_mask_output);
        fprintf('  Size: %.1f MB\n', file_info.bytes/1024/1024);
        fprintf('  Dimensions: %d x %d x %d\n', size(mask_data));
        fprintf('  Labels: %d (range: %d-%d)\n', length(unique_labels), min(unique_labels), max(unique_labels));
        
    catch ME
        fprintf('✗ Error validating parcellation mask: %s\n', ME.message);
    end
else
    fprintf('✗ Parcellation mask file missing\n');
end

% Check atlas labels
if isfile(atlas_labels_file)
    try
        labels_data = load(atlas_labels_file);
        if isfield(labels_data, 'atlas_labels')
            file_info = dir(atlas_labels_file);
            fprintf('✓ Atlas Labels: %s\n', atlas_labels_file);
            fprintf('  Size: %.1f KB\n', file_info.bytes/1024);
            fprintf('  Label count: %d\n', labels_data.atlas_labels.map.Count);
            fprintf('  Atlas type: %s\n', labels_data.atlas_labels.atlas_type);
        else
            fprintf('✗ Invalid atlas labels file structure\n');
        end
    catch ME
        fprintf('✗ Error validating atlas labels: %s\n', ME.message);
    end
else
    fprintf('✗ Atlas labels file missing or not created\n');
end

fprintf('==========================================\n');
end