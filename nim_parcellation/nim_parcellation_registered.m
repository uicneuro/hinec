function nim = nim_parcellation_registered(nim, registration_data, parcellation_mask_file)
% nim_parcellation_registered: Enhanced parcellation using proper registration
%
% This function performs brain parcellation using proper multi-modal registration
% instead of relying only on image header information. It transforms atlas data
% from MNI space through T1 space to DTI space using computed transformations.
%
% Arguments:
%   nim - NIM structure containing DTI data
%   registration_data - Registration data from nim_registration
%   parcellation_mask_file - Output file for parcellation mask
%
% Returns:
%   nim - Updated nim structure with enhanced parcellation_mask field

fprintf('=== Enhanced Parcellation with Registration ===\n');

if ~isstruct(registration_data) || ~isfield(registration_data, 'transforms')
    error('Invalid registration data. Please run nim_registration first.');
end

% Get atlas information from registration options
if isfield(registration_data, 'options') && isfield(registration_data.options, 'atlas_type')
    atlas_type = registration_data.options.atlas_type;
else
    atlas_type = 'JHU-tract';  % Default
end

fprintf('Using atlas: %s\n', atlas_type);
fprintf('Registration method: %s\n', registration_data.options.registration_method);

% Define transformation chain: MNI -> T1 -> DTI
transform_chain = {'mni_to_t1', 't1_to_dti'};

%% Step 1: Get MNI atlas
fprintf('Step 1: Loading MNI atlas...\n');
mni_atlas_file = get_mni_atlas_file(atlas_type);

if ~isfile(mni_atlas_file)
    error('MNI atlas file not found: %s', mni_atlas_file);
end

fprintf('  Atlas file: %s\n', mni_atlas_file);

%% Step 2: Transform atlas from MNI to DTI space
fprintf('Step 2: Transforming atlas MNI -> T1 -> DTI...\n');

% Set up transformation options
transform_options = struct();
transform_options.interpolation = 'nearest';  % Important for label preservation
transform_options.output_file = parcellation_mask_file;
transform_options.method = registration_data.options.registration_method;

% Apply transformation chain
try
    output_file = nim_apply_transforms(mni_atlas_file, registration_data, transform_chain, transform_options);
    
    if ~isfile(output_file)
        error('Atlas transformation failed - output file not created');
    end
    
    fprintf('  ✓ Atlas transformed to DTI space: %s\n', output_file);
    
catch ME
    error('Atlas transformation failed: %s', ME.message);
end

%% Step 3: Load and validate transformed atlas
fprintf('Step 3: Loading transformed atlas...\n');

try
    % Load the transformed parcellation mask
    V_parcel = spm_vol(parcellation_mask_file);
    parcellation_data = spm_read_vols(V_parcel);
    
    % Validate dimensions match DTI data
    if ~isequal(size(parcellation_data), size(nim.FA))
        warning('Parcellation dimensions do not match DTI data');
        fprintf('  DTI dimensions: %s\n', mat2str(size(nim.FA)));
        fprintf('  Parcellation dimensions: %s\n', mat2str(size(parcellation_data)));
        
        % Attempt to resize parcellation to match DTI
        fprintf('  Attempting to resize parcellation...\n');
        parcellation_data = resize_parcellation_to_dti(parcellation_data, nim.FA);
    end
    
    % Store in nim structure
    nim.parcellation_mask = parcellation_data;
    
    % Get unique labels and statistics
    unique_labels = unique(parcellation_data(:));
    num_labels = length(unique_labels) - 1;  % Exclude background (0)
    num_voxels = sum(parcellation_data(:) > 0);
    
    fprintf('  ✓ Parcellation loaded successfully\n');
    fprintf('  Number of regions: %d\n', num_labels);
    fprintf('  Number of parcellated voxels: %d\n', num_voxels);
    fprintf('  Coverage: %.1f%% of brain\n', 100 * num_voxels / sum(nim.FA(:) > 0.1));
    
catch ME
    error('Failed to load transformed atlas: %s', ME.message);
end

%% Step 4: Quality validation
fprintf('Step 4: Validating parcellation quality...\n');

% Check for reasonable overlap with high FA regions
high_fa_mask = nim.FA > 0.3;  % White matter regions
wm_labels = parcellation_data(high_fa_mask);
wm_coverage = sum(wm_labels > 0) / sum(high_fa_mask(:));

fprintf('  White matter coverage: %.1f%%\n', wm_coverage * 100);

if wm_coverage > 0.7
    fprintf('  ✓ Good white matter coverage\n');
elseif wm_coverage > 0.5
    fprintf('  ⚠ Moderate white matter coverage\n');
else
    fprintf('  ❌ Poor white matter coverage - check registration quality\n');
end

% Check for anatomically reasonable distribution
label_volumes = histcounts(parcellation_data(parcellation_data > 0), unique_labels(unique_labels > 0));
mean_volume = mean(label_volumes);
std_volume = std(label_volumes);

fprintf('  Mean region volume: %.0f voxels\n', mean_volume);
fprintf('  Volume variability: %.1f%%\n', 100 * std_volume / mean_volume);

%% Step 5: Store registration information in nim structure
fprintf('Step 5: Storing registration metadata...\n');

% Store enhanced parcellation metadata
nim.parcellation_info = struct();
nim.parcellation_info.method = 'registration_based';
nim.parcellation_info.atlas_type = atlas_type;
nim.parcellation_info.transform_chain = transform_chain;
nim.parcellation_info.registration_quality = [];

% Include registration quality metrics if available
if isfield(registration_data, 'quality_metrics')
    nim.parcellation_info.registration_quality = registration_data.quality_metrics;
end

% Store transformation information for future use
nim.parcellation_info.transforms = struct();
nim.parcellation_info.transforms.mni_atlas_file = mni_atlas_file;
nim.parcellation_info.transforms.registration_file = registration_data.output_file;

fprintf('  ✓ Registration metadata stored\n');

fprintf('=== Enhanced Parcellation Complete ===\n');
fprintf('Parcellation mask: %s\n', parcellation_mask_file);

% Generate quality report
generate_parcellation_quality_report(nim, registration_data, parcellation_mask_file);

end

function mni_atlas_file = get_mni_atlas_file(atlas_type)
% Get the path to MNI atlas file based on type

fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSL not found. Please set FSLDIR environment variable.');
end

switch lower(atlas_type)
    case 'jhu'
        mni_atlas_file = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');
        
    case 'jhu-tract'
        mni_atlas_file = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz');
        
    case 'harvardoxford'
        mni_atlas_file = fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz');
        
    otherwise
        warning('Unknown atlas type "%s", defaulting to JHU-tract', atlas_type);
        mni_atlas_file = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz');
end

end

function resized_parcellation = resize_parcellation_to_dti(parcellation_data, fa_data)
% Resize parcellation data to match DTI dimensions using nearest neighbor

fprintf('    Resizing parcellation to match DTI...\n');

dti_size = size(fa_data);
parcel_size = size(parcellation_data);

% Create coordinate grids
[X_parcel, Y_parcel, Z_parcel] = ndgrid(1:parcel_size(1), 1:parcel_size(2), 1:parcel_size(3));
[X_dti, Y_dti, Z_dti] = ndgrid(1:dti_size(1), 1:dti_size(2), 1:dti_size(3));

% Scale coordinates
scale_x = parcel_size(1) / dti_size(1);
scale_y = parcel_size(2) / dti_size(2);
scale_z = parcel_size(3) / dti_size(3);

X_dti_scaled = X_dti * scale_x;
Y_dti_scaled = Y_dti * scale_y;
Z_dti_scaled = Z_dti * scale_z;

% Interpolate using nearest neighbor (preserves integer labels)
resized_parcellation = interp3(X_parcel, Y_parcel, Z_parcel, parcellation_data, ...
    X_dti_scaled, Y_dti_scaled, Z_dti_scaled, 'nearest', 0);

% Convert to integer labels
resized_parcellation = round(resized_parcellation);

fprintf('    ✓ Parcellation resized from %s to %s\n', ...
    mat2str(parcel_size), mat2str(dti_size));

end

function generate_parcellation_quality_report(nim, registration_data, parcellation_mask_file)
% Generate quality report for enhanced parcellation

[output_dir, file_name, ~] = fileparts(parcellation_mask_file);
report_file = fullfile(output_dir, [file_name '_quality_report.txt']);

try
    fid = fopen(report_file, 'w');
    
    fprintf(fid, 'HINEC Enhanced Parcellation Quality Report\n');
    fprintf(fid, '==========================================\n\n');
    fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));
    
    % Atlas information
    fprintf(fid, 'Atlas Information:\n');
    fprintf(fid, '  Type: %s\n', nim.parcellation_info.atlas_type);
    fprintf(fid, '  Method: %s\n', nim.parcellation_info.method);
    fprintf(fid, '  Transform chain: %s\n', strjoin(nim.parcellation_info.transform_chain, ' -> '));
    fprintf(fid, '\n');
    
    % Parcellation statistics
    unique_labels = unique(nim.parcellation_mask(:));
    num_labels = length(unique_labels) - 1;
    num_voxels = sum(nim.parcellation_mask(:) > 0);
    
    fprintf(fid, 'Parcellation Statistics:\n');
    fprintf(fid, '  Number of regions: %d\n', num_labels);
    fprintf(fid, '  Total parcellated voxels: %d\n', num_voxels);
    fprintf(fid, '  Brain coverage: %.1f%%\n', 100 * num_voxels / sum(nim.FA(:) > 0.1));
    fprintf(fid, '\n');
    
    % White matter analysis
    high_fa_mask = nim.FA > 0.3;
    wm_labels = nim.parcellation_mask(high_fa_mask);
    wm_coverage = sum(wm_labels > 0) / sum(high_fa_mask(:));
    
    fprintf(fid, 'White Matter Analysis:\n');
    fprintf(fid, '  WM coverage: %.1f%%\n', wm_coverage * 100);
    if wm_coverage > 0.7
        fprintf(fid, '  Quality: Good\n');
    elseif wm_coverage > 0.5
        fprintf(fid, '  Quality: Moderate\n');
    else
        fprintf(fid, '  Quality: Poor\n');
    end
    fprintf(fid, '\n');
    
    % Registration quality
    if isfield(nim.parcellation_info, 'registration_quality') && ...
       ~isempty(nim.parcellation_info.registration_quality)
        
        fprintf(fid, 'Registration Quality:\n');
        qual = nim.parcellation_info.registration_quality;
        
        if isfield(qual, 'dti_t1_nmi')
            fprintf(fid, '  DTI->T1 NMI: %.4f\n', qual.dti_t1_nmi);
        end
        
        if isfield(qual, 't1_mni_nmi')
            fprintf(fid, '  T1->MNI NMI: %.4f\n', qual.t1_mni_nmi);
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, 'Files:\n');
    fprintf(fid, '  Parcellation mask: %s\n', parcellation_mask_file);
    fprintf(fid, '  Registration data: %s\n', registration_data.output_file);
    
    fclose(fid);
    
    fprintf('  ✓ Quality report: %s\n', report_file);
    
catch ME
    warning('Failed to generate quality report: %s', ME.message);
    if exist('fid', 'var') && fid > 0
        fclose(fid);
    end
end

end