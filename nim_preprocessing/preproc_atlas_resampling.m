function [parcellation_mask_output, atlas_labels_file] = preproc_atlas_resampling(reference_file, output_dir, file_prefix, atlas_type, options)
% preproc_atlas_resampling: Register atlas to DWI space using proper registration chain
%
% Arguments:
%   reference_file - Reference DWI file for dimensions
%   output_dir - Directory for output files
%   file_prefix - Prefix for output files
%   atlas_type - Atlas type ('HarvardOxford', 'JHU', 'JHU-tract')
%   options - Preprocessing options structure
%
% Returns:
%   parcellation_mask_output - Path to final parcellation mask
%   atlas_labels_file - Path to atlas labels file

fprintf('Step 8: Atlas registration to DWI space...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output paths
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');
atlas_labels_file = [file_prefix '_atlas_labels.mat'];

% Select atlas based on type
atlas_path = get_atlas_path(fsl_path, atlas_type);
fprintf('Using atlas: %s\n', atlas_path);

% Verify atlas file exists
if ~isfile(atlas_path)
    error('Atlas file not found: %s', atlas_path);
end

% Check if T1-based registration is available and enabled
if isfield(options, 'use_t1_registration') && options.use_t1_registration && isfield(options, 't1_available') && options.t1_available
    fprintf('Using T1-based atlas registration...\n');
    parcellation_mask_output = register_atlas_with_t1(atlas_path, reference_file, output_dir, file_prefix, options);
else
    fprintf('T1 not available, using direct registration...\n');
    parcellation_mask_output = register_atlas_direct(atlas_path, reference_file, output_dir);
end

% Load and save atlas labels
fprintf('Loading atlas labels for %s...\n', atlas_type);
try
    % Add the path to nim_utils if needed
    if ~exist('nim_load_atlas_labels', 'file')
        addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'nim_utils'));
    end
    
    % Load the atlas labels
    atlas_labels = nim_load_atlas_labels(atlas_type);
    fprintf('Loaded %d parcellation labels\n', atlas_labels.map.Count);
    
    % Add atlas type information to the structure
    atlas_labels.atlas_type = atlas_type;
    
    % For JHU atlases, also add which specific variant was used
    if strcmpi(atlas_type, 'JHU') || strcmpi(atlas_type, 'JHU-tract')
        [~, atlas_name, ~] = fileparts(atlas_path);
        atlas_labels.atlas_variant = atlas_name;
    end
    
    % Save the atlas labels
    save(atlas_labels_file, 'atlas_labels');
    
    % Get file size for reporting
    file_info = dir(atlas_labels_file);
    fprintf('✓ Atlas labels saved: %s (%.1f MB)\n', atlas_labels_file, file_info.bytes/1024/1024);
    
catch ex
    warning(ex.identifier, 'Error loading atlas labels: %s', ex.message);
    atlas_labels_file = ''; % Return empty if failed
end

% Validate atlas registration
validate_atlas_registration(parcellation_mask_output, reference_file);

% Get file size for reporting
file_info = dir(parcellation_mask_output);
fprintf('✓ Final parcellation mask: %s (%.1f MB)\n', parcellation_mask_output, file_info.bytes/1024/1024);

end

function parcellation_mask_output = register_atlas_with_t1(atlas_path, dwi_ref_file, output_dir, file_prefix, options)
% Register atlas using T1-based registration chain: MNI→T1→DWI

fprintf('Performing T1-based atlas registration (MNI→T1→DWI)...\n');

fsl_path = getenv('FSLDIR');

% Step 1: Create T1 brain extraction if not already done
t1_brain_file = fullfile(output_dir, [file_prefix '_T1_brain.nii.gz']);
t1_brain_mask_file = fullfile(output_dir, [file_prefix '_T1_brain_mask.nii.gz']);

if ~isfile(t1_brain_file) || ~isfile(t1_brain_mask_file)
    fprintf('Performing T1 brain extraction...\n');
    [t1_brain_file, t1_brain_mask_file] = preproc_t1_brain_extraction(options.t1_file, output_dir);
end

% Step 2: Create DWI reference if not already done
dwi_ref_final = fullfile(output_dir, [file_prefix '_dwi_ref.nii.gz']);
if ~isfile(dwi_ref_final)
    fprintf('Creating DWI reference volume...\n');
    brain_mask_file = fullfile(output_dir, [file_prefix '_M.nii.gz']);
    dwi_ref_final = preproc_create_dwi_reference(dwi_ref_file, brain_mask_file, output_dir, file_prefix);
end
dwi_ref_file = dwi_ref_final;

% Step 3: T1-DWI registration if not already done
t1_to_dwi_mat = fullfile(output_dir, [file_prefix '_T1_to_dwi.mat']);
if ~isfile(t1_to_dwi_mat)
    fprintf('Performing T1-DWI registration...\n');
    t1_to_dwi_mat = preproc_t1_dwi_registration(dwi_ref_file, options.t1_file, t1_brain_file, output_dir, file_prefix);
end

% Step 4: T1-MNI registration if not already done
mni_to_t1_warp = fullfile(output_dir, [file_prefix '_MNI_to_T1_warp.nii.gz']);
if ~isfile(mni_to_t1_warp)
    fprintf('Performing T1-MNI registration...\n');
    mni_to_t1_warp = preproc_t1_mni_registration(options.t1_file, t1_brain_file, output_dir, file_prefix);
end

% Step 5: Apply composite transformation to atlas
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');

fprintf('Applying composite transformation (MNI→T1→DWI) to atlas...\n');
cmd_applywarp = sprintf('%s/bin/applywarp -i %s -r %s --warp=%s --postmat=%s --interp=nn -o %s', ...
    fsl_path, atlas_path, dwi_ref_file, mni_to_t1_warp, t1_to_dwi_mat, parcellation_mask_output);

fprintf('Running: %s\n', cmd_applywarp);
[status, cmdout] = system(cmd_applywarp);

if status ~= 0
    error('Atlas transformation failed: %s', cmdout);
end

if ~isfile(parcellation_mask_output)
    error('Atlas transformation failed: output file not found at %s', parcellation_mask_output);
end

% Step 6: Copy geometry to ensure consistency
fprintf('Ensuring geometry consistency...\n');
cmd_copygeom = sprintf('%s/bin/fslcpgeom %s %s', fsl_path, dwi_ref_file, parcellation_mask_output);
[status, ~] = system(cmd_copygeom);

if status ~= 0
    warning('Could not copy geometry information');
end

fprintf('✓ T1-based atlas registration completed\n');

end

function parcellation_mask_output = register_atlas_direct(atlas_path, reference_file, output_dir)
% Fallback: Direct atlas registration (improved version of original method)

fprintf('Performing direct atlas registration...\n');

fsl_path = getenv('FSLDIR');
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');

% Use improved FLIRT registration without problematic -usesqform
cmd_resample = sprintf('%s/bin/flirt -in %s -ref %s -out %s -interp nearestneighbour', ...
    fsl_path, atlas_path, reference_file, parcellation_mask_output);

fprintf('Running: %s\n', cmd_resample);
[status, cmdout] = system(cmd_resample);

if status ~= 0
    error('Direct atlas registration failed: %s', cmdout);
end

if ~isfile(parcellation_mask_output)
    error('Direct atlas registration failed: output file not found at %s', parcellation_mask_output);
end

fprintf('✓ Direct atlas registration completed\n');

end

function validate_atlas_registration(atlas_file, reference_file)
% Validate atlas registration quality

fprintf('Validating atlas registration...\n');

fsl_path = getenv('FSLDIR');

try
    % Check label value range
    cmd_range = sprintf('%s/bin/fslstats %s -R', fsl_path, atlas_file);
    [status, range_output] = system(cmd_range);
    
    if status == 0
        range_values = str2num(range_output);
        if length(range_values) >= 2
            min_val = range_values(1);
            max_val = range_values(2);
            
            fprintf('Atlas label range: %.0f to %.0f\n', min_val, max_val);
            
            % Check if values are integers (as expected for labels)
            if abs(min_val - round(min_val)) > 0.01 || abs(max_val - round(max_val)) > 0.01
                warning('Atlas contains non-integer values - may indicate interpolation issues');
            end
            
            % Check if range is reasonable for anatomical atlas
            if max_val > 200
                warning('Atlas has unusually high label values (max: %.0f)', max_val);
            end
            
            if min_val < 0
                warning('Atlas contains negative values');
            end
        end
    end
    
    % Check voxel count
    cmd_voxels = sprintf('%s/bin/fslstats %s -V', fsl_path, atlas_file);
    [status, voxel_output] = system(cmd_voxels);
    
    if status == 0
        voxel_info = str2num(voxel_output);
        if length(voxel_info) >= 1
            nonzero_voxels = voxel_info(1);
            fprintf('Atlas coverage: %d non-zero voxels\n', nonzero_voxels);
            
            if nonzero_voxels < 1000
                warning('Atlas has very few non-zero voxels (%d) - registration may have failed', nonzero_voxels);
            end
        end
    end
    
catch ME
    warning(ME.identifier, 'Atlas validation failed: %s', ME.message);
end

fprintf('✓ Atlas validation completed\n');

end

function atlas_path = get_atlas_path(fsl_path, atlas_type)
% Get the full path to the atlas file based on type
switch lower(atlas_type)
    case 'jhu'
        % JHU White-Matter labels atlas
        atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');
        
    case 'jhu-tract'
        % JHU White-Matter Tract atlas (maxprob)
        atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz');
        
    case 'harvardoxford'
        % Harvard-Oxford cortical atlas (default)
        atlas_path = fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz');
        
    otherwise
        % Default to Harvard-Oxford if unknown type
        warning('Unknown atlas type "%s", defaulting to HarvardOxford', atlas_type);
        atlas_path = fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz');
end
end