function [parcellation_mask_output, atlas_labels_file] = preproc_atlas_resampling(reference_file, output_dir, file_prefix, atlas_type)
% preproc_atlas_resampling: Resample atlas to match DWI dimensions and load labels
%
% Arguments:
%   reference_file - Reference DWI file for dimensions
%   output_dir - Directory for output files
%   file_prefix - Prefix for output files
%   atlas_type - Atlas type ('HarvardOxford', 'JHU', 'JHU-tract')
%
% Returns:
%   parcellation_mask_output - Path to final parcellation mask
%   atlas_labels_file - Path to atlas labels file

fprintf('Step 4: Resampling atlas and loading labels...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Define output paths
resampled_parcellation_mask = fullfile(output_dir, 'resampled_parcellation_mask.nii.gz');
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');
atlas_labels_file = file_prefix + "_atlas_labels.mat";

% Select atlas based on type
atlas_path = get_atlas_path(fsl_path, atlas_type);
fprintf('Using atlas: %s\n', atlas_path);

% Verify atlas file exists
if ~isfile(atlas_path)
    error('Atlas file not found: %s', atlas_path);
end

% Resample atlas to match DWI dimensions using FLIRT
cmd_resample = sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyxfm -usesqform -interp nearestneighbour', ...
    fsl_path, atlas_path, reference_file, resampled_parcellation_mask);

fprintf('Running: %s\n', cmd_resample);
[status, cmdout] = system(cmd_resample);

if status ~= 0
    error('Error in resampling the parcellation atlas with flirt: %s', cmdout);
end

% Verify resampled atlas was created
if ~isfile(resampled_parcellation_mask)
    error('Atlas resampling failed: output file not found at %s', resampled_parcellation_mask);
end

% Copy resampled atlas to final location
try
    copyfile(resampled_parcellation_mask, parcellation_mask_output);
    fprintf('✓ Parcellation mask created: %s\n', parcellation_mask_output);
catch ME
    error('Failed to copy parcellation mask: %s', ME.message);
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

% Get file size for reporting
file_info = dir(parcellation_mask_output);
fprintf('✓ Final parcellation mask: %s (%.1f MB)\n', parcellation_mask_output, file_info.bytes/1024/1024);

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