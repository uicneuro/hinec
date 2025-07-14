function nim_preprocessing(file_prefix, run_eddy, atlas_type)
% nim_preprocessing: Preprocess DWI data with optional eddy current correction
%
% Arguments:
% file_prefix - The prefix for the file paths
% run_eddy - Boolean flag to indicate whether to run eddy current correction (true/false)
% atlas_type - Optional parameter to specify atlas type ('HarvardOxford' or 'JHU')
%              Default is 'HarvardOxford'

if nargin < 3
  atlas_type = 'HarvardOxford';
end

% Ensure the necessary FSL environment variables are set
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
  error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

% Source the FSL environment to ensure dependencies like fslpython are available
system('source $FSLDIR/etc/fslconf/fsl.sh');

% Define the suffixes for different files
suffix_raw = '_raw.nii.gz';
suffix_bvec = '.bvec';
suffix_bval = '.bval';
suffix_brain_mask = '_M.nii.gz';  % This is the final name for the mask file
suffix_processed = '.nii.gz';
suffix_eddy_corrected = '_eddy_corrected.nii.gz';
suffix_acqp = '_acqp.txt';
suffix_index = '_index.txt';
suffix_atlas_labels = '_atlas_labels.mat';

% Construct file paths based on the file prefix and suffixes
dwi_file = file_prefix + suffix_raw;
bvec_file = file_prefix + suffix_bvec;
bval_file = file_prefix + suffix_bval;
brain_mask_file = file_prefix + suffix_brain_mask;
output_file = file_prefix + suffix_processed;
output_dir = fileparts(output_file);
atlas_labels_file = file_prefix + suffix_atlas_labels;

% Intermediate file paths
b0_file = fullfile(output_dir, 'b0.nii.gz');
bet_mask_file = fullfile(output_dir, 'nodif_brain_mask.nii.gz');
eddy_corrected_file = file_prefix + suffix_eddy_corrected;
dti_output_prefix = fullfile(output_dir, 'dti');
acqp_file = file_prefix + suffix_acqp;
index_file = file_prefix + suffix_index;

% Step 1: Extract b0 image (first volume) from DWI data
cmd_extract_b0 = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, b0_file);
[status, cmdout] = system(cmd_extract_b0);
if status ~= 0
  error('Error in fslroi: %s', cmdout);
end

% Step 2: Create brain mask using FSL's BET tool
cmd_bet = sprintf('%s/bin/bet %s %s -m', fsl_path, b0_file, fullfile(output_dir, 'nodif_brain'));
[status, cmdout] = system(cmd_bet);
if status ~= 0
  error('Error in bet: %s', cmdout);
end

% Step 2.5: Rename the mask file to the final mask filename
if isfile(bet_mask_file)
  movefile(bet_mask_file, brain_mask_file);
else
  error('BET mask file not found: %s', bet_mask_file);
end

% Step 3: Create a copy of the raw file and rename it as the final product
copyfile(dwi_file, output_file);

% Step 3.5: Eddy Current Correction
% if run_eddy
% Use FSL's eddy tool with the necessary parameters
% cmd_eddy = sprintf('%s/bin/eddy --imain=%s --mask=%s --bvecs=%s --bvals=%s --out=%s --acqp=%s --index=%s', ...
%   fsl_path, dwi_file, brain_mask_file, bvec_file, bval_file, eddy_corrected_file, acqp_file, index_file);
% [status, cmdout] = system(cmd_eddy);
% if status ~= 0
%   error('Error in eddy: %s', cmdout);
% end
% end

% Step 4: Resample the selected atlas to match the dimensions of the final product
if strcmpi(atlas_type, 'JHU')
  % JHU White-Matter atlas
  atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');
  % You can choose one of the following JHU atlases:
  % atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');
  % atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz');
  % atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-WhiteMatter-labels-1mm.nii.gz');
elseif strcmpi(atlas_type, 'JHU-tract')
  % JHU White-Matter Tract atlas
  atlas_path = fullfile(fsl_path, 'data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz');
else
  % Default: Harvard-Oxford atlas
  atlas_path = fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz');
end
resampled_parcellation_mask = fullfile(output_dir, 'resampled_parcellation_mask.nii.gz');

% Use flirt with nearest neighbor interpolation to preserve labels
cmd_resample = sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyxfm -usesqform -interp nearestneighbour', ...
  fsl_path, atlas_path, output_file, resampled_parcellation_mask);
[status, cmdout] = system(cmd_resample);
if status ~= 0
  error('Error in resampling the parcellation atlas with flirt: %s', cmdout);
end

% Step 5: (Optional) Use flirt to apply the resampled parcellation mask to the final product
% Update output filename to include atlas information
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');
copyfile(resampled_parcellation_mask, parcellation_mask_output);

% Step 6: Load atlas labels and save them
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
  
  % For JHU atlas, also add which specific JHU atlas was used
  if strcmpi(atlas_type, 'JHU')
    [~, atlas_name, ~] = fileparts(atlas_path);
    atlas_labels.atlas_variant = atlas_name;
  elseif strcmpi(atlas_type, 'JHU-tract')
    [~, atlas_name, ~] = fileparts(atlas_path);
    atlas_labels.atlas_variant = atlas_name;
  end
  
  % Save the atlas labels
  save(atlas_labels_file, 'atlas_labels');
  
  disp(['Atlas labels saved to: ' atlas_labels_file]);
catch ex
  warning(ex.identifier, 'Error loading atlas labels: %s', ex.message);
end

% Display completion message
disp('DWI preprocessing complete.');
end
