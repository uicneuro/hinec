function nim_preprocessing(file_prefix)
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

% Construct file paths based on the file prefix and suffixes
dwi_file = file_prefix + suffix_raw;
bvec_file = file_prefix + suffix_bvec;
bval_file = file_prefix + suffix_bval;
brain_mask_file = file_prefix + suffix_brain_mask;
output_file = file_prefix + suffix_processed;
output_dir = fileparts(output_file);

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

% Step 4: Resample the Harvard-Oxford atlas to match the dimensions of the final product
atlas_path = fullfile(fsl_path, 'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz');
resampled_parcellation_mask = fullfile(output_dir, 'resampled_parcellation_mask.nii.gz');

% Use flirt with nearest neighbor interpolation to preserve labels
cmd_resample = sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyxfm -usesqform -interp nearestneighbour', ...
  fsl_path, atlas_path, output_file, resampled_parcellation_mask);
[status, cmdout] = system(cmd_resample);
if status ~= 0
  error('Error in resampling the parcellation atlas with flirt: %s', cmdout);
end

% Step 5: (Optional) Use flirt to apply the resampled parcellation mask to the final product
parcellation_mask_output = fullfile(output_dir, 'parcellation_mask.nii.gz');
copyfile(resampled_parcellation_mask, parcellation_mask_output);

% Display completion message
disp('DWI preprocessing complete.');
end
