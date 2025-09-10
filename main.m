function main(imgpath, nimpath)
arguments
    % Path to NiFTI image file
    imgpath string
    
    % Path to save processed .mat file (must end in `.mat`)
    nimpath string
end

% include folders to path
addpath('nim_preprocessing/');
addpath('nim_plots');
addpath('nim_utils');
addpath('nim_calculation');
addpath('nim_parcellation');
addpath('nim_tractography');
addpath('nim_tests');  % Add diagnostic and test functions
addpath(genpath('spm12'));
addpath('utils');
addpath('nifti_sample');
addpath('bfgs');


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
        
        % Set up preprocessing options for comprehensive processing
        preproc_options = struct();
        preproc_options.run_denoising = true;
        preproc_options.denoise_method = 'dwidenoise';  % Use the best denoising method
        preproc_options.run_motion_correction = true;
        preproc_options.run_eddy = true;  % Enable eddy correction for better quality
        preproc_options.improve_mask = true;  % Critical for good tractography
        preproc_options.atlas_type = 'JHU-tract';
        
        % Preprocess the raw data with enhanced pipeline
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


% load NiFTI image
nim = nim_read(imgpath);
% Calculate diffusion tensor
nim = nim_dt_spd(nim);
nim = nim_eig(nim);
nim = nim_fa(nim);
% Calculate parcellation using the fixed loader
[output_dir, ~, ~] = fileparts(imgpath);
parcellation_mask_file = fullfile(output_dir, 'parcellation_mask.nii.gz');
nim = nim_parcellation(nim, parcellation_mask_file);
% Store parcellation mask file path for reference
nim.parcellation_mask_file = parcellation_mask_file;
% Load parcellation labels
nim = nim_load_labels(nim);

% Brain mask improvement using FA data (final step)
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

nim_save(nim, nimpath);

end_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC END: %s\n", string(end_time));
end
