function main(imgpath, nimpath)
arguments
    % Path to NiFTI image file
    imgpath string
    
    % Path to save processed .mat file (must end in `.mat`)
    nimpath string
end

% include folders to path
addpaths();


% Define the file extensions and suffixes
img_file = imgpath + ".nii.gz";
raw_file = imgpath + "_raw.nii.gz";

% Check if the processed NIfTI file exists
if isfile(img_file)
    fprintf("Found processed NIfTI image: %s\n", img_file);
else
    fprintf("Processed NIfTI image not found: %s\n", img_file);
    
    % Check if the raw data exists
    if isfile(raw_file)
        fprintf("Found raw data: %s\n", raw_file);
        fprintf("Starting preprocessing...\n");
        
        % Preprocess the raw data
        nim_preprocessing(imgpath);
        
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
%nim = nim_eig(nim);
nim = nim_fa(nim);
% Calculate parcellation
nim = nim_parcellation(nim, 'nifti_sample/parcellation_sample/parcellation_mask.nii.gz');
nim_save(nim, nimpath);

end_time = datetime('now', 'Format', 'yyyy-MM-dd hh:mm:ss');
fprintf("HINEC END: %s\n", string(end_time));
end
