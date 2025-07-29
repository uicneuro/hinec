function nim = nim_parcellation(nim, parcellation_file)
% nim_parcellation: Load and process parcellation mask for DTI analysis
%
% Arguments:
%   nim - The nim structure containing DTI data
%   parcellation_file - Path to the parcellation NIfTI file
%
% Returns:
%   nim - Updated nim structure with parcellation_mask field

arguments
    nim
    parcellation_file {mustBeFile} % Path to the parcellation NIfTI file
end

fprintf("Loading parcellation data from: %s\n", parcellation_file);

%% Load parcellation data using multiple methods
try
    % Method 1: Try MATLAB's built-in niftiread (if available)
    if exist('niftiread', 'file')
        parcellation_data = niftiread(parcellation_file);
        fprintf("Loaded parcellation using niftiread\n");
    else
        % Method 2: Use SPM's nifti class
        parcellation_nii = nifti(parcellation_file);
        parcellation_data = parcellation_nii.dat(:,:,:);
        fprintf("Loaded parcellation using SPM nifti\n");
    end
catch ME
    % Method 3: Fallback to SPM vol functions
    try
        fprintf("Falling back to spm_vol method...\n");
        V = spm_vol(parcellation_file);
        parcellation_data = spm_read_vols(V);
        fprintf("Loaded parcellation using spm_read_vols\n");
    catch ME2
        error('All parcellation loading methods failed. Error: %s', ME.message);
    end
end

%% Process parcellation data
% Ensure integer format
parcellation_data = round(parcellation_data);
parcellation_data = int32(parcellation_data);

% Get expected dimensions from nim data
if isfield(nim, 'FA')
    expected_dims = size(nim.FA);
elseif isfield(nim, 'img')
    expected_dims = size(nim.img);
    expected_dims = expected_dims(1:3); % Take first 3 dimensions
elseif isfield(nim, 'hdr') && isfield(nim.hdr, 'ImageSize')
    expected_dims = nim.hdr.ImageSize(1:3);
else
    error('Cannot determine expected dimensions from nim structure');
end

%% Handle dimension mismatches
parc_dims = size(parcellation_data);
fprintf("Parcellation dimensions: %d x %d x %d\n", parc_dims);
fprintf("Expected dimensions: %d x %d x %d\n", expected_dims);

if ~isequal(parc_dims, expected_dims)
    fprintf("Dimension mismatch detected, attempting to fix...\n");
    
    % Try common dimension fixes
    if isequal(parc_dims, [expected_dims(2), expected_dims(1), expected_dims(3)])
        fprintf("Applying X-Y transpose fix...\n");
        parcellation_data = permute(parcellation_data, [2, 1, 3]);
    elseif isequal(parc_dims, [expected_dims(1), expected_dims(3), expected_dims(2)])
        fprintf("Applying Y-Z swap fix...\n");
        parcellation_data = permute(parcellation_data, [1, 3, 2]);
    else
        % Try resampling if dimensions are close
        dim_ratios = parc_dims ./ expected_dims;
        if all(dim_ratios > 0.8 & dim_ratios < 1.2)
            fprintf("Resampling parcellation to match expected dimensions...\n");
            [X, Y, Z] = meshgrid(1:parc_dims(2), 1:parc_dims(1), 1:parc_dims(3));
            [Xq, Yq, Zq] = meshgrid(linspace(1, parc_dims(2), expected_dims(2)), ...
                                   linspace(1, parc_dims(1), expected_dims(1)), ...
                                   linspace(1, parc_dims(3), expected_dims(3)));
            
            parcellation_data = interp3(X, Y, Z, double(parcellation_data), Xq, Yq, Zq, 'nearest');
            parcellation_data = int32(round(parcellation_data));
        else
            error('Cannot resolve dimension mismatch: parcellation %dx%dx%d vs expected %dx%dx%d', ...
                  parc_dims, expected_dims);
        end
    end
end

%% Final validation
final_dims = size(parcellation_data);
if ~isequal(final_dims, expected_dims)
    error('Final dimension check failed: parcellation %dx%dx%d vs expected %dx%dx%d', ...
          final_dims, expected_dims);
end

%% Add parcellation data to nim structure
nim.parcellation_mask = parcellation_data;

%% Report statistics
unique_labels = unique(parcellation_data(:));
unique_labels = unique_labels(unique_labels > 0); % Exclude background (0)

fprintf("Parcellation statistics:\n");
fprintf("  Final dimensions: %d x %d x %d\n", size(parcellation_data));
fprintf("  Unique labels: %d (excluding background)\n", length(unique_labels));
if ~isempty(unique_labels)
    fprintf("  Label range: %d to %d\n", min(unique_labels), max(unique_labels));
end
fprintf("  Background voxels: %.1f%%\n", 100*sum(parcellation_data(:) == 0)/numel(parcellation_data));

fprintf("Parcellation data successfully loaded and processed.\n");

end