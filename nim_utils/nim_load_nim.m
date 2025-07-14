function nim = nim_load_nim(file_prefix)
% nim_load_nim: Load NIM structure with data, mask and atlas labels
%
% Arguments:
%   file_prefix - The prefix for the file paths
%
% Returns:
%   nim - Structure containing all necessary data for analysis and visualization

% Define the suffixes for different files
suffix_processed = '.nii.gz';
suffix_brain_mask = '_M.nii.gz';
suffix_atlas_labels = '_atlas_labels.mat';

% Construct file paths
nii_file = file_prefix + suffix_processed;
mask_file = file_prefix + suffix_brain_mask;
parcellation_mask_file = fullfile(fileparts(file_prefix), 'parcellation_mask.nii.gz');
atlas_labels_file = file_prefix + suffix_atlas_labels;

% Load the main data
if ~exist(nii_file, 'file')
  error('NII file not found: %s', nii_file);
end
nim = load_nii(nii_file);

% Load the brain mask
if exist(mask_file, 'file')
  mask_nii = load_nii(mask_file);
  nim.mask = mask_nii.img;
else
  warning('Brain mask file not found: %s. Creating default mask.', mask_file);
  nim.mask = ones(nim.hdr.dime.dim(2:4));
end

% Load the parcellation mask
if exist(parcellation_mask_file, 'file')
  parcellation_nii = load_nii(parcellation_mask_file);
  nim.parcellation_mask = parcellation_nii.img;
else
  warning('Parcellation mask file not found: %s', parcellation_mask_file);
  nim.parcellation_mask = zeros(nim.hdr.dime.dim(2:4));
end

% Load atlas labels if available
if exist(atlas_labels_file, 'file')
  try
    labels_data = load(atlas_labels_file);
    nim.atlas_labels = labels_data.atlas_labels;
    disp('Atlas labels loaded successfully.');
  catch ex
    warning(ex.identifier, 'Error loading atlas labels: %s', ex.message);
  end
else
  warning('Atlas labels file not found: %s', atlas_labels_file);
end

disp('NIM data loaded successfully.');
end
