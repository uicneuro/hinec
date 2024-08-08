function nim = nim_parcellation(nim, parcellation_file)
arguments
  nim
  parcellation_file {mustBeFile} % Path to the parcellation NIfTI file
end

disp("Loading parcellation data...");

% Read the parcellation NIfTI file
parcellation_nii = load_nii(parcellation_file);
parcellation_data = parcellation_nii.img;

% Ensure the parcellation data is in integer format
parcellation_data = round(parcellation_data); % Round if needed, then convert to integers
parcellation_data = int32(parcellation_data); % Convert to 32-bit integers

% Get the dimensions of the nim data
Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

% Check if the dimensions match
if ~isequal(size(parcellation_data), [Nvox_x, Nvox_y, Nvox_z])
  error('Parcellation dimensions do not match nim data dimensions');
end

% Add the parcellation data to the nim struct
nim.parcellation_mask = parcellation_data;

disp("Parcellation data successfully added to nim structure.");
end
