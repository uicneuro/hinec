function nim = nim_parcellation(nim, parcellation_file, varargin)
% nim_parcellation: Load or generate parcellation mask for DTI analysis
%
% If parcellation file exists, loads it. If not, generates it from MNI atlas.
%
% Arguments:
%   nim - The nim structure containing DTI data
%   parcellation_file - Path to the parcellation NIfTI file
%   options - Optional structure with:
%       .dwi_file - DWI file path (required for generation)
%       .atlas_type - 'jhu', 'HarvardOxford', 'aal' (default: 'jhu')
%       .force_regenerate - Regenerate even if file exists (default: false)
%
% Returns:
%   nim - Updated nim structure with parcellation_mask field
%
% Examples:
%   % Load existing parcellation
%   nim = nim_parcellation(nim, 'parcellation_mask.nii.gz');
%
%   % Generate if missing
%   opts.dwi_file = 'sub-MR256.nii.gz';
%   nim = nim_parcellation(nim, 'parcellation_mask.nii.gz', opts);

% Parse options
if nargin > 2 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

if ~isfield(options, 'atlas_type')
    options.atlas_type = 'jhu';
end
if ~isfield(options, 'force_regenerate')
    options.force_regenerate = false;
end

%% Generate parcellation if it doesn't exist
if ~isfile(parcellation_file) || options.force_regenerate
    if ~isfield(options, 'dwi_file') || isempty(options.dwi_file)
        error('Parcellation file not found and no dwi_file provided for generation: %s', parcellation_file);
    end

    fprintf("Parcellation file not found. Generating from MNI atlas...\n");
    generate_parcellation_from_mni(options.dwi_file, parcellation_file, options.atlas_type);
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

%% Helper Functions

function generate_parcellation_from_mni(dwi_file, parcellation_file, atlas_type)
% Generate parcellation by registering MNI atlas to DWI space
    fprintf('=== Generating Parcellation from MNI Atlas ===\n');
    fprintf('Atlas: %s\n', atlas_type);

    % Validate FSL
    fsl_path = getenv('FSLDIR');
    if isempty(fsl_path)
        error('FSLDIR not set. Please ensure FSL is installed.');
    end

    [output_dir, ~, ~] = fileparts(parcellation_file);
    if isempty(output_dir)
        output_dir = pwd;
    end

    % Setup paths
    b0_file = fullfile(output_dir, 'temp_b0_for_registration.nii.gz');
    mni_template = fullfile(fsl_path, 'data', 'standard', 'MNI152_T1_1mm.nii.gz');
    [atlas_file, labels_xml] = select_atlas(atlas_type, fsl_path);

    fprintf('Step 1: Extracting b0 reference...\n');
    extract_b0_for_registration(dwi_file, b0_file, fsl_path);

    fprintf('Step 2: Registering MNI → DWI (2-5 min)...\n');
    transform_mat = register_mni_to_dwi(mni_template, b0_file, output_dir, fsl_path);

    fprintf('Step 3: Transforming atlas to DWI space...\n');
    apply_transform_to_atlas(atlas_file, b0_file, parcellation_file, transform_mat, fsl_path);

    % Copy labels if available
    if ~isempty(labels_xml) && isfile(labels_xml)
        labels_file = fullfile(output_dir, 'atlas_labels.xml');
        copyfile(labels_xml, labels_file);
        fprintf('✓ Labels file: %s\n', labels_file);
    end

    % Cleanup
    cleanup_temp_files({b0_file, transform_mat});
    fprintf('✓ Parcellation generated: %s\n', parcellation_file);
end

function [atlas_file, labels_xml] = select_atlas(atlas_type, fsl_path)
    switch lower(atlas_type)
        case 'jhu'
            atlas_file = fullfile(fsl_path, 'data', 'atlases', 'JHU', 'JHU-ICBM-labels-1mm.nii.gz');
            labels_xml = fullfile(fsl_path, 'data', 'atlases', 'JHU-labels.xml');
        case 'harvardoxford'
            atlas_file = fullfile(fsl_path, 'data', 'atlases', 'HarvardOxford', ...
                'HarvardOxford-cort-maxprob-thr25-1mm.nii.gz');
            labels_xml = fullfile(fsl_path, 'data', 'atlases', 'HarvardOxford-Cortical.xml');
        case 'aal'
            atlas_file = fullfile(fsl_path, 'data', 'atlases', 'AAL', 'AAL.nii.gz');
            labels_xml = fullfile(fsl_path, 'data', 'atlases', 'AAL.xml');
        otherwise
            error('Unsupported atlas: %s', atlas_type);
    end
    if ~isfile(atlas_file)
        error('Atlas file not found: %s', atlas_file);
    end
end

function extract_b0_for_registration(dwi_file, b0_file, fsl_path)
    bval_file = strrep(dwi_file, '.nii.gz', '.bval');
    if isfile(bval_file)
        bval_data = load(bval_file);
        b0_indices = find(bval_data < 50);
        if isempty(b0_indices), b0_indices = 1; end
        idx = b0_indices(1) - 1;
        cmd = sprintf('%s/bin/fslroi %s %s %d 1', fsl_path, dwi_file, b0_file, idx);
    else
        warning('No .bval file. Using first volume.');
        cmd = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, b0_file);
    end
    [status, cmdout] = system(cmd);
    if status ~= 0, error('Failed to extract b0: %s', cmdout); end
end

function transform_mat = register_mni_to_dwi(mni_template, b0_file, output_dir, fsl_path)
    transform_mat = fullfile(output_dir, 'mni_to_dwi.mat');
    cmd = sprintf(['%s/bin/flirt -in %s -ref %s -omat %s ', ...
                   '-bins 256 -cost corratio ', ...
                   '-searchrx -90 90 -searchry -90 90 -searchrz -90 90 ', ...
                   '-dof 12 -interp trilinear'], ...
                  fsl_path, mni_template, b0_file, transform_mat);
    [status, cmdout] = system(cmd);
    if status ~= 0, error('Registration failed: %s', cmdout); end
end

function apply_transform_to_atlas(atlas_file, b0_file, parcellation_file, transform_mat, fsl_path)
    cmd = sprintf('%s/bin/flirt -in %s -ref %s -out %s -applyxfm -init %s -interp nearestneighbour', ...
                  fsl_path, atlas_file, b0_file, parcellation_file, transform_mat);
    [status, cmdout] = system(cmd);
    if status ~= 0, error('Atlas transformation failed: %s', cmdout); end
end

function cleanup_temp_files(file_list)
    for i = 1:length(file_list)
        if isfile(file_list{i}), delete(file_list{i}); end
    end
end