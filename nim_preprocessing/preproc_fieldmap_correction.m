function corrected_file = preproc_fieldmap_correction(dwi_file, fieldmap_file, output_prefix, varargin)
% PREPROC_FIELDMAP_CORRECTION: Apply susceptibility distortion correction using field maps
%
% Usage:
%   corrected_file = preproc_fieldmap_correction(dwi_file, fieldmap_file, output_prefix)
%   corrected_file = preproc_fieldmap_correction(..., 'TE', echo_time, 'dwell_time', dt)
%
% Arguments:
%   dwi_file      - Input DWI file to correct
%   fieldmap_file - Field map file (Hz or rad/s)
%   output_prefix - Prefix for output files
%
% Optional Parameters:
%   'TE'          - Echo time difference for field map (seconds) [default: auto-detect]
%   'dwell_time'  - Effective echo spacing (seconds) [default: 0.00058]
%   'phase_dir'   - Phase encoding direction: 'x', 'y', 'z', 'x-', 'y-', 'z-' [default: 'y']
%   'units'       - Field map units: 'Hz' or 'rad/s' [default: auto-detect from filename]
%   'smooth_sigma'- Gaussian smoothing sigma for field map (mm) [default: 2.0]
%   'mask_file'   - Brain mask file [default: auto-generate from DWI]
%
% Output:
%   corrected_file - Path to distortion-corrected DWI file

%% Parse input arguments
p = inputParser;
addRequired(p, 'dwi_file', @(x) ischar(x) || isstring(x));
addRequired(p, 'fieldmap_file', @(x) ischar(x) || isstring(x));
addRequired(p, 'output_prefix', @(x) ischar(x) || isstring(x));
addParameter(p, 'TE', [], @isnumeric);
addParameter(p, 'dwell_time', 0.00058, @isnumeric);
addParameter(p, 'phase_dir', 'y', @(x) ismember(x, {'x', 'y', 'z', 'x-', 'y-', 'z-'}));
addParameter(p, 'units', 'auto', @(x) ismember(x, {'Hz', 'rad/s', 'auto'}));
addParameter(p, 'smooth_sigma', 2.0, @isnumeric);
addParameter(p, 'mask_file', '', @(x) ischar(x) || isstring(x));

parse(p, dwi_file, fieldmap_file, output_prefix, varargin{:});
opts = p.Results;

fprintf('=== Field Map Distortion Correction ===\n');
fprintf('DWI file: %s\n', dwi_file);
fprintf('Field map: %s\n', fieldmap_file);
fprintf('Output prefix: %s\n', output_prefix);

%% Check FSL installation
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end

%% Verify input files
if ~exist(dwi_file, 'file')
    error('DWI file not found: %s', dwi_file);
end
if ~exist(fieldmap_file, 'file')
    error('Field map file not found: %s', fieldmap_file);
end

%% Auto-detect field map units from filename
if strcmp(opts.units, 'auto')
    [~, fname, ~] = fileparts(fieldmap_file);
    if contains(lower(fname), 'hz')
        opts.units = 'Hz';
        fprintf('Auto-detected field map units: Hz\n');
    elseif contains(lower(fname), 'rad')
        opts.units = 'rad/s';
        fprintf('Auto-detected field map units: rad/s\n');
    else
        opts.units = 'Hz';  % Default assumption
        fprintf('Could not auto-detect units, assuming Hz\n');
    end
end

%% Generate output filenames
corrected_file = [output_prefix '_fieldmap_corrected.nii.gz'];
fmap_processed = [output_prefix '_fieldmap_processed.nii.gz'];
fmap_unwrapped = [output_prefix '_fieldmap_unwrapped.nii.gz'];

fprintf('Field map units: %s\n', opts.units);
fprintf('Phase encoding direction: %s\n', opts.phase_dir);
fprintf('Dwell time: %.6f s\n', opts.dwell_time);
fprintf('Smoothing sigma: %.1f mm\n', opts.smooth_sigma);

%% Create or use brain mask
if isempty(opts.mask_file) || ~exist(opts.mask_file, 'file')
    fprintf('Creating brain mask from DWI data...\n');
    temp_b0 = [output_prefix '_temp_b0.nii.gz'];
    temp_mask = [output_prefix '_temp_mask.nii.gz'];

    % Extract first volume (typically b=0)
    cmd = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, temp_b0);
    [status, result] = system(cmd);
    if status ~= 0
        error('Failed to extract b0 volume: %s', result);
    end

    % Create brain mask
    cmd = sprintf('%s/bin/bet %s %s -m -f 0.3', fsl_path, temp_b0, strrep(temp_mask, '.nii.gz', ''));
    [status, result] = system(cmd);
    if status ~= 0
        error('Failed to create brain mask: %s', result);
    end

    mask_file = temp_mask;

    % Clean up temporary b0
    if exist(temp_b0, 'file')
        delete(temp_b0);
    end
else
    mask_file = opts.mask_file;
    fprintf('Using provided brain mask: %s\n', mask_file);
end

%% Process field map
fprintf('Processing field map...\n');

% Convert units if needed
if strcmp(opts.units, 'rad/s')
    fprintf('Converting field map from rad/s to Hz...\n');
    cmd = sprintf('%s/bin/fslmaths %s -div 6.28318 %s', fsl_path, fieldmap_file, fmap_processed);
    [status, result] = system(cmd);
    if status ~= 0
        error('Failed to convert field map units: %s', result);
    end
else
    % Copy field map to processing location
    copyfile(fieldmap_file, fmap_processed);
end

% Apply brain mask to field map
fprintf('Masking field map...\n');
cmd = sprintf('%s/bin/fslmaths %s -mas %s %s', fsl_path, fmap_processed, mask_file, fmap_processed);
[status, result] = system(cmd);
if status ~= 0
    error('Failed to mask field map: %s', result);
end

% Smooth field map to reduce noise
if opts.smooth_sigma > 0
    fprintf('Smoothing field map (sigma=%.1f mm)...\n', opts.smooth_sigma);
    cmd = sprintf('%s/bin/fslmaths %s -s %.2f %s', fsl_path, fmap_processed, opts.smooth_sigma, fmap_processed);
    [status, result] = system(cmd);
    if status ~= 0
        error('Failed to smooth field map: %s', result);
    end
end

% Skip phase unwrapping for Hz field maps (they are already unwrapped frequency maps)
fprintf('Skipping phase unwrapping (Hz field maps are already unwrapped)...\n');
copyfile(fmap_processed, fmap_unwrapped);

%% Apply distortion correction using FUGUE
fprintf('Applying susceptibility distortion correction...\n');

% Convert phase encoding direction to FUGUE format
switch opts.phase_dir
    case 'x'
        fugue_dir = 'x';
    case 'x-'
        fugue_dir = 'x-';
    case 'y'
        fugue_dir = 'y';
    case 'y-'
        fugue_dir = 'y-';
    case 'z'
        fugue_dir = 'z';
    case 'z-'
        fugue_dir = 'z-';
    otherwise
        error('Unsupported phase encoding direction: %s', opts.phase_dir);
end

% Apply FUGUE correction
cmd = sprintf('%s/bin/fugue -i %s --loadfmap=%s --dwell=%.6f --unwarpdir=%s -u %s', ...
    fsl_path, dwi_file, fmap_unwrapped, opts.dwell_time, fugue_dir, corrected_file);

fprintf('Running FUGUE command: %s\n', cmd);
[status, result] = system(cmd);

if status ~= 0
    error('FUGUE distortion correction failed: %s', result);
end

%% Verify output
if ~exist(corrected_file, 'file')
    error('Distortion correction failed - output file not created');
end

% Get file size info
file_info = dir(corrected_file);
fprintf('✓ Distortion correction complete\n');
fprintf('  Output file: %s (%.1f MB)\n', corrected_file, file_info.bytes/1024/1024);

%% Clean up temporary files
temp_files = {fmap_processed, fmap_unwrapped};
if exist([output_prefix '_temp_mask.nii.gz'], 'file')
    temp_files{end+1} = [output_prefix '_temp_mask.nii.gz'];
end

for i = 1:length(temp_files)
    if exist(temp_files{i}, 'file')
        delete(temp_files{i});
    end
end

fprintf('Field map distortion correction completed successfully!\n');
end