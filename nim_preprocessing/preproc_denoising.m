function denoised_file = preproc_denoising(dwi_file, file_prefix, method)
% preproc_denoising: Apply denoising to DWI data using various methods
%
% Arguments:
%   dwi_file - Path to the DWI file
%   file_prefix - Prefix for output files
%   method - Denoising method ('dwidenoise', 'nlmeans', 'gaussian') [default: 'dwidenoise']
%
% Returns:
%   denoised_file - Path to the denoised DWI file

if nargin < 3
    method = 'dwidenoise';
end

fprintf('Step: Denoising using %s method...\n', method);

% Define output file path
denoised_file = [strrep(file_prefix, '_raw', '') '_denoised.nii.gz'];

% Verify input file exists
if ~isfile(dwi_file)
    error('DWI file not found: %s', dwi_file);
end

switch lower(method)
    case 'dwidenoise'
        denoised_file = denoise_mrtrix(dwi_file, denoised_file, file_prefix);
        
    case 'nlmeans'
        denoised_file = denoise_nlmeans(dwi_file, denoised_file);
        
    case 'gaussian'
        denoised_file = denoise_gaussian(dwi_file, denoised_file);
        
    otherwise
        error('Unknown denoising method: %s. Supported methods: dwidenoise, nlmeans, gaussian', method);
end

% Verify the output file was created
if ~isfile(denoised_file)
    error('Denoising failed: output file not found at %s', denoised_file);
end

% Get file size for reporting
file_info = dir(denoised_file);
fprintf('✓ Denoising completed: %s (%.1f MB)\n', denoised_file, file_info.bytes/1024/1024);

end

function denoised_file = denoise_mrtrix(dwi_file, denoised_file, file_prefix)
% Denoising using MRtrix3's dwidenoise (MP-PCA method)

fprintf('Using MRtrix3 dwidenoise (MP-PCA method)...\n');

% Check if MRtrix3 is available
[status, ~] = system('which dwidenoise');
if status ~= 0
    fprintf('⚠ MRtrix3 not found, falling back to Gaussian smoothing\n');
    denoised_file = denoise_gaussian(dwi_file, denoised_file);
    return;
end

% Define noise map output
noise_map_file = [strrep(file_prefix, '_raw', '') '_noise_map.nii.gz'];

% Run MRtrix3 dwidenoise
cmd_denoise = sprintf('dwidenoise %s %s -noise %s -force', ...
                     dwi_file, denoised_file, noise_map_file);

fprintf('Command: %s\n', cmd_denoise);

tic;
[status, cmdout] = system(cmd_denoise);
elapsed_time = toc;

if status ~= 0
    fprintf('⚠ MRtrix3 dwidenoise failed, falling back to FSL-based denoising\n');
    denoised_file = denoise_fsl_susan(dwi_file, denoised_file);
    return;
end

fprintf('  Processing time: %.1f seconds\n', elapsed_time);

% Report on noise characteristics if noise map was created
if isfile(noise_map_file)
    fprintf('✓ Noise map saved: %s\n', noise_map_file);
    
    % Load and analyze noise map
    try
        % Simple noise analysis - load first volume to get statistics
        [status, noise_stats] = system(sprintf('fslstats %s -M -S', noise_map_file));
        if status == 0
            noise_values = str2num(noise_stats);
            if length(noise_values) >= 2
                fprintf('  Mean noise level: %.3f\n', noise_values(1));
                fprintf('  Noise std: %.3f\n', noise_values(2));
            end
        end
    catch
        fprintf('  Noise map created but unable to analyze statistics\n');
    end
end

end

function denoised_file = denoise_fsl_susan(dwi_file, denoised_file)
% Fallback denoising using FSL's SUSAN filter

fprintf('Using FSL SUSAN denoising...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('Neither MRtrix3 nor FSL found for denoising');
end

% Calculate brightness threshold for SUSAN (typically 10-25% of image max)
[status, max_val_str] = system(sprintf('%s/bin/fslstats %s -R', fsl_path, dwi_file));
if status == 0
    range_vals = str2num(max_val_str);
    if length(range_vals) >= 2
        brightness_threshold = range_vals(2) * 0.15; % 15% of max value
    else
        brightness_threshold = 100; % Default fallback
    end
else
    brightness_threshold = 100; % Default fallback
end

% SUSAN parameters
sigma = 1.0; % Spatial kernel size (mm)

% Run FSL SUSAN
cmd_susan = sprintf('%s/bin/susan %s %.1f %.1f 3 1 1 %s', ...
                   fsl_path, dwi_file, brightness_threshold, sigma, denoised_file);

fprintf('Command: %s\n', cmd_susan);
fprintf('Brightness threshold: %.1f\n', brightness_threshold);

tic;
[status, cmdout] = system(cmd_susan);
elapsed_time = toc;

if status ~= 0
    error('FSL SUSAN denoising failed: %s', cmdout);
end

fprintf('  Processing time: %.1f seconds\n', elapsed_time);

end

function denoised_file = denoise_nlmeans(dwi_file, denoised_file)
% Non-local means denoising (requires additional toolbox or external tool)

fprintf('Using Non-local means denoising...\n');

% Check if we have Image Processing Toolbox for imnlmfilt
if license('test', 'image_toolbox')
    fprintf('Using MATLAB Image Processing Toolbox...\n');
    denoise_matlab_nlmeans(dwi_file, denoised_file);
else
    fprintf('⚠ Image Processing Toolbox not available, falling back to Gaussian smoothing\n');
    denoised_file = denoise_gaussian(dwi_file, denoised_file);
end

end

function denoise_matlab_nlmeans(dwi_file, denoised_file)
% Use MATLAB's built-in non-local means (requires Image Processing Toolbox)

% Load the DWI data using SPM
V = spm_vol(dwi_file);
dwi_data = spm_read_vols(V);

fprintf('Applying non-local means to %d volumes...\n', size(dwi_data, 4));

% Apply non-local means denoising to each volume
denoised_data = zeros(size(dwi_data));

for vol = 1:size(dwi_data, 4)
    if mod(vol, 10) == 0
        fprintf('  Processing volume %d/%d\n', vol, size(dwi_data, 4));
    end
    
    % Get current volume and normalize for better denoising
    current_vol = dwi_data(:,:,:,vol);
    
    % Apply slice-wise denoising (more efficient for 3D volumes)
    denoised_vol = zeros(size(current_vol));
    for slice = 1:size(current_vol, 3)
        current_slice = current_vol(:,:,slice);
        if max(current_slice(:)) > 0 % Only process non-zero slices
            denoised_slice = imnlmfilt(current_slice, 'DegreeOfSmoothing', 0.1);
            denoised_vol(:,:,slice) = denoised_slice;
        end
    end
    
    denoised_data(:,:,:,vol) = denoised_vol;
end

% Save the denoised data
V_out = V(1);
V_out.fname = denoised_file;
V_out.private.dat.fname = denoised_file;

spm_write_vol(V_out, denoised_data(:,:,:,1));

% For 4D data, we need to save differently
if size(denoised_data, 4) > 1
    % Create 4D NIfTI
    V_4d = V;
    for i = 1:length(V_4d)
        V_4d(i).fname = denoised_file;
        V_4d(i).private.dat.fname = denoised_file;
    end
    spm_write_vol(V_4d, denoised_data);
end

end

function denoised_file = denoise_gaussian(dwi_file, denoised_file)
% Simple Gaussian smoothing denoising

fprintf('Using Gaussian smoothing denoising...\n');

% Ensure FSL is available
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSL not found for Gaussian smoothing');
end

% Gaussian kernel FWHM (mm) - conservative to preserve resolution
fwhm = 1.0;

% Run FSL's fslmaths with Gaussian smoothing
cmd_smooth = sprintf('%s/bin/fslmaths %s -s %.1f %s', ...
                    fsl_path, dwi_file, fwhm, denoised_file);

fprintf('Command: %s\n', cmd_smooth);
fprintf('FWHM: %.1f mm\n', fwhm);

tic;
[status, cmdout] = system(cmd_smooth);
elapsed_time = toc;

if status ~= 0
    error('Gaussian smoothing failed: %s', cmdout);
end

fprintf('  Processing time: %.1f seconds\n', elapsed_time);

end