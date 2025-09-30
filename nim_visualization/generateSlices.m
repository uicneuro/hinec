function cache_dir = generateSlices(tracks_file, nim_file, output_dir, options)
% GENERATESLICES: Server-side slice image generation for distributed viewing
%
% This function generates all pre-computed slice images on a high-resource
% server. The output can be transferred to any computer for instant viewing
% with the Python GUI (no MATLAB required on local machine).
%
% WORKFLOW:
%   1. SERVER: Run this function to generate slice images
%   2. TRANSFER: Copy output_dir to local computer
%   3. LOCAL: Run Python viewer to view slices
%
% Usage:
%   cache_dir = generateSlices(tracks_file, nim_file)
%   cache_dir = generateSlices(tracks_file, nim_file, output_dir)
%   cache_dir = generateSlices(tracks_file, nim_file, output_dir, options)
%
% Arguments:
%   tracks_file - Path to tracks .mat file
%   nim_file    - Path to nim .mat file with FA data
%   output_dir  - Output directory for slice images (optional)
%                 Default: './tractography_cache'
%   options     - Rendering options structure (optional)
%
% Options:
%   .tolerance         - Slice thickness in voxels (default: 2)
%   .color_mode        - 'direction' or 'uniform' (default: 'direction')
%   .show_anatomy      - Show FA background (default: true)
%   .alpha             - FA transparency (default: 0.6)
%   .image_resolution  - [width, height] pixels (default: [1024, 1024])
%   .image_format      - 'png' or 'jpg' (default: 'png')
%   .compression_level - PNG: 0-9, JPG: 0-100 (default: 6 or 90)
%   .parallel_workers  - Number of workers (default: auto)
%
% Returns:
%   cache_dir - Full path to generated cache directory
%
% Output Structure:
%   output_dir/
%   ├── datasets/
%   │   └── {dataset_hash}/
%   │       ├── metadata.json
%   │       └── parameters/
%   │           └── {param_hash}/
%   │               ├── config.json
%   │               ├── axial/*.png (or .jpg)
%   │               ├── sagittal/*.png
%   │               └── coronal/*.png
%
% Transfer Instructions:
%   1. Zip the output directory:
%      tar -czf slices.tar.gz output_dir/
%   2. Transfer to local computer (scp, rsync, USB, etc.)
%   3. Extract on local computer
%   4. Run Python viewer:
%      python FastTractographyViewer.py output_dir/
%
% Example:
%   % On high-resource server:
%   cache = generateSlices('tracks.mat', 'nim.mat', '/export/slices');
%
%   % Transfer /export/slices to local computer
%
%   % On local computer (no MATLAB needed):
%   % python FastTractographyViewer.py /path/to/slices

fprintf('╔═══════════════════════════════════════════════════════════╗\n');
fprintf('║   TRACTOGRAPHY SLICE GENERATOR (Server-Side Processing)  ║\n');
fprintf('╚═══════════════════════════════════════════════════════════╝\n\n');

%% Input validation
if nargin < 2
    error('generateSlices:InvalidInput', ...
          'Must provide tracks_file and nim_file');
end

if nargin < 3 || isempty(output_dir)
    output_dir = './tractography_cache';
    fprintf('Using default output directory: %s\n', output_dir);
end

if nargin < 4
    options = struct();
end

% Validate input files
fprintf('Validating input files...\n');
if ~exist(tracks_file, 'file')
    error('generateSlices:FileNotFound', ...
          'Tracks file not found: %s', tracks_file);
end

if ~exist(nim_file, 'file')
    error('generateSlices:FileNotFound', ...
          'NIM file not found: %s', nim_file);
end

fprintf('  ✓ Tracks: %s\n', tracks_file);
fprintf('  ✓ NIM:    %s\n', nim_file);

%% Set default options
options = setDefaultGenerationOptions(options);

fprintf('\nGeneration Settings:\n');
fprintf('  Image format:    %s\n', upper(options.image_format));
fprintf('  Resolution:      %dx%d pixels\n', options.image_resolution);
fprintf('  Slice thickness: %.1f voxels\n', options.tolerance);
fprintf('  Color mode:      %s\n', options.color_mode);
fprintf('  Parallel workers: %d\n', options.parallel_workers);
fprintf('  Show anatomy:    %s\n', mat2str(options.show_anatomy));

%% Create output directory
fprintf('\nPreparing output directory...\n');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('  ✓ Created: %s\n', output_dir);
else
    fprintf('  ✓ Using existing: %s\n', output_dir);
end

%% Generate slice cache using optimized pipeline
fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              GENERATING SLICE IMAGES                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

try
    success = generateTractographySliceCache(tracks_file, nim_file, output_dir, options);

    if ~success
        error('generateSlices:GenerationFailed', ...
              'Slice generation failed');
    end

catch ME
    fprintf('\n✗ ERROR: %s\n', ME.message);
    fprintf('  Location: %s\n', ME.stack(1).name);
    fprintf('  Line: %d\n', ME.stack(1).line);
    rethrow(ME);
end

%% Get cache directory structure
cache_manager = TractographyCacheManager(output_dir);
dataset_hash = cache_manager.generateDatasetHash(tracks_file, nim_file);
param_hash = cache_manager.generateParameterHash(options);

cache_dir = fullfile(output_dir, 'datasets', dataset_hash, 'parameters', param_hash);

%% Validate output
fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              VALIDATING OUTPUT                             ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

validation = validateGeneratedCache(cache_dir, options);

if validation.valid
    fprintf('✓ Validation PASSED\n');
    fprintf('  Total images: %d\n', validation.total_images);
    fprintf('  Axial:        %d slices\n', validation.counts.axial);
    fprintf('  Sagittal:     %d slices\n', validation.counts.sagittal);
    fprintf('  Coronal:      %d slices\n', validation.counts.coronal);
    fprintf('  Total size:   %.2f MB\n', validation.total_size_mb);
else
    warning('generateSlices:ValidationFailed', ...
            'Output validation failed: %s', validation.error);
end

%% Display transfer instructions
fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║              TRANSFER INSTRUCTIONS                         ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

fprintf('Cache generated successfully!\n\n');
fprintf('OUTPUT LOCATION:\n');
fprintf('  %s\n\n', output_dir);

fprintf('TRANSFER TO LOCAL COMPUTER:\n');
fprintf('  Method 1 - Compress and transfer:\n');
fprintf('    tar -czf slices.tar.gz "%s"\n', output_dir);
fprintf('    scp slices.tar.gz user@local:/destination/\n\n');

fprintf('  Method 2 - Direct sync:\n');
fprintf('    rsync -avz "%s" user@local:/destination/\n\n', output_dir);

fprintf('  Method 3 - USB/Network drive:\n');
fprintf('    Copy entire directory to portable storage\n\n');

fprintf('VIEW ON LOCAL COMPUTER (no MATLAB required):\n');
fprintf('  python FastTractographyViewer.py "%s"\n\n', output_dir);

fprintf('REQUIREMENTS ON LOCAL COMPUTER:\n');
fprintf('  • Python 3.7+\n');
fprintf('  • pip install pillow numpy\n');
fprintf('  • FastTractographyViewer.py script\n\n');

fprintf('═══════════════════════════════════════════════════════════\n');
fprintf('Generation complete! Ready for transfer.\n');
fprintf('═══════════════════════════════════════════════════════════\n\n');

end


function options = setDefaultGenerationOptions(options)
% Set default options for slice generation

default_opts = struct();
default_opts.tolerance = 2;
default_opts.color_mode = 'direction';
default_opts.show_anatomy = true;
default_opts.alpha = 0.6;
default_opts.image_resolution = [1024, 1024];
default_opts.image_format = 'png';
default_opts.compression_level = [];  % Auto-set based on format
default_opts.parallel_workers = [];   % Auto-detect

% Merge with provided options
field_names = fieldnames(default_opts);
for i = 1:length(field_names)
    field = field_names{i};
    if ~isfield(options, field) || isempty(options.(field))
        options.(field) = default_opts.(field);
    end
end

% Set compression level based on format
if isempty(options.compression_level)
    if strcmp(options.image_format, 'png')
        options.compression_level = 6;  % PNG compression 0-9
    else
        options.compression_level = 90; % JPG quality 0-100
    end
end

% Auto-detect parallel workers
if isempty(options.parallel_workers)
    try
        pool = gcp('nocreate');
        if isempty(pool)
            options.parallel_workers = min(4, feature('numcores'));
        else
            options.parallel_workers = pool.NumWorkers;
        end
    catch
        options.parallel_workers = 1;
    end
end

end


function validation = validateGeneratedCache(cache_dir, options)
% Validate generated cache structure and contents

validation = struct();
validation.valid = false;
validation.total_images = 0;
validation.counts = struct('axial', 0, 'sagittal', 0, 'coronal', 0);
validation.total_size_mb = 0;
validation.error = '';

try
    % Check if cache directory exists
    if ~exist(cache_dir, 'dir')
        validation.error = 'Cache directory not found';
        return;
    end

    % Check for required subdirectories
    orientations = {'axial', 'sagittal', 'coronal'};
    for i = 1:length(orientations)
        orient = orientations{i};
        orient_dir = fullfile(cache_dir, orient);

        if ~exist(orient_dir, 'dir')
            validation.error = sprintf('Missing %s directory', orient);
            return;
        end

        % Count images in this orientation
        image_ext = ['*.' options.image_format];
        image_files = dir(fullfile(orient_dir, image_ext));
        count = length(image_files);

        validation.counts.(orient) = count;
        validation.total_images = validation.total_images + count;

        % Sum file sizes
        for j = 1:length(image_files)
            validation.total_size_mb = validation.total_size_mb + image_files(j).bytes / (1024*1024);
        end

        if count == 0
            validation.error = sprintf('No images found in %s directory', orient);
            return;
        end
    end

    % Check for metadata files
    if ~exist(fullfile(cache_dir, 'config.json'), 'file')
        validation.error = 'Missing config.json';
        return;
    end

    % All checks passed
    validation.valid = true;

catch ME
    validation.error = ME.message;
end

end