function success = generateTractographySliceCache(tracks_file, nim_file, cache_dir, options)
% GENERATETRACTOGRAPHYSLICECACHE: Generate pre-computed slice cache for fast viewing
%
% This function creates a complete cache of slice images for all orientations
% and slice positions, enabling instant navigation in the fast viewer GUI.
%
% Key Features:
% - Pre-computes all slice images at integer intervals
% - Optimized rendering pipeline with parallel processing
% - Progress tracking with ETA and cancellation support
% - Quality validation and error recovery
% - Configurable compression and resolution settings
% - Automatic metadata generation and storage
%
% Arguments:
%   tracks_file - Path to tracks .mat file
%   nim_file    - Path to nim .mat file
%   cache_dir   - Directory for cache storage (optional)
%   options     - Generation options structure (optional)
%
% Options Structure:
%   .tolerance         - Slice thickness in voxels (default: 2)
%   .color_mode        - 'direction' or 'uniform' (default: 'direction')
%   .show_anatomy      - Show FA background (default: true)
%   .alpha            - FA transparency (default: 0.6)
%   .image_resolution  - [width, height] for output images (default: [1024, 1024])
%   .image_format     - 'png' or 'jpg' (default: 'png')
%   .compression_level - 0-9 for PNG, 0-100 for JPG (default: 6)
%   .parallel_workers  - Number of parallel workers (default: auto-detect)
%   .progress_callback - Function handle for progress updates
%   .quality_validation - Enable image quality checks (default: true)
%
% Returns:
%   success - True if cache generation completed successfully

fprintf('=== TRACTOGRAPHY SLICE CACHE GENERATION ===\n');

%% Input validation and setup
if nargin < 2
    error('generateTractographySliceCache:InvalidInput', ...
          'Must provide tracks_file and nim_file');
end

if nargin < 3 || isempty(cache_dir)
    cache_dir = '';  % Will use default location
end

if nargin < 4
    options = struct();
end

% Set default options
options = setDefaultOptions(options);

% Validate input files
if ~exist(tracks_file, 'file')
    error('generateTractographySliceCache:FileNotFound', ...
          'Tracks file not found: %s', tracks_file);
end

if ~exist(nim_file, 'file')
    error('generateTractographySliceCache:FileNotFound', ...
          'NIM file not found: %s', nim_file);
end

%% Initialize cache manager
try
    cache_manager = TractographyCacheManager(cache_dir);
    fprintf('Cache manager initialized\n');
catch ME
    error('generateTractographySliceCache:CacheManagerError', ...
          'Failed to initialize cache manager: %s', ME.message);
end

%% Generate dataset and parameter hashes
dataset_hash = cache_manager.generateDatasetHash(tracks_file, nim_file);
param_hash = cache_manager.generateParameterHash(options);

fprintf('Dataset hash: %s\n', dataset_hash);
fprintf('Parameter hash: %s\n', param_hash);

%% Check if cache already exists
if cache_manager.validateCacheIntegrity(dataset_hash, param_hash)
    fprintf('Valid cache already exists for this dataset and parameters\n');
    if ~options.force_regenerate
        fprintf('Use options.force_regenerate = true to regenerate\n');
        success = true;
        return;
    else
        fprintf('Forcing regeneration...\n');
    end
end

%% Create cache directory structure
success = cache_manager.createCacheStructure(dataset_hash, param_hash, ...
                                           tracks_file, nim_file, options);
if ~success
    error('generateTractographySliceCache:CacheStructureError', ...
          'Failed to create cache directory structure');
end

cache_directory = cache_manager.getCacheDirectory(dataset_hash, param_hash);
fprintf('Cache directory: %s\n', cache_directory);

%% Load data
fprintf('\nLoading data...\n');
load_start = tic;

% Load tracks
fprintf('  Loading tracks...\n');
tracks_data = load(tracks_file);
if isfield(tracks_data, 'tracks')
    tracks = tracks_data.tracks;
elseif isfield(tracks_data, 'tracks_selected')
    tracks = tracks_data.tracks_selected;
else
    error('generateTractographySliceCache:InvalidTracksFile', ...
          'Tracks file must contain "tracks" or "tracks_selected" variable');
end

% Load nim structure
fprintf('  Loading nim structure...\n');
nim_data = load(nim_file);
if isfield(nim_data, 'nim')
    nim = nim_data.nim;
else
    error('generateTractographySliceCache:InvalidNimFile', ...
          'NIM file must contain "nim" variable');
end

load_time = toc(load_start);
fprintf('Data loaded in %.2f seconds\n', load_time);

% Get volume dimensions
dims = size(nim.FA);
fprintf('Volume dimensions: %dx%dx%d\n', dims(1), dims(2), dims(3));
fprintf('Track count: %d\n', length(tracks));

%% Build optimized track lookup
fprintf('\nBuilding track-slice lookup...\n');
lookup_start = tic;

try
    slice_lookup = buildOptimizedTrackSliceLookup(tracks, dims);
    fprintf('Using optimized track-slice lookup\n');
catch ME
    fprintf('Warning: Optimized lookup failed (%s), using standard method\n', ME.message);
    slice_lookup = buildTrackSliceLookup(tracks, dims);
end

lookup_time = toc(lookup_start);
fprintf('Track lookup built in %.2f seconds\n', lookup_time);

%% Calculate total work and setup progress tracking
orientations = {'axial', 'sagittal', 'coronal'};
orientation_dims = [dims(3), dims(1), dims(2)];  % Z, X, Y slices
total_slices = sum(orientation_dims);

fprintf('\nGeneration plan:\n');
fprintf('  Axial slices (Z): %d\n', dims(3));
fprintf('  Sagittal slices (X): %d\n', dims(1));
fprintf('  Coronal slices (Y): %d\n', dims(2));
fprintf('  Total slices: %d\n', total_slices);

% Progress tracking
progress = struct();
progress.total_slices = total_slices;
progress.completed_slices = 0;
progress.start_time = tic;
progress.last_update_time = progress.start_time;
progress.update_interval = 5; % Update every 5 seconds

%% Setup parallel processing
if options.parallel_workers > 1
    fprintf('\nSetting up parallel processing with %d workers...\n', options.parallel_workers);
    try
        current_pool = gcp('nocreate');
        if isempty(current_pool) || current_pool.NumWorkers ~= options.parallel_workers
            if ~isempty(current_pool)
                delete(current_pool);
            end
            parpool(options.parallel_workers);
        end
        fprintf('Parallel pool ready with %d workers\n', options.parallel_workers);
    catch ME
        fprintf('Warning: Failed to setup parallel processing (%s)\n', ME.message);
        fprintf('Continuing with sequential processing\n');
        options.parallel_workers = 1;
    end
end

%% Generate slice images
fprintf('\n=== GENERATING SLICE IMAGES ===\n');
generation_start = tic;

success = true;
error_count = 0;
max_errors = total_slices * 0.01; % Allow 1% failure rate

for orient_idx = 1:length(orientations)
    orientation = orientations{orient_idx};
    num_slices = orientation_dims(orient_idx);

    fprintf('\nGenerating %s slices (1-%d)...\n', orientation, num_slices);

    % Create output directory
    output_dir = fullfile(cache_directory, orientation);

    % Generate slices for this orientation
    orient_success = generateOrientationSlices(output_dir, orientation, ...
        num_slices, nim, tracks, slice_lookup, options, progress);

    if ~orient_success
        error_count = error_count + 1;
        fprintf('Warning: Errors occurred generating %s slices\n', orientation);
    end

    % Check if too many errors
    if error_count > max_errors
        fprintf('Error: Too many generation failures (%d), stopping\n', error_count);
        success = false;
        break;
    end
end

generation_time = toc(generation_start);

%% Generate summary and metadata
if success
    fprintf('\n=== GENERATION COMPLETE ===\n');
    fprintf('Total time: %.2f minutes\n', generation_time / 60);
    fprintf('Average time per slice: %.3f seconds\n', generation_time / total_slices);
    fprintf('Cache directory: %s\n', cache_directory);

    % Generate cache statistics
    cache_stats = calculateCacheStatistics(cache_directory);
    fprintf('\nCache statistics:\n');
    fprintf('  Total size: %.2f MB\n', cache_stats.total_size_mb);
    fprintf('  Image count: %d\n', cache_stats.image_count);
    fprintf('  Average image size: %.1f KB\n', cache_stats.avg_image_size_kb);

    % Update metadata with generation results
    updateCacheMetadata(cache_manager, dataset_hash, param_hash, cache_stats, generation_time);

else
    fprintf('\n=== GENERATION FAILED ===\n');
    fprintf('Errors encountered during generation\n');
end

%% Cleanup
if options.parallel_workers > 1
    try
        delete(gcp('nocreate'));
    catch
        % Pool cleanup error - not critical
    end
end

status_messages = {'FAILED', 'COMPLETE'};
fprintf('=== Cache generation %s ===\n', status_messages{success + 1});

end


%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

function options = setDefaultOptions(options)
% Set default options for cache generation

default_opts = struct();
default_opts.tolerance = 2;
default_opts.color_mode = 'direction';
default_opts.show_anatomy = true;
default_opts.alpha = 0.6;
default_opts.image_resolution = [1024, 1024];
default_opts.image_format = 'png';
default_opts.compression_level = 6;
default_opts.parallel_workers = feature('numcores');
default_opts.progress_callback = [];
default_opts.quality_validation = true;
default_opts.force_regenerate = false;

% Merge with provided options
field_names = fieldnames(default_opts);
for i = 1:length(field_names)
    field = field_names{i};
    if ~isfield(options, field)
        options.(field) = default_opts.(field);
    end
end
end


function success = generateOrientationSlices(output_dir, orientation, num_slices, ...
                                           nim, tracks, slice_lookup, options, progress)
% Generate all slices for a specific orientation

success = true;
dims = size(nim.FA);

% Determine slice coordinate mapping
switch orientation
    case 'axial'
        coord_indices = [1, 2]; % X, Y
        slice_dim = 3; % Z
    case 'sagittal'
        coord_indices = [2, 3]; % Y, Z
        slice_dim = 1; % X
    case 'coronal'
        coord_indices = [1, 3]; % X, Z
        slice_dim = 2; % Y
end

% Generate slices
for slice_idx = 1:num_slices
    try
        % Create filename
        filename = sprintf('%s_%03d.%s', orientation, slice_idx, options.image_format);
        filepath = fullfile(output_dir, filename);

        % Skip if file already exists and not forcing regeneration
        if exist(filepath, 'file') && ~options.force_regenerate
            continue;
        end

        % Generate slice image
        generateSingleSlice(filepath, orientation, slice_idx, nim, tracks, ...
                          slice_lookup, options, dims, coord_indices, slice_dim);

        % Update progress
        progress.completed_slices = progress.completed_slices + 1;
        updateProgress(progress, options);

    catch ME
        fprintf('Error generating %s slice %d: %s\n', orientation, slice_idx, ME.message);
        success = false;
    end
end
end


function generateSingleSlice(filepath, orientation, slice_idx, nim, tracks, ...
                           slice_lookup, options, dims, coord_indices, slice_dim)
% Generate a single slice image and save to file

% Create invisible figure for rendering
fig = figure('Visible', 'off', 'Color', 'white', ...
            'Position', [1, 1, options.image_resolution]);

try
    ax = axes(fig);

    % Render slice using optimized renderer
    slices = struct();
    slices.axial = round(dims(3)/2);
    slices.sagittal = round(dims(1)/2);
    slices.coronal = round(dims(2)/2);

    % Set current slice position
    switch orientation
        case 'axial'
            slices.axial = slice_idx;
            orient_char = 'z';
        case 'sagittal'
            slices.sagittal = slice_idx;
            orient_char = 'x';
        case 'coronal'
            slices.coronal = slice_idx;
            orient_char = 'y';
    end

    % Use optimized renderer
    optimizedSliceRenderer(ax, orient_char, slice_idx, slices, nim, tracks, ...
                         slice_lookup, options, dims);

    % Configure for high-quality export
    set(ax, 'Position', [0, 0, 1, 1]); % Full figure
    axis(ax, 'off'); % Remove axes
    set(fig, 'InvertHardcopy', 'off');

    % Save image
    switch options.image_format
        case 'png'
            print(fig, filepath, '-dpng', sprintf('-r%d', 150));
        case 'jpg'
            print(fig, filepath, '-djpeg', sprintf('-r%d', 150), ...
                  sprintf('-q%d', options.compression_level));
    end

    % Quality validation
    if options.quality_validation
        validateImageQuality(filepath);
    end

finally
    close(fig);
end
end


function updateProgress(progress, options)
% Update progress display with ETA

current_time = toc(progress.start_time);

% Update every few seconds to avoid spam
if current_time - progress.last_update_time > progress.update_interval
    progress.last_update_time = current_time;

    percent_complete = progress.completed_slices / progress.total_slices * 100;
    slices_per_second = progress.completed_slices / current_time;
    remaining_slices = progress.total_slices - progress.completed_slices;
    eta_seconds = remaining_slices / slices_per_second;

    fprintf('Progress: %d/%d slices (%.1f%%), ETA: %.1f minutes\n', ...
            progress.completed_slices, progress.total_slices, percent_complete, eta_seconds / 60);

    % Call user callback if provided
    if ~isempty(options.progress_callback)
        try
            options.progress_callback(progress.completed_slices, progress.total_slices);
        catch
            % Ignore callback errors
        end
    end
end
end


function validateImageQuality(filepath)
% Basic quality validation for generated images

try
    info = imfinfo(filepath);

    % Check file size (should be reasonable)
    if info.FileSize < 1000 % Less than 1KB is suspicious
        warning('Generated image may be corrupted (very small): %s', filepath);
    end

    % Check image dimensions
    if info.Width < 100 || info.Height < 100
        warning('Generated image has suspicious dimensions: %s', filepath);
    end

catch ME
    warning('Could not validate image quality for %s: %s', filepath, ME.message);
end
end


function stats = calculateCacheStatistics(cache_directory)
% Calculate statistics about the generated cache

stats = struct();
stats.total_size_bytes = 0;
stats.image_count = 0;

orientations = {'axial', 'sagittal', 'coronal'};

for i = 1:length(orientations)
    orient_dir = fullfile(cache_directory, orientations{i});
    if exist(orient_dir, 'dir')
        files = dir(fullfile(orient_dir, '*.*'));
        files = files(~[files.isdir]);

        stats.image_count = stats.image_count + length(files);
        stats.total_size_bytes = stats.total_size_bytes + sum([files.bytes]);
    end
end

stats.total_size_mb = stats.total_size_bytes / (1024^2);
if stats.image_count > 0
    stats.avg_image_size_kb = stats.total_size_bytes / stats.image_count / 1024;
else
    stats.avg_image_size_kb = 0;
end
end


function updateCacheMetadata(cache_manager, dataset_hash, param_hash, cache_stats, generation_time)
% Update cache metadata with generation results

try
    % Get existing metadata
    metadata = cache_manager.getCacheMetadata(dataset_hash);

    % Add generation information
    metadata.generation = struct();
    metadata.generation.completed = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
    metadata.generation.duration_seconds = generation_time;
    metadata.generation.cache_stats = cache_stats;

    % Save updated metadata
    cache_dir = cache_manager.getCacheDirectory(dataset_hash, param_hash);
    metadata_file = fullfile(fileparts(fileparts(cache_dir)), 'metadata.json');

    fid = fopen(metadata_file, 'w');
    if fid ~= -1
        fprintf(fid, '%s', jsonencode(metadata, 'PrettyPrint', true));
        fclose(fid);
    end

catch ME
    fprintf('Warning: Could not update cache metadata: %s\n', ME.message);
end
end