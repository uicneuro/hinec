function visualizeTractography(tracks_file, nim_file, varargin)
% VISUALIZETRACTOGRAPHY: Complete tractography visualization with multiple modes and export
%
% This single function replaces ALL previous visualization files:
% - visualizeTractography.m (whole brain)
% - visualizeTractographyRegion.m (single region)
% - visualizeTractographyAllRegions.m (grid layout)
% - nim_plot_tractography_region.m (core plotting)
%
% Usage:
%   visualizeTractography(tracks_file, nim_file)                           % Whole brain
%   visualizeTractography(tracks_file, nim_file, 'region', 5)              % Single region
%   visualizeTractography(tracks_file, nim_file, 'mode', 'grid')           % All regions grid
%   visualizeTractography(tracks_file, nim_file, 'mode', 'sequential')     % Sequential regions
%   visualizeTractography(..., 'export_dir', 'figures/')                  % With image export
%
% Arguments:
%   tracks_file - Path to tracks .mat file (REQUIRED)
%   nim_file - Path to nim .mat file (REQUIRED)
%
% Parameters:
%   'mode' - 'whole', 'region', 'grid', 'sequential' (default: 'whole')
%   'region' - Region ID for single region mode (alternative to mode='region')
%   'export_dir' - Directory for image export (auto-creates subdirectories)
%   'export_format' - 'png', 'pdf', 'eps', 'fig' (default: 'png')
%   'export_dpi' - Export resolution (default: 300)
%
% Track Filtering:
%   'filter_mode' - 'pass_through', 'start_in', 'end_in', 'connect_to' (default: 'pass_through')
%   'min_overlap' - Minimum region overlap ratio 0-1 (default: 0.05)
%   'min_region_points' - Minimum in-region track points for pass-through (default: 3)
%   'max_tracks' - Max tracks to display (default: 2000 whole, inf others)
%
% Visualization:
%   'color_mode' - 'direction', 'fa', 'uniform', 'region' (default: 'direction')
%   'show_anatomy' - Show FA background (default: true)
%   'show_region' - Show region overlay (default: true)
%   'line_width' - Track line width (default: 1.2)
%   'region_alpha' - Region transparency (default: 0.3)
%
% Grid Layout:
%   'grid_cols' - Grid columns (default: auto)
%   'subplot_size' - [width, height] (default: [400, 350])
%   'show_all_tracks' - Force show ALL tracks in grid (default: true)
%
% Examples:
%   % Basic whole brain
%   visualizeTractography('tracks.mat', 'data.mat')
%
%   % Single region with export
%   visualizeTractography('tracks.mat', 'data.mat', 'region', 5, 'export_dir', 'figs/')
%
%   % All regions grid
%   visualizeTractography('tracks.mat', 'data.mat', 'mode', 'grid', 'filter_mode', 'start_in')
%
%   % Sequential with export
%   visualizeTractography('tracks.mat', 'data.mat', 'mode', 'sequential', 'export_dir', 'output/')

%% INPUT VALIDATION AND PARSING
if nargin < 2
    error('visualizeTractography requires 2 arguments: tracks_file, nim_file');
end

if ~ischar(tracks_file) && ~isstring(tracks_file)
    error('tracks_file must be a string path to .mat file');
end

if ~ischar(nim_file) && ~isstring(nim_file)
    error('nim_file must be a string path to .mat file');
end

% Parse input parameters with comprehensive options
p = inputParser;
addRequired(p, 'tracks_file');
addRequired(p, 'nim_file');

% Mode selection - support both 'mode' and direct 'region' for compatibility
addParameter(p, 'mode', 'whole', @(x) ismember(x, {'whole', 'region', 'grid', 'sequential'}));
addParameter(p, 'region', [], @(x) isempty(x) || (isnumeric(x) && x > 0));

% Export parameters
addParameter(p, 'export_dir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'export_format', 'png', @(x) ismember(x, {'png', 'pdf', 'eps', 'fig'}));
addParameter(p, 'export_dpi', 300, @(x) isnumeric(x) && x > 0);

% Track filtering parameters
addParameter(p, 'filter_mode', 'pass_through', @(x) ismember(x, {'pass_through', 'start_in', 'end_in', 'connect_to'}));
addParameter(p, 'min_overlap', 0.05, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'min_region_points', 3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'max_tracks', [], @(x) isempty(x) || (isnumeric(x) && x > 0));

% Visualization parameters
addParameter(p, 'color_mode', 'direction', @(x) ismember(x, {'direction', 'fa', 'uniform', 'region'}));
addParameter(p, 'show_anatomy', true, @islogical);
addParameter(p, 'show_region', true, @islogical);
addParameter(p, 'line_width', 1.2, @(x) isnumeric(x) && x > 0);
addParameter(p, 'region_alpha', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);

% Grid layout parameters
addParameter(p, 'grid_cols', 0, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'subplot_size', [400, 350], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'show_all_tracks', true, @islogical);

% Legacy compatibility
addParameter(p, 'show_region_overlays', [], @(x) isempty(x) || islogical(x));
addParameter(p, 'show_background', [], @(x) isempty(x) || islogical(x));

% Export behavior
addParameter(p, 'silent_export', true, @islogical); % Don't display figures when exporting

parse(p, tracks_file, nim_file, varargin{:});
options = p.Results;

% Handle legacy parameter names
if ~isempty(options.show_region_overlays)
    options.show_region = options.show_region_overlays;
end
if ~isempty(options.show_background)
    options.show_anatomy = options.show_background;
end

% Auto-detect mode from 'region' parameter
if ~isempty(options.region)
    options.mode = 'region';
    options.region_id = options.region;
else
    options.region_id = [];
end

% Auto-enable silent export when export directory is specified
if ~isempty(options.export_dir) && options.silent_export
    options.display_figures = false;
else
    options.display_figures = true;
end

% Validate mode-specific requirements
if strcmp(options.mode, 'region') && isempty(options.region_id)
    error('region_id is required when mode is "region". Use: visualizeTractography(..., "region", 5)');
end

% Normalize thresholds
options.min_region_points = max(0, round(options.min_region_points));
options.min_overlap = max(0, min(1, options.min_overlap));

% Set default max_tracks based on mode
if isempty(options.max_tracks)
    switch options.mode
        case 'whole'
            options.max_tracks = 2000;  % Performance limit
        case {'region', 'grid', 'sequential'}
            options.max_tracks = inf;   % Show all tracks
    end
end

%% ADD PATHS AND INITIALIZE
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');
addpath('nim_parcellation');

fprintf('=== UNIFIED TRACTOGRAPHY VISUALIZATION ===\n');
fprintf('Mode: %s\n', options.mode);
if strcmp(options.mode, 'region')
    fprintf('Region ID: %d\n', options.region_id);
end
if ~isempty(options.export_dir)
    fprintf('Export mode: Silent export enabled (no figure display or pauses)\n');
    fprintf('Export directory: %s\n', options.export_dir);
end

%% LOAD AND VALIDATE DATA
[tracks, nim] = load_and_validate_data(tracks_file, nim_file, options);

%% EXECUTE VISUALIZATION MODE
fig_handle = [];
export_filename = '';

switch options.mode
    case 'whole'
        [fig_handle, export_filename] = visualize_whole_brain(tracks, nim, options);

    case 'region'
        [fig_handle, export_filename] = visualize_single_region(tracks, nim, options);

    case 'grid'
        [fig_handle, export_filename] = visualize_all_regions_grid(tracks, nim, options);

    case 'sequential'
        [fig_handle, export_filename] = visualize_all_regions_sequential(tracks, nim, options);
end

%% EXPORT FIGURE IF REQUESTED
if ~isempty(options.export_dir) && ~isempty(fig_handle)
    export_figure(fig_handle, export_filename, options);
end

fprintf('Visualization complete!\n');
fprintf('==========================================\n');
end


%% ============================================================================
%% CORE DATA LOADING AND VALIDATION
%% ============================================================================

function [tracks, nim] = load_and_validate_data(tracks_file, nim_file, options)
% Comprehensive data loading with validation

fprintf('Loading and validating data...\n');

%% Load tracks data
fprintf('  Loading tracks from %s...\n', tracks_file);
if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s\nPlease run tractography first.', tracks_file);
end

try
    track_data = load(tracks_file);
    if ~isfield(track_data, 'tracks')
        error('Invalid tracks file: missing "tracks" field');
    end
    tracks = track_data.tracks;

    if isempty(tracks)
        error('No tracks found in tracks file. Please run tractography first.');
    end

    fprintf('    ✓ Loaded %d tracks\n', length(tracks));
catch ME
    error('Failed to load tracks file: %s', ME.message);
end

%% Load nim data
fprintf('  Loading nim data from %s...\n', nim_file);
if ~exist(nim_file, 'file')
    error('Nim file not found: %s\nPlease run main() first.', nim_file);
end

try
    nim_data = load(nim_file);
    if ~isfield(nim_data, 'nim')
        error('Invalid nim file: missing "nim" field');
    end
    nim = nim_data.nim;

    % Validate required fields
    if ~isfield(nim, 'FA')
        error('nim structure missing required field: FA');
    end

    fprintf('    ✓ Loaded nim data with FA dimensions: %s\n', mat2str(size(nim.FA)));
catch ME
    error('Failed to load nim file: %s', ME.message);
end

%% Validate parcellation for region modes
if ismember(options.mode, {'region', 'grid', 'sequential'})
    if ~isfield(nim, 'parcellation_mask')
        error('nim structure must contain parcellation_mask for region-based visualization. Please run parcellation first.');
    end

    unique_regions = unique(nim.parcellation_mask(:));
    unique_regions = unique_regions(unique_regions > 0);

    if isempty(unique_regions)
        error('No regions found in parcellation mask');
    end

    fprintf('    ✓ Found %d regions in parcellation\n', length(unique_regions));

    % Validate specific region if specified
    if strcmp(options.mode, 'region')
        max_region = max(unique_regions);
        if options.region_id > max_region
            error('Region %d does not exist. Available regions: 1-%d', options.region_id, max_region);
        end

        if sum(nim.parcellation_mask(:) == options.region_id) == 0
            error('Region %d contains no voxels in the parcellation mask', options.region_id);
        end
    end
end

fprintf('  Data validation complete ✓\n');
end


%% ============================================================================
%% WHOLE BRAIN VISUALIZATION (replaces original visualizeTractography.m)
%% ============================================================================

function [fig_handle, export_filename] = visualize_whole_brain(tracks, nim, options)
% Complete whole brain visualization with comprehensive layout

fprintf('Creating whole brain visualization...\n');

% Create figure with comprehensive layout
if options.display_figures
    fig_handle = figure('Name', 'Tractography Visualization', 'Position', [100, 100, 1400, 900], 'Color', 'w');
else
    fig_handle = figure('Name', 'Tractography Visualization', 'Position', [100, 100, 1400, 900], 'Color', 'w', 'Visible', 'off');
end

% Main 3D tracks view (takes up most space)
subplot(2,3,[1,2,4,5]);
hold on;

% Show anatomical background if requested
if options.show_anatomy
    plot_anatomical_background(nim, 0.15);
end

% Plot tracks with direction-based coloring
max_tracks_display = min(options.max_tracks, length(tracks));
if max_tracks_display < length(tracks)
    fprintf('  Limiting display to %d tracks for performance\n', max_tracks_display);
    track_indices = round(linspace(1, length(tracks), max_tracks_display));
    display_tracks = tracks(track_indices);
else
    display_tracks = tracks;
end

plot_tracks_with_color(display_tracks, nim, options);

% Set up main plot
dims = size(nim.FA);
xlim([1 dims(1)]); ylim([1 dims(2)]); zlim([1 dims(3)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Tractography Results (%d tracks, showing %d)', length(tracks), length(display_tracks)));
grid on; axis equal; view(45, 30);

% Add color legend
add_color_legend(options.color_mode);

% Track length histogram
subplot(2,3,3);
plot_track_length_histogram(tracks);

% Seed distribution
subplot(2,3,6);
plot_seed_distribution(tracks);

% Print summary statistics
print_track_summary(tracks, options);

hold off;

export_filename = 'tractography_whole';
fprintf('  ✓ Whole brain visualization complete\n');
end


%% ============================================================================
%% SINGLE REGION VISUALIZATION (replaces visualizeTractographyRegion.m)
%% ============================================================================

function [fig_handle, export_filename] = visualize_single_region(tracks, nim, options)
% Complete single region visualization with detailed analysis

region_id = options.region_id;
fprintf('=== Visualizing Region %d Tractography ===\n', region_id);

%% Display region information
display_region_info(nim, region_id);

%% Filter tracks for the region
fprintf('Filtering tracks for region %d (%s mode)...\n', region_id, options.filter_mode);
filtered_tracks = filter_tracks_by_region(tracks, nim.parcellation_mask, region_id, options);

if isempty(filtered_tracks)
    warning('No tracks found for region %d with current filter settings', region_id);
    if strcmp(options.filter_mode, 'pass_through')
        fprintf('  Hint: lower "min_overlap" (currently %.2f) or "min_region_points" (currently %d) to capture thin tracts.\n', ...
                options.min_overlap, options.min_region_points);
    end
    fig_handle = [];
    export_filename = '';
    return;
end

fprintf('Found %d tracks related to region %d\n', length(filtered_tracks), region_id);

%% Limit tracks if specified
if ~isinf(options.max_tracks) && length(filtered_tracks) > options.max_tracks
    fprintf('Limiting display to %d tracks for performance\n', options.max_tracks);
    track_indices = round(linspace(1, length(filtered_tracks), options.max_tracks));
    filtered_tracks = filtered_tracks(track_indices);
end

%% Get region name and create figure
region_name = get_region_name(nim, region_id);
figure_title = sprintf('Region %d Tractography', region_id);
if ~isempty(region_name) && ~strcmp(region_name, sprintf('Region %d', region_id))
    figure_title = sprintf('Region %d: %s', region_id, region_name);
end

if options.display_figures
    fig_handle = figure('Name', figure_title, 'Color', 'w', 'Position', [100, 100, 1200, 800]);
else
    fig_handle = figure('Name', figure_title, 'Color', 'w', 'Position', [100, 100, 1200, 800], 'Visible', 'off');
end
hold on;

%% Show region overlay if requested
if options.show_region
    plot_region_overlay(nim.parcellation_mask, region_id, options.region_alpha);
end

%% Show anatomical background
if options.show_anatomy
    plot_anatomical_background(nim, 0.1);
end

%% Plot filtered tracks
plot_tracks_with_color(filtered_tracks, nim, options);

%% Set up axes and formatting
dims = size(nim.FA);
xlim([1 dims(1)]); ylim([1 dims(2)]); zlim([1 dims(3)]);
xlabel('X (voxels)'); ylabel('Y (voxels)'); zlabel('Z (voxels)');
title(figure_title, 'FontSize', 14, 'FontWeight', 'bold');
axis equal; grid on; view(45, 30);
camlight; lighting gouraud;

%% Add legend and statistics
add_color_legend(options.color_mode);
display_track_statistics(filtered_tracks, region_id, region_name, options);

hold off;

export_filename = sprintf('tractography_region-%02d', region_id);
fprintf('Region %d visualization complete!\n', region_id);
end


%% ============================================================================
%% GRID LAYOUT VISUALIZATION (replaces visualizeTractographyAllRegions.m)
%% ============================================================================

function [fig_handle, export_filename] = visualize_all_regions_grid(tracks, nim, options)
% Complete grid layout for all regions with comprehensive display

fprintf('=== Visualizing ALL Regions Tractography in Grid Layout ===\n');

%% Get all available regions
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0);

if isempty(unique_regions)
    error('No regions found in parcellation mask');
end

fprintf('Found %d regions to visualize: %s\n', length(unique_regions), mat2str(unique_regions'));

%% Calculate optimal grid layout
[grid_rows, grid_cols] = calculate_optimal_grid(length(unique_regions), options.grid_cols);
fprintf('Using %dx%d grid layout for %d regions\n', grid_rows, grid_cols, length(unique_regions));

%% Create figure
figure_width = grid_cols * options.subplot_size(1) + 100;
figure_height = grid_rows * options.subplot_size(2) + 150;
figure_title = sprintf('All Regions Tractography Grid (%d regions)', length(unique_regions));
if options.display_figures
    fig_handle = figure('Name', figure_title, 'Color', 'w', 'Position', [50, 50, figure_width, figure_height]);
else
    fig_handle = figure('Name', figure_title, 'Color', 'w', 'Position', [50, 50, figure_width, figure_height], 'Visible', 'off');
end

%% Process each region in grid layout
fprintf('Processing regions and creating grid subplots...\n');
region_info = struct('id', {}, 'name', {}, 'track_count', {}, 'subplot_idx', {});
total_tracks_displayed = 0;

for r = 1:length(unique_regions)
    region_id = unique_regions(r);

    % Create subplot for this region
    subplot(grid_rows, grid_cols, r);
    hold on;

    % Get region name
    region_name = get_region_name(nim, region_id);

    % Filter tracks for this region - show ALL tracks if requested
    region_options = options;
    if options.show_all_tracks
        region_options.max_tracks = inf;
    end

    filtered_tracks = filter_tracks_by_region(tracks, nim.parcellation_mask, region_id, region_options);

    if ~isempty(filtered_tracks)
        % Show anatomical background if requested
        if options.show_anatomy
            plot_anatomical_background_subplot(nim);
        end

        % Show region overlay if requested
        if options.show_region
            plot_region_overlay_subplot(nim.parcellation_mask, region_id, options.region_alpha);
        end

        % Plot ALL tracks for this region
        plot_all_region_tracks(filtered_tracks, nim, options);

        total_tracks_displayed = total_tracks_displayed + length(filtered_tracks);
        if options.show_all_tracks
            fprintf('Region %d (%s): %d tracks (ALL displayed)\n', region_id, region_name, length(filtered_tracks));
        else
            fprintf('Region %d (%s): %d tracks\n', region_id, region_name, length(filtered_tracks));
        end
    else
        fprintf('Region %d (%s): No tracks found\n', region_id, region_name);
    end

    % Store region info
    region_info(end+1).id = region_id;
    region_info(end).name = region_name;
    region_info(end).track_count = length(filtered_tracks);
    region_info(end).subplot_idx = r;

    % Set up individual subplot
    setup_region_subplot(nim, region_id, region_name, length(filtered_tracks));
    hold off;
end

%% Add overall figure title and statistics
title_text = sprintf('%s - Total: %d tracks from %d regions', figure_title, total_tracks_displayed, length(unique_regions));
if options.show_all_tracks
    title_text = [title_text ' (ALL TRACKS SHOWN)'];
end
sgtitle(title_text, 'FontSize', 16, 'FontWeight', 'bold');

% Display comprehensive statistics
display_grid_statistics(region_info, options, total_tracks_displayed);

export_filename = 'tractography_grid';
fprintf('Grid visualization complete - %d total tracks displayed!\n', total_tracks_displayed);
end


%% ============================================================================
%% SEQUENTIAL VISUALIZATION (replaces visualizeAllRegionsSequentially.m)
%% ============================================================================

function [fig_handle, export_filename] = visualize_all_regions_sequential(tracks, nim, options)
% Sequential visualization of all regions with individual detailed views

fprintf('=== Visualizing ALL Regions Sequentially ===\n');

%% Get all available regions
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0);

if isempty(unique_regions)
    error('No regions found in parcellation mask');
end

fprintf('Will visualize %d regions sequentially...\n', length(unique_regions));

%% Process each region sequentially
for r = 1:length(unique_regions)
    region_id = unique_regions(r);

    fprintf('\n=== Region %d/%d ===\n', r, length(unique_regions));

    % Create individual region visualization
    temp_options = options;
    temp_options.mode = 'region';
    temp_options.region_id = region_id;

    [region_fig, ~] = visualize_single_region(tracks, nim, temp_options);

    % Export individual region if export directory specified
    if ~isempty(options.export_dir) && ~isempty(region_fig)
        individual_export_filename = sprintf('sequential_region-%02d', region_id);
        temp_export_options = options;
        temp_export_options.mode = 'sequential';
        export_figure(region_fig, individual_export_filename, temp_export_options);
    end

    % Pause for user interaction (only if not exporting)
    if r < length(unique_regions) && isempty(options.export_dir)
        fprintf('Press any key to continue to next region (or Ctrl+C to stop)...\n');
        pause;
    elseif r < length(unique_regions) && ~isempty(options.export_dir)
        fprintf('Processing region %d/%d (export mode - no pause)\n', r, length(unique_regions));
    end
end

% Return the last figure
fig_handle = gcf;
export_filename = 'tractography_sequential_complete';

fprintf('Sequential visualization complete for all %d regions!\n', length(unique_regions));
end


%% ============================================================================
%% CORE PLOTTING FUNCTIONS
%% ============================================================================

function plot_tracks_with_color(tracks, nim, options)
% Plot tracks with specified coloring mode

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 4  % Require at least 4 points to eliminate gray straight lines
        continue;
    end

    % Get track color based on mode
    track_color = get_track_color(track, nim, options.color_mode);

    % Skip tracks with invalid color (empty means skip)
    if isempty(track_color)
        continue;
    end

    % Plot track
    plot3(track(:,1), track(:,2), track(:,3), ...
          'Color', track_color, 'LineWidth', options.line_width);
end
end


function color = get_track_color(track, nim, color_mode)
% Get color for a track based on specified mode

switch color_mode
    case 'direction'
        % Direction-based coloring (RGB from normalized direction vector)
        directions = diff(track);
        if ~isempty(directions)
            avg_direction = mean(directions, 1);
            avg_direction_norm = norm(avg_direction);
            if avg_direction_norm > 0
                color = abs(avg_direction) / avg_direction_norm;
            else
                color = [];  % Return empty to skip tracks with no meaningful direction
            end
        else
            color = [];  % Return empty to skip tracks with no direction data
        end

    case 'fa'
        % FA-based coloring
        fa_values = sample_fa_along_track(track, nim.FA);
        avg_fa = mean(fa_values);
        cmap = hot(256);
        fa_idx = round(avg_fa * 255) + 1;
        fa_idx = max(1, min(256, fa_idx));
        color = cmap(fa_idx, :);

    case 'uniform'
        % Single color for all tracks
        color = [0.2, 0.6, 0.8];

    case 'region'
        % Region-based coloring
        if isfield(nim, 'parcellation_mask')
            track_labels = get_track_parcellation_labels(track, nim.parcellation_mask);
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                mode_label = mode(valid_labels);
                rng(mode_label); % Consistent color for same region
                color = rand(1, 3) * 0.8 + 0.2;
            else
                color = [];  % Return empty to skip tracks with no meaningful direction
            end
        else
            color = [];  % Return empty to skip tracks with no direction data
        end

    otherwise
        color = [0.2, 0.6, 0.8];
end
end


function fa_values = sample_fa_along_track(track, fa_volume)
% Sample FA values along a track
fa_values = zeros(size(track, 1), 1);

for i = 1:size(track, 1)
    pos = track(i, :);
    if all(pos >= 1) && all(pos <= size(fa_volume))
        try
            fa_values(i) = interp3(fa_volume, pos(2), pos(1), pos(3), 'linear', 0);
        catch
            % Fallback to nearest neighbor
            pos_round = round(pos);
            if all(pos_round >= 1) && all(pos_round <= size(fa_volume))
                fa_values(i) = fa_volume(pos_round(1), pos_round(2), pos_round(3));
            end
        end
    end
end
end


function plot_anatomical_background(nim, alpha_value)
% Plot FA slices as anatomical background

if nargin < 2
    alpha_value = 0.15;
end

dims = size(nim.FA);
slice_step = round(dims(3) / 5);
slices = slice_step:slice_step:dims(3)-slice_step;

for s = slices
    fa_slice = nim.FA(:, :, s);
    [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
    Z = ones(size(X)) * s;
    surf(X, Y, Z, fa_slice, 'EdgeColor', 'none', 'FaceAlpha', alpha_value);
end
colormap(gray);
end


function plot_anatomical_background_subplot(nim)
% Plot FA background for subplot (lighter)
dims = size(nim.FA);
slice_step = max(1, round(dims(3) / 6));
slices = slice_step:slice_step:dims(3)-slice_step;

for s = slices
    fa_slice = nim.FA(:, :, s);
    [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
    Z = ones(size(X)) * s;
    surf(X, Y, Z, fa_slice, 'EdgeColor', 'none', 'FaceAlpha', 0.05);
end
colormap(gray);
end


function plot_region_overlay(parcellation_mask, region_id, alpha_value)
% Plot the parcellation region as a 3D overlay

region_mask = parcellation_mask == region_id;

if sum(region_mask(:)) > 100 % Only if region has sufficient voxels
    try
        % Smooth the mask slightly for better visualization
        region_mask_smooth = smooth3(double(region_mask), 'box', 3);

        % Create isosurface
        [faces, vertices] = isosurface(region_mask_smooth, 0.5);

        if ~isempty(faces)
            % COORDINATE FIX: Reorder vertices to match track coordinate system
            vertices_fixed = vertices;
            vertices_fixed(:,1) = vertices(:,2); % X = old Y
            vertices_fixed(:,2) = vertices(:,1); % Y = old X
            % Z stays the same

            % Plot region surface
            patch('Vertices', vertices_fixed, 'Faces', faces, ...
                  'FaceColor', 'red', 'EdgeColor', 'none', ...
                  'FaceAlpha', alpha_value, 'SpecularStrength', 0.1);
        end
    catch
        % Fallback to scatter plot of region voxels
        [i, j, k] = ind2sub(size(region_mask), find(region_mask));
        if ~isempty(i)
            scatter3(j, i, k, 10, 'red', 'filled', 'MarkerFaceAlpha', alpha_value);
        end
    end
end
end


function plot_region_overlay_subplot(parcellation_mask, region_id, alpha_value)
% Plot region overlay for subplot
region_mask = parcellation_mask == region_id;

if sum(region_mask(:)) > 50
    try
        region_mask_smooth = smooth3(double(region_mask), 'box', 3);
        [faces, vertices] = isosurface(region_mask_smooth, 0.5);

        if ~isempty(faces)
            vertices_fixed = vertices;
            vertices_fixed(:,1) = vertices(:,2);
            vertices_fixed(:,2) = vertices(:,1);

            patch('Vertices', vertices_fixed, 'Faces', faces, ...
                  'FaceColor', [1.0, 0.2, 0.2], 'EdgeColor', 'none', ...
                  'FaceAlpha', alpha_value, 'SpecularStrength', 0.1);
        end
    catch
        [i, j, k] = ind2sub(size(region_mask), find(region_mask));
        if ~isempty(i)
            scatter3(j, i, k, 8, 'red', 'filled', 'MarkerFaceAlpha', alpha_value);
        end
    end
end
end


function plot_all_region_tracks(tracks, ~, options)
% Plot ALL tracks for a region with direction-based coloring
fprintf('  Plotting %d tracks (ALL)...\n', length(tracks));

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) >= 4  % Require at least 4 points for meaningful tracks
        % Direction-based color
        directions = diff(track);
        if ~isempty(directions)
            avg_direction = mean(directions, 1);
            avg_direction_norm = norm(avg_direction);
            if avg_direction_norm > 0
                color = abs(avg_direction) / avg_direction_norm;
            else
                color = [];  % Return empty to skip tracks with no meaningful direction
            end
        else
            color = [];  % Return empty to skip tracks with no direction data
        end

        % Skip tracks with invalid color (empty means skip)
        if ~isempty(color)
            plot3(track(:,1), track(:,2), track(:,3), ...
                  'Color', color, 'LineWidth', options.line_width);
        end
    end
end
end


%% ============================================================================
%% TRACK FILTERING FUNCTIONS
%% ============================================================================

function filtered_tracks = filter_tracks_by_region(tracks, parcellation_mask, region_id, options)
% Filter tracks based on their relationship to the specified region

filtered_tracks = cell(length(tracks), 1);
track_count = 0;

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 10  % Require at least 10 points to reduce spurious tracks
        continue;
    end

    % Get parcellation labels along the track
    track_labels = get_track_parcellation_labels(track, parcellation_mask);
    if isempty(track_labels)
        continue;
    end

    include_track = false;
    valid_labels = track_labels(track_labels > 0);

    switch options.filter_mode
        case 'pass_through'
            include_track = any(track_labels == region_id);

        case 'start_in'
            if ~isempty(valid_labels)
                include_track = valid_labels(1) == region_id;
            end

        case 'end_in'
            if ~isempty(valid_labels)
                include_track = valid_labels(end) == region_id;
            end

        case 'connect_to'
            if ~isempty(valid_labels)
                unique_regions = unique(valid_labels);
                include_track = ismember(region_id, unique_regions) && numel(unique_regions) > 1;
            end
    end

    if include_track
        region_points = sum(track_labels == region_id);
        if region_points == 0
            include_track = false;
        else
            total_points = numel(track_labels);
            overlap_ratio = region_points / total_points;
            starts_or_ends_in_region = ~isempty(valid_labels) && ...
                (valid_labels(1) == region_id || valid_labels(end) == region_id);

            if strcmp(options.filter_mode, 'pass_through')
                meets_overlap = (options.min_overlap <= 0) || (overlap_ratio >= options.min_overlap);
                meets_points = (options.min_region_points <= 0) || (region_points >= options.min_region_points);
                include_track = starts_or_ends_in_region || meets_overlap || meets_points;
            elseif options.min_overlap > 0
                include_track = overlap_ratio >= options.min_overlap;
            end
        end
    end

    if include_track
        track_count = track_count + 1;
        filtered_tracks{track_count} = track;
    end
end

filtered_tracks = filtered_tracks(1:track_count);
end


function track_labels = get_track_parcellation_labels(track, parcellation_mask)
% Get parcellation labels for each point along the track
track_labels = zeros(size(track, 1), 1);

mask_dims = size(parcellation_mask);

for i = 1:size(track, 1)
    % Convert from continuous coordinates to discrete voxel indices
    pos_continuous = track(i, :);

    % Ensure we have exactly 3 coordinates
    if length(pos_continuous) ~= 3
        continue;
    end

    % Round to nearest voxel and clamp to valid range
    pos = round(pos_continuous);
    pos = max(1, min(pos, mask_dims));  % Clamp to valid indices

    % Additional bounds check (should be redundant now)
    if pos(1) >= 1 && pos(1) <= mask_dims(1) && ...
       pos(2) >= 1 && pos(2) <= mask_dims(2) && ...
       pos(3) >= 1 && pos(3) <= mask_dims(3)
        track_labels(i) = parcellation_mask(pos(1), pos(2), pos(3));
    end
end
end


%% ============================================================================
%% UTILITY AND HELPER FUNCTIONS
%% ============================================================================

function region_name = get_region_name(nim, region_id)
% Get the name of a region if available
region_name = '';

% Try new direct access arrays first (preferred)
if isfield(nim, 'labels') && region_id <= length(nim.labels) && ~isempty(nim.labels{region_id})
    region_name = nim.labels{region_id};
elseif isfield(nim, 'region_names') && region_id <= length(nim.region_names) && strlength(nim.region_names(region_id)) > 0
    region_name = char(nim.region_names(region_id));
% Fall back to old map access for compatibility
elseif isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
    try
        if nim.atlas_labels.map.isKey(region_id)
            region_name = nim.atlas_labels.map(region_id);
        end
    catch
        % Continue without name
    end
end

if isempty(region_name)
    region_name = sprintf('Region %d', region_id);
end
end


function [grid_rows, grid_cols] = calculate_optimal_grid(num_regions, preferred_cols)
% Calculate optimal grid layout for displaying regions
if preferred_cols > 0
    grid_cols = preferred_cols;
    grid_rows = ceil(num_regions / grid_cols);
else
    % Auto-calculate optimal grid dimensions
    if num_regions <= 4
        grid_rows = 2; grid_cols = 2;
    elseif num_regions <= 6
        grid_rows = 2; grid_cols = 3;
    elseif num_regions <= 9
        grid_rows = 3; grid_cols = 3;
    elseif num_regions <= 12
        grid_rows = 3; grid_cols = 4;
    elseif num_regions <= 16
        grid_rows = 4; grid_cols = 4;
    elseif num_regions <= 20
        grid_rows = 4; grid_cols = 5;
    elseif num_regions <= 25
        grid_rows = 5; grid_cols = 5;
    elseif num_regions <= 30
        grid_rows = 5; grid_cols = 6;
    else
        % For very large numbers, use square-ish layout
        grid_cols = ceil(sqrt(num_regions));
        grid_rows = ceil(num_regions / grid_cols);
    end
end
end


function setup_region_subplot(nim, region_id, region_name, track_count)
% Set up individual region subplot with proper formatting
dims = size(nim.FA);
xlim([1 dims(1)]); ylim([1 dims(2)]); zlim([1 dims(3)]);
axis equal; grid on; view(45, 30);

% Truncate long region names for title
display_name = region_name;
if length(display_name) > 25
    display_name = [display_name(1:22) '...'];
end

% Create informative title
if track_count > 0
    title_str = sprintf('R%d: %s\n%d tracks', region_id, display_name, track_count);
    title_color = 'black';
else
    title_str = sprintf('R%d: %s\nNo tracks', region_id, display_name);
    title_color = [0.6, 0.6, 0.6];
end

title(title_str, 'FontSize', 10, 'FontWeight', 'bold', 'Color', title_color);

% Add small axis labels
xlabel('X', 'FontSize', 8);
ylabel('Y', 'FontSize', 8);
zlabel('Z', 'FontSize', 8);

% Add subtle lighting
camlight; lighting gouraud;
end


function add_color_legend(color_mode)
% Add legend explaining the visualization
switch color_mode
    case 'direction'
        legend_text = 'Colors: Red=L-R, Green=A-P, Blue=S-I';
    case 'fa'
        legend_text = 'Colors: Hot colormap (FA values)';
    case 'uniform'
        legend_text = 'All tracks same color';
    case 'region'
        legend_text = 'Colors: Region-based';
    otherwise
        legend_text = '';
end

if ~isempty(legend_text)
    text(0.02, 0.98, legend_text, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
         'FontSize', 9, 'EdgeColor', 'black');
end
end


%% ============================================================================
%% STATISTICAL ANALYSIS AND DISPLAY FUNCTIONS
%% ============================================================================

function plot_track_length_histogram(tracks)
% Plot track length distribution
if isfield(tracks, 'track_lengths') && ~isempty(tracks.track_lengths)
    % Use pre-computed lengths if available
    track_lengths = tracks.track_lengths;
else
    % Calculate lengths from tracks
    track_lengths = cellfun(@(x) size(x, 1), tracks);
end

if ~isempty(track_lengths)
    histogram(track_lengths, 20);
    xlabel('Track Length (points)');
    ylabel('Count');
    title('Track Length Distribution');

    % Add statistics
    mean_len = mean(track_lengths);
    hold on;
    xline(mean_len, 'r--', sprintf('Mean: %.1f', mean_len), 'LineWidth', 2);
    grid on;
else
    text(0.5, 0.5, 'No length data available', 'HorizontalAlignment', 'center');
    title('Track Lengths');
end
end


function plot_seed_distribution(tracks)
% Plot where tracks start (seed distribution)
if ~isempty(tracks)
    start_points = zeros(length(tracks), 3);
    for i = 1:length(tracks)
        if ~isempty(tracks{i}) && size(tracks{i}, 1) > 0
            start_points(i, :) = tracks{i}(1, :);
        end
    end

    % Remove invalid points
    valid_idx = all(start_points > 0, 2);
    start_points = start_points(valid_idx, :);

    % Sample points for visualization
    max_points = 1000;
    if size(start_points, 1) > max_points
        idx = randperm(size(start_points, 1), max_points);
        start_points = start_points(idx, :);
    end

    scatter3(start_points(:,1), start_points(:,2), start_points(:,3), ...
             10, 'filled', 'MarkerFaceAlpha', 0.6);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Seed Point Distribution');
    grid on; axis equal; view(3);
else
    text(0.5, 0.5, 'No tracks to display', 'HorizontalAlignment', 'center');
    title('Seed Points');
end
end


function print_track_summary(tracks, options)
% Print summary statistics (for whole brain mode)
fprintf('\n=== TRACK SUMMARY ===\n');
fprintf('Number of tracks: %d\n', length(tracks));

if ~isempty(tracks)
    track_lengths = cellfun(@(x) size(x, 1), tracks);
    fprintf('Track length (points): Mean=%.1f, Min=%d, Max=%d\n', ...
            mean(track_lengths), min(track_lengths), max(track_lengths));
    fprintf('Estimated length (mm): Mean=%.1f, Min=%.1f, Max=%.1f\n', ...
            mean(track_lengths), min(track_lengths), max(track_lengths));
end

fprintf('Color mode: %s\n', options.color_mode);
fprintf('Max tracks displayed: %d\n', options.max_tracks);
fprintf('====================\n');
end


function display_region_info(nim, region_id)
% Display information about the selected region
fprintf('\n=== REGION INFORMATION ===\n');
fprintf('Region ID: %d\n', region_id);

% Get region name
region_name = get_region_name(nim, region_id);
if ~strcmp(region_name, sprintf('Region %d', region_id))
    fprintf('Region name: %s\n', region_name);
end

% Region statistics
region_voxels = sum(nim.parcellation_mask(:) == region_id);
total_voxels = sum(nim.parcellation_mask(:) > 0);
region_percentage = (region_voxels / total_voxels) * 100;

fprintf('Region size: %d voxels (%.1f%% of brain)\n', region_voxels, region_percentage);

% Available regions info
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0);
fprintf('Total regions available: %d (IDs: %d-%d)\n', ...
        length(unique_regions), min(unique_regions), max(unique_regions));

fprintf('========================\n\n');
end


function display_track_statistics(tracks, region_id, region_name, options)
% Display statistics about the filtered tracks
fprintf('\n=== TRACK STATISTICS FOR REGION %d ===\n', region_id);
if ~strcmp(region_name, sprintf('Region %d', region_id))
    fprintf('Region name: %s\n', region_name);
end

fprintf('Filter mode: %s\n', options.filter_mode);
fprintf('Minimum overlap: %.1f%%\n', options.min_overlap * 100);
if strcmp(options.filter_mode, 'pass_through')
    fprintf('Min in-region points: %d\n', options.min_region_points);
end
fprintf('Number of tracks: %d\n', length(tracks));

if ~isempty(tracks)
    track_lengths = cellfun(@(x) size(x, 1), tracks);
    fprintf('Track length (points): Mean=%.1f, Min=%d, Max=%d\n', ...
            mean(track_lengths), min(track_lengths), max(track_lengths));

    % Estimate physical lengths (assuming 1mm spacing)
    fprintf('Estimated track length (mm): Mean=%.1f, Min=%.1f, Max=%.1f\n', ...
            mean(track_lengths), min(track_lengths), max(track_lengths));
end

fprintf('=====================================\n');
end


function display_grid_statistics(region_info, options, total_tracks)
% Display comprehensive statistics for grid visualization
fprintf('\n=== GRID VISUALIZATION STATISTICS ===\n');
fprintf('Filter mode: %s\n', options.filter_mode);
fprintf('Minimum overlap: %.1f%%\n', options.min_overlap * 100);
if strcmp(options.filter_mode, 'pass_through')
    fprintf('Min in-region points: %d\n', options.min_region_points);
end
if options.show_all_tracks
    fprintf('Track limits: NONE (ALL tracks displayed)\n');
else
    fprintf('Max tracks per region: %d\n', options.max_tracks);
end

if ~isempty(region_info)
    fprintf('\n--- REGION BREAKDOWN ---\n');
    fprintf('Region ID | Region Name                    | Tracks     | Subplot\n');
    fprintf('----------|--------------------------------|------------|--------\n');

    regions_with_tracks = 0;
    track_counts = [];

    for i = 1:length(region_info)
        info = region_info(i);

        % Truncate long names for display
        display_name = info.name;
        if length(display_name) > 30
            display_name = [display_name(1:27) '...'];
        end

        fprintf('%8d | %-30s | %10d | %7d\n', ...
                info.id, display_name, info.track_count, info.subplot_idx);

        if info.track_count > 0
            regions_with_tracks = regions_with_tracks + 1;
            track_counts(end+1) = info.track_count;
        end
    end

    fprintf('\nTOTAL: %d regions, %d with tracks, %d tracks displayed\n', ...
            length(region_info), regions_with_tracks, total_tracks);

    if ~isempty(track_counts)
        fprintf('Track distribution: Mean=%.1f, Min=%d, Max=%d\n', ...
                mean(track_counts), min(track_counts), max(track_counts));
        fprintf('Standard deviation: %.1f tracks\n', std(track_counts));
    end
else
    fprintf('No regions found.\n');
end

if options.show_all_tracks
    fprintf('\nIMPORTANT: ALL available tracks are displayed (no limits applied)\n');
end
fprintf('=====================================\n');
end


%% ============================================================================
%% IMAGE EXPORT FUNCTIONALITY
%% ============================================================================

function export_figure(fig_handle, base_filename, options)
% Export figure with automatic directory management and naming

if isempty(options.export_dir)
    return;
end

fprintf('Exporting figure...\n');

% Create main export directory
export_dir = char(options.export_dir);
if ~exist(export_dir, 'dir')
    mkdir(export_dir);
    fprintf('  Created export directory: %s\n', export_dir);
end

% Create mode-specific subdirectories
mode_dir = fullfile(export_dir, options.mode);
if ~exist(mode_dir, 'dir')
    mkdir(mode_dir);
    fprintf('  Created mode directory: %s\n', mode_dir);
end

% Create metadata directory
metadata_dir = fullfile(export_dir, 'metadata');
if ~exist(metadata_dir, 'dir')
    mkdir(metadata_dir);
end

% Generate full filename
full_filename = fullfile(mode_dir, [base_filename '.' options.export_format]);

% Export based on format
try
    switch options.export_format
        case 'png'
            print(fig_handle, full_filename, '-dpng', sprintf('-r%d', options.export_dpi));
        case 'pdf'
            print(fig_handle, full_filename, '-dpdf', '-bestfit');
        case 'eps'
            print(fig_handle, full_filename, '-depsc', sprintf('-r%d', options.export_dpi));
        case 'fig'
            savefig(fig_handle, full_filename);
    end

    fprintf('  ✓ Exported: %s\n', full_filename);

    % Log export info
    log_export_info(full_filename, options, metadata_dir);

catch ME
    warning('EXPORT:Failed', 'Failed to export figure: %s', ME.message);
end
end


function log_export_info(filename, options, metadata_dir)
% Log export information for reproducibility

log_file = fullfile(metadata_dir, 'export_log.txt');

% Create log entry
log_entry = sprintf('[%s] %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), filename);
log_entry = [log_entry sprintf('  Mode: %s, Format: %s, DPI: %d\n', ...
                               options.mode, options.export_format, options.export_dpi)];
log_entry = [log_entry sprintf('  Filter: %s (overlap %.2f, minPts %d), Color: %s, Max tracks: %d\n', ...
                               options.filter_mode, options.min_overlap, options.min_region_points, ...
                               options.color_mode, options.max_tracks)];
if strcmp(options.mode, 'region')
    log_entry = [log_entry sprintf('  Region ID: %d\n', options.region_id)];
end
log_entry = [log_entry newline];

% Append to log file
fid = fopen(log_file, 'a');
if fid > 0
    fprintf(fid, '%s', log_entry);
    fclose(fid);
end
end
