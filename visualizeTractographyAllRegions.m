function visualizeTractographyAllRegions(tracks_file, nim_file, varargin)
% visualizeTractographyAllRegions: Visualize ALL tractography tracks for each brain region in grid layout
%
% Usage:
%   visualizeTractographyAllRegions('tractography_results/tracks_*.mat', 'sample_parcellated.mat')
%   visualizeTractographyAllRegions(tracks_file, nim_file, 'filter_mode', 'start_in')
%
% Arguments:
%   tracks_file - Path to tracks .mat file (REQUIRED)
%   nim_file - Path to .mat file with nim structure (REQUIRED)
%
% Optional name-value pairs:
%   'filter_mode' - 'pass_through', 'start_in', 'end_in' (default: 'pass_through')
%   'min_overlap' - Minimum overlap with region (0-1, default: 0.1)
%   'show_region_overlays' - Show region boundaries (true/false, default: true)
%   'show_background' - Show FA background (true/false, default: true)
%   'region_alpha' - Transparency of region overlays (default: 0.3)
%   'line_width' - Width of track lines (default: 1.2)
%   'grid_cols' - Number of columns in grid (default: auto-calculated)
%   'subplot_size' - Size of each subplot (default: auto-calculated)
%   'show_all_tracks' - Force display of ALL tracks (default: true)

% Validate required arguments
if nargin < 2
    error('visualizeTractographyAllRegions requires 2 arguments: tracks_file, nim_file');
end

if ~ischar(tracks_file) && ~isstring(tracks_file)
    error('tracks_file must be a string path to .mat file');
end

if ~ischar(nim_file) && ~isstring(nim_file)
    error('nim_file must be a string path to .mat file');
end

% Parse optional arguments
p = inputParser;
addParameter(p, 'filter_mode', 'pass_through', @ischar);
addParameter(p, 'min_overlap', 0.1, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'show_region_overlays', true, @islogical);
addParameter(p, 'show_background', true, @islogical);
addParameter(p, 'region_alpha', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'line_width', 1.2, @(x) isnumeric(x) && x > 0);
addParameter(p, 'grid_cols', 0, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'subplot_size', [400, 350], @(x) isnumeric(x) && length(x) == 2);
addParameter(p, 'show_all_tracks', true, @islogical);

parse(p, varargin{:});
options = p.Results;

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

fprintf('=== Visualizing ALL Regions Tractography in Grid Layout ===\n');

%% Load data files

% Load tracks data
fprintf('Loading tracks from %s...\n', tracks_file);
if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s\nPlease run tractography first.', tracks_file);
end
track_data = load(tracks_file);
tracks = track_data.tracks;

% Load nim data
fprintf('Loading nim data from %s...\n', nim_file);
if ~exist(nim_file, 'file')
    error('Nim file not found: %s\nPlease run main() first.', nim_file);
end
nim_data = load(nim_file);
nim = nim_data.nim;

%% Validate data
if ~isfield(nim, 'parcellation_mask')
    error('nim structure must contain parcellation_mask. Please run parcellation first.');
end

if isempty(tracks)
    error('No tracks found in tracks file. Please run tractography first.');
end

%% Get all available regions
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0); % Exclude background

if isempty(unique_regions)
    error('No regions found in parcellation mask');
end

fprintf('Found %d regions to visualize: %s\n', length(unique_regions), ...
        mat2str(unique_regions'));

%% Calculate optimal grid layout
[grid_rows, grid_cols] = calculate_optimal_grid(length(unique_regions), options.grid_cols);
fprintf('Using %dx%d grid layout for %d regions\n', grid_rows, grid_cols, length(unique_regions));

%% Calculate figure size based on grid and subplot size
figure_width = grid_cols * options.subplot_size(1) + 100;
figure_height = grid_rows * options.subplot_size(2) + 150;
figure_title = sprintf('All Regions Tractography Grid (%d regions)', length(unique_regions));
figure('Name', figure_title, 'Color', 'w', 'Position', [50, 50, figure_width, figure_height]);

%% Process each region in grid layout
fprintf('Processing regions and creating grid subplots...\n');
region_info = struct('id', {}, 'name', {}, 'track_count', {}, 'subplot_idx', {});
total_tracks_displayed = 0;

for r = 1:length(unique_regions)
    region_id = unique_regions(r);
    
    % Create subplot for this region
    subplot_idx = r;
    subplot(grid_rows, grid_cols, subplot_idx);
    hold on;
    
    % Get region name
    region_name = get_region_name(nim, region_id);
    
    % Filter tracks for this region - NO LIMITS, show ALL tracks
    region_options = struct();
    region_options.filter_mode = options.filter_mode;
    region_options.min_overlap = options.min_overlap;
    region_options.max_tracks = inf; % Show ALL tracks
    
    filtered_tracks = filter_tracks_by_region(tracks, nim.parcellation_mask, region_id, region_options);
    
    if ~isempty(filtered_tracks)
        % Show anatomical background for this region if requested
        if options.show_background
            plot_anatomical_background_subplot(nim, region_id);
        end
        
        % Show region overlay if requested
        if options.show_region_overlays
            plot_region_overlay_subplot(nim.parcellation_mask, region_id, options.region_alpha);
        end
        
        % Plot ALL tracks for this region with direction-based coloring
        plot_all_region_tracks(filtered_tracks, nim, options.line_width);
        
        % Store region info
        region_info(end+1).id = region_id;
        region_info(end).name = region_name;
        region_info(end).track_count = length(filtered_tracks);
        region_info(end).subplot_idx = subplot_idx;
        
        total_tracks_displayed = total_tracks_displayed + length(filtered_tracks);
        
        fprintf('Region %d (%s): %d tracks (ALL displayed)\n', region_id, region_name, length(filtered_tracks));
    else
        fprintf('Region %d (%s): No tracks found\n', region_id, region_name);
        % Still store info for empty regions
        region_info(end+1).id = region_id;
        region_info(end).name = region_name;
        region_info(end).track_count = 0;
        region_info(end).subplot_idx = subplot_idx;
    end
    
    % Set up individual subplot
    setup_region_subplot(nim, region_id, region_name, length(filtered_tracks));
    
    hold off;
end

%% Add overall figure title and statistics
sgtitle(sprintf('%s - Total: %d tracks from %d regions (ALL TRACKS SHOWN)', ...
        figure_title, total_tracks_displayed, length(unique_regions)), ...
        'FontSize', 16, 'FontWeight', 'bold');

% Display comprehensive statistics
display_grid_statistics(region_info, options, total_tracks_displayed);

fprintf('Grid visualization complete - ALL tracks displayed!\n');
fprintf('========================================\n');
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


function plot_anatomical_background_subplot(nim, region_id)
% Plot FA slices as anatomical background for individual subplot
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


function plot_region_overlay_subplot(parcellation_mask, region_id, alpha_value)
% Plot the parcellation region as a 3D overlay for individual subplot
region_mask = parcellation_mask == region_id;

if sum(region_mask(:)) > 50 % Only if region has sufficient voxels
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
            
            % Plot region surface with region-specific color
            region_color = [1.0, 0.2, 0.2]; % Red for region boundary
            patch('Vertices', vertices_fixed, 'Faces', faces, ...
                  'FaceColor', region_color, 'EdgeColor', 'none', ...
                  'FaceAlpha', alpha_value, 'SpecularStrength', 0.1);
        end
    catch
        % Fallback to scatter plot of region voxels
        [i, j, k] = ind2sub(size(region_mask), find(region_mask));
        if ~isempty(i)
            scatter3(j, i, k, 8, 'red', 'filled', 'MarkerFaceAlpha', alpha_value);
        end
    end
end
end


function plot_all_region_tracks(tracks, nim, line_width)
% Plot ALL tracks for a region with direction-based coloring - NO LIMITS
fprintf('  Plotting %d tracks (ALL)...\n', length(tracks));

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) >= 2
        % Direction-based color
        directions = diff(track);
        if ~isempty(directions)
            avg_direction = mean(directions, 1);
            avg_direction_norm = norm(avg_direction);
            if avg_direction_norm > 0
                color = abs(avg_direction) / avg_direction_norm;
            else
                color = [0.5, 0.5, 0.5];
            end
        else
            color = [0.5, 0.5, 0.5];
        end
        
        plot3(track(:,1), track(:,2), track(:,3), ...
              'Color', color, 'LineWidth', line_width);
    end
end
end


function setup_region_subplot(nim, region_id, region_name, track_count)
% Set up individual region subplot with proper formatting
dims = size(nim.FA);
xlim([1 dims(1)]); ylim([1 dims(2)]); zlim([1 dims(3)]);
axis equal;
grid on;
view(45, 30);

% Truncate long region names for title
display_name = region_name;
if length(display_name) > 25
    display_name = [display_name(1:22) '...'];
end

% Create informative title
if track_count > 0
    title_str = sprintf('R%d: %s\n%d tracks (ALL)', region_id, display_name, track_count);
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


function display_grid_statistics(region_info, options, total_tracks)
% Display comprehensive statistics for grid visualization
fprintf('\n=== GRID VISUALIZATION STATISTICS (ALL TRACKS SHOWN) ===\n');
fprintf('Filter mode: %s\n', options.filter_mode);
fprintf('Minimum overlap: %.1f%%\n', options.min_overlap * 100);
fprintf('Track limits: NONE (ALL tracks displayed)\n');

if ~isempty(region_info)
    fprintf('\n--- COMPLETE REGION BREAKDOWN ---\n');
    fprintf('Region ID | Region Name                    | ALL Tracks | Subplot\n');
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
    
    fprintf('\nTOTAL: %d regions, %d with tracks, %d tracks displayed (ALL)\n', ...
            length(region_info), regions_with_tracks, total_tracks);
    
    if ~isempty(track_counts)
        fprintf('Track distribution: Mean=%.1f, Min=%d, Max=%d\n', ...
                mean(track_counts), min(track_counts), max(track_counts));
        fprintf('Standard deviation: %.1f tracks\n', std(track_counts));
    end
else
    fprintf('No regions found.\n');
end

fprintf('\nIMPORTANT: ALL available tracks are displayed (no limits applied)\n');
fprintf('=======================================================\n');
end


function colormap_regions = generate_region_colormap(num_regions)
% Generate distinct colors for each region using HSV colorspace
if num_regions <= 12
    % Use predefined distinct colors for small number of regions
    distinct_colors = [
        1.0, 0.0, 0.0;  % Red
        0.0, 1.0, 0.0;  % Green
        0.0, 0.0, 1.0;  % Blue
        1.0, 1.0, 0.0;  % Yellow
        1.0, 0.0, 1.0;  % Magenta
        0.0, 1.0, 1.0;  % Cyan
        1.0, 0.5, 0.0;  % Orange
        0.5, 0.0, 1.0;  % Purple
        0.0, 0.5, 0.0;  % Dark Green
        0.5, 0.5, 0.0;  % Olive
        0.8, 0.2, 0.6;  % Pink
        0.2, 0.8, 0.6;  % Teal
    ];
    colormap_regions = distinct_colors(1:num_regions, :);
else
    % Generate colors using HSV for larger number of regions
    hue_values = linspace(0, 1, num_regions + 1);
    hue_values = hue_values(1:end-1); % Remove duplicate
    saturation = 0.8 + 0.2 * rand(num_regions, 1); % Random saturation 0.8-1.0
    value = 0.7 + 0.3 * rand(num_regions, 1);       % Random value 0.7-1.0
    
    hsv_colors = [hue_values', saturation, value];
    colormap_regions = hsv2rgb(hsv_colors);
end
end


% Old single-window functions removed - now using grid layout system


% Include helper functions from nim_plot_tractography_region.m
function filtered_tracks = filter_tracks_by_region(tracks, parcellation_mask, region_id, options)
% Filter tracks based on their relationship to the specified region
% MODIFIED: No artificial limits - shows ALL qualifying tracks
filtered_tracks = cell(length(tracks), 1);
track_count = 0;

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 2
        continue;
    end
    
    % Get parcellation labels along the track
    track_labels = get_track_parcellation_labels(track, parcellation_mask);
    
    % Apply filtering based on mode
    include_track = false;
    
    switch options.filter_mode
        case 'pass_through'
            % Track passes through the region at any point
            include_track = any(track_labels == region_id);
            
        case 'start_in'
            % Track starts in the region
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                include_track = valid_labels(1) == region_id;
            end
            
        case 'end_in'
            % Track ends in the region
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                include_track = valid_labels(end) == region_id;
            end
            
        case 'connect_to'
            % Track connects the region to any other region
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                unique_regions = unique(valid_labels);
                include_track = ismember(region_id, unique_regions) && length(unique_regions) > 1;
            end
    end
    
    % Apply minimum overlap constraint
    if include_track && options.min_overlap > 0
        region_points = sum(track_labels == region_id);
        total_points = length(track_labels);
        overlap_ratio = region_points / total_points;
        include_track = overlap_ratio >= options.min_overlap;
    end
    
    if include_track
        track_count = track_count + 1;
        filtered_tracks{track_count} = track;
        
        % REMOVED: No max tracks limit - show ALL qualifying tracks
    end
end

% Trim to actual size
filtered_tracks = filtered_tracks(1:track_count);
end


function track_labels = get_track_parcellation_labels(track, parcellation_mask)
% Get parcellation labels for each point along the track
track_labels = zeros(size(track, 1), 1);

for i = 1:size(track, 1)
    pos = round(track(i, :));
    
    % Check bounds
    if all(pos >= 1) && all(pos <= size(parcellation_mask))
        track_labels(i) = parcellation_mask(pos(1), pos(2), pos(3));
    end
end
end


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
        % Continue without name if there's an error
    end
end

if isempty(region_name)
    region_name = sprintf('Region %d', region_id);
end
end