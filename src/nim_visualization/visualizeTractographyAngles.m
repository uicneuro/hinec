function visualizeTractographyAngles(varargin)
% VISUALIZETRACTOGRAPHYANGLES: Generate tractography figures from multiple viewing angles
%
% Automatically generates and exports 3D tractography visualizations from
% standard anatomical viewing angles (axial, sagittal, coronal, oblique).
%
% Usage:
%   % Auto-detect files from run directory
%   visualizeTractographyAngles('hinec_runs/run_20251125_*')
%   visualizeTractographyAngles('hinec_runs/run_20251125_*', 'output_dir', './figures')
%
%   % Traditional usage with explicit file paths
%   visualizeTractographyAngles(tracks_file, nim_file)
%   visualizeTractographyAngles(tracks_file, nim_file, 'output_dir', './figures')
%
% Arguments:
%   EITHER:
%     run_dir - Path to run directory (e.g., 'hinec_runs/run_20251125_*') - auto-detects files
%   OR:
%     tracks_file - Path to tracks .mat file
%     nim_file - Path to nim .mat file
%
% Parameters:
%   'output_dir' - Directory for exported figures (default: './tractography_angles')
%   'format' - Export format: 'png', 'pdf', 'eps', 'fig' (default: 'png')
%   'dpi' - Export resolution (default: 300)
%   'max_tracks' - Maximum tracks to display (default: 5000)
%   'line_width' - Track line width (default: 0.8)
%   'color_mode' - 'direction', 'fa', 'uniform' (default: 'direction')
%   'show_anatomy' - Show FA background slices (default: true)
%   'show_region' - Show region overlay when filtering (default: true)
%   'region_alpha' - Region overlay transparency (default: 0.3)
%   'fig_size' - Figure size [width, height] (default: [1200, 1000])
%   'angles' - Cell array of angle names to generate (default: all)
%              Options: 'superior', 'inferior', 'anterior', 'posterior',
%                       'left', 'right', 'oblique_left', 'oblique_right'
%
% Region Filtering:
%   'region' - Region ID or array of IDs to filter tracks (default: [] = all tracks)
%   'filter_mode' - 'pass_through', 'start_in', 'end_in', 'connect_to' (default: 'pass_through')
%   'min_overlap' - Minimum region overlap ratio 0-1 (default: 0.05)
%   'min_region_points' - Minimum in-region track points (default: 3)
%
% Standard Anatomical Views:
%   Superior (top-down)    - Looking down from above the head
%   Inferior (bottom-up)   - Looking up from below the head
%   Anterior (front)       - Looking at the face
%   Posterior (back)       - Looking at the back of the head
%   Left (left side)       - Looking at left ear
%   Right (right side)     - Looking at right ear
%   Oblique Left           - 3D view from front-left
%   Oblique Right          - 3D view from front-right
%
% Examples:
%   % Generate all angles from run directory
%   visualizeTractographyAngles('hinec_runs/run_*')
%
%   % Generate specific angles only
%   visualizeTractographyAngles('hinec_runs/run_*', 'angles', {'superior', 'anterior', 'left'})
%
%   % Custom output directory and format
%   visualizeTractographyAngles('hinec_runs/run_*', 'output_dir', './figs', 'format', 'pdf')
%
%   % Filter by region (single or multiple)
%   visualizeTractographyAngles('hinec_runs/run_*', 'region', 5)
%   visualizeTractographyAngles('hinec_runs/run_*', 'region', [5, 10, 15], 'filter_mode', 'start_in')

%% INPUT VALIDATION AND PARSING
if nargin < 1
    error('visualizeTractographyAngles requires at least 1 argument: run_dir OR (tracks_file, nim_file)');
end

% Check if first argument is a run directory or a tracks file
first_arg = varargin{1};
if ~ischar(first_arg) && ~isstring(first_arg)
    error('First argument must be a string path');
end
first_arg = char(first_arg);

% Detect if this is a run directory (contains 'run_' or is a directory with output/tractography)
[tracks_file, nim_file, remaining_args] = resolve_input_paths(first_arg, varargin);

% Parse input parameters
p = inputParser;
addRequired(p, 'tracks_file');
addRequired(p, 'nim_file');

% Output parameters
addParameter(p, 'output_dir', './tractography_angles', @(x) ischar(x) || isstring(x));
addParameter(p, 'format', 'png', @(x) ismember(x, {'png', 'pdf', 'eps', 'fig'}));
addParameter(p, 'dpi', 300, @(x) isnumeric(x) && x > 0);

% Visualization parameters
addParameter(p, 'max_tracks', 5000, @(x) isnumeric(x) && x > 0);
addParameter(p, 'line_width', 0.8, @(x) isnumeric(x) && x > 0);
addParameter(p, 'color_mode', 'direction', @(x) ismember(x, {'direction', 'fa', 'uniform'}));
addParameter(p, 'show_anatomy', true, @islogical);
addParameter(p, 'show_region', true, @islogical);
addParameter(p, 'region_alpha', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'fig_size', [1200, 1000], @(x) isnumeric(x) && length(x) == 2);

% Region filtering parameters
addParameter(p, 'region', [], @(x) isempty(x) || (isnumeric(x) && all(x > 0)));
addParameter(p, 'filter_mode', 'pass_through', @(x) ismember(x, {'pass_through', 'start_in', 'end_in', 'connect_to'}));
addParameter(p, 'min_overlap', 0.05, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'min_region_points', 3, @(x) isnumeric(x) && x >= 0);

% Angle selection
all_angles = {'superior', 'inferior', 'anterior', 'posterior', 'left', 'right', 'oblique_left', 'oblique_right'};
addParameter(p, 'angles', all_angles, @(x) iscell(x) || ischar(x) || isstring(x));

parse(p, tracks_file, nim_file, remaining_args{:});
options = p.Results;

% Ensure angles is a cell array
if ischar(options.angles) || isstring(options.angles)
    options.angles = {char(options.angles)};
end

%% SETUP
fprintf('=== MULTI-ANGLE TRACTOGRAPHY VISUALIZATION ===\n');
fprintf('Output directory: %s\n', options.output_dir);
fprintf('Format: %s @ %d DPI\n', options.format, options.dpi);
fprintf('Angles to generate: %s\n', strjoin(options.angles, ', '));

% Process region filtering settings
if ~isempty(options.region)
    options.region_ids = options.region;
    options.multi_region = length(options.region) > 1;
    if options.multi_region
        fprintf('Region filter: %s (regions %s)\n', options.filter_mode, mat2str(options.region_ids));
    else
        fprintf('Region filter: %s (region %d)\n', options.filter_mode, options.region_ids(1));
    end
else
    options.region_ids = [];
    options.multi_region = false;
end

% Create output directory
if ~exist(options.output_dir, 'dir')
    mkdir(options.output_dir);
    fprintf('Created output directory: %s\n', options.output_dir);
end

%% LOAD DATA
fprintf('\nLoading data...\n');
[tracks, nim] = load_data(tracks_file, nim_file);
fprintf('Loaded %d tracks\n', length(tracks));

%% FILTER BY REGION IF SPECIFIED
display_tracks = tracks;
if ~isempty(options.region_ids)
    % Validate parcellation exists
    if ~isfield(nim, 'parcellation_mask')
        error('Region filtering requires parcellation_mask in nim structure');
    end

    fprintf('Filtering tracks by region(s)...\n');
    filtered_tracks = {};

    for r = 1:length(options.region_ids)
        region_id = options.region_ids(r);
        region_tracks = filter_tracks_by_region(tracks, nim.parcellation_mask, region_id, options);
        if ~isempty(region_tracks)
            filtered_tracks = [filtered_tracks; region_tracks(:)];
            fprintf('  Region %d: %d tracks\n', region_id, length(region_tracks));
        else
            fprintf('  Region %d: No tracks found\n', region_id);
        end
    end

    if isempty(filtered_tracks)
        warning('No tracks found for specified region(s). Showing all tracks instead.');
    else
        display_tracks = filtered_tracks;
        fprintf('Total filtered tracks: %d\n', length(display_tracks));
    end
end

% Limit tracks for performance
if length(display_tracks) > options.max_tracks
    fprintf('Limiting to %d tracks for visualization\n', options.max_tracks);
    track_indices = round(linspace(1, length(display_tracks), options.max_tracks));
    display_tracks = display_tracks(track_indices);
end

tracks = display_tracks;

%% DEFINE VIEWING ANGLES
% Format: [azimuth, elevation] in degrees
% Azimuth: rotation around Z-axis (0 = +X, 90 = +Y)
% Elevation: angle above XY plane
angle_definitions = struct();
angle_definitions.superior = struct('az', 0, 'el', 90, 'name', 'Superior (Top)', 'desc', 'Looking down from above');
angle_definitions.inferior = struct('az', 0, 'el', -90, 'name', 'Inferior (Bottom)', 'desc', 'Looking up from below');
angle_definitions.anterior = struct('az', 180, 'el', 0, 'name', 'Anterior (Front)', 'desc', 'Looking at the face');
angle_definitions.posterior = struct('az', 0, 'el', 0, 'name', 'Posterior (Back)', 'desc', 'Looking at back of head');
angle_definitions.left = struct('az', 90, 'el', 0, 'name', 'Left Side', 'desc', 'Looking at left ear');
angle_definitions.right = struct('az', -90, 'el', 0, 'name', 'Right Side', 'desc', 'Looking at right ear');
angle_definitions.oblique_left = struct('az', 135, 'el', 30, 'name', 'Oblique Left', 'desc', '3D view from front-left');
angle_definitions.oblique_right = struct('az', -135, 'el', 30, 'name', 'Oblique Right', 'desc', '3D view from front-right');

%% GENERATE FIGURES FOR EACH ANGLE
dims = size(nim.FA);
num_angles = length(options.angles);
exported_files = cell(num_angles, 1);

fprintf('\nGenerating %d views...\n', num_angles);

for i = 1:num_angles
    angle_name = options.angles{i};

    if ~isfield(angle_definitions, angle_name)
        warning('Unknown angle: %s, skipping', angle_name);
        continue;
    end

    angle_def = angle_definitions.(angle_name);
    fprintf('\n[%d/%d] %s (%s)...\n', i, num_angles, angle_def.name, angle_def.desc);

    % Create figure (invisible for batch processing)
    fig = figure('Name', angle_def.name, ...
                 'Position', [100, 100, options.fig_size(1), options.fig_size(2)], ...
                 'Color', 'w', ...
                 'Visible', 'off');

    hold on;

    % Plot anatomical background if requested
    if options.show_anatomy
        plot_anatomical_background(nim, 0.12);
    end

    % Plot region overlay if filtering by region
    if options.show_region && ~isempty(options.region_ids) && isfield(nim, 'parcellation_mask')
        for r = 1:length(options.region_ids)
            plot_region_overlay(nim.parcellation_mask, options.region_ids(r), options.region_alpha);
        end
    end

    % Plot tracks with direction-based coloring
    plot_tracks_with_color(tracks, nim, options);

    % Set up axes
    xlim([1 dims(1)]);
    ylim([1 dims(2)]);
    zlim([1 dims(3)]);
    xlabel('X (L-R)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y (A-P)', 'FontSize', 12, 'FontWeight', 'bold');
    zlabel('Z (S-I)', 'FontSize', 12, 'FontWeight', 'bold');

    % Set viewing angle
    view(angle_def.az, angle_def.el);

    % Title with angle info (include region if filtering)
    if ~isempty(options.region_ids)
        if options.multi_region
            region_str = sprintf('Regions %s', mat2str(options.region_ids));
        else
            region_str = sprintf('Region %d', options.region_ids(1));
        end
        title_str = sprintf('%s View | %s\n%d tracks | %s coloring', ...
                            angle_def.name, region_str, length(tracks), options.color_mode);
    else
        title_str = sprintf('%s View\n%d tracks | %s coloring', ...
                            angle_def.name, length(tracks), options.color_mode);
    end
    title(title_str, 'FontSize', 14, 'FontWeight', 'bold');

    % Formatting
    axis equal;
    grid on;
    box on;

    % Add lighting for 3D effect (except for pure top/bottom views)
    if abs(angle_def.el) < 85
        camlight('headlight');
        lighting gouraud;
    end

    % Add color legend
    add_color_legend(options.color_mode);

    hold off;

    % Export figure (include region in filename if filtering)
    if ~isempty(options.region_ids)
        if options.multi_region
            region_suffix = sprintf('_regions-%s', strrep(mat2str(options.region_ids), ' ', '-'));
        else
            region_suffix = sprintf('_region-%02d', options.region_ids(1));
        end
        output_filename = fullfile(options.output_dir, sprintf('tractography_%s%s.%s', angle_name, region_suffix, options.format));
    else
        output_filename = fullfile(options.output_dir, sprintf('tractography_%s.%s', angle_name, options.format));
    end
    export_figure(fig, output_filename, options);
    exported_files{i} = output_filename;

    fprintf('  Exported: %s\n', output_filename);

    % Close figure to free memory
    close(fig);
end

%% SUMMARY
fprintf('\n=== EXPORT COMPLETE ===\n');
fprintf('Output directory: %s\n', options.output_dir);
fprintf('Files generated:\n');
for i = 1:length(exported_files)
    if ~isempty(exported_files{i})
        fprintf('  - %s\n', exported_files{i});
    end
end
fprintf('===========================\n');

end


%% ============================================================================
%% HELPER FUNCTIONS
%% ============================================================================

function [tracks_file, nim_file, remaining_args] = resolve_input_paths(first_arg, all_args)
% RESOLVE_INPUT_PATHS: Detect input format and resolve file paths
% (Same logic as visualizeTractography.m)

is_run_dir = false;

% Handle wildcard patterns
if contains(first_arg, '*')
    matches = dir(first_arg);
    if ~isempty(matches)
        [~, idx] = sort({matches.name});
        newest = matches(idx(end));
        first_arg = fullfile(newest.folder, newest.name);
        fprintf('Wildcard matched: %s\n', first_arg);
    else
        error('No directories match pattern: %s', first_arg);
    end
end

% Check if it's a directory with the expected structure
if isfolder(first_arg)
    output_dir = fullfile(first_arg, 'output');
    tractography_dir = fullfile(first_arg, 'tractography');

    if isfolder(output_dir) && isfolder(tractography_dir)
        is_run_dir = true;
    end
end

if is_run_dir
    fprintf('=== Auto-detecting files from run directory ===\n');
    fprintf('Run directory: %s\n', first_arg);

    % Find nim file in output/
    output_dir = fullfile(first_arg, 'output');
    nim_files = dir(fullfile(output_dir, '*.mat'));

    if isempty(nim_files)
        error('No .mat files found in output directory: %s', output_dir);
    elseif length(nim_files) > 1
        [~, idx] = max([nim_files.bytes]);
        nim_file = fullfile(output_dir, nim_files(idx).name);
    else
        nim_file = fullfile(output_dir, nim_files(1).name);
    end
    fprintf('  nim file: %s\n', nim_file);

    % Find tracks file in tractography/
    tractography_dir = fullfile(first_arg, 'tractography');
    tracks_files = dir(fullfile(tractography_dir, 'tracks*.mat'));

    if isempty(tracks_files)
        error('No tracks*.mat files found in tractography directory: %s', tractography_dir);
    elseif length(tracks_files) > 1
        [~, idx] = max([tracks_files.datenum]);
        tracks_file = fullfile(tractography_dir, tracks_files(idx).name);
    else
        tracks_file = fullfile(tractography_dir, tracks_files(1).name);
    end
    fprintf('  tracks file: %s\n', tracks_file);
    fprintf('================================================\n\n');

    remaining_args = all_args(2:end);
else
    if length(all_args) < 2
        error(['First argument is not a valid run directory.\n' ...
               'Either provide a run directory: visualizeTractographyAngles(''hinec_runs/run_*'')\n' ...
               'Or provide both files: visualizeTractographyAngles(tracks_file, nim_file)']);
    end

    tracks_file = char(all_args{1});
    nim_file = char(all_args{2});
    remaining_args = all_args(3:end);
end

tracks_file = char(tracks_file);
nim_file = char(nim_file);
end


function [tracks, nim] = load_data(tracks_file, nim_file)
% Load tracks and nim data

if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s', tracks_file);
end
if ~exist(nim_file, 'file')
    error('Nim file not found: %s', nim_file);
end

track_data = load(tracks_file);
if ~isfield(track_data, 'tracks')
    error('Invalid tracks file: missing "tracks" field');
end
tracks = track_data.tracks;

nim_data = load(nim_file);
if ~isfield(nim_data, 'nim')
    error('Invalid nim file: missing "nim" field');
end
nim = nim_data.nim;

if ~isfield(nim, 'FA')
    error('nim structure missing required field: FA');
end
end


function plot_tracks_with_color(tracks, nim, options)
% Plot tracks with specified coloring mode

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 4
        continue;
    end

    track_color = get_track_color(track, nim, options.color_mode);

    if isempty(track_color)
        continue;
    end

    plot3(track(:,1), track(:,2), track(:,3), ...
          'Color', track_color, 'LineWidth', options.line_width);
end
end


function color = get_track_color(track, nim, color_mode)
% Get color for a track based on specified mode

switch color_mode
    case 'direction'
        directions = diff(track);
        if ~isempty(directions)
            avg_direction = mean(directions, 1);
            avg_direction_norm = norm(avg_direction);
            if avg_direction_norm > 0
                color = abs(avg_direction) / avg_direction_norm;
            else
                color = [];
            end
        else
            color = [];
        end

    case 'fa'
        fa_values = sample_fa_along_track(track, nim.FA);
        avg_fa = mean(fa_values);
        cmap = hot(256);
        fa_idx = round(avg_fa * 255) + 1;
        fa_idx = max(1, min(256, fa_idx));
        color = cmap(fa_idx, :);

    case 'uniform'
        color = [0.2, 0.6, 0.8];

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
    [Y_plot, X_plot] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
    Z = ones(size(X_plot)) * s;
    surf(X_plot, Y_plot, Z, fa_slice, 'EdgeColor', 'none', 'FaceAlpha', alpha_value);
end
colormap(gray);
end


function add_color_legend(color_mode)
% Add legend explaining the color scheme

switch color_mode
    case 'direction'
        legend_text = 'Colors: Red=L-R, Green=A-P, Blue=S-I';
    case 'fa'
        legend_text = 'Colors: Hot colormap (FA values)';
    case 'uniform'
        legend_text = 'All tracks: uniform color';
    otherwise
        legend_text = '';
end

if ~isempty(legend_text)
    annotation('textbox', [0.02, 0.02, 0.3, 0.04], ...
               'String', legend_text, ...
               'FitBoxToText', 'on', ...
               'BackgroundColor', 'white', ...
               'EdgeColor', 'black', ...
               'FontSize', 10);
end
end


function export_figure(fig, filename, options)
% Export figure to file

switch options.format
    case 'png'
        print(fig, filename, '-dpng', sprintf('-r%d', options.dpi));
    case 'pdf'
        print(fig, filename, '-dpdf', '-bestfit');
    case 'eps'
        print(fig, filename, '-depsc', sprintf('-r%d', options.dpi));
    case 'fig'
        savefig(fig, filename);
end
end


function filtered_tracks = filter_tracks_by_region(tracks, parcellation_mask, region_id, options)
% Filter tracks based on their relationship to the specified region
% (Ported from visualizeTractography.m)

filtered_tracks = cell(length(tracks), 1);
track_count = 0;

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 10  % Require at least 10 points
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
    pos_continuous = track(i, :);

    if length(pos_continuous) ~= 3
        continue;
    end

    pos = round(pos_continuous);
    pos = max(1, min(pos, mask_dims));

    if pos(1) >= 1 && pos(1) <= mask_dims(1) && ...
       pos(2) >= 1 && pos(2) <= mask_dims(2) && ...
       pos(3) >= 1 && pos(3) <= mask_dims(3)
        track_labels(i) = parcellation_mask(pos(1), pos(2), pos(3));
    end
end
end


function plot_region_overlay(parcellation_mask, region_id, alpha_value)
% Plot the parcellation region as a 3D overlay
% COORDINATE FIX: isosurface returns (col, row, slice), tracks use (row, col, slice)

region_mask = parcellation_mask == region_id;

if sum(region_mask(:)) > 100
    try
        region_mask_smooth = smooth3(double(region_mask), 'box', 3);
        [faces, vertices] = isosurface(region_mask_smooth, 0.5);

        if ~isempty(faces)
            % Swap X and Y to match track coordinates
            vertices_fixed = vertices;
            vertices_fixed(:,1) = vertices(:,2);  % X_plot = row
            vertices_fixed(:,2) = vertices(:,1);  % Y_plot = col

            patch('Vertices', vertices_fixed, 'Faces', faces, ...
                  'FaceColor', 'red', 'EdgeColor', 'none', ...
                  'FaceAlpha', alpha_value, 'SpecularStrength', 0.1);
        end
    catch
        % Fallback to scatter
        [row_idx, col_idx, slice_idx] = ind2sub(size(region_mask), find(region_mask));
        if ~isempty(row_idx)
            scatter3(row_idx, col_idx, slice_idx, 10, 'red', 'filled', 'MarkerFaceAlpha', alpha_value);
        end
    end
end
end
