function nim_plot_tractography_region(tracks, nim, region_id, varargin)
% nim_plot_tractography_region: Visualize tractography tracks for a specific brain region
%
% Usage:
%   nim_plot_tractography_region(tracks, nim, region_id)
%   nim_plot_tractography_region(tracks, nim, region_id, options)
%
% Arguments:
%   tracks - Cell array of fiber tracks (from tractography)
%   nim - NIM structure with parcellation data
%   region_id - Parcellation region index (integer) to visualize
%   options - Optional struct with visualization parameters
%
% Options:
%   filter_mode - How to filter tracks ('pass_through', 'start_in', 'end_in', 'connect_to')
%   min_overlap - Minimum percentage of track that must be in region (0-1)
%   show_region - Show the parcellation region as overlay (true/false)
%   region_alpha - Transparency of region overlay (0-1)
%   track_color - Color mode ('direction', 'fa', 'uniform', 'region')
%   max_tracks - Maximum number of tracks to display for performance
%   line_width - Width of track lines
%   show_stats - Display track statistics (true/false)

% Parse input arguments
if nargin > 3 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default options
if ~isfield(options, 'filter_mode')
    options.filter_mode = 'pass_through';
end
if ~isfield(options, 'min_overlap')
    options.min_overlap = 0.1; % 10% of track must be in region
end
if ~isfield(options, 'show_region')
    options.show_region = true;
end
if ~isfield(options, 'region_alpha')
    options.region_alpha = 0.3;
end
if ~isfield(options, 'track_color')
    options.track_color = 'direction';
end
if ~isfield(options, 'max_tracks')
    options.max_tracks = inf;
end
if ~isfield(options, 'line_width')
    options.line_width = 1.2;
end
if ~isfield(options, 'show_stats')
    options.show_stats = true;
end

% Validate inputs
if ~isfield(nim, 'parcellation_mask')
    error('nim structure must contain parcellation_mask field');
end

if region_id < 1 || region_id > max(nim.parcellation_mask(:))
    error('region_id must be between 1 and %d', max(nim.parcellation_mask(:)));
end

% Check if region exists
region_voxels = sum(nim.parcellation_mask(:) == region_id);
if region_voxels == 0
    error('Region %d does not exist in parcellation mask', region_id);
end

fprintf('Filtering tracks for region %d...\n', region_id);

% Filter tracks based on region
filtered_tracks = filter_tracks_by_region(tracks, nim.parcellation_mask, region_id, options);

if isempty(filtered_tracks)
    warning('No tracks found for region %d with current filter settings', region_id);
    return;
end

fprintf('Found %d tracks related to region %d\n', length(filtered_tracks), region_id);

% Limit tracks if specified
if ~isinf(options.max_tracks) && length(filtered_tracks) > options.max_tracks
    fprintf('Limiting display to %d tracks for performance\n', options.max_tracks);
    track_indices = round(linspace(1, length(filtered_tracks), options.max_tracks));
    filtered_tracks = filtered_tracks(track_indices);
end


% Get region name if available
region_name = get_region_name(nim, region_id);

% Create figure
figure_title = sprintf('Region %d Tractography', region_id);
if ~isempty(region_name)
    figure_title = sprintf('Region %d: %s', region_id, region_name);
end

figure('Name', figure_title, 'Color', 'w', 'Position', [100, 100, 1200, 800]);
hold on;

% Show region overlay if requested
if options.show_region
    plot_region_overlay(nim.parcellation_mask, region_id, options.region_alpha);
end

% Show anatomical background (FA slices)
plot_anatomical_background(nim);

% Plot filtered tracks
plot_region_tracks(filtered_tracks, nim, options);

% Set up axes and labels
axis equal;
grid on;
xlabel('X (voxels)'); ylabel('Y (voxels)'); zlabel('Z (voxels)');
title(figure_title, 'FontSize', 14, 'FontWeight', 'bold');

% Set reasonable view
dims = size(nim.FA);
xlim([1 dims(1)]); ylim([1 dims(2)]); zlim([1 dims(3)]);
view(45, 30);
camlight; lighting gouraud;

% Add legend and information
add_visualization_legend(options.track_color, region_name);

% Display statistics if requested
if options.show_stats
    display_track_statistics(filtered_tracks, region_id, region_name, options);
end

hold off;
fprintf('Visualization complete for region %d\n', region_id);
end

function filtered_tracks = filter_tracks_by_region(tracks, parcellation_mask, region_id, options)
% Filter tracks based on their relationship to the specified region
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

function plot_region_overlay(parcellation_mask, region_id, alpha_value)
% Plot the parcellation region as a 3D overlay
region_mask = parcellation_mask == region_id;

% Create isosurface of the region
if sum(region_mask(:)) > 100 % Only if region has sufficient voxels
    try
        % Smooth the mask slightly for better visualization
        region_mask_smooth = smooth3(double(region_mask), 'box', 3);
        
        % Create isosurface
        [faces, vertices] = isosurface(region_mask_smooth, 0.5);
        
        if ~isempty(faces)
            % COORDINATE FIX: Reorder vertices to match track coordinate system
            % MATLAB isosurface uses (Y,X,Z) ordering, but tracks use (X,Y,Z)
            % So we need to swap X and Y coordinates
            vertices_fixed = vertices;
            vertices_fixed(:,1) = vertices(:,2); % X = old Y
            vertices_fixed(:,2) = vertices(:,1); % Y = old X
            % Z stays the same
            
            % Plot region surface with fixed coordinates
            patch('Vertices', vertices_fixed, 'Faces', faces, ...
                  'FaceColor', 'red', 'EdgeColor', 'none', ...
                  'FaceAlpha', alpha_value, 'SpecularStrength', 0.1);
        end
    catch
        % Fallback to scatter plot of region voxels
        [i, j, k] = ind2sub(size(region_mask), find(region_mask));
        if ~isempty(i)
            % COORDINATE FIX: Match track coordinate system
            % ind2sub returns (i,j,k) but tracks use (x,y,z) where x=j, y=i, z=k
            scatter3(j, i, k, 10, 'red', 'filled', 'MarkerFaceAlpha', alpha_value);
        end
    end
end
end

function plot_anatomical_background(nim)
% Plot FA slices as anatomical background
dims = size(nim.FA);
slice_step = max(1, round(dims(3) / 8));
slices = slice_step:slice_step:dims(3)-slice_step;

for s = slices
    fa_slice = nim.FA(:, :, s);
    % COORDINATE FIX: Create meshgrid to match track coordinate system
    % Tracks use (X,Y,Z) where X corresponds to columns, Y to rows
    [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
    Z = ones(size(X)) * s;
    
    % Use fa_slice directly (not transposed) to match coordinate system
    surf(X, Y, Z, fa_slice, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
end
colormap(gray);
end

function plot_region_tracks(tracks, nim, options)
% Plot the filtered tracks with specified coloring
for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 2
        continue;
    end
    
    % Get track color based on mode
    track_color = get_track_color(track, nim, options.track_color);
    
    % Plot track
    plot3(track(:,1), track(:,2), track(:,3), ...
          'Color', track_color, 'LineWidth', options.line_width);
end
end

function color = get_track_color(track, nim, color_mode)
% Get color for a track based on the specified mode
switch color_mode
    case 'direction'
        % Direction-based coloring (RGB from normalized direction vector)
        directions = diff(track);
        if ~isempty(directions)
            avg_direction = mean(directions, 1);
            avg_direction = avg_direction / norm(avg_direction);
            color = abs(avg_direction);
        else
            color = [0.5, 0.5, 0.5];
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
        % Random color based on most common region
        if isfield(nim, 'parcellation_mask')
            track_labels = get_track_parcellation_labels(track, nim.parcellation_mask);
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                mode_label = mode(valid_labels);
                rng(mode_label); % Consistent color for same region
                color = rand(1, 3) * 0.8 + 0.2;
            else
                color = [0.5, 0.5, 0.5];
            end
        else
            color = [0.5, 0.5, 0.5];
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

function add_visualization_legend(color_mode, region_name)
% Add legend explaining the visualization
legend_text = {};

switch color_mode
    case 'direction'
        legend_text{end+1} = 'Colors: Red=L-R, Green=A-P, Blue=S-I';
    case 'fa'
        legend_text{end+1} = 'Colors: Hot colormap (FA values)';
    case 'uniform'
        legend_text{end+1} = 'All tracks same color';
    case 'region'
        legend_text{end+1} = 'Colors: Region-based';
end

if ~isempty(region_name) && ~strcmp(region_name, sprintf('Region %d', 1))
    legend_text{end+1} = sprintf('Region: %s', region_name);
end

if ~isempty(legend_text)
    text(0.02, 0.98, legend_text, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'BackgroundColor', 'white', ...
         'FontSize', 9, 'EdgeColor', 'black');
end
end

function display_track_statistics(tracks, region_id, region_name, options)
% Display statistics about the filtered tracks
fprintf('\n=== TRACK STATISTICS FOR REGION %d ===\n', region_id);
if ~isempty(region_name) && ~strcmp(region_name, sprintf('Region %d', region_id))
    fprintf('Region name: %s\n', region_name);
end

fprintf('Filter mode: %s\n', options.filter_mode);
fprintf('Minimum overlap: %.1f%%\n', options.min_overlap * 100);
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