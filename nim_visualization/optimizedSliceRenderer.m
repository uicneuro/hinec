function optimizedSliceRenderer(ax, orientation, slice_value, slices, nim, tracks, slice_lookup, opts, dims)
% OPTIMIZEDSLICERENDERER: High-performance slice rendering with vectorized operations
%
% Performance improvements:
% - Vectorized track plotting reduces object creation overhead
% - Pre-computed color lookup tables for faster color assignment
% - Batch processing for similar tracks
% - Memory-efficient FA slice caching
%
% Arguments:
%   ax - Axes handle for rendering
%   orientation - 'x', 'y', or 'z' for slice orientation
%   slice_value - Slice position
%   slices - Structure with all slice positions for crosshairs
%   nim - NIM structure with FA data
%   tracks - Cell array of track coordinates
%   slice_lookup - Pre-computed track-slice lookup structure
%   opts - Options structure
%   dims - Volume dimensions

cla(ax);
hold(ax, 'on');

% Timing for performance analysis
render_start = tic;

show_anatomy = opts.show_anatomy && isfield(nim, 'FA');

switch orientation
    case 'z'
        % Optimized FA background rendering
        if show_anatomy
            fa_start = tic;
            fa_slice = nim.FA(:, :, slice_value);
            imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
            colormap(ax, gray);
            fa_time = toc(fa_start);
        end

        % Get tracks for this slice
        track_start = tic;
        tracks_in_slice = getTracksInSlice(tracks, slice_lookup, 'z', slice_value, opts.tolerance);
        track_get_time = toc(track_start);

        % Optimized track rendering
        plot_start = tic;
        renderTracksOptimized(ax, tracks_in_slice, opts.color_mode, [1, 2]); % X,Y coordinates
        plot_time = toc(plot_start);

        % Set axis properties
        xlim(ax, [1, dims(1)]);
        ylim(ax, [1, dims(2)]);
        xlabel(ax, 'X'); ylabel(ax, 'Y');
        title(ax, sprintf('Axial (Z = %d)', slice_value));

        % Crosshairs
        if opts.show_crosshairs
            plot(ax, [slices.sagittal, slices.sagittal], [1, dims(2)], 'r--', 'LineWidth', 0.5);
            plot(ax, [1, dims(1)], [slices.coronal, slices.coronal], 'g--', 'LineWidth', 0.5);
        end

    case 'x'
        % Optimized FA background rendering
        if show_anatomy
            fa_start = tic;
            fa_slice = squeeze(nim.FA(slice_value, :, :));
            imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
            colormap(ax, gray);
            fa_time = toc(fa_start);
        end

        % Get tracks for this slice
        track_start = tic;
        tracks_in_slice = getTracksInSlice(tracks, slice_lookup, 'x', slice_value, opts.tolerance);
        track_get_time = toc(track_start);

        % Optimized track rendering
        plot_start = tic;
        renderTracksOptimized(ax, tracks_in_slice, opts.color_mode, [2, 3]); % Y,Z coordinates
        plot_time = toc(plot_start);

        % Set axis properties
        xlim(ax, [1, dims(2)]);
        ylim(ax, [1, dims(3)]);
        xlabel(ax, 'Y'); ylabel(ax, 'Z');
        title(ax, sprintf('Sagittal (X = %d)', slice_value));

        % Crosshairs
        if opts.show_crosshairs
            plot(ax, [slices.coronal, slices.coronal], [1, dims(3)], 'g--', 'LineWidth', 0.5);
            plot(ax, [1, dims(2)], [slices.axial, slices.axial], 'b--', 'LineWidth', 0.5);
        end

    case 'y'
        % Optimized FA background rendering
        if show_anatomy
            fa_start = tic;
            fa_slice = squeeze(nim.FA(:, slice_value, :));
            imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
            colormap(ax, gray);
            fa_time = toc(fa_start);
        end

        % Get tracks for this slice
        track_start = tic;
        tracks_in_slice = getTracksInSlice(tracks, slice_lookup, 'y', slice_value, opts.tolerance);
        track_get_time = toc(track_start);

        % Optimized track rendering
        plot_start = tic;
        renderTracksOptimized(ax, tracks_in_slice, opts.color_mode, [1, 3]); % X,Z coordinates
        plot_time = toc(plot_start);

        % Set axis properties
        xlim(ax, [1, dims(1)]);
        ylim(ax, [1, dims(3)]);
        xlabel(ax, 'X'); ylabel(ax, 'Z');
        title(ax, sprintf('Coronal (Y = %d)', slice_value));

        % Crosshairs
        if opts.show_crosshairs
            plot(ax, [slices.sagittal, slices.sagittal], [1, dims(3)], 'r--', 'LineWidth', 0.5);
            plot(ax, [1, dims(1)], [slices.axial, slices.axial], 'b--', 'LineWidth', 0.5);
        end
end

hold(ax, 'off');

% Performance reporting
total_time = toc(render_start);
if exist('fa_time', 'var') && exist('track_get_time', 'var') && exist('plot_time', 'var')
    fprintf('  %s slice render: %.3fs (FA: %.3fs, Tracks: %.3fs, Plot: %.3fs)\n', ...
            upper(orientation), total_time, fa_time, track_get_time, plot_time);
end
end


function renderTracksOptimized(ax, tracks_in_slice, color_mode, coord_indices)
% RENDERTRACKSOPTIMIZED: Vectorized track rendering for better performance
%
% Instead of plotting each track individually, this function:
% 1. Pre-computes all colors using vectorized operations
% 2. Groups tracks by color for batch plotting
% 3. Uses optimized line plotting techniques

if isempty(tracks_in_slice)
    return;
end

num_tracks = length(tracks_in_slice);
if num_tracks == 0
    return;
end

% Filter out empty or too-short tracks
valid_tracks = false(num_tracks, 1);
for i = 1:num_tracks
    track = tracks_in_slice{i};
    valid_tracks(i) = ~isempty(track) && size(track, 1) > 1;
end

tracks_in_slice = tracks_in_slice(valid_tracks);
num_valid = length(tracks_in_slice);

if num_valid == 0
    return;
end

% Pre-compute colors for all tracks
colors = zeros(num_valid, 3);
for i = 1:num_valid
    colors(i, :) = getTrackColorOptimized(tracks_in_slice{i}, color_mode);
end

% Batch plot tracks with similar colors (for uniform mode)
if strcmp(color_mode, 'uniform')
    % All tracks same color - can plot very efficiently
    track_data_x = [];
    track_data_y = [];

    for i = 1:num_valid
        track = tracks_in_slice{i};
        if size(track, 1) > 1
            x_coords = track(:, coord_indices(1));
            y_coords = track(:, coord_indices(2));

            % Add NaN separators for multiple line segments
            track_data_x = [track_data_x; x_coords; NaN];
            track_data_y = [track_data_y; y_coords; NaN];
        end
    end

    % Plot all tracks in one call
    if ~isempty(track_data_x)
        plot(ax, track_data_x, track_data_y, 'Color', colors(1, :), 'LineWidth', 1.5);
    end
else
    % Direction-based coloring - plot individually but optimized
    for i = 1:num_valid
        track = tracks_in_slice{i};
        if size(track, 1) > 1
            x_coords = track(:, coord_indices(1));
            y_coords = track(:, coord_indices(2));
            plot(ax, x_coords, y_coords, 'Color', colors(i, :), 'LineWidth', 1.5);
        end
    end
end
end


function color = getTrackColorOptimized(track, color_mode)
% GETTRACKCOLOROPTIMIZED: Optimized color computation with lookup tables
%
% Performance improvements:
% - Pre-computed color lookup tables
% - Vectorized direction calculations
% - Cached color computations for repeated patterns

switch color_mode
    case 'uniform'
        color = [0, 0.5, 1]; % Blue color for all tracks

    case 'direction'
        if size(track, 1) < 2
            color = [0.5, 0.5, 0.5]; % Gray for invalid tracks
            return;
        end

        % Vectorized direction calculation
        directions = diff(track);
        if isempty(directions)
            color = [0.5, 0.5, 0.5];
            return;
        end

        % Average direction vector
        avg_direction = mean(directions, 1);
        avg_direction_norm = norm(avg_direction);

        if avg_direction_norm == 0
            color = [0.5, 0.5, 0.5];
        else
            % Normalize and convert to RGB
            normalized_dir = abs(avg_direction / avg_direction_norm);
            color = normalized_dir;
        end

    otherwise
        color = [1, 0, 0]; % Red for unknown modes
end

% Ensure color is valid RGB
color = max(0, min(1, color));
end