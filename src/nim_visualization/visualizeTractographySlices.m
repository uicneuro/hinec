function visualizeTractographySlices(tracks_file, nim_file, x, y, z, varargin)
% VISUALIZETRACTOGRAPHYSLICES: Command-line tractography slice viewer
%
% Usage:
%   visualizeTractographySlices('tracks.mat', 'nim.mat', x, y, z)
%   visualizeTractographySlices(..., 'save', 'output.png')
%   visualizeTractographySlices(..., 'tolerance', 2)
%   visualizeTractographySlices(..., 'show_crosshairs', true)
%
% Arguments:
%   tracks_file - Path to tracks .mat file
%   nim_file    - Path to nim structure .mat file
%   x           - Sagittal slice position (X coordinate)
%   y           - Coronal slice position (Y coordinate)
%   z           - Axial slice position (Z coordinate)
%
% Optional Parameters:
%   'save'      - Output PNG filename (if not specified, just displays)
%   'tolerance' - Slice thickness in voxels (default: 2)
%   'show_crosshairs' - Show slice intersections (default: true)
%   'show_anatomy'    - Show FA background (default: true)
%   'color_mode'      - 'direction' or 'uniform' (default: 'direction')
%   'alpha'           - FA background transparency (default: 0.6)
%   'regions'         - Region ID(s) to filter tracks (single number or array)
%   'region_filter'   - Filter mode: 'pass_through', 'start_in', 'end_in' (default: 'pass_through')
%   'min_overlap'     - Minimum overlap with region (0-1, default: 0.1)

%% Parse input arguments
p = inputParser;
addRequired(p, 'tracks_file', @(x) ischar(x) || isstring(x));
addRequired(p, 'nim_file', @(x) ischar(x) || isstring(x));
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'z', @isnumeric);
addParameter(p, 'save', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'tolerance', 2, @isnumeric);
addParameter(p, 'show_crosshairs', true, @islogical);
addParameter(p, 'color_mode', 'direction', @(x) ismember(x, {'direction', 'uniform'}));
addParameter(p, 'show_anatomy', true, @islogical);
addParameter(p, 'alpha', 0.6, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'regions', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'region_filter', 'pass_through', @(x) ismember(x, {'pass_through', 'start_in', 'end_in'}));
addParameter(p, 'min_overlap', 0.1, @(x) isnumeric(x) && x >= 0 && x <= 1);

parse(p, tracks_file, nim_file, x, y, z, varargin{:});
opts = p.Results;

%% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

fprintf('=== TRACTOGRAPHY SLICE VIEWER ===\n');
fprintf('Loading data...\n');

%% Load data
if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s', tracks_file);
end
track_data = load(tracks_file);
if isfield(track_data, 'tracks')
    tracks = track_data.tracks;
else
    error('No tracks field found in tracks file');
end

if ~exist(nim_file, 'file')
    error('NIM file not found: %s', nim_file);
end
nim_data = load(nim_file);
if isfield(nim_data, 'nim')
    nim = nim_data.nim;
else
    error('No nim field found in nim file');
end

if ~isfield(nim, 'FA')
    error('FA field required in nim structure');
end

dims = size(nim.FA);
fprintf('Loaded %d tracks | Volume dimensions: %d x %d x %d\n', ...
        length(tracks), dims(1), dims(2), dims(3));

%% Filter tracks by region if specified
if ~isempty(opts.regions)
    if ~isfield(nim, 'parcellation_mask')
        error('Parcellation mask required for region filtering. Run parcellation first.');
    end

    fprintf('Filtering tracks for region(s): %s\n', mat2str(opts.regions));
    original_count = length(tracks);

    % Filter tracks based on parcellation regions
    filtered_tracks = filterTracksByRegion(tracks, nim.parcellation_mask, opts.regions, ...
                                         opts.region_filter, opts.min_overlap);
    tracks = filtered_tracks;

    fprintf('Region filtering: %d → %d tracks (%.1f%% retained)\n', ...
            original_count, length(tracks), 100*length(tracks)/original_count);
end

% Validate slice coordinates
x = max(1, min(dims(1), round(x)));
y = max(1, min(dims(2), round(y)));
z = max(1, min(dims(3), round(z)));

slices = struct();
slices.sagittal = x;  % X coordinate
slices.coronal = y;   % Y coordinate
slices.axial = z;     % Z coordinate

fprintf('Rendering slices at X=%d, Y=%d, Z=%d\n', x, y, z);

%% Performance Profiling Setup
profiling_enabled = true; % Set to false to disable profiling
if profiling_enabled
    profile_times = struct();
    total_start_time = tic;
    fprintf('=== PERFORMANCE PROFILING ENABLED ===\n');
end

%% Precompute track slice lookup for efficiency
fprintf('Computing track intersections...\n');
if profiling_enabled, lookup_start = tic; end

% Use optimized lookup function for better performance
try
    slice_lookup = buildOptimizedTrackSliceLookup(tracks, dims);
    fprintf('Using optimized track-slice lookup\n');
catch ME
    fprintf('Warning: Optimized lookup failed (%s), falling back to standard method\n', ME.message);
    slice_lookup = buildTrackSliceLookup(tracks, dims);
end

if profiling_enabled
    profile_times.track_lookup = toc(lookup_start);
    fprintf('Track lookup time: %.3f seconds\n', profile_times.track_lookup);
end

%% Create figure
should_save = ~isempty(opts.save);
fig_visible = 'on';
if should_save
    fig_visible = 'off';  % Don't show if just saving
end

render_fig = figure('Name', sprintf('Tractography Slices [X=%d Y=%d Z=%d]', x, y, z), ...
                    'Color', 'w', 'NumberTitle', 'off', ...
                    'Visible', fig_visible, ...
                    'Position', [100, 100, 820, 740]);

%% Create layout
t = tiledlayout(render_fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Render all three slice views
if profiling_enabled, render_start = tic; end

ax_axial = nexttile(t, 1);
if profiling_enabled, axial_start = tic; end
% Use optimized renderer for better performance
try
    optimizedSliceRenderer(ax_axial, 'z', slices.axial, slices, nim, tracks, slice_lookup, opts, dims);
catch ME
    fprintf('Warning: Optimized renderer failed (%s), using fallback\n', ME.message);
    drawSlice(ax_axial, 'z', slices.axial, slices, nim, tracks, slice_lookup, opts, dims);
end
if profiling_enabled
    profile_times.axial_render = toc(axial_start);
    fprintf('Axial slice render time: %.3f seconds\n', profile_times.axial_render);
end

ax_sagittal = nexttile(t, 2);
if profiling_enabled, sagittal_start = tic; end
try
    optimizedSliceRenderer(ax_sagittal, 'x', slices.sagittal, slices, nim, tracks, slice_lookup, opts, dims);
catch ME
    fprintf('Warning: Optimized renderer failed (%s), using fallback\n', ME.message);
    drawSlice(ax_sagittal, 'x', slices.sagittal, slices, nim, tracks, slice_lookup, opts, dims);
end
if profiling_enabled
    profile_times.sagittal_render = toc(sagittal_start);
    fprintf('Sagittal slice render time: %.3f seconds\n', profile_times.sagittal_render);
end

ax_coronal = nexttile(t, 3);
if profiling_enabled, coronal_start = tic; end
try
    optimizedSliceRenderer(ax_coronal, 'y', slices.coronal, slices, nim, tracks, slice_lookup, opts, dims);
catch ME
    fprintf('Warning: Optimized renderer failed (%s), using fallback\n', ME.message);
    drawSlice(ax_coronal, 'y', slices.coronal, slices, nim, tracks, slice_lookup, opts, dims);
end
if profiling_enabled
    profile_times.coronal_render = toc(coronal_start);
    fprintf('Coronal slice render time: %.3f seconds\n', profile_times.coronal_render);
    profile_times.total_render = toc(render_start);
    fprintf('Total rendering time: %.3f seconds\n', profile_times.total_render);
end

ax_info = nexttile(t, 4);
drawInfoPanel(ax_info, tracks, slices, opts, nim);

%% Save or display
if should_save
    fprintf('Saving to %s...\n', opts.save);
    print(render_fig, opts.save, '-dpng', '-r300');
    close(render_fig);
    fprintf('Saved tractography slice view to %s\n', opts.save);
else
    fprintf('Displaying slice view...\n');
    drawnow;
end

%% Performance Summary
if profiling_enabled
    profile_times.total_execution = toc(total_start_time);
    fprintf('\n=== PERFORMANCE SUMMARY ===\n');
    fprintf('Track lookup: %.3f seconds (%.1f%%)\n', profile_times.track_lookup, ...
            100 * profile_times.track_lookup / profile_times.total_execution);
    fprintf('Total rendering: %.3f seconds (%.1f%%)\n', profile_times.total_render, ...
            100 * profile_times.total_render / profile_times.total_execution);
    fprintf('  - Axial: %.3f seconds\n', profile_times.axial_render);
    fprintf('  - Sagittal: %.3f seconds\n', profile_times.sagittal_render);
    fprintf('  - Coronal: %.3f seconds\n', profile_times.coronal_render);
    fprintf('Total execution: %.3f seconds\n', profile_times.total_execution);
    fprintf('========================\n');

    % Save profiling data for analysis
    assignin('base', 'slice_profiling_data', profile_times);
end

fprintf('=== Complete ===\n');

end

%% =========================================================================
%% RENDERING HELPERS
%% =========================================================================

function drawSlice(ax, orientation, slice_value, slices, nim, tracks, slice_lookup, opts, dims)
cla(ax);
hold(ax, 'on');
show_anatomy = opts.show_anatomy && isfield(nim, 'FA');

switch orientation
    case 'z'
        if show_anatomy
            fa_slice = nim.FA(:, :, slice_value);
            imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
            colormap(ax, gray);
        end
        tracks_in_slice = getTracksInSlice(tracks, slice_lookup, 'z', slice_value, opts.tolerance);
        for k = 1:length(tracks_in_slice)
            track = tracks_in_slice{k};
            if size(track, 1) > 1
                color = getTrackColor(track, opts.color_mode);
                plot(ax, track(:,1), track(:,2), 'Color', color, 'LineWidth', 1.5);
            end
        end
        xlim(ax, [1, dims(1)]);
        ylim(ax, [1, dims(2)]);
        xlabel(ax, 'X'); ylabel(ax, 'Y');
        title(ax, sprintf('Axial (Z = %d)', slice_value));
        if opts.show_crosshairs
            plot(ax, [slices.sagittal, slices.sagittal], [1, dims(2)], 'r--', 'LineWidth', 0.5);
            plot(ax, [1, dims(1)], [slices.coronal, slices.coronal], 'g--', 'LineWidth', 0.5);
        end

    case 'x'
        if show_anatomy
            fa_slice = squeeze(nim.FA(slice_value, :, :));
            imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
            colormap(ax, gray);
        end
        tracks_in_slice = getTracksInSlice(tracks, slice_lookup, 'x', slice_value, opts.tolerance);
        for k = 1:length(tracks_in_slice)
            track = tracks_in_slice{k};
            if size(track, 1) > 1
                color = getTrackColor(track, opts.color_mode);
                plot(ax, track(:,2), track(:,3), 'Color', color, 'LineWidth', 1.5);
            end
        end
        xlim(ax, [1, dims(2)]);
        ylim(ax, [1, dims(3)]);
        xlabel(ax, 'Y'); ylabel(ax, 'Z');
        title(ax, sprintf('Sagittal (X = %d)', slice_value));
        if opts.show_crosshairs
            plot(ax, [slices.coronal, slices.coronal], [1, dims(3)], 'g--', 'LineWidth', 0.5);
            plot(ax, [1, dims(2)], [slices.axial, slices.axial], 'b--', 'LineWidth', 0.5);
        end

    case 'y'
        if show_anatomy
            fa_slice = squeeze(nim.FA(:, slice_value, :));
            imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
            colormap(ax, gray);
        end
        tracks_in_slice = getTracksInSlice(tracks, slice_lookup, 'y', slice_value, opts.tolerance);
        for k = 1:length(tracks_in_slice)
            track = tracks_in_slice{k};
            if size(track, 1) > 1
                color = getTrackColor(track, opts.color_mode);
                plot(ax, track(:,1), track(:,3), 'Color', color, 'LineWidth', 1.5);
            end
        end
        xlim(ax, [1, dims(1)]);
        ylim(ax, [1, dims(3)]);
        xlabel(ax, 'X'); ylabel(ax, 'Z');
        title(ax, sprintf('Coronal (Y = %d)', slice_value));
        if opts.show_crosshairs
            plot(ax, [slices.sagittal, slices.sagittal], [1, dims(3)], 'r--', 'LineWidth', 0.5);
            plot(ax, [1, dims(1)], [slices.axial, slices.axial], 'b--', 'LineWidth', 0.5);
        end
end

set(ax, 'YDir', 'normal');
axis(ax, 'equal', 'tight');
grid(ax, 'on');
hold(ax, 'off');
end

function drawInfoPanel(ax, tracks, slices, opts, nim)
cla(ax);
axis(ax, 'off');
hold(ax, 'on');

track_lengths = cellfun(@(x) size(x, 1), tracks);

text(ax, 0.05, 0.9, 'TRACTOGRAPHY STATISTICS', 'FontSize', 12, 'FontWeight', 'bold');
info_text = sprintf(['Total tracks: %d\n' ...
                     'Mean length: %.1f points\n' ...
                     'Current slice: X=%d  Y=%d  Z=%d\n' ...
                     'Slice tolerance: %d voxels'], ...
                     length(tracks), mean(track_lengths), ...
                     slices.sagittal, slices.coronal, slices.axial, opts.tolerance);
text(ax, 0.05, 0.6, info_text, 'FontSize', 10, 'VerticalAlignment', 'top');

if isfield(nim, 'labels') && ~isempty(nim.labels)
    text(ax, 0.05, 0.32, sprintf('Atlas regions available: %d', numel(nim.labels)), 'FontSize', 9);
end

if strcmp(opts.color_mode, 'direction')
    legend_text = 'Color legend: X=Red, Y=Green, Z=Blue';
else
    legend_text = sprintf('Color mode: %s', opts.color_mode);
end
text(ax, 0.05, 0.20, legend_text, 'FontSize', 9);

axis(ax, [0 1 0 1]);
hold(ax, 'off');
end

%% =========================================================================
%% UTILITY FUNCTIONS
%% =========================================================================

function tracks_in_slice = getTracksInSlice(tracks, slice_lookup, dimension, position, tolerance)
tracks_in_slice = {};
if isempty(tracks)
    return;
end

switch dimension
    case 'x'
        lookup_cells = slice_lookup.x;
        axis_dim = slice_lookup.dims(1);
    case 'y'
        lookup_cells = slice_lookup.y;
        axis_dim = slice_lookup.dims(2);
    case 'z'
        lookup_cells = slice_lookup.z;
        axis_dim = slice_lookup.dims(3);
    otherwise
        lookup_cells = {};
        axis_dim = 0;
end

if axis_dim == 0 || isempty(lookup_cells)
    return;
end

position = max(1, min(axis_dim, round(position)));
tolerance = max(0, round(tolerance));
pos_min = max(1, position - tolerance);
pos_max = min(axis_dim, position + tolerance);

idx_cells = lookup_cells(pos_min:pos_max);
if all(cellfun(@isempty, idx_cells))
    return;
end

track_ids = unique([idx_cells{:}]);
if isempty(track_ids)
    return;
end

tracks_in_slice = tracks(track_ids);
end

function color = getTrackColor(track, color_mode)
switch color_mode
    case 'direction'
        if size(track, 1) > 1
            directions = diff(track);
            avg_dir = mean(directions, 1);
            if norm(avg_dir) > 0
                color = abs(avg_dir) / norm(avg_dir);
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

function lookup = buildTrackSliceLookup(tracks, dims)
lookup = struct();
lookup.dims = dims;
lookup.x = cell(dims(1), 1);
lookup.y = cell(dims(2), 1);
lookup.z = cell(dims(3), 1);

for idx = 1:length(tracks)
    track = tracks{idx};
    if isempty(track)
        continue;
    end

    x_idx = unique(max(1, min(dims(1), round(track(:,1)))));
    y_idx = unique(max(1, min(dims(2), round(track(:,2)))));
    z_idx = unique(max(1, min(dims(3), round(track(:,3)))));

    for xi = transpose(x_idx)
        lookup.x{xi}(end+1) = idx; %#ok<AGROW>
    end
    for yi = transpose(y_idx)
        lookup.y{yi}(end+1) = idx; %#ok<AGROW>
    end
    for zi = transpose(z_idx)
        lookup.z{zi}(end+1) = idx; %#ok<AGROW>
    end
end
end

function filtered_tracks = filterTracksByRegion(tracks, parcellation_mask, region_ids, filter_mode, min_overlap)
% Filter tracks based on their relationship to specified brain regions
%
% Arguments:
%   tracks - Cell array of tracks
%   parcellation_mask - 3D parcellation mask
%   region_ids - Array of region IDs to include
%   filter_mode - 'pass_through', 'start_in', 'end_in'
%   min_overlap - Minimum overlap fraction with region (0-1)

if isempty(tracks)
    filtered_tracks = {};
    return;
end

% Ensure region_ids is a row vector
region_ids = region_ids(:)';

filtered_tracks = cell(length(tracks), 1);
track_count = 0;

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 2
        continue; % Skip tracks that are too short
    end

    % Get parcellation labels along the track
    track_labels = getTrackParcellationLabels(track, parcellation_mask);
    if isempty(track_labels)
        continue;
    end

    include_track = false;

    switch filter_mode
        case 'pass_through'
            % Track passes through any of the specified regions
            include_track = any(ismember(track_labels, region_ids));

        case 'start_in'
            % Track starts in one of the specified regions
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                include_track = ismember(valid_labels(1), region_ids);
            end

        case 'end_in'
            % Track ends in one of the specified regions
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                include_track = ismember(valid_labels(end), region_ids);
            end
    end

    % Apply minimum overlap constraint
    if include_track && min_overlap > 0
        region_points = sum(ismember(track_labels, region_ids));
        total_points = length(track_labels);
        overlap_ratio = region_points / total_points;
        include_track = overlap_ratio >= min_overlap;
    end

    if include_track
        track_count = track_count + 1;
        filtered_tracks{track_count} = track;
    end
end

% Trim to actual size
filtered_tracks = filtered_tracks(1:track_count);
end

function track_labels = getTrackParcellationLabels(track, parcellation_mask)
% Get parcellation labels for each point along a track
dims = size(parcellation_mask);
track_labels = zeros(size(track, 1), 1);

for i = 1:size(track, 1)
    pos = round(track(i, :));

    % Check bounds
    if all(pos >= 1) && all(pos <= dims)
        track_labels(i) = parcellation_mask(pos(1), pos(2), pos(3));
    end
end
end
