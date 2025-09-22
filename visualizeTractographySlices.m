function visualizeTractographySlices(tracks_file, nim_file, varargin)
% VISUALIZETRACTOGRAPHYSLICES: Interactive 2D slice viewer for tractography
%
% Displays tractography in three orthogonal planes (axial, sagittal, coronal)
% with interactive sliders for navigating through brain slices.
%
% Usage:
%   visualizeTractographySlices('tracks.mat', 'data.mat')
%   visualizeTractographySlices(..., 'tolerance', 2)
%   visualizeTractographySlices(..., 'show_crosshairs', true)
%
% Arguments:
%   tracks_file - Path to tracks .mat file
%   nim_file - Path to nim .mat file
%
% Optional Parameters:
%   'tolerance' - Slice thickness in voxels (default: 2)
%   'show_crosshairs' - Show crosshair lines (default: true)
%   'color_mode' - Track coloring ('direction', 'uniform') (default: 'direction')
%   'show_anatomy' - Show FA background (default: true)
%   'alpha' - FA background transparency (default: 0.6)

%% Parse input arguments
p = inputParser;
addRequired(p, 'tracks_file', @(x) ischar(x) || isstring(x));
addRequired(p, 'nim_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'tolerance', 2, @isnumeric);
addParameter(p, 'show_crosshairs', true, @islogical);
addParameter(p, 'color_mode', 'direction', @(x) ismember(x, {'direction', 'uniform'}));
addParameter(p, 'show_anatomy', true, @islogical);
addParameter(p, 'alpha', 0.6, @(x) isnumeric(x) && x >= 0 && x <= 1);

parse(p, tracks_file, nim_file, varargin{:});
opts = p.Results;

%% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

fprintf('=== TRACTOGRAPHY SLICE VIEWER ===\n');
fprintf('Loading data...\n');

%% Load data
% Load tracks
if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s', tracks_file);
end
track_data = load(tracks_file);
if isfield(track_data, 'tracks')
    tracks = track_data.tracks;
else
    error('No tracks field found in tracks file');
end

% Load nim data
if ~exist(nim_file, 'file')
    error('NIM file not found: %s', nim_file);
end
nim_data = load(nim_file);
if isfield(nim_data, 'nim')
    nim = nim_data.nim;
else
    error('No nim field found in nim file');
end

% Validate required fields
if ~isfield(nim, 'FA')
    error('FA field required in nim structure');
end

fprintf('Loaded %d tracks\n', length(tracks));
dims = size(nim.FA);
fprintf('Volume dimensions: %d x %d x %d\n', dims(1), dims(2), dims(3));

%% Initialize figure and UI
fig = figure('Name', 'Tractography Slice Viewer', ...
             'Position', [100, 50, 1400, 900], ...
             'Color', 'w', ...
             'CloseRequestFcn', @closeFigure);

% Store data in figure for callbacks
setappdata(fig, 'tracks', tracks);
setappdata(fig, 'nim', nim);
setappdata(fig, 'dims', dims);
setappdata(fig, 'opts', opts);

% Initialize slice positions (middle of each dimension)
slices = struct();
slices.axial = round(dims(3) / 2);      % Z slice
slices.sagittal = round(dims(1) / 2);   % X slice
slices.coronal = round(dims(2) / 2);    % Y slice
setappdata(fig, 'slices', slices);

% Initialize crosshair lines storage
crosshairs = struct();
setappdata(fig, 'crosshairs', crosshairs);

%% Create UI layout
% Main viewing area (top 75% of figure)
main_panel = uipanel('Parent', fig, ...
                     'Position', [0, 0.25, 1, 0.75], ...
                     'BackgroundColor', 'w', ...
                     'BorderType', 'none');

% Control panel (bottom 25% of figure)
control_panel = uipanel('Parent', fig, ...
                       'Position', [0, 0, 1, 0.25], ...
                       'BackgroundColor', [0.95, 0.95, 0.95], ...
                       'BorderType', 'line');

%% Create subplot axes
% Axial view (top-left)
ax_axial = axes('Parent', main_panel, ...
                'Position', [0.05, 0.52, 0.4, 0.4]);
title(ax_axial, sprintf('Axial (Top View) - Z = %d', slices.axial));
xlabel(ax_axial, 'X'); ylabel(ax_axial, 'Y');
axis(ax_axial, 'equal', 'tight');
hold(ax_axial, 'on');

% Sagittal view (top-right)
ax_sagittal = axes('Parent', main_panel, ...
                   'Position', [0.55, 0.52, 0.4, 0.4]);
title(ax_sagittal, sprintf('Sagittal (Side View) - X = %d', slices.sagittal));
xlabel(ax_sagittal, 'Y'); ylabel(ax_sagittal, 'Z');
axis(ax_sagittal, 'equal', 'tight');
hold(ax_sagittal, 'on');

% Coronal view (bottom-left)
ax_coronal = axes('Parent', main_panel, ...
                  'Position', [0.05, 0.05, 0.4, 0.4]);
title(ax_coronal, sprintf('Coronal (Front View) - Y = %d', slices.coronal));
xlabel(ax_coronal, 'X'); ylabel(ax_coronal, 'Z');
axis(ax_coronal, 'equal', 'tight');
hold(ax_coronal, 'on');

% Info panel (bottom-right)
ax_info = axes('Parent', main_panel, ...
               'Position', [0.55, 0.05, 0.4, 0.4]);
axis(ax_info, 'off');
hold(ax_info, 'on');

% Store axes handles
setappdata(fig, 'ax_axial', ax_axial);
setappdata(fig, 'ax_sagittal', ax_sagittal);
setappdata(fig, 'ax_coronal', ax_coronal);
setappdata(fig, 'ax_info', ax_info);

%% Create slider controls
% Axial slider
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', 'Axial (Z):', ...
          'Position', [20, 140, 80, 20], ...
          'FontSize', 10, ...
          'HorizontalAlignment', 'left');

slider_axial = uicontrol('Parent', control_panel, ...
                        'Style', 'slider', ...
                        'Position', [100, 140, 300, 20], ...
                        'Min', 1, 'Max', dims(3), ...
                        'Value', slices.axial, ...
                        'SliderStep', [1/(dims(3)-1), 10/(dims(3)-1)], ...
                        'Callback', @updateAxialSlice);

text_axial = uicontrol('Parent', control_panel, ...
                      'Style', 'text', ...
                      'Position', [410, 140, 50, 20], ...
                      'String', num2str(slices.axial), ...
                      'FontSize', 10);

% Sagittal slider
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', 'Sagittal (X):', ...
          'Position', [20, 100, 80, 20], ...
          'FontSize', 10, ...
          'HorizontalAlignment', 'left');

slider_sagittal = uicontrol('Parent', control_panel, ...
                           'Style', 'slider', ...
                           'Position', [100, 100, 300, 20], ...
                           'Min', 1, 'Max', dims(1), ...
                           'Value', slices.sagittal, ...
                           'SliderStep', [1/(dims(1)-1), 10/(dims(1)-1)], ...
                           'Callback', @updateSagittalSlice);

text_sagittal = uicontrol('Parent', control_panel, ...
                         'Style', 'text', ...
                         'Position', [410, 100, 50, 20], ...
                         'String', num2str(slices.sagittal), ...
                         'FontSize', 10);

% Coronal slider
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', 'Coronal (Y):', ...
          'Position', [20, 60, 80, 20], ...
          'FontSize', 10, ...
          'HorizontalAlignment', 'left');

slider_coronal = uicontrol('Parent', control_panel, ...
                          'Style', 'slider', ...
                          'Position', [100, 60, 300, 20], ...
                          'Min', 1, 'Max', dims(2), ...
                          'Value', slices.coronal, ...
                          'SliderStep', [1/(dims(2)-1), 10/(dims(2)-1)], ...
                          'Callback', @updateCoronalSlice);

text_coronal = uicontrol('Parent', control_panel, ...
                        'Style', 'text', ...
                        'Position', [410, 60, 50, 20], ...
                        'String', num2str(slices.coronal), ...
                        'FontSize', 10);

% Store slider handles
setappdata(fig, 'slider_axial', slider_axial);
setappdata(fig, 'slider_sagittal', slider_sagittal);
setappdata(fig, 'slider_coronal', slider_coronal);
setappdata(fig, 'text_axial', text_axial);
setappdata(fig, 'text_sagittal', text_sagittal);
setappdata(fig, 'text_coronal', text_coronal);

%% Add control buttons
% Tolerance control
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'String', 'Slice Thickness:', ...
          'Position', [500, 140, 100, 20], ...
          'FontSize', 10, ...
          'HorizontalAlignment', 'left');

tolerance_spinner = uicontrol('Parent', control_panel, ...
                             'Style', 'edit', ...
                             'Position', [600, 140, 40, 20], ...
                             'String', num2str(opts.tolerance), ...
                             'Callback', @updateTolerance);

% Crosshairs toggle
checkbox_crosshairs = uicontrol('Parent', control_panel, ...
                               'Style', 'checkbox', ...
                               'Position', [500, 100, 150, 20], ...
                               'String', 'Show Crosshairs', ...
                               'Value', opts.show_crosshairs, ...
                               'Callback', @toggleCrosshairs);

% Anatomy toggle
checkbox_anatomy = uicontrol('Parent', control_panel, ...
                            'Style', 'checkbox', ...
                            'Position', [500, 60, 150, 20], ...
                            'String', 'Show Anatomy', ...
                            'Value', opts.show_anatomy, ...
                            'Callback', @toggleAnatomy);

% Reset button
uicontrol('Parent', control_panel, ...
          'Style', 'pushbutton', ...
          'Position', [700, 100, 100, 30], ...
          'String', 'Reset View', ...
          'FontSize', 10, ...
          'Callback', @resetView);

% Export button
uicontrol('Parent', control_panel, ...
          'Style', 'pushbutton', ...
          'Position', [700, 60, 100, 30], ...
          'String', 'Export View', ...
          'FontSize', 10, ...
          'Callback', @exportView);

% Instructions text
uicontrol('Parent', control_panel, ...
          'Style', 'text', ...
          'Position', [850, 20, 500, 140], ...
          'String', sprintf(['INSTRUCTIONS:\n' ...
                            '• Use sliders to navigate through slices\n' ...
                            '• Crosshairs show current position across views\n' ...
                            '• Adjust slice thickness to see more tracks\n' ...
                            '• Toggle anatomy to show/hide FA background\n' ...
                            '• Export saves current view as image']), ...
          'FontSize', 9, ...
          'HorizontalAlignment', 'left', ...
          'BackgroundColor', [0.95, 0.95, 0.95]);

%% Initial rendering
updateAllViews(fig);
updateInfo(fig);

fprintf('Slice viewer ready!\n');
fprintf('Use sliders to navigate through slices\n');

end


%% ============================================================================
%% CALLBACK FUNCTIONS
%% ============================================================================

function updateAxialSlice(hObject, ~)
% Update axial slice position
fig = ancestor(hObject, 'figure');
slices = getappdata(fig, 'slices');
slices.axial = round(get(hObject, 'Value'));
setappdata(fig, 'slices', slices);

% Update display
text_handle = getappdata(fig, 'text_axial');
set(text_handle, 'String', num2str(slices.axial));

% Redraw
updateAxialView(fig);
updateCrosshairs(fig);
end

function updateSagittalSlice(hObject, ~)
% Update sagittal slice position
fig = ancestor(hObject, 'figure');
slices = getappdata(fig, 'slices');
slices.sagittal = round(get(hObject, 'Value'));
setappdata(fig, 'slices', slices);

% Update display
text_handle = getappdata(fig, 'text_sagittal');
set(text_handle, 'String', num2str(slices.sagittal));

% Redraw
updateSagittalView(fig);
updateCrosshairs(fig);
end

function updateCoronalSlice(hObject, ~)
% Update coronal slice position
fig = ancestor(hObject, 'figure');
slices = getappdata(fig, 'slices');
slices.coronal = round(get(hObject, 'Value'));
setappdata(fig, 'slices', slices);

% Update display
text_handle = getappdata(fig, 'text_coronal');
set(text_handle, 'String', num2str(slices.coronal));

% Redraw
updateCoronalView(fig);
updateCrosshairs(fig);
end

function updateTolerance(hObject, ~)
% Update slice tolerance
fig = ancestor(hObject, 'figure');
opts = getappdata(fig, 'opts');
new_tol = str2double(get(hObject, 'String'));
if ~isnan(new_tol) && new_tol > 0 && new_tol <= 10
    opts.tolerance = new_tol;
    setappdata(fig, 'opts', opts);
    updateAllViews(fig);
else
    set(hObject, 'String', num2str(opts.tolerance));
end
end

function toggleCrosshairs(hObject, ~)
% Toggle crosshairs visibility
fig = ancestor(hObject, 'figure');
opts = getappdata(fig, 'opts');
opts.show_crosshairs = get(hObject, 'Value');
setappdata(fig, 'opts', opts);
updateCrosshairs(fig);
end

function toggleAnatomy(hObject, ~)
% Toggle anatomy visibility
fig = ancestor(hObject, 'figure');
opts = getappdata(fig, 'opts');
opts.show_anatomy = get(hObject, 'Value');
setappdata(fig, 'opts', opts);
updateAllViews(fig);
end

function resetView(~, ~)
% Reset to middle slices
fig = gcf;
dims = getappdata(fig, 'dims');
slices = struct();
slices.axial = round(dims(3) / 2);
slices.sagittal = round(dims(1) / 2);
slices.coronal = round(dims(2) / 2);
setappdata(fig, 'slices', slices);

% Update sliders
set(getappdata(fig, 'slider_axial'), 'Value', slices.axial);
set(getappdata(fig, 'slider_sagittal'), 'Value', slices.sagittal);
set(getappdata(fig, 'slider_coronal'), 'Value', slices.coronal);
set(getappdata(fig, 'text_axial'), 'String', num2str(slices.axial));
set(getappdata(fig, 'text_sagittal'), 'String', num2str(slices.sagittal));
set(getappdata(fig, 'text_coronal'), 'String', num2str(slices.coronal));

updateAllViews(fig);
end

function exportView(~, ~)
% Export current view
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('tractography_slices_%s.png', timestamp);
print(gcf, filename, '-dpng', '-r300');
fprintf('Exported view to: %s\n', filename);
end

function closeFigure(hObject, ~)
% Clean up when closing
delete(hObject);
end


%% ============================================================================
%% VIEW UPDATE FUNCTIONS
%% ============================================================================

function updateAllViews(fig)
% Update all three views
updateAxialView(fig);
updateSagittalView(fig);
updateCoronalView(fig);
updateCrosshairs(fig);
updateInfo(fig);
end

function updateAxialView(fig)
% Update axial (top) view - XY plane at Z slice
ax = getappdata(fig, 'ax_axial');
cla(ax);
hold(ax, 'on');

nim = getappdata(fig, 'nim');
tracks = getappdata(fig, 'tracks');
slices = getappdata(fig, 'slices');
opts = getappdata(fig, 'opts');
dims = getappdata(fig, 'dims');

% Show anatomy slice
if opts.show_anatomy && isfield(nim, 'FA')
    fa_slice = nim.FA(:, :, slices.axial);
    imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
    colormap(ax, gray);
end

% Get and plot tracks in slice
tracks_in_slice = getTracksInSlice(tracks, 'z', slices.axial, opts.tolerance);
for i = 1:length(tracks_in_slice)
    track = tracks_in_slice{i};
    if size(track, 1) > 1
        color = getTrackColor(track, opts.color_mode);
        plot(ax, track(:,1), track(:,2), 'Color', color, 'LineWidth', 1.5);
    end
end

% Set view properties
xlim(ax, [1, dims(1)]);
ylim(ax, [1, dims(2)]);
title(ax, sprintf('Axial (Top View) - Z = %d', slices.axial));
xlabel(ax, 'X'); ylabel(ax, 'Y');
axis(ax, 'equal', 'tight');
grid(ax, 'on');
hold(ax, 'off');
end

function updateSagittalView(fig)
% Update sagittal (side) view - YZ plane at X slice
ax = getappdata(fig, 'ax_sagittal');
cla(ax);
hold(ax, 'on');

nim = getappdata(fig, 'nim');
tracks = getappdata(fig, 'tracks');
slices = getappdata(fig, 'slices');
opts = getappdata(fig, 'opts');
dims = getappdata(fig, 'dims');

% Show anatomy slice
if opts.show_anatomy && isfield(nim, 'FA')
    fa_slice = squeeze(nim.FA(slices.sagittal, :, :));
    imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
    colormap(ax, gray);
end

% Get and plot tracks in slice
tracks_in_slice = getTracksInSlice(tracks, 'x', slices.sagittal, opts.tolerance);
for i = 1:length(tracks_in_slice)
    track = tracks_in_slice{i};
    if size(track, 1) > 1
        color = getTrackColor(track, opts.color_mode);
        plot(ax, track(:,2), track(:,3), 'Color', color, 'LineWidth', 1.5);
    end
end

% Set view properties
xlim(ax, [1, dims(2)]);
ylim(ax, [1, dims(3)]);
title(ax, sprintf('Sagittal (Side View) - X = %d', slices.sagittal));
xlabel(ax, 'Y'); ylabel(ax, 'Z');
axis(ax, 'equal', 'tight');
grid(ax, 'on');
hold(ax, 'off');
end

function updateCoronalView(fig)
% Update coronal (front) view - XZ plane at Y slice
ax = getappdata(fig, 'ax_coronal');
cla(ax);
hold(ax, 'on');

nim = getappdata(fig, 'nim');
tracks = getappdata(fig, 'tracks');
slices = getappdata(fig, 'slices');
opts = getappdata(fig, 'opts');
dims = getappdata(fig, 'dims');

% Show anatomy slice
if opts.show_anatomy && isfield(nim, 'FA')
    fa_slice = squeeze(nim.FA(:, slices.coronal, :));
    imagesc(ax, fa_slice', 'AlphaData', opts.alpha);
    colormap(ax, gray);
end

% Get and plot tracks in slice
tracks_in_slice = getTracksInSlice(tracks, 'y', slices.coronal, opts.tolerance);
for i = 1:length(tracks_in_slice)
    track = tracks_in_slice{i};
    if size(track, 1) > 1
        color = getTrackColor(track, opts.color_mode);
        plot(ax, track(:,1), track(:,3), 'Color', color, 'LineWidth', 1.5);
    end
end

% Set view properties
xlim(ax, [1, dims(1)]);
ylim(ax, [1, dims(3)]);
title(ax, sprintf('Coronal (Front View) - Y = %d', slices.coronal));
xlabel(ax, 'X'); ylabel(ax, 'Z');
axis(ax, 'equal', 'tight');
grid(ax, 'on');
hold(ax, 'off');
end

function updateCrosshairs(fig)
% Update crosshair lines across all views
opts = getappdata(fig, 'opts');
if ~opts.show_crosshairs
    % Remove existing crosshairs
    removeCrosshairs(fig);
    return;
end

slices = getappdata(fig, 'slices');
dims = getappdata(fig, 'dims');

% Remove old crosshairs
removeCrosshairs(fig);

% Axial view crosshairs (show sagittal and coronal positions)
ax = getappdata(fig, 'ax_axial');
hold(ax, 'on');
h1 = plot(ax, [slices.sagittal, slices.sagittal], [1, dims(2)], 'r--', 'LineWidth', 0.5);
h2 = plot(ax, [1, dims(1)], [slices.coronal, slices.coronal], 'g--', 'LineWidth', 0.5);
hold(ax, 'off');

% Sagittal view crosshairs (show coronal and axial positions)
ax = getappdata(fig, 'ax_sagittal');
hold(ax, 'on');
h3 = plot(ax, [slices.coronal, slices.coronal], [1, dims(3)], 'g--', 'LineWidth', 0.5);
h4 = plot(ax, [1, dims(2)], [slices.axial, slices.axial], 'b--', 'LineWidth', 0.5);
hold(ax, 'off');

% Coronal view crosshairs (show sagittal and axial positions)
ax = getappdata(fig, 'ax_coronal');
hold(ax, 'on');
h5 = plot(ax, [slices.sagittal, slices.sagittal], [1, dims(3)], 'r--', 'LineWidth', 0.5);
h6 = plot(ax, [1, dims(1)], [slices.axial, slices.axial], 'b--', 'LineWidth', 0.5);
hold(ax, 'off');

% Store crosshair handles
crosshairs = struct('h1', h1, 'h2', h2, 'h3', h3, 'h4', h4, 'h5', h5, 'h6', h6);
setappdata(fig, 'crosshairs', crosshairs);
end

function removeCrosshairs(fig)
% Remove existing crosshair lines
crosshairs = getappdata(fig, 'crosshairs');
if ~isempty(crosshairs)
    fields = fieldnames(crosshairs);
    for i = 1:length(fields)
        h = crosshairs.(fields{i});
        if ishandle(h)
            delete(h);
        end
    end
end
setappdata(fig, 'crosshairs', struct());
end

function updateInfo(fig)
% Update information panel
ax = getappdata(fig, 'ax_info');
cla(ax);

tracks = getappdata(fig, 'tracks');
slices = getappdata(fig, 'slices');
opts = getappdata(fig, 'opts');

% Calculate statistics
total_tracks = length(tracks);
track_lengths = cellfun(@(x) size(x, 1), tracks);

% Display info
text(ax, 0.1, 0.9, 'TRACTOGRAPHY STATISTICS', ...
     'FontSize', 12, 'FontWeight', 'bold');

info_text = sprintf(['Total tracks: %d\n' ...
                    'Mean length: %.1f points\n' ...
                    'Current position: [%d, %d, %d]\n' ...
                    'Slice tolerance: %d voxels\n'], ...
                    total_tracks, mean(track_lengths), ...
                    slices.sagittal, slices.coronal, slices.axial, ...
                    opts.tolerance);

text(ax, 0.1, 0.7, info_text, 'FontSize', 10, ...
     'VerticalAlignment', 'top');

% Color legend
text(ax, 0.1, 0.3, 'COLOR LEGEND:', 'FontSize', 10, 'FontWeight', 'bold');
if strcmp(opts.color_mode, 'direction')
    text(ax, 0.1, 0.2, 'Red: Left-Right (X)', 'FontSize', 9, 'Color', [1, 0, 0]);
    text(ax, 0.1, 0.1, 'Green: Anterior-Posterior (Y)', 'FontSize', 9, 'Color', [0, 0.7, 0]);
    text(ax, 0.1, 0.0, 'Blue: Superior-Inferior (Z)', 'FontSize', 9, 'Color', [0, 0, 1]);
else
    text(ax, 0.1, 0.2, 'Uniform coloring', 'FontSize', 9);
end

axis(ax, 'off');
xlim(ax, [0, 1]);
ylim(ax, [0, 1]);
end


%% ============================================================================
%% UTILITY FUNCTIONS
%% ============================================================================

function tracks_in_slice = getTracksInSlice(tracks, dimension, position, tolerance)
% Extract tracks that pass through a slice
tracks_in_slice = {};
count = 0;

for i = 1:length(tracks)
    track = tracks{i};
    if isempty(track)
        continue;
    end

    % Get relevant coordinate based on dimension
    switch dimension
        case 'x'
            coords = track(:, 1);
        case 'y'
            coords = track(:, 2);
        case 'z'
            coords = track(:, 3);
        otherwise
            continue;
    end

    % Check if track passes through slice (with tolerance)
    min_pos = position - tolerance;
    max_pos = position + tolerance;
    in_slice = (coords >= min_pos) & (coords <= max_pos);

    if any(in_slice)
        % Extract the portion of track in/near slice
        count = count + 1;
        tracks_in_slice{count} = track;
    end
end
end

function color = getTrackColor(track, color_mode)
% Get color for a track based on coloring mode
switch color_mode
    case 'direction'
        % Color by average direction
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

    case 'uniform'
        % Single color for all tracks
        color = [0.2, 0.6, 0.8];

    otherwise
        color = [0.5, 0.5, 0.5];
end
end