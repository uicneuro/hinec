function visualizeTractographySlices(tracks_file, nim_file, varargin)
% VISUALIZETRACTOGRAPHYSLICES: Interactive slice controller for tractography
%
% Usage:
%   visualizeTractographySlices('tracks.mat', 'nim.mat')
%   visualizeTractographySlices(..., 'tolerance', 2)
%   visualizeTractographySlices(..., 'show_crosshairs', true)
%
% The controller window lets you adjust axial, sagittal, and coronal slice
% positions. Press "Render View" to open a dedicated figure showing the
% selected slices with the current options.

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

fprintf('=== TRACTOGRAPHY SLICE CONTROLLER ===\n');
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

fprintf('Precomputing track slice lookup tables (one-time)...\n');
slice_lookup = buildTrackSliceLookup(tracks, dims);
fprintf('Lookup ready. Interactive rendering will reuse cached indices.\n');

%% Initialize controller UI
control_fig = figure('Name', 'Tractography Slice Controls', ...
                     'Position', [100, 100, 630, 260], ...
                     'Color', 'w', ...
                     'NumberTitle', 'off', ...
                     'Resize', 'off', ...
                     'CloseRequestFcn', @closeFigure);

% Store data for callbacks
setappdata(control_fig, 'tracks', tracks);
setappdata(control_fig, 'nim', nim);
setappdata(control_fig, 'dims', dims);
setappdata(control_fig, 'opts', opts);
setappdata(control_fig, 'slice_lookup', slice_lookup);

slices = struct();
slices.axial = round(dims(3) / 2);
slices.sagittal = round(dims(1) / 2);
slices.coronal = round(dims(2) / 2);
setappdata(control_fig, 'slices', slices);

control_panel = uipanel('Parent', control_fig, ...
                        'Position', [0, 0, 1, 1], ...
                        'BackgroundColor', [0.95, 0.95, 0.95], ...
                        'BorderType', 'line');

%% Slider controls
uicontrol('Parent', control_panel, 'Style', 'text', ...
          'String', 'Axial (Z):', ...
          'Position', [20, 190, 80, 18], ...
          'FontSize', 10, 'HorizontalAlignment', 'left');
slider_axial = uicontrol('Parent', control_panel, 'Style', 'slider', ...
          'Position', [110, 190, 260, 20], ...
          'Min', 1, 'Max', dims(3), 'Value', slices.axial, ...
          'SliderStep', [1/(dims(3)-1), min(1, 10/(dims(3)-1))], ...
          'Callback', @updateAxialSlice);
text_axial = uicontrol('Parent', control_panel, 'Style', 'text', ...
          'Position', [380, 190, 50, 20], 'String', num2str(slices.axial), ...
          'FontSize', 10);

uicontrol('Parent', control_panel, 'Style', 'text', ...
          'String', 'Sagittal (X):', ...
          'Position', [20, 150, 80, 18], ...
          'FontSize', 10, 'HorizontalAlignment', 'left');
slider_sagittal = uicontrol('Parent', control_panel, 'Style', 'slider', ...
          'Position', [110, 150, 260, 20], ...
          'Min', 1, 'Max', dims(1), 'Value', slices.sagittal, ...
          'SliderStep', [1/(dims(1)-1), min(1, 10/(dims(1)-1))], ...
          'Callback', @updateSagittalSlice);
text_sagittal = uicontrol('Parent', control_panel, 'Style', 'text', ...
          'Position', [380, 150, 50, 20], 'String', num2str(slices.sagittal), ...
          'FontSize', 10);

uicontrol('Parent', control_panel, 'Style', 'text', ...
          'String', 'Coronal (Y):', ...
          'Position', [20, 110, 80, 18], ...
          'FontSize', 10, 'HorizontalAlignment', 'left');
slider_coronal = uicontrol('Parent', control_panel, 'Style', 'slider', ...
          'Position', [110, 110, 260, 20], ...
          'Min', 1, 'Max', dims(2), 'Value', slices.coronal, ...
          'SliderStep', [1/(dims(2)-1), min(1, 10/(dims(2)-1))], ...
          'Callback', @updateCoronalSlice);
text_coronal = uicontrol('Parent', control_panel, 'Style', 'text', ...
          'Position', [380, 110, 50, 20], 'String', num2str(slices.coronal), ...
          'FontSize', 10);

setappdata(control_fig, 'slider_axial', slider_axial);
setappdata(control_fig, 'slider_sagittal', slider_sagittal);
setappdata(control_fig, 'slider_coronal', slider_coronal);
setappdata(control_fig, 'text_axial', text_axial);
setappdata(control_fig, 'text_sagittal', text_sagittal);
setappdata(control_fig, 'text_coronal', text_coronal);

%% Additional controls
uicontrol('Parent', control_panel, 'Style', 'text', ...
          'String', 'Slice thickness (voxels):', ...
          'Position', [20, 70, 150, 20], ...
          'FontSize', 10, 'HorizontalAlignment', 'left');

tolerance_spinner = uicontrol('Parent', control_panel, 'Style', 'edit', ...
          'Position', [170, 70, 60, 22], ...
          'String', num2str(opts.tolerance), ...
          'Callback', @updateTolerance);

checkbox_crosshairs = uicontrol('Parent', control_panel, 'Style', 'checkbox', ...
          'Position', [250, 70, 150, 22], ...
          'String', 'Show crosshairs', ...
          'Value', opts.show_crosshairs, ...
          'Callback', @toggleCrosshairs);

checkbox_anatomy = uicontrol('Parent', control_panel, 'Style', 'checkbox', ...
          'Position', [250, 40, 150, 22], ...
          'String', 'Show anatomy', ...
          'Value', opts.show_anatomy, ...
          'Callback', @toggleAnatomy);

uicontrol('Parent', control_panel, 'Style', 'pushbutton', ...
          'Position', [430, 190, 170, 32], ...
          'String', 'Render View', ...
          'FontSize', 11, 'FontWeight', 'bold', ...
          'Callback', @renderCurrentView);

uicontrol('Parent', control_panel, 'Style', 'pushbutton', ...
          'Position', [430, 145, 170, 30], ...
          'String', 'Reset to Center', ...
          'FontSize', 10, ...
          'Callback', @resetView);

uicontrol('Parent', control_panel, 'Style', 'pushbutton', ...
          'Position', [430, 105, 170, 30], ...
          'String', 'Render && Export PNG', ...
          'FontSize', 10, ...
          'Callback', @exportView);

uicontrol('Parent', control_panel, 'Style', 'text', ...
          'Position', [20, 10, 580, 24], ...
          'String', 'Adjust slices, then press "Render View" to open a figure. Use "Render && Export" to save a snapshot.', ...
          'FontSize', 9, 'HorizontalAlignment', 'left', ...
          'BackgroundColor', [0.95, 0.95, 0.95]);

fprintf('Controller ready. Adjust sliders and press "Render View" when ready.\n');

end

%% =========================================================================
%% CALLBACKS
%% =========================================================================

function updateAxialSlice(hObject, ~)
fig = ancestor(hObject, 'figure');
dims = getappdata(fig, 'dims');
slices = getappdata(fig, 'slices');
value = round(get(hObject, 'Value'));
value = min(max(value, 1), dims(3));
slices.axial = value;
setappdata(fig, 'slices', slices);
set(getappdata(fig, 'text_axial'), 'String', num2str(slices.axial));
end

function updateSagittalSlice(hObject, ~)
fig = ancestor(hObject, 'figure');
dims = getappdata(fig, 'dims');
slices = getappdata(fig, 'slices');
value = round(get(hObject, 'Value'));
value = min(max(value, 1), dims(1));
slices.sagittal = value;
setappdata(fig, 'slices', slices);
set(getappdata(fig, 'text_sagittal'), 'String', num2str(slices.sagittal));
end

function updateCoronalSlice(hObject, ~)
fig = ancestor(hObject, 'figure');
dims = getappdata(fig, 'dims');
slices = getappdata(fig, 'slices');
value = round(get(hObject, 'Value'));
value = min(max(value, 1), dims(2));
slices.coronal = value;
setappdata(fig, 'slices', slices);
set(getappdata(fig, 'text_coronal'), 'String', num2str(slices.coronal));
end

function updateTolerance(hObject, ~)
fig = ancestor(hObject, 'figure');
opts = getappdata(fig, 'opts');
new_tol = str2double(get(hObject, 'String'));
if isnan(new_tol) || new_tol <= 0
    set(hObject, 'String', num2str(opts.tolerance));
    return;
end
opts.tolerance = max(0, round(new_tol));
setappdata(fig, 'opts', opts);
set(hObject, 'String', num2str(opts.tolerance));
end

function toggleCrosshairs(hObject, ~)
fig = ancestor(hObject, 'figure');
opts = getappdata(fig, 'opts');
opts.show_crosshairs = logical(get(hObject, 'Value'));
setappdata(fig, 'opts', opts);
end

function toggleAnatomy(hObject, ~)
fig = ancestor(hObject, 'figure');
opts = getappdata(fig, 'opts');
opts.show_anatomy = logical(get(hObject, 'Value'));
setappdata(fig, 'opts', opts);
end

function resetView(hObject, ~)
fig = ancestor(hObject, 'figure');
dims = getappdata(fig, 'dims');
slices = struct();
slices.axial = round(dims(3) / 2);
slices.sagittal = round(dims(1) / 2);
slices.coronal = round(dims(2) / 2);
setappdata(fig, 'slices', slices);

set(getappdata(fig, 'slider_axial'), 'Value', slices.axial);
set(getappdata(fig, 'slider_sagittal'), 'Value', slices.sagittal);
set(getappdata(fig, 'slider_coronal'), 'Value', slices.coronal);

set(getappdata(fig, 'text_axial'), 'String', num2str(slices.axial));
set(getappdata(fig, 'text_sagittal'), 'String', num2str(slices.sagittal));
set(getappdata(fig, 'text_coronal'), 'String', num2str(slices.coronal));
end

function renderCurrentView(hObject, ~)
fig = ancestor(hObject, 'figure');
renderSlices(fig, 'display');
end

function exportView(hObject, ~)
fig = ancestor(hObject, 'figure');
render_fig = renderSlices(fig, 'export');
if isempty(render_fig)
    return;
end
filename = sprintf('tractography_slices_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
print(render_fig, filename, '-dpng', '-r300');
close(render_fig);
fprintf('Exported render to %s\n', filename);
end

function closeFigure(fig, ~)
if ishghandle(fig)
    delete(fig);
end
end

%% =========================================================================
%% RENDERING HELPERS
%% =========================================================================

function render_fig = renderSlices(control_fig, mode)
if nargin < 2
    mode = 'display';
end

tracks = getappdata(control_fig, 'tracks');
nim = getappdata(control_fig, 'nim');
dims = getappdata(control_fig, 'dims');
slices = getappdata(control_fig, 'slices');
opts = getappdata(control_fig, 'opts');
slice_lookup = getappdata(control_fig, 'slice_lookup');

if isempty(tracks)
    warning('No tracks available to render.');
    render_fig = [];
    return;
end

fig_visible = 'on';
if strcmp(mode, 'export')
    fig_visible = 'off';
end

render_fig = figure('Name', sprintf('Slice Render [X=%d Y=%d Z=%d]', ...
                                    slices.sagittal, slices.coronal, slices.axial), ...
                    'Color', 'w', 'NumberTitle', 'off', ...
                    'Visible', fig_visible, ...
                    'Position', [760, 120, 820, 740]);

t = tiledlayout(render_fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_axial = nexttile(t, 1);
drawSlice(ax_axial, 'z', slices.axial, slices, nim, tracks, slice_lookup, opts, dims);

ax_sagittal = nexttile(t, 2);
drawSlice(ax_sagittal, 'x', slices.sagittal, slices, nim, tracks, slice_lookup, opts, dims);

ax_coronal = nexttile(t, 3);
drawSlice(ax_coronal, 'y', slices.coronal, slices, nim, tracks, slice_lookup, opts, dims);

ax_info = nexttile(t, 4);
drawInfoPanel(ax_info, tracks, slices, opts, nim);

drawnow;
end

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
