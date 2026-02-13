function nim_plot_tractography(tracks, nim, varargin)
% nim_plot_tractography: Visualize tractography results
%
% Arguments:
%   tracks - Cell array of fiber tracks
%   nim - NIM structure
%   options - Visualization options (optional struct)

% Parse input arguments
if nargin > 2 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default values
if ~isfield(options, 'color_mode')
    options.color_mode = "direction";
end
if ~isfield(options, 'show_anatomy')
    options.show_anatomy = true;
end
if ~isfield(options, 'track_subset')
    options.track_subset = [];
end
if ~isfield(options, 'line_width')
    options.line_width = 1;
end
if ~isfield(options, 'alpha')
    options.alpha = 0.8;
end

% Set default track subset
if isempty(options.track_subset)
    options.track_subset = 1:length(tracks);
end

% Limit tracks for performance
max_tracks = 2000;
if length(options.track_subset) > max_tracks
    fprintf('Limiting display to %d tracks for performance\n', max_tracks);
    options.track_subset = options.track_subset(1:max_tracks);
end

figure('Name', 'Tractography Visualization', 'Color', 'w');
hold on;

% Show anatomical background if requested
if options.show_anatomy
    plot_anatomy_slice(nim);
end

% Plot tracks
for i = options.track_subset
    if i > length(tracks)
        continue;
    end
    
    track = tracks{i};
    if size(track, 1) < 2
        continue;
    end
    
    % Get color based on mode
    switch options.color_mode
        case "direction"
            track_color = get_direction_color(track);
        case "fa"
            track_color = get_fa_color(track, nim);
        case "parcellation"
            track_color = get_parcellation_color(track, nim);
        otherwise
            track_color = [0.2, 0.6, 0.8]; % Default blue
    end
    
    % Plot track
    plot3(track(:,1), track(:,2), track(:,3), ...
          'Color', [track_color, options.alpha], ...
          'LineWidth', options.line_width);
end

% Formatting
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Tractography (%s coloring) - %d tracks', ...
              options.color_mode, length(options.track_subset)));

% Set view
view(3);
camlight; lighting gouraud;

hold off;
end

function plot_anatomy_slice(nim)
% Plot anatomical slice as background
slice_idx = round(size(nim.FA, 3) / 2);
fa_slice = nim.FA(:, :, slice_idx);

[X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
Z = ones(size(X)) * slice_idx;

surf(Y, X, Z, fa_slice', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
colormap(gray);
end

function color = get_direction_color(track)
% Color based on local direction
if size(track, 1) < 2
    color = [0.5, 0.5, 0.5];
    return;
end

% Calculate average direction
directions = diff(track);
avg_direction = mean(directions, 1);
avg_direction = avg_direction / norm(avg_direction);

% Map direction to RGB (absolute values)
color = abs(avg_direction);
color = color / max(color); % Normalize
end

function color = get_fa_color(track, nim)
% Color based on FA values along track
if size(track, 1) < 2
    color = [0.5, 0.5, 0.5];
    return;
end

% Sample FA values along track
fa_values = zeros(size(track, 1), 1);
for i = 1:size(track, 1)
    pos = track(i, :);
    if all(pos >= 1) && all(pos <= size(nim.FA))
        fa_values(i) = interp3(nim.FA, pos(2), pos(1), pos(3), 'linear', 0);
    end
end

% Map FA to hot colormap
avg_fa = mean(fa_values);
cmap = hot(256);
fa_idx = round(avg_fa * 255) + 1;
fa_idx = max(1, min(256, fa_idx));
color = cmap(fa_idx, :);
end

function color = get_parcellation_color(track, nim)
% Color based on parcellation regions
if ~isfield(nim, 'parcellation_mask')
    color = [0.5, 0.5, 0.5];
    return;
end

if size(track, 1) < 2
    color = [0.5, 0.5, 0.5];
    return;
end

% Sample parcellation values along track
parcel_values = zeros(size(track, 1), 1);
for i = 1:size(track, 1)
    pos = round(track(i, :));
    if all(pos >= 1) && all(pos <= size(nim.parcellation_mask))
        parcel_values(i) = nim.parcellation_mask(pos(1), pos(2), pos(3));
    end
end

% Get most common parcellation label
parcel_labels = parcel_values(parcel_values > 0);
if isempty(parcel_labels)
    color = [0.5, 0.5, 0.5];
    return;
end

mode_label = mode(parcel_labels);

% Generate color based on label
rng(mode_label); % Seed for consistent colors
color = rand(1, 3);
color = color * 0.8 + 0.2; % Ensure not too dark
end 