function visualizeTractography(tracks_file, nim_file)
% visualizeTractography: Visualize saved tractography results
%
% Usage:
%   visualizeTractography('tractography_results/tracks_2024-01-01_12-00-00.mat', 'sample_parcellated.mat')
%   visualizeTractography() % Uses most recent tracks file

if nargin < 1
    % Find most recent tracks file
    tracks_dir = 'tractography_results';
    if ~exist(tracks_dir, 'dir')
        error('No tractography results directory found. Run tractography first.');
    end
    
    files = dir(fullfile(tracks_dir, 'tracks_*.mat'));
    if isempty(files)
        error('No tracks files found in %s', tracks_dir);
    end
    
    % Get most recent file
    [~, idx] = max([files.datenum]);
    tracks_file = fullfile(tracks_dir, files(idx).name);
    fprintf('Using most recent tracks file: %s\n', tracks_file);
end

if nargin < 2
    nim_file = 'sample_parcellated.mat';
end

% Load tracks data
fprintf('Loading tracks from %s...\n', tracks_file);
track_data = load(tracks_file);
tracks = track_data.tracks;

% Load nim data for anatomical context
fprintf('Loading anatomical data from %s...\n', nim_file);
nim_data = load(nim_file);
nim = nim_data.nim;

% Create comprehensive visualization
figure('Name', 'Tractography Visualization', 'Position', [100, 100, 1400, 900]);

% Plot 1: Main 3D tracks view
subplot(2,3,[1,2,4,5]);
plot_tracks_3d(tracks, nim);

% Plot 2: Track length histogram
subplot(2,3,3);
plot_track_lengths(track_data);

% Plot 3: Seed distribution
subplot(2,3,6);
plot_seed_distribution(tracks);

% Print summary
print_track_summary(track_data);

fprintf('Visualization complete!\n');
end

function plot_tracks_3d(tracks, nim)
% Main 3D visualization of tracks with anatomy

hold on;

% Show anatomical background (FA slices)
dims = size(nim.FA);
slice_step = round(dims(3) / 5);
slices = slice_step:slice_step:dims(3)-slice_step;

for s = slices
    fa_slice = nim.FA(:, :, s);
    [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
    Z = ones(size(X)) * s;
    surf(X, Y, Z, fa_slice', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end
colormap(gray);

% Plot tracks with direction-based coloring
max_tracks = min(1000, length(tracks));
track_indices = round(linspace(1, length(tracks), max_tracks));

for i = 1:max_tracks
    track = tracks{track_indices(i)};
    if size(track, 1) > 2
        % Direction-based color
        directions = diff(track);
        avg_dir = mean(directions, 1);
        avg_dir_norm = norm(avg_dir);
        if avg_dir_norm > 0
            color = abs(avg_dir) / avg_dir_norm;
            plot3(track(:,1), track(:,2), track(:,3), ...
                  'Color', color, 'LineWidth', 1.2);
        end
    end
end

% Formatting
xlim([1 dims(1)]); ylim([1 dims(2)]); zlim([1 dims(3)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('Tractography Results (%d tracks, showing %d)', length(tracks), max_tracks));
grid on; axis equal; view(45, 30);

% Add color legend
text(0.02, 0.98, 'Colors: Red=L-R, Green=A-P, Blue=S-I', ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'FontSize', 8);
end

function plot_track_lengths(track_data)
% Plot track length distribution

if isfield(track_data, 'track_lengths') && ~isempty(track_data.track_lengths)
    histogram(track_data.track_lengths, 20);
    xlabel('Track Length (mm)');
    ylabel('Count');
    title('Track Length Distribution');
    
    % Add statistics
    mean_len = mean(track_data.track_lengths);
    hold on;
    xline(mean_len, 'r--', sprintf('Mean: %.1f mm', mean_len), 'LineWidth', 2);
    
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
        start_points(i, :) = tracks{i}(1, :);
    end
    
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

function print_track_summary(track_data)
% Print summary statistics

fprintf('\n=== TRACK SUMMARY ===\n');
if isfield(track_data, 'track_stats')
    stats = track_data.track_stats;
    fprintf('Number of tracks: %d\n', stats.num_tracks);
    if stats.num_tracks > 0
        fprintf('Mean length: %.2f mm\n', stats.mean_length);
        fprintf('Length range: %.2f - %.2f mm\n', stats.min_length, stats.max_length);
        fprintf('Total length: %.2f mm\n', stats.total_length);
    end
else
    fprintf('Number of tracks: %d\n', length(track_data.tracks));
end

if isfield(track_data, 'options')
    opts = track_data.options;
    fprintf('\nTracking parameters:\n');
    fprintf('- FA threshold: %.2f\n', opts.fa_threshold);
    fprintf('- Step size: %.2f\n', opts.step_size);
    fprintf('- Angle threshold: %.1f°\n', opts.angle_thresh);
end
fprintf('====================\n');
end 