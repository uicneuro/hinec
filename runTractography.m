function runTractography(data_path)
% runTractography: Simple entry point for DTI tractography
%
% Usage: runTractography('sample_parcellated.mat')

if nargin < 1
    data_path = 'sample_parcellated.mat';
end

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

fprintf('=== HINEC Tractography Pipeline ===\n');

%% Load data
fprintf('Loading data from %s...\n', data_path);
if ~exist(data_path, 'file')
    error('Data file not found: %s', data_path);
end
load(data_path, 'nim');

%% Check required fields
if ~isfield(nim, 'evec')
    error('Eigenvectors not found. Please run main() first to generate DTI data.');
end
if ~isfield(nim, 'FA')
    error('FA not found. Please run main() first to generate DTI data.');
end

%% Set tractography parameters
fprintf('Setting up tractography parameters...\n');
options = struct();
options.seed_density = 1;
options.step_size = 0.5;
options.fa_threshold = 0.15;  % Lower threshold to get more tracks
options.angle_thresh = 60;    % More permissive angle
options.max_steps = 1000;     % Reasonable step limit
options.order = 1;
options.interp_method = 'linear';

% Create brain-only seed mask with lower FA threshold
seed_mask = nim.FA > 0.1;  % Lower FA threshold for more seeds
if isfield(nim, 'parcellation_mask')
    brain_mask = nim.parcellation_mask > 0;
    seed_mask = seed_mask & brain_mask;
    fprintf('Seed mask restricted to brain tissue (parcellation > 0)\n');
end
options.seed_mask = seed_mask;

%% Run tractography
fprintf('Running standard tractography...\n');
tic;
tracks = nim_tractography_standard(nim, options);
elapsed_time = toc;

fprintf('Tractography completed in %.1f seconds\n', elapsed_time);
fprintf('Generated %d tracks\n', length(tracks));

if isempty(tracks)
    error('No tracks generated! Check FA threshold and seed mask.');
end

%% Compute statistics
track_lengths = cellfun(@(x) size(x, 1), tracks);
fprintf('\nTrack Statistics:\n');
fprintf('  Mean length: %.1f points\n', mean(track_lengths));
fprintf('  Max length: %d points\n', max(track_lengths));
fprintf('  Min length: %d points\n', min(track_lengths));

%% Visualize results
fprintf('Creating visualization...\n');

% Create main tractography figure
figure('Name', 'HINEC Tractography Results', 'Position', [100, 100, 1200, 800]);

% Plot 1: FA background with tracks
subplot(2,2,1);
visualize_tractography_with_anatomy(tracks, nim);
title('Tractography with FA Background');

% Plot 2: Tracks only (direction colored)
subplot(2,2,2);
visualize_tracks_only(tracks, nim);
title('Direction-Colored Tracks');

% Plot 3: FA map
subplot(2,2,3);
slice_idx = round(size(nim.FA, 3)/2);
imagesc(nim.FA(:,:,slice_idx)');
colormap(gca, gray);
colorbar;
title(sprintf('FA Map - Slice %d', slice_idx));
axis equal; axis tight;

% Plot 4: Seed distribution
subplot(2,2,4);
visualize_seed_distribution(options.seed_mask);
title('Seed Point Distribution');

%% Save results
output_dir = 'tractography_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

save(fullfile(output_dir, 'tracks_standard.mat'), 'tracks', 'options', 'elapsed_time');
fprintf('\nResults saved to %s/\n', output_dir);
fprintf('=== Tractography Complete ===\n');
end

function visualize_tractography_with_anatomy(tracks, nim)
% Visualize tracks with anatomical background
hold on;

% Show FA slice as background - use multiple slices for better context
dims = size(nim.FA);
slice_step = round(dims(3) / 4);
slices = [slice_step, 2*slice_step, 3*slice_step];

for s = slices
    if s <= dims(3)
        fa_slice = nim.FA(:, :, s);
        [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
        Z = ones(size(X)) * s;
        surf(X, Y, Z, fa_slice', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
end
colormap(gray);

% Plot tracks throughout the volume
max_tracks = min(800, length(tracks));
for i = 1:max_tracks
    track = tracks{i};
    if size(track, 1) > 3  % Show shorter tracks too
        % Color by direction
        directions = diff(track);
        if size(directions, 1) > 0
            avg_dir = mean(directions, 1);
            avg_dir_norm = norm(avg_dir);
            if avg_dir_norm > 0
                color = abs(avg_dir) / avg_dir_norm;
                plot3(track(:,1), track(:,2), track(:,3), 'Color', color, 'LineWidth', 1.2);
            end
        end
    end
end

% Set proper view and limits
xlim([1 dims(1)]);
ylim([1 dims(2)]);
zlim([1 dims(3)]);
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 30);  % Better viewing angle
end

function visualize_tracks_only(tracks, nim)
% Visualize tracks without background
hold on;

% Plot all tracks with direction coloring
max_tracks = min(1500, length(tracks));
colors_used = [];

for i = 1:max_tracks
    track = tracks{i};
    if size(track, 1) > 3
        % Direction-based coloring
        directions = diff(track);
        if size(directions, 1) > 0
            avg_dir = mean(directions, 1);
            avg_dir_norm = norm(avg_dir);
            if avg_dir_norm > 0
                color = abs(avg_dir) / avg_dir_norm;
                plot3(track(:,1), track(:,2), track(:,3), 'Color', color, 'LineWidth', 1.0);
                colors_used = [colors_used; color];
            end
        end
    end
end

% Set proper axis limits based on brain size
dims = size(nim.FA);
xlim([1 dims(1)]);
ylim([1 dims(2)]);
zlim([1 dims(3)]);

axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 30);  % Better viewing angle

% Add color explanation
if ~isempty(colors_used)
    text(0.02, 0.98, 'Red: L-R, Green: A-P, Blue: S-I', ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'white', 'FontSize', 8);
end
end

function visualize_seed_distribution(seed_mask)
% Show where seeds are distributed
[x, y, z] = ind2sub(size(seed_mask), find(seed_mask));

% Sample points for visualization
max_points = 2000;
if length(x) > max_points
    idx = randperm(length(x), max_points);
    x = x(idx);
    y = y(idx);
    z = z(idx);
end

scatter3(x, y, z, 2, 'filled', 'MarkerFaceAlpha', 0.6);
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
end