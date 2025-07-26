function test_improved_tractography()
% Test script to validate improved tractography implementation

fprintf('=== Testing Improved Tractography ===\n');

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

%% Load sample data
data_path = 'sample_parcellated.mat';
if ~exist(data_path, 'file')
    error('Sample data not found: %s. Please run main() first.', data_path);
end

fprintf('Loading data from %s...\n', data_path);
load(data_path, 'nim');

%% Test 1: Standard tractography with improved parameters
fprintf('\n--- Test 1: Standard Tractography ---\n');
options_std = struct();
options_std.seed_density = 1;
options_std.step_size = 0.2;
options_std.fa_threshold = 0.1;
options_std.termination_fa = 0.05;
options_std.angle_thresh = 60;
options_std.max_steps = 3000;
options_std.min_length = 8;
options_std.interp_method = 'linear';

% Create seed mask (reduced for testing)
seed_mask = nim.FA > 0.08;
if isfield(nim, 'parcellation_mask')
    brain_mask = nim.parcellation_mask > 0;
    seed_mask = seed_mask & brain_mask;
end

% Use full seed mask for complete tractography
options_std.seed_mask = seed_mask;

fprintf('Testing standard method with %d seed points...\n', sum(options_std.seed_mask(:)));

tic;
tracks_std = nim_tractography_standard(nim, options_std);
time_std = toc;

%% Test 2: High-order tractography
fprintf('\n--- Test 2: High-Order Tractography ---\n');
options_ho = options_std;
options_ho.order = 3;
options_ho.interp_method = 'spline';

tic;
tracks_ho = nim_tractography_highorder(nim, options_ho);
time_ho = toc;

%% Analyze results
fprintf('\n=== Results Comparison ===\n');

% Standard method results
if ~isempty(tracks_std)
    lengths_std = cellfun(@(x) calculate_track_length(x, options_std.step_size), tracks_std);
    fprintf('Standard Method:\n');
    fprintf('  Tracks generated: %d\n', length(tracks_std));
    fprintf('  Mean length: %.1f mm (%.1f points)\n', mean(lengths_std), mean(cellfun(@(x) size(x,1), tracks_std)));
    fprintf('  Max length: %.1f mm (%.1f points)\n', max(lengths_std), max(cellfun(@(x) size(x,1), tracks_std)));
    fprintf('  Min length: %.1f mm (%.1f points)\n', min(lengths_std), min(cellfun(@(x) size(x,1), tracks_std)));
    fprintf('  Processing time: %.2f seconds\n', time_std);
else
    fprintf('Standard Method: No tracks generated!\n');
end

% High-order method results
if ~isempty(tracks_ho)
    lengths_ho = cellfun(@(x) calculate_track_length(x, options_ho.step_size), tracks_ho);
    fprintf('\nHigh-Order Method:\n');
    fprintf('  Tracks generated: %d\n', length(tracks_ho));
    fprintf('  Mean length: %.1f mm (%.1f points)\n', mean(lengths_ho), mean(cellfun(@(x) size(x,1), tracks_ho)));
    fprintf('  Max length: %.1f mm (%.1f points)\n', max(lengths_ho), max(cellfun(@(x) size(x,1), tracks_ho)));
    fprintf('  Min length: %.1f mm (%.1f points)\n', min(lengths_ho), min(cellfun(@(x) size(x,1), tracks_ho)));
    fprintf('  Processing time: %.2f seconds\n', time_ho);
else
    fprintf('High-Order Method: No tracks generated!\n');
end

%% Visualization comparison
if ~isempty(tracks_std) || ~isempty(tracks_ho)
    fprintf('\nCreating comparison visualization...\n');
    
    figure('Name', 'Tractography Comparison', 'Position', [100, 100, 1400, 600]);
    
    if ~isempty(tracks_std)
        subplot(1, 2, 1);
        visualize_tracks_sample(tracks_std, nim, 'Standard Method');
        title(sprintf('Standard: %d tracks, %.1f mm avg', ...
                     length(tracks_std), mean(lengths_std)));
    end
    
    if ~isempty(tracks_ho)
        subplot(1, 2, 2);
        visualize_tracks_sample(tracks_ho, nim, 'High-Order Method');
        title(sprintf('High-Order: %d tracks, %.1f mm avg', ...
                     length(tracks_ho), mean(lengths_ho)));
    end
    
    % Display the figure (no saving needed)
    fprintf('Comparison visualization displayed.\n');
end

fprintf('\n=== Test Complete ===\n');
end

function track_length = calculate_track_length(track, step_size)
% Calculate actual track length in mm
track_length = 0;
if size(track, 1) > 1
    for i = 2:size(track, 1)
        track_length = track_length + norm(track(i,:) - track(i-1,:)) * step_size;
    end
end
end

function visualize_tracks_sample(tracks, nim, method_name)
% Visualize tracks in 3D with FA background
hold on;

% Show FA background as semi-transparent slices
dims = size(nim.FA);
slice_step = round(dims(3) / 4);
slices = [slice_step, 2*slice_step, 3*slice_step];

for s = slices
    if s <= dims(3) && s > 0
        fa_slice = nim.FA(:, :, s);
        [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
        Z = ones(size(X)) * s;
        surf(X, Y, Z, fa_slice', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    end
end
colormap(gray);

% Plot sample of tracks in 3D
max_tracks = min(800, length(tracks));
if length(tracks) > max_tracks
    track_indices = randperm(length(tracks), max_tracks);
else
    track_indices = 1:length(tracks);
end

for i = 1:length(track_indices)
    track = tracks{track_indices(i)};
    if size(track, 1) > 3
        % Color by direction (RGB = |xyz direction|)
        directions = diff(track);
        if size(directions, 1) > 0
            avg_dir = mean(directions, 1);
            avg_dir_norm = norm(avg_dir);
            if avg_dir_norm > 0
                color = abs(avg_dir) / avg_dir_norm;
                plot3(track(:,1), track(:,2), track(:,3), 'Color', color, 'LineWidth', 1.0);
            end
        end
    end
end

% Set proper 3D view
xlim([1 dims(1)]);
ylim([1 dims(2)]);
zlim([1 dims(3)]);
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(45, 30);
title(sprintf('%s - 3D Tractography', method_name));

% Add color legend
text(0.02, 0.02, 0.98, 'Color: Red=L-R, Green=A-P, Blue=S-I', ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'FontSize', 8);
end