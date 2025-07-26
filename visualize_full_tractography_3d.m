function visualize_full_tractography_3d(data_path, method)
% Visualize full 3D tractography results
%
% Usage: 
%   visualize_full_tractography_3d('sample_parcellated.mat', 'standard')
%   visualize_full_tractography_3d('sample_parcellated.mat', 'highorder')

if nargin < 2
    method = 'standard';
end

fprintf('=== Full 3D Tractography Visualization ===\n');

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

%% Load data
fprintf('Loading data from %s...\n', data_path);
if ~exist(data_path, 'file')
    error('Data file not found: %s', data_path);
end
load(data_path, 'nim');

%% Set tractography parameters for full run
fprintf('Setting up parameters for %s method...\n', method);
options = struct();
options.seed_density = 1;
options.step_size = 0.2;
options.fa_threshold = 0.1;
options.termination_fa = 0.05;
options.angle_thresh = 60;
options.max_steps = 5000;
options.min_length = 10;

if strcmp(method, 'highorder')
    options.order = 3;
    options.interp_method = 'spline';
else
    options.order = 1;
    options.interp_method = 'linear';
end

% Create seed mask
seed_mask = nim.FA > 0.08;
if isfield(nim, 'parcellation_mask')
    brain_mask = nim.parcellation_mask > 0;
    seed_mask = seed_mask & brain_mask;
    fprintf('Applied brain mask from parcellation\n');
end
options.seed_mask = seed_mask;

fprintf('Seed points: %d\n', sum(options.seed_mask(:)));

%% Run tractography
fprintf('Running %s tractography...\n', method);
tic;
if strcmp(method, 'highorder')
    tracks = nim_tractography_highorder(nim, options);
else
    tracks = nim_tractography_standard(nim, options);
end
elapsed_time = toc;

if isempty(tracks)
    error('No tracks generated!');
end

%% Calculate statistics
track_lengths = cellfun(@(x) calculate_track_length_mm(x, options.step_size), tracks);
fprintf('\n=== Tractography Results ===\n');
fprintf('Total tracks: %d\n', length(tracks));
fprintf('Mean length: %.1f mm (%.1f points)\n', mean(track_lengths), mean(cellfun(@(x) size(x,1), tracks)));
fprintf('Max length: %.1f mm (%.1f points)\n', max(track_lengths), max(cellfun(@(x) size(x,1), tracks)));
fprintf('Min length: %.1f mm (%.1f points)\n', min(track_lengths), min(cellfun(@(x) size(x,1), tracks)));
fprintf('Processing time: %.1f seconds\n', elapsed_time);

%% Create comprehensive 3D visualization
fprintf('Creating 3D visualization...\n');

figure('Name', sprintf('%s Tractography - 3D View', upper(method)), 'Position', [100, 100, 1200, 900]);

% Main 3D plot
subplot(2,2,[1,2]);
plot_3d_tractography_with_anatomy(tracks, nim, options);
title(sprintf('%s Method: %d Tracks (%.1f mm avg length)', upper(method), length(tracks), mean(track_lengths)));

% FA slice view
subplot(2,2,3);
slice_idx = round(size(nim.FA, 3)/2);
imagesc(nim.FA(:,:,slice_idx)');
colormap(gca, gray);
colorbar;
title(sprintf('FA Map - Axial Slice %d', slice_idx));
axis equal; axis tight;
xlabel('X'); ylabel('Y');

% Track length histogram
subplot(2,2,4);
histogram(track_lengths, 30, 'FaceColor', [0.3 0.6 0.9]);
xlabel('Track Length (mm)');
ylabel('Number of Tracks');
title('Track Length Distribution');
grid on;

fprintf('3D visualization complete.\n');
fprintf('=== Done ===\n');
end

function plot_3d_tractography_with_anatomy(tracks, nim, options)
% Plot 3D tractography with anatomical context
hold on;

% Show FA background as semi-transparent slices
dims = size(nim.FA);
slice_step = max(1, round(dims(3) / 6));
slices = slice_step:slice_step:dims(3);

fprintf('Rendering %d FA background slices...\n', length(slices));
for s = slices
    if s <= dims(3)
        fa_slice = nim.FA(:, :, s);
        fa_slice(fa_slice < 0.1) = NaN;  % Hide low FA areas
        [X, Y] = meshgrid(1:dims(2), 1:dims(1));
        Z = ones(size(X)) * s;
        surf(X, Y, Z, fa_slice, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
end

% Plot all tracks with direction-based coloring
fprintf('Rendering %d fiber tracks...\n', length(tracks));
for i = 1:length(tracks)
    if mod(i, 1000) == 0
        fprintf('  Track %d/%d\n', i, length(tracks));
    end
    
    track = tracks{i};
    if size(track, 1) > 3
        % Direction-based RGB coloring (absolute values normalized)
        directions = diff(track);
        if size(directions, 1) > 0
            avg_dir = mean(directions, 1);
            avg_dir_norm = norm(avg_dir);
            if avg_dir_norm > 0
                % RGB color: Red=L-R, Green=A-P, Blue=S-I
                color = abs(avg_dir) / avg_dir_norm;
                plot3(track(:,2), track(:,1), track(:,3), 'Color', color, 'LineWidth', 0.8);
            end
        end
    end
end

% Set proper 3D visualization
xlim([1 dims(2)]);
ylim([1 dims(1)]);
zlim([1 dims(3)]);
axis equal;
grid on;
xlabel('Y (Anterior-Posterior)');
ylabel('X (Left-Right)');
zlabel('Z (Superior-Inferior)');
view(45, 20);

% Add color legend
text(0.02, 0.02, 0.98, 'Fiber Color: Red=L-R, Green=A-P, Blue=S-I', ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'BackgroundColor', 'white', 'FontSize', 10, 'EdgeColor', 'black');

% Add lighting for better 3D effect
lighting gouraud;
light('Position', [1 1 1]);
end

function length_mm = calculate_track_length_mm(track, step_size)
% Calculate track length in millimeters
length_mm = 0;
if size(track, 1) > 1
    for i = 2:size(track, 1)
        length_mm = length_mm + norm(track(i,:) - track(i-1,:)) * step_size;
    end
end
end