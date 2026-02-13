function visualize_tractography_example(output_dir)
% visualize_tractography_example: Synthetic tractography workflow demonstration
%
% Creates 2x2 panel demonstration using SYNTHETIC data showing:
% 1. Seed placement in white matter (high FA region)
% 2. Tract propagation following principal directions
% 3. FA threshold termination at tissue boundary
% 4. Final complete fiber tract
%
% Features:
% - Realistic synthetic tensor field with coherent flow
% - Proper tissue segmentation (WM, GM, CSF)
% - Large grid for detailed visualization
%
% Arguments:
%   output_dir - Directory to save figures (default: 'presentation_figures')

if nargin < 1
    output_dir = 'presentation_figures';
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('=== Generating Synthetic Tractography Demonstration ===\n');

%% Parameters
grid_size = 32;  % Large grid for detailed visualization
step_size = 0.3;
fa_threshold = 0.15;
angle_threshold = 60;  % degrees

%% Generate synthetic tissue masks and tensor field
fprintf('Generating synthetic brain tissue...\n');
[fa_map, tissue_mask, directions] = generate_synthetic_brain(grid_size);

%% Find optimal seed location and generate tract
fprintf('Computing demonstration tract...\n');
[seed, tract] = compute_demo_tract(fa_map, tissue_mask, directions, ...
    step_size, fa_threshold, angle_threshold, grid_size);

fprintf('Tract length: %d points\n', size(tract, 1));

%% Create 2x2 figure
fig = figure('Position', [100, 100, 1200, 1200], 'Color', 'w');

% Define consistent styling
arrow_scale = 0.35;
seed_size = 14;
tract_width = 2.5;

%% Panel 1: Seed Placement (top-left)
subplot(2, 2, 1);
draw_base_panel(grid_size, fa_map, tissue_mask, directions, arrow_scale);
title('{\bf1. Seed Placement}', 'FontSize', 14);

% Draw seed point
plot(seed(1), seed(2), 'o', 'MarkerSize', seed_size, 'LineWidth', 3, ...
    'MarkerEdgeColor', [0.8 0 0], 'MarkerFaceColor', [1 0.3 0.3]);
text(seed(1), seed(2) - 1.2, 'SEED', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0 0]);

% Add annotation
text(2, grid_size - 1, 'Place seed in WM (high FA)', ...
    'FontSize', 9, 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');

%% Panel 2: Tract Propagation (top-right)
subplot(2, 2, 2);
draw_base_panel(grid_size, fa_map, tissue_mask, directions, arrow_scale);
title('{\bf2. Tract Propagation}', 'FontSize', 14);

% Show partial tract (first 50%)
n_show = max(3, round(0.5 * size(tract, 1)));
plot(tract(1:n_show, 1), tract(1:n_show, 2), '-', ...
    'Color', [0 0.4 0.8], 'LineWidth', tract_width);
plot(tract(1:n_show, 1), tract(1:n_show, 2), 'o', ...
    'MarkerSize', 6, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', [0 0.3 0.7], 'MarkerFaceColor', [0.3 0.6 1]);

% Seed marker
plot(seed(1), seed(2), 'o', 'MarkerSize', seed_size, 'LineWidth', 3, ...
    'MarkerEdgeColor', [0 0.6 0], 'MarkerFaceColor', [0.3 0.8 0.3]);

% Show current direction arrow at propagation front
front_idx = n_show;
front_pos = tract(front_idx, :);
if front_idx > 1
    dir_vec = tract(front_idx, :) - tract(front_idx-1, :);
    dir_vec = dir_vec / norm(dir_vec) * 1.5;
    quiver(front_pos(1), front_pos(2), dir_vec(1), dir_vec(2), 0, ...
        'Color', [0.8 0.4 0], 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
end

text(2, grid_size - 1, 'Follow principal eigenvector', ...
    'FontSize', 9, 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');

%% Panel 3: FA Threshold Check (bottom-left)
subplot(2, 2, 3);
draw_base_panel(grid_size, fa_map, tissue_mask, directions, arrow_scale);
title('{\bf3. Termination Check}', 'FontSize', 14);

% Full tract
plot(tract(:, 1), tract(:, 2), '-', 'Color', [0 0.4 0.8], 'LineWidth', tract_width);

% Seed (green - start)
plot(seed(1), seed(2), 'o', 'MarkerSize', seed_size, 'LineWidth', 3, ...
    'MarkerEdgeColor', [0 0.6 0], 'MarkerFaceColor', [0.3 0.8 0.3]);

% Termination point (red X)
end_pt = tract(end, :);
plot(end_pt(1), end_pt(2), 'x', 'MarkerSize', 16, 'LineWidth', 4, ...
    'Color', [0.8 0 0]);
text(end_pt(1) + 1.5, end_pt(2), 'STOP', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.8 0 0], ...
    'BackgroundColor', [1 1 1 0.9]);

% Show FA value at termination
end_row = round(end_pt(2));
end_col = round(end_pt(1));
end_row = max(1, min(grid_size, end_row));
end_col = max(1, min(grid_size, end_col));
fa_at_end = fa_map(end_row, end_col);

text(2, grid_size - 1, sprintf('Stop when FA < %.2f', fa_threshold), ...
    'FontSize', 9, 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');
text(2, grid_size - 3, sprintf('FA at end: %.2f', fa_at_end), ...
    'FontSize', 9, 'BackgroundColor', [1 0.9 0.9 0.8], 'EdgeColor', [0.8 0 0]);

%% Panel 4: Final Result (bottom-right)
subplot(2, 2, 4);
draw_base_panel(grid_size, fa_map, tissue_mask, directions, arrow_scale);
title('{\bf4. Complete Fiber Tract}', 'FontSize', 14);

% Draw tract with gradient coloring
n_pts = size(tract, 1);
colors = [linspace(0.3, 0, n_pts)', linspace(0.7, 0.3, n_pts)', linspace(0.3, 0.8, n_pts)'];
for i = 1:n_pts-1
    plot(tract(i:i+1, 1), tract(i:i+1, 2), '-', ...
        'Color', colors(i, :), 'LineWidth', tract_width + 1);
end

% Start marker (green)
plot(seed(1), seed(2), 'o', 'MarkerSize', seed_size, 'LineWidth', 3, ...
    'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', [0.2 0.8 0.2]);
text(seed(1), seed(2) - 1.2, 'START', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', [0 0.5 0]);

% End marker (red)
plot(end_pt(1), end_pt(2), 'o', 'MarkerSize', seed_size, 'LineWidth', 3, ...
    'MarkerEdgeColor', [0.6 0 0], 'MarkerFaceColor', [1 0.3 0.3]);
text(end_pt(1), end_pt(2) + 1.5, 'END', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.6 0 0]);

% Statistics
text(2, grid_size - 1, sprintf('Tract length: %d steps', n_pts), ...
    'FontSize', 9, 'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');

%% Add legend for tissue types
annotation('textbox', [0.35, 0.02, 0.3, 0.04], ...
    'String', '\color[rgb]{0.9,0.9,0.95}■\color{black} CSF   \color[rgb]{0.6,0.6,0.65}■\color{black} GM   \color[rgb]{1,1,1}■\color{black} WM', ...
    'FontSize', 11, 'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none', 'BackgroundColor', [0.95 0.95 0.95]);

%% Save figure
sgtitle('Tractography Algorithm: Step-by-Step', 'FontSize', 16, 'FontWeight', 'bold');

saveas(fig, fullfile(output_dir, 'tractography_workflow_2x2.png'));
saveas(fig, fullfile(output_dir, 'tractography_workflow_2x2.fig'));
fprintf('Saved: %s\n', fullfile(output_dir, 'tractography_workflow_2x2.png'));

% Also save high-res version
print(fig, fullfile(output_dir, 'tractography_workflow_2x2_hires.png'), '-dpng', '-r300');
fprintf('Saved high-res: %s\n', fullfile(output_dir, 'tractography_workflow_2x2_hires.png'));

close(fig);
fprintf('Tractography demonstration complete\n\n');

end


%% ========== HELPER FUNCTIONS ==========

function [fa_map, tissue_mask, directions] = generate_synthetic_brain(grid_size)
% Generate realistic synthetic brain tissue with coherent tensor field
%
% Tissue mask values:
%   0 = Background/outside brain
%   1 = CSF (cerebrospinal fluid) - very low FA
%   2 = Gray matter - low FA
%   3 = White matter - high FA with coherent directions

% Initialize
fa_map = zeros(grid_size, grid_size);
tissue_mask = zeros(grid_size, grid_size);
directions = zeros(grid_size, grid_size, 2);

center = grid_size / 2;

% Create layered structure (like a coronal slice)
for row = 1:grid_size
    for col = 1:grid_size
        % Distance from center
        dist_from_center = sqrt((row - center)^2 + (col - center)^2);

        % Elliptical brain boundary
        brain_radius_x = grid_size * 0.42;
        brain_radius_y = grid_size * 0.38;
        normalized_dist = sqrt(((row - center)/brain_radius_y)^2 + ((col - center)/brain_radius_x)^2);

        if normalized_dist > 1.0
            % Outside brain
            tissue_mask(row, col) = 0;
            fa_map(row, col) = 0;
            directions(row, col, :) = [0, 0];
        elseif normalized_dist > 0.85
            % Cortical gray matter (outer layer)
            tissue_mask(row, col) = 2;
            fa_map(row, col) = 0.08 + 0.06 * rand();
            % Random low-coherence directions in GM
            angle = 2 * pi * rand();
            directions(row, col, :) = [cos(angle), sin(angle)];
        elseif normalized_dist < 0.15
            % Central CSF (ventricles)
            tissue_mask(row, col) = 1;
            fa_map(row, col) = 0.02 + 0.03 * rand();
            directions(row, col, :) = [0, 0];
        else
            % White matter with coherent fiber bundles
            tissue_mask(row, col) = 3;

            % Create flowing direction field (like corpus callosum + corona radiata)
            % Multiple fiber populations for realism

            % Horizontal component (corpus callosum-like)
            horiz_weight = exp(-((row - center)^2) / (grid_size * 3)^2);

            % Vertical component (corona radiata-like)
            vert_weight = exp(-((col - center)^2) / (grid_size * 4)^2);

            % Curved fiber bundle (arcuate-like)
            angle_to_center = atan2(row - center, col - center);
            curve_weight = 0.3 * sin(2 * angle_to_center);

            % Combine directions with smooth transitions
            base_angle = atan2(row - center * 0.7, col - center);

            % Add controlled randomness for realism
            noise = 0.15 * randn();

            % Direction varies smoothly across the field
            if abs(row - center) < grid_size * 0.2
                % More horizontal in middle (corpus callosum region)
                dir_angle = 0.1 * (row - center) / grid_size + noise;
            else
                % More radial toward edges (corona radiata)
                dir_angle = base_angle * 0.6 + noise;
            end

            directions(row, col, :) = [cos(dir_angle), sin(dir_angle)];

            % FA varies: higher in dense WM bundles, lower near boundaries
            dist_to_gm = (0.85 - normalized_dist) / 0.7;  % 0 at GM boundary, 1 at center
            dist_to_csf = (normalized_dist - 0.15) / 0.7;  % 0 at CSF, 1 at edge

            boundary_factor = min(dist_to_gm, dist_to_csf);
            boundary_factor = max(0, min(1, boundary_factor * 2));  % Sharper transition

            base_fa = 0.45 + 0.25 * boundary_factor;
            fa_map(row, col) = base_fa + 0.08 * randn();
            fa_map(row, col) = max(0.2, min(0.85, fa_map(row, col)));
        end
    end
end

% Smooth the FA map slightly for realism
fa_map = imgaussfilt(fa_map, 0.8);

% Ensure directions are normalized
for row = 1:grid_size
    for col = 1:grid_size
        dir_vec = squeeze(directions(row, col, :));
        if norm(dir_vec) > 0.01
            directions(row, col, :) = dir_vec / norm(dir_vec);
        end
    end
end

end


function [seed, tract] = compute_demo_tract(fa_map, tissue_mask, directions, ...
    step_size, fa_threshold, angle_threshold, grid_size)
% Find optimal seed and compute demonstration tract

cos_thresh = cosd(angle_threshold);

% Find good seed location in white matter with high FA
best_seed = [];
best_tract = [];
best_length = 0;

% Try multiple seeds to find a good demonstration tract
for attempt = 1:100
    % Random seed in white matter region
    wm_indices = find(tissue_mask == 3 & fa_map > 0.4);
    if isempty(wm_indices)
        wm_indices = find(tissue_mask == 3 & fa_map > 0.25);
    end

    idx = wm_indices(randi(length(wm_indices)));
    [seed_row, seed_col] = ind2sub([grid_size, grid_size], idx);
    test_seed = [seed_col, seed_row];  % [x, y] format

    % Track from this seed
    test_tract = track_fiber(test_seed, fa_map, directions, ...
        step_size, fa_threshold, cos_thresh, grid_size);

    % Keep if longer
    if size(test_tract, 1) > best_length
        best_length = size(test_tract, 1);
        best_seed = test_seed;
        best_tract = test_tract;
    end

    % Good enough?
    if best_length >= 25
        break;
    end
end

seed = best_seed;
tract = best_tract;

if isempty(tract) || size(tract, 1) < 5
    % Fallback: create a simple curved tract for demonstration
    warning('Could not find good tract, using fallback');
    seed = [grid_size/2, grid_size/3];
    t = linspace(0, 1, 30)';
    tract = [seed(1) + 8*t - 4*t.^2, seed(2) + 15*t];
end

end


function tract = track_fiber(seed, fa_map, directions, step_size, fa_threshold, cos_thresh, grid_size)
% Simple Euler fiber tracking

max_steps = 200;
tract = zeros(max_steps, 2);
tract(1, :) = seed;
n_pts = 1;

pos = seed;
prev_dir = [0, 0];

for step = 1:max_steps
    % Get current voxel
    col = round(pos(1));
    row = round(pos(2));

    % Boundary check
    if col < 1 || col > grid_size || row < 1 || row > grid_size
        break;
    end

    % FA check
    fa_val = fa_map(row, col);
    if fa_val < fa_threshold
        break;
    end

    % Get interpolated direction
    dir_vec = interpolate_direction(pos, directions, grid_size);
    if norm(dir_vec) < 0.01
        break;
    end
    dir_vec = dir_vec / norm(dir_vec);

    % Ensure consistent direction (don't flip)
    if n_pts > 1 && dot(dir_vec, prev_dir) < 0
        dir_vec = -dir_vec;
    end

    % Angle check
    if n_pts > 1 && dot(dir_vec, prev_dir) < cos_thresh
        break;
    end

    % Take step
    new_pos = pos + step_size * dir_vec;

    % Store
    n_pts = n_pts + 1;
    tract(n_pts, :) = new_pos;

    pos = new_pos;
    prev_dir = dir_vec;
end

tract = tract(1:n_pts, :);

end


function dir_vec = interpolate_direction(pos, directions, grid_size)
% Bilinear interpolation of direction field

x = pos(1);
y = pos(2);

% Get corner indices
x0 = floor(x); x1 = x0 + 1;
y0 = floor(y); y1 = y0 + 1;

% Clamp to valid range
x0 = max(1, min(grid_size, x0));
x1 = max(1, min(grid_size, x1));
y0 = max(1, min(grid_size, y0));
y1 = max(1, min(grid_size, y1));

% Interpolation weights
wx = x - floor(x);
wy = y - floor(y);

% Get corner directions
d00 = squeeze(directions(y0, x0, :))';
d01 = squeeze(directions(y1, x0, :))';
d10 = squeeze(directions(y0, x1, :))';
d11 = squeeze(directions(y1, x1, :))';

% Handle direction ambiguity (eigenvectors can flip sign)
if dot(d01, d00) < 0, d01 = -d01; end
if dot(d10, d00) < 0, d10 = -d10; end
if dot(d11, d00) < 0, d11 = -d11; end

% Bilinear interpolation
dir_vec = (1-wx)*(1-wy)*d00 + (1-wx)*wy*d01 + wx*(1-wy)*d10 + wx*wy*d11;

end


function draw_base_panel(grid_size, fa_map, tissue_mask, directions, arrow_scale)
% Draw the base visualization panel with tissue coloring and direction arrows

hold on;
axis equal;
xlim([0.5, grid_size + 0.5]);
ylim([0.5, grid_size + 0.5]);
set(gca, 'YDir', 'reverse', 'XTick', [], 'YTick', []);
box on;

% Create RGB image based on tissue type
rgb_img = zeros(grid_size, grid_size, 3);
for row = 1:grid_size
    for col = 1:grid_size
        switch tissue_mask(row, col)
            case 0  % Background
                rgb_img(row, col, :) = [0.2, 0.2, 0.2];
            case 1  % CSF
                rgb_img(row, col, :) = [0.85, 0.85, 0.92];
            case 2  % Gray matter
                rgb_img(row, col, :) = [0.55, 0.55, 0.6];
            case 3  % White matter - shade by FA
                fa_val = fa_map(row, col);
                brightness = 0.7 + 0.3 * fa_val;
                rgb_img(row, col, :) = [brightness, brightness, brightness];
        end
    end
end

% Display tissue image
imagesc([0.5, grid_size + 0.5], [0.5, grid_size + 0.5], rgb_img);

% Draw subtle grid lines
grid_spacing = 4;  % Draw every 4th line for cleaner look
for i = 0.5:grid_spacing:(grid_size + 0.5)
    plot([0.5, grid_size + 0.5], [i, i], '-', 'Color', [0.4 0.4 0.4 0.3], 'LineWidth', 0.5);
    plot([i, i], [0.5, grid_size + 0.5], '-', 'Color', [0.4 0.4 0.4 0.3], 'LineWidth', 0.5);
end

% Draw direction arrows (subsampled for clarity)
arrow_spacing = 2;
for row = 2:arrow_spacing:grid_size-1
    for col = 2:arrow_spacing:grid_size-1
        % Only show arrows in white matter
        if tissue_mask(row, col) == 3 && fa_map(row, col) > 0.2
            dir_vec = squeeze(directions(row, col, :));
            if norm(dir_vec) > 0.01
                dir_vec = dir_vec / norm(dir_vec) * arrow_scale;

                % Color by FA value
                fa_val = fa_map(row, col);
                arrow_color = [0.2, 0.5 + 0.3*fa_val, 0.2];

                quiver(col, row, dir_vec(1), dir_vec(2), 0, ...
                    'Color', arrow_color, 'LineWidth', 1, ...
                    'MaxHeadSize', 0.8, 'AutoScale', 'off');
            end
        end
    end
end

end
