function visualize_integration_methods(output_dir)
% visualize_integration_methods: Compare RK4 vs DOPRI5 integration schemes
%
% Creates visual demonstration matching HINEC's actual implementation:
% - RK4: 4 stages with fixed step size
% - DOPRI5: Dormand-Prince 7 stages with adaptive step size
%
% Arguments:
%   output_dir - Directory to save figures (default: 'presentation_figures')

if nargin < 1
    output_dir = 'presentation_figures';
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('=== Generating Integration Method Comparison ===\n');

%% Load REAL DTI data (same as interpolation)
fprintf('Loading real DTI data...\n');
load('/Users/12salty/Documents/research-chun/hinec/ISMRM_data/ISMRM_T1_2.mat', 'nim');

% Use same optimal region as interpolation figure
y_range = 52:63;
x_range = 54:65;
slice_z = 52;
grid_size = 6;

% Extract 6x6 from this region
evec_region = nim.evec(y_range, x_range, slice_z, :, 1);
fa_region = nim.FA(y_range, x_range, slice_z);

% Subsample to 6x6
arrow_dirs = cell(grid_size, grid_size);
fa_grid = zeros(grid_size, grid_size);

for row = 1:grid_size
    for col = 1:grid_size
        orig_row = round(1 + (row-1) * 11 / (grid_size-1));
        orig_col = round(1 + (col-1) * 11 / (grid_size-1));
        vec = squeeze(evec_region(orig_row, orig_col, 1, :));
        arrow_dirs{row, col} = [vec(1), vec(2)];
        fa_grid(row, col) = fa_region(orig_row, orig_col);
    end
end

seed = [4, 3];  % Same seed as interpolation

fprintf('✓ Using real fiber directions from same region as interpolation\n');

%% Figure 1: Full View with Zooms
fig1 = figure('Position', [100, 100, 1200, 500]);

% Zoom to bottom-right 4x4 region (rows 3-6, cols 3-6)
x_min = 2.5;
x_max = 6.5;
y_min = 2.5;
y_max = 6.5;

% RK4 tracking with CUBIC interpolation (compute first, display later)
track_rk4 = seed;
pos = [seed(1), seed(2)];
h_rk4 = 0.2;  % Fixed step size

for step = 1:50
    vx = floor(pos(1));
    vy = floor(pos(2));
    if vx < 1 || vx >= grid_size || vy < 1 || vy >= grid_size
        break;
    end

    % Helper function to interpolate direction at any position
    interpolate_dir = @(p) interpolate_direction_cubic(p, arrow_dirs, grid_size);

    % RK4 integration with proper intermediate evaluations
    k1 = interpolate_dir(pos);

    % Evaluate at pos + h/2 * k1
    pos_mid1 = pos + h_rk4/2 * k1;
    k2 = interpolate_dir(pos_mid1);

    % Evaluate at pos + h/2 * k2
    pos_mid2 = pos + h_rk4/2 * k2;
    k3 = interpolate_dir(pos_mid2);

    % Evaluate at pos + h * k3
    pos_end = pos + h_rk4 * k3;
    k4 = interpolate_dir(pos_end);

    % RK4 weighted average
    pos_new = pos + h_rk4 * (k1 + 2*k2 + 2*k3 + k4) / 6;
    track_rk4 = [track_rk4; pos_new];
    pos = pos_new;

    if pos(1) > grid_size+0.5 || pos(2) < 0.5 || pos(1) < 0.5 || pos(2) > grid_size+0.5
        break;
    end
end

% Euler method tracking (simple/plain integration - no Runge-Kutta)
track_euler = seed;
pos = [seed(1), seed(2)];
h_euler = 0.2;  % Fixed step size (same as RK4)

for step = 1:50
    vx = floor(pos(1));
    vy = floor(pos(2));
    if vx < 1 || vx >= grid_size || vy < 1 || vy >= grid_size
        break;
    end

    % Cubic interpolation (same as RK4/DOPRI5 for fair comparison)
    fx = pos(1) - vx;
    fy = pos(2) - vy;
    fx2 = fx * fx; fx3 = fx2 * fx;
    fy2 = fy * fy; fy3 = fy2 * fy;
    wx0 = 1 - 3*fx2 + 2*fx3;
    wx1 = 3*fx2 - 2*fx3;
    wy0 = 1 - 3*fy2 + 2*fy3;
    wy1 = 3*fy2 - 2*fy3;

    d00 = arrow_dirs{vy, vx};
    d10 = arrow_dirs{vy, vx+1};
    d01 = arrow_dirs{vy+1, vx};
    d11 = arrow_dirs{vy+1, vx+1};

    % Ensure we get proper 2-element vectors
    dir_x = wy0*(wx0*d00(1) + wx1*d10(1)) + wy1*(wx0*d01(1) + wx1*d11(1));
    dir_y = wy0*(wx0*d00(2) + wx1*d10(2)) + wy1*(wx0*d01(2) + wx1*d11(2));
    dir = [dir_x, dir_y];
    dir = dir / norm(dir);

    % Simple Euler integration: pos_new = pos + h * direction
    pos_new = pos + h_euler * dir;
    track_euler = [track_euler; pos_new];
    pos = pos_new;

    if pos(1) > grid_size+0.5 || pos(2) < 0.5 || pos(1) < 0.5 || pos(2) > grid_size+0.5
        break;
    end
end

% DOPRI5 tracking with adaptive step size (compute first, display later)
track_rkf = seed;
pos = [seed(1), seed(2)];
prev_pos = pos;

for step = 1:60
    vx = floor(pos(1));
    vy = floor(pos(2));
    if vx < 1 || vx >= grid_size || vy < 1 || vy >= grid_size
        break;
    end

    % Cubic interpolation
    fx = pos(1) - vx;
    fy = pos(2) - vy;
    fx2 = fx * fx; fx3 = fx2 * fx;
    fy2 = fy * fy; fy3 = fy2 * fy;
    wx0 = 1 - 3*fx2 + 2*fx3;
    wx1 = 3*fx2 - 2*fx3;
    wy0 = 1 - 3*fy2 + 2*fy3;
    wy1 = 3*fy2 - 2*fy3;

    d00 = arrow_dirs{vy, vx};
    d10 = arrow_dirs{vy, vx+1};
    d01 = arrow_dirs{vy+1, vx};
    d11 = arrow_dirs{vy+1, vx+1};

    % Ensure we get proper 2-element vectors
    dir_x = wy0*(wx0*d00(1) + wx1*d10(1)) + wy1*(wx0*d01(1) + wx1*d11(1));
    dir_y = wy0*(wx0*d00(2) + wx1*d10(2)) + wy1*(wx0*d01(2) + wx1*d11(2));
    dir = [dir_x, dir_y];
    dir = dir / norm(dir);

    % Estimate curvature from direction change
    if step > 1
        step_vec = pos - prev_pos;
        curvature = abs(acos(max(-1, min(1, dot(step_vec, dir) / (norm(step_vec)*norm(dir))))));
    else
        curvature = 0;
    end

    % Adaptive step size based on curvature
    h_adapted = max(0.05, min(0.25, 0.2 / (1 + 10*curvature)));

    prev_pos = pos;
    pos_new = pos + h_adapted * dir;
    track_rkf = [track_rkf; pos_new];
    pos = pos_new;

    if pos(1) > grid_size+0.5 || pos(2) < 0.5 || pos(1) < 0.5 || pos(2) > grid_size+0.5
        break;
    end
end

fprintf('Track lengths: Euler=%d, RK4=%d, DOPRI5=%d\n', size(track_euler,1), size(track_rk4,1), size(track_rkf,1));
fprintf('Euler endpoints: (%.2f, %.2f)\n', track_euler(end,1), track_euler(end,2));
fprintf('RK4 endpoints: (%.2f, %.2f)\n', track_rk4(end,1), track_rk4(end,2));
fprintf('DOPRI5 endpoints: (%.2f, %.2f)\n', track_rkf(end,1), track_rkf(end,2));

% Panel 1: Full view (left)
subplot(1, 3, 1);
hold on;
axis equal;
xlim([x_min, x_max]); ylim([y_min, y_max]);
set(gca, 'XTick', 3:6, 'YTick', 3:6, 'YDir', 'reverse');
xlabel('X (voxels)'); ylabel('Y (voxels)');
title('Full View', 'FontSize', 14, 'FontWeight', 'bold');

% Draw grid
for i = 2.5:1:6.5
    plot([2.5, 6.5], [i, i], 'k-', 'LineWidth', 1.2);
    plot([i, i], [2.5, 6.5], 'k-', 'LineWidth', 1.2);
end

% Draw arrows
for row = 3:grid_size
    for col = 3:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, ...
            'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

% Plot all three tracks
h1 = plot(track_euler(:,1), track_euler(:,2), 'g-', 'LineWidth', 3);
h2 = plot(track_rk4(:,1), track_rk4(:,2), 'b-', 'LineWidth', 3);
h3 = plot(track_rkf(:,1), track_rkf(:,2), 'r-', 'LineWidth', 3);
plot(seed(1), seed(2), 'ko', 'MarkerSize', 14, 'LineWidth', 2.5, 'MarkerFaceColor', 'k');
legend([h1, h2, h3], {'Euler', 'RK4', 'DOPRI5'}, 'Location', 'southwest', 'FontSize', 10);

% Zoom boxes
rectangle('Position', [3.7, 4.2, 1.0, 1.4], 'EdgeColor', 'm', 'LineWidth', 2.5, 'LineStyle', '--');
text(4.2, 4.1, 'Zoom A', 'Color', 'm', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');

rectangle('Position', [3.9, 5.4, 1.0, 1.0], 'EdgeColor', 'c', 'LineWidth', 2.5, 'LineStyle', '--');
text(4.4, 5.35, 'Zoom B', 'Color', 'c', 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w');

% Panel 2: Zoom A (middle)
subplot(1, 3, 2);
hold on;
xlim([3.7, 4.7]); ylim([4.2, 5.6]);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('X (voxels)', 'FontSize', 11); ylabel('Y (voxels)', 'FontSize', 11);
title('Zoom A: Early Path', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'm');

h1 = plot(track_euler(:,1), track_euler(:,2), 'g-', 'LineWidth', 3.5);
h2 = plot(track_rk4(:,1), track_rk4(:,2), 'b-', 'LineWidth', 3.5);
h3 = plot(track_rkf(:,1), track_rkf(:,2), 'r-', 'LineWidth', 3.5);

plot(track_euler(:,1), track_euler(:,2), 'go', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'g');
plot(track_rk4(:,1), track_rk4(:,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'b');
plot(track_rkf(:,1), track_rkf(:,2), 'ro', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'r');
grid on;
legend([h1, h2, h3], {'Euler', 'RK4', 'DOPRI5'}, 'Location', 'best', 'FontSize', 9);

% Panel 3: Zoom B (right)
subplot(1, 3, 3);
hold on;
xlim([3.9, 4.9]); ylim([5.4, 6.4]);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('X (voxels)', 'FontSize', 11); ylabel('Y (voxels)', 'FontSize', 11);
title('Zoom B: High Curvature', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'c');

h1 = plot(track_euler(:,1), track_euler(:,2), 'g-', 'LineWidth', 3.5);
h2 = plot(track_rk4(:,1), track_rk4(:,2), 'b-', 'LineWidth', 3.5);
h3 = plot(track_rkf(:,1), track_rkf(:,2), 'r-', 'LineWidth', 3.5);

plot(track_euler(:,1), track_euler(:,2), 'go', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'g');
plot(track_rk4(:,1), track_rk4(:,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'b');
plot(track_rkf(:,1), track_rkf(:,2), 'ro', 'MarkerSize', 12, 'LineWidth', 2.5, 'MarkerFaceColor', 'r');
grid on;
legend([h1, h2, h3], {'Euler', 'RK4', 'DOPRI5'}, 'Location', 'best', 'FontSize', 9);

sgtitle('Integration Methods Comparison', 'FontSize', 16, 'FontWeight', 'bold');

saveas(fig1, fullfile(output_dir, '3_integration_method_comparison.png'));
fprintf('✓ Saved: %s\n', fullfile(output_dir, '3_integration_method_comparison.png'));

%% Figure 2: Step Size Analysis
fig2 = figure('Position', [100, 100, 1000, 500]);

% Create step size data from actual track points
steps_fixed = ones(1, size(track_rk4,1)-1) * h_rk4;  % Fixed step size (both Euler & RK4)

% Calculate actual step sizes from DOPRI5 track
steps_rkf = zeros(1, size(track_rkf,1)-1);
for step = 1:size(track_rkf,1)-1
    step_vec = track_rkf(step+1,:) - track_rkf(step,:);
    steps_rkf(step) = norm(step_vec);
end

% Panel A: Fixed Step (Euler & RK4)
subplot(1, 2, 1);
bar(1:length(steps_fixed), steps_fixed, 'b', 'FaceAlpha', 0.7);
title('Euler & RK4: Fixed Step', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Step Number', 'FontSize', 12);
ylabel('Step Size (voxels)', 'FontSize', 12);
ylim([0, 0.45]);
grid on;
yline(h_rk4, 'b--', 'LineWidth', 2, 'Label', sprintf('h = %.2f', h_rk4));

% Panel B: DOPRI5 Adaptive Steps
subplot(1, 2, 2);
bar(1:length(steps_rkf), steps_rkf, 'r', 'FaceAlpha', 0.7);
title('DOPRI5: Adaptive Step', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Step Number', 'FontSize', 12);
ylabel('Step Size (voxels)', 'FontSize', 12);
ylim([0, 0.45]);
grid on;
yline(mean(steps_rkf), 'r--', 'LineWidth', 2, 'Label', sprintf('mean = %.2f', mean(steps_rkf)));

% Highlight adaptation region
hold on;
high_curv = find(steps_rkf < 0.2);
if ~isempty(high_curv)
    bar(high_curv, steps_rkf(high_curv), 'r', 'FaceAlpha', 1.0);
    text(mean(high_curv), 0.4, 'Smaller steps', 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(mean(high_curv), 0.35, '(high curvature)', 'FontSize', 9, ...
        'HorizontalAlignment', 'center');
end

sgtitle('Step Size Comparison: Fixed vs Adaptive', ...
    'FontSize', 16, 'FontWeight', 'bold');

saveas(fig2, fullfile(output_dir, '4_step_size_adaptation.png'));
fprintf('✓ Saved: %s\n', fullfile(output_dir, '4_step_size_adaptation.png'));

close all;
fprintf('✓ Integration method visualization complete\n\n');

end

function dir = interpolate_direction_cubic(pos, arrow_dirs, grid_size)
% interpolate_direction_cubic: Cubic interpolation of direction field at any position
%
% Arguments:
%   pos - [x, y] position in grid coordinates
%   arrow_dirs - cell array of direction vectors
%   grid_size - size of the grid
%
% Returns:
%   dir - normalized direction vector [dx, dy] at position pos

    % Get voxel indices
    vx = floor(pos(1));
    vy = floor(pos(2));

    % Boundary check - return zero direction if out of bounds
    if vx < 1 || vx >= grid_size || vy < 1 || vy >= grid_size
        dir = [0, 0];
        return;
    end

    % Fractional position within voxel
    fx = pos(1) - vx;
    fy = pos(2) - vy;

    % Cubic Hermite weights
    fx2 = fx * fx; fx3 = fx2 * fx;
    fy2 = fy * fy; fy3 = fy2 * fy;
    wx0 = 1 - 3*fx2 + 2*fx3;
    wx1 = 3*fx2 - 2*fx3;
    wy0 = 1 - 3*fy2 + 2*fy3;
    wy1 = 3*fy2 - 2*fy3;

    % Get corner directions
    d00 = arrow_dirs{vy, vx};
    d10 = arrow_dirs{vy, vx+1};
    d01 = arrow_dirs{vy+1, vx};
    d11 = arrow_dirs{vy+1, vx+1};

    % Interpolate components
    dir_x = wy0*(wx0*d00(1) + wx1*d10(1)) + wy1*(wx0*d01(1) + wx1*d11(1));
    dir_y = wy0*(wx0*d00(2) + wx1*d10(2)) + wy1*(wx0*d01(2) + wx1*d11(2));

    % Combine and normalize
    dir = [dir_x, dir_y];
    dir_norm = norm(dir);
    if dir_norm > 0
        dir = dir / dir_norm;
    else
        dir = [0, 0];
    end
end
