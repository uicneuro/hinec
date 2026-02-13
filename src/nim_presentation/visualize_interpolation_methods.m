function visualize_interpolation_methods(output_dir)
% visualize_interpolation_methods: Show interpolation value with dynamic directions
%
% Creates large grid visualizations showing:
% - Figure 1: FACT vs Interpolation comparison
% - Figure 2: How interpolation creates smooth paths
%
% Arguments:
%   output_dir - Directory to save figures (default: 'presentation_figures')

if nargin < 1
    output_dir = 'presentation_figures';
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('=== Generating Interpolation Method Visualization ===\n');

%% Load REAL DTI data from sample
fprintf('Loading real DTI data...\n');
load('/Users/12salty/Documents/research-chun/hinec/ISMRM_data/ISMRM_T1_2.mat', 'nim');

%% Use fixed optimal region
fprintf('Using optimal region found: [52:63, 54:65, z=52], seed [4, 3]\n');

grid_size = 6;  % Display grid size

% Fixed optimal region
y_range = 52:63;
x_range = 54:65;
slice_z = 52;
best_seed_local = [4, 3];

% Extract final 6x6 from best region
evec_region = nim.evec(y_range, x_range, slice_z, :, 1);
fa_region = nim.FA(y_range, x_range, slice_z);

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

seed = best_seed_local;

%% Figure 1: FACT vs Interpolation Comparison - ZOOMED to bottom-right 4x4
fig1 = figure('Position', [100, 100, 1400, 500]);

% Zoom to bottom-right 4x4 region (rows 3-6, cols 3-6)
x_min = 2.5;
x_max = 6.5;
y_min = 2.5;
y_max = 6.5;

% Panel A: FACT (no interpolation - discrete jumps)
subplot(1, 3, 1);
hold on;
axis equal;
xlim([x_min, x_max]);
ylim([y_min, y_max]);
set(gca, 'XTick', 3:6, 'YTick', 3:6, 'YDir', 'reverse');
title('FACT (No Interpolation)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Draw grid for 4x4 region
for i = 2.5:1:6.5
    plot([2.5, 6.5], [i, i], 'k-', 'LineWidth', 1.5);
    plot([i, i], [2.5, 6.5], 'k-', 'LineWidth', 1.5);
end

% Draw direction arrows (only bottom-right 4x4)
for row = 3:grid_size
    for col = 3:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'k', 'LineWidth', 2.5, 'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

% FACT tracking: follow direction until voxel boundary, then switch to new voxel's direction
% Use optimally found seed location
track_fact = seed;
pos = [seed(1), seed(2)];

for step = 1:150
    % Get direction at current voxel center
    voxel = round(pos);
    if voxel(1) < 1 || voxel(1) > grid_size || voxel(2) < 1 || voxel(2) > grid_size
        break;
    end

    dir = arrow_dirs{voxel(2), voxel(1)};

    % Move in this direction until hitting voxel boundary
    step_size = 0.02;  % Very small for detail
    moved = false;
    for substep = 1:300
        pos_new = pos + step_size * dir;

        % Check if we crossed into new voxel
        new_voxel = round(pos_new);
        if ~isequal(new_voxel, voxel)
            % Hit boundary - record this point and switch direction
            track_fact = [track_fact; pos_new];
            pos = pos_new;
            moved = true;
            break;
        end
        pos = pos_new;
    end

    if ~moved || pos(1) > grid_size+0.5 || pos(2) < 0.5
        break;
    end
end

plot(track_fact(:,1), track_fact(:,2), 'r-', 'LineWidth', 4);
plot(track_fact(1,1), track_fact(1,2), 'go', 'MarkerSize', 14, ...
    'LineWidth', 2, 'MarkerFaceColor', 'g');

% Panel B: Trilinear Interpolation
subplot(1, 3, 2);
hold on;
axis equal;
xlim([x_min, x_max]);
ylim([y_min, y_max]);
set(gca, 'XTick', 3:6, 'YTick', 3:6, 'YDir', 'reverse');
title('Trilinear Interpolation', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Draw grid for 4x4 region (lighter)
for i = 2.5:1:6.5
    plot([2.5, 6.5], [i, i], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    plot([i, i], [2.5, 6.5], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
end

% Draw direction arrows (only bottom-right 4x4, lighter)
for row = 3:grid_size
    for col = 3:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, ...
            'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

% Trilinear interpolation: linear interpolation (simpler than cubic)
track_interp = seed;  % Same seed as FACT
pos = [seed(1), seed(2)];
step_size = 0.15;  % Medium step size
for step = 1:80
    % Get integer voxel indices
    vx = floor(pos(1));
    vy = floor(pos(2));

    if vx < 1 || vx >= grid_size || vy < 1 || vy >= grid_size
        break;
    end

    % Fractional position within voxel (0 to 1)
    fx = pos(1) - vx;
    fy = pos(2) - vy;

    % Bilinear interpolation (2D trilinear)
    d00 = arrow_dirs{vy, vx};
    d10 = arrow_dirs{vy, vx+1};
    d01 = arrow_dirs{vy+1, vx};
    d11 = arrow_dirs{vy+1, vx+1};

    % Linear blend
    dir_bottom = (1-fx)*d00 + fx*d10;
    dir_top = (1-fx)*d01 + fx*d11;
    dir = (1-fy)*dir_bottom + fy*dir_top;

    dir = dir / norm(dir);  % Normalize
    pos_new = pos + step_size * dir;
    track_interp = [track_interp; pos_new];
    pos = pos_new;

    if pos(1) > grid_size+0.5 || pos(2) < 0.5 || pos(1) < 0.5
        break;
    end
end

plot(track_interp(:,1), track_interp(:,2), 'b-', 'LineWidth', 4);
plot(track_interp(1,1), track_interp(1,2), 'go', 'MarkerSize', 14, ...
    'LineWidth', 2, 'MarkerFaceColor', 'g');

% Panel C: Cubic Interpolation
subplot(1, 3, 3);
hold on;
axis equal;
xlim([x_min, x_max]);
ylim([y_min, y_max]);
set(gca, 'XTick', 3:6, 'YTick', 3:6, 'YDir', 'reverse');
title('Cubic Interpolation', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Draw grid for 4x4 region (lighter)
for i = 2.5:1:6.5
    plot([2.5, 6.5], [i, i], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    plot([i, i], [2.5, 6.5], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
end

% Draw direction arrows (only bottom-right 4x4, lighter)
for row = 3:grid_size
    for col = 3:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, ...
            'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

% Cubic interpolation: VERY SMOOTH with cubic Hermite spline interpolation
track_cubic = seed;  % Same seed as FACT
pos = [seed(1), seed(2)];
step_size = 0.08;  % Very small steps for ultra-smooth path
for step = 1:150
    % Get integer voxel indices
    vx = floor(pos(1));
    vy = floor(pos(2));

    if vx < 1 || vx >= grid_size || vy < 1 || vy >= grid_size
        break;
    end

    % Fractional position within voxel (0 to 1)
    fx = pos(1) - vx;
    fy = pos(2) - vy;

    % Cubic Hermite interpolation weights
    fx2 = fx * fx;
    fx3 = fx2 * fx;
    fy2 = fy * fy;
    fy3 = fy2 * fy;

    % Cubic basis functions (smoother than linear)
    wx0 = 1 - 3*fx2 + 2*fx3;
    wx1 = 3*fx2 - 2*fx3;
    wy0 = 1 - 3*fy2 + 2*fy3;
    wy1 = 3*fy2 - 2*fy3;

    % Get 4 neighboring directions
    d00 = arrow_dirs{vy, vx};
    d10 = arrow_dirs{vy, vx+1};
    d01 = arrow_dirs{vy+1, vx};
    d11 = arrow_dirs{vy+1, vx+1};

    % Cubic blend (smoother than bilinear)
    dir_bottom = wx0*d00 + wx1*d10;
    dir_top = wx0*d01 + wx1*d11;
    dir = wy0*dir_bottom + wy1*dir_top;

    dir = dir / norm(dir);  % Normalize
    pos_new = pos + step_size * dir;
    track_cubic = [track_cubic; pos_new];
    pos = pos_new;

    if pos(1) > grid_size+0.5 || pos(2) < 0.5 || pos(1) < 0.5
        break;
    end
end

plot(track_cubic(:,1), track_cubic(:,2), 'm-', 'LineWidth', 4);
plot(track_cubic(1,1), track_cubic(1,2), 'go', 'MarkerSize', 14, ...
    'LineWidth', 2, 'MarkerFaceColor', 'g');

sgtitle('Interpolation Methods: FACT vs Trilinear vs Cubic', ...
    'FontSize', 16, 'FontWeight', 'bold');

saveas(fig1, fullfile(output_dir, '1_interpolation_comparison.png'));
fprintf('✓ Saved: %s\n', fullfile(output_dir, '1_interpolation_comparison.png'));

%% Figure 2: Interpolation Detail
fig2 = figure('Position', [100, 100, 1400, 500]);

% Panel A: Discrete directions
subplot(1, 3, 1);
hold on;
axis equal;
xlim([0.5, grid_size+0.5]);
ylim([0.5, grid_size+0.5]);
set(gca, 'XTick', 1:grid_size, 'YTick', 1:grid_size, 'YDir', 'reverse');
title('Discrete Voxel Directions', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Draw grid
for i = 0.5:1:grid_size+0.5
    plot([0.5, grid_size+0.5], [i, i], 'k-', 'LineWidth', 1.5);
    plot([i, i], [0.5, grid_size+0.5], 'k-', 'LineWidth', 1.5);
end

% Draw all direction arrows
for row = 1:grid_size
    for col = 1:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'k', 'LineWidth', 2.5, 'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

plot(seed(1), seed(2), 'ro', 'MarkerSize', 14, 'LineWidth', 3, 'MarkerFaceColor', 'r');
text(seed(1), seed(2)-0.5, 'Seed', 'HorizontalAlignment', 'center', ...
    'FontSize', 11, 'FontWeight', 'bold');

% Panel B: Smooth path
subplot(1, 3, 2);
hold on;
axis equal;
xlim([0.5, grid_size+0.5]);
ylim([0.5, grid_size+0.5]);
set(gca, 'XTick', 1:grid_size, 'YTick', 1:grid_size, 'YDir', 'reverse');
title('Interpolated Track Path', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Draw grid (lighter)
for i = 0.5:1:grid_size+0.5
    plot([0.5, grid_size+0.5], [i, i], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
    plot([i, i], [0.5, grid_size+0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
end

% Draw direction arrows (lighter)
for row = 1:grid_size
    for col = 1:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, ...
            'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

plot(track_interp(:,1), track_interp(:,2), 'r-', 'LineWidth', 5);
plot(track_interp(1,1), track_interp(1,2), 'ro', 'MarkerSize', 14, ...
    'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(track_interp(end,1), track_interp(end,2), 'ro', 'MarkerSize', 12, ...
    'LineWidth', 2, 'MarkerFaceColor', 'r');

% Panel C: Combined
subplot(1, 3, 3);
hold on;
axis equal;
xlim([0.5, grid_size+0.5]);
ylim([0.5, grid_size+0.5]);
set(gca, 'XTick', 1:grid_size, 'YTick', 1:grid_size, 'YDir', 'reverse');
title('Combined View', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Draw grid
for i = 0.5:1:grid_size+0.5
    plot([0.5, grid_size+0.5], [i, i], 'k-', 'LineWidth', 1.5);
    plot([i, i], [0.5, grid_size+0.5], 'k-', 'LineWidth', 1.5);
end

% Draw direction arrows
for row = 1:grid_size
    for col = 1:grid_size
        dir = arrow_dirs{row, col};
        dir_norm = dir / norm(dir) * 0.35;
        quiver(col, row, dir_norm(1), dir_norm(2), 0, ...
            'k', 'LineWidth', 2, 'MaxHeadSize', 0.7, 'AutoScale', 'off');
    end
end

% Draw smooth track
plot(track_interp(:,1), track_interp(:,2), 'r-', 'LineWidth', 5);
plot(track_interp(1,1), track_interp(1,2), 'ro', 'MarkerSize', 14, ...
    'LineWidth', 2, 'MarkerFaceColor', 'r');

% Highlight intermediate points
for i = 2:size(track_interp,1)-1
    plot(track_interp(i,1), track_interp(i,2), 'ro', 'MarkerSize', 6, ...
        'MarkerFaceColor', 'r');
end

sgtitle('Tractography: Discrete Directions → Continuous Path via Interpolation', ...
    'FontSize', 16, 'FontWeight', 'bold');

saveas(fig2, fullfile(output_dir, '2_interpolation_field_comparison.png'));
fprintf('✓ Saved: %s\n', fullfile(output_dir, '2_interpolation_field_comparison.png'));

close all;
fprintf('✓ Interpolation visualization complete\n\n');

end
