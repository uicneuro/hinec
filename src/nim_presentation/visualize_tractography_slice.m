function visualize_tractography_slice(output_dir)
% visualize_tractography_slice: Create ACT mask visualization
%
% Shows tissue classification (ACT mask) from real brain data
% ACT = Anatomically Constrained Tractography
%
% Arguments:
%   output_dir - Directory to save figures (default: 'presentation_figures')

if nargin < 1
    output_dir = 'presentation_figures';
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('=== Generating ACT Mask Visualization ===\n');

%% Load NIM data
nim_file = '/Users/12salty/Documents/research-chun/hinec/ISMRM_data/ISMRM_T1_2.mat';

fprintf('Loading NIM data...\n');
if ~exist(nim_file, 'file')
    error('NIM file not found: %s', nim_file);
end

% Load nim data
nim_data = load(nim_file, 'nim');
nim = nim_data.nim;

fprintf('✓ NIM data loaded (size: %dx%dx%d)\n', size(nim.FA, 1), size(nim.FA, 2), size(nim.FA, 3));

%% Select optimal slice for visualization
fprintf('Selecting representative slice...\n');

vol_size = size(nim.FA);

% Choose a central slice with good brain tissue coverage
best_slice = round(vol_size(3) / 2);  % Middle slice

fprintf('✓ Selected slice z=%d\n', best_slice);

%% Create visualization
fig = figure('Position', [100, 100, 1400, 500]);

% Extract slice data
fa_slice = nim.FA(:,:,best_slice);

% Panel A: FA Map
subplot(1, 3, 1);
imagesc(fa_slice);
colormap(gca, gray);
axis equal tight;
set(gca, 'YDir', 'normal');
title('Fractional Anisotropy (FA)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');
colorbar;
caxis([0 1]);

% Panel B: ACT-like tissue classification
subplot(1, 3, 2);

% Create ACT-like mask: WM (high FA), GM (medium FA), CSF (low FA)
act_mask = zeros(size(fa_slice));
act_mask(fa_slice > 0.5) = 3;   % White matter (high FA)
act_mask(fa_slice > 0.2 & fa_slice <= 0.5) = 2;  % Gray matter (medium FA)
act_mask(fa_slice > 0.05 & fa_slice <= 0.2) = 1;  % CSF (low FA)

% Create color map: CSF=blue, GM=gray, WM=white
cmap = [0 0 0;      % Background (black)
        0.3 0.6 1;    % CSF (light blue)
        0.6 0.6 0.6;  % Gray matter (gray)
        1 1 1];     % White matter (white)

imagesc(act_mask);
colormap(gca, cmap);
axis equal tight;
set(gca, 'YDir', 'normal');
title('ACT Tissue Mask', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');

% Manual legend with patches
hold on;
% Create invisible patches for legend
h1 = patch(NaN, NaN, cmap(2,:), 'EdgeColor', 'none');
h2 = patch(NaN, NaN, cmap(3,:), 'EdgeColor', 'none');
h3 = patch(NaN, NaN, cmap(4,:), 'EdgeColor', 'none');
legend([h1, h2, h3], {'CSF', 'Gray Matter', 'White Matter'}, 'Location', 'northeast');

% Panel C: Combined view with tissue boundaries
subplot(1, 3, 3);
imagesc(fa_slice);
colormap(gca, gray);
axis equal tight;
set(gca, 'YDir', 'normal');
title('FA with Tissue Boundaries', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('X (voxels)'); ylabel('Y (voxels)');
hold on;

% Overlay tissue boundaries
contour(act_mask, [0.5, 1.5, 2.5], 'LineWidth', 2, 'LineColor', 'r');

sgtitle('ACT (Anatomically Constrained Tractography) Mask Example', 'FontSize', 16, 'FontWeight', 'bold');

%% Save figure
saveas(fig, fullfile(output_dir, '7_tractography_slice_visualization.png'));
fprintf('✓ Saved: %s\n', fullfile(output_dir, '7_tractography_slice_visualization.png'));

close(fig);
fprintf('✓ ACT mask visualization complete\n\n');

end
