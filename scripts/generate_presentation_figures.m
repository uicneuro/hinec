function generate_presentation_figures(output_dir)
% generate_presentation_figures: Main executor for all presentation figures
%
% Generates all methodological figures for HINEC academic presentation:
% 1-2. Interpolation methods (FACT vs Trilinear vs Cubic)
% 3-4. Integration methods (RK4 vs RKF45)
% 5-6. Tractography workflow
%
% Usage:
%   generate_presentation_figures()                    % Uses default 'presentation_figures'
%   generate_presentation_figures('my_figures')        % Uses custom directory
%
% Output:
%   Saves 6 high-quality PNG figures to output_dir

if nargin < 1
    output_dir = 'presentation_figures';
end

% Add presentation code to path
addpath('src/nim_presentation');

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  HINEC Presentation Figure Generator                         ║\n');
fprintf('║  Academic Poster - Methodology Figures                       ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

start_time = tic;

%% 1. Interpolation Methods
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('SECTION 1: Interpolation Methods\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('\n');

try
    visualize_interpolation_methods(output_dir);
    fprintf('✓ Interpolation figures generated successfully\n\n');
catch ME
    fprintf('✗ Error generating interpolation figures:\n');
    fprintf('  %s\n\n', ME.message);
end

%% 2. Integration Methods
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('SECTION 2: Integration Methods (RK4 vs RKF45)\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('\n');

try
    visualize_integration_methods(output_dir);
    fprintf('✓ Integration figures generated successfully\n\n');
catch ME
    fprintf('✗ Error generating integration figures:\n');
    fprintf('  %s\n\n', ME.message);
end

%% 3. Tractography Example
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('SECTION 3: Tractography Step-by-Step Example\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('\n');

try
    visualize_tractography_example(output_dir);
    fprintf('✓ Tractography figures generated successfully\n\n');
catch ME
    fprintf('✗ Error generating tractography figures:\n');
    fprintf('  %s\n\n', ME.message);
end

%% 4. Real Tractography Slice
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('SECTION 4: Real Tractography Slice Visualization\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('\n');

try
    visualize_tractography_slice(output_dir);
    fprintf('✓ Tractography slice visualization generated successfully\n\n');
catch ME
    fprintf('✗ Error generating tractography slice visualization:\n');
    fprintf('  %s\n\n', ME.message);
end

%% Summary
elapsed_time = toc(start_time);

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  Generation Complete                                         ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('Output directory: %s\n', fullfile(pwd, output_dir));
fprintf('Total execution time: %.1f seconds\n\n', elapsed_time);

fprintf('Generated figures:\n');
fprintf('  1. 1_interpolation_comparison.png\n');
fprintf('     → FACT vs Trilinear vs Cubic (3 panels)\n\n');

fprintf('  2. 2_interpolation_field_comparison.png\n');
fprintf('     → Direction field detail (3 panels)\n\n');

fprintf('  3. 3_integration_method_comparison.png\n');
fprintf('     → RK4 vs RKF45 visual comparison\n\n');

fprintf('  4. 4_step_size_adaptation.png\n');
fprintf('     → Fixed vs adaptive step size\n\n');

fprintf('  5. 5_tractography_step_by_step.png\n');
fprintf('     → Complete workflow (6 panels)\n\n');

fprintf('  6. 6_tractography_3d_visualization.png\n');
fprintf('     → 3D multi-track visualization\n\n');

fprintf('  7. 7_tractography_slice_visualization.png\n');
fprintf('     → Real tractography results with tissue classification\n\n');

fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('USAGE IN PRESENTATION:\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
fprintf('\n');
fprintf('Methodology Section:\n');
fprintf('  • Figures 1-2: Interpolation methods\n');
fprintf('  • Figures 3-4: Integration methods\n');
fprintf('  • Figures 5-6: Tractography workflow\n');
fprintf('  • Figure 7: Real tractography results\n\n');

fprintf('Recommended Layout:\n');
fprintf('  Top row:    Figures 1, 3, 5\n');
fprintf('  Bottom row: Figures 2, 4, 6\n\n');

fprintf('Figure Quality:\n');
fprintf('  • Resolution: High-DPI suitable for poster printing\n');
fprintf('  • Format: PNG with transparent backgrounds where applicable\n');
fprintf('  • Colors: Optimized for both screen and print\n\n');

fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

end
