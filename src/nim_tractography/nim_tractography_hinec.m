function tracks = nim_tractography_hinec(data_path, varargin)
% nim_tractography_hinec: High-Order Deterministic Tractography with ACT
%
% HIGH-ORDER TRACTOGRAPHY IMPLEMENTATION:
% - Trilinear interpolation of diffusion direction fields at sub-voxel positions
% - Runge-Kutta 4th order (RK4) numerical integration for smooth trajectories
% - Anatomically Constrained Tractography (ACT) using tissue segmentation
% - Continuous tracking with fixed step sizes (not voxel boundary jumps)
% - Gray matter termination and CSF avoidance for biological plausibility
%
% FOUR CORE ENHANCEMENTS:
% 1. INTERPOLATION: Smooth direction estimation between voxels
% 2. HIGH-ORDER INTEGRATION: RK2/RK4/RKF45 numerical methods
% 3. ADAPTIVE STEPPING: RKF45 with automatic error control (optional)
% 4. ACT CONSTRAINTS: Anatomical validity using WM/GM/CSF masks
%
% RKF45 ADAPTIVE INTEGRATION (integration_order=5, adaptive_step=true):
% - Dormand-Prince embedded RK pair (5th/4th order)
% - 7-stage integration with shared k_i evaluations
% - Automatic step size control based on local error estimate
% - Error tolerance: rkf_tolerance (default 0.01 voxels)
% - Step bounds: [step_min, step_max] (default [0.01, 1.0] voxels)
% - Safety factor: rkf_safety (default 0.9)
% - Provides superior accuracy with minimal user intervention
%
% Arguments:
%   data_path - Path to .mat file containing nim structure or nim structure itself
%   options - Structure containing tractography parameters (optional struct)
%
% Returns:
%   tracks - Cell array of fiber tracks (each track is Nx3 matrix)
%
% HIGH-ORDER TRACTOGRAPHY PARAMETERS:
%   step_size - Fixed integration step size in voxel units (default: 0.2)
%   interp_method - Interpolation method: 'trilinear' or 'none' (default: 'trilinear')
%   integration_order - 1=Euler, 2=RK2, 4=RK4, 5=RKF45 (default: 4)
%   adaptive_step - Enable RKF adaptive step sizing (default: false)
%   rkf_tolerance - Error tolerance for RKF in voxel units (default: 0.01)
%   rkf_safety - Safety factor for step adjustment (default: 0.9)
%   step_min - Minimum step size for RKF (default: 0.01 voxels)
%   step_max - Maximum step size for RKF (default: 1.0 voxels)
%   termination_fa - FA threshold for termination (default: 0.05)
%   angle_thresh - Maximum angle change between directions (default: 60°)
%   max_steps - Maximum integration steps (default: 5000)
%   seed_density - Seeds per voxel with flexible positioning (default: 1)
%
% ACT PARAMETERS:
%   act_enabled - Enable anatomically constrained tracking (default: true)
%   wm_mask - White matter mask for seeding and propagation
%   gm_mask - Gray matter mask for valid termination points
%   csf_mask - CSF mask for invalid termination (avoid CSF entry)
%
% USAGE EXAMPLES:
%   tracks = nim_tractography_hinec('data.mat'); % High-order with ACT
%   options.integration_order = 2; % Use RK2 instead of RK4
%   tracks = nim_tractography_hinec('data.mat', options);
%
%   % RKF adaptive step size (recommended for high accuracy)
%   options.integration_order = 5;
%   options.adaptive_step = true;
%   options.rkf_tolerance = 0.01;  % 0.01 voxel error tolerance
%   tracks = nim_tractography_hinec('data.mat', options);
%
% REFERENCES:
%   Basser et al. (2000). In vivo fiber tractography using DT-MRI data.
%   Magnetic Resonance in Medicine, 44(4), 625-632.
%
%   Smith et al. (2012). Anatomically-constrained tractography.
%   NeuroImage, 62(3), 1924-1938.
%
%   Dormand & Prince (1980). A family of embedded Runge-Kutta formulae.
%   Journal of Computational and Applied Mathematics, 6(1), 19-26.

% Parse input arguments
if nargin > 1 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default value
if ~isfield(options, 'seed_density')
    options.seed_density = 1;
end
if ~isfield(options, 'seed_strategy')
    options.seed_strategy = "uniform";
end
if ~isfield(options, 'step_size')
    options.step_size = 0.2;  % HINEC: Smaller step size for RK4 accuracy
end
if ~isfield(options, 'fa_threshold')
    options.fa_threshold = 0.1;
end
if ~isfield(options, 'angle_thresh')
    options.angle_thresh = 60;
end
if ~isfield(options, 'max_steps')
    options.max_steps = 5000;
end
if ~isfield(options, 'min_length')
    options.min_length = 10;
end
if ~isfield(options, 'termination_fa')
    options.termination_fa = 0.05;
end
if ~isfield(options, 'order')
    options.order = 1;  % Kept for backward compatibility
end
if ~isfield(options, 'interp_method')
    options.interp_method = "trilinear";  % HINEC: Use trilinear interpolation
end
if ~isfield(options, 'integration_order')
    options.integration_order = 4;  % HINEC: RK4 integration by default
end
if ~isfield(options, 'adaptive_step')
    options.adaptive_step = false;  % RKF: Disable adaptive stepping by default
end
if ~isfield(options, 'rkf_tolerance')
    options.rkf_tolerance = 0.01;  % RKF: Error tolerance in voxel units
end
if ~isfield(options, 'rkf_safety')
    options.rkf_safety = 0.9;  % RKF: Safety factor for step adjustment
end
if ~isfield(options, 'step_min')
    options.step_min = 0.01;  % RKF: Minimum step size (voxels)
end
if ~isfield(options, 'step_max')
    options.step_max = 1.0;  % RKF: Maximum step size (voxels)
end
if ~isfield(options, 'act_enabled')
    options.act_enabled = true;  % HINEC: Enable ACT by default
end
if ~isfield(options, 'wm_mask')
    options.wm_mask = [];  % HINEC: White matter mask for ACT
end
if ~isfield(options, 'gm_mask')
    options.gm_mask = [];  % HINEC: Gray matter mask for ACT
end
if ~isfield(options, 'csf_mask')
    options.csf_mask = [];  % HINEC: CSF mask for ACT
end
if ~isfield(options, 'seed_mask')
    options.seed_mask = [];
end
if ~isfield(options, 'enable_diagnostics')
    options.enable_diagnostics = true;  % Enable timing diagnostics
end
if ~isfield(options, 'use_fa_seed_filter')
    options.use_fa_seed_filter = false;
end
if ~isfield(options, 'inferior_slice_fraction')
    options.inferior_slice_fraction = 0.1;
end

options.seed_strategy = lower(string(options.seed_strategy));

% Initialize timing diagnostics
if options.enable_diagnostics
    timing = struct();
    timing.total_start = tic;
end

% Load data if path is provided
if ischar(data_path) || isstring(data_path)
    fprintf('Loading data from %s...\n', data_path);
    data = load(data_path);
    nim = data.nim;
else
    nim = data_path;
end

% Verify required fields
if ~isfield(nim, 'evec')
    error('Eigenvectors not found in nim structure. Please run nim_eig() first.');
end
if ~isfield(nim, 'FA')
    error('FA values not found in nim structure. Please run nim_fa() first.');
end

fprintf('Starting HINEC High-Order Tractography...\n');
fprintf('Parameters: step=%.2f, FA_thresh=%.2f, angle_thresh=%.1f\u00b0, integration_order=%d\n', ...
    options.step_size, options.fa_threshold, options.angle_thresh, options.integration_order);
fprintf('Integration: %s', get_integration_method_name(options.integration_order));
if options.integration_order == 5 && options.adaptive_step
    fprintf(' (adaptive, tol=%.4f)', options.rkf_tolerance);
end
fprintf('\n');
fprintf('Interpolation: %s, ACT: %s\n', options.interp_method, ...
    string(options.act_enabled && (~isempty(options.wm_mask) || ~isempty(options.gm_mask))));

% Get image dimensions
dims = size(nim.FA);

% HINEC: Pre-extract eigenvector components for efficient interpolation
if options.enable_diagnostics
    timing.precompute_start = tic;
end

fprintf('HINEC: Pre-extracting eigenvector components for interpolation...\n');
% Verify eigenvectors are available
if size(nim.evec, 4) ~= 3 || size(nim.evec, 5) ~= 3
    error('HINEC requires eigenvectors in format [h x w x d x 3 x 3]');
end

% Pre-extract primary eigenvector components for efficient interpolation
% This avoids repeated 5D array access during tracking
nim.v1_x = squeeze(nim.evec(:,:,:,1,1));  % X component of primary eigenvector
nim.v1_y = squeeze(nim.evec(:,:,:,2,1));  % Y component of primary eigenvector
nim.v1_z = squeeze(nim.evec(:,:,:,3,1));  % Z component of primary eigenvector

% Verify primary eigenvector at center voxel
center_idx = round(dims/2);
if isfield(nim, 'eval')
    center_eigenvals = squeeze(nim.eval(center_idx(1), center_idx(2), center_idx(3), :));
    if center_eigenvals(1) < center_eigenvals(2) || center_eigenvals(1) < center_eigenvals(3)
        warning('HINEC: Primary eigenvector may not correspond to largest eigenvalue!');
    end
    fprintf('HINEC: Primary eigenvalue at center: %.4f\n', center_eigenvals(1));
end

fprintf('HINEC: Extracted v1_x, v1_y, v1_z components (size: %dx%dx%d)\n', dims(1), dims(2), dims(3));

% HINEC OPTIMIZATION: Pre-create griddedInterpolant objects for fast repeated interpolation
% griddedInterpolant is 2-5x faster than interp3 for repeated queries on the same grid
fprintf('HINEC: Creating griddedInterpolant objects for fast interpolation...\n');

% Determine interpolation method
if strcmp(options.interp_method, 'cubic')
    interp_method = 'cubic';
else
    interp_method = 'linear';  % Default for 'trilinear' or any other
end
fprintf('HINEC: Interpolation method: %s\n', interp_method);

% Create grid vectors for griddedInterpolant
grid_vectors = {1:dims(1), 1:dims(2), 1:dims(3)};

% Pre-create interpolant objects (much faster than repeated interp3 calls)
nim.FA_interp = griddedInterpolant(grid_vectors, nim.FA, interp_method, 'none');
nim.v1_x_interp = griddedInterpolant(grid_vectors, nim.v1_x, interp_method, 'none');
nim.v1_y_interp = griddedInterpolant(grid_vectors, nim.v1_y, interp_method, 'none');
nim.v1_z_interp = griddedInterpolant(grid_vectors, nim.v1_z, interp_method, 'none');

fprintf('HINEC: griddedInterpolant objects created successfully\n');

if options.enable_diagnostics
    timing.precompute_time = toc(timing.precompute_start);
    fprintf('HINEC: Eigenvector extraction + interpolant setup took: %.2f seconds\n', timing.precompute_time);
end

% Verify seed mask was provided by caller
if isempty(options.seed_mask)
    error(['No seed mask provided. Seeding strategy must be configured in runTractography.m\n' ...
           'Example: options.seed_mask = brain_mask > 0.5;']);
end

% Ensure seed mask is logical
options.seed_mask = logical(options.seed_mask > 0);

options.seed_mask = logical(options.seed_mask);
fprintf('Pre-computing dilated brain mask for boundary checking...\n');
nim.dilated_brain_mask = imdilate(options.seed_mask, ones(3,3,3));

% Generate seed points
if options.enable_diagnostics
    timing.seed_start = tic;
end

[seed_points, seed_info] = generate_seed_points_hinec(options.seed_mask, options, dims);
fprintf('Seed strategy: %s\n', char(options.seed_strategy));
fprintf('Seeding layout: %s\n', seed_info.description);
fprintf('Seeds per seeded voxel: %.2f\n', seed_info.seeds_per_voxel);
if ~isnan(seed_info.voxel_spacing)
    fprintf('Approximate inter-seed spacing: %.3f voxels\n', seed_info.voxel_spacing);
end
fprintf('Generated %d seed points\n', size(seed_points, 1));

if options.enable_diagnostics
    timing.seed_time = toc(timing.seed_start);
    fprintf('Seed generation took: %.2f seconds\n', timing.seed_time);
end

% Print diagnostics
fprintf('=== TRACTOGRAPHY DIAGNOSTICS ===\n');
fprintf('Volume dimensions: %d x %d x %d\n', dims);
fprintf('Seed mask voxels: %d\n', sum(options.seed_mask(:)));
fprintf('Total seeds to process: %d\n', size(seed_points, 1));
fprintf('Estimated tracks: %d\n', size(seed_points, 1) * 2);
fprintf('==============================\n');

% Pre-allocate tracks
tracks = cell(size(seed_points, 1) * 2, 1);
track_count = 0;

% Convert angle threshold to cosine for efficiency
cos_angle_thresh = cos(deg2rad(options.angle_thresh));

% Initialize timing for tracking
if options.enable_diagnostics
    timing.tracking_start = tic;
end

% PARFOR PARALLELIZATION: Process seeds in parallel for massive speedup
% Each seed is independent - embarrassingly parallel problem
num_seeds = size(seed_points, 1);
fprintf('Processing %d seeds using parallel workers...\n', num_seeds);

% Pre-allocate sliced output arrays for parfor compatibility
% parfor requires sliced variables (one element per iteration)
all_tracks = cell(num_seeds, 1);           % Store combined tracks
term_fwd_all = cell(num_seeds, 1);         % Forward termination reasons
term_bwd_all = cell(num_seeds, 1);         % Backward termination reasons
track_valid = false(num_seeds, 1);         % Track validity flags
step_counts = zeros(num_seeds, 1);         % Step counts per seed
boundary_times = zeros(num_seeds, 1);      % Boundary check times per seed

% Check if Parallel Computing Toolbox is available
use_parfor = ~isempty(ver('parallel'));
if use_parfor
    % Ensure parallel pool is started
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('Starting parallel pool...\n');
        pool = parpool('local');
    end
    fprintf('Using %d parallel workers\n', pool.NumWorkers);
else
    fprintf('WARNING: Parallel Computing Toolbox not available. Using serial processing.\n');
end

% Progress tracking for parfor (using DataQueue)
if use_parfor
    progress_queue = parallel.pool.DataQueue;
    progress_count = 0;
    afterEach(progress_queue, @(~) fprintf('.'));
end

% Main parallel loop
tracking_start = tic;

if use_parfor
    parfor i = 1:num_seeds
        seed = seed_points(i, :);

        % Track in both directions
        [track_forward, step_timing_fwd, term_fwd] = track_fiber_hinec(nim, seed, +1, options, cos_angle_thresh);
        [track_backward, step_timing_bwd, term_bwd] = track_fiber_hinec(nim, seed, -1, options, cos_angle_thresh);

        % Store termination reasons
        term_fwd_all{i} = term_fwd;
        term_bwd_all{i} = term_bwd;

        % Store timing info
        step_counts(i) = step_timing_fwd.step_count + step_timing_bwd.step_count;
        boundary_times(i) = step_timing_fwd.boundary_time + step_timing_bwd.boundary_time;

        % Check if tracks were generated
        if isempty(track_forward) && isempty(track_backward)
            track_valid(i) = false;
            all_tracks{i} = [];
        else
            % Combine tracks: backward (flipped) + seed + forward
            if size(track_backward, 1) > 1
                track_backward = flipud(track_backward(2:end, :));
            else
                track_backward = [];
            end

            if size(track_forward, 1) > 1
                track_forward = track_forward(2:end, :);
            else
                track_forward = [];
            end

            % Combine into one continuous track
            combined_track = [track_backward; seed; track_forward];

            if size(combined_track, 1) > 1
                all_tracks{i} = combined_track;
                track_valid(i) = true;
            else
                all_tracks{i} = [];
                track_valid(i) = false;
            end
        end

        % Send progress update (every 100 seeds)
        if mod(i, 100) == 0
            send(progress_queue, i);
        end
    end
else
    % Serial fallback (no Parallel Computing Toolbox)
    for i = 1:num_seeds
        if mod(i, 100) == 0
            elapsed = toc(tracking_start);
            rate = i / elapsed;
            eta = (num_seeds - i) / rate;
            fprintf('\n%d/%d (%.1f seeds/s, ETA: %.1f min) ', i, num_seeds, rate, eta/60);
        end

        seed = seed_points(i, :);

        % Track in both directions
        [track_forward, step_timing_fwd, term_fwd] = track_fiber_hinec(nim, seed, +1, options, cos_angle_thresh);
        [track_backward, step_timing_bwd, term_bwd] = track_fiber_hinec(nim, seed, -1, options, cos_angle_thresh);

        % Store termination reasons
        term_fwd_all{i} = term_fwd;
        term_bwd_all{i} = term_bwd;

        % Store timing info
        step_counts(i) = step_timing_fwd.step_count + step_timing_bwd.step_count;
        boundary_times(i) = step_timing_fwd.boundary_time + step_timing_bwd.boundary_time;

        % Check if tracks were generated
        if isempty(track_forward) && isempty(track_backward)
            track_valid(i) = false;
            all_tracks{i} = [];
        else
            % Combine tracks: backward (flipped) + seed + forward
            if size(track_backward, 1) > 1
                track_backward = flipud(track_backward(2:end, :));
            else
                track_backward = [];
            end

            if size(track_forward, 1) > 1
                track_forward = track_forward(2:end, :);
            else
                track_forward = [];
            end

            % Combine into one continuous track
            combined_track = [track_backward; seed; track_forward];

            if size(combined_track, 1) > 1
                all_tracks{i} = combined_track;
                track_valid(i) = true;
            else
                all_tracks{i} = [];
                track_valid(i) = false;
            end
        end
    end
end

fprintf('\nParallel tracking completed.\n');

% Aggregate results from parallel workers
% Extract valid tracks
tracks = all_tracks(track_valid);
track_count = sum(track_valid);

% Count termination reasons (aggregate from all workers)
failure_reasons = struct();
failure_reasons.no_initial_direction = sum(cellfun(@isempty, term_fwd_all) & cellfun(@isempty, term_bwd_all));
failure_reasons.immediate_fa_fail = 0;
failure_reasons.no_boundary_exit = 0;
failure_reasons.short_tracks = 0;
failure_reasons.successful = track_count;

% HINEC ACT-specific termination reasons
failure_reasons.csf_termination = sum(strcmp(term_fwd_all, 'csf')) + sum(strcmp(term_bwd_all, 'csf'));
failure_reasons.gm_termination = sum(strcmp(term_fwd_all, 'gm')) + sum(strcmp(term_bwd_all, 'gm'));
failure_reasons.outside_termination = sum(strcmp(term_fwd_all, 'outside')) + sum(strcmp(term_bwd_all, 'outside'));
failure_reasons.fa_termination = sum(strcmp(term_fwd_all, 'fa')) + sum(strcmp(term_bwd_all, 'fa'));
failure_reasons.angle_termination = sum(strcmp(term_fwd_all, 'angle')) + sum(strcmp(term_bwd_all, 'angle'));

% Aggregate timing statistics
if options.enable_diagnostics
    timing.step_count = sum(step_counts);
    timing.boundary_time = sum(boundary_times);
    timing.interpolation_time = 0;  % Not tracked per-worker for simplicity
    timing.rkf_rejections = 0;      % Not tracked in parallel mode
    timing.rkf_retries = 0;         % Not tracked in parallel mode
end

% Print final timing report
if options.enable_diagnostics
    timing.tracking_time = toc(timing.tracking_start);
    timing.total_time = toc(timing.total_start);
    
    fprintf('\n\n=== HINEC TIMING REPORT ===\n');
    fprintf('Integration method: Order %d (%s)\n', options.integration_order, ...
        get_integration_method_name(options.integration_order));
    if options.integration_order == 5 && options.adaptive_step
        fprintf('Adaptive stepping: ENABLED (tolerance=%.4f voxels)\n', options.rkf_tolerance);
    end
    fprintf('Total time: %.2f seconds\n', timing.total_time);
    fprintf('Eigenvector extraction: %.2f seconds (%.1f%%)\n', timing.precompute_time, 100*timing.precompute_time/timing.total_time);
    fprintf('Seed generation: %.2f seconds (%.1f%%)\n', timing.seed_time, 100*timing.seed_time/timing.total_time);
    fprintf('HINEC tracking: %.2f seconds (%.1f%%)\n', timing.tracking_time, 100*timing.tracking_time/timing.total_time);
    fprintf('  - Interpolation + integration overhead included\n');
    fprintf('  - Brain boundary checks: %.2f seconds (%.1f%% of tracking)\n', timing.boundary_time, 100*timing.boundary_time/timing.tracking_time);
    if options.integration_order == 5 && options.adaptive_step && timing.rkf_rejections > 0
        fprintf('  - RKF step rejections: %d (%.1f%% of steps)\n', timing.rkf_rejections, 100*timing.rkf_rejections/timing.step_count);
        fprintf('  - RKF retry attempts: %d\n', timing.rkf_retries);
    end
    fprintf('Total integration steps: %d\n', timing.step_count);
    fprintf('Average steps per track: %.1f\n', timing.step_count / (size(seed_points, 1) * 2));
    fprintf('Integration steps per second: %.1f\n', timing.step_count / timing.tracking_time);
    fprintf('========================\n');
end

total_attempts = size(seed_points, 1) * 2;
success_rate = (track_count / total_attempts) * 100;
fprintf('\nGenerated %d valid tracks (filtered from %d total attempts - %.1f%% success rate)\n', track_count, total_attempts, success_rate);

% Detailed failure analysis
fprintf('\n=== TERMINATION ANALYSIS ===\n');
fprintf('No initial direction: %d (%.1f%%)\n', failure_reasons.no_initial_direction, 100*failure_reasons.no_initial_direction/total_attempts);
fprintf('Successful tracks: %d (%.1f%%)\n', failure_reasons.successful, 100*failure_reasons.successful/total_attempts);

% ACT-specific terminations
act_enabled = ~isempty(options.wm_mask) && ~isempty(options.gm_mask) && ~isempty(options.csf_mask);
if act_enabled
    fprintf('\n--- ACT Termination Statistics ---\n');
    fprintf('GM terminations (valid):   %d (%.1f%%)\n', failure_reasons.gm_termination, 100*failure_reasons.gm_termination/total_attempts);
    fprintf('CSF terminations (invalid): %d (%.1f%%)\n', failure_reasons.csf_termination, 100*failure_reasons.csf_termination/total_attempts);
    fprintf('Outside brain:             %d (%.1f%%)\n', failure_reasons.outside_termination, 100*failure_reasons.outside_termination/total_attempts);

    % ACT effectiveness
    total_act_terminations = failure_reasons.gm_termination + failure_reasons.csf_termination + failure_reasons.outside_termination;
    if total_act_terminations > 0
        gm_ratio = 100 * failure_reasons.gm_termination / total_act_terminations;
        fprintf('ACT effectiveness: %.1f%% valid GM terminations\n', gm_ratio);
    end
end

% Standard termination reasons
fprintf('\n--- Standard Termination Reasons ---\n');
fprintf('FA threshold:    %d (%.1f%%)\n', failure_reasons.fa_termination, 100*failure_reasons.fa_termination/total_attempts);
fprintf('Angle threshold: %d (%.1f%%)\n', failure_reasons.angle_termination, 100*failure_reasons.angle_termination/total_attempts);

fprintf('============================\n');

if success_rate < 10
    fprintf('⚠️  WARNING: Extremely low success rate! Check algorithm parameters.\n');
end

% Print summary
fprintf('\n========================================\n');
fprintf('HINEC HIGH-ORDER TRACTOGRAPHY COMPLETE\n');
fprintf('========================================\n');
fprintf('Algorithm: Interpolation + RK4 + ACT\n');
if act_enabled
    fprintf('ACT Status: ENABLED (WM/GM/CSF masks active)\n');
    fprintf('ACT Effectiveness: %.1f%% valid GM terminations\n', track_stats.act_effectiveness_pct);
else
    fprintf('ACT Status: DISABLED (no tissue masks)\n');
end
fprintf('Tracks Generated: %d\n', track_count);
fprintf('Mean Track Length: %.2f mm\n', track_stats.mean_length);
fprintf('Output File: %s\n', output_file);
fprintf('========================================\n');

end

function [seed_points, seed_info] = generate_seed_points_hinec(seed_mask, options, dims)
% HINEC: Generate seed points for high-order tractography
% Default strategy places a uniform lattice of seeds inside each seeded voxel

mask_idx = find(seed_mask);
if isempty(mask_idx)
    seed_points = zeros(0, 3);
    seed_info = struct('description', 'uniform (empty mask)', ...
                       'seeds_per_voxel', 0, ...
                       'voxel_spacing', NaN);
    return;
end

[x, y, z] = ind2sub(dims, mask_idx);
base_voxels = [x, y, z];
num_voxels = size(base_voxels, 1);

strategy = lower(string(options.seed_strategy));

density = max(1, options.seed_density);

if strcmp(strategy, "random")
    % Backwards-compatible stochastic seeding inside voxels
    if density <= 1
        offsets = (rand(num_voxels, 3) - 0.5) * 0.4;
        seed_points = base_voxels + offsets;
        seeds_per_voxel = 1;
    else
        density_int = max(1, round(density));
        seed_points = zeros(num_voxels * density_int, 3);
        idx = 1;
        for i = 1:num_voxels
            for j = 1:density_int
                offset = (rand(1, 3) - 0.5) * 0.8;
                seed_points(idx, :) = base_voxels(i, :) + offset;
                idx = idx + 1;
            end
        end
        seeds_per_voxel = density_int;
    end

    seed_info = struct('description', 'random jitter within seeded voxels', ...
                       'seeds_per_voxel', double(seeds_per_voxel), ...
                       'voxel_spacing', NaN);
    return;
end

% Deterministic lattice seeding inside each voxel
per_axis = max(1, ceil(density^(1/3)));
axis_edges = linspace(-0.5, 0.5, per_axis + 1);
axis_offsets = (axis_edges(1:end-1) + axis_edges(2:end)) / 2;
[ox, oy, oz] = ndgrid(axis_offsets, axis_offsets, axis_offsets);
offsets = [ox(:), oy(:), oz(:)];
num_offsets = size(offsets, 1);

seed_points = zeros(num_voxels * num_offsets, 3);
idx = 1;
for i = 1:num_voxels
    voxel_center = base_voxels(i, :);
    for j = 1:num_offsets
        seed_points(idx, :) = voxel_center + offsets(j, :);
        idx = idx + 1;
    end
end

voxel_spacing = 1 / per_axis;
description = sprintf('uniform grid (%d×%d×%d sub-voxel lattice)', per_axis, per_axis, per_axis);

seed_info = struct('description', description, ...
                   'seeds_per_voxel', double(num_offsets), ...
                   'voxel_spacing', voxel_spacing);
end

function [track, step_timing, termination_reason] = track_fiber_hinec(nim, seed, direction, options, cos_angle_thresh)
% HINEC: High-order deterministic tractography with interpolation, RK4, and ACT
%
% HINEC ALGORITHM ENHANCEMENTS:
% - Trilinear interpolation of diffusion direction fields (sub-voxel precision)
% - Runge-Kutta 4th order numerical integration (higher accuracy)
% - Anatomically Constrained Tractography (tissue-based termination)
% - Continuous tracking with fixed step size (no voxel boundary jumps)
%
% Arguments:
%   nim - NIM structure with diffusion tensors
%   seed - Starting position
%   direction - Initial tracking direction (+1 or -1)
%   options - HINEC tracking parameters
%   cos_angle_thresh - Cosine of maximum angle change
%
% Returns:
%   track - Array of positions along fiber track

% Initialize position and streamline
current_pos = seed;
track = zeros(options.max_steps + 1, 3);
track(1, :) = current_pos;
track_length = 1;

% Initialize termination reason tracking
if nargout > 2
    termination_reason = 'unknown';
end

% Initialize timing
if nargout > 1
    step_timing = struct();
    step_timing.interpolation_time = 0;  % Included in RK4 timing
    step_timing.boundary_time = 0;
    step_timing.step_count = 0;
end

% Get initial direction using interpolation
[dir_vec, fa_val] = interpolate_direction_trilinear(nim, current_pos, options);
if isempty(dir_vec) || fa_val < options.termination_fa
    track = track(1:track_length, :);
    if nargout > 2
        termination_reason = 'no_direction';
    end
    return;
end

% Apply direction flip for bidirectional tracking
dir_vec = dir_vec * direction;

% Pre-compute frequently used values
dims = size(nim.FA);
has_parcellation = isfield(nim, 'dilated_brain_mask');

% HINEC algorithm: continuous tracking with RK4 integration
while true
    if nargout > 1
        step_timing.step_count = step_timing.step_count + 1;
    end

    % Check termination criteria BEFORE advancing
    if fa_val < options.termination_fa
        if nargout > 2
            termination_reason = 'fa';
        end
        break;
    end

    % Check angle constraint (only after first step)
    if track_length > 1
        prev_pos = track(track_length-1, :);
        current_step = track(track_length, :) - prev_pos;
        if norm(current_step) > 1e-6
            current_step = current_step / norm(current_step);
            if dot(dir_vec, current_step) < cos_angle_thresh
                if nargout > 2
                    termination_reason = 'angle';
                end
                break;
            end
        end
    end

    % Check maximum steps
    if track_length >= options.max_steps
        if nargout > 2
            termination_reason = 'max_steps';
        end
        break;
    end

    % HINEC: Advance position using selected integration method (Euler/RK2/RK4/RKF45)
    % For RKF adaptive stepping, may need to retry with smaller step
    max_retries = 5;
    retry_count = 0;
    step_accepted = false;

    while ~step_accepted && retry_count < max_retries
        [next_pos, step_accepted, new_step_size] = advance_position(nim, current_pos, dir_vec, options);

        if ~step_accepted
            % Step rejected - reduce step size and retry
            retry_count = retry_count + 1;
            options.step_size = new_step_size;

            % Track RKF statistics (only for RKF adaptive mode)
            if options.enable_diagnostics && options.integration_order == 5 && options.adaptive_step
                step_timing.rkf_rejections = step_timing.rkf_rejections + 1;
                step_timing.rkf_retries = step_timing.rkf_retries + 1;
            end
        else
            % Step accepted - update step size for next iteration
            options.step_size = new_step_size;
        end
    end

    % If step still not accepted after max retries, terminate tracking
    if ~step_accepted
        if nargout > 2
            termination_reason = 'rkf_failure';
        end
        break;
    end

    % Check if next position is valid (within volume bounds)
    if any(next_pos < 1) || any(next_pos > dims)
        if nargout > 2
            termination_reason = 'outside';
        end
        break;
    end

    % Brain tissue check at next position
    if has_parcellation
        if nargout > 1
            boundary_tic = tic;
        end

        next_voxel = round(next_pos);
        if all(next_voxel >= 1) && all(next_voxel <= dims)
            if ~nim.dilated_brain_mask(next_voxel(1), next_voxel(2), next_voxel(3))
                if nargout > 2
                    termination_reason = 'outside';
                end
                break;
            end
        end

        if nargout > 1
            step_timing.boundary_time = step_timing.boundary_time + toc(boundary_tic);
        end
    end

    % HINEC ACT: Check tissue type at next position
    tissue_type = check_tissue_type(next_pos, options, dims);

    % ACT termination logic
    if strcmp(tissue_type, 'CSF')
        % Entered CSF - invalid termination, discard track
        % Set track_length to 0 to signal invalid track
        track_length = 0;
        if nargout > 2
            termination_reason = 'csf';
        end
        break;

    elseif strcmp(tissue_type, 'GM')
        % Reached gray matter - valid termination point
        % Add final position in GM and stop tracking
        track_length = track_length + 1;
        track(track_length, :) = next_pos;
        if nargout > 2
            termination_reason = 'gm';
        end
        break;

    elseif strcmp(tissue_type, 'OUTSIDE')
        % Left brain volume - stop tracking
        if nargout > 2
            termination_reason = 'outside';
        end
        break;

    elseif strcmp(tissue_type, 'WM') || strcmp(tissue_type, 'UNKNOWN')
        % White matter or no ACT - continue tracking normally
        % Add position to streamline
        track_length = track_length + 1;
        track(track_length, :) = next_pos;

    else
        % Unknown tissue type - stop tracking
        if nargout > 2
            termination_reason = 'unknown';
        end
        break;
    end

    % If track was invalidated by CSF, exit immediately
    if track_length == 0
        break;
    end

    % Move to next position and update direction
    current_pos = next_pos;

    % Get new direction using interpolation
    [new_dir, fa_val] = interpolate_direction_trilinear(nim, current_pos, options);
    if isempty(new_dir)
        if nargout > 2
            termination_reason = 'no_direction';
        end
        break;
    end

    % Ensure new direction is oriented consistently with current direction
    if dot(dir_vec, new_dir) < 0
        new_dir = -new_dir;
    end
    dir_vec = new_dir;
end

% Trim track array to actual length
track = track(1:track_length, :);
end

function name = get_integration_method_name(order)
% Get human-readable name for integration method
switch order
    case 1
        name = 'Euler';
    case 2
        name = 'RK2/Midpoint';
    case 4
        name = 'RK4';
    case 5
        name = 'RKF45 (Dormand-Prince)';
    otherwise
        name = sprintf('Unknown (order %d)', order);
end
end


function initial_dir = get_initial_direction_hinec(nim, pos, options)
% HINEC: Get initial direction from seed using interpolation
[initial_dir, fa_val] = interpolate_direction_trilinear(nim, pos, options);
if isempty(initial_dir) || fa_val < options.termination_fa
    initial_dir = [];
end
end


function [new_pos, error_est, success] = rkf45_integration_step(nim, pos, dir_vec, options)
% HINEC: Runge-Kutta-Fehlberg 4(5) adaptive integration step
%
% RKF45 (Dormand-Prince) EMBEDDED RK PAIR:
% - 7 stages with shared k_i evaluations
% - 5th-order solution (primary, used for advancement)
% - 4th-order solution (embedded, used for error estimation)
% - Adaptive step size based on local error estimate
%
% DORMAND-PRINCE RK5(4)7M COEFFICIENTS:
% Butcher tableau with c_i, a_ij, b_hat_i (5th order), b_i (4th order)
%
% Arguments:
%   nim - NIM structure with interpolation data
%   pos - Current position
%   dir_vec - Current direction (k1)
%   options - Contains step_size (h), rkf_tolerance, rkf_safety
%
% Returns:
%   new_pos - Next position using 5th-order RKF solution
%   error_est - Local error estimate (norm of difference)
%   success - Boolean indicating if error is within tolerance

% Dormand-Prince RK5(4)7M coefficients
% c_i coefficients (node positions)
c = [0; 1/5; 3/10; 4/5; 8/9; 1; 1];

% a_ij coefficients (Butcher tableau rows)
a = zeros(7, 7);
a(2, 1) = 1/5;
a(3, 1:2) = [3/40, 9/40];
a(4, 1:3) = [44/45, -56/15, 32/9];
a(5, 1:4) = [19372/6561, -25360/2187, 64448/6561, -212/729];
a(6, 1:5) = [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656];
a(7, 1:6) = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84];

% b_hat_i weights (5th-order solution - PRIMARY)
b_hat = [35/384; 0; 500/1113; 125/192; -2187/6784; 11/84; 0];

% b_i weights (4th-order solution - for error estimation)
b = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];

h = options.step_size;

% Stage 1: k1 = f(x_n, y_n) - already provided as dir_vec
k1 = dir_vec;

% Stage 2: k2 = f(x_n + c_2*h, y_n + h*(a_21*k1))
pos_k2 = pos + h * (a(2,1) * k1);
[k2, ~] = interpolate_direction_trilinear(nim, pos_k2, options);
if isempty(k2), k2 = k1; else, if dot(k2, k1) < 0, k2 = -k2; end; end

% Stage 3: k3 = f(x_n + c_3*h, y_n + h*(a_31*k1 + a_32*k2))
pos_k3 = pos + h * (a(3,1)*k1 + a(3,2)*k2);
[k3, ~] = interpolate_direction_trilinear(nim, pos_k3, options);
if isempty(k3), k3 = k2; else, if dot(k3, k1) < 0, k3 = -k3; end; end

% Stage 4: k4
pos_k4 = pos + h * (a(4,1)*k1 + a(4,2)*k2 + a(4,3)*k3);
[k4, ~] = interpolate_direction_trilinear(nim, pos_k4, options);
if isempty(k4), k4 = k3; else, if dot(k4, k1) < 0, k4 = -k4; end; end

% Stage 5: k5
pos_k5 = pos + h * (a(5,1)*k1 + a(5,2)*k2 + a(5,3)*k3 + a(5,4)*k4);
[k5, ~] = interpolate_direction_trilinear(nim, pos_k5, options);
if isempty(k5), k5 = k4; else, if dot(k5, k1) < 0, k5 = -k5; end; end

% Stage 6: k6
pos_k6 = pos + h * (a(6,1)*k1 + a(6,2)*k2 + a(6,3)*k3 + a(6,4)*k4 + a(6,5)*k5);
[k6, ~] = interpolate_direction_trilinear(nim, pos_k6, options);
if isempty(k6), k6 = k5; else, if dot(k6, k1) < 0, k6 = -k6; end; end

% Stage 7: k7
pos_k7 = pos + h * (a(7,1)*k1 + a(7,2)*k2 + a(7,3)*k3 + a(7,4)*k4 + a(7,5)*k5 + a(7,6)*k6);
[k7, ~] = interpolate_direction_trilinear(nim, pos_k7, options);
if isempty(k7), k7 = k6; else, if dot(k7, k1) < 0, k7 = -k7; end; end

% Compute 5th-order solution (PRIMARY - used for advancement)
y_hat = pos + h * (b_hat(1)*k1 + b_hat(2)*k2 + b_hat(3)*k3 + b_hat(4)*k4 + ...
                   b_hat(5)*k5 + b_hat(6)*k6 + b_hat(7)*k7);

% Compute 4th-order solution (for error estimation)
y_tilde = pos + h * (b(1)*k1 + b(2)*k2 + b(3)*k3 + b(4)*k4 + ...
                     b(5)*k5 + b(6)*k6 + b(7)*k7);

% Error estimate: difference between 5th and 4th order solutions
error_vec = y_hat - y_tilde;
error_est = norm(error_vec);

% Check if error is within tolerance
success = error_est <= options.rkf_tolerance;

% Return 5th-order solution (higher accuracy)
new_pos = y_hat;
end


function new_pos = rk4_integration_step(nim, pos, dir_vec, options)
% HINEC: Runge-Kutta 4th order integration step
%
% RK4 MATHEMATICAL FORMULATION:
% k1 = v(r_n)
% k2 = v(r_n + 0.5*h*k1)
% k3 = v(r_n + 0.5*h*k2)
% k4 = v(r_n + h*k3)
% r_{n+1} = r_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
%
% Arguments:
%   nim - NIM structure with interpolation data
%   pos - Current position
%   dir_vec - Current direction (k1)
%   options - Contains step_size (h)
%
% Returns:
%   new_pos - Next position using RK4 integration

h = options.step_size;

% k1: direction at current position (already provided as dir_vec)
k1 = dir_vec;

% k2: direction at position + 0.5*h*k1
pos_k2 = pos + 0.5 * h * k1;
[k2, ~] = interpolate_direction_trilinear(nim, pos_k2, options);
if isempty(k2)
    % Fallback to k1 if interpolation fails
    k2 = k1;
else
    % Ensure direction consistency
    if dot(k2, dir_vec) < 0
        k2 = -k2;
    end
end

% k3: direction at position + 0.5*h*k2
pos_k3 = pos + 0.5 * h * k2;
[k3, ~] = interpolate_direction_trilinear(nim, pos_k3, options);
if isempty(k3)
    % Fallback to k2 if interpolation fails
    k3 = k2;
else
    % Ensure direction consistency
    if dot(k3, dir_vec) < 0
        k3 = -k3;
    end
end

% k4: direction at position + h*k3
pos_k4 = pos + h * k3;
[k4, ~] = interpolate_direction_trilinear(nim, pos_k4, options);
if isempty(k4)
    % Fallback to k3 if interpolation fails
    k4 = k3;
else
    % Ensure direction consistency
    if dot(k4, dir_vec) < 0
        k4 = -k4;
    end
end

% RK4 weighted combination
new_pos = pos + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end


function [new_pos, step_accepted, new_step_size] = advance_position(nim, pos, dir_vec, options)
% HINEC: Advance position using selected integration method
%
% INTEGRATION METHODS:
% Order 1 (Euler): r_{n+1} = r_n + h*v(r_n)
% Order 2 (RK2):   Midpoint method
% Order 4 (RK4):   Full Runge-Kutta 4th order
% Order 5 (RKF45): Adaptive Runge-Kutta-Fehlberg with error control
%
% Arguments:
%   nim - NIM structure
%   pos - Current position
%   dir_vec - Current direction
%   options - Contains integration_order, step_size, adaptive_step
%
% Returns:
%   new_pos - Next position
%   step_accepted - Boolean (true for fixed methods, adaptive for RKF)
%   new_step_size - Updated step size (adaptive for RKF, unchanged otherwise)

h = options.step_size;
step_accepted = true;  % Default: always accept for non-adaptive methods
new_step_size = h;     % Default: no change

switch options.integration_order
    case 1
        % Euler integration (first order)
        new_pos = pos + h * dir_vec;

    case 2
        % RK2 / Midpoint method (second order)
        k1 = dir_vec;
        mid_pos = pos + 0.5 * h * k1;

        [k2, ~] = interpolate_direction_trilinear(nim, mid_pos, options);
        if ~isempty(k2)
            % Ensure direction consistency
            if dot(k2, dir_vec) < 0
                k2 = -k2;
            end
            new_pos = pos + h * k2;
        else
            % Fallback to Euler if interpolation fails
            new_pos = pos + h * k1;
        end

    case 4
        % RK4 integration (fourth order)
        new_pos = rk4_integration_step(nim, pos, dir_vec, options);

    case 5
        % RKF45 integration (fifth order with adaptive step size)
        if options.adaptive_step
            % Adaptive RKF with step size control
            [new_pos, error_est, success] = rkf45_integration_step(nim, pos, dir_vec, options);

            % Compute new step size based on error estimate
            if error_est > 1e-10  % Avoid division by zero
                % h_new = safety * h * (tol / err)^(1/5)
                % Exponent 1/5 because local error scales as h^5 for 4th-order method
                scale_factor = (options.rkf_tolerance / error_est)^(1/5);
                new_step_size = options.rkf_safety * h * scale_factor;
            else
                % Error is negligible - allow step increase
                new_step_size = options.rkf_safety * h * 2.0;
            end

            % Apply step size bounds
            new_step_size = max(options.step_min, min(options.step_max, new_step_size));

            % Limit growth rate (max 2x increase per step for stability)
            new_step_size = min(new_step_size, 2.0 * h);

            % Accept or reject step based on error tolerance
            step_accepted = success;

            if ~step_accepted
                % Step rejected - use smaller step for retry
                % Return original position (caller will retry)
                new_pos = pos;
            end
        else
            % Fixed step RKF (no adaptivity)
            [new_pos, ~, ~] = rkf45_integration_step(nim, pos, dir_vec, options);
        end

    otherwise
        % Default to Euler for unknown orders
        warning('Unknown integration order %d, using Euler', options.integration_order);
        new_pos = pos + h * dir_vec;
end
end


function [direction, fa_value] = interpolate_direction_trilinear(nim, pos, options)
% HINEC: Fast interpolation of direction and FA using pre-created griddedInterpolant
%
% OPTIMIZED INTERPOLATION:
% - Uses pre-created griddedInterpolant objects (2-5x faster than interp3)
% - Supports both 'linear' (trilinear) and 'cubic' interpolation
% - Interpolates primary eigenvector components at continuous positions
% - Provides smooth direction transitions across voxel boundaries
%
% Arguments:
%   nim - NIM structure with pre-created interpolant objects:
%         nim.FA_interp, nim.v1_x_interp, nim.v1_y_interp, nim.v1_z_interp
%   pos - Current position (continuous, sub-voxel precision)
%   options - Tractography parameters
%
% Returns:
%   direction - Interpolated primary eigenvector [3x1]
%   fa_value - Interpolated FA value

direction = [];
fa_value = 0;

dims = size(nim.FA);

% Boundary check: need buffer for interpolation
% griddedInterpolant with 'none' extrapolation returns NaN outside bounds
% Use slightly tighter bounds for cubic interpolation
if strcmp(options.interp_method, 'cubic')
    margin = 1.5;  % Cubic needs more margin
else
    margin = 1.1;  % Linear needs less margin
end

if any(pos < margin) || pos(1) > dims(1)-margin+1 || ...
   pos(2) > dims(2)-margin+1 || pos(3) > dims(3)-margin+1
    return;
end

% Interpolate FA value first using pre-created interpolant (fast termination check)
% griddedInterpolant uses (X, Y, Z) ordering directly - no coordinate swap needed
try
    fa_value = nim.FA_interp(pos(1), pos(2), pos(3));
catch
    return;
end

% Check for NaN (outside bounds) or too low FA
if isnan(fa_value) || fa_value < options.termination_fa
    fa_value = 0;
    return;
end

% Interpolate primary eigenvector components using pre-created interpolants
try
    v_x = nim.v1_x_interp(pos(1), pos(2), pos(3));
    v_y = nim.v1_y_interp(pos(1), pos(2), pos(3));
    v_z = nim.v1_z_interp(pos(1), pos(2), pos(3));

    % Check for NaN values (outside interpolation domain)
    if isnan(v_x) || isnan(v_y) || isnan(v_z)
        return;
    end

    direction = [v_x, v_y, v_z];

    % Validate and normalize
    dir_norm = norm(direction);
    if dir_norm > 1e-6 && ~any(isinf(direction))
        direction = direction / dir_norm;
    else
        direction = [];
    end
catch
    % Interpolation failed - return empty
    direction = [];
end
end

function tissue_type = check_tissue_type(pos, options, dims)
% HINEC ACT: Check tissue type at given position
% Returns: 'WM', 'GM', 'CSF', 'OUTSIDE', or 'UNKNOWN'
%
% ACT Logic:
%   - WM: Valid for tracking (continue)
%   - GM: Valid termination point (stop tracking, keep track)
%   - CSF: Invalid termination point (stop tracking, discard track)
%   - OUTSIDE: Outside brain volume (stop tracking, discard track)
%   - UNKNOWN: No tissue masks available (use FA-based termination only)

tissue_type = 'UNKNOWN';

% Check if ACT is enabled (tissue masks available)
if isempty(options.wm_mask) || isempty(options.gm_mask) || isempty(options.csf_mask)
    return;
end

% Round position to nearest voxel for tissue lookup
voxel_pos = round(pos);

% Check bounds
if any(voxel_pos < 1) || voxel_pos(1) > dims(1) || ...
   voxel_pos(2) > dims(2) || voxel_pos(3) > dims(3)
    tissue_type = 'OUTSIDE';
    return;
end

% Query tissue masks at voxel position (binary masks)
% Check in priority order: CSF > GM > WM
% (CSF and GM are termination conditions, WM allows continuation)

try
    % Convert to linear index for efficient access
    linear_idx = sub2ind(dims, voxel_pos(1), voxel_pos(2), voxel_pos(3));

    % Check CSF first (highest priority termination)
    if options.csf_mask(linear_idx) > 0.5
        tissue_type = 'CSF';
        return;
    end

    % Check GM (valid termination)
    if options.gm_mask(linear_idx) > 0.5
        tissue_type = 'GM';
        return;
    end

    % Check WM (continue tracking)
    if options.wm_mask(linear_idx) > 0.5
        tissue_type = 'WM';
        return;
    end

    % Position not classified in any tissue mask
    tissue_type = 'UNKNOWN';

catch ME
    % Error accessing masks - treat as outside
    tissue_type = 'OUTSIDE';
end

end
