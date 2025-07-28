function tracks = nim_tractography_standard(data_path, varargin)
% nim_tractography_standard: Standard deterministic tractography (OPTIMIZED)
%
% Arguments:
%   data_path - Path to .mat file containing nim structure or nim structure itself
%   options - Structure containing tractography parameters (optional struct)
%
% Returns:
%   tracks - Cell array of fiber tracks (each track is Nx3 matrix)

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
if ~isfield(options, 'step_size')
    options.step_size = 0.2;
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
    options.order = 1;
end
if ~isfield(options, 'interp_method')
    options.interp_method = "linear";
end
if ~isfield(options, 'seed_mask')
    options.seed_mask = [];
end
if ~isfield(options, 'enable_diagnostics')
    options.enable_diagnostics = true;  % Enable timing diagnostics
end

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

fprintf('Starting standard tractography...\n');
fprintf('Parameters: step=%.2f, FA_thresh=%.2f, angle_thresh=%.1f\n', ...
        options.step_size, options.fa_threshold, options.angle_thresh);

% Get image dimensions
dims = size(nim.FA);

% OPTIMIZATION 1: Pre-compute eigenvector components
if options.enable_diagnostics
    timing.precompute_start = tic;
end

fprintf('Pre-computing eigenvector components...\n');
% FIX: Ensure we're extracting PRIMARY eigenvector (largest eigenvalue)
nim.v1_x = squeeze(nim.evec(:,:,:,1,1));  % First component of first eigenvector
nim.v1_y = squeeze(nim.evec(:,:,:,2,1));  % Second component of first eigenvector  
nim.v1_z = squeeze(nim.evec(:,:,:,3,1));  % Third component of first eigenvector

% Verify eigenvector extraction at center voxel
center_idx = round(dims/2);
if isfield(nim, 'eval')
    center_eigenvals = squeeze(nim.eval(center_idx(1), center_idx(2), center_idx(3), :));
    if center_eigenvals(1) < center_eigenvals(2) || center_eigenvals(1) < center_eigenvals(3)
        warning('Primary eigenvector may not correspond to largest eigenvalue!');
    end
end

if options.enable_diagnostics
    timing.precompute_time = toc(timing.precompute_start);
    fprintf('Pre-computation took: %.2f seconds\n', timing.precompute_time);
end

% Create seed mask if not provided
if isempty(options.seed_mask)
    options.seed_mask = nim.FA > options.fa_threshold;
    
    % Exclude non-brain areas if parcellation is available
    if isfield(nim, 'parcellation_mask')
        brain_mask = nim.parcellation_mask > 0;
        
        % FIX: Exclude bottom slices to avoid inferior brain artifacts
        z_exclude = max(1, round(dims(3) * 0.1)); % Exclude bottom 10% of slices
        brain_mask(:, :, 1:z_exclude) = 0;
        fprintf('Excluded bottom %d slices to avoid susceptibility artifacts\n', z_exclude);
        
        options.seed_mask = options.seed_mask & brain_mask;
        fprintf('Applied brain mask from parcellation (excluding label 0)\n');
        
        % OPTIMIZATION 2: Pre-compute dilated brain mask for boundary checking
        fprintf('Pre-computing dilated brain mask...\n');
        nim.dilated_brain_mask = imdilate(brain_mask, ones(3,3,3));
    end
end

% Generate seed points
if options.enable_diagnostics
    timing.seed_start = tic;
end

seed_points = generate_seed_points_standard_optimized(options.seed_mask, options.seed_density, dims);
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
    timing.interpolation_time = 0;
    timing.boundary_time = 0;
    timing.step_count = 0;
end

% Process each seed point
fprintf('Processing seeds: ');
last_report_time = tic;

for i = 1:size(seed_points, 1)
    % Progress reporting with time estimate
    if mod(i, 10) == 0
        elapsed = toc(timing.tracking_start);
        rate = i / elapsed;
        eta = (size(seed_points, 1) - i) / rate;
        fprintf('\n%d/%d (%.1f seeds/s, ETA: %.1f min) ', i, size(seed_points, 1), rate, eta/60);
        m = memory;
        fprintf('\nMemory: %.1f GB used', m.MemUsedMATLAB/1e9);
    end
    
    seed = seed_points(i, :);
    
    % Track in both directions
    for direction = [-1, 1]
        if options.enable_diagnostics
            [track, step_timing] = track_fiber_standard_optimized(nim, seed, direction, options, cos_angle_thresh);
            timing.interpolation_time = timing.interpolation_time + step_timing.interpolation_time;
            timing.boundary_time = timing.boundary_time + step_timing.boundary_time;
            timing.step_count = timing.step_count + step_timing.step_count;
        else
            track = track_fiber_standard_optimized(nim, seed, direction, options, cos_angle_thresh);
        end
        
        % Calculate track length in mm (CORRECTED)
        if size(track, 1) > 1
            % Vectorized length calculation
            track_length_mm = sum(vecnorm(diff(track), 2, 2));
            
            % Only keep tracks above minimum length
            if track_length_mm >= options.min_length
                track_count = track_count + 1;
                tracks{track_count} = track;
            end
        end
    end
end

% Trim tracks array
tracks = tracks(1:track_count);

% Print final timing report
if options.enable_diagnostics
    timing.tracking_time = toc(timing.tracking_start);
    timing.total_time = toc(timing.total_start);
    
    fprintf('\n\n=== TIMING REPORT ===\n');
    fprintf('Total time: %.2f seconds\n', timing.total_time);
    fprintf('Pre-computation: %.2f seconds (%.1f%%)\n', timing.precompute_time, 100*timing.precompute_time/timing.total_time);
    fprintf('Seed generation: %.2f seconds (%.1f%%)\n', timing.seed_time, 100*timing.seed_time/timing.total_time);
    fprintf('Tracking: %.2f seconds (%.1f%%)\n', timing.tracking_time, 100*timing.tracking_time/timing.total_time);
    fprintf('  - Interpolation: %.2f seconds (%.1f%% of tracking)\n', timing.interpolation_time, 100*timing.interpolation_time/timing.tracking_time);
    fprintf('  - Boundary checks: %.2f seconds (%.1f%% of tracking)\n', timing.boundary_time, 100*timing.boundary_time/timing.tracking_time);
    fprintf('Total steps processed: %d\n', timing.step_count);
    fprintf('Average steps per track: %.1f\n', timing.step_count / (size(seed_points, 1) * 2));
    fprintf('Steps per second: %.1f\n', timing.step_count / timing.tracking_time);
    fprintf('====================\n');
end

fprintf('\nGenerated %d valid tracks\n', track_count);

% SAVE RESULTS AUTOMATICALLY
output_dir = 'tractography_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save tracks with metadata
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
output_file = fullfile(output_dir, sprintf('tracks_%s.mat', timestamp));

% Calculate track statistics
if track_count > 0
    track_lengths = zeros(track_count, 1);
    for i = 1:track_count
        if size(tracks{i}, 1) > 1
            track_lengths(i) = sum(vecnorm(diff(tracks{i}), 2, 2));
        end
    end
    
    track_stats = struct();
    track_stats.num_tracks = track_count;
    track_stats.mean_length = mean(track_lengths);
    track_stats.median_length = median(track_lengths);
    track_stats.max_length = max(track_lengths);
    track_stats.min_length = min(track_lengths);
    track_stats.total_length = sum(track_lengths);
else
    track_lengths = [];
    track_stats = struct('num_tracks', 0);
end

% Save everything
save(output_file, 'tracks', 'options', 'track_stats', 'track_lengths', 'dims');
fprintf('Results saved to: %s\n', output_file);

end

function seed_points = generate_seed_points_standard_optimized(seed_mask, density, dims)
% OPTIMIZED: Generate seed points with pre-allocation
[x, y, z] = ind2sub(dims, find(seed_mask));
base_seeds = [x, y, z];
n_base = size(base_seeds, 1);

if density <= 1
    seed_points = base_seeds;
else
    % Pre-allocate for efficiency
    n_total = n_base * density;
    seed_points = zeros(n_total, 3);
    
    idx = 1;
    for i = 1:n_base
        for j = 1:density
            offset = (rand(1, 3) - 0.5) * 0.8;
            seed_points(idx, :) = base_seeds(i, :) + offset;
            idx = idx + 1;
        end
    end
end
end

function [track, step_timing] = track_fiber_standard_optimized(nim, seed, direction, options, cos_angle_thresh)
% OPTIMIZED: Track with improved boundary checking and diagnostics
track = seed;
current_pos = seed;
prev_direction = [];

% Initialize timing
if nargout > 1
    step_timing = struct();
    step_timing.interpolation_time = 0;
    step_timing.boundary_time = 0;
    step_timing.step_count = 0;
end

% Get initial direction
initial_dir = get_initial_direction_optimized(nim, seed, options);
if isempty(initial_dir)
    return;
end

% Apply direction flip for bidirectional tracking
initial_dir = initial_dir * direction;
prev_direction = initial_dir;

% Pre-compute frequently used values
dims = size(nim.FA);
has_parcellation = isfield(nim, 'dilated_brain_mask');

for step = 1:options.max_steps
    % Interpolation timing
    if nargout > 1
        interp_tic = tic;
    end
    
    % Get interpolated direction at current position
    [dir_vec, fa_val] = interpolate_direction_standard_optimized(nim, current_pos, options);
    
    if nargout > 1
        step_timing.interpolation_time = step_timing.interpolation_time + toc(interp_tic);
        step_timing.step_count = step_timing.step_count + 1;
    end
    
    % Termination criteria
    if isempty(dir_vec) || fa_val < options.termination_fa
        break;
    end
    
    % Ensure consistent direction (flip if needed)
    if ~isempty(prev_direction)
        if dot(dir_vec, prev_direction) < 0
            dir_vec = -dir_vec;
        end
        
        % Check curvature constraint
        if dot(dir_vec, prev_direction) < cos_angle_thresh
            break;
        end
    end
    
    % Step integration (simplified for performance testing)
    current_pos = current_pos + dir_vec * options.step_size;
    
    % Check bounds
    if any(current_pos < 1.5) || any(current_pos > dims - 0.5)
        break;
    end
    
    % OPTIMIZED: Simplified brain tissue check
    if has_parcellation
        if nargout > 1
            boundary_tic = tic;
        end
        
        pos_int = round(current_pos);
        if all(pos_int >= 1) && all(pos_int <= dims)
            % Use pre-computed dilated mask
            if ~nim.dilated_brain_mask(pos_int(1), pos_int(2), pos_int(3))
                break;
            end
        end
        
        if nargout > 1
            step_timing.boundary_time = step_timing.boundary_time + toc(boundary_tic);
        end
    end
    
    track = [track; current_pos];
    prev_direction = dir_vec;
end
end

function initial_dir = get_initial_direction_optimized(nim, pos, options)
[initial_dir, fa_val] = interpolate_direction_standard_optimized(nim, pos, options);
if isempty(initial_dir) || fa_val < options.termination_fa
    initial_dir = [];
end
end

function [direction, fa_value] = interpolate_direction_standard_optimized(nim, pos, options)
% OPTIMIZED: Use pre-computed eigenvector components
direction = [];
fa_value = 0;

% Check bounds
dims = size(nim.FA);
if any(pos < 1.1) || any(pos > dims - 0.1)
    return;
end

% Get FA value
try
    fa_value = interp3(nim.FA, pos(2), pos(1), pos(3), 'linear', 0);
catch
    return;
end

if fa_value < options.termination_fa
    return;
end

% OPTIMIZED: Use pre-computed components
try
    x_comp = interp3(nim.v1_x, pos(2), pos(1), pos(3), 'linear', 0);
    y_comp = interp3(nim.v1_y, pos(2), pos(1), pos(3), 'linear', 0);
    z_comp = interp3(nim.v1_z, pos(2), pos(1), pos(3), 'linear', 0);
    
    direction = [x_comp, y_comp, z_comp];
    
    % Normalize
    dir_norm = norm(direction);
    if dir_norm > 1e-6 && ~any(isnan(direction)) && ~any(isinf(direction))
        direction = direction / dir_norm;
    else
        direction = [];
    end
catch
    direction = [];
end
end