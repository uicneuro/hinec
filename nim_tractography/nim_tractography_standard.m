function tracks = nim_tractography_standard(data_path, varargin)
% nim_tractography_standard: FACT (Fiber Assignment by Continuous Tracking) Tractography
%
% FACT ALGORITHM IMPLEMENTATION:
% - Uses discrete voxel tensors (NO interpolation)
% - Voxel-to-voxel stepping (NO sub-voxel positioning)
% - Primary eigenvector direction from nearest voxel center
% - True deterministic tractography without smoothing
%
% DIFFERENCES FROM INTERPOLATED METHODS:
% - Direction sampling: Voxel center only (vs. linear/spline interpolation)
% - Position updates: Voxel-to-voxel (vs. sub-voxel positioning)
% - Computational cost: Lower (no interpolation calculations)
% - Track characteristics: More angular, discrete steps
%
% Arguments:
%   data_path - Path to .mat file containing nim structure or nim structure itself
%   options - Structure containing tractography parameters (optional struct)
%
% Returns:
%   tracks - Cell array of fiber tracks (each track is Nx3 matrix)
%
% FACT PARAMETERS:
%   termination_fa - FA threshold for termination (default: 0.05)
%   angle_thresh - Maximum angle change between voxels (default: 60°)
%   max_steps - Maximum tracking steps (default: 5000)
%   step_size - Should be 1.0 for pure FACT (voxel-to-voxel)
%
% USAGE EXAMPLES:
%   tracks = nim_tractography_standard('data.mat'); % Basic FACT
%   options.termination_fa = 0.02; % Lower threshold
%   tracks = nim_tractography_standard('data.mat', options);
%
% REFERENCES:
%   Mori et al. (1999). Three-dimensional tracking of axonal projections
%   in the brain by magnetic resonance imaging. Ann Neurol, 45(2), 265-269.

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
    options.step_size = 1.0;  % FACT: Use 1.0 for voxel-to-voxel stepping
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
    options.interp_method = "none";  % FACT: No interpolation used
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

fprintf('Starting FACT tractography...\n');
fprintf('FACT Parameters: step=%.2f (voxel units), FA_thresh=%.2f, angle_thresh=%.1f\u00b0\n', ...
    options.step_size, options.fa_threshold, options.angle_thresh);
fprintf('FACT Mode: Discrete voxel tensors, no interpolation\n');

% Get image dimensions
dims = size(nim.FA);

% FACT: Verify eigenvector structure (no pre-computation needed)
if options.enable_diagnostics
    timing.precompute_start = tic;
end

fprintf('FACT: Verifying eigenvector structure...\n');
% Verify eigenvectors are available for FACT algorithm
if size(nim.evec, 4) ~= 3 || size(nim.evec, 5) ~= 3
    error('FACT requires eigenvectors in format [h x w x d x 3 x 3]');
end

% Verify primary eigenvector at center voxel
center_idx = round(dims/2);
if isfield(nim, 'eval')
    center_eigenvals = squeeze(nim.eval(center_idx(1), center_idx(2), center_idx(3), :));
    if center_eigenvals(1) < center_eigenvals(2) || center_eigenvals(1) < center_eigenvals(3)
        warning('FACT: Primary eigenvector may not correspond to largest eigenvalue!');
    end
    fprintf('FACT: Primary eigenvalue at center: %.4f\n', center_eigenvals(1));
end

if options.enable_diagnostics
    timing.precompute_time = toc(timing.precompute_start);
    fprintf('FACT verification took: %.2f seconds\n', timing.precompute_time);
end

% Create seed mask if not provided
if isempty(options.seed_mask)
    options.seed_mask = nim.FA > options.fa_threshold;
    
    % Priority 1: Use the preprocessed brain mask if available
    brain_mask = [];
    if isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
        brain_mask = nim.mask > 0.5;
        fprintf('Using preprocessed brain mask from nim.mask\n');
    elseif isfield(nim, 'parcellation_mask')
        % Fallback: Use parcellation mask if no preprocessed brain mask
        brain_mask = nim.parcellation_mask > 0;
        fprintf('Using parcellation mask as brain mask (fallback)\n');
    end
    
    if ~isempty(brain_mask)
        % FIX: Exclude bottom slices to avoid inferior brain artifacts
        z_exclude = max(1, round(dims(3) * 0.1)); % Exclude bottom 10% of slices
        brain_mask(:, :, 1:z_exclude) = 0;
        fprintf('Excluded bottom %d slices to avoid susceptibility artifacts\n', z_exclude);
        
        options.seed_mask = options.seed_mask & brain_mask;
        fprintf('Applied brain mask constraint to seed generation\n');
        
        % OPTIMIZATION 2: Pre-compute dilated brain mask for boundary checking
        fprintf('Pre-computing dilated brain mask...\n');
        nim.dilated_brain_mask = imdilate(brain_mask, ones(3,3,3));
    else
        fprintf('⚠ WARNING: No brain mask found - using FA-only seed mask\n');
    end
end

% Generate seed points
if options.enable_diagnostics
    timing.seed_start = tic;
end

seed_points = generate_seed_points_fact(options.seed_mask, options.seed_density, dims);
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
        % m = memory;
        % fprintf('\nMemory: %.1f GB used', m.MemUsedMATLAB/1e9);
    end
    
    seed = seed_points(i, :);
    
    % Track in both directions
    for direction = [-1, 1]
        if options.enable_diagnostics
            [track, step_timing] = track_fiber_fact(nim, seed, direction, options, cos_angle_thresh);
            timing.interpolation_time = timing.interpolation_time + step_timing.interpolation_time;
            timing.boundary_time = timing.boundary_time + step_timing.boundary_time;
            timing.step_count = timing.step_count + step_timing.step_count;
        else
            track = track_fiber_fact(nim, seed, direction, options, cos_angle_thresh);
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
    
    fprintf('\n\n=== FACT TIMING REPORT ===\n');
    fprintf('Total time: %.2f seconds\n', timing.total_time);
    fprintf('Verification: %.2f seconds (%.1f%%)\n', timing.precompute_time, 100*timing.precompute_time/timing.total_time);
    fprintf('Seed generation: %.2f seconds (%.1f%%)\n', timing.seed_time, 100*timing.seed_time/timing.total_time);
    fprintf('FACT tracking: %.2f seconds (%.1f%%)\n', timing.tracking_time, 100*timing.tracking_time/timing.total_time);
    fprintf('  - Voxel access: 0.00 seconds (no interpolation)\n');
    fprintf('  - Boundary checks: %.2f seconds (%.1f%% of tracking)\n', timing.boundary_time, 100*timing.boundary_time/timing.tracking_time);
    fprintf('Total voxel steps: %d\n', timing.step_count);
    fprintf('Average voxels per track: %.1f\n', timing.step_count / (size(seed_points, 1) * 2));
    fprintf('Voxels per second: %.1f\n', timing.step_count / timing.tracking_time);
    fprintf('========================\n');
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

% Save everything - use v7.3 format for large variables (>2GB support)
fprintf('Saving results (using MAT v7.3 for large file support)...\n');
save(output_file, 'tracks', 'options', 'track_stats', 'track_lengths', 'dims', '-v7.3');
fprintf('Results saved to: %s\n', output_file);

end

function seed_points = generate_seed_points_fact(seed_mask, density, dims)
% FACT: Generate seed points for voxel-center based tracking
% FACT uses only voxel centers, no sub-voxel seeding
[x, y, z] = ind2sub(dims, find(seed_mask));
base_seeds = [x, y, z];
n_base = size(base_seeds, 1);

if density <= 1
    % FACT: Use only voxel centers (no random offsets)
    seed_points = base_seeds;
else
    % FACT: Multiple seeds per voxel at voxel center only
    % Generate multiple instances of the same voxel center
    seed_points = repmat(base_seeds, density, 1);

    % Optional: Add very small perturbation to avoid identical tracks
    % but keep within voxel center for FACT
    if density > 1
        n_total = size(seed_points, 1);
        % Add tiny random offset (±0.01 voxels) to break ties
        perturbation = (rand(n_total, 3) - 0.5) * 0.02;
        seed_points = seed_points + perturbation;
    end
end
end

function [track, step_timing] = track_fiber_fact(nim, seed, direction, options, cos_angle_thresh)
% FACT: Track fiber using discrete voxel tensors (NO interpolation)
%
% FACT TRACKING ALGORITHM:
% - Start from seed voxel center
% - Get direction from current voxel tensor only
% - Step to adjacent voxel in direction of primary eigenvector
% - No sub-voxel positioning or interpolation
%
% Arguments:
%   nim - NIM structure with diffusion tensors
%   seed - Starting position (will be rounded to voxel center)
%   direction - Initial tracking direction (+1 or -1)
%   options - FACT tracking parameters
%   cos_angle_thresh - Cosine of maximum angle change
%
% Returns:
%   track - Array of voxel positions along fiber track

% FACT: Start from voxel center
current_pos = round(seed);
% Pre-allocate track array for performance
track = zeros(options.max_steps + 1, 3);
track(1, :) = current_pos;
track_length = 1;
prev_direction = [];

% Initialize timing
if nargout > 1
    step_timing = struct();
    step_timing.interpolation_time = 0;  % Not used in FACT
    step_timing.boundary_time = 0;
    step_timing.step_count = 0;
end

% Get initial direction from seed voxel
initial_dir = get_initial_direction_fact(nim, current_pos, options);
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
    % FACT: Get direction from current voxel (no interpolation)
    if nargout > 1
        step_timing.step_count = step_timing.step_count + 1;
    end

    [dir_vec, fa_val] = get_voxel_direction_fact(nim, current_pos, options);

    % Termination criteria
    if isempty(dir_vec) || fa_val < options.termination_fa
        break;
    end

    % Ensure consistent direction (flip if needed)
    if ~isempty(prev_direction)
        if dot(dir_vec, prev_direction) < 0
            dir_vec = -dir_vec;
        end

        % Check curvature constraint between voxels
        if dot(dir_vec, prev_direction) < cos_angle_thresh
            break;
        end
    end

    % FACT: Step in direction of primary eigenvector using step size
    % Use proper step size (default 1.0 for FACT) in direction of eigenvector
    next_pos = current_pos + dir_vec * options.step_size;

    % FACT: Round to nearest voxel center for discrete tracking
    next_pos = round(next_pos);

    % Check bounds
    if any(next_pos < 1) || any(next_pos > dims)
        break;
    end

    % Brain tissue check
    if has_parcellation
        if nargout > 1
            boundary_tic = tic;
        end

        if all(next_pos >= 1) && all(next_pos <= dims)
            % Use pre-computed dilated mask
            if ~nim.dilated_brain_mask(next_pos(1), next_pos(2), next_pos(3))
                break;
            end
        end

        if nargout > 1
            step_timing.boundary_time = step_timing.boundary_time + toc(boundary_tic);
        end
    end

    % FACT: Move to next voxel center
    current_pos = next_pos;

    % Avoid adding identical positions (can happen with rounding)
    if ~isequal(current_pos, track(track_length, :))
        track_length = track_length + 1;
        track(track_length, :) = current_pos;
    end

    prev_direction = dir_vec;
end

% Trim track array to actual length
track = track(1:track_length, :);
end

function initial_dir = get_initial_direction_fact(nim, pos, options)
% FACT: Get initial direction from seed voxel (no interpolation)
[initial_dir, fa_val] = get_voxel_direction_fact(nim, pos, options);
if isempty(initial_dir) || fa_val < options.termination_fa
    initial_dir = [];
end
end

function [direction, fa_value] = get_voxel_direction_fact(nim, pos, options)
% FACT: Get direction from nearest voxel tensor (NO interpolation)
%
% FACT ALGORITHM:
% - Round position to nearest voxel center
% - Extract primary eigenvector directly from voxel tensor
% - No smoothing or interpolation between adjacent voxels
%
% Arguments:
%   nim - NIM structure with diffusion tensors
%   pos - Current position (will be rounded to nearest voxel)
%   options - Tractography parameters
%
% Returns:
%   direction - Primary eigenvector from nearest voxel [3x1]
%   fa_value - FA value at voxel center

direction = [];
fa_value = 0;

% FACT: Round to nearest voxel center (no sub-voxel positioning)
voxel_idx = round(pos);
dims = size(nim.FA);

% Check bounds
if any(voxel_idx < 1) || any(voxel_idx > dims)
    return;
end

% FACT: Get FA value directly from voxel center (no interpolation)
fa_value = nim.FA(voxel_idx(1), voxel_idx(2), voxel_idx(3));

if fa_value < options.termination_fa
    return;
end

% FACT: Extract primary eigenvector directly from voxel tensor
try
    % Primary eigenvector (largest eigenvalue) from voxel center
    direction = squeeze(nim.evec(voxel_idx(1), voxel_idx(2), voxel_idx(3), :, 1))';

    % Validate and normalize direction
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