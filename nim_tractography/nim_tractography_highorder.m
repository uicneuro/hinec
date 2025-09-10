function tracks = nim_tractography_highorder(data_path, varargin)
% nim_tractography_highorder: High-order deterministic tractography
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

% Set default values
if ~isfield(options, 'seed_density')
    options.seed_density = 1;
end
if ~isfield(options, 'step_size')
    options.step_size = 0.2;  % Smaller step size for better tracking
end
if ~isfield(options, 'fa_threshold')
    options.fa_threshold = 0.1;  % Lower threshold to continue in lower FA areas
end
if ~isfield(options, 'angle_thresh')
    options.angle_thresh = 60;  % More permissive angle threshold
end
if ~isfield(options, 'max_steps')
    options.max_steps = 5000;  % Higher step limit for longer tracks
end
if ~isfield(options, 'min_length')
    options.min_length = 10;  % Minimum track length in mm
end
if ~isfield(options, 'termination_fa')
    options.termination_fa = 0.05;  % Lower termination threshold
end
if ~isfield(options, 'order')
    options.order = 3;
end
if ~isfield(options, 'interp_method')
    options.interp_method = "spline";
end
if ~isfield(options, 'seed_mask')
    options.seed_mask = [];
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

fprintf('Starting high-order tractography...\n');
fprintf('Parameters: order=%d, step=%.2f, FA_thresh=%.2f\n', ...
        options.order, options.step_size, options.fa_threshold);

% Get image dimensions
dims = size(nim.FA);
voxel_size = [1, 1, 1]; % Default voxel size, could be extracted from header

% Create seed mask if not provided
if isempty(options.seed_mask)
    options.seed_mask = nim.FA > options.fa_threshold;
    
    % Priority 1: Use the preprocessed brain mask if available
    if isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
        brain_mask = nim.mask > 0.5;
        options.seed_mask = options.seed_mask & brain_mask;
        fprintf('Applied brain mask from nim.mask (preprocessed)\n');
    elseif isfield(nim, 'parcellation_mask')
        % Fallback: Use parcellation mask if no preprocessed brain mask
        brain_mask = nim.parcellation_mask > 0;
        options.seed_mask = options.seed_mask & brain_mask;
        fprintf('Applied brain mask from parcellation (fallback)\n');
    else
        fprintf('⚠ WARNING: No brain mask found - using FA-only seed mask\n');
    end
end

% Generate seed points
seed_points = generate_seed_points(options.seed_mask, options.seed_density, dims);
fprintf('Generated %d seed points\n', size(seed_points, 1));

% Initialize tracks
tracks = {};
track_count = 0;

% Convert angle threshold to cosine for efficiency
cos_angle_thresh = cos(deg2rad(options.angle_thresh));

% Process each seed point
for i = 1:size(seed_points, 1)
    if mod(i, 1000) == 0
        fprintf('Processing seed %d/%d\n', i, size(seed_points, 1));
    end
    
    seed = seed_points(i, :);
    
    % Track in both directions
    for direction = [-1, 1]
        track = track_fiber(nim, seed, direction, options, cos_angle_thresh, voxel_size);
        
        % Calculate track length in mm
        track_length_mm = 0;
        if size(track, 1) > 1
            for t = 2:size(track, 1)
                track_length_mm = track_length_mm + norm(track(t,:) - track(t-1,:)) * options.step_size;
            end
        end
        
        % Only keep tracks above minimum length
        if track_length_mm >= options.min_length
            track_count = track_count + 1;
            tracks{track_count} = track;
        end
    end
end

fprintf('Generated %d valid tracks\n', length(tracks));
end

function seed_points = generate_seed_points(seed_mask, density, dims)
% Generate seed points based on mask and density
[x, y, z] = ind2sub(dims, find(seed_mask));
seed_points = [x, y, z];

% If density > 1, add sub-voxel seeds
if density > 1
    n_seeds = size(seed_points, 1);
    all_seeds = zeros(n_seeds * density^3, 3);
    
    idx = 1;
    for i = 1:n_seeds
        base_pos = seed_points(i, :);
        for dx = 0:density-1
            for dy = 0:density-1
                for dz = 0:density-1
                    offset = [dx, dy, dz] / density - 0.5 + 0.5/density;
                    all_seeds(idx, :) = base_pos + offset;
                    idx = idx + 1;
                end
            end
        end
    end
    seed_points = all_seeds;
end
end

function track = track_fiber(nim, seed, direction, options, cos_angle_thresh, voxel_size)
% Track a single fiber from seed point using improved high-order FACT
track = seed;
current_pos = seed;
prev_direction = [];

% Get initial direction
initial_dir = get_initial_direction_ho(nim, seed, options);
if isempty(initial_dir)
    return;
end

% Apply direction flip for bidirectional tracking
initial_dir = initial_dir * direction;
prev_direction = initial_dir;

for step = 1:options.max_steps
    % Get interpolated direction at current position
    [dir_vec, fa_val] = interpolate_direction(nim, current_pos, options);
    
    % More lenient termination criteria
    if isempty(dir_vec) || fa_val < options.termination_fa
        break;
    end
    
    % Ensure consistent direction (flip if needed)
    if ~isempty(prev_direction)
        if dot(dir_vec, prev_direction) < 0
            dir_vec = -dir_vec;
        end
        
        % Check curvature constraint
        curvature_check = dot(dir_vec, prev_direction);
        if curvature_check < cos_angle_thresh
            break;
        end
    end
    
    % High-order integration (RK4 for spline, RK2 for linear)
    if strcmp(options.interp_method, "spline") && options.order > 1
        % RK4 integration for high-order
        current_pos = rk4_step(nim, current_pos, dir_vec, options);
    else
        % RK2 integration (midpoint method)
        k1 = dir_vec * options.step_size;
        mid_pos = current_pos + 0.5 * k1;
        
        % Get direction at midpoint
        [mid_dir, mid_fa] = interpolate_direction(nim, mid_pos, options);
        if ~isempty(mid_dir) && mid_fa >= options.termination_fa
            % Ensure direction consistency at midpoint
            if dot(mid_dir, dir_vec) < 0
                mid_dir = -mid_dir;
            end
            k2 = mid_dir * options.step_size;
            current_pos = current_pos + k2;
        else
            current_pos = current_pos + k1;
        end
    end
    
    % Check bounds with buffer
    dims = size(nim.FA);
    if any(current_pos < 1.5) || any(current_pos > dims - 0.5)
        break;
    end
    
    % Brain tissue boundary check - use proper brain mask
    pos_int = round(current_pos);
    if all(pos_int >= 1) && all(pos_int <= dims)
        % Priority 1: Use preprocessed brain mask
        if isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
            if nim.mask(pos_int(1), pos_int(2), pos_int(3)) <= 0.5
                break;
            end
        elseif isfield(nim, 'parcellation_mask')
            % Fallback: Use parcellation mask with boundary check
            if nim.parcellation_mask(pos_int(1), pos_int(2), pos_int(3)) <= 0
                % Check boundary neighbors
                boundary_ok = false;
                for dx = -1:1
                    for dy = -1:1
                        for dz = -1:1
                            check_pos = pos_int + [dx, dy, dz];
                            if all(check_pos >= 1) && all(check_pos <= size(nim.parcellation_mask))
                                if nim.parcellation_mask(check_pos(1), check_pos(2), check_pos(3)) > 0
                                    boundary_ok = true;
                                    break;
                                end
                            end
                        end
                        if boundary_ok, break; end
                    end
                    if boundary_ok, break; end
                end
                if ~boundary_ok
                    break;
                end
            end
        end
    end
    
    track = [track; current_pos];
    prev_direction = dir_vec;
end
end

function initial_dir = get_initial_direction_ho(nim, pos, options)
% Get initial direction at seed point for high-order method
[initial_dir, fa_val] = interpolate_direction(nim, pos, options);
if isempty(initial_dir) || fa_val < options.termination_fa
    initial_dir = [];
end
end

function new_pos = rk4_step(nim, pos, dir_vec, options)
% RK4 integration step for high-order tracking
h = options.step_size;

k1 = h * dir_vec;
[k2_dir, ~] = interpolate_direction(nim, pos + 0.5*k1, options);
if ~isempty(k2_dir)
    if dot(k2_dir, dir_vec) < 0, k2_dir = -k2_dir; end
    k2 = h * k2_dir;
else
    k2 = k1;
end

[k3_dir, ~] = interpolate_direction(nim, pos + 0.5*k2, options);
if ~isempty(k3_dir)
    if dot(k3_dir, dir_vec) < 0, k3_dir = -k3_dir; end
    k3 = h * k3_dir;
else
    k3 = k2;
end

[k4_dir, ~] = interpolate_direction(nim, pos + k3, options);
if ~isempty(k4_dir)
    if dot(k4_dir, dir_vec) < 0, k4_dir = -k4_dir; end
    k4 = h * k4_dir;
else
    k4 = k3;
end

new_pos = pos + (k1 + 2*k2 + 2*k3 + k4) / 6;
end

function [direction, fa_value] = interpolate_direction(nim, pos, options)
% Improved interpolation with better termination criteria
direction = [];
fa_value = 0;

% Check if position is within bounds with buffer
dims = size(nim.FA);
if any(pos < 1.1) || any(pos > dims - 0.1)
    return;
end

% Get FA value at position with bounds checking
try
    fa_value = interp3(nim.FA, pos(2), pos(1), pos(3), 'linear', 0);
catch
    return;
end

% More lenient FA check
if fa_value < options.termination_fa
    return;
end

% Get primary eigenvector (V1) at position
if strcmp(options.interp_method, "spline") && options.order > 1
    % High-order spline interpolation
    direction = interpolate_eigenvector_spline(nim, pos, options.order);
else
    % Linear interpolation
    direction = interpolate_eigenvector_linear(nim, pos);
end

% Normalize direction and check validity
if ~isempty(direction)
    dir_norm = norm(direction);
    if dir_norm > 1e-6  % More lenient normalization check
        direction = direction / dir_norm;
    else
        direction = [];
    end
end
end

function dir_vec = interpolate_eigenvector_linear(nim, pos)
% Improved linear interpolation of eigenvector
try
    % Check bounds more carefully
    dims = size(nim.evec);
    if any(pos < 1.01) || pos(1) > dims(1)-0.01 || pos(2) > dims(2)-0.01 || pos(3) > dims(3)-0.01
        dir_vec = [];
        return;
    end
    
    % Extract primary eigenvector components
    v1_x = squeeze(nim.evec(:,:,:,1,1));
    v1_y = squeeze(nim.evec(:,:,:,2,1));
    v1_z = squeeze(nim.evec(:,:,:,3,1));
    
    % Robust trilinear interpolation
    x_comp = interp3(v1_x, pos(2), pos(1), pos(3), 'linear', 0);
    y_comp = interp3(v1_y, pos(2), pos(1), pos(3), 'linear', 0);
    z_comp = interp3(v1_z, pos(2), pos(1), pos(3), 'linear', 0);
    
    dir_vec = [x_comp, y_comp, z_comp];
    
    % Improved validity check
    dir_norm = norm(dir_vec);
    if dir_norm < 1e-6 || any(isnan(dir_vec)) || any(isinf(dir_vec))
        dir_vec = [];
    end
    
catch ME
    fprintf('Warning: Linear eigenvector interpolation failed at position [%.2f, %.2f, %.2f]: %s\n', ...
            pos(1), pos(2), pos(3), ME.message);
    dir_vec = [];
end
end

function dir_vec = interpolate_eigenvector_spline(nim, pos, order)
% Improved high-order spline interpolation of eigenvector
try
    % Check bounds first
    dims = size(nim.evec);
    radius = max(2, ceil(order/2));
    
    if pos(1) < radius || pos(1) > dims(1)-radius || ...
       pos(2) < radius || pos(2) > dims(2)-radius || ...
       pos(3) < radius || pos(3) > dims(3)-radius
        % Near boundary, use linear interpolation
        dir_vec = interpolate_eigenvector_linear(nim, pos);
        return;
    end
    
    % Get neighborhood for spline interpolation
    x_range = max(1, floor(pos(1)) - radius) : min(dims(1), ceil(pos(1)) + radius);
    y_range = max(1, floor(pos(2)) - radius) : min(dims(2), ceil(pos(2)) + radius);
    z_range = max(1, floor(pos(3)) - radius) : min(dims(3), ceil(pos(3)) + radius);
    
    % Ensure we have enough points for spline
    if length(x_range) < 4 || length(y_range) < 4 || length(z_range) < 4
        dir_vec = interpolate_eigenvector_linear(nim, pos);
        return;
    end
    
    % Extract primary eigenvector components in neighborhood
    v1_local = nim.evec(x_range, y_range, z_range, :, 1);
    
    % Create coordinate grids
    [Y, X, Z] = meshgrid(y_range, x_range, z_range);
    
    % Spline interpolation with error handling
    x_comp = interp3(Y, X, Z, squeeze(v1_local(:,:,:,1)), pos(2), pos(1), pos(3), 'spline', 0);
    y_comp = interp3(Y, X, Z, squeeze(v1_local(:,:,:,2)), pos(2), pos(1), pos(3), 'spline', 0);
    z_comp = interp3(Y, X, Z, squeeze(v1_local(:,:,:,3)), pos(2), pos(1), pos(3), 'spline', 0);
    
    dir_vec = [x_comp, y_comp, z_comp];
    
    % Check validity
    dir_norm = norm(dir_vec);
    if dir_norm < 1e-6 || any(isnan(dir_vec)) || any(isinf(dir_vec))
        % Fallback to linear if spline fails
        dir_vec = interpolate_eigenvector_linear(nim, pos);
    end
    
catch ME
    fprintf('Warning: Spline interpolation failed, using linear fallback: %s\n', ME.message);
    % Fallback to linear interpolation
    dir_vec = interpolate_eigenvector_linear(nim, pos);
end
end 