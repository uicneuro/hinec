function tracks = nim_tractography_standard(data_path, varargin)
% nim_tractography_standard: Standard deterministic tractography
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
    options.order = 1;
end
if ~isfield(options, 'interp_method')
    options.interp_method = "linear";
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

fprintf('Starting standard tractography...\n');
fprintf('Parameters: step=%.2f, FA_thresh=%.2f, angle_thresh=%.1f\n', ...
        options.step_size, options.fa_threshold, options.angle_thresh);

% Get image dimensions
dims = size(nim.FA);

% Create seed mask if not provided
if isempty(options.seed_mask)
    options.seed_mask = nim.FA > options.fa_threshold;
    
    % Exclude non-brain areas (parcellation label 0) if parcellation is available
    if isfield(nim, 'parcellation_mask')
        brain_mask = nim.parcellation_mask > 0;
        options.seed_mask = options.seed_mask & brain_mask;
        fprintf('Applied brain mask from parcellation (excluding label 0)\n');
    end
end

% Generate seed points
seed_points = generate_seed_points_standard(options.seed_mask, options.seed_density, dims);
fprintf('Generated %d seed points\n', size(seed_points, 1));

% Initialize tracks
tracks = {};
track_count = 0;

% Convert angle threshold to cosine for efficiency
cos_angle_thresh = cos(deg2rad(options.angle_thresh));

% Process each seed point
fprintf('Processing seeds: ');
for i = 1:size(seed_points, 1)
    if mod(i, 500) == 0
        fprintf('%d/%d ', i, size(seed_points, 1));
        if mod(i, 2500) == 0
            fprintf('\n               ');
        end
    end
    
    seed = seed_points(i, :);
    
    % Track in both directions
    for direction = [-1, 1]
        track = track_fiber_standard(nim, seed, direction, options, cos_angle_thresh);
        
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

function seed_points = generate_seed_points_standard(seed_mask, density, dims)
% Generate seed points based on mask and density (simplified version)
[x, y, z] = ind2sub(dims, find(seed_mask));
seed_points = [x, y, z];

% For standard tractography, use simpler seeding
if density > 1
    % Add random sub-voxel offsets
    n_seeds = size(seed_points, 1);
    all_seeds = [];
    for i = 1:n_seeds
        for j = 1:density
            offset = (rand(1, 3) - 0.5) * 0.8; % Random offset within voxel
            all_seeds = [all_seeds; seed_points(i, :) + offset];
        end
    end
    seed_points = all_seeds;
end
end

function track = track_fiber_standard(nim, seed, direction, options, cos_angle_thresh)
% Track a single fiber from seed point using improved FACT algorithm
track = seed;
current_pos = seed;
prev_direction = [];

% Get initial direction
initial_dir = get_initial_direction(nim, seed, options);
if isempty(initial_dir)
    return;
end

% Apply direction flip for bidirectional tracking
initial_dir = initial_dir * direction;
prev_direction = initial_dir;

for step = 1:options.max_steps
    % Get interpolated direction at current position
    [dir_vec, fa_val] = interpolate_direction_standard(nim, current_pos, options);
    
    % More lenient termination criteria
    if isempty(dir_vec) || fa_val < options.termination_fa
        break;
    end
    
    % Ensure consistent direction (flip if needed)
    if ~isempty(prev_direction)
        if dot(dir_vec, prev_direction) < 0
            dir_vec = -dir_vec;
        end
        
        % Check curvature constraint - use dot product directly
        curvature_check = dot(dir_vec, prev_direction);
        if curvature_check < cos_angle_thresh
            break;
        end
    end
    
    % Improved step integration with RK2 (midpoint method)
    k1 = dir_vec * options.step_size;
    mid_pos = current_pos + 0.5 * k1;
    
    % Get direction at midpoint
    [mid_dir, mid_fa] = interpolate_direction_standard(nim, mid_pos, options);
    if ~isempty(mid_dir) && mid_fa >= options.termination_fa
        % Ensure direction consistency at midpoint
        if dot(mid_dir, dir_vec) < 0
            mid_dir = -mid_dir;
        end
        k2 = mid_dir * options.step_size;
        current_pos = current_pos + k2;  % Use midpoint slope
    else
        % Fall back to Euler if midpoint fails
        current_pos = current_pos + k1;
    end
    
    % Check bounds with small buffer
    dims = size(nim.FA);
    if any(current_pos < 1.5) || any(current_pos > dims - 0.5)
        break;
    end
    
    % More lenient brain tissue check
    if isfield(nim, 'parcellation_mask')
        pos_int = round(current_pos);
        if all(pos_int >= 1) && all(pos_int <= size(nim.parcellation_mask))
            % Allow tracking in very low parcellation values (near boundaries)
            if nim.parcellation_mask(pos_int(1), pos_int(2), pos_int(3)) <= 0
                % Check if we're just at boundary - look at neighbors
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

function initial_dir = get_initial_direction(nim, pos, options)
% Get initial direction at seed point
[initial_dir, fa_val] = interpolate_direction_standard(nim, pos, options);
if isempty(initial_dir) || fa_val < options.termination_fa
    initial_dir = [];
end
end

function [direction, fa_value] = interpolate_direction_standard(nim, pos, options)
% Improved linear interpolation of direction vector
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

% More lenient FA check - don't return if above termination threshold
if fa_value < options.termination_fa
    return;
end

% Get primary eigenvector using improved linear interpolation
direction = interpolate_eigenvector_linear_standard(nim, pos);

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

function dir_vec = interpolate_eigenvector_linear_standard(nim, pos)
% Improved linear interpolation of eigenvector with better handling
try
    % Check bounds more carefully
    dims = size(nim.evec);
    if any(pos < 1.01) || pos(1) > dims(1)-0.01 || pos(2) > dims(2)-0.01 || pos(3) > dims(3)-0.01
        dir_vec = [];
        return;
    end
    
    % Extract primary eigenvector components (first column of evec)
    v1_x = squeeze(nim.evec(:,:,:,1,1));
    v1_y = squeeze(nim.evec(:,:,:,2,1));
    v1_z = squeeze(nim.evec(:,:,:,3,1));
    
    % Robust trilinear interpolation with error handling
    x_comp = interp3(v1_x, pos(2), pos(1), pos(3), 'linear', 0);
    y_comp = interp3(v1_y, pos(2), pos(1), pos(3), 'linear', 0);
    z_comp = interp3(v1_z, pos(2), pos(1), pos(3), 'linear', 0);
    
    dir_vec = [x_comp, y_comp, z_comp];
    
    % Improved consistency check
    dir_norm = norm(dir_vec);
    if dir_norm < 1e-6 || any(isnan(dir_vec)) || any(isinf(dir_vec))
        dir_vec = [];
    end
    
catch ME
    fprintf('Warning: Eigenvector interpolation failed at position [%.2f, %.2f, %.2f]: %s\n', ...
            pos(1), pos(2), pos(3), ME.message);
    dir_vec = [];
end
end 