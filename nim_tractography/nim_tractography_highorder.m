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
    options.step_size = 0.5;
end
if ~isfield(options, 'fa_threshold')
    options.fa_threshold = 0.2;
end
if ~isfield(options, 'angle_thresh')
    options.angle_thresh = 45;
end
if ~isfield(options, 'max_steps')
    options.max_steps = 2000;
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
    
    % Exclude non-brain areas (parcellation label 0) if parcellation is available
    if isfield(nim, 'parcellation_mask')
        brain_mask = nim.parcellation_mask > 0;
        options.seed_mask = options.seed_mask & brain_mask;
        fprintf('Applied brain mask from parcellation (excluding label 0)\n');
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
        
        if size(track, 1) > 2 % Minimum track length
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
% Track a single fiber from seed point
track = seed;
current_pos = seed;
prev_direction = [];

for step = 1:options.max_steps
    % Get interpolated direction at current position
    [dir_vec, fa_val] = interpolate_direction(nim, current_pos, options);
    
    if isempty(dir_vec) || fa_val < options.fa_threshold
        break;
    end
    
    % Ensure consistent direction
    if ~isempty(prev_direction)
        if dot(dir_vec, prev_direction) < 0
            dir_vec = -dir_vec;
        end
        
        % Check angle constraint
        if dot(dir_vec, prev_direction) < cos_angle_thresh
            break;
        end
    end
    
    % Apply initial direction
    dir_vec = dir_vec * direction;
    direction = 1; % Only apply initial direction once
    
    % Take step
    current_pos = current_pos + dir_vec * options.step_size;
    
    % Check bounds
    if any(current_pos < 1) || any(current_pos > size(nim.FA))
        break;
    end
    
    % Check if we're still in brain tissue (parcellation > 0)
    if isfield(nim, 'parcellation_mask')
        pos_int = round(current_pos);
        if all(pos_int >= 1) && all(pos_int <= size(nim.parcellation_mask))
            if nim.parcellation_mask(pos_int(1), pos_int(2), pos_int(3)) == 0
                break; % Stop tracking in non-brain areas
            end
        end
    end
    
    track = [track; current_pos];
    prev_direction = dir_vec;
end
end

function [direction, fa_value] = interpolate_direction(nim, pos, options)
% Interpolate direction vector at given position using high-order interpolation
direction = [];
fa_value = 0;

% Check if position is within bounds
dims = size(nim.FA);
if any(pos < 1) || any(pos > dims)
    return;
end

% Get FA value at position
fa_value = interp3(nim.FA, pos(2), pos(1), pos(3), 'linear', 0);

if fa_value < options.fa_threshold
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

% Normalize direction
if ~isempty(direction)
    direction = direction / norm(direction);
end
end

function dir_vec = interpolate_eigenvector_linear(nim, pos)
% Linear interpolation of eigenvector
try
    % Extract primary eigenvector (first column of evec)
    v1_x = squeeze(nim.evec(:,:,:,1,1));
    v1_y = squeeze(nim.evec(:,:,:,2,1));
    v1_z = squeeze(nim.evec(:,:,:,3,1));
    
    % Interpolate each component
    x_comp = interp3(v1_x, pos(2), pos(1), pos(3), 'linear', 0);
    y_comp = interp3(v1_y, pos(2), pos(1), pos(3), 'linear', 0);
    z_comp = interp3(v1_z, pos(2), pos(1), pos(3), 'linear', 0);
    
    dir_vec = [x_comp, y_comp, z_comp];
catch
    dir_vec = [];
end
end

function dir_vec = interpolate_eigenvector_spline(nim, pos, order)
% High-order spline interpolation of eigenvector
try
    % Get neighborhood for interpolation
    radius = ceil(order/2);
    x_range = max(1, floor(pos(1)) - radius) : min(size(nim.evec,1), ceil(pos(1)) + radius);
    y_range = max(1, floor(pos(2)) - radius) : min(size(nim.evec,2), ceil(pos(2)) + radius);
    z_range = max(1, floor(pos(3)) - radius) : min(size(nim.evec,3), ceil(pos(3)) + radius);
    
    % Extract primary eigenvector components in neighborhood
    v1_local = nim.evec(x_range, y_range, z_range, :, 1);
    
    % Use spline interpolation
    [X, Y, Z] = meshgrid(y_range, x_range, z_range);
    
    x_comp = interp3(X, Y, Z, squeeze(v1_local(:,:,:,1)), pos(2), pos(1), pos(3), 'spline', 0);
    y_comp = interp3(X, Y, Z, squeeze(v1_local(:,:,:,2)), pos(2), pos(1), pos(3), 'spline', 0);
    z_comp = interp3(X, Y, Z, squeeze(v1_local(:,:,:,3)), pos(2), pos(1), pos(3), 'spline', 0);
    
    dir_vec = [x_comp, y_comp, z_comp];
catch
    % Fallback to linear interpolation
    dir_vec = interpolate_eigenvector_linear(nim, pos);
end
end 