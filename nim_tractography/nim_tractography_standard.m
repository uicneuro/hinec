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
for i = 1:size(seed_points, 1)
    if mod(i, 1000) == 0
        fprintf('Processing seed %d/%d\n', i, size(seed_points, 1));
    end
    
    seed = seed_points(i, :);
    
    % Track in both directions
    for direction = [-1, 1]
        track = track_fiber_standard(nim, seed, direction, options, cos_angle_thresh);
        
        if size(track, 1) > 2 % Minimum track length
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
% Track a single fiber from seed point using standard linear interpolation
track = seed;
current_pos = seed;
prev_direction = [];

for step = 1:options.max_steps
    % Get interpolated direction at current position (linear only)
    [dir_vec, fa_val] = interpolate_direction_standard(nim, current_pos, options);
    
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
    
    % Take step using Euler integration
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

function [direction, fa_value] = interpolate_direction_standard(nim, pos, options)
% Standard linear interpolation of direction vector
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

% Get primary eigenvector using linear interpolation only
direction = interpolate_eigenvector_linear_standard(nim, pos);

% Normalize direction
if ~isempty(direction)
    direction = direction / norm(direction);
end
end

function dir_vec = interpolate_eigenvector_linear_standard(nim, pos)
% Standard linear interpolation of eigenvector (simplified)
try
    % Extract primary eigenvector (first column of evec)
    v1_x = squeeze(nim.evec(:,:,:,1,1));
    v1_y = squeeze(nim.evec(:,:,:,2,1));
    v1_z = squeeze(nim.evec(:,:,:,3,1));
    
    % Simple trilinear interpolation
    x_comp = interp3(v1_x, pos(2), pos(1), pos(3), 'linear', 0);
    y_comp = interp3(v1_y, pos(2), pos(1), pos(3), 'linear', 0);
    z_comp = interp3(v1_z, pos(2), pos(1), pos(3), 'linear', 0);
    
    dir_vec = [x_comp, y_comp, z_comp];
    
    % Basic consistency check
    if norm(dir_vec) < 0.1
        dir_vec = [];
    end
catch
    dir_vec = [];
end
end 