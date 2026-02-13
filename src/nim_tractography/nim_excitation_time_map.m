function [time_map, velocity_field] = nim_excitation_time_map(nim, seed_points, varargin)
% nim_excitation_time_map: Compute excitation time map based on DTI
%
% Arguments:
%   nim - NIM structure with FA and eigenvectors
%   seed_points - Nx3 matrix of seed coordinates
%   options - Options for time map computation (optional struct)
%
% Returns:
%   time_map - 3D map of excitation times
%   velocity_field - 3D velocity field

% Parse input arguments
if nargin > 2 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default values
if ~isfield(options, 'conduction_speed')
    options.conduction_speed = 1.0;
end
if ~isfield(options, 'fa_scaling')
    options.fa_scaling = true;
end
if ~isfield(options, 'max_time')
    options.max_time = 100;
end
if ~isfield(options, 'method')
    options.method = "fast_marching";
end

fprintf('Computing excitation time map...\n');
fprintf('Seed points: %d\n', size(seed_points, 1));
fprintf('Conduction speed: %.2f\n', options.conduction_speed);

% Get image dimensions
dims = size(nim.FA);

% Initialize time map
time_map = inf(dims);
velocity_field = zeros([dims, 3]);

% Compute velocity field based on FA and eigenvectors
velocity_field = compute_velocity_field(nim, options);

% Set seed points to time zero
for i = 1:size(seed_points, 1)
    seed = round(seed_points(i, :));
    if all(seed >= 1) && all(seed <= dims)
        time_map(seed(1), seed(2), seed(3)) = 0;
    end
end

% Compute time map using chosen method
switch options.method
    case "fast_marching"
        time_map = fast_marching_time_map(time_map, velocity_field, options);
    case "dijkstra"
        time_map = dijkstra_time_map(time_map, velocity_field, options);
    otherwise
        error('Unknown method: %s', options.method);
end

% Cap maximum time
time_map(time_map > options.max_time) = options.max_time;

fprintf('Excitation time map completed\n');
end

function velocity_field = compute_velocity_field(nim, options)
% Compute velocity field based on FA and eigenvectors
dims = size(nim.FA);
velocity_field = zeros([dims, 3]);

if ~isfield(nim, 'evec')
    warning('Eigenvectors not found. Using isotropic velocity.');
    velocity_field = repmat(options.conduction_speed, [dims, 3]);
    return;
end

for x = 1:dims(1)
    for y = 1:dims(2)
        for z = 1:dims(3)
            fa_val = nim.FA(x, y, z);
            
            if fa_val > 0.1 % Only process voxels with meaningful FA
                % Get primary eigenvector
                v1 = squeeze(nim.evec(x, y, z, :, 1));
                
                % Scale velocity by FA if requested
                if options.fa_scaling
                    speed = options.conduction_speed * (0.2 + 0.8 * fa_val);
                else
                    speed = options.conduction_speed;
                end
                
                % Velocity in direction of primary eigenvector
                velocity_field(x, y, z, :) = speed * v1;
            else
                % Low FA: isotropic slow conduction
                velocity_field(x, y, z, :) = options.conduction_speed * 0.1;
            end
        end
    end
end
end

function time_map = fast_marching_time_map(time_map, velocity_field, options)
% Fast marching method for time map computation
dims = size(time_map);

% Priority queue simulation using simple sorted list
active_points = [];

% Find initial points (seeds)
[seed_x, seed_y, seed_z] = ind2sub(dims, find(time_map == 0));
for i = 1:length(seed_x)
    active_points = [active_points; seed_x(i), seed_y(i), seed_z(i), 0];
end

% Neighbor offsets (6-connectivity)
neighbors = [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1];

iteration = 0;
while ~isempty(active_points) && iteration < 1e6
    iteration = iteration + 1;
    
    if mod(iteration, 10000) == 0
        fprintf('Fast marching iteration %d, active points: %d\n', iteration, size(active_points, 1));
    end
    
    % Sort by time and get point with minimum time
    [~, min_idx] = min(active_points(:, 4));
    current_point = active_points(min_idx, :);
    active_points(min_idx, :) = [];
    
    x = current_point(1);
    y = current_point(2);
    z = current_point(3);
    current_time = current_point(4);
    
    % Process neighbors
    for n = 1:size(neighbors, 1)
        nx = x + neighbors(n, 1);
        ny = y + neighbors(n, 2);
        nz = z + neighbors(n, 3);
        
        % Check bounds
        if nx < 1 || nx > dims(1) || ny < 1 || ny > dims(2) || nz < 1 || nz > dims(3)
            continue;
        end
        
        % Calculate travel time to neighbor
        vel = squeeze(velocity_field(x, y, z, :));
        if norm(vel) > 0
            step_time = 1 / norm(vel);
        else
            step_time = 10; % Very slow for zero velocity
        end
        
        new_time = current_time + step_time;
        
        % Update if faster path found
        if new_time < time_map(nx, ny, nz)
            time_map(nx, ny, nz) = new_time;
            
            % Add to active points if not already there
            existing = find(active_points(:, 1) == nx & active_points(:, 2) == ny & active_points(:, 3) == nz);
            if isempty(existing)
                active_points = [active_points; nx, ny, nz, new_time];
            else
                active_points(existing, 4) = new_time;
            end
        end
    end
end

fprintf('Fast marching completed after %d iterations\n', iteration);
end

function time_map = dijkstra_time_map(time_map, velocity_field, options)
% Dijkstra's algorithm for time map computation
dims = size(time_map);

% Convert to linear indices for easier handling
visited = false(dims);
[seed_x, seed_y, seed_z] = ind2sub(dims, find(time_map == 0));

% Simple priority queue using sorted array
queue = [];
for i = 1:length(seed_x)
    idx = sub2ind(dims, seed_x(i), seed_y(i), seed_z(i));
    queue = [queue; 0, idx];
end

% Neighbor offsets
neighbors = [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1];

iteration = 0;
while ~isempty(queue) && iteration < 1e6
    iteration = iteration + 1;
    
    if mod(iteration, 10000) == 0
        fprintf('Dijkstra iteration %d, queue size: %d\n', iteration, size(queue, 1));
    end
    
    % Get minimum time point
    [~, min_idx] = min(queue(:, 1));
    current_time = queue(min_idx, 1);
    current_idx = queue(min_idx, 2);
    queue(min_idx, :) = [];
    
    % Skip if already visited
    if visited(current_idx)
        continue;
    end
    visited(current_idx) = true;
    
    [x, y, z] = ind2sub(dims, current_idx);
    
    % Process neighbors
    for n = 1:size(neighbors, 1)
        nx = x + neighbors(n, 1);
        ny = y + neighbors(n, 2);
        nz = z + neighbors(n, 3);
        
        % Check bounds
        if nx < 1 || nx > dims(1) || ny < 1 || ny > dims(2) || nz < 1 || nz > dims(3)
            continue;
        end
        
        neighbor_idx = sub2ind(dims, nx, ny, nz);
        
        % Skip if already visited
        if visited(neighbor_idx)
            continue;
        end
        
        % Calculate travel time
        vel = squeeze(velocity_field(x, y, z, :));
        if norm(vel) > 0
            step_time = 1 / norm(vel);
        else
            step_time = 10;
        end
        
        new_time = current_time + step_time;
        
        % Update time map and add to queue
        if new_time < time_map(nx, ny, nz)
            time_map(nx, ny, nz) = new_time;
            queue = [queue; new_time, neighbor_idx];
        end
    end
end

fprintf('Dijkstra completed after %d iterations\n', iteration);
end 