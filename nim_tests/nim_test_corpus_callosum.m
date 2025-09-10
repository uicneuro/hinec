function test_tracks = nim_test_corpus_callosum(nim)
% nim_test_corpus_callosum: Test tracking from corpus callosum
%
% This function tests fiber tracking from the corpus callosum region
% to validate that tracking works in a known high-anisotropy area
%
% Arguments:
%   nim - Structure containing processed DTI data
%
% Returns:
%   test_tracks - Cell array containing test fiber tracks

fprintf('=== CORPUS CALLOSUM TEST ===\n');

% Verify required fields
if ~isfield(nim, 'evec') || ~isfield(nim, 'FA')
    fprintf('ERROR: Required fields (evec, FA) not found\n');
    test_tracks = {};
    return;
end

dims = size(nim.FA);
fprintf('Volume dimensions: %d x %d x %d\n', dims);

% Define corpus callosum region (middle sagittal, upper part)
cc_center = [round(dims(1)/2), round(dims(2)/2), round(dims(3)*0.6)];
cc_size = 5;  % 5x5x5 region

x_range = max(1, cc_center(1)-cc_size):min(dims(1), cc_center(1)+cc_size);
y_range = max(1, cc_center(2)-cc_size):min(dims(2), cc_center(2)+cc_size);
z_range = max(1, cc_center(3)-cc_size):min(dims(3), cc_center(3)+cc_size);

fprintf('Testing corpus callosum region around [%d, %d, %d]\n', cc_center);

% Find highest FA voxel in region
fa_region = nim.FA(x_range, y_range, z_range);
[max_fa, max_idx] = max(fa_region(:));
[rel_x, rel_y, rel_z] = ind2sub(size(fa_region), max_idx);
seed_point = [x_range(rel_x), y_range(rel_y), z_range(rel_z)];

fprintf('Seed point: [%d, %d, %d] with FA = %.4f\n', seed_point, max_fa);

if max_fa < 0.3
    fprintf('WARNING: Low FA in corpus callosum region (%.4f < 0.3)\n', max_fa);
end

% Set conservative tracking parameters for test
options = struct();
options.step_size = 0.5;
options.fa_threshold = 0.25;
options.termination_fa = 0.2;
options.angle_thresh = 45;
options.max_steps = 200;
options.min_length = 15;

% Simple tracking test
test_tracks = {};
track_count = 0;

% Track in both directions
for direction = [-1, 1]
    track = simple_track_from_seed(nim, seed_point, direction, options);
    
    if size(track, 1) > 1
        track_length = sum(vecnorm(diff(track), 2, 2));
        if track_length >= options.min_length
            track_count = track_count + 1;
            test_tracks{track_count} = track;
            fprintf('Direction %d: Track length = %.2f mm (%d points)\n', ...
                direction, track_length, size(track, 1));
        end
    end
end

if track_count == 0
    fprintf('ERROR: No valid tracks generated from corpus callosum!\n');
    fprintf('This indicates serious tracking issues.\n');
else
    fprintf('SUCCESS: Generated %d tracks from corpus callosum\n', track_count);
end

fprintf('============================\n');

end

function track = simple_track_from_seed(nim, seed, direction, options)
% Simple tracking function for testing
track = seed;
current_pos = seed;
prev_direction = [];

dims = size(nim.FA);
cos_angle_thresh = cos(deg2rad(options.angle_thresh));

% Get initial direction
fa_val = nim.FA(round(seed(1)), round(seed(2)), round(seed(3)));
if fa_val < options.termination_fa
    return;
end

initial_dir = squeeze(nim.evec(round(seed(1)), round(seed(2)), round(seed(3)), :, 1))';
initial_dir = initial_dir * direction;
prev_direction = initial_dir;

for step = 1:options.max_steps
    % Get direction at current position
    pos_int = round(current_pos);
    
    % Check bounds
    if any(pos_int < 1) || any(pos_int > dims)
        break;
    end
    
    fa_val = nim.FA(pos_int(1), pos_int(2), pos_int(3));
    
    % Check termination
    if fa_val < options.termination_fa
        break;
    end
    
    % Get eigenvector
    dir_vec = squeeze(nim.evec(pos_int(1), pos_int(2), pos_int(3), :, 1))';
    
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
    
    % Step forward
    current_pos = current_pos + dir_vec * options.step_size;
    track = [track; current_pos];
    prev_direction = dir_vec;
end

end 