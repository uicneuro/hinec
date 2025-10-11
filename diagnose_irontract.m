% diagnose_irontract.m - Diagnostic script for IronTract submission issues
%
% This script helps identify why binary mask contains 0 voxels

fprintf('=== IronTract Diagnostic Tool ===\n\n');

%% Load data
fprintf('1. Loading data files...\n');
data_path = 'ironTract.mat';
injection_file = 'ironTract/sub-MR243_injection.nii.gz';
tracks_file = 'tractography_results/tracks_standard.mat';

if ~exist(data_path, 'file')
    error('Data file not found: %s', data_path);
end
if ~exist(injection_file, 'file')
    error('Injection file not found: %s', injection_file);
end
if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s', tracks_file);
end

load(data_path, 'nim');
load(tracks_file, 'tracks');
injection_mask = niftiread(injection_file);
injection_info = niftiinfo(injection_file);

fprintf('  ✓ Loaded nim structure\n');
fprintf('  ✓ Loaded %d tracks\n', length(tracks));
fprintf('  ✓ Loaded injection mask\n\n');

%% Check dimensions
fprintf('2. Checking dimensions...\n');
nim_dims = [nim.xdim, nim.ydim, nim.zdim];
inj_dims = size(injection_mask);

fprintf('  NIM dimensions: %d x %d x %d\n', nim_dims);
fprintf('  Injection dimensions: %d x %d x %d\n', inj_dims);

if isequal(nim_dims, inj_dims)
    fprintf('  ✓ Dimensions match\n\n');
else
    fprintf('  ❌ DIMENSION MISMATCH!\n\n');
end

%% Analyze injection mask
fprintf('3. Analyzing injection site...\n');
injection_voxels = sum(injection_mask(:) > 0);
injection_percentage = 100 * injection_voxels / prod(inj_dims);
fprintf('  Injection site voxels: %d (%.3f%% of volume)\n', injection_voxels, injection_percentage);

if injection_voxels == 0
    fprintf('  ❌ CRITICAL: Injection mask is empty!\n\n');
    return;
end

% Find injection site extent
[inj_x, inj_y, inj_z] = ind2sub(inj_dims, find(injection_mask > 0));
inj_range_x = [min(inj_x), max(inj_x)];
inj_range_y = [min(inj_y), max(inj_y)];
inj_range_z = [min(inj_z), max(inj_z)];

fprintf('  Injection X range: [%d, %d]\n', inj_range_x);
fprintf('  Injection Y range: [%d, %d]\n', inj_range_y);
fprintf('  Injection Z range: [%d, %d]\n', inj_range_z);
fprintf('  Injection center: [%.1f, %.1f, %.1f]\n', mean(inj_x), mean(inj_y), mean(inj_z));
fprintf('\n');

%% Analyze track coordinates
fprintf('4. Analyzing track coordinate system...\n');
if isempty(tracks)
    fprintf('  ❌ No tracks to analyze!\n\n');
    return;
end

% Sample first 10 tracks
num_sample = min(10, length(tracks));
fprintf('  Sampling first %d tracks...\n', num_sample);

all_coords = [];
for i = 1:num_sample
    all_coords = [all_coords; tracks{i}]; %#ok<AGROW>
end

track_range_x = [min(all_coords(:,1)), max(all_coords(:,1))];
track_range_y = [min(all_coords(:,2)), max(all_coords(:,2))];
track_range_z = [min(all_coords(:,3)), max(all_coords(:,3))];

fprintf('  Track X range: [%.2f, %.2f]\n', track_range_x);
fprintf('  Track Y range: [%.2f, %.2f]\n', track_range_y);
fprintf('  Track Z range: [%.2f, %.2f]\n', track_range_z);
fprintf('  Track center: [%.1f, %.1f, %.1f]\n', mean(all_coords(:,1)), mean(all_coords(:,2)), mean(all_coords(:,3)));
fprintf('\n');

%% Check coordinate system alignment
fprintf('5. Checking coordinate system alignment...\n');

% Check if coordinates are in voxel indices (1-based) or world coordinates
if track_range_x(1) < 1 || track_range_y(1) < 1 || track_range_z(1) < 1
    fprintf('  ⚠️  WARNING: Track coordinates include values < 1\n');
    fprintf('      Tracks may be in world coordinates, not voxel indices!\n');
end

if track_range_x(2) > nim_dims(1) || track_range_y(2) > nim_dims(2) || track_range_z(2) > nim_dims(3)
    fprintf('  ⚠️  WARNING: Track coordinates exceed image dimensions\n');
    fprintf('      Expected range: [1, %d] x [1, %d] x [1, %d]\n', nim_dims);
end

% Check overlap between track and injection ranges
x_overlap = (track_range_x(2) >= inj_range_x(1)) && (track_range_x(1) <= inj_range_x(2));
y_overlap = (track_range_y(2) >= inj_range_y(1)) && (track_range_y(1) <= inj_range_y(2));
z_overlap = (track_range_z(2) >= inj_range_z(1)) && (track_range_z(1) <= inj_range_z(2));

fprintf('  X-axis overlap: %s\n', iif(x_overlap, 'YES', 'NO'));
fprintf('  Y-axis overlap: %s\n', iif(y_overlap, 'YES', 'NO'));
fprintf('  Z-axis overlap: %s\n', iif(z_overlap, 'YES', 'NO'));

if x_overlap && y_overlap && z_overlap
    fprintf('  ✓ Track and injection ranges overlap in all dimensions\n\n');
else
    fprintf('  ❌ CRITICAL: Track and injection ranges DO NOT overlap!\n');
    fprintf('      This explains why no tracks pass through injection site.\n\n');
end

%% Test filtering with current implementation
fprintf('6. Testing track filtering...\n');
filtered_count = 0;
for i = 1:length(tracks)
    track = tracks{i};
    track_voxels = round(track);

    intersects = false;
    for j = 1:size(track_voxels, 1)
        pos = track_voxels(j, :);

        if all(pos >= 1) && all(pos <= inj_dims)
            if injection_mask(pos(1), pos(2), pos(3)) > 0
                intersects = true;
                break;
            end
        end
    end

    if intersects
        filtered_count = filtered_count + 1;
    end
end

fprintf('  Tracks passing through injection site: %d / %d (%.2f%%)\n', ...
    filtered_count, length(tracks), 100*filtered_count/length(tracks));

if filtered_count == 0
    fprintf('  ❌ CONFIRMED: No tracks intersect injection site\n\n');
else
    fprintf('  ✓ Some tracks intersect injection site\n\n');
end

%% Check NIfTI header info
fprintf('7. Checking NIfTI header information...\n');
fprintf('  Injection file header:\n');
fprintf('    Datatype: %s\n', injection_info.Datatype);
fprintf('    Transform:\n');
disp(injection_info.Transform.T);
fprintf('\n');

if isfield(nim, 'hdr')
    fprintf('  NIM file header:\n');
    fprintf('    Transform:\n');
    disp(nim.hdr.Transform.T);
end

%% Recommendations
fprintf('\n=== DIAGNOSIS SUMMARY ===\n');
if filtered_count == 0
    fprintf('❌ PROBLEM IDENTIFIED: No tracks pass through injection site\n\n');
    fprintf('Possible causes:\n');
    if ~(x_overlap && y_overlap && z_overlap)
        fprintf('  1. COORDINATE SYSTEM MISMATCH (most likely)\n');
        fprintf('     - Tracks and injection mask are in different coordinate spaces\n');
        fprintf('     - May need to transform coordinates or check voxel-to-world mapping\n\n');
    end

    fprintf('  2. INSUFFICIENT TRACK COVERAGE\n');
    fprintf('     - Injection site may be in region with low FA\n');
    fprintf('     - Try lowering FA thresholds or increasing seed density\n\n');

    fprintf('  3. DATA PROCESSING ISSUE\n');
    fprintf('     - Check if injection mask and nim structure are from same subject\n');
    fprintf('     - Verify preprocessing steps are compatible\n\n');
else
    fprintf('✓ Tracks are passing through injection site\n');
    fprintf('  The niftiwrite error is likely a data type issue.\n\n');
end

function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
