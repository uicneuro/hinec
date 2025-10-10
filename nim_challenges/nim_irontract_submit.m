function nim_irontract_submit(nim_file, injection_file, output_dir, varargin)
% nim_irontract_submit: Generate IronTract Challenge submission files
%
% Creates binary NIfTI voxel masks from tractography results for IronTract
% Challenge evaluation. Automatically filters streamlines by injection site
% and generates multiple volumes with varying parameters for ROC analysis.
%
% Arguments:
%   nim_file - Path to .mat file with nim structure OR nim structure itself
%   injection_file - Path to injection site binary mask (.nii.gz)
%   output_dir - Directory to save submission volumes
%   options - Optional structure with tractography parameters
%
% Options:
%   angle_thresholds - Array of angle thresholds for parameter sweep
%                     (default: [30, 45, 60, 75, 90])
%   tracks_file - Pre-computed tracks .mat file (skips tractography)
%   base_options - Base tractography options structure
%
% Output:
%   Saves numbered NIfTI volumes: submission_001.nii.gz, submission_002.nii.gz, ...
%   Each volume is a binary mask where 1 = voxel visited by streamlines
%   Headers match the injection site mask for proper alignment
%
% Usage Examples:
%   % Basic usage with default parameters
%   nim_irontract_submit('data.mat', 'injection.nii.gz', 'submissions/');
%
%   % Custom angle sweep
%   options.angle_thresholds = [20, 40, 60, 80];
%   nim_irontract_submit('data.mat', 'injection.nii.gz', 'out/', options);
%
%   % Use pre-computed tracks
%   options.tracks_file = 'tracks.mat';
%   nim_irontract_submit('data.mat', 'injection.nii.gz', 'out/', options);
%
% Reference:
%   Maffei et al. (2022). Insights from the IronTract challenge. NeuroImage.

% Parse input arguments
if nargin > 3 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set defaults
if ~isfield(options, 'angle_thresholds')
    options.angle_thresholds = [30, 45, 60, 75, 90];
end
if ~isfield(options, 'tracks_file')
    options.tracks_file = '';
end
if ~isfield(options, 'base_options')
    options.base_options = struct();
end

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Load nim structure if path provided
if ischar(nim_file) || isstring(nim_file)
    fprintf('Loading nim structure from %s...\n', nim_file);
    data = load(nim_file);
    nim = data.nim;
else
    nim = nim_file;
end

% Load injection site mask
fprintf('Loading injection site mask from %s...\n', injection_file);
injection_mask = niftiread(injection_file);
injection_info = niftiinfo(injection_file);

% Get dimensions from injection mask
dims = size(injection_mask);
fprintf('Injection mask dimensions: %d x %d x %d\n', dims);

% Verify dimensions match nim structure
if ~isequal(dims, [nim.xdim, nim.ydim, nim.zdim])
    warning('Dimension mismatch between nim structure and injection mask!');
    fprintf('NIM dims: %d x %d x %d\n', nim.xdim, nim.ydim, nim.zdim);
    fprintf('Injection dims: %d x %d x %d\n', dims);
end

% Parameter sweep: generate multiple volumes for ROC curve
num_params = length(options.angle_thresholds);
fprintf('\n=== IronTract Submission Generation ===\n');
fprintf('Generating %d volumes with angle thresholds: %s\n', ...
    num_params, mat2str(options.angle_thresholds));

for i = 1:num_params
    angle = options.angle_thresholds(i);
    fprintf('\n--- Volume %d/%d: Angle threshold = %.1f° ---\n', i, num_params, angle);

    % Set tractography options
    tract_opts = options.base_options;
    tract_opts.angle_thresh = angle;

    % Run tractography or load pre-computed tracks
    if ~isempty(options.tracks_file)
        fprintf('Loading pre-computed tracks from %s...\n', options.tracks_file);
        track_data = load(options.tracks_file);
        tracks = track_data.tracks;
    else
        fprintf('Running tractography with angle threshold %.1f°...\n', angle);
        tracks = nim_tractography_standard(nim, tract_opts);
    end

    fprintf('Generated %d total tracks\n', length(tracks));

    % Filter tracks by injection site
    fprintf('Filtering tracks by injection site...\n');
    filtered_tracks = filter_tracks_by_injection(tracks, injection_mask);
    fprintf('Filtered to %d tracks passing through injection site (%.1f%%)\n', ...
        length(filtered_tracks), 100*length(filtered_tracks)/length(tracks));

    % Convert to binary voxel mask
    fprintf('Converting tracks to voxel mask...\n');
    voxel_mask = tracks_to_voxel_mask(filtered_tracks, dims);
    num_voxels = sum(voxel_mask(:));
    fprintf('Binary mask contains %d voxels (%.2f%% of volume)\n', ...
        num_voxels, 100*num_voxels/prod(dims));

    % Save with matching header
    output_file = fullfile(output_dir, sprintf('submission_%03d.nii.gz', i));
    fprintf('Saving to %s...\n', output_file);
    niftiwrite(uint8(voxel_mask), output_file, injection_info, 'Compressed', true);
end

fprintf('\n=== Submission Generation Complete ===\n');
fprintf('Generated %d volumes in %s\n', num_params, output_dir);
fprintf('Ready for IronTract Challenge submission!\n');
end

%% Helper Functions

function filtered_tracks = filter_tracks_by_injection(tracks, injection_mask)
% Filter tracks that pass through injection site
%
% Arguments:
%   tracks - Cell array of tracks (each Nx3 matrix)
%   injection_mask - Binary 3D mask of injection site
%
% Returns:
%   filtered_tracks - Cell array of tracks intersecting injection site

dims = size(injection_mask);
filtered_tracks = {};

for i = 1:length(tracks)
    track = tracks{i};

    % Round track coordinates to voxel indices
    track_voxels = round(track);

    % Check each point along the track
    intersects = false;
    for j = 1:size(track_voxels, 1)
        pos = track_voxels(j, :);

        % Verify position is within bounds
        if all(pos >= 1) && all(pos <= dims)
            % Check if this voxel is in the injection site
            if injection_mask(pos(1), pos(2), pos(3)) > 0
                intersects = true;
                break;
            end
        end
    end

    % Keep track if it intersects injection site
    if intersects
        filtered_tracks{end+1} = track; %#ok<AGROW>
    end
end

% Convert to column cell array
filtered_tracks = filtered_tracks(:);
end

function voxel_mask = tracks_to_voxel_mask(tracks, dims)
% Convert tracks to binary voxel mask
%
% Arguments:
%   tracks - Cell array of tracks (each Nx3 matrix)
%   dims - Dimensions of output volume [x, y, z]
%
% Returns:
%   voxel_mask - Binary 3D array where 1 = voxel visited by streamlines

% Initialize empty mask
voxel_mask = zeros(dims);

% Mark all voxels visited by any track
for i = 1:length(tracks)
    track = tracks{i};

    % Round track coordinates to voxel indices
    track_voxels = round(track);

    % Mark each voxel along the track
    for j = 1:size(track_voxels, 1)
        pos = track_voxels(j, :);

        % Verify position is within bounds
        if all(pos >= 1) && all(pos <= dims)
            voxel_mask(pos(1), pos(2), pos(3)) = 1;
        end
    end
end
end
