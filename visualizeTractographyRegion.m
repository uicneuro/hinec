function visualizeTractographyRegion(region_id, tracks_file, nim_file, varargin)
% visualizeTractographyRegion: Visualize tracks for a specific brain region
%
% Usage:
%   visualizeTractographyRegion(5, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat')
%   visualizeTractographyRegion(5, 'tractography_results/tracks_*.mat', 'sample_parcellated.mat', 'start_in')
%
% Arguments:
%   region_id - Parcellation region index (integer)
%   tracks_file - Path to tracks .mat file (REQUIRED)
%   nim_file - Path to .mat file with nim structure (REQUIRED)
%   filter_mode - 'pass_through' (default), 'start_in', 'end_in', 'connect_to'
%   
% Optional name-value pairs:
%   'min_overlap' - Minimum overlap with region (0-1, default: 0.1)
%   'show_region' - Show region overlay (true/false, default: true)
%   'color_mode' - 'direction', 'fa', 'uniform', 'region' (default: 'direction')
%   'max_tracks' - Maximum tracks to display (default: unlimited)

% Validate required arguments
if nargin < 3
    error('visualizeTractographyRegion requires 3 arguments: region_id, tracks_file, nim_file');
end

if ~isnumeric(region_id) || region_id < 1
    error('region_id must be a positive integer');
end

if ~ischar(tracks_file) && ~isstring(tracks_file)
    error('tracks_file must be a string path to .mat file');
end

if ~ischar(nim_file) && ~isstring(nim_file)
    error('nim_file must be a string path to .mat file');
end

% Parse optional arguments
p = inputParser;
addOptional(p, 'filter_mode', 'pass_through', @ischar);
addParameter(p, 'min_overlap', 0.1, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'show_region', true, @islogical);
addParameter(p, 'color_mode', 'direction', @ischar);
addParameter(p, 'max_tracks', inf, @(x) isnumeric(x) && x > 0);
addParameter(p, 'region_alpha', 0.3, @(x) isnumeric(x) && x >= 0 && x <= 1);

parse(p, varargin{:});

% Extract parsed values
filter_mode = p.Results.filter_mode;

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

fprintf('=== Visualizing Region %d Tractography ===\n', region_id);

%% Load data files

% Load tracks data
fprintf('Loading tracks from %s...\n', tracks_file);
if ~exist(tracks_file, 'file')
    error('Tracks file not found: %s\nPlease run tractography first.', tracks_file);
end
track_data = load(tracks_file);
tracks = track_data.tracks;

% Load nim data
fprintf('Loading nim data from %s...\n', nim_file);
if ~exist(nim_file, 'file')
    error('Nim file not found: %s\nPlease run main() first.', nim_file);
end
nim_data = load(nim_file);
nim = nim_data.nim;

%% Validate data
if ~isfield(nim, 'parcellation_mask')
    error('nim structure must contain parcellation_mask. Please run parcellation first.');
end

if isempty(tracks)
    error('No tracks found in tracks file. Please run tractography first.');
end

% Check if region exists
max_region = max(nim.parcellation_mask(:));
if region_id > max_region
    error('Region %d does not exist. Available regions: 1-%d', region_id, max_region);
end

if sum(nim.parcellation_mask(:) == region_id) == 0
    error('Region %d contains no voxels in the parcellation mask', region_id);
end

%% Set up visualization options
options = struct();
options.filter_mode = filter_mode;
options.min_overlap = p.Results.min_overlap;
options.show_region = p.Results.show_region;
options.track_color = p.Results.color_mode;
options.max_tracks = p.Results.max_tracks;
options.region_alpha = p.Results.region_alpha;
options.show_stats = true;

%% Display region information
display_region_info(nim, region_id);

%% Visualize
nim_plot_tractography_region(tracks, nim, region_id, options);

fprintf('Region %d visualization complete!\n', region_id);
end


function display_region_info(nim, region_id)
% Display information about the selected region
fprintf('\n=== REGION INFORMATION ===\n');
fprintf('Region ID: %d\n', region_id);

% Get region name if available
region_name = '';

% Try new direct access arrays first (preferred)
if isfield(nim, 'labels') && region_id <= length(nim.labels) && ~isempty(nim.labels{region_id})
    region_name = nim.labels{region_id};
    fprintf('Region name: %s\n', region_name);
elseif isfield(nim, 'region_names') && region_id <= length(nim.region_names) && strlength(nim.region_names(region_id)) > 0
    region_name = char(nim.region_names(region_id));
    fprintf('Region name: %s\n', region_name);
% Fall back to old map access for compatibility
elseif isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
    try
        if nim.atlas_labels.map.isKey(region_id)
            region_name = nim.atlas_labels.map(region_id);
            fprintf('Region name: %s\n', region_name);
        end
    catch
        % Continue without name
    end
end

% Region statistics
region_voxels = sum(nim.parcellation_mask(:) == region_id);
total_voxels = sum(nim.parcellation_mask(:) > 0);
region_percentage = (region_voxels / total_voxels) * 100;

fprintf('Region size: %d voxels (%.1f%% of brain)\n', region_voxels, region_percentage);

% Available regions info
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0);
fprintf('Total regions available: %d (IDs: %d-%d)\n', ...
        length(unique_regions), min(unique_regions), max(unique_regions));

fprintf('========================\n\n');
end