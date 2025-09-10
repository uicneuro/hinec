function visualizeAllRegionsSequentially(tracks_file, nim_file, varargin)
% visualizeAllRegionsSequentially: Run visualizeTractographyRegion for all regions one by one
%
% Usage:
%   visualizeAllRegionsSequentially('tractography_results/tracks_*.mat', 'sample_parcellated.mat')
%   visualizeAllRegionsSequentially(tracks_file, nim_file, 'filter_mode', 'start_in', 'pause_time', 3)
%
% Arguments:
%   tracks_file - Path to tracks .mat file (REQUIRED)
%   nim_file - Path to .mat file with nim structure (REQUIRED)
%
% Optional name-value pairs:
%   'pause_time' - Seconds to pause between regions (default: 2)
%   'filter_mode' - 'pass_through', 'start_in', 'end_in' (default: 'pass_through')
%   'close_previous' - Close previous figure before showing next (default: true)
%   'start_region' - Start from specific region ID (default: first available)
%   'end_region' - Stop at specific region ID (default: last available)
%   'show_empty_regions' - Show regions with no tracks (default: false)
%   'auto_advance' - Auto-advance without user input (default: true)
%   Additional parameters are passed to visualizeTractographyRegion

% Validate required arguments
if nargin < 2
    error('visualizeAllRegionsSequentially requires 2 arguments: tracks_file, nim_file');
end

if ~ischar(tracks_file) && ~isstring(tracks_file)
    error('tracks_file must be a string path to .mat file');
end

if ~ischar(nim_file) && ~isstring(nim_file)
    error('nim_file must be a string path to .mat file');
end

% Parse optional arguments
p = inputParser;
addParameter(p, 'pause_time', 2, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'filter_mode', 'pass_through', @ischar);
addParameter(p, 'close_previous', true, @islogical);
addParameter(p, 'start_region', 0, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'end_region', inf, @(x) isnumeric(x) && x > 0);
addParameter(p, 'show_empty_regions', false, @islogical);
addParameter(p, 'auto_advance', true, @islogical);

% Parse known parameters and keep unknown ones for visualizeTractographyRegion
[known_args, unknown_args] = parse_known_args(varargin, p);
parse(p, known_args{:});
options = p.Results;

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');

fprintf('=== Sequential Region Visualization ===\n');

%% Load nim data to get available regions
fprintf('Loading nim data from %s...\n', nim_file);
if ~exist(nim_file, 'file')
    error('Nim file not found: %s\nPlease run main() first.', nim_file);
end
nim_data = load(nim_file);
nim = nim_data.nim;

if ~isfield(nim, 'parcellation_mask')
    error('nim structure must contain parcellation_mask. Please run parcellation first.');
end

%% Get all available regions
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0); % Exclude background

if isempty(unique_regions)
    error('No regions found in parcellation mask');
end

% Apply region range filters
if options.start_region > 0
    unique_regions = unique_regions(unique_regions >= options.start_region);
end
if ~isinf(options.end_region)
    unique_regions = unique_regions(unique_regions <= options.end_region);
end

fprintf('Found %d regions to visualize: %s\n', length(unique_regions), ...
        mat2str(unique_regions'));

if options.auto_advance
    fprintf('Auto-advancing every %.1f seconds. Press Ctrl+C to stop.\n', options.pause_time);
else
    fprintf('Manual mode: Press any key to advance to next region.\n');
end

%% Pre-check regions for tracks (if not showing empty regions)
if ~options.show_empty_regions
    fprintf('Pre-checking regions for tracks...\n');
    
    % Load tracks for pre-checking
    if ~exist(tracks_file, 'file')
        error('Tracks file not found: %s\nPlease run tractography first.', tracks_file);
    end
    track_data = load(tracks_file);
    tracks = track_data.tracks;
    
    regions_with_tracks = [];
    for r = 1:length(unique_regions)
        region_id = unique_regions(r);
        
        % Quick check if region has tracks
        region_options = struct();
        region_options.filter_mode = options.filter_mode;
        region_options.min_overlap = 0.1;
        region_options.max_tracks = 1; % Just check if any exist
        
        try
            filtered_tracks = filter_tracks_by_region(tracks, nim.parcellation_mask, region_id, region_options);
            if ~isempty(filtered_tracks)
                regions_with_tracks(end+1) = region_id;
            end
        catch
            % Skip regions that cause errors
        end
    end
    
    unique_regions = regions_with_tracks;
    fprintf('Found %d regions with tracks: %s\n', length(unique_regions), ...
            mat2str(unique_regions));
end

if isempty(unique_regions)
    fprintf('No regions with tracks found.\n');
    return;
end

%% Sequential visualization
fprintf('\n=== Starting Sequential Visualization ===\n');
success_count = 0;
error_count = 0;

for r = 1:length(unique_regions)
    region_id = unique_regions(r);
    
    % Get region name for display
    region_name = get_region_name(nim, region_id);
    
    fprintf('\n--- Region %d/%d: ID=%d (%s) ---\n', ...
            r, length(unique_regions), region_id, region_name);
    
    try
        % Close previous figure if requested
        if options.close_previous && r > 1
            close(gcf);
        end
        
        % Call visualizeTractographyRegion with all parameters
        if isempty(unknown_args)
            visualizeTractographyRegion(region_id, tracks_file, nim_file, options.filter_mode);
        else
            visualizeTractographyRegion(region_id, tracks_file, nim_file, options.filter_mode, unknown_args{:});
        end
        
        success_count = success_count + 1;
        fprintf('✓ Region %d visualization successful\n', region_id);
        
        % Handle advancing to next region
        if r < length(unique_regions) % Not the last region
            if options.auto_advance
                fprintf('Auto-advancing in %.1f seconds...\n', options.pause_time);
                pause(options.pause_time);
            else
                fprintf('Press any key for next region (Ctrl+C to stop)...\n');
                pause;
            end
        end
        
    catch ME
        error_count = error_count + 1;
        fprintf('✗ Region %d visualization failed: %s\n', region_id, ME.message);
        
        % Brief pause even on error before continuing
        if r < length(unique_regions)
            pause(0.5);
        end
    end
end

%% Final summary
fprintf('\n=== Sequential Visualization Complete ===\n');
fprintf('Total regions processed: %d\n', length(unique_regions));
fprintf('Successful visualizations: %d\n', success_count);
if error_count > 0
    fprintf('Failed visualizations: %d\n', error_count);
end
fprintf('Success rate: %.1f%%\n', (success_count / length(unique_regions)) * 100);

if success_count > 0
    fprintf('\nAll successful region visualizations are complete!\n');
    if ~options.close_previous
        fprintf('Multiple figure windows may be open.\n');
    end
end

fprintf('==========================================\n');
end


function [known_args, unknown_args] = parse_known_args(args, parser)
% Separate known arguments from unknown ones
known_params = parser.Parameters;
known_args = {};
unknown_args = {};

i = 1;
while i <= length(args)
    if i < length(args) && ischar(args{i}) && any(strcmp(args{i}, known_params))
        % This is a known parameter
        known_args{end+1} = args{i};
        known_args{end+1} = args{i+1};
        i = i + 2;
    else
        % This is an unknown parameter or value
        unknown_args{end+1} = args{i};
        i = i + 1;
    end
end
end


function filtered_tracks = filter_tracks_by_region(tracks, parcellation_mask, region_id, options)
% Quick filter function for pre-checking regions
filtered_tracks = {};
track_count = 0;

for i = 1:length(tracks)
    track = tracks{i};
    if size(track, 1) < 2
        continue;
    end
    
    % Get parcellation labels along the track
    track_labels = get_track_parcellation_labels(track, parcellation_mask);
    
    % Apply filtering based on mode
    include_track = false;
    
    switch options.filter_mode
        case 'pass_through'
            include_track = any(track_labels == region_id);
        case 'start_in'
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                include_track = valid_labels(1) == region_id;
            end
        case 'end_in'
            valid_labels = track_labels(track_labels > 0);
            if ~isempty(valid_labels)
                include_track = valid_labels(end) == region_id;
            end
    end
    
    % Apply minimum overlap constraint
    if include_track && options.min_overlap > 0
        region_points = sum(track_labels == region_id);
        total_points = length(track_labels);
        overlap_ratio = region_points / total_points;
        include_track = overlap_ratio >= options.min_overlap;
    end
    
    if include_track
        track_count = track_count + 1;
        filtered_tracks{track_count} = track;
        
        % Early exit if we just need to check existence
        if track_count >= options.max_tracks
            break;
        end
    end
end
end


function track_labels = get_track_parcellation_labels(track, parcellation_mask)
% Get parcellation labels for each point along the track
track_labels = zeros(size(track, 1), 1);

for i = 1:size(track, 1)
    pos = round(track(i, :));
    
    % Check bounds
    if all(pos >= 1) && all(pos <= size(parcellation_mask))
        track_labels(i) = parcellation_mask(pos(1), pos(2), pos(3));
    end
end
end


function region_name = get_region_name(nim, region_id)
% Get the name of a region if available
region_name = '';

% Try new direct access arrays first (preferred)
if isfield(nim, 'labels') && region_id <= length(nim.labels) && ~isempty(nim.labels{region_id})
    region_name = nim.labels{region_id};
elseif isfield(nim, 'region_names') && region_id <= length(nim.region_names) && strlength(nim.region_names(region_id)) > 0
    region_name = char(nim.region_names(region_id));
% Fall back to old map access for compatibility
elseif isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
    try
        if nim.atlas_labels.map.isKey(region_id)
            region_name = nim.atlas_labels.map(region_id);
        end
    catch
        % Continue without name if there's an error
    end
end

if isempty(region_name)
    region_name = sprintf('Region %d', region_id);
end
end