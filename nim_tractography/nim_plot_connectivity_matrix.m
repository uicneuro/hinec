function connectivity_matrix = nim_plot_connectivity_matrix(tracks, nim, varargin)
% nim_plot_connectivity_matrix: Compute and visualize connectivity matrix
%
% Arguments:
%   tracks - Cell array of fiber tracks
%   nim - NIM structure with parcellation
%   options - Options for connectivity analysis (optional struct)
%
% Returns:
%   connectivity_matrix - NxN matrix of connections between regions

% Parse input arguments
if nargin > 2 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default values
if ~isfield(options, 'min_track_length')
    options.min_track_length = 10;
end
if ~isfield(options, 'normalize')
    options.normalize = true;
end
if ~isfield(options, 'symmetric')
    options.symmetric = true;
end

if ~isfield(nim, 'parcellation_mask')
    error('Parcellation mask not found in nim structure');
end

% Get unique parcellation labels
parcel_labels = unique(nim.parcellation_mask(:));
parcel_labels = parcel_labels(parcel_labels > 0); % Remove background
n_regions = length(parcel_labels);

fprintf('Computing connectivity matrix for %d regions...\n', n_regions);

% Initialize connectivity matrix
connectivity_matrix = zeros(n_regions, n_regions);

% Process each track
valid_tracks = 0;
for i = 1:length(tracks)
    track = tracks{i};
    
    % Skip short tracks
    if size(track, 1) < options.min_track_length
        continue;
    end
    
    % Get parcellation labels along track
    track_labels = get_track_labels(track, nim.parcellation_mask);
    
    % Find start and end regions
    [start_region, end_region] = get_track_endpoints(track_labels, parcel_labels);
    
    if ~isempty(start_region) && ~isempty(end_region) && start_region ~= end_region
        % Add connection
        start_idx = find(parcel_labels == start_region);
        end_idx = find(parcel_labels == end_region);
        
        connectivity_matrix(start_idx, end_idx) = connectivity_matrix(start_idx, end_idx) + 1;
        
        if options.symmetric
            connectivity_matrix(end_idx, start_idx) = connectivity_matrix(end_idx, start_idx) + 1;
        end
        
        valid_tracks = valid_tracks + 1;
    end
end

fprintf('Used %d valid tracks for connectivity\n', valid_tracks);

% Normalize if requested
if options.normalize
    max_connections = max(connectivity_matrix(:));
    if max_connections > 0
        connectivity_matrix = connectivity_matrix / max_connections;
    end
end

% Visualize connectivity matrix
figure('Name', 'Connectivity Matrix', 'Color', 'w');

% Main connectivity matrix
subplot(2, 2, [1, 3]);
imagesc(connectivity_matrix);
colormap(hot);
colorbar;
title('Region-to-Region Connectivity Matrix');
xlabel('Target Region');
ylabel('Source Region');

% Add region labels if available
if isfield(nim, 'atlas_labels') && length(nim.atlas_labels) >= n_regions
    % Create abbreviated labels
    region_names = cell(n_regions, 1);
    for i = 1:n_regions
        if iscell(nim.atlas_labels)
            full_name = nim.atlas_labels{parcel_labels(i)};
        else
            full_name = sprintf('Region_%d', parcel_labels(i));
        end
        if length(full_name) > 10
            region_names{i} = full_name(1:10);
        else
            region_names{i} = full_name;
        end
    end
    
    set(gca, 'XTick', 1:n_regions, 'XTickLabel', region_names, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:n_regions, 'YTickLabel', region_names);
end

% Connection strength histogram
subplot(2, 2, 2);
connection_strengths = connectivity_matrix(connectivity_matrix > 0);
if ~isempty(connection_strengths)
    histogram(connection_strengths, 20);
    title('Distribution of Connection Strengths');
    xlabel('Connection Strength');
    ylabel('Frequency');
end

% Network metrics
subplot(2, 2, 4);
node_strengths = sum(connectivity_matrix, 2);
bar(node_strengths);
title('Node Strengths (Total Connections)');
xlabel('Region');
ylabel('Total Connections');

fprintf('Connectivity analysis complete\n');
end

function track_labels = get_track_labels(track, parcellation_mask)
% Get parcellation labels along a track
track_labels = zeros(size(track, 1), 1);

for i = 1:size(track, 1)
    pos = round(track(i, :));
    
    % Check bounds
    if all(pos >= 1) && all(pos <= size(parcellation_mask))
        track_labels(i) = parcellation_mask(pos(1), pos(2), pos(3));
    end
end
end

function [start_region, end_region] = get_track_endpoints(track_labels, parcel_labels)
% Find start and end regions of a track
start_region = [];
end_region = [];

% Find first non-zero label
valid_labels = track_labels(track_labels > 0);
if isempty(valid_labels)
    return;
end

start_region = valid_labels(1);

% Find last non-zero label
end_region = valid_labels(end);

% Ensure labels are in our parcellation
if ~ismember(start_region, parcel_labels)
    start_region = [];
end
if ~ismember(end_region, parcel_labels)
    end_region = [];
end
end 