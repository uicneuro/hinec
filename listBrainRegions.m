function region_info = listBrainRegions(data_file)
% listBrainRegions: Display available brain regions from parcellation
%
% Usage:
%   listBrainRegions('my_data.mat')       % List regions from specific .mat file
%   listBrainRegions()                    % List regions from sample_parcellated.mat
%   regions = listBrainRegions('data.mat'); % Return region info as struct array
%
% Arguments:
%   data_file - Path to .mat file containing nim structure with parcellation_mask
%
% Returns:
%   region_info - Struct array with fields: id, name, voxel_count, percentage

% Parse input
if nargin < 1
    data_file = 'sample_parcellated.mat';
end

% Add paths
addpath('nim_utils');

fprintf('=== BRAIN REGIONS AVAILABLE ===\n');

% Load nim data
if ~exist(data_file, 'file')
    error('Data file not found: %s', data_file);
end

fprintf('Loading data from %s...\n', data_file);
data = load(data_file);
nim = data.nim;

if ~isfield(nim, 'parcellation_mask')
    error('nim structure does not contain parcellation_mask field');
end

% Get unique regions (excluding background = 0)
unique_regions = unique(nim.parcellation_mask(:));
unique_regions = unique_regions(unique_regions > 0);

if isempty(unique_regions)
    error('No brain regions found in parcellation mask');
end

% Calculate total brain voxels
total_brain_voxels = sum(nim.parcellation_mask(:) > 0);

% Initialize region info
region_info = struct('id', {}, 'name', {}, 'voxel_count', {}, 'percentage', {});

fprintf('\nRegion ID | Region Name                    | Voxels | %% Brain\n');
fprintf('----------|--------------------------------|--------|--------\n');

for i = 1:length(unique_regions)
    region_id = unique_regions(i);
    
    % Count voxels
    voxel_count = sum(nim.parcellation_mask(:) == region_id);
    percentage = (voxel_count / total_brain_voxels) * 100;
    
    % Get region name
    region_name = 'Unknown';
    
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
            % Keep default name
        end
    end
    
    % Truncate long names for display
    display_name = region_name;
    if length(display_name) > 30
        display_name = [display_name(1:27) '...'];
    end
    
    fprintf('%8d | %-30s | %6d | %6.1f\n', ...
            region_id, display_name, voxel_count, percentage);
    
    % Store in return structure
    region_info(end+1).id = region_id;
    region_info(end).name = region_name;
    region_info(end).voxel_count = voxel_count;
    region_info(end).percentage = percentage;
end

fprintf('\nTotal regions: %d\n', length(unique_regions));
fprintf('Total brain voxels: %d\n', total_brain_voxels);

% Show usage examples
% fprintf('\n=== USAGE EXAMPLES ===\n');
% fprintf('Visualize region 5:\n');
% fprintf('  visualizeTractographyRegion(5)\n\n');
% 
% fprintf('Show tracks starting in region 12:\n');
% fprintf('  visualizeTractographyRegion(12, ''start_in'')\n\n');
% 
% fprintf('Show tracks connecting to region 8:\n');
% fprintf('  visualizeTractographyRegion(8, ''connect_to'')\n\n');
% 
% fprintf('Advanced usage with custom settings:\n');
% fprintf('  visualizeTractographyRegion(5, ''pass_through'', ...\n');
% fprintf('                             ''color_mode'', ''fa'', ...\n');
% fprintf('                             ''max_tracks'', 500, ...\n');
% fprintf('                             ''show_region'', true)\n\n');
% 
% fprintf('Available filter modes:\n');
% fprintf('  - ''pass_through'': Tracks that pass through the region (default)\n');
% fprintf('  - ''start_in'': Tracks that start in the region\n');
% fprintf('  - ''end_in'': Tracks that end in the region\n');
% fprintf('  - ''connect_to'': Tracks that connect the region to other regions\n\n');
% 
% fprintf('Available color modes:\n');
% fprintf('  - ''direction'': Color by fiber direction (default)\n');
% fprintf('  - ''fa'': Color by FA values along track\n');
% fprintf('  - ''uniform'': Single color for all tracks\n');
% fprintf('  - ''region'': Color by parcellation region\n');
% 
% fprintf('=============================\n');

if nargout == 0
    clear region_info;
end
end