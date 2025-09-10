function fig = nim_plotparcellation(nim, varargin)
% NIM_PLOTPARCELLATION Plot all parcels in a single figure with different colors
%
% Arguments:
%   nim - The nim structure containing parcellation data
%   Optional name-value pairs:
%     'ParcelIDs' - Array of parcel IDs to include (default: all parcels)
%     'Keyword' - String to filter parcels by name (case-insensitive) (default: '')
%     'Colormap' - Colormap to use (default: 'jet')
%     'Alpha' - Transparency value (0-1) for arrows (default: 0.7)
%     'ScaleFactor' - Scaling factor for arrows (default: 1.0)
%     'ShowLegend' - Whether to show the legend (default: true)
%     'ShowAxes' - Whether to show the axes (default: true)
%     'TitleText' - Custom title text (default: 'Brain Parcellation')
%
% Returns:
%   fig - Handle to the created figure

% Parse inputs
p = inputParser;
addRequired(p, 'nim');
addParameter(p, 'ParcelIDs', [], @isnumeric);
addParameter(p, 'Keyword', '', @ischar);
addParameter(p, 'Colormap', 'jet', @ischar);
addParameter(p, 'Alpha', 0.7, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'ScaleFactor', 1.0, @isnumeric);
addParameter(p, 'ShowLegend', true, @islogical);
addParameter(p, 'ShowAxes', true, @islogical);
addParameter(p, 'TitleText', 'Brain Parcellation', @ischar);
parse(p, nim, varargin{:});

params = p.Results;

% Get dimensions
Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

% Extract the parcellation mask and the brain mask
mask = nim.parcellation_mask;
brain_mask = nim.mask > 0;

% Determine unique parcel IDs (excluding background usually labeled as 0)
unique_parcels = unique(mask(brain_mask & mask > 0));

% If specific parcels were requested, filter to only those
if ~isempty(params.ParcelIDs)
  unique_parcels = intersect(unique_parcels, params.ParcelIDs);
end

% If a keyword was specified, filter parcels by keyword in their name
keyword_filtered_parcels = [];
if ~isempty(params.Keyword) && isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
  keyword_filter = lower(params.Keyword); % Convert to lowercase for case-insensitive comparison
  
  for i = 1:length(unique_parcels)
    parcel_id = unique_parcels(i);
    try
      if nim.atlas_labels.map.isKey(parcel_id)
        parcel_name = nim.atlas_labels.map(parcel_id);
        % Case-insensitive search for keyword in parcel name
        if contains(lower(parcel_name), keyword_filter)
          keyword_filtered_parcels = [keyword_filtered_parcels, parcel_id];
        end
      end
    catch
      % If error occurs, skip this parcel
    end
  end
  
  % If we found matching parcels, use only those
  if ~isempty(keyword_filtered_parcels)
    unique_parcels = keyword_filtered_parcels;
    fprintf('Found %d parcels matching keyword "%s"\n', length(unique_parcels), params.Keyword);
  else
    warning('No parcels matching keyword "%s" were found.', params.Keyword);
  end
end

% Create the figure
fig = figure('Name', 'Parcellation Visualization');
hold on;

% Generate a set of distinct colors manually instead of using colormap
% This ensures we get different colors for each parcel
num_parcels = length(unique_parcels);
if strcmp(params.Colormap, 'jet')
  % Use a preset of distinct colors for better visualization
  distinct_colors = [
    1 0 0;      % Red
    0 1 0;      % Green
    0 0 1;      % Blue
    1 1 0;      % Yellow
    1 0 1;      % Magenta
    0 1 1;      % Cyan
    0.5 0 0;    % Dark red
    0 0.5 0;    % Dark green
    0 0 0.5;    % Dark blue
    0.5 0.5 0;  % Olive
    0.5 0 0.5;  % Purple
    0 0.5 0.5;  % Teal
    1 0.5 0;    % Orange
    0.5 0 1;    % Violet
    1 0 0.5;    % Pink
    0.5 1 0;    % Lime
    0.5 0.5 0.5;% Gray
    0.7 0.3 0.3;% Brown
    0 0.7 0.7;  % Turquoise
    0.9 0.9 0;  % Gold
    ];
  
  % If we need more colors than in our preset, generate them
  if num_parcels > size(distinct_colors, 1)
    % Generate additional colors using HSV colorspace for maximum distinction
    hsv_colors = hsv(num_parcels);
    distinct_colors = hsv_colors;
  end
  
  fprintf('Using %d distinct colors for %d parcels\n', size(distinct_colors, 1), num_parcels);
else
  % If user specified a different colormap, use that
  distinct_colors = eval([params.Colormap '(num_parcels)']);
end

% Prepare for legend
legend_entries = cell(length(unique_parcels), 1);

% Vertex points
[X,~] = zwuni(Nvox_x);     % [-1 1]
[Y,~] = zwuni(Nvox_y);
[Z,~] = zwuni(Nvox_z);

X = X .* floor(Nvox_x / 2);  % Scale. [-Nvox/2 Nvox/2]
Y = Y .* floor(Nvox_y / 2);
Z = Z .* floor(Nvox_z / 2);

% Center points
Xc = 0.5 .* (X(2:end) + X(1:end-1));
Yc = 0.5 .* (Y(2:end) + Y(1:end-1));
Zc = 0.5 .* (Z(2:end) + Z(1:end-1));

% Vectors for each voxel
Nvox = nim.hdr.ImageSize(1:3);
Vx = reshape(nim.evec(:, :, :, 1, 1), Nvox);
Vy = reshape(nim.evec(:, :, :, 1, 2), Nvox);
Vz = reshape(nim.evec(:, :, :, 1, 3), Nvox);

% Voxel Vertices
[XXc, YYc, ZZc] = meshgrid(Xc, Yc, Zc);

% Plot each parcel with a different color
for i = 1:length(unique_parcels)
  parcel_id = unique_parcels(i);
  
  % Calculate color for this parcel - ensure we don't exceed array bounds
  color_idx = min(i, size(distinct_colors, 1));
  color = distinct_colors(color_idx, :);
  
  fprintf('Parcel %d: Using color [%.2f, %.2f, %.2f]\n', parcel_id, color(1), color(2), color(3));
  
  % Find voxels belonging to this parcel
  [indx, indy, indz] = ind2sub(size(mask), find(mask == parcel_id & brain_mask));
  
  if isempty(indx)
    continue; % Skip if no voxels found
  end
  
  % Convert subscripts to linear indices for consistent indexing
  linear_indices = sub2ind(size(mask), indx, indy, indz);
  
  XXcind = XXc(linear_indices);
  YYcind = YYc(linear_indices);
  ZZcind = ZZc(linear_indices);
  Vxind = Vx(linear_indices);
  Vyind = Vy(linear_indices);
  Vzind = Vz(linear_indices);
  
  % Scale vectors by the scale factor
  Vxind = Vxind * params.ScaleFactor;
  Vyind = Vyind * params.ScaleFactor;
  Vzind = Vzind * params.ScaleFactor;
  
  % Get parcel name if available
  parcel_name = '';
  if isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
    try
      if nim.atlas_labels.map.isKey(parcel_id)
        parcel_name = nim.atlas_labels.map(parcel_id);
      end
    catch
      % If any error occurs, continue without the name
    end
  end
  
  % Create legend entry
  if ~isempty(parcel_name)
    legend_entries{i} = sprintf('Parcel %d: %s', parcel_id, parcel_name);
  else
    legend_entries{i} = sprintf('Parcel %d', parcel_id);
  end
  
  % Plot the parcel
  h = quiver3(XXcind, YYcind, ZZcind, Vxind, Vyind, Vzind, 'AutoScale', 'off', 'Color', color);
  
  % Apply transparency - simplified approach
  try
    h.Color(4) = params.Alpha;  % Try to set transparency directly
  catch
    % If that fails, don't worry about transparency
  end
end

% Set axis properties
axis([-Nvox_x/2, Nvox_x/2, -Nvox_y/2, Nvox_y/2, -Nvox_z/2, Nvox_z/2]);
xlabel('X');
ylabel('Y');
zlabel('Z');
title(params.TitleText);

% Show legend if requested
if params.ShowLegend
  % Limit legend to a reasonable number of entries if there are too many
  if length(legend_entries) > 20
    warning('Too many parcels for legend. Showing only first 20.');
    legend_entries = legend_entries(1:20);
    legend_entries{end} = '...more parcels';
  end
  legend(legend_entries, 'Location', 'eastoutside');
end

% Show/hide axes if requested
if ~params.ShowAxes
  axis off;
end

% Set up a nice view
view(3); % 3D view
grid on;
box on;
rotate3d on; % Allow rotation

hold off;
end

% Helper function (assuming this is defined elsewhere in your code)
function [x, xp] = zwuni(N)
% Computes uniformly distributed points in [-1,1]
x = linspace(-1, 1, N+1);
xp = x;
end