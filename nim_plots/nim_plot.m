function nim_plot(nim, varargin)
% nim_plot: Unified plotting function for diffusion tensor data
%
% Usage:
%   nim_plot(nim)                                    % Plot single voxel subset
%   nim_plot('filename.mat')                         % Load and plot
%   nim_plot(nim, 'mode', 'parcel', 'parcel_id', 5) % Plot ONE specific parcel
%   nim_plot(nim, 'mode', 'parcels')                % Plot ALL parcels (separate figures)
%   nim_plot(nim, 'indx', 1:10, 'indy', 1:10)      % Custom indices
%
% Options:
%   mode              - 'single'(default), 'parcel', 'parcels'
%   indx, indy, indz  - Voxel indices to plot (mode='single')
%   parcel_id         - Parcel ID to plot (mode='parcel')
%   figindex          - Figure number (default: 1)
%   downsample_factor - Downsample by factor (default: 1)
%   show_figure       - Create new figure (default: true)
%   show_progress     - Show pause/progress (default: true)
%
% Arrow Colors (DTI Standard):
%   Red   - Left/Right diffusion (X-direction)
%   Green - Anterior/Posterior diffusion (Y-direction)
%   Blue  - Superior/Inferior diffusion (Z-direction)
%   Mixed colors indicate diagonal diffusion directions

% Parse input arguments
p = inputParser;
addRequired(p, 'nim');
addParameter(p, 'mode', 'single');
addParameter(p, 'indx', []);
addParameter(p, 'indy', []);
addParameter(p, 'indz', []);
addParameter(p, 'parcel_id', []);
addParameter(p, 'figindex', 1);
addParameter(p, 'downsample_factor', 1);
addParameter(p, 'show_figure', true);
addParameter(p, 'show_progress', true);

parse(p, nim, varargin{:});
opts = p.Results;

% Ensure utils directory is in path
if ~exist('zwuni', 'file')
    addpath('utils');
end

% Load nim structure if filename provided
if ischar(opts.nim) || isstring(opts.nim)
    fprintf('Loading nim data from %s...\n', opts.nim);
    nim_data = load(opts.nim);
    nim = nim_data.nim;
else
    nim = opts.nim;
end

% Route to appropriate plotting mode
switch opts.mode
    case 'single'
        plot_single_region(nim, opts);
    case 'parcel'
        if isempty(opts.parcel_id)
            error('parcel_id must be specified for mode="parcel"');
        end
        plot_single_parcel(nim, opts);
    case 'parcels'
        plot_all_parcels(nim, opts);
    otherwise
        error('Unknown plotting mode: %s. Use "single", "parcel", or "parcels"', opts.mode);
end

end

%% Helper Functions

function plot_single_region(nim, opts)
% Plot a single region with specified indices

% Set default indices if not provided
if isempty(opts.indx)
    opts.indx = 1:nim.hdr.ImageSize(1);
end
if isempty(opts.indy)
    opts.indy = 1:nim.hdr.ImageSize(2);
end
if isempty(opts.indz)
    opts.indz = 1:nim.hdr.ImageSize(3);
end

% Get image dimensions
Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

% Apply downsampling if specified
indx = opts.indx;
indy = opts.indy;
indz = opts.indz;

if opts.downsample_factor > 1
    indx = indx(1:opts.downsample_factor:end);
    indy = indy(1:opts.downsample_factor:end);
    indz = indz(1:opts.downsample_factor:end);
end

% Create spatial coordinates
[X,~] = zwuni(Nvox_x);
[Y,~] = zwuni(Nvox_y);
[Z,~] = zwuni(Nvox_z);

X = X .* floor(Nvox_x/2);
Y = Y .* floor(Nvox_y/2);
Z = Z .* floor(Nvox_z/2);

Xc = 0.5.*(X(2:end) + X(1:end-1));
Yc = 0.5.*(Y(2:end) + Y(1:end-1));
Zc = 0.5.*(Z(2:end) + Z(1:end-1));

% Extract eigenvectors
Nvox = nim.hdr.ImageSize(1:3);
Vx = reshape(nim.evec(:, :, :, 1, 1), Nvox);
Vy = reshape(nim.evec(:, :, :, 1, 2), Nvox);
Vz = reshape(nim.evec(:, :, :, 1, 3), Nvox);

[XXc,YYc,ZZc] = meshgrid(Xc,Yc,Zc);

XXcind = XXc(indx,indy,indz);
YYcind = YYc(indx,indy,indz);
ZZcind = ZZc(indx,indy,indz);
Vxind = Vx(indx,indy,indz);
Vyind = Vy(indx,indy,indz);
Vzind = Vz(indx,indy,indz);

% Colors: RGB represents eigenvector direction (Red=L/R, Green=A/P, Blue=S/I)
colors = vector_to_color(Vxind, Vyind, Vzind);


if opts.show_figure
    figure(opts.figindex);
    clf; % Clear the figure
end

hold on
% Plot colored arrows - each arrow shows principal diffusion direction
arrow_count = 0;
for i = 1:numel(Vxind)
  if ~isnan(Vxind(i)) && ~isnan(Vyind(i)) && ~isnan(Vzind(i))
    quiver3(XXcind(i), YYcind(i), ZZcind(i), Vxind(i)*3, Vyind(i)*3, Vzind(i)*3, ...
            'Color', colors(i, :), 'LineWidth', 1.5, 'AutoScale', 'off', ...
            'MaxHeadSize', 0.8);
    arrow_count = arrow_count + 1;
  end
end
hold off

axis([-Nvox_x/2, Nvox_x/2, -Nvox_y/2, Nvox_y/2, -Nvox_z/2, Nvox_z/2])
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Diffusion Tensor Principal Eigenvectors (DTI Color Coding)');
grid on;
view(45, 30);  % Corner view to show all xyz axes

% Add legend for DTI colors
h1 = plot3(NaN,NaN,NaN,'Color',[1 0 0],'LineWidth',3); hold on;
h2 = plot3(NaN,NaN,NaN,'Color',[0 1 0],'LineWidth',3);
h3 = plot3(NaN,NaN,NaN,'Color',[0 0 1],'LineWidth',3);
legend([h1 h2 h3], {'Red: Left-Right (X)', 'Green: Anterior-Posterior (Y)', 'Blue: Superior-Inferior (Z)'}, 'Location', 'northeast');

drawnow
end


function plot_single_parcel(nim, opts)
% Plot vectors for a specific brain parcel (based on original nim_plotparcel)

parcel_id = opts.parcel_id;

% Get image dimensions
Nvox_x = nim.hdr.ImageSize(1);
Nvox_y = nim.hdr.ImageSize(2);
Nvox_z = nim.hdr.ImageSize(3);

% Extract the parcellation mask and the brain mask (ORIGINAL METHOD)
mask = nim.parcellation_mask;
brain_mask = nim.mask > 0;

% Find the indices of the voxels belonging to the specified parcel within the brain mask
[indx, indy, indz] = ind2sub(size(mask), find(mask == parcel_id & brain_mask));

if isempty(indx)
    fprintf('No voxels found for parcel %d\n', parcel_id);
    return;
end

% Vertex points (ORIGINAL METHOD)
[X,~] = zwuni(Nvox_x);
[Y,~] = zwuni(Nvox_y);
[Z,~] = zwuni(Nvox_z);

X = X .* floor(Nvox_x / 2);
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

% Convert subscripts to linear indices for consistent indexing (ORIGINAL METHOD)
linear_indices = sub2ind(size(mask), indx, indy, indz);

XXcind = XXc(linear_indices);
YYcind = YYc(linear_indices);
ZZcind = ZZc(linear_indices);
Vxind = Vx(linear_indices);
Vyind = Vy(linear_indices);
Vzind = Vz(linear_indices);

% Apply downsampling if needed
if opts.downsample_factor > 1
    step = opts.downsample_factor;
    idx = 1:step:length(XXcind);
    XXcind = XXcind(idx);
    YYcind = YYcind(idx);
    ZZcind = ZZcind(idx);
    Vxind = Vxind(idx);
    Vyind = Vyind(idx);
    Vzind = Vzind(idx);
end

% Colors: RGB represents eigenvector direction (Red=L/R, Green=A/P, Blue=S/I)
colors = vector_to_color(Vxind, Vyind, Vzind);

% Get parcel name
parcel_name = get_parcel_name(nim, parcel_id);

figure(opts.figindex);
hold on;
quiver3(XXcind, YYcind, ZZcind, Vxind, Vyind, Vzind, 'AutoScale', 'off');
for i = 1:numel(Vxind)
    quiver3(XXcind(i), YYcind(i), ZZcind(i), Vxind(i), Vyind(i), Vzind(i), 'Color', colors(i, :), 'LineWidth', 2);
end

axis([-Nvox_x / 2, Nvox_x / 2, -Nvox_y / 2, Nvox_y / 2, -Nvox_z / 2, Nvox_z / 2]);
xlabel('X');
ylabel('Y');
zlabel('Z');
view(45, 30);  % Corner view to show all xyz axes

% Create figure title
if ~isempty(parcel_name)
    title(['Parcel ' num2str(parcel_id) ': ' parcel_name]);
else
    title(['Parcel ' num2str(parcel_id)]);
end

drawnow
end

function plot_all_parcels(nim, opts)
% Plot all parcels sequentially in separate figures

num_parcels = max(nim.parcellation_mask(:));
figindex = opts.figindex;

% Use aggressive downsampling for all parcels to avoid memory issues
if opts.downsample_factor < 3
    opts.downsample_factor = 3;
end

for parcel_id = 1:num_parcels
  if sum(nim.parcellation_mask(:) == parcel_id) > 0
    sub_opts = opts;
    sub_opts.parcel_id = parcel_id;
    sub_opts.figindex = figindex;

    plot_single_parcel(nim, sub_opts);

    figindex = figindex + 1;
  end
end
end

function parcel_name = get_parcel_name(nim, parcel_id)
% Extract parcel name from nim structure
parcel_name = '';

if isfield(nim, 'labels') && parcel_id <= length(nim.labels) && ~isempty(nim.labels{parcel_id})
    parcel_name = nim.labels{parcel_id};
elseif isfield(nim, 'region_names') && parcel_id <= length(nim.region_names) && strlength(nim.region_names(parcel_id)) > 0
    parcel_name = char(nim.region_names(parcel_id));
elseif isfield(nim, 'atlas_labels') && isfield(nim.atlas_labels, 'map')
    try
        if nim.atlas_labels.map.isKey(parcel_id)
            parcel_name = nim.atlas_labels.map(parcel_id);
        end
    catch
        % Continue without name
    end
end
end