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

% Debug: Show parsed arguments
if ischar(opts.nim)
    fprintf('DEBUG: nim argument is a string: "%s"\n', opts.nim);
end
fprintf('DEBUG: mode = "%s"\n', opts.mode);
if ~isempty(opts.parcel_id)
    fprintf('DEBUG: parcel_id = %s (class: %s)\n', mat2str(opts.parcel_id), class(opts.parcel_id));
end

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
% Plot vectors for a specific brain parcel

parcel_id = opts.parcel_id;

% Debug: Check parcel_id
fprintf('parcel_id type: %s, size: [%s], value: %s\n', ...
        class(parcel_id), num2str(size(parcel_id)), mat2str(parcel_id));

% Ensure parcel_id is a scalar
if ~isscalar(parcel_id)
    error('parcel_id must be a scalar, got size %s', mat2str(size(parcel_id)));
end

% Get dimensions from FA volume
dims = size(nim.FA);
Nvox_x = dims(1);
Nvox_y = dims(2);
Nvox_z = dims(3);

fprintf('Plotting parcel %d (volume: %dx%dx%d)\n', parcel_id, Nvox_x, Nvox_y, Nvox_z);

% Extract the parcellation mask and the brain mask
mask = nim.parcellation_mask;
brain_mask = nim.mask > 0;

% Debug: verify sizes match
fprintf('  mask size: [%s]\n', num2str(size(mask)));
fprintf('  brain_mask size: [%s]\n', num2str(size(brain_mask)));
if ~isequal(size(mask), size(brain_mask))
    error('Dimension mismatch: mask is %s but brain_mask is %s', ...
          mat2str(size(mask)), mat2str(size(brain_mask)));
end

% Find voxels belonging to this parcel
linear_indices = find(mask == parcel_id & brain_mask);

if isempty(linear_indices)
    fprintf('No voxels found for parcel %d\n', parcel_id);
    return;
end

fprintf('Found %d voxels in parcel %d\n', length(linear_indices), parcel_id);

% Create simple voxel coordinate grids (1-based indexing)
[X_grid, Y_grid, Z_grid] = meshgrid(1:Nvox_x, 1:Nvox_y, 1:Nvox_z);
X_grid = permute(X_grid, [2, 1, 3]);
Y_grid = permute(Y_grid, [2, 1, 3]);
Z_grid = permute(Z_grid, [2, 1, 3]);

% Extract primary eigenvector components
Vx = nim.evec(:, :, :, 1, 1);
Vy = nim.evec(:, :, :, 1, 2);
Vz = nim.evec(:, :, :, 1, 3);

% Extract coordinates and vectors for this parcel
X_coords = X_grid(linear_indices);
Y_coords = Y_grid(linear_indices);
Z_coords = Z_grid(linear_indices);
Vx_coords = Vx(linear_indices);
Vy_coords = Vy(linear_indices);
Vz_coords = Vz(linear_indices);

% Apply downsampling if needed
if opts.downsample_factor > 1
    step = opts.downsample_factor;
    idx = 1:step:length(X_coords);
    X_coords = X_coords(idx);
    Y_coords = Y_coords(idx);
    Z_coords = Z_coords(idx);
    Vx_coords = Vx_coords(idx);
    Vy_coords = Vy_coords(idx);
    Vz_coords = Vz_coords(idx);
end

% Scale vectors for visibility
scale_factor = 5.0;
Vx_coords = Vx_coords * scale_factor;
Vy_coords = Vy_coords * scale_factor;
Vz_coords = Vz_coords * scale_factor;

% Colors: RGB represents eigenvector direction (Red=L/R, Green=A/P, Blue=S/I)
colors = vector_to_color(Vx_coords/scale_factor, Vy_coords/scale_factor, Vz_coords/scale_factor);

% Get parcel name
parcel_name = get_parcel_name(nim, parcel_id);

figure(opts.figindex);
clf;
hold on;

% Plot all arrows with individual colors
for i = 1:length(X_coords)
    quiver3(X_coords(i), Y_coords(i), Z_coords(i), ...
            Vx_coords(i), Vy_coords(i), Vz_coords(i), ...
            'Color', colors(i, :), 'LineWidth', 1.5, 'AutoScale', 'off', 'MaxHeadSize', 0.5);
end

% Set proper axis limits - fixed for all parcels to show full brain volume
xlim([1, Nvox_x]);
ylim([1, Nvox_y]);
zlim([1, Nvox_z]);
xlabel('X (voxels)');
ylabel('Y (voxels)');
zlabel('Z (voxels)');
grid on;
box on;
view(45, 30);  % Corner view to show all xyz axes
daspect([1 1 1]);  % Equal data unit lengths in all directions

% Create figure title
if ~isempty(parcel_name)
    title(sprintf('Parcel %d: %s (%d voxels)', parcel_id, parcel_name, length(X_coords)));
else
    title(sprintf('Parcel %d (%d voxels)', parcel_id, length(X_coords)));
end

% Add DTI color legend
h1 = plot3(NaN,NaN,NaN,'Color',[1 0 0],'LineWidth',3);
h2 = plot3(NaN,NaN,NaN,'Color',[0 1 0],'LineWidth',3);
h3 = plot3(NaN,NaN,NaN,'Color',[0 0 1],'LineWidth',3);
legend([h1 h2 h3], {'Red: Left-Right (X)', 'Green: Anterior-Posterior (Y)', 'Blue: Superior-Inferior (Z)'}, 'Location', 'northeast');

hold off;
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