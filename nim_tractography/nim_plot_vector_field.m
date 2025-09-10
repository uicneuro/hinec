function nim_plot_vector_field(nim, varargin)
% nim_plot_vector_field: Visualize eigenvector field on a slice
%
% Arguments:
%   nim - NIM structure with eigenvectors
%   options - Visualization options (optional struct)

% Parse input arguments
if nargin > 1 && isstruct(varargin{1})
    options = varargin{1};
else
    options = struct();
end

% Set default values
if ~isfield(options, 'slice')
    options.slice = [];
end
if ~isfield(options, 'downsample')
    options.downsample = 2;
end
if ~isfield(options, 'scale')
    options.scale = 1;
end
if ~isfield(options, 'axis_view')
    options.axis_view = "axial";
end

% Set default slice to middle
if isempty(options.slice)
    switch options.axis_view
        case "axial"
            options.slice = round(size(nim.FA, 3) / 2);
        case "coronal"
            options.slice = round(size(nim.FA, 2) / 2);
        case "sagittal"
            options.slice = round(size(nim.FA, 1) / 2);
    end
end

% Check if eigenvectors exist
if ~isfield(nim, 'evec')
    warning('Eigenvectors not found. Please run nim_eig() first.');
    return;
end

figure('Name', 'Vector Field Visualization', 'Color', 'w');

% Extract slice data
switch options.axis_view
    case "axial"
        fa_slice = nim.FA(:, :, options.slice)';
        evec_slice = squeeze(nim.evec(:, :, options.slice, :, 1)); % Primary eigenvector
        [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
        
    case "coronal"
        fa_slice = squeeze(nim.FA(:, options.slice, :))';
        evec_slice = squeeze(nim.evec(:, options.slice, :, :, 1));
        [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
        
    case "sagittal"
        fa_slice = squeeze(nim.FA(options.slice, :, :))';
        evec_slice = squeeze(nim.evec(options.slice, :, :, :, 1));
        [X, Y] = meshgrid(1:size(fa_slice, 2), 1:size(fa_slice, 1));
end

% Display FA background
imagesc(fa_slice);
colormap(gray);
hold on;

% Downsample for vector display
step = options.downsample;
X_down = X(1:step:end, 1:step:end);
Y_down = Y(1:step:end, 1:step:end);

% Get vector components for the slice
switch options.axis_view
    case "axial"
        U = evec_slice(1:step:end, 1:step:end, 1)'; % X component
        V = evec_slice(1:step:end, 1:step:end, 2)'; % Y component
        fa_down = fa_slice(1:step:end, 1:step:end);
        
    case "coronal"
        U = evec_slice(1:step:end, 1:step:end, 1)'; % X component
        V = evec_slice(1:step:end, 1:step:end, 3)'; % Z component
        fa_down = fa_slice(1:step:end, 1:step:end);
        
    case "sagittal"
        U = evec_slice(1:step:end, 1:step:end, 2)'; % Y component
        V = evec_slice(1:step:end, 1:step:end, 3)'; % Z component
        fa_down = fa_slice(1:step:end, 1:step:end);
end

% Mask vectors by FA threshold
fa_thresh = 0.2;
mask = fa_down > fa_thresh;
U(~mask) = 0;
V(~mask) = 0;

% Scale vectors
scale_factor = options.scale * step * 0.4;
U = U * scale_factor;
V = V * scale_factor;

% Plot vectors
quiver(X_down, Y_down, U, V, 0, 'r', 'LineWidth', 1.2);

% Formatting
axis equal;
axis tight;
title(sprintf('Primary Eigenvector Field - %s Slice %d', ...
              capitalize(options.axis_view), options.slice));
xlabel('X'); ylabel('Y');

% Add colorbar for FA
c = colorbar;
c.Label.String = 'Fractional Anisotropy';
caxis([0 1]);

hold off;
end

function str = capitalize(str)
% Capitalize first letter
str = char(str);
str(1) = upper(str(1));
end 