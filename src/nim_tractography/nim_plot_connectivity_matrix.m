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
% Region membership that lets a voxel belong to SEVERAL regions.
%
% This used to read `unique(nim.parcellation_mask(:))` and then take one label
% per track point. A label volume has one owner per voxel, but anatomy does not:
% on the ISMRM bundle masks 84.8% of labelled voxels lie in two or more bundles
% and one lies in ten. Under the old lookup a streamline ending in a voxel shared
% by the corpus callosum and the cingulum was credited to whichever region won
% the tie, so edges were attributed to the wrong pair and the regions that lose
% ties (CC_u_shaped keeps 1% of its voxels in the label volume) were effectively
% absent from the matrix.
R = nim_region_lookup(nim);
parcel_labels = R.ids(:);
n_regions = numel(parcel_labels);
if R.n_labelled > 0
    fprintf('Connectivity: %d regions, %d of %d labelled voxels shared by >1 region (%.1f%%)\n', ...
        n_regions, R.n_multi, R.n_labelled, 100*R.n_multi/R.n_labelled);
end
if ~all(R.exact)
    warning('nim_plot_connectivity_matrix:approxRegions', ...
        ['%d of %d regions have no true mask and fall back to the label volume, ' ...
         'so their extent is whatever no other region claimed.'], nnz(~R.exact), n_regions);
end

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
    
    % Every region containing each ENDPOINT, not one label per endpoint.
    a = R.at(track(1, :));
    b = R.at(track(end, :));

    if ~isempty(a) && ~isempty(b)
        % A streamline whose endpoint sits in overlapping regions is evidence for
        % each pair it could join. Crediting only one would pick a winner
        % arbitrarily; crediting all of them is what the masks actually say.
        added = false;
        for ra = a(:)'
            for rb = b(:)'
                if ra == rb, continue; end
                start_idx = find(parcel_labels == ra, 1);
                end_idx   = find(parcel_labels == rb, 1);
                if isempty(start_idx) || isempty(end_idx), continue; end
                connectivity_matrix(start_idx, end_idx) = connectivity_matrix(start_idx, end_idx) + 1;
                if options.symmetric
                    connectivity_matrix(end_idx, start_idx) = connectivity_matrix(end_idx, start_idx) + 1;
                end
                added = true;
            end
        end
        if added
            valid_tracks = valid_tracks + 1;
        end
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
