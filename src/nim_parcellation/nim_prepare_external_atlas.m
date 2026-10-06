function [parcellation, atlas_labels] = nim_prepare_external_atlas(source, reference, destination, labels_file)
% Validate a DWI-aligned integer atlas and optionally stage a compressed copy.
% No registration, axis swapping, interpolation, or index offsets are applied.
% labels_file is an optional TSV with index and name columns. Every positive
% voxel label must have a name when a TSV is supplied; unused TSV rows are OK.
if nargin < 3, destination = ''; end
if nargin < 4, labels_file = ''; end
source = char(source); reference = char(reference);
destination = char(destination); labels_file = char(labels_file);
if ~isfile(source)
    error('nim:externalAtlasMissing', 'Configured atlas_file not found: %s', source);
end
if ~(endsWith(source, '.nii', 'IgnoreCase', true) || endsWith(source, '.nii.gz', 'IgnoreCase', true))
    error('nim:externalAtlasFormat', 'atlas_file must be .nii or .nii.gz.');
end
ref = niftiinfo(reference); info = niftiinfo(source);
if numel(info.ImageSize) ~= 3 || numel(ref.ImageSize) < 3 || ...
        ~isequal(info.ImageSize, ref.ImageSize(1:3)) || ...
        ~strcmp(info.SpaceUnits, ref.SpaceUnits) || ...
        any(abs(info.PixelDimensions-ref.PixelDimensions(1:3)) > 1e-5) || ...
        any(abs(info.Transform.T-ref.Transform.T) > 1e-5, 'all')
    error('nim:externalAtlasGrid', ...
        ['External atlas must match the DWI dimensions, voxel spacing, units and affine. ' ...
         'Register/resample labels to DWI space with nearest-neighbour interpolation before use.']);
end
data = double(niftiread(source));
if any(~isfinite(data) | data < 0 | data ~= round(data) | data > double(intmax('int32')), 'all')
    error('nim:externalAtlasValues', 'Atlas values must be finite, nonnegative int32-range integers (0 = background).');
end
ids = unique(data(data > 0));
if isempty(ids)
    error('nim:externalAtlasEmpty', 'External atlas has no positive region labels.');
end
map = containers.Map('KeyType', 'double', 'ValueType', 'char');
if isempty(labels_file)
    for id = ids', map(id) = sprintf('Region_%d', id); end
else
    if ~isfile(labels_file)
        error('nim:externalAtlasLabels', 'Configured atlas_labels_file not found: %s', labels_file);
    end
    opts = detectImportOptions(labels_file, 'FileType', 'text', 'Delimiter', '\t', ...
        'VariableNamingRule', 'preserve');
    if ~all(ismember({'index', 'name'}, opts.VariableNames))
        error('nim:externalAtlasLabels', 'Atlas label TSV must contain index and name columns.');
    end
    opts = setvartype(opts, {'index', 'name'}, 'string');
    table_data = readtable(labels_file, opts);
    indices = str2double(table_data.index);
    names = strtrim(table_data.name);
    if isempty(indices) || any(~isfinite(indices) | indices < 0 | indices ~= round(indices) | ...
            indices > double(intmax('int32'))) || numel(unique(indices)) ~= numel(indices) || ...
            any(ismissing(names) | strlength(names) == 0) || ...
            numel(unique(lower(names))) ~= numel(names)
        error('nim:externalAtlasLabels', 'Atlas TSV requires unique integer indices and nonempty, case-insensitively unique names.');
    end
    if ~all(ismember(ids, indices))
        error('nim:externalAtlasLabels', 'Atlas TSV must name every positive label present in the image.');
    end
    for k = 1:numel(indices), map(indices(k)) = char(names(k)); end
end
parcellation = int32(data);
[~, attr] = fileattrib(source);
atlas_labels = struct('map', map, 'atlas_type', 'external', ...
    'atlas_variant', attr.Name, 'source', attr.Name, 'labels_source', labels_file, ...
    'coordinate_policy', 'DWI-grid exact match; no registration or resampling');
if ~isempty(labels_file)
    [~, attr] = fileattrib(labels_file); atlas_labels.labels_source = attr.Name;
end
if ~isempty(destination)
    if ~endsWith(destination, '.nii.gz', 'IgnoreCase', true)
        error('nim:externalAtlasDestination', 'Staged atlas destination must end in .nii.gz.');
    end
    [same_exists, dst] = fileattrib(destination);
    if same_exists && strcmp(atlas_labels.source, dst.Name), return; end
    if endsWith(source, '.nii.gz', 'IgnoreCase', true)
        copyfile(source, destination);
    else
        scratch = tempname; mkdir(scratch);
        cleanup = onCleanup(@()rmdir(scratch, 's')); %#ok<NASGU>
        files = gzip(source, scratch);
        copyfile(files{1}, destination);
    end
end
end
