function labels = nim_atlas_label_map(nim)
% nim_atlas_label_map: Return the atlas index -> region name map for a nim.
%
%   labels = nim_atlas_label_map(nim);   % containers.Map, keys are numeric
%
% Resolution order:
%   1. a map embedded in the nim (nim.atlas_labels.map) - written by the
%      parcellation step and travelling with the data
%   2. the sidecar <prefix>_atlas_labels.mat next to the nim, if recorded
%   3. nim_load_atlas_labels(atlas_type), which reads the FSL atlas XML
%
% Errors with an actionable message rather than returning an empty map, since
% every caller needs real names to resolve an ROI request.

    labels = [];

    % 1. embedded
    if isfield(nim, 'atlas_labels')
        al = nim.atlas_labels;
        if isa(al, 'containers.Map')
            labels = al;
        elseif isstruct(al) && isfield(al, 'map') && isa(al.map, 'containers.Map')
            labels = al.map;
        end
    end

    % 2. sidecar recorded on the nim
    if isempty(labels) && isfield(nim, 'atlas_labels_file') && ~isempty(nim.atlas_labels_file) ...
            && isfile(nim.atlas_labels_file)
        S = load(nim.atlas_labels_file);
        labels = map_from_struct(S);
    end

    % 3. rebuild from the FSL atlas XML
    if isempty(labels)
        atlas_type = 'jhu';
        if isfield(nim, 'atlas_type') && ~isempty(nim.atlas_type)
            atlas_type = char(nim.atlas_type);
        end
        try
            L = nim_load_atlas_labels(atlas_type);
            labels = map_from_struct(L);
        catch ME
            error('nim_atlas_label_map:unavailable', ...
                ['Could not obtain atlas label names (%s). ROIs can still be given by ' ...
                 'numeric index, but names require the label map. Original error: %s'], ...
                atlas_type, ME.message);
        end
    end

    if isempty(labels) || ~isa(labels, 'containers.Map') || labels.Count == 0
        error('nim_atlas_label_map:unavailable', ...
            'No usable atlas label map found for this nim.');
    end
end

% =========================================================================
function m = map_from_struct(S)
% Accept the several shapes these label files have been saved in.
    m = [];
    if isa(S, 'containers.Map'), m = S; return; end
    if ~isstruct(S), return; end

    if isfield(S, 'map') && isa(S.map, 'containers.Map'), m = S.map; return; end
    if isfield(S, 'atlas_labels')
        m = map_from_struct(S.atlas_labels); if ~isempty(m), return; end
    end
    % index/name array form
    if isfield(S, 'index') && isfield(S, 'name')
        m = containers.Map('KeyType','double','ValueType','char');
        for i = 1:numel(S)
            m(double(S(i).index)) = char(string(S(i).name));
        end
        return;
    end
    % single-field wrapper
    f = fieldnames(S);
    if numel(f) == 1
        m = map_from_struct(S.(f{1}));
    end
end
