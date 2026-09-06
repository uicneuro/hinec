function [mask, info] = nim_roi_mask(nim, roi_spec, dilate)
% nim_roi_mask: Build a binary mask from a list of atlas regions.
%
%   [mask, info] = nim_roi_mask(nim, {41, 42})
%   [mask, info] = nim_roi_mask(nim, {'SLF_L', 'Uncinate fasciculus R'}, 1)
%
% roi_spec is a cell array mixing atlas INDICES (numeric) and NAMES (char),
% freely. Names are matched against the atlas label map, case-insensitively,
% by exact match first, then by unique substring, then by short alias
% (e.g. SLF_L, CST_R, CC_genu).
%
% dilate (default 0) grows the result by that many voxels with a spherical
% structuring element. The JHU parcellation is thin - 9,641 labelled voxels in
% the ISMRM data, and Uncinate is 23-24 - so dilation is often necessary to get
% a usable seed count, not merely cosmetic.
%
% Returns info with per-label voxel counts at each stage, so an empty or tiny
% result can be diagnosed rather than guessed at.

    if nargin < 3 || isempty(dilate), dilate = 0; end

    if ~isfield(nim, 'parcellation_mask') || isempty(nim.parcellation_mask)
        error('nim_roi_mask:noParcellation', ...
            ['This nim has no parcellation_mask, so ROI selection is impossible. ' ...
             'Re-run main() with a parcellation step, or drop the roi setting.']);
    end
    P = nim.parcellation_mask;

    if ischar(roi_spec) || isstring(roi_spec), roi_spec = cellstr(roi_spec); end
    if ~iscell(roi_spec), roi_spec = {roi_spec}; end

    labels = nim_atlas_label_map(nim);

    % Names present only in roi_masks - endpoint and inclusion ROIs, which are
    % addressable regions but are deliberately not parcellation labels - resolve
    % directly, without going through the label map.
    extra = containers.Map('KeyType','char','ValueType','double');
    if isfield(nim, 'roi_masks') && isa(nim.roi_masks, 'containers.Map')
        lab_names = values(labels);
        for kk = keys(nim.roi_masks)
            if ~any(strcmpi(kk{1}, lab_names)), extra(lower(kk{1})) = 0; end
        end
    end

    idx = zeros(1, numel(roi_spec));
    names = cell(1, numel(roi_spec));
    for i = 1:numel(roi_spec)
        s_i = roi_spec{i};
        if (ischar(s_i) || isstring(s_i)) && isKey(extra, lower(strtrim(char(s_i))))
            kk = keys(nim.roi_masks);
            hit = kk{find(strcmpi(strtrim(char(s_i)), kk), 1)};
            idx(i) = 0; names{i} = hit;
            continue;
        end
        [idx(i), names{i}] = resolve_one(s_i, labels);
    end

    % --- build ------------------------------------------------------------
    % Prefer nim.roi_masks when it exists. A parcellation LABEL VOLUME gives each
    % voxel one owner, so where regions genuinely overlap - the temporal stem
    % belongs to the uncinate and the ILF, midline voxels to the corpus callosum
    % and whatever crosses it - the label volume can only keep one, and asking
    % for the loser returns a hollowed-out region. roi_masks holds the regions as
    % they were defined, overlaps intact, so a request for UF_right returns all
    % of UF_right rather than the part no larger structure claimed.
    have_true_masks = isfield(nim, 'roi_masks') && isa(nim.roi_masks, 'containers.Map');

    mask = false(size(P));
    per_label = zeros(1, numel(idx));
    for i = 1:numel(idx)
        m = [];
        if have_true_masks && isKey(nim.roi_masks, names{i})
            m = logical(nim.roi_masks(names{i}));
        end
        if isempty(m), m = (P == idx(i)); end
        per_label(i) = sum(m(:));
        mask = mask | m;
    end
    raw_count = sum(mask(:));

    if dilate > 0
        mask = imdilate(mask, strel('sphere', double(dilate)));
    end
    dilated_count = sum(mask(:));

    info = struct('labels', idx, 'names', {names}, 'per_label_voxels', per_label, ...
                  'raw_voxels', raw_count, 'dilated_voxels', dilated_count, ...
                  'dilate', dilate);

    % --- report -----------------------------------------------------------
    fprintf('  ROI selection (%d region(s)):\n', numel(idx));
    for i = 1:numel(idx)
        fprintf('    [%2d] %-58s %6d vox\n', idx(i), truncate(names{i}, 58), per_label(i));
        if per_label(i) == 0
            warning('nim_roi_mask:emptyLabel', ...
                'Atlas label %d ("%s") has no voxels in this parcellation.', idx(i), names{i});
        end
    end
    if dilate > 0
        fprintf('    dilate %d vox: %d -> %d voxels\n', dilate, raw_count, dilated_count);
    else
        fprintf('    total: %d voxels\n', raw_count);
    end
end

% =========================================================================
function [id, name] = resolve_one(spec, labels)
    if isnumeric(spec) && isscalar(spec)
        id = round(spec);
        if isKey(labels, id)
            name = labels(id);
        else
            name = '<not in label map>';
            warning('nim_roi_mask:unknownIndex', ...
                'Atlas index %d is not in the label map; using it anyway.', id);
        end
        return;
    end

    if ~(ischar(spec) || isstring(spec))
        error('nim_roi_mask:badSpec', 'ROI entries must be an atlas index or a name.');
    end
    q = strtrim(char(spec));

    ids = cell2mat(keys(labels));
    nm  = values(labels);

    % 1. exact, case-insensitive
    hit = find(strcmpi(nm, q), 1);
    if ~isempty(hit), id = ids(hit); name = nm{hit}; return; end

    % 2. short alias (SLF_L, CST_R, CC_genu, ...)
    [aid, aname] = resolve_alias(q, ids, nm);
    if ~isempty(aid), id = aid; name = aname; return; end

    % 3. unique case-insensitive substring
    hits = find(contains(lower(nm), lower(q)));
    if numel(hits) == 1
        id = ids(hits); name = nm{hits}; return;
    elseif numel(hits) > 1
        cand = cell(1, min(numel(hits), 6));
        for k = 1:numel(cand), cand{k} = sprintf('%d="%s"', ids(hits(k)), nm{hits(k)}); end
        error('nim_roi_mask:ambiguousName', ...
            'ROI name "%s" is ambiguous. Candidates: %s%s', q, strjoin(cand, ', '), ...
            repmat(' ...', 1, numel(hits) > 6));
    end

    error('nim_roi_mask:unknownName', ...
        ['ROI "%s" does not match any atlas region. Use an index, the full name, ' ...
         'or a short alias like SLF_L / CST_R / CC_genu.'], q);
end

% -------------------------------------------------------------------------
function [id, name] = resolve_alias(q, ids, nm)
% Short aliases for the JHU regions people actually type.
    id = []; name = '';
    t = upper(strrep(strtrim(q), ' ', '_'));

    side = '';
    if endsWith(t, '_L') || endsWith(t, '_LEFT'), side = 'L'; t = regexprep(t, '_(L|LEFT)$', ''); end
    if endsWith(t, '_R') || endsWith(t, '_RIGHT'), side = 'R'; t = regexprep(t, '_(R|RIGHT)$', ''); end

    base = containers.Map( ...
        {'SLF','CST','UF','ILF','IFOF','SFOF','ICP','SCP','MCP','CP','EC','ALIC','PLIC', ...
         'RLIC','ACR','SCR','PCR','PTR','SS','TAP','CGC','CGH','FX','CC_GENU','CC_BODY','CC_SPLENIUM','ML'}, ...
        {'Superior longitudinal fasciculus', 'Corticospinal tract', 'Uncinate fasciculus', ...
         'Inferior longitidinal fasciculus', 'Inferior fronto-occipital fasciculus', ...
         'Superior fronto-occipital fasciculus', 'Inferior cerebellar peduncle', ...
         'Superior cerebellar peduncle', 'Middle cerebellar peduncle', 'Cerebral peduncle', ...
         'External capsule', 'Anterior limb of internal capsule', ...
         'Posterior limb of internal capsule', 'Retrolenticular part of internal capsule', ...
         'Anterior corona radiata', 'Superior corona radiata', 'Posterior corona radiata', ...
         'Posterior thalamic radiation', 'Sagittal stratum', 'Tapetum', ...
         'Cingulum (cingulate gyrus)', 'Cingulum (hippocampus)', 'Fornix', ...
         'Genu of corpus callosum', 'Body of corpus callosum', 'Splenium of corpus callosum', ...
         'Medial lemniscus'});

    if ~isKey(base, t), return; end
    stem = base(t);

    for k = 1:numel(nm)
        n = nm{k};
        if ~contains(lower(n), lower(stem)), continue; end
        if isempty(side)
            % Unsided alias: accept only an unsided region name.
            if ~endsWith(n, ' R') && ~endsWith(n, ' L')
                id = ids(k); name = n; return;
            end
        else
            if endsWith(n, [' ' side]), id = ids(k); name = n; return; end
        end
    end
end

% -------------------------------------------------------------------------
function s = truncate(s, n)
    s = char(s);
    if numel(s) > n, s = [s(1:n-3) '...']; end
end
