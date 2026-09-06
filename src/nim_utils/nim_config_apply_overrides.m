function [config, applied] = nim_config_apply_overrides(config, sets)
% nim_config_apply_overrides: Apply CLI --set overrides to a loaded config.
%
%   config = nim_config_apply_overrides(config, {'integrator.step=0.05', ...})
%
% Every configurable parameter is reachable, not just tractography ones.
% Accepted key forms, in resolution order:
%
%   tractography.integrator.step   full canonical path
%   preprocessing.run_eddy         full canonical path (preprocessing reachable)
%   integrator.step                'tractography.' assumed
%   upsample                       bare leaf name, if unique across the schema
%   step_size                      legacy alias, mapped with a warning
%
% A bare leaf that is ambiguous (e.g. `method`, `fa_min`) is an error listing the
% candidate full paths. Unknown keys are an error, never a silent no-op - the
% previous behaviour assigned them blindly to config.tractography, where nothing
% read them.
%
% Values are parsed like YAML scalars, plus lists: [41,42] or 41,42 or SLF_L,SLF_R

    applied = {};
    if nargin < 2 || isempty(sets), return; end
    if ischar(sets) || isstring(sets), sets = cellstr(sets); end

    S = nim_config_schema();

    for i = 1:numel(sets)
        kv = strtrim(char(sets{i}));
        if isempty(kv), continue; end
        eq = find(kv == '=', 1);
        if isempty(eq)
            error('nim_config_apply_overrides:badFormat', ...
                '--set expects key=value, got "%s".', kv);
        end
        key = strtrim(kv(1:eq-1));
        val_str = strtrim(kv(eq+1:end));

        [path, entry] = resolve_key(key, S);
        value = parse_override_value(val_str, entry);
        value = check_override(value, entry, key);

        config = set_path(config, path, value);
        applied{end+1} = sprintf('%s = %s', path, val_str); %#ok<AGROW>
        fprintf('  override: %s = %s\n', path, val_str);
    end
end

% =========================================================================
function [path, entry] = resolve_key(key, S)
    paths = {S.path};

    % 1. exact canonical path
    idx = find(strcmp(paths, key), 1);
    if ~isempty(idx), path = paths{idx}; entry = S(idx); return; end

    % 2. 'tractography.' assumed
    cand = ['tractography.' key];
    idx = find(strcmp(paths, cand), 1);
    if ~isempty(idx), path = paths{idx}; entry = S(idx); return; end

    % 3. bare leaf name, if unique
    hits = [];
    for i = 1:numel(S)
        p = split(S(i).path, '.');
        if strcmp(p{end}, key), hits(end+1) = i; end %#ok<AGROW>
    end
    if numel(hits) == 1
        path = S(hits).path; entry = S(hits); return;
    elseif numel(hits) > 1
        cands = cell(1, numel(hits));
        for j = 1:numel(hits), cands{j} = S(hits(j)).path; end
        error('nim_config_apply_overrides:ambiguousKey', ...
            '--set "%s" is ambiguous; use a full path. Candidates: %s', ...
            key, strjoin(cands, ', '));
    end

    % 4. legacy alias
    for i = 1:numel(S)
        if any(strcmp(S(i).legacy, key))
            warning('nim_config_apply_overrides:legacyKey', ...
                '--set "%s" is deprecated; use "%s".', key, S(i).path);
            path = S(i).path; entry = S(i); return;
        end
    end

    error('nim_config_apply_overrides:unknownKey', ...
        '--set "%s" is not a known parameter.', key);
end

% =========================================================================
function value = parse_override_value(str, entry)
    if strcmp(entry.type, 'list')
        s = str;
        if startsWith(s, '[') && endsWith(s, ']'), s = s(2:end-1); end
        s = strtrim(s);
        if isempty(s), value = {}; return; end
        parts = strsplit(s, ',');
        value = cell(1, numel(parts));
        for i = 1:numel(parts)
            value{i} = parse_scalar(strtrim(parts{i}));
        end
        return;
    end
    value = parse_scalar(str);
end

function value = parse_scalar(str)
    if numel(str) >= 2 && ((startsWith(str,'''') && endsWith(str,'''')) || ...
                           (startsWith(str,'"')  && endsWith(str,'"')))
        value = str(2:end-1); return;
    end
    if strcmpi(str, 'true'),  value = true;  return; end
    if strcmpi(str, 'false'), value = false; return; end
    num = str2double(str);
    if ~isnan(num) && ~isempty(regexp(str, '^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$', 'once'))
        value = num; return;
    end
    value = str;
end

% =========================================================================
function v = check_override(v, entry, key)
    switch entry.type
        case 'numeric'
            if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
                error('nim_config_apply_overrides:badType', ...
                    '--set %s: expected a number.', key);
            end
            if ~isempty(entry.range) && (v < entry.range(1) || v > entry.range(2))
                error('nim_config_apply_overrides:outOfRange', ...
                    '--set %s: must be in [%g, %g], got %g.', key, entry.range(1), entry.range(2), v);
            end
        case 'logical'
            if isnumeric(v) && isscalar(v) && (v == 0 || v == 1), v = logical(v); end
            if ~islogical(v) || ~isscalar(v)
                error('nim_config_apply_overrides:badType', ...
                    '--set %s: expected true or false.', key);
            end
        case 'string'
            if ~(ischar(v) || isstring(v))
                error('nim_config_apply_overrides:badType', '--set %s: expected a string.', key);
            end
            v = lower(char(v));
            if ~isempty(entry.allowed) && ~any(strcmpi(v, entry.allowed))
                error('nim_config_apply_overrides:badValue', ...
                    '--set %s: must be one of {%s}, got "%s".', key, strjoin(entry.allowed, ', '), v);
            end
        case 'list'
            if ~iscell(v), v = {v}; end
    end
end

function s = set_path(s, path, value)
    p = split(path, '.');
    s = set_rec(s, p, value);
end

function s = set_rec(s, p, value)
    if numel(p) == 1
        s.(p{1}) = value;
    else
        if ~isfield(s, p{1}) || ~isstruct(s.(p{1})), s.(p{1}) = struct(); end
        s.(p{1}) = set_rec(s.(p{1}), p(2:end), value);
    end
end
