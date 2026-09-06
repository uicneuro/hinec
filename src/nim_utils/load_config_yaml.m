function config = load_config_yaml(config_file)
% load_config_yaml: Load, migrate, validate and default a HINEC config.
%
%   config = load_config_yaml('config/ismrm2015.yml')
%
% Returns the CANONICAL NESTED config (section.group.key), fully defaulted.
% Trackers do not read this directly - runTractography converts it with
% nim_config_to_options, which is the single place legacy option names live.
%
% Pipeline:
%   1. nim_yaml_parse    - indentation-aware parse, rejects >2 levels of nesting
%   2. migrate_legacy    - old flat keys -> canonical paths, with warnings
%   3. validate/default  - unknown-key rejection, types, enums, ranges
%
% The schema (nim_config_schema) is the single source of truth; retired keys
% are declared in nim_config_retired.

    if ~isfile(config_file)
        error('load_config_yaml:notFound', 'Configuration file not found: %s', config_file);
    end

    fprintf('Loading configuration from: %s\n', config_file);

    raw = nim_yaml_parse(config_file);
    S   = nim_config_schema();
    R   = nim_config_retired();

    [canon, provided, from_legacy] = migrate_legacy(raw, S, R, config_file);
    config = validate_and_default(canon, provided, from_legacy, S, config_file);

    fprintf('Configuration loaded successfully\n');
end

% =========================================================================
function [out, provided, from_legacy] = migrate_legacy(raw, S, R, src)
% Move legacy flat keys onto canonical paths; drop retired keys with a warning.
% `provided` records canonical paths the user actually supplied (as opposed to
% defaulted), so combination checks only fire on explicit settings.

    out = struct();
    provided = {};
    from_legacy = {};

    % Map legacy flat name -> canonical path (simple 1:1 renames).
    legacy_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for i = 1:numel(S)
        for j = 1:numel(S(i).legacy)
            legacy_map(S(i).legacy{j}) = S(i).path;
        end
    end
    retired_map = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for i = 1:numel(R)
        retired_map(R(i).key) = R(i).reason;
    end

    % Canonical group names per section, so we can tell a group from a key.
    groups = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(S)
        p = split(S(i).path, '.');
        if numel(p) == 3
            sec = p{1};
            if ~isKey(groups, sec), groups(sec) = {}; end
            g = groups(sec);
            if ~any(strcmp(g, p{2})), g{end+1} = p{2}; groups(sec) = g; end %#ok<AGROW>
        end
    end

    % Deferred transforms that need more than a rename.
    pending = struct('integration_order', [], 'max_steps', []);

    secs = fieldnames(raw);
    for si = 1:numel(secs)
        sec = secs{si};
        if ~any(strcmp(sec, {'preprocessing', 'tractography'}))
            error('load_config_yaml:unknownSection', ...
                '%s: unknown top-level section "%s" (expected preprocessing or tractography).', src, sec);
        end
        body = raw.(sec);
        keys = fieldnames(body);
        for ki = 1:numel(keys)
            k = keys{ki};
            v = body.(k);

            % A nested group: keys inside are already canonical.
            if isstruct(v)
                gk = fieldnames(v);
                for gi = 1:numel(gk)
                    path = sprintf('%s.%s.%s', sec, k, gk{gi});
                    out = set_path(out, path, v.(gk{gi}));
                    provided{end+1} = path; %#ok<AGROW>
                end
                continue;
            end

            % Transforming migrations: deferred because they need another key's
            % value (max_steps needs step) or a value mapping (integration_order).
            if strcmp(k, 'max_steps')
                pending.max_steps = v;
                continue;
            end
            if strcmp(k, 'integration_order')
                pending.integration_order = v;
                continue;
            end

            % Retired: warn and drop.
            if isKey(retired_map, k)
                warning('load_config_yaml:retiredKey', ...
                    '%s: "%s" is retired and has been ignored. %s', src, k, retired_map(k));
                continue;
            end

            % Simple legacy rename.
            if isKey(legacy_map, k)
                path = legacy_map(k);
                % Value-level migration: interp_method 'none' never actually
                % disabled interpolation - hinec.m:237 maps anything that is not
                % 'cubic' to 'linear'. Preserve behaviour, name it honestly.
                if strcmp(k, 'interp_method') && (ischar(v) || isstring(v)) && strcmpi(char(v), 'none')
                    warning('load_config_yaml:legacyKey', ...
                        ['%s: "interp_method: none" never disabled interpolation - the tracker ' ...
                         'mapped any non-cubic value to linear. Migrated to ' ...
                         '"interpolation.method: trilinear" (identical behaviour).'], src);
                    v = 'trilinear';
                else
                    warning('load_config_yaml:legacyKey', ...
                        '%s: "%s" is deprecated; use "%s".', src, k, path);
                end
                out = set_path(out, path, v);
                provided{end+1} = path;      %#ok<AGROW>
                from_legacy{end+1} = path;   %#ok<AGROW>
                continue;
            end

            % Already canonical at section level (e.g. tractography.algorithm)?
            path = sprintf('%s.%s', sec, k);
            if any(strcmp({S.path}, path))
                out = set_path(out, path, v);
                provided{end+1} = path; %#ok<AGROW>
                continue;
            end

            % A bare group name with a non-struct value.
            if isKey(groups, sec) && any(strcmp(groups(sec), k))
                error('load_config_yaml:groupNotScalar', ...
                    '%s: "%s.%s" is a group and cannot take a value directly.', src, sec, k);
            end

            error('load_config_yaml:unknownKey', ...
                '%s: unknown key "%s.%s".%s', src, sec, k, suggest(k, S));
        end
    end

    % --- deferred: integration_order -> integrator.method --------------------
    if ~isempty(pending.integration_order)
        io = pending.integration_order;
        switch io
            case 1, m = 'euler';
            case 2, m = 'rk2';
            case 4, m = 'rk4';
            case 5, m = 'rkf45';
            otherwise
                error('load_config_yaml:badIntegrationOrder', ...
                    '%s: integration_order must be 1, 2, 4 or 5, got %g.', src, io);
        end
        warning('load_config_yaml:legacyKey', ...
            ['%s: "integration_order: %g" is deprecated; use "integrator.method: %s". ' ...
             '(It named a METHOD, not an order - 5 meant RKF45, a 4(5) pair of order 4.)'], ...
            src, io, m);
        p = 'tractography.integrator.method';
        if ~any(strcmp(provided, p))
            out = set_path(out, p, m);
            provided{end+1} = p;
        end
    end

    % --- deferred: max_steps -> termination.max_arc --------------------------
    if ~isempty(pending.max_steps)
        step = get_path(out, 'tractography.integrator.step', []);
        if isempty(step)
            step = schema_default(S, 'tractography.integrator.step');
        end
        max_arc = pending.max_steps * step;
        warning('load_config_yaml:legacyKey', ...
            ['%s: "max_steps: %g" is retired; converted to "termination.max_arc: %g" ' ...
             '(= max_steps x step %g). max_arc is step-invariant; max_steps is now derived.'], ...
            src, pending.max_steps, max_arc, step);
        p = 'tractography.termination.max_arc';
        if ~any(strcmp(provided, p))
            out = set_path(out, p, max_arc);
            provided{end+1} = p;
        end
    end
end

% =========================================================================
function config = validate_and_default(canon, provided, from_legacy, S, src)
% Type/enum/range check every provided value, then fill defaults.

    config = canon;

    for i = 1:numel(S)
        path = S(i).path;
        if any(strcmp(provided, path))
            v = get_path(config, path, []);
            v = check_value(v, S(i), src);
            config = set_path(config, path, v);
        else
            config = set_path(config, path, S(i).default);
        end
    end

    % --- combination checks, only on EXPLICITLY provided keys ---------------
    method = get_path(config, 'tractography.integrator.method', 'rk4');
    if ~strcmp(method, 'rkf45')
        adaptive_only = {'tractography.integrator.tolerance', ...
                         'tractography.integrator.step_min', ...
                         'tractography.integrator.step_max', ...
                         'tractography.integrator.safety', ...
                         'tractography.integrator.adaptive'};
        for i = 1:numel(adaptive_only)
            if any(strcmp(provided, adaptive_only{i}))
                if any(strcmp(from_legacy, adaptive_only{i}))
                    % Came from a legacy key in an old config, where setting an
                    % inert adaptive param was tolerated. Warn, do not break it.
                    warning('load_config_yaml:adaptiveWithFixedStep', ...
                        ['%s: "%s" applies only to integrator.method: rkf45, but the method ' ...
                         'is "%s", so it has no effect. Remove it from the config.'], ...
                        src, adaptive_only{i}, method);
                else
                    error('load_config_yaml:adaptiveWithFixedStep', ...
                        ['%s: "%s" applies only to integrator.method: rkf45, but the method ' ...
                         'is "%s". Remove it, or switch method to rkf45.'], ...
                        src, adaptive_only{i}, method);
                end
            end
        end
    end

    smin = get_path(config, 'tractography.integrator.step_min', 0);
    smax = get_path(config, 'tractography.integrator.step_max', inf);
    if strcmp(method, 'rkf45') && smin >= smax
        error('load_config_yaml:badStepBounds', ...
            '%s: integrator.step_min (%g) must be < step_max (%g).', src, smin, smax);
    end

    fprintf('  Parameter validation passed\n');
end

% =========================================================================
function v = check_value(v, entry, src)
    switch entry.type
        case 'numeric'
            if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v)
                error('load_config_yaml:badType', '%s: "%s" must be a finite number, got "%s".', ...
                    src, entry.path, val2str(v));
            end
            if ~isempty(entry.range) && (v < entry.range(1) || v > entry.range(2))
                error('load_config_yaml:outOfRange', '%s: "%s" must be in [%g, %g], got %g.', ...
                    src, entry.path, entry.range(1), entry.range(2), v);
            end
        case 'logical'
            if islogical(v) && isscalar(v), return; end
            if isnumeric(v) && isscalar(v) && (v == 0 || v == 1), v = logical(v); return; end
            error('load_config_yaml:badType', '%s: "%s" must be true or false, got "%s".', ...
                src, entry.path, val2str(v));
        case 'string'
            if ~(ischar(v) || isstring(v))
                error('load_config_yaml:badType', '%s: "%s" must be a string, got "%s".', ...
                    src, entry.path, val2str(v));
            end
            v = char(v);
            if ~isempty(entry.allowed) && ~any(strcmpi(v, entry.allowed))
                error('load_config_yaml:badValue', '%s: "%s" must be one of {%s}, got "%s".', ...
                    src, entry.path, strjoin(entry.allowed, ', '), v);
            end
            v = lower(v);
        case 'list'
            if isempty(v), v = {}; return; end
            if ~iscell(v), v = {v}; end
    end
end

% =========================================================================
function s = suggest(key, S)
% "did you mean" against canonical leaf names and legacy aliases.
    cands = {};
    for i = 1:numel(S)
        p = split(S(i).path, '.');
        cands{end+1} = p{end};   %#ok<AGROW>
        cands = [cands, S(i).legacy]; %#ok<AGROW>
    end
    cands = unique(cands);
    best = ''; bestd = inf;
    for i = 1:numel(cands)
        d = edit_distance(lower(key), lower(cands{i}));
        if d < bestd, bestd = d; best = cands{i}; end
    end
    if bestd <= max(2, floor(numel(key)/3)) && ~isempty(best)
        % Report the full canonical path for the suggestion when we have one.
        full = best;
        for i = 1:numel(S)
            p = split(S(i).path, '.');
            if strcmp(p{end}, best) || any(strcmp(S(i).legacy, best))
                full = S(i).path; break;
            end
        end
        s = sprintf(' Did you mean "%s"?', full);
    else
        s = '';
    end
end

function d = edit_distance(a, b)
    m = numel(a); n = numel(b);
    D = zeros(m+1, n+1);
    D(:,1) = (0:m)'; D(1,:) = 0:n;
    for i = 1:m
        for j = 1:n
            c = double(a(i) ~= b(j));
            D(i+1,j+1) = min([D(i,j+1)+1, D(i+1,j)+1, D(i,j)+c]);
        end
    end
    d = D(m+1, n+1);
end

function s = val2str(v)
    try
        if ischar(v) || isstring(v), s = char(v);
        elseif islogical(v), s = mat2str(v);
        elseif isnumeric(v), s = mat2str(v);
        elseif iscell(v), s = sprintf('<list of %d>', numel(v));
        else, s = class(v);
        end
    catch
        s = class(v);
    end
end

function d = schema_default(S, path)
    d = [];
    for i = 1:numel(S)
        if strcmp(S(i).path, path), d = S(i).default; return; end
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

function v = get_path(s, path, dflt)
    p = split(path, '.');
    v = dflt;
    cur = s;
    for i = 1:numel(p)
        if ~isstruct(cur) || ~isfield(cur, p{i}), return; end
        cur = cur.(p{i});
    end
    v = cur;
end
