function nim_config_write(config, out_file, opts)
% nim_config_write: Write a canonical config as nested YAML.
%
%   nim_config_write(config, 'config/foo.yml')
%   nim_config_write(config, 'run/config.yml', struct('all', true, 'header', {{'line one'}}))
%
% By default only values that DIFFER from the schema default are written, which
% keeps configs minimal and makes the intent of each file obvious. Pass
% opts.all = true to emit every key (useful for run manifests, where a complete
% record matters more than brevity).
%
% opts.header : cell array of comment lines placed at the top of the file.

    if nargin < 3, opts = struct(); end
    emit_all = isfield(opts, 'all') && opts.all;
    header = {};
    if isfield(opts, 'header'), header = opts.header; end

    S = nim_config_schema();

    fid = fopen(out_file, 'w');
    if fid == -1
        error('nim_config_write:cannotOpen', 'Cannot write: %s', out_file);
    end
    cleaner = onCleanup(@() fclose(fid));

    for i = 1:numel(header)
        fprintf(fid, '# %s\n', header{i});
    end
    if ~isempty(header), fprintf(fid, '\n'); end

    for sec = {'preprocessing', 'tractography'}
        section = sec{1};
        lines = collect_section(config, S, section, emit_all);
        if isempty(lines), continue; end
        fprintf(fid, '%s:\n', section);
        for i = 1:numel(lines)
            fprintf(fid, '%s\n', lines{i});
        end
        fprintf(fid, '\n');
    end
end

% =========================================================================
function lines = collect_section(config, S, section, emit_all)
    lines = {};
    prefix = [section '.'];

    % Adaptive-only keys are inert unless the integrator is rkf45, and emitting
    % them would produce a config the loader rejects. Skip them.
    method = get_path(config, 'tractography.integrator.method', 'rk4');
    skip = {};
    if ~strcmp(method, 'rkf45')
        skip = {'tractography.integrator.tolerance', ...
                'tractography.integrator.step_min', ...
                'tractography.integrator.step_max', ...
                'tractography.integrator.safety', ...
                'tractography.integrator.adaptive'};
    end

    % Section-level keys first, then groups in schema order.
    order = {}; groups = {};
    for i = 1:numel(S)
        if ~startsWith(S(i).path, prefix), continue; end
        p = split(S(i).path, '.');
        if numel(p) == 2
            order{end+1} = S(i).path; %#ok<AGROW>
        else
            if ~any(strcmp(groups, p{2})), groups{end+1} = p{2}; end %#ok<AGROW>
        end
    end

    for i = 1:numel(order)
        e = entry_for(S, order{i});
        v = get_path(config, order{i}, e.default);
        if emit_all || ~is_default(v, e.default)
            lines{end+1} = sprintf('  %s: %s', leaf(order{i}), fmt(v)); %#ok<AGROW>
        end
    end

    for gi = 1:numel(groups)
        g = groups{gi};
        gl = {};
        for i = 1:numel(S)
            p = split(S(i).path, '.');
            if numel(p) ~= 3 || ~strcmp(p{1}, section) || ~strcmp(p{2}, g), continue; end
            if any(strcmp(skip, S(i).path)), continue; end
            e = S(i);
            v = get_path(config, S(i).path, e.default);
            if emit_all || ~is_default(v, e.default)
                gl{end+1} = sprintf('    %s: %s', p{3}, fmt(v)); %#ok<AGROW>
            end
        end
        if ~isempty(gl)
            lines{end+1} = sprintf('  %s:', g); %#ok<AGROW>
            lines = [lines, gl];
        end
    end
end

function e = entry_for(S, path)
    e = [];
    for i = 1:numel(S)
        if strcmp(S(i).path, path), e = S(i); return; end
    end
end

function s = leaf(path)
    p = split(path, '.');
    s = p{end};
end

function tf = is_default(v, d)
    try
        if iscell(v) && iscell(d), tf = isempty(v) && isempty(d); return; end
        if (ischar(v) || isstring(v)) && (ischar(d) || isstring(d))
            tf = strcmp(char(v), char(d)); return;
        end
        if islogical(v) || islogical(d), tf = (logical(v) == logical(d)); return; end
        if isnumeric(v) && isnumeric(d), tf = (v == d); return; end
        tf = isequal(v, d);
    catch
        tf = false;
    end
end

function s = fmt(v)
    if islogical(v)
        if v, s = 'true'; else, s = 'false'; end
    elseif isnumeric(v) && isscalar(v)
        if v == floor(v) && abs(v) < 1e15
            s = sprintf('%d', v);
        else
            s = strtrim(sprintf('%.10g', v));
        end
    elseif ischar(v) || isstring(v)
        s = char(v);
        if isempty(s) || any(isspace(s)) || any(s == ':' ) || any(s == '#')
            s = ['''' s ''''];
        end
    elseif iscell(v)
        parts = cell(1, numel(v));
        for i = 1:numel(v)
            parts{i} = fmt(v{i});
        end
        s = ['[' strjoin(parts, ', ') ']'];
    else
        s = '';
    end
end

function v = get_path(s, path, dflt)
    p = split(path, '.');
    v = dflt; cur = s;
    for i = 1:numel(p)
        if ~isstruct(cur) || ~isfield(cur, p{i}), return; end
        cur = cur.(p{i});
    end
    v = cur;
end
