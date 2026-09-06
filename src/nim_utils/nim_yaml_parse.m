function raw = nim_yaml_parse(filename)
% nim_yaml_parse: Indentation-aware YAML reader for the HINEC config subset.
%
% Replaces the previous two-level, indentation-BLIND parser, which treated any
% line ending in ':' as a new top-level section regardless of indentation. That
% silently misparsed nested YAML instead of rejecting it.
%
% Supports exactly:
%   section:              (column 0)
%     key: value          (one level in)
%     group:              (one level in, no value)
%       key: value        (two levels in)
%
% A third level of nesting is a hard error - nesting is deliberately kept mild.
%
% Values: numbers, booleans, null/~, quoted and unquoted strings, and inline
% lists [a, b, c]. Inline '#' comments are stripped outside quotes.
%
% Returns a nested struct mirroring the file. No defaults, no validation - see
% load_config_yaml for those.

    fid = fopen(filename, 'r');
    if fid == -1
        error('nim_yaml_parse:cannotOpen', 'Cannot open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));

    raw = struct();
    stack = {};       % stack{d} = field name at depth d
    indents = [];     % indents(d) = column at which depth d starts
    line_num = 0;

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line_num = line_num + 1;

        stripped = strip_comment(line);
        if isempty(strtrim(stripped)), continue; end

        indent = numel(stripped) - numel(strtrim_left(stripped));
        content = strtrim(stripped);

        % Determine depth from the indent stack.
        while ~isempty(indents) && indent < indents(end)
            indents(end) = []; %#ok<AGROW>
            stack(end)   = []; %#ok<AGROW>
        end
        if isempty(indents)
            if indent ~= 0
                error('nim_yaml_parse:badIndent', ...
                    '%s line %d: unexpected indentation (expected a top-level section).', ...
                    filename, line_num);
            end
            depth = 0;
        elseif indent == indents(end)
            depth = numel(indents) - 1;
            stack(end) = []; %#ok<AGROW>
            indents(end) = []; %#ok<AGROW>
        else % indent > indents(end)
            depth = numel(indents);
        end

        if depth > 2
            error('nim_yaml_parse:tooDeep', ...
                ['%s line %d: nesting deeper than section.group.key is not supported ' ...
                 '(config nesting is deliberately limited to two levels below a section).'], ...
                filename, line_num);
        end

        colon = find(content == ':', 1);
        if isempty(colon)
            error('nim_yaml_parse:noColon', '%s line %d: expected "key: value" or "group:", got "%s".', ...
                filename, line_num, content);
        end

        key = strtrim(content(1:colon-1));
        val_str = strtrim(content(colon+1:end));

        if ~isvarname(key)
            error('nim_yaml_parse:badKey', '%s line %d: "%s" is not a valid key name.', ...
                filename, line_num, key);
        end

        stack{end+1} = key;      %#ok<AGROW>
        indents(end+1) = indent; %#ok<AGROW>

        if isempty(val_str)
            % A group / section header. Create it so empty groups still exist.
            raw = assign_path(raw, stack, struct(), true);
        else
            if depth == 0
                error('nim_yaml_parse:valueAtTop', ...
                    '%s line %d: "%s" is at the top level but has a value; only sections may appear there.', ...
                    filename, line_num, key);
            end
            raw = assign_path(raw, stack, parse_value(val_str, filename, line_num), false);
            % A leaf never has children; pop it so siblings resolve correctly.
            stack(end) = [];   %#ok<AGROW>
            indents(end) = []; %#ok<AGROW>
        end
    end
end

% -------------------------------------------------------------------------
function s = strtrim_left(s)
    idx = find(~isspace(s), 1);
    if isempty(idx), s = ''; else, s = s(idx:end); end
end

function out = strip_comment(line)
% Remove an inline # comment, respecting quotes and bracketed lists.
    out = line;
    in_s = false; in_d = false;
    for i = 1:numel(line)
        c = line(i);
        if c == '''' && ~in_d, in_s = ~in_s;
        elseif c == '"' && ~in_s, in_d = ~in_d;
        elseif c == '#' && ~in_s && ~in_d
            out = line(1:i-1);
            return;
        end
    end
end

function s = assign_path(s, path, value, is_group)
% Set s.(path{1}).(path{2})... = value. Groups never clobber existing content.
    if numel(path) == 1
        f = path{1};
        if is_group && isfield(s, f) && isstruct(s.(f)), return; end
        s.(f) = value;
    else
        f = path{1};
        if ~isfield(s, f) || ~isstruct(s.(f)), s.(f) = struct(); end
        s.(f) = assign_path(s.(f), path(2:end), value, is_group);
    end
end

function value = parse_value(str, filename, line_num)
% Scalar or inline-list value.
    if startsWith(str, '[')
        if ~endsWith(str, ']')
            error('nim_yaml_parse:badList', '%s line %d: unterminated list "%s".', ...
                filename, line_num, str);
        end
        inner = strtrim(str(2:end-1));
        if isempty(inner), value = {}; return; end
        parts = split_list(inner);
        value = cell(1, numel(parts));
        for i = 1:numel(parts)
            value{i} = parse_scalar(strtrim(parts{i}));
        end
        return;
    end
    value = parse_scalar(str);
end

function parts = split_list(inner)
% Split on commas that are outside quotes.
    parts = {}; cur = ''; in_s = false; in_d = false;
    for i = 1:numel(inner)
        c = inner(i);
        if c == '''' && ~in_d, in_s = ~in_s; cur(end+1) = c; %#ok<AGROW>
        elseif c == '"' && ~in_s, in_d = ~in_d; cur(end+1) = c; %#ok<AGROW>
        elseif c == ',' && ~in_s && ~in_d
            parts{end+1} = cur; cur = ''; %#ok<AGROW>
        else
            cur(end+1) = c; %#ok<AGROW>
        end
    end
    parts{end+1} = cur;
end

function value = parse_scalar(str)
    if (startsWith(str, '''') && endsWith(str, '''') && numel(str) >= 2) || ...
       (startsWith(str, '"')  && endsWith(str, '"')  && numel(str) >= 2)
        value = str(2:end-1);
        return;
    end
    if strcmpi(str, 'true'),  value = true;  return; end
    if strcmpi(str, 'false'), value = false; return; end
    if strcmpi(str, 'null') || strcmp(str, '~'), value = []; return; end
    num = str2double(str);
    if ~isnan(num) && ~isempty(regexp(str, '^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$', 'once'))
        value = num;
        return;
    end
    value = str;
end
