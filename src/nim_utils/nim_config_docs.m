function md = nim_config_docs(out_file)
% nim_config_docs: Generate the config reference from the schema.
%
%   nim_config_docs('docs/YAML_CONFIG.md')   % write
%   md = nim_config_docs();                  % return the markdown
%
% The reference is GENERATED, so it cannot drift from the code. Every documented
% key exists in nim_config_schema, and every schema key is documented, by
% construction. tests/unit/TestConfigDocs.m regenerates and diffs to enforce it.
%
% Covers all three tractography algorithms in one table, marking which of them
% actually reads each parameter - previously the docs described only hinec, under
% the wrong key names.

    S = nim_config_schema();
    R = nim_config_retired();
    L = {};

    a = @(varargin) assignin('caller', 'L', [evalin('caller','L'), {sprintf(varargin{:})}]);

    L{end+1} = '# HINEC Configuration Reference';
    L{end+1} = '';
    L{end+1} = '<!-- GENERATED FILE - do not edit by hand.';
    L{end+1} = '     Regenerate with:  nim_config_docs(''docs/YAML_CONFIG.md'')';
    L{end+1} = '     Source of truth:  src/nim_utils/nim_config_schema.m -->';
    L{end+1} = '';
    L{end+1} = 'Configuration is YAML with **two levels of nesting** below a section';
    L{end+1} = '(`section` -> `group` -> `key`). A third level is a parse error, deliberately:';
    L{end+1} = 'deeper nesting was judged more confusing than helpful.';
    L{end+1} = '';
    L{end+1} = 'Every key is optional. Anything you omit takes the default below, so a working';
    L{end+1} = 'config can be very short.';
    L{end+1} = '';
    L{end+1} = '## Minimal example';
    L{end+1} = '';
    L{end+1} = '```yaml';
    L{end+1} = 'tractography:';
    L{end+1} = '  algorithm: hinec';
    L{end+1} = '  seeding:';
    L{end+1} = '    roi: [SLF_L, SLF_R]';
    L{end+1} = '```';
    L{end+1} = '';
    L{end+1} = '## Overriding from the command line';
    L{end+1} = '';
    L{end+1} = 'Every parameter in this reference can be overridden with `--set`, on both';
    L{end+1} = '`bin/run_hinec.sh` and `bin/run_tractography.sh`:';
    L{end+1} = '';
    L{end+1} = '```bash';
    L{end+1} = './bin/run_tractography.sh hinec_dti --set tractography.integrator.step=0.05';
    L{end+1} = './bin/run_tractography.sh hinec_dti --set integrator.step=0.05    # section assumed';
    L{end+1} = './bin/run_tractography.sh hinec_dti --set upsample=2              # bare leaf, if unique';
    L{end+1} = './bin/run_hinec.sh data/x x.mat config/y.yml --set preprocessing.run_eddy=false';
    L{end+1} = './bin/run_tractography.sh hinec_dti --set seeding.roi=[41,42]     # lists';
    L{end+1} = '```';
    L{end+1} = '';
    L{end+1} = 'Overrides are checked against the schema. An unknown or mistyped key is an';
    L{end+1} = 'error, never a silent no-op. A bare leaf name that is ambiguous (`method`,';
    L{end+1} = '`fa_min`) is rejected with the candidate full paths.';
    L{end+1} = '';
    L{end+1} = '## Parameters';
    L{end+1} = '';
    L{end+1} = 'The **Applies to** column says which tracking algorithms actually read the';
    L{end+1} = 'parameter. A key marked `hinec` is ignored by `standard` and `mmf`.';
    L{end+1} = '';

    % ---- group the schema by section, then group ---------------------------
    sections = {'tractography', 'preprocessing'};
    for si = 1:numel(sections)
        sec = sections{si};
        L{end+1} = sprintf('### `%s`', sec);
        L{end+1} = '';

        % section-level keys
        rows = {};
        for i = 1:numel(S)
            p = split(S(i).path, '.');
            if numel(p) == 2 && strcmp(p{1}, sec)
                rows{end+1} = row_for(S(i)); %#ok<AGROW>
            end
        end
        if ~isempty(rows)
            L = [L, table_head(), rows];
            L{end+1} = '';
        end

        % groups
        seen = {};
        for i = 1:numel(S)
            p = split(S(i).path, '.');
            if numel(p) ~= 3 || ~strcmp(p{1}, sec), continue; end
            if any(strcmp(seen, p{2})), continue; end
            seen{end+1} = p{2}; %#ok<AGROW>

            L{end+1} = sprintf('#### `%s.%s`', sec, p{2});
            L{end+1} = '';
            grows = {};
            for j = 1:numel(S)
                q = split(S(j).path, '.');
                if numel(q) == 3 && strcmp(q{1}, sec) && strcmp(q{2}, p{2})
                    grows{end+1} = row_for(S(j)); %#ok<AGROW>
                end
            end
            L = [L, table_head(), grows];
            L{end+1} = '';
        end
    end

    % ---- migration ---------------------------------------------------------
    L{end+1} = '## Migrating from the old flat config';
    L{end+1} = '';
    L{end+1} = 'Old configs still load; each superseded key produces a deprecation warning';
    L{end+1} = 'naming its replacement.';
    L{end+1} = '';
    L{end+1} = '| Old key | Replacement |';
    L{end+1} = '|---|---|';
    for i = 1:numel(S)
        for j = 1:numel(S(i).legacy)
            L{end+1} = sprintf('| `%s` | `%s` |', S(i).legacy{j}, S(i).path); %#ok<AGROW>
        end
    end
    L{end+1} = '| `integration_order: 1\|2\|4\|5` | `integrator.method: euler\|rk2\|rk4\|rkf45` |';
    L{end+1} = '| `max_steps` | `termination.max_arc` (converted as `max_steps x step`) |';
    L{end+1} = '';
    L{end+1} = '`integration_order` was a **method selector wearing a numeric-order name**:';
    L{end+1} = 'the value `5` selected RKF45. That value was not wrong numerically - the';
    L{end+1} = 'implementation uses Dormand-Prince coefficients and advances on the 5th-order';
    L{end+1} = 'solution, using the embedded 4th-order one only for error control - but a';
    L{end+1} = 'method selector should not be spelled as a number. `integrator.method` names';
    L{end+1} = 'the method instead.';
    L{end+1} = '';
    L{end+1} = '`max_steps` counted integration steps, so halving the step size silently';
    L{end+1} = 'halved how far a track could travel. `termination.max_arc` is an arc length';
    L{end+1} = 'in voxels and is step-invariant; `max_steps` is now derived as';
    L{end+1} = '`ceil(max_arc / step)`.';
    L{end+1} = '';

    % ---- retired -----------------------------------------------------------
    L{end+1} = '## Retired keys';
    L{end+1} = '';
    L{end+1} = 'These were accepted by the old loader but read by no tracker. They are';
    L{end+1} = 'ignored with a warning rather than silently accepted.';
    L{end+1} = '';
    L{end+1} = '| Key | Why it is gone |';
    L{end+1} = '|---|---|';
    for i = 1:numel(R)
        L{end+1} = sprintf('| `%s` | %s |', R(i).key, one_line(R(i).reason)); %#ok<AGROW>
    end
    L{end+1} = '';

    md = strjoin(L, newline);
    if nargin >= 1 && ~isempty(out_file)
        fid = fopen(out_file, 'w');
        if fid == -1, error('nim_config_docs:cannotWrite', 'Cannot write %s', out_file); end
        fwrite(fid, md); fclose(fid);
        fprintf('Wrote %s (%d lines)\n', out_file, numel(L));
    end
end

% =========================================================================
function h = table_head()
    h = {'| Key | Type | Default | Applies to | Description |', '|---|---|---|---|---|'};
end

function s = row_for(e)
    p = split(e.path, '.');
    s = sprintf('| `%s` | %s | `%s` | %s | %s |', ...
        p{end}, e.type, fmt_default(e.default), fmt_algos(e.algos), one_line(e.help));
end

function s = fmt_algos(algos)
    if isempty(algos), s = '-'; return; end
    if any(strcmp(algos, 'all')), s = 'all'; return; end
    if any(strcmp(algos, 'n/a')), s = '-'; return; end
    s = strjoin(algos, ', ');
end

function s = fmt_default(d)
    if islogical(d)
        if d, s = 'true'; else, s = 'false'; end
    elseif isnumeric(d) && isscalar(d)
        if d == floor(d), s = sprintf('%d', d); else, s = strtrim(sprintf('%.10g', d)); end
    elseif ischar(d) || isstring(d)
        s = char(d);
    elseif iscell(d)
        if isempty(d), s = '[]'; else, s = '[...]'; end
    else
        s = '';
    end
end

function s = one_line(t)
    s = regexprep(char(t), '\s+', ' ');
    s = strrep(s, '|', '\|');
end
