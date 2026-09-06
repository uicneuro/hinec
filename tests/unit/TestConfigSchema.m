classdef TestConfigSchema < matlab.unittest.TestCase
    % Tests for the Phase 0 config schema rework: parser, migration,
    % validation, CLI overrides, and documentation drift.

    properties
        Root
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fullfile(here, '..', '..');
            addpath(fullfile(tc.Root, 'src', 'nim_utils'));
        end
    end

    methods (Test)

        % ---------------------------------------------------------- schema
        function schemaLoadsAndIsWellFormed(tc)
            S = nim_config_schema();
            tc.verifyGreaterThan(numel(S), 30);
            for i = 1:numel(S)
                tc.verifyNotEmpty(S(i).path);
                tc.verifyTrue(any(strcmp(S(i).type, {'numeric','string','logical','list'})), ...
                    sprintf('%s has an unknown type "%s"', S(i).path, S(i).type));
                tc.verifyNotEmpty(S(i).help, sprintf('%s has no help text', S(i).path));
                n = numel(split(S(i).path, '.'));
                tc.verifyTrue(n == 2 || n == 3, ...
                    sprintf('%s must be section.key or section.group.key', S(i).path));
            end
        end

        function legacyAliasesAreUnique(tc)
            % A duplicated alias would make migration ambiguous.
            S = nim_config_schema();
            all_leg = {};
            for i = 1:numel(S), all_leg = [all_leg, S(i).legacy]; end %#ok<AGROW>
            tc.verifyEqual(numel(unique(all_leg)), numel(all_leg), ...
                'Duplicate legacy alias across schema entries.');
        end

        function retiredKeysAreNotAlsoSchemaKeys(tc)
            S = nim_config_schema(); R = nim_config_retired();
            leaves = {};
            for i = 1:numel(S)
                p = split(S(i).path, '.'); leaves{end+1} = p{end}; %#ok<AGROW>
                leaves = [leaves, S(i).legacy]; %#ok<AGROW>
            end
            for i = 1:numel(R)
                tc.verifyFalse(any(strcmp(leaves, R(i).key)), ...
                    sprintf('"%s" is both retired and live in the schema.', R(i).key));
            end
        end


        function schemaDefaultsMatchTrackerDefaults(tc)
            % Regression guard. The schema now OWNS defaults for keys the old
            % loader never set, where the tracker's own `if ~isfield` default used
            % to apply. If the two disagree, a config that omits the key silently
            % changes behaviour. This caught sel_power (schema said 16, the hinec
            % tracker defaults to 0) - a config omitting it would have flipped from
            % plain interpolation to trajectory-dependent selection.
            S = nim_config_schema();
            expected = { ...
                'tractography.field',                   'dti';  % hinec.m:113
                'tractography.seeding.fa_min',          0.05;   % runTractography.m:162
                'tractography.csd.lmax',                6;      % runTractography.m:282
                'tractography.csd.n_iter',              50;
                'tractography.csd.peak_thresh',         0.5;
                'tractography.csd.peak_min_sep',        45;
                'tractography.csd.max_peaks',           3;
                'tractography.mmf.anchor',              0;      % mmf_connframe.m:29
                'tractography.mmf.frame_sel_power',     16};    % mmf_connframe.m:25
            for i = 1:size(expected, 1)
                path = expected{i,1}; want = expected{i,2};
                idx = find(strcmp({S.path}, path), 1);
                tc.assertNotEmpty(idx, sprintf('%s missing from schema', path));
                got = S(idx).default;
                if ischar(want)
                    tc.verifyEqual(char(got), want, ...
                        sprintf('%s default must match the tracker default', path));
                else
                    tc.verifyEqual(double(got), want, ...
                        sprintf('%s default must match the tracker default', path));
                end
            end
        end

        function selPowerIsFullyRemovedFromHinec(tc)
            % sel_power was a direction-steering term with no justification for
            % DTI. HINEC is interpolation + integration only; the key is retired
            % and the tracker must not reference it at all.
            S = nim_config_schema();
            tc.verifyEmpty(find(contains({S.path}, 'sel_power') & ...
                                ~contains({S.path}, 'frame_sel_power'), 1), ...
                'sel_power must not be in the schema.');
            R = nim_config_retired();
            tc.verifyTrue(any(strcmp({R.key}, 'sel_power')), ...
                'sel_power must be listed as retired.');
            src = fileread(fullfile(tc.Root, 'src', 'nim_tractography', ...
                                    'nim_tractography_hinec.m'));
            tc.verifyFalse(contains(src, 'sel_power'), ...
                'The hinec tracker still references sel_power.');
            tc.verifyFalse(contains(src, 'interp_e1_traj_h'), ...
                'The sel_power DTI interpolation path is still present.');
        end

        function mmfKeepsItsOwnFrameSelectivity(tc)
            % frame_sel_power is mmf's and is unaffected by removing sel_power.
            S = nim_config_schema();
            i = find(strcmp({S.path}, 'tractography.mmf.frame_sel_power'), 1);
            tc.assertNotEmpty(i);
            tc.verifyEqual(S(i).algos, {'mmf'});
            tc.verifyEqual(double(S(i).default), 16);
        end

        % ---------------------------------------------------------- parser
        function allShippedConfigsLoad(tc)
            w = warning('off','all'); c = onCleanup(@() warning(w));
            f = dir(fullfile(tc.Root, 'config', '*.yml'));
            tc.verifyGreaterThan(numel(f), 0);
            for i = 1:numel(f)
                cfg = load_config_yaml(fullfile(tc.Root, 'config', f(i).name));
                tc.verifyTrue(isfield(cfg, 'tractography'));
                tc.verifyTrue(isfield(cfg.tractography, 'integrator'));
            end
        end

        function nestingDeeperThanTwoLevelsIsRejected(tc)
            p = [tempname '.yml'];
            fid = fopen(p,'w');
            fprintf(fid, 'tractography:\n  integrator:\n    sub:\n      key: 1\n');
            fclose(fid);
            c = onCleanup(@() delete(p));
            tc.verifyError(@() nim_yaml_parse(p), 'nim_yaml_parse:tooDeep');
        end

        function inlineListsParse(tc)
            p = [tempname '.yml'];
            fid = fopen(p,'w');
            fprintf(fid, 'tractography:\n  seeding:\n    roi: [41, 42, ''Uncinate fasciculus R'']\n');
            fclose(fid);
            c = onCleanup(@() delete(p));
            r = nim_yaml_parse(p);
            tc.verifyEqual(numel(r.tractography.seeding.roi), 3);
            tc.verifyEqual(r.tractography.seeding.roi{1}, 41);
            tc.verifyEqual(r.tractography.seeding.roi{3}, 'Uncinate fasciculus R');
        end

        % ------------------------------------------------------ validation
        function unknownKeyIsRejected(tc)
            p = [tempname '.yml'];
            fid = fopen(p,'w'); fprintf(fid, 'tractography:\n  nonsense_key: 1\n'); fclose(fid);
            c = onCleanup(@() delete(p));
            tc.verifyError(@() load_config_yaml(p), 'load_config_yaml:unknownKey');
        end

        function adaptiveKeyWithFixedStepMethodIsRejected(tc)
            p = [tempname '.yml'];
            fid = fopen(p,'w');
            fprintf(fid, 'tractography:\n  integrator:\n    method: rk4\n    tolerance: 0.01\n');
            fclose(fid);
            c = onCleanup(@() delete(p));
            tc.verifyError(@() load_config_yaml(p), 'load_config_yaml:adaptiveWithFixedStep');
        end

        function outOfRangeValueIsRejected(tc)
            p = [tempname '.yml'];
            fid = fopen(p,'w'); fprintf(fid, 'tractography:\n  seeding:\n    density: -3\n'); fclose(fid);
            c = onCleanup(@() delete(p));
            tc.verifyError(@() load_config_yaml(p), 'load_config_yaml:outOfRange');
        end

        % ------------------------------------------------------- migration
        function legacyIntegrationOrderMigrates(tc)
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            p = [tempname '.yml'];
            fid = fopen(p,'w');
            fprintf(fid, 'tractography:\n  integration_order: 5\n  adaptive_step: true\n  step_size: 0.2\n');
            fclose(fid);
            c = onCleanup(@() delete(p));
            cfg = load_config_yaml(p);
            tc.verifyEqual(cfg.tractography.integrator.method, 'rkf45');
            o = nim_config_to_options(cfg);
            tc.verifyEqual(o.integration_order, 5);
            tc.verifyTrue(o.adaptive_step);
        end

        function maxStepsMigratesToArcLengthAndRoundTrips(tc)
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            p = [tempname '.yml'];
            fid = fopen(p,'w');
            fprintf(fid, 'tractography:\n  step_size: 0.2\n  max_steps: 3000\n');
            fclose(fid);
            c = onCleanup(@() delete(p));
            cfg = load_config_yaml(p);
            tc.verifyEqual(cfg.tractography.termination.max_arc, 600, 'AbsTol', 1e-9);
            o = nim_config_to_options(cfg);
            tc.verifyEqual(o.max_steps, 3000);   % exact round trip
        end

        function maxStepsIsDerivedFromStepSoRefiningDoesNotTruncate(tc)
            % The failure mode this rework exists to remove: halving the step
            % must NOT halve how far a track may travel.
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            cfg = load_config_yaml(fullfile(tc.Root, 'config', 'ismrm2015.yml'));
            arc = cfg.tractography.termination.max_arc;
            steps = [0.5 0.25 0.125 0.0625];
            for s = steps
                c2 = nim_config_apply_overrides(cfg, {sprintf('integrator.step=%g', s)});
                o = nim_config_to_options(c2);
                tc.verifyEqual(o.max_steps, ceil(arc / s));
                tc.verifyEqual(o.max_arc, arc, 'AbsTol', 1e-9);
            end
        end

        % -------------------------------------------------------- overrides
        function overridesReachEveryScope(tc)
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            cfg = load_config_yaml(fullfile(tc.Root, 'config', 'ismrm2015.yml'));
            c2 = nim_config_apply_overrides(cfg, { ...
                'tractography.integrator.step=0.05', ...
                'interpolation.method=cubic', ...
                'interpolation.upsample=2', ...
                'preprocessing.run_eddy=false', ...
                'seeding.roi=[41,42]'});
            tc.verifyEqual(c2.tractography.integrator.step, 0.05);
            tc.verifyEqual(c2.tractography.interpolation.method, 'cubic');
            tc.verifyEqual(c2.tractography.interpolation.upsample, 2);
            tc.verifyFalse(c2.preprocessing.run_eddy);
            tc.verifyEqual(numel(c2.tractography.seeding.roi), 2);
        end

        function ambiguousBareLeafIsRejected(tc)
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            cfg = load_config_yaml(fullfile(tc.Root, 'config', 'ismrm2015.yml'));
            tc.verifyError(@() nim_config_apply_overrides(cfg, {'method=rk4'}), ...
                'nim_config_apply_overrides:ambiguousKey');
            tc.verifyError(@() nim_config_apply_overrides(cfg, {'fa_min=0.2'}), ...
                'nim_config_apply_overrides:ambiguousKey');
        end

        function unknownOverrideIsRejected(tc)
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            cfg = load_config_yaml(fullfile(tc.Root, 'config', 'ismrm2015.yml'));
            tc.verifyError(@() nim_config_apply_overrides(cfg, {'not_a_param=1'}), ...
                'nim_config_apply_overrides:unknownKey');
        end

        % ------------------------------------------------------------ docs
        function documentationMatchesSchema(tc)
            % Guards the drift that made the old docs describe integration_order
            % while the code had moved on, and omit sel_power and mmf entirely.
            doc_path = fullfile(tc.Root, 'docs', 'YAML_CONFIG.md');
            tc.assertTrue(isfile(doc_path), 'docs/YAML_CONFIG.md is missing.');
            on_disk = fileread(doc_path);
            generated = nim_config_docs();
            tc.verifyEqual(strtrim(strrep(on_disk, sprintf('\r'), '')), ...
                           strtrim(strrep(generated, sprintf('\r'), '')), ...
                ['docs/YAML_CONFIG.md is out of date. Regenerate with ' ...
                 'nim_config_docs(''docs/YAML_CONFIG.md'').']);
        end

        function everySchemaKeyAppearsInTheDocs(tc)
            md = nim_config_docs();
            S = nim_config_schema();
            for i = 1:numel(S)
                p = split(S(i).path, '.');
                tc.verifyTrue(contains(md, sprintf('`%s`', p{end})), ...
                    sprintf('%s is not documented.', S(i).path));
            end
        end

        % ------------------------------------------------- round-trip write
        function writeThenReadIsStable(tc)
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            f = dir(fullfile(tc.Root, 'config', '*.yml'));
            for i = 1:numel(f)
                src = fullfile(tc.Root, 'config', f(i).name);
                c1 = load_config_yaml(src);
                p = [tempname '.yml'];
                nim_config_write(c1, p);
                c2 = load_config_yaml(p);
                delete(p);
                o1 = nim_config_to_options(c1);
                o2 = nim_config_to_options(c2);
                tc.verifyEqual(o2, o1, sprintf('%s does not survive a write/read round trip.', f(i).name));
            end
        end
    end
end
