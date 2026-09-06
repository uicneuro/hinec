classdef TestConfigYaml < matlab.unittest.TestCase
% TestConfigYaml: Unit tests for load_config_yaml

    properties
        ConfigDir
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
            testCase.ConfigDir = fullfile(projRoot, 'config');
        end
    end

    methods (Test)
        function testAllPresetsLoad(testCase)
            % All 5 YAML config files should load without error
            presets = {'hinec_default.yml', 'hinec_dti.yml', ...
                       'hinec_csd.yml', 'standard_dti.yml', 'irontract.yml'};
            for i = 1:numel(presets)
                cfg = load_config_yaml(fullfile(testCase.ConfigDir, presets{i}));
                testCase.verifyTrue(isstruct(cfg), ...
                    sprintf('%s should load as struct', presets{i}));
            end
        end

        function testDefaultsApplied(testCase)
            % Default config should have all expected fields
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            testCase.verifyTrue(isfield(cfg, 'preprocessing'));
            testCase.verifyTrue(isfield(cfg, 'tractography'));
            testCase.verifyTrue(isfield(cfg.tractography, 'integrator'));
            testCase.verifyTrue(isfield(cfg.tractography.integrator, 'step'));
            testCase.verifyTrue(isfield(cfg.tractography.seeding, 'density'));
            testCase.verifyTrue(isfield(cfg.tractography.termination, 'fa_min'));
        end

        function testRequiredFields(testCase)
            % Tractography section must contain critical fields
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            % Canonical groups, not the old flat keys.
            required = {'algorithm', 'field', 'act', 'integrator', ...
                        'interpolation', 'seeding', 'termination', 'filter'};
            for i = 1:numel(required)
                testCase.verifyTrue(isfield(cfg.tractography, required{i}), ...
                    sprintf('Missing required field: %s', required{i}));
            end
        end

        function testInvalidFileThrows(testCase)
            threw = false;
            try
                load_config_yaml('/nonexistent/file.yml');
            catch
                threw = true;
            end
            testCase.verifyTrue(threw, 'Should throw error for nonexistent file');
        end

        function testParameterRanges(testCase)
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            testCase.verifyGreaterThan(cfg.tractography.integrator.step, 0);
            testCase.verifyGreaterThan(cfg.tractography.termination.angle_max, 0);
            testCase.verifyLessThanOrEqual(cfg.tractography.termination.angle_max, 180);
            testCase.verifyGreaterThanOrEqual(cfg.tractography.termination.fa_min, 0);
            testCase.verifyLessThanOrEqual(cfg.tractography.termination.fa_min, 1);
        end

        function testIntegrationOrderValid(testCase)
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            valid_methods = {'euler', 'rk2', 'rk4', 'rkf45'};
            testCase.verifyTrue(ismember(cfg.tractography.integrator.method, valid_methods));
        end

        function testPreprocessingDefaults(testCase)
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            testCase.verifyTrue(isfield(cfg.preprocessing, 'run_denoising'));
            testCase.verifyTrue(isfield(cfg.preprocessing, 'atlas_type'));
        end

        function testNumericParsing(testCase)
            % Ensure numeric values are parsed as numbers, not strings
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            testCase.verifyClass(cfg.tractography.integrator.step, 'double');
            testCase.verifyClass(cfg.tractography.seeding.density, 'double');
        end
    end
end
