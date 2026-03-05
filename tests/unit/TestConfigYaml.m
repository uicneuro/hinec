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
            presets = {'hinec_default.yml', 'fast_exploration.yml', ...
                       'high_precision.yml', 'standard_fact.yml', 'irontract.yml'};
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
            testCase.verifyTrue(isfield(cfg.tractography, 'step_size'));
            testCase.verifyTrue(isfield(cfg.tractography, 'seed_density'));
            testCase.verifyTrue(isfield(cfg.tractography, 'termination_fa'));
        end

        function testRequiredFields(testCase)
            % Tractography section must contain critical fields
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            required = {'algorithm', 'integration_order', 'step_size', ...
                        'seed_density', 'termination_fa', 'angle_thresh', ...
                        'max_steps', 'min_length'};
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
            testCase.verifyGreaterThan(cfg.tractography.step_size, 0);
            testCase.verifyGreaterThan(cfg.tractography.angle_thresh, 0);
            testCase.verifyLessThanOrEqual(cfg.tractography.angle_thresh, 180);
            testCase.verifyGreaterThanOrEqual(cfg.tractography.termination_fa, 0);
            testCase.verifyLessThanOrEqual(cfg.tractography.termination_fa, 1);
        end

        function testIntegrationOrderValid(testCase)
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            valid_orders = [1, 2, 4, 5];
            testCase.verifyTrue(ismember(cfg.tractography.integration_order, valid_orders));
        end

        function testPreprocessingDefaults(testCase)
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            testCase.verifyTrue(isfield(cfg.preprocessing, 'run_denoising'));
            testCase.verifyTrue(isfield(cfg.preprocessing, 'atlas_type'));
        end

        function testNumericParsing(testCase)
            % Ensure numeric values are parsed as numbers, not strings
            cfg = load_config_yaml(fullfile(testCase.ConfigDir, 'hinec_default.yml'));
            testCase.verifyClass(cfg.tractography.step_size, 'double');
            testCase.verifyClass(cfg.tractography.seed_density, 'double');
        end
    end
end
