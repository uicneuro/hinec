classdef TestSeedClMask < matlab.unittest.TestCase
    % seeding.cl_min: seeds only where the voxel tensor is linear,
    % CL = (lambda_1 - lambda_2)/lambda_1 > cl_min (nim_seed_cl_mask).

    methods (TestClassSetup)
        function setup(tc) %#ok<MANU>
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, '..', '..', 'src', 'nim_tractography'));
            addpath(fullfile(here, '..', '..', 'src', 'nim_utils'));
        end
    end

    methods (Test)
        function keepsLinearDropsPancake(tc)
            ev = zeros(3, 1, 1, 3);
            ev(1, 1, 1, :) = [1.7 0.3 0.3];     % single bundle: CL = 0.824
            ev(2, 1, 1, :) = [1.0 0.9 0.3];     % crossing as a pancake: CL = 0.1
            ev(3, 1, 1, :) = [0 0 0];           % no fit: CL = 0
            [keep, cl] = nim_seed_cl_mask(struct('eval', ev), 0.15);
            tc.verifyEqual(cl(:)', [1.4/1.7, 0.1, 0], 'AbsTol', 1e-12);
            tc.verifyEqual(keep(:)', [true false false]);
        end

        function thresholdIsStrict(tc)
            ev = reshape([1 0.8 0.1], 1, 1, 1, 3);   % CL = 0.2 exactly
            tc.verifyFalse(nim_seed_cl_mask(struct('eval', ev), 0.2));
            tc.verifyTrue(nim_seed_cl_mask(struct('eval', ev), 0.19));
        end

        function needsEval(tc)
            tc.verifyError(@() nim_seed_cl_mask(struct('FA', 1), 0.1), 'nim_seed_cl_mask:noEval');
        end

        function schemaDefaultIsOff(tc)
            c = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti_hold.yml'));
            tc.verifyEqual(c.tractography.seeding.cl_min, 0);
            o = nim_config_to_options(c);
            tc.verifyEqual(o.seed_cl_min, 0);
        end
    end
end
