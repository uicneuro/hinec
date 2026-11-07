classdef TestSeedCoherenceMask < matlab.unittest.TestCase
    % seeding.coherence_min: seeds only where the principal direction agrees with
    % its 3 x 3 x 3 neighbours (nim_seed_coherence_mask).

    methods (TestClassSetup)
        function setup(tc) %#ok<MANU>
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, '..', '..', 'src', 'nim_tractography'));
            addpath(fullfile(here, '..', '..', 'src', 'nim_utils'));
        end
    end

    methods (Test)
        function uniformFieldIsOne(tc)
            nim = field(repmat(reshape([0 0 1], 1, 1, 1, 3), 3, 3, 3));
            [~, coh] = nim_seed_coherence_mask(nim, 0.5);
            tc.verifyEqual(coh, ones(3, 3, 3), 'AbsTol', 1e-12);
        end

        function signDoesNotMatter(tc)
            v = repmat(reshape([1 0 0], 1, 1, 1, 3), 3, 3, 3);
            v(1:2:end, :, :, :) = -v(1:2:end, :, :, :);
            [~, coh] = nim_seed_coherence_mask(field(v), 0.5);
            tc.verifyEqual(coh, ones(3, 3, 3), 'AbsTol', 1e-12);
        end

        function oddVoxelOutIsLowAndDropped(tc)
            % centre points along z, every neighbour along x: M = diag(26, 0, 1)/27
            v = repmat(reshape([1 0 0], 1, 1, 1, 3), 3, 3, 3);
            v(2, 2, 2, :) = [0 0 1];
            [keep, coh] = nim_seed_coherence_mask(field(v), 0.9);
            tc.verifyEqual(coh(2, 2, 2), 26 / 27, 'AbsTol', 1e-12);
            tc.verifyTrue(keep(2, 2, 2));
            % three orthogonal directions mixed equally around a voxel: coherence 1/3
            v = zeros(3, 1, 1, 3);
            v(1, 1, 1, :) = [1 0 0]; v(2, 1, 1, :) = [0 1 0]; v(3, 1, 1, :) = [0 0 1];
            [keep, coh] = nim_seed_coherence_mask(field(v), 0.5);
            tc.verifyEqual(coh(2), 1 / 3, 'AbsTol', 1e-12);
            tc.verifyFalse(keep(2));
        end

        function noTensorIsDropped(tc)
            nim = field(repmat(reshape([0 1 0], 1, 1, 1, 3), 3, 1, 1));
            nim.eval(2, 1, 1, :) = 0;
            [keep, coh] = nim_seed_coherence_mask(nim, 0.1);
            tc.verifyEqual(coh(2), 0);
            tc.verifyEqual(keep(:)', [true false true]);
        end

        function needsEig(tc)
            tc.verifyError(@() nim_seed_coherence_mask(struct('FA', 1), 0.5), 'nim_seed_coherence_mask:noEig');
        end

        function schemaDefaultIsOff(tc)
            c = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti_hold.yml'));
            tc.verifyEqual(c.tractography.seeding.coherence_min, 0);
            tc.verifyEqual(nim_config_to_options(c).seed_coherence_min, 0);
        end
    end
end

function nim = field(v1)
sz = size(v1);
nim.evec = zeros([sz(1:3) 3 3]);
nim.evec(:, :, :, :, 1) = v1;
nim.eval = repmat(reshape([1.5 0.3 0.2], 1, 1, 1, 3), sz(1), sz(2), sz(3));
end
