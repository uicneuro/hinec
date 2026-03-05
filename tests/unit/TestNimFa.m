classdef TestNimFa < matlab.unittest.TestCase
% TestNimFa: Unit tests for nim_fa (fractional anisotropy calculation)

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        function testIsotropicFA(testCase)
            % Equal eigenvalues -> FA = 0
            nim = makeMinimalNim(2);
            nim.eval = repmat(reshape([1 1 1], 1, 1, 1, 3), [2 2 2 1]);
            nim.mask = ones(2, 2, 2);
            result = nim_fa(nim);
            testCase.verifyEqual(result.FA(1,1,1), 0, 'AbsTol', 1e-10);
        end

        function testHighAnisotropy(testCase)
            % One dominant eigenvalue -> FA close to 1
            nim = makeMinimalNim(2);
            nim.eval = repmat(reshape([1 0 0], 1, 1, 1, 3), [2 2 2 1]);
            nim.mask = ones(2, 2, 2);
            result = nim_fa(nim);
            testCase.verifyEqual(result.FA(1,1,1), 1, 'AbsTol', 1e-10);
        end

        function testKnownFA(testCase)
            % lambda = [0.0017, 0.0003, 0.0003] -> FA ~ 0.8
            nim = makeMinimalNim(2);
            nim.eval = repmat(reshape([0.0017 0.0003 0.0003], 1, 1, 1, 3), [2 2 2 1]);
            nim.mask = ones(2, 2, 2);
            result = nim_fa(nim);
            testCase.verifyEqual(result.FA(1,1,1), 0.8, 'AbsTol', 0.05);
        end

        function testMaskRespected(testCase)
            % FA should be zero outside the mask
            nim = makeMinimalNim(2);
            nim.eval = repmat(reshape([0.0017 0.0003 0.0003], 1, 1, 1, 3), [2 2 2 1]);
            nim.mask = zeros(2, 2, 2);
            nim.mask(1, 1, 1) = 1;
            result = nim_fa(nim);
            testCase.verifyGreaterThan(result.FA(1,1,1), 0);
            testCase.verifyEqual(result.FA(2,2,2), 0);
        end

        function testZeroEigenvalues(testCase)
            % All-zero eigenvalues should not cause divide-by-zero
            nim = makeMinimalNim(2);
            nim.eval = zeros(2, 2, 2, 3);
            nim.mask = ones(2, 2, 2);
            result = nim_fa(nim);
            testCase.verifyEqual(result.FA(1,1,1), 0);
        end

        function testFARange(testCase)
            % FA must always be in [0, 1]
            nim = makeMinimalNim(2);
            nim.eval = repmat(reshape([3 0.5 0.1], 1, 1, 1, 3), [2 2 2 1]);
            nim.mask = ones(2, 2, 2);
            result = nim_fa(nim);
            testCase.verifyGreaterThanOrEqual(result.FA(:), 0);
            testCase.verifyLessThanOrEqual(result.FA(:), 1);
        end

        function testOutputDimensions(testCase)
            nim = makeMinimalNim(4);
            nim.eval = zeros(4, 4, 4, 3);
            nim.mask = ones(4, 4, 4);
            result = nim_fa(nim);
            testCase.verifySize(result.FA, [4 4 4]);
        end
    end
end

function nim = makeMinimalNim(dim)
    nim.hdr.ImageSize = [dim dim dim 7];
    nim.xdim = dim;
    nim.ydim = dim;
    nim.zdim = dim;
end
