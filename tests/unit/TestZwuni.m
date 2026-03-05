classdef TestZwuni < matlab.unittest.TestCase
% TestZwuni: Unit tests for zwuni (uniform nodes and weights on [-1,1])

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        function testEndpoints(testCase)
            [z, ~] = zwuni(4);
            testCase.verifyEqual(z(1), -1, 'AbsTol', 1e-14);
            testCase.verifyEqual(z(end), 1, 'AbsTol', 1e-14);
        end

        function testNodeCount(testCase)
            % zwuni(N) returns N+1 nodes
            for N = [2, 5, 10]
                [z, w] = zwuni(N);
                testCase.verifyLength(z, N + 1);
                testCase.verifyLength(w, N + 1);
            end
        end

        function testUniformSpacing(testCase)
            % Spacing between consecutive nodes should be constant
            [z, ~] = zwuni(8);
            diffs = diff(z);
            testCase.verifyEqual(diffs, repmat(diffs(1), size(diffs)), 'AbsTol', 1e-14);
        end

        function testWeightsSum(testCase)
            % Weights must sum to 2 (interval length)
            for N = [3, 6, 12]
                [~, w] = zwuni(N);
                testCase.verifyEqual(sum(w), 2, 'AbsTol', 1e-14);
            end
        end

        function testHalfEndpointWeights(testCase)
            % Endpoint weights are half the interior weights
            [~, w] = zwuni(4);
            interior_w = w(2);
            testCase.verifyEqual(w(1), 0.5 * interior_w, 'AbsTol', 1e-14);
            testCase.verifyEqual(w(end), 0.5 * interior_w, 'AbsTol', 1e-14);
        end

        function testColumnVectors(testCase)
            [z, w] = zwuni(5);
            testCase.verifyTrue(iscolumn(z), 'z must be a column vector');
            testCase.verifyTrue(iscolumn(w), 'w must be a column vector');
        end
    end
end
