classdef TestZwgll < matlab.unittest.TestCase
% TestZwgll: Unit tests for zwgll (Gauss-Lobatto-Legendre nodes/weights)

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        function testEndpoints(testCase)
            % First and last nodes must be -1 and 1
            for p = [2, 3, 5, 8]
                [z, ~] = zwgll(p);
                testCase.verifyEqual(z(1), -1, 'AbsTol', 1e-14);
                testCase.verifyEqual(z(end), 1, 'AbsTol', 1e-14);
            end
        end

        function testNodeCount(testCase)
            % zwgll(p) returns p+1 nodes
            for p = [1, 2, 4, 7]
                [z, w] = zwgll(p);
                testCase.verifyLength(z, p + 1);
                testCase.verifyLength(w, p + 1);
            end
        end

        function testSorted(testCase)
            % Nodes must be strictly increasing
            [z, ~] = zwgll(6);
            testCase.verifyTrue(all(diff(z) > 0), 'Nodes must be sorted ascending');
        end

        function testWeightsPositive(testCase)
            % All weights must be positive
            [~, w] = zwgll(5);
            testCase.verifyTrue(all(w > 0), 'All weights must be positive');
        end

        function testWeightsSum(testCase)
            % Weights on [-1,1] must sum to 2
            for p = [2, 4, 6, 10]
                [~, w] = zwgll(p);
                testCase.verifyEqual(sum(w), 2, 'AbsTol', 1e-12);
            end
        end

        function testP1(testCase)
            % p=1: two nodes at -1 and 1
            [z, ~] = zwgll(1);
            testCase.verifyEqual(z, [-1; 1], 'AbsTol', 1e-14);
        end

        function testP2(testCase)
            % p=2: three nodes at -1, 0, 1
            [z, ~] = zwgll(2);
            testCase.verifyEqual(z, [-1; 0; 1], 'AbsTol', 1e-14);
        end
    end
end
