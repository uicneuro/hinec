classdef TestNimReshapeD < matlab.unittest.TestCase
% TestNimReshapeD: Unit tests for nim_reshape_d
%
% nim_reshape_d maps a 6-element vector [d1,d2,d3,d4,d5,d6] to:
%   D = [ d1 d4 d6 ;
%         d4 d2 d5 ;
%         d6 d5 d3 ]

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        function testIdentity(testCase)
            % [1 1 1 0 0 0] should produce eye(3)
            D = nim_reshape_d([1 1 1 0 0 0]);
            testCase.verifyEqual(D, eye(3));
        end

        function testSymmetry(testCase)
            % Output matrix must always be symmetric
            d = [2.5, 1.3, 0.7, 0.4, -0.1, 0.3];
            D = nim_reshape_d(d);
            testCase.verifyEqual(D, D', 'Output must be symmetric');
        end

        function testKnownValues(testCase)
            % Verify exact placement of each element
            d = [1 2 3 4 5 6];
            D = nim_reshape_d(d);
            expected = [1 4 6; 4 2 5; 6 5 3];
            testCase.verifyEqual(D, expected);
        end

        function testOutputSize(testCase)
            D = nim_reshape_d([1 2 3 4 5 6]);
            testCase.verifySize(D, [3 3]);
        end

        function testZeroInput(testCase)
            D = nim_reshape_d([0 0 0 0 0 0]);
            testCase.verifyEqual(D, zeros(3));
        end
    end
end
