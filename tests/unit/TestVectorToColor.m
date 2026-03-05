classdef TestVectorToColor < matlab.unittest.TestCase
% TestVectorToColor: Unit tests for vector_to_color

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        function testAxisDirections(testCase)
            % X-axis -> red, Y-axis -> green, Z-axis -> blue
            colors = vector_to_color([1; 0; 0], [0; 1; 0], [0; 0; 1]);
            testCase.verifyEqual(colors(1, :), [1 0 0], 'AbsTol', 1e-10);
            testCase.verifyEqual(colors(2, :), [0 1 0], 'AbsTol', 1e-10);
            testCase.verifyEqual(colors(3, :), [0 0 1], 'AbsTol', 1e-10);
        end

        function testNormalization(testCase)
            % (2,0,0) should normalize to red [1,0,0]
            colors = vector_to_color(2, 0, 0);
            testCase.verifyEqual(colors, [1 0 0], 'AbsTol', 1e-10);
        end

        function testZeroVector(testCase)
            % Zero vector should map to [0,0,0] (NaN handled)
            colors = vector_to_color(0, 0, 0);
            testCase.verifyEqual(colors, [0 0 0]);
        end

        function testNegativeComponents(testCase)
            % Negative directions produce same color as positive (abs used)
            c_pos = vector_to_color(1, 0, 0);
            c_neg = vector_to_color(-1, 0, 0);
            testCase.verifyEqual(c_pos, c_neg, 'AbsTol', 1e-10);
        end

        function testMultipleVectors(testCase)
            Vx = [1; 0; 0];
            Vy = [0; 1; 0];
            Vz = [0; 0; 1];
            colors = vector_to_color(Vx, Vy, Vz);
            testCase.verifySize(colors, [3, 3]);
        end

        function testOutputRange(testCase)
            % All color values must be in [0, 1]
            rng(42);
            Vx = randn(10, 1);
            Vy = randn(10, 1);
            Vz = randn(10, 1);
            colors = vector_to_color(Vx, Vy, Vz);
            testCase.verifyGreaterThanOrEqual(colors, 0);
            testCase.verifyLessThanOrEqual(colors, 1);
        end

        function testDiagonalVector(testCase)
            % (1,1,1) should produce equal RGB components
            colors = vector_to_color(1, 1, 1);
            testCase.verifyEqual(colors(1), colors(2), 'AbsTol', 1e-10);
            testCase.verifyEqual(colors(2), colors(3), 'AbsTol', 1e-10);
        end
    end
end
