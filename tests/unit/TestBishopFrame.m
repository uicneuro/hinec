classdef TestBishopFrame < matlab.unittest.TestCase
% TestBishopFrame: Unit tests for the MMF Phase-1 frame helpers
% (mmf_reference_axis_frame, mmf_bishop_update, mmf_gram_schmidt).

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        % --- reference-axis projection (§2.2) --------------------------------
        function testRefAxisOrthonormal(testCase)
            rng(1);
            for k = 1:50
                e1 = randn(1, 3); e1 = e1 / norm(e1);
                [e2, e3] = mmf_reference_axis_frame(e1);
                testCase.verifyEqual(norm(e2), 1, 'AbsTol', 1e-12);
                testCase.verifyEqual(norm(e3), 1, 'AbsTol', 1e-12);
                testCase.verifyEqual(dot(e1, e2), 0, 'AbsTol', 1e-12);
                testCase.verifyEqual(dot(e1, e3), 0, 'AbsTol', 1e-12);
                testCase.verifyEqual(dot(e2, e3), 0, 'AbsTol', 1e-12);
                % right-handed: e3 == e1 x e2
                testCase.verifyEqual(e3, cross(e1, e2), 'AbsTol', 1e-12);
            end
        end

        function testRefAxisPicksLeastAlignedAxis(testCase)
            % e1 mostly along x -> chosen reference axis is NOT x, so e2 has a
            % large component along the least-aligned axis.
            e1 = [0.98, 0.14, 0.14]; e1 = e1 / norm(e1);
            [e2, ~] = mmf_reference_axis_frame(e1);
            testCase.verifyLessThan(abs(e2(1)), max(abs(e2(2)), abs(e2(3))));
        end

        % --- Bishop transport (§2.3) -----------------------------------------
        function testBishopStraightSegmentNoRotation(testCase)
            % Parallel tangents -> e2 only re-orthogonalized (essentially fixed).
            t0 = [1 0 0]; t1 = [1 0 0];
            e2 = [0 1 0];
            [e2n, e3n] = mmf_bishop_update(t0, t1, e2, 1e-6);
            testCase.verifyEqual(e2n, [0 1 0], 'AbsTol', 1e-12);
            testCase.verifyEqual(e3n, cross(t1, e2n), 'AbsTol', 1e-12);
        end

        function testBishopOrthonormalAfterRotation(testCase)
            rng(3);
            for k = 1:50
                t0 = randn(1,3); t0 = t0/norm(t0);
                t1 = randn(1,3); t1 = t1/norm(t1);
                [e2, ~] = mmf_reference_axis_frame(t0);
                [e2n, e3n] = mmf_bishop_update(t0, t1, e2, 1e-6);
                testCase.verifyEqual(norm(e2n), 1, 'AbsTol', 1e-10);
                testCase.verifyEqual(norm(e3n), 1, 'AbsTol', 1e-10);
                % Frame must be orthogonal to the NEW tangent.
                testCase.verifyEqual(dot(e2n, t1), 0, 'AbsTol', 1e-10);
                testCase.verifyEqual(dot(e3n, t1), 0, 'AbsTol', 1e-10);
                testCase.verifyEqual(dot(e2n, e3n), 0, 'AbsTol', 1e-10);
            end
        end

        function testBishopMinimalRotation(testCase)
            % 90-degree turn in the xy-plane: e2 initially = z stays ~ z
            % (rotation axis is z, e2 along z is invariant).
            t0 = [1 0 0]; t1 = [0 1 0];
            e2 = [0 0 1];
            [e2n, ~] = mmf_bishop_update(t0, t1, e2, 1e-6);
            testCase.verifyEqual(e2n, [0 0 1], 'AbsTol', 1e-10);
        end

        % --- Gram-Schmidt ----------------------------------------------------
        function testGramSchmidtReorthonormalizes(testCase)
            e1 = [1 0.02 -0.01];
            e2 = [0.03 1 0.02];
            e3 = [0 0 1];
            [q1, q2, q3] = mmf_gram_schmidt(e1, e2, e3);
            M = [q1; q2; q3];
            testCase.verifyEqual(M * M.', eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(q3, cross(q1, q2), 'AbsTol', 1e-12);
        end
    end
end
