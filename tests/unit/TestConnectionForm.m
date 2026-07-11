classdef TestConnectionForm < matlab.unittest.TestCase
% TestConnectionForm: Unit tests for nim_connection_form (the connection 1-form
% field omega). Uses an analytic frame field with a hand-computable connection.
%
% Analytic oracle -- frame rotating about z by angle theta = c * x1:
%   e1 = ( cos(c x), sin(c x), 0 )
%   e2 = (-sin(c x), cos(c x), 0 )
%   e3 = ( 0, 0, 1 )
% Then de1/dx1 = c e2, and the only nonzero form is
%   omega12(e_k) = c * e_k^1     (component of e_k along axis 1),
% while omega13 = omega23 = 0. (See §2.4; derived in the test docstring.)

    methods (TestClassSetup)
        function addPaths(testCase) %#ok<MANU>
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
        end
    end

    methods (Test)
        function testRotationAboutZOracle(testCase)
            c = 0.1;                      % small so central-difference error is tiny
            [frames, X] = makeRotZFrames([30 24 8], c);
            wijk = nim_connection_form(frames, [], struct());

            % Interior nodes only (finite differences degrade at the boundary).
            in1 = 3:28; in2 = 3:22; in3 = 3:6;

            % Extract e_k^1 (component along axis 1) for each k.
            e1x = frames(:,:,:,1,1);
            e2x = frames(:,:,:,1,2);
            e3x = frames(:,:,:,1,3);

            w12_e1 = wijk(:,:,:,1,2,1);
            w12_e2 = wijk(:,:,:,1,2,2);
            w12_e3 = wijk(:,:,:,1,2,3);

            % omega12(e_k) should equal c * e_k^1.
            testCase.verifyLessThan(relErr(w12_e1(in1,in2,in3), c*e1x(in1,in2,in3)), 1e-2);
            testCase.verifyLessThan(relErr(w12_e2(in1,in2,in3), c*e2x(in1,in2,in3)), 1e-2);
            testCase.verifyLessThan(maxAbs(w12_e3(in1,in2,in3) - c*e3x(in1,in2,in3)), 1e-2);

            % omega13 and omega23 vanish for this field.
            for k = 1:3
                testCase.verifyLessThan(maxAbs(wijk(in1,in2,in3,1,3,k)), 1e-2);
                testCase.verifyLessThan(maxAbs(wijk(in1,in2,in3,2,3,k)), 1e-2);
            end
        end

        function testAntisymmetryAndDiagonal(testCase)
            [frames, ~] = makeRotZFrames([20 20 6], 0.15);
            wijk = nim_connection_form(frames, [], struct());
            for k = 1:3
                % Diagonal is exactly zero.
                testCase.verifyEqual(maxAbs(wijk(:,:,:,1,1,k)), 0, 'AbsTol', 1e-12);
                testCase.verifyEqual(maxAbs(wijk(:,:,:,2,2,k)), 0, 'AbsTol', 1e-12);
                testCase.verifyEqual(maxAbs(wijk(:,:,:,3,3,k)), 0, 'AbsTol', 1e-12);
                % omega_ij = -omega_ji.
                testCase.verifyLessThan(maxAbs(wijk(:,:,:,1,2,k) + wijk(:,:,:,2,1,k)), 1e-12);
                testCase.verifyLessThan(maxAbs(wijk(:,:,:,1,3,k) + wijk(:,:,:,3,1,k)), 1e-12);
                testCase.verifyLessThan(maxAbs(wijk(:,:,:,2,3,k) + wijk(:,:,:,3,2,k)), 1e-12);
            end
        end

        function testOutputShape(testCase)
            [frames, ~] = makeRotZFrames([12 10 5], 0.1);
            wijk = nim_connection_form(frames, [], struct());
            testCase.verifySize(wijk, [12 10 5 3 3 3]);
        end
    end
end

% ---- helpers ---------------------------------------------------------------
function [frames, X] = makeRotZFrames(dims, c)
% Build the analytic rotation-about-z frame field on a dims grid.
[X, ~, ~] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
th = c * X;
frames = zeros([dims, 3, 3]);
% e1
frames(:,:,:,1,1) = cos(th);  frames(:,:,:,2,1) = sin(th);  frames(:,:,:,3,1) = 0;
% e2
frames(:,:,:,1,2) = -sin(th); frames(:,:,:,2,2) = cos(th);  frames(:,:,:,3,2) = 0;
% e3
frames(:,:,:,1,3) = 0;        frames(:,:,:,2,3) = 0;        frames(:,:,:,3,3) = 1;
end

function e = relErr(a, b)
% Relative Frobenius error, guarded against a zero reference.
den = norm(b(:));
if den < 1e-12, den = 1; end
e = norm(a(:) - b(:)) / den;
end

function m = maxAbs(a)
m = max(abs(a(:)));
end
