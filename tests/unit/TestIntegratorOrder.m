classdef TestIntegratorOrder < matlab.unittest.TestCase
% Order verification for the integrators themselves, on a synthetic field where
% the answer is known independently of any dMRI data.
%
% This exists because the ladder on real ISMRM data reports RK4 at order 2.02 -
% indistinguishable from RK2 in both order AND error constant - and that could
% mean either "the RK4 implementation is wrong" or "the measured field is too
% rough for high order to pay". Only a smooth field separates the two.
%
% The field is a rigid rotation: analytic, unit-norm everywhere, tensor
% eigenvalues a uniform 3/1/1 so the principal eigenvector is well conditioned.
% Every step size integrates the SAME arc, so comparing endpoints is a clean
% global-error measure with no termination and no windowing.

    properties
        Nim
        Seed
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            R = 15; C = 30.5; n = 60; nz = 9;
            [X,Y,~] = ndgrid(1:n,1:n,1:nz);
            dx = X-C; dy = Y-C; r = sqrt(dx.^2+dy.^2); r(r<1e-9) = 1e-9;
            evec = zeros(n,n,nz,3,3);
            evec(:,:,:,1,1) = -dy./r; evec(:,:,:,2,1) = dx./r;
            evec(:,:,:,1,2) =  dx./r; evec(:,:,:,2,2) = dy./r; evec(:,:,:,3,3) = 1;
            FA = zeros(n,n,nz); FA(r>4 & r<27) = 1;
            tc.Nim = struct('FA', FA, 'evec', evec, ...
                'eval', repmat(reshape([3 1 1],1,1,1,3),n,n,nz), 'mask', FA>0);
            tc.Seed = false(size(FA)); tc.Seed(round(C+R), round(C), 5) = true;
        end
    end

    methods
        function e = endpointAfterArc(tc, order, h, A)
            % max_steps = N+1 gives exactly N steps: track_length starts at 1 for
            % the seed and the loop stops at track_length >= max_steps. Using
            % round(A/h) instead gives every h a different arc, A-h, and the
            % endpoints then differ by one step length - a clean O(h) signal that
            % reads as first-order convergence for every method.
            o = struct('seed_mask', tc.Seed, 'seed_density', 1, 'seed_strategy', "uniform", ...
                'step_size', h, 'angle_thresh', 0, 'integration_order', order, ...
                'interp_method', "cubic", 'termination_fa', 0.5, 'seed_fa_threshold', 0.5, ...
                'min_length', 0, 'max_steps', round(A/h)+1, 'enable_diagnostics', false, ...
                'act_enabled', false, 'wm_mask', [], 'gm_mask', [], 'csf_mask', []);
            evalc('T = nim_tractography_hinec(tc.Nim, o);');
            tr = []; for i = 1:numel(T), if size(T{i},1) > size(tr,1), tr = T{i}; end, end
            e = tr(end,:);
        end
        function [err, p] = ladder(tc, order, H, A)
            ref = tc.endpointAfterArc(order, 1e-4, A);
            err = arrayfun(@(h) norm(tc.endpointAfterArc(order,h,A) - ref), H);
            c = polyfit(log(H), log(err), 1); p = c(1);
        end
    end

    methods (Test)

        function eulerIsFirstOrder(tc)
            [~, p] = tc.ladder(1, [0.4 0.2 0.1 0.05], 8);
            tc.verifyGreaterThan(p, 0.90); tc.verifyLessThan(p, 1.15);
        end

        function rk2IsSecondOrder(tc)
            [~, p] = tc.ladder(2, [0.4 0.2 0.1 0.05], 8);
            tc.verifyGreaterThan(p, 1.85); tc.verifyLessThan(p, 2.20);
        end

        function rk4BeatsRk2ByOrdersOfMagnitudeOnASmoothField(tc)
            % The claim under test: RK4's collapse to RK2's accuracy on real data
            % is a property of the DATA, not a broken tableau. On a smooth field
            % the two must separate enormously - measured ~750x at h = 0.4.
            H = [0.4 0.2 0.1];
            e2 = tc.ladder(2, H, 8);
            e4 = tc.ladder(4, H, 8);
            tc.verifyGreaterThan(e2(1)/e4(1), 100, ...
                ['RK4 is not materially more accurate than RK2 even on a smooth ' ...
                 'analytic field - that points at the RK4 implementation, not the data.']);
        end

        function rk4OutranksRk2OnASmoothField(tc)
            % Deliberately a COMPARISON, not an absolute threshold. RK4's observed
            % order here measures ~2.4 with visibly noisy local slopes (3.6, 1.0,
            % 3.2, 2.9), because even this field reaches the integrator through a
            % cubic-convolution interpolant, which is only C1 - short of the C4
            % that order 4 requires. Pinning a number would be pinning noise. What
            % is robust, and what this guards, is that RK4 converges strictly
            % faster than RK2 when the field is smooth, which is exactly what it
            % fails to do on real data (2.02 vs 2.01).
            [~, p2] = tc.ladder(2, [0.4 0.2 0.1 0.05], 8);
            [~, p4] = tc.ladder(4, [0.4 0.2 0.1 0.05], 8);
            tc.verifyGreaterThan(p4, p2 + 0.25, sprintf( ...
                ['RK4 (p=%.2f) did not outconverge RK2 (p=%.2f) even on a smooth ' ...
                 'field - that would point at the tableau, not the data.'], p4, p2));
        end
    end
end
