classdef TestAngleCriterionBehaviour < matlab.unittest.TestCase
% BEHAVIOURAL tests of the angle termination criterion, running the actual
% tracker. The companion TestAngleLimit covers nim_angle_limit as a function and
% then checks the trackers by GREPPING THEIR SOURCE, which proves only that some
% text is present - it would pass against a tracker that never calls the helper
% on a live streamline. These tests integrate real streamlines instead.
%
% The field is a rigid rotation about the volume centre, so v1 is everywhere
% tangent to a circle and a streamline seeded at radius R follows that circle
% exactly. Its turning rate is therefore KNOWN in closed form:
%
%       57.2958 / R   degrees per voxel of arc
%
% which makes the criterion's firing threshold analytically predictable rather
% than a number read off a previous run.

    properties (Constant)
        R  = 15;                       % seed radius, voxels
        C  = 30.5;                     % centre of the volume, voxels
    end
    properties
        Nim
        Rate                            % true turning rate, deg per voxel of arc
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            tc.Rate = 180/pi / tc.R;
            tc.Nim  = tc.buildRotationField();
        end
    end

    methods
        function nim = buildRotationField(tc)
            % v1(x,y,z) = tangent to the circle about (C,C). FA is 1 inside a wide
            % annulus so nothing terminates on FA before the angle criterion or
            % the arc cap gets its chance.
            n = 60; nz = 9;
            [X, Y, ~] = ndgrid(1:n, 1:n, 1:nz);
            dx = X - tc.C; dy = Y - tc.C;
            r  = sqrt(dx.^2 + dy.^2); r(r < 1e-9) = 1e-9;
            vx = -dy ./ r; vy = dx ./ r; vz = zeros(size(r));

            evec = zeros(n, n, nz, 3, 3);
            evec(:,:,:,1,1) = vx; evec(:,:,:,2,1) = vy; evec(:,:,:,3,1) = vz;
            % v2, v3 arbitrary but orthonormal - the tracker only reads column 1.
            evec(:,:,:,1,2) = vy; evec(:,:,:,2,2) = -vx;
            evec(:,:,:,3,3) = 1;

            FA = zeros(n, n, nz);
            FA(r > 4 & r < 27) = 1;

            nim = struct('FA', FA, 'evec', evec, ...
                         'eval', repmat(reshape([3 1 1], 1,1,1,3), n, n, nz), ...
                         'mask', FA > 0);
        end

        function [tracks, opts] = track(tc, angle_max, step, order, min_arc)
            if nargin < 5 || isempty(min_arc), min_arc = 0; end
            seed = false(size(tc.Nim.FA));
            seed(round(tc.C + tc.R), round(tc.C), 5) = true;
            opts = struct( ...
                'seed_mask', seed, 'seed_density', 1, 'seed_strategy', "uniform", ...
                'step_size', step, 'angle_thresh', angle_max, ...
                'integration_order', order, 'interp_method', "trilinear", ...
                'termination_fa', 0.5, 'seed_fa_threshold', 0.5, ...
                'min_length', min_arc, 'max_steps', ceil(60 / step), ...
                'enable_diagnostics', false, 'act_enabled', false, ...
                'wm_mask', [], 'gm_mask', [], 'csf_mask', []);
            evalc('tracks = nim_tractography_hinec(tc.Nim, opts);');
        end

        function L = arc(~, tracks)
            L = 0;
            for i = 1:numel(tracks)
                p = tracks{i};
                if size(p,1) > 1, L = max(L, sum(sqrt(sum(diff(p,1,1).^2, 2)))); end
            end
        end
    end

    methods (Test)

        function firesOnlyWhenBudgetIsBelowTheTrueTurningRate(tc)
            % The field turns at 3.82 deg per voxel of arc. A budget below that
            % must stop the streamline almost at once; a budget above it must not
            % interfere at all. Anything else means the criterion is not measuring
            % turning per unit arc.
            tight = tc.arc(tc.track(tc.Rate * 0.5, 0.25, 4));
            loose = tc.arc(tc.track(tc.Rate * 2.0, 0.25, 4));
            tc.verifyLessThan(tight, 5, ...
                'A budget below the true turning rate should terminate immediately.');
            tc.verifyGreaterThan(loose, 55, ...
                'A budget above the true turning rate should not terminate at all.');
        end

        function budgetIsStepInvariant(tc)
            % The point of expressing the limit as a RATE. The same physical
            % constraint must produce the same outcome across an 8x change in
            % step, which is the regression that matters for a step-size ladder:
            % a per-STEP angle would let the loose case survive at fine steps and
            % kill it at coarse ones.
            for h = [0.4 0.2 0.1 0.05]
                tc.verifyGreaterThan(tc.arc(tc.track(tc.Rate * 2.0, h, 4)), 55, ...
                    sprintf('loose budget wrongly terminated at h = %g', h));
                tc.verifyLessThan(tc.arc(tc.track(tc.Rate * 0.5, h, 4)), 5, ...
                    sprintf('tight budget failed to terminate at h = %g', h));
            end
        end

        function budgetIsIndependentOfIntegrationMethod(tc)
            % THE defect this whole exercise was about. The budget used to be
            % scaled by the realised step CHORD; Euler and RK2 advance by
            % h*(unit vector) so their chord is exactly h, while RK4 averages four
            % stage vectors and its chord is shorter. That handed RK4 a tighter
            % angle limit than Euler at the same nominal step - inside an
            % experiment built to compare the two. All three must now agree.
            for order = [1 2 4]
                tc.verifyGreaterThan(tc.arc(tc.track(tc.Rate * 2.0, 0.25, order)), 55, ...
                    sprintf('loose budget wrongly terminated for integration_order %d', order));
                tc.verifyLessThan(tc.arc(tc.track(tc.Rate * 0.5, 0.25, order)), 5, ...
                    sprintf('tight budget failed to terminate for integration_order %d', order));
            end
        end

        function zeroDisablesTheCriterionInTheTracker(tc)
            % 0 must mean "no angle criterion", not "no turning allowed". Read
            % naively it is a 0 degree budget, which would cut every streamline at
            % its second point.
            tc.verifyGreaterThan(tc.arc(tc.track(0, 0.25, 4)), 55, ...
                'angle_max = 0 must disable the criterion, not forbid all turning.');
        end

        function factTrackerHonoursTheSameUnits(tc)
            % The angle fix reaches nim_tractography_standard too, and until now
            % that was asserted only by grepping its source. FACT is piecewise
            % constant per voxel, so one transition is budgeted as one voxel of
            % arc: a rate below the field's true turning rate must stop it, a rate
            % above must not. Same field, same closed-form turning rate as the
            % streamline tests.
            seed = false(size(tc.Nim.FA));
            seed(round(tc.C + tc.R), round(tc.C), 5) = true;
            base = struct('seed_mask', seed, 'seed_density', 1, 'seed_strategy', "uniform", ...
                'step_size', 0.5, 'termination_fa', 0.5, 'seed_fa_threshold', 0.5, ...
                'min_length', 0, 'max_steps', 400, 'enable_diagnostics', false, ...
                'act_enabled', false, 'wm_mask', [], 'gm_mask', [], 'csf_mask', []);

            tight = base; tight.angle_thresh = tc.Rate * 0.5;
            loose = base; loose.angle_thresh = tc.Rate * 4;
            evalc('tt = nim_tractography_standard(tc.Nim, tight);');
            evalc('tl = nim_tractography_standard(tc.Nim, loose);');
            tc.verifyLessThan(tc.arc(tt), tc.arc(tl), ...
                ['FACT ignored the angle rate: a budget below the true turning rate ' ...
                 'produced a track no shorter than one well above it.']);
            tc.verifyGreaterThan(tc.arc(tl), 20, 'setup: loose budget should track freely');
        end

        function minArcIsActuallyEnforced(tc)
            % min_arc was defaulted at the top of hinec and standard and then
            % never read - the standard tracker even carried a comment saying
            % "Save ALL generated tracks - no filters at all", which was true.
            % A schema-validated key that every config sets was honoured by one
            % of the three trackers.
            long_run = tc.arc(tc.track(0, 0.25, 4, 0));
            tc.assertGreaterThan(long_run, 20, 'setup: expected a long streamline');

            % A floor below the track's length keeps it...
            kept = tc.track(0, 0.25, 4, long_run * 0.5);
            tc.verifyNotEmpty(kept, 'min_arc below the track length wrongly discarded it.');

            % ... and a floor above it must discard the track entirely.
            dropped = tc.track(0, 0.25, 4, long_run * 2);
            tc.verifyEmpty(dropped, ...
                'min_arc above the track length did not discard it - the key is inert.');
        end

        function budgetAboveNinetyDegreesIsInert(tc)
            % Tangents are sign-aligned, so a measured turn never exceeds 90 deg.
            % At h = 0.5 a rate of 225 gives a 112 deg budget: it cannot fire, and
            % must behave exactly like the disabled criterion. This is the trap the
            % coarsest rung of the step ladder was sitting in.
            inert    = tc.arc(tc.track(225, 0.5, 4));
            disabled = tc.arc(tc.track(0,   0.5, 4));
            tc.verifyEqual(inert, disabled, 'AbsTol', 1e-9, ...
                'A budget above the 90 deg ceiling must be indistinguishable from disabled.');
        end
    end
end
