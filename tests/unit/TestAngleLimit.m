classdef TestAngleLimit < matlab.unittest.TestCase
% The angle criterion, as a unit. Every case here is a bug that was live in the
% tracker: the budget scaled by the realised chord, the turn measured against
% that chord, no ceiling at 90 degrees, and three trackers reading the same
% config key in three different units.

    properties
        Root
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fileparts(fileparts(here));
            addpath(genpath(tc.Root));
        end
    end

    methods (Test)

        function budgetScalesWithArc(tc)
            % The whole point of a RATE: halve the step, halve the budget, so the
            % same curvature bound is enforced at every step size.
            for h = [0.5 0.25 0.125 0.0625]
                tc.verifyEqual(nim_angle_limit(60, h), 60*h, 'AbsTol', 1e-12);
            end
        end

        function ceilingIsNinetyDegrees(tc)
            % Tangents are sign-aligned, so a measured turn is at most 90 degrees.
            % A budget above that cannot fire and must not be reported as if it
            % could - clamping is what makes the inertness visible.
            tc.verifyEqual(nim_angle_limit(225, 0.5), 90, 'AbsTol', 1e-12);
            tc.verifyEqual(nim_angle_limit(3600, 1.0), 90, 'AbsTol', 1e-12);
            % ... and 225 at the default step 0.2 is the classic 45 deg/step,
            % which is where that default came from and is NOT clamped.
            tc.verifyEqual(nim_angle_limit(225, 0.2), 45, 'AbsTol', 1e-12);
        end

        function zeroDisablesRatherThanForbids(tc)
            % A rate of 0 must mean "no angle criterion", not "no turning
            % allowed". Read naively, 0 gives a 0 degree budget and terminates
            % every streamline at its second point.
            [deg, c] = nim_angle_limit(0, 0.5);
            tc.verifyEqual(deg, Inf);
            tc.verifyEqual(c, -1);
            % dot(u,v) >= -1 for unit vectors, so the test `dot < c` never fires.
            tc.verifyFalse(-1 < c);
        end

        function trackersAgreeOnUnits(tc)
            % angle_max is degrees per VOXEL OF ARC in all three trackers. The
            % regression this guards: the config rework changed the units of the
            % key and only hinec was taught the new meaning, leaving FACT and MMF
            % reading a per-voxel rate as a per-step angle - cos(225 deg) is
            % negative, so their angle criterion was silently dead.
            for f = {'nim_tractography_hinec', 'nim_tractography_standard', ...
                     'nim_tractography_mmf_connframe'}
                src = fileread(fullfile(tc.Root, 'src', 'nim_tractography', [f{1} '.m']));
                tc.verifyTrue(contains(src, 'nim_angle_limit'), ...
                    sprintf('%s must budget its angle criterion through nim_angle_limit', f{1}));
                tc.verifyFalse(contains(src, 'cos(deg2rad(options.angle_thresh))'), ...
                    sprintf('%s still reads angle_thresh as a per-step angle', f{1}));
            end
        end

        function hinecBudgetsByNominalStepNotChord(tc)
            % The defect that biased the method comparison. Euler and RK2 advance
            % by h*(unit vector) so their chord is exactly h; RK4 averages four
            % stage vectors and its chord fell to 0.25*h on real data. Budgeting
            % by the chord therefore imposed a 4x tighter angle limit on RK4 alone
            % - inside an experiment whose purpose was to compare RK4 against
            % Euler and RK2.
            src = fileread(fullfile(tc.Root, 'src', 'nim_tractography', ...
                                    'nim_tractography_hinec.m'));
            tc.verifyTrue(contains(src, 'nim_angle_limit(angle_rate, last_step_arc)'), ...
                'hinec must budget the turn by the nominal step arc.');
            tc.verifyTrue(contains(src, 'dot(prev_dir, dir_vec)'), ...
                'hinec must measure the turn between consecutive TANGENTS.');
            tc.verifyFalse(contains(src, 'angle_rate * step_len'), ...
                'hinec still scales the angle budget by the realised chord length.');
        end

        function configsDeclareAnActiveBudget(tc)
            % A config whose angle_max * step lands at or above 90 has no angle
            % criterion at all. That is legal only when it is asked for (0).
            cfgs = dir(fullfile(tc.Root, 'config', '*.yml'));
            for i = 1:numel(cfgs)
                c = load_config_yaml(fullfile(cfgs(i).folder, cfgs(i).name));
                if ~isfield(c, 'tractography') || ~isfield(c.tractography, 'termination')
                    continue;
                end
                a = c.tractography.termination.angle_max;
                if a == 0, continue; end   % deliberately disabled
                h = c.tractography.integrator.step;
                tc.verifyLessThan(a * h, 90, sprintf( ...
                    ['%s: angle_max %g at step %g gives a %.1f deg budget, at or above ' ...
                     'the 90 deg ceiling - the criterion cannot fire. Use 0 if that is ' ...
                     'intended.'], cfgs(i).name, a, h, a*h));
            end
        end
    end
end
