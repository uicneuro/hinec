classdef TestTrackResample < matlab.unittest.TestCase
    % Output decimation: resampling saved streamlines to fixed arc-length
    % spacing, so file size stops scaling as 1/integration_step.

    properties
        Root
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fullfile(here, '..', '..');
            addpath(fullfile(tc.Root, 'src', 'nim_utils'));
            addpath(fullfile(tc.Root, 'src', 'nim_tractography'));
        end
    end

    methods (Test)

        function offByDefaultReturnsInputUnchanged(tc)
            t = [zeros(101,1), zeros(101,1), (0:0.1:10)'];
            tc.verifyEqual(nim_resample_track_arc(t, 0), t);
            tc.verifyEqual(nim_resample_track_arc(t, []), t);
        end

        function spacingIsHonouredOnAStraightLine(tc)
            t = [zeros(101,1), zeros(101,1), (0:0.1:10)'];   % arc length 10
            for a = [0.5 1 2.5]
                r = nim_resample_track_arc(t, a);
                tc.verifyEqual(size(r,1), round(10/a) + 1, ...
                    sprintf('arc_step %.2f gave %d points', a, size(r,1)));
                d = sqrt(sum(diff(r,1,1).^2, 2));
                tc.verifyEqual(d, repmat(a, numel(d), 1), 'AbsTol', 1e-9);
            end
        end

        function endpointsArePreservedExactly(tc)
            % Endpoints carry the connectivity that scoring is built on; they
            % must never be interpolated away.
            th = linspace(0, pi, 400)';
            c = [10*cos(th), 10*sin(th), zeros(size(th))];
            r = nim_resample_track_arc(c, 0.5);
            tc.verifyEqual(r(1,:),   c(1,:));
            tc.verifyEqual(r(end,:), c(end,:));
        end

        function arcLengthIsEssentiallyPreserved(tc)
            th = linspace(0, pi, 400)';
            c = [10*cos(th), 10*sin(th), zeros(size(th))];
            Lc = sum(sqrt(sum(diff(c,1,1).^2, 2)));
            r  = nim_resample_track_arc(c, 0.5);
            Lr = sum(sqrt(sum(diff(r,1,1).^2, 2)));
            % Chord-vs-arc loss only; must stay far below a voxel.
            tc.verifyLessThan(abs(Lc - Lr) / Lc, 1e-3);
        end

        function trackShorterThanOneStepKeepsBothEndpoints(tc)
            t = [0 0 0; 0 0 0.1];
            r = nim_resample_track_arc(t, 5);
            tc.verifyEqual(size(r,1), 2);
            tc.verifyEqual(r, t);
        end

        function degenerateAndDuplicatePointsAreHandled(tc)
            tc.verifyEqual(size(nim_resample_track_arc([1 2 3], 0.5), 1), 1);   % single point
            tc.verifyEqual(size(nim_resample_track_arc(zeros(0,3), 0.5), 1), 0); % empty
            dup = [0 0 0; 0 0 0; 0 0 1; 0 0 1; 0 0 2];   % zero-length segments
            r = nim_resample_track_arc(dup, 0.5);
            tc.verifyEqual(size(r,1), 5);
            tc.verifyEqual(r(1,:),   dup(1,:));
            tc.verifyEqual(r(end,:), dup(end,:));
        end

        function cellArraysAreResampledElementwise(tc)
            t = [zeros(101,1), zeros(101,1), (0:0.1:10)'];
            out = nim_resample_track_arc({t, t}, 1);
            tc.verifyClass(out, 'cell');
            tc.verifyEqual(numel(out), 2);
            tc.verifyEqual(size(out{1},1), 11);
        end

        function decouplesStorageFromIntegrationStep(tc)
            % The point of the feature: the same physical curve sampled at two
            % very different integration steps must yield the SAME stored size.
            L = 10;
            fine   = [zeros(1001,1), zeros(1001,1), linspace(0,L,1001)'];  % step 0.01
            coarse = [zeros(51,1),   zeros(51,1),   linspace(0,L,51)'];    % step 0.2
            rf = nim_resample_track_arc(fine,   0.5);
            rc = nim_resample_track_arc(coarse, 0.5);
            tc.verifyEqual(size(rf,1), size(rc,1));
            tc.verifyEqual(rf, rc, 'AbsTol', 1e-9);
        end

        function configKeyIsWiredThrough(tc)
            S = nim_config_schema();
            tc.verifyTrue(any(strcmp({S.path}, 'tractography.output.arc_step')));
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            cfg = load_config_yaml(fullfile(tc.Root, 'config', 'hinec_dti.yml'));
            tc.verifyEqual(nim_config_to_options(cfg).output_arc_step, 0);  % off by default
            c2 = nim_config_apply_overrides(cfg, {'output.arc_step=0.5'});
            tc.verifyEqual(nim_config_to_options(c2).output_arc_step, 0.5);
        end
    end
end
