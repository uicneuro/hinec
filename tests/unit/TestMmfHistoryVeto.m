classdef TestMmfHistoryVeto < matlab.unittest.TestCase
% The history veto of the mmf hold scheme (advanceFrame_ml 'VetoK'/'VetoTol', tractography.mmf.veto_arc /
% veto_tol). A crossing in which the tract's bundle is the minor one is fitted as a LINEAR tensor along
% the other bundle: the confidence CL is high, so the zone never opens, and re-anchoring follows the wrong
% e1. The veto compares the data's e1 at a step's end with the heading veto_arc back; above veto_tol the
% tract returns there and coasts on its ring.

    methods (TestMethodSetup)
        function setup(tc) %#ok<MANU>
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
        end
    end

    methods
        function o = opts(~, nim, seed, scheme)
            cfg = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti.yml'));
            o = nim_config_to_options(cfg);
            o.seed_mask = false(size(nim.FA)); o.seed_mask(seed(1), seed(2), seed(3)) = true; o.seed_density = 8;
            o.wm_mask = []; o.gm_mask = []; o.csf_mask = [];
            o.min_length = 0; o.mmf_connection = 'nodal'; o.mmf_hold = scheme; o.mmf_reanchor_tol = 0.2;
            o.step_size = 0.2; o.max_arc = 100; o.max_steps = 500; o.angle_thresh = 0;
            o.mmf_hold_min = 0.15; o.mmf_exit_tol = 1.2;
        end
        function nim = hiddenCrossing(~)
            % make_synthetic_nim('crossing') with the crossing fitted as a linear tensor along band B (y):
            % CL = (l1 - l2)/l1 = 0.82 there, far above any hold_min, and e1 90 degrees from band A.
            nim = make_synthetic_nim('crossing');
            lam = [1.7e-3 3e-4 3e-4];
            for x = 40:46
                for y = 27:33
                    nim.eval(x, y, :, :) = repmat(reshape(lam, 1, 1, 1, 3), [1 1 size(nim.FA, 3) 1]);
                    nim.FA(x, y, :) = nim.FA(1, 30, 1);
                end
            end
        end
        function nim = gradualCrossing(~)
            % Band A of make_synthetic_nim('crossing') with e1 turning from x towards y over x = 38..40 at
            % 30 degrees per voxel (seed 19 of ISMRM 2015 turned at ~50 per voxel) and along y over x = 41..42,
            % linear tensors throughout (CL = 0.82): re-anchoring follows the turn a few degrees per step and no
            % node is in the zone. From x = 43 band A runs along x again (linear: the pancakes of the fixture are gone). The data are wrong over 5 voxels; a
            % coast much longer than that drifts with whatever curvature the ring holds (4h coasts on a circle). A turn as slow as a real bundle's
            % bend (the 95th percentile of the ISMRM bundles is 16 degrees per voxel) cannot be told from one
            % by the heading history, and is not what this fixture tests.
            nim = make_synthetic_nim('crossing');
            lam = [1.7e-3 3e-4 3e-4];
            for x = 38:46
                t = deg2rad(min(30 * (x - 37), 90)) * (x <= 42);
                for y = 27:33
                    nim.evec(x, y, :, :, 1) = repmat(reshape([cos(t) sin(t) 0], 1, 1, 1, 3), [1 1 size(nim.FA, 3) 1]);
                    nim.evec(x, y, :, :, 2) = repmat(reshape([-sin(t) cos(t) 0], 1, 1, 1, 3), [1 1 size(nim.FA, 3) 1]);
                    nim.eval(x, y, :, :) = repmat(reshape(lam, 1, 1, 1, 3), [1 1 size(nim.FA, 3) 1]);
                    nim.FA(x, y, :) = nim.FA(1, 30, 1);
                end
            end
        end
        function c = countOf(~, meta, reason)
            j = strcmp(meta.termination.reason, reason);
            if any(j), c = meta.termination.count(j); else, c = 0; end
        end
    end

    methods (Test)
        function theHoldAloneDoesNotSeeAHiddenCrossing(tc)
            % The failure the veto is for: no node of the crossing is in the zone, so neither scheme
            % coasts, and no tract reaches band A beyond the crossing.
            nim = tc.hiddenCrossing();
            for scheme = {'hold', 'march'}
                o = tc.opts(nim, [30 30 3], scheme{1});
                evalc('[T, meta] = nim_tractography_mmf_ml(nim, o);');
                tc.verifyEqual(meta.mmf.coasted_fraction, 0, scheme{1});
                tc.verifyFalse(any(cellfun(@(t) max(t(:, 1)) > 60 && max(abs(t(:, 2) - 30)) < 1, T)), scheme{1});
            end
        end
        function theVetoCarriesTheTractThroughAHiddenCrossing(tc)
            % With the veto (2 voxels back, 30 degrees) every tract returns to the node before the turn,
            % coasts straight along band A through the crossing, re-enters band A and leaves the volume.
            nim = tc.hiddenCrossing();
            for scheme = {'hold', 'march'}
                o = tc.opts(nim, [30 30 3], scheme{1}); o.mmf_veto_arc = 2; o.mmf_veto_tol = deg2rad(30);
                evalc('[T, meta] = nim_tractography_mmf_ml(nim, o);');
                tc.assertNumElements(T, 8, scheme{1});
                tc.verifyEqual(meta.termination.reason, {'outside'}, scheme{1});
                tc.verifyTrue(all(cellfun(@(t) max(t(:, 1)), T) > 79), scheme{1});
                tc.verifyLessThan(max(cellfun(@(t) max(abs(t(:, 2) - 30)), T)), 0.5, scheme{1});
                tc.verifyGreaterThan(meta.mmf.veto_returns, 0, scheme{1});
                tc.verifyGreaterThan(meta.mmf.coasted_fraction, 0, scheme{1});
                evalc('nim_check_tracker_output(T, meta);');
            end
        end
        function aVetoedCoastPassesAGradualTurn(tc)
            % The coast that follows a return must reach the node where the veto fired before the exit test
            % may end it: there the data still agree with the heading they are about to leave, and an end
            % there re-anchors onto the turn one step at a time (seed 19 of ISMRM 2015, 2026-11-06).
            nim = tc.gradualCrossing();
            o = tc.opts(nim, [30 30 3], 'hold');
            evalc('[T0, ~] = nim_tractography_mmf_ml(nim, o);');
            tc.verifyFalse(any(cellfun(@(t) max(t(:, 1)) > 60 && max(abs(t(:, 2) - 30)) < 1, T0)));
            for scheme = {'hold', 'march'}
                o = tc.opts(nim, [30 30 3], scheme{1}); o.mmf_veto_arc = 2; o.mmf_veto_tol = deg2rad(30);
                o.trace = true; o.trace_max = 0;
                evalc('[T, meta] = nim_tractography_mmf_ml(nim, o);');
                % Every tract passes the crossing on band A and none turns into band B (x 40..46, along y).
                % What follows the crossing is not the veto's: the coast ends where the field turns back and
                % re-anchors to an e1 a few degrees off, and adaptive re-anchoring (reanchor_tol 0.2 rad) never
                % corrects a heading error below 11.5 degrees on a straight bundle, so tracts slide 1-3.5 voxels
                % sideways downstream.
                tc.verifyTrue(all(cellfun(@(t) max(t(:, 1)), T) > 50), scheme{1});
                inB = @(t) any(t(:, 1) >= 40 & t(:, 1) <= 46 & abs(t(:, 2) - 30) > 4);
                tc.verifyFalse(any(cellfun(inB, T)), scheme{1});
                f = meta.trace(1).forward;
                tc.verifyFalse(any(f.coasted(2:end) & f.reanchored(2:end)), scheme{1});   % a coasting node is never re-anchored
            end
        end
        function theTraceStaysConsistentAfterAReturn(tc)
            % A return truncates every per-node and per-step record to the node returned to.
            nim = tc.hiddenCrossing();
            o = tc.opts(nim, [30 30 3], 'hold'); o.mmf_veto_arc = 2; o.mmf_veto_tol = deg2rad(30);
            o.trace = true; o.trace_max = 0;
            evalc('[~, meta] = nim_tractography_mmf_ml(nim, o);');
            f = meta.trace(1).forward;
            n = size(f.pos, 1);
            tc.verifyEqual(numel(f.coasted), n); tc.verifyEqual(numel(f.confidence), n);
            tc.verifyEqual(numel(f.kappa), n); tc.verifyEqual(numel(f.s), n);
            d = sqrt(sum(diff(f.pos) .^ 2, 2));
            tc.verifyLessThan(max(abs(d - 0.2)), 1e-6);              % consecutive nodes one step apart: nothing left of the turn
            x = f.pos(f.coasted, 1);
            tc.verifyLessThan(min(x), 40); tc.verifyGreaterThan(max(x), 46);   % the coast spans the crossing
        end
        function aCurvedBundleIsNotVetoed(tc)
            % The circles of radius ~15 about (40, 40) of make_synthetic_nim('arc') turn by 2/15 rad = 7.6 degrees per
            % 2 voxels of arc: below veto_tol, so the veto changes nothing.
            nim = make_synthetic_nim('arc');
            o = tc.opts(nim, [55 40 4], 'hold'); o.max_arc = 40;
            evalc('[T0, m0] = nim_tractography_mmf_ml(nim, o);');
            o.mmf_veto_arc = 2; o.mmf_veto_tol = deg2rad(30);
            evalc('[T1, m1] = nim_tractography_mmf_ml(nim, o);');
            tc.verifyEqual(m1.mmf.veto_returns, 0);
            tc.verifyEqual(T1, T0);
            tc.verifyEqual(m0.mmf.veto_arc, 0);
        end
        function offIsTheDefaultAndChangesNothing(tc)
            c = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti_hold.yml'));
            tc.verifyEqual(c.tractography.mmf.veto_arc, 0);
            nim = make_synthetic_nim('crossing');
            o = tc.opts(nim, [30 30 3], 'hold');
            evalc('[T0, m0] = nim_tractography_mmf_ml(nim, o);');
            o.mmf_veto_tol = 0.1;                                    % the tolerance alone does nothing
            evalc('[T1, m1] = nim_tractography_mmf_ml(nim, o);');
            tc.verifyEqual(T1, T0); tc.verifyEqual(m1.mmf.veto_returns, 0);
            tc.verifyEqual(m1.termination, m0.termination);
        end
        function theVetoNeedsTheHold(tc)
            nim = make_synthetic_nim('crossing');
            o = tc.opts(nim, [30 30 3], 'none'); o.mmf_veto_arc = 2;
            tc.verifyError(@() nim_tractography_mmf_ml(nim, o), 'mmf:hold');
            dom = struct('xlim', [0 1], 'ylim', [0 1]);
            tc.verifyError(@() advanceFrame_ml(@(X) zeros(size(X, 1), 2, 2, 2), [], dom, [0.5 0.5], [1 0], 0.1, 1, ...
                'Frame0', eye(2), 'VetoK', 3), 'advanceFrame:badOptions');
        end
    end
end
