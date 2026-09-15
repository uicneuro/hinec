classdef TestACT < matlab.unittest.TestCase
% Anatomically constrained tracking: the two policy rules, and the trace that
% makes them auditable.
%
% ACT was previously unusable, and not because of the masks. Measured against the
% ISMRM ground truth on Cingulum_right:
%
%   * Discarding a streamline on ANY contact with a CSF voxel deletes 27% of the
%     GROUND TRUTH. Only 3.5% of ground-truth bundle voxels are labelled CSF, but
%     streamlines are long enough that a quarter of them clip one.
%   * Grey matter appears in the INTERIOR of 68% of ground-truth streamlines, so
%     terminating on GM contact truncates two thirds of them mid-bundle - while
%     64% of ground-truth ENDPOINTS are in grey matter, so it is exactly where a
%     streamline should be allowed to end.
%
% So CSF is tolerated up to a budget and then TRUNCATES (what came before the
% excursion is still valid), and GM marks a valid place to stop without forcing
% one.

    methods (Static)
        function [nim, opts, masks] = corridor()
        % A straight WM corridor along +y with a 2-voxel CSF speck and a 7-voxel
        % CSF slab further along.
            d = [20 60 20];
            ev = zeros([d 3 3]); ev(:,:,:,2,1) = 1;      % fibres along +y
            nim = struct('FA', 0.6*ones(d), 'evec', ev, ...
                         'eval', repmat(reshape([3 1 1],1,1,1,3), [d 1]), 'mask', true(d));
            wm = true(d); gm = false(d); cs = false(d);
            cs(:,25:26,:) = true; wm(:,25:26,:) = false;   % 2-voxel speck
            cs(:,40:46,:) = true; wm(:,40:46,:) = false;   % 7-voxel slab
            sm = false(d); sm(10,5,10) = true;
            masks = struct('wm', wm, 'gm', gm, 'csf', cs);
            opts = struct('step_size',0.5,'max_steps',200,'termination_fa',0.1, ...
                'integration_order',4,'interp_method','linear','field','dti', ...
                'enable_diagnostics',false,'adaptive_step',false, ...
                'seed_mask',sm,'seed_density',1,'seed_strategy','uniform', ...
                'min_arc',0,'max_arc',200,'trace',true,'trace_max',1, ...
                'wm_mask',wm,'gm_mask',gm,'csf_mask',cs,'act_csf_run',3);
        end
        function y = furthest(T)
            y = 0; for i = 1:numel(T), y = max(y, max(T{i}(:,2))); end
        end
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
        end
    end

    methods (Test)

        function briefCsfContactIsToleratedAndSustainedCsfTruncates(tc)
            [nim, o] = TestACT.corridor();
            o.act_csf_run = 3;                     % voxels of arc
            [T, ~] = evalc_track(nim, o);
            tc.assertNotEmpty(T, 'ACT removed every streamline.');
            y = TestACT.furthest(T);
            tc.verifyGreaterThan(y, 30, ...
                ['the 2-voxel CSF speck at y=25 must be tolerated under a 3-voxel ' ...
                 'budget; the streamline stopped at it instead.']);
            tc.verifyLessThan(y, 41, ...
                'the 7-voxel CSF slab at y=40 must stop the streamline.');
        end

        function theBudgetIsArcNotSteps(tc)
            % Halving the step size must not change where tracking stops. A
            % budget counted in STEPS would halve in physical terms.
            [nim, o] = TestACT.corridor();
            o.act_csf_run = 3;
            [Ta, ~] = evalc_track(nim, o);
            o.step_size = 0.25; o.max_steps = 800;
            [Tb, ~] = evalc_track(nim, o);
            ya = TestACT.furthest(Ta); yb = TestACT.furthest(Tb);
            tc.verifyLessThan(abs(ya - yb), 1.0, ...
                sprintf(['the CSF budget must be arc-based: halving the step moved the ' ...
                         'stopping point from y=%.1f to y=%.1f.'], ya, yb));
        end

        function truncationKeepsThePrefixInsteadOfDiscarding(tc)
            % The old policy set track_length to 0 on CSF contact, throwing away
            % every valid voxel the streamline had already traversed.
            [nim, o] = TestACT.corridor();
            o.act_csf_run = 1;                     % stop at the very first speck
            [T, ~] = evalc_track(nim, o);
            tc.assertNotEmpty(T, ...
                'a streamline that reaches CSF must keep the valid part before it.');
            tc.verifyGreaterThan(TestACT.furthest(T), 15, ...
                'the surviving prefix should run from the seed up to the CSF.');
        end

        function greyMatterDoesNotStopTracking(tc)
            % GM in the middle of a bundle must not truncate the streamline.
            [nim, o] = TestACT.corridor();
            o.gm_mask(:, 30:31, :) = true;         % a GM band mid-corridor
            o.wm_mask(:, 30:31, :) = false;
            [T, ~] = evalc_track(nim, o);
            tc.assertNotEmpty(T);
            tc.verifyGreaterThan(TestACT.furthest(T), 33, ...
                ['grey matter mid-bundle must not terminate tracking - it marks a ' ...
                 'valid place to END, and 68% of ground-truth streamlines cross it.']);
        end

        function everyActDecisionIsTraced(tc)
            [nim, o] = TestACT.corridor();
            [~, M] = evalc_track(nim, o);
            tc.assertTrue(isfield(M, 'trace') && ~isempty(M.trace), 'no trace recorded');
            tr = M.trace(1).forward;
            tc.assertTrue(isfield(tr, 'tissue') && isfield(tr, 'act'), ...
                'the trace must record the tissue class and the ACT decision per step.');
            acts = tr.act(~cellfun(@isempty, tr.act));
            tc.verifyTrue(any(strcmp(acts, 'csf_tolerated')), ...
                'tolerating a brief CSF contact must be visible in the trace.');
            tc.verifyTrue(any(strcmp(acts, 'truncate')), ...
                'the truncation decision must be visible in the trace.');
            tiss = tr.tissue(~cellfun(@isempty, tr.tissue));
            tc.verifyTrue(any(strcmp(tiss, 'CSF')) && any(strcmp(tiss, 'WM')), ...
                'the tissue class must be recorded per step.');
        end
    end
end

function [T, M] = evalc_track(nim, o)
    evalc('[T, M] = nim_tractography_hinec(nim, o);');
end
