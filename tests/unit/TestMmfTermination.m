classdef TestMmfTermination < matlab.unittest.TestCase
% Why the MMF tracker stopped, and the trace that reports it.
%
% valid_point collapsed four unrelated outcomes into one false - leaving the
% volume, leaving the propagation mask, landing outside the interpolation
% domain, and FA falling under the termination floor - and the caller labelled
% every one of them 'outside'. A traced run therefore reported 100% of
% streamlines ending 'outside' and the FA floor looked as though it were never
% applied. It was; the label hid it. The fix is label-only: tracks before and
% after are bit-identical.

    methods (Static)
        function [nim, o] = corridor()
        % Straight +y fibre field whose FA falls off a cliff partway along, so
        % the streamline must stop on the FA floor well inside the volume.
            d = [20 60 20];
            ev = zeros([d 3 3]); ev(:,:,:,2,1) = 1;
            fa = 0.6*ones(d); fa(:, 40:end, :) = 0.01;      % cliff at y = 40
            nim = struct('FA', fa, 'evec', ev, ...
                         'eval', repmat(reshape([3 1 1],1,1,1,3), [d 1]), 'mask', true(d));
            sm = false(d); sm(10,5,10) = true;
            o = struct('algorithm','mmf','field','dti','step_size',0.5,'max_steps',200, ...
                'termination_fa',0.1,'integration_order',4,'interp_method','cubic', ...
                'enable_diagnostics',false,'adaptive_step',false, ...
                'seed_mask',sm,'seed_density',1,'seed_strategy','uniform', ...
                'min_length',0,'max_arc',200,'angle_thresh',225, ...
                'wm_mask',[],'gm_mask',[],'csf_mask',[], ...
                'mmf_anchor',0.25,'trace',true,'trace_max',1);
        end
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
        end
    end

    methods (Test)

        function fallingBelowTheFaFloorIsReportedAsFa(tc)
            [nim, o] = TestMmfTermination.corridor();
            [T, info] = evalc_mmf(nim, o);
            tc.assumeNotEmpty(T, 'the tracker produced nothing to diagnose');
            tc.assertTrue(isfield(info,'trace') && ~isempty(info.trace), ...
                'the MMF tracer must record a trace when debug.trace is set');
            r = {info.trace(1).forward.termination, info.trace(1).backward.termination};
            tc.verifyTrue(any(strcmp(r,'fa')), ...
                sprintf(['stopping at an FA cliff must be reported as ''fa'', not ''outside'' ' ...
                         '- got {%s}.'], strjoin(r, ', ')));
        end

        function leavingTheVolumeIsStillOutside(tc)
            % The other branch must keep its own name, or the fix has merely
            % moved the conflation rather than removed it.
            [nim, o] = TestMmfTermination.corridor();
            nim.FA(:) = 0.6;                       % no cliff: run to the edge
            [~, info] = evalc_mmf(nim, o);
            tc.assertTrue(isfield(info,'trace') && ~isempty(info.trace));
            r = {info.trace(1).forward.termination, info.trace(1).backward.termination};
            tc.verifyTrue(any(strcmp(r,'outside')), ...
                sprintf('running out of the volume must still be ''outside'' - got {%s}.', ...
                        strjoin(r, ', ')));
        end

        function theTraceRecordsDriftAgainstTheField(tc)
            % Drift - the angle between the carried frame and the measured e1 -
            % is the diagnostic MMF needs and hinec has no counterpart for. On a
            % uniform field with the anchor on it must be near zero.
            [nim, o] = TestMmfTermination.corridor();
            [~, info] = evalc_mmf(nim, o);
            t = info.trace(1).forward;
            tc.assertTrue(isfield(t,'drift') && isfield(t,'aeff') && isfield(t,'kappa'), ...
                'the trace must record drift, the applied anchor and the curvature');
            dr = t.drift(~isnan(t.drift));
            tc.assumeNotEmpty(dr, 'no drift samples recorded');
            tc.verifyLessThan(median(dr), 5, ...
                'on a uniform straight field the carried frame must track the field closely');
        end
    end
end

function [T, info] = evalc_mmf(nim, o)
    % runTractography step 3: the connection geometry is built per run, before
    % the tracker, which now asserts it arrived rather than building it itself.
    % Mirror that order here or the test is not testing the pipeline.
    evalc('nim = nim_mmf_geometry(nim, o);');
    evalc('[T, info] = nim_tractography_mmf_connframe(nim, o);');
end
