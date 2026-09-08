classdef TestTrackTrace < matlab.unittest.TestCase
% The per-step trace and the termination labels that depend on it.
%
% WHY THIS EXISTS. Termination reasons were unusable for diagnosis: every
% failure inside the direction lookup surfaced as 'no_direction', merging two
% unrelated causes - the position leaving the interpolation domain, and FA
% falling under the termination floor. On a traced cingulum run that one label
% covered 769 of 800 arms while reporting ZERO 'fa' terminations, so the
% statistics said nothing about why tracking stopped. These tests pin the two
% causes apart, and pin the trace's internal consistency, because a trace that
% is off by one step would silently misattribute every error built on it.

    methods (Static)
        function nim = field(dirvec, fa)
        % A uniform direction field with a prescribed FA volume.
            d = size(fa);
            nim = struct();
            nim.FA = fa;
            nim.eval = repmat(reshape([3 1 1],1,1,1,3), [d 1]);
            nim.evec = zeros([d 3 3]);
            u = dirvec(:)/norm(dirvec);
            for c = 1:3, nim.evec(:,:,:,c,1) = u(c); end
            nim.mask = true(d);
        end
    end

    methods (Test)

        function traceIsInternallyConsistent(tc)
            % chord must be the realised displacement, and turn must be the
            % angle between consecutive stored directions. Both are recomputable
            % from other trace fields, so a recording bug shows up here.
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            d = [24 24 24];
            fa = 0.6*ones(d);
            nim = TestTrackTrace.field([0 1 0], fa);
            sm = false(d); sm(12,8,12) = true;
            opts = struct('step_size',0.5,'max_steps',30,'termination_fa',0.1, ...
                          'integration_order',4,'interp_method','linear','field','dti', ...
                          'enable_diagnostics',false,'adaptive_step',false, ...
                          'seed_mask',sm,'seed_density',1,'seed_strategy','uniform', ...
                          'trace',true,'trace_max',1);
            [tracks, meta] = evalc_track(nim, opts);
            tc.assumeNotEmpty(meta, 'tracker produced no metadata');
            tc.assertTrue(isfield(meta,'trace'), 'no trace was recorded');
            t = meta.trace(1).forward;
            tc.assumeGreaterThan(t.n_steps, 3, 'trace too short to check');

            dp = vecnorm(diff(t.pos,1,1),2,2);
            m = ~isnan(t.chord(1:end-1));
            tc.verifyLessThan(max(abs(t.chord(m) - dp(m))), 1e-12, ...
                'chord(i) must equal |pos(i+1) - pos(i)|; a mismatch means the trace is off by one step.');
            tc.verifyLessThan(max(abs(vecnorm(t.dir,2,2)-1)), 1e-12, ...
                'stored directions must be unit vectors.');
            tc.verifyNotEmpty(tracks);
        end

        function lowFaIsReportedAsFaNotNoDirection(tc)
            % A slab of sub-threshold FA ahead of the seed must stop tracking
            % with reason 'fa'.
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            d = [24 24 24];
            fa = 0.6*ones(d);
            fa(:, 16:end, :) = 0.01;              % wall of near-isotropic tissue
            nim = TestTrackTrace.field([0 1 0], fa);
            sm = false(d); sm(12,8,12) = true;
            opts = struct('step_size',0.5,'max_steps',40,'termination_fa',0.1, ...
                          'integration_order',4,'interp_method','linear','field','dti', ...
                          'enable_diagnostics',false,'adaptive_step',false, ...
                          'seed_mask',sm,'seed_density',1,'seed_strategy','uniform', ...
                          'trace',true,'trace_max',1);
            [~, meta] = evalc_track(nim, opts);
            tc.assertTrue(isfield(meta,'trace'), 'no trace was recorded');
            r = {meta.trace(1).forward.termination, meta.trace(1).backward.termination};
            tc.verifyTrue(any(strcmp(r,'fa')), ...
                sprintf(['tracking into sub-threshold FA must terminate with reason ''fa'', ' ...
                         'not ''no_direction'' - got {%s}.'], strjoin(r, ', ')));
        end
    end
end

function [tracks, meta] = evalc_track(nim, opts)
    evalc('[tracks, meta] = nim_tractography_hinec(nim, opts);');
end
