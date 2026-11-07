classdef TestAdvanceFrameMl < matlab.unittest.TestCase
% advanceFrame_ml and nim_tractography_mmf_ml are copies of Sehun Chun's advanceFrame and
% nim_tractography_mmf with two additions: heading-aware fields (MMF-ml) and the history veto.
% With neither in use they must take exactly the steps of the originals, which stay unchanged.
% Fixture: the rigid-rotation phantom of TestMmfTracker (circles about the axis, FA 1 in an
% annulus), with a band of low linear shape CL = (l1 - l2)/l1 = 0.17 for the hold scheme to coast
% through.

    properties
        Nim
        Opts
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(fileparts(here));
            addpath(genpath(root));
            C = 30.5; n = 60; nz = 9;
            [X, Y, ~] = ndgrid(1:n, 1:n, 1:nz);
            dx = X - C; dy = Y - C; r = sqrt(dx.^2 + dy.^2); r(r < 1e-9) = 1e-9;
            evec = zeros(n, n, nz, 3, 3);
            evec(:,:,:,1,1) = -dy./r; evec(:,:,:,2,1) = dx./r;
            evec(:,:,:,1,2) =  dx./r; evec(:,:,:,2,2) = dy./r; evec(:,:,:,3,3) = 1;
            FA = zeros(n, n, nz); FA(r > 4 & r < 27) = 1;
            ev = repmat(reshape([3 1 1], 1, 1, 1, 3), n, n, nz);
            band = Y > C + 3 & Y < C + 7;                                % the zone
            ev(repmat(band, [1 1 1 3])) = repmat([1.2 1 1], nnz(band), 1);
            tc.Nim = struct('FA', FA, 'evec', evec, 'eval', ev, 'mask', FA > 0);
            seed = false(size(FA)); seed(round(C + 15), round(C) - 4 : round(C) + 1, 5) = true;
            cfg = load_config_yaml(fullfile(root, 'config', 'mmf_dti.yml'));
            tc.Opts = nim_config_to_options(cfg);
            tc.Opts.seed_mask = seed;
            tc.Opts.wm_mask = []; tc.Opts.gm_mask = []; tc.Opts.csf_mask = [];
            tc.Opts.step_size = 0.25; tc.Opts.max_arc = 60; tc.Opts.max_steps = 240;
            tc.Opts.min_length = 0; tc.Opts.seed_density = 1;
            tc.Opts.mmf_reanchor_tol = 0.2;
        end
    end

    methods (Test)
        function theTrackerCopyMatchesTheOriginal(tc)
            % MMF-4a, MMF-4h and MMF-4s on field: dti, veto off: identical tracks and metadata counts.
            setenv('HINEC_DISABLE_PARFOR', '1'); c = onCleanup(@() setenv('HINEC_DISABLE_PARFOR', ''));
            for hold = {'none', 'hold', 'march'}
                o = tc.Opts; o.mmf_hold = hold{1}; o.mmf_veto_arc = 0;
                evalc('[T0, m0] = nim_tractography_mmf(tc.Nim, o);');
                evalc('[T1, m1] = nim_tractography_mmf_ml(tc.Nim, o);');
                tc.verifyEqual(numel(T1), numel(T0), hold{1});
                tc.verifyTrue(isequal(T1, T0), ['tracks differ with mmf.hold ' hold{1}]);
                tc.verifyEqual(m1.mmf.coasted_fraction, m0.mmf.coasted_fraction, hold{1});
                tc.verifyEqual(m1.mmf.reanchored_fraction, m0.mmf.reanchored_fraction, hold{1});
                if ~strcmp(hold{1}, 'none')
                    tc.verifyGreaterThan(m0.mmf.coasted_fraction, 0, 'the band should be coasted');
                end
            end
        end

        function theSchemeCopyMatchesTheOriginal(tc)
            % advanceFrame vs advanceFrame_ml on an analytic planar field (circles of radius 1 about
            % the origin), with and without the hold scheme: identical tracts, frames and info.
            w  = @(X) circleW(X);
            mf = @(X) circleFrames(X);
            dom = struct('xlim', [-2 2], 'ylim', [-2 2]);
            seeds = [0.5 0; 0.8 0; 1.2 0];
            args = {{}, {'Mode', 'adaptive', 'ReanchorTol', 0.2}, ...
                    {'Mode', 'adaptive', 'ReanchorTol', 0.2, 'Confidence', @(X) 0.1 + (X(:, 2) < 0.3), ...
                     'HoldMin', 0.5, 'ExitTol', 1.2}};
            for a = 1:numel(args)
                [t0, f0, i0] = advanceFrame(w, mf, dom, seeds, [0 1], 0.05, 3, args{a}{:});
                [t1, f1, i1] = advanceFrame_ml(w, mf, dom, seeds, [0 1], 0.05, 3, args{a}{:});
                tc.verifyTrue(isequal(t1, t0) && isequal(f1, f0), sprintf('option set %d differs', a));
                for fn = fieldnames(i0)'     % every output of the original (NaN marks no data in both)
                    x0 = i0.(fn{1}); x1 = i1.(fn{1});
                    if isstruct(x0)          % the option record: the copy adds the veto's own options
                        x1 = rmfield(x1, setdiff(fieldnames(x1), fieldnames(x0)));
                    end
                    tc.verifyTrue(isequaln(x1, x0), sprintf('info.%s, option set %d', fn{1}, a));
                end
            end
        end
    end
end

function W = circleW(X)
% connection of the circles about the origin in the polar frame (e1 = tangent, e2 = inward normal):
% w_12 along e1 is 1/r, so W(:, 1, 2, :) = e1 / r and W(:, 2, 1, :) = -e1 / r.
r = max(vecnorm(X, 2, 2), 1e-9); t = [-X(:, 2), X(:, 1)] ./ r;
W = zeros(size(X, 1), 2, 2, 2);
W(:, 1, 2, :) = reshape(t ./ r, [], 1, 1, 2);
W(:, 2, 1, :) = -W(:, 1, 2, :);
end

function E = circleFrames(X)
r = max(vecnorm(X, 2, 2), 1e-9); t = [-X(:, 2), X(:, 1)] ./ r;
E = zeros(size(X, 1), 2, 2);
E(:, 1, :) = reshape(t, [], 1, 2); E(:, 2, :) = reshape(-X ./ r, [], 1, 2);
end
