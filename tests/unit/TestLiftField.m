classdef TestLiftField < matlab.unittest.TestCase
    % The lifts of crossing bundles on M = Omega x RP^2 (nim_lift_field) and the MMF tracker on them
    % (field: lift, nim_lift_connection). Fixture: make_lift_crossing_nim, the arc field of the MMF notes
    % (Sec. 9, bundle A) crossed by the logarithmic spiral of pitch alpha (bundle B) in a sector band,
    % noise-free so that the assertions test the construction, not a noise draw. The construction tests
    % switch the prior of few lifts off (lift.sparsity 0): it pulls the lobes of a crossing toward each other
    % by a few degrees, the price of removing false lifts on real data, which its own test measures.

    methods (TestClassSetup)
        function setup(tc) %#ok<MANU>
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
        end
    end

    methods
        function [nim, t] = crossing(~, alpha, fmax)
            [nim, t] = make_lift_crossing_nim('alpha', alpha, 'fmax', fmax, 'dims', [80 52 5], 'snr', Inf);
        end
        function a = liftError(~, L, u, sel)
            % angle from u to the voxel's nearest lift, over the voxels sel
            d = abs(sum(double(L.U) .* reshape(u, [size(sel) 1 3]), 5)); d(~L.ok) = 0;
            a = acosd(min(max(d, [], 4), 1)); a = a(sel);
        end
        function o = opts(~, nim, seeds, field)
            cfg = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti.yml'));
            o = nim_config_to_options(cfg);
            o.seed_mask = false(size(nim.FA));
            for s = 1:size(seeds, 1), p = round(seeds(s, :)); o.seed_mask(p(1), p(2), p(3)) = true; end
            o.seed_density = 1; o.wm_mask = []; o.gm_mask = []; o.csf_mask = [];
            o.min_length = 0; o.mmf_connection = 'nodal'; o.mmf_hold = 'none'; o.mmf_reanchor_tol = 0.2;
            o.step_size = 0.2; o.max_arc = 120; o.max_steps = 2000; o.termination_fa = 0.05; o.angle_thresh = 0;
            o.field = field;
        end
        function d = exitOffset(~, T, c, seed)
            % distance from the exact circle where the tract leaves the band (phi = 140 deg); Inf if it never does
            dd = cellfun(@(P) min(vecnorm(P - round(seed), 2, 2)), T); [~, k] = min(dd); P = T{k};
            r0 = hypot(round(seed(1)) - c(1), round(seed(2)) - c(2));
            ph = atan2d(P(:, 2) - c(2), P(:, 1) - c(1)); j = find(ph >= 140, 1);
            d = Inf; if ~isempty(j), d = abs(hypot(P(j, 1) - c(1), P(j, 2) - c(2)) - r0); end
        end
    end

    methods (Test)
        function aSingleBundleHasOneLiftAlongIt(tc)
            [nim, t] = tc.crossing(60, 0);
            [~, L] = evalc('nim_lift_field(nim, struct(''workers'', 0, ''sparsity'', 0))');
            m = logical(nim.mask);
            tc.verifyGreaterThan(mean(L.n(m) == 1), 0.99);
            tc.verifyLessThan(median(tc.liftError(L, t.uA, m)), 2);
        end

        function aCrossingHasOneLiftPerBundle(tc)
            [nim, t] = tc.crossing(60, 0.5);
            [~, L] = evalc('nim_lift_field(nim, struct(''workers'', 0, ''sparsity'', 0))');
            band = logical(nim.mask) & t.fB > 0.3;
            tc.verifyGreaterThan(mean(L.n(band) >= 2), 0.9);
            tc.verifyLessThan(median(tc.liftError(L, t.uA, band)), 4);
            tc.verifyLessThan(median(tc.liftError(L, t.uB, band)), 6);
            v1 = reshape(nim.evec(:, :, :, :, 1), [size(band) 3]);
            e1 = acosd(min(abs(sum(v1 .* t.uA, 4)), 1));
            tc.verifyGreaterThan(median(e1(band)), 10);          % the single-valued field is the mixture there
        end

        function theLiftedTrackerStaysOnItsBundle(tc)
            [nim, t] = tc.crossing(60, 0.7);
            c = t.c; seeds = [c(1) + [18 26]' * cosd(40), c(2) + [18 26]' * sind(40), [3 3]'];
            o = tc.opts(nim, seeds, 'dti');
            evalc('[T0, ~] = nim_tractography_mmf_ml(nim, o);');
            o.field = 'lift'; o.lift_workers = 0; o.lift_sparsity = 0;
            evalc('nimL = nim_field(nim, o, ''''); [T1, ~] = nim_tractography_mmf_ml(nimL, o);');
            for s = 1:2
                d0 = tc.exitOffset(T0, c, seeds(s, :));
                tc.verifyTrue(isfinite(d0));                                   % the dti tract does get through ...
                tc.verifyGreaterThan(d0, 2);                                   % ... carried off A's circle by v1
                tc.verifyLessThan(tc.exitOffset(T1, c, seeds(s, :)), 1);      % the lift does not
            end
        end

        function theStreamlineTrackerFollowsTheLiftsToo(tc)
            % algorithm: hinec reads the lifts as a peak field (the lift nearest the incoming tangent at every
            % stencil voxel, as for FOD peaks) and stays on bundle A's circle through the crossing, where the
            % same tracker on v1 is carried off it
            [nim, t] = tc.crossing(60, 0.7);
            c = t.c; seeds = [c(1) + [18 26]' * cosd(40), c(2) + [18 26]' * sind(40), [3 3]'];
            seed = false(size(nim.FA));
            for s = 1:2, p = round(seeds(s, :)); seed(p(1), p(2), p(3)) = true; end
            o = struct('field', 'dti', 'seed_mask', seed, 'seed_density', 1, 'seed_strategy', 'uniform', ...
                'step_size', 0.2, 'integration_order', 4, 'interp_method', 'trilinear', 'termination_fa', 0.05, ...
                'fa_threshold', 0.05, 'angle_thresh', 225, 'min_length', 0, 'max_steps', 2000, ...
                'enable_diagnostics', false, 'act_enabled', false, 'wm_mask', [], 'gm_mask', [], 'csf_mask', []);
            prior = getenv('HINEC_DISABLE_PARFOR'); cleanup = onCleanup(@() setenv('HINEC_DISABLE_PARFOR', prior)); %#ok<NASGU>
            setenv('HINEC_DISABLE_PARFOR', '1');
            evalc('[T0, ~] = nim_tractography_hinec(nim, o);');
            o.field = 'lift'; o.lift_workers = 0; o.lift_sparsity = 0;
            evalc('nimL = nim_field(nim, o, ''''); [T1, ~] = nim_tractography_hinec(nimL, o);');
            for s = 1:2
                d0 = tc.exitOffset(T0, c, seeds(s, :));
                tc.verifyTrue(isfinite(d0));
                tc.verifyGreaterThan(d0, 2);                                   % v1 carries it off A's circle
                tc.verifyLessThan(tc.exitOffset(T1, c, seeds(s, :)), 1);      % the lift does not
            end
        end

        function theFewLiftPriorRemovesFalseLiftsFromFreeWater(tc)
            % one bundle with 30% free water at SNR 20: without the prior the stick sum builds part of the
            % isotropic signal out of sticks, which become false lifts (no smoothing along the fibres, the
            % default: one lift in 36% of the voxels at mu = 0, 98.5% at mu = 0.1)
            [nim, t] = make_lift_crossing_nim('alpha', 60, 'fmax', 0, 'dims', [80 52 5], 'snr', 20, 'fiso', 0.3, 'seed', 3);
            m = logical(nim.mask);
            [~, L0] = evalc('nim_lift_field(nim, struct(''workers'', 0, ''sparsity'', 0))');
            [~, L1] = evalc('nim_lift_field(nim, struct(''workers'', 0, ''sparsity'', 0.1))');
            tc.verifyLessThan(mean(L0.n(m) == 1), 0.5);
            tc.verifyGreaterThan(mean(L1.n(m) == 1), 0.97);
            tc.verifyLessThan(median(tc.liftError(L1, t.uA, m)), median(tc.liftError(L0, t.uA, m)));
        end

        function theHoldCoastsWhereTheTractsOwnLiftIsMissing(tc)
            % bundle A's lift deleted over a strip of the crossing (phi 86..94 deg): its continuation is lost
            % there. Without the hold the tract has nothing to read and stops; with the hold, whose confidence
            % for field: lift is the read weight of the tract's own lift, it coasts on its held connection
            % and comes out on its circle.
            [nim, t] = tc.crossing(60, 0.5);
            c = t.c; seeds = [c(1) + [18 26]' * cosd(40), c(2) + [18 26]' * sind(40), [3 3]'];
            o = tc.opts(nim, seeds, 'lift'); o.lift_workers = 0; o.lift_sparsity = 0;
            evalc('nimL = nim_field(nim, o, '''');');
            [X, Y, ~] = ndgrid(1:size(nim.FA, 1), 1:size(nim.FA, 2), 1:size(nim.FA, 3));
            strip = abs(atan2d(Y - c(2), X - c(1)) - 90) <= 4;
            for s = 1:size(nimL.lift_ok, 4)
                u = reshape(nimL.lift_U(:, :, :, s, :), [size(strip) 3]);
                isA = abs(sum(double(u) .* t.uA, 4)) >= cosd(20);
                nimL.lift_ok(:, :, :, s) = nimL.lift_ok(:, :, :, s) & ~(strip & isA);
            end
            evalc('[T0, ~] = nim_tractography_mmf_ml(nimL, o);');
            o.mmf_hold = 'hold'; o.mmf_hold_min = 0.25; o.mmf_exit_tol = 1.2; o.lift_read_min = 0.125;
            evalc('[T1, ~] = nim_tractography_mmf_ml(nimL, o);');
            for s = 1:2
                tc.verifyEqual(tc.exitOffset(T0, c, seeds(s, :)), Inf);       % stops in the strip
                tc.verifyLessThan(tc.exitOffset(T1, c, seeds(s, :)), 1);      % coasts through it
            end
        end

        function theSheetReadRecoversTheCurvatureOfOneBundle(tc)
            % the connection of the local-sheet read on an exact one-bundle field (the arcs of the notes,
            % u_A = phi-hat, kappa = 1 / rho): the tangent and kappa = |omega(e1)| = 1 / rho, for one query point as
            % for many (a batch of one is the case MATLAB's logical indexing gets wrong), and on a single sheet the
            % same read as the corner read. (The lifts of nim_lift_field are not exact: their refinement on the
            % dictionary of directions makes the field turn in small steps, kappa * rho 0.73-0.97 noise-free.)
            [nim, t] = tc.crossing(60, 0);
            U = reshape(t.uA, [size(nim.FA) 1 3]); ok = reshape(logical(nim.mask), [size(nim.FA) 1]);
            c = t.c; rho = [16 20 24 28]'; ph = 60;
            X = [c(1) + rho * cosd(ph), c(2) + rho * sind(ph), 3 * ones(4, 1)]; H = repmat([-sind(ph) cosd(ph) 0], 4, 1);
            kap = @(W) sqrt(sum(squeeze(sum(W(:, 1, 2, :) .* reshape(H(1:size(W, 1), :), [], 1, 1, 3), 4)) .^ 2 + ...
                                squeeze(sum(W(:, 1, 3, :) .* reshape(H(1:size(W, 1), :), [], 1, 1, 3), 4)) .^ 2, 2));
            for interp = {'sheet_linear', 'sheet_cubic'}
                [w, mf] = nim_lift_connection(U, ok, 0.125, 30, 'fd4', interp{1});
                k = kap(w(X, H)); E = mf(X, H);
                tc.verifyTrue(all(isfinite(k)), interp{1});
                tc.verifyLessThan(max(abs(k .* rho - 1)), 0.01, interp{1});            % kappa within 1% of 1 / rho
                tc.verifyLessThan(max(acosd(min(abs(sum(squeeze(E(:, 1, :)) .* H, 2)), 1))), 1, interp{1});
                k1 = kap(w(X(1, :), H(1, :)));                                           % a batch of one point
                tc.verifyEqual(k1, k(1), 'AbsTol', 1e-12, interp{1});
            end
            w0 = nim_lift_connection(U, ok, 0.125, 30, 'fd4', 'corner');
            w1 = nim_lift_connection(U, ok, 0.125, 30, 'fd4', 'sheet_linear');
            tc.verifyEqual(kap(w1(X, H)), kap(w0(X, H)), 'RelTol', 1e-6);                % one sheet: the same read
        end

        function theCurvedKernelReachesTheLiftField(tc)
            % every lift.* key of the schema must reach the options nim_field passes to nim_lift_field
            cfg = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti.yml'));
            cfg.tractography.lift.smooth = 2; cfg.tractography.lift.smooth_kernel = 'curved';
            cfg.tractography.lift.smooth_kappa = [0.2 0.4]; cfg.tractography.lift.smooth_normals = 6;
            o = nim_config_to_options(cfg);
            for k = fieldnames(cfg.tractography.lift)'
                tc.verifyTrue(isfield(o, ['lift_' k{1}]), ['lift.' k{1} ' does not reach the options']);
            end
            tc.verifyEqual(char(string(o.lift_smooth_kernel)), 'curved');
            tc.verifyEqual(o.lift_smooth_kappa(:)', [0.2 0.4]);
            tc.verifyEqual(o.lift_smooth_normals, 6);
        end

        function dtiIsTheDefaultField(tc)
            cfg = load_config_yaml(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'mmf_dti.yml'));
            o = nim_config_to_options(cfg);
            tc.verifyEqual(lower(char(string(o.field))), 'dti');
        end
    end
end
