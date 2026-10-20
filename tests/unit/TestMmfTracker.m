classdef TestMmfTracker < matlab.unittest.TestCase
% Moving frames tractography (nim_tractography_mmf, advanceFrame) on fields
% whose fibres are known exactly, plus the tracker contract, the option
% plumbing and the termination logic.
%
% The rigid-rotation phantom of TestIntegratorOrder (circles of every radius,
% FA 1 in an annulus) is where the exponential step is at its best: along a
% circle the curvature is constant, so once the interpolated dyadic reproduces
% the tangent field the step is exact for any step size and the residual is
% the interpolation error of the field (spline 5e-6 voxel over 60 voxels of
% arc) or, with the cubic and trilinear kernels, the error of integrating a
% connection that is discontinuous at the voxel faces (7e-4, 1.6e-2). The
% helix fixture of make_synthetic_nim has torsion, where the twist is dropped
% and the step is exact only to the order of the integrator. The straight
% field is the degenerate case (kappa = 0, a line).

    properties
        Nim
        Opts
        R0      % radius of the circle through the seed voxel centre
        C = 30.5
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            C = tc.C; n = 60; nz = 9;
            [X, Y, ~] = ndgrid(1:n, 1:n, 1:nz);
            dx = X - C; dy = Y - C; r = sqrt(dx.^2 + dy.^2); r(r < 1e-9) = 1e-9;
            evec = zeros(n, n, nz, 3, 3);
            evec(:,:,:,1,1) = -dy./r; evec(:,:,:,2,1) = dx./r;
            evec(:,:,:,1,2) =  dx./r; evec(:,:,:,2,2) = dy./r; evec(:,:,:,3,3) = 1;
            FA = zeros(n, n, nz); FA(r > 4 & r < 27) = 1;
            tc.Nim = struct('FA', FA, 'evec', evec, ...
                'eval', repmat(reshape([3 1 1], 1, 1, 1, 3), n, n, nz), 'mask', FA > 0);
            seed = false(size(FA)); seed(round(C + 15), round(C), 5) = true;
            tc.R0 = hypot(round(C + 15) - C, round(C) - C);   % 15.508: the seed voxel centre is off the radius-15 circle

            cfg = load_config_yaml(fullfile(fileparts(fileparts(here)), 'config', 'mmf_dti.yml'));
            tc.Opts = nim_config_to_options(cfg);
            tc.Opts.seed_mask = seed;
            tc.Opts.wm_mask = []; tc.Opts.gm_mask = []; tc.Opts.csf_mask = [];
            tc.Opts.step_size = 0.25; tc.Opts.max_arc = 60; tc.Opts.max_steps = 240;
            tc.Opts.min_length = 0; tc.Opts.seed_density = 1;
        end
    end

    methods
        function err = circleError(tc, T, R0)
            % largest distance of the track's nodes from the circle of radius R0
            err = max(abs(sqrt(sum((T(:, 1:2) - tc.C).^2, 2)) - R0));
        end
        function nim = straightNim(~, dims)
            nim = struct('FA', ones(dims), 'mask', true(dims), 'evec', zeros([dims 3 3]));
            nim.evec(:,:,:,3,1) = 1; nim.evec(:,:,:,1,2) = 1; nim.evec(:,:,:,2,3) = 1;   % fibres along z
        end
        function c = countOf(~, meta, reason)
            j = strcmp(meta.termination.reason, reason);
            if any(j), c = meta.termination.count(j); else, c = 0; end
        end
    end

    methods (Test)
        function followsTheCircleAndHonoursTheContract(tc)
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, tc.Opts);');
            evalc('nim_check_tracker_output(T, meta);');
            tc.assertNumElements(T, 1);
            tc.verifyEqual(meta.n_seeds, 1);
            tc.verifyEqual(meta.seed_index, 1);
            tc.verifyEqual(size(meta.seed_points), [1 3]);
            tr = T{1};
            % The seed is on the track once and joins the two halves.
            seedrow = find(all(abs(tr - meta.seed_points) < 1e-12, 2));
            tc.verifyNumElements(seedrow, 1);
            tc.verifyEqual(seedrow, 241);                       % 240 backward nodes before it
            % Both halves ran to max_arc: 60/0.25 = 240 steps each.
            tc.verifyEqual(size(tr, 1), 481);
            tc.verifyEqual(meta.termination.reason, {'max_arc'});
            tc.verifyEqual(meta.termination.count, 2);
            % Consecutive points are the chord of an arc of length h, i.e. h to
            % O(h^3 kappa^2) = 1e-6 here: the nodes are the integration nodes.
            steps = sqrt(sum(diff(tr, 1, 1).^2, 2));
            tc.verifyEqual(steps, 0.25 * ones(size(steps)), 'AbsTol', 1e-4);
            % The field is planar: z never moves.
            tc.verifyLessThan(max(abs(tr(:, 3) - 5)), 1e-9);
            % On the circle through the seed the exponential step is exact up to
            % the interpolation error of the dyadic field: 5.3e-6 voxel with the
            % spline kernel of mmf_dti.yml (0.48 for the Euler template over
            % the same 60 voxels of arc).
            tc.verifyLessThan(tc.circleError(tr, tc.R0), 2e-5);
            tc.verifyEqual(meta.mmf.order, 4);
            tc.verifyEqual(meta.mmf.mode, 'adaptive');
            tc.verifyEqual(meta.mmf.seeds_alive, 1);
            tc.verifyEmpty(meta.mmf.seeds_skipped.reason);
        end

        function signAlignmentAcrossTheAxisFlipLocus(tc)
            % principal_axis returns the axis with a sign that flips across the
            % diagonals of this phantom (dx*dy > 0); the connection must align
            % the signs of its central difference there. A seed on that
            % diagonal runs the whole circle across the flip locus.
            o = tc.Opts; o.seed_mask = false(size(tc.Nim.FA)); o.seed_mask(41, 41, 5) = true;
            evalc('[T, ~] = nim_tractography_mmf(tc.Nim, o);');
            tc.assertNumElements(T, 1);
            tc.verifyLessThan(tc.circleError(T{1}, hypot(41 - tc.C, 41 - tc.C)), 2e-5);
        end

        function straightFieldIsALine(tc)
            dims = [20 20 60];
            nim = tc.straightNim(dims);
            o = tc.Opts; o.seed_mask = false(dims); o.seed_mask(10, 11, 30) = true;
            o.max_arc = 25; o.max_steps = 100;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.assertNumElements(T, 1);
            tr = T{1};
            tc.verifyLessThan(max(abs(tr(:, 1) - 10)), 1e-12);
            tc.verifyLessThan(max(abs(tr(:, 2) - 11)), 1e-12);
            tc.verifyEqual(sqrt(sum(diff(tr, 1, 1).^2, 2)), 0.25 * ones(size(tr, 1) - 1, 1), 'AbsTol', 1e-12);
            tc.verifyEqual(meta.termination.reason, {'max_arc'});
            tc.verifyEqual(meta.termination.count, 2);
        end

        function kernelsAndModesAreDistinguishable(tc)
            % Spline: every mode and order lands within the interpolation error.
            for mode = {'transport', 'reanchor', 'adaptive'}
                for integ = {'euler', 'rk2', 'rk4'}
                    o = tc.Opts; o.mmf_mode = mode{1}; o.integrator = integ{1};
                    evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
                    tc.assertNumElements(T, 1, sprintf('%s/%s gave no track', mode{1}, integ{1}));
                    tc.verifyLessThan(tc.circleError(T{1}, tc.R0), 2e-5, sprintf('%s/%s left the circle', mode{1}, integ{1}));
                    tc.verifyEqual(meta.mmf.mode, mode{1});
                    tc.verifyEqual(meta.mmf.order, 2^(find(strcmp(integ{1}, {'euler', 'rk2', 'rk4'})) - 1));
                end
            end
            % Cubic and trilinear: the connection is discontinuous at the voxel
            % faces, so transport integrates a kinked right-hand side (cubic
            % 7e-4, trilinear 1.6e-2 with CF4; the measured hierarchy spline <
            % cubic < trilinear), the orders separate (cubic Lie-Euler 3.7e-3
            % against 7e-4), and re-anchoring to the data after every step
            % brings the track back to the interpolated circle (trilinear
            % 1.4e-4). The adaptive mode keeps the transported heading here
            % (nothing drifts by 1e-3 rad) and coincides with transport.
            o = tc.Opts; o.mmf_mode = 'transport';
            e = struct();
            for kern = {'spline', 'cubic', 'trilinear'}
                for integ = {'euler', 'rk4'}
                    o.interp_method = kern{1}; o.integrator = integ{1};
                    evalc('[T, ~] = nim_tractography_mmf(tc.Nim, o);');
                    e.(sprintf('%s_%s', kern{1}, integ{1})) = tc.circleError(T{1}, tc.R0);
                end
            end
            tc.verifyLessThan(e.spline_rk4, e.cubic_rk4);
            tc.verifyLessThan(e.cubic_rk4, e.trilinear_rk4);
            tc.verifyLessThan(e.cubic_rk4, 2e-3);
            tc.verifyGreaterThan(e.cubic_euler, 2 * e.cubic_rk4);
            tc.verifyLessThan(e.trilinear_rk4, 5e-2);
            o.interp_method = 'trilinear'; o.integrator = 'rk4';
            o.mmf_mode = 'reanchor';  evalc('[Tr, mr] = nim_tractography_mmf(tc.Nim, o);');
            o.mmf_mode = 'adaptive';  evalc('[Ta, ma] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyLessThan(tc.circleError(Tr{1}, tc.R0), 1e-3);
            tc.verifyGreaterThan(e.trilinear_rk4, 5 * tc.circleError(Tr{1}, tc.R0));
            tc.verifyEqual(mr.mmf.reanchored_fraction, 1);
            tc.verifyEqual(ma.mmf.reanchored_fraction, 0);
            o.mmf_mode = 'transport'; evalc('[Tt, ~] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(Ta{1}, Tt{1}, 'AbsTol', 1e-12);
        end

        function helixIsFollowedToTheIntegratorOrder(tc)
            % Constant curvature AND torsion (make_synthetic_nim 'helix': radius
            % a = 15, pitch b = 6, kappa = a/(a^2+b^2)). The tracker drops the
            % twist ('Twist' 'zero'), so Omega is not constant along the fibre
            % in the carried frame and the step is exact only to the order of
            % the integrator: measured 1.0e-3 / 3.5e-4 / 5.0e-6 voxel for
            % orders 1 / 2 / 4 over 40 voxels of arc each way (spline).
            nim = make_synthetic_nim('helix');
            o = tc.Opts; o.seed_mask = false(size(nim.FA)); o.seed_mask(45, 30, 30) = true;
            o.max_arc = 40; o.max_steps = 160; o.trace = true; o.trace_max = 1;
            err = zeros(1, 3); integs = {'euler', 'rk2', 'rk4'};
            for j = 1:3
                o.integrator = integs{j};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 1, [integs{j} ' gave no track']);
                rho = sqrt(sum((T{1}(:, 1:2) - 30).^2, 2));
                err(j) = max(abs(rho - nim.truth.a));
                tc.verifyGreaterThan(max(T{1}(:, 3)) - min(T{1}(:, 3)), 25, 'the track does not advance along the helix axis');
            end
            tc.verifyLessThan(err(3), 2e-5);
            tc.verifyGreaterThan(err(2), 5 * err(3));
            tc.verifyGreaterThan(err(1), 2 * err(2));
            % |kappa| returned with the CF4 track is the curvature of the helix.
            tc.verifyEqual(median(meta.trace(1).forward.kappa(2:end)), nim.truth.kappa, 'AbsTol', 1e-6);
        end

        function adaptiveReanchorsOnlyWhenTheHeadingDrifts(tc)
            % On the exact-ish circle the carried heading agrees with the data
            % to far better than the default tolerance, so nothing re-anchors;
            % with tolerance 0 every step does.
            evalc('[~, m1] = nim_tractography_mmf(tc.Nim, tc.Opts);');
            tc.verifyEqual(m1.mmf.reanchored_fraction, 0);
            o = tc.Opts; o.mmf_reanchor_tol = 0;
            evalc('[T2, m2] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(m2.mmf.reanchored_fraction, 1);
            tc.verifyLessThan(tc.circleError(T2{1}, tc.R0), 2e-5);
        end

        function rejectsWhatItCannotDo(tc)
            o = tc.Opts; o.integrator = 'rkf45';
            tc.verifyError(@() nim_tractography_mmf(tc.Nim, o), 'mmf:integrator');
            o = tc.Opts; o.field = 'csd';
            tc.verifyError(@() nim_tractography_mmf(tc.Nim, o), 'mmf:field');
        end

        function stopsAtTheMaskAndReportsWhy(tc)
            % Both halves leave the allowed half-plane y > 30 and both say so,
            % with every integrator: the scheme examines the stages of a step,
            % and a step end it never examined (orders 1, 2) is dropped by the
            % tracker when it lies outside tissue.
            nim = tc.Nim; nim.mask(:, 1:30, :) = false;
            o = tc.Opts; o.propagation_mask = nim.mask;
            for integ = {'euler', 'rk2', 'rk4'}
                o.integrator = integ{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 1);
                tc.verifyTrue(all(round(T{1}(:, 2)) >= 31), [integ{1} ': a node lies in the blocked half-plane']);
                tc.verifyEqual(meta.termination.reason, {'outside_mask'}, integ{1});
                tc.verifyEqual(meta.termination.count, 2, integ{1});
            end
        end

        function propagationMaskTakesPrecedenceOverTheBrainMask(tc)
            nim = tc.Nim; nim.mask(:, 1:30, :) = false;
            o = tc.Opts; o.propagation_mask = tc.Nim.mask;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.assertNumElements(T, 1);
            tc.verifyEqual(size(T{1}, 1), 481);
            tc.verifyEqual(meta.termination.reason, {'max_arc'});
        end

        function reportsVolumeAndFaStops(tc)
            dims = [20 20 60];
            nim = tc.straightNim(dims);
            o = tc.Opts; o.seed_mask = false(dims); o.seed_mask(10, 11, 30) = true;
            o.max_arc = 100; o.max_steps = 400;
            for integ = {'euler', 'rk4'}
                o.integrator = integ{1};
                % Straight to both faces of the volume: 'outside' twice, every
                % node inside the volume.
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.verifyEqual(meta.termination.reason, {'outside'}, integ{1});
                tc.verifyEqual(meta.termination.count, 2, integ{1});
                tc.verifyTrue(all(T{1}(:, 3) >= 1 & T{1}(:, 3) <= 60), integ{1});
                % An FA cliff at z = 46: the forward half stops with 'fa' and no
                % node lies where the interpolated FA is below the floor (0.1,
                % i.e. beyond z = 45.9).
                nf = nim; nf.FA(:, :, 46:end) = 0;
                evalc('[T, meta] = nim_tractography_mmf(nf, o);');
                tc.verifyEqual(tc.countOf(meta, 'fa'), 1, integ{1});
                tc.verifyEqual(tc.countOf(meta, 'outside'), 1, integ{1});
                tc.verifyLessThan(max(T{1}(:, 3)), 45.9, integ{1});
            end
        end

        function angleBudgetStopsASharpTurn(tc)
            % angle_max is degrees per voxel of arc; with the circle turning at
            % 57.3/15.5 = 3.7 deg per voxel, a budget of 2 deg/voxel stops the
            % tract at its seed (the first step would exceed it), a budget of
            % 10 lets it run.
            o = tc.Opts; o.angle_thresh = 2;
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEmpty(T);
            tc.verifyEqual(meta.termination.reason, {'angle'});
            tc.verifyEqual(meta.termination.count, 2);
            tc.verifyEqual(meta.mmf.max_turn_rad, deg2rad(0.5));
            o.angle_thresh = 10;
            evalc('[T, ~] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyNumElements(T, 1);
        end

        function reanchoringCannotBypassTheAngleBudget(tc)
            % Fibres along y for y <= 30 and along x beyond: the interpolated
            % principal axis is piecewise constant with a swap at y = 30.5, so
            % the connection is zero except within gradient_delta of the swap
            % plane and a step whose stages straddle it (h = 0.4 from a voxel
            % centre: nodes 30.4 -> 30.8, stages at 30.4, 30.6, 30.8) sees no
            % curvature at all. The data's e1 at the end of that step is 90
            % degrees from the carried one; re-anchoring to it would put a 90
            % degree kink into the track, past any budget. Every mode must stop
            % with 'angle' instead: the budget bounds the turn of the heading as
            % the data see it (the data's e1 at the end of a step against the
            % heading at its start), which the re-anchoring modes would realise
            % and which transport, without re-anchoring, would run through.
            dims = [60 60 6];
            nim = struct('FA', ones(dims), 'mask', true(dims), 'evec', zeros([dims 3 3]));
            nim.evec(:, 1:30, :, 2, 1) = 1;   nim.evec(:, 1:30, :, 1, 2) = 1;
            nim.evec(:, 31:end, :, 1, 1) = 1; nim.evec(:, 31:end, :, 2, 2) = 1;
            nim.evec(:, :, :, 3, 3) = 1;
            o = tc.Opts; o.seed_mask = false(dims); o.seed_mask(30, 20, 3) = true;
            o.step_size = 0.4; o.max_arc = 40; o.max_steps = 100; o.angle_thresh = 100;   % 40 deg per step
            unit = @(v) v ./ sqrt(sum(v.^2, 2));
            kink = @(t) max([0; acos(min(1, max(-1, sum(unit(diff(t(1:end-1, :), 1, 1)) .* unit(diff(t(2:end, :), 1, 1)), 2))))]);
            for mode = {'adaptive', 'reanchor', 'transport'}
                o.mmf_mode = mode{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 1, mode{1});
                tc.verifyEqual(tc.countOf(meta, 'angle'), 1, mode{1});
                tc.verifyEqual(tc.countOf(meta, 'outside'), 1, mode{1});    % the backward half, at y = 1
                tc.verifyEqual(max(T{1}(:, 2)), 30.4, 'AbsTol', 1e-9, mode{1});
                tc.verifyLessThan(kink(T{1}), 1e-9, mode{1});
            end
        end

        function traceIsRecordedPerStep(tc)
            o = tc.Opts; o.trace = true; o.trace_max = 1;
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
            tc.assertTrue(isfield(meta, 'trace'));
            tr = T{1}; seedrow = 241;
            for fn = {'direction', 'pos', 'dir', 'fa', 'h', 'chord', 'turn', 'flip', 'drift', 'reanchored', 'kappa', 's', 'n_steps', 'termination'}
                tc.verifyTrue(isfield(meta.trace(1).forward, fn{1}), ['missing field ' fn{1}]);
            end
            f = meta.trace(1).forward; b = meta.trace(1).backward;
            tc.verifyEqual(meta.trace(1).seed_index, 1);
            tc.verifyEqual(f.direction, 'forward'); tc.verifyEqual(b.direction, 'backward');
            % The records are the two halves of the track, seed first.
            tc.verifyEqual(f.n_steps, 240); tc.verifyEqual(b.n_steps, 240);
            tc.verifyEqual(f.pos, tr(seedrow:end, :));
            tc.verifyEqual(b.pos, flipud(tr(1:seedrow, :)));
            tc.verifyEqual(f.termination, 'max_arc'); tc.verifyEqual(b.termination, 'max_arc');
            % The carried direction is the unit tangent; the step bookkeeping
            % is what nim_track_diagnose reads: nominal h, realised chord,
            % interpolated FA, turn in degrees, no sign flips.
            tc.verifyEqual(sqrt(sum(f.dir.^2, 2)), ones(241, 1), 'AbsTol', 1e-12);
            tc.verifyEqual(f.h, [0.25 * ones(240, 1); NaN]);
            tc.verifyEqual(f.chord(1:240), sqrt(sum(diff(f.pos, 1, 1).^2, 2)));
            tc.verifyTrue(isnan(f.chord(end)) && isnan(f.turn(end)));
            tc.verifyEqual(f.fa, ones(241, 1), 'AbsTol', 1e-12);
            tc.verifyEqual(median(f.turn(1:240)), rad2deg(0.25 / tc.R0), 'AbsTol', 1e-3);
            tc.verifyTrue(islogical(f.flip) && ~any(f.flip));
            tc.verifyTrue(islogical(f.reanchored) && ~any(f.reanchored));
            tc.verifyEqual(f.s, 0.25 * (0:240)', 'AbsTol', 1e-12);
            % |kappa| returned with the tract is the curvature of the circle,
            % 1/R0 (measured error 1e-8); the heading never drifts (4e-8 rad).
            tc.verifyEqual(median(f.kappa(2:end)), 1/tc.R0, 'AbsTol', 1e-6);
            tc.verifyLessThan(max(f.drift), 1e-6);
        end

        function seedsAreTheLatticeOfTheOtherTrackersAndChunksDoNotMatter(tc)
            % 4001 seeds over three slices cross the chunk boundary (4000):
            % tracks are the same as from a single seed, the contract holds.
            FA = tc.Nim.FA; n = size(FA, 1);
            cand = find(FA(:, :, 4:6) > 0); [i, j, k] = ind2sub([n n 3], cand(1:4001));
            sm = false(size(FA)); sm(sub2ind(size(FA), i, j, k + 3)) = true;
            o = tc.Opts; o.seed_mask = sm; o.max_arc = 5; o.max_steps = 20;
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
            evalc('nim_check_tracker_output(T, meta);');
            tc.verifyEqual(meta.n_seeds, 4001);
            tc.verifyEqual(numel(T), 4001);
            tc.verifyEqual(meta.seed_index, 1:4001);
            sl = find(sm);
            for s = [1 4000 4001]
                o1 = o; o1.seed_mask = false(size(FA)); o1.seed_mask(sl(s)) = true;
                evalc('[T1, ~] = nim_tractography_mmf(tc.Nim, o1);');
                tc.verifyEqual(T{s}, T1{1}, sprintf('seed %d differs from its single-seed run', s));
            end
            % Density 8 on two voxels: the sub-voxel lattice of nim_seed_offsets,
            % seed_points aligned with the tracks, and min_length applied to the
            % arc length: on the circle through (46, 31) a mask at y > 40 cuts
            % the forward half after ~10 voxels (the backward half runs its 60
            % before coming round), ~70 voxels in all; the circle through
            % (36, 31) runs the full 120.
            nim = tc.Nim; nim.mask(:, 41:end, :) = false;
            o = tc.Opts; o.seed_density = 8; o.propagation_mask = nim.mask;
            o.seed_mask = false(size(FA)); o.seed_mask(46, 31, 5) = true; o.seed_mask(36, 31, 5) = true;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(meta.n_seeds, 16);
            tc.verifyEqual(numel(T), 16);
            off = nim_seed_offsets(8);
            expected = [repmat([46 31 5], 8, 1) + off; repmat([36 31 5], 8, 1) + off];
            tc.verifyEqual(sortrows(meta.seed_points), sortrows(expected), 'AbsTol', 1e-12);
            for t = 1:16
                tc.verifyTrue(any(all(abs(T{t} - meta.seed_points(t, :)) < 1e-12, 2)), 'seed_points row is not on its track');
            end
            arc = cellfun(@(t) sum(sqrt(sum(diff(t, 1, 1).^2, 2))), T);
            tc.verifyEqual(nnz(arc > 100), 8); tc.verifyEqual(nnz(arc < 80), 8);
            o.min_length = 100;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(numel(T), 8);
            tc.verifyEqual(size(meta.seed_index), [1 8]);
            tc.verifyTrue(all(abs(meta.seed_points(:, 1) - 36) <= 0.25 + 1e-12));
            tc.verifyEqual(meta.n_seeds, 16);
            % Random seeding jitters as hinec does: +-0.2 voxel at density 1,
            % +-0.4 above.
            o = tc.Opts; o.seed_strategy = 'random'; o.max_arc = 2; o.max_steps = 8;
            o.seed_mask = false(size(FA)); o.seed_mask(46, 31, 5) = true; o.seed_mask(45, 31, 5) = true;
            rng(1);
            evalc('[~, m1] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(m1.n_seeds, 2);
            tc.verifyLessThanOrEqual(max(abs(m1.seed_points(:) - round(m1.seed_points(:)))), 0.2);
            o.seed_density = 3;
            evalc('[~, m3] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(m3.n_seeds, 6);
            d = max(abs(m3.seed_points(:) - round(m3.seed_points(:))));
            tc.verifyLessThanOrEqual(d, 0.4); tc.verifyGreaterThan(d, 0.2);
        end

        function legacyOptionsWithoutMaxArcStillRun(tc)
            % runTractography(data, 'mmf') hands over the legacy option struct,
            % which has max_steps but no max_arc.
            o = rmfield(tc.Opts, 'max_arc'); o.max_steps = 240;
            evalc('[T, ~] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(size(T{1}, 1), 481);
        end

        function deadSeedGivesNoTrackButIsAccountedFor(tc)
            o = tc.Opts; o.seed_mask = false(size(tc.Nim.FA)); o.seed_mask(30, 31, 5) = true;   % FA 0, outside the mask
            o.trace = true; o.trace_max = 0;
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
            evalc('nim_check_tracker_output(T, meta);');
            tc.verifyEqual(size(T), [0 1]);
            tc.verifyEqual(meta.n_seeds, 1);
            tc.verifyEqual(size(meta.seed_index), [1 0]);
            tc.verifyEqual(size(meta.seed_points), [0 3]);
            tc.verifyEqual(meta.mmf.seeds_alive, 0);
            tc.verifyEqual(meta.mmf.seeds_skipped.reason, {'outside_mask'});
            tc.verifyEqual(meta.mmf.seeds_skipped.count, 1);
            tc.verifyEmpty(meta.termination.reason); tc.verifyEmpty(meta.termination.count);
            tc.verifyEqual(meta.trace(1).forward.n_steps, 0);
            tc.verifyEqual(meta.trace(1).forward.termination, 'outside_mask');
            tc.verifyEqual(meta.trace(1).backward.pos, [30 31 5]);
        end

        function crossingStopsEveryModeOnTheAngleBudget(tc)
            % make_synthetic_nim('crossing'): a band along x (y in [27, 33]) is
            % crossed by a band along y (x in [40, 46]) whose direction the
            % crossing voxels carry. Seeds at sub-voxel phases (density 8), so
            % that no stage lands within gradient_delta of the swap plane
            % x = 39.5: the connection sees no curvature there, and only the
            % data's e1 at the end of a step, 90 degrees from the heading, can
            % stop the tract. hinec stops on its node-to-node angle; every
            % mode of the frame tracker must stop as well.
            nim = make_synthetic_nim('crossing');
            o = tc.Opts; o.seed_mask = false(size(nim.FA)); o.seed_mask(10, 30, 3) = true;
            o.seed_density = 8; o.step_size = 0.2; o.max_arc = 100; o.max_steps = 500;
            o.angle_thresh = 225;                                   % 45 deg per step
            for mode = {'adaptive', 'reanchor', 'transport'}
                o.mmf_mode = mode{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.verifyEqual(meta.n_seeds, 8);
                tc.assertNumElements(T, 8, mode{1});
                tc.verifyEqual(tc.countOf(meta, 'angle'), 8, mode{1});   % every forward half, at x ~ 39.5
                tc.verifyTrue(all(cellfun(@(t) max(t(:, 1)), T) < 40), mode{1});
                tc.verifyTrue(all(cellfun(@(t) max(abs(t(:, 2) - 30)), T) < 1), mode{1});   % never turned into the other band
            end
        end

        function disabledAngleBudgetKeepsASafetyCeiling(tc)
            % angle_max 0 switches the criterion off (MaxTurn Inf) in every
            % tracker. A stage landing within gradient_delta of a jump of e1
            % then meets |kappa| ~ 1/delta, and the exponential step would
            % rotate the frame by an arbitrary angle and retrace the track; the
            % tracker keeps a ceiling of pi per step there and reports the stop
            % as 'no_direction' (the data are singular), not as 'angle'. Here
            % the step from the voxel centre with h = 0.25 puts a CF4 stage on
            % the swap plane y = 30.5 of the two-bundle field.
            dims = [60 60 6];
            nim = struct('FA', ones(dims), 'mask', true(dims), 'evec', zeros([dims 3 3]));
            nim.evec(:, 1:30, :, 2, 1) = 1;   nim.evec(:, 1:30, :, 1, 2) = 1;
            nim.evec(:, 31:end, :, 1, 1) = 1; nim.evec(:, 31:end, :, 2, 2) = 1;
            nim.evec(:, :, :, 3, 3) = 1;
            o = tc.Opts; o.seed_mask = false(dims); o.seed_mask(30, 20, 3) = true;
            o.step_size = 0.25; o.max_arc = 40; o.max_steps = 160; o.angle_thresh = 0;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.assertNumElements(T, 1);
            tc.verifyEqual(meta.mmf.max_turn_rad, Inf);
            tc.verifyEqual(tc.countOf(meta, 'angle'), 0);
            tc.verifyEqual(tc.countOf(meta, 'no_direction'), 1);
            tc.verifyEqual(tc.countOf(meta, 'outside'), 1);
            tc.verifyLessThan(max(T{1}(:, 2)), 30.5);
            tc.verifyTrue(all(diff(T{1}(:, 2)) > 0), 'the track retraces itself');
        end

        function nodesOnAVolumeFaceAreNotDead(tc)
            % The difference stencil is one-sided at a face of the volume, so a
            % seed on the face x = 1 is tracked along it.
            dims = [20 20 60];
            nim = tc.straightNim(dims);
            o = tc.Opts; o.seed_mask = false(dims); o.seed_mask(1, 11, 30) = true;
            o.max_arc = 25; o.max_steps = 100;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(meta.mmf.seeds_alive, 1);
            tc.assertNumElements(T, 1);
            tc.verifyEqual(size(T{1}, 1), 201);
            tc.verifyTrue(all(T{1}(:, 1) == 1));
            tc.verifyEqual(meta.termination.reason, {'max_arc'});
        end

        function aNonFiniteVoxelDoesNotPoisonTheSplineField(tc)
            % griddedInterpolant 'spline' is global along grid lines: one NaN
            % voxel anywhere would make the dyadic field NaN everywhere.
            nim = tc.Nim; nim.evec(1, 1, 1, :, :) = NaN;      % far outside the annulus
            evalc('[T, meta] = nim_tractography_mmf(nim, tc.Opts);');
            tc.verifyEqual(meta.mmf.seeds_alive, 1);
            tc.assertNumElements(T, 1);
            tc.verifyLessThan(tc.circleError(T{1}, tc.R0), 2e-5);
        end
    end
end
