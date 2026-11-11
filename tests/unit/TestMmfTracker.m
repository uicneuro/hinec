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
            % The tests below pin the 'derivative' connection (central difference of the
            % interpolated field), whose tolerances they were written for; the 'nodal'
            % connection of the notes, the default, has its own tests at the end.
            tc.Opts.mmf_connection = 'derivative';
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

        function angleMaxStopsASharpTurn(tc)
            % angle_max is degrees per voxel of arc; with the circle turning at
            % 57.3/15.5 = 3.7 deg per voxel, an angle_max of 2 deg/voxel stops
            % the tract at its seed (the first step would exceed it), an
            % angle_max of 10 lets it run.
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

        function reanchoringCannotBypassAngleMax(tc)
            % Fibres along y for y <= 30 and along x beyond: the interpolated
            % principal axis is piecewise constant with a swap at y = 30.5, so
            % the connection is zero except within gradient_delta of the swap
            % plane and a step whose stages straddle it (h = 0.4 from a voxel
            % centre: nodes 30.4 -> 30.8, stages at 30.4, 30.6, 30.8) sees no
            % curvature at all. The data's e1 at the end of that step is 90
            % degrees from the carried one; re-anchoring to it would put a 90
            % degree kink into the track, past any angle_max. Every mode must stop
            % with 'angle' instead: angle_max bounds the turn of the heading as
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

        function crossingStopsEveryModeOnAngleMax(tc)
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

        function disabledAngleMaxKeepsASafetyCeiling(tc)
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
    methods (Test)
        function nodalConnectionIsTheDefault(tc)
            here = fileparts(mfilename('fullpath'));
            cfg = load_config_yaml(fullfile(fileparts(fileparts(here)), 'config', 'mmf_dti.yml'));
            o = nim_config_to_options(cfg);
            tc.verifyEqual(o.mmf_connection, 'nodal');
            tc.verifyEqual(o.mmf_max_turn, 0);
        end
        function nodalConnectionFollowsTheCircle(tc)
            % Connection computed at the voxel nodes (fourth-order differences of the dyadic)
            % and interpolated: the circle is kept to the node-level accuracy, 2.6e-4 voxel
            % over 60 voxels of arc with the spline kernel, and the curvature returned along
            % the tract is 1/R0 to 0.1%.
            o = tc.Opts; o.mmf_connection = 'nodal'; o.interp_method = 'spline';
            o.trace = true; o.trace_max = 0;
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
            tc.assertNumElements(T, 1);
            tc.verifyEqual(meta.mmf.connection, 'nodal');
            tc.verifyLessThan(tc.circleError(T{1}, tc.R0), 1e-3);
            k = meta.trace(1).forward.kappa; k = k(isfinite(k));
            tc.verifyLessThan(max(abs(k * tc.R0 - 1)), 2e-3);
        end
        function nodalConnectionIsSignFree(tc)
            % Flipping the arbitrary sign of v1 in half of the voxels changes nothing.
            o = tc.Opts; o.mmf_connection = 'nodal';
            evalc('T1 = nim_tractography_mmf(tc.Nim, o);');
            nim = tc.Nim; flip = rand(size(nim.FA)) > 0.5;
            for c = 1:3
                v = nim.evec(:, :, :, c, 1); v(flip) = -v(flip); nim.evec(:, :, :, c, 1) = v;
            end
            evalc('T2 = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(size(T2{1}), size(T1{1}));
            tc.verifyLessThan(max(abs(T2{1}(:) - T1{1}(:))), 1e-9);
        end
        function maxTurnOptionReplacesAngleMax(tc)
            % mmf.max_turn (radians per step) replaces the angle_max mapping: a max_turn of 0.01 rad
            % stops the circle (h*kappa = 0.25/15.5 = 0.016 rad) at its seed.
            o = tc.Opts; o.mmf_connection = 'nodal'; o.mmf_max_turn = 0.01;
            evalc('[T, meta] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyTrue(isempty(T) || size(T{1}, 1) <= 3);
            tc.verifyEqual(meta.mmf.max_turn_rad, 0.01);
        end
    end
    methods
        function o = holdOpts(tc, scheme)
            % The options of config/mmf_dti_hold.yml and mmf_dti_march.yml on the phantoms.
            o = tc.Opts; o.mmf_connection = 'nodal'; o.mmf_hold = scheme; o.mmf_reanchor_tol = 0.2;
            o.step_size = 0.2; o.max_arc = 40; o.max_steps = 200;
        end
        function nim = degenerate(~, nim, zone)
            % The voxels of zone become crossings fitted as pancakes: CL = 0 and the
            % principal axis the second eigenvector, 90 degrees from the fibres.
            pancake = [2 2 1];
            for c = 1:3
                v1 = nim.evec(:, :, :, c, 1); v2 = nim.evec(:, :, :, c, 2);
                t = v1; v1(zone) = v2(zone); v2(zone) = t(zone);
                nim.evec(:, :, :, c, 1) = v1; nim.evec(:, :, :, c, 2) = v2;
                e = nim.eval(:, :, :, c); e(zone) = pancake(c); nim.eval(:, :, :, c) = e;
            end
        end
        function [nim, r, phi] = spiralNim(tc, alpha)
            % The phantom with the fibres at the angle alpha from the radial direction:
            % logarithmic spirals about the centre, kappa = sin(alpha)/r, sigma = cos(alpha)/r
            % (the circles of the class for alpha = pi/2).
            dims = size(tc.Nim.FA);
            [X, Y, ~] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
            dx = X - tc.C; dy = Y - tc.C; r = sqrt(dx.^2 + dy.^2); phi = atan2(dy, dx);
            e1x = (cos(alpha) * dx - sin(alpha) * dy) ./ r;
            e1y = (cos(alpha) * dy + sin(alpha) * dx) ./ r;
            nim = tc.Nim;
            nim.evec(:, :, :, 1, 1) = e1x;  nim.evec(:, :, :, 2, 1) = e1y;
            nim.evec(:, :, :, 1, 2) = -e1y; nim.evec(:, :, :, 2, 2) = e1x;
        end
        function e = spiralError(tc, T, p, alpha)
            % Largest distance of the track's nodes from the spiral through p: the gap
            % along the circle about the centre, whose normal component is cos(alpha) of it.
            invariant = @(P) atan2(P(:, 2) - tc.C, P(:, 1) - tc.C) - tan(alpha) * log(hypot(P(:, 1) - tc.C, P(:, 2) - tc.C));
            d = invariant(T) - invariant(p);
            e = max(abs(atan2(sin(d), cos(d))) .* hypot(T(:, 1) - tc.C, T(:, 2) - tc.C)) * cos(alpha);
        end
    end
    methods (Test)
        function holdChangesNothingWithoutAZone(tc)
            % mmf_dti sets it off (the schema default is march). On data with CL above the floor everywhere the zone is
            % empty, so the scheme, its fill of the zone and its stencils change nothing.
            tc.verifyEqual(tc.Opts.mmf_hold, 'none');
            o = tc.holdOpts('none');
            evalc('[T0, m0] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(m0.mmf.hold, 'none');
            tc.verifyEqual(m0.mmf.coasted_fraction, 0);
            for scheme = {'hold', 'march'}
                o.mmf_hold = scheme{1};
                evalc('[T1, m1] = nim_tractography_mmf(tc.Nim, o);');
                tc.verifyEqual(T1{1}, T0{1});
                tc.verifyEqual(m1.mmf.hold, scheme{1});
                tc.verifyEqual(m1.mmf.coasted_fraction, 0);
            end
        end
        function theCrossingOfTheFixtureIsCrossed(tc)
            % make_synthetic_nim('crossing') at the sub-voxel phases of the angle_max test: every
            % mode reads the swap of the axis and stops at x ~ 39.5 without the scheme. With
            % it the tract coasts straight through the pancakes, re-anchors on band A beyond
            % them and leaves the volume, never turning into band B.
            nim = make_synthetic_nim('crossing');
            o = tc.holdOpts('none');
            o.seed_mask = false(size(nim.FA)); o.seed_mask(10, 30, 3) = true; o.seed_density = 8;
            o.max_arc = 100; o.max_steps = 500; o.angle_thresh = 225;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'angle'), 8);
            tc.verifyTrue(all(cellfun(@(t) max(t(:, 1)), T) < 40));
            for scheme = {'hold', 'march'}
                o.mmf_hold = scheme{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 8, scheme{1});
                tc.verifyEqual(meta.termination.reason, {'outside'}, scheme{1});
                tc.verifyTrue(all(cellfun(@(t) max(t(:, 1)), T) > 79), scheme{1});
                tc.verifyLessThan(max(cellfun(@(t) max(abs(t(:, 2) - 30)), T)), 0.3, scheme{1});
                tc.verifyLessThan(max(cellfun(@(t) max(abs(diff(t(:, 2)) ./ diff(t(:, 1)))), T)), 1e-3, scheme{1});
                tc.verifyGreaterThan(meta.mmf.coasted_fraction, 0.05, scheme{1});
                evalc('nim_check_tracker_output(T, meta);');
            end
        end
        function theHoldStaysOnTheCircleThroughADegenerateSector(tc)
            % Circles, kappa constant: the case the hold is exact on, and sigma = 0, so
            % that the march is the hold. The sector turns the axis by 90 degrees over about
            % 4 to 8 voxels of arc. The seeds are 32 sub-voxel phases at four radii; what is
            % left (0.04 voxel median, 0.10 at worst) is the staircase of the voxelized edge
            % of the zone.
            [nim, ~, phi] = tc.spiralNim(pi/2);
            nim = tc.degenerate(nim, phi > deg2rad(60) & phi < deg2rad(80));
            seed = false(size(nim.FA));
            seed(sub2ind(size(seed), round(tc.C + (12:3:21)), round(tc.C) * ones(1, 4), 5 * ones(1, 4))) = true;
            o = tc.holdOpts('none'); o.seed_mask = seed; o.seed_density = 8;
            evalc('[~, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'angle'), 32);               % the half that reaches the sector first
            for scheme = {'hold', 'march'}
                o.mmf_hold = scheme{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 32, scheme{1});
                tc.verifyEqual(meta.termination.reason, {'max_arc'}, scheme{1});
                e = zeros(32, 1);
                for k = 1:32
                    e(k) = tc.circleError(T{k}, hypot(meta.seed_points(k, 1) - tc.C, meta.seed_points(k, 2) - tc.C));
                end
                tc.verifyLessThan(median(e), 0.06, scheme{1});
                tc.verifyLessThan(max(e), 0.15, scheme{1});
            end
        end
        function theMarchFollowsASpiralWhereTheHoldDoesNot(tc)
            % kappa/sigma constant: the march is exact and the hold is not, through an annulus
            % of 3 voxels crossed over 6 voxels of arc at alpha = 60 degrees, outwards from
            % seeds inside it, where the fan opens, and inwards from seeds beyond it, where
            % it closes: 0.02 voxel median and 0.05, 0.07 at worst for the march, against
            % 0.21 and 0.25 outwards and 0.09 and 0.17 inwards for the hold.
            alpha = deg2rad(60);
            [nim, r] = tc.spiralNim(alpha);
            nim = tc.degenerate(nim, r > 15 & r < 18);
            o = tc.holdOpts('hold'); o.seed_density = 8; o.max_arc = 26;
            radii = {8:2:12, 20:2:24};
            e = struct();
            for way = 1:2
                seed = false(size(nim.FA));
                seed(sub2ind(size(seed), round(tc.C + radii{way}), round(tc.C) * ones(1, 3), 5 * ones(1, 3))) = true;
                o.seed_mask = seed;
                for scheme = {'hold', 'march'}
                    o.mmf_hold = scheme{1};
                    evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                    tc.assertNumElements(T, 24, scheme{1});
                    e.(scheme{1}){way} = arrayfun(@(k) tc.spiralError(T{k}, meta.seed_points(k, :), alpha), (1:24)');
                end
                tc.verifyLessThan(median(e.march{way}), 0.03);
                tc.verifyLessThan(max(e.march{way}), 0.1);
                tc.verifyGreaterThan(median(e.hold{way}), 4 * median(e.march{way}));
            end
        end
        function theMarchKeepsTheSignsOfTheAxialData(tc)
            % The fill of the zone and the nodal tensor are free of the arbitrary sign of v1
            % in half of the voxels. (The direction of travel in the block of the march is
            % pinned by the spiral, crossed both ways.)
            alpha = deg2rad(60);
            [nim, r] = tc.spiralNim(alpha);
            nim = tc.degenerate(nim, r > 15 & r < 18);
            seed = false(size(nim.FA)); seed(round(tc.C + 10), round(tc.C), 5) = true;
            o = tc.holdOpts('march'); o.seed_mask = seed; o.max_arc = 26;
            evalc('T1 = nim_tractography_mmf(nim, o);');
            flip = rand(size(nim.FA)) > 0.5;
            for c = 1:3
                v = nim.evec(:, :, :, c, 1); v(flip) = -v(flip); nim.evec(:, :, :, c, 1) = v;
            end
            evalc('T2 = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(size(T2{1}), size(T1{1}));
            tc.verifyLessThan(max(abs(T2{1}(:) - T1{1}(:))), 1e-9);
        end
        function theMarchTurnsTheCurvatureOnAHelix(tc)
            % The helix of the fixture (radius 15, pitch 6) through a degenerate sector of
            % 40 degrees, about 11 voxels of arc: the twist is not a datum, so the held
            % curvature keeps its direction and the coast is a circle; the transverse block of
            % the connection turns it by the torsion, the law of the march in 3-D. At the end of
            % the coast the tract is 0.35 from its helix with the hold and 0.04 with the march.
            nim = make_synthetic_nim('helix');
            d = size(nim.FA);
            [X, Y, ~] = ndgrid(1:d(1), 1:d(2), 1:d(3));
            ctr = d(1:2) / 2;
            phi = atan2(Y - ctr(2), X - ctr(1));
            for c = 1:3                                                  % an orthonormal completion of the fixture's e1
                nim.evec(:, :, :, c, 2) = 0; nim.evec(:, :, :, c, 3) = 0;
            end
            rad = cat(4, cos(phi), sin(phi), zeros(d));
            e1  = nim.evec(:, :, :, :, 1);
            nim.evec(:, :, :, :, 2) = rad;
            nim.evec(:, :, :, :, 3) = cross(e1, rad, 4);
            nim = tc.degenerate(nim, phi > deg2rad(20) & phi < deg2rad(60));
            a = nim.truth.a; b = nim.truth.b;
            seed = false(d); seed(round(ctr(1) + a), round(ctr(2)), 30) = true;
            o = tc.holdOpts('hold'); o.seed_mask = seed; o.max_arc = 40;
            o.seed_density = 1; o.trace = true; o.trace_max = 0;
            gap = struct();
            for scheme = {'hold', 'march'}
                o.mmf_hold = scheme{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 1, scheme{1});
                f = meta.trace(1).forward;
                tc.assertGreaterThan(nnz(f.coasted), 20, scheme{1});
                last = find(f.coasted, 1, 'last');
                p  = f.pos(last, :);  p0 = f.pos(1, :);
                ph = atan2(p(2) - ctr(2), p(1) - ctr(1)) - atan2(p0(2) - ctr(2), p0(1) - ctr(1));
                gap.(scheme{1}) = hypot(hypot(p(1) - ctr(1), p(2) - ctr(2)) - hypot(p0(1) - ctr(1), p0(2) - ctr(2)), ...
                                        (p(3) - p0(3)) - b * ph);
            end
            tc.verifyLessThan(gap.march, 0.2 * gap.hold);
            tc.verifyLessThan(gap.march, 0.1);
        end
        function aCoastThatLeavesTheBundleStops(tc)
            % Beyond the pancakes band A runs 30 degrees off, 0.52 rad: within the
            % limit per step from angle_max, and above an exit_tol of 0.2.
            nim = make_synthetic_nim('crossing');
            nim.evec(47:end, 27:33, :, 1, 1) = cosd(30); nim.evec(47:end, 27:33, :, 2, 1) = sind(30);
            nim.evec(47:end, 27:33, :, 1, 2) = -sind(30); nim.evec(47:end, 27:33, :, 2, 2) = cosd(30);
            o = tc.holdOpts('hold'); o.seed_mask = false(size(nim.FA)); o.seed_mask(30, 30, 3) = true;
            o.seed_density = 8; o.max_arc = 100; o.max_steps = 500; o.angle_thresh = 225; o.mmf_exit_tol = 0.2;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'exit'), 8);
            tc.verifyTrue(all(cellfun(@(t) max(t(:, 1)), T) < 50));
            tc.verifyEqual(meta.mmf.exit_tol, 0.2);
            o.mmf_exit_tol = 1;                                          % the same tract with a tolerance above the 0.52 rad
            evalc('[~, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'exit'), 0);
        end
        function aCoastLongerThanTheCapStops(tc)
            nim = make_synthetic_nim('crossing');
            o = tc.holdOpts('hold'); o.seed_mask = false(size(nim.FA)); o.seed_mask(30, 30, 3) = true;
            o.max_arc = 100; o.max_steps = 500; o.angle_thresh = 225; o.mmf_coast_max = 2; o.trace = true; o.trace_max = 0;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'coast'), 1);
            f = meta.trace(1).forward;
            tc.verifyEqual(f.termination, 'coast');
            tc.verifyLessThan(nnz(f.coasted) * o.step_size, 2.5);
            tc.verifyLessThan(max(T{1}(:, 1)), 43);
            tc.verifyEqual(meta.mmf.coast_max, 2);
        end
        function theMarchStopsAtTheFocusOfAFan(tc)
            % Rays through a point in 3-D, kappa = 0, B = I/r. Inwards the fan closes and the
            % density the march carries, 1/u = (r_0/r)^2, reaches ten times its value at
            % r/r_0 = 0.32: 'focus'. The hold has no fan and runs into the hole around the
            % centre, the FA floor.
            d = [41 41 41]; ctr = 21;
            [X, Y, Z] = ndgrid(1:d(1), 1:d(2), 1:d(3));
            x = X - ctr; y = Y - ctr; z = Z - ctr; r = sqrt(x.^2 + y.^2 + z.^2); r(r == 0) = 1;
            nim = struct('FA', double(r > 1.5 & r < 19.5), 'mask', r < 19.5, 'evec', zeros([d 3 3]), ...
                         'eval', repmat(reshape([3 1 1], 1, 1, 1, 3), d));
            e1 = cat(4, x ./ r, y ./ r, z ./ r);
            ref = zeros([d 3]); ref(:, :, :, 1) = 1; flipx = abs(x ./ r) > 0.8; ref(:, :, :, 1) = ~flipx; ref(:, :, :, 2) = flipx;
            e2 = ref - sum(ref .* e1, 4) .* e1; e2 = e2 ./ sqrt(sum(e2 .^ 2, 4));
            nim.evec(:, :, :, :, 1) = e1; nim.evec(:, :, :, :, 2) = e2; nim.evec(:, :, :, :, 3) = cross(e1, e2, 4);
            nim = tc.degenerate(nim, r < 11);
            seed = false(d); seed(ctr + 9, ctr + 9, ctr + 5) = true;     % r = 13.7
            o = tc.holdOpts('march'); o.seed_mask = seed; o.max_arc = 30; o.max_steps = 150; o.mmf_coast_max = 20;
            o.trace = true; o.trace_max = 0;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.assertNumElements(T, 1);
            tc.verifyEqual(tc.countOf(meta, 'focus'), 1);
            rec = meta.trace(1).forward;
            if ~strcmp(rec.termination, 'focus'), rec = meta.trace(1).backward; end
            k = find(rec.coasted, 1);
            rn = @(i) norm(rec.pos(i, :) - ctr);
            tc.verifyGreaterThan(rn(size(rec.pos, 1)) / rn(k), 0.27);
            tc.verifyLessThan(rn(size(rec.pos, 1)) / rn(k), 0.42);
            o.mmf_hold = 'hold';
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'focus'), 0);
            tc.verifyEqual(tc.countOf(meta, 'fa'), 1);
            tc.verifyLessThan(min(sqrt(sum((T{1} - ctr) .^ 2, 2))), 2.5);
        end
        function theFaFloorAndTheMaskStayStopsInTheZone(tc)
            % A slab across the sector below the FA floor, then outside the mask: the coast
            % does not read the data, but a tract that goes out of tissue stops there, and
            % for that reason.
            [nim, ~, phi] = tc.spiralNim(pi/2);
            nim = tc.degenerate(nim, phi > deg2rad(60) & phi < deg2rad(80));
            cut = phi > deg2rad(68) & phi < deg2rad(72);
            seed = false(size(nim.FA)); seed(round(tc.C + 15), round(tc.C), 5) = true;
            o = tc.holdOpts('hold'); o.seed_mask = seed; o.termination_fa = 0.5;
            n1 = nim; n1.FA(cut) = 0.1;
            evalc('[T, meta] = nim_tractography_mmf(n1, o);');
            tc.verifyEqual(tc.countOf(meta, 'fa'), 1);                   % the half that reaches the sector
            FAI = griddedInterpolant(n1.FA, 'linear');
            tc.verifyGreaterThanOrEqual(min(FAI(T{1})), o.termination_fa);
            n2 = nim; n2.mask(cut) = false; n2.FA = tc.Nim.FA;
            evalc('[T, meta] = nim_tractography_mmf(n2, o);');
            tc.verifyEqual(tc.countOf(meta, 'outside_mask'), 1);
            vox = sub2ind(size(n2.mask), round(T{1}(:, 1)), round(T{1}(:, 2)), round(T{1}(:, 3)));
            tc.verifyTrue(all(n2.mask(vox)));
        end
        function theEdgeOfTheTissueIsNotAZone(tc)
            % A straight bundle of CL 0.22, just above the floor, between fitted, isotropic
            % voxels, at sub-voxel phases beside them. Whether those voxels are below the FA
            % floor or outside the mask, their low CL must reach the tract neither through
            % the erosion nor through the interpolation in the cell, or it coasts on perfect
            % data. The zone is empty in both.
            d = [30 30 80];
            nim = struct('FA', 0.05 * ones(d), 'mask', true(d), 'evec', zeros([d 3 3]), 'eval', zeros([d 3]));
            nim.evec(:, :, :, 3, 1) = 1; nim.evec(:, :, :, 1, 2) = 1; nim.evec(:, :, :, 2, 3) = 1;
            in = false(d); in(13:17, 13:17, :) = true;
            bundle = [1 0.78 0.7]; grey = [0.8 0.76 0.72];
            for c = 1:3
                e = grey(c) * ones(d); e(in) = bundle(c); nim.eval(:, :, :, c) = e;
            end
            nim.FA(in) = 0.8;
            seed = false(d); seed(13:17, 13:17, 40) = true;
            o = tc.holdOpts('hold'); o.seed_mask = seed; o.seed_density = 8;
            o.termination_fa = 0.1; o.max_arc = 100; o.max_steps = 500;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyNumElements(T, 200);
            tc.verifyEqual(meta.mmf.coasted_fraction, 0);
            tc.verifyEqual(meta.mmf.zone_fraction, 0);
            tc.verifyEqual(meta.termination.reason, {'outside'});
            nim.mask = in; nim.FA(:) = 0.8;                              % the same voxels above the floor, outside the mask
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(meta.mmf.coasted_fraction, 0);
            tc.verifyEqual(meta.mmf.zone_fraction, 0);
        end
        function theRingTurnsWithTheFrame(tc)
            % Circles in planes tilted from every plane of the completion of frames_at, and
            % two sectors a few trusted nodes apart. A re-anchoring replaces the normal plane
            % by the data's arbitrary one, so the ring of the second coast has to turn with
            % it: the mean of two gauges would turn the curvature out of the circle's plane.
            d = [64 64 64]; c = [32.5 32.5 32.5]; n = [0.6 0 0.8];
            [X, Y, Z] = ndgrid(1:d(1), 1:d(2), 1:d(3)); P = [X(:) Y(:) Z(:)] - c;
            a = cross(n, [0.3 0.5 0.8]); a = a / norm(a); b = cross(n, a);
            h = P * n'; rho = P - h * n; r = sqrt(sum(rho .^ 2, 2)); r(r == 0) = 1; er = rho ./ r;
            e1 = cross(repmat(n, size(P, 1), 1), er, 2);
            th = atan2(P * b', P * a');
            FA = reshape(double(abs(h) < 7 & r > 8 & r < 29), d);
            z = abs(h) < 7 & ((th > deg2rad(40) & th < deg2rad(55)) | (th > deg2rad(59) & th < deg2rad(74)));
            ev1 = e1; ev1(z, :) = er(z, :); ev2 = er; ev2(z, :) = e1(z, :);
            ev = zeros([d 3 3]);
            for q = 1:3
                ev(:, :, :, q, 1) = reshape(ev1(:, q), d);
                ev(:, :, :, q, 2) = reshape(ev2(:, q), d);
                ev(:, :, :, q, 3) = n(q);
            end
            L = repmat([3 1 1], size(P, 1), 1); L(z, :) = repmat([2 2 1], nnz(z), 1);
            nim = struct('FA', FA, 'mask', FA > 0, 'evec', ev, 'eval', reshape(L, [d 3]));
            seed = false(d);
            for rr = 15:2:21
                q = round(c + rr * a); seed(q(1), q(2), q(3)) = true;
            end
            o = tc.holdOpts('hold'); o.seed_mask = seed; o.seed_density = 8;
            o.mmf_coast_max = 20; o.max_arc = 60; o.max_steps = 300;
            bound = struct('hold', 1e-6, 'march', 0.1);                  % the block of the march reads the voxels
            for scheme = {'hold', 'march'}
                o.mmf_hold = scheme{1};
                evalc('[T, meta] = nim_tractography_mmf(nim, o);');
                tc.assertNumElements(T, 32, scheme{1});
                tc.verifyEqual(tc.countOf(meta, 'exit'), 0, scheme{1});
                out = zeros(32, 1);
                for k = 1:32
                    out(k) = max(abs((T{k} - meta.seed_points(k, :)) * n'));
                end
                tc.verifyLessThan(max(out), bound.(scheme{1}), scheme{1});   % 0.4 to 0.6 with the ring left behind
            end
        end
        function theExitDelayLengthensTheCoast(tc)
            nim = make_synthetic_nim('crossing');
            o = tc.holdOpts('hold'); o.seed_mask = false(size(nim.FA)); o.seed_mask(30, 30, 3) = true;
            o.max_arc = 100; o.max_steps = 500; o.angle_thresh = 225; o.trace = true; o.trace_max = 0;
            evalc('[~, m0] = nim_tractography_mmf(nim, o);');
            o.mmf_exit_delay = 2;
            evalc('[~, m2] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(nnz(m2.trace(1).forward.coasted) - nnz(m0.trace(1).forward.coasted), 10);
            tc.verifyEqual(m2.mmf.exit_delay, 2);
            % It is no way out of the tissue: here the band ends two voxels beyond the crossing.
            nim.FA(49:end, :, :) = 0; nim.mask(49:end, :, :) = 0;
            o.mmf_exit_delay = 6;
            evalc('[T, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'outside_mask'), 1);
            vox = sub2ind(size(nim.mask), round(T{1}(:, 1)), round(T{1}(:, 2)), round(T{1}(:, 3)));
            tc.verifyTrue(all(nim.mask(vox) > 0));
        end
        function theTraceRecordsTheCoast(tc)
            nim = make_synthetic_nim('crossing');
            o = tc.holdOpts('march'); o.seed_mask = false(size(nim.FA)); o.seed_mask(30, 30, 3) = true;
            o.max_arc = 100; o.max_steps = 500; o.angle_thresh = 225; o.trace = true; o.trace_max = 0;
            evalc('[~, meta] = nim_tractography_mmf(nim, o);');
            f = meta.trace(1).forward;
            tc.verifyEqual(size(f.coasted), [size(f.pos, 1) 1]);
            tc.verifyEqual(size(f.confidence), [size(f.pos, 1) 1]);
            tc.verifyTrue(islogical(f.coasted));
            x = f.pos(f.coasted, 1);
            tc.verifyGreaterThan(min(x), 38.5); tc.verifyLessThan(max(x), 48);
            tc.verifyEqual(diff(find(f.coasted)), ones(nnz(f.coasted) - 1, 1));
            tc.verifyTrue(all(f.confidence(f.coasted) < 0.2));
            tc.verifyTrue(all(f.confidence(~f.coasted) >= 0.2));
            o.mmf_hold = 'none';
            evalc('[~, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyFalse(any(meta.trace(1).forward.coasted));
            tc.verifyTrue(all(isnan(meta.trace(1).forward.confidence)));
        end
        function theSchemeNeedsTheEigenvaluesAndAnIntegerRing(tc)
            o = tc.holdOpts('hold');
            tc.verifyError(@() nim_tractography_mmf(rmfield(tc.Nim, 'eval'), o), 'mmf:hold');
            o.mmf_hold = 'zero';
            tc.verifyError(@() nim_tractography_mmf(tc.Nim, o), 'mmf:hold');
            o.mmf_hold = 'march'; o.mmf_hold_k = 2.5;
            tc.verifyError(@() nim_tractography_mmf(tc.Nim, o), 'mmf:hold');
        end
        function aZoneThatIsTheTissueItselfIsReported(tc)
            % 1000 voxels of tissue, 1000 below the floor, n of the tissue pancakes: the share
            % is of the tissue, and the warning starts above 40%.
            d = [10 10 20];
            nim = tc.straightNim(d);
            nim.FA(1001:end) = 0.01;
            nim.eval = repmat(reshape([3 1 1], 1, 1, 1, 3), d);
            o = tc.holdOpts('hold'); o.termination_fa = 0.1; o.max_arc = 2; o.max_steps = 10;
            o.seed_mask = false(d); o.seed_mask(5, 5, 8) = true;
            for n = [400 401]
                lastwarn('');
                z = false(d); z(1:n) = true;
                m = tc.degenerate(nim, z);
                evalc('[~, meta] = nim_tractography_mmf(m, o);');
                [~, id] = lastwarn;
                tc.verifyEqual(meta.mmf.zone_fraction, n / 1000, 'AbsTol', 1e-12);
                tc.verifyEqual(strcmp(id, 'mmf:holdZone'), n > 400);
            end
        end
        function theFloorAndTheRingReachTheScheme(tc)
            % On the spiral a coast from one trusted node is not one from three, and a floor
            % above the CL of the phantom makes everything a zone: a straight coast to the cap.
            alpha = deg2rad(60);
            [nim, r] = tc.spiralNim(alpha);
            nim = tc.degenerate(nim, r > 15 & r < 18);
            seed = false(size(nim.FA)); seed(round(tc.C + 10), round(tc.C), 5) = true;
            o = tc.holdOpts('hold'); o.seed_mask = seed; o.max_arc = 26;
            evalc('T3 = nim_tractography_mmf(nim, o);');
            o.mmf_hold_k = 1;
            evalc('T1 = nim_tractography_mmf(nim, o);');
            tc.verifyGreaterThan(max(abs(T1{1}(end, :) - T3{1}(end, :))), 1e-3);
            o.mmf_hold_k = 3; o.mmf_hold_min = 0.7;
            evalc('[Tz, meta] = nim_tractography_mmf(nim, o);');
            tc.verifyEqual(tc.countOf(meta, 'coast'), 2);
            tc.verifyLessThan(sum(sqrt(sum(diff(Tz{1}) .^ 2, 2))), 26);
            tc.verifyEqual(meta.mmf.hold_min, 0.7);
        end
        function algorithmMmfAloneIsMMF4s(tc)
            % The schema defaults of tractography.mmf are those of MMF-4s (hold: march,
            % reanchor_tol 0.2 rad): a config that only selects the algorithm runs it, and a
            % tracker handed neither option does too. mmf_dti.yml sets hold: none and
            % reanchor_tol 0.001 and stays MMF-4a.
            here = fileparts(mfilename('fullpath'));
            cfgdir = fullfile(fileparts(fileparts(here)), 'config');
            p = [tempname '.yml'];
            fid = fopen(p, 'w'); fprintf(fid, 'tractography:\n  algorithm: mmf\n'); fclose(fid);
            c = onCleanup(@() delete(p));
            evalc('d = nim_config_to_options(load_config_yaml(p));');
            tc.verifyEqual(d.mmf_hold, 'march');
            tc.verifyEqual(d.mmf_reanchor_tol, 0.2);
            evalc('m = nim_config_to_options(load_config_yaml(fullfile(cfgdir, ''mmf_dti_march.yml'')));');
            tc.verifyEqual(d.mmf_hold, m.mmf_hold);
            tc.verifyEqual(d.mmf_reanchor_tol, m.mmf_reanchor_tol);
            evalc('a = nim_config_to_options(load_config_yaml(fullfile(cfgdir, ''mmf_dti.yml'')));');
            tc.verifyEqual(a.mmf_hold, 'none');
            tc.verifyEqual(a.mmf_reanchor_tol, 1e-3);
            o = rmfield(tc.holdOpts('none'), {'mmf_hold', 'mmf_reanchor_tol'});   % the tracker's own defaults
            evalc('[~, meta] = nim_tractography_mmf(tc.Nim, o);');
            tc.verifyEqual(meta.mmf.hold, 'march');
            tc.verifyEqual(meta.mmf.reanchor_tol, 0.2);
        end
        function thePresetsAreMmfDtiWithTheSchemeOn(tc)
            here = fileparts(mfilename('fullpath'));
            cfgdir = fullfile(fileparts(fileparts(here)), 'config');
            evalc('base = nim_config_to_options(load_config_yaml(fullfile(cfgdir, ''mmf_dti.yml'')));');
            for preset = {'hold', 'march'}
                evalc('p = nim_config_to_options(load_config_yaml(fullfile(cfgdir, [''mmf_dti_'' preset{1} ''.yml''])));');
                tc.verifyEqual(p.mmf_hold, preset{1});
                tc.verifyEqual(p.mmf_reanchor_tol, 0.2);
                tc.verifyEqual(p.angle_thresh, 0);                       % the exit test takes the place of angle_max
                p.mmf_hold = 'none'; p.mmf_reanchor_tol = base.mmf_reanchor_tol; p.angle_thresh = base.angle_thresh;
                tc.verifyEqual(p, base);
            end
        end
    end
end
