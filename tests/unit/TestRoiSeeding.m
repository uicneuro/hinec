classdef TestRoiSeeding < matlab.unittest.TestCase
    % ROI seeding, include/exclude track filtering, and the ROI <-> whole-brain
    % identity that ROI scoring rests on (plans/NEXT_STEPS.md A0): the tracks of
    % an ROI-seeded run are EXACTLY the tracks of a whole-brain run whose seeds
    % lie in the ROI - same count, bit-identical polylines - for every tracker.
    % Uses small self-contained nims so the tests stay fast and do not depend on
    % the 260 MB ISMRM nim being present.

    properties
        Root
        Nim
        TrackNim   % phantom carrying a direction field, for the tracking tests
        SeedW      % whole-brain seed mask on that phantom
        SeedR      % a proper subset of it (the "ROI")
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fullfile(here, '..', '..');
            addpath(fullfile(tc.Root, 'src', 'nim_utils'));
            addpath(fullfile(tc.Root, 'src', 'nim_tractography'));
            addpath(fullfile(tc.Root, 'src', 'nim_calculation'));

            % 20x20x20 volume with three labelled blocks and a label map.
            d = [20 20 20];
            P = zeros(d);
            P(3:7,   3:7,   3:7)   = 41;   % "Superior longitudinal fasciculus R"
            P(12:16, 3:7,   3:7)   = 42;   % "Superior longitudinal fasciculus L"
            P(3:5,   12:14, 12:14) = 7;    % "Corticospinal tract R"

            m = containers.Map('KeyType','double','ValueType','char');
            m(41) = 'Superior longitudinal fasciculus R';
            m(42) = 'Superior longitudinal fasciculus L';
            m(7)  = 'Corticospinal tract R';
            m(3)  = 'Genu of corpus callosum';

            n = struct();
            n.parcellation_mask = P;
            n.FA   = 0.5 * ones(d);
            n.mask = ones(d);
            n.atlas_labels = struct('map', m, 'atlas_type', 'jhu');
            tc.Nim = n;

            % Separate phantom for the ROI <-> whole-brain identity (A0). It
            % needs a direction field, which the parcellation phantom above has
            % no use for.
            [tc.TrackNim, tc.SeedW, tc.SeedR] = TestRoiSeeding.trackingPhantom();
        end
    end

    methods (Test)

        % ------------------------------------------------------- resolution
        function resolvesByIndex(tc)
            [mask, info] = nim_roi_mask(tc.Nim, {41});
            tc.verifyEqual(info.labels, 41);
            tc.verifyEqual(sum(mask(:)), 125);        % 5^3
        end

        function resolvesByFullName(tc)
            [~, info] = nim_roi_mask(tc.Nim, {'Superior longitudinal fasciculus R'});
            tc.verifyEqual(info.labels, 41);
        end

        function resolvesByShortAlias(tc)
            [~, i1] = nim_roi_mask(tc.Nim, {'SLF_R'});
            tc.verifyEqual(i1.labels, 41);
            [~, i2] = nim_roi_mask(tc.Nim, {'SLF_L'});
            tc.verifyEqual(i2.labels, 42);
            [~, i3] = nim_roi_mask(tc.Nim, {'CST_R'});
            tc.verifyEqual(i3.labels, 7);
        end

        function resolvesMixedIndexAndName(tc)
            % The requirement was that indices and names mix freely in one list.
            [mask, info] = nim_roi_mask(tc.Nim, {41, 'SLF_L', 'Corticospinal tract R'});
            tc.verifyEqual(sort(info.labels), [7 41 42]);
            tc.verifyEqual(sum(mask(:)), 125 + 125 + 27);
        end

        function dilationGrowsTheMask(tc)
            [~, i0] = nim_roi_mask(tc.Nim, {41}, 0);
            [~, i1] = nim_roi_mask(tc.Nim, {41}, 1);
            tc.verifyGreaterThan(i1.dilated_voxels, i0.dilated_voxels);
            tc.verifyEqual(i0.dilated_voxels, i0.raw_voxels);
        end

        % ------------------------------------------------------------ errors
        function unknownNameErrors(tc)
            tc.verifyError(@() nim_roi_mask(tc.Nim, {'NotARegion'}), ...
                'nim_roi_mask:unknownName');
        end

        function ambiguousNameErrorsWithCandidates(tc)
            n = tc.Nim;
            m = n.atlas_labels.map;
            m(35) = 'Cingulum (cingulate gyrus) R';
            m(37) = 'Cingulum (hippocampus) R';
            n.atlas_labels.map = m;
            tc.verifyError(@() nim_roi_mask(n, {'Cingulum'}), 'nim_roi_mask:ambiguousName');
        end

        function missingParcellationErrors(tc)
            n = rmfield(tc.Nim, 'parcellation_mask');
            tc.verifyError(@() nim_roi_mask(n, {41}), 'nim_roi_mask:noParcellation');
        end

        % --------------------------------------------------------- filtering
        function filterIsNoOpWhenUnset(tc)
            tracks = {[3 3 3; 4 4 4; 5 5 5], [12 3 3; 13 4 4]};
            [out, st] = nim_filter_tracks_roi(tracks, tc.Nim, struct());
            tc.verifyFalse(st.applied);
            tc.verifyEqual(numel(out), 2);
        end

        function includeRoiKeepsOnlyTracksTouchingIt(tc)
            inside  = [3 3 3; 4 4 4; 5 5 5];      % inside label 41
            outside = [12 3 3; 13 4 4; 14 5 5];   % inside label 42, not 41
            opts = struct('include_roi', {{41}}, 'roi_filter_mode', 'all');
            [out, st] = nim_filter_tracks_roi({inside, outside}, tc.Nim, opts);
            tc.verifyTrue(st.applied);
            tc.verifyEqual(numel(out), 1);
            tc.verifyEqual(out{1}, inside);
            tc.verifyEqual(st.n_dropped_include, 1);
        end

        function excludeRoiDropsTracksTouchingIt(tc)
            inside  = [3 3 3; 4 4 4];
            outside = [12 3 3; 13 4 4];
            opts = struct('exclude_roi', {{41}});
            [out, st] = nim_filter_tracks_roi({inside, outside}, tc.Nim, opts);
            tc.verifyEqual(numel(out), 1);
            tc.verifyEqual(out{1}, outside);
            tc.verifyEqual(st.n_dropped_exclude, 1);
        end

        function includeModeAllVersusAny(tc)
            % A track through 41 only. mode 'all' over {41,42} must reject it;
            % mode 'any' must keep it.
            tr = {[3 3 3; 4 4 4; 5 5 5]};
            all_opts = struct('include_roi', {{41, 42}}, 'roi_filter_mode', 'all');
            any_opts = struct('include_roi', {{41, 42}}, 'roi_filter_mode', 'any');
            tc.verifyEqual(numel(nim_filter_tracks_roi(tr, tc.Nim, all_opts)), 0);
            tc.verifyEqual(numel(nim_filter_tracks_roi(tr, tc.Nim, any_opts)), 1);
        end

        function tracksOutsideVolumeAreDropped(tc)
            tr = {[999 999 999; 1000 1000 1000]};
            opts = struct('include_roi', {{41}});
            tc.verifyEqual(numel(nim_filter_tracks_roi(tr, tc.Nim, opts)), 0);
        end

        % ------------------------------------------------------- seed count
        function seedDensityIsHonouredExactly(tc)
            % Regression: the old inline lattice used per_axis = ceil(d^(1/3))
            % and returned per_axis^3 offsets, so 2..7 all gave 8 and 9..26 gave
            % 27. seed_density: 4 placed 8 seeds per voxel.
            for d = [1 2 3 4 5 6 7 8 9 12 16 27 64]
                off = nim_seed_offsets(d);
                tc.verifyEqual(size(off, 1), d, ...
                    sprintf('seed_density %d produced %d offsets', d, size(off,1)));
                tc.verifyTrue(all(off(:) >= -0.5 & off(:) <= 0.5), ...
                    'seed offsets must lie inside the voxel');
            end
        end

        function cubeDensitiesKeepTheOriginalLattice(tc)
            % Perfect cubes must be byte-identical to the previous behaviour, so
            % every shipped config (which uses 1 or 8) is unchanged.
            for d = [1 8 27 64]
                per = max(1, ceil(d^(1/3) - 1e-12));
                e = linspace(-0.5, 0.5, per + 1);
                c = (e(1:end-1) + e(2:end)) / 2;
                [ox, oy, oz] = ndgrid(c, c, c);
                tc.verifyEqual(nim_seed_offsets(d), [ox(:), oy(:), oz(:)], ...
                    sprintf('cube density %d changed', d));
            end
        end

        function seedOffsetsAreDeterministic(tc)
            % The convergence ladder compares streamline i across runs, which
            % requires identical seed placement every time. No RNG anywhere.
            for d = [4 5 12]
                tc.verifyEqual(nim_seed_offsets(d), nim_seed_offsets(d));
            end
        end

        function nonCubeSeedsAreSpreadNotClustered(tc)
            % Taking the first N lattice points would clump them against one
            % face; the farthest-point subset keeps them spread.
            off = nim_seed_offsets(4);
            D = pdist2(off, off); D(1:size(D,1)+1:end) = inf;
            tc.verifyGreaterThanOrEqual(min(D(:)), 0.4, ...
                '4 seeds in a voxel should be well separated');
        end

        % ------------------------------------ ROI <-> whole-brain identity
        function hinecRoiTracksAreExactlyTheWholeBrainSubset(tc)
            % THE result ROI seeding rests on (NEXT_STEPS A0): seeds are a
            % deterministic lattice (nim_seed_offsets; density 1 = the voxel
            % centre) and the trackers keep no cross-seed state, so seeding a
            % subset R of the whole-brain mask W must reproduce EXACTLY the
            % whole-brain streamlines whose seed lies in R - same count,
            % bit-identical polylines. That is what makes an ROI run's
            % bundle_wise row equal the whole-brain run's (verified on real
            % data: UF_left, VS 1492 / TP 9195 / FP 8098 / FN 3221 both ways),
            % and it is the gate for every Phase B edit. hinec is the spine, so
            % it gets the direct check: pair the two runs by SEED, using the
            % per-track seed the tracker reports.
            o = TestRoiSeeding.trackerOptions(tc.Root, 'hinec', {});
            [TW, MW] = TestRoiSeeding.runTracker(tc.TrackNim, o, tc.SeedW);
            TR       = TestRoiSeeding.runTracker(tc.TrackNim, o, tc.SeedR);
            tc.assertNotEmpty(TR, 'the ROI run produced no tracks - nothing was tested');

            inR = TestRoiSeeding.seedsIn(MW.seed_points, tc.SeedR);
            tc.verifyEqual(numel(TR), nnz(inR), sprintf( ...
                ['ROI-seeded run produced %d tracks, but %d of the %d whole-brain ' ...
                 'tracks were seeded inside the ROI.'], numel(TR), nnz(inR), numel(TW)));

            expected = TW(inR);
            seeds    = MW.seed_points(inR, :);
            for i = 1:min(numel(TR), numel(expected))
                if ~isequal(TR{i}, expected{i})
                    tc.verifyFail(sprintf( ...
                        ['ROI track %d (seed [%s]) differs from its whole-brain ' ...
                         'counterpart: %s'], i, num2str(seeds(i,:), '%.4f '), ...
                        TestRoiSeeding.trackDiff(TR{i}, expected{i})));
                    break
                end
            end
        end

        function theIdentityCheckCanFail(tc)
            % Negative control. The assertion above is worth having only if a
            % settings difference between the two runs breaks it - which is
            % exactly the failure it exists to catch: the "ROI vs whole-brain
            % disagree" confusion was a whole-brain run carrying an extra
            % --set that the ROI runs did not. A different integration step
            % must make the paired polylines differ.
            oW = TestRoiSeeding.trackerOptions(tc.Root, 'hinec', {});
            oR = TestRoiSeeding.trackerOptions(tc.Root, 'hinec', {'integrator.step=0.4'});
            [TW, MW] = TestRoiSeeding.runTracker(tc.TrackNim, oW, tc.SeedW);
            TR       = TestRoiSeeding.runTracker(tc.TrackNim, oR, tc.SeedR);
            inR = TestRoiSeeding.seedsIn(MW.seed_points, tc.SeedR);
            tc.verifyFalse(isequal(TR, TW(inR)), ...
                ['two runs with DIFFERENT integrator steps produced identical ' ...
                 'polylines - the identity assertion cannot fail, so it tests nothing.']);
        end

        function standardRoiTracksAreExactlyTheWholeBrainSubset(tc)
            % FACT reports no per-track seed, so the same identity is checked as
            % a partition (see verifyPartitionIdentity).
            o = TestRoiSeeding.trackerOptions(tc.Root, 'standard', {});
            tc.verifyPartitionIdentity(tc.TrackNim, o);
        end

        function mmfRoiTracksAreExactlyTheWholeBrainSubset(tc)
            o = TestRoiSeeding.trackerOptions(tc.Root, 'mmf', {});
            % The connection-form tracer CONSUMES the moving-frame geometry;
            % building it is runTractography's step 3, so build it here.
            w = warning('off','all'); cw = onCleanup(@() warning(w)); %#ok<NASGU>
            nim = tc.TrackNim;
            evalc('nim = nim_mmf_geometry(nim, o);');
            tc.verifyPartitionIdentity(nim, o);
        end

        % ------------------------------------------------------------ config
        function roiKeysAreInTheSchemaAndReachableFromCli(tc)
            S = nim_config_schema();
            for p = {'tractography.seeding.roi', 'tractography.seeding.roi_dilate', ...
                     'tractography.filter.include_roi', 'tractography.filter.exclude_roi', ...
                     'tractography.filter.mode', 'tractography.filter.roi_dilate'}
                tc.verifyTrue(any(strcmp({S.path}, p{1})), sprintf('%s missing from schema', p{1}));
            end
            w = warning('off','all'); cw = onCleanup(@() warning(w));
            cfg = load_config_yaml(fullfile(tc.Root, 'config', 'hinec_dti.yml'));
            c2 = nim_config_apply_overrides(cfg, ...
                {'seeding.roi=[41,42]', 'seeding.roi_dilate=1', 'filter.include_roi=SLF_R'});
            tc.verifyEqual(numel(c2.tractography.seeding.roi), 2);
            tc.verifyEqual(c2.tractography.seeding.roi_dilate, 1);
            o = nim_config_to_options(c2);
            tc.verifyEqual(numel(o.seed_roi), 2);
            tc.verifyEqual(o.include_roi{1}, 'SLF_R');
        end
    end

    methods (Access = private)
        function verifyPartitionIdentity(tc, nim, o)
            % Identity check for the trackers that report no per-track seed.
            % Seeds are enumerated in find(seed_mask) order, so the seeds of a
            % run on R and a run on its complement Rc = W \ R are complementary
            % ordered subsequences of the whole-brain run's. Hence the two runs
            % must reproduce the whole-brain tracks exactly, track for track and
            % in order - which is the A0 identity plus the statement that no
            % whole-brain track is left over.
            Rc = tc.SeedW & ~tc.SeedR;
            TW = TestRoiSeeding.runTracker(nim, o, tc.SeedW);
            TR = TestRoiSeeding.runTracker(nim, o, tc.SeedR);
            TC = TestRoiSeeding.runTracker(nim, o, Rc);
            tc.assertNotEmpty(TR, 'the ROI run produced no tracks - nothing was tested');
            tc.verifyEqual(numel(TW), numel(TR) + numel(TC), sprintf( ...
                ['%s: whole-brain produced %d tracks; the two halves of the same ' ...
                 'seed mask produced %d + %d - seeding does not partition.'], ...
                o.algorithm, numel(TW), numel(TR), numel(TC)));
            i = 1; j = 1;
            for k = 1:numel(TW)
                if i <= numel(TR) && isequal(TW{k}, TR{i})
                    i = i + 1;
                elseif j <= numel(TC) && isequal(TW{k}, TC{j})
                    j = j + 1;
                else
                    ref = TR{min(i, numel(TR))};
                    tc.verifyFail(sprintf( ...
                        ['%s: whole-brain track %d matches neither the next ROI track ' ...
                         '(%d of %d) nor the next complement track (%d of %d). Against ' ...
                         'the ROI one: %s'], o.algorithm, k, i, numel(TR), j, numel(TC), ...
                        TestRoiSeeding.trackDiff(TW{k}, ref)));
                    return
                end
            end
            tc.verifyEqual([i-1, j-1], [numel(TR), numel(TC)], sprintf( ...
                '%s: not every ROI / complement track appeared in the whole-brain run.', ...
                o.algorithm));
        end
    end

    methods (Static, Access = private)

        function [nim, W, R] = trackingPhantom()
            % A circular direction field inside an annulus: curved, so a changed
            % integrator setting visibly moves the polylines (the negative
            % control needs that), and small enough to track three times per
            % tracker in a couple of seconds.
            d = [20 20 20]; c = 10.5;
            [X, Y, ~] = ndgrid(1:d(1), 1:d(2), 1:d(3));
            dx = X - c; dy = Y - c; r = sqrt(dx.^2 + dy.^2); r(r < 1e-9) = 1e-9;
            evec = zeros([d 3 3]);
            evec(:,:,:,1,1) = -dy ./ r;  evec(:,:,:,2,1) = dx ./ r;   % e1 circular
            evec(:,:,:,1,2) =  dx ./ r;  evec(:,:,:,2,2) = dy ./ r;   % e2 radial
            evec(:,:,:,3,3) = 1;
            FA = 0.02 * ones(d); FA(r > 3 & r < 9) = 0.6;             % trackable annulus
            nim = struct('FA', FA, 'evec', evec, ...
                'eval', repmat(reshape([3 1 1],1,1,1,3), [d 1]), 'mask', true(d));
            W = false(d); W(:,:,10) = true; W = W & (FA > 0.4);  % 224 voxels, one slice
            R = false(d); R(:,1:10,:) = true; R = R & W;         % half of them
        end

        function o = trackerOptions(root, algorithm, extra)
            % Options come through the SHIPPED surface - load_config_yaml ->
            % nim_config_apply_overrides -> nim_config_to_options - so what is
            % pinned is the option path run_tractography.sh actually uses.
            sets = [{sprintf('algorithm=%s', algorithm), 'seeding.density=1', ...
                     'integrator.method=rk4', 'integrator.step=0.5', ...
                     'termination.angle_max=60', 'termination.max_arc=25', ...
                     'termination.min_arc=0', 'termination.fa_min=0.1', ...
                     'diagnostics=false'}, extra];
            w = warning('off','all'); cw = onCleanup(@() warning(w)); %#ok<NASGU>
            evalc(['cfg = load_config_yaml(fullfile(root, ''config'', ''hinec_dti.yml''));' ...
                   'cfg = nim_config_apply_overrides(cfg, sets);']);
            o = nim_config_to_options(cfg);
            o.wm_mask = []; o.gm_mask = []; o.csf_mask = [];
            o.enable_diagnostics = false;
        end

        function [tracks, meta] = runTracker(nim, o, seed_mask)
            o.seed_mask = seed_mask;
            meta = struct();
            switch char(o.algorithm)
                case 'hinec',    evalc('[tracks, meta] = nim_tractography_hinec(nim, o);');
                case 'mmf',      evalc('[tracks, meta] = nim_tractography_mmf_connframe(nim, o);');
                case 'standard', evalc('tracks = nim_tractography_standard(nim, o);');
                otherwise, error('TestRoiSeeding:algorithm', 'unknown algorithm %s', o.algorithm);
            end
        end

        function tf = seedsIn(pts, M)
            % seed_density 1 puts one seed at each voxel CENTRE, so rounding
            % recovers the voxel that produced it exactly.
            d = size(M); p = round(pts);
            ok = all(p >= 1, 2) & all(p <= d, 2);
            tf = false(size(pts, 1), 1);
            tf(ok) = M(sub2ind(d, p(ok,1), p(ok,2), p(ok,3)));
        end

        function s = trackDiff(a, b)
            if size(a, 1) ~= size(b, 1)
                s = sprintf('lengths differ, %d vs %d points', size(a,1), size(b,1));
                return
            end
            k = find(any(a ~= b, 2), 1);
            if isempty(k), s = 'identical'; return; end
            s = sprintf('first differing point %d of %d: [%s] vs [%s]', ...
                k, size(a,1), num2str(a(k,:), '%.6f '), num2str(b(k,:), '%.6f '));
        end
    end

end
