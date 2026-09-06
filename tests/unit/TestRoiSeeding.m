classdef TestRoiSeeding < matlab.unittest.TestCase
    % Phase 1: ROI seeding and include/exclude track filtering.
    % Uses a small self-contained nim so the tests stay fast and do not depend
    % on the 260 MB ISMRM nim being present.

    properties
        Root
        Nim
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            tc.Root = fullfile(here, '..', '..');
            addpath(fullfile(tc.Root, 'src', 'nim_utils'));
            addpath(fullfile(tc.Root, 'src', 'nim_tractography'));

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
end
