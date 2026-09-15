classdef TestRegionOverlap < matlab.unittest.TestCase
% Regions that share voxels, and the consumers that used to lose them.
%
% A parcellation LABEL VOLUME names one owner per voxel. Anatomy does not: on the
% ISMRM bundle masks 84.8% of labelled voxels belong to two or more bundles, one
% belongs to ten, and the label volume retains a median 43% of each region -
% CC_u_shaped keeps 1275 of its 106502 voxels. Any code that reads a region as
% `parcellation_mask == id`, lists regions with `unique(parcellation_mask)`, or
% takes one label per streamline point is wrong by that much, silently and always
% in the same direction.

    methods (Static)
        function nim = overlapping()
        % big (64 vox) fully covered by two smaller halves, so smallest-wins
        % leaves big with NOTHING in the label volume.
            d = [8 8 8];
            big = false(d); big(2:5, 2:5, 2:5) = true;          % 64
            a   = false(d); a(2:5, 2:5, 2:3) = true;            % 32, lower half
            b   = false(d); b(2:5, 2:5, 4:5) = true;            % 32, upper half
            R = containers.Map({'big','a','b'}, {big, a, b});
            L = zeros(d, 'uint16');
            L(big) = 1; L(a) = 2; L(b) = 3;                     % smaller painted last
            map = containers.Map({1,2,3}, {'big','a','b'});
            ret = containers.Map({'big','a','b'}, {0, 1, 1});
            nim = struct('FA', zeros(d), 'parcellation_mask', L, ...
                'atlas_labels', struct('map', map), 'roi_masks', R, ...
                'roi_overlap', struct('retained', ret, ...
                    'region_sizes', containers.Map({'big','a','b'},{64,32,32})));
        end
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
        end
    end

    methods (Test)

        function regionMaskReturnsTheWholeRegion(tc)
            nim = TestRegionOverlap.overlapping();
            tc.verifyEqual(nnz(nim.parcellation_mask == 1), 0, ...
                'fixture is wrong: big should be erased from the label volume');
            m = nim_region_mask(nim, 1);
            tc.verifyEqual(nnz(m), 64, ...
                'nim_region_mask returned the label-volume remnant, not the region.');
            tc.verifyEqual(nnz(nim_region_mask(nim, 'big')), 64, ...
                'lookup by name must agree with lookup by index.');
        end

        function regionIdsListsARegionTheLabelVolumeErased(tc)
            nim = TestRegionOverlap.overlapping();
            u = unique(nim.parcellation_mask(:)); u = u(u > 0);
            tc.verifyFalse(ismember(1, u), 'fixture is wrong: big should be absent');
            ids = nim_region_ids(nim);
            tc.verifyEqual(sort(ids(:))', [1 2 3], ...
                'a region with no surviving voxels must still be listed.');
        end

        function aSharedVoxelBelongsToEveryRegionThatContainsIt(tc)
            nim = TestRegionOverlap.overlapping();
            R = nim_region_lookup(nim);
            ids = R.at([3 3 3]);                 % inside big AND a
            tc.verifyEqual(sort(ids(:))', [1 2], ...
                'a voxel shared by two regions must report both.');
            tc.verifyEqual(R.n_multi, 64, ...
                'every voxel of big is shared, so all 64 are multi-region.');
        end

        function trackMembershipCountsSharedPoints(tc)
            nim = TestRegionOverlap.overlapping();
            R = nim_region_lookup(nim);
            track = [repmat([3 3 2], 4, 1); repmat([3 3 5], 4, 1)];  % in a, then in b
            touched = R.touched(track);
            tc.verifyEqual(sort(touched(:))', [1 2 3], ...
                'the track runs through big, a and b, so all three must be reported.');
            in_reg = nim_track_membership(track, nim_region_mask(nim, 1), R.any);
            tc.verifyEqual(sum(in_reg), 8, ...
                'every point is inside big, even though big owns no label voxels.');
        end

        function connectivityCreditsBothEndsOfAnOverlap(tc)
            % A streamline ending in voxels that each belong to two regions is
            % evidence for every pair it could join. Crediting only the label
            % winner attributes the edge to one arbitrary pair.
            nim = TestRegionOverlap.overlapping();
            track = [repmat([3 3 2], 6, 1); repmat([3 3 5], 6, 1)];
            evalc('M = nim_plot_connectivity_matrix({track}, nim, ''plot'', false);');
            tc.verifyEqual(size(M), [3 3]);
            ids = nim_region_ids(nim);
            ia = find(ids == 2); ib = find(ids == 3);
            tc.verifyGreaterThan(M(ia, ib), 0, ...
                'the a-b connection implied by the endpoints was not counted.');
            ibig = find(ids == 1);
            tc.verifyGreaterThan(M(ibig, ia) + M(ia, ibig), 0, ...
                ['big contains both endpoints and shares them with a and b, so the ' ...
                 'pairs involving big must be credited too.']);
        end
    end
end
