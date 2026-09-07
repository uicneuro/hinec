classdef TestGateLengthCriteria < matlab.unittest.TestCase
% The scorer criteria that nim_filter_tracks_roi used to ignore: any_mask and
% the length family. Two things are easy to get backwards and are pinned here.
%
% UNITS. scilpy calls sft.to_rasmm() before measuring, so every length is in
% MILLIMETRES. On this data a voxel is 2 mm, so a gate written in voxel units is
% wrong by exactly a factor of two - which passes casual inspection because the
% numbers still look plausible.
%
% NET vs TOTAL. length_x is |sum(dx)| and length_x_abs is sum(|dx|). For a track
% that doubles back they differ maximally: the first is zero, the second is the
% whole journey. Reading the scilpy docstring quickly gets these swapped.

    properties
        Nim
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            d = [40 40 10];
            blob = false(d); blob(13:17, 18:22, 4:6) = true;
            masks = containers.Map({'blob'}, {blob});
            tc.Nim = struct('FA', zeros(d), 'parcellation_mask', zeros(d), ...
                            'atlas_labels', struct('map', containers.Map(1,'unused')), ...
                            'roi_masks', masks, ...
                            'hdr', struct('PixelDimensions', [2 2 2 1]));
        end
    end

    methods
        function [kept, stats] = filt(tc, tracks, opts)
            evalc('[kept, stats] = nim_filter_tracks_roi(tracks, tc.Nim, opts);');
        end
        function p = seg(~, a, b, n)
            p = [linspace(a(1),b(1),n)', linspace(a(2),b(2),n)', linspace(a(3),b(3),n)'];
        end
    end

    methods (Test)

        function lengthIsMillimetresNotVoxels(tc)
            % 10 voxel steps at 2 mm = 20 mm. A gate of [0 15] must REJECT it;
            % it would wrongly accept if the code measured in voxels (10 < 15).
            s = tc.seg([5 20 5], [15 20 5], 41);
            tc.verifyEmpty(tc.filt({s}, struct('length', [0 15])), ...
                'length must be in mm: a 10-voxel track is 20 mm, not 10.');
            tc.verifyNumElements(tc.filt({s}, struct('length', [0 25])), 1);
        end

        function netCancelsOnADoubleBackButTotalDoesNot(tc)
            out  = tc.seg([5 20 5], [15 20 5], 41);
            back = tc.seg([15 20 5], [5 20 5], 41);
            loop = [out; back(2:end, :)];          % net x = 0, total x = 40 mm

            % net: the loop reads 0 mm and passes a tight gate the straight one fails
            tc.verifyNumElements(tc.filt({loop}, struct('length_x', [0 5])), 1);
            tc.verifyEmpty(tc.filt({out}, struct('length_x', [0 5])), ...
                'straight 20 mm track must fail a net-x gate of [0 5].');

            % total travel: exactly the other way round
            tc.verifyNumElements(tc.filt({out},  struct('length_x_abs', [0 30])), 1);
            tc.verifyEmpty(tc.filt({loop}, struct('length_x_abs', [0 30])), ...
                'doubling back is 40 mm of travel and must fail [0 30].');
        end

        function anyInNeedsOnlyASinglePointOfContact(tc)
            through = tc.seg([5 20 5], [30 20 5], 80);   % crosses the blob
            beside  = tc.seg([5 30 5], [30 30 5], 80);   % parallel, never enters
            kept = tc.filt({through, beside}, struct('any_in', {{'blob'}}));
            tc.verifyNumElements(kept, 1);
            tc.verifyEqual(kept{1}(1,2), 20, 'the wrong track survived any_in.');
        end

        function lengthWithoutVoxelSizeIsRefusedNotGuessed(tc)
            n = tc.Nim; n = rmfield(n, 'hdr');
            s = tc.seg([5 20 5], [15 20 5], 41);
            tc.verifyError(@() nim_filter_tracks_roi({s}, n, struct('length', [0 25])), ...
                'nim_filter_tracks_roi:noVoxelSize');
        end

        function noCriteriaStillReturnsEverythingUntouched(tc)
            s = tc.seg([5 20 5], [15 20 5], 41);
            [kept, stats] = tc.filt({s}, struct());
            tc.verifyNumElements(kept, 1);
            tc.verifyFalse(stats.applied);
        end
    end
end
