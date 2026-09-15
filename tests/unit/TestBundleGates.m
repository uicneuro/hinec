classdef TestBundleGates < matlab.unittest.TestCase
% endpoints_in and contained_in - the two predicates that make up the ISMRM 2015
% bundle definition (head + tail + all_mask).
%
% The distinction being pinned here is the one that makes them separate keys:
% include_roi asks whether a track PASSES THROUGH a region, endpoints_in asks
% where it STOPS. A streamline that runs straight through both endpoint regions
% and carries on is a waypoint match and an endpoint failure, and conflating the
% two is how a bundle fills up with streamlines that merely transit it.

    properties
        Nim
    end

    methods (TestMethodSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            addpath(genpath(fileparts(fileparts(here))));
            d = [40 40 10];
            head = false(d); head(5:8,   18:22, 4:6) = true;
            tail = false(d); tail(32:35, 18:22, 4:6) = true;
            corr = false(d); corr(3:37,  16:24, 3:7) = true;   % containment corridor
            masks = containers.Map({'head','tail','corridor'}, {head, tail, corr});
            tc.Nim = struct('FA', zeros(d), 'parcellation_mask', zeros(d), ...
                            'atlas_labels', struct('map', containers.Map(1,'unused')), ...
                            'roi_masks', masks);
        end

        function r = noop(~), r = []; end
    end

    methods
        function [kept, stats] = filt(tc, tracks, opts)
            evalc('[kept, stats] = nim_filter_tracks_roi(tracks, tc.Nim, opts);');
        end
        function p = line3(~, a, b, n)
            p = [linspace(a(1),b(1),n)', linspace(a(2),b(2),n)', linspace(a(3),b(3),n)'];
        end
    end

    methods (Test)

        function endpointPairIsKeptEitherWayRound(tc)
            fwd = tc.line3([6 20 5], [34 20 5], 60);
            rev = flipud(fwd);
            [kept, st] = tc.filt({fwd, rev}, struct('endpoints_in', {{'head','tail'}}));
            tc.verifyEqual(numel(kept), 2, 'head->tail and tail->head must both count.');
            tc.verifyEqual(st.n_dropped_endpoints, 0);
        end

        function trackEndingOutsideBothIsDropped(tc)
            % Passes through head AND tail, then keeps going and stops elsewhere.
            % A waypoint test would accept this; an endpoint test must not.
            through = tc.line3([6 20 5], [39 20 5], 70);
            [kept, st] = tc.filt({through}, struct('endpoints_in', {{'head','tail'}}));
            tc.verifyEmpty(kept, 'A track that only transits the endpoint regions was kept.');
            tc.verifyEqual(st.n_dropped_endpoints, 1);

            % ... and the waypoint test does accept it, which is the difference.
            kept2 = tc.filt({through}, struct('include_roi', {{'head','tail'}}, ...
                                              'roi_filter_mode', 'all'));
            tc.verifyEqual(numel(kept2), 1, ...
                'include_roi is a waypoint test and should accept a transiting track.');
        end

        function bothEndsInTheSameRegionIsDropped(tc)
            loop = tc.line3([6 19 5], [7 21 5], 20);
            kept = tc.filt({loop}, struct('endpoints_in', {{'head','tail'}}));
            tc.verifyEmpty(kept, 'Both endpoints in head must not satisfy a head/tail pair.');
        end

        function aSingleExcursionIsNotStrictlyContained(tc)
            % The SCORER's rule is all-or-nothing, and that rule is still what
            % n_strict_contained reports - it is the scorer-comparable number.
            % What changed is that the streamline is no longer thrown away for
            % it: the excursion is cut and the valid part kept. Both facts are
            % asserted here so neither can drift.
            inside  = tc.line3([6 20 5], [34 20 5], 60);
            strayer = inside;
            strayer(30, 2) = 30;                 % one point outside the corridor
            [kept, st] = tc.filt({inside, strayer}, struct('contained_in', {{'corridor'}}));
            tc.verifyEqual(st.n_strict_contained, 1, ...
                'one excursion must disqualify a track from being STRICTLY contained.');
            tc.verifyEqual(numel(kept), 2, ...
                'but the straying track must survive in truncated form, not be discarded.');
            tc.verifyEqual(st.n_truncated_contained, 1);
            tc.verifyLessThan(size(kept{2},1), size(inside,1), ...
                'the surviving remnant must be shorter than the original.');
        end

        function gatesCombine(tc)
            good    = tc.line3([6 20 5], [34 20 5], 60);
            wrongEnd= tc.line3([6 20 5], [20 20 5], 40);   % inside corridor, wrong endpoint
            outside = good; outside(30,2) = 30;            % right endpoints, leaves corridor
            [kept, st] = tc.filt({good, wrongEnd, outside}, ...
                struct('endpoints_in', {{'head','tail'}}, 'contained_in', {{'corridor'}}));
            % 'outside' is now truncated rather than dropped, and its truncated
            % form no longer reaches the tail, so the endpoint gate removes it.
            tc.verifyEqual(numel(kept), 1, 'Only the track satisfying BOTH gates should survive.');
            tc.verifyEqual(st.n_strict_contained, 2, ...
                'good and wrongEnd lie wholly inside; only outside strays.');
            tc.verifyEqual(st.n_dropped_endpoints, 2);
        end

        function containmentTruncatesInsteadOfDiscarding(tc)
            % The scorer's rule is all-or-nothing: every point inside, or the
            % streamline is not that bundle. That is right for SCORING and wrong
            % for keeping data - a track that follows the corridor for most of its
            % length and strays briefly used to lose all of it. On Cingulum_right
            % containment alone rejects 80% of produced streamlines, and it
            % rejects hardest the further a tracker explores.
            %
            % So the excursion is cut and the longest inside-run kept, while the
            % strict count stays available for scorer-comparable reporting.
            inside_a = tc.line3([6 20 5], [20 20 5], 40);    % inside the corridor
            outside  = tc.line3([20 20 5], [20 30 5], 20);   % strays out (y>24)
            inside_b = tc.line3([20 30 5], [30 30 5], 20);
            strayer  = [inside_a; outside(2:end,:); inside_b(2:end,:)];
            clean    = tc.line3([6 20 5], [34 20 5], 60);

            [kept, st] = tc.filt({clean, strayer}, struct('contained_in', {{'corridor'}}));
            tc.verifyNumElements(kept, 2, ...
                'the straying track must be kept in truncated form, not discarded.');
            tc.verifyEqual(st.n_strict_contained, 1, ...
                'exactly one track lies wholly inside and must be reported as strict.');
            tc.verifyEqual(st.n_truncated_contained, 1, ...
                'the straying track must be reported as truncated.');

            % and what survives really is the inside portion
            lens = cellfun(@(p) size(p,1), kept);
            [~, shorter] = min(lens);
            v = round(kept{shorter});
            m = tc.Nim.roi_masks('corridor');
            idx = sub2ind(size(tc.Nim.FA), v(:,1), v(:,2), v(:,3));
            tc.verifyTrue(all(m(idx)), ...
                'every point of the truncated track must lie inside the corridor.');
        end
    end
end
