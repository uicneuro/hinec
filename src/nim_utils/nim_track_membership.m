function [in_region, in_any] = nim_track_membership(track, region_mask, any_mask)
% nim_track_membership: Per-point membership of a track in a region.
%
%   [in_region, in_any] = nim_track_membership(track, region_mask, any_mask)
%
% Returns logical vectors, one entry per track point: whether the point lies in
% this region, and whether it lies in any region at all.
%
% This replaces reading a single label per point out of the parcellation volume.
% A label volume names one owner per voxel, so a point in a voxel shared between
% two regions was attributed to whichever won the tie - and on the ISMRM bundle
% masks 84.8% of labelled voxels are shared. Testing the region's own mask counts
% a track for every region it genuinely enters.
%
% any_mask is optional; pass [] if the "inside any region" test is not needed.

    dims = size(region_mask);
    n = size(track, 1);
    in_region = false(n, 1);
    in_any    = false(n, 1);
    v = round(track);
    ok = all(v >= 1, 2) & v(:,1) <= dims(1) & v(:,2) <= dims(2) & v(:,3) <= dims(3);
    if ~any(ok), return; end
    li = sub2ind(dims, v(ok,1), v(ok,2), v(ok,3));
    in_region(ok) = region_mask(li);
    if nargin >= 3 && ~isempty(any_mask)
        in_any(ok) = any_mask(li);
    else
        in_any = in_region;
    end
end
