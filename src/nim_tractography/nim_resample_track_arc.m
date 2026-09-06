function out = nim_resample_track_arc(tracks, arc_step)
% nim_resample_track_arc: Resample streamlines to fixed arc-length spacing.
%
%   tracks = nim_resample_track_arc(tracks, 0.5)   % a point every 0.5 voxels
%
% Decouples the STORED representation from the integration step. Every
% integration step is otherwise kept, so storage grows as 1/step: a full-brain
% run at step 0.01 would write ~60 GB. Integration still runs at full
% resolution; only the saved polyline is coarsened.
%
% Endpoints are always preserved exactly - they carry the anatomical
% connectivity that scoring is built on, so they must not be interpolated away.
% A track shorter than one arc_step degenerates to its two endpoints.
%
% arc_step <= 0 (or empty) returns the input untouched.
%
% Accepts a single N-by-3 track or a cell array of them.

    if nargin < 2 || isempty(arc_step) || arc_step <= 0
        out = tracks;
        return;
    end

    if ~iscell(tracks)
        out = resample_one(tracks, arc_step);
        return;
    end

    out = cell(size(tracks));
    for i = 1:numel(tracks)
        out{i} = resample_one(tracks{i}, arc_step);
    end
end

% =========================================================================
function t = resample_one(t, arc_step)
    if isempty(t) || size(t, 1) < 2
        return;
    end

    seg = sqrt(sum(diff(t, 1, 1).^2, 2));
    s = [0; cumsum(seg)];
    L = s(end);

    if L <= 0
        t = t([1 end], :);
        return;
    end
    if L <= arc_step
        t = [t(1, :); t(end, :)];
        return;
    end

    % Strictly increasing knots for interp1: drop zero-length segments.
    keep = [true; seg > 0];
    s = s(keep);
    p = t(keep, :);
    if numel(s) < 2
        t = t([1 end], :);
        return;
    end

    q = (0:arc_step:L)';
    if q(end) < L - 1e-12
        q(end+1) = L;      % always land exactly on the far endpoint
    else
        q(end) = L;
    end

    t = interp1(s, p, q, 'linear');

    % Guard against floating-point drift at the ends.
    t(1, :)   = p(1, :);
    t(end, :) = p(end, :);
end
