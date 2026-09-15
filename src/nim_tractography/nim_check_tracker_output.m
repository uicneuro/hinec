function nim_check_tracker_output(tracks, meta)
% nim_check_tracker_output: Assert a tracker's return values obey the contract.
%
%   nim_check_tracker_output(tracks, meta)
%
% Run by runTractography right after the tracker returns, so that a tracker
% that breaks the contract is caught at the boundary with a message naming the
% rule, not three functions later with an index error. The contract:
%
%   tracks  cell {T x 1}; tracks{t} is an N_t x 3 real double, N_t >= 2,
%           voxel-space coordinates, finite.
%   meta    struct. May be empty. If it carries seed_index it must be 1 x T
%           (one entry per KEPT track); if it carries seed_points, T x 3.
%
% Logs a one-paragraph summary of what came back, so the run log shows the
% output shape next to the input description from nim_describe_tracker_input.

if ~iscell(tracks)
    error('tracker:outputContract', 'tracks must be a cell array, got %s.', class(tracks));
end
tracks = tracks(:);
T = numel(tracks);

n_pts = zeros(T, 1);
for t = 1:T
    tr = tracks{t};
    if ~isnumeric(tr) || ~isreal(tr) || ndims(tr) ~= 2 || size(tr, 2) ~= 3 %#ok<ISMAT>
        error('tracker:outputContract', ...
            'tracks{%d} must be an N x 3 real numeric matrix of voxel coordinates, got %s %s.', ...
            t, class(tr), mat2str(size(tr)));
    end
    if size(tr, 1) < 2
        error('tracker:outputContract', 'tracks{%d} has %d point(s); a track needs at least 2.', t, size(tr, 1));
    end
    if ~all(isfinite(tr(:)))
        error('tracker:outputContract', 'tracks{%d} contains NaN or Inf coordinates.', t);
    end
    n_pts(t) = size(tr, 1);
end

if ~isstruct(meta)
    error('tracker:outputContract', 'meta must be a struct (empty struct() is fine), got %s.', class(meta));
end
if isfield(meta, 'seed_index') && numel(meta.seed_index) ~= T
    error('tracker:outputContract', ...
        'meta.seed_index has %d entries but %d tracks were returned; it must index the KEPT tracks.', ...
        numel(meta.seed_index), T);
end
if isfield(meta, 'seed_points') && ~isequal(size(meta.seed_points), [T 3])
    error('tracker:outputContract', ...
        'meta.seed_points is %s; it must be T x 3 with T = %d kept tracks.', mat2str(size(meta.seed_points)), T);
end

fprintf('=== Tracker output (contract satisfied) ===\n');
if T == 0
    fprintf('tracks: cell {0 x 1} - nothing survived\n');
else
    fprintf('tracks: cell {%d x 1}, each N x 3 voxel coordinates; N: min %d, median %d, max %d\n', ...
        T, min(n_pts), round(median(n_pts)), max(n_pts));
    fprintf('  tracks{1}(1,:)   = %s   (first point of the first track)\n', mat2str(tracks{1}(1, :), 5));
    fprintf('  tracks{1}(end,:) = %s   (last point)\n', mat2str(tracks{1}(end, :), 5));
end
mf = fieldnames(meta);
if isempty(mf)
    fprintf('meta:   empty struct\n');
else
    fprintf('meta:   %s\n', strjoin(mf', ', '));
    if isfield(meta, 'n_seeds')
        fprintf('  %d of %d seeds produced a track that passed min_arc\n', T, meta.n_seeds);
    end
end
fprintf('===========================================\n');
end
