function bench_spline_callchain(nim_file)
% bench_spline_callchain: reproduce the call chain of nim_tractography_hinec on synthetic
% positions to time one RK4-like step (4 lookups of six dyadic interpolants + FA) by kernel:
%   client builds the interpolants into a struct -> parfor on a one-worker pool ->
%   per half-track function receives the struct by value -> 4 lookups per step,
%   each a function that receives the struct by value and queries ONE point.
% Variant 'warm': the same, but every interpolant is queried once inside the
% half-track function before its steps (the cache then lives in that workspace).
if nargin < 1 || isempty(nim_file), nim_file = 'data/ismrm2015/ismrm2015.mat'; end
setenv('HINEC_MAX_WORKERS', '1'); delete(gcp('nocreate')); parpool('local', 1);
S = load(nim_file, 'nim'); nim = S.nim; dims = size(nim.FA);
v = reshape(nim.evec(:, :, :, :, 1), [dims 3]); v(~isfinite(v)) = 0;
D = {v(:,:,:,1).^2, v(:,:,:,2).^2, v(:,:,:,3).^2, v(:,:,:,1).*v(:,:,:,2), v(:,:,:,1).*v(:,:,:,3), v(:,:,:,2).*v(:,:,:,3)};
gv = {1:dims(1), 1:dims(2), 1:dims(3)};
nsteps = 50;
pos0 = dims / 2;
for kern = {'spline', 'cubic'}
    s = struct();
    for q = 1:6, s.(sprintf('G%d', q)) = griddedInterpolant(gv, D{q}, kern{1}, 'none'); end
    s.FA = griddedInterpolant(gv, double(nim.FA), 'linear', 'none');
    for warm = [false true]
        t = zeros(1, 3);
        parfor r = 1:3
            t(r) = half_track(s, pos0 + r, nsteps, warm);
        end
        fprintf('%-6s warm=%d  per step %.3f ms  (3 half-tracks of %d steps: %s ms/step)\n', ...
            kern{1}, warm, 1e3 * mean(t) / nsteps, nsteps, mat2str(round(1e3 * t / nsteps, 2)));
    end
    % Variant 'client-warm': query every interpolant once in the CLIENT, before the struct
    % is sent to the worker; the half-track function does not warm.
    s = warm_all(s, pos0);
    t = zeros(1, 3);
    parfor r = 1:3
        t(r) = half_track(s, pos0 + r, nsteps, false);
    end
    fprintf('%-6s client-warm  per step %.3f ms  (%s ms/step)\n', kern{1}, 1e3 * mean(t) / nsteps, ...
        mat2str(round(1e3 * t / nsteps, 2)));
    % Variant 'pool-constant': a parallel.pool.Constant whose value is built AND warmed on
    % each worker once (function-handle form); every iteration reads C.Value.
    C = parallel.pool.Constant(@() warm_all(s, pos0));
    t = zeros(1, 3);
    parfor r = 1:3
        t(r) = half_track(C.Value, pos0 + r, nsteps, false);
    end
    fprintf('%-6s pool-constant per step %.3f ms  (%s ms/step)\n', kern{1}, 1e3 * mean(t) / nsteps, ...
        mat2str(round(1e3 * t / nsteps, 2)));
end
delete(gcp('nocreate'));
end

function t = half_track(s, pos, nsteps, warm)
if warm, lookup(s, pos); s = warm_all(s, pos); end
t0 = tic;
for n = 1:nsteps
    for stage = 1:4, lookup(s, pos + 0.01 * stage); end
    pos = pos + [0.25 0 0] * 0.1;
end
t = toc(t0);
end

function s = warm_all(s, pos)
f = fieldnames(s);
for q = 1:numel(f), s.(f{q})(pos(1), pos(2), pos(3)); end
end

function d = lookup(s, pos)
fa = s.FA(pos(1), pos(2), pos(3)); %#ok<NASGU>
d = zeros(1, 6);
for q = 1:6, d(q) = s.(sprintf('G%d', q))(pos(1), pos(2), pos(3)); end
end
