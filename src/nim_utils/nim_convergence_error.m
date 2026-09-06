function E = nim_convergence_error(test_file, ref_file, opts)
% nim_convergence_error: Per-streamline geometric error against a reference run.
%
%   E = nim_convergence_error('run_h0.125/.../tracks.mat', 'run_ref/.../tracks.mat')
%   E = nim_convergence_error(test, ref, struct('prefix_arc', 20))
%
% Streamlines are matched by SEED, not by position in the cell array. The tracker
% compacts out seeds that produced no track, so track index != seed index: if a
% seed succeeds at one step size and fails at another, comparing "streamline i"
% would silently compare different fibres. track_meta.seed_index (written by
% runTractography) makes the correspondence explicit.
%
% Two error families are reported, and the distinction matters:
%
%   PREFIX  - compared over +/- opts.prefix_arc voxels of arc EITHER SIDE OF THE
%             SEED (default 20), on tracks long enough in both directions in BOTH
%             runs. Anchoring at the seed matters: a track is stored as
%             [reversed backward half; forward half], so its first point is a
%             TERMINATION point, and measuring a prefix from there reintroduces
%             exactly the termination difference this metric exists to exclude.
%             Use this to fit the observed order.
%
%   WHOLE   - compared over the full common arc length, plus endpoint distance
%             and arc-length difference. This is dominated by TERMINATION, which
%             is discontinuous in the step size: a track that stops one step
%             earlier moves its endpoint macroscopically. Use it to characterise
%             termination sensitivity, NOT to fit an order.
%
% Statistics are reported as median and p95, never mean: crossing regions make
% the error distribution heavy-tailed and a mean is dominated by a handful of
% diverged streamlines.

    if nargin < 3, opts = struct(); end
    prefix_arc = getf(opts, 'prefix_arc', 20);      % voxels
    n_samp     = getf(opts, 'n_samples', 64);       % resample points per curve

    [T, Tm, h_test] = load_tracks(test_file);
    [R, Rm, h_ref]  = load_tracks(ref_file);

    if isempty(Tm) || isempty(Rm)
        error('nim_convergence_error:noMeta', ...
            ['One of these runs has no track_meta.seed_index, so streamlines cannot be ' ...
             'matched by seed. Re-run with the current runTractography (older runs did ' ...
             'not record it).']);
    end

    % ---- match by seed ---------------------------------------------------
    [common, it, ir] = intersect(Tm.seed_index, Rm.seed_index, 'stable');
    E = struct();
    E.n_test = numel(T); E.n_ref = numel(R); E.n_matched = numel(common);
    E.matched_frac = numel(common) / max(numel(R), 1);
    if isempty(common)
        error('nim_convergence_error:noOverlap', 'No seeds in common between the two runs.');
    end

    % ---- per-streamline errors -------------------------------------------
    np = numel(common);
    pre_mean = nan(np,1); pre_max = nan(np,1);
    who_mean = nan(np,1); who_max = nan(np,1);
    endp     = nan(np,1); dlen    = nan(np,1);

    seedxyz_all = Tm.seed_points(it, :);
    for k = 1:np
        a = T{it(k)};  b = R{ir(k)};
        if size(a,1) < 2 || size(b,1) < 2, continue; end
        seed_xyz = seedxyz_all(k, :);

        sa = arclen(a); sb = arclen(b);
        dlen(k) = sa(end) - sb(end);

        % --- prefix: fixed arc EITHER SIDE OF THE SEED --------------------
        % Anchored at the seed, not at index 1. A track is
        % [reversed backward half; forward half], so index 1 is the far end of
        % the backward half - a TERMINATION point. Measuring from there put the
        % termination difference straight into the metric the prefix exists to
        % keep it out of. Locate the seed in each track and walk outward.
        ka = seed_pos_index(a, seed_xyz);
        kb = seed_pos_index(b, seed_xyz);

        % NOMINAL arc, not the cumulative chord length.
        %
        % Summing chords underestimates true arc by h^2/(24R^2) per unit length,
        % so two tracks stored at different step sizes end up on parameters that
        % disagree by O(h^2) - and comparing them "at matched arc" then reports
        % that disagreement as positional error. Measured here it was 1.2x the
        % predicted h^2 mismatch at every step size, and it was the second of two
        % O(h^2) artefacts pinning RK4 at order 2.
        %
        % The integrator solves dx/ds = v(x) with |v| = 1, so s IS arc length and
        % each step advances exactly h of it. The point n steps from the seed is
        % at arc n*h by construction - exact, with nothing estimated from the
        % geometry. Fall back to chord arc only when the step size was not
        % recorded with the run.
        if ~isempty(h_test)
            oa = ((1:size(a,1))' - ka) * h_test;
        else
            oa = sa - sa(ka);
        end
        if ~isempty(h_ref)
            ob = ((1:size(b,1))' - kb) * h_ref;
        else
            ob = sb - sb(kb);
        end
        if max(oa) >= prefix_arc && -min(oa) >= prefix_arc && ...
           max(ob) >= prefix_arc && -min(ob) >= prefix_arc
            % Compare the test track's OWN STORED POINTS against the reference.
            % Do not resample the test curve.
            %
            % This function used to resample BOTH curves linearly at 64 matched
            % arc positions, and that silently destroyed the measurement it
            % exists to make. A track stored every h voxels has a linear-
            % interpolation chord error of about h^2/(8R) for a curve of radius
            % R - an O(h^2) floor that has nothing to do with the integrator.
            % On this data (R ~ 9 voxels) that floor is 0.00347 at h = 0.5 and
            % the ladder measured 0.00352 for RK4: a ratio of 1.01, holding to
            % within 2% at every step size. Euler's own error is 37x above the
            % floor so it was reported correctly at order 1; RK2's is the same
            % order as the floor; and RK4's is well below it, so RK4 and RK2 both
            % reported the FLOOR - which is why their error constants agreed to
            % 6% and both read as order 2.
            %
            % The test points are exact integrator output, so they need no
            % interpolation at all. Only the reference is interpolated, and it is
            % dense (h_ref = 0.0078125 here) and interpolated with a spline, so
            % its own contribution is ~h_ref^4 and negligible.
            sel = abs(oa) <= prefix_arc;
            if nnz(sel) >= 4
                pa = a(sel, :);
                pb = resample_spline(b, ob, oa(sel));
                d  = sqrt(sum((pa-pb).^2, 2));
                pre_mean(k) = mean(d); pre_max(k) = max(d);
            end
        end

        % --- whole: over the common arc length ---------------------------
        L = min(sa(end), sb(end));
        if L > 0
            q  = linspace(0, L, n_samp)';
            pa = resample_at(a, sa, q);
            pb = resample_at(b, sb, q);
            d  = sqrt(sum((pa-pb).^2, 2));
            who_mean(k) = mean(d); who_max(k) = max(d);
        end
        endp(k) = norm(a(end,:) - b(end,:));
    end

    % Per-seed errors, not just their summary. A convergence ladder must be fitted
    % on ONE population: which streamlines are long enough to measure varies with
    % the step size (here by ~6% across the ladder), and if the seeds that only
    % qualify at fine h happen to be easier ones, the median falls for a reason
    % that has nothing to do with the integrator and the observed order is
    % inflated. Callers should intersect these seed lists across every level and
    % re-reduce on the common set.
    E.seed_index      = common(:)';
    E.prefix_per_seed = pre_mean(:)';

    E.prefix_arc      = prefix_arc;
    E.n_prefix_valid  = sum(~isnan(pre_mean));
    E.prefix          = stats(pre_mean);
    E.prefix_max      = stats(pre_max);
    E.whole           = stats(who_mean);
    E.whole_max       = stats(who_max);
    E.endpoint        = stats(endp);
    E.length_diff     = stats(dlen);
end

% =========================================================================
function [tracks, meta, step] = load_tracks(f)
    if isfolder(f)
        d = dir(fullfile(f, '**', 'tracks*.mat'));
        if isempty(d), error('nim_convergence_error:notFound', 'No tracks*.mat under %s', f); end
        [~,i] = max([d.datenum]);
        f = fullfile(d(i).folder, d(i).name);
    end
    S = load(f);
    tracks = S.tracks;
    step = [];
    if isfield(S, 'options') && isfield(S.options, 'step_size')
        step = S.options.step_size;
    end
    meta = [];
    if isfield(S, 'track_meta') && isfield(S.track_meta, 'seed_index')
        meta = S.track_meta;
    end
end

function k = seed_pos_index(p, seed_xyz)
% Index of the track point nearest the seed. The seed sits where the backward and
% forward halves were joined, so this recovers the true origin of integration.
    [~, k] = min(sum((p - seed_xyz).^2, 2));
end

function s = arclen(p)
    s = [0; cumsum(sqrt(sum(diff(p,1,1).^2, 2)))];
end

function q = resample_spline(p, s, at)
% Position along a DENSE polyline at given arc lengths, splined. Used only for
% the reference curve, whose point spacing is far finer than any test track, so
% the interpolation error here is orders below the quantity being measured.
    keep = [true; diff(s) > 0];
    q = interp1(s(keep), p(keep,:), at, 'spline', 'extrap');
end

function q = resample_at(p, s, at)
% Position along the polyline at given arc lengths. Duplicate knots are dropped
% so interp1 has a strictly increasing grid.
    keep = [true; diff(s) > 0];
    q = interp1(s(keep), p(keep,:), at, 'linear', 'extrap');
end

function st = stats(v)
    v = v(~isnan(v));
    if isempty(v)
        st = struct('n',0,'median',NaN,'p95',NaN,'max',NaN);
        return;
    end
    st = struct('n', numel(v), 'median', median(v), ...
                'p95', prctile(v,95), 'max', max(v));
end

function v = getf(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
