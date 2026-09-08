function R = nim_track_diagnose(run_dir, nim, gt_file, opts)
% nim_track_diagnose: Attribute streamline error to the stage that caused it.
%
%   R = nim_track_diagnose(run_dir, nim, 'bundles/Cingulum_right.trk')
%
% Consumes the PER-STEP trace written by tractography.debug.trace and answers
% one question: when a streamline stops following the true tract, WHICH PART of
% the tracking equation is responsible?
%
% THE DECOMPOSITION. At each step the tracker turns a true fibre direction into
% a direction it actually steps along, through three stages. Each stage can be
% measured separately because each has an observable input and output:
%
%   d_true  --(tensor fit: model + acquisition)-->  v1 at the nearest voxel
%   v1      --(interpolation of the dyadic)------>  d_used, the traced direction
%   d_used  --(Runge-Kutta stage averaging)------>  the realised chord
%
% giving three angles per step:
%
%   e_model  = angle(d_true, v1_raw)   what a single tensor cannot represent
%   e_interp = angle(v1_raw, d_used)   what the interpolation kernel changed
%   e_total  = angle(d_true, d_used)   the error the streamline actually carries
%
% and a fourth, scalar, term for the integrator: chord/h, the fraction of the
% requested arc the step actually advanced. RK4 averages four stage vectors, so
% when they disagree the chord falls short of h; chord/h near 1 means the step
% was effectively straight and the integrator contributed nothing.
%
% All comparisons use abs(dot(.)) because v1 and the ground-truth tangent are
% LINE fields - their sign is arbitrary, and a signed comparison reports 180
% degrees of error for two vectors describing the same orientation.
%
% ATTRIBUTION. A departure is the first step whose e_total exceeds opts.thresh
% (default 45 deg) and stays above it for at least opts.persist voxels of arc -
% the persistence requirement matters, because roughly half of bare threshold
% crossings are the tangent wobbling across the line and recover on their own.
% Each departure is then charged to a stage:
%
%   model        e_model  > thresh   the tensor at that voxel is already wrong
%   interpolation e_model <= thresh but e_interp > thresh
%   sign          a flip was applied at that step and neither of the above
%   other        neither stage is individually over threshold
%
% Steps where the ground truth has no support are left unattributed rather than
% guessed: off-bundle there is no true direction to compare against.

    if nargin < 4, opts = struct(); end
    thresh  = getf(opts, 'thresh', 45);
    persist = getf(opts, 'persist', 2);     % voxels of arc
    dims = size(nim.FA);

    [T, meta, topts] = load_run(run_dir);
    if ~isfield(meta, 'trace') || isempty(meta.trace)
        error('nim_track_diagnose:noTrace', ...
            ['This run has no per-step trace. Re-run with debug.trace: true - the trace ' ...
             'cannot be reconstructed afterwards, because output.arc_step decimates the ' ...
             'saved polyline and the direction actually used is never stored.']);
    end

    [GD, GC] = gt_direction_field(gt_file, nim, dims);
    [V1, IN] = field_interpolants(nim, topts, dims);
    [CL, CP] = tensor_shape(nim, dims);

    n = numel(meta.trace);
    rows = struct('e_model',{},'e_interp',{},'e_total',{},'chord_ratio',{}, ...
                  'fa',{},'cl',{},'cp',{},'flip',{},'arc',{},'pos',{},'seed',{});
    dep = struct('cause',{},'e_model',{},'e_interp',{},'cl',{},'cp',{},'fa',{}, ...
                 'arc',{},'pos',{},'seed',{},'termination',{});
    nsteps = 0; nsupported = 0;

    for k = 1:n
        for arm = {'forward','backward'}
            tr = meta.trace(k).(arm{1});
            if isempty(tr) || tr.n_steps < 3, continue; end
            arc = [0; cumsum(tr.chord(1:end-1))];
            arc(isnan(arc)) = 0;

            em = nan(tr.n_steps,1); ei = em; et = em; cr = em;
            cl = em; cp = em;
            for i = 1:tr.n_steps
                p = tr.pos(i,:); d = tr.dir(i,:);
                if any(isnan(p)) || any(isnan(d)), continue; end
                v = min(max(round(p),1),dims); ix = sub2ind(dims,v(1),v(2),v(3));
                cl(i) = CL(ix); cp(i) = CP(ix);
                raw = [V1{1}(ix) V1{2}(ix) V1{3}(ix)];
                nr = norm(raw);
                if nr > 0
                    raw = raw/nr;
                    ei(i) = ang(raw, d);
                end
                if GC(ix) > 0
                    g = GD(ix,:);
                    et(i) = ang(g, d);
                    if nr > 0, em(i) = ang(g, raw); end
                end
                if ~isnan(tr.h(i)) && tr.h(i) > 0, cr(i) = tr.chord(i)/tr.h(i); end
            end
            nsteps = nsteps + tr.n_steps;
            nsupported = nsupported + nnz(~isnan(et));

            % first PERSISTENT departure
            bad = et > thresh;
            kd = [];
            for i = 1:tr.n_steps
                if ~bad(i) || isnan(et(i)), continue; end
                j = find(arc >= arc(i) + persist, 1, 'first');
                if isempty(j), break; end
                w = i:j; w = w(~isnan(et(w)));
                if ~isempty(w) && all(et(w) > thresh), kd = i; break; end
            end
            if ~isempty(kd)
                if em(kd) > thresh,        c = 'model';
                elseif ei(kd) > thresh,    c = 'interpolation';
                elseif tr.flip(kd),        c = 'sign';
                else,                      c = 'other';
                end
                dep(end+1) = struct('cause',c,'e_model',em(kd),'e_interp',ei(kd), ...
                    'cl',cl(kd),'cp',cp(kd),'fa',tr.fa(kd),'arc',arc(kd), ...
                    'pos',tr.pos(kd,:),'seed',meta.trace(k).seed_index, ...
                    'termination',tr.termination); %#ok<AGROW>
            end
            good = ~isnan(et);
            if any(good)
                rows(end+1) = struct('e_model',em(good),'e_interp',ei(good), ...
                    'e_total',et(good),'chord_ratio',cr(good),'fa',tr.fa(good), ...
                    'cl',cl(good),'cp',cp(good),'flip',tr.flip(good), ...
                    'arc',arc(good),'pos',tr.pos(good,:), ...
                    'seed',meta.trace(k).seed_index); %#ok<AGROW>
            end
        end
    end

    R = struct();
    R.n_arms_traced = n*2;
    R.n_steps = nsteps;
    R.n_steps_on_gt = nsupported;
    R.steps = rows;
    R.departures = dep;
    R.thresh = thresh; R.persist = persist;
    R.summary = summarise(rows, dep, thresh);
end

% =========================================================================
function S = summarise(rows, dep, thresh)
    cat_ = @(f) vertcat(rows.(f));
    S = struct();
    if isempty(rows), return; end
    em = cat_('e_model'); ei = cat_('e_interp'); et = cat_('e_total'); cr = cat_('chord_ratio');
    S.median_e_model  = median(em,'omitnan');
    S.median_e_interp = median(ei,'omitnan');
    S.median_e_total  = median(et,'omitnan');
    S.p95_e_model     = prctile(em(~isnan(em)),95);
    S.p95_e_interp    = prctile(ei(~isnan(ei)),95);
    S.frac_model_over  = mean(em > thresh,'omitnan');
    S.frac_interp_over = mean(ei > thresh,'omitnan');
    S.median_chord_ratio = median(cr,'omitnan');
    S.n_departures = numel(dep);
    if ~isempty(dep)
        c = {dep.cause};
        for u = {'model','interpolation','sign','other'}
            S.(['dep_' u{1}]) = sum(strcmp(c,u{1}));
        end
    end
end

function [GD, GC] = gt_direction_field(gt_file, nim, dims)
% Per-voxel ground-truth orientation, from the average tangent DYADIC. Averaging
% the tangents themselves would cancel: they are a line field with arbitrary
% per-point sign.
    NV = prod(dims);
    G = nim_read_trk(gt_file, gt_ref(nim));
    A = zeros(NV,6); GC = zeros(NV,1);
    for t = 1:numel(G)
        g = G{t}; if size(g,1) < 3, continue; end
        d = g(3:end,:) - g(1:end-2,:); d = d./max(vecnorm(d,2,2),eps);
        q = round(g(2:end-1,:));
        m = all(q>=1,2)&q(:,1)<=dims(1)&q(:,2)<=dims(2)&q(:,3)<=dims(3);
        q = q(m,:); d = d(m,:); if isempty(q), continue; end
        ix = sub2ind(dims,q(:,1),q(:,2),q(:,3));
        A(:,1)=A(:,1)+accumarray(ix,d(:,1).^2,[NV 1]); A(:,2)=A(:,2)+accumarray(ix,d(:,2).^2,[NV 1]);
        A(:,3)=A(:,3)+accumarray(ix,d(:,3).^2,[NV 1]); A(:,4)=A(:,4)+accumarray(ix,d(:,1).*d(:,2),[NV 1]);
        A(:,5)=A(:,5)+accumarray(ix,d(:,1).*d(:,3),[NV 1]); A(:,6)=A(:,6)+accumarray(ix,d(:,2).*d(:,3),[NV 1]);
        GC = GC + accumarray(ix,1,[NV 1]);
    end
    GD = zeros(NV,3); gi = find(GC>0);
    for i = 1:numel(gi)
        j = gi(i);
        GD(j,:) = nim_principal_dir(A(j,1),A(j,2),A(j,3),A(j,4),A(j,5),A(j,6));
    end
end

function f = gt_ref(nim)
    f = 'data/ismrm2015/ismrm2015_dwi_ref.nii.gz';
    if isfield(nim,'dwi_ref_file') && ~isempty(nim.dwi_ref_file), f = nim.dwi_ref_file; end
end

function [V1, IN] = field_interpolants(nim, topts, dims)
    V1 = {squeeze(nim.evec(:,:,:,1,1)), squeeze(nim.evec(:,:,:,2,1)), squeeze(nim.evec(:,:,:,3,1))};
    IN = topts;  % kept for provenance; the raw field is what the comparison needs
end

function [CL, CP] = tensor_shape(nim, dims)
    L = sort(reshape(nim.eval,[],3), 2, 'descend');
    s = max(L(:,1), eps);
    CL = reshape((L(:,1)-L(:,2))./s, dims);
    CP = reshape((L(:,2)-L(:,3))./s, dims);
end

function [T, meta, opts] = load_run(run_dir)
    if isfolder(run_dir)
        d = dir(fullfile(run_dir,'**','tracks*.mat'));
        if isempty(d), error('nim_track_diagnose:notFound','No tracks*.mat under %s', run_dir); end
        [~,i] = max([d.datenum]); f = fullfile(d(i).folder, d(i).name);
    else
        f = run_dir;
    end
    S = load(f);
    T = S.tracks; meta = []; opts = struct();
    if isfield(S,'track_meta'), meta = S.track_meta; end
    if isfield(S,'options'), opts = S.options; end
end

function a = ang(u, v)
    a = acosd(min(abs(dot(u(:), v(:))), 1));    % LINE fields: sign is meaningless
end

function v = getf(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
