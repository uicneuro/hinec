function nim = nim_field(nim, options, nim_path)
%NIM_FIELD  Provision the DIRECTION FIELD a tractography run needs (step 2).
%
%   nim = nim_field(nim, options, nim_path)
%
%   THE BOUNDARY RULE THIS FILE EXISTS TO ENFORCE:
%     * the nim on disk is the DATASET (img, evec/eval, FA, masks, parcellation)
%       and nothing else - main.m never reads config.tractography;
%     * everything that depends on config.tractography (the field here, the MMF
%       geometry in step 3) is built PER RUN by runTractography;
%     * a per-run product is cached to a sidecar next to the nim only when it
%       costs more than a minute to rebuild.
%
%   Measured on data/ismrm2015/ismrm2015.mat (90x108x90, 33 volumes):
%     field: csd   nim_csd              122 s   -> cached by dataset and CSD settings
%     field: dwi   nim_mmf_from_dwi     286 s   -> cached by dataset and fit settings
%                                        (45 s per-voxel fit + 240 s of the 3
%                                         continuity sweeps)
%     geometry     nim_mmf_geometry     4.9 s (dti) / 8.1 s (csd)  -> never cached
%
%   Cases:
%     'dti' (default) - nothing to do; the principal eigenvector is dataset (a).
%     'csd'           - FOD peaks (nim.peaks/npeaks/peak_w, optional fod_sh) via
%                       nim_csd or a supplied csd_peaks_file. Needed by ANY
%                       tracker running field=csd (hinec AND mmf), so it is
%                       provisioned before the dispatch. Supplied peaks use
%                       the dataset's voxel-axis direction coordinates.
%     'dwi'           - joint frame+curvature fit to the raw DW signal
%                       (nim.mmf_e1_dwi / mmf_kappa_dwi) via nim_mmf_from_dwi,
%                       consumed by nim_mmf_geometry in step 3.
%
%   nim_path is the path of the nim on disk; the sidecars are written beside it.

if nargin < 2 || isempty(options), options = struct(); end
if nargin < 3, nim_path = ''; end

fld = 'dti';
if isfield(options, 'field') && ~isempty(options.field)
    fld = lower(char(string(options.field)));
end

switch fld
case 'csd'
    nim = provision_csd(nim, options, nim_path);
case 'dwi'
    nim = provision_dwi(nim, options, nim_path);
otherwise
    fprintf('field=%s: principal eigenvector from the nim; nothing to provision\n', fld);
end
end

% ---------------------------------------------------------------------------

function nim = provision_csd(nim, options, nim_path)
% CSD FOD peaks are needed by ANY tracker running field=csd (hinec AND mmf), so
% provision them BEFORE the algorithm dispatch. The cached sidecar includes every
% CSD setting, so a requested configuration never inherits stale embedded peaks.

source_file = getfield_default(options, 'csd_peaks_file', '');
if ~isempty(source_file)
    source_file = char(string(source_file));
    if ~isfile(source_file)
        error('nim_field:missingPeaksFile', 'CSD peaks file not found: %s', source_file);
    end
    source = load(source_file, 'peaks', 'npeaks', 'peak_w');
    required = {'peaks', 'npeaks', 'peak_w'};
    if ~all(isfield(source, required))
        error('nim_field:invalidPeaksFile', ...
            'CSD peaks file must contain peaks, npeaks, and peak_w.');
    end
    dims = size(nim.FA);
    npeak_slots = size(source.peaks, 4);
    if ~isequal(size(source.peaks, [1 2 3]), dims) || ...
            size(source.peaks, 5) ~= 3 || npeak_slots < 1 || ...
            ~isequal(size(source.npeaks), dims) || ...
            ~isequal(size(source.peak_w, [1 2 3]), dims) || ...
            size(source.peak_w, 4) ~= npeak_slots || ...
            any(source.npeaks(:) < 0) || ...
            any(source.npeaks(:) > npeak_slots) || ...
            any(source.npeaks(:) ~= floor(source.npeaks(:)))
        error('nim_field:invalidPeaksFile', ...
            'CSD peaks file dimensions or values do not match the dataset.');
    end
    active = reshape(1:npeak_slots, [1 1 1 npeak_slots]) <= source.npeaks;
    active_vectors = repmat(active, [1 1 1 1 3]);
    if any(~isfinite(source.peaks(active_vectors))) || ...
            any(~isfinite(source.peak_w(active))) || ...
            any(source.peak_w(active) < 0)
        error('nim_field:invalidPeaksFile', ...
            'CSD peaks file has invalid values in active peak slots.');
    end
    % MRtrix peak exporters may encode unused slots as NaN. They are never
    % selected (npeaks excludes them), so normalize them to zero on import.
    source.peaks(~active_vectors) = 0;
    source.peak_w(~active) = 0;
    nim.peaks = source.peaks;
    nim.npeaks = source.npeaks;
    nim.peak_w = source.peak_w;
    fprintf('field=csd: using supplied deterministic peak field from %s\n', source_file);
    return;
end

% The peak field changes with every option below. A dataset-only cache name
% silently reused peaks fitted under a different configuration.
csd_opts = struct('lmax', 4, 'n_iter', 50, 'peak_thresh', 0.2, ...
                  'peak_min_sep', 45, 'max_peaks', 3);
csd_keys = {'lmax', 'n_iter', 'peak_thresh', 'peak_min_sep', 'max_peaks'};
for ci = 1:numel(csd_keys)
    ck = ['csd_' csd_keys{ci}];
    if isfield(options, ck) && ~isempty(options.(ck))
        csd_opts.(csd_keys{ci}) = options.(ck);
    end
end
csd_cache = field_cache_path(nim_path, 'csd', csd_opts, csd_keys);
if ~isempty(csd_cache) && isfile(csd_cache)
    fprintf('field=csd: loading cached CSD FOD from %s\n', csd_cache);
    Sc = load(csd_cache);
    nim.peaks = Sc.peaks; nim.npeaks = Sc.npeaks; nim.peak_w = Sc.peak_w;
    if isfield(Sc, 'fod_sh'), nim.fod_sh = Sc.fod_sh; end
else
    fprintf('field=csd: computing CSD FOD peaks (nim_csd)...\n');
    nim = nim_csd(nim, csd_opts);
    if ~isempty(csd_cache)
        try
            peaks = nim.peaks; npeaks = nim.npeaks; peak_w = nim.peak_w; %#ok<NASGU>
            if isfield(nim, 'fod_sh')
                fod_sh = nim.fod_sh; %#ok<NASGU>
                save(csd_cache, 'peaks', 'npeaks', 'peak_w', 'fod_sh', '-v7.3');
            else
                save(csd_cache, 'peaks', 'npeaks', 'peak_w', '-v7.3');
            end
            fprintf('  cached CSD FOD -> %s\n', csd_cache);
        catch
            % non-fatal: proceed without caching
        end
    end
end
end

% ---------------------------------------------------------------------------

function nim = provision_dwi(nim, options, nim_path)
% Direct route: e1 and the curvature vector are fitted jointly to the raw DW
% signal, so the curvature was never differenced out of an already-collapsed
% direction field. nim_mmf_geometry (step 3) only CONSUMES these two fields.
%
% MEASURED at 286 s on data/ismrm2015/ismrm2015.mat (81328 voxels, 32 gradients:
% 45 s of per-voxel fitting, 240 s for the 3 continuity sweeps), which is over the
% one-minute line, so it gets a parameter-specific sidecar next to the CSD
% one. Embedded fields cannot prove which fit settings produced them.

dwi_opts = struct('termination_fa', getfield_default(options,'termination_fa',0.08), ...
                  'continuity_sweeps', getfield_default(options,'continuity_sweeps',3), ...
                  'continuity_weight', getfield_default(options,'continuity_weight',1.0));
dwi_cache = field_cache_path(nim_path, 'dwi', dwi_opts, fieldnames(dwi_opts));
if ~isempty(dwi_cache) && isfile(dwi_cache)
    fprintf('field=dwi: loading cached DW-fitted frame from %s\n', dwi_cache);
    Sd = load(dwi_cache);
    nim.mmf_e1_dwi = Sd.mmf_e1_dwi; nim.mmf_kappa_dwi = Sd.mmf_kappa_dwi;
    if isfield(Sd, 'mmf_dwi_resid'), nim.mmf_dwi_resid = Sd.mmf_dwi_resid; end
    return;
end

fprintf('field=dwi: fitting frame + curvature to the raw DW signal (nim_mmf_from_dwi)...\n');
t_dwi = tic;
nim = nim_mmf_from_dwi(nim, options);
fprintf('  nim_mmf_from_dwi: %.1f s\n', toc(t_dwi));
if ~isempty(dwi_cache)
    try
        mmf_e1_dwi = nim.mmf_e1_dwi; mmf_kappa_dwi = nim.mmf_kappa_dwi; %#ok<NASGU>
        if isfield(nim, 'mmf_dwi_resid')
            mmf_dwi_resid = nim.mmf_dwi_resid; %#ok<NASGU>
            save(dwi_cache, 'mmf_e1_dwi', 'mmf_kappa_dwi', 'mmf_dwi_resid', '-v7.3');
        else
            save(dwi_cache, 'mmf_e1_dwi', 'mmf_kappa_dwi', '-v7.3');
        end
        fprintf('  cached DW-fitted frame -> %s\n', dwi_cache);
    catch
        % non-fatal: proceed without caching
    end
end
end

function path = field_cache_path(nim_path, field, values, keys)
% IEEE-754 hex keeps distinct numeric settings from sharing a cache name.
path = '';
if isempty(nim_path), return; end
parts = cell(1,numel(keys));
for k = 1:numel(keys)
    key = keys{k};
    parts{k} = [key '_' num2hex(double(values.(key)))];
end
tag = strjoin(parts,'_');
path = regexprep(nim_path,'\.mat$',['_' field '_' tag '.mat']);
end

function value = getfield_default(s,key,fallback)
if isfield(s,key) && ~isempty(s.(key)), value=s.(key); else, value=fallback; end
end
