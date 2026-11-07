function nim = nim_field(nim, options, nim_path)
%NIM_FIELD  Provision the DIRECTION FIELD a tractography run needs (step 2).
%
%   nim = nim_field(nim, options, nim_path)
%
%   THE BOUNDARY RULE THIS FILE EXISTS TO ENFORCE:
%     * the nim on disk is the DATASET (img, evec/eval, FA, masks, parcellation)
%       and nothing else - main.m never reads config.tractography;
%     * everything that depends on config.tractography (the field here) is
%       built PER RUN by runTractography;
%     * a per-run product is cached to a sidecar next to the nim only when it
%       costs more than a minute to rebuild.
%
%   Measured on data/ismrm2015/ismrm2015.mat (90x108x90, 33 volumes):
%     field: csd   nim_csd              122 s   -> cached by dataset and CSD settings
%
%   Cases:
%     'dti' (default) - nothing to do; the principal eigenvector is dataset (a).
%     'csd'           - FOD peaks (nim.peaks/npeaks/peak_w, optional fod_sh) via
%                       nim_csd or a supplied csd_peaks_file. Needed by ANY
%                       tracker running field=csd (hinec and stitching), so it
%                       is provisioned before the dispatch. Supplied peaks use
%                       the dataset's voxel-axis direction coordinates.
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
case 'lift'
    nim = provision_lift(nim, options, nim_path);
otherwise
    fprintf('field=%s: principal eigenvector from the nim; nothing to provision\n', fld);
end
end

% ---------------------------------------------------------------------------

function nim = provision_csd(nim, options, nim_path)
% CSD FOD peaks are needed by ANY tracker running field=csd (hinec and stitching), so
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

function nim = provision_lift(nim, options, nim_path)
% The lifts (nim_lift_field) for field: lift. Every setting that changes them is in the cache name.
keys = {'smooth', 'bridge', 'rel', 'mass_min', 'max_lifts', 'sep', 'ndirs', 'kappa', 'ref_fa_pct', 'fit', 'iso', 'sparsity'};
lo = struct('smooth', 0, 'bridge', 2, 'rel', 0.1, 'mass_min', 0.1, 'max_lifts', 3, 'sep', 15, ...
            'ndirs', 300, 'kappa', 20, 'ref_fa_pct', 90, 'fit', false, 'iso', 'auto', 'sparsity', 0.1);
for k = 1:numel(keys)
    ok_ = ['lift_' keys{k}];
    if isfield(options, ok_) && ~isempty(options.(ok_)), lo.(keys{k}) = options.(ok_); end
end
if lo.smooth > 0 && isfield(options, 'lift_smooth_kernel') && strcmpi(options.lift_smooth_kernel, 'curved')
    % the curved kernel's settings enter the cache name only when it is used, so the caches of the straight
    % kernel keep their names
    lo.smooth_kernel = 'curved';
    for k = {'smooth_kappa', 'smooth_normals'}
        if isfield(options, ['lift_' k{1}]) && ~isempty(options.(['lift_' k{1}])), lo.(k{1}) = options.(['lift_' k{1}]); end
    end
end
if isfield(options, 'lift_continuation') && ~isempty(options.lift_continuation) && options.lift_continuation > 0
    lo.continuation = options.lift_continuation;   % in the cache name only when on (the old caches keep their names)
end
cache = lift_cache_path(nim_path, lo);
if isfield(options, 'lift_workers') && ~isempty(options.lift_workers), lo.workers = options.lift_workers; end   % not a setting of the lifts
if ~isempty(cache) && isfile(cache)
    fprintf('field=lift: loading cached lifts from %s\n', cache);
    Sc = load(cache, 'U', 'ok', 'mass', 'lambda');
else
    fprintf('field=lift: building the lifts (nim_lift_field)...\n');
    L = nim_lift_field(nim, lo);
    Sc = struct('U', L.U, 'ok', L.ok, 'mass', L.mass, 'lambda', L.lambda);
    if ~isempty(cache)
        try
            U = L.U; ok = L.ok; mass = L.mass; lambda = L.lambda; %#ok<NASGU>
            save(cache, 'U', 'ok', 'mass', 'lambda', '-v7.3');
            fprintf('  cached lifts -> %s\n', cache);
        catch
            % non-fatal: proceed without caching
        end
    end
end
nim.lift_U = Sc.U; nim.lift_ok = Sc.ok; nim.lift_mass = Sc.mass; nim.lift_lambda = Sc.lambda;
end

function path = lift_cache_path(nim_path, lo)
% <nim>_lift_<hash>.mat, the hash (MD5, 16 hex digits) of every lift setting, the code version of
% nim_lift_field and the nim file itself (size and modification time): a re-preprocessed dataset or a
% changed builder never reuses stale lifts, and the name stays short (a per-key name exceeds NAME_MAX).
path = '';
if isempty(nim_path) || ~isfile(nim_path), return; end
fi = dir(nim_path); bi = dir(which('nim_lift_field'));
key = sprintf('%s|%d|%.10f|%d|%.10f', jsonencode(lo), fi.bytes, fi.datenum, bi.bytes, bi.datenum);
md = java.security.MessageDigest.getInstance('MD5');
h = reshape(dec2hex(typecast(md.digest(uint8(key)), 'uint8'))', 1, []);
path = regexprep(nim_path, '\.mat$', ['_lift_' lower(h(1:16)) '.mat']);
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
