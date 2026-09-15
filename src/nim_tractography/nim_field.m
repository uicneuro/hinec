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
%     field: csd   nim_csd              122 s   -> CACHED to <nim>_csd.mat
%     field: dwi   nim_mmf_from_dwi     286 s   -> CACHED to <nim>_dwi.mat
%                                        (45 s per-voxel fit + 240 s of the 3
%                                         continuity sweeps)
%     geometry     nim_mmf_geometry     4.9 s (dti) / 8.1 s (csd)  -> never cached
%
%   Cases:
%     'dti' (default) - nothing to do; the principal eigenvector is dataset (a).
%     'csd'           - FOD peaks (nim.peaks/npeaks/peak_w, optional fod_sh) via
%                       nim_csd. Needed by ANY tracker running field=csd (hinec
%                       AND mmf), so it is provisioned before the dispatch.
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
% provision them BEFORE the algorithm dispatch. Compute with nim_csd when the config
% sets field=csd, cached next to the source nim (<source>_csd.mat) so it is computed
% once per preprocessed dataset and reused by every tractography config.
if isfield(nim, 'peaks'), return; end

csd_cache = '';
if ~isempty(nim_path)
    csd_cache = regexprep(nim_path, '\.mat$', '_csd.mat');
end
if ~isempty(csd_cache) && isfile(csd_cache)
    fprintf('field=csd: loading cached CSD FOD from %s\n', csd_cache);
    Sc = load(csd_cache);
    nim.peaks = Sc.peaks; nim.npeaks = Sc.npeaks; nim.peak_w = Sc.peak_w;
    if isfield(Sc, 'fod_sh'), nim.fod_sh = Sc.fod_sh; end
else
    fprintf('field=csd: computing CSD FOD peaks (nim_csd)...\n');
    % lmax 4 and peak_thresh 0.2 are set from this acquisition, not convention;
    % see nim_config_schema for the measurements behind both.
    csd_opts = struct('lmax', 4, 'n_iter', 50, 'peak_thresh', 0.2, ...
                      'peak_min_sep', 45, 'max_peaks', 3);
    csd_keys = {'lmax', 'n_iter', 'peak_thresh', 'peak_min_sep', 'max_peaks'};
    for ci = 1:numel(csd_keys)
        ck = ['csd_' csd_keys{ci}];
        if isfield(options, ck) && ~isempty(options.(ck))
            csd_opts.(csd_keys{ci}) = options.(ck);
        end
    end
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
% one-minute line, so it gets a sidecar next to the CSD one - same style, same
% rule: computed once per preprocessed dataset, reused by every dwi config.
if isfield(nim, 'mmf_kappa_dwi') && ~isempty(nim.mmf_kappa_dwi) ...
        && isfield(nim, 'mmf_e1_dwi') && ~isempty(nim.mmf_e1_dwi)
    return;
end

dwi_cache = '';
if ~isempty(nim_path)
    dwi_cache = regexprep(nim_path, '\.mat$', '_dwi.mat');
end
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
