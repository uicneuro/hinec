function txt = nim_describe_tracker_input(nim, options, algorithm, out_file)
% nim_describe_tracker_input: Log EXACTLY what a tracker is handed at step 5.
%
%   txt = nim_describe_tracker_input(nim, options, algorithm)
%   txt = nim_describe_tracker_input(nim, options, algorithm, out_file)
%
% runTractography calls this immediately before dispatching to the tracker, so
% the log of every run (and <run_dir>/tractography/tracker_input.txt) records
% the two arguments the tracker receives - every nim field and every option,
% with its class, size, where it came from and what it means - followed by the
% contract the tracker must satisfy on the way back.
%
% This exists so that someone writing a NEW tracker can read one run log and
% know the data structure they will be given, without reverse-engineering
% runTractography, nim_field and nim_config_to_options. The same tables, with
% layouts and worked examples, are in docs/TRACKER_INTERFACE.md; this function
% is the executable, always-current version of that document.
%
% Anything present that is not in the dictionaries below is printed as
% "(undocumented)" rather than hidden: an undocumented field in the hand-off is
% a defect in the boundary, and the log should say so.

if nargin < 3 || isempty(algorithm), algorithm = ''; end
if nargin < 4, out_file = ''; end

L = {};

% ---------------------------------------------------------------- the call
L = put(L, '=== Tracker input (step 5 hand-off) ===');
L = put(L, 'Call: [tracks, meta] = nim_tractography_%s(nim, options)', lower(char(string(algorithm))));
L = put(L, 'Coordinates: voxel space, 1-based (voxel (1,1,1) is centred at [1 1 1]);');
L = put(L, '             positions and step sizes are in VOXELS, angles in degrees.');
L = put(L, '');

% ---------------------------------------------------------------- nim
L = put(L, '--- nim: fields present (%d) ---', numel(fieldnames(nim)));
[groups, dict] = nim_dictionary();
present = fieldnames(nim);
seen = false(size(present));
for g = 1:numel(groups)
    names = groups(g).fields;
    rows = {};
    for k = 1:numel(names)
        idx = find(strcmp(present, names{k}), 1);
        if isempty(idx), continue; end
        seen(idx) = true;
        rows(end+1, :) = {names{k}, describe_value(nim.(names{k})), dict.(names{k})}; %#ok<AGROW>
    end
    if isempty(rows), continue; end
    L = put(L, '[%s]', groups(g).title);
    L = put_rows(L, rows);
end
extra = present(~seen);
if ~isempty(extra)
    L = put(L, '[other - undocumented, not read by any shipped tracker]');
    rows = cell(numel(extra), 3);
    for k = 1:numel(extra)
        rows(k, :) = {extra{k}, describe_value(nim.(extra{k})), '(undocumented)'};
    end
    L = put_rows(L, rows);
end
L = put(L, '');

% ---------------------------------------------------------------- options
L = put(L, '--- options: fields present (%d) ---', numel(fieldnames(options)));
L = put(L, 'Flat option name <- config key it came from : meaning. Values as handed over.');
[ogroups, odict] = options_dictionary();
present = fieldnames(options);
seen = false(size(present));
for g = 1:numel(ogroups)
    names = ogroups(g).fields;
    rows = {};
    for k = 1:numel(names)
        idx = find(strcmp(present, names{k}), 1);
        if isempty(idx), continue; end
        seen(idx) = true;
        d = odict.(names{k});
        rows(end+1, :) = {names{k}, format_value(options.(names{k})), sprintf('<- %s : %s', d{1}, d{2})}; %#ok<AGROW>
    end
    if isempty(rows), continue; end
    L = put(L, '[%s]', ogroups(g).title);
    L = put_rows(L, rows);
end
extra = present(~seen);
if ~isempty(extra)
    L = put(L, '[other - undocumented]');
    rows = cell(numel(extra), 3);
    for k = 1:numel(extra)
        rows(k, :) = {extra{k}, format_value(options.(extra{k})), '(undocumented)'};
    end
    L = put_rows(L, rows);
end
L = put(L, '');

% ---------------------------------------------------------------- output contract
L = put(L, '--- output contract (what the tracker must return) ---');
L = put(L, 'tracks : cell {T x 1}. tracks{t} is an N_t x 3 double of voxel coordinates,');
L = put(L, '         one row per point, ordered [backward half reversed; seed; forward half]');
L = put(L, '         so that the polyline reads end-to-end. Keep a track only if its');
L = put(L, '         chord length sum(|diff|) >= options.min_length; drop the rest.');
L = put(L, 'meta   : struct (may be empty). Recognised fields:');
L = put(L, '         seed_index  1 x T  index into the seed list of each KEPT track');
L = put(L, '         seed_points T x 3  the seed position of each kept track');
L = put(L, '         n_seeds     scalar number of seeds attempted');
L = put(L, '         trace       optional per-step record (see debug.trace)');
L = put(L, 'runTractography then applies the ROI filter (step 6), decimates to');
L = put(L, 'output.arc_step, and saves tracks/options/elapsed_time/algorithm/track_meta.');
L = put(L, '=======================================');

txt = strjoin(L, newline);
fprintf('%s\n', txt);

if ~isempty(out_file)
    fid = fopen(out_file, 'w');
    if fid > 0
        fprintf(fid, '%s\n', txt);
        fclose(fid);
        fprintf('Tracker input description written to %s\n', out_file);
    else
        warning('nim_describe_tracker_input:write', 'Could not write %s', out_file);
    end
end
end

% =========================================================================
function L = put(L, varargin)
    L{end+1, 1} = sprintf(varargin{:});
end

function L = put_rows(L, rows)
% Three aligned columns: name, class/size or value, meaning.
    w1 = max(cellfun(@length, rows(:, 1)));
    w2 = max(cellfun(@length, rows(:, 2)));
    for k = 1:size(rows, 1)
        L{end+1, 1} = sprintf('  %-*s  %-*s  %s', w1, rows{k, 1}, w2, rows{k, 2}, rows{k, 3});
    end
end

function s = describe_value(v)
    sz = size(v);
    s = sprintf('%s [%s]', class(v), strjoin(arrayfun(@(d) sprintf('%d', d), sz, 'uni', 0), ' '));
    if isa(v, 'containers.Map'), s = sprintf('containers.Map (%d keys)', v.Count); end
end

function s = format_value(v)
% Options are mostly scalars and short strings; print those verbatim so the log
% is a record of the run. Anything bigger is summarised by class and size.
    if ischar(v) || isstring(v)
        s = sprintf('''%s''', char(string(v)));
    elseif islogical(v) && isscalar(v)
        if v, s = 'true'; else, s = 'false'; end
    elseif isnumeric(v) && isempty(v)
        s = '[]';
    elseif isnumeric(v) && numel(v) <= 6
        s = mat2str(v, 6);
    elseif iscell(v)
        if isempty(v)
            s = '{}';
        elseif all(cellfun(@(c) ischar(c) || isstring(c), v)) && numel(v) <= 8
            s = ['{' strjoin(cellfun(@(c) char(string(c)), v(:)', 'uni', 0), ', ') '}'];
        else
            s = sprintf('cell {%d}', numel(v));
        end
    elseif isstruct(v)
        s = sprintf('struct (%d fields)', numel(fieldnames(v)));
    else
        s = describe_value(v);
    end
end

% =========================================================================
function [groups, d] = nim_dictionary()
% What each nim field is. Layouts: X Y Z are the volume dimensions; the trailing
% dimensions are documented per field. Only fields actually read by a shipped
% tracker or needed to interpret one are described.
    d = struct();
    % dataset - written by main.m, on disk, independent of config.tractography
    d.hdr        = 'NIfTI header (niftiinfo): ImageSize, PixelDimensions, Transform (voxel->world affine)';
    d.xdim       = 'volume size along X'; d.ydim = 'along Y'; d.zdim = 'along Z';
    d.size3      = 'xdim*ydim*zdim';
    d.img        = '[X Y Z V] raw DWI signal, V volumes (b0 + diffusion weighted)';
    d.bval       = '[V 1] b-value of each volume';
    d.bvec       = '[V 3] gradient direction of each volume (unit)';
    d.size_b0    = 'number of b0 volumes'; d.size_bi = 'number of diffusion-weighted volumes';
    d.img_b0     = '[X Y Z] mean b0 image'; d.img_bi = '[X Y Z size_bi] diffusion-weighted volumes only';
    d.thrsh_b0   = 'b-value below which a volume counts as b0';
    d.mask       = '[X Y Z] brain mask, 1 inside. Tracking stops on leaving it';
    d.DT         = '[X Y Z 6] diffusion tensor as [Dxx Dyy Dzz Dxy Dyz Dxz] (nim_reshape_d)';
    d.evec       = '[X Y Z 3 3] eigenvectors: evec(x,y,z,:,k) is the k-th, k=1 principal. Sign is ARBITRARY (a line field)';
    d.eval       = '[X Y Z 3] eigenvalues, descending. [0 0 0] where the fit was skipped (mask edge)';
    d.FA         = '[X Y Z] fractional anisotropy in [0,1]; the termination criterion';
    d.parcellation_mask      = '[X Y Z] region label per voxel (0 = none); used for ROI seeding / filtering';
    d.parcellation_mask_file = 'path of the label volume';
    d.atlas_labels           = 'label index -> name for parcellation_mask';
    d.parcellation_mask_jhu  = '[X Y Z] JHU atlas labels (kept when a bundle-ROI parcellation replaced the atlas)';
    d.atlas_labels_jhu       = 'label index -> name for parcellation_mask_jhu';
    d.atlas_type             = 'which atlas/parcellation source is active';
    d.roi_masks  = 'name -> [X Y Z] logical bundle masks (ISMRM scoring ROIs), read by nim_roi_mask';
    d.roi_source = 'directory the roi_masks came from';
    d.wm_mask    = '[X Y Z] white-matter probability/mask (ACT)';
    d.gm_mask    = '[X Y Z] grey-matter mask (ACT: a track ending here is KEPT)';
    d.csf_mask   = '[X Y Z] CSF mask (ACT: a track entering here is DISCARDED)';
    d.wm_mask_file = 'path'; d.gm_mask_file = 'path'; d.csf_mask_file = 'path';
    % direction field - built per run by nim_field (step 2)
    d.peaks      = '[X Y Z P 3] CSD FOD peak directions, unit, largest first; P = csd.max_peaks';
    d.npeaks     = '[X Y Z] how many of the P peaks are valid at each voxel';
    d.peak_w     = '[X Y Z P] FOD amplitude of each peak';
    d.fod_sh     = '[X Y Z C] FOD spherical-harmonic coefficients (csd.lmax)';
    d.mmf_e1_dwi    = '[X Y Z 3] frame e1 fitted directly to the DW signal (field: dwi)';
    d.mmf_kappa_dwi = '[X Y Z 3] curvature vector fitted with it (field: dwi)';
    d.mmf_dwi_resid = '[X Y Z] fit residual (field: dwi)';
    % geometry - built per run by nim_mmf_geometry (step 3), algorithm: mmf only
    d.mmf_frames   = '[X Y Z 3 3] moving frame; mmf_frames(x,y,z,:,k) = e_k, e1 = tangent';
    d.mmf_kappa    = '[X Y Z 3] connection curvature vector (Eq 6-9)';
    d.mmf_tau      = '[X Y Z] torsion';
    d.mmf_peakdirs = '[X Y Z P 3] per-peak tangents (csd: one pathway per peak)';
    d.mmf_kappa_p  = '[X Y Z P 3] per-peak curvature vectors';
    d.mmf_npeaks   = '[X Y Z] valid peaks per voxel (copy of npeaks)';
    d.mmf_multi    = 'true when the per-peak geometry is present';

    groups = struct('title', {}, 'fields', {});
    groups(end+1) = struct('title', 'dataset - from main.m, on disk, config-independent', ...
        'fields', {{'hdr','xdim','ydim','zdim','size3','img','bval','bvec','size_b0','size_bi', ...
                    'img_b0','img_bi','thrsh_b0','mask','DT','evec','eval','FA', ...
                    'parcellation_mask','parcellation_mask_file','atlas_labels', ...
                    'parcellation_mask_jhu','atlas_labels_jhu','atlas_type','roi_masks','roi_source', ...
                    'wm_mask','gm_mask','csf_mask','wm_mask_file','gm_mask_file','csf_mask_file'}});
    groups(end+1) = struct('title', 'direction field - nim_field, step 2, per run (field: csd | dwi)', ...
        'fields', {{'peaks','npeaks','peak_w','fod_sh','mmf_e1_dwi','mmf_kappa_dwi','mmf_dwi_resid'}});
    groups(end+1) = struct('title', 'geometry - nim_mmf_geometry, step 3, per run (algorithm: mmf only)', ...
        'fields', {{'mmf_frames','mmf_kappa','mmf_tau','mmf_peakdirs','mmf_kappa_p','mmf_npeaks','mmf_multi'}});
end

function [groups, d] = options_dictionary()
% Flat option name -> {config key, meaning}. The flat names are what
% nim_config_to_options emits and the trackers read; the config keys are the
% schema paths (nim_config_schema) so the log can be traced back to the YAML.
    d = struct();
    d.algorithm  = {'tractography.algorithm', 'which tracker runTractography dispatches to'};
    d.field      = {'tractography.field', 'direction source: dti (nim.evec) | csd (nim.peaks) | dwi (nim.mmf_e1_dwi)'};
    d.integrator        = {'tractography.integrator.method', 'stepping scheme name: euler | rk2 | rk4 | rkf45'};
    d.integration_order = {'tractography.integrator.method', 'same choice as a number: 1 | 2 | 4 | 5 (legacy selector, NOT an order claim)'};
    d.step_size    = {'tractography.integrator.step', 'step h in voxels (initial step for rkf45)'};
    d.adaptive_step = {'tractography.integrator.adaptive', 'rkf45 only: adapt h to the error estimate'};
    d.rkf_tolerance = {'tractography.integrator.tolerance', 'rkf45 only: local error tolerance in voxels'};
    d.rkf_tol       = {'tractography.integrator.tolerance', 'same value under the name mmf reads'};
    d.step_min      = {'tractography.integrator.step_min', 'rkf45 only: smallest h'};
    d.step_max      = {'tractography.integrator.step_max', 'rkf45 only: largest h'};
    d.rkf_safety    = {'tractography.integrator.safety', 'rkf45 only: safety factor on the step update'};
    d.interp_method = {'tractography.interpolation.method', 'kernel for the direction field: trilinear (C0) | cubic (C1) | spline (C2)'};
    d.upsample      = {'tractography.interpolation.upsample', 'sample the field on a grid of spacing 1/upsample before interpolating'};
    d.seed_density      = {'tractography.seeding.density', 'seeds per seeded voxel (nim_seed_offsets places them)'};
    d.seed_strategy     = {'tractography.seeding.strategy', 'uniform (deterministic lattice) | random (jittered)'};
    d.seed_fa_threshold = {'tractography.seeding.fa_min', 'voxels below this FA were excluded from seed_mask'};
    d.seed_roi          = {'tractography.seeding.roi', 'atlas regions the seeds were restricted to ({} = whole brain)'};
    d.seed_roi_dilate   = {'tractography.seeding.roi_dilate', 'dilation applied to those regions'};
    d.seed_mask     = {'(built in runTractography step 4)', '[X Y Z] logical: seed inside these voxels. THE seed input - the tracker generates seed points from it'};
    d.seed_roi_info = {'(built in runTractography step 4)', 'voxel counts of the ROI seed mask at each masking stage'};
    d.termination_fa = {'tractography.termination.fa_min', 'stop when interpolated FA drops below this'};
    d.angle_thresh   = {'tractography.termination.angle_max', 'max turn in DEGREES PER VOXEL OF ARC (min radius 57.3/angle_max voxels)'};
    d.min_length     = {'tractography.termination.min_arc', 'discard tracks whose chord length is below this (voxels)'};
    d.max_arc        = {'tractography.termination.max_arc', 'stop a half-track after this arc length (voxels)'};
    d.max_steps      = {'(derived) ceil(max_arc / step)', 'the same limit as a step count'};
    d.act_enabled = {'tractography.act', 'ACT requested. Trackers decide from whether the masks below are non-empty'};
    d.wm_mask  = {'(nim.wm_mask, step 4b)', '[X Y Z] or [] - white matter'};
    d.gm_mask  = {'(nim.gm_mask, step 4b)', '[X Y Z] or [] - grey matter: a half-track ending here is kept'};
    d.csf_mask = {'(nim.csf_mask, step 4b)', '[X Y Z] or [] - CSF: a half-track entering here is discarded'};
    d.include_roi       = {'tractography.filter.include_roi', 'step 6, after tracking: keep tracks touching these regions'};
    d.exclude_roi       = {'tractography.filter.exclude_roi', 'step 6: discard tracks touching these'};
    d.roi_filter_mode   = {'tractography.filter.mode', 'step 6: all | any of the include regions'};
    d.roi_filter_dilate = {'tractography.filter.roi_dilate', 'step 6: dilation of the filter masks'};
    d.endpoints_in = {'tractography.filter.endpoints_in', 'step 6: endpoint test, one end in each of two regions'};
    d.contained_in = {'tractography.filter.contained_in', 'step 6: every point inside these regions'};
    d.any_in       = {'tractography.filter.any_in', 'step 6: at least one point inside these regions'};
    d.length       = {'tractography.filter.length', 'step 6: [min max] track length in mm'};
    d.length_x = {'tractography.filter.length_x', 'step 6: [min max] net x displacement, mm'};
    d.length_y = {'tractography.filter.length_y', 'step 6: net y, mm'};
    d.length_z = {'tractography.filter.length_z', 'step 6: net z, mm'};
    d.length_x_abs = {'tractography.filter.length_x_abs', 'step 6: [min max] total x travel, mm'};
    d.length_y_abs = {'tractography.filter.length_y_abs', 'step 6: total y, mm'};
    d.length_z_abs = {'tractography.filter.length_z_abs', 'step 6: total z, mm'};
    d.csd_lmax         = {'tractography.csd.lmax', 'used by nim_field (step 2) to build nim.peaks; not read by trackers'};
    d.csd_max_peaks    = {'tractography.csd.max_peaks', 'P, the peak dimension of nim.peaks'};
    d.csd_peak_thresh  = {'tractography.csd.peak_thresh', 'relative amplitude below which a peak is dropped'};
    d.csd_peak_min_sep = {'tractography.csd.peak_min_sep', 'minimum angle between peaks (degrees)'};
    d.csd_n_iter       = {'tractography.csd.n_iter', 'CSD iterations'};
    d.mmf_anchor = {'tractography.mmf.anchor', 'mmf only: 0 = pure connection-form evolution (Eq 10-11); >0 re-anchors the carried frame to the field'};
    d.trace     = {'tractography.debug.trace', 'record a per-step trace (position, direction, FA, turn, stop reason) into meta.trace'};
    d.trace_max = {'tractography.debug.trace_max', 'how many seeds to trace, sampled evenly (0 = all)'};
    d.output_arc_step   = {'tractography.output.arc_step', 'step 7: resample saved tracks to this arc spacing (0 = keep every step)'};
    d.enable_diagnostics = {'tractography.diagnostics', 'extra per-run statistics in the log'};

    groups = struct('title', {}, 'fields', {});
    groups(end+1) = struct('title', 'identity', 'fields', {{'algorithm','field'}});
    groups(end+1) = struct('title', 'integrator', 'fields', {{'integrator','integration_order','step_size', ...
        'adaptive_step','rkf_tolerance','rkf_tol','step_min','step_max','rkf_safety'}});
    groups(end+1) = struct('title', 'interpolation', 'fields', {{'interp_method','upsample'}});
    groups(end+1) = struct('title', 'seeding', 'fields', {{'seed_mask','seed_density','seed_strategy', ...
        'seed_fa_threshold','seed_roi','seed_roi_dilate','seed_roi_info'}});
    groups(end+1) = struct('title', 'termination', 'fields', {{'termination_fa','angle_thresh','min_length','max_arc','max_steps'}});
    groups(end+1) = struct('title', 'ACT', 'fields', {{'act_enabled','wm_mask','gm_mask','csf_mask'}});
    groups(end+1) = struct('title', 'field construction (consumed in step 2, not by the tracker)', 'fields', ...
        {{'csd_lmax','csd_max_peaks','csd_peak_thresh','csd_peak_min_sep','csd_n_iter'}});
    groups(end+1) = struct('title', 'mmf', 'fields', {{'mmf_anchor'}});
    groups(end+1) = struct('title', 'post-processing (steps 6-7, not read by the tracker)', 'fields', ...
        {{'include_roi','exclude_roi','roi_filter_mode','roi_filter_dilate','endpoints_in','contained_in','any_in', ...
          'length','length_x','length_y','length_z','length_x_abs','length_y_abs','length_z_abs','output_arc_step'}});
    groups(end+1) = struct('title', 'debug', 'fields', {{'trace','trace_max','enable_diagnostics'}});
end
