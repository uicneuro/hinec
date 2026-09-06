function runTractography(data_path, varargin)
% runTractography: Entry point for DTI tractography with algorithm selection
%
% Usage:
%   % Standard FACT tractography (default)
%   runTractography('sample_parcellated.mat')
%   runTractography('sample_parcellated.mat', 'standard')
%
%   % HINEC high-order tractography (interpolation + RK4 + ACT)
%   runTractography('sample_parcellated.mat', 'hinec')
%
%   % IronTract Challenge mode (generates tracks + submission files)
%   runTractography('data.mat', 'IronTract', 'injection.nii.gz', 'submissions/')
%
% Arguments:
%   data_path - Path to .mat file with nim structure
%   algorithm - 'standard' (FACT) or 'hinec' (high-order) [default: 'standard']
%   'IronTract' - Optional flag to enable IronTract Challenge submission
%   injection_file - Path to injection site mask (required if IronTract enabled)
%   output_dir - Output directory for submissions (required if IronTract enabled)
%
% Examples:
%   % Standard FACT tractography
%   runTractography('processed.mat', 'standard');
%
%   % HINEC high-order tractography
%   runTractography('processed.mat', 'hinec');
%
%   % IronTract Challenge with standard algorithm
%   main('ironTract/sub-MR243', 'ironTract.mat');
%   runTractography('ironTract.mat', 'IronTract', 'ironTract/injection.nii.gz', 'ironTract_submissions/');

if nargin < 1
    data_path = 'sample_parcellated.mat';
end

% Parse arguments
algorithm = 'standard';  % Default to standard FACT
enable_irontract = false;
injection_file = '';
irontract_output_dir = '';
config = struct();
use_config = false;
run_info = struct();

if nargin >= 2
    arg1 = varargin{1};

    % Check if first argument is YAML config structure
    if isstruct(arg1) && isfield(arg1, 'tractography')
        config = arg1;
        use_config = true;
        fprintf('Using YAML configuration for tractography\n');

        % Extract algorithm from config
        if isfield(config.tractography, 'algorithm')
            algorithm = lower(config.tractography.algorithm);
        end

        % Check for run_info in second varargin
        if nargin >= 3 && isstruct(varargin{2}) && isfield(varargin{2}, 'run_dir')
            run_info = varargin{2};
            fprintf('Using run directory: %s\n', run_info.run_dir);
        end
    elseif isstruct(arg1) && isfield(arg1, 'run_dir')
        % First argument is run_info structure
        run_info = arg1;
        fprintf('Using run directory: %s\n', run_info.run_dir);
    elseif strcmpi(arg1, 'standard') || strcmpi(arg1, 'hinec') || strcmpi(arg1, 'mmf')
        % Legacy: algorithm selection string
        algorithm = lower(arg1);
    elseif strcmpi(arg1, 'IronTract')
        % IronTract mode
        enable_irontract = true;
        if nargin < 4
            error('IronTract mode requires: runTractography(data_path, ''IronTract'', injection_file, output_dir)');
        end
        injection_file = varargin{2};
        irontract_output_dir = varargin{3};
    else
        error('Invalid argument. Use YAML config, ''standard'', ''hinec'', or ''IronTract''.');
    end
end

% Check if using run directory organization
use_run_dir = ~isempty(fieldnames(run_info));

% Add necessary paths
addpath('src/nim_tractography');
addpath('src/nim_utils');
addpath('src/nim_plots');
addpath('src/nim_challenges');
addpath('src/nim_calculation');

fprintf('=== HINEC Tractography Pipeline ===\n');

%% Load data
fprintf('Loading data from %s...\n', data_path);
if ~exist(data_path, 'file')
    error('Data file not found: %s', data_path);
end
load(data_path, 'nim');

%% Check required fields
if ~isfield(nim, 'evec')
    error('Eigenvectors not found. Please run main() first to generate DTI data.');
end
if ~isfield(nim, 'FA')
    error('FA not found. Please run main() first to generate DTI data.');
end

%% Set tractography parameters
if use_config
    fprintf('Loading tractography parameters from YAML config...\n');
    % The config is the canonical nested schema; nim_config_to_options is the
    % single place the legacy flat option names are produced for the trackers.
    options = nim_config_to_options(config);

    % Display key parameters
    fprintf('  Algorithm: %s\n', algorithm);
    fprintf('  Integrator: %s\n', config.tractography.integrator.method);
    fprintf('  Step size: %.4g voxels\n', options.step_size);
    fprintf('  Max arc length: %g voxels (max_steps = %d, derived)\n', ...
        options.max_arc, options.max_steps);
    fprintf('  Seed density: %d seeds/voxel\n', options.seed_density);
    fprintf('  Termination FA: %.2f\n', options.termination_fa);
    fprintf('  Angle threshold: %.1f deg\n', options.angle_thresh);

    if options.adaptive_step
        fprintf('  RKF45 adaptive: enabled (tol=%.4g)\n', options.rkf_tolerance);
    end
else
    fprintf('Setting up default tractography parameters...\n');
    % Use hardcoded defaults (legacy behavior)
    options = struct();
    options.seed_density = 4;
    options.step_size = 0.5;
    options.fa_threshold = 0.15;
    options.termination_fa = 0.15;
    options.angle_thresh = 35;
    options.max_steps = 1000;
    options.min_length = 35;
    options.order = 1;
    options.interp_method = 'none';
end

%% CENTRALIZED SEEDING STRATEGY
% All seeding decisions happen here - nim_tractography_standard.m only executes
fprintf('\n=== Configuring Seeding Strategy ===\n');

% Strategy priority (best to worst):
% 1. Preprocessed brain mask (nim.mask) - most accurate
% 2. Expanded parcellation mask - good for labeled regions
% 3. FA-threshold fallback - basic but works everywhere

seed_mask = [];
seeding_strategy = '';

% Seed FA threshold: voxels below this FA are excluded from seeding.
% Default 0.05 (just excludes CSF). Set seed_fa_threshold: 0.2 in the config
% to seed only from white matter (matched WM-seeding comparisons).
if isfield(options, 'seed_fa_threshold') && ~isempty(options.seed_fa_threshold)
    seed_fa_thr = options.seed_fa_threshold;
else
    seed_fa_thr = 0.05;
end

% Strategy 0: explicit ROI (highest priority). When seeding.roi is set, seeds go
% only in those atlas regions - the cheapest way to interrogate one bundle, and
% the basis of the convergence ladder's fixed seed set.
if isfield(options, 'seed_roi') && ~isempty(options.seed_roi)
    fprintf('Strategy: explicit ROI seeding\n');
    seeding_strategy = 'roi';

    dil = 0;
    if isfield(options, 'seed_roi_dilate') && ~isempty(options.seed_roi_dilate)
        dil = options.seed_roi_dilate;
    end
    [roi_mask, roi_info] = nim_roi_mask(nim, options.seed_roi, dil);

    % Constrain to brain and FA exactly as the other strategies do, so switching
    % to ROI seeding changes WHERE seeds go and nothing else.
    if isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
        roi_mask = roi_mask & (nim.mask > 0.5);
        fprintf('  Intersected with brain mask: %d voxels\n', sum(roi_mask(:)));
    end
    n_before_fa = sum(roi_mask(:));
    roi_mask = roi_mask & (nim.FA > seed_fa_thr);
    fprintf('  FA > %.2f: %d -> %d voxels\n', seed_fa_thr, n_before_fa, sum(roi_mask(:)));

    if ~any(roi_mask(:))
        error('runTractography:emptyRoiSeedMask', ...
            ['ROI seed mask is empty after masking. Voxels: %d raw -> %d dilated -> ' ...
             '%d in brain -> 0 after FA > %.2f. Try a larger seeding.roi_dilate, a ' ...
             'lower seeding.fa_min, or check the region names.'], ...
            roi_info.raw_voxels, roi_info.dilated_voxels, n_before_fa, seed_fa_thr);
    end

    seed_mask = roi_mask;
    options.seed_roi_info = roi_info;

% Strategy 1: Preprocessed brain mask (BEST)
elseif isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
    brain_mask = nim.mask > 0.5;
    fprintf('Strategy: Preprocessed brain mask\n');
    seeding_strategy = 'brain_mask';

    % Refine with FA threshold to exclude CSF (and optionally GM if WM-seeding)
    brain_mask = brain_mask & (nim.FA > seed_fa_thr);
    fprintf('  Refined with FA > %.2f for seeding\n', seed_fa_thr);

    seed_mask = brain_mask;

% Strategy 2: Expanded parcellation mask
elseif isfield(nim, 'parcellation_mask') && any(nim.parcellation_mask(:) > 0)
    fprintf('Strategy: Expanded parcellation mask\n');
    parcel_mask = nim.parcellation_mask > 0;

    % Aggressive expansion to capture surrounding white matter
    se = strel('sphere', 3);  % 3-voxel dilation for better coverage
    brain_mask = imdilate(parcel_mask, se);

    % Constrain to reasonable FA to avoid CSF (very low threshold)
    brain_mask = brain_mask & (nim.FA > 0.05);
    fprintf('  Dilated parcellation by 3 voxels\n');
    fprintf('  Constrained with FA > 0.05\n');

    seeding_strategy = 'expanded_parcellation';
    seed_mask = brain_mask;

% Strategy 3: FA-based fallback (CONSERVATIVE)
else
    fprintf('⚠ Strategy: FA-threshold fallback (no brain mask available)\n');
    % Use low FA threshold for better coverage
    brain_mask = nim.FA > 0.10;  % Lower threshold = more coverage
    fprintf('  Using FA > 0.10 for brain boundary\n');
    fprintf('  WARNING: This misses low-anisotropy regions (fornix, cingulum)\n');

    seeding_strategy = 'fa_threshold';
    seed_mask = brain_mask;
end

% Quality check
seed_voxel_count = sum(seed_mask(:));
if seed_voxel_count == 0
    error('Seed mask is empty! Cannot proceed with tractography.');
end

fprintf('\n✓ Seeding strategy: %s\n', seeding_strategy);
fprintf('✓ Seed voxels: %d (%.1f%% of volume)\n', ...
    seed_voxel_count, 100 * seed_voxel_count / numel(nim.FA));
options.seed_mask = seed_mask;

%% ACT Configuration: Add tissue masks if available
fprintf('\n=== Anatomically Constrained Tractography (ACT) Configuration ===\n');
act_enabled = ~isfield(options,'act_enabled') || ~isequal(options.act_enabled, 0) && ~isequal(options.act_enabled, false);
if act_enabled && isfield(nim, 'wm_mask') && isfield(nim, 'gm_mask') && isfield(nim, 'csf_mask')
    % All three tissue masks are available - enable full ACT
    options.wm_mask = nim.wm_mask;
    options.gm_mask = nim.gm_mask;
    options.csf_mask = nim.csf_mask;

    wm_voxels = sum(nim.wm_mask(:) > 0.5);
    gm_voxels = sum(nim.gm_mask(:) > 0.5);
    csf_voxels = sum(nim.csf_mask(:) > 0.5);

    fprintf('✓ ACT enabled with tissue segmentation masks:\n');
    fprintf('  WM mask:  %d voxels (%.1f%% of volume)\n', ...
            wm_voxels, 100 * wm_voxels / numel(nim.FA));
    fprintf('  GM mask:  %d voxels (%.1f%% of volume)\n', ...
            gm_voxels, 100 * gm_voxels / numel(nim.FA));
    fprintf('  CSF mask: %d voxels (%.1f%% of volume)\n', ...
            csf_voxels, 100 * csf_voxels / numel(nim.FA));
    fprintf('  ACT will constrain seeding to WM and terminate at GM/CSF boundaries\n');
else
    % ACT off. Hand the tracker empty masks either way - it decides whether ACT
    % is active purely from whether it received them - but say WHICH of the two
    % reasons applies, because the remedies are different and telling someone to
    % regenerate masks they already have sends them in a circle.
    options.wm_mask = [];
    options.gm_mask = [];
    options.csf_mask = [];

    have_masks = isfield(nim, 'wm_mask') && isfield(nim, 'gm_mask') && isfield(nim, 'csf_mask');
    if have_masks
        fprintf('⚠ ACT disabled by config (tractography.act: false)\n');
        fprintf('  Tissue masks ARE present in the nim (WM %d, GM %d, CSF %d voxels)\n', ...
            sum(nim.wm_mask(:) > 0.5), sum(nim.gm_mask(:) > 0.5), sum(nim.csf_mask(:) > 0.5));
        fprintf('  Tractography will use FA-based termination only\n');
        fprintf('  To enable ACT, set  act: true  under tractography: in the config\n');
    else
        fprintf('⚠ ACT disabled: tissue masks not present in the nim\n');
        fprintf('  Tractography will use FA-based termination only\n');
        fprintf('  To enable ACT, run main() to generate tissue masks first\n');
    end
end

% Show detailed seeding statistics
total_seed_locations = sum(seed_mask(:));
estimated_seeds = total_seed_locations * options.seed_density;
fprintf('\n=== Seeding Statistics ===\n');
fprintf('Seed voxels: %d (%.1f%% of volume)\n', ...
    total_seed_locations, 100 * total_seed_locations / numel(nim.FA));
fprintf('Seeds per voxel: %d (density parameter)\n', options.seed_density);
fprintf('Estimated total seeds: %d\n', estimated_seeds);
fprintf('Expected tracks: ~%d (bidirectional from each seed)\n', estimated_seeds * 2);
fprintf('==========================\n');

%% Run tractography (algorithm-dependent)
% Generate timestamp for output filename
timestamp = datestr(now, 'yyyy-mm-dd_HH_MM_SS');


% CSD FOD peaks are needed by ANY tracker running field=csd (hinec AND mmf), so
% provision them BEFORE the algorithm dispatch. Compute with nim_csd when the config
% sets field=csd, cached next to the source nim (<source>_csd.mat) so it is computed
% once per preprocessed dataset and reused by every tractography config.
fld = 'dti';
if isfield(options, 'field') && ~isempty(options.field)
    fld = lower(char(string(options.field)));
end
if strcmp(fld, 'csd') && ~isfield(nim, 'peaks')
    csd_cache = regexprep(data_path, '\.mat$', '_csd.mat');
    if isfile(csd_cache)
        fprintf('field=csd: loading cached CSD FOD from %s\n', csd_cache);
        Sc = load(csd_cache);
        nim.peaks = Sc.peaks; nim.npeaks = Sc.npeaks; nim.peak_w = Sc.peak_w;
        if isfield(Sc, 'fod_sh'), nim.fod_sh = Sc.fod_sh; end
    else
        fprintf('field=csd: computing CSD FOD peaks (nim_csd)...\n');
        csd_opts = struct('lmax', 6, 'n_iter', 50, 'peak_thresh', 0.5, ...
                          'peak_min_sep', 45, 'max_peaks', 3);
        csd_keys = {'lmax', 'n_iter', 'peak_thresh', 'peak_min_sep', 'max_peaks'};
        for ci = 1:numel(csd_keys)
            ck = ['csd_' csd_keys{ci}];
            if isfield(options, ck) && ~isempty(options.(ck))
                csd_opts.(csd_keys{ci}) = options.(ck);
            end
        end
        nim = nim_csd(nim, csd_opts);
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

track_meta = struct();   % populated by hinec; empty for other trackers
if strcmpi(algorithm, 'mmf')
    fprintf('Running MMF connection-frame tractography (moving frames + connection 1-form)...\n');
    tic;
    tracks = nim_tractography_mmf_connframe(nim, options);
    elapsed_time = toc;
    output_filename = sprintf('tracks_mmf_%s.mat', timestamp);
elseif strcmpi(algorithm, 'hinec')
    fprintf('Running HINEC high-order tractography (interpolation + RK4 + ACT)...\n');
    tic;
    [tracks, track_meta] = nim_tractography_hinec(nim, options);
    elapsed_time = toc;
    output_filename = sprintf('tracks_hinec_%s.mat', timestamp);
else
    fprintf('Running standard FACT tractography...\n');
    tic;
    tracks = nim_tractography_standard(nim, options);
    elapsed_time = toc;
    output_filename = sprintf('tracks_standard_%s.mat', timestamp);
end

fprintf('Tractography completed in %.1f seconds\n', elapsed_time);
fprintf('Generated %d tracks\n', length(tracks));

if isempty(tracks)
    error('No tracks generated! Check FA threshold and seed mask.');
end

%% ROI filtering (include / exclude waypoints), applied before saving.
% No-op unless filter.include_roi or filter.exclude_roi is set.
[tracks, roi_filter_stats] = nim_filter_tracks_roi(tracks, nim, options);
% Keep per-track metadata aligned with the surviving tracks, otherwise
% track_meta.seed_index would point at the pre-filter ordering.
if roi_filter_stats.applied && isfield(track_meta, 'seed_index')
    track_meta.seed_index  = track_meta.seed_index(roi_filter_stats.keep);
    track_meta.seed_points = track_meta.seed_points(roi_filter_stats.keep, :);
end
if roi_filter_stats.applied && isempty(tracks)
    error('runTractography:emptyAfterRoiFilter', ...
        ['ROI filtering removed all %d tracks. Loosen filter.include_roi, raise ' ...
         'filter.roi_filter_dilate, or check the region names.'], roi_filter_stats.n_in);
end

%% Output decimation. Applied AFTER ROI filtering, because filtering tests which
% voxels a track visits and needs the full-resolution polyline to do that; a
% decimated track could skip a voxel it actually passed through.
if isfield(options, 'output_arc_step') && ~isempty(options.output_arc_step) ...
        && options.output_arc_step > 0
    pts_before = sum(cellfun(@(x) size(x, 1), tracks));
    tracks = nim_resample_track_arc(tracks, options.output_arc_step);
    pts_after = sum(cellfun(@(x) size(x, 1), tracks));
    fprintf('\nOutput decimation: arc_step %.4g voxels, %d -> %d points (%.1f%% of original)\n', ...
        options.output_arc_step, pts_before, pts_after, 100 * pts_after / max(pts_before, 1));
end

%% Compute statistics
track_lengths = cellfun(@(x) size(x, 1), tracks);
fprintf('\nTrack Statistics:\n');
fprintf('  Mean length: %.1f points\n', mean(track_lengths));
fprintf('  Max length: %d points\n', max(track_lengths));
fprintf('  Min length: %d points\n', min(track_lengths));

%% Save results
if use_run_dir
    output_dir = run_info.tractography_dir;
    fprintf('\nUsing run directory for tractography output\n');
else
    output_dir = 'tractography_results';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
end

save(fullfile(output_dir, output_filename), 'tracks', 'options', 'elapsed_time', 'algorithm', 'track_meta', '-v7.3');
fprintf('\nResults saved to %s/%s\n', output_dir, output_filename);
fprintf('Algorithm used: %s\n', algorithm);

% Save diagnostic information if using run directory
if use_run_dir
    diagnostics_file = fullfile(run_info.diagnostics_dir, 'track_statistics.txt');
    fid = fopen(diagnostics_file, 'w');
    fprintf(fid, 'HINEC Tractography Statistics\n');
    fprintf(fid, '==============================\n\n');
    fprintf(fid, 'Run ID: %s\n', run_info.run_id);
    fprintf(fid, 'Timestamp: %s\n', timestamp);
    fprintf(fid, 'Algorithm: %s\n\n', algorithm);
    fprintf(fid, 'Track Statistics:\n');
    fprintf(fid, '  Total tracks: %d\n', length(tracks));
    fprintf(fid, '  Mean length: %.1f points\n', mean(track_lengths));
    fprintf(fid, '  Max length: %d points\n', max(track_lengths));
    fprintf(fid, '  Min length: %d points\n', min(track_lengths));
    fprintf(fid, '  Computation time: %.1f seconds\n', elapsed_time);
    fprintf(fid, '\nSeeding Statistics:\n');
    fprintf(fid, '  Seed voxels: %d\n', total_seed_locations);
    fprintf(fid, '  Seeds per voxel: %d\n', options.seed_density);
    fprintf(fid, '  Estimated total seeds: %d\n', estimated_seeds);
    fclose(fid);
    fprintf('Diagnostics saved to %s\n', diagnostics_file);
end

%% IronTract Challenge Submission (if enabled)
if enable_irontract
    fprintf('\n=== IronTract Challenge Submission ===\n');
    fprintf('Generating submission files...\n');

    % Set up IronTract options
    irontract_opts = struct();
    irontract_opts.angle_thresholds = [30, 45, 60, 75, 90];
    % DON'T pass tracks_file - force regeneration with injection site seeds
    irontract_opts.tracks_file = '';  % Empty = regenerate with injection-site seeding
    irontract_opts.base_options = options;

    % NOTE: IronTract has special seeding requirements
    % nim_irontract_submit.m will MODIFY the seed_mask to include injection site
    % This is handled internally by adding: seed_mask | injection_mask
    % See nim_challenges/nim_irontract_submit.m lines 102-128 for details

    % Generate submission files
    nim_irontract_submit(data_path, injection_file, irontract_output_dir, irontract_opts);

    fprintf('IronTract submissions ready in %s/\n', irontract_output_dir);
end

fprintf('=== Tractography Complete ===\n');
end

function val = get_opt(s, field, default)
% Safe struct field accessor with a default.
if isfield(s, field) && ~isempty(s.(field))
    val = s.(field);
else
    val = default;
end
end