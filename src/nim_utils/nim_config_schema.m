function S = nim_config_schema()
% nim_config_schema: Single source of truth for the HINEC configuration surface.
%
% Every configurable parameter is declared here exactly once: its canonical
% (nested) path, default, type, permitted values, legacy aliases, which
% tractography algorithms actually read it, and its one-line description.
%
% Everything else is derived from this:
%   - load_config_yaml  : defaults, validation, unknown-key rejection, migration
%   - nim_config_to_options : flattening to the option names the trackers read
%   - the config reference doc : generated, so docs cannot drift from code
%
% Entry fields:
%   path    - canonical dotted path, e.g. 'tractography.integrator.method'
%   default - default value (used when the key is absent)
%   type    - 'numeric' | 'string' | 'logical' | 'list'
%   allowed - cell array of permitted values for enums, {} otherwise
%   range   - [lo hi] inclusive bounds for numerics, [] otherwise
%   legacy  - cell array of superseded FLAT key names that map here
%   algos   - algorithms that read it: any of 'hinec','standard','mmf', or
%             'all'; 'n/a' marks preprocessing keys
%   help    - one-line description (used to generate documentation)
%
% Nesting is deliberately limited to two levels below a section
% (section.group.key). A third level is a parse error - see load_config_yaml.

% NOTE: cell-valued fields (default, allowed, legacy, algos) must be wrapped in
% braces when handed to struct(), otherwise struct() expands them into a struct
% ARRAY (or, for {}, an empty struct) instead of storing the cell as one value.
e = @(path, default, type, allowed, range, legacy, algos, help) struct( ...
    'path', path, 'default', {default}, 'type', type, 'allowed', {allowed}, ...
    'range', range, 'legacy', {legacy}, 'algos', {algos}, 'help', help);

S = [
% ---------------------------------------------------------------- preprocessing
e('preprocessing.run_denoising',        true,          'logical', {}, [],       {}, {'n/a'}, 'Run denoising before tensor fitting.')
e('preprocessing.denoise_method',       'dwidenoise',  'string',  {'dwidenoise','nlmeans','gaussian'}, [], {}, {'n/a'}, 'Denoising backend. Falls back to gaussian when MRtrix3 is absent.')
e('preprocessing.run_motion_correction',true,          'logical', {}, [],       {}, {'n/a'}, 'Run motion correction (FSL).')
e('preprocessing.run_eddy',             true,          'logical', {}, [],       {}, {'n/a'}, 'Run eddy-current correction (FSL).')
e('preprocessing.improve_mask',         true,          'logical', {}, [],       {}, {'n/a'}, 'Refine the brain mask using FA.')
e('preprocessing.atlas_type',           'jhu',         'string',  {'jhu','harvardoxford','jhu-tract'}, [], {}, {'n/a'}, 'Atlas used for parcellation.')
e('preprocessing.t1_available',         false,         'logical', {}, [],       {}, {'n/a'}, 'A T1 volume is present alongside the DWI.')
e('preprocessing.use_t1_registration',  false,         'logical', {}, [],       {}, {'n/a'}, 'Register T1 to DWI space.')
e('preprocessing.register_to_mni',      false,         'logical', {}, [],       {}, {'n/a'}, 'Register to MNI space.')

% ------------------------------------------------------------------ tractography
e('tractography.algorithm',             'hinec',       'string',  {'hinec','standard','mmf'}, [], {}, {'all'}, 'Tracking algorithm.')
e('tractography.field',                 'dti',         'string',  {'dti','csd'}, [],  {}, {'hinec','mmf'}, 'Direction source: DTI principal eigenvector, or CSD FOD peaks.')
e('tractography.act',                   false,         'logical', {}, [],       {'act_enabled'}, {'hinec'}, 'Anatomically constrained tracking using WM/GM/CSF masks.')

% --- integrator -------------------------------------------------------------
e('tractography.integrator.method',     'rk4',         'string',  {'euler','rk2','rk4','rkf45'}, [], {'integrator'}, {'hinec','mmf'}, 'Numerical stepping scheme. NOTE: this is a method name, not an order - RKF45 is a 4(5) embedded pair.')
e('tractography.integrator.step',       0.2,           'numeric', {}, [1e-6 10], {'step_size'}, {'all'}, 'Integration step in voxels. Fixed step, or the initial step for rkf45.')
e('tractography.integrator.tolerance',  0.01,          'numeric', {}, [1e-12 1], {'rkf_tolerance','rkf_tol'}, {'hinec','mmf'}, 'Adaptive error tolerance in voxels. rkf45 ONLY.')
e('tractography.integrator.step_min',   0.01,          'numeric', {}, [1e-9 10], {'step_min'}, {'hinec','mmf'}, 'Minimum adaptive step. rkf45 ONLY.')
e('tractography.integrator.step_max',   1.0,           'numeric', {}, [1e-6 10], {'step_max'}, {'hinec','mmf'}, 'Maximum adaptive step. rkf45 ONLY.')
e('tractography.integrator.safety',     0.9,           'numeric', {}, [0.01 1],  {'rkf_safety'}, {'hinec'}, 'Adaptive step safety factor. rkf45 ONLY.')
e('tractography.integrator.adaptive',   true,          'logical', {}, [],        {'adaptive_step'}, {'hinec'}, 'Adaptive step-size control. rkf45 ONLY. rkf45 with adaptive:false is fixed-step RKF45 (a real third mode, kept for completeness).')

% --- interpolation ----------------------------------------------------------
e('tractography.interpolation.method',  'trilinear',   'string',  {'trilinear','cubic','spline'}, [], {'interp_method'}, {'hinec','mmf'}, 'Spatial interpolation kernel for the direction field. These differ in SMOOTHNESS, which caps how much of an integrator formal order is reachable: a Runge-Kutta method of order p needs a right-hand side with p continuous derivatives. trilinear is C0 (kinked at every voxel face), cubic is C1 (Keys cubic convolution, NOT a spline - its second derivative jumps), spline is C2 (a genuine cubic spline). Measured here: RK4 reaches 2.02 on trilinear and 2.02 on cubic, identical, which is what a cap below C2 predicts.')
e('tractography.interpolation.upsample',  1,           'numeric', {}, [0.05 8],  {}, {'hinec','mmf'}, 'Spatial sampling factor for the direction field: the field is sampled on a grid of spacing 1/upsample voxels before the interpolants are built. 1 = the acquisition grid. Above 1 refines toward the continuous field the samples imply; BELOW 1 coarsens, discarding spatial information, which is how the space axis of a convergence study is swept. The coordinate frame is unchanged, so positions, step sizes and lengths stay in native voxel units and runs at different factors are directly comparable. Note the u -> infinity limit is the native-resolution interpolant, not ground-truth anatomy.')

% --- seeding ----------------------------------------------------------------
e('tractography.seeding.density',       8,             'numeric', {}, [1 1000],  {'seed_density'}, {'all'}, 'Seeds per seeded voxel - honoured exactly. Placed on a deterministic sub-voxel lattice (perfect cubes) or a spread subset of one (other counts), so runs are reproducible.')
e('tractography.seeding.strategy',      'uniform',     'string',  {'uniform','random'}, [], {'seed_strategy'}, {'all'}, 'Seed placement. uniform = deterministic lattice (reproducible); random = jittered.')
e('tractography.seeding.fa_min',        0.05,          'numeric', {}, [0 1],     {'seed_fa_threshold'}, {'all'}, 'Minimum FA for a voxel to be seeded. Default 0.05 excludes CSF only.')
e('tractography.seeding.roi',           {},            'list',    {}, [],        {}, {'all'}, 'Restrict seeding to these atlas regions. JHU indices and/or names, freely mixed. Empty = whole brain mask.')
e('tractography.seeding.roi_dilate',    0,             'numeric', {}, [0 10],    {}, {'all'}, 'Dilate the seed ROI by this many voxels before seeding.')

% --- termination ------------------------------------------------------------
e('tractography.termination.fa_min',    0.10,          'numeric', {}, [0 1],     {'termination_fa'}, {'all'}, 'Stop tracking below this FA. (NOT the legacy fa_threshold - see nim_config_retired.)')
e('tractography.termination.angle_max', 225,           'numeric', {}, [0 3600],  {}, {'all'}, 'Maximum turn in DEGREES PER VOXEL OF ARC, i.e. a minimum radius of curvature R = 57.3/angle_max voxels. The default 225 is the classic 45 deg/step evaluated at the default step of 0.2. Step-invariant: the budget for one step is angle_max x the NOMINAL step arc (never the realised chord - that would make the budget method-dependent), so refining the step does not loosen the constraint. CEILING: tangents are sign-aligned because v1 is a line field, so a measured turn never exceeds 90 deg - any angle_max above 90/step is INERT, not merely loose, and 225 goes inert for any step >= 0.4. 0 disables the criterion outright, which is what a control run should use.')
e('tractography.termination.max_arc',   400,           'numeric', {}, [1 100000],{}, {'all'}, 'Maximum track arc length in voxels. Step-invariant: max_steps is derived as ceil(max_arc/step).')
e('tractography.termination.min_arc',   15,            'numeric', {}, [0 100000],{'min_length'}, {'all'}, 'Discard tracks shorter than this arc length in voxels.')

% --- filter -----------------------------------------------------------------
e('tractography.filter.include_roi',    {},            'list',    {}, [],        {}, {'all'}, 'Keep only tracks touching these regions (waypoint).')
e('tractography.filter.exclude_roi',    {},            'list',    {}, [],        {}, {'all'}, 'Discard tracks touching these regions.')
e('tractography.filter.mode',           'all',         'string',  {'all','any'}, [], {}, {'all'}, 'Whether a track must touch ALL include regions or ANY of them.')
e('tractography.filter.roi_dilate',     0,             'numeric', {}, [0 10],    {}, {'all'}, 'Dilate the include/exclude masks by this many voxels before testing.')
e('tractography.filter.endpoints_in',   {},            'list',    {}, [],        {}, {'all'}, 'Two regions; keep a track only if one END lands in the first and the other END in the second, either way round. This is an ENDPOINT test, not a waypoint test - include_roi asks whether a track passes through a region, this asks where it stops. It is half of how the ISMRM 2015 scorer defines a bundle (head + tail).')
e('tractography.filter.contained_in',   {},            'list',    {}, [],        {}, {'all'}, 'Keep a track only if EVERY point lies inside these regions. The other half of the ISMRM bundle definition (all_mask): a streamline that wanders outside the corridor is not that bundle, however it ends.')

% --- output -----------------------------------------------------------------
e('tractography.output.arc_step',       0,             'numeric', {}, [0 100],   {}, {'all'}, 'Resample saved streamlines to this arc-length spacing in voxels. 0 = store every integration step. Decouples file size from step size; integration accuracy is unaffected.')

% --- csd (field: csd only) --------------------------------------------------
e('tractography.csd.lmax',              6,             'numeric', {}, [2 16],    {'csd_lmax'}, {'hinec','mmf'}, 'Spherical harmonic order for CSD.')
e('tractography.csd.max_peaks',         3,             'numeric', {}, [1 10],    {'csd_max_peaks'}, {'hinec','mmf'}, 'Maximum FOD peaks retained per voxel.')
e('tractography.csd.peak_thresh',       0.5,           'numeric', {}, [0 1],     {'csd_peak_thresh'}, {'hinec','mmf'}, 'Relative amplitude threshold for accepting an FOD peak.')
e('tractography.csd.peak_min_sep',      45,            'numeric', {}, [0 180],   {'csd_peak_min_sep'}, {'hinec','mmf'}, 'Minimum angular separation between FOD peaks (degrees).')
e('tractography.csd.n_iter',            50,            'numeric', {}, [1 1000],  {}, {'hinec','mmf'}, 'CSD deconvolution iterations.')

% --- mmf (algorithm: mmf only) ----------------------------------------------
e('tractography.mmf.anchor',            0,             'numeric', {}, [0 1],     {'mmf_anchor'}, {'mmf'}, 'Re-anchor strength of e1 toward the field tangent. 0 = pure connection-form evolution.')
e('tractography.mmf.frame_sel_power',   16,            'numeric', {}, [0 256],   {'frame_sel_power'}, {'mmf'}, 'Directional selectivity used when building the moving frame.')

% --- diagnostics ------------------------------------------------------------
e('tractography.diagnostics',           true,          'logical', {}, [],        {'enable_diagnostics'}, {'all'}, 'Write per-run diagnostic reports.')
];

end
