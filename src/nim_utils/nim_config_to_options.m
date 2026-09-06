function options = nim_config_to_options(config)
% nim_config_to_options: Canonical nested config -> flat tracker option struct.
%
% This is the ONE place where the legacy flat option names live. Trackers keep
% reading `options.step_size`, `options.integration_order`, etc.; the config
% surface users see is the canonical nested schema. When the trackers are later
% updated to read canonical names, this function shrinks and then disappears.
%
%   options = nim_config_to_options(load_config_yaml('config/ismrm2015.yml'))

    if ~isstruct(config) || ~isfield(config, 'tractography')
        error('nim_config_to_options:badInput', 'Expected a config struct with a tractography section.');
    end
    t = config.tractography;

    options = struct();

    % --- identity -----------------------------------------------------------
    options.algorithm = t.algorithm;
    options.field     = t.field;

    % --- integrator ---------------------------------------------------------
    % integration_order is a METHOD selector in the trackers, despite the name.
    switch t.integrator.method
        case 'euler', options.integration_order = 1;
        case 'rk2',   options.integration_order = 2;
        case 'rk4',   options.integration_order = 4;
        case 'rkf45', options.integration_order = 5;
        otherwise
            error('nim_config_to_options:badMethod', 'Unknown integrator.method "%s".', t.integrator.method);
    end
    options.step_size = t.integrator.step;

    % Adaptivity applies only to rkf45. Fixed-step methods must report false,
    % which is exactly what the old configs carried.
    options.adaptive_step = strcmp(t.integrator.method, 'rkf45') && t.integrator.adaptive;

    options.rkf_tolerance = t.integrator.tolerance;
    options.rkf_tol       = t.integrator.tolerance;   % mmf spells it this way
    options.step_min      = t.integrator.step_min;
    options.step_max      = t.integrator.step_max;
    options.rkf_safety    = t.integrator.safety;
    options.integrator    = t.integrator.method;      % mmf reads a method name

    % --- interpolation ------------------------------------------------------
    options.interp_method = t.interpolation.method;
    options.upsample      = t.interpolation.upsample;
    % NOTE: no sel_power. It was a direction-steering term, removed from hinec -
    % tracking is interpolation + integration only (see nim_config_retired).

    % --- seeding ------------------------------------------------------------
    options.seed_density       = t.seeding.density;
    options.seed_strategy      = t.seeding.strategy;
    options.seed_fa_threshold  = t.seeding.fa_min;
    options.seed_roi           = t.seeding.roi;
    options.seed_roi_dilate    = t.seeding.roi_dilate;

    % --- termination --------------------------------------------------------
    options.termination_fa = t.termination.fa_min;
    options.angle_thresh   = t.termination.angle_max;
    options.min_length     = t.termination.min_arc;

    % An angle budget above 90 degrees cannot fire. Tangents are sign-aligned
    % because v1 is a line field, so the measured turn lives in [0, 90] and a
    % larger budget is not a loose criterion - it is an absent one. Say so, loudly:
    % a silently inert termination criterion looks exactly like a satisfied one,
    % and it is worst at the coarse end of a step-size ladder, where it switches
    % itself off and then back on as the ladder refines.
    if t.termination.angle_max > 0
        [budget_deg, ~] = nim_angle_limit(t.termination.angle_max, t.integrator.step);
        if budget_deg >= 90
            warning('nim_config:angleMaxInert', ...
                ['angle_max = %g deg/voxel at step = %g gives a %.1f deg budget, which is >= the ' ...
                 '90 deg ceiling: the angle criterion CANNOT fire and is effectively disabled. ' ...
                 'Use angle_max <= %.0f to make it active at this step, or angle_max: 0 to ' ...
                 'disable it deliberately.'], ...
                t.termination.angle_max, t.integrator.step, ...
                t.termination.angle_max * t.integrator.step, floor(90 / t.integrator.step));
        end
    end

    % max_steps is DERIVED from arc length, so refining the step can no longer
    % silently truncate tracks. This is the round-trip inverse of the legacy
    % migration (max_arc = max_steps * step), so existing configs are unchanged.
    options.max_steps = ceil(t.termination.max_arc / t.integrator.step);
    options.max_arc   = t.termination.max_arc;

    % NOTE: fa_threshold is deliberately NOT emitted. It was functionally dead
    % (see nim_config_retired) and every original config set it to 0.1, which is
    % also every tracker's internal default - so omitting it is behaviour-identical.

    % --- filtering ----------------------------------------------------------
    options.include_roi = t.filter.include_roi;
    options.exclude_roi = t.filter.exclude_roi;
    options.roi_filter_mode = t.filter.mode;
    options.roi_filter_dilate = t.filter.roi_dilate;
    options.endpoints_in      = t.filter.endpoints_in;
    options.contained_in      = t.filter.contained_in;

    % --- output -------------------------------------------------------------
    options.output_arc_step = t.output.arc_step;

    % --- act / diagnostics --------------------------------------------------
    options.act_enabled        = t.act;
    options.enable_diagnostics = t.diagnostics;

    % --- csd ----------------------------------------------------------------
    options.csd_lmax         = t.csd.lmax;
    options.csd_max_peaks    = t.csd.max_peaks;
    options.csd_peak_thresh  = t.csd.peak_thresh;
    options.csd_peak_min_sep = t.csd.peak_min_sep;
    options.csd_n_iter       = t.csd.n_iter;

    % --- mmf ----------------------------------------------------------------
    options.mmf_anchor      = t.mmf.anchor;
    options.frame_sel_power = t.mmf.frame_sel_power;
end
