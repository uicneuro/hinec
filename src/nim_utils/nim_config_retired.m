function R = nim_config_retired()
% nim_config_retired: Keys that were previously accepted and now do nothing.
%
% These are NOT migrated to a new path - they had no functional effect, so
% silently mapping them somewhere would change behaviour. The loader warns and
% drops them, naming the reason, so a user who had one in a config learns that
% it never did what its name suggested.
%
% Each entry: key, reason.

r = @(key, reason) struct('key', key, 'reason', reason);

R = [
% --- dead schema keys: defaulted AND validated in the old loader, read by no
%     tracker. Verified by grep across src/nim_tractography/ (0 references).
r('gate_power',         'Never read by any tracker. Defaulted and validated by the old loader but had no effect.')
r('crossing_cp',        'Never read by any tracker.')
r('curv_beta',          'Never read by any tracker.')
r('crossing_detect',    'Never read by any tracker.')
r('swing_ratio_max',    'Never read by any tracker.')
r('transport_gate',     'Never read by any tracker.')
r('transport_strength', 'Never read by any tracker.')
r('bishop_eps',         'Referenced only in a comment; no tracker reads it as an option.')

% --- fa_threshold: the subtle one. Verified per tracker:
%       hinec    (:188) - fprintf only
%       standard (:120) - fprintf only
%       mmf      (:86)  - seed-mask fallback, guarded by `if isempty(seed_mask)`,
%                         which never fires because runTractography.m:215 always
%                         supplies a seed mask
%       highorder(:77)  - same fallback pattern
%     So it is functionally dead on every path reachable from runTractography.
%     Termination is governed by termination_fa -> tractography.termination.fa_min;
%     seeding by seed_fa_threshold -> tractography.seeding.fa_min.
r('fa_threshold',       ['Functionally dead: printed by hinec/standard, and only a ' ...
                         'seed-mask fallback in mmf/highorder that never fires because ' ...
                         'runTractography always supplies a seed mask. Use ' ...
                         'termination.fa_min to stop tracking, or seeding.fa_min to ' ...
                         'restrict seeding. NOTE: in 5 of the 12 original configs ' ...
                         'fa_threshold and termination_fa held DIFFERENT values, so ' ...
                         'they were never interchangeable.'])

r('order',              'Legacy backward-compatibility key, read by no tracker as an option.')

% --- sel_power: removed from the hinec tracker entirely.
%     It re-weighted each stencil voxel by (alignment with the incoming
%     direction)^sel_power, so at sel_power=16 a voxel 45 deg off the current
%     heading contributed a weight of ~0.004 - near winner-take-all selection
%     biased toward where the track was already going. For DTI there is exactly
%     one principal eigenvector per voxel, so there is no ambiguity for it to
%     resolve; it simply bent a single-valued field toward the current heading,
%     which is self-reinforcing (it suppresses genuine sharp turns and locks in
%     error) and makes the ODE direction-dependent, dx/ds = v(x, dx/ds), so
%     classical RK order theory no longer applies.
%     HINEC tracking is now interpolation + integration only.
%     CSD still picks the FOD peak nearest the incoming tangent - that is
%     structurally required to reduce a multi-valued field to one direction, and
%     is not a tunable steering term. The alignment exponent on top of it is gone.
r('sel_power',          ['Removed. HINEC is now pure interpolation + integration. It biased ' ...
                         'interpolation toward the incoming direction (weight = alignment^sel_power), ' ...
                         'which has no justification for DTI (one eigenvector per voxel, nothing to ' ...
                         'disambiguate) and made the ODE direction-dependent. For CSD, nearest-peak ' ...
                         'selection is retained because it is structural, but the alignment exponent is gone.'])
];

% NOTE - deliberately NOT listed here, because they are MIGRATED (transformed),
% not dropped. They are handled by migrate_legacy() in load_config_yaml:
%   integration_order -> tractography.integrator.method  (1|2|4|5 -> euler|rk2|rk4|rkf45)
%   adaptive_step     -> tractography.integrator.adaptive (simple rename; NOT
%                        implied by the method, because rkf45 with adaptive:false
%                        is a genuine fixed-step RKF45 mode - hinec.m:1170)
%   max_steps         -> tractography.termination.max_arc  (max_arc = max_steps * step_size)
% Simple 1:1 renames are declared in the `legacy` field of nim_config_schema.

end
