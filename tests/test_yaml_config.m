%% test_yaml_config.m - Test YAML configuration system
% Verifies that YAML config loading and validation works correctly

fprintf('=== HINEC YAML Configuration System Test ===\n\n');

% Add paths
addpath('src/nim_utils');

%% Test 1: Load default config
fprintf('Test 1: Loading default configuration...\n');
try
    config = load_config_yaml('config/hinec_default.yml');
    fprintf('  ✓ Default config loaded successfully\n');
    fprintf('    Algorithm: %s\n', config.tractography.algorithm);
    fprintf('    Integrator: %s\n', config.tractography.integrator.method);
    fprintf('    Step size: %.2f\n', config.tractography.integrator.step);
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    return;
end

%% Test 2: Load an RKF45 tracker config
fprintf('\nTest 2: Loading RKF45 tracker configuration (hinec_dti)...\n');
try
    config_hp = load_config_yaml('config/hinec_dti.yml');
    fprintf('  ✓ hinec_dti config loaded\n');
    fprintf('    Integrator: %s\n', config_hp.tractography.integrator.method);
    fprintf('    Adaptive step: %d\n', config_hp.tractography.integrator.adaptive);
    fprintf('    RKF tolerance: %.4f\n', config_hp.tractography.integrator.tolerance);
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    return;
end

%% Test 3: Validate parameter constraints
fprintf('\nTest 3: Parameter validation...\n');
try
    % This should pass validation
    config = load_config_yaml('config/hinec_default.yml');
    fprintf('  ✓ Default parameters pass validation\n');

    % Check key parameters are within expected ranges
    assert(config.tractography.integrator.step > 0, 'Step size validation');
    assert(config.tractography.termination.angle_max > 0 && config.tractography.termination.angle_max <= 180, 'Angle threshold validation');
    assert(config.tractography.seeding.density >= 1, 'Seed density validation');
    fprintf('  ✓ All parameter ranges validated\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    return;
end

%% Test 4: Check all preset configs exist
fprintf('\nTest 4: Checking preset configurations...\n');
presets = {
    'config/hinec_default.yml', 'Default';
    'config/hinec_dti.yml',     'HINEC DTI';
    'config/hinec_csd.yml',     'HINEC CSD';
    'config/mmf_dti.yml',       'MMF DTI';
    'config/standard_dti.yml',  'Standard FACT';
    'config/irontract.yml',     'IronTract';
    'config/ismrm2015.yml',     'ISMRM-2015'
};

all_exist = true;
for i = 1:size(presets, 1)
    if isfile(presets{i, 1})
        fprintf('  ✓ %s: %s\n', presets{i, 2}, presets{i, 1});
    else
        fprintf('  ✗ MISSING: %s\n', presets{i, 1});
        all_exist = false;
    end
end

if all_exist
    fprintf('  ✓ All preset configs found\n');
else
    fprintf('  ✗ Some preset configs missing\n');
    return;
end

%% Test 5: RKF45 parameter validation
fprintf('\nTest 5: RKF45 parameter validation...\n');
try
    config_rkf = load_config_yaml('config/hinec_dti.yml');

    % Verify RKF45 parameters are present (nested under integrator: since the
    % schema; the flat rkf_* names are legacy aliases that migrate on load)
    integ = config_rkf.tractography.integrator;
    assert(isfield(integ, 'tolerance'), 'RKF tolerance field');
    assert(isfield(integ, 'safety'),    'RKF safety field');
    assert(isfield(integ, 'step_min'),  'Step min field');
    assert(isfield(integ, 'step_max'),  'Step max field');

    % Verify constraints
    assert(integ.tolerance > 0, 'RKF tolerance positive');
    assert(integ.step_min < integ.step_max, 'Step bounds');
    assert(integ.safety > 0 && integ.safety <= 1, 'Safety factor range');

    fprintf('  ✓ RKF45 parameters validated\n');
catch ME
    fprintf('  ✗ FAILED: %s\n', ME.message);
    return;
end

%% Test 6: Compare configs
fprintf('\nTest 6: Comparing configurations...\n');
config_default = load_config_yaml('config/hinec_default.yml');
config_fast = load_config_yaml('config/standard_dti.yml');
config_precise = load_config_yaml('config/hinec_dti.yml');

fprintf('\n  Parameter Comparison:\n');
fprintf('  %-20s %10s %10s %10s\n', 'Parameter', 'Default', 'Fast', 'Precise');
fprintf('  %-20s %10s %10s %10s\n', repmat('-', 1, 20), repmat('-', 1, 10), repmat('-', 1, 10), repmat('-', 1, 10));
t_def = config_default.tractography; t_fast = config_fast.tractography; t_prec = config_precise.tractography;
fprintf('  %-20s %10s %10s %10s\n', 'integrator.method', t_def.integrator.method, t_fast.integrator.method, t_prec.integrator.method);
fprintf('  %-20s %10.2f %10.2f %10.2f\n', 'integrator.step', t_def.integrator.step, t_fast.integrator.step, t_prec.integrator.step);
fprintf('  %-20s %10d %10d %10d\n', 'seeding.density', t_def.seeding.density, t_fast.seeding.density, t_prec.seeding.density);
fprintf('  %-20s %10.2f %10.2f %10.2f\n', 'termination.fa_min', t_def.termination.fa_min, t_fast.termination.fa_min, t_prec.termination.fa_min);

%% Summary
fprintf('\n========================================\n');
fprintf('ALL TESTS PASSED ✓\n');
fprintf('========================================\n');
fprintf('\nYAML configuration system is ready to use!\n');
fprintf('\nQuick start:\n');
fprintf('  ./run_hinec.sh data/subject processed.mat\n');
fprintf('  ./run_hinec.sh data/subject processed.mat config/ismrm2015.yml\n');
fprintf('\nSee docs/YAML_CONFIG.md for complete documentation.\n');
