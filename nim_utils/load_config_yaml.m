function config = load_config_yaml(config_file)
% load_config_yaml: Load HINEC configuration from YAML file
%
% Arguments:
%   config_file - Path to YAML configuration file
%
% Returns:
%   config - Structure with configuration parameters
%
% Usage:
%   config = load_config_yaml('config/hinec_default.yml');
%   config = load_config_yaml('my_experiment.yml');

    if nargin < 1 || isempty(config_file)
        % Use default config
        config_file = 'config/hinec_default.yml';
    end

    if ~isfile(config_file)
        error('Configuration file not found: %s', config_file);
    end

    fprintf('Loading configuration from: %s\n', config_file);

    % Read YAML file
    config = read_yaml_file(config_file);

    % Validate configuration
    config = validate_config(config);

    fprintf('Configuration loaded successfully\n');
end

function config = read_yaml_file(filename)
% Read and parse YAML file (simple parser for MATLAB)
%
% This is a basic YAML parser that handles the HINEC config format.
% For complex YAML features, consider using YAML toolbox:
% https://github.com/ewiger/yamlmatlab

    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    config = struct();
    current_section = '';
    line_num = 0;

    try
        while ~feof(fid)
            line = fgetl(fid);
            line_num = line_num + 1;

            % Skip empty lines and comments
            line = strtrim(line);
            if isempty(line) || startsWith(line, '#')
                continue;
            end

            % Check for section header (e.g., "preprocessing:")
            if endsWith(line, ':') && ~contains(line, ' ')
                current_section = line(1:end-1);
                config.(current_section) = struct();
                continue;
            end

            % Parse key-value pair
            if contains(line, ':')
                % Split on first colon
                colon_idx = find(line == ':', 1);
                key = strtrim(line(1:colon_idx-1));
                value_str = strtrim(line(colon_idx+1:end));

                % Remove inline comments
                comment_idx = strfind(value_str, '#');
                if ~isempty(comment_idx)
                    value_str = strtrim(value_str(1:comment_idx(1)-1));
                end

                % Skip if commented out or empty
                if isempty(value_str)
                    continue;
                end

                % Parse value
                value = parse_yaml_value(value_str);

                % Store in appropriate section
                if isempty(current_section)
                    config.(key) = value;
                else
                    config.(current_section).(key) = value;
                end
            end
        end
    catch ME
        fclose(fid);
        error('Error parsing YAML at line %d: %s', line_num, ME.message);
    end

    fclose(fid);
end

function value = parse_yaml_value(str)
% Parse YAML value from string

    % Remove quotes if present
    if (startsWith(str, '''') && endsWith(str, '''')) || ...
       (startsWith(str, '"') && endsWith(str, '"'))
        value = str(2:end-1);
        return;
    end

    % Boolean values
    if strcmpi(str, 'true')
        value = true;
        return;
    elseif strcmpi(str, 'false')
        value = false;
        return;
    end

    % Null/empty
    if strcmpi(str, 'null') || strcmpi(str, '~')
        value = [];
        return;
    end

    % Try to parse as number
    num = str2double(str);
    if ~isnan(num)
        value = num;
        return;
    end

    % Default to string
    value = str;
end

function config = validate_config(config)
% Validate and set defaults for missing fields

    % Preprocessing defaults
    if ~isfield(config, 'preprocessing')
        config.preprocessing = struct();
    end
    config.preprocessing = set_defaults(config.preprocessing, struct(...
        'run_denoising', true, ...
        'denoise_method', 'dwidenoise', ...
        'run_motion_correction', true, ...
        'run_eddy', true, ...
        'improve_mask', true, ...
        'atlas_type', 'jhu', ...
        't1_available', false, ...
        'use_t1_registration', false, ...
        'register_to_mni', false));

    % Tractography defaults
    if ~isfield(config, 'tractography')
        config.tractography = struct();
    end
    config.tractography = set_defaults(config.tractography, struct(...
        'algorithm', 'hinec', ...
        'integration_order', 4, ...
        'interp_method', 'trilinear', ...
        'step_size', 0.2, ...
        'adaptive_step', false, ...
        'rkf_tolerance', 0.01, ...
        'rkf_safety', 0.9, ...
        'step_min', 0.01, ...
        'step_max', 1.0, ...
        'seed_density', 4, ...
        'seed_strategy', 'uniform', ...
        'termination_fa', 0.15, ...
        'angle_thresh', 35, ...
        'max_steps', 1000, ...
        'min_length', 35, ...
        'act_enabled', true, ...
        'enable_diagnostics', true));

    % Validate tractography parameters
    validate_tractography_params(config.tractography);
end

function s = set_defaults(s, defaults)
% Set default values for missing fields in structure

    fields = fieldnames(defaults);
    for i = 1:length(fields)
        field = fields{i};
        if ~isfield(s, field)
            s.(field) = defaults.(field);
        end
    end
end

function validate_tractography_params(params)
% Validate tractography parameters

    % Step size
    if params.step_size <= 0
        error('step_size must be positive, got: %.4f', params.step_size);
    end

    % Angle threshold
    if params.angle_thresh <= 0 || params.angle_thresh > 180
        error('angle_thresh must be in (0, 180], got: %.1f', params.angle_thresh);
    end

    % Integration order
    valid_orders = [1, 2, 4, 5];
    if ~ismember(params.integration_order, valid_orders)
        error('integration_order must be 1, 2, 4, or 5, got: %d', params.integration_order);
    end

    % RKF45-specific validation
    if params.integration_order == 5
        if params.rkf_tolerance <= 0
            error('rkf_tolerance must be positive, got: %.4f', params.rkf_tolerance);
        end
        if params.step_min <= 0 || params.step_min >= params.step_max
            error('Invalid step bounds: step_min=%.4f, step_max=%.4f', ...
                  params.step_min, params.step_max);
        end
        if params.rkf_safety <= 0 || params.rkf_safety > 1
            error('rkf_safety must be in (0, 1], got: %.2f', params.rkf_safety);
        end
    end

    % Seed density
    if params.seed_density < 1 || params.seed_density > 20
        warning('Unusual seed_density: %d (recommended: 1-8)', params.seed_density);
    end

    % FA threshold
    if params.termination_fa < 0 || params.termination_fa > 1
        error('termination_fa must be in [0, 1], got: %.2f', params.termination_fa);
    end

    fprintf('✓ Parameter validation passed\n');
end
