function run_info = create_run_directory(config_file, varargin)
% create_run_directory: Create organized run directory with proper structure
%
% Creates a timestamped directory for a HINEC pipeline run with subdirectories
% for logs, intermediate files, output data, and tractography results.
%
% Arguments:
%   config_file - Path to YAML configuration file used for this run
%   varargin    - Optional name-value pairs:
%                 'base_dir'   - Base directory for runs (default: 'hinec_runs')
%                 'run_name'   - Custom run name (default: auto-generated)
%                 'description' - Run description for metadata
%
% Returns:
%   run_info - Structure with fields:
%              .run_dir        - Full path to run directory
%              .config_file    - Path to copied config file
%              .logs_dir       - Path to logs subdirectory
%              .intermediate_dir - Path to intermediate files
%              .output_dir     - Path to output files
%              .tractography_dir - Path to tractography results
%              .diagnostics_dir  - Path to diagnostics
%              .run_id         - Unique run identifier
%              .timestamp      - Run timestamp string
%
% Example:
%   run_info = create_run_directory('config/high_precision.yml');
%   run_info = create_run_directory('config/irontract.yml', 'run_name', 'my_experiment');
%   run_info = create_run_directory('config/fast_exploration.yml', 'base_dir', 'experiments');

    % Convert config_file to char array for consistent handling
    config_file = char(config_file);

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'base_dir', 'hinec_runs', @ischar);
    addParameter(p, 'run_name', '', @ischar);
    addParameter(p, 'description', '', @ischar);
    parse(p, varargin{:});

    base_dir = char(p.Results.base_dir);
    custom_run_name = char(p.Results.run_name);
    description = char(p.Results.description);

    % Create base directory if it doesn't exist
    if ~exist(base_dir, 'dir')
        mkdir(base_dir);
        fprintf('Created base directory: %s\n', base_dir);
    end

    % Generate timestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    % Extract config preset name from config file
    [~, config_name, ~] = fileparts(config_file);

    % Create run directory name
    if isempty(custom_run_name)
        run_name = sprintf('run_%s_%s', timestamp, config_name);
    else
        run_name = custom_run_name;
    end

    run_dir = fullfile(base_dir, run_name);

    % Check if directory already exists
    if exist(run_dir, 'dir')
        warning('Run directory already exists: %s', run_dir);
        % Append counter to make it unique
        counter = 1;
        while exist(sprintf('%s_%d', run_dir, counter), 'dir')
            counter = counter + 1;
        end
        run_dir = sprintf('%s_%d', run_dir, counter);
        fprintf('Using unique name: %s\n', run_dir);
    end

    % Create run directory
    mkdir(run_dir);
    fprintf('Created run directory: %s\n', run_dir);

    % Create subdirectories
    logs_dir = fullfile(run_dir, 'logs');
    intermediate_dir = fullfile(run_dir, 'intermediate');
    output_dir = fullfile(run_dir, 'output');
    tractography_dir = fullfile(run_dir, 'tractography');
    diagnostics_dir = fullfile(tractography_dir, 'diagnostics');

    mkdir(logs_dir);
    mkdir(intermediate_dir);
    mkdir(output_dir);
    mkdir(tractography_dir);
    mkdir(diagnostics_dir);

    fprintf('Created subdirectories:\n');
    fprintf('  logs/         - Pipeline logs\n');
    fprintf('  intermediate/ - Preprocessing artifacts\n');
    fprintf('  output/       - Final processed data\n');
    fprintf('  tractography/ - Tractography results\n');
    fprintf('  tractography/diagnostics/ - Diagnostic reports\n');

    % Copy configuration file to run directory
    config_copy = fullfile(run_dir, 'config.yml');
    copyfile(config_file, config_copy);
    fprintf('Copied config: %s -> %s\n', config_file, config_copy);

    % Create run metadata file
    run_info_file = fullfile(run_dir, 'run_info.txt');
    fid = fopen(run_info_file, 'w');

    fprintf(fid, 'HINEC Pipeline Run Information\n');
    fprintf(fid, '==============================\n\n');
    fprintf(fid, 'Run ID: %s\n', run_name);
    fprintf(fid, 'Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, 'Configuration: %s\n', config_file);
    fprintf(fid, 'Config Preset: %s\n\n', config_name);

    if ~isempty(description)
        fprintf(fid, 'Description: %s\n\n', description);
    end

    % System information
    fprintf(fid, 'System Information:\n');
    fprintf(fid, '  MATLAB Version: %s\n', version);
    fprintf(fid, '  Platform: %s\n', computer);
    fprintf(fid, '  Working Directory: %s\n\n', pwd);

    % Git information (if available)
    try
        [status, git_commit] = system('git rev-parse HEAD 2>/dev/null');
        if status == 0
            fprintf(fid, 'Git Information:\n');
            fprintf(fid, '  Commit: %s', git_commit);
            [~, git_branch] = system('git rev-parse --abbrev-ref HEAD 2>/dev/null');
            fprintf(fid, '  Branch: %s', git_branch);
            [~, git_status] = system('git status --porcelain 2>/dev/null');
            if ~isempty(strtrim(git_status))
                fprintf(fid, '  Status: MODIFIED (uncommitted changes)\n');
            else
                fprintf(fid, '  Status: CLEAN\n');
            end
            fprintf(fid, '\n');
        end
    catch
        % Git not available, skip
    end

    fprintf(fid, 'Directory Structure:\n');
    fprintf(fid, '  Run Directory: %s\n', run_dir);
    fprintf(fid, '  Logs: %s\n', logs_dir);
    fprintf(fid, '  Intermediate: %s\n', intermediate_dir);
    fprintf(fid, '  Output: %s\n', output_dir);
    fprintf(fid, '  Tractography: %s\n', tractography_dir);
    fprintf(fid, '  Diagnostics: %s\n', diagnostics_dir);

    fclose(fid);
    fprintf('Created run metadata: %s\n', run_info_file);

    % Update 'latest' symlink (platform-dependent)
    latest_link = fullfile(base_dir, 'latest');
    try
        if isunix
            % Remove old symlink if exists
            if exist(latest_link, 'dir') || exist(latest_link, 'file')
                delete(latest_link);
            end
            % Create new symlink
            system(sprintf('ln -s "%s" "%s"', run_name, latest_link));
            fprintf('Updated symlink: latest -> %s\n', run_name);
        else
            % Windows: Create a text file with the path instead
            fid = fopen(latest_link, 'w');
            fprintf(fid, '%s', run_dir);
            fclose(fid);
            fprintf('Updated latest marker: %s\n', latest_link);
        end
    catch ME
        warning('Could not create latest link: %s', ME.message);
    end

    % Build return structure
    run_info = struct();
    run_info.run_dir = run_dir;
    run_info.config_file = config_copy;
    run_info.logs_dir = logs_dir;
    run_info.intermediate_dir = intermediate_dir;
    run_info.output_dir = output_dir;
    run_info.tractography_dir = tractography_dir;
    run_info.diagnostics_dir = diagnostics_dir;
    run_info.run_id = run_name;
    run_info.timestamp = timestamp;

    fprintf('\n✅ Run directory created successfully\n');
    fprintf('Run ID: %s\n', run_name);
    fprintf('Location: %s\n\n', run_dir);
end
