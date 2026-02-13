function launchFastViewer(tracks_file, nim_file, options)
% LAUNCHFASTVIEWER: MATLAB bridge to Python FastTractographyViewer
%
% This function provides seamless integration between MATLAB tractography
% processing and the Python-based fast slice viewer. It automatically
% generates the cache if needed and launches the viewer.
%
% Usage:
%   launchFastViewer(tracks_file, nim_file)
%   launchFastViewer(tracks_file, nim_file, options)
%
% Arguments:
%   tracks_file - Path to tracks .mat file
%   nim_file    - Path to nim .mat file
%   options     - Optional parameters structure
%
% Options:
%   .cache_dir         - Custom cache directory (default: auto)
%   .force_regenerate  - Force cache regeneration (default: false)
%   .python_path       - Python executable path (default: 'python')
%   .viewer_path       - Path to FastTractographyViewer.py (default: auto)
%   .auto_launch       - Launch viewer after cache generation (default: true)
%   .background_mode   - Run cache generation in background (default: false)
%
% Example:
%   % Basic usage
%   launchFastViewer('tracks.mat', 'nim.mat');
%
%   % With custom options
%   opts.force_regenerate = true;
%   opts.cache_dir = '/custom/cache/path';
%   launchFastViewer('tracks.mat', 'nim.mat', opts);

fprintf('=== FAST TRACTOGRAPHY VIEWER LAUNCHER ===\n');

%% Input validation
if nargin < 2
    error('launchFastViewer:InvalidInput', ...
          'Must provide tracks_file and nim_file');
end

if nargin < 3
    options = struct();
end

% Set default options
options = setDefaultLauncherOptions(options);

% Validate input files
if ~exist(tracks_file, 'file')
    error('launchFastViewer:FileNotFound', ...
          'Tracks file not found: %s', tracks_file);
end

if ~exist(nim_file, 'file')
    error('launchFastViewer:FileNotFound', ...
          'NIM file not found: %s', nim_file);
end

%% Determine cache directory
if isempty(options.cache_dir)
    % Use default cache location based on platform
    if ispc
        cache_root = fullfile(getenv('LOCALAPPDATA'), 'TractographyCache');
    elseif ismac
        cache_root = fullfile(getenv('HOME'), 'Library', 'Caches', 'TractographyCache');
    else
        cache_root = fullfile(getenv('HOME'), '.cache', 'TractographyCache');
    end
    options.cache_dir = cache_root;
end

fprintf('Cache directory: %s\n', options.cache_dir);

%% Check if cache exists and is valid
cache_manager = TractographyCacheManager(options.cache_dir);
dataset_hash = cache_manager.generateDatasetHash(tracks_file, nim_file);
param_hash = cache_manager.generateParameterHash(options);

cache_exists = cache_manager.validateCacheIntegrity(dataset_hash, param_hash);

if cache_exists && ~options.force_regenerate
    fprintf('Valid cache found, skipping generation\n');
else
    fprintf('Generating cache (this may take several minutes)...\n');

    %% Generate cache
    if options.background_mode
        generateCacheBackground(tracks_file, nim_file, options.cache_dir, options);
    else
        success = generateTractographySliceCache(tracks_file, nim_file, options.cache_dir, options);
        if ~success
            error('launchFastViewer:CacheGenerationFailed', ...
                  'Failed to generate cache');
        end
    end
end

%% Launch Python viewer
if options.auto_launch
    launchPythonViewer(options);
else
    fprintf('Cache ready. Launch viewer manually with:\n');
    fprintf('  %s "%s" "%s"\n', options.python_path, options.viewer_path, options.cache_dir);
end

end


function options = setDefaultLauncherOptions(options)
% Set default options for the launcher

% Get script directory for default viewer path
script_dir = fileparts(mfilename('fullpath'));

default_opts = struct();
default_opts.cache_dir = '';
default_opts.force_regenerate = false;
default_opts.python_path = 'python';
default_opts.viewer_path = fullfile(script_dir, 'FastTractographyViewer.py');
default_opts.auto_launch = true;
default_opts.background_mode = false;

% Merge with provided options
field_names = fieldnames(default_opts);
for i = 1:length(field_names)
    field = field_names{i};
    if ~isfield(options, field)
        options.(field) = default_opts.(field);
    end
end

% Validate Python path
if strcmp(options.python_path, 'python')
    % Try to find Python executable
    if ispc
        [status, ~] = system('python --version');
        if status ~= 0
            [status, ~] = system('py --version');
            if status == 0
                options.python_path = 'py';
            end
        end
    end
end

% Validate viewer path
if ~exist(options.viewer_path, 'file')
    warning('launchFastViewer:ViewerNotFound', ...
            'FastTractographyViewer.py not found at: %s', options.viewer_path);
end

end


function generateCacheBackground(tracks_file, nim_file, cache_dir, options)
% Generate cache in background using parallel processing

fprintf('Starting background cache generation...\n');

% Create a timer to monitor progress
progress_timer = timer('ExecutionMode', 'fixedRate', 'Period', 5, ...
                      'TimerFcn', @(~,~) fprintf('Cache generation in progress...\n'));

try
    start(progress_timer);

    % Run cache generation
    success = generateTractographySliceCache(tracks_file, nim_file, cache_dir, options);

    stop(progress_timer);
    delete(progress_timer);

    if success
        fprintf('Background cache generation completed successfully\n');
    else
        fprintf('Background cache generation failed\n');
    end

catch ME
    stop(progress_timer);
    delete(progress_timer);
    rethrow(ME);
end

end


function launchPythonViewer(options)
% Launch the Python viewer application

fprintf('Launching Fast Tractography Viewer...\n');

% Construct command
cmd = sprintf('"%s" "%s" "%s"', options.python_path, options.viewer_path, options.cache_dir);

fprintf('Executing: %s\n', cmd);

try
    if ispc
        % On Windows, use system with proper handling
        system(cmd);
    else
        % On Unix systems, run in background
        system([cmd ' &']);
    end

    fprintf('Viewer launched successfully\n');

catch ME
    fprintf('Failed to launch viewer: %s\n', ME.message);
    fprintf('You can launch manually with:\n');
    fprintf('  %s "%s" "%s"\n', options.python_path, options.viewer_path, options.cache_dir);
end

end


function checkPythonDependencies(python_path)
% Check if required Python dependencies are installed

fprintf('Checking Python dependencies...\n');

required_packages = {'tkinter', 'PIL', 'numpy', 'scipy', 'matplotlib'};

for i = 1:length(required_packages)
    package = required_packages{i};

    % Try to import the package
    import_cmd = sprintf('"%s" -c "import %s; print(''%s: OK'')"', python_path, package, package);

    [status, output] = system(import_cmd);

    if status == 0
        fprintf('  %s\n', strtrim(output));
    else
        fprintf('  %s: MISSING\n', package);
        fprintf('    Install with: pip install %s\n', getDependencyInstallName(package));
    end
end

end


function install_name = getDependencyInstallName(package)
% Get the correct pip install name for a package

switch package
    case 'PIL'
        install_name = 'pillow';
    case 'tkinter'
        install_name = 'tk';  % Usually built-in
    otherwise
        install_name = package;
end

end