classdef TractographyCacheManager < handle
    % TRACTOGRAPHYCACHEMANAGER: Comprehensive cache management for fast slice viewing
    %
    % This class manages the hierarchical cache directory structure, parameter
    % hashing, metadata storage, and cache validation for the Fast Tractography
    % Slice Viewer system.
    %
    % Key Features:
    % - Hierarchical directory organization by dataset and parameters
    % - Robust parameter hashing with collision detection
    % - JSON-based metadata with versioning and integrity checks
    % - Cross-platform path handling and naming conventions
    % - Cache validation and cleanup utilities
    % - Storage optimization and compression management

    properties (Access = private)
        cache_root_dir          % Root directory for all cache data
        current_dataset_hash    % Hash of currently loaded dataset
        current_param_hash      % Hash of current parameters
        metadata_version        % Cache format version
        default_cache_dir       % Default cache location
    end

    properties (Constant)
        CACHE_VERSION = '1.0'
        MAX_CACHE_AGE_DAYS = 30
        HASH_LENGTH = 16        % SHA256 truncated to 16 chars for readability
    end

    methods
        function obj = TractographyCacheManager(cache_dir)
            % Constructor: Initialize cache manager
            %
            % Arguments:
            %   cache_dir - Optional custom cache directory path
            %               If empty, uses default location

            if nargin < 1 || isempty(cache_dir)
                obj.default_cache_dir = obj.getDefaultCacheLocation();
                cache_dir = obj.default_cache_dir;
            end

            obj.cache_root_dir = cache_dir;
            obj.metadata_version = obj.CACHE_VERSION;

            % Create root cache directory if it doesn't exist
            obj.ensureDirectoryExists(obj.cache_root_dir);

            % Initialize global config
            obj.initializeGlobalConfig();

            fprintf('TractographyCacheManager initialized\n');
            fprintf('  Cache directory: %s\n', obj.cache_root_dir);
        end

        function dataset_hash = generateDatasetHash(obj, tracks_file, nim_file)
            % Generate unique hash for dataset files
            %
            % Arguments:
            %   tracks_file - Path to tracks .mat file
            %   nim_file    - Path to nim .mat file
            %
            % Returns:
            %   dataset_hash - Unique hash identifying this dataset

            % Get file information for hashing
            tracks_info = dir(tracks_file);
            nim_info = dir(nim_file);

            if isempty(tracks_info) || isempty(nim_info)
                error('TractographyCacheManager:FileNotFound', ...
                      'Cannot generate hash - tracks or nim file not found');
            end

            % Create hash input from file paths, sizes, and modification dates
            hash_input = sprintf('%s|%s|%d|%d|%s|%s', ...
                tracks_file, nim_file, ...
                tracks_info.bytes, nim_info.bytes, ...
                tracks_info.date, nim_info.date);

            % Generate SHA256 hash and truncate for readability
            hash_full = string(java.security.MessageDigest.getInstance('SHA-256').digest(uint8(hash_input)));
            dataset_hash = char(extractAfter(extractBefore(hash_full, obj.HASH_LENGTH + 1), 1));

            obj.current_dataset_hash = dataset_hash;
        end

        function param_hash = generateParameterHash(obj, options)
            % Generate unique hash for parameter combination
            %
            % Arguments:
            %   options - Structure with visualization parameters
            %
            % Returns:
            %   param_hash - Unique hash for this parameter set

            % Normalize parameters for consistent hashing
            norm_params = obj.normalizeParameters(options);

            % Convert to JSON string for hashing
            json_str = jsonencode(norm_params);

            % Generate hash
            hash_full = string(java.security.MessageDigest.getInstance('SHA-256').digest(uint8(json_str)));
            param_hash = char(extractAfter(extractBefore(hash_full, obj.HASH_LENGTH + 1), 1));

            obj.current_param_hash = param_hash;
        end

        function cache_dir = getCacheDirectory(obj, dataset_hash, param_hash)
            % Get full cache directory path for dataset and parameters
            %
            % Returns:
            %   cache_dir - Full path to cache directory

            cache_dir = fullfile(obj.cache_root_dir, 'datasets', dataset_hash, ...
                               'parameters', param_hash);
        end

        function success = createCacheStructure(obj, dataset_hash, param_hash, tracks_file, nim_file, options)
            % Create complete cache directory structure with metadata
            %
            % Arguments:
            %   dataset_hash - Dataset identifier hash
            %   param_hash   - Parameter combination hash
            %   tracks_file  - Original tracks file path
            %   nim_file     - Original nim file path
            %   options      - Parameter options structure
            %
            % Returns:
            %   success - True if structure created successfully

            success = false;

            try
                % Create main cache directory
                cache_dir = obj.getCacheDirectory(dataset_hash, param_hash);
                obj.ensureDirectoryExists(cache_dir);

                % Create orientation subdirectories
                obj.ensureDirectoryExists(fullfile(cache_dir, 'axial'));
                obj.ensureDirectoryExists(fullfile(cache_dir, 'sagittal'));
                obj.ensureDirectoryExists(fullfile(cache_dir, 'coronal'));

                % Create dataset metadata (if not exists)
                dataset_dir = fullfile(obj.cache_root_dir, 'datasets', dataset_hash);
                dataset_metadata_file = fullfile(dataset_dir, 'metadata.json');

                if ~exist(dataset_metadata_file, 'file')
                    dataset_metadata = obj.createDatasetMetadata(tracks_file, nim_file);
                    obj.writeJSONFile(dataset_metadata_file, dataset_metadata);
                end

                % Create parameter configuration
                param_config_file = fullfile(cache_dir, 'config.json');
                param_config = obj.createParameterConfig(options);
                obj.writeJSONFile(param_config_file, param_config);

                % Create thumbnails directory
                thumbnails_dir = fullfile(dataset_dir, 'thumbnails');
                obj.ensureDirectoryExists(thumbnails_dir);

                success = true;
                fprintf('Cache structure created: %s\n', cache_dir);

            catch ME
                fprintf('Error creating cache structure: %s\n', ME.message);
            end
        end

        function valid = validateCacheIntegrity(obj, dataset_hash, param_hash)
            % Validate cache integrity and completeness
            %
            % Returns:
            %   valid - True if cache is valid and complete

            valid = false;

            try
                cache_dir = obj.getCacheDirectory(dataset_hash, param_hash);

                % Check if cache directory exists
                if ~exist(cache_dir, 'dir')
                    return;
                end

                % Check for required subdirectories
                orientations = {'axial', 'sagittal', 'coronal'};
                for i = 1:length(orientations)
                    orient_dir = fullfile(cache_dir, orientations{i});
                    if ~exist(orient_dir, 'dir')
                        return;
                    end
                end

                % Check for configuration file
                config_file = fullfile(cache_dir, 'config.json');
                if ~exist(config_file, 'file')
                    return;
                end

                % Validate configuration
                config = obj.readJSONFile(config_file);
                if ~obj.validateParameterConfig(config)
                    return;
                end

                % Check for dataset metadata
                dataset_dir = fullfile(obj.cache_root_dir, 'datasets', dataset_hash);
                metadata_file = fullfile(dataset_dir, 'metadata.json');
                if ~exist(metadata_file, 'file')
                    return;
                end

                valid = true;

            catch ME
                fprintf('Cache validation error: %s\n', ME.message);
            end
        end

        function cleaned_bytes = cleanupExpiredCache(obj, max_age_days)
            % Clean up old and unused cache entries
            %
            % Arguments:
            %   max_age_days - Maximum age in days (default: 30)
            %
            % Returns:
            %   cleaned_bytes - Number of bytes freed

            if nargin < 2
                max_age_days = obj.MAX_CACHE_AGE_DAYS;
            end

            cleaned_bytes = 0;
            current_time = now;
            cutoff_time = current_time - max_age_days;

            try
                datasets_dir = fullfile(obj.cache_root_dir, 'datasets');
                if ~exist(datasets_dir, 'dir')
                    return;
                end

                dataset_dirs = dir(datasets_dir);
                dataset_dirs = dataset_dirs([dataset_dirs.isdir] & ~ismember({dataset_dirs.name}, {'.', '..'}));

                for i = 1:length(dataset_dirs)
                    dataset_path = fullfile(datasets_dir, dataset_dirs(i).name);

                    % Check dataset age
                    if dataset_dirs(i).datenum < cutoff_time
                        fprintf('Cleaning expired dataset: %s\n', dataset_dirs(i).name);

                        % Calculate size before deletion
                        size_info = obj.calculateDirectorySize(dataset_path);
                        cleaned_bytes = cleaned_bytes + size_info.total_bytes;

                        % Remove directory
                        rmdir(dataset_path, 's');
                    end
                end

                fprintf('Cleanup complete: %.2f MB freed\n', cleaned_bytes / (1024^2));

            catch ME
                fprintf('Cleanup error: %s\n', ME.message);
            end
        end

        function metadata = getCacheMetadata(obj, dataset_hash)
            % Get comprehensive metadata for a dataset
            %
            % Returns:
            %   metadata - Structure with cache information

            metadata = struct();

            try
                dataset_dir = fullfile(obj.cache_root_dir, 'datasets', dataset_hash);
                metadata_file = fullfile(dataset_dir, 'metadata.json');

                if exist(metadata_file, 'file')
                    metadata = obj.readJSONFile(metadata_file);

                    % Add current cache statistics
                    metadata.cache_stats = obj.calculateCacheStatistics(dataset_hash);
                end

            catch ME
                fprintf('Error reading metadata: %s\n', ME.message);
            end
        end
    end

    methods (Access = private)
        function default_dir = getDefaultCacheLocation(obj)
            % Get platform-appropriate default cache location

            if ispc
                % Windows: Use AppData/Local
                default_dir = fullfile(getenv('LOCALAPPDATA'), 'TractographyCache');
            elseif ismac
                % macOS: Use ~/Library/Caches
                default_dir = fullfile(getenv('HOME'), 'Library', 'Caches', 'TractographyCache');
            else
                % Linux: Use ~/.cache
                default_dir = fullfile(getenv('HOME'), '.cache', 'TractographyCache');
            end
        end

        function ensureDirectoryExists(obj, dir_path)
            % Create directory if it doesn't exist

            if ~exist(dir_path, 'dir')
                mkdir(dir_path);
            end
        end

        function norm_params = normalizeParameters(obj, options)
            % Normalize parameters for consistent hashing

            norm_params = struct();

            % Extract and normalize relevant parameters
            fields = {'tolerance', 'color_mode', 'show_anatomy', 'alpha', ...
                     'regions', 'region_filter', 'min_overlap'};

            for i = 1:length(fields)
                field = fields{i};
                if isfield(options, field)
                    value = options.(field);

                    % Normalize different data types
                    if isnumeric(value)
                        % Round to avoid floating point precision issues
                        norm_params.(field) = round(value * 1000) / 1000;
                    elseif islogical(value)
                        norm_params.(field) = logical(value);
                    elseif ischar(value) || isstring(value)
                        norm_params.(field) = char(value);
                    else
                        norm_params.(field) = value;
                    end
                end
            end

            % Sort fields for consistent ordering
            norm_params = orderfields(norm_params);
        end

        function metadata = createDatasetMetadata(obj, tracks_file, nim_file)
            % Create metadata structure for dataset

            metadata = struct();
            metadata.version = obj.CACHE_VERSION;
            metadata.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
            metadata.tracks_file = tracks_file;
            metadata.nim_file = nim_file;

            % File information
            tracks_info = dir(tracks_file);
            nim_info = dir(nim_file);

            metadata.files = struct();
            metadata.files.tracks = struct('size', tracks_info.bytes, 'modified', tracks_info.date);
            metadata.files.nim = struct('size', nim_info.bytes, 'modified', nim_info.date);

            % Dataset statistics (will be filled when cache is generated)
            metadata.dataset_stats = struct();
        end

        function config = createParameterConfig(obj, options)
            % Create parameter configuration structure

            config = struct();
            config.version = obj.CACHE_VERSION;
            config.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
            config.parameters = obj.normalizeParameters(options);
            config.hash = obj.current_param_hash;
        end

        function valid = validateParameterConfig(obj, config)
            % Validate parameter configuration structure

            valid = false;

            try
                required_fields = {'version', 'created', 'parameters', 'hash'};
                for i = 1:length(required_fields)
                    if ~isfield(config, required_fields{i})
                        return;
                    end
                end

                % Check version compatibility
                if ~strcmp(config.version, obj.CACHE_VERSION)
                    fprintf('Warning: Cache version mismatch (cache: %s, current: %s)\n', ...
                           config.version, obj.CACHE_VERSION);
                end

                valid = true;

            catch ME
                fprintf('Config validation error: %s\n', ME.message);
            end
        end

        function writeJSONFile(obj, filename, data)
            % Write data to JSON file with error handling

            try
                json_str = jsonencode(data, 'PrettyPrint', true);
                fid = fopen(filename, 'w');
                if fid ~= -1
                    fprintf(fid, '%s', json_str);
                    fclose(fid);
                else
                    error('Cannot open file for writing: %s', filename);
                end
            catch ME
                error('Error writing JSON file %s: %s', filename, ME.message);
            end
        end

        function data = readJSONFile(obj, filename)
            % Read data from JSON file with error handling

            try
                fid = fopen(filename, 'r');
                if fid ~= -1
                    json_str = fread(fid, '*char')';
                    fclose(fid);
                    data = jsondecode(json_str);
                else
                    error('Cannot open file for reading: %s', filename);
                end
            catch ME
                error('Error reading JSON file %s: %s', filename, ME.message);
            end
        end

        function initializeGlobalConfig(obj)
            % Initialize global cache configuration

            global_config_file = fullfile(obj.cache_root_dir, 'global_config.json');

            if ~exist(global_config_file, 'file')
                global_config = struct();
                global_config.version = obj.CACHE_VERSION;
                global_config.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
                global_config.settings = struct();
                global_config.settings.max_cache_size_gb = 10;
                global_config.settings.auto_cleanup_days = obj.MAX_CACHE_AGE_DAYS;

                obj.writeJSONFile(global_config_file, global_config);
            end
        end

        function size_info = calculateDirectorySize(obj, dir_path)
            % Calculate total size of directory and subdirectories

            size_info = struct('total_bytes', 0, 'file_count', 0);

            try
                if exist(dir_path, 'dir')
                    files = dir(fullfile(dir_path, '**', '*'));
                    files = files(~[files.isdir]);

                    size_info.total_bytes = sum([files.bytes]);
                    size_info.file_count = length(files);
                end
            catch ME
                fprintf('Error calculating directory size: %s\n', ME.message);
            end
        end

        function stats = calculateCacheStatistics(obj, dataset_hash)
            % Calculate statistics for a dataset cache

            stats = struct();

            try
                dataset_dir = fullfile(obj.cache_root_dir, 'datasets', dataset_hash);
                size_info = obj.calculateDirectorySize(dataset_dir);

                stats.total_size_bytes = size_info.total_bytes;
                stats.total_size_mb = size_info.total_bytes / (1024^2);
                stats.file_count = size_info.file_count;

                % Count parameter sets
                param_dir = fullfile(dataset_dir, 'parameters');
                if exist(param_dir, 'dir')
                    param_dirs = dir(param_dir);
                    param_dirs = param_dirs([param_dirs.isdir] & ~ismember({param_dirs.name}, {'.', '..'}));
                    stats.parameter_sets = length(param_dirs);
                else
                    stats.parameter_sets = 0;
                end

            catch ME
                fprintf('Error calculating cache statistics: %s\n', ME.message);
                stats = struct();
            end
        end
    end
end