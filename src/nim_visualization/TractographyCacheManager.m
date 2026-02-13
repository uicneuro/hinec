classdef TractographyCacheManager < handle
    % TRACTOGRAPHYCACHEMANAGER: Simple cache management for slice viewing
    %
    % Creates a simple, portable directory structure:
    %   output_dir/
    %   ├── metadata.json      (all info about this tractography)
    %   ├── axial/*.png
    %   ├── sagittal/*.png
    %   └── coronal/*.png

    properties (Access = private)
        output_dir              % Output directory for this cache
        metadata                % Metadata structure
    end

    properties (Constant)
        CACHE_VERSION = '2.0'
    end

    methods
        function obj = TractographyCacheManager(output_dir)
            % Constructor: Initialize cache manager
            %
            % Arguments:
            %   output_dir - Directory where cache will be created

            if nargin < 1 || isempty(output_dir)
                error('TractographyCacheManager:InvalidInput', ...
                      'Must provide output_dir');
            end

            obj.output_dir = output_dir;
            obj.metadata = struct();

            % Create output directory if it doesn't exist
            if ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end

            fprintf('Cache manager initialized\n');
            fprintf('  Output directory: %s\n', obj.output_dir);
        end

        function initializeCache(obj, tracks_file, nim_file, options)
            % Initialize cache with metadata
            %
            % Arguments:
            %   tracks_file - Path to tracks .mat file
            %   nim_file    - Path to nim .mat file
            %   options     - Generation options

            % Load data files to get information
            tracks_data = load(tracks_file);
            nim_data = load(nim_file);

            % Get file info
            tracks_info = dir(tracks_file);
            nim_info = dir(nim_file);

            % Build metadata
            obj.metadata = struct();
            obj.metadata.version = obj.CACHE_VERSION;
            obj.metadata.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');

            % Source files
            obj.metadata.source = struct();
            obj.metadata.source.tracks_file = tracks_file;
            obj.metadata.source.nim_file = nim_file;
            obj.metadata.source.tracks_size_bytes = tracks_info.bytes;
            obj.metadata.source.nim_size_bytes = nim_info.bytes;
            obj.metadata.source.tracks_modified = tracks_info.date;
            obj.metadata.source.nim_modified = nim_info.date;

            % Track information
            if isfield(tracks_data, 'tracks')
                obj.metadata.tracks = struct();
                obj.metadata.tracks.count = length(tracks_data.tracks);
                obj.metadata.tracks.total_points = sum(cellfun(@(t) size(t,1), tracks_data.tracks));
            end

            % Volume dimensions
            if isfield(nim_data, 'nim') && isfield(nim_data.nim, 'FA')
                obj.metadata.volume = struct();
                obj.metadata.volume.dimensions = size(nim_data.nim.FA);
            end

            % Rendering parameters
            obj.metadata.parameters = struct();
            obj.metadata.parameters.tolerance = options.tolerance;
            obj.metadata.parameters.color_mode = options.color_mode;
            obj.metadata.parameters.show_anatomy = options.show_anatomy;
            obj.metadata.parameters.alpha = options.alpha;
            obj.metadata.parameters.image_resolution = options.image_resolution;
            obj.metadata.parameters.image_format = options.image_format;
            obj.metadata.parameters.compression_level = options.compression_level;

            % Region filtering (if available)
            if isfield(tracks_data, 'region_filter')
                obj.metadata.region_filter = tracks_data.region_filter;
            end

            % Save metadata
            obj.saveMetadata();
        end

        function createOrientationDirs(obj)
            % Create subdirectories for each orientation

            orientations = {'axial', 'sagittal', 'coronal'};
            for i = 1:length(orientations)
                orient_dir = fullfile(obj.output_dir, orientations{i});
                if ~exist(orient_dir, 'dir')
                    mkdir(orient_dir);
                end
            end
        end

        function saveMetadata(obj)
            % Save metadata to JSON file

            metadata_file = fullfile(obj.output_dir, 'metadata.json');

            % Convert to JSON
            json_str = jsonencode(obj.metadata);

            % Write to file
            fid = fopen(metadata_file, 'w');
            if fid == -1
                error('TractographyCacheManager:FileError', ...
                      'Could not open metadata file for writing');
            end

            fprintf(fid, '%s', json_str);
            fclose(fid);

            fprintf('Metadata saved: %s\n', metadata_file);
        end

        function metadata = loadMetadata(obj)
            % Load metadata from JSON file

            metadata_file = fullfile(obj.output_dir, 'metadata.json');

            if ~exist(metadata_file, 'file')
                error('TractographyCacheManager:FileNotFound', ...
                      'Metadata file not found: %s', metadata_file);
            end

            % Read JSON file
            fid = fopen(metadata_file, 'r');
            json_str = fread(fid, '*char')';
            fclose(fid);

            % Parse JSON
            metadata = jsondecode(json_str);
            obj.metadata = metadata;
        end

        function is_valid = validateCache(obj)
            % Validate cache structure and completeness

            is_valid = false;

            % Check metadata exists
            metadata_file = fullfile(obj.output_dir, 'metadata.json');
            if ~exist(metadata_file, 'file')
                fprintf('Validation failed: metadata.json not found\n');
                return;
            end

            % Check orientation directories exist
            orientations = {'axial', 'sagittal', 'coronal'};
            for i = 1:length(orientations)
                orient_dir = fullfile(obj.output_dir, orientations{i});
                if ~exist(orient_dir, 'dir')
                    fprintf('Validation failed: %s directory not found\n', orientations{i});
                    return;
                end

                % Check for images
                image_files = dir(fullfile(orient_dir, '*.png'));
                if isempty(image_files)
                    image_files = dir(fullfile(orient_dir, '*.jpg'));
                end

                if isempty(image_files)
                    fprintf('Validation failed: No images in %s directory\n', orientations{i});
                    return;
                end
            end

            is_valid = true;
            fprintf('Cache validation passed\n');
        end

        function stats = getCacheStats(obj)
            % Get statistics about the cache

            stats = struct();
            stats.output_dir = obj.output_dir;

            orientations = {'axial', 'sagittal', 'coronal'};
            total_size = 0;
            total_images = 0;

            for i = 1:length(orientations)
                orient = orientations{i};
                orient_dir = fullfile(obj.output_dir, orient);

                if exist(orient_dir, 'dir')
                    % Count images
                    png_files = dir(fullfile(orient_dir, '*.png'));
                    jpg_files = dir(fullfile(orient_dir, '*.jpg'));
                    image_files = [png_files; jpg_files];

                    count = length(image_files);
                    size_bytes = sum([image_files.bytes]);

                    stats.(orient) = struct('count', count, 'size_mb', size_bytes / (1024*1024));

                    total_images = total_images + count;
                    total_size = total_size + size_bytes;
                end
            end

            stats.total_images = total_images;
            stats.total_size_mb = total_size / (1024*1024);
        end
    end
end