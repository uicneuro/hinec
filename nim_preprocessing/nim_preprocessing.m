function nim_preprocessing(file_prefix, varargin)
% nim_preprocessing: Enhanced modular DWI preprocessing pipeline
%
% Arguments:
%   file_prefix - The prefix for the file paths
%   options - Structure with preprocessing options (optional):
%     .run_denoising - Boolean flag for denoising (default: true)
%     .denoise_method - 'dwidenoise', 'nlmeans', 'gaussian' (default: 'dwidenoise')
%     .run_motion_correction - Boolean flag for motion correction (default: true)
%     .run_eddy - Boolean flag for eddy current correction (default: true)
%     .improve_mask - Boolean flag for mask improvement (default: true)
%     .atlas_type - Atlas type (default: 'HarvardOxford')
%     .phase_encoding_direction - Phase encoding axis for eddy (e.g. 'j-', default: "")
%     .total_readout_time - Readout time in seconds for each acquisition (default: [])
%     .eddy_index_vector - Optional acquisition indices per volume (default: [])
%
% Legacy usage (backward compatibility):
%   nim_preprocessing(file_prefix, run_eddy, atlas_type)

fprintf('HINEC DWI Preprocessing Pipeline\n');
fprintf('================================\n');

% Ensure file_prefix is a character array (handle string input)
file_prefix = char(file_prefix);

% Handle backward compatibility and set default parameters
if length(varargin) == 0
    options = struct();
elseif length(varargin) == 1
    if isstruct(varargin{1})
        % New usage: nim_preprocessing(file_prefix, options)
        options = varargin{1};
    else
        % Legacy usage: nim_preprocessing(file_prefix, run_eddy)
        run_eddy = varargin{1};
        options = struct();
        options.run_eddy = run_eddy;
        options.atlas_type = 'HarvardOxford';
        options.run_denoising = true;
        options.run_motion_correction = true;
        options.improve_mask = true;
        options.denoise_method = 'dwidenoise';
    end
elseif length(varargin) == 2
    % Legacy usage: nim_preprocessing(file_prefix, run_eddy, atlas_type)
    run_eddy = varargin{1};
    atlas_type = varargin{2};
    options = struct();
    options.run_eddy = run_eddy;
    options.atlas_type = atlas_type;
    options.run_denoising = true;
    options.run_motion_correction = true;
    options.improve_mask = true;
    options.denoise_method = 'dwidenoise';
else
    error('Too many input arguments');
end

% Set default options
default_options = struct(...
    'run_denoising', true, ...
    'denoise_method', 'dwidenoise', ...
    'run_motion_correction', true, ...
    'run_fieldmap_correction', true, ...
    'fieldmap_file', '', ...
    'fieldmap_units', 'auto', ...
    'phase_encoding_dir', 'y', ...
    'dwell_time', 0.00058, ...
    'eddy_method', 'auto', ...
    'run_eddy', true, ...
    'create_wm_mask', true, ...
    'improve_mask', true, ...
    'atlas_type', 'HarvardOxford', ...
    'phase_encoding_direction', "", ...
    'total_readout_time', [], ...
    'eddy_index_vector', [] ...
);

% Merge user options with defaults
option_fields = fieldnames(default_options);
for i = 1:length(option_fields)
    field = option_fields{i};
    if ~isfield(options, field)
        options.(field) = default_options.(field);
    end
end

%% Auto-detect field map file if not provided
if options.run_fieldmap_correction && isempty(options.fieldmap_file)
    fprintf('Auto-detecting field map file...\n');

    % Get directory from file_prefix
    [file_dir, file_base, ~] = fileparts(file_prefix);
    if isempty(file_dir)
        file_dir = pwd;
    end

    % Look for field map files following {name}_ naming convention
    fmap_patterns = {
        fullfile(file_dir, [file_base '_fmap_Hz.nii.gz']),
        fullfile(file_dir, [file_base '_fmap_RadPerSec.nii.gz'])
    };

    for i = 1:length(fmap_patterns)
        if exist(fmap_patterns{i}, 'file')
            options.fieldmap_file = fmap_patterns{i};
            fprintf('Found field map: %s\n', options.fieldmap_file);
            break;
        end
    end

    if isempty(options.fieldmap_file)
        fprintf('⚠ No field map found. Field map correction will be skipped.\n');
        options.run_fieldmap_correction = false;
    end
end

% Display pipeline information
fprintf('File prefix: %s\n', file_prefix);
fprintf('Configuration:\n');
fprintf('  Denoising: %s (%s)\n', char(string(options.run_denoising)), options.denoise_method);
fprintf('  Motion correction: %s\n', char(string(options.run_motion_correction)));
fprintf('  Field map correction: %s\n', char(string(options.run_fieldmap_correction)));
if options.run_fieldmap_correction
    fprintf('    Field map file: %s\n', options.fieldmap_file);
    fprintf('    Phase encoding: %s\n', options.phase_encoding_dir);
end
fprintf('  Eddy correction: %s (%s)\n', char(string(options.run_eddy)), options.eddy_method);
fprintf('  Create WM mask: %s\n', char(string(options.create_wm_mask)));
fprintf('  Mask improvement: %s\n', char(string(options.improve_mask)));
fprintf('  Atlas type: %s\n', options.atlas_type);
fprintf('-------------------------------------------\n');

% Initialize preprocessing report
preprocessing_report = struct();
preprocessing_report.start_time = datetime('now');
preprocessing_report.file_prefix = char(file_prefix);
preprocessing_report.options = options;
preprocessing_report.steps_completed = {};
preprocessing_report.errors = {};
preprocessing_report.warnings = {};

% Source FSL environment
fprintf('Initializing FSL environment...\n');
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end
system('source $FSLDIR/etc/fslconf/fsl.sh');

% Define file paths
dwi_raw_file = [file_prefix '_raw.nii.gz'];
bvec_file = [file_prefix '.bvec'];
bval_file = [file_prefix '.bval'];
output_file = [file_prefix '.nii.gz'];
output_dir = fileparts(output_file);
if isempty(output_dir)
    output_dir = pwd;
end

% Verify input files exist
input_files = {dwi_raw_file, bvec_file, bval_file};
input_descriptions = {'raw DWI', 'b-vectors', 'b-values'};

fprintf('Verifying input files...\n');
for i = 1:length(input_files)
    fprintf('  Checking: %s\n', input_files{i});
    if ~isfile(input_files{i})
        error('Required input file not found: %s (%s)', char(input_files{i}), char(input_descriptions{i}));
    end
    file_info = dir(input_files{i});
    fprintf('  ✓ %s: %s (%.1f MB)\n', input_descriptions{i}, input_files{i}, file_info.bytes/1024/1024);
end

preprocessing_report.input_files = containers.Map(input_descriptions, input_files);
current_dwi_file = dwi_raw_file;
current_bvec_file = bvec_file;

try
    %% Step 1: Extract b0 volume
    fprintf('\n=== Step 1: B0 Extraction ===\n');
    step_start = tic;
    
    b0_file = preproc_extract_b0(current_dwi_file, output_dir);
    preprocessing_report.b0_file = b0_file;
    preprocessing_report.steps_completed{end+1} = 'b0_extraction';
    
    step_time = toc(step_start);
    fprintf('Step 1 completed in %.1f seconds\n', step_time);
    
    %% Step 2: Initial brain extraction
    fprintf('\n=== Step 2: Initial Brain Extraction ===\n');
    step_start = tic;
    
    initial_brain_mask_file = [file_prefix '_M_initial.nii.gz'];
    brain_mask_file = preproc_brain_extraction(b0_file, output_dir, initial_brain_mask_file);
    preprocessing_report.initial_brain_mask_file = brain_mask_file;
    preprocessing_report.steps_completed{end+1} = 'initial_brain_extraction';
    
    step_time = toc(step_start);
    fprintf('Step 2 completed in %.1f seconds\n', step_time);
    
    %% Step 3: Denoising (if enabled)
    if options.run_denoising
        fprintf('\n=== Step 3: Denoising ===\n');
        step_start = tic;
        
        denoised_file = preproc_denoising(current_dwi_file, file_prefix, options.denoise_method);
        current_dwi_file = denoised_file;
        preprocessing_report.denoised_file = denoised_file;
        preprocessing_report.steps_completed{end+1} = 'denoising';
        
        step_time = toc(step_start);
        fprintf('Step 3 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 3: Denoising (SKIPPED) ===\n');
    end

    %% Step 4: Field Map Distortion Correction (if enabled)
    if options.run_fieldmap_correction
        fprintf('\n=== Step 4: Field Map Distortion Correction ===\n');
        step_start = tic;

        try
            fieldmap_corrected_file = preproc_fieldmap_correction(current_dwi_file, ...
                options.fieldmap_file, file_prefix, ...
                'phase_dir', options.phase_encoding_dir, ...
                'dwell_time', options.dwell_time, ...
                'units', options.fieldmap_units, ...
                'mask_file', brain_mask_file);

            current_dwi_file = fieldmap_corrected_file;
            preprocessing_report.fieldmap_corrected_file = fieldmap_corrected_file;
            preprocessing_report.steps_completed{end+1} = 'fieldmap_correction';

            fprintf('✓ Field map correction successful\n');
        catch ME
            warning('HINEC:FieldMapFailed', 'Field map correction failed: %s', ME.message);
            preprocessing_report.warnings{end+1} = sprintf('Field map correction failed: %s', ME.message);
        end

        step_time = toc(step_start);
        fprintf('Step 4 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 4: Field Map Distortion Correction (SKIPPED) ===\n');
    end
    
    %% Step 5: Motion correction (if enabled)
    if options.run_motion_correction
        fprintf('\n=== Step 5: Motion Correction ===\n');
        step_start = tic;
        
        motion_corrected_file = preproc_motion_correction(current_dwi_file, current_bvec_file, bval_file, file_prefix);
        current_dwi_file = motion_corrected_file;
        
        % Update b-vectors to motion-corrected ones
        motion_corrected_bvec = strrep(current_bvec_file, '.bvec', '_motion_corrected.bvec');
        if isfile(motion_corrected_bvec)
            current_bvec_file = motion_corrected_bvec;
        end
        
        preprocessing_report.motion_corrected_file = motion_corrected_file;
        preprocessing_report.motion_corrected_bvec = current_bvec_file;
        preprocessing_report.steps_completed{end+1} = 'motion_correction';
        
        step_time = toc(step_start);
        fprintf('Step 5 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 5: Motion Correction (SKIPPED) ===\n');
    end

    %% Step 6: Enhanced Eddy current correction (if enabled)
    if options.run_eddy
        fprintf('\n=== Step 6: Enhanced Eddy Current Correction ===\n');
        step_start = tic;

        % Determine best eddy correction method
        if strcmp(options.eddy_method, 'auto')
            % Check if proper eddy files exist
            acqp_file = [file_prefix '_acqp.txt'];
            index_file = [file_prefix '_index.txt'];

            if exist(acqp_file, 'file') && exist(index_file, 'file')
                eddy_method = 'eddy';
            else
                eddy_method = 'eddy_correct';
            end
        else
            eddy_method = options.eddy_method;
            acqp_file = [file_prefix '_acqp.txt'];
            index_file = [file_prefix '_index.txt'];
        end

        fprintf('Using eddy correction method: %s\n', eddy_method);

        try
            if strcmp(eddy_method, 'eddy_correct')
                % Use simple eddy_correct as fallback
                fprintf('Applying basic eddy current correction...\n');
                eddy_corrected_file = [file_prefix '_eddy_corrected.nii.gz'];

                % Set PATH to include FSL bin directory for eddy_correct script
                old_path = getenv('PATH');
                new_path = sprintf('%s/bin:%s', fsl_path, old_path);
                setenv('PATH', new_path);

                cmd = sprintf('%s/bin/eddy_correct %s %s 0', fsl_path, current_dwi_file, eddy_corrected_file);
                [status, result] = system(cmd);

                % Restore original PATH
                setenv('PATH', old_path);

                if status == 0 && exist(eddy_corrected_file, 'file')
                    current_dwi_file = eddy_corrected_file;
                    preprocessing_report.eddy_corrected_file = eddy_corrected_file;
                    preprocessing_report.eddy_method_used = 'eddy_correct';
                    preprocessing_report.eddy_corrected = true;
                    preprocessing_report.steps_completed{end+1} = 'eddy_correction';
                    fprintf('✓ Basic eddy correction successful\n');
                else
                    error('eddy_correct failed: %s', result);
                end
            else
                % Use standard eddy with parameter files
                % Create parameter files from preprocessing options when needed
                missing_params_warning = '';
                if ~isfile(acqp_file) || ~isfile(index_file)
                    fprintf('Creating eddy parameter files from preprocessing options...\n');
                    [created, creation_message] = create_eddy_parameter_files(acqp_file, index_file, current_dwi_file, options);
                    if created
                        fprintf('  ✓ %s\n', creation_message);
                    else
                        missing_params_warning = creation_message;
                    end
                else
                    fprintf('Using existing eddy parameter files\n');
                end

                % Check if we can run advanced eddy
                can_run_eddy = isfile(acqp_file) && isfile(index_file) && isempty(missing_params_warning);

                if can_run_eddy
                    eddy_corrected_file = preproc_eddy_correction(current_dwi_file, brain_mask_file, ...
                        current_bvec_file, bval_file, acqp_file, index_file, file_prefix);

                    if ~isempty(eddy_corrected_file) && exist(eddy_corrected_file, 'file')
                        current_dwi_file = eddy_corrected_file;

                        % Update b-vectors to eddy-corrected ones
                        eddy_corrected_bvec = strrep(eddy_corrected_file, '.nii.gz', '.eddy_rotated_bvecs');
                        if exist(eddy_corrected_bvec, 'file')
                            current_bvec_file = eddy_corrected_bvec;
                        end

                        preprocessing_report.eddy_corrected_file = eddy_corrected_file;
                        preprocessing_report.eddy_corrected_bvec = current_bvec_file;
                        preprocessing_report.eddy_method_used = 'eddy';
                        preprocessing_report.eddy_corrected = true;
                        preprocessing_report.steps_completed{end+1} = 'eddy_correction';
                        fprintf('✓ Advanced eddy correction successful\n');
                    else
                        error('Advanced eddy correction failed');
                    end
                else
                    if isempty(missing_params_warning)
                        missing_params_warning = 'Eddy correction skipped: parameter files not available';
                    end
                    fprintf('⚠ %s\n', missing_params_warning);
                    preprocessing_report.warnings{end+1} = missing_params_warning;
                    preprocessing_report.eddy_corrected = false;
                end
            end

        catch ME
            warning('HINEC:EddyFailed', 'Eddy correction failed: %s', ME.message);
            preprocessing_report.warnings{end+1} = sprintf('Eddy correction failed: %s', ME.message);
            preprocessing_report.eddy_corrected = false;
        end
        
        step_time = toc(step_start);
        fprintf('Step 6 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 6: Eddy Current Correction (SKIPPED) ===\n');
        preprocessing_report.eddy_corrected = false;
    end

    %% Step 7: Enhanced Brain Mask and White Matter Segmentation
    if options.create_wm_mask
        fprintf('\n=== Step 7: Enhanced Brain Mask and White Matter Segmentation ===\n');
        step_start = tic;

        fprintf('Creating white matter mask for improved tractography...\n');
        wm_mask_file = [file_prefix '_WM_mask.nii.gz'];

        % Create a preliminary white matter mask using simple thresholding
        % This will be refined after DTI calculation in main.m
        temp_fa_file = [file_prefix '_temp_fa.nii.gz'];

        % Run quick DTI to get preliminary FA
        fprintf('Computing preliminary FA for white matter segmentation...\n');
        try
            % Basic DTI calculation for FA estimation
            final_mask_file = [file_prefix '_M.nii.gz'];
            copyfile(brain_mask_file, final_mask_file);

            cmd = sprintf('%s/bin/dtifit -k %s -o %s_temp_dti -m %s -r %s -b %s', ...
                fsl_path, current_dwi_file, file_prefix, final_mask_file, current_bvec_file, bval_file);
            [status, result] = system(cmd);

            if status == 0
                temp_fa_file = [file_prefix '_temp_dti_FA.nii.gz'];
                if exist(temp_fa_file, 'file')
                    % Create white matter mask (FA > 0.2, eroded to avoid boundaries)
                    cmd = sprintf('%s/bin/fslmaths %s -thr 0.2 -bin %s', fsl_path, temp_fa_file, wm_mask_file);
                    system(cmd);

                    % Erode to remove boundary voxels
                    cmd = sprintf('%s/bin/fslmaths %s -ero -ero %s', fsl_path, wm_mask_file, wm_mask_file);
                    system(cmd);

                    preprocessing_report.wm_mask_file = wm_mask_file;
                    preprocessing_report.steps_completed{end+1} = 'white_matter_segmentation';
                    fprintf('✓ White matter mask created for enhanced tractography\n');

                    % Clean up temporary DTI files
                    temp_files = dir([file_prefix '_temp_dti*']);
                    for i = 1:length(temp_files)
                        delete(fullfile(temp_files(i).folder, temp_files(i).name));
                    end
                end
            end
        catch
            fprintf('⚠ Could not create white matter mask, will use brain mask\n');
            preprocessing_report.warnings{end+1} = 'White matter mask creation failed';
        end

        step_time = toc(step_start);
        fprintf('Step 7 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 7: Enhanced Segmentation (SKIPPED) ===\n');
        final_mask_file = [file_prefix '_M.nii.gz'];
        copyfile(brain_mask_file, final_mask_file);
    end

    %% Step 8: Copy processed data to final location
    fprintf('\n=== Step 8: Copy Processed Data to Final Location ===\n');
    step_start = tic;
    
    % Define final output bvec file path
    final_bvec_file = [file_prefix '.bvec'];

    % Copy current processed DWI to final location
    copyfile(current_dwi_file, output_file);

    % Copy final b-vectors
    copyfile(current_bvec_file, final_bvec_file);
    
    preprocessing_report.steps_completed{end+1} = 'copy_final_data';
    
    step_time = toc(step_start);
    fprintf('Step 8 completed in %.1f seconds\n', step_time);

    %% Step 9: Atlas processing
    fprintf('\n=== Step 9: Atlas Processing ===\n');
    step_start = tic;
    
    [parcellation_mask_output, atlas_labels_file] = preproc_atlas_resampling(...
        output_file, output_dir, file_prefix, options.atlas_type);
    
    preprocessing_report.parcellation_mask = parcellation_mask_output;
    preprocessing_report.atlas_labels_file = atlas_labels_file;
    preprocessing_report.steps_completed{end+1} = 'atlas_processing';
    
    step_time = toc(step_start);
    fprintf('Step 9 completed in %.1f seconds\n', step_time);

    %% Step 10: Finalization
    fprintf('\n=== Step 10: Finalization ===\n');
    step_start = tic;
    
    % Copy files to final locations with standard names
final_output_bval = [file_prefix '.bval'];   % Copy of b-values
final_output_mask = [file_prefix '_M.nii.gz'];  % Final brain mask

% B-vectors are already at final_bvec_file, which is the final output location

% Copy b-values
if ~strcmp(bval_file, final_output_bval)
    copyfile(bval_file, final_output_bval);
end

% Copy final brain mask
if ~strcmp(brain_mask_file, final_output_mask)
    copyfile(brain_mask_file, final_output_mask);
end
    
    % Define files to keep for cleanup
    final_files = {
        output_file;                               % Final processed DWI
        final_bvec_file;                         % Final b-vectors
        final_output_bval;                         % B-values
        final_output_mask;                         % Final brain mask
        parcellation_mask_output;                  % Parcellation
        atlas_labels_file                         % Atlas labels
    };

    % Add enhanced files if they exist
    if isfield(preprocessing_report, 'wm_mask_file') && exist(preprocessing_report.wm_mask_file, 'file')
        final_files{end+1} = preprocessing_report.wm_mask_file;
    end
    
    cleanup_report = preproc_cleanup(output_dir, file_prefix, final_files);
    preprocessing_report.cleanup_report = cleanup_report;
    preprocessing_report.steps_completed{end+1} = 'finalization';
    
    % Update final file assignments in report
preprocessing_report.final_dwi_file = output_file;
preprocessing_report.final_bvec_file = final_bvec_file;
preprocessing_report.final_bval_file = final_output_bval;
    preprocessing_report.final_mask_file = final_output_mask;
    
    step_time = toc(step_start);
    fprintf('Step 10 completed in %.1f seconds\n', step_time);

    % Finalize report
    preprocessing_report.end_time = datetime('now');
    preprocessing_report.total_duration = preprocessing_report.end_time - preprocessing_report.start_time;
    preprocessing_report.success = true;

    % Store report for potential future use
    report_file = [file_prefix '_preprocessing_report.mat'];
    save(report_file, 'preprocessing_report');

    fprintf('\n========================================\n');
    fprintf('🎉 ENHANCED PREPROCESSING COMPLETE 🎉\n');
    fprintf('========================================\n');
    fprintf('✅ Enhanced preprocessing successful\n');
    fprintf('⏱ Total processing time: %s\n', char(preprocessing_report.total_duration));
    fprintf('📁 Output directory: %s\n', output_dir);

    % Display enhanced features applied
    enhanced_features = {};
    if options.run_fieldmap_correction && isfield(preprocessing_report, 'fieldmap_corrected_file')
        enhanced_features{end+1} = 'Field map distortion correction';
    end
    if isfield(preprocessing_report, 'eddy_method_used') && strcmp(preprocessing_report.eddy_method_used, 'eddy_correct')
        enhanced_features{end+1} = 'Alternative eddy current correction';
    end
    if isfield(preprocessing_report, 'wm_mask_file')
        enhanced_features{end+1} = 'White matter mask for better tractography';
    end

    if ~isempty(enhanced_features)
        fprintf('\n🚀 ENHANCED FEATURES APPLIED:\n');
        for i = 1:length(enhanced_features)
            fprintf('  ✓ %s\n', enhanced_features{i});
        end
    end

    fprintf('\n🧠 Atlas used: %s\n', options.atlas_type);
    fprintf('📊 Report saved to: %s\n', report_file);
    
    % Display processing summary
    fprintf('\nPROCESSING SUMMARY:\n');
    for i = 1:length(preprocessing_report.steps_completed)
        fprintf('  ✓ %s\n', preprocessing_report.steps_completed{i});
    end
    
    if ~isempty(preprocessing_report.warnings)
        fprintf('\nWARNINGS:\n');
        for i = 1:length(preprocessing_report.warnings)
            fprintf('  ⚠ %s\n', preprocessing_report.warnings{i});
        end
    end
    
    % Display final files
    fprintf('\nFINAL OUTPUT FILES:\n');
    final_outputs = {
    output_file, 'Processed DWI data';
    final_bvec_file, 'Final b-vectors';
    final_output_bval, 'B-values';
    final_output_mask, 'Brain mask';
    parcellation_mask_output, 'Parcellation mask';
    atlas_labels_file, 'Atlas labels'
};
    
    for i = 1:size(final_outputs, 1)
        file_path = final_outputs{i, 1};
        description = final_outputs{i, 2};
        if isfile(file_path)
            file_info = dir(file_path);
            fprintf('  ✓ %s: %s (%.1f MB)\n', description, file_path, file_info.bytes/1024/1024);
        else
            fprintf('  ✗ %s: %s (MISSING)\n', description, file_path);
        end
    end
    
    fprintf('\n🚀 Your data is now ready for DTI analysis and tractography! 🚀\n');
    fprintf('Next steps:\n');
    fprintf('  1. Run main(''%s'', ''output.mat'') for full DTI analysis\n', char(output_file));
    fprintf('  2. Use runTractography(''output.mat'') for tractography analysis\n');
    
catch ME
    % Record comprehensive error information
    error_info = struct();
    error_info.message = ME.message;
    error_info.identifier = ME.identifier;
    error_info.stack = ME.stack;
    error_info.timestamp = datetime('now');
    
    preprocessing_report.errors{end+1} = error_info;
    preprocessing_report.success = false;
    preprocessing_report.end_time = datetime('now');
    
    % Save error report
    error_report_file = [file_prefix '_preprocessing_error_report.mat'];
    save(error_report_file, 'preprocessing_report');
    
    fprintf('\n❌ ENHANCED PREPROCESSING FAILED ❌\n');
    fprintf('Error in %s (line %d): %s\n', ME.stack(1).name, ME.stack(1).line, ME.message);
    fprintf('Error report saved to: %s\n', error_report_file);
    
    if ~isempty(preprocessing_report.steps_completed)
        fprintf('\nSteps completed before failure:\n');
        for i = 1:length(preprocessing_report.steps_completed)
            fprintf('  ✓ %s\n', preprocessing_report.steps_completed{i});
        end
    end
    
    % Re-throw the error
    rethrow(ME);
end

end

function [success, message] = create_eddy_parameter_files(acqp_file, index_file, current_dwi_file, options)
% Helper: generate eddy parameter files from preprocessing options
success = false;
message = '';

phase_info = options.phase_encoding_direction;
readout_info = options.total_readout_time;
index_vector = options.eddy_index_vector;

% Normalise phase encoding specification to a cell array of strings
if isstring(phase_info)
    phase_list = cellstr(phase_info(:));
elseif ischar(phase_info)
    phase_list = {phase_info};
elseif iscell(phase_info)
    phase_list = phase_info(:);
else
    phase_list = {};
end
phase_list = cellfun(@(c) strtrim(char(c)), phase_list, 'UniformOutput', false);
phase_list = phase_list(~cellfun(@isempty, phase_list));

% Normalise readout times to a numeric row vector
if isempty(readout_info)
    readout_values = [];
elseif iscell(readout_info)
    readout_values = cellfun(@double, readout_info);
else
    readout_values = double(readout_info);
end
readout_values = readout_values(:)';

if isempty(phase_list) || isempty(readout_values)
    message = 'Provide options.phase_encoding_direction and options.total_readout_time to generate eddy parameter files.';
    return;
end

if numel(readout_values) == 1 && numel(phase_list) > 1
    readout_values = repmat(readout_values, 1, numel(phase_list));
end

if numel(readout_values) ~= numel(phase_list)
    message = 'Number of entries in total_readout_time must match phase_encoding_direction.';
    return;
end

if any(~isfinite(readout_values) | readout_values <= 0)
    message = 'total_readout_time values must be positive, finite numbers (seconds).';
    return;
end

% Write acqp file with one line per acquisition
try
    fid = fopen(acqp_file, 'w');
    if fid == -1
        message = sprintf('Could not create acqp file: %s', acqp_file);
        return;
    end
    cleanup_acqp = onCleanup(@() fclose(fid));
    for i = 1:numel(phase_list)
        pe_vector = phase_encoding_direction_to_vector(phase_list{i});
        fprintf(fid, '%s %.6f\n', pe_vector, readout_values(i));
    end
    clear cleanup_acqp; % closes file
catch ME
    message = sprintf('Failed to write acqp file: %s', ME.message);
    return;
end

% Count diffusion volumes to build index file
[status, vol_output] = system(sprintf('fslnvols %s', current_dwi_file));
if status ~= 0
    message = sprintf('Could not count volumes in %s: %s', current_dwi_file, strtrim(vol_output));
    return;
end

num_vols = str2double(strtrim(vol_output));
if isnan(num_vols) || num_vols <= 0
    message = sprintf('Invalid volume count returned by fslnvols: %s', vol_output);
    return;
end

if isempty(index_vector)
    index_vector = ones(1, num_vols);
else
    index_vector = double(index_vector(:)');
    if numel(index_vector) ~= num_vols
        message = sprintf('Provided eddy_index_vector has %d entries but dataset has %d volumes.', numel(index_vector), num_vols);
        return;
    end
    if any(~isfinite(index_vector)) || any(index_vector < 1) || any(abs(index_vector - round(index_vector)) > 0)
        message = 'eddy_index_vector must contain positive integer acquisition indices.';
        return;
    end
end

if max(index_vector) > numel(phase_list)
    message = sprintf('eddy_index_vector references acquisition %d, but only %d acquisition lines provided.', max(index_vector), numel(phase_list));
    return;
end

try
    fid = fopen(index_file, 'w');
    if fid == -1
        message = sprintf('Could not create index file: %s', index_file);
        return;
    end
    cleanup_index = onCleanup(@() fclose(fid));
    fprintf(fid, '%d ', index_vector);
    fprintf(fid, '\n');
    clear cleanup_index;
catch ME
    message = sprintf('Failed to write index file: %s', ME.message);
    return;
end

success = true;
message = 'Created acqp/index parameter files using preprocessing options.';
end

function pe_vector = phase_encoding_direction_to_vector(direction)
% Convert a BIDS-style phase encoding direction (e.g., 'j-') to FSL vector
direction = lower(strtrim(direction));
if isempty(direction)
    error('Phase encoding direction cannot be empty.');
end

sign = 1;
if direction(end) == '-'
    sign = -1;
    direction = direction(1:end-1);
elseif direction(end) == '+'
    direction = direction(1:end-1);
end

switch direction
    case {'i', 'x'}
        pe_vector = sprintf('%d 0 0', sign);
    case {'j', 'y'}
        pe_vector = sprintf('0 %d 0', sign);
    case {'k', 'z'}
        pe_vector = sprintf('0 0 %d', sign);
    otherwise
        error('Unsupported phase encoding axis: %s', direction);
end
end
