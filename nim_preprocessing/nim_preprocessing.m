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
    'run_eddy', true, ...
    'improve_mask', true, ...
    'atlas_type', 'HarvardOxford' ...
);

% Merge user options with defaults
option_fields = fieldnames(default_options);
for i = 1:length(option_fields)
    field = option_fields{i};
    if ~isfield(options, field)
        options.(field) = default_options.(field);
    end
end

% Display pipeline information
fprintf('File prefix: %s\n', file_prefix);
fprintf('Configuration:\n');
fprintf('  Denoising: %s (%s)\n', char(string(options.run_denoising)), options.denoise_method);
fprintf('  Motion correction: %s\n', char(string(options.run_motion_correction)));
fprintf('  Eddy correction: %s\n', char(string(options.run_eddy)));
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
    
    %% Step 4: Motion correction (if enabled)
    if options.run_motion_correction
        fprintf('\n=== Step 4: Motion Correction ===\n');
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
        fprintf('Step 4 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 4: Motion Correction (SKIPPED) ===\n');
    end
    
    %% Step 5: Eddy current correction (if enabled)
    if options.run_eddy
        fprintf('\n=== Step 5: Eddy Current Correction ===\n');
        step_start = tic;
        
        % Define eddy parameter files
        acqp_file = [file_prefix '_acqp.txt'];
        index_file = [file_prefix '_index.txt'];
        json_file = [file_prefix '.json'];
        
        % Create parameter files from JSON if they don't exist
        if ~isfile(acqp_file) || ~isfile(index_file)
            fprintf('Creating eddy parameter files from JSON metadata...\n');
            
            if isfile(json_file)
                try
                    % Read JSON file
                    json_text = fileread(json_file);
                    json_data = jsondecode(json_text);
                    
                    % Extract phase encoding direction and readout time
                    if isfield(json_data, 'PhaseEncodingDirection') && isfield(json_data, 'TotalReadoutTime')
                        phase_dir = json_data.PhaseEncodingDirection;
                        readout_time = json_data.TotalReadoutTime;
                        
                        fprintf('  Phase encoding direction: %s\n', phase_dir);
                        fprintf('  Total readout time: %.6f seconds\n', readout_time);
                        
                        % Convert phase encoding direction to FSL format
                        % i = x-direction, j = y-direction, k = z-direction
                        % + = positive direction, - = negative direction
                        switch phase_dir
                            case 'i'
                                pe_vector = '1 0 0';
                            case 'i-'
                                pe_vector = '-1 0 0';
                            case 'j'
                                pe_vector = '0 1 0';
                            case 'j-'
                                pe_vector = '0 -1 0';
                            case 'k'
                                pe_vector = '0 0 1';
                            case 'k-'
                                pe_vector = '0 0 -1';
                            otherwise
                                error('Unknown phase encoding direction: %s', phase_dir);
                        end
                        
                        % Create acqp.txt file
                        acqp_content = sprintf('%s %.6f', pe_vector, readout_time);
                        fid = fopen(acqp_file, 'w');
                        if fid == -1
                            error('Could not create acqp file: %s', acqp_file);
                        end
                        fprintf(fid, '%s\n', acqp_content);
                        fclose(fid);
                        fprintf('  ✓ Created %s\n', acqp_file);
                        
                        % Count volumes and create index.txt
                        [status, vol_output] = system(sprintf('fslnvols %s', current_dwi_file));
                        if status == 0
                            num_vols = str2double(strtrim(vol_output));
                            if isnan(num_vols) || num_vols <= 0
                                error('Invalid volume count: %s', vol_output);
                            end
                            
                            % Create index file (all volumes use acquisition 1)
                            index_content = repmat('1 ', 1, num_vols);
                            fid = fopen(index_file, 'w');
                            if fid == -1
                                error('Could not create index file: %s', index_file);
                            end
                            fprintf(fid, '%s\n', strtrim(index_content));
                            fclose(fid);
                            fprintf('  ✓ Created %s for %d volumes\n', index_file, num_vols);
                        else
                            error('Could not count volumes in %s: %s', current_dwi_file, vol_output);
                        end
                        
                    else
                        error('JSON file missing required fields: PhaseEncodingDirection and/or TotalReadoutTime');
                    end
                    
                catch ME
                    fprintf('⚠ Failed to create parameter files from JSON: %s\n', ME.message);
                    fprintf('  You may need to create %s and %s manually\n', acqp_file, index_file);
                    preprocessing_report.warnings{end+1} = sprintf('Failed to create eddy parameter files: %s', ME.message);
                    preprocessing_report.eddy_corrected = false;
                    step_time = toc(step_start);
                    fprintf('Step 5 completed in %.1f seconds\n', step_time);
                    return;
                end
            else
                fprintf('⚠ JSON file not found: %s\n', json_file);
                fprintf('  Cannot create eddy parameter files automatically\n');
                preprocessing_report.warnings{end+1} = 'JSON file not found for eddy parameter creation';
                preprocessing_report.eddy_corrected = false;
                step_time = toc(step_start);
                fprintf('Step 5 completed in %.1f seconds\n', step_time);
                return;
            end
        else
            fprintf('Using existing eddy parameter files\n');
        end
        
        % Now run eddy correction with parameter files
        if isfile(acqp_file) && isfile(index_file)
            try
                eddy_corrected_file = preproc_eddy_correction(current_dwi_file, brain_mask_file, ...
                    current_bvec_file, bval_file, acqp_file, index_file, file_prefix);
                
                if ~isempty(eddy_corrected_file) && isfile(eddy_corrected_file)
                    current_dwi_file = eddy_corrected_file;
                    
                    % Update b-vectors to eddy-corrected ones
                    eddy_corrected_bvec = strrep(eddy_corrected_file, '.nii.gz', '.eddy_rotated_bvecs');
                    if isfile(eddy_corrected_bvec)
                        current_bvec_file = eddy_corrected_bvec;
                    end
                    
                    preprocessing_report.eddy_corrected_file = eddy_corrected_file;
                    preprocessing_report.eddy_corrected_bvec = current_bvec_file;
                    preprocessing_report.eddy_corrected = true;
                    preprocessing_report.steps_completed{end+1} = 'eddy_correction';
                else
                    preprocessing_report.warnings{end+1} = 'Eddy correction failed, continuing without it';
                    preprocessing_report.eddy_corrected = false;
                end
            catch ME
                preprocessing_report.warnings{end+1} = sprintf('Eddy correction error: %s', ME.message);
                preprocessing_report.eddy_corrected = false;
                fprintf('⚠ Eddy correction failed: %s\n', ME.message);
            end
        else
            fprintf('⚠ Could not create or find eddy parameter files, skipping eddy correction\n');
            preprocessing_report.warnings{end+1} = 'Could not create eddy parameter files';
            preprocessing_report.eddy_corrected = false;
        end
        
        step_time = toc(step_start);
        fprintf('Step 5 completed in %.1f seconds\n', step_time);
    else
        fprintf('\n=== Step 5: Eddy Current Correction (SKIPPED) ===\n');
        preprocessing_report.eddy_corrected = false;
    end
    
    %% Step 6: Copy processed data to final location
    fprintf('\n=== Step 6: Copy Processed Data to Final Location ===\n');
    step_start = tic;
    
    % Define final output bvec file path
    final_bvec_file = [file_prefix '.bvec'];

    % Copy current processed DWI to final location
    copyfile(current_dwi_file, output_file);

    % Copy final b-vectors
    copyfile(current_bvec_file, final_bvec_file);
    
    preprocessing_report.steps_completed{end+1} = 'copy_final_data';
    
    step_time = toc(step_start);
    fprintf('Step 6 completed in %.1f seconds\n', step_time);
    
    %% Step 7: Copy brain mask to final location
    fprintf('\n=== Step 7: Copy Brain Mask to Final Location ===\n');
    step_start = tic;
    
    % Copy initial mask to final location (improvement will happen in main.m with FA data)
    final_mask_file = [file_prefix '_M.nii.gz'];
    copyfile(brain_mask_file, final_mask_file);
    brain_mask_file = final_mask_file;
    
    preprocessing_report.steps_completed{end+1} = 'copy_brain_mask';
    
    step_time = toc(step_start);
    fprintf('Step 7 completed in %.1f seconds\n', step_time);
    
    %% Step 8: Atlas processing
    fprintf('\n=== Step 8: Atlas Processing ===\n');
    step_start = tic;
    
    [parcellation_mask_output, atlas_labels_file] = preproc_atlas_resampling(...
        output_file, output_dir, file_prefix, options.atlas_type);
    
    preprocessing_report.parcellation_mask = parcellation_mask_output;
    preprocessing_report.atlas_labels_file = atlas_labels_file;
    preprocessing_report.steps_completed{end+1} = 'atlas_processing';
    
    step_time = toc(step_start);
    fprintf('Step 8 completed in %.1f seconds\n', step_time);
    
    %% Step 9: Finalization
    fprintf('\n=== Step 9: Finalization ===\n');
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
    
    cleanup_report = preproc_cleanup(output_dir, file_prefix, final_files);
    preprocessing_report.cleanup_report = cleanup_report;
    preprocessing_report.steps_completed{end+1} = 'finalization';
    
    % Update final file assignments in report
preprocessing_report.final_dwi_file = output_file;
preprocessing_report.final_bvec_file = final_bvec_file;
preprocessing_report.final_bval_file = final_output_bval;
    preprocessing_report.final_mask_file = final_output_mask;
    
    step_time = toc(step_start);
    fprintf('Step 9 completed in %.1f seconds\n', step_time);
    
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
    fprintf('✅ All requested steps completed successfully\n');
    fprintf('⏱ Total processing time: %s\n', char(preprocessing_report.total_duration));
    fprintf('📁 Output directory: %s\n', output_dir);
    fprintf('🧠 Atlas used: %s\n', options.atlas_type);
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