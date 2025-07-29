function nim_preprocessing(file_prefix, run_eddy, atlas_type)
% nim_preprocessing: Modular DWI preprocessing pipeline
%
% Arguments:
%   file_prefix - The prefix for the file paths
%   run_eddy - Boolean flag for eddy current correction (optional, default: false)
%   atlas_type - Atlas type ('HarvardOxford', 'JHU', 'JHU-tract', default: 'HarvardOxford')

fprintf('HINEC DWI Preprocessing Pipeline\n');
fprintf('================================\n');

% Set default parameters
if nargin < 3
    atlas_type = 'HarvardOxford';
end
if nargin < 2
    run_eddy = false;
end

% Display pipeline information
fprintf('File prefix: %s\n', file_prefix);
fprintf('Atlas type: %s\n', atlas_type);
fprintf('Eddy correction: %s\n', char(string(run_eddy)));
fprintf('--------------------------------\n');

% Initialize preprocessing report
preprocessing_report = struct();
preprocessing_report.start_time = datetime('now');
preprocessing_report.file_prefix = char(file_prefix);
preprocessing_report.atlas_type = atlas_type;
preprocessing_report.eddy_enabled = run_eddy;
preprocessing_report.steps_completed = {};
preprocessing_report.errors = {};

% Source FSL environment
fprintf('Initializing FSL environment...\n');
fsl_path = getenv('FSLDIR');
if isempty(fsl_path)
    error('FSLDIR environment variable is not set. Please ensure FSL is installed and configured.');
end
system('source $FSLDIR/etc/fslconf/fsl.sh');

% Define file paths
dwi_file = file_prefix + "_raw.nii.gz";
bvec_file = file_prefix + ".bvec";
bval_file = file_prefix + ".bval";
brain_mask_file = file_prefix + "_M.nii.gz";
output_file = file_prefix + ".nii.gz";
output_dir = fileparts(output_file);

% Verify input files exist
input_files = {dwi_file, bvec_file, bval_file};
input_descriptions = {'raw DWI', 'b-vectors', 'b-values'};

fprintf('Verifying input files...\n');
for i = 1:length(input_files)
    if ~isfile(input_files{i})
        error('Required input file not found: %s (%s)', input_files{i}, input_descriptions{i});
    end
    file_info = dir(input_files{i});
    fprintf('  ✓ %s: %s (%.1f MB)\n', input_descriptions{i}, input_files{i}, file_info.bytes/1024/1024);
end

preprocessing_report.input_files = containers.Map(input_descriptions, input_files);

try
    %% Step 1: Extract b0 volume
    fprintf('\n--- Step 1: B0 Extraction ---\n');
    step_start = tic;
    
    b0_file = preproc_extract_b0(dwi_file, output_dir);
    preprocessing_report.b0_file = b0_file;
    preprocessing_report.steps_completed{end+1} = 'b0_extraction';
    
    step_time = toc(step_start);
    fprintf('Step 1 completed in %.1f seconds\n', step_time);
    
    %% Step 2: Brain extraction
    fprintf('\n--- Step 2: Brain Extraction ---\n');
    step_start = tic;
    
    brain_mask_file = preproc_brain_extraction(b0_file, output_dir, brain_mask_file);
    preprocessing_report.brain_mask_file = brain_mask_file;
    preprocessing_report.steps_completed{end+1} = 'brain_extraction';
    
    step_time = toc(step_start);
    fprintf('Step 2 completed in %.1f seconds\n', step_time);
    
    %% Step 3: Copy raw data as processed (or eddy correction if enabled)
    fprintf('\n--- Step 3: DWI Processing ---\n');
    step_start = tic;
    
    if run_eddy
        fprintf('Eddy current correction requested...\n');
        
        % Check for required eddy files
        acqp_file = file_prefix + "_acqp.txt";
        index_file = file_prefix + "_index.txt";
        
        if isfile(acqp_file) && isfile(index_file)
            eddy_corrected_file = preproc_eddy_correction(dwi_file, brain_mask_file, ...
                bvec_file, bval_file, acqp_file, index_file, file_prefix);
            
            if ~isempty(eddy_corrected_file) && isfile(eddy_corrected_file)
                copyfile(eddy_corrected_file, output_file);
                preprocessing_report.eddy_corrected = true;
                preprocessing_report.eddy_file = eddy_corrected_file;
            else
                fprintf('Eddy correction failed, using raw data\n');
                copyfile(dwi_file, output_file);
                preprocessing_report.eddy_corrected = false;
            end
        else
            fprintf('Missing eddy files (%s, %s), using raw data\n', acqp_file, index_file);
            copyfile(dwi_file, output_file);
            preprocessing_report.eddy_corrected = false;
        end
    else
        fprintf('Using raw DWI data (no eddy correction)\n');
        copyfile(dwi_file, output_file);
        preprocessing_report.eddy_corrected = false;
    end
    
    preprocessing_report.output_file = output_file;
    preprocessing_report.steps_completed{end+1} = 'dwi_processing';
    
    step_time = toc(step_start);
    fprintf('Step 3 completed in %.1f seconds\n', step_time);
    
    %% Step 4: Atlas resampling and label loading
    fprintf('\n--- Step 4: Atlas Processing ---\n');
    step_start = tic;
    
    [parcellation_mask_output, atlas_labels_file] = preproc_atlas_resampling_fixed(...
        output_file, output_dir, file_prefix, atlas_type);
    
    preprocessing_report.parcellation_mask = parcellation_mask_output;
    preprocessing_report.atlas_labels_file = atlas_labels_file;
    preprocessing_report.steps_completed{end+1} = 'atlas_processing';
    
    step_time = toc(step_start);
    fprintf('Step 4 completed in %.1f seconds\n', step_time);
    
    %% Step 5: Cleanup intermediate files
    fprintf('\n--- Step 5: Cleanup ---\n');
    step_start = tic;
    
    final_files = {
        output_file;                                       % Processed DWI
        brain_mask_file;                                   % Brain mask
        parcellation_mask_output;                          % Parcellation
        atlas_labels_file                                  % Atlas labels
    };
    
    cleanup_report = preproc_cleanup(output_dir, file_prefix, final_files);
    preprocessing_report.cleanup_report = cleanup_report;
    preprocessing_report.steps_completed{end+1} = 'cleanup';
    
    step_time = toc(step_start);
    fprintf('Step 5 completed in %.1f seconds\n', step_time);
    
    % Finalize report
    preprocessing_report.end_time = datetime('now');
    preprocessing_report.total_duration = preprocessing_report.end_time - preprocessing_report.start_time;
    preprocessing_report.success = true;
    
    % Store report for potential future use
    report_file = file_prefix + "_preprocessing_report.mat";
    save(report_file, 'preprocessing_report');
    
    fprintf('\n=== PREPROCESSING COMPLETE ===\n');
    fprintf('✅ All steps completed successfully\n');
    fprintf('⏱ Total processing time: %s\n', char(preprocessing_report.total_duration));
    fprintf('📁 Output directory: %s\n', output_dir);
    fprintf('🧠 Atlas used: %s\n', atlas_type);
    fprintf('📊 Report saved to: %s\n', report_file);
    
    % Display final summary
    fprintf('\nFINAL FILES:\n');
    if isfield(preprocessing_report, 'cleanup_report') && isfield(preprocessing_report.cleanup_report, 'files_kept')
        for i = 1:length(preprocessing_report.cleanup_report.files_kept)
            file_record = preprocessing_report.cleanup_report.files_kept{i};
            if file_record.exists
                fprintf('  ✓ %s (%.1f MB)\n', file_record.path, file_record.size_mb);
            else
                fprintf('  ✗ %s (MISSING)\n', file_record.path);
            end
        end
    end
    
    fprintf('\nPreprocessing pipeline completed successfully! 🎉\n');
    
catch ME
    % Record error information
    error_info = struct();
    error_info.message = ME.message;
    error_info.identifier = ME.identifier;
    error_info.stack = ME.stack;
    error_info.timestamp = datetime('now');
    
    preprocessing_report.errors{end+1} = error_info;
    preprocessing_report.success = false;
    
    fprintf('\n❌ PREPROCESSING PIPELINE FAILED ❌\n');
    fprintf('Error in %s (line %d): %s\n', ME.stack(1).name, ME.stack(1).line, ME.message);
    
    % Re-throw the error
    rethrow(ME);
end

end