function test_enhanced_preprocessing()
% TEST_ENHANCED_PREPROCESSING: Test script for field map-based preprocessing
%
% This script demonstrates how to use the enhanced preprocessing pipeline
% with field map correction for the ISMRM dataset.

fprintf('=== Testing Enhanced Preprocessing Pipeline ===\n');

%% Setup paths
addpath('nim_preprocessing');
addpath('nim_utils');

%% Test with ISMRM data
data_prefix = 'ISMRM/ISMRM';  % Adjust path as needed

% Check if required files exist
required_files = {
    [data_prefix '_raw.nii.gz']; ...
    [data_prefix '.bvec']; ...
    [data_prefix '.bval']; ...
    [data_prefix '_fmap_Hz.nii.gz']  % Field map file
};

fprintf('Checking required files...\n');
all_files_exist = true;
for i = 1:length(required_files)
    if exist(required_files{i}, 'file')
        fprintf('  ✓ %s\n', required_files{i});
    else
        fprintf('  ✗ %s (MISSING)\n', required_files{i});
        all_files_exist = false;
    end
end

if ~all_files_exist
    error('Required files missing. Please ensure ISMRM data and field maps are available.');
end

%% Configure enhanced preprocessing options
options = struct();

% Field map correction settings
options.run_fieldmap_correction = true;
options.fieldmap_file = 'ISMRM/ISMRM_fmap_Hz.nii.gz';  % Use Hz field map
options.fieldmap_units = 'Hz';
options.phase_encoding_dir = 'y';  % Adjust based on your acquisition
options.dwell_time = 0.00058;  % Adjust based on your acquisition parameters

% Enhanced eddy correction
options.eddy_method = 'eddy_correct';  % Use fallback method since no acqp/index files

% Standard preprocessing options
options.run_denoising = true;
options.denoise_method = 'dwidenoise';
options.run_motion_correction = true;
options.create_wm_mask = true;  % Create white matter mask for better tractography
options.atlas_type = 'HarvardOxford';

%% Display configuration
fprintf('\nEnhanced Preprocessing Configuration:\n');
fprintf('  Field map correction: %s\n', char(string(options.run_fieldmap_correction)));
fprintf('  Field map file: %s\n', options.fieldmap_file);
fprintf('  Phase encoding: %s\n', options.phase_encoding_dir);
fprintf('  Eddy method: %s\n', options.eddy_method);
fprintf('  Create WM mask: %s\n', char(string(options.create_wm_mask)));

%% Run enhanced preprocessing
fprintf('\nStarting enhanced preprocessing...\n');
try
    nim_preprocessing(data_prefix, options);

    fprintf('\n✅ Enhanced preprocessing completed successfully!\n');

    %% Verify output quality
    fprintf('\nVerifying output quality...\n');

    % Check if enhanced files were created
    output_files = {
        [data_prefix '.nii.gz']; ...
        [data_prefix '_M.nii.gz']; ...
        [data_prefix '_WM_mask.nii.gz']; ...
        [data_prefix '_preprocessing_report.mat']
    };

    for i = 1:length(output_files)
        if exist(output_files{i}, 'file')
            file_info = dir(output_files{i});
            fprintf('  ✓ %s (%.1f MB)\n', output_files{i}, file_info.bytes/1024/1024);
        else
            fprintf('  ⚠ %s (not created)\n', output_files{i});
        end
    end

    %% Load and display preprocessing report
    report_file = [data_prefix '_preprocessing_report.mat'];
    if exist(report_file, 'file')
        fprintf('\nLoading preprocessing report...\n');
        load(report_file, 'preprocessing_report');

        fprintf('Enhanced features applied:\n');
        if isfield(preprocessing_report, 'enhanced_features')
            for i = 1:length(preprocessing_report.enhanced_features)
                fprintf('  • %s\n', preprocessing_report.enhanced_features{i});
            end
        end

        if ~isempty(preprocessing_report.warnings)
            fprintf('Warnings:\n');
            for i = 1:length(preprocessing_report.warnings)
                fprintf('  ⚠ %s\n', preprocessing_report.warnings{i});
            end
        end
    end

    %% Suggest next steps
    fprintf('\n🚀 Next Steps:\n');
    fprintf('1. Run DTI analysis:\n');
    fprintf('   main(''%s.nii.gz'', ''%s_processed.mat'')\n', data_prefix, data_prefix);
    fprintf('\n2. Run improved tractography:\n');
    fprintf('   runTractography(''%s_processed.mat'')\n', data_prefix);
    fprintf('\n3. Compare with original results using slice viewer:\n');
    fprintf('   python tractography_slice_gui.py\n');

    fprintf('\n✨ Expected improvements:\n');
    fprintf('  • Reduced edge artifacts from field map correction\n');
    fprintf('  • Better fiber coherence from eddy correction\n');
    fprintf('  • Improved seeding from white matter mask\n');
    fprintf('  • Cleaner tractography boundaries\n');

catch ME
    fprintf('\n❌ Enhanced preprocessing failed: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d): %s\n', ME.stack(i).name, ME.stack(i).line, ME.stack(i).file);
    end

    fprintf('\n🔧 Troubleshooting suggestions:\n');
    fprintf('1. Check FSL installation: echo $FSLDIR\n');
    fprintf('2. Verify field map file format and units\n');
    fprintf('3. Adjust phase encoding direction if needed\n');
    fprintf('4. Check dwell time parameter for your scanner\n');

    rethrow(ME);
end

end