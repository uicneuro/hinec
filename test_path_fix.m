function test_path_fix()
% TEST_PATH_FIX: Verify the path construction fix works correctly

fprintf('Testing path construction fix...\n');

% Test case: ISMRM/ISMRM file prefix (the problematic case)
file_prefix = 'ISMRM/ISMRM';
[~, filename_only] = fileparts(file_prefix);
output_file = [file_prefix '.nii.gz'];
output_dir = fileparts(output_file);
if isempty(output_dir)
    output_dir = pwd;
end

% Fixed path construction
t1_to_dwi_mat = fullfile(output_dir, [filename_only '_T1_to_dwi.mat']);
t1_to_dwi_output = fullfile(output_dir, [filename_only '_T1_to_dwi']);

fprintf('\nTest Case: file_prefix = ''%s''\n', file_prefix);
fprintf('filename_only: %s\n', filename_only);
fprintf('output_dir: %s\n', output_dir);
fprintf('t1_to_dwi_mat: %s\n', t1_to_dwi_mat);
fprintf('t1_to_dwi_output: %s\n', t1_to_dwi_output);

% Verify no duplicate paths
if contains(t1_to_dwi_output, 'ISMRM/ISMRM/ISMRM')
    fprintf('❌ FAILED: Still contains duplicate path\n');
else
    fprintf('✅ PASSED: No duplicate path detected\n');
end

% Test case 2: Simple prefix (should still work)
file_prefix2 = 'sample';
[~, filename_only2] = fileparts(file_prefix2);
output_file2 = [file_prefix2 '.nii.gz'];
output_dir2 = fileparts(output_file2);
if isempty(output_dir2)
    output_dir2 = pwd;
end

t1_to_dwi_output2 = fullfile(output_dir2, [filename_only2 '_T1_to_dwi']);

fprintf('\nTest Case 2: file_prefix = ''%s''\n', file_prefix2);
fprintf('t1_to_dwi_output: %s\n', t1_to_dwi_output2);

if strcmp(t1_to_dwi_output2, 'sample_T1_to_dwi')
    fprintf('✅ PASSED: Simple prefix works correctly\n');
else
    fprintf('❌ FAILED: Simple prefix broken\n');
end

fprintf('\nPath construction fix verification complete.\n');
end