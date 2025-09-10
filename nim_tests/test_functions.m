% Simple test script for tractography diagnostics
%
% This script ONLY tests if the diagnostic functions work properly
% It does NOT run full tractography (use runTractography.m for that)

fprintf('=== QUICK DIAGNOSTIC TESTS ===\n');

% Load test data
nim_file = 'sample_parcellated.mat';
if ~exist(nim_file, 'file')
    fprintf('ERROR: Test file not found: %s\n', nim_file);
    return;
end

data = load(nim_file);
nim = data.nim;

% Test 1: Diagnostic check (fast)
fprintf('\nTest 1: Running diagnostic check...\n');
nim_diagnostic_check(nim);

% Test 2: Corpus callosum test (fast)
fprintf('\nTest 2: Testing corpus callosum...\n');
cc_tracks = nim_test_corpus_callosum(nim);

% Summary
fprintf('\n=== TEST RESULTS ===\n');
if length(cc_tracks) > 0
    fprintf('✓ Diagnostics PASSED - %d corpus callosum tracks found\n', length(cc_tracks));
    fprintf('✓ Tractography functions are working\n');
else
    fprintf('✗ Diagnostics FAILED - no corpus callosum tracks\n');
    fprintf('✗ Check your data quality or tensor calculation\n');
end

fprintf('\nTo run full tractography: runTractography(''%s'')\n', nim_file);
fprintf('==============================\n'); 