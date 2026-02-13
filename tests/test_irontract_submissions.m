% test_irontract_submissions.m - Quick validation of IronTract submissions
%
% Tests submission files to verify they're valid and non-empty

fprintf('=== IronTract Submission Validator ===\n\n');

submission_dir = 'irontract-challenge/submissions';

% Find all submission files
files = dir(fullfile(submission_dir, 'submission_*.nii.gz'));
if isempty(files)
    error('No submission files found in %s', submission_dir);
end

fprintf('Found %d submission files\n\n', length(files));

% Test each submission
all_valid = true;
results = struct();

for i = 1:length(files)
    filepath = fullfile(submission_dir, files(i).name);
    fprintf('Testing %s...\n', files(i).name);

    try
        % Read NIfTI file
        data = niftiread(filepath);
        info = niftiinfo(filepath);

        % Get statistics
        dims = size(data);
        num_voxels = sum(data(:) > 0);
        total_voxels = prod(dims);
        percentage = 100 * num_voxels / total_voxels;

        % Check validity
        is_valid = true;
        issues = {};

        if num_voxels == 0
            is_valid = false;
            issues{end+1} = 'ZERO VOXELS';
        end

        if any(dims ~= [96, 112, 96])
            is_valid = false;
            issues{end+1} = sprintf('WRONG DIMS: %dx%dx%d', dims);
        end

        if ~strcmp(info.Datatype, 'single') && ~strcmp(info.Datatype, 'uint8')
            issues{end+1} = sprintf('Unusual datatype: %s', info.Datatype);
        end

        % Store results
        results(i).filename = files(i).name;
        results(i).dims = dims;
        results(i).datatype = info.Datatype;
        results(i).num_voxels = num_voxels;
        results(i).percentage = percentage;
        results(i).valid = is_valid;
        results(i).issues = issues;

        % Print results
        if is_valid
            fprintf('  ✓ Valid: %d voxels (%.3f%% coverage)\n', num_voxels, percentage);
        else
            fprintf('  ❌ INVALID: %s\n', strjoin(issues, ', '));
            all_valid = false;
        end
        fprintf('  Dimensions: %d x %d x %d\n', dims);
        fprintf('  Datatype: %s\n', info.Datatype);
        fprintf('  File size: %.1f KB\n', files(i).bytes/1024);

    catch ME
        fprintf('  ❌ ERROR: %s\n', ME.message);
        all_valid = false;
        results(i).valid = false;
        results(i).issues = {ME.message};
    end

    fprintf('\n');
end

% Summary
fprintf('=== SUMMARY ===\n');
fprintf('Total submissions: %d\n', length(files));

valid_count = sum([results.valid]);
fprintf('Valid: %d\n', valid_count);
fprintf('Invalid: %d\n', length(files) - valid_count);

if all_valid
    fprintf('\n✓ ALL SUBMISSIONS VALID!\n');
    fprintf('\nNext steps:\n');
    fprintf('1. Submit files to IronTract challenge\n');
    fprintf('2. Compare with ground truth (if available)\n');
else
    fprintf('\n❌ SOME SUBMISSIONS INVALID\n');
    fprintf('Please regenerate submissions with fixed pipeline\n');
end

% Show voxel distribution
fprintf('\n=== Voxel Coverage Distribution ===\n');
for i = 1:length(results)
    if results(i).valid
        fprintf('%s: %6d voxels (%.3f%%)\n', ...
            results(i).filename, results(i).num_voxels, results(i).percentage);
    end
end

fprintf('\nExpected pattern:\n');
fprintf('  - submission_001 (strictest): Fewest voxels\n');
fprintf('  - submission_005 (loosest):   Most voxels\n');
fprintf('  - Progressive increase from 001 to 005\n');

% Check for monotonic increase
if all_valid && length(results) == 5
    voxel_counts = [results.num_voxels];
    is_increasing = all(diff(voxel_counts) >= 0);

    fprintf('\nMonotonic increase check: ');
    if is_increasing
        fprintf('✓ PASS (voxels increase as angle threshold increases)\n');
    else
        fprintf('⚠️  WARNING (voxels should increase from 001→005)\n');
    end
end
