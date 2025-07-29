function cleanup_report = preproc_cleanup(output_dir, file_prefix, keep_files)
% preproc_cleanup: Clean up intermediate files and report final outputs
%
% Arguments:
%   output_dir - Directory containing files to clean up
%   file_prefix - Prefix used for files
%   keep_files - Cell array of files to keep (optional)
%
% Returns:
%   cleanup_report - Structure with cleanup statistics

fprintf('Step 7: Cleaning up intermediate files...\n');

% Default files to keep if not specified
if nargin < 3 || isempty(keep_files)
    keep_files = {
        file_prefix + ".nii.gz";                           % Processed DWI
        file_prefix + "_M.nii.gz";                         % Brain mask
        fullfile(output_dir, 'parcellation_mask.nii.gz');  % Parcellation
        file_prefix + "_atlas_labels.mat"                  % Atlas labels
    };
end

% Standard intermediate files that can be removed
intermediate_files = {
    fullfile(output_dir, 'b0.nii.gz');                          % Extracted b0 volume
    fullfile(output_dir, 'nodif_brain.nii.gz');                 % BET brain volume
    fullfile(output_dir, 'nodif_brain_mask.nii.gz');            % BET mask (before rename)
    fullfile(output_dir, 'resampled_parcellation_mask.nii.gz')  % Intermediate atlas
};

% Initialize cleanup report
cleanup_report = struct();
cleanup_report.files_removed = {};
cleanup_report.files_kept = {};
cleanup_report.removal_errors = {};
cleanup_report.files_removed_count = 0;
cleanup_report.space_freed_mb = 0;

% Remove intermediate files
for i = 1:length(intermediate_files)
    file_path = intermediate_files{i};
    
    if isfile(file_path)
        try
            % Get file size before deletion
            file_info = dir(file_path);
            file_size_mb = file_info.bytes / (1024 * 1024);
            
            % Delete the file
            delete(file_path);
            
            % Record successful removal
            cleanup_report.files_removed{end+1} = file_path;
            cleanup_report.files_removed_count = cleanup_report.files_removed_count + 1;
            cleanup_report.space_freed_mb = cleanup_report.space_freed_mb + file_size_mb;
            
            fprintf('  Removed: %s (%.1f MB)\n', file_path, file_size_mb);
            
        catch ME
            % Record removal error
            error_info = struct('file', file_path, 'error', ME.message);
            cleanup_report.removal_errors{end+1} = error_info;
            warning('Could not remove intermediate file %s: %s', file_path, ME.message);
        end
    end
end

fprintf('Cleanup complete: %d files removed, %.1f MB freed\n', ...
        cleanup_report.files_removed_count, cleanup_report.space_freed_mb);

% Verify and report final files
fprintf('\nFinal preprocessed files:\n');
total_size_mb = 0;

for i = 1:length(keep_files)
    file_path = keep_files{i};
    
    if isfile(file_path)
        file_info = dir(file_path);
        file_size_mb = file_info.bytes / (1024 * 1024);
        total_size_mb = total_size_mb + file_size_mb;
        
        % Record kept file
        file_record = struct('path', file_path, 'size_mb', file_size_mb, 'exists', true);
        cleanup_report.files_kept{end+1} = file_record;
        
        fprintf('  ✓ %s (%.1f MB)\n', file_path, file_size_mb);
    else
        % Record missing file
        file_record = struct('path', file_path, 'size_mb', 0, 'exists', false);
        cleanup_report.files_kept{end+1} = file_record;
        
        fprintf('  ✗ MISSING: %s\n', file_path);
    end
end

fprintf('Total final files size: %.1f MB\n', total_size_mb);
cleanup_report.total_final_size_mb = total_size_mb;

% Add summary statistics
cleanup_report.summary = struct();
cleanup_report.summary.intermediate_files_removed = cleanup_report.files_removed_count;
cleanup_report.summary.space_freed_mb = cleanup_report.space_freed_mb;
cleanup_report.summary.final_files_count = length(keep_files);
cleanup_report.summary.final_files_size_mb = total_size_mb;
cleanup_report.summary.removal_errors_count = length(cleanup_report.removal_errors);

% Check for any missing final files
missing_final_files = 0;
for i = 1:length(cleanup_report.files_kept)
    if ~cleanup_report.files_kept{i}.exists
        missing_final_files = missing_final_files + 1;
    end
end
cleanup_report.summary.missing_final_files = missing_final_files;

if missing_final_files > 0
    warning('⚠ %d expected final files are missing!', missing_final_files);
end

fprintf('Cleanup report complete.\n');

end