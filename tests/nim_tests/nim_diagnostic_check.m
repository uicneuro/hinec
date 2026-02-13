function nim_diagnostic_check(nim)
% nim_diagnostic_check: Quick diagnostic check for tensor quality
%
% This function performs the diagnostic checks mentioned in the analysis:
% 1. Verify eigenvector orientations  
% 2. Check tensor validity
% 3. Look for data quality issues
%
% Arguments:
%   nim - Structure containing diffusion tensor data

fprintf('=== DIAGNOSTIC CHECK ===\n');

% Check if required fields exist
if ~isfield(nim, 'DT')
    fprintf('ERROR: Diffusion tensors not found. Run nim_dt_spd() first.\n');
    return;
end

if ~isfield(nim, 'evec') || ~isfield(nim, 'eval')
    fprintf('ERROR: Eigenvectors/eigenvalues not found. Run nim_eig() first.\n');
    return;
end

if ~isfield(nim, 'FA')
    fprintf('ERROR: FA values not found. Run nim_fa() first.\n');
    return;
end

dims = size(nim.FA);
fprintf('Volume dimensions: %d x %d x %d\n', dims);

% Check center voxel for eigenvector orientation
center_coords = round(dims/2);
x = center_coords(1); y = center_coords(2); z = center_coords(3);

fprintf('\nChecking center voxel [%d, %d, %d]:\n', x, y, z);
fprintf('FA: %.4f\n', nim.FA(x,y,z));

% Check eigenvalues
eigenvals = squeeze(nim.eval(x,y,z,:));
fprintf('Eigenvalues: [%.6f, %.6f, %.6f]\n', eigenvals(1), eigenvals(2), eigenvals(3));

% Verify eigenvalue sorting (largest should be first)
is_sorted = eigenvals(1) >= eigenvals(2) && eigenvals(2) >= eigenvals(3);
fprintf('Eigenvalues properly sorted: %s\n', mat2str(is_sorted));

if ~is_sorted
    fprintf('WARNING: Eigenvalues not properly sorted!\n');
end

% Check primary eigenvector
primary_evec = squeeze(nim.evec(x,y,z,:,1));
evec_norm = norm(primary_evec);
fprintf('Primary eigenvector: [%.4f, %.4f, %.4f] (norm: %.4f)\n', ...
    primary_evec(1), primary_evec(2), primary_evec(3), evec_norm);

if abs(evec_norm - 1.0) > 0.01
    fprintf('WARNING: Eigenvector not normalized!\n');
end

% Check for invalid tensors
invalid_count = 0;
total_voxels = numel(nim.DT(:,:,:,1));

for i = 1:dims(1)
    for j = 1:dims(2)
        for k = 1:dims(3)
            D = nim_reshape_d(nim.DT(i,j,k,:));
            if any(isnan(D(:))) || any(isinf(D(:)))
                invalid_count = invalid_count + 1;
            end
        end
    end
end

fprintf('\nTensor validation:\n');
fprintf('Total voxels: %d\n', total_voxels);
fprintf('Invalid tensors: %d (%.2f%%)\n', invalid_count, 100*invalid_count/total_voxels);

if invalid_count > total_voxels * 0.05
    fprintf('WARNING: High rate of invalid tensors (>5%%)!\n');
end

% Check FA distribution
fa_values = nim.FA(nim.FA > 0);
fprintf('\nFA statistics:\n');
fprintf('Mean FA: %.4f\n', mean(fa_values));
fprintf('Max FA: %.4f\n', max(fa_values));
fprintf('Voxels with FA > 0.2: %d (%.1f%%)\n', sum(nim.FA(:) > 0.2), 100*sum(nim.FA(:) > 0.2)/total_voxels);

if max(fa_values) > 0.9
    fprintf('WARNING: Suspiciously high FA values detected (>0.9) - possible CSF contamination\n');
end

fprintf('========================\n');
end 