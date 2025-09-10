function registration_data = extract_reference_volumes(registration_data, options)
% extract_reference_volumes: Extract reference volumes for registration
%
% This function extracts b0 volume from DWI data for registration purposes

fprintf('Extracting reference volumes...\n');

dwi_file = registration_data.input.dwi_file;
output_dir = registration_data.output_dir;

% Extract b0 volume from DWI
fprintf('  Extracting b0 volume from DWI...\n');

[dwi_dir, dwi_name, dwi_ext] = fileparts(dwi_file);
if strcmp(dwi_ext, '.gz')
    [~, dwi_name, dwi_ext2] = fileparts(dwi_name);
    dwi_ext = [dwi_ext2 dwi_ext];
end

b0_file = fullfile(output_dir, [dwi_name '_b0' dwi_ext]);

% Use FSL to extract first volume (assumed to be b0 or low b-value)
fsl_path = getenv('FSLDIR');
if ~isempty(fsl_path)
    cmd_extract = sprintf('%s/bin/fslroi %s %s 0 1', fsl_path, dwi_file, b0_file);
    [status, cmdout] = system(cmd_extract);
    
    if status ~= 0
        error('Failed to extract b0 volume: %s', cmdout);
    end
    
    fprintf('    ✓ B0 volume extracted: %s\n', b0_file);
else
    % Fallback: copy entire DWI file (not ideal but works)
    warning('FSL not found, using entire DWI file as reference');
    copyfile(dwi_file, b0_file);
end

% Store reference volume information
registration_data.reference_volumes = struct();
registration_data.reference_volumes.b0_file = b0_file;

fprintf('  ✓ Reference volumes extracted\n');

end

function registration_data = compute_registration_quality(registration_data, options)
% compute_registration_quality: Compute quality metrics for registration

fprintf('Computing registration quality metrics...\n');

registration_data.quality_metrics = struct();

% DTI to T1 quality
if isfield(registration_data.registered_images, 'b0_in_t1')
    fprintf('  Computing DTI->T1 quality...\n');
    
    try
        % Load registered b0 and T1 images
        if isfile(registration_data.registered_images.b0_in_t1)
            V_b0_reg = spm_vol(registration_data.registered_images.b0_in_t1);
            V_t1 = spm_vol(registration_data.input.t1_file);
            
            % Sample both images at same locations
            [X, Y, Z] = ndgrid(1:2:V_t1.dim(1), 1:2:V_t1.dim(2), 1:2:V_t1.dim(3));
            coords = [X(:), Y(:), Z(:)]';
            
            % Get image intensities
            b0_vals = spm_sample_vol(V_b0_reg, coords(1,:), coords(2,:), coords(3,:), 1);
            t1_vals = spm_sample_vol(V_t1, coords(1,:), coords(2,:), coords(3,:), 1);
            
            % Remove NaN values
            valid_idx = ~isnan(b0_vals) & ~isnan(t1_vals) & b0_vals > 0 & t1_vals > 0;
            b0_vals = b0_vals(valid_idx);
            t1_vals = t1_vals(valid_idx);
            
            if length(b0_vals) > 100
                % Compute normalized mutual information
                nmi = compute_normalized_mutual_information(b0_vals, t1_vals);
                registration_data.quality_metrics.dti_t1_nmi = nmi;
                
                fprintf('    DTI->T1 NMI: %.4f\n', nmi);
                
                if nmi > 0.3
                    fprintf('    ✓ Good registration quality\n');
                elseif nmi > 0.2
                    fprintf('    ⚠ Moderate registration quality\n');
                else
                    fprintf('    ❌ Poor registration quality\n');
                end
            end
        end
    catch ME
        warning('Failed to compute DTI->T1 quality: %s', ME.message);
    end
end

% T1 to MNI quality
if isfield(registration_data.registered_images, 't1_in_mni') && options.register_to_mni
    fprintf('  Computing T1->MNI quality...\n');
    
    try
        if isfile(registration_data.registered_images.t1_in_mni)
            V_t1_reg = spm_vol(registration_data.registered_images.t1_in_mni);
            V_mni = spm_vol(registration_data.input.mni_template);
            
            % Sample both images
            [X, Y, Z] = ndgrid(1:3:V_mni.dim(1), 1:3:V_mni.dim(2), 1:3:V_mni.dim(3));
            coords = [X(:), Y(:), Z(:)]';
            
            t1_vals = spm_sample_vol(V_t1_reg, coords(1,:), coords(2,:), coords(3,:), 1);
            mni_vals = spm_sample_vol(V_mni, coords(1,:), coords(2,:), coords(3,:), 1);
            
            % Remove NaN values
            valid_idx = ~isnan(t1_vals) & ~isnan(mni_vals) & t1_vals > 0 & mni_vals > 0;
            t1_vals = t1_vals(valid_idx);
            mni_vals = mni_vals(valid_idx);
            
            if length(t1_vals) > 100
                % Compute normalized mutual information
                nmi = compute_normalized_mutual_information(t1_vals, mni_vals);
                registration_data.quality_metrics.t1_mni_nmi = nmi;
                
                fprintf('    T1->MNI NMI: %.4f\n', nmi);
                
                if nmi > 0.4
                    fprintf('    ✓ Good registration quality\n');
                elseif nmi > 0.3
                    fprintf('    ⚠ Moderate registration quality\n');
                else
                    fprintf('    ❌ Poor registration quality\n');
                end
            end
        end
    catch ME
        warning('Failed to compute T1->MNI quality: %s', ME.message);
    end
end

fprintf('  ✓ Quality metrics computed\n');

end

function nmi = compute_normalized_mutual_information(img1, img2)
% Compute normalized mutual information between two image vectors

% Number of bins for histogram
nbins = 64;

% Normalize intensities to [0, nbins-1]
img1_norm = round((img1 - min(img1)) / (max(img1) - min(img1)) * (nbins - 1)) + 1;
img2_norm = round((img2 - min(img2)) / (max(img2) - min(img2)) * (nbins - 1)) + 1;

% Ensure valid range
img1_norm = max(1, min(nbins, img1_norm));
img2_norm = max(1, min(nbins, img2_norm));

% Compute joint histogram
joint_hist = accumarray([img1_norm(:), img2_norm(:)], 1, [nbins, nbins]);
joint_hist = joint_hist / sum(joint_hist(:));

% Compute marginal histograms
hist1 = sum(joint_hist, 2);
hist2 = sum(joint_hist, 1);

% Compute entropies
H1 = -sum(hist1(hist1 > 0) .* log2(hist1(hist1 > 0)));
H2 = -sum(hist2(hist2 > 0) .* log2(hist2(hist2 > 0)));
H12 = -sum(joint_hist(joint_hist > 0) .* log2(joint_hist(joint_hist > 0)));

% Compute normalized mutual information
nmi = (H1 + H2) / H12;

end

function save_registration_data(registration_data, options)
% save_registration_data: Save registration data to file

output_file = registration_data.output_file;

fprintf('Saving registration data...\n');

try
    save(output_file, 'registration_data', '-v7.3');
    
    file_info = dir(output_file);
    fprintf('  ✓ Registration data saved: %s (%.1f MB)\n', output_file, file_info.bytes/1024/1024);
    
catch ME
    error('Failed to save registration data: %s', ME.message);
end

end

function generate_registration_report(registration_data, options)
% generate_registration_report: Generate HTML report of registration results

fprintf('Generating registration report...\n');

report_file = [registration_data.output_prefix '_registration_report.html'];

try
    fid = fopen(report_file, 'w');
    
    % HTML header
    fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n');
    fprintf(fid, '<title>HINEC Registration Report</title>\n');
    fprintf(fid, '<style>\n');
    fprintf(fid, 'body { font-family: Arial, sans-serif; margin: 40px; }\n');
    fprintf(fid, 'h1 { color: #2E86AB; }\n');
    fprintf(fid, 'h2 { color: #A23B72; margin-top: 30px; }\n');
    fprintf(fid, 'table { border-collapse: collapse; width: 100%%; margin: 20px 0; }\n');
    fprintf(fid, 'th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }\n');
    fprintf(fid, 'th { background-color: #f2f2f2; }\n');
    fprintf(fid, '.quality-good { color: #28a745; font-weight: bold; }\n');
    fprintf(fid, '.quality-moderate { color: #ffc107; font-weight: bold; }\n');
    fprintf(fid, '.quality-poor { color: #dc3545; font-weight: bold; }\n');
    fprintf(fid, '</style>\n</head>\n<body>\n');
    
    % Report content
    fprintf(fid, '<h1>HINEC Registration Report</h1>\n');
    fprintf(fid, '<p><strong>Generated:</strong> %s</p>\n', char(registration_data.created));
    
    % Input files
    fprintf(fid, '<h2>Input Files</h2>\n');
    fprintf(fid, '<table>\n');
    fprintf(fid, '<tr><th>File Type</th><th>Path</th></tr>\n');
    fprintf(fid, '<tr><td>DWI</td><td>%s</td></tr>\n', registration_data.input.dwi_file);
    fprintf(fid, '<tr><td>T1</td><td>%s</td></tr>\n', registration_data.input.t1_file);
    if options.register_to_mni
        fprintf(fid, '<tr><td>MNI Template</td><td>%s</td></tr>\n', registration_data.input.mni_template);
    end
    fprintf(fid, '</table>\n');
    
    % Registration parameters
    fprintf(fid, '<h2>Registration Parameters</h2>\n');
    fprintf(fid, '<table>\n');
    fprintf(fid, '<tr><th>Parameter</th><th>Value</th></tr>\n');
    fprintf(fid, '<tr><td>Method</td><td>%s</td></tr>\n', options.registration_method);
    fprintf(fid, '<tr><td>DTI->T1 DOF</td><td>%d</td></tr>\n', options.dti_reg_dof);
    if options.register_to_mni
        fprintf(fid, '<tr><td>T1->MNI Type</td><td>%s</td></tr>\n', options.t1_mni_reg_type);
    end
    fprintf(fid, '</table>\n');
    
    % Quality metrics
    if isfield(registration_data, 'quality_metrics')
        fprintf(fid, '<h2>Registration Quality</h2>\n');
        fprintf(fid, '<table>\n');
        fprintf(fid, '<tr><th>Registration</th><th>NMI Score</th><th>Quality</th></tr>\n');
        
        if isfield(registration_data.quality_metrics, 'dti_t1_nmi')
            nmi = registration_data.quality_metrics.dti_t1_nmi;
            if nmi > 0.3
                quality_class = 'quality-good';
                quality_text = 'Good';
            elseif nmi > 0.2
                quality_class = 'quality-moderate';
                quality_text = 'Moderate';
            else
                quality_class = 'quality-poor';
                quality_text = 'Poor';
            end
            fprintf(fid, '<tr><td>DTI → T1</td><td>%.4f</td><td class="%s">%s</td></tr>\n', ...
                nmi, quality_class, quality_text);
        end
        
        if isfield(registration_data.quality_metrics, 't1_mni_nmi')
            nmi = registration_data.quality_metrics.t1_mni_nmi;
            if nmi > 0.4
                quality_class = 'quality-good';
                quality_text = 'Good';
            elseif nmi > 0.3
                quality_class = 'quality-moderate';
                quality_text = 'Moderate';
            else
                quality_class = 'quality-poor';
                quality_text = 'Poor';
            end
            fprintf(fid, '<tr><td>T1 → MNI</td><td>%.4f</td><td class="%s">%s</td></tr>\n', ...
                nmi, quality_class, quality_text);
        end
        
        fprintf(fid, '</table>\n');
    end
    
    % Output files
    fprintf(fid, '<h2>Output Files</h2>\n');
    fprintf(fid, '<table>\n');
    fprintf(fid, '<tr><th>Description</th><th>Path</th></tr>\n');
    
    if isfield(registration_data.registered_images, 'b0_in_t1')
        fprintf(fid, '<tr><td>B0 in T1 space</td><td>%s</td></tr>\n', ...
            registration_data.registered_images.b0_in_t1);
    end
    
    if isfield(registration_data.registered_images, 't1_in_dti')
        fprintf(fid, '<tr><td>T1 in DTI space</td><td>%s</td></tr>\n', ...
            registration_data.registered_images.t1_in_dti);
    end
    
    if isfield(registration_data.registered_images, 't1_in_mni')
        fprintf(fid, '<tr><td>T1 in MNI space</td><td>%s</td></tr>\n', ...
            registration_data.registered_images.t1_in_mni);
    end
    
    fprintf(fid, '<tr><td>Registration data</td><td>%s</td></tr>\n', registration_data.output_file);
    
    fprintf(fid, '</table>\n');
    
    % Footer
    fprintf(fid, '<hr>\n');
    fprintf(fid, '<p><em>Generated by HINEC Registration Pipeline</em></p>\n');
    fprintf(fid, '</body>\n</html>\n');
    
    fclose(fid);
    
    fprintf('  ✓ Registration report: %s\n', report_file);
    
catch ME
    warning('Failed to generate registration report: %s', ME.message);
    if exist('fid', 'var') && fid > 0
        fclose(fid);
    end
end

end