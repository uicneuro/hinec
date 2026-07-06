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
        warning(ME.identifier, 'Failed to compute DTI->T1 quality: %s', ME.message);
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
        warning(ME.identifier, 'Failed to compute T1->MNI quality: %s', ME.message);
    end
end

fprintf('  ✓ Quality metrics computed\n');

end
