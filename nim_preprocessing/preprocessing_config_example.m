function options = preprocessing_config_example(profile)
% preprocessing_config_example: Generate preprocessing configuration examples
%
% Arguments:
%   profile - Configuration profile:
%     'quick' - Fast processing, minimal corrections
%     'standard' - Balanced processing (recommended)
%     'comprehensive' - All corrections enabled (slowest but best quality)
%     'tractography_optimized' - Optimized for tractography quality
%
% Returns:
%   options - Structure with preprocessing options

if nargin < 1
    profile = 'standard';
end

switch lower(profile)
    case 'quick'
        % Fast processing with minimal corrections
        options = struct();
        options.run_denoising = false;
        options.denoise_method = 'gaussian';
        options.run_motion_correction = false;
        options.run_eddy = false;
        options.improve_mask = false;
        options.atlas_type = 'HarvardOxford';
        
        fprintf('Quick processing profile:\n');
        fprintf('  - No denoising\n');
        fprintf('  - No motion correction\n');
        fprintf('  - No eddy correction\n');
        fprintf('  - Basic brain mask\n');
        fprintf('  - Processing time: ~2-5 minutes\n');
        
    case 'standard'
        % Balanced processing (recommended for most users)
        options = struct();
        options.run_denoising = true;
        options.denoise_method = 'dwidenoise';
        options.run_motion_correction = true;
        options.run_eddy = false;  % Often requires special parameter files
        options.improve_mask = true;
        options.atlas_type = 'HarvardOxford';
        
        fprintf('Standard processing profile (recommended):\n');
        fprintf('  - MRtrix3 MP-PCA denoising\n');
        fprintf('  - Motion correction with b-vector rotation\n');
        fprintf('  - Eddy correction (if parameter files available)\n');
        fprintf('  - Improved brain masking\n');
        fprintf('  - Processing time: ~10-20 minutes\n');
        
    case 'comprehensive'
        % All corrections enabled (best quality)
        options = struct();
        options.run_denoising = true;
        options.denoise_method = 'dwidenoise';
        options.run_motion_correction = true;
        options.run_eddy = true;
        options.improve_mask = true;
        options.atlas_type = 'JHU-tract';
        
        fprintf('Comprehensive processing profile:\n');
        fprintf('  - MRtrix3 MP-PCA denoising\n');
        fprintf('  - Motion correction with b-vector rotation\n');
        fprintf('  - Eddy current correction (requires acqp.txt and index.txt)\n');
        fprintf('  - Advanced brain mask improvement\n');
        fprintf('  - JHU tractography atlas\n');
        fprintf('  - Processing time: ~1-3 hours (depending on eddy)\n');
        
    case 'tractography_optimized'
        % Optimized specifically for tractography quality
        options = struct();
        options.run_denoising = true;
        options.denoise_method = 'dwidenoise';
        options.run_motion_correction = true;
        options.run_eddy = true;
        options.improve_mask = true;  % Critical for tractography
        options.atlas_type = 'JHU-tract';
        
        fprintf('Tractography-optimized processing profile:\n');
        fprintf('  - MP-PCA denoising for clean diffusion estimates\n');
        fprintf('  - Motion correction to prevent wobbly streamlines\n');
        fprintf('  - Eddy correction for accurate tensor estimation\n');
        fprintf('  - Advanced mask improvement to prevent tracts outside brain\n');
        fprintf('  - JHU tractography atlas for white matter analysis\n');
        fprintf('  - Processing time: ~1-3 hours\n');
        fprintf('  - ⭐ RECOMMENDED for fixing tractography issues\n');
        
    otherwise
        error('Unknown profile: %s. Available: quick, standard, comprehensive, tractography_optimized', profile);
end

fprintf('\nUsage:\n');
fprintf('  options = preprocessing_config_example(''%s'');\n', profile);
fprintf('  nim_preprocessing(''path/to/your/data'', options);\n\n');

end