function create_test_phantom(output_dir)
% create_test_phantom: Generate synthetic diffusion MRI phantom for testing
%
% Creates a 16x16x16 volume with known fiber orientations:
%   Left half (x <= 8):  fibers along X axis (FA ~ 0.8)
%   Right half (x > 8):  fibers along Z axis (FA ~ 0.8)
%   Background:          low signal (FA ~ 0)
%
% Generates files:
%   phantom.nii.gz              - 4D diffusion data (16x16x16x7)
%   phantom.bval                - b-values
%   phantom.bvec                - gradient directions
%   phantom_M.nii.gz            - binary brain mask
%   phantom_raw.nii.gz          - copy for preprocessing tests
%   phantom_acqp.txt            - acquisition parameters (for eddy)
%   phantom_index.txt           - volume index (for eddy)
%   parcellation_mask.nii.gz    - 4-region parcellation
%
% Arguments:
%   output_dir - Directory to write test files

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    dims = [16, 16, 16];
    n_directions = 6;
    n_volumes = n_directions + 1;

    % Diffusion parameters
    S0 = 1000;
    b_value = 1000;
    lambda_par = 0.0017;
    lambda_perp = 0.0003;

    % 6 non-collinear gradient directions
    grad_dirs = [
        1, 0, 0;
        0, 1, 0;
        0, 0, 1;
        1, 1, 0;
        0, 1, 1;
        1, 0, 1
    ];
    grad_dirs = grad_dirs ./ vecnorm(grad_dirs, 2, 2);

    % Brain mask (sphere)
    [X, Y, Z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
    center = (dims + 1) / 2;
    radius = 6;
    mask = sqrt((X - center(1)).^2 + (Y - center(2)).^2 + (Z - center(3)).^2) <= radius;

    % Generate diffusion signal: S(b,g) = S0 * exp(-b * g' * D * g)
    rng(42);
    img = zeros([dims, n_volumes]);

    for x = 1:dims(1)
        for y = 1:dims(2)
            for z = 1:dims(3)
                if mask(x, y, z)
                    if x <= dims(1) / 2
                        D = diag([lambda_par, lambda_perp, lambda_perp]);
                    else
                        D = diag([lambda_perp, lambda_perp, lambda_par]);
                    end

                    img(x, y, z, 1) = S0;
                    for i = 1:n_directions
                        g = grad_dirs(i, :);
                        img(x, y, z, i + 1) = S0 * exp(-b_value * (g * D * g'));
                    end
                else
                    img(x, y, z, :) = 10;
                end
            end
        end
    end

    % Add noise
    img = img + 5 * abs(randn(size(img)));
    img = max(img, 1);

    % --- Write files ---

    % 4D diffusion data
    niftiwrite(int16(img), fullfile(output_dir, 'phantom'), 'Compressed', true);

    % Copy as "raw" for preprocessing tests
    niftiwrite(int16(img), fullfile(output_dir, 'phantom_raw'), 'Compressed', true);

    % Brain mask
    niftiwrite(int16(mask), fullfile(output_dir, 'phantom_M'), 'Compressed', true);

    % Parcellation mask (4 quadrants within brain)
    parc = zeros(dims, 'int32');
    parc(1:8, 1:8, :) = 1;
    parc(9:16, 1:8, :) = 2;
    parc(1:8, 9:16, :) = 3;
    parc(9:16, 9:16, :) = 4;
    parc = parc .* int32(mask);
    niftiwrite(parc, fullfile(output_dir, 'parcellation_mask'), 'Compressed', true);

    % b-values (standard format: single line)
    fid = fopen(fullfile(output_dir, 'phantom.bval'), 'w');
    fprintf(fid, '0');
    for i = 1:n_directions
        fprintf(fid, ' %d', b_value);
    end
    fclose(fid);

    % b-vectors (standard format: 3 rows of components)
    all_bvecs = [0 0 0; grad_dirs]';
    fid = fopen(fullfile(output_dir, 'phantom.bvec'), 'w');
    for row = 1:3
        fprintf(fid, '%.6f', all_bvecs(row, 1));
        for col = 2:n_volumes
            fprintf(fid, ' %.6f', all_bvecs(row, col));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);

    % Acquisition parameters (for eddy)
    fid = fopen(fullfile(output_dir, 'phantom_acqp.txt'), 'w');
    fprintf(fid, '0 1 0 0.05\n');
    fclose(fid);

    % Volume index (for eddy)
    fid = fopen(fullfile(output_dir, 'phantom_index.txt'), 'w');
    for i = 1:n_volumes
        fprintf(fid, '1 ');
    end
    fclose(fid);

    fprintf('Test phantom created in: %s\n', output_dir);
end
