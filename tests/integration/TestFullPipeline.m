classdef TestFullPipeline < matlab.unittest.TestCase
% TestFullPipeline: Integration tests for the full HINEC processing pipeline
%
% Runs the complete pipeline once in TestClassSetup using synthetic phantom
% data, then individual tests verify aspects of the stored results.
%
% Pipeline: nim_read -> nim_dt_spd -> nim_eig -> nim_fa -> nim_parcellation -> tractography

    properties (ClassSetupParameter)
    end

    properties
        nim           % nim structure after full pipeline
        nimRead       % nim structure after nim_read only
        nimDT         % nim structure after nim_dt_spd
        nimEig        % nim structure after nim_eig
        nimFA         % nim structure after nim_fa
        nimParc       % nim structure after parcellation
        tracks        % tractography output
        fixtureDir    % path to phantom fixture files
        origDir       % original working directory
    end

    methods (TestClassSetup)
        function runPipeline(testCase)
            % Add project paths
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
            addpath(genpath(fullfile(projRoot, 'lib', 'bfgs')));
            addpath(fullfile(projRoot, 'tests', 'fixtures'));

            % Create phantom in temp directory
            testCase.fixtureDir = tempname;
            mkdir(testCase.fixtureDir);
            create_test_phantom(testCase.fixtureDir);

            phantomFile = fullfile(testCase.fixtureDir, 'phantom.nii.gz');
            parcFile = fullfile(testCase.fixtureDir, 'parcellation_mask.nii.gz');

            % Step 1: Read
            testCase.nimRead = nim_read(phantomFile);

            % Step 2: Diffusion tensor (SPD-constrained)
            testCase.nimDT = nim_dt_spd(testCase.nimRead);

            % Step 3: Eigendecomposition
            testCase.nimEig = nim_eig(testCase.nimDT);

            % Step 4: Fractional anisotropy
            testCase.nimFA = nim_fa(testCase.nimEig);

            % Step 5: Parcellation
            testCase.nimParc = nim_parcellation(testCase.nimFA, parcFile);

            % Step 6: Tractography (run from temp dir to contain output files)
            testCase.origDir = pwd;
            tractDir = tempname;
            mkdir(tractDir);
            cd(tractDir);

            options = struct();
            options.seed_density = 1;
            options.seed_mask = testCase.nimParc.mask;
            options.termination_fa = 0.1;
            options.fa_threshold = 0.1;
            options.angle_thresh = 60;
            options.max_steps = 500;
            options.min_length = 3;
            options.enable_diagnostics = false;
            testCase.tracks = nim_tractography_standard(testCase.nimParc, options);

            % Store final nim
            testCase.nim = testCase.nimParc;
        end
    end

    methods (TestClassTeardown)
        function cleanup(testCase)
            if ~isempty(testCase.origDir)
                cd(testCase.origDir);
            end
            if exist(testCase.fixtureDir, 'dir')
                rmdir(testCase.fixtureDir, 's');
            end
        end
    end

    methods (Test)
        function testNimRead(testCase)
            nim = testCase.nimRead;
            % 16x16x16 volume with 7 volumes (1 b0 + 6 DWI)
            testCase.verifyEqual(nim.xdim, 16);
            testCase.verifyEqual(nim.ydim, 16);
            testCase.verifyEqual(nim.zdim, 16);
            testCase.verifyEqual(nim.hdr.ImageSize(4), 7);
            testCase.verifyLength(nim.bval, 7);
            testCase.verifySize(nim.bvec, [7, 3]);
            testCase.verifySize(nim.mask, [16, 16, 16]);
        end

        function testDiffusionTensor(testCase)
            nim = testCase.nimDT;
            % DT is stored as 6-element vector per voxel
            testCase.verifySize(nim.DT, [16, 16, 16, 6]);
            % No NaN or Inf within mask
            mask = nim.mask == 1;
            dt_masked = nim.DT(repmat(mask, [1 1 1 6]));
            testCase.verifyTrue(all(isfinite(dt_masked(:))), 'DT should have no NaN/Inf in mask');
            % Eigenvalues should be non-negative within mask
            eval_masked = nim.eval(repmat(mask, [1 1 1 3]));
            testCase.verifyGreaterThanOrEqual(min(eval_masked(:)), -1e-6);
        end

        function testEigenvaluesSorted(testCase)
            nim = testCase.nimEig;
            mask = nim.mask == 1;
            % eval(:,:,:,1) >= eval(:,:,:,2) >= eval(:,:,:,3) within mask
            e1 = nim.eval(:,:,:,1);
            e2 = nim.eval(:,:,:,2);
            e3 = nim.eval(:,:,:,3);
            testCase.verifyTrue(all(e1(mask) >= e2(mask) - 1e-10), ...
                'First eigenvalue must be >= second within mask');
            testCase.verifyTrue(all(e2(mask) >= e3(mask) - 1e-10), ...
                'Second eigenvalue must be >= third within mask');
        end

        function testFractionalAnisotropy(testCase)
            nim = testCase.nimFA;
            % FA must be in [0, 1]
            testCase.verifyGreaterThanOrEqual(min(nim.FA(:)), 0);
            testCase.verifyLessThanOrEqual(max(nim.FA(:)), 1);
            % Mean FA in mask should be approximately 0.8
            mask = nim.mask == 1;
            mean_fa = mean(nim.FA(mask));
            testCase.verifyEqual(mean_fa, 0.8, 'AbsTol', 0.15, ...
                'Mean FA in mask should be ~0.8 for phantom with lambda_par=0.0017, lambda_perp=0.0003');
        end

        function testFAValues(testCase)
            nim = testCase.nimFA;
            % FA outside mask should be approximately 0
            mask = nim.mask == 1;
            fa_outside = nim.FA(~mask);
            testCase.verifyEqual(max(fa_outside(:)), 0, 'AbsTol', 1e-10, ...
                'FA outside mask should be 0');
        end

        function testEigenvectorDirections(testCase)
            nim = testCase.nimEig;
            mask = nim.mask == 1;
            % Left half (x<=8): primary eigenvector should align with [1,0,0]
            % Right half (x>8): primary eigenvector should align with [0,0,1]
            left_dots = [];
            right_dots = [];
            for x = 1:16
                for y = 1:16
                    for z = 1:16
                        if mask(x, y, z)
                            evec1 = squeeze(nim.evec(x, y, z, :, 1));
                            if x <= 8
                                left_dots(end+1) = abs(dot(evec1, [1; 0; 0])); %#ok<AGROW>
                            else
                                right_dots(end+1) = abs(dot(evec1, [1; 0; 0])); %#ok<AGROW>
                            end
                        end
                    end
                end
            end
            % Left half primary evec should align with X (dot product > 0.8)
            testCase.verifyGreaterThan(mean(left_dots), 0.8, ...
                'Left half primary eigenvector should align with X axis');
            % Right half primary evec should NOT align with X (dot < 0.3 = aligned with Z)
            testCase.verifyLessThan(mean(right_dots), 0.3, ...
                'Right half primary eigenvector should align with Z axis, not X');
        end

        function testParcellation(testCase)
            nim = testCase.nimParc;
            testCase.verifyTrue(isfield(nim, 'parcellation_mask'));
            % Should have 4 unique labels plus background (0)
            labels = unique(nim.parcellation_mask(:));
            testCase.verifyTrue(ismember(0, labels), 'Should have background label 0');
            non_zero_labels = labels(labels > 0);
            testCase.verifyEqual(numel(non_zero_labels), 4, ...
                'Should have exactly 4 parcellation regions');
        end

        function testTractography(testCase)
            tracks = testCase.tracks;
            testCase.verifyTrue(iscell(tracks), 'Tracks must be a cell array');
            testCase.verifyGreaterThan(numel(tracks), 0, 'Should produce at least some tracks');
            % Each track should be Nx3
            for i = 1:min(10, numel(tracks))
                testCase.verifySize(tracks{i}, [NaN, 3], ...
                    'Each track should be Nx3');
                testCase.verifyGreaterThanOrEqual(size(tracks{i}, 1), 2, ...
                    'Each track should have at least 2 points');
            end
        end

        function testTractographyBounds(testCase)
            tracks = testCase.tracks;
            % All coordinates should be within volume bounds [1, 16]
            for i = 1:numel(tracks)
                coords = tracks{i};
                testCase.verifyGreaterThanOrEqual(min(coords(:)), 0.5, ...
                    'Track coordinates must be >= 0.5');
                testCase.verifyLessThanOrEqual(max(coords(:)), 16.5, ...
                    'Track coordinates must be <= 16.5');
            end
        end
    end
end
