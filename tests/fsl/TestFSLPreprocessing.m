classdef TestFSLPreprocessing < matlab.unittest.TestCase
% TestFSLPreprocessing: Tests for FSL-dependent preprocessing functions
%
% All tests use assumeNotEmpty(getenv('FSLDIR')) to skip gracefully
% when FSL is not installed (e.g., on local dev machines).

    properties
        FixtureDir   % temp directory with phantom files
    end

    methods (TestClassSetup)
        function setupFixtures(testCase)
            projRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
            addpath(genpath(fullfile(projRoot, 'src')));
            addpath(genpath(fullfile(projRoot, 'lib', 'bfgs')));
            addpath(fullfile(projRoot, 'tests', 'fixtures'));

            testCase.FixtureDir = tempname;
            mkdir(testCase.FixtureDir);
            create_test_phantom(testCase.FixtureDir);
        end
    end

    methods (TestClassTeardown)
        function cleanup(testCase)
            if exist(testCase.FixtureDir, 'dir')
                rmdir(testCase.FixtureDir, 's');
            end
        end
    end

    methods (Test)
        function testFSLInstalled(testCase)
            fsldir = getenv('FSLDIR');
            testCase.assumeNotEmpty(fsldir, 'FSLDIR not set - skipping FSL tests');
            testCase.verifyTrue(isfolder(fsldir), 'FSLDIR should point to a valid directory');
        end

        function testFslinfo(testCase)
            testCase.assumeNotEmpty(getenv('FSLDIR'), 'FSLDIR not set');
            phantomFile = fullfile(testCase.FixtureDir, 'phantom.nii.gz');
            [status, ~] = system(sprintf('fslinfo %s', phantomFile));
            testCase.verifyEqual(status, 0, 'fslinfo should succeed on phantom data');
        end

        function testFslroi(testCase)
            testCase.assumeNotEmpty(getenv('FSLDIR'), 'FSLDIR not set');
            phantomFile = fullfile(testCase.FixtureDir, 'phantom.nii.gz');
            b0File = fullfile(testCase.FixtureDir, 'phantom_b0.nii.gz');
            % Extract first volume (b0)
            [status, ~] = system(sprintf('fslroi %s %s 0 1', phantomFile, b0File));
            testCase.verifyEqual(status, 0, 'fslroi should succeed');
            testCase.verifyTrue(isfile(b0File), 'Extracted b0 file should exist');
        end

        function testBet(testCase)
            testCase.assumeNotEmpty(getenv('FSLDIR'), 'FSLDIR not set');
            phantomFile = fullfile(testCase.FixtureDir, 'phantom.nii.gz');
            b0File = fullfile(testCase.FixtureDir, 'phantom_bet_b0.nii.gz');
            betOutput = fullfile(testCase.FixtureDir, 'phantom_bet');
            maskFile = fullfile(testCase.FixtureDir, 'phantom_bet_mask.nii.gz');

            % Extract b0 first
            system(sprintf('fslroi %s %s 0 1', phantomFile, b0File));

            % Run BET
            [status, ~] = system(sprintf('bet %s %s -m -f 0.3', b0File, betOutput));
            testCase.verifyEqual(status, 0, 'BET should succeed');
            testCase.verifyTrue(isfile(maskFile), 'BET mask file should exist');
        end

        function testBrainExtraction(testCase)
            testCase.assumeNotEmpty(getenv('FSLDIR'), 'FSLDIR not set');
            phantomFile = fullfile(testCase.FixtureDir, 'phantom.nii.gz');
            outputDir = fullfile(testCase.FixtureDir, 'bet_output');
            mkdir(outputDir);
            maskFile = fullfile(outputDir, 'brain_mask.nii.gz');

            % Call preproc_brain_extraction
            result = preproc_brain_extraction(phantomFile, outputDir, maskFile);
            testCase.verifyTrue(isfile(result), 'Brain mask should be created');

            % Verify mask is binary and 3D
            maskData = niftiread(result);
            testCase.verifyTrue(all(maskData(:) == 0 | maskData(:) == 1), ...
                'Mask should be binary');
            testCase.verifySize(maskData, [16 16 16]);
        end
    end
end
