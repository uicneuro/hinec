function runTractography(data_path, varargin)
% runTractography: Simple entry point for DTI tractography
%
% Usage:
%   % Standard tractography (generates tracks only)
%   runTractography('sample_parcellated.mat')
%
%   % IronTract Challenge mode (generates tracks + submission files)
%   runTractography('data.mat', 'IronTract', 'injection.nii.gz', 'submissions/')
%
% Arguments:
%   data_path - Path to .mat file with nim structure
%   'IronTract' - Optional flag to enable IronTract Challenge submission
%   injection_file - Path to injection site mask (required if IronTract enabled)
%   output_dir - Output directory for submissions (required if IronTract enabled)
%
% Examples:
%   % Process IronTract data
%   main('ironTract/sub-MR243', 'ironTract.mat');
%   runTractography('ironTract.mat', 'IronTract', 'ironTract/injection.nii.gz', 'ironTract_submissions/');

if nargin < 1
    data_path = 'sample_parcellated.mat';
end

% Parse optional IronTract arguments
enable_irontract = false;
injection_file = '';
irontract_output_dir = '';

if nargin >= 2 && strcmpi(varargin{1}, 'IronTract')
    enable_irontract = true;
    if nargin < 4
        error('IronTract mode requires: runTractography(data_path, ''IronTract'', injection_file, output_dir)');
    end
    injection_file = varargin{2};
    irontract_output_dir = varargin{3};
end

% Add necessary paths
addpath('nim_tractography');
addpath('nim_utils');
addpath('nim_plots');
addpath('nim_challenges');

fprintf('=== HINEC Tractography Pipeline ===\n');

%% Load data
fprintf('Loading data from %s...\n', data_path);
if ~exist(data_path, 'file')
    error('Data file not found: %s', data_path);
end
load(data_path, 'nim');

%% Check required fields
if ~isfield(nim, 'evec')
    error('Eigenvectors not found. Please run main() first to generate DTI data.');
end
if ~isfield(nim, 'FA')
    error('FA not found. Please run main() first to generate DTI data.');
end

%% Set FACT tractography parameters
fprintf('Setting up FACT tractography parameters...\n');
options = struct();
options.seed_density = 4;           % VERY DENSE: 8 seeds per voxel for complete uniform coverage
options.step_size = 0.5;            % FACT step size Δ for Euler method: r_{i+1} = r_i + e1(r_i)*Δ
options.fa_threshold = 0.15;        % Low FA threshold - not used for seeding, only for tracking quality
options.termination_fa = 0.15;      % Low termination threshold for extensive tracking
options.angle_thresh = 35;          % Maximum turning angle in degrees - slightly more permissive for FACT
options.max_steps = 1000;           % Maximum Euler steps per track
options.min_length = 35;            % Minimum track length in mm - filters out short spurious tracks
options.order = 1;                  % FACT uses first-order integration
options.interp_method = 'none';     % FACT samples nearest voxel tensor (no interpolation)

% ENHANCED SEEDING STRATEGY - Use brain mask for comprehensive coverage
[data_dir, data_name, ~] = fileparts(data_path);
if isempty(data_dir)
    data_dir = pwd;
end

% White matter masking disabled - use brain mask for seeding
seed_mask = [];

% Use brain mask for comprehensive seeding (preserves parcellation regions)
if isempty(seed_mask)
    % Strategy 1: Use preprocessed brain mask (best option)
    if isfield(nim, 'mask') && ~isempty(nim.mask) && any(nim.mask(:) > 0)
        brain_mask = nim.mask > 0.5;
        fprintf('Using preprocessed brain mask for comprehensive seeding\n');
        seeding_strategy = 'brain_mask';
    % Strategy 2: Expand parcellation mask to include surrounding brain tissue
    elseif isfield(nim, 'parcellation_mask') && any(nim.parcellation_mask(:) > 0)
        % Create brain-like mask from parcellation + expansion
        parcel_mask = nim.parcellation_mask > 0;
        % Dilate parcellation to capture surrounding brain tissue
        se = strel('sphere', 2);  % 2-voxel dilation
        brain_mask = imdilate(parcel_mask, se);
        % Constrain to reasonable FA values to avoid CSF
        brain_mask = brain_mask & (nim.FA > 0.05);
        fprintf('Using expanded parcellation mask for brain seeding\n');
        seeding_strategy = 'expanded_parcellation';
    % Strategy 3: FA-based fallback (conservative)
    else
        brain_mask = nim.FA > 0.1;  % Conservative FA threshold
        fprintf('⚠ Using FA-based brain boundary (FA > 0.1)\n');
        seeding_strategy = 'fa_based_brain';
    end

    seed_mask = brain_mask;
    fprintf('✓ Using %s for seeding (%d voxels)\n', seeding_strategy, sum(seed_mask(:)));
end
options.seed_mask = seed_mask;

% Show enhanced seeding statistics
total_seed_locations = sum(seed_mask(:));
estimated_seeds = total_seed_locations * options.seed_density;
fprintf('ENHANCED SEEDING STRATEGY: %s\n', seeding_strategy);
fprintf('  Seed locations: %d\n', total_seed_locations);
fprintf('  Seeds per voxel: %d\n', options.seed_density);
fprintf('  Estimated total seeds: %d\n', estimated_seeds);
fprintf('  Quality improvement: Boundary erosion reduces edge artifacts\n');

%% Run FACT tractography
fprintf('Running FACT tractography...\n');
tic;
tracks = nim_tractography_standard(nim, options);
elapsed_time = toc;

fprintf('Tractography completed in %.1f seconds\n', elapsed_time);
fprintf('Generated %d tracks\n', length(tracks));

if isempty(tracks)
    error('No tracks generated! Check FA threshold and seed mask.');
end

%% Compute statistics
track_lengths = cellfun(@(x) size(x, 1), tracks);
fprintf('\nTrack Statistics:\n');
fprintf('  Mean length: %.1f points\n', mean(track_lengths));
fprintf('  Max length: %d points\n', max(track_lengths));
fprintf('  Min length: %d points\n', min(track_lengths));

%% Save results
output_dir = 'tractography_results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

save(fullfile(output_dir, 'tracks_standard.mat'), 'tracks', 'options', 'elapsed_time');
fprintf('\nResults saved to %s/\n', output_dir);

%% IronTract Challenge Submission (if enabled)
if enable_irontract
    fprintf('\n=== IronTract Challenge Submission ===\n');
    fprintf('Generating submission files...\n');

    % Set up IronTract options
    irontract_opts = struct();
    irontract_opts.angle_thresholds = [30, 45, 60, 75, 90];
    % DON'T pass tracks_file - force regeneration with injection site seeds
    irontract_opts.tracks_file = '';  % Empty = regenerate with injection-site seeding
    irontract_opts.base_options = options;

    % Generate submission files
    nim_irontract_submit(data_path, injection_file, irontract_output_dir, irontract_opts);

    fprintf('IronTract submissions ready in %s/\n', irontract_output_dir);
end

fprintf('=== Tractography Complete ===\n');
end