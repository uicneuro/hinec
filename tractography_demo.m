%% Simple Tractography Demo
% This script demonstrates basic tractography functionality

clear; close all; clc;

fprintf('=== HINEC Tractography Demo ===\n');
fprintf('Running simplified tractography pipeline...\n\n');

%% Run tractography
runTractography('sample_parcellated.mat');

fprintf('\nDemo completed!\n');
fprintf('Check the generated figure for results.\n'); 