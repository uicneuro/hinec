#!/usr/bin/env bash
# run_hinec.sh - Unified HINEC preprocessing and tractography pipeline
#
# Runs complete HINEC workflow:
#   1. Preprocessing (main.m): Raw data → DTI → FA → tissue segmentation
#   2. Tractography (nim_tractography_hinec): High-order fiber tracking
#
# Usage:
#   ./run_hinec.sh <data_prefix> <output_mat>
#
# Arguments:
#   data_prefix - Path prefix for input NIfTI files (e.g., 'sample' for 'sample_raw.nii.gz')
#   output_mat  - Path to output .mat file (e.g., 'processed.mat')
#
# Example:
#   ./run_hinec.sh nifti_sample/sample sample_hinec.mat
#
# Input files expected:
#   <data_prefix>_raw.nii.gz  - Raw diffusion data
#   <data_prefix>.bval        - b-values
#   <data_prefix>.bvec        - b-vectors
#   <data_prefix>_T1.nii.gz   - T1 anatomical (optional)
#
# Output files:
#   <output_mat>              - Processed nim structure
#   tractography_results/tracks_hinec.mat - HINEC tractography results

set -euo pipefail

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <data_prefix> <output_mat>" >&2
    echo "Example: $0 nifti_sample/sample sample_hinec.mat" >&2
    exit 1
fi

data_prefix=$1
output_mat=$2

# Validate input files
raw_file="${data_prefix}_raw.nii.gz"
if [[ ! -f "$raw_file" ]]; then
    echo "Error: Raw diffusion file not found: $raw_file" >&2
    exit 1
fi

bval_file="${data_prefix}.bval"
if [[ ! -f "$bval_file" ]]; then
    echo "Error: b-values file not found: $bval_file" >&2
    exit 1
fi

bvec_file="${data_prefix}.bvec"
if [[ ! -f "$bvec_file" ]]; then
    echo "Error: b-vectors file not found: $bvec_file" >&2
    exit 1
fi

# Check for optional T1
t1_file="${data_prefix}_T1.nii.gz"
has_t1=false
if [[ -f "$t1_file" ]]; then
    has_t1=true
    printf 'Found T1 anatomical: %s\n' "$t1_file"
fi

# Create timestamped log file
timestamp=$(date +"%Y-%m-%d_%H_%M_%S")
log_file="hinec_${timestamp}.log"

# Helper function to escape MATLAB strings
escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

# Build MATLAB command
printf '\n========================================\n'
printf 'HINEC UNIFIED PIPELINE\n'
printf '========================================\n'
printf 'Input: %s\n' "$raw_file"
printf 'Output: %s\n' "$output_mat"
printf 'Log: %s\n' "$log_file"
printf '========================================\n\n'

# Build combined MATLAB command for both phases
printf 'Starting background execution...\n'

# Combine both phases into single MATLAB command
matlab_command="addpath(genpath('.')); "
matlab_command+="fprintf('\\n========================================\\n'); "
matlab_command+="fprintf('PHASE 1: Preprocessing and DTI calculation\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="main('$(escape_matlab_string "$data_prefix")', '$(escape_matlab_string "$output_mat")'); "
matlab_command+="fprintf('\\n✅ Phase 1 complete: %s created\\n\\n', '$(escape_matlab_string "$output_mat")'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="fprintf('PHASE 2: HINEC High-order Tractography\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="runTractography('$(escape_matlab_string "$output_mat")', 'hinec'); "
matlab_command+="fprintf('\\n========================================\\n'); "
matlab_command+="fprintf('HINEC PIPELINE COMPLETE\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="fprintf('Processed data: %s\\n', '$(escape_matlab_string "$output_mat")'); "
matlab_command+="fprintf('Tractography: tractography_results/tracks_hinec.mat\\n'); "
matlab_command+="fprintf('========================================\\n');"

# Run in background with nohup
nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &
pid=$!

printf '\n========================================\n'
printf 'HINEC PIPELINE LAUNCHED IN BACKGROUND\n'
printf '========================================\n'
printf 'Process ID (PID): %d\n' "$pid"
printf 'Log file: %s\n' "$log_file"
printf 'Output .mat: %s\n' "$output_mat"
printf 'Tractography: tractography_results/tracks_hinec_YYYY-MM-DD_HH_MM_SS.mat\n'
printf '              (timestamped filename will be generated)\n'
printf '========================================\n'
printf '\nMonitor progress:\n'
printf '  tail -f %s\n' "$log_file"
printf '\nCheck if running:\n'
printf '  ps -p %d\n' "$pid"
printf '\nKill if needed:\n'
printf '  kill %d\n' "$pid"
printf '\nView results when complete:\n'
printf '  ls -lt tractography_results/\n'
printf '========================================\n'
