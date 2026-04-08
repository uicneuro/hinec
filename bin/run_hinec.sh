#!/usr/bin/env bash
# run_hinec.sh - Unified HINEC preprocessing and tractography pipeline with YAML config
#
# Runs complete HINEC workflow:
#   1. Preprocessing (main.m): Raw data → DTI → FA → tissue segmentation
#   2. Tractography (nim_tractography_hinec): High-order fiber tracking
#
# Usage:
#   ./run_hinec.sh <data_prefix> <output_mat> [config_file]
#
# Arguments:
#   data_prefix - Path prefix for input NIfTI files (e.g., 'sample' for 'sample_raw.nii.gz')
#   output_mat  - Path to output .mat file (e.g., 'processed.mat')
#   config_file - Optional YAML configuration file (default: config/hinec_default.yml)
#
# Examples:
#   # Use default configuration
#   ./run_hinec.sh data/parcellation_sample/sample sample_hinec.mat
#
#   # Use custom configuration
#   ./run_hinec.sh data/parcellation_sample/sample sample_hinec.mat config/high_precision.yml
#
# Input files expected:
#   <data_prefix>_raw.nii.gz  - Raw diffusion data
#   <data_prefix>.bval        - b-values
#   <data_prefix>.bvec        - b-vectors
#   <data_prefix>_T1.nii.gz   - T1 anatomical (optional)
#
# Output files:
#   <output_mat>              - Processed nim structure
#   tractography_results/tracks_hinec_*.mat - HINEC tractography results

set -euo pipefail

# Parse arguments
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <data_prefix> <output_mat> [config_file]" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  $0 data/parcellation_sample/sample sample_hinec.mat" >&2
    echo "  $0 data/parcellation_sample/sample sample_hinec.mat config/high_precision.yml" >&2
    exit 1
fi

data_prefix=$1
output_mat=$2
config_file="${3:-config/hinec_default.yml}"  # Default config if not specified

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

# Validate config file
if [[ ! -f "$config_file" ]]; then
    echo "Error: Configuration file not found: $config_file" >&2
    echo "Available configs:" >&2
    ls -1 config/*.yml 2>/dev/null || echo "  No config files found in config/" >&2
    exit 1
fi

# Check for optional T1
t1_file="${data_prefix}_T1.nii.gz"
has_t1=false
if [[ -f "$t1_file" ]]; then
    has_t1=true
    printf 'Found T1 anatomical: %s\n' "$t1_file"
fi

# Helper function to escape MATLAB strings
escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

# Create run directory from shell (matches MATLAB create_run_directory naming)
timestamp=$(date +"%Y%m%d_%H%M%S")
config_name=$(basename "$config_file" .yml)
run_name="run_${timestamp}_${config_name}"
run_dir="hinec_runs/${run_name}"

mkdir -p "${run_dir}/logs" "${run_dir}/intermediate" "${run_dir}/output" "${run_dir}/tractography/diagnostics"
cp "$config_file" "${run_dir}/config.yml"

log_file="${run_dir}/logs/pipeline.log"

# Build MATLAB command
printf '\n========================================\n'
printf 'HINEC UNIFIED PIPELINE (YAML CONFIG)\n'
printf '========================================\n'
printf 'Input: %s\n' "$raw_file"
printf 'Output: %s\n' "$output_mat"
printf 'Config: %s\n' "$config_file"
printf 'Run dir: %s\n' "$run_dir"
printf 'Log: %s\n' "$log_file"
printf '========================================\n\n'

# Build combined MATLAB command for both phases
printf 'Starting background execution...\n'

# Combine both phases into single MATLAB command
# Pass pre-created run directory to MATLAB instead of creating it there
matlab_command="addpath(genpath('.')); "
matlab_command+="fprintf('\\n========================================\\n'); "
matlab_command+="fprintf('LOADING CONFIGURATION\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="config = load_config_yaml('$(escape_matlab_string "$config_file")'); "
matlab_command+="fprintf('\\n========================================\\n'); "
matlab_command+="fprintf('USING RUN DIRECTORY\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="run_info = create_run_directory('$(escape_matlab_string "$config_file")', 'run_name', '$(escape_matlab_string "$run_name")'); "
matlab_command+="fprintf('\\n========================================\\n'); "
matlab_command+="fprintf('PHASE 1: Preprocessing and DTI calculation\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="main('$(escape_matlab_string "$data_prefix")', '$(escape_matlab_string "$output_mat")', config, run_info); "
matlab_command+="fprintf('\\n✅ Phase 1 complete\\n\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="fprintf('PHASE 2: Tractography\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="[~, output_mat_name, output_mat_ext] = fileparts('$(escape_matlab_string "$output_mat")'); "
matlab_command+="output_mat_path = fullfile(run_info.output_dir, [output_mat_name output_mat_ext]); "
matlab_command+="runTractography(output_mat_path, config, run_info); "
matlab_command+="fprintf('\\n========================================\\n'); "
matlab_command+="fprintf('HINEC PIPELINE COMPLETE\\n'); "
matlab_command+="fprintf('========================================\\n'); "
matlab_command+="fprintf('Run directory: %s\\n', run_info.run_dir); "
matlab_command+="fprintf('  Output: %s\\n', run_info.output_dir); "
matlab_command+="fprintf('  Tractography: %s\\n', run_info.tractography_dir); "
matlab_command+="fprintf('  Logs: %s\\n', run_info.logs_dir); "
matlab_command+="fprintf('Configuration: %s\\n', '$(escape_matlab_string "$config_file")'); "
matlab_command+="fprintf('========================================\\n');"

# Run in background with nohup, log goes directly into run directory
nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &
pid=$!

printf '\n========================================\n'
printf 'HINEC PIPELINE LAUNCHED IN BACKGROUND\n'
printf '========================================\n'
printf 'Process ID (PID): %d\n' "$pid"
printf 'Run directory: %s\n' "$run_dir"
printf 'Config file: %s\n' "$config_file"
printf '\nDirectory structure:\n'
printf '  %s/\n' "$run_dir"
printf '    ├── config.yml          (copy of config used)\n'
printf '    ├── logs/pipeline.log   (live log output)\n'
printf '    ├── intermediate/       (preprocessing artifacts)\n'
printf '    ├── output/             (processed .mat file)\n'
printf '    ├── tractography/       (tracks + diagnostics)\n'
printf '    └── figures/            (visualization output)\n'
printf '========================================\n'
printf '\nMonitor progress:\n'
printf '  tail -f %s\n' "$log_file"
printf '\nCheck if running:\n'
printf '  ps -p %d\n' "$pid"
printf '\nKill if needed:\n'
printf '  kill %d\n' "$pid"
printf '\nView results when complete:\n'
printf '  ls -lt hinec_runs/latest/\n'
printf '  ls -lt %s/\n' "$run_dir"
printf '\nEdit configuration:\n'
printf '  nano %s\n' "$config_file"
printf '========================================\n'
