#!/usr/bin/env bash
# Wrapper to launch generateSlices via MATLAB batch mode with timestamped logging.
# Runs slice generation in background for distributed viewing workflow.
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <tracks_file> <nim_file> [output_dir]" >&2
    echo "" >&2
    echo "Generate pre-computed slice images for fast viewing." >&2
    echo "" >&2
    echo "Arguments:" >&2
    echo "  tracks_file  - Path to tracks .mat file" >&2
    echo "  nim_file     - Path to nim .mat file" >&2
    echo "  output_dir   - Output directory for slice cache (optional)" >&2
    echo "" >&2
    echo "Example:" >&2
    echo "  $0 tracks.mat nim.mat /export/slices" >&2
    exit 1
fi

tracks_file=$1
nim_file=$2
output_dir=""
if [[ $# -eq 3 ]]; then
    output_dir=$3
fi

timestamp=$(date +"%Y-%m-%d_%H_%M_%S")
log_file="generateSlices_${timestamp}.log"

escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

# Build MATLAB command
matlab_command="addpath(genpath('.'));"
if [[ -n "$output_dir" ]]; then
    matlab_command+=" generateSlices('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")', '$(escape_matlab_string "$output_dir")');"
else
    matlab_command+=" generateSlices('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")');"
fi

# Launch in background with nohup
nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &

printf '╔═══════════════════════════════════════════════════════════╗\n'
printf '║   Slice Generation Launched in Background                ║\n'
printf '╚═══════════════════════════════════════════════════════════╝\n'
printf '\n'
printf 'Tracks file:  %s\n' "$tracks_file"
printf 'NIM file:     %s\n' "$nim_file"
if [[ -n "$output_dir" ]]; then
    printf 'Output dir:   %s\n' "$output_dir"
else
    printf 'Output dir:   ./tractography_cache (default)\n'
fi
printf 'Log file:     %s\n' "$log_file"
printf 'PID:          %d\n' "$!"
printf '\n'
printf 'Monitor progress:\n'
printf '  tail -f %s\n' "$log_file"
printf '\n'
printf 'Check if still running:\n'
printf '  ps -p %d\n' "$!"
printf '\n'
printf 'This will take 5-30 minutes depending on data size.\n'
printf '═══════════════════════════════════════════════════════════\n'