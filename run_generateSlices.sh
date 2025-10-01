#!/usr/bin/env bash
# Wrapper to launch generateSlices via MATLAB batch mode with timestamped logging.
# Runs slice generation in background for distributed viewing workflow.
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "Usage: $0 <tracks_file> <nim_file> [output_dir] [regions]" >&2
    echo "" >&2
    echo "Generate pre-computed slice images for fast viewing." >&2
    echo "" >&2
    echo "Arguments:" >&2
    echo "  tracks_file  - Path to tracks .mat file" >&2
    echo "  nim_file     - Path to nim .mat file" >&2
    echo "  output_dir   - Output directory for slice cache (optional)" >&2
    echo "  regions      - Region IDs to filter (optional, e.g., '1,2,3' or '5:10')" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  $0 tracks.mat nim.mat /export/slices" >&2
    echo "  $0 tracks.mat nim.mat /export/slices '1,2,3'" >&2
    echo "  $0 tracks.mat nim.mat '' '5:10'" >&2
    exit 1
fi

tracks_file=$1
nim_file=$2
output_dir=""
regions=""

if [[ $# -ge 3 && -n "$3" ]]; then
    output_dir=$3
fi

if [[ $# -eq 4 ]]; then
    regions=$4
fi

timestamp=$(date +"%Y-%m-%d_%H_%M_%S")
log_file="generateSlices_${timestamp}.log"

escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

# Build MATLAB command
matlab_command="addpath(genpath('.'));"

# Build options structure if regions specified
if [[ -n "$regions" ]]; then
    # Convert region string to MATLAB array format
    # '1,2,3' -> [1,2,3] or '5:10' -> [5:10]
    region_array="[${regions}]"
    matlab_command+=" opts = struct(); opts.regions = ${region_array};"

    if [[ -n "$output_dir" ]]; then
        matlab_command+=" generateSlices('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")', '$(escape_matlab_string "$output_dir")', opts);"
    else
        matlab_command+=" generateSlices('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")', '', opts);"
    fi
else
    if [[ -n "$output_dir" ]]; then
        matlab_command+=" generateSlices('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")', '$(escape_matlab_string "$output_dir")');"
    else
        matlab_command+=" generateSlices('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")');"
    fi
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
if [[ -n "$regions" ]]; then
    printf 'Regions:      %s\n' "$regions"
else
    printf 'Regions:      all (no filtering)\n'
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