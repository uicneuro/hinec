#!/usr/bin/env bash
# Lightweight wrapper to launch runTractography via MATLAB batch mode with timestamped logging.
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 <input_mat_file> [output_mat_file]" >&2
    exit 1
fi

mat_file=$1
output_arg=""
if [[ $# -eq 2 ]]; then
    output_arg=$2
fi

timestamp=$(date +"%Y-%m-%d_%H_%M_%S")
log_file="output_${timestamp}.log"

escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

matlab_command="addpath(genpath('.'));"
if [[ -n "$output_arg" ]]; then
    matlab_command+=" runTractography('$(escape_matlab_string "$mat_file")', '$(escape_matlab_string "$output_arg")');"
else
    matlab_command+=" runTractography('$(escape_matlab_string "$mat_file")');"
fi

nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &

printf 'Launched tractography for %s\n' "$mat_file"
printf 'Logging to %s\n' "$log_file"
printf 'PID: %d\n' "$!"
