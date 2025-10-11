#!/usr/bin/env bash
# Lightweight wrapper to launch runTractography via MATLAB batch mode with timestamped logging.
#
# Usage:
#   Standard mode:
#     ./run_tractography.sh <input_mat_file>
#
#   IronTract mode:
#     ./run_tractography.sh <input_mat_file> IronTract <injection_mask> <output_dir>
#
# Examples:
#   ./run_tractography.sh sample_parcellated.mat
#   ./run_tractography.sh ironTract.mat IronTract ironTract/injection.nii.gz ironTract_submissions/

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_mat_file> [IronTract <injection_file> <output_dir>]" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  Standard mode:  $0 data.mat" >&2
    echo "  IronTract mode: $0 data.mat IronTract injection.nii.gz submissions/" >&2
    exit 1
fi

mat_file=$1
mode=""
injection_file=""
output_dir=""

# Parse IronTract mode arguments
if [[ $# -eq 4 && "$2" == "IronTract" ]]; then
    mode="IronTract"
    injection_file=$3
    output_dir=$4
elif [[ $# -gt 1 ]]; then
    echo "Error: Invalid arguments. For IronTract mode, use: $0 <mat_file> IronTract <injection_file> <output_dir>" >&2
    exit 1
fi

timestamp=$(date +"%Y-%m-%d_%H_%M_%S")
log_file="tractography_${timestamp}.log"

escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

matlab_command="addpath(genpath('.'));"

if [[ "$mode" == "IronTract" ]]; then
    matlab_command+=" runTractography('$(escape_matlab_string "$mat_file")', 'IronTract', '$(escape_matlab_string "$injection_file")', '$(escape_matlab_string "$output_dir")');"
    printf 'Launching IronTract tractography...\n'
    printf '  Input MAT file: %s\n' "$mat_file"
    printf '  Injection mask: %s\n' "$injection_file"
    printf '  Output dir: %s\n' "$output_dir"
else
    matlab_command+=" runTractography('$(escape_matlab_string "$mat_file")');"
    printf 'Launching standard tractography...\n'
    printf '  Input MAT file: %s\n' "$mat_file"
fi

nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &

printf 'Logging to: %s\n' "$log_file"
printf 'PID: %d\n' "$!"
printf '\nTo monitor progress: tail -f %s\n' "$log_file"
