#!/usr/bin/env bash
# Lightweight wrapper to launch runTractography via MATLAB batch mode with timestamped logging.
#
# Usage:
#   Standard FACT tractography (default):
#     ./run_tractography.sh <input_mat_file>
#     ./run_tractography.sh <input_mat_file> standard
#
#   HINEC high-order tractography:
#     ./run_tractography.sh <input_mat_file> hinec
#
#   IronTract Challenge mode:
#     ./run_tractography.sh <input_mat_file> IronTract <injection_mask> <output_dir>
#
# Examples:
#   ./run_tractography.sh sample_parcellated.mat
#   ./run_tractography.sh sample_parcellated.mat standard
#   ./run_tractography.sh sample_parcellated.mat hinec
#   ./run_tractography.sh ironTract.mat IronTract ironTract/injection.nii.gz ironTract_submissions/

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input_mat_file> [algorithm]" >&2
    echo "       $0 <input_mat_file> IronTract <injection_file> <output_dir>" >&2
    echo "" >&2
    echo "Algorithms: standard (default), hinec" >&2
    echo "" >&2
    echo "Examples:" >&2
    echo "  Standard FACT:  $0 data.mat" >&2
    echo "  HINEC mode:     $0 data.mat hinec" >&2
    echo "  IronTract mode: $0 data.mat IronTract injection.nii.gz submissions/" >&2
    exit 1
fi

mat_file=$1
algorithm="standard"  # Default algorithm
mode=""
injection_file=""
output_dir=""

# Parse arguments
if [[ $# -eq 4 && "$2" == "IronTract" ]]; then
    # IronTract mode
    mode="IronTract"
    injection_file=$3
    output_dir=$4
elif [[ $# -eq 2 ]]; then
    # Algorithm selection
    if [[ "$2" == "standard" || "$2" == "hinec" ]]; then
        algorithm=$2
    else
        echo "Error: Invalid algorithm '$2'. Use 'standard' or 'hinec'." >&2
        exit 1
    fi
elif [[ $# -gt 1 ]]; then
    echo "Error: Invalid arguments." >&2
    echo "Usage: $0 <input_mat_file> [algorithm]" >&2
    echo "       $0 <input_mat_file> IronTract <injection_file> <output_dir>" >&2
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
    matlab_command+=" runTractography('$(escape_matlab_string "$mat_file")', '$(escape_matlab_string "$algorithm")');"
    if [[ "$algorithm" == "hinec" ]]; then
        printf 'Launching HINEC high-order tractography...\n'
        printf '  Algorithm: HINEC (interpolation + RK4 + ACT)\n'
    else
        printf 'Launching standard FACT tractography...\n'
        printf '  Algorithm: FACT\n'
    fi
    printf '  Input MAT file: %s\n' "$mat_file"
fi

nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &

printf 'Logging to: %s\n' "$log_file"
printf 'PID: %d\n' "$!"
printf '\nTo monitor progress: tail -f %s\n' "$log_file"
