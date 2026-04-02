#!/usr/bin/env bash
# run_visualization.sh - Export tractography visualization figures (headless)
#
# Generates all 8 anatomical viewing angles by default (superior, inferior,
# anterior, posterior, left, right, oblique_left, oblique_right).
# Optionally generates single-mode views (whole, grid, region, sequential).
#
# Usage:
#   ./bin/run_visualization.sh <run_dir> <output_dir> [format] [region] [dpi]
#   ./bin/run_visualization.sh <tracks_file> <nim_file> <output_dir> [format] [region] [dpi]
#
# Arguments:
#   run_dir      Path to a HINEC run directory (e.g., hinec_runs/run_20260330_*)
#   tracks_file  Path to tracks .mat file
#   nim_file     Path to nim .mat file
#   output_dir   Output directory for images (REQUIRED)
#   format       png | pdf | eps (default: png)
#   region       Region ID(s), e.g. '5' or '5,10,15' (default: all tracks)
#   dpi          Export resolution (default: 300)
#
# Examples:
#   # All 8 angles, whole brain
#   ./bin/run_visualization.sh hinec_runs/run_20260330_*/ experiments/tractography_figures/FACT
#
#   # All 8 angles, filtered by regions
#   ./bin/run_visualization.sh hinec_runs/run_20260330_*/ figures/corpus_callosum png '23,24,25,26,27,28'
#
#   # PDF at high DPI
#   ./bin/run_visualization.sh hinec_runs/run_20260330_*/ figures/FACT pdf '' 600
#
#   # Explicit files
#   ./bin/run_visualization.sh tracks.mat nim.mat figures/euler png

set -euo pipefail

# ============================================================
# Usage
# ============================================================
usage() {
    printf 'Usage:\n'
    printf '  %s <run_dir>                 <output_dir> [format] [region] [dpi]\n' "$0"
    printf '  %s <tracks_file> <nim_file>  <output_dir> [format] [region] [dpi]\n' "$0"
    printf '\nGenerates all 8 anatomical viewing angles (superior, inferior, anterior,\n'
    printf 'posterior, left, right, oblique_left, oblique_right).\n'
    printf '\nFormats: png (default), pdf, eps\n'
    printf '\nExamples:\n'
    printf '  %s hinec_runs/run_20260330_*/ experiments/tractography_figures/FACT\n' "$0"
    printf '  %s hinec_runs/run_20260330_*/ figures/FACT png "23,24,25"\n' "$0"
    printf '  %s tracks.mat nim.mat figures/euler\n' "$0"
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
fi

# ============================================================
# Helpers
# ============================================================
escape_matlab_string() {
    local str=$1
    printf "%s" "${str//\'/''}"
}

resolve_run_dir() {
    local candidate="$1"
    # Wildcard expansion: pick the newest match
    if [[ "$candidate" == *"*"* ]]; then
        local match
        match=$(ls -dt $candidate 2>/dev/null | head -1) || true
        if [[ -z "$match" ]]; then
            printf 'Error: No directories matching pattern: %s\n' "$candidate" >&2
            exit 1
        fi
        candidate="$match"
    fi
    if [[ -d "$candidate" ]] && [[ -d "${candidate}/tractography" || -d "${candidate}/output" ]]; then
        printf '%s' "$candidate"
    else
        printf ''
    fi
}

# ============================================================
# Detect input mode: run directory vs explicit files
# ============================================================
input_mode=""
run_dir=""
tracks_file=""
nim_file=""

resolved=$(resolve_run_dir "$1")
if [[ -n "$resolved" ]]; then
    input_mode="run_dir"
    run_dir="$resolved"
    output_dir="$2"
    shift 2
else
    if [[ $# -lt 3 ]]; then
        printf 'Error: "%s" is not a run directory. Provide <tracks_file> <nim_file> <output_dir>.\n' "$1" >&2
        usage
    fi
    tracks_file="$1"
    nim_file="$2"
    output_dir="$3"
    shift 3

    if [[ ! -f "$tracks_file" ]]; then
        printf 'Error: Tracks file not found: %s\n' "$tracks_file" >&2
        exit 1
    fi
    if [[ ! -f "$nim_file" ]]; then
        printf 'Error: Nim file not found: %s\n' "$nim_file" >&2
        exit 1
    fi
    input_mode="explicit"
fi

# ============================================================
# Parse optional arguments
# ============================================================
format="${1:-png}"
region="${2:-}"
dpi="${3:-300}"

# Validate format
case "$format" in
    png|pdf|eps) ;;
    *)
        printf 'Error: Invalid format "%s". Use: png, pdf, eps\n' "$format" >&2
        exit 1
        ;;
esac

# ============================================================
# Build MATLAB command
# ============================================================
region_arg=""
if [[ -n "$region" ]]; then
    region_clean="${region//[\[\]]/}"
    region_arg=", 'region', [${region_clean}]"
fi

matlab_command="addpath(genpath('.')); "

if [[ "$input_mode" == "run_dir" ]]; then
    matlab_command+="visualizeTractographyAngles('$(escape_matlab_string "$run_dir")', "
    matlab_command+="'output_dir', '$(escape_matlab_string "$output_dir")', "
    matlab_command+="'format', '${format}', "
    matlab_command+="'dpi', ${dpi}"
    matlab_command+="${region_arg});"
else
    matlab_command+="visualizeTractographyAngles('$(escape_matlab_string "$tracks_file")', '$(escape_matlab_string "$nim_file")', "
    matlab_command+="'output_dir', '$(escape_matlab_string "$output_dir")', "
    matlab_command+="'format', '${format}', "
    matlab_command+="'dpi', ${dpi}"
    matlab_command+="${region_arg});"
fi

# ============================================================
# Determine log file location
# ============================================================
timestamp=$(date +"%Y-%m-%d_%H_%M_%S")
if [[ "$input_mode" == "run_dir" ]]; then
    mkdir -p "${run_dir}/logs"
    log_file="${run_dir}/logs/visualization_${timestamp}.log"
else
    log_file="visualization_${timestamp}.log"
fi

# ============================================================
# Launch
# ============================================================
printf '\n========================================\n'
printf ' HINEC Visualization Export\n'
printf '========================================\n'
if [[ "$input_mode" == "run_dir" ]]; then
    printf '  Run dir:    %s\n' "$run_dir"
else
    printf '  Tracks:     %s\n' "$tracks_file"
    printf '  Nim:        %s\n' "$nim_file"
fi
printf '  Output dir: %s\n' "$output_dir"
printf '  Format:     %s\n' "$format"
printf '  DPI:        %s\n' "$dpi"
if [[ -n "$region" ]]; then
    printf '  Region:     [%s]\n' "$region"
else
    printf '  Region:     all (whole brain)\n'
fi
printf '  Log:        %s\n' "$log_file"
printf '  Views:      8 angles (superior, inferior, anterior, posterior,\n'
printf '              left, right, oblique_left, oblique_right)\n'
printf '========================================\n\n'

nohup matlab -batch "$matlab_command" > "$log_file" 2>&1 &
pid=$!

printf 'Launched (PID: %d)\n\n' "$pid"
printf 'Monitor:\n'
printf '  tail -f %s\n\n' "$log_file"
printf 'Output:\n'
printf '  %s/tractography_*.%s\n' "$output_dir" "$format"
printf '========================================\n'
