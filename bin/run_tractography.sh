#!/usr/bin/env bash
# run_tractography.sh — Tractography-ONLY on already-preprocessed HINEC data.
#
# Companion to run_hinec.sh. Split so you preprocess ONCE and iterate tractography
# cheaply:
#   1. run_hinec.sh        → preprocess a dataset once → a run dir with
#                            output/<nim>.mat + intermediate/ (DWI + masks).
#   2. run_tractography.sh → run any number of tractography CONFIGS on that
#                            processed nim WITHOUT re-preprocessing, each into its
#                            own scorable run dir. The preprocessed intermediate/
#                            is REUSED (symlinked) so run_ismrm_scoring.sh works.
#
# Usage (config-driven — the intended path):
#   ./bin/run_tractography.sh <config> [--score] [--source <dir_or_mat>]
#     <config>          a YAML config: 'config/<name>.yml' or just '<name>'.
#     --score           after tracking, run bin/run_ismrm_scoring.sh on the run dir.
#     --source <s>      the preprocessed source: a run_hinec run dir (uses its
#                       output/*.mat + intermediate/) or a processed .mat.
#                       DEFAULT (no --source): data/ismrm2015/ismrm2015.mat, the
#                       canonical nim in the DATA dir. Run dirs hold tractography
#                       OUTPUT; the preprocessed nim is a stable input and is not
#                       read from inside one.
#
# Examples:
#   ./bin/run_tractography.sh hinec_csd --score
#   ./bin/run_tractography.sh config/hinec_dti.yml --source hinec_runs/run_20260627_054713_rk4_step005
#
# Legacy forms (unchanged — no run dir / no scoring):
#   ./bin/run_tractography.sh <input.mat> standard|hinec
#   ./bin/run_tractography.sh <input.mat> IronTract <injection.nii.gz> <output_dir>

set -uo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
escape_matlab_string() { printf "%s" "${1//\'/\'\'}"; }

if [[ $# -lt 1 ]]; then sed -n '2,30p' "$0" >&2; exit 1; fi

# --- separate flags from positionals ---------------------------------------
# Parallel workers. The tracker caps parpool at 8 unless HINEC_MAX_WORKERS says
# otherwise, which leaves most of a large machine idle. Default to half the
# cores (capped at 64) so a co-tenant still has room; override by exporting
# HINEC_MAX_WORKERS before running.
if [[ -z "${HINEC_MAX_WORKERS:-}" ]]; then
    _ncpu=$(nproc 2>/dev/null || echo 8)
    _w=$(( _ncpu / 2 )); (( _w < 1 )) && _w=1; (( _w > 64 )) && _w=64
    export HINEC_MAX_WORKERS="$_w"
fi

POS=(); do_score=false; src_override=""; sets=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --score|-s) do_score=true; shift;;
        --source)   src_override="${2:?--source needs a path}"; shift 2;;
        --set)      sets+=("${2:?--set needs key=value}"); shift 2;;   # DSE param override
        *)          POS+=("$1"); shift;;
    esac
done
a1="${POS[0]:-}"; a2="${POS[1]:-}"

# --- legacy forms: <input.mat> standard|hinec | IronTract ------------------
if [[ "$a2" == "standard" || "$a2" == "hinec" ]]; then
    mcmd="addpath(genpath('.')); runTractography('$(escape_matlab_string "$a1")', '$(escape_matlab_string "$a2")');"
    ts=$(date +"%Y-%m-%d_%H_%M_%S"); log="tractography_${ts}.log"
    nohup matlab -batch "$mcmd" > "$log" 2>&1 &
    printf 'Legacy tractography (%s) launched. Log: %s  PID: %d\n' "$a2" "$log" "$!"; exit 0
fi
if [[ "$a2" == "IronTract" ]]; then
    [[ ${#POS[@]} -eq 4 ]] || { echo "IronTract: $0 <input.mat> IronTract <injection> <output_dir>" >&2; exit 1; }
    mcmd="addpath(genpath('.')); runTractography('$(escape_matlab_string "$a1")', 'IronTract', '$(escape_matlab_string "${POS[2]}")', '$(escape_matlab_string "${POS[3]}")');"
    ts=$(date +"%Y-%m-%d_%H_%M_%S"); log="tractography_${ts}.log"
    nohup matlab -batch "$mcmd" > "$log" 2>&1 &
    printf 'IronTract tractography launched. Log: %s  PID: %d\n' "$log" "$!"; exit 0
fi

# --- config-driven, reuse-preprocessing, scorable --------------------------
# Resolve config (accept 'config/x.yml', 'x.yml', or bare 'x')
cfg="$a1"; [[ "$cfg" != *.yml ]] && cfg="config/${cfg}.yml"
[[ -f "$cfg" ]] || { echo "Error: config not found: $cfg" >&2; echo "Available:" >&2; ls -1 config/*.yml >&2; exit 1; }
cfg_name=$(basename "$cfg" .yml)

# Resolve the preprocessed SOURCE nim. Precedence: --source > 2nd positional >
# DEFAULT = the canonical nim in the DATA dir (NOT inside a hinec_runs/ run dir).
# The preprocessed nim is a clean, stable input; run dirs hold tractography OUTPUT.
DEFAULT_NIM="data/ismrm2015/ismrm2015.mat"
src="${src_override:-$a2}"
[[ -z "$src" ]] && src="$DEFAULT_NIM"

# From <src> derive the nim .mat and the intermediate DWI reference for scoring.
if [[ -f "$src" ]]; then
    src_mat="$src"
    # References live next to the nim (data dir) or in a sibling intermediate/ (run dir).
    d="$(cd "$(dirname "$src")" && pwd)"
    if [[ -d "$d/../intermediate" ]]; then src_inter="$d/../intermediate"; else src_inter="$d"; fi
elif [[ -d "$src" ]]; then
    src_mat=$(ls -t "$src"/output/*.mat 2>/dev/null | head -1 || ls -t "$src"/*.mat 2>/dev/null | head -1 || true)
    src_inter="$src/intermediate"; [[ -d "$src_inter" ]] || src_inter="$src"
else
    echo "Error: source nim not found: $src" >&2
    echo "Preprocess first (./bin/run_hinec.sh writes the nim to the data dir), or pass --source <nim.mat|run_dir>." >&2
    exit 1
fi
[[ -n "${src_mat:-}" && -f "$src_mat" ]] || { echo "Error: no processed .mat at source: $src" >&2; exit 1; }
src_mat="$(cd "$(dirname "$src_mat")" && pwd)/$(basename "$src_mat")"

# DSE overrides. Each --set key=value is handed to nim_config_apply_overrides,
# which resolves it against the schema (full path, tractography-assumed, bare
# leaf, or legacy alias), type-checks it, and rejects unknown keys instead of
# silently assigning a field nothing reads. Every parameter is reachable this
# way, including preprocessing.* and nested paths.
override_mcmd=""; name_suffix=""
if [[ ${#sets[@]} -gt 0 ]]; then
    set_list=""
    for kv in "${sets[@]}"; do
        [[ -n "$set_list" ]] && set_list+=", "
        set_list+="'$(escape_matlab_string "$kv")'"
        key="${kv%%=*}"; val="${kv#*=}"
        name_suffix+="_$(printf '%s%s' "${key##*.}" "$val" | tr -c 'A-Za-z0-9' '_' | sed 's/__*/_/g; s/_$//')"
    done
    override_mcmd="config = nim_config_apply_overrides(config, {${set_list}}); "
fi

# New run dir (named by config [+ overrides], matching run_hinec convention)
timestamp=$(date +"%Y%m%d_%H%M%S")
run_name="run_${timestamp}_${cfg_name}${name_suffix}"
run_dir="hinec_runs/${run_name}"
mkdir -p "${run_dir}/logs" "${run_dir}/intermediate" "${run_dir}/output" "${run_dir}/tractography/diagnostics"
cp "$cfg" "${run_dir}/config.yml"
[[ ${#sets[@]} -gt 0 ]] && printf '%s\n' "${sets[@]}" > "${run_dir}/overrides.txt"

# Reuse preprocessing SAFELY: COPY (not symlink) the DWI reference into the run
# dir. Scoring only needs this one file (its voxel->world affine). Copying freezes
# it, so the run dir is a self-contained snapshot — a later reprocess of the data
# dir can NEVER silently change an old run's scoring reference. Provenance recorded.
base=$(basename "$src_mat" .mat)
ref="$src_inter/${base}.nii.gz"
[[ -e "$ref" ]] || ref="data/ismrm2015/${base}.nii.gz"
[[ -e "$ref" ]] || ref="data/ismrm2015/ismrm2015.nii.gz"
if [[ -e "$ref" ]]; then
    cp -L "$ref" "${run_dir}/intermediate/$(basename "$ref")"   # -L: copy the real file, frozen
    linked=1
else
    linked=0; echo "Warning: no DWI reference found to copy; scoring will fail." >&2
fi
{ echo "source_nim:    $src_mat"
  echo "dwi_reference: $ref (copied, frozen at run time)"
  echo "created:       $(date -Is)"; } > "${run_dir}/intermediate/SOURCE.txt"

log_file="${run_dir}/logs/pipeline.log"
printf '========================================\nHINEC TRACTOGRAPHY-ONLY (config-driven)\n========================================\n'
printf 'Source nim:   %s\nConfig:       %s\nRun dir:      %s\nDWI ref:      copied (frozen)\nOverrides:    %s\nScore after:  %s\n' \
    "$src_mat" "$cfg" "$run_dir" "${sets[*]:-(none)}" "$do_score"
printf '========================================\n'

mcmd="addpath(genpath('.')); "
mcmd+="config = load_config_yaml('$(escape_matlab_string "$cfg")'); "
mcmd+="${override_mcmd}"
mcmd+="run_info = create_run_directory('$(escape_matlab_string "$cfg")', 'run_name', '$(escape_matlab_string "$run_name")'); "
mcmd+="runTractography('$(escape_matlab_string "$src_mat")', config, run_info); "
mcmd+="fprintf('\\n=== tractography done: %s ===\\n', '$(escape_matlab_string "$cfg_name")');"

echo "[$(date -Is)] running tractography -> ${run_dir}"
matlab -batch "$mcmd" > "$log_file" 2>&1
rc=$?
echo "[$(date -Is)] matlab exit ${rc} (log: $log_file)"
[[ $rc -ne 0 ]] && { echo "tractography failed; see $log_file" >&2; exit $rc; }

if $do_score; then
    echo "[$(date -Is)] scoring ${run_dir}"
    ./bin/run_ismrm_scoring.sh "$run_dir" > "${run_dir}/logs/scoring.log" 2>&1; sc=$?
    echo "[$(date -Is)] scoring exit ${sc}"
    [[ -f "${run_dir}/scoring/renauld2023/results.json" ]] && echo "results: ${run_dir}/scoring/renauld2023/results.json"
    exit $sc
fi
echo "done: ${run_dir}   (score with: ./bin/run_ismrm_scoring.sh ${run_dir})"
