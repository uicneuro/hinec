#!/usr/bin/env bash
# run_assessment.sh — convergence / repeatability / consistency of ANY tracker.
#
# Compares tractograms on four separate tiers (yield, space, matched paths, MMF
# geometry) and writes hinec_runs/assess_<ts>_<study>[_name]/ with report.html,
# summary.json and pairs.csv. See docs/TRACT_ASSESSMENT.md.
#
# A) Assess runs that already exist (any run dir or tracks .mat):
#   ./bin/run_assessment.sh convergence   --param integrator.step --runs <run> <run> ... [--reference <run>]
#   ./bin/run_assessment.sh consistency   --runs <run> <run> ... [--labels a b ...]
#   ./bin/run_assessment.sh repeatability --group subj1=<run>,<run> --group subj2=<run>,<run>
#   ./bin/run_assessment.sh compare <runA> <runB>
#
# B) Generate the runs with run_tractography.sh, then assess them:
#   ./bin/run_assessment.sh convergence --config hinec_dti --sweep integrator.step=0.4,0.2,0.1,0.05,0.025 \
#         [--source <nim>] [--set key=value ...]
#   ./bin/run_assessment.sh consistency --config hinec_dti --sweep interpolation.method=trilinear,cubic,spline
#   ./bin/run_assessment.sh consistency --configs hinec_dti mmf_dti standard_dti [--set seeding.roi=...]
#   Every generated run shares --source and every --set, so the swept knob (or
#   the config) is the only difference. Seeding must be deterministic
#   (seeding.strategy: uniform, the default) for seed-matched comparisons.
#
# Repeatability needs INDEPENDENT repeats — separate acquisitions preprocessed
# with run_hinec.sh. It is not generated here: re-running the same input
# reproduces the same tractogram (seeding.strategy=random uses MATLAB's default
# startup seed, so it is not a source of independent repeats either).
#
# Any other flag (--name, --out, --grid-mm, --max-tracks, --max-pairs,
# --sensitivity, --workers, --image, --finer) is passed to scripts/assess_tractography.py.
#
# Resources: generated trackings use HINEC_MAX_WORKERS=4 MATLAB workers and the
# analysis uses 4 Python processes unless you raise them. Check the seed count
# first (an ROI such as Fornix at the default density is ~88k seeds per run).

set -uo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

if [[ $# -lt 1 ]]; then sed -n '2,31p' "$0" >&2; exit 1; fi

# Python with h5py/nibabel/scipy/matplotlib: HINEC_PYTHON > ISMRM venv > fsl > python3
pick_python() {
    for py in "${HINEC_PYTHON:-}" "${ISMRM_VENV:-$HOME/venvs/ismrm}/bin/python" "$HOME/fsl/bin/python" python3; do
        [[ -n "$py" ]] && command -v "$py" >/dev/null 2>&1 \
            && "$py" -c 'import h5py, nibabel, scipy, matplotlib' >/dev/null 2>&1 && { echo "$py"; return; }
    done
    echo "Error: no Python with h5py, nibabel, scipy, matplotlib (set HINEC_PYTHON)" >&2; exit 1
}
PY=$(pick_python) || exit 1

study="$1"; shift
config=""; configs=(); sweep=""; source_arg=(); sets=(); pass=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)  config="${2:?--config needs a name}"; shift 2;;
        --configs) shift; while [[ $# -gt 0 && "$1" != --* ]]; do configs+=("$1"); shift; done;;
        --sweep)   sweep="${2:?--sweep needs key=v1,v2,...}"; shift 2;;
        --source)  source_arg=(--source "${2:?--source needs a path}"); shift 2;;
        --set)     sets+=(--set "${2:?--set needs key=value}"); shift 2;;
        *)         pass+=("$1"); shift;;
    esac
done

# --- mode A: analysis only ---------------------------------------------------
if [[ -z "$config" && ${#configs[@]} -eq 0 ]]; then
    exec "$PY" scripts/assess_tractography.py "$study" "${pass[@]}"
fi

# --- mode B: generate runs, then analyse --------------------------------------
# run_tractography.sh defaults to HALF THE CORES of MATLAB workers (64 here),
# which overloads the shared machine. A sweep runs several trackings back to
# back, so cap it at 4 unless the caller explicitly asks for more.
export HINEC_MAX_WORKERS="${HINEC_MAX_WORKERS:-4}"
echo "[assess] MATLAB workers per tracking run: ${HINEC_MAX_WORKERS} (set HINEC_MAX_WORKERS to change)" >&2
[[ "$study" == convergence || "$study" == consistency ]] \
    || { echo "Error: --config/--configs generation supports convergence and consistency only" >&2; exit 1; }

track() {   # run_tractography.sh "$@" ; echo the run dir it created
    local out rc
    out=$(./bin/run_tractography.sh "$@" 2>&1 | tee /dev/stderr); rc=$?   # pipefail: tracker status
    [[ $rc -ne 0 ]] && return $rc
    local dir; dir=$(printf '%s\n' "$out" | sed -n 's/^Run dir: *//p' | tail -1)
    [[ -d "$dir" ]] || { echo "Error: could not find run dir in run_tractography output" >&2; return 1; }
    echo "$dir" > "$RUNLIST_ITEM"
}

runs=(); labels=(); values=()
RUNLIST_ITEM=$(mktemp); trap 'rm -f "$RUNLIST_ITEM"' EXIT
if [[ -n "$sweep" ]]; then
    [[ -n "$config" ]] || { echo "Error: --sweep needs --config" >&2; exit 1; }
    key="${sweep%%=*}"; IFS=',' read -r -a vals <<< "${sweep#*=}"
    for v in "${vals[@]}"; do
        echo "[assess] tracking ${config} with ${key}=${v}" >&2
        track "$config" "${source_arg[@]}" "${sets[@]}" --set "${key}=${v}" >/dev/null || exit 1
        runs+=("$(cat "$RUNLIST_ITEM")"); labels+=("${key##*.}=${v}"); values+=("$v")
    done
else
    [[ "$study" == consistency ]] || { echo "Error: convergence needs --config with --sweep" >&2; exit 1; }
    for c in "${configs[@]}"; do
        echo "[assess] tracking ${c}" >&2
        track "$c" "${source_arg[@]}" "${sets[@]}" >/dev/null || exit 1
        runs+=("$(cat "$RUNLIST_ITEM")"); labels+=("$(basename "$c" .yml)")
    done
fi

if [[ "$study" == convergence ]]; then
    exec "$PY" scripts/assess_tractography.py convergence --param "$key" --runs "${runs[@]}" \
        --values "${values[@]}" "${pass[@]}"
fi
exec "$PY" scripts/assess_tractography.py consistency --runs "${runs[@]}" --labels "${labels[@]}" "${pass[@]}"
