#!/usr/bin/env bash
# run_sift.sh — SIFT-filter an existing tractography run's tractogram (nim_sift),
# writing a new scorable run dir. No re-tracking: it prunes streamlines so their
# density matches the FOD fiber density, cutting overreach (OR) for precision.
#
# Usage:
#   ./bin/run_sift.sh <run_dir> [--score] [--set sift_key=val ...] [--csd <csd.mat>]
#     <run_dir>   an existing tractography run dir (with tractography/tracks*.mat).
#     --score     score the filtered result with bin/run_ismrm_scoring.sh.
#     --set k=v   nim_sift options (sift_batch_frac, sift_n_iter, sift_min_keep).
#     --csd <f>   FOD source (peak_w). Default: data/ismrm2015/ismrm2015_csd.mat.
#
# Requires a CSD FOD field (nim.peak_w). Output: hinec_runs/<run>_sift/.

set -uo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$REPO_ROOT"
escape_matlab_string() { printf "%s" "${1//\'/\'\'}"; }

[[ $# -ge 1 ]] || { sed -n '2,15p' "$0" >&2; exit 1; }
POS=(); do_score=false; sets=(); csd="data/ismrm2015/ismrm2015_csd.mat"; nim="data/ismrm2015/ismrm2015.mat"
src_rundir=""
while [[ $# -gt 0 ]]; do case "$1" in
    --score|-s) do_score=true; shift;;
    --set) sets+=("${2:?}"); shift 2;;
    --csd) csd="${2:?}"; shift 2;;
    --nim) nim="${2:?}"; shift 2;;
    *) POS+=("$1"); shift;;
esac; done
src_rundir="${POS[0]:?need a run dir}"

tracks_mat=$(ls -t "$src_rundir"/tractography/tracks*.mat 2>/dev/null | head -1 || true)
[[ -n "$tracks_mat" ]] || { echo "Error: no tracks*.mat in $src_rundir/tractography/" >&2; exit 1; }
[[ -f "$nim" ]] || { echo "Error: nim not found: $nim (pass --nim)" >&2; exit 1; }
[[ -f "$csd" ]] || { echo "Error: CSD FOD not found: $csd (pass --csd)" >&2; exit 1; }

# override snippet + run-name tag
override=""; suffix=""
for kv in ${sets[@]+"${sets[@]}"}; do
    key="${kv%%=*}"; val="${kv#*=}"
    if [[ "$val" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then mv="$val"; else mv="'$(escape_matlab_string "$val")'"; fi
    override+="opts.$(escape_matlab_string "$key") = $mv; "; suffix+="_$(printf '%s%s' "$key" "$val" | tr -c 'A-Za-z0-9' '_' | sed 's/__*/_/g; s/_$//')"
done

base_run=$(basename "$src_rundir")
timestamp=$(date +"%Y%m%d_%H%M%S")
run_name="run_${timestamp}_${base_run#run_*_}_sift${suffix}"
run_dir="hinec_runs/${run_name}"
mkdir -p "${run_dir}/logs" "${run_dir}/intermediate" "${run_dir}/tractography/diagnostics"

# freeze the DWI reference (copy) for self-contained scoring
ref=$(ls "$src_rundir"/intermediate/*.nii.gz 2>/dev/null | grep -v _M | grep -v parcellation | grep -v mask | head -1 || true)
[[ -n "$ref" ]] || ref="data/ismrm2015/ismrm2015.nii.gz"
[[ -e "$ref" ]] && cp -L "$ref" "${run_dir}/intermediate/$(basename "$ref")"
{ echo "sift_of:  $src_rundir"; echo "tracks:   $tracks_mat"; echo "csd_fod:  $csd"; echo "created:  $(date -Is)"; } > "${run_dir}/intermediate/SOURCE.txt"

ts=$(date +"%Y-%m-%d_%H_%M_%S")
mcmd="addpath(genpath('.')); "
mcmd+="L = load('$(escape_matlab_string "$tracks_mat")', 'tracks'); tracks = L.tracks; "
mcmd+="S = load('$(escape_matlab_string "$nim")', 'nim'); nim = S.nim; "
mcmd+="C = load('$(escape_matlab_string "$csd")'); nim.peak_w = C.peak_w; nim.peaks = C.peaks; nim.npeaks = C.npeaks; "
mcmd+="opts = struct(); ${override}"
mcmd+="[tracks, sift_info] = nim_sift(tracks, nim, opts); "
mcmd+="options = opts; algorithm = 'sift'; elapsed_time = 0; "
mcmd+="save('${run_dir}/tractography/tracks_sift_${ts}.mat', 'tracks', 'sift_info', 'options', 'algorithm', 'elapsed_time', '-v7.3'); "
mcmd+="fprintf('\\nSIFT wrote %d streamlines -> ${run_dir}\\n', numel(tracks));"

echo "SIFT-filtering ${src_rundir} -> ${run_dir}"
matlab -batch "$mcmd" > "${run_dir}/logs/sift.log" 2>&1
rc=$?
echo "matlab exit ${rc} (log: ${run_dir}/logs/sift.log)"
[[ $rc -ne 0 ]] && { echo "SIFT failed; see log" >&2; tail -5 "${run_dir}/logs/sift.log" >&2; exit $rc; }
grep -a "SIFT:" "${run_dir}/logs/sift.log" || true

if $do_score; then
    ./bin/run_ismrm_scoring.sh "$run_dir" > "${run_dir}/logs/scoring.log" 2>&1
    echo "scored: ${run_dir}/scoring/renauld2023/results.json"
fi
echo "done: ${run_dir}"
