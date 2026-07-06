#!/usr/bin/env bash
# run_ismrm_scoring.sh - Score a HINEC tractography run against ISMRM 2015 ground truth
#
# Converts a HINEC tracks .mat to TRK (RAS world coords), runs scilpy's Renauld
# 2023 scorer against data/ismrm2015/scoring_data_Renauld2023/, and (if the
# original 2015 scoring data is also present) runs the dedicated 2015 scorer
# for a cross-check. Writes everything under <run_dir>/scoring/.
#
# Usage:
#   ./bin/run_ismrm_scoring.sh <hinec_run_dir>
#
# Example:
#   ./bin/run_ismrm_scoring.sh hinec_runs/run_20260605_*_ismrm2015/
#
# Assumes:
#   - ~/venvs/ismrm exists with scilpy + nibabel + ismrm_2015_tractography_challenge_scoring
#   - data/ismrm2015/scoring_data_Renauld2023/ is populated (gt_config + bundles + ROI + t1)
#   - data/ismrm2015/scoring_data_2015/ is populated (optional - 2015 dedicated scorer)

set -uo pipefail
# (intentionally not -e: we want Step 3 to run even if Step 2 fails)

VENV="${ISMRM_VENV:-$HOME/venvs/ismrm}"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCORING_RENAULD="${REPO_ROOT}/data/ismrm2015/scoring_data_Renauld2023"
SCORING_2015="${REPO_ROOT}/data/ismrm2015/scoring_data_2015"
DEDICATED_SCORER="${REPO_ROOT}/../ismrm_2015_tractography_challenge_scoring/scripts/score_tractogram.py"
# Fall back to the install path if the cloned repo moved
[[ -f "$DEDICATED_SCORER" ]] || DEDICATED_SCORER="$HOME/ismrm_2015_tractography_challenge_scoring/scripts/score_tractogram.py"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <hinec_run_dir>" >&2
    echo "" >&2
    echo "Scores a HINEC run against ISMRM 2015 ground truth (Renauld 2023 method," >&2
    echo "plus the dedicated 2015 scorer if scoring_data_2015/ is present)." >&2
    exit 1
fi

RUN_DIR="$(cd "$1" && pwd)"
[[ -d "$RUN_DIR" ]] || { echo "Error: not a directory: $1" >&2; exit 1; }

# Locate the latest tracks .mat in the run's tractography/ dir
TRACKS_MAT=$(ls -t "$RUN_DIR"/tractography/tracks_*.mat 2>/dev/null | head -1 || true)
[[ -n "$TRACKS_MAT" ]] || { echo "Error: no tracks_*.mat in $RUN_DIR/tractography/" >&2; exit 1; }

# Locate the preprocessed DWI NIfTI (its affine is what hinec_to_trk uses)
DWI_NII=$(ls "$RUN_DIR"/intermediate/*.nii.gz 2>/dev/null | grep -v _M.nii.gz | grep -v parcellation | head -1 || true)
[[ -n "$DWI_NII" ]] || { echo "Error: no preprocessed DWI in $RUN_DIR/intermediate/" >&2; exit 1; }

SCORING_DIR="$RUN_DIR/scoring"
mkdir -p "$SCORING_DIR"
TRACKS_TRK="$SCORING_DIR/tracks.trk"

echo "========================================"
echo "ISMRM 2015 Scoring"
echo "========================================"
echo "Run dir:    $RUN_DIR"
echo "Tracks:     $TRACKS_MAT"
echo "DWI ref:    $DWI_NII"
echo "Output:     $SCORING_DIR"
echo "Venv:       $VENV"
echo "========================================"

# Step 1: Convert HINEC .mat to TRK in RAS world coordinates. Use the DWI
# for voxel->world (correct anatomical placement) but attach the scoring T1
# as the saved TRK reference so scilpy's space-compatibility check against
# the ROI masks (which live in scoring T1 space) passes.
echo ""
echo "[1/3] Converting tracks to TRK (RAS world space)..."
SCORING_T1="$SCORING_RENAULD/t1.nii.gz"
TRK_REF_ARGS=()
if [[ -f "$SCORING_T1" ]]; then
    TRK_REF_ARGS=(--trk-ref "$SCORING_T1")
fi
"$VENV/bin/python" "$REPO_ROOT/scripts/hinec_to_trk.py" \
    "$TRACKS_MAT" "$DWI_NII" "$TRACKS_TRK" --space rasmm "${TRK_REF_ARGS[@]}"

# Step 2: Renauld 2023 scoring (scilpy)
echo ""
echo "[2/3] Renauld 2023 scoring (scilpy)..."
if [[ -d "$SCORING_RENAULD" ]] && [[ -f "$SCORING_RENAULD/config_file_segmentation.json" ]]; then
    RENAULD_OUT="$SCORING_DIR/renauld2023"
    rm -rf "$RENAULD_OUT"
    # Prefer the merged config (segmentation rules + per-bundle gt_mask) so we
    # get geometric dice/OL/OR/F1 scores per bundle, not just VB/VS counts.
    RENAULD_CFG="$SCORING_RENAULD/config_file_merged.json"
    [[ -f "$RENAULD_CFG" ]] || RENAULD_CFG="$SCORING_RENAULD/config_file_segmentation.json"
    "$VENV/bin/scil_tractogram_segment_with_ROI_and_score" \
        "$TRACKS_TRK" \
        "$RENAULD_CFG" \
        "$RENAULD_OUT" \
        --gt_dir "$SCORING_RENAULD" \
        --reference "$SCORING_RENAULD/t1.nii.gz" \
        --no_bbox_check -f
    echo "  → $RENAULD_OUT/results.json"
else
    echo "  Skipped: $SCORING_RENAULD missing or incomplete."
fi

# Step 3: Original 2015 dedicated scorer (optional)
echo ""
echo "[3/3] Original 2015 dedicated scorer..."
if [[ -d "$SCORING_2015" ]] && [[ -f "$SCORING_2015/gt_bundles_attributes.json" ]] && [[ -f "$DEDICATED_SCORER" ]]; then
    LEGACY_OUT="$SCORING_DIR/dedicated2015"
    rm -rf "$LEGACY_OUT"
    mkdir -p "$LEGACY_OUT"
    "$VENV/bin/python" "$DEDICATED_SCORER" \
        "$TRACKS_TRK" "$SCORING_2015" "$LEGACY_OUT" \
        --compute_ic_ib -f -v
    echo "  → $LEGACY_OUT/scores/"
else
    echo "  Skipped: original 2015 scoring data or dedicated scorer not found."
    echo "  (Expected $SCORING_2015/gt_bundles_attributes.json and $DEDICATED_SCORER)"
fi

echo ""
echo "========================================"
echo "Scoring complete"
echo "========================================"
echo "Renauld 2023 metrics: $SCORING_DIR/renauld2023/results.json"
[[ -d "$SCORING_DIR/dedicated2015" ]] && echo "2015 dedicated metrics: $SCORING_DIR/dedicated2015/scores/"
