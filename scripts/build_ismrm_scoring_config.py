#!/usr/bin/env python3
"""build_ismrm_scoring_config.py — Build the merged scoring config for the
ISMRM 2015 (Renauld 2023) scorer.

The scoring package ships two configs:
  - config_file_segmentation.json : per-bundle ROI/length segmentation rules
  - config_file_tractometry.json  : per-bundle ground-truth mask (gt_mask)

`scil_tractogram_segment_with_ROI_and_score` needs BOTH in one file to compute
per-bundle Dice/overlap/overreach. This script merges the gt_mask from the
tractometry config onto each segmentation entry (with the CC sub-bundles and
ICP parts inheriting their parent bundle's mask), writing config_file_merged.json.

Usage:
    python scripts/build_ismrm_scoring_config.py \
        [--scoring-dir data/ismrm2015/scoring_data_Renauld2023]
"""
import argparse
import json
from pathlib import Path


def build(scoring_dir: Path) -> Path:
    seg = json.loads((scoring_dir / "config_file_segmentation.json").read_text())
    tract = json.loads((scoring_dir / "config_file_tractometry.json").read_text())

    matched = 0
    # direct name matches
    for name in seg:
        if name in tract and "gt_mask" in tract[name]:
            seg[name]["gt_mask"] = tract[name]["gt_mask"]
            matched += 1
    # CC sub-bundles inherit the whole-CC mask
    if "CC" in tract and "gt_mask" in tract["CC"]:
        for name in seg:
            if name.startswith("CC") and "gt_mask" not in seg[name]:
                seg[name]["gt_mask"] = tract["CC"]["gt_mask"]
                matched += 1
    # ICP_<side>_partN inherit the whole-side mask
    for side in ("left", "right"):
        whole = f"ICP_{side}"
        if whole in tract and "gt_mask" in tract[whole]:
            for name in seg:
                if name.startswith(whole + "_part") and "gt_mask" not in seg[name]:
                    seg[name]["gt_mask"] = tract[whole]["gt_mask"]
                    matched += 1

    with_mask = sum(1 for b in seg.values() if "gt_mask" in b)
    out = scoring_dir / "config_file_merged.json"
    out.write_text(json.dumps(seg, indent=2))
    print(f"Bundles with gt_mask: {with_mask}/{len(seg)} (matched {matched})")
    print(f"Wrote {out}")
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scoring-dir", default="data/ismrm2015/scoring_data_Renauld2023",
                   help="directory containing the two provided config files")
    args = p.parse_args()
    build(Path(args.scoring_dir))
