#!/usr/bin/env python3
"""ismrm_report_scores.py — Report the ONE number a scored run actually produced.

Reads <run_dir>/scoring/renauld2023/results.json and decides what to print from
how the run was seeded:

  whole brain (seeding.roi empty)  -> the headline metrics (mean_f1, VB, ...)
  bundle      (seeding.roi = B)    -> B's `bundle_wise` row, also written to
                                      <run_dir>/scoring/renauld2023/bundle.txt

`mean_f1` of an ROI-seeded run is meaningless: the scorer averages over all 26
bundles, and only those whose corridor overlaps the seeded ROI can receive any
streamline at all. The seeded bundle's row is the comparable quantity — because
seeds are a deterministic lattice and the scorer credits a streamline to B only
if its whole length lies inside B's all_mask, that row equals the same row of a
whole-brain run with identical settings. Every OTHER row of an ROI run is only a
subset of its whole-brain counterpart, so it is not comparable either.

The seeding ROI is read from the run dir itself: `overrides.txt`
(`seeding.roi=<name>`, written by run_tractography.sh --set) takes precedence
over `tractography.seeding.roi` in `config.yml`.

Usage:
    python scripts/ismrm_report_scores.py <run_dir>
"""
import json
import sys
from pathlib import Path

RESULTS_REL = Path("scoring") / "renauld2023" / "results.json"


def read_overrides_roi(run_dir: Path):
    """seeding.roi=... from overrides.txt (run_tractography.sh --set)."""
    f = run_dir / "overrides.txt"
    if not f.exists():
        return None
    roi = None
    for line in f.read_text().splitlines():
        line = line.strip()
        if not line.startswith("seeding.roi="):
            continue
        roi = line.split("=", 1)[1].strip()
    if roi is None:
        return None
    return parse_roi_value(roi)


def parse_roi_value(raw: str):
    """'[]' / '[UF_left, CA]' / 'UF_left' / 'UF_left,CA' -> list of names."""
    raw = raw.strip()
    if raw.startswith("[") and raw.endswith("]"):
        raw = raw[1:-1]
    return [p.strip().strip("'\"") for p in raw.split(",") if p.strip()]


def read_config_roi(run_dir: Path):
    """tractography.seeding.roi from the run's config.yml.

    A hand-rolled two-level scan rather than a YAML dependency: the configs are
    section -> group -> key with 2-space indent (config/README.md), and this
    script has to run inside the scoring venv, which has no pyyaml.
    """
    f = run_dir / "config.yml"
    if not f.exists():
        return []
    section = group = None
    collecting = False
    roi = []
    for line in f.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())
        stripped = line.strip()
        if collecting:
            # block-list continuation:  roi:\n    - UF_left
            if stripped.startswith("-"):
                roi += parse_roi_value(stripped[1:])
                continue
            collecting = False
        if indent == 0:
            section, group = stripped.split(":", 1)[0], None
            continue
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        key, value = key.strip(), value.split(" #", 1)[0].strip()
        if indent == 2:
            group = key if value == "" else None
            continue
        if indent >= 4 and section == "tractography" and group == "seeding" and key == "roi":
            if value == "":
                collecting = True
            else:
                roi = parse_roi_value(value)
    return roi


def bundle_row(name, row):
    return (
        f"  {name:<22} {row['VS']:>7d} {row['TP']:>8d} {row['FP']:>8d} "
        f"{row['FN']:>8d} {row['OL']:>8.3f} {row['OR_pct_gt']:>8.3f} {row['f1']:>8.3f}"
    )


BUNDLE_HEADER = (
    f"  {'bundle':<22} {'VS':>7} {'TP':>8} {'FP':>8} {'FN':>8} "
    f"{'OL':>8} {'OR_gt':>8} {'f1':>8}"
)


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    run_dir = Path(argv[1]).resolve()
    results = run_dir / RESULTS_REL
    if not results.exists():
        print(f"No {RESULTS_REL} in {run_dir} - nothing to report.", file=sys.stderr)
        return 1
    with open(results) as fp:
        d = json.load(fp)

    roi = read_overrides_roi(run_dir)
    if not roi:
        roi = read_config_roi(run_dir)

    bw = d.get("bundle_wise", {})
    lines = []
    if roi:
        lines.append(f"Seeded in ROI: {', '.join(roi)}")
        lines.append(BUNDLE_HEADER)
        for name in roi:
            # seeding.roi also accepts the scorer's endpoint gates
            # (<bundle>_head / <bundle>_tail); the scored bundle is the stem.
            key = name
            if key not in bw:
                for suffix in ("_head", "_tail"):
                    if key.endswith(suffix) and key[: -len(suffix)] in bw:
                        key = key[: -len(suffix)]
                        break
            if key in bw:
                lines.append(bundle_row(key, bw[key]))
            else:
                lines.append(f"  {name:<22} (not a scored bundle - no bundle_wise row)")
        out = run_dir / RESULTS_REL.parent / "bundle.txt"
        out.write_text("\n".join(lines) + "\n")
        nonzero = sum(1 for v in bw.values() if v.get("VS", 0) > 0)
        print("\n".join(lines))
        print(
            f"NOTE: mean_f1 ({d.get('mean_f1', float('nan')):.4f}) is NOT comparable for an "
            f"ROI run - only {nonzero} of {len(bw)} bundles got any streamline (the rest are "
            "empty by construction). The seeded row above is the number: it equals the same "
            "row of a whole-brain run with identical settings, while every other row here is "
            "only a subset of its whole-brain counterpart."
        )
        print(f"  -> {out}")
    else:
        print("Whole brain (seeding.roi empty)")
        print(
            f"  mean_f1 {d.get('mean_f1', 0):.4f}   mean_OL {d.get('mean_OL', 0):.4f}   "
            f"mean_OR_gt {d.get('mean_OR_gt', 0):.4f}   VB {d.get('VB', 0)}"
        )
        print(
            f"  streamlines {d.get('total_streamlines', 0)}   VS {d.get('VS', 0)}   "
            f"IS {d.get('IS', 0)} ({d.get('IS_ratio', 0):.3f})"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
