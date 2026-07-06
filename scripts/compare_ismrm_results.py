#!/usr/bin/env python3
"""compare_ismrm_results.py — Tabulate scoring results across HINEC configs.

Reads <run_dir>/scoring/renauld2023/results.json from each provided run dir
and prints a side-by-side comparison of headline metrics and per-bundle F1.

Usage:
    python compare_ismrm_results.py <run_dir1> <run_dir2> ...
"""
import json
import sys
from pathlib import Path


def load(run_dir: Path):
    f = run_dir / "scoring" / "renauld2023" / "results.json"
    if not f.exists():
        return None
    with open(f) as fp:
        return json.load(fp)


def label_from_dir(d: Path) -> str:
    # Extract config name from run_<TS>_<config>
    name = d.name
    parts = name.split("_", 2)
    return parts[2] if len(parts) >= 3 else name


def main(run_dirs):
    results = {}
    for d in run_dirs:
        d = Path(d).resolve()
        r = load(d)
        if r is None:
            print(f"[skip] {d}: no results.json", file=sys.stderr)
            continue
        results[label_from_dir(d)] = r

    if not results:
        print("No results to compare.", file=sys.stderr)
        sys.exit(1)

    labels = list(results.keys())
    headline_keys = [
        ("VB", "{:.0f}"),
        ("VS", "{:.0f}"),
        ("VS_ratio", "{:.3f}"),
        ("IS_ratio", "{:.3f}"),
        ("mean_OL", "{:.3f}"),
        ("mean_OR_gt", "{:.3f}"),
        ("mean_OR_vs", "{:.3f}"),
        ("mean_f1", "{:.3f}"),
    ]

    print("\n=== Headline metrics ===")
    print(f"{'metric':<15}" + "".join(f"{l:>18}" for l in labels))
    for k, fmt in headline_keys:
        row = f"{k:<15}"
        for l in labels:
            v = results[l].get(k, "—")
            row += f"{(fmt.format(v) if isinstance(v, (int, float)) else str(v)):>18}"
        print(row)

    # Per-bundle F1 across configs
    bundles = set()
    for r in results.values():
        for b, info in r.get("bundle_wise", {}).items():
            if isinstance(info, dict) and "f1" in info:
                bundles.add(b)
    bundles = sorted(bundles)

    print("\n=== Per-bundle F1 ===")
    print(f"{'bundle':<28}" + "".join(f"{l:>18}" for l in labels))
    for b in bundles:
        row = f"{b:<28}"
        for l in labels:
            info = results[l]["bundle_wise"].get(b, {})
            v = info.get("f1") if isinstance(info, dict) else None
            row += f"{(f'{v:.3f}' if isinstance(v, (int, float)) else '—'):>18}"
        print(row)

    print("\n=== Per-bundle VS (streamline counts) ===")
    print(f"{'bundle':<28}" + "".join(f"{l:>18}" for l in labels))
    for b in bundles:
        row = f"{b:<28}"
        for l in labels:
            info = results[l]["bundle_wise"].get(b, {})
            v = info.get("VS") if isinstance(info, dict) else None
            row += f"{(str(v) if v is not None else '—'):>18}"
        print(row)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
