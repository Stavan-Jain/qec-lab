"""A36 — summarize the triage ledgers into the note's final table."""
import collections
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data" / "a36"


def main() -> None:
    rows = []
    for ax in ("x", "y"):
        p = DATA / f"triage_{ax}_results.jsonl"
        if p.exists():
            for line in p.read_text().splitlines():
                rows.append(json.loads(line))
    # dedupe identical (group, A, B, axis) rows (duplicate-line races
    # append twice; verdicts are deterministic)
    seen = set()
    uniq = []
    for r in rows:
        key = (r["group"], r["A"], r["B"], r["axis"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r)
    by = collections.defaultdict(
        lambda: {"rows": 0, "pass_cells": 0, "cells": 0, "lights": [],
                 "errs": 0, "all_light": 0})
    for r in uniq:
        d = by[(r["group"], r["axis"], r["W"])]
        d["rows"] += 1
        if "error" in r:
            d["errs"] += 1
            continue
        d["pass_cells"] += r.get("n_pass", 0)
        d["cells"] += r.get("n_cells", 0)
        if r.get("all_light"):
            d["all_light"] += 1
        d["lights"].append(r.get("n_lights"))
    for g in sorted(by, key=str):
        d = by[g]
        ls = [x for x in d["lights"] if x is not None]
        print(f"{g[0]:8s} {g[1]} W={g[2]:2d}  rows={d['rows']:3d} "
              f"errs={d['errs']} all_light={d['all_light']:3d} "
              f"pass={d['pass_cells']}/{d['cells']}  "
              f"lights min/max={min(ls) if ls else '-'}"
              f"/{max(ls) if ls else '-'}")


if __name__ == "__main__":
    main()
