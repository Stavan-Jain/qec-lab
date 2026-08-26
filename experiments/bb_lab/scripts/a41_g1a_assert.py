#!/usr/bin/env python3
"""A41 G1a — falsify-first assert: the regenerated w = 5 structural
census must reproduce the Entry 15 banked table EXACTLY (A5_goal2_log
Entry 15, 2026-07-12; the sweep's own output was never banked — data/
is gitignored and the original jsonl lived in the main checkout only).

Reads the per-frame summary JSON lines from the census log and
hard-asserts every recorded field of every banked frame.  Tail frames
(Entry 15's "in flight", never recorded) are reported as NEW.

Usage: uv run python scripts/a41_g1a_assert.py [--log data/a41/t42_census.log]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Entry 15's census table, verbatim (frame -> poolA poolB zdA zdB d2
# kpos members).  Z2/Z3-axis frames are Lemma-E empty and were skipped
# there; the five [41,60] frames re-measured in this session's chunk
# runs are included where Entry 15 lists them.
BANKED = {
    "Z6xZ7":  (792, 408, 192, 168, 0, 0, 0),
    "Z5xZ9":  (2388, 618, 1092, 138, 0, 0, 0),
    "Z7xZ7":  (2052, 2052, 1512, 1512, 144, 0, 0),
    "Z5xZ10": (3424, 544, 416, 96, 0, 0, 0),
    "Z6xZ9":  (3996, 1080, 2832, 708, 1440, 216, 18),
    "Z5xZ11": (8060, 1030, 0, 0, 0, 0, 0),
    "Z6xZ10": (4298, 848, 626, 352, 336, 64, 8),
    "Z7xZ9":  (8838, 5004, 3852, 1962, 34956, 540, 15),
    "Z5xZ13": (20140, 1644, 64, 108, 0, 0, 0),
    "Z6xZ11": (12800, 1800, 0, 0, 0, 0, 0),
    "Z5xZ14": (23872, 1632, 5968, 96, 7008, 0, 0),
    "Z7xZ10": (11760, 4920, 528, 1344, 9984, 0, 0),
    "Z5xZ15": (42300, 2460, 30540, 1500, 691008, 67296, 2103),
    "Z7xZ11": (25650, 8340, 0, 0, 0, 0, 0),
    "Z6xZ13": (31252, 2832, 24, 816, 960, 0, 0),
}
FIELDS = ("poolA", "poolB", "zdA", "zdB", "d2_pairs", "kpos_pairs",
          "members")
BANKED_TOTAL_MEMBERS = 2144  # 2103 + 18 + 8 + 15


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", type=str, default="data/a41/t42_census.log")
    args = ap.parse_args()

    seen = {}
    for line in Path(args.log).read_text().splitlines():
        line = line.strip()
        if not line.startswith("{"):
            continue
        row = json.loads(line)
        if "poolA" in row:
            seen[row["frame"]] = row

    failures = []
    new_frames = []
    for frame, banked in BANKED.items():
        if frame not in seen:
            failures.append(f"{frame}: MISSING from the regenerated census")
            continue
        row = seen[frame]
        for f, want in zip(FIELDS, banked):
            got = row[f]
            if got != want:
                failures.append(f"{frame}.{f}: banked {want}, got {got}")
    for frame, row in seen.items():
        if frame not in BANKED:
            new_frames.append(
                f"{frame}: members={row['members']} (pools {row['poolA']}/"
                f"{row['poolB']}, kpos={row['kpos_pairs']}) — NEW, no "
                f"Entry-15 reference (was 'in flight')")

    banked_members = sum(seen[f]["members"] for f in BANKED if f in seen)
    print(f"banked frames matched: "
          f"{len(BANKED) - sum(1 for x in failures if 'MISSING' in x)}"
          f"/{len(BANKED)}; banked-scope member total: {banked_members} "
          f"(Entry 15: {BANKED_TOTAL_MEMBERS})")
    for nf in new_frames:
        print("NEW:", nf)
    if failures:
        print("\nGATE RED — deviations from the banked table:")
        for f in failures:
            print("  ", f)
        sys.exit(1)
    assert banked_members == BANKED_TOTAL_MEMBERS, (
        f"member total {banked_members} != banked {BANKED_TOTAL_MEMBERS}")
    print("GATE GREEN — every banked Entry-15 census number reproduced "
          "exactly.")


if __name__ == "__main__":
    main()
