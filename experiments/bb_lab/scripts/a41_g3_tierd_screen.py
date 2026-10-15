#!/usr/bin/env python3
"""A41 G3 — Tier-D screen over the w = 5 class population (S1 gate).

For each census member and each doubling axis, run the A14 forward
battery's FREE tiers (k-gate + S0/S1+/S2 via
`a17_corpus_battery.process_cell(run_s4=False)`) on the literal-lift
Z2 cover at floor 2*(2w) = 20.  No SAT anywhere; verdicts are
CHEAP-PASS (SF@20 survives the cheap tiers -> S4-eligible),
CHEAP-REJECT (an explicit coset element < 20: SF@20 false for this
presentation), or K-GATE-FAIL ((R) fails on that axis).

Claim discipline: the floor is 2*(class value) = 20, NOT 2*d(member)
— member exact distances are not individually certified (Entry 15's
450@10 are witness-certified upper bounds + the conjectured class
floor).  A verdict here is an SF@20 statement for the stored
presentation, never a "doubles"/"does not double" claim (A11:
presentation-sensitivity; SF sufficient-not-necessary).

Usage (from experiments/bb_lab):
    uv run python scripts/a41_g3_tierd_screen.py --limit 20   # pricing
    uv run python scripts/a41_g3_tierd_screen.py --jobs 6
Output: data/a41/g3_tierd_screen.jsonl (+ summary on stdout).
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IN = ROOT / "data/a41/t42_w5_census.jsonl"
DEFAULT_OUT = ROOT / "data/a41/g3_tierd_screen.jsonl"
CLASS_D = 10  # 2w at w = 5; floor = 2*CLASS_D = 20


def support_of(poly_str: str, ell: int, m: int):
    from bb_lab.group import AbelianGroup
    from bb_lab.poly import Poly
    G = AbelianGroup((ell, m))
    return [tuple(g) for g in Poly.from_string(poly_str, G).support]


def run_cell(args):
    member_id, A_str, B_str, ell, m, k, axis = args
    from a17_corpus_battery import process_cell
    row = {"instance_id": member_id, "code_id": f"w5:{member_id}",
           "group_struct": f"Z{ell}xZ{m}", "ell": ell, "m": m,
           "n": 2 * ell * m, "k": k, "d_exact": CLASS_D,
           "A_poly": A_str, "B_poly": B_str}
    # process_cell parses A_poly/B_poly via a14's parse_poly; feed it
    # support-list literals instead by monkey-friendly path: the a14
    # parser accepts the canonical "x^i*y^j" sums the census stores.
    t0 = time.time()
    try:
        out = process_cell(row, axis, conf_budget=0, run_s4=False,
                           run_ladder=False)
    except Exception as e:  # surface, never silently drop a cell
        out = {"instance_id": member_id, "axis": axis,
               "status": "ERROR", "error": repr(e)}
    out["member_id"] = member_id
    out["axis_order"] = ell if axis == "x" else m
    out["axis_odd"] = (out["axis_order"] % 2 == 1)
    out["secs"] = round(time.time() - t0, 2)
    # drop the bulky per-class detail from the row (keep counts)
    pc = out.pop("per_class", None)
    if pc:
        out["n_classes"] = len(pc)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", type=str, default=str(DEFAULT_IN))
    ap.add_argument("--frames", type=str, default=None,
                    help="comma list like 5x15 to filter; default all")
    ap.add_argument("--limit", type=int, default=None,
                    help="first N members (pricing runs)")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", type=str, default=str(DEFAULT_OUT))
    args = ap.parse_args()

    want = None
    if args.frames:
        want = set()
        for tok in args.frames.split(","):
            a, b = tok.lower().split("x")
            want.add((int(a), int(b)))

    members = []
    with open(args.census) as fh:
        for line in fh:
            row = json.loads(line)
            ell, m = row["frame"]
            if want and (ell, m) not in want:
                continue
            members.append(row)
    if args.limit:
        members = members[: args.limit]

    # resume: skip (member_id, axis) cells already in the output
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if outp.exists():
        with open(outp) as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                    done.add((r.get("member_id"), r.get("axis")))
                except json.JSONDecodeError:
                    pass

    cells = []
    for i, row in enumerate(members):
        ell, m = row["frame"]
        mid = f"Z{ell}xZ{m}#{i}:{row['A']}|{row['B']}"
        for axis in ("x", "y"):
            if (mid, axis) in done:
                continue
            cells.append((mid, row["A"], row["B"], ell, m, row["k"], axis))

    print(f"members: {len(members)}, cells to run: {len(cells)} "
          f"(resumed past {len(done)})", flush=True)
    t0 = time.time()
    results = []

    def emit(r):
        results.append(r)
        with open(outp, "a") as fh:
            fh.write(json.dumps(r) + "\n")

    if args.jobs > 1 and len(cells) > 1:
        ctx = mp.get_context("spawn")
        with ctx.Pool(args.jobs) as pool:
            for r in pool.imap_unordered(run_cell, cells):
                emit(r)
                if len(results) % 25 == 0:
                    print(f"  {len(results)}/{len(cells)} "
                          f"({time.time() - t0:.0f}s)", flush=True)
    else:
        for c in cells:
            emit(run_cell(c))

    hist = Counter((r.get("status"), r.get("axis_odd")) for r in results)
    by_frame = Counter((r.get("group", "?"), r.get("axis"),
                        r.get("status")) for r in results)
    print("\nverdict histogram (status, axis_odd):")
    for kk, v in sorted(hist.items(), key=lambda kv: str(kv[0])):
        print(f"  {kk}: {v}")
    print("by frame/axis:")
    for kk, v in sorted(by_frame.items(), key=lambda kv: str(kv[0])):
        print(f"  {kk}: {v}")
    print(f"total {len(results)} cells in {time.time() - t0:.0f}s "
          f"-> {outp}")


if __name__ == "__main__":
    main()
