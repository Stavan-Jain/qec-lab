"""A28 — general fibering engine on curated non-docket safe-floor targets.

Targets (safe floor at 2·d_base on the BASE code, per doubling axis):

* a8-z6z14: the A8 open core. Base [[168,12,6]] over Z6xZ14,
  A = 1 + y + x^3*y^3, B = 1 + x + x^2*y^7; x-cover [[336,12,12]] is
  SAT-exact d = 12 (A8 §3) but its confined/safe floor ("every base
  1-cycle in a nonzero Smith class has weight ≥ 12", A8 §4.3) is the
  declared open core. Engine target: floor 12, axis x.
* by90: Bravyi-360 tower bottom [[90,8,8]] over Z30xZ3 (A19 coords),
  A = x^9 + y + y^2, B = 1 + x^25 + x^26; its x-doubling rung freezes
  at 12 < 16 (A14 §13), so floor 16 should REFUTE — a negative control
  with a known answer.

Usage: uv run python scripts/a28_extra_targets.py [--only NAME]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))

from bb_lab.fibering import best_frames, safe_floor_certify  # noqa: E402
from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402

OUT = LAB / "data" / "a28" / "extra_targets.jsonl"

TARGETS = [
    {
        "name": "a8-z6z14-x",
        "orders": (6, 14),
        "A": "1 + y + x^3*y^3",
        "B": "1 + x + x^2*y^7",
        "axis": 0,
        "floor": 12,
        "expect": "open (A8 §4.3)",
    },
    {
        "name": "by90-x",
        "orders": (30, 3),
        "A": "x^9 + y + y^2",
        "B": "1 + x^25 + x^26",
        "axis": 0,
        "floor": 16,
        "expect": "refute (rung freezes at 12, A14 §13)",
    },
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("a") as fh:
        for tgt in TARGETS:
            if args.only and tgt["name"] != args.only:
                continue
            G = AbelianGroup(tgt["orders"])
            A = Poly.from_string(tgt["A"], G)
            B = Poly.from_string(tgt["B"], G)
            print(f"== {tgt['name']} floor {tgt['floor']} "
                  f"(expect: {tgt['expect']}) ==", flush=True)
            scores = best_frames(A, B, tgt["floor"])
            for s in scores:
                print("  ", s, flush=True)
            rec = {"name": tgt["name"], "floor": tgt["floor"],
                   "expect": tgt["expect"],
                   "frames": [
                       {"z": list(s.z), "q": s.q, "S": s.S, "mode": s.mode,
                        "n_masks": s.n_masks, "margin": s.margin,
                        "link_shift": s.link_shift, "feasible": s.feasible}
                       for s in scores]}
            feas = [s for s in scores if s.feasible]
            if not feas:
                rec["engine_status"] = "NO-FEASIBLE-FIBER"
                print("  NO FEASIBLE FIBER", flush=True)
            else:
                t0 = time.time()
                try:
                    rep = safe_floor_certify(
                        A, B, axis=tgt["axis"], target=tgt["floor"]
                    )
                    rec["seconds"] = round(time.time() - t0, 1)
                    rec["per_class"] = [r.summary() for r in rep.per_class]
                    rec["notes"] = rep.notes
                    rec["engine_status"] = (
                        "CERTIFIED" if rep.certified
                        else "REFUTED" if rep.refuted else "UNDECIDED"
                    )
                    if rep.refuted:
                        rec["witness_weights"] = sorted({
                            v["weight"] for r in rep.per_class
                            for v in r.violations
                        })
                    print(rep.summary(), flush=True)
                    print(f"  [{rec['seconds']}s]", flush=True)
                except (AssertionError, ValueError) as e:
                    rec["engine_status"] = f"ERROR: {e}"
                    print("  ERROR:", e, flush=True)
            fh.write(json.dumps(rec) + "\n")
            fh.flush()


if __name__ == "__main__":
    main()
