"""A28 — run the general fibering engine over the A17 safe-floor docket.

For every distinct (instance, axis) row of data/a17/docket_decision.jsonl:
enumerate fiber frames, score feasibility, and for the best feasible frame
run the certified safe-floor sweep at the docket floor (2d). Stored SAT
verdicts are used only to ASSERT agreement afterwards:

  * SF-CERTIFIED rows must come out CERTIFIED (engine agreement = the
    docket's CMS UNSAT@2d−2 replaced by an analytic sweep), and
  * UNKNOWN rows are new decisions (certified or refuted-with-witness).

Also records, per row: the (R) sanity check (k̃ = k for the doubling
cover along the row's axis) and the full fiber scoreboard.

Usage:  uv run python scripts/a28_docket_engine.py [--only ID8] [--out PATH]
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

from bb_lab.fibering import (  # noqa: E402
    FiberFrame,
    best_frames,
    kernel_basis,
    safe_floor_certify,
)
from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402

DOCKET = LAB / "data" / "a17" / "docket_decision.jsonl"
OUT = LAB / "data" / "a28" / "docket_engine.jsonl"


def parse_group(label: str) -> AbelianGroup:
    orders = tuple(int(p[1:]) for p in label.split("x"))
    return AbelianGroup(orders)


def k_of(A: Poly, B: Poly) -> int:
    return 2 * kernel_basis(A, B).shape[0]


def double_axis(G: AbelianGroup, axis: int) -> AbelianGroup:
    orders = list(G.orders)
    orders[axis] *= 2
    return AbelianGroup(tuple(orders))


def run_row(row: dict) -> dict:
    G = parse_group(row["group"])
    A = Poly.from_string(row["A"], G)
    B = Poly.from_string(row["B"], G)
    axis = 0 if row["axis"] == "x" else 1
    target = int(row["floor"])
    out: dict = {
        "id": row["instance_id"],
        "axis": row["axis"],
        "group": row["group"],
        "n": row["n"],
        "k": row["k"],
        "d_base": row["d_base"],
        "floor": target,
        "docket_status": row["status"],
    }
    # (R) sanity: k preserved under the axis-doubling cover
    k_base = k_of(A, B)
    Gc = double_axis(G, axis)
    k_cover = k_of(
        Poly.from_support(A.support, Gc), Poly.from_support(B.support, Gc)
    )
    out["k_base"] = k_base
    out["k_cover"] = k_cover
    out["R_holds"] = bool(k_base == k_cover)
    assert k_base == row["k"], (k_base, row["k"])

    scores = best_frames(A, B, target)
    out["frames"] = [
        {
            "z": list(s.z), "q": s.q, "S": s.S, "mode": s.mode,
            "maxact": s.maxact, "n_masks": s.n_masks,
            "r_delta": s.r_delta, "margin": s.margin,
            "link_shift": s.link_shift, "feasible": s.feasible,
            "note": s.note,
        }
        for s in scores
    ]
    feas = [s for s in scores if s.feasible]
    if not feas:
        out["engine_status"] = "NO-FEASIBLE-FIBER"
        return out
    pick = feas[0]
    out["picked_fiber"] = list(pick.z)
    t0 = time.time()
    try:
        rep = safe_floor_certify(A, B, axis=axis, target=target, z=pick.z)
    except AssertionError as e:
        out["engine_status"] = f"ENGINE-ASSERT: {e}"
        return out
    out["sweep_seconds"] = round(time.time() - t0, 1)
    out["per_class"] = [r.summary() for r in rep.per_class]
    out["notes"] = rep.notes
    if rep.certified:
        out["engine_status"] = "CERTIFIED"
    elif rep.refuted:
        out["engine_status"] = "REFUTED"
        wts = sorted(
            {v["weight"] for r in rep.per_class for v in r.violations}
        )
        out["witness_weights"] = wts
        # persist one verified witness for the record
        for r in rep.per_class:
            if r.violations:
                v0 = r.violations[0]
                out["witness_u"] = [int(i) for i in np.flatnonzero(v0["u"])]
                out["witness_v"] = [int(i) for i in np.flatnonzero(v0["v"])]
                break
    else:
        out["engine_status"] = "UNDECIDED"
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", default=None, help="instance id prefix filter")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    rows = [json.loads(l) for l in DOCKET.read_text().splitlines() if l]
    seen = set()
    todo = []
    for r in rows:
        key = (r["instance_id"], r["axis"])
        if key in seen:
            continue
        seen.add(key)
        if args.only and not r["instance_id"].startswith(args.only):
            continue
        todo.append(r)
    print(f"{len(todo)} distinct (id, axis) rows", flush=True)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    results = []
    agree = new_decided = 0
    with outp.open("w") as fh:
        for r in todo:
            t0 = time.time()
            res = run_row(r)
            res["row_seconds"] = round(time.time() - t0, 1)
            results.append(res)
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            tag = ""
            if r["status"] == "SF-CERTIFIED":
                ok = res.get("engine_status") == "CERTIFIED"
                agree += ok
                tag = "AGREE" if ok else "DISAGREE-OR-UNREACHED"
            elif res.get("engine_status") in ("CERTIFIED", "REFUTED"):
                new_decided += 1
                tag = "NEW-DECISION"
            print(
                f"[{res['id'][:8]}:{res['axis']}] {res['group']} "
                f"floor {res['floor']} docket={res['docket_status']} -> "
                f"engine={res.get('engine_status')} "
                f"({res['row_seconds']}s) {tag}",
                flush=True,
            )
    print(
        f"\nagreements on SF-CERTIFIED: {agree}; new decisions: {new_decided}",
        flush=True,
    )


if __name__ == "__main__":
    main()
