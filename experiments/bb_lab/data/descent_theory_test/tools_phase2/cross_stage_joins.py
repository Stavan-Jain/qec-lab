"""Stage-C cross-stage synthesis: join Stage-A certified floors with
verified minimum-weight witnesses (coordinator request, 2026-08-17).

Two provenance families, kept separate:

  A. SWEEP JOIN (batch-2 rows): Stage-A floor F + order-144 sweep
     d_ub = F with `d_ub_witness_verified`.  The sweep stored the bool,
     not the vector, but its sampler is deterministic
     (l1_distance_ub, n_samples=100000, seed=3 — driver.py verbatim), so
     the witness is RE-DERIVED here and re-verified from scratch
     (weight == F, HZ v = 0, v not in rowspace HX via
     verify_witness_in_nontrivial_coset).  d = F EXACT, certificate
     tier, cross-stage sandwich.

  B. CORPUS-UB CANDIDATES (batch-1 rows with floor == corpus d_ub):
     the corpus witness vector is not banked either; a FRESH L1 hunt
     (60k samples over several seeds) looks for a weight-F witness and
     verifies it.  Found -> d = F exact (floor + freshly-verified
     witness).  Not found -> stays a floor (recorded honestly).

Floors are never derived from witnesses; witnesses never tighten floors.
Output: ../phase3_joins.jsonl
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DTT = HERE.parent
LAB = DTT.parent.parent
SWEEP = Path("/Users/stavanjain/Code/qec-lab/.claude/worktrees/"
             "quantum-distance-presentation-73cb27/experiments/bb_lab/"
             "data/order144_sweep/results.jsonl")
sys.path.insert(0, str(LAB / "src"))

from bb_lab.group import AbelianGroup            # noqa: E402
from bb_lab.poly import Poly                      # noqa: E402
from bb_lab.checks import bb_check_matrices       # noqa: E402
from bb_lab.l1_sampling import (                  # noqa: E402
    l1_distance_ub, verify_witness_in_nontrivial_coset)


def load(p) -> list[dict]:
    return [json.loads(line) for line in open(p)]


def main() -> None:
    by = {}
    for r in load(DTT / "phase2_results.jsonl"):
        if r.get("kind") == "closure":
            by[(r["instance_id"], r["q"])] = r
    floors: dict[str, int] = {}
    strat: dict[str, str] = {}
    for (iid, _q), r in by.items():
        e = r["eval"]
        if e["outcome"] == "CERTIFIED_FLOOR" and e.get("floor"):
            floors[iid] = max(floors.get(iid, 0), e["floor"])
            strat[iid] = r.get("stratum")
    sweep = {s["instance_id"]: s for s in load(SWEEP)}
    preds1 = {r["instance_id"]: r for r in load(DTT / "predictions.jsonl")}
    preds2 = {r["instance_id"]: r
              for r in load(DTT / "predictions_batch2.jsonl")}

    out_p = DTT / "phase3_joins.jsonl"
    results = []

    def emit(rec):
        rec["ts"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        results.append(rec)
        with out_p.open("a") as f:
            f.write(json.dumps(rec) + "\n")
        print(json.dumps({k: rec[k] for k in
                          ("instance_id", "family", "floor", "witness_w",
                           "witness_verified", "join")}, sort_keys=True),
              flush=True)

    # ---------------- family A: sweep joins
    for iid, F in sorted(floors.items()):
        s = sweep.get(iid)
        if s is None:
            continue
        p = preds2[iid]
        rec = {"instance_id": iid, "family": "A-sweep-join",
               "stratum": strat[iid], "n": p["n"], "k": p["k"],
               "group": p["group"], "floor": F,
               "floor_leg": "Stage-A counting-invariant tower census "
                            "(phase2_results.jsonl)",
               "sweep_d_ub": s.get("d_ub"),
               "sweep_witness_verified_flag": s.get(
                   "d_ub_witness_verified")}
        if s.get("d_ub") != F:
            rec.update({"join": False, "witness_w": None,
                        "witness_verified": None,
                        "note": f"floor {F} < sweep d_ub {s.get('d_ub')}: "
                                "no join; remains a floor"})
            emit(rec)
            continue
        G = AbelianGroup((p["ell"], p["m"]))
        checks = bb_check_matrices(Poly.from_string(p["A_poly"], G),
                                   Poly.from_string(p["B_poly"], G))
        res = l1_distance_ub(checks, n_samples=100_000, seed=3)
        w = int(res.distance_ub)
        wv = bool(verify_witness_in_nontrivial_coset(checks, res.witness))
        ww = int(res.witness.sum())
        assert ww == w
        if w < F:
            rec.update({"join": False, "witness_w": w,
                        "witness_verified": wv,
                        "ALARM": "witness BELOW certified floor — "
                                 "STOP AND INVESTIGATE"})
        else:
            rec.update({
                "join": bool(w == F and wv),
                "witness_w": w, "witness_verified": wv,
                "witness_support": [int(x) for x in
                                    np.nonzero(res.witness)[0]],
                "witness_leg": "re-derived l1_distance_ub(n=100000, "
                               "seed=3) [sweep driver params], F2 "
                               "re-verified here",
                "d_exact": F if (w == F and wv) else None,
                "trust": "certificate tier (cross-stage sandwich: "
                         "counting-invariant floor + verified witness)"})
        emit(rec)

    # ---------------- family B: batch-1 corpus-ub candidates
    for iid, F in sorted(floors.items()):
        p = preds1.get(iid)
        if p is None or strat[iid] == "anchor":
            continue
        if p.get("d_ub") != F:
            continue
        G = AbelianGroup((p["ell"], p["m"]))
        checks = bb_check_matrices(Poly.from_string(p["A_poly"], G),
                                   Poly.from_string(p["B_poly"], G))
        best_w, best_sup, wv = None, None, False
        for seed in (3, 11, 20260817):
            res = l1_distance_ub(checks, n_samples=60_000, seed=seed)
            w = int(res.distance_ub)
            if best_w is None or w < best_w:
                best_w = w
                best_sup = res.witness.copy()
            if w == F:
                break
        if best_w is not None and best_w < F:
            emit({"instance_id": iid, "family": "B-corpus-ub",
                  "stratum": strat[iid], "floor": F, "witness_w": best_w,
                  "witness_verified": None, "join": False,
                  "ALARM": "fresh witness BELOW certified floor — "
                           "STOP AND INVESTIGATE"})
            continue
        if best_w == F:
            wv = bool(verify_witness_in_nontrivial_coset(checks, best_sup))
        emit({"instance_id": iid, "family": "B-corpus-ub",
              "stratum": strat[iid], "n": p["n"], "k": p["k"],
              "group": p.get("group"), "floor": F,
              "floor_leg": "Stage-A counting-invariant tower census",
              "witness_w": best_w if best_w == F else None,
              "witness_verified": wv if best_w == F else None,
              "witness_support": ([int(x) for x in np.nonzero(best_sup)[0]]
                                  if best_w == F and wv else None),
              "witness_leg": ("fresh l1_distance_ub 60k x <=3 seeds, F2 "
                              "verified here" if best_w == F else
                              f"no weight-{F} witness re-found (best "
                              f"fresh {best_w}); corpus d_ub stands "
                              "un-reverified"),
              "join": bool(best_w == F and wv),
              "d_exact": F if (best_w == F and wv) else None,
              "trust": ("exact (certificate floor + freshly verified "
                        "witness)" if best_w == F and wv else
                        "floor only (certificate tier)")})

    a_joins = sum(1 for r in results
                  if r["family"] == "A-sweep-join" and r.get("join"))
    b_joins = sum(1 for r in results
                  if r["family"] == "B-corpus-ub" and r.get("join"))
    alarms = [r for r in results if "ALARM" in r]
    print(f"\nSUMMARY: A-sweep joins {a_joins}, B-corpus joins {b_joins}, "
          f"no-join rows {len(results)-a_joins-b_joins}, "
          f"ALARMS {len(alarms)}", flush=True)


if __name__ == "__main__":
    main()
