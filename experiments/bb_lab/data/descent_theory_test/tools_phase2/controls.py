"""Phase-3 negative/soundness controls (PROTOCOL + tasking):

  C-A  planted-logical FOUND (tower-census engine family): run a known
       code at W = its known d; the engine must FIND a logical of that
       weight (counterexample lane); re-run with that support planted and
       assert the planted-detection path fires.
  C-B  known-false floor must NOT certify: "certify d >= known_d + 2"
       must end COUNTEREXAMPLE, never CERTIFIED_FLOOR (the by90-style
       control, run on the gross anchor and pair72).
  C-C  SAT-lane witness sanity: x_distance's witness verified in a
       nontrivial coset (arm-B engine family).

Emits kind=control records into phase2_results.jsonl.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DTT = HERE.parent
LAB = DTT.parent.parent
RESULTS = DTT / "phase2_results.jsonl"
PY = sys.executable
sys.path.insert(0, str(LAB / "src"))

GROSS = {"ell": 12, "m": 6, "A": "x^3 + y + y^2", "B": "y^3 + x + x^2",
         "folds": [[1, 3], [0, 6], [0, 3]], "d_known": 12,
         "name": "gross [[144,12,12]] (anchor 6f9e8285)"}
PAIR72 = {"ell": 6, "m": 6, "A": "x^2 + y + y^3", "B": "1 + x + y^2",
          "folds": [[0, 3], [1, 3]], "d_known": 8,
          "name": "pair72 [[72,4,8]] (A35 docket)"}


def emit(rec: dict) -> None:
    rec["kind"] = "control"
    rec["ts"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with RESULTS.open("a") as f:
        f.write(json.dumps(rec) + "\n")


def run_worker(spec: dict, tag: str) -> dict:
    wd = HERE / "work" / f"ctrl_{tag}"
    wd.mkdir(parents=True, exist_ok=True)
    spec = dict(spec)
    spec.update({"threads": 4, "workdir": str(wd), "budget_s": 600,
                 "tag": tag})
    sp = wd / "spec.json"
    sp.write_text(json.dumps(spec))
    p = subprocess.run(["nice", "-n", "10", PY,
                        str(HERE / "dtt_close_worker.py"), str(sp)],
                       capture_output=True, text=True, timeout=700,
                       cwd=str(LAB))
    for line in p.stdout.splitlines():
        if line.startswith("CLOSE_JSON "):
            return json.loads(line[len("CLOSE_JSON "):])
    return {"outcome": "ERROR", "reason": p.stderr[-300:]}


def main() -> None:
    # ---- C-B on both codes: false floor W = d_known must counterexample
    for code in (GROSS, PAIR72):
        spec = {k: code[k] for k in ("ell", "m", "A", "B", "folds")}
        spec["W"] = code["d_known"]
        r = run_worker(spec, f"falsefloor_{code['ell']}x{code['m']}")
        ok = (r.get("outcome") == "COUNTEREXAMPLE"
              and r["counterexample"]["weight"] == code["d_known"])
        emit({"instance_id": f"control_{code['ell']}x{code['m']}",
              "control": "false-floor-must-not-certify",
              "target": code["name"], "W": spec["W"],
              "outcome": r.get("outcome"),
              "counterexample_weight": (r.get("counterexample") or {}
                                        ).get("weight"),
              "passed": bool(ok),
              "wall_s": r.get("wall_s")})
        print(f"C-B {code['name']}: {r.get('outcome')} "
              f"w={(r.get('counterexample') or {}).get('weight')} "
              f"passed={ok}", flush=True)
        if code is GROSS and ok:
            cex_support = r["counterexample"]["support"]
            # ---- C-A planted: re-run with the found logical planted
            spec2 = dict(spec)
            spec2["planted_support"] = cex_support
            r2 = run_worker(spec2, "planted_gross")
            ok2 = bool(r2.get("planted_found"))
            emit({"instance_id": "control_planted_gross",
                  "control": "planted-logical-found",
                  "target": code["name"], "W": spec["W"],
                  "planted_weight": len(cex_support),
                  "outcome": r2.get("outcome"),
                  "planted_found": r2.get("planted_found"),
                  "passed": ok2, "wall_s": r2.get("wall_s")})
            print(f"C-A planted gross: found={r2.get('planted_found')} "
                  f"passed={ok2}", flush=True)

    # ---- C-C SAT-lane witness verification (arm-B family)
    from bb_lab.group import AbelianGroup
    from bb_lab.poly import Poly
    from bb_lab.checks import bb_check_matrices
    from bb_lab.sat_distance import x_distance
    from bb_lab.l1_sampling import verify_witness_in_nontrivial_coset
    G = AbelianGroup((6, 6))
    checks = bb_check_matrices(Poly.from_string(PAIR72["A"], G),
                               Poly.from_string(PAIR72["B"], G))
    res = x_distance(checks, weight_upper_bound=8)
    wv = bool(verify_witness_in_nontrivial_coset(checks, res.witness))
    ok3 = res.distance == 8 and wv
    emit({"instance_id": "control_sat_pair72",
          "control": "sat-lane-witness-verified",
          "target": PAIR72["name"],
          "d_sat": int(res.distance), "witness_verified": wv,
          "passed": bool(ok3)})
    print(f"C-C sat pair72: d={res.distance} witness_verified={wv} "
          f"passed={ok3}", flush=True)


if __name__ == "__main__":
    main()
