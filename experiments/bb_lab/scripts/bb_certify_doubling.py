"""CLI for the doubling-aware certification front-end (bb_lab.doubling_certify).

Feed it a BB code — typically a suspected doubling cover like the
[[360,4,?]] instances — and it detects the base, certifies the template
inputs, and assembles d = 2*d_base with a certificate bundle.  With
--tandem it then runs the composed Tandem lane: the certified floor goes
in as -init-lb, so the solver only has to FIND the witness (the proof
phase is deleted), giving an independent cross-check in minutes.

Examples:
  uv run python scripts/bb_certify_doubling.py Z30xZ6 \
      "1 + y + x" "y^4 + x + x^11*y^2"
  uv run python scripts/bb_certify_doubling.py Z15xZ12 \
      "1 + y + x" "y^4 + x^8*y^2 + x^13" --tandem --budget 2400
Output: verdict JSON to stdout tail + data/certify_runs/<stamp>/verdict.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.doubling_certify import certify  # noqa: E402


def _progress(stage: str, **kw) -> None:
    detail = " ".join(f"{k}={v}" for k, v in kw.items())
    print(f"[{time.strftime('%H:%M:%S')}] {stage} {detail}".rstrip(),
          flush=True)


def _seed_phase_file(verdict: dict, checks, workdir: Path) -> Path | None:
    """Write the certified witness as a -phase-file (pure bias).

    Recomputes the diagonal-lift witness through the fast path, then
    validates it with the SAME verify_witness the solver path uses — so
    a layout mismatch turns into a skipped seed, never an unsound one."""
    import time as _t

    import numpy as np

    from bb_lab.doubling_certify import (
        BaseTools, RungEngine, d_side_exact, detect, witness_lift,
    )
    from bb_lab.group import AbelianGroup
    from bb_lab.maxsat_distance import verify_witness
    from bb_lab.poly import Poly

    orders = tuple(verdict["input"]["orders"])
    G = AbelianGroup(orders)
    A = Poly.from_string(verdict["input"]["A"], G)
    B = Poly.from_string(verdict["input"]["B"], G)
    cands = [c for c in detect(G, A, B) if c.R_holds]
    if not cands or cands[0].shift_A or cands[0].shift_B:
        print("  [seed] skipped (no shift-free candidate)", flush=True)
        return None
    cand = cands[0]
    Gb = AbelianGroup(cand.base_group)
    bt = BaseTools(Gb, Poly.from_string(cand.base_A, Gb),
                   Poly.from_string(cand.base_B, Gb))
    side = d_side_exact(bt, "X", 15, _t.monotonic() + 300, 8, workdir)
    eng = RungEngine(bt, cand.axis, side["d"])
    wit = witness_lift(eng, side["witness"])
    if wit is None:
        print("  [seed] skipped (no lift)", flush=True)
        return None
    for w in (wit, np.concatenate([wit[len(wit) // 2:],
                                   wit[: len(wit) // 2]])):
        try:
            verify_witness(checks, w, int(w.sum()))
        except Exception:
            continue
        p = workdir / "phase_seed.txt"
        p.write_text("".join("1" if b else "0" for b in w))
        print(f"  [seed] witness weight {int(w.sum())} verified; "
              f"phase file {p}", flush=True)
        return p
    print("  [seed] skipped (witness failed verify in both layouts)",
          flush=True)
    return None


def run_tandem_lane(verdict: dict, args, phase_seed: bool = False) -> dict:
    """The composed witness lane: Tandem with the certified -init-lb."""
    from bb_lab.webui import analysis, solver
    from bb_lab.webui.solver import Job

    floor = verdict["distance"]["floor"]
    report, checks = analysis.analyse(
        "x".join(str(o) for o in verdict["input"]["orders"]),
        verdict["input"]["A"], verdict["input"]["B"], lookup_corpus=False,
    )
    backend = solver.probe(args.binary or solver.DEFAULT_TANDEM)
    if not backend.available:
        return {"ran": False, "error": backend.error}
    flags = dict(verdict["tandem"]["suggested_flags"])
    if args.tandem_cpu:
        flags["-cpu-lim"] = args.tandem_cpu
    seeded = False
    if phase_seed:
        pf = _seed_phase_file(verdict, checks,
                              Path(verdict.get("workdir", ".")))
        if pf is not None:
            flags["-phase-file"] = str(pf)
            seeded = True
    job = Job("tandem-crosscheck", "tandem")

    import threading

    def follow() -> None:
        seq = 0
        while not job.finished:
            for ev in job.since(seq, timeout=5.0):
                seq = ev.seq
                if ev.kind == "incumbent":
                    print(f"  tandem incumbent: {ev.payload['cost']}",
                          flush=True)
                elif ev.kind == "stage":
                    print(f"  tandem {ev.payload['stage']}", flush=True)

    t = threading.Thread(target=follow, daemon=True)
    t.start()
    try:
        solver.run_tandem(
            job, checks, backend.path, flags,
            backend=backend, premises=report.premises,
            acknowledged=verdict["tandem"]["acknowledge"],
        )
    except Exception as e:
        return {"ran": True, "error": f"{type(e).__name__}: {e}"}
    t.join(timeout=1)
    res = job.result or {}
    agree = res.get("distance") == verdict["distance"]["value"]
    return {"ran": True, "distance": res.get("distance"),
            "seconds": res.get("solver_seconds"),
            "method": res.get("method"), "agrees": agree,
            "seeded": seeded}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("group", help="e.g. Z30xZ6")
    ap.add_argument("A")
    ap.add_argument("B")
    ap.add_argument("--budget", type=float, default=2400.0)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--tandem", action="store_true",
                    help="after certifying, run the Tandem witness lane")
    ap.add_argument("--tandem-cpu", type=int, default=900)
    ap.add_argument("--phase-seed", action="store_true",
                    help="seed Tandem's phases with the certified witness")
    ap.add_argument("--binary", default=None, help="tandem binary path")
    args = ap.parse_args()

    orders = tuple(int(p[1:]) for p in args.group.split("x"))
    verdict = certify(orders, args.A, args.B, budget_s=args.budget,
                      threads=args.threads, progress=_progress)

    if args.tandem and verdict.get("distance", {}).get("floor"):
        print("[tandem] launching the composed witness lane "
              f"(-init-lb={verdict['distance']['floor']})", flush=True)
        verdict["tandem"]["crosscheck"] = run_tandem_lane(
            verdict, args, phase_seed=args.phase_seed)

    wd = Path(verdict.get("workdir", "."))
    outp = wd / "verdict.json"

    def scrub(o):
        import numpy as np
        if isinstance(o, dict):
            return {k: scrub(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [scrub(v) for v in o]
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.integer,)):
            return int(o)
        return o

    outp.write_text(json.dumps(scrub(verdict), indent=1))
    d = verdict.get("distance", {})
    print(json.dumps({
        "status": verdict.get("status"),
        "distance": {k: d.get(k) for k in
                     ("value", "floor", "upper", "d_base", "statement")},
        "tier": verdict.get("tier"),
        "wall_s": verdict.get("wall_s"),
        "tandem": {k: v for k, v in verdict.get("tandem", {}).items()
                   if k != "note"},
        "verdict_file": str(outp),
    }, indent=1))


if __name__ == "__main__":
    main()
