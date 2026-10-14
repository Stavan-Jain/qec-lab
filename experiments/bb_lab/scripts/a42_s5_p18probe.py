#!/usr/bin/env python3
"""A42 S5 — falsify-first witness hunt on the period-18 cylinder (the
r = 3 member period): is there a class-nontrivial compact X-cycle of
weight <= 34 (< 36 = 2p; parity excludes 35)?

SOLVER TIER ONLY: a SAT witness would REFUTE Theorem W's r = 3 row at
the conjectured value (and be re-verified end to end); UNSAT / timeout
is an observation, never a floor.  The pure sector is certificate-empty
below 36 (§2.10), so any witness must be MIXED (nonzero barren
content) — exactly the HM_9 (Z_9-fibre) residual.

Run: cd experiments/bb_lab && uv run python scripts/a42_s5_p18probe.py
     [Wcols] [wmax] [time_limit_s]
Output: data/a42/s5_p18probe.json (+ .log)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib.util  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FB = _load("a42_s4_fiber")
DATA = LAB / "data" / "a42"
LOGF = (DATA / "s5_p18probe.log").open("a")


def log(s):
    line = f"[{time.strftime('%H:%M:%S')}] {s}"
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


def main():
    Wc = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    wmax = int(sys.argv[2]) if len(sys.argv) > 2 else 34
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 1800.0
    p = 18
    log(f"p={p} window {Wc} columns, wmax {wmax}, time limit {tl}s")
    pr = FB.FibreProbe(p, Wc, time_limit=tl)
    log(f"window: n={pr.cw.n} vars, H {pr.H.shape}, class functionals "
        f"{pr.L.shape[0]} (dim Z {pr.dimZ}, dim B {pr.dimB})")
    t0 = time.time()
    status, v, dt = pr.solve(wmax, nontrivial=True)
    rec = {"p": p, "Wcols": Wc, "wmax": wmax, "status": status,
           "solve_s": dt, "time_limit": tl}
    if status == "sat":
        w = int(v.sum())
        prof = pr.profile(v)
        cells = [(blk, c, j) for blk in (0, 1) for c in range(Wc)
                 for j in range(p) if v[pr.cw.vid(blk, c, j)]]
        rec.update({"weight": w, "fibre_profile_Z3": prof,
                    "cells": cells})
        log(f"WITNESS weight {w} (verified nontrivial by CylWindow.verify) "
            f"Z3-fibre profile {prof}")
    else:
        log(f"{status} after {dt}s (observation only)")
    (DATA / "s5_p18probe.json").write_text(json.dumps(rec, indent=1))
    log(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
