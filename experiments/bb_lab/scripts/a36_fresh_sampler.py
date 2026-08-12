"""A36 Mode 2 — fresh-code construction at a target parameter point.

When every corpus code's presentation orbit at a point is safe-floor
dead (the A36 Mode-1 outcome so far), the remaining supply is NEW codes.
This samples weight-3 pairs on the target frame and pushes them through
the same funnel, with two extra gates because d(base) is unknown:

  G0  k-gate: rank gives k = target (fast, ~ms)          [exact]
  G1  d-probe: budgeted SAT "is there an X-logical of weight
      <= d_target - 2?"  SAT -> kill (d too small).  UNDET/UNSAT ->
      proceed (certify()'s d_base stage is the exact decider).
  S0/T1  the a36 seam screens at floor 2*d_target (presentation cell =
      the sampled pair itself; its orbit can be swept later if the pair
      is close).
  T1.5/certify  via a36_certify_runner (exact, budget-ruled).

Emits surviving candidates as cover specs to
data/a36/fresh_<tag>_candidates.jsonl, ranked by S0.

Usage:
    uv run python scripts/a36_fresh_sampler.py --orders 7,9 --k 12 \
        --d 10 [--samples 20000] [--t1-cap 30] [--seed 1]
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pysat.card import CardEnc, EncType          # noqa: E402
from pysat.formula import CNF, IDPool            # noqa: E402
from pysat.solvers import Cadical195             # noqa: E402

from a36_orbit_screen import AxisProblem, poly_str  # noqa: E402
from bb_lab.checks import bb_check_matrices      # noqa: E402
from bb_lab.diffset_predicates import (          # noqa: E402
    difference_sets_disjoint, is_sidon,
)
from bb_lab.fibering import kernel_basis         # noqa: E402
from bb_lab.group import AbelianGroup            # noqa: E402
from bb_lab.poly import Poly                     # noqa: E402
from bb_lab.sat_distance import _xor_chain, find_logical_z  # noqa: E402

DATA_DIR = LAB_ROOT / "data" / "a36"


def d_probe_light_logical(A: Poly, B: Poly, cap: int,
                          conf_budget: int = 60_000) -> str:
    """Budgeted SAT: exists nontrivial X-logical of weight <= cap?
    Returns SAT / UNSAT / UNDET.  SAT kills the candidate (d < target).

    X-logicals: v in ker H_Z, v not in rowspace(H_X).  Encode ker H_Z
    by XOR rows; exclude stabilizers by requiring nonzero pairing with
    some Z-logical (v . lz = 1 for at least one basis Z-logical) —
    sound and complete: v is a nontrivial X-logical iff it pairs
    nontrivially with H_1^Z."""
    ch = bb_check_matrices(A, B)
    HZ = ch.H_Z.astype(np.uint8)
    LZ = find_logical_z(ch)
    n = HZ.shape[1]
    pool = IDPool()
    qv = [pool.id() for _ in range(n)]
    cnf = CNF()
    for row in HZ:
        idxs = np.flatnonzero(row)
        out = _xor_chain((qv[i] for i in idxs), pool, cnf)
        if out is not None:
            cnf.append([-out])
    pair_lits = []
    for lz in LZ:
        idxs = np.flatnonzero(lz)
        out = _xor_chain((qv[i] for i in idxs), pool, cnf)
        if out is not None:
            pair_lits.append(out)
    cnf.append(pair_lits)  # at least one nontrivial pairing
    cnf.extend(CardEnc.atmost(lits=qv, bound=cap, vpool=pool,
                              encoding=EncType.seqcounter).clauses)
    with Cadical195(bootstrap_with=cnf.clauses) as solver:
        solver.conf_budget(conf_budget)
        res = solver.solve_limited()
    return "UNDET" if res is None else ("SAT" if res else "UNSAT")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", type=str, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--d", type=int, required=True)
    ap.add_argument("--axis", type=str, default="x")
    ap.add_argument("--samples", type=int, default=20000)
    ap.add_argument("--t1-cap", type=int, default=30)
    ap.add_argument("--s0-keep", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    ell, m = (int(t) for t in args.orders.split(","))
    G = AbelianGroup((ell, m))
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    tag = f"Z{ell}xZ{m}_k{args.k}_d{args.d}_{args.axis}"

    # ---- G0: sample weight-3 pairs, k-gate by kernel dimension -------
    elems = [(e, f) for e in range(ell) for f in range(m)]
    kept: list[tuple[list, list]] = []
    seen: set = set()
    trials = 0
    while trials < args.samples:
        trials += 1
        # translation-normalised: first monomial of each poly = 1
        ia = rng.choice(len(elems), size=2, replace=False)
        ib = rng.choice(len(elems), size=2, replace=False)
        Asup = [(0, 0)] + [elems[i] for i in ia]
        Bsup = [(0, 0)] + [elems[i] for i in ib]
        if len(set(Asup)) < 3 or len(set(Bsup)) < 3:
            continue
        key = (frozenset(Asup), frozenset(Bsup))
        if key in seen:
            continue
        seen.add(key)
        A = Poly.from_support(Asup, G)
        B = Poly.from_support(Bsup, G)
        if 2 * kernel_basis(A, B).shape[0] != args.k:
            continue
        kept.append((Asup, Bsup))
    print(f"G0: {trials} trials -> {len(kept)} k={args.k} pairs "
          f"({time.time() - t0:.0f}s)", flush=True)

    # ---- S0 + diffset annotation -------------------------------------
    scored = []
    target = 2 * args.d
    for Asup, Bsup in kept:
        apr = AxisProblem((ell, m), Asup, Bsup,
                          0 if args.axis == "x" else 1, args.d)
        zg = apr.kernel_grids(apr.A0, apr.B0)
        WA = apr.seam_weights(apr.A0, zg)[:, 0]
        WB = apr.seam_weights(apr.B0, zg)[:, 0]
        s0 = int((WA + WB).min())
        if s0 < target:
            continue
        A = Poly.from_support([tuple(t) for t in Asup], G)
        B = Poly.from_support([tuple(t) for t in Bsup], G)
        d1d2 = bool(is_sidon(A) and is_sidon(B)
                    and difference_sets_disjoint(A, B))
        scored.append({"A": poly_str(sorted(Asup)),
                       "B": poly_str(sorted(Bsup)),
                       "s0": s0, "d1d2": d1d2})
    scored.sort(key=lambda r: (-r["s0"], not r["d1d2"]))
    hist = Counter(r["s0"] for r in scored)
    print(f"S0 >= {target}: {len(scored)} pairs; histogram "
          f"{dict(sorted(hist.items()))} ({time.time() - t0:.0f}s)",
          flush=True)
    scored = scored[: args.s0_keep]

    # ---- G1 (d-probe) + T1 on the S0-ranked head ----------------------
    out_path = DATA_DIR / f"fresh_{tag}_candidates.jsonl"
    n_t1 = n_dkill = 0
    kept_final = []
    with out_path.open("w") as fh:
        for r in scored:
            if n_t1 >= args.t1_cap:
                break
            from a36_orbit_screen import parse_poly
            Asup = parse_poly(r["A"])
            Bsup = parse_poly(r["B"])
            A = Poly.from_support([tuple(t) for t in Asup], G)
            B = Poly.from_support([tuple(t) for t in Bsup], G)
            dp = d_probe_light_logical(A, B, args.d - 2)
            if dp == "SAT":
                n_dkill += 1
                continue
            apr = AxisProblem((ell, m), Asup, Bsup,
                              0 if args.axis == "x" else 1, args.d)
            n_t1 += 1
            t1 = apr.t1_cell(apr.A0, apr.B0)
            status = ("SF-PASS" if t1["pass"]
                      else f"kill@{t1.get('kill_weight')}")
            print(f"  [{n_t1}/{args.t1_cap}] s0={r['s0']} "
                  f"d1d2={r['d1d2']} dprobe={dp} {status} "
                  f"({t1['wall_s']}s) A={r['A']} B={r['B']}", flush=True)
            if t1["pass"]:
                rec = {**r, "d_probe": dp, **t1,
                       **apr.cover_spec(apr.A0, apr.B0)}
                kept_final.append(rec)
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
    print(f"DONE: {n_t1} T1 runs, {n_dkill} d-probe kills, "
          f"{len(kept_final)} candidates -> {out_path} "
          f"({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
