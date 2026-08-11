"""A32 Part 9: d(BY) >= 12 solver-free — the calculus applied one rung down.

The d(C) = 24 assembly's last SAT-tier input is d(BY) = 12 (A19 floor@11
SAT ladders), consumed by the b = 0 tau-stratum.  The identical tower
trisection re-derives it at budget 10 (weights are even, Lemma 2):

A nontrivial BY X-logical u with |u| <= 10 would decompose over GB as
u = s0'(beta') + tau'(gamma), |u| = |beta'| + 2 m2, with exactly one of:

  [beta'] != 0 : beta' is a weight-8/10 member of a nonzero GB class
                 (census-complete, Part 5a: 2 + 36 orbits; d(GB) = 8) and
                 m2 <= 1 / 0.  ANY lift is a counterexample.
  beta' = 0    : u = tau'(gamma), gamma a GB logical outside im p_x*,
                 |u| = 2|gamma| >= 16 > 10 — dead (census: min 8).
  beta' stab   : |beta'| in {6, 10} (no weight-8 GB stabilizers), caps
                 m2 <= 2 / 0.  A lift that is NOT a BY-stabilizer is a
                 counterexample.

So: enumerate every fiber; zero logical lifts  =>  d(BY) >= 12.
(The = 12 side is the explicit banked weight-12 vector, re-verified.)

45 fibers total.  Output: data/a32/dby_floor_summary.json.
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

from a30_rung_pass import i2v  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_subclosures import enumerate_lifts  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a32"


def main():
    t0 = time.monotonic()
    GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck_x = TS.Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                     lambda e, s: (e[0] + 15 * s, e[1]))

    # fibers: 38 logical orbits (census 5a) + hexagon + 6 w10 stab orbits
    log_recs = [json.loads(x) for x in
                (DATA / "gb_census_logical.jsonl").open()]
    stab_recs = [json.loads(x) for x in
                 (DATA / "gb_census_stab.jsonl").open()]
    fibers = []  # (support, w, kind, cap)
    seen = set()
    for r in log_recs:
        if r["canon"] not in seen:
            seen.add(r["canon"])
            cap = (10 - r["w"]) // 2
            fibers.append((r["support"], r["w"], "logical", cap))
    seen = set()
    for r in stab_recs:
        if r["w"] in (6, 10) and r["canon"] not in seen:
            seen.add(r["canon"])
            cap = (10 - r["w"]) // 2
            fibers.append((r["support"], r["w"], "stab", cap))
    assert len(fibers) == 45, len(fibers)

    counterexamples = 0
    n_lifts = 0
    stab_lifts = []
    for supp, w, kind, cap in fibers:
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[supp] = 1
        lifts = enumerate_lifts(deck_x, beta, cap=cap)
        for v0c, m2 in lifts.items():
            u = deck_x.lift(i2v(v0c, GB.n), beta)
            wu = int(u.sum())
            assert wu == w + 2 * m2 <= 10 and BY.is_cycle(u)
            n_lifts += 1
            if kind == "logical" or not BY.is_stab(u):
                counterexamples += 1
                print(f"  !! weight-{wu} BY LOGICAL over |beta'|={w} "
                      f"{kind} — d(BY) < 12 ?!")
            else:
                stab_lifts.append((w, m2, wu))
    # beta' = 0 branch: analytic (2|gamma| >= 16); nothing to run.
    # witness side: the banked weight-12 logical
    b12 = json.loads(next(iter(
        (MAIN / "a24" / "band12_census.jsonl").open())))
    wvec = np.zeros(BY.n, dtype=np.uint8)
    wvec[b12["b_support"]] = 1
    assert BY.is_cycle(wvec) and not BY.is_stab(wvec) \
        and int(wvec.sum()) == 12
    out = {"fibers": len(fibers), "lifts": n_lifts,
           "counterexamples": counterexamples,
           "stab_lifts": stab_lifts,
           "witness12_verified": True,
           "wall_s": round(time.monotonic() - t0, 1)}
    verdict = "d(BY) = 12 SOLVER-FREE" if counterexamples == 0 else \
        "REFUTED — d(BY) < 12 (investigate!)"
    print(f"[{out['wall_s']}s] {verdict}: {len(fibers)} fibers, "
          f"{n_lifts} lifts ({len(stab_lifts)} stabilizer, "
          f"{counterexamples} logical), witness-12 verified")
    (DATA / "dby_floor_summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
