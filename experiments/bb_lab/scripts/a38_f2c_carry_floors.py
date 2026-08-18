"""A38 S2, front F2c (time-boxed): analytic carry floors — can fiber
emptiness be certified without MITM enumeration?

The object: the carry system E v0 = RHS b over a shadow b, with the
window |v0 off supp(b)| <= cap.  Banked ground truth: the a32 sector-C
production layers (the m2-hists ARE the enumeration verdicts):

  sectorAC_C18to22.jsonl.gz — 76,954 fibers at caps 2/1/0 (empties
      319/1,733 @ |b|=18, 5,635/10,602 @ 20, 55,555/64,619 @ 22)
  sectorAC_C14to16.jsonl    — the 397-fiber caps-4/3 layer

CANDIDATE FLOOR (the only structure-free one available): every
off-support solution set X satisfies XOR_{j in X} red[j] = rhs_res
(on-support-reduced), so

  m2 >= GREEDY(rhs_res) := min k s.t. the k largest reduced-column
        weights sum to >= wt(rhs_res)   [>= ceil(wt/max_w)]

Certification: GREEDY > cap => fiber empty, no enumeration.

STRUCTURAL NOTE (the honest reduction, F2b language): solutions form
v0p + Z(base), so the fiber minimum m2 is the distance from v0p (off
supp b) to the PUNCTURED base cycle code — a structured coset-distance
problem.  The charter names the self-similarity risk explicitly; this
probe measures only the structure-free bound and records the gap.

Time-box: measure the bound's certification power on ALL banked
sector-C fibers; if it is no stronger than enumeration where
enumeration is already cheap, record the negative with the measured
numbers and stop (charter F2c).

Output: data/a38/f2c_carry_floors.json
"""

from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import (  # noqa: E402
    AxisDeck, TowerCode, reduce_int, rref_ints, v2i,
)

DATA = LAB / "data"
OUT = DATA / "a38"


def bound_for(deck: AxisDeck, beta: np.ndarray,
              E_cols: list[int]) -> tuple[int, int]:
    """(greedy lower bound on m2, wt(rhs_res)) for the fiber over beta."""
    n = deck.base.n
    rhs = (deck.RHS @ beta) % 2
    bsupp = np.nonzero(beta)[0]
    bcols = [E_cols[int(j)] for j in bsupp]
    bb, bp = rref_ints(bcols)
    rhs_res = reduce_int(v2i(rhs), bb, bp)
    if rhs_res == 0:
        return 0, 0
    wt = bin(rhs_res).count("1")
    bmask = v2i(beta)
    red_w = sorted(
        (bin(reduce_int(E_cols[j], bb, bp)).count("1")
         for j in range(n) if not (bmask >> j) & 1),
        reverse=True)
    acc = 0
    for k, w in enumerate(red_w, start=1):
        acc += w
        if acc >= wt:
            return k, wt
    return len(red_w) + 1, wt        # unreachable: certainly empty


def main() -> None:
    t0 = time.monotonic()
    out: dict = {}
    GB = TowerCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TowerCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck = AxisDeck(BY, GB, 0)
    E_cols = [v2i(deck.E[:, j]) for j in range(GB.n)]
    print(f"[{time.monotonic()-t0:5.1f}s] frames built")

    # ------- the full GB stab-census orbit battery (caps 0..8 spanned):
    # ground truth recomputed by enumeration per orbit rep (this battery
    # CONTAINS the banked 397-fiber caps-4/3 layer and the 28 deep
    # fibers); the banked heavy-band rows carry no supports, so their
    # aggregate rates serve as context, not as the per-fiber battery.
    from bb_lab.tower import enumerate_lifts_deep
    reps: dict[int, dict] = {}
    for line in (DATA / "a32" / "gb_census_stab.jsonl").open():
        r = json.loads(line)
        reps.setdefault(r["canon"], r)
    per_band: dict[int, dict] = {}
    n_done = 0
    for r in sorted(reps.values(), key=lambda r: (r["w"], r["canon"])):
        wb = r["w"]
        cap = (22 - wb) // 2
        if cap < 0:
            continue
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[r["support"]] = 1
        k, wt = bound_for(deck, beta, E_cols)
        lifts = enumerate_lifts_deep(deck, beta, cap=min(cap, 8))
        band = per_band.setdefault(wb, {
            "orbit_fibers": 0, "empty": 0, "cap": cap, "bound_hist": {},
            "certified_empty": 0, "wt_rhs_max": 0, "min_m2_hist": {}})
        band["orbit_fibers"] += 1
        band["wt_rhs_max"] = max(band["wt_rhs_max"], wt)
        band["bound_hist"][str(k)] = band["bound_hist"].get(str(k), 0) + 1
        if not lifts:
            band["empty"] += 1
            if k > cap:
                band["certified_empty"] += 1
        else:
            mm = min(lifts.values())
            band["min_m2_hist"][str(mm)] = \
                band["min_m2_hist"].get(str(mm), 0) + 1
            assert k <= mm, \
                f"UNSOUND: bound {k} > enumerated min m2 {mm} (|b|={wb})"
        n_done += 1
        if n_done % 100 == 0:
            print(f"[{time.monotonic()-t0:5.1f}s]   ... {n_done} orbit "
                  f"fibers")
    print(f"[{time.monotonic()-t0:5.1f}s] GB stab-orbit battery "
          f"({n_done} fibers, caps 0-8; soundness asserted against "
          f"every enumerated minimum):")
    tot_e = tot_c = 0
    for wb in sorted(per_band):
        b = per_band[wb]
        rate = (b["certified_empty"] / b["empty"]) if b["empty"] else None
        tot_e += b["empty"]
        tot_c += b["certified_empty"]
        print(f"    |b|={wb} cap={b['cap']}: {b['empty']}/"
              f"{b['orbit_fibers']} empty; bound-certified "
              f"{b['certified_empty']}"
              + (f" ({100*rate:.1f}% of empties)" if rate is not None
                 else "")
              + f"; bound hist {b['bound_hist']}")
        b["certified_rate_of_empties"] = \
            round(rate, 4) if rate is not None else None
    print(f"[{time.monotonic()-t0:5.1f}s] TOTAL: {tot_c}/{tot_e} "
          f"empties certified by the greedy syndrome-weight bound")
    out["gb_orbit_battery"] = {str(k): v for k, v in per_band.items()}
    out["total"] = {"empties": tot_e, "certified": tot_c}

    # banked heavy-band aggregate rates (context; per-element, not
    # per-orbit — read-only reproduction of the S1-gate numbers)
    band_ctx: dict[int, tuple[int, int]] = {}
    with gzip.open(DATA / "a32" / "sectorAC_C18to22.jsonl.gz") as f:
        for line in f:
            r = json.loads(line)
            e, t = band_ctx.get(r["wbeta"], (0, 0))
            band_ctx[r["wbeta"]] = (e + (0 if r.get("m2_hist") else 1),
                                    t + 1)
    assert band_ctx == {18: (319, 1733), 20: (5635, 10602),
                        22: (55555, 64619)}
    out["banked_heavy_context"] = {
        str(w): {"empty": e, "fibers": t}
        for w, (e, t) in sorted(band_ctx.items())}
    print(f"[{time.monotonic()-t0:5.1f}s] banked heavy-band context "
          f"rates reproduced (319/1733, 5635/10602, 55555/64619)")

    out["wall_s"] = round(time.monotonic() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "f2c_carry_floors.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {OUT / 'f2c_carry_floors.json'}")


if __name__ == "__main__":
    main()
