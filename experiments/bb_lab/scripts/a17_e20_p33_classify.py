#!/usr/bin/env python3
"""A17 E20 — classify the 555 P-33 terminal patterns.

The abstract enumeration (a17_e20_overlay_engine.py --mode p33) is
COMPLETE: 555 injective patterns survive the frame-free d ∈ {1,2}
kill discipline. This classifier applies the CLASS-frame reductions
(G = Z_ℓ × Z_m, 4∤ℓ, 4∤m — rank 2, 2-adic valuation ≤ 1):

  R1 (2-power kill): a difference-class form forced with d = 2^k
     torsion has order | gcd(2^k, 2·odd) = 2 ⟹ 2-torsion ⟹ D1-dead.
  R2 (torsion reduction): d-torsion ⟹ d̃-torsion, d̃ = 2^{min(v₂,1)}
     · odd(d).
  R3 (confinement cardinality): if ALL FOUR a-units (or b-units)
     are torsion-forced with combined reduced exponent d̃, then
     A ⊆ a₀ + T_d̃ and dA ⊆ T_d̃; a Sidon 5-set needs 20 distinct
     nonzero differences, so |T_d̃| ≥ 21. On rank-2 frames
     |T_d̃| = gcd(d̃,ℓ)·gcd(d̃,m) ≤ d̃², so d̃ ≤ 4 is DEAD, and each
     surviving d̃ pins the frame divisibility (the ZONE).
  R4 (zone Sidon census): for a zone whose torsion subgroup is
     forced ≅ a specific small group T (e.g. Z₅² for d̃ = 5 — the
     only rank-2 option with |T₅| ≥ 21), exhaustively census Sidon
     5-sets in T; an empty census empties the zone.

Output: the family table — per-terminal verdicts and the residue.

Usage: uv run python scripts/a17_e20_p33_classify.py \
    --table data/a17/e20_p33_table.json \
    --out data/a17/e20_p33_families.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a17_e20_overlay_engine import (ZLattice, avec, bvec,
                                    p33_grid_forms, static_battery)


def v2(n: int) -> int:
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def reduce_torsion(d: int) -> int:
    """R2: order | d ⟹ order | d̃ on 4∤ frames."""
    odd = d >> v2(d)
    return odd * (2 if d % 2 == 0 else 1)


def sidon5_census(mods: tuple) -> int:
    """Number of Sidon 5-sets (20 distinct differences) in
    Z_mods[0] × ... (small groups only)."""
    from itertools import product
    elems = list(product(*(range(m) for m in mods)))
    n = len(elems)

    def sub(x, y):
        return tuple((a - b) % m for a, b, m in zip(x, y, mods))

    count = 0
    idx = {e: i for i, e in enumerate(elems)}
    for comb in combinations(range(1, n), 4):
        S = [elems[0 - 0]] + [elems[i] for i in comb]  # anchor 0
        diffs = set()
        ok = True
        for x, y in combinations(S, 2):
            d1, d2 = sub(x, y), sub(y, x)
            if d1 in diffs or d2 in diffs or d1 == d2:
                ok = False
                break
            diffs.add(d1)
            diffs.add(d2)
        if ok:
            count += 1
    return count


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()
    data = json.loads(Path(args.table).read_text())
    terms = data["terminal_patterns"]
    stat = static_battery()

    verdicts = Counter()
    residue = []
    zone_demands = Counter()
    for t in terms:
        lat = ZLattice()
        lat.rows = [r[:] for r in t["lattice"]]
        gforms, gkinds, r, g = p33_grid_forms(t["alpha"], t["beta"])
        allf = stat.forms + gforms
        allk = stat.kinds + gkinds
        # R1: any diff-class form forced at a pure 2-power
        killed = None
        for f, k in zip(allf, allk):
            if k != "diff":
                continue
            d = lat.forced_denom(f)
            if d is not None and (d & (d - 1)) == 0:
                killed = f"R1-2power(d={d})"
                break
        if killed:
            verdicts[killed.split("(")[0]] += 1
            continue
        # unit torsion profile
        da = [lat.forced_denom(avec(i)) for i in range(1, 5)]
        db = [lat.forced_denom(bvec(i)) for i in range(1, 5)]
        zones = []
        for side, ds in (("A", da), ("B", db)):
            if all(d is not None for d in ds):
                dt = reduce_torsion(lcm(*ds))
                if dt <= 4:
                    killed = f"R3-confined(side={side},d~={dt})"
                    break
                zones.append((side, dt))
        if killed:
            verdicts[killed.split("(")[0]] += 1
            continue
        if zones:
            for side, dt in zones:
                zone_demands[dt] += 1
            verdicts["ZONE"] += 1
            residue.append({"kind": "zone",
                            "zones": zones, "rank": t["rank"],
                            "alpha": t["alpha"], "beta": t["beta"],
                            "da": da, "db": db})
        else:
            verdicts["FREE"] += 1
            residue.append({"kind": "free", "rank": t["rank"],
                            "alpha": t["alpha"], "beta": t["beta"],
                            "da": da, "db": db})

    out = {"terminals": len(terms), "verdicts": dict(verdicts),
           "zone_demand_hist": dict(sorted(zone_demands.items()))}

    # R4: Sidon-existence in the small forced zones
    zone_census = {}
    demands = sorted(zone_demands)
    if any(d % 5 == 0 for d in demands):
        zone_census["Z5xZ5"] = sidon5_census((5, 5))
    if any(d % 7 == 0 for d in demands):
        zone_census["Z7xZ7"] = sidon5_census((7, 7))
    out["zone_sidon_census"] = zone_census
    out["residue"] = residue
    out["secs"] = round(time.time() - t0, 1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items()
                      if k != "residue"}, indent=1))
    print(f"\nresidue families: "
          f"{sum(1 for x in residue if x['kind'] == 'zone')} zone + "
          f"{sum(1 for x in residue if x['kind'] == 'free')} free")


if __name__ == "__main__":
    main()
