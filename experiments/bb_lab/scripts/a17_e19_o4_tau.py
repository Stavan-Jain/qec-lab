#!/usr/bin/env python3
"""A17 E19 — O4 REDONE with the relative-translation degree of
freedom (τ-corrected shift-grid census).

The original O4 (a17_e19_odd_census.py) compared canon triangle
REPRESENTATIVES at a fixed alignment — but the 9-shift demand
s_ic = r_c − g_i ∈ A − B is invariant only under JOINT translation
of (u_L, u_R); the relative offset τ is a free parameter. Corrected
census: for every dA-triangle class Tg, dB-triangle class Tr, and
τ ∈ G, test u_L = Tg, u_R = Tr + τ:
  shift_pass:   all 9 shifts in A − B (witnesses forced, D2
                bijection |A − B| = 25);
  witness_pass: forced witnesses α row- AND column-injective and β
                row- and column-injective (the P-33 demand; the
                α-column/β-row parts are D2-automatic and asserted).
P-33 predicts witness_pass = 0 over the FULL τ population.

Usage: uv run python scripts/a17_e19_o4_tau.py \
    --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl \
    --out data/a17/e19_o4_tau.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a17_e18_k4_census import diffs
from a17_e19_odd_census import triangles
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()

    out = {"members": 0, "class_pairs": 0, "tau_trials": 0,
           "shift_pass": 0, "witness_pass": 0, "auto_violations": 0,
           "examples": []}
    rows = []
    for path in args.members.split(","):
        rows += [json.loads(l) for l in open(path) if '"A"' in l]
    for rmem in rows:
        G = AbelianGroup(tuple(rmem["frame"]))
        A = sorted(Poly.from_string(rmem["A"], G).support)
        B = sorted(Poly.from_string(rmem["B"], G).support)
        dAs, dBs = diffs(frozenset(A), G), diffs(frozenset(B), G)
        S_map = {}
        for a in A:
            for b in B:
                S_map[G.sub(a, b)] = (a, b)   # r − g = a − b
        assert len(S_map) == 25, "D2 bijection"
        triA = triangles(dAs, G)
        triB = triangles(dBs, G)
        elems = list(G)
        for Tg in triA:
            for Tr in triB:
                out["class_pairs"] += 1
                for tau in elems:
                    out["tau_trials"] += 1
                    wit = {}
                    ok = True
                    for i in range(3):
                        for c in range(3):
                            w = S_map.get(G.sub(G.add(Tr[c], tau),
                                                Tg[i]))
                            if w is None:
                                ok = False
                                break
                            wit[(i, c)] = w
                        if not ok:
                            break
                    if not ok:
                        continue
                    out["shift_pass"] += 1
                    # D2-automatic parts (must always hold)
                    for c in range(3):
                        if len({wit[(i, c)][0]
                                for i in range(3)}) < 3:
                            out["auto_violations"] += 1
                    for i in range(3):
                        if len({wit[(i, c)][1]
                                for c in range(3)}) < 3:
                            out["auto_violations"] += 1
                    rows_ok = all(len({wit[(i, c)][0]
                                       for c in range(3)}) == 3
                                  for i in range(3))
                    cols_ok = all(len({wit[(i, c)][1]
                                       for i in range(3)}) == 3
                                  for c in range(3))
                    if rows_ok and cols_ok:
                        out["witness_pass"] += 1
                        if len(out["examples"]) < 5:
                            out["examples"].append({
                                "frame": rmem["frame"],
                                "A": rmem["A"], "B": rmem["B"],
                                "Tg": [list(t) for t in Tg],
                                "Tr": [list(G.add(t, tau))
                                       for t in Tr]})
        out["members"] += 1
    out["secs"] = round(time.time() - t0, 1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    ok = out["witness_pass"] == 0 and out["auto_violations"] == 0
    print(f"\nO4-tau: {'CLEAN' if ok else 'RESIDUALS'} "
          f"(shift_pass={out['shift_pass']} over "
          f"{out['tau_trials']} tau-trials)")


if __name__ == "__main__":
    main()
