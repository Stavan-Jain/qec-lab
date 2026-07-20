#!/usr/bin/env python3
"""Post-classifier for a17_c_matching_table.py residuals: test each
saved OTHER pattern against ALL vertex-bijections pi: T->B and both
signs — pattern is (anti-)translation iff the lattice forces
b_{pi(i)} - s*t_i constant over i. Genuine exotics = none expected."""
import itertools
import json
import sys
sys.path.insert(0, "scripts")
from a17_c_matching_table import (Lattice, T_EDGES, B_EDGES, tvec, bvec,
                                  vsub, vneg, edge_diff_T, edge_diff_B)

raw = json.load(open("data/a17/e17_c_table_raw.json"))
patterns = raw["OTHER_patterns"]
print(f"{len(patterns)} OTHER patterns to reclassify")

genuine = []
tally = {"translation(perm)": 0, "anti-translation(perm)": 0, "EXOTIC": 0}
for pat in patterns:
    lat = Lattice()
    for te_s, (be, sgn) in pat.items():
        te = tuple(int(x) for x in te_s.strip("()").split(","))
        bd = edge_diff_B(tuple(be))
        rel = vsub(edge_diff_T(te), bd if sgn > 0 else vneg(bd))
        lat.add(rel)
    found = None
    for pi in itertools.permutations(range(5)):
        for s in (+1, -1):
            base = vsub(bvec(pi[0]), [s * x for x in tvec(0)])
            if all(lat.contains(vsub(vsub(bvec(pi[i]),
                                          [s * x for x in tvec(i)]), base))
                   for i in range(1, 5)):
                found = ("translation(perm)" if s > 0
                         else "anti-translation(perm)", pi)
                break
        if found:
            break
    if found:
        tally[found[0]] += 1
    else:
        tally["EXOTIC"] += 1
        genuine.append(pat)

print(json.dumps(tally, indent=1))
if genuine:
    print("GENUINE EXOTICS:")
    for g in genuine:
        print(json.dumps(g))
print("(C) TABLE " + ("COMPLETE" if not genuine else "INCOMPLETE"))
