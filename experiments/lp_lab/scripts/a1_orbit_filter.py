"""Post-filter sweep JSONL: rebuild each promising a=2 code, compute the number of
X-check orbits of Aut(Tanner graph) with nauty, and report the best (k, d_ub) among
codes with ≥ 2 orbits (certified non-2BGA) vs 1 orbit (possibly 2BGA)."""
import json, sys, glob
from collections import defaultdict
from lp_lab.fastlp import Ab, lp_ab_rows, lp_AB_rows, to_code
from lp_lab.structure import tanner_automorphisms

files = sys.argv[1:-1]; dmin = int(sys.argv[-1])
for fn in files:
    G = None; best = defaultdict(dict); counts = defaultdict(int)
    for line in open(fn):
        r = json.loads(line)
        if r["d_ub"] < dmin: continue
        orders = tuple(r["orders"]) if "orders" in r else (r["ell"],)
        if G is None or G.orders != orders: G = Ab(orders)
        mode = r.get("mode", "ab2")
        if mode == "ab2":
            A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["q"])}
            rX, rZ = lp_ab_rows(G, A, 2, tuple(r["b"]))
        elif mode == "AB2":
            A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["q"])}
            B = {(0, 0): tuple(r["r"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["s"])}
            rX, rZ = lp_AB_rows(G, A, B, 2)
        elif mode == "ab3c":
            A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 1): tuple(r["q"]), (1, 2): (0,), (2, 2): tuple(r["r"]), (2, 0): (0,)}
            rX, rZ = lp_ab_rows(G, A, 3, tuple(r["b"]))
        elif mode == "ab3":
            e = r["e"]
            A = {(0, 0): (0,), (0, 1): (0,), (0, 2): (0,), (1, 0): (0,), (2, 0): (0,), (1, 1): (e[0],), (1, 2): (e[1],), (2, 1): (e[2],), (2, 2): (e[3],)}
            rX, rZ = lp_ab_rows(G, A, 3, tuple(r["b"]))
        code = to_code(rX, rZ, r["n"])
        aut = tanner_automorphisms(code.H_X, code.H_Z)
        ox = aut["orbits_X"]; cls = "non2BGA" if ox >= 2 else "transitive"
        counts[cls] += 1
        key = (cls, r["k"])
        if key not in best or r["d_ub"] > best[key]["d_ub"] or (r["d_ub"] == best[key]["d_ub"] and aut["order"] < best[key]["aut"]):
            best[key] = dict(r, aut=aut["order"], orbitsX=ox)
    print(f"== {fn}: {dict(counts)}")
    for key in sorted(best, key=lambda t: (t[0], -t[1])):
        r = best[key]
        desc = {k: v for k, v in r.items() if k in ("p", "q", "r", "s", "b", "e")}
        print(f"   {key[0]:10s} k={key[1]:3d} d_ub={r['d_ub']:3d} fom={r['fom']:6.2f} |Aut|={r['aut']:.0f} orbitsX={r['orbitsX']} {json.dumps(desc)}")
