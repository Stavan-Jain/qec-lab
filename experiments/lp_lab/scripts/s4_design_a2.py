"""Design-based search: a=2 codes LP([[p,1],[1,q]], b) with maximal k = 2·dim ker(b).

For each abelian group: table of b-classes (δ_b = dim ker b, D_b = exact min weight of
ker b, fom cap δ_b D_b²/(2N)); for b with cap ≥ --cap enumerate all (p,q) with
pq ≡ 1 mod (b*) (⇒ k = 2δ_b), ISD-screen, classify the best by nauty orbits.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
from lp_lab.fastlp import Ab, lp_ab_rows, k_from_rows, connected, logical_reps, to_code
from lp_lab.design import b_table, max_k_designs
from lp_lab import f2
from lp_lab.distance import isd_min_weight


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", required=True)
    ap.add_argument("--delta-min", type=int, default=3); ap.add_argument("--delta-max", type=int, default=12)
    ap.add_argument("--cap", type=float, default=9.0); ap.add_argument("--trials", type=int, default=60)
    ap.add_argument("--dstop", type=int, default=13); ap.add_argument("--brute-cap", type=int, default=20)
    ap.add_argument("--max-designs", type=int, default=0, help="cap designs per b (0 = all)")
    args = ap.parse_args()
    orders = tuple(int(x) for x in args.orders.split(","))
    G = Ab(orders); N = G.N; n = 4 * N
    tag = "x".join(map(str, orders))
    out = Path(f"data/design_a2_{tag}.jsonl"); out.parent.mkdir(exist_ok=True)
    t0 = time.time()
    tab = b_table(G, delta_range=(args.delta_min, args.delta_max), cap=args.brute_cap)
    print(f"G={tag} n={n}: {len(tab)} b-classes with δ∈[{args.delta_min},{args.delta_max}]", flush=True)
    for bi in tab:
        print(f"   b={G.poly_str(bi.b):28s} δ={bi.delta:2d} D={bi.D} cap={bi.fom_cap:.2f}", flush=True)
    results = []
    with out.open("w") as fh:
        for bi in tab:
            if not (bi.fom_cap == bi.fom_cap) or bi.fom_cap < args.cap:
                continue
            t1 = time.time(); cnt = 0; best = None
            for p, q in max_k_designs(G, bi.b):
                if args.max_designs and cnt >= args.max_designs: break
                A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
                rX, rZ = lp_ab_rows(G, A, 2, bi.b)
                k = k_from_rows(rX, rZ, n, True)
                assert k == 2 * bi.delta, (k, bi)
                if not connected(rX, rZ, n): continue
                cnt += 1
                LZ = logical_reps(rZ, rX, n)
                r = isd_min_weight(f2.matrix_from_rows(rZ, n), LZ, n, trials=args.trials, seed=1, stop_at=args.dstop)
                rec = {"orders": list(orders), "n": n, "k": k, "d_ub": r.weight, "p": list(p), "q": list(q), "b": list(bi.b),
                       "delta": bi.delta, "D_b": bi.D, "fom": round(k * r.weight ** 2 / n, 2)}
                results.append(rec); fh.write(json.dumps(rec) + "\n"); fh.flush()
                if best is None or r.weight > best["d_ub"]: best = rec
            print(f"  b={G.poly_str(bi.b)}: {cnt} designs, best d_ub={best['d_ub'] if best else None} (D_b={bi.D}) t={time.time()-t1:.0f}s", flush=True)
    print(f"done {tag}: {len(results)} designs, {time.time()-t0:.0f}s")
    # summary with orbits for the top codes per (k)
    from lp_lab.structure import tanner_automorphisms
    results.sort(key=lambda r: (-r["fom"], -r["d_ub"]))
    shown = 0; seen = set()
    for r in results:
        key = (r["k"], r["d_ub"], tuple(r["b"]))
        if key in seen: continue
        seen.add(key)
        A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["q"])}
        rX, rZ = lp_ab_rows(G, A, 2, tuple(r["b"]))
        aut = tanner_automorphisms(*[m for m in (to_code(rX, rZ, n).H_X, to_code(rX, rZ, n).H_Z)])
        print(f"  [[{n},{r['k']},≤{r['d_ub']}]] fom≤{r['fom']} b={G.poly_str(tuple(r['b']))} p={G.poly_str(tuple(r['p']))} q={G.poly_str(tuple(r['q']))} |Aut|={aut['order']:.0f} X-orbits={aut['orbits_X']}")
        shown += 1
        if shown >= 12: break


if __name__ == "__main__":
    main()
