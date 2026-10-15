"""Design search for LP(A,B), A=[[p,1],[1,q]], B=[[r,1],[1,s]] over cyclic-or-abelian G (n = 8|G|).

k = 2·dim R/(1+p*q*, 1+rs)  (left kernel of [A⊗I | I⊗B†] is z00 ∈ Ann(1+p*q*) ∩ Ann(1+rs)).
Constituent cap: ker(I⊗B) ⊇ {(u0, r u0) : u0 ∈ Ann(1+rs)} ⇒ d ≤ D_B := min |u0| + |r u0|.
Step 1: table over (r,s) of c_B = 1+rs: δ = dim Ann(c_B), D_B (exact if δ ≤ cap).
Step 2: for top c_B: all (p,q) with p*q* ≡ 1 mod (c_B)  [⇒ 1+p*q* ∈ (c_B) ⇒ k = 2δ]; ISD screen.
"""
from __future__ import annotations
import argparse, itertools, json, time
from pathlib import Path
from lp_lab.fastlp import Ab, lp_AB_rows, k_from_rows, connected, logical_reps, canon_poly
from lp_lab.design import poly_int, mul, kernel_code, min_weight_brute, max_k_designs, translate
from lp_lab import f2
from lp_lab.distance import isd_min_weight


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", required=True); ap.add_argument("--delta-min", type=int, default=3); ap.add_argument("--delta-max", type=int, default=10)
    ap.add_argument("--top", type=int, default=6, help="number of c_B classes to design for"); ap.add_argument("--trials", type=int, default=60)
    ap.add_argument("--dstop", type=int, default=13); ap.add_argument("--max-designs", type=int, default=4000)
    args = ap.parse_args()
    orders = tuple(int(x) for x in args.orders.split(",")); G = Ab(orders); N = G.N; n = 8 * N
    tag = "x".join(map(str, orders)); t0 = time.time()
    perms = G.units()
    # Step 1: table of c_B = 1 + r s over normalised (r, s): r = 1 + x^{r2}, s 2-subset; dedupe by canonical c
    seen = {}
    for r2 in range(1, N):
        r = (0, r2)
        for s in itertools.combinations(range(N), 2):
            c = mul(G, poly_int(r), poly_int(s)) ^ 1
            csupp = tuple(i for i in range(N) if (c >> i) & 1)
            if not csupp: continue
            key = canon_poly(G, csupp, perms)
            if key in seen: continue
            K = kernel_code(G, c); delta = len(K)
            if delta < args.delta_min or delta > args.delta_max: seen[key] = None; continue
            # D_B = min over u0 ∈ Ann(c) of |u0| + |r u0|  (exact if delta small)
            if delta <= 16:
                best = None; v = 0
                for i in range(1, 1 << delta):
                    j = (i & -i).bit_length() - 1; v ^= K[j]
                    w = v.bit_count() + mul(G, poly_int(r), v).bit_count()
                    if best is None or w < best: best = w
            else:
                best = None
            seen[key] = {"r": r, "s": s, "c": csupp, "delta": delta, "D_B": best, "cap": (2 * delta * best * best / n) if best else None}
    tab = [v for v in seen.values() if v is not None]
    tab.sort(key=lambda t: -(t["cap"] or -1))
    print(f"G={tag} n={n}: {len(tab)} c_B classes with δ∈[{args.delta_min},{args.delta_max}]", flush=True)
    for t in tab[:15]:
        print(f"   c_B={G.poly_str(t['c']):34s} (r={G.poly_str(t['r'])}, s={G.poly_str(t['s'])}) δ={t['delta']:2d} D_B={t['D_B']} k={2*t['delta']} cap={t['cap'] and round(t['cap'],2)}", flush=True)
    out = Path(f"data/design_AB2_{tag}.jsonl"); fh = out.open("w"); results = []
    for t in tab[: args.top]:
        if not t["cap"] or t["cap"] < 9: continue
        r, s = t["r"], t["s"]; B = {(0, 0): r, (0, 1): (0,), (1, 0): (0,), (1, 1): s}
        c_int = poly_int(t["c"])
        # need 1 + p*q* ∈ (c_B): i.e. (pq)* ≡ 1 mod (c_B) ⇔ pq ≡ 1 mod (c_B*)  → modulus = conj(c_B)
        cconj = 0
        for g in t["c"]: cconj |= 1 << G.neg[g]
        cnt = 0; best = None; t1 = time.time()
        for p, q in max_k_designs(G, t["c"], modulus=cconj):
            if cnt >= args.max_designs: break
            A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
            rX, rZ = lp_AB_rows(G, A, B, 2)
            k = k_from_rows(rX, rZ, n, symmetric=False)
            if k != 2 * t["delta"]:
                print("   k mismatch", k, 2 * t["delta"], p, q, flush=True); continue
            if not connected(rX, rZ, n): continue
            cnt += 1
            LZ = logical_reps(rZ, rX, n)
            res = isd_min_weight(f2.matrix_from_rows(rZ, n), LZ, n, trials=args.trials, seed=1, stop_at=args.dstop)
            dX = res.weight
            dZ = None
            if dX > args.dstop:
                LX = logical_reps(rX, rZ, n)
                resZ = isd_min_weight(f2.matrix_from_rows(rX, n), LX, n, trials=args.trials, seed=2, stop_at=args.dstop)
                dZ = resZ.weight
            d = dX if dZ is None else min(dX, dZ)
            rec = {"mode": "AB2", "orders": list(orders), "n": n, "k": k, "d_ub": d, "dX_ub": dX, "dZ_ub": dZ, "p": list(p), "q": list(q), "r": list(r), "s": list(s), "cB": list(t["c"]), "delta": t["delta"], "D_B": t["D_B"], "fom": round(k * d * d / n, 2)}
            results.append(rec); fh.write(json.dumps(rec) + "\n"); fh.flush()
            if best is None or d > best["d_ub"]: best = rec
        print(f"  c_B={G.poly_str(t['c'])}: {cnt} designs, best d_ub={best['d_ub'] if best else None} (D_B={t['D_B']}) t={time.time()-t1:.0f}s", flush=True)
    fh.close(); print(f"done {tag} {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
