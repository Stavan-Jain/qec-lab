"""General exhaustive/sampled sweeps over lifted-product classes.

Modes:
  ab2 : LP(A,b), a=2, A=[[p,1],[1,q]] over any abelian G   (n = 4N)
  AB2 : LP(A,B), both 2×2, A=[[p,1],[1,q]], B=[[r,1],[1,s]]  (n = 8N)
  ab3 : LP(A,b), a=3, A all-monomial with first row/col = 1  (n = 6N)
Pipeline: k → keep k ≥ kmin → connectivity → ISD d_X (and d_Z if A not
symmetric) with early stop below --dstop → JSONL + ranked summary.
"""
from __future__ import annotations
import argparse, itertools, json, math, random, time
from pathlib import Path
from lp_lab.fastlp import Ab, canon_poly, lp_ab_rows, lp_AB_rows, k_from_rows, connected, logical_reps
from lp_lab import f2
from lp_lab.distance import isd_min_weight


def b_classes(G: Ab):
    perms = G.units()
    seen = set()
    N = G.N
    for b1 in range(1, N):
        for b2 in range(b1 + 1, N):
            seen.add(canon_poly(G, (0, b1, b2), perms))
    return sorted(seen)


def screen(rowsX, rowsZ, n, trials, dstop, symmetric):
    LZ = logical_reps(rowsZ, rowsX, n)   # Z-logical reps: ker(H_X)/rs(H_Z)
    HZ = f2.matrix_from_rows(rowsZ, n)
    rX = isd_min_weight(HZ, LZ, n, trials=trials, seed=1, stop_at=dstop)
    dX = rX.weight
    if symmetric:
        return dX, dX
    LX = logical_reps(rowsX, rowsZ, n)
    HX = f2.matrix_from_rows(rowsX, n)
    rZ = isd_min_weight(HX, LX, n, trials=trials, seed=2, stop_at=min(dstop, dX) if dstop else dX)
    return dX, rZ.weight


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["ab2", "AB2", "ab3", "ab3c"], required=True)
    ap.add_argument("--orders", type=str, required=True, help="e.g. 6,6 or 18")
    ap.add_argument("--kmin", type=int, default=8)
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument("--dstop", type=int, default=9, help="stop ISD once a logical of weight ≤ dstop is found")
    ap.add_argument("--sample", type=float, default=1.0, help="fraction of the A-space to sample")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    orders = tuple(int(x) for x in args.orders.split(","))
    G = Ab(orders); N = G.N
    rng = random.Random(args.seed)
    tag = f"{args.mode}_" + "x".join(map(str, orders))
    out = Path(args.out or f"data/sweep_{tag}.jsonl"); out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); total = kept = 0; results = []
    fh = out.open("w")

    def record(k, dX, dZ, desc, n):
        d = min(dX, dZ)
        rec = {"mode": args.mode, "orders": list(orders), "n": n, "k": k, "dX_ub": dX, "dZ_ub": dZ, "d_ub": d,
               "fom": round(k * d * d / n, 3), **desc}
        results.append(rec); fh.write(json.dumps(rec) + "\n"); fh.flush()

    if args.mode == "ab2":
        n = 4 * N
        Bs = b_classes(G)
        print(f"{tag}: n={n}, {len(Bs)} b-classes × {N-1} p × {N*(N-1)//2} q", flush=True)
        for b in Bs:
            for p2 in range(1, N):
                p = (0, p2)
                for q in itertools.combinations(range(N), 2):
                    if args.sample < 1 and rng.random() > args.sample: continue
                    total += 1
                    A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
                    rowsX, rowsZ = lp_ab_rows(G, A, 2, b)
                    k = k_from_rows(rowsX, rowsZ, n, symmetric=True)
                    if k < args.kmin or not connected(rowsX, rowsZ, n): continue
                    dX, dZ = screen(rowsX, rowsZ, n, args.trials, args.dstop, True)
                    kept += 1
                    record(k, dX, dZ, {"p": list(p), "q": list(q), "b": list(b)}, n)
            print(f"  b={b} ({G.poly_str(b)}): total={total} kept={kept} t={time.time()-t0:.0f}s", flush=True)
    elif args.mode == "AB2":
        n = 8 * N
        print(f"{tag}: n={n}, ({N-1} × {N*(N-1)//2})² A,B pairs (sample {args.sample})", flush=True)
        pairs = [((0, p2), q) for p2 in range(1, N) for q in itertools.combinations(range(N), 2)]
        for ia, (p, q) in enumerate(pairs):
            A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
            for ib in range(ia, len(pairs)):        # LP(A,B) ≅ LP(B,A) up to relabeling: enumerate unordered
                r, s = pairs[ib]
                if args.sample < 1 and rng.random() > args.sample: continue
                total += 1
                B = {(0, 0): r, (0, 1): (0,), (1, 0): (0,), (1, 1): s}
                rowsX, rowsZ = lp_AB_rows(G, A, B, 2)
                k = k_from_rows(rowsX, rowsZ, n, symmetric=False)
                if k < args.kmin or not connected(rowsX, rowsZ, n): continue
                dX, dZ = screen(rowsX, rowsZ, n, args.trials, args.dstop, False)
                kept += 1
                record(k, dX, dZ, {"p": list(p), "q": list(q), "r": list(r), "s": list(s)}, n)
            if ia % 50 == 0:
                print(f"  A#{ia}/{len(pairs)}: total={total} kept={kept} t={time.time()-t0:.0f}s", flush=True)
    elif args.mode == "ab3":
        n = 6 * N
        Bs = b_classes(G)
        print(f"{tag}: n={n}, {len(Bs)} b-classes × {N**4} A (sample {args.sample})", flush=True)
        for b in Bs:
            for e in itertools.product(range(N), repeat=4):
                if args.sample < 1 and rng.random() > args.sample: continue
                total += 1
                A = {(0, 0): (0,), (0, 1): (0,), (0, 2): (0,), (1, 0): (0,), (2, 0): (0,),
                     (1, 1): (e[0],), (1, 2): (e[1],), (2, 1): (e[2],), (2, 2): (e[3],)}
                rowsX, rowsZ = lp_ab_rows(G, A, 3, b)
                k = k_from_rows(rowsX, rowsZ, n, symmetric=False)
                if k < args.kmin or not connected(rowsX, rowsZ, n): continue
                dX, dZ = screen(rowsX, rowsZ, n, args.trials, args.dstop, False)
                kept += 1
                record(k, dX, dZ, {"e": list(e), "b": list(b)}, n)
            print(f"  b={b} ({G.poly_str(b)}): total={total} kept={kept} t={time.time()-t0:.0f}s", flush=True)
    elif args.mode == "ab3c":
        n = 6 * N
        Bs = b_classes(G)
        print(f"{tag}: n={n}, {len(Bs)} b-classes × {(N-1)**2 * N*(N-1)//2} A (sample {args.sample})", flush=True)
        for b in Bs:
            for p2 in range(1, N):
                for q2 in range(1, N):
                    for r in itertools.combinations(range(N), 2):
                        if args.sample < 1 and rng.random() > args.sample: continue
                        total += 1
                        A = {(0, 0): (0, p2), (0, 1): (0,), (1, 1): (0, q2), (1, 2): (0,), (2, 2): r, (2, 0): (0,)}
                        rowsX, rowsZ = lp_ab_rows(G, A, 3, b)
                        k = k_from_rows(rowsX, rowsZ, n, symmetric=False)
                        if k < args.kmin or not connected(rowsX, rowsZ, n): continue
                        dX, dZ = screen(rowsX, rowsZ, n, args.trials, args.dstop, False)
                        kept += 1
                        record(k, dX, dZ, {"p": [0, p2], "q": [0, q2], "r": list(r), "b": list(b)}, n)
            print(f"  b={b} ({G.poly_str(b)}): total={total} kept={kept} t={time.time()-t0:.0f}s", flush=True)
    fh.close()
    print(f"done {tag}: {total} codes, {kept} kept, {time.time()-t0:.0f}s")
    best = {}
    for r in results:
        key = r["k"]
        if key not in best or r["d_ub"] > best[key]["d_ub"]:
            best[key] = r
    for k in sorted(best):
        r = best[k]
        print(f"  k={k:3d} d_ub={r['d_ub']:3d} (dX {r['dX_ub']}, dZ {r['dZ_ub']}) fom={r['fom']:6.2f} {json.dumps({kk: v for kk, v in r.items() if kk in ('p','q','r','s','b','e')})}")


if __name__ == "__main__":
    main()
