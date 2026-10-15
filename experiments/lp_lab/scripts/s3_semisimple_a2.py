"""Algebraic targeting for a=2 LP([[p,1],[1,q]], b) over Z_N with N odd:
k = 2·deg gcd(b, pq+1, x^N+1)  (semisimple ring).  Verified numerically on
random instances at startup.  Enumerate all normalised triples with k ≥ kmin,
then ISD-screen.  Polynomials are Python ints (bit i ↔ x^i)."""
from __future__ import annotations
import argparse, itertools, json, math, random, time
from pathlib import Path
from lp_lab.fastlp import Ab, lp_ab_rows, k_from_rows, connected, logical_reps
from lp_lab import f2
from lp_lab.distance import isd_min_weight


def pdeg(a): return a.bit_length() - 1
def pmod(a, m):
    dm = pdeg(m)
    while a and pdeg(a) >= dm:
        a ^= m << (pdeg(a) - dm)
    return a
def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a
def pmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return r
def poly_of(support): 
    v = 0
    for s in support: v ^= 1 << s
    return v


def k_alg(N, p, q, b):
    F = (1 << N) | 1
    det = pmod(pmul(poly_of(p), poly_of(q)), F) ^ 1
    g = pgcd(pgcd(poly_of(b), det), F)
    return 2 * pdeg(g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True); ap.add_argument("--kmin", type=int, default=12)
    ap.add_argument("--trials", type=int, default=40); ap.add_argument("--dstop", type=int, default=12)
    ap.add_argument("--sample", type=float, default=1.0)
    args = ap.parse_args()
    N = args.N; assert N % 2 == 1
    G = Ab((N,)); n = 4 * N
    rng = random.Random(0)
    # verify formula
    for _ in range(30):
        p = (0, rng.randrange(1, N)); q = tuple(rng.sample(range(N), 2)); b = (0,) + tuple(rng.sample(range(1, N), 2))
        A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
        rX, rZ = lp_ab_rows(G, A, 2, b)
        assert k_from_rows(rX, rZ, n, True) == k_alg(N, p, q, b), (p, q, b)
    print(f"k formula verified (N={N}); n={n}", flush=True)
    perms = G.units(); seen = set()
    from lp_lab.fastlp import canon_poly
    for b1 in range(1, N):
        for b2 in range(b1 + 1, N):
            seen.add(canon_poly(G, (0, b1, b2), perms))
    Bs = sorted(seen)
    out = Path(f"data/sweep_semisimple_a2_N{N}.jsonl"); out.parent.mkdir(exist_ok=True)
    t0 = time.time(); total = kept = 0; results = []
    with out.open("w") as fh:
        for b in Bs:
            for p2 in range(1, N):
                p = (0, p2)
                for q in itertools.combinations(range(N), 2):
                    if args.sample < 1 and rng.random() > args.sample: continue
                    total += 1
                    k = k_alg(N, p, q, b)
                    if k < args.kmin: continue
                    A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
                    rX, rZ = lp_ab_rows(G, A, 2, b)
                    if not connected(rX, rZ, n): continue
                    LZ = logical_reps(rZ, rX, n)
                    r = isd_min_weight(f2.matrix_from_rows(rZ, n), LZ, n, trials=args.trials, seed=1, stop_at=args.dstop)
                    kept += 1
                    rec = {"N": N, "n": n, "k": k, "d_ub": r.weight, "p": list(p), "q": list(q), "b": list(b), "fom": round(k * r.weight ** 2 / n, 2)}
                    results.append(rec); fh.write(json.dumps(rec) + "\n")
            print(f"  b={b}: total={total} kept={kept} t={time.time()-t0:.0f}s", flush=True)
    print(f"done N={N}: total={total} kept={kept} {time.time()-t0:.0f}s")
    best = {}
    for r in results:
        if r["k"] not in best or r["d_ub"] > best[r["k"]]["d_ub"]: best[r["k"]] = r
    for k in sorted(best):
        r = best[k]; print(f"  k={k:3d} d_ub={r['d_ub']:3d} fom={r['fom']:6.2f} p={r['p']} q={r['q']} b={r['b']}")


if __name__ == "__main__":
    main()
