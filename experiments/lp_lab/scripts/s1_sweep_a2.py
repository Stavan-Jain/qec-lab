"""Exhaustive sweep of the a=2 class LP(A,b) over Z_ℓ, A=[[p,1],[1,q]].

Normalizations (all code equivalences, see notes/LP1_plan.md):
  b ∋ 0, canonical under shifts and units of Z_ℓ;  p = {0, p2};  q any 2-subset.
Pipeline: k (bit-packed rank) → keep k ≥ kmin → ISD upper bound on d_X
(= d_Z by A-symmetry) → connectivity check → JSONL + ranked summary.
"""
from __future__ import annotations
import argparse, itertools, json, math, time, sys
from pathlib import Path
from lp_lab.fastlp2 import hx_rows, hz_rows, k_of
from lp_lab import f2
from lp_lab.distance import isd_min_weight


def units(ell):
    return [u for u in range(1, ell) if math.gcd(u, ell) == 1]


def canon_b(bset, ell, U):
    best = None
    for u in U:
        s = sorted((u * x) % ell for x in bset)
        for sh in s:
            t = tuple(sorted((x - sh) % ell for x in s))
            if best is None or t < best:
                best = t
    return best


def b_classes(ell):
    U = units(ell)
    seen = set()
    for b1 in range(1, ell):
        for b2 in range(b1 + 1, ell):
            c = canon_b((0, b1, b2), ell, U)
            seen.add(c)
    return sorted(seen)


def connected(ell, rowsX, rowsZ, n):
    # union-find over qubits via checks
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for r in rowsX + rowsZ:
        bits = []
        while r:
            lb = r & -r; bits.append(lb.bit_length() - 1); r ^= lb
        for x in bits[1:]:
            union(bits[0], x)
    return len({find(i) for i in range(n)}) == 1


def logical_reps(rowsX, rowsZ, n):
    """Z-logical reps: ker(H_X)/rs(H_Z) as ints (for X-distance nontriviality)."""
    ker = f2.nullspace(rowsX, n)
    acc = dict(f2.echelon(rowsZ))
    reps = []
    for v in ker:
        r = f2.reduce(v, acc)
        if r:
            reps.append(v); acc[r & -r] = r
    return reps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ell", type=int, required=True)
    ap.add_argument("--kmin", type=int, default=6)
    ap.add_argument("--trials", type=int, default=40)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--max-b", type=int, default=None)
    args = ap.parse_args()
    ell = args.ell; n = 4 * ell
    out = Path(args.out or f"data/sweep_a2_ell{ell}.jsonl")
    out.parent.mkdir(parents=True, exist_ok=True)
    Bs = b_classes(ell)
    if args.max_b: Bs = Bs[: args.max_b]
    print(f"ell={ell} n={n}: {len(Bs)} b-classes, {ell-1} p, {ell*(ell-1)//2} q", flush=True)
    t0 = time.time(); total = 0; kept = 0
    results = []
    with out.open("w") as fh:
        for b in Bs:
            for p2 in range(1, ell):
                p = (0, p2)
                for q in itertools.combinations(range(ell), 2):
                    total += 1
                    k = k_of(ell, p, q, b)
                    if k < args.kmin:
                        continue
                    rowsX = hx_rows(ell, p, q, b); rowsZ = hz_rows(ell, p, q, b)
                    if not connected(ell, rowsX, rowsZ, n):
                        continue
                    from lp_lab.f2 import matrix_from_rows
                    LZ = logical_reps(rowsX, rowsZ, n)
                    HZ = matrix_from_rows(rowsZ, n)
                    r = isd_min_weight(HZ, LZ, n, trials=args.trials, seed=1)
                    kept += 1
                    rec = {"ell": ell, "n": n, "k": k, "d_ub": r.weight, "p": list(p), "q": list(q), "b": list(b),
                           "fom": round(k * r.weight ** 2 / n, 3)}
                    results.append(rec); fh.write(json.dumps(rec) + "\n")
            print(f"  b={b}: total={total} kept={kept} elapsed={time.time()-t0:.0f}s", flush=True)
    results.sort(key=lambda r: (-r["d_ub"], -r["k"]))
    print(f"done: {total} codes, {kept} with k>={args.kmin} & connected, {time.time()-t0:.0f}s")
    best_by_k = {}
    for r in results:
        if r["k"] not in best_by_k or r["d_ub"] > best_by_k[r["k"]]["d_ub"]:
            best_by_k[r["k"]] = r
    for k in sorted(best_by_k):
        r = best_by_k[k]; print(f"  k={k:3d}  d_ub={r['d_ub']:3d}  fom={r['fom']:6.2f}  p={r['p']} q={r['q']} b={r['b']}")


if __name__ == "__main__":
    main()
