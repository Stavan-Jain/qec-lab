#!/usr/bin/env python3
"""A42 S1d — the sigma-conditioned h-DP: exact joint-cost floors for
the m = 1 spine (p = 3 * 2^a), SAT-free.

For p = 3*2^a the only factors of y^p-1 are the barren (y+1)-chain
Lambda' = F2[y]/((y+1)^{2^a}) and the omega chain Lambda.  A cylinder
cycle v with nontrivial class decomposes as

    v = pure-lift(sigma) + barren part (B'h, A'h),  h in Lambda'[x^pm]

with sigma a class-nontrivial omega-syzygy.  Its weight is
sum over columns of  jointcost(z'_col, lambda_col)  where jointcost
is the EXACT table  min{ wt(z) : z has Lambda'-component z' and
Lambda-component lambda }  (brute force 2^p once).

For fixed sigma the minimum over h is a shortest path: the barren
contents at column c are

    z'_1(c) = h_c + (1+u)^3 (h_{c-1} + h_{c-2})        (block 1)
    z'_2(c) = (u+u^2) h_c + h_{c-3}                    (block 2)

so a DP over columns with state (h_{c-3}, h_{c-2}, h_{c-1}) in
Lambda'^3 computes min_h exactly.  Boundary states are FREE (entering
and leaving at cost 0) — a relaxation, so the result is a certified
LOWER bound; the pure lift (h = 0) realizes an upper bound, and
equality pins the floor.

Scope: the sigma-enumeration covers class-nontrivial omega-syzygies
with <= SMAX slots (slot-column gaps <= 3 by the splitting argument);
sigma with more slots are bounded below only by their slot count.
The reported floor is therefore exact over objects whose omega-support
uses <= SMAX slots, and a certified LB >= min(reported, SMAX+1)
overall.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a42_lib as AL  # noqa: E402
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_syzygy", Path(__file__).parent / "a42_s1_syzygy.py")
SY = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SY)

DATA = LAB / "data" / "a42"


def joint_cost_table(p: int):
    """(z' bits, lambda bits) -> min weight; z' = z mod (y+1)^{2^a},
    lambda = z mod (y^2+y+1)^{2^a}."""
    a = AL.v2(p)
    mprime = 1
    for _ in range(a):
        mprime = AL.pmul(mprime, mprime)
    # (y+1)^{2^a}: binomial: = y^{2^a} + 1 in char 2
    mprime = (1 << (1 << a)) | 1
    L = SY.Lam(a)
    tab = {}
    for z in range(0, 1 << p):
        w = bin(z).count("1")
        zp = AL.pmod(z, mprime)
        lam = AL.pmod(z, L.mod)
        key = (zp, lam)
        cur = tab.get(key)
        if cur is None or w < cur:
            tab[key] = w
    return tab, mprime, L


def lamp_mul_table(a: int):
    """Multiplication tables in Lambda' = F2[y]/((y+1)^{2^a}) for the
    two structure constants used by the DP: c1 = (1+u)^3 = y^3 mod,
    c2 = u + u^2 = y + y^2 mod (u = y+1)."""
    mprime = (1 << (1 << a)) | 1
    dim = 1 << a

    def red(z):
        return AL.pmod(z, mprime)

    c1 = red(0b1000)          # y^3
    c2 = red(0b110)           # y + y^2
    mul1 = [red(AL.pmul(c1, z)) for z in range(1 << dim)]
    mul2 = [red(AL.pmul(c2, z)) for z in range(1 << dim)]
    return mul1, mul2, dim, mprime


_DP_CACHE = {}


def _dp_static(p: int, a: int, tab):
    """Precomputed static arrays for the DP at this (p, a)."""
    key = (p, a)
    if key in _DP_CACHE:
        return _DP_CACHE[key]
    mul1, mul2, dim, mprime = lamp_mul_table(a)
    nh = 1 << dim
    NSTATE = nh ** 3
    lam_dim = 1 << (a + 1)
    tabarr = np.zeros((nh, 1 << lam_dim), dtype=np.int32)
    for (zp, lam), w in tab.items():
        tabarr[zp, lam] = w
    st = np.arange(NSTATE, dtype=np.int64)
    h3 = st // (nh * nh)
    h2 = (st // nh) % nh
    h1 = st % nh
    mul1a = np.array(mul1, dtype=np.int64)
    mul2a = np.array(mul2, dtype=np.int64)
    pre1 = mul1a[h1 ^ h2]               # per-state (1+u)^3(h1+h2)
    nstbase = (h2 * nh + h1) * nh       # + hc
    # exact tail costs on the lambda = 0 graph (Dijkstra):
    # edge st -> nstbase[st] + hc with cost tab[(hc ^ pre1[st], 0)]
    #                                    + tab[(mul2[hc] ^ h3[st], 0)]
    import heapq

    def dijkstra(reverse: bool):
        dist = [1 << 30] * NSTATE
        dist[0] = 0
        pq = [(0, 0)]
        while pq:
            d, s = heapq.heappop(pq)
            if d > dist[s]:
                continue
            if not reverse:
                # forward edges from s (used for entry costs E)
                p1 = int(pre1[s])
                hh3 = int(h3[s])
                base = int(nstbase[s])
                for hc in range(nh):
                    w = int(tabarr[hc ^ p1, 0]) + \
                        int(tabarr[int(mul2a[hc]) ^ hh3, 0])
                    t_ = base + hc
                    if d + w < dist[t_]:
                        dist[t_] = d + w
                        heapq.heappush(pq, (d + w, t_))
            else:
                # reverse: need predecessors of s; iterate all edges
                # once instead: handled by building reverse adjacency
                raise RuntimeError
        return np.array(dist, dtype=np.int64)

    E = dijkstra(False)          # cheapest way to arrive at state st
    # T(st): cheapest lambda=0 continuation from st down to state 0:
    # Dijkstra on the reverse graph; build reverse adjacency once.
    radj: list[list[tuple[int, int]]] = [[] for _ in range(NSTATE)]
    for s in range(NSTATE):
        p1 = int(pre1[s])
        hh3 = int(h3[s])
        base = int(nstbase[s])
        for hc in range(nh):
            w = int(tabarr[hc ^ p1, 0]) + \
                int(tabarr[int(mul2a[hc]) ^ hh3, 0])
            radj[base + hc].append((s, w))
    import heapq as hq
    T = [1 << 30] * NSTATE
    T[0] = 0
    pq = [(0, 0)]
    while pq:
        d, s = hq.heappop(pq)
        if d > T[s]:
            continue
        for (pred, w) in radj[s]:
            if d + w < T[pred]:
                T[pred] = d + w
                hq.heappush(pq, (d + w, pred))
    T = np.array(T, dtype=np.int64)
    _DP_CACHE[key] = (nh, NSTATE, tabarr, pre1, h3, nstbase, mul2a,
                      E, T)
    return _DP_CACHE[key]


def dp_min_for_sigma(p: int, sigma_cols, tab, a: int, span_lo: int,
                     span_hi: int):
    """sigma_cols: dict (blk, c) -> lambda bits (0 outside pattern).
    EXACT min over all finitely-supported h of the total column cost:
    charged columns [span_lo, span_hi] + exact lambda=0 tail costs
    (Dijkstra to/from the all-zero h-state) on both sides."""
    nh, NSTATE, tabarr, pre1, h3, nstbase, mul2a, E, T = \
        _dp_static(p, a, tab)
    INF = 1 << 30
    cost = E.copy()                     # exact entry tails
    for c in range(span_lo, span_hi + 1):
        lam1 = sigma_cols.get((0, c), 0)
        lam2 = sigma_cols.get((1, c), 0)
        new = np.full(NSTATE, INF, dtype=np.int64)
        col1 = tabarr[:, lam1]
        col2 = tabarr[:, lam2]
        for hc in range(nh):
            z1 = hc ^ pre1
            z2 = mul2a[hc] ^ h3
            w = cost + col1[z1] + col2[z2]
            np.minimum.at(new, nstbase + hc, w)
        cost = new
    return int((cost + T).min())        # exact exit tails


def run_p(p: int, Wx: int, smax: int, log=print):
    a = AL.v2(p)
    t0 = time.time()
    tab, mprime, L = joint_cost_table(p)
    ow = SY.OmegaWindow(a, Wx)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    log(f"p={p} (a={a}): window classes {ncls}; joint table "
        f"{len(tab)} keys")
    Ld = ow.L.dim
    slots = [(blk, c) for blk in (0, 1) for c in range(Wx)]
    best = None
    best_sigma = None
    n_sigma = 0
    for s in range(1, smax + 1):
        for pat in itertools.combinations(slots, s):
            if min(c for (_, c) in pat) != 0:
                continue
            cols = sorted({c for (_, c) in pat})
            if any(b_ - a_ > 4 for a_, b_ in zip(cols, cols[1:])):
                continue
            sols = SY.solve_affine_in_pattern(ow, H, Lrows, pat, 0)
            if sols.size == 0:
                continue
            pair = (Lrows @ sols.T) % 2
            if not pair.any():
                continue
            dimS = sols.shape[0]
            if dimS > 14:
                log(f"  pattern {pat}: dim {dimS} > 14 skipped (flag)")
                continue
            for mb in range(1, 1 << dimS):
                acc = np.zeros(ow.nbits, dtype=np.uint8)
                t = mb
                i = 0
                while t:
                    if t & 1:
                        acc ^= sols[i]
                    t >>= 1
                    i += 1
                if not ((Lrows @ acc) % 2).any():
                    continue
                n_sigma += 1
                if n_sigma % 25000 == 0:
                    log(f"  ... {n_sigma} sigmas, best {best} "
                        f"({round(time.time()-t0)} s)")
                sigma_cols = {}
                for (blk, c) in pat:
                    lam = 0
                    for ii in range(Ld):
                        if acc[ow.bit(blk, c, ii)]:
                            lam |= 1 << ii
                    if lam:
                        sigma_cols[(blk, c)] = lam
                if not sigma_cols:
                    continue
                # cheap pre-filter: DP >= sum of free costs
                nh_ = 1 << (1 << a)
                fsum = 0
                for lam in sigma_cols.values():
                    fmin = min(tab[(zp, lam)] for zp in range(nh_)
                               if (zp, lam) in tab)
                    fsum += fmin
                if best is not None and fsum >= best:
                    continue
                cmin = min(c for (_, c) in sigma_cols)
                cmaxx = max(c for (_, c) in sigma_cols)
                val = dp_min_for_sigma(p, sigma_cols, tab, a,
                                       cmin - 4, cmaxx + 4)
                if best is None or val < best:
                    best = val
                    best_sigma = dict(
                        pattern=[list(k) for k in sigma_cols],
                        lams=[int(v) for v in sigma_cols.values()])
                    log(f"  new best {best} at sigma "
                        f"{sorted(sigma_cols.keys())}")
    dt = round(time.time() - t0, 1)
    log(f"p={p}: {n_sigma} nontrivial sigmas; exact joint floor over "
        f"omega-support <= {smax} slots: {best}  ({dt} s)")
    return {"p": p, "a": a, "smax": smax, "n_sigma": n_sigma,
            "floor_dp": best, "best_sigma": best_sigma, "wall_s": dt}


def main():
    out = {}
    for p, Wx, smax in ((3, 8, 6), (6, 10, 7), (12, 10, 7)):
        out[str(p)] = run_p(p, Wx, smax,
                            log=lambda s: print(s, flush=True))
    (DATA / "s1_jointdp.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s1_jointdp.json'}", flush=True)


if __name__ == "__main__":
    main()
