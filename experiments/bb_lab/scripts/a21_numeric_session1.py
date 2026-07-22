#!/usr/bin/env python
"""A21 session 1: numeric ground truth for f2a6f17e's weight-<=6 cycle census.

Code: bb_neigh_z5z15_f2a6f17e — G = Z5 x Z15, A = 1 + y + x,
B = x*y^6 + x*y^10 + x^2*y^12, [[150,8,8]].

REPO convention throughout (QECLean BBChainComplex.lean):
  cycle:    B*u_L = A*u_R          (∂₁(u_L,u_R) = B⋆u_L + A⋆u_R = 0)
  boundary: (u_L, u_R) = (A⋆f, B⋆f)
  generator column ∂₂δ_g = (A+g, B+g).

Outputs (printed, no files): k, V_A/V_B orbit table, μ(Ann A), μ(Ann B),
exhaustive census of cycles at every split of total weight ≤ 6, boundary
verdicts, family classification.  This is the falsify-first ground truth
for the analytic split map (notes/A21_analytic_base_floor.md).
"""
from __future__ import annotations

import numpy as np
from itertools import combinations

L, M = 5, 15
N = L * M  # 75

A_SUPP = [(0, 0), (0, 1), (1, 0)]        # 1 + y + x
B_SUPP = [(1, 6), (1, 10), (2, 12)]      # xy^6 + xy^10 + x^2 y^12


def idx(gx: int, gy: int) -> int:
    return (gx % L) * M + (gy % M)


def conv_matrix(supp) -> np.ndarray:
    """C[g,h] = P(g-h), so (C f)(g) = sum_h P(g-h) f(h) = (P * f)(g)."""
    C = np.zeros((N, N), dtype=np.uint8)
    for g in range(N):
        gx, gy = divmod(g, M)
        for (px, py) in supp:
            C[g, idx(gx - px, gy - py)] ^= 1
    return C


def rref(Min: np.ndarray):
    """GF(2) RREF.  Returns (R, E, piv) with R = E @ Min (mod 2) in RREF."""
    A = Min.copy().astype(np.uint8)
    rows, cols = A.shape
    E = np.eye(rows, dtype=np.uint8)
    piv: list[int] = []
    r = 0
    for c in range(cols):
        pr = next((rr for rr in range(r, rows) if A[rr, c]), None)
        if pr is None:
            continue
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
            E[[r, pr]] = E[[pr, r]]
        for rr in range(rows):
            if rr != r and A[rr, c]:
                A[rr] ^= A[r]
                E[rr] ^= E[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return A, E, piv


class Solver:
    """Solve M x = c over GF(2); solution set = x0 + ker."""

    def __init__(self, Mmat: np.ndarray):
        self.M = Mmat
        self.R, self.E, self.piv = rref(Mmat)
        self.rank = len(self.piv)
        self.rows, self.cols = Mmat.shape
        # kernel basis
        free = [c for c in range(self.cols) if c not in set(self.piv)]
        ker = []
        for fc in free:
            v = np.zeros(self.cols, dtype=np.uint8)
            v[fc] = 1
            for i, pc in enumerate(self.piv):
                if self.R[i, fc]:
                    v[pc] = 1
            ker.append(v)
        self.ker = np.array(ker, dtype=np.uint8) if ker else \
            np.zeros((0, self.cols), dtype=np.uint8)

    def solve(self, c: np.ndarray):
        y = (self.E @ c.astype(np.int64)) % 2
        if y[self.rank:].any():
            return None
        x = np.zeros(self.cols, dtype=np.uint8)
        x[self.piv] = y[: self.rank].astype(np.uint8)
        return x

    def solve_batch(self, Cmat: np.ndarray):
        """Cmat: (n, rows).  Returns (mask solvable, X particular (n, cols))."""
        Y = (Cmat.astype(np.int64) @ self.E.T) % 2
        mask = ~Y[:, self.rank:].any(axis=1)
        X = np.zeros((Cmat.shape[0], self.cols), dtype=np.uint8)
        X[:, self.piv] = Y[:, : self.rank].astype(np.uint8)
        return mask, X

    def all_kernel_elements(self):
        d = self.ker.shape[0]
        out = np.zeros((1 << d, self.cols), dtype=np.uint8)
        for i in range(1 << d):
            v = np.zeros(self.cols, dtype=np.uint8)
            for b in range(d):
                if i >> b & 1:
                    v ^= self.ker[b]
            out[i] = v
        return out


# ─── GF(16) Fourier ───────────────────────────────────────────────────
EXPT = [1] * 15
for i in range(1, 15):
    v = EXPT[i - 1] << 1
    if v & 16:
        v ^= 0b10011  # t^4 = t + 1
    EXPT[i] = v


def fhat(supp, i, j) -> int:
    """P̂ at character (i,j): sum of ω5^(i·ax) ω15^(j·ay), ω15 = t, ω5 = t³."""
    acc = 0
    for (ax, ay) in supp:
        acc ^= EXPT[(3 * i * ax + j * ay) % 15]
    return acc


def orbits():
    seen, orbs = set(), []
    for i in range(L):
        for j in range(M):
            if (i, j) in seen:
                continue
            o, ci, cj = [], i, j
            while (ci, cj) not in o:
                o.append((ci, cj))
                seen.add((ci, cj))
                ci, cj = (2 * ci) % L, (2 * cj) % M
            orbs.append(o)
    return orbs


def wt(v: np.ndarray) -> int:
    return int(v.sum())


def classify_side(v: np.ndarray, supp, tag: str):
    """Is v a translate of supp?  Returns g or None."""
    s = frozenset(np.flatnonzero(v))
    for g in range(N):
        gx, gy = divmod(g, M)
        t = frozenset(idx(px + gx, py + gy) for (px, py) in supp)
        if t == s:
            return g
    return None


def main() -> None:
    CA, CB = conv_matrix(A_SUPP), conv_matrix(B_SUPP)
    solA, solB = Solver(CA), Solver(CB)          # for A⋆x = c, B⋆x = c
    d2 = np.vstack([CA, CB])                     # ∂₂ : f ↦ (A⋆f, B⋆f)
    sol2 = Solver(d2)
    d1 = np.hstack([CB, CA])                     # ∂₁ : (u_L,u_R) ↦ B⋆u_L+A⋆u_R
    _, _, piv1 = rref(d1)
    rank1, rank2 = len(piv1), sol2.rank
    k = (2 * N - rank1) - rank2
    print(f"ranks: ∂₁ {rank1}, ∂₂ {rank2};  k = {k}  (expect 8)")

    # Fourier table
    print("\n─── vanishing orbits ───")
    VA, VB = set(), set()
    for o in orbits():
        i, j = o[0]
        a, b = fhat(A_SUPP, i, j), fhat(B_SUPP, i, j)
        za, zb = a == 0, b == 0
        if za:
            VA.update(o)
        if zb:
            VB.update(o)
        if za or zb:
            print(f"  orbit rep ({i:2d},{j:2d}) size {len(o)}: "
                  f"Â={'0' if za else hex(a)}, B̂={'0' if zb else hex(b)}"
                  f"{'   *** SHARED ***' if za and zb else ''}")
    print(f"|V_A| = {len(VA)}, |V_B| = {len(VB)}, "
          f"|V_A ∩ V_B| = {len(VA & VB)}  (k = 2·|shared| ⟹ {2*len(VA & VB)})")
    print(f"dim Ann(A) = {solA.ker.shape[0]} (=|V_A|?), "
          f"dim Ann(B) = {solB.ker.shape[0]}")

    # Annihilator weight tables
    annA = solA.all_kernel_elements()
    annB = solB.all_kernel_elements()
    wA = np.sort(annA.sum(axis=1))[1:]  # drop 0
    wB = np.sort(annB.sum(axis=1))[1:]
    print(f"\nμ(Ann A) = {wA.min() if wA.size else '∞'}; "
          f"weights {np.unique(wA, return_counts=True)}")
    print(f"μ(Ann B) = {wB.min() if wB.size else '∞'}; "
          f"weights {np.unique(wB, return_counts=True)}")

    results = []  # (split, u_L, u_R, boundary?)

    def record(split, uL, uR):
        c1 = (CB @ uL.astype(np.int64) + CA @ uR.astype(np.int64)) % 2
        assert not c1.any(), f"not a cycle at split {split}"
        bd = sol2.solve(np.concatenate([uL, uR])) is not None
        results.append((split, uL.copy(), uR.copy(), bd))

    # ── one-sided: (a,0) needs u_L ∈ Ann(B); (0,b) needs u_R ∈ Ann(A)
    for v in annB:
        if 0 < wt(v) <= 6:
            record((wt(v), 0), v, np.zeros(N, dtype=np.uint8))
    for v in annA:
        if 0 < wt(v) <= 6:
            record((0, wt(v)), np.zeros(N, dtype=np.uint8), v)

    # ── (1,b): B+g = A⋆u_R, b ∈ {1,3,5}
    for g in range(N):
        uL = np.zeros(N, dtype=np.uint8)
        uL[g] = 1
        c = (CB @ uL.astype(np.int64)) % 2
        x0 = solA.solve(c)
        if x0 is None:
            continue
        cand = (annA ^ x0)
        for v in cand:
            if wt(v) in (1, 3, 5):
                record((1, wt(v)), uL, v)

    # ── (a,1): B⋆u_L = A+r, a ∈ {3,5}  ((1,1) covered above)
    for r in range(N):
        uR = np.zeros(N, dtype=np.uint8)
        uR[r] = 1
        c = (CA @ uR.astype(np.int64)) % 2
        x0 = solB.solve(c)
        if x0 is None:
            continue
        cand = (annB ^ x0)
        for v in cand:
            if wt(v) in (3, 5):
                record((wt(v), 1), v, uR)

    # ── (2,b): u_L = pair, b ∈ {2,4}
    pairs = list(combinations(range(N), 2))
    Cmat = np.zeros((len(pairs), N), dtype=np.uint8)
    for i, (g1, g2) in enumerate(pairs):
        e = np.zeros(N, dtype=np.uint8)
        e[g1] = e[g2] = 1
        Cmat[i] = (CB @ e.astype(np.int64)) % 2
    mask, X = solA.solve_batch(Cmat)
    for i in np.flatnonzero(mask):
        g1, g2 = pairs[i]
        uL = np.zeros(N, dtype=np.uint8)
        uL[g1] = uL[g2] = 1
        for v in (annA ^ X[i]):
            if wt(v) in (2, 4):
                record((2, wt(v)), uL, v)

    # ── (a,2): u_R = pair, a = 4  ((2,2) covered above)
    Cmat = np.zeros((len(pairs), N), dtype=np.uint8)
    for i, (g1, g2) in enumerate(pairs):
        e = np.zeros(N, dtype=np.uint8)
        e[g1] = e[g2] = 1
        Cmat[i] = (CA @ e.astype(np.int64)) % 2
    mask, X = solB.solve_batch(Cmat)
    for i in np.flatnonzero(mask):
        g1, g2 = pairs[i]
        uR = np.zeros(N, dtype=np.uint8)
        uR[g1] = uR[g2] = 1
        for v in (annB ^ X[i]):
            if wt(v) == 4:
                record((4, wt(v)), v, uR)

    # ── (3,3): all weight-3 u_L
    trips = list(combinations(range(N), 3))
    Cmat = np.zeros((len(trips), N), dtype=np.uint8)
    base = np.zeros(N, dtype=np.uint8)
    for i, t in enumerate(trips):
        base[:] = 0
        base[list(t)] = 1
        Cmat[i] = (CB @ base.astype(np.int64)) % 2
    mask, X = solA.solve_batch(Cmat)
    n33 = 0
    for i in np.flatnonzero(mask):
        uL = np.zeros(N, dtype=np.uint8)
        uL[list(trips[i])] = 1
        for v in (annA ^ X[i]):
            if wt(v) == 3:
                record((3, 3), uL, v)
                n33 += 1
    print(f"\n(3,3) candidates: {mask.sum()} solvable u_L of "
          f"{len(trips)}; {n33} with weight-3 u_R")

    # ── report ───────────────────────────────────────────────────────
    print("\n─── census: all cycles of total weight ≤ 6 ───")
    from collections import Counter
    tally = Counter()
    nonboundary = []
    for split, uL, uR, bd in results:
        tally[(split, bd)] += 1
        if not bd:
            nonboundary.append((split, uL, uR))
    for (split, bd), n in sorted(tally.items()):
        print(f"  split {split}: {n} cycles, "
              f"{'ALL BOUNDARIES' if bd else '*** NON-BOUNDARY ***'}")
    print(f"\nnon-boundary cycles of weight ≤ 6: {len(nonboundary)} "
          f"(target claim: 0)")

    # family classification of the (3,3) census
    print("\n─── (3,3) family classification ───")
    fam = Counter()
    exotic = []
    for split, uL, uR, bd in results:
        if split != (3, 3):
            continue
        gA = classify_side(uL, A_SUPP, "L")
        gB = classify_side(uR, B_SUPP, "R")
        if gA is not None and gB is not None and gA == gB:
            fam["generator column (A+g, B+g)"] += 1
        else:
            fam[f"exotic (uL A-translate: {gA is not None}, "
                f"uR B-translate: {gB is not None})"] += 1
            exotic.append((uL, uR, gA, gB, bd))
    for f, n in fam.items():
        print(f"  {f}: {n}")
    for uL, uR, gA, gB, bd in exotic[:10]:
        print(f"    exotic: uL supp {sorted(np.flatnonzero(uL))} (gA={gA}), "
              f"uR supp {sorted(np.flatnonzero(uR))} (gB={gB}), boundary={bd}")

    # ── the (1,5)/(5,1)/(2,4)/(4,2) coset tables for the analytic attack
    print("\n─── coset structure for the analytic attack ───")
    # w0 = A^{-1}B representative: A⋆w0 = B  (g = 0 case of (1,·))
    cB = np.zeros(N, dtype=np.int64)
    for (px, py) in B_SUPP:
        cB[idx(px, py)] = 1
    x0 = solA.solve(cB % 2)
    if x0 is not None:
        ws = sorted(wt(x0 ^ a) for a in annA)
        print(f"A⋆w = B solvable; coset weight profile {ws}")
    else:
        print("A⋆w = B NOT solvable")
    cA = np.zeros(N, dtype=np.int64)
    for (px, py) in A_SUPP:
        cA[idx(px, py)] = 1
    y0 = solB.solve(cA % 2)
    if y0 is not None:
        ws = sorted(wt(y0 ^ b) for b in annB)
        print(f"B⋆w = A solvable; coset weight profile {ws}")
    else:
        print("B⋆w = A NOT solvable")


if __name__ == "__main__":
    main()
