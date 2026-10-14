#!/usr/bin/env python3
"""A42 S1b — the omega-syzygy structure over Lambda_a[x^pm] and the
exact per-column cost tables: the two factors of the conjectured
compact-floor law  floor(p) = min over patterns of sum of slot costs.

Objects (plain-pair convention; the X-sector is its antipode image
with identical invariants):

  Lambda_a = F_2[y]/((y^2+y+1)^{2^a})   (the omega chain ring),
  system row (A~, B^) over Lambda_a[x^pm],
  syzygies  Syz = {(f,g): A~ f + B^ g = 0},
  trivial   T   = {(B^ t, A~ t)},
  classes   H_a = Syz / T   (== the cylinder omega-homology at
                             period p = 3*2^a; the odd dilation m
                             changes only the cost tables).

Computed here, exactly:
  1. S-table: for each a <= AMAX and each pure level l (class in
     pi^l * H realizable with all slot contents of valuation >= l),
     the minimum number of (block, column) slots of a representative
     — by exhaustive slot-pattern search with linear-algebra
     feasibility per pattern.
  2. cost tables: for p in PS, the exact minimum weight of an element
     of F_2[y]/(y^p-1) with prescribed omega-chain content class:
     brute force over all 2^p contents for p <= 14 (bucketed by
     omega-component), giving c_nu(p) = min weight with omega
     valuation exactly nu, both with free and with zero complementary
     content.
  3. the pure-omega pattern floor: min over slot patterns and
     class-nontrivial syzygies of  sum_slots c_{nu(slot)}(p) — with
     the complementary factors' content set to ZERO (upper-structure
     floor; the mixed-content refinement is S1c).

Everything is small exact linear algebra over F_2; no SAT.
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

from bb_lab.linalg import nullspace_f2  # noqa: E402
import a42_lib as AL  # noqa: E402

DATA = LAB / "data" / "a42"
DATA.mkdir(parents=True, exist_ok=True)

OMEGA2 = 0b111  # y^2+y+1


class Lam:
    """Lambda_a = F_2[y]/((y^2+y+1)^{2^a}); elements = ints (y-polys)."""

    def __init__(self, a: int):
        self.a = a
        m = OMEGA2
        for _ in range(a):
            m = AL.pmul(m, m)
        self.mod = m
        self.dim = AL.pdeg(m)          # = 2^{a+1}
        self.pi = OMEGA2

    def red(self, z: int) -> int:
        return AL.pmod(z, self.mod)

    def mul(self, x: int, y: int) -> int:
        return AL.pmod(AL.pmul(x, y), self.mod)

    def val(self, z: int) -> int:
        """pi-adic valuation (2^a for z = 0)."""
        if z == 0:
            return 1 << self.a
        v = 0
        while True:
            q, r = AL.pdivmod(z, OMEGA2)
            if r != 0:
                return v
            z = q
            v += 1

    def basis(self):
        return [1 << i for i in range(self.dim)]

    def to_bits(self, z):
        return [(z >> i) & 1 for i in range(self.dim)]

    def from_bits(self, bits):
        return sum((b & 1) << i for i, b in enumerate(bits))


# supports of A~ and B^ as (dx, y-poly-int)
def pair_polys(L: Lam):
    # A~ = y + y^2 + x^3 : columns {0: y+y^2, 3: 1}
    At = {0: L.red(0b110), 3: 1}
    # B^ = 1 + (x + x^2) y^3 : {0: 1, 1: y^3, 2: y^3}
    y3 = L.red(0b1000)
    Bh = {0: 1, 1: y3, 2: y3}
    return At, Bh


class OmegaWindow:
    """Window [0, Wx) of the Lambda_a[x^pm] system."""

    def __init__(self, a: int, Wx: int):
        self.L = Lam(a)
        self.Wx = Wx
        self.At, self.Bh = pair_polys(self.L)
        self.nslots = 2 * Wx                    # (block, col)
        self.nbits = self.nslots * self.L.dim   # F2 variables

    def bit(self, blk, c, i):
        return (blk * self.Wx + c) * self.L.dim + i

    def build_syz_matrix(self) -> np.ndarray:
        """F2 matrix of the map (f,g) -> A~ f + B^ g on the window
        (equations at all columns touched)."""
        L, Wx = self.L, self.Wx
        rows = []
        for c in range(-1, Wx + 4):
            # equation column c: sum over terms
            cols = {}
            for dx, coef in self.At.items():
                cc = c - dx
                if 0 <= cc < Wx:
                    cols.setdefault((0, cc), 0)
                    cols[(0, cc)] ^= coef
            for dx, coef in self.Bh.items():
                cc = c - dx
                if 0 <= cc < Wx:
                    cols.setdefault((1, cc), 0)
                    cols[(1, cc)] ^= coef
            if not cols:
                continue
            # Lambda-valued equation -> L.dim F2 rows
            for out_i in range(L.dim):
                row = np.zeros(self.nbits, dtype=np.uint8)
                any_on = False
                for (blk, cc), coef in cols.items():
                    if coef == 0:
                        continue
                    for in_i in range(L.dim):
                        prod = L.mul(coef, 1 << in_i)
                        if (prod >> out_i) & 1:
                            row[self.bit(blk, cc, in_i)] ^= 1
                            any_on = True
                if any_on:
                    rows.append(row)
        return np.array(rows, dtype=np.uint8)

    def build_trivial(self) -> np.ndarray:
        """Windowed trivial syzygies (B^ t, A~ t), t in a slack window,
        restricted to window-supported images (unit-edge exact)."""
        L, Wx = self.L, self.Wx
        SL = 4
        lo, hi = -SL - 0, Wx + SL
        # image x-range: t col + dx, dx in [0,3]
        ilo, ihi = lo + 0, hi + 3
        WA = ihi - ilo
        nam = 2 * WA * L.dim

        def amb(blk, c, i):
            return (blk * WA + (c - ilo)) * L.dim + i

        rows = []
        for tc in range(lo, hi):
            for ti in range(L.dim):
                row = np.zeros(nam, dtype=np.uint8)
                for dx, coef in self.Bh.items():
                    prod = L.mul(coef, 1 << ti)
                    for oi in range(L.dim):
                        if (prod >> oi) & 1:
                            row[amb(0, tc + dx, oi)] ^= 1
                for dx, coef in self.At.items():
                    prod = L.mul(coef, 1 << ti)
                    for oi in range(L.dim):
                        if (prod >> oi) & 1:
                            row[amb(1, tc + dx, oi)] ^= 1
                rows.append(row)
        D = np.array(rows, dtype=np.uint8)
        win_idx = np.zeros(self.nbits, dtype=np.int64)
        for blk in (0, 1):
            for c in range(Wx):
                for i in range(L.dim):
                    win_idx[self.bit(blk, c, i)] = amb(blk, c, i)
        mask = np.zeros(nam, dtype=bool)
        mask[win_idx] = True
        K = nullspace_f2(D[:, ~mask].T)
        if K.size == 0:
            return np.zeros((0, self.nbits), dtype=np.uint8)
        B = (K @ D[:, win_idx]) % 2
        B = B[B.any(axis=1)]
        return B.astype(np.uint8)

    def functionals(self):
        H = self.build_syz_matrix()
        T = self.build_trivial()
        Z = nullspace_f2(H)
        dimZ = Z.shape[0]
        dimT = int(np.linalg.matrix_rank(T.astype(np.float64))) if False \
            else rank_f2(T)
        N = nullspace_f2(T) if T.size else np.eye(self.nbits,
                                                 dtype=np.uint8)
        P = (N @ Z.T) % 2
        picked, basis = [], []
        for i in range(P.shape[0]):
            v = P[i].copy()
            for (bv, bp) in basis:
                if v[bp]:
                    v ^= bv
            nz = np.flatnonzero(v)
            if nz.size:
                basis.append((v, nz[0]))
                picked.append(i)
        Lrows = N[picked]
        return H, T, Z, Lrows, dimZ, len(basis)


def rank_f2(M: np.ndarray) -> int:
    M = (M % 2).astype(np.uint8).copy()
    if M.size == 0:
        return 0
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def solve_affine_in_pattern(ow: OmegaWindow, H, Lrows, pattern,
                            min_val: int = 0):
    """All syzygies supported in the slot pattern, content valuation
    >= min_val per slot; returns (particular-free basis of the
    solution space restricted to pattern bits, nontrivial-detection
    functionals restricted).  Solution space = nullspace of H
    restricted to allowed bits."""
    L = ow.L
    allowed = np.zeros(ow.nbits, dtype=bool)
    for (blk, c) in pattern:
        # basis of pi^{min_val} Lambda inside Lambda (as bit rows)
        for i in range(L.dim):
            allowed[ow.bit(blk, c, i)] = True
    idx = np.flatnonzero(allowed)
    Hr = H[:, idx]
    if min_val > 0:
        # substitute variables: content = pi^v * u; build basis matrix
        # of pi^v-multiples per slot and compose
        piv_ = 1
        for _ in range(min_val):
            piv_ = L.mul(piv_, L.pi)
        cols = []
        for (blk, c) in pattern:
            for i in range(L.dim):
                z = L.mul(piv_, 1 << i)
                col = np.zeros(ow.nbits, dtype=np.uint8)
                for oi in range(L.dim):
                    if (z >> oi) & 1:
                        col[ow.bit(blk, c, oi)] = 1
                cols.append(col)
        Bmat = np.array(cols, dtype=np.uint8).T  # nbits x nvars
        Hr = (H @ Bmat) % 2
        Zr = nullspace_f2(Hr)
        sols = (Zr @ Bmat.T) % 2   # rows = syzygies in full coords
    else:
        Zr = nullspace_f2(Hr)
        sols = np.zeros((Zr.shape[0], ow.nbits), dtype=np.uint8)
        sols[:, idx] = Zr
    if sols.size == 0:
        return sols
    keep = sols[sols.any(axis=1)]
    return keep


def s_table(a: int, Wx: int, smax: int, log=print):
    """Min slot-count of a class-nontrivial syzygy with all slot
    valuations >= l, for each l = 0..2^a-1."""
    ow = OmegaWindow(a, Wx)
    H, T, Z, Lrows, dimZ, dimT = ow.functionals()
    log(f"  a={a} Wx={Wx}: dim syz={dimZ} dim triv={dimT} "
        f"classes={dimZ-dimT}")
    res = {}
    slots = [(blk, c) for blk in (0, 1) for c in range(Wx)]
    for lev in range(1 << a):
        found = None
        for s in range(1, smax + 1):
            hits = 0
            # patterns up to x-translation: anchor min col = 0
            for pat in itertools.combinations(slots, s):
                if min(c for (_, c) in pat) != 0:
                    continue
                sols = solve_affine_in_pattern(ow, H, Lrows, pat,
                                               min_val=lev)
                if sols.size == 0:
                    continue
                # any nontrivial?
                pair = (Lrows @ sols.T) % 2
                if pair.any():
                    hits += 1
                    found = (s, pat)
                    break
            if found:
                break
        res[lev] = found[0] if found else None
        log(f"    level >= {lev}: min slots = {res[lev]}"
            + (f" (pattern {found[1]})" if found else " (none <= smax)"))
    return res, (dimZ, dimT)


def cost_table(p: int, log=print):
    """Exact per-column cost: for each omega-content valuation class nu
    (and zero content), the min weight of an element of F_2[y]/(y^p-1)
    whose omega-chain component has valuation exactly nu — (i) with
    the complementary content free, (ii) forced zero.
    Brute force 2^p (p <= 14)."""
    assert p <= 14
    assert p % 3 == 0
    a = AL.v2(p)
    L = Lam(a)
    # omega-chain component of z: z mod (y^2+y+1)^{2^a} after reducing
    # mod y^p-1 — the CRT projection is just reduction of the poly.
    ymod = (1 << p) | 1  # y^p - 1 = y^p + 1
    # complementary content zero <=> z is a multiple of
    # (y^p-1)/((y^2+y+1)^{2^a})
    cof = AL.pdivmod(ymod, L.mod)[0]
    free_min = {}
    pure_min = {}
    for z in range(1, 1 << p):
        w = bin(z).count("1")
        lam = AL.pmod(z, L.mod)
        nu = L.val(lam)
        if nu == (1 << a):
            nu = None  # omega content zero
        cur = free_min.get(nu)
        if cur is None or w < cur:
            free_min[nu] = w
        if AL.pmod(z, cof) == 0:
            cur = pure_min.get(nu)
            if cur is None or w < cur:
                pure_min[nu] = w
    log(f"  p={p} (a={a}): free {free_min} | pure {pure_min}")
    return free_min, pure_min


def main():
    t0 = time.time()
    out = {}
    print("== S-tables (omega-syzygy min slots) ==", flush=True)
    stabs = {}
    for a, Wx, smax in ((0, 8, 5), (1, 10, 7), (2, 12, 8)):
        res, dims = s_table(a, Wx, smax,
                            log=lambda s: print(s, flush=True))
        stabs[a] = {"levels": {str(k): v for k, v in res.items()},
                    "dims": dims}
    out["s_tables"] = stabs

    print("\n== cost tables ==", flush=True)
    costs = {}
    for p in (3, 6, 12, 9):
        free_min, pure_min = cost_table(
            p, log=lambda s: print(s, flush=True))
        costs[p] = {"free": {str(k): v for k, v in free_min.items()},
                    "pure": {str(k): v for k, v in pure_min.items()}}
    out["cost_tables"] = costs

    print("\n== pure-omega pattern floors (S x c synthesis) ==",
          flush=True)
    # floor_pred(p) = min over levels l of S_l(a) * pure-cost(l)(p)
    synth = {}
    for p in (3, 6, 12):
        a = AL.v2(p)
        terms = {}
        for lev in range(1 << a):
            S = stabs[a]["levels"].get(str(lev))
            c = costs[p]["pure"].get(str(lev)) if p in costs else None
            c = costs[p]["pure"].get(str(lev))
            if S is not None and c is not None:
                terms[lev] = (S, c, S * c)
        floor = min(v[2] for v in terms.values()) if terms else None
        synth[p] = {"terms": {str(k): v for k, v in terms.items()},
                    "floor_pred_pure": floor}
        print(f"  p={p}: terms level->(S, c, S*c): {terms} => "
              f"pure-pattern floor {floor}", flush=True)
    out["synthesis"] = synth

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s1_syzygy.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s1_syzygy.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
