"""A20 V7: complete analytic enumeration of Y4's light boundaries, |b| <= 19.

Independent, solver-free re-derivation of the boundary census of
Y4 = BB(Z18 x Z4, A = 1+x+x^14 y, B = 1+x y^2+x^2 y^3): the COMPLETE list of
G-translation classes of nonzero b in im d2 with |b| <= 19.  Closes the
band-18 exhaustiveness gap of the SAT census (data/a20/m_census_classes.jsonl:
bands <= 16 UNSAT-closed at 469 classes; band 18 UNTERMINATED at 1,186).

Every ingredient is a finite algebraic fact (A22-V7 grade):

  I1  CRT fibering (a20_fibering, re-derived here from D2 = H_X^T itself):
      z = x^2, base Z2 x Z4 (8 sites), R = F2[z]/(z^9-1) ~= F2 x GF4 x GF64;
      512-triple bijection with exact site-weight table W(eps,d4,d64).
  I2  Component operators extracted from D2 by 8 site-delta applications
      (fiber value 1 has CRT triple (1,1,1), a unit in every component);
      verified against D2 on random f AND on all census rows.
  I3  im d2 = F2^8 x GF4^6 x GF64^8 exactly: eps free (A_eps invertible,
      v_eps = C_eps u_eps), delta4 = the 6-dim GF4 joint image (the 16-elt
      kernel of f4 -> (u4,v4) = ker d2, cf. v7_lever0: ker A* = ker B* =
      ker d2), delta64 free (A64, B64 invertible, v64 = C64 u64).
  I4  Weight parity: W(e,d4,xi) == e (mod 2) for all xi ==> per-site excess
      over W_min(e,d4) is EVEN; |b| always even; B0 always even.
  I5  Outer bound: B0(h,pair) = sum of 16 site W_min's <= |b| for ANY d64
      data ==> cells with B0 > 18 contain no light boundary.
  I6  Inner bound: |argmin(e,d4)| = 3^{W_min(e,d4)}; for any information set
      I of the graph code {(w, C64 w)} (8 independent coordinates of 16),
      any light solution restricted to I has total excess <= E = 18 - B0.
      Min-cost I found by matroid greedy => per-cell enumeration is
      3^{cost(I)} x (even-excess deviation terms), exhaustive by I5+I6.
      Additional complete strategy: two one-sided passes at halved budget
      (min(exc_u, exc_v) <= 2 floor(E/4) since both even, sum <= E).
  I7  Reconstruction: unique fiber polynomial per CRT triple (I1) glues each
      solution to its b vector; canonicalization over the 72 G-translations.

Cross-validation: classes at w <= 16 must equal the SAT census's 469 exactly
(those bands are UNSAT-closed, so ANY discrepancy is a bug); band 18 must
contain the census's 1,186; the surplus is the census's incompleteness.

Usage (from experiments/bb_lab):
    uv run python scripts/a20_v7_completeness.py [--max-cells N]
Writes: data/a20/v7_complete_classes.jsonl   (the definitive class list)
        data/a20/v7_band18_new.jsonl         (classes the SAT census missed)
        data/a20/v7_summary.json
"""
import argparse
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a20_fibering as fib                    # noqa: E402  (read-only reuse)
from bb_lab.checks import bb_check_matrices   # noqa: E402
from bb_lab.group import AbelianGroup         # noqa: E402
from bb_lab.linalg import rank_f2             # noqa: E402
from bb_lab.poly import Poly                  # noqa: E402

OUT = LAB / "data" / "a20"
CENSUS = OUT / "m_census_classes.jsonl"
Y4 = {"frame": (18, 4), "A": "1 + x + x^14*y", "B": "1 + x*y^2 + x^2*y^3"}
LX, LY, N2 = 18, 4, 72
SITES, SIDX = fib.SITES, fib.SIDX             # 8 base sites (s, y)
BUDGET = 18                                   # |b| <= 19 and even => <= 18

# ------------------------------------------------------------ GF(64) tables
def _build_gf64():
    for g in range(2, 64):
        v, seen = 1, set()
        for _ in range(63):
            v = fib.gf64_mul(v, g)
            seen.add(v)
        if len(seen) == 63:
            EXP = np.zeros(126, dtype=np.int64)
            LOG = np.zeros(64, dtype=np.int64)
            v = 1
            for k in range(63):
                EXP[k], LOG[v] = v, k
                v = fib.gf64_mul(v, g)
            EXP[63:] = EXP[:63]
            return EXP, LOG, g
    raise RuntimeError("no primitive element?!")


EXP64, LOG64, GEN64 = _build_gf64()


def g64inv(a):
    assert a
    return int(EXP64[(63 - LOG64[a]) % 63])


def g64_scale(vec, s):
    """s * vec, vec np int array."""
    if s == 0:
        return np.zeros_like(vec)
    out = np.zeros_like(vec)
    nz = vec != 0
    out[nz] = EXP64[(LOG64[vec[nz]] + LOG64[s]) % 63]
    return out


def g64_matvec(M, x):
    out = np.zeros(M.shape[0], dtype=np.int64)
    for j in range(M.shape[1]):
        if x[j]:
            out ^= g64_scale(M[:, j], int(x[j]))
    return out


def g64_mat_inv(M):
    n = M.shape[0]
    aug = np.concatenate([M.astype(np.int64), np.eye(n, dtype=np.int64)], axis=1)
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if aug[i, c]), None)
        assert piv is not None, "matrix singular"
        aug[[r, piv]] = aug[[piv, r]]
        aug[r] = g64_scale(aug[r], g64inv(int(aug[r, c])))
        for i in range(n):
            if i != r and aug[i, c]:
                aug[i] ^= g64_scale(aug[r], int(aug[i, c]))
        r += 1
    return aug[:, n:]


def g64_batch(M, X):
    """M (8x8) applied to batch X (N,8) -> (N,8), all GF64."""
    N = X.shape[0]
    out = np.zeros((N, 8), dtype=np.int16)
    for i in range(8):
        acc = np.zeros(N, dtype=np.int16)
        for j in range(8):
            m = int(M[i, j])
            if m == 0:
                continue
            xj = X[:, j]
            nz = xj != 0
            if nz.any():
                prod = np.zeros(N, dtype=np.int16)
                prod[nz] = EXP64[(LOG64[xj[nz]] + LOG64[m]) % 63]
                acc ^= prod
        out[:, i] = acc
    return out


# -------------------------------------------------------------- GF(4) tools
MUL4 = np.array([[fib.gf4_mul(a, b) for b in range(4)] for a in range(4)],
                dtype=np.int64)
INV4 = [0, 1, 3, 2]  # 1->1, w->w^2, w^2->w  (w*w^2 = 1)


def gf4_rref(rows):
    """RREF of a list of GF4 row vectors; returns (rank, reduced rows)."""
    m = [r.astype(np.int64).copy() for r in rows]
    nrows, ncols = len(m), len(m[0])
    r, out = 0, []
    for c in range(ncols):
        piv = next((i for i in range(r, nrows) if m[i][c]), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        m[r] = MUL4[INV4[m[r][c]]][m[r]]
        for i in range(nrows):
            if i != r and m[i][c]:
                m[i] = m[i] ^ MUL4[m[i][c]][m[r]]
        r += 1
        if r == nrows:
            break
    return r, [row for row in m if row.any()]


# ---------------------------------------------------- component extraction
def extract_components():
    """Build D2 from the census's own H_X and read off the six component
    operators (I2).  Returns dict with everything downstream needs."""
    G = AbelianGroup(Y4["frame"])
    checks = bb_check_matrices(Poly.from_string(Y4["A"], G),
                               Poly.from_string(Y4["B"], G))
    HX = checks.H_X % 2
    D2 = HX.T.copy() % 2                      # b = D2 f, 144 x 72

    AE = np.zeros((8, 8), dtype=np.uint8)
    BE = np.zeros((8, 8), dtype=np.uint8)
    A4 = np.zeros((8, 8), dtype=np.int64)
    B4 = np.zeros((8, 8), dtype=np.int64)
    A64 = np.zeros((8, 8), dtype=np.int64)
    B64 = np.zeros((8, 8), dtype=np.int64)
    for sig, (ss, yy) in enumerate(SITES):
        f = np.zeros(N2, dtype=np.uint8)
        a = 0 if ss == 0 else 9               # fiber exponent j=0, s-part ss
        f[a * LY + yy] = 1
        b = (D2 @ f) % 2
        tu = fib.site_triples(b[:N2])
        tv = fib.site_triples(b[N2:])
        for st in SITES:
            i = SIDX[st]
            AE[i, sig], A4[i, sig], A64[i, sig] = tu[st]
            BE[i, sig], B4[i, sig], B64[i, sig] = tv[st]
    return dict(D2=D2, AE=AE, BE=BE, A4=A4, B4=B4, A64=A64, B64=B64)


def f2_mat_inv(M):
    n = M.shape[0]
    aug = np.concatenate([M.astype(np.uint8) % 2, np.eye(n, dtype=np.uint8)], 1)
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if aug[i, c]), None)
        assert piv is not None, "F2 matrix singular"
        aug[[r, piv]] = aug[[piv, r]]
        for i in range(n):
            if i != r and aug[i, c]:
                aug[i] ^= aug[r]
        r += 1
    return aug[:, n:]


def gf4_matvec(M, x):
    out = np.zeros(M.shape[0], dtype=np.int64)
    for j in range(M.shape[1]):
        if x[j]:
            out ^= MUL4[int(x[j])][M[:, j]]
    return out


# --------------------------------------------------------------- main flow
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-cells", type=int, default=0,
                    help="debug: stop after N survivor cells (0 = full run)")
    args = ap.parse_args()
    t00 = time.time()

    # ---- I1: weight table, parity, argmin structure
    W = fib.build_weight_table()              # (e,d4,d64) -> weight, 512 keys
    WARR = np.zeros((2, 4, 64), dtype=np.int64)
    for (e, d4, d64), w in W.items():
        WARR[e, d4, d64] = w
    for e in range(2):
        for d4 in range(4):
            assert (WARR[e, d4] % 2 == e).all(), "I4 parity violated"
    WMIN = WARR.min(axis=2)                   # (2,4)
    expect_min = np.array([[0, 2, 2, 2], [3, 1, 1, 1]])
    assert (WMIN == expect_min).all(), WMIN
    LEVELS, LEVELCNT = {}, {}
    for e in range(2):
        for d4 in range(4):
            exc = WARR[e, d4] - WMIN[e, d4]
            LEVELS[(e, d4)] = {int(v): np.flatnonzero(exc == v).astype(np.int16)
                               for v in sorted(set(exc.tolist()))}
            cnt = np.zeros(BUDGET + 1, dtype=np.float64)
            for v, arr in LEVELS[(e, d4)].items():
                if v <= BUDGET:
                    cnt[v] = len(arr)
            LEVELCNT[(e, d4)] = cnt
            assert len(LEVELS[(e, d4)][0]) == 3 ** WMIN[e, d4], \
                "argmin size != 3^cost"
    print(f"I1/I4 PASS: 512-triple weight table; parity; argmin sizes 3^cost "
          f"({time.time()-t00:.1f}s)")

    # ---- I2: extract + verify component operators
    C = extract_components()
    D2 = C["D2"]
    assert rank_f2(C["AE"].copy()) == 8 and rank_f2(C["BE"].copy()) == 8
    AEinv = f2_mat_inv(C["AE"])
    CE = (C["BE"] @ AEinv) % 2                # v_eps = CE u_eps
    ce_elt = CE[:, SIDX[(0, 0)]]
    A64inv = g64_mat_inv(C["A64"])
    C64 = np.zeros((8, 8), dtype=np.int64)    # v64 = C64 u64
    for j in range(8):
        C64[:, j] = g64_matvec(C["B64"], A64inv[:, j])
    C64INV = g64_mat_inv(C64)
    c64_elt = C64[:, SIDX[(0, 0)]]
    print(f"I2: C_eps support {int(ce_elt.sum())}/8 (expect 5); "
          f"C64 support {int((c64_elt != 0).sum())}/8 (expect 8)")
    assert int(ce_elt.sum()) == 5 and int((c64_elt != 0).sum()) == 8

    # delta4: joint image basis (I3)
    J = [np.concatenate([C["A4"][:, j], C["B4"][:, j]]) for j in range(8)]
    rk, img_basis = gf4_rref(J)
    assert rk == 6, f"joint delta4 rank {rk} != 6"
    PAIRS = np.zeros((4096, 16), dtype=np.int64)
    for idx in range(4096):
        v = np.zeros(16, dtype=np.int64)
        ii = idx
        for bvec in img_basis:
            d = ii & 3
            ii >>= 2
            if d:
                v ^= MUL4[d][bvec]
        PAIRS[idx] = v
    PBY = {PAIRS[i].tobytes(): i for i in range(4096)}
    assert len(PBY) == 4096, "image parametrization not injective"
    print("I3: delta4 joint image = 6-dim GF4 (4096 pairs), eps/delta64 free")

    # verify components against D2 on random f
    rng = np.random.default_rng(7)
    for _ in range(50):
        f = rng.integers(0, 2, N2).astype(np.uint8)
        tf = fib.site_triples(f)
        fe = np.array([tf[st][0] for st in SITES], dtype=np.uint8)
        f4 = np.array([tf[st][1] for st in SITES], dtype=np.int64)
        f64 = np.array([tf[st][2] for st in SITES], dtype=np.int64)
        b = (D2 @ f) % 2
        tu, tv = fib.site_triples(b[:N2]), fib.site_triples(b[N2:])
        assert all(tu[st] == (int((C["AE"] @ fe)[SIDX[st]] % 2),
                              int(gf4_matvec(C["A4"], f4)[SIDX[st]]),
                              int(g64_matvec(C["A64"], f64)[SIDX[st]]))
                   for st in SITES)
        assert all(tv[st] == (int((C["BE"] @ fe)[SIDX[st]] % 2),
                              int(gf4_matvec(C["B4"], f4)[SIDX[st]]),
                              int(g64_matvec(C["B64"], f64)[SIDX[st]]))
                   for st in SITES)
    print("I2 PASS: component operators == D2 on 50 random f")

    # ---- census decomposition + full-relation verification
    census = [json.loads(l) for l in CENSUS.read_text().splitlines()]
    census = [r for r in census if "w" in r]
    assert len(census) == 1655
    for r in census:
        b = np.zeros(2 * N2, dtype=np.uint8)
        b[r["b_support"]] = 1
        tu, tv = fib.site_triples(b[:N2]), fib.site_triples(b[N2:])
        h = np.array([tu[st][0] for st in SITES], dtype=np.uint8)
        u4 = np.array([tu[st][1] for st in SITES], dtype=np.int64)
        u64 = np.array([tu[st][2] for st in SITES], dtype=np.int64)
        ve = np.array([tv[st][0] for st in SITES], dtype=np.uint8)
        v4 = np.array([tv[st][1] for st in SITES], dtype=np.int64)
        v64 = np.array([tv[st][2] for st in SITES], dtype=np.int64)
        assert (ve == (CE @ h) % 2).all(), "census: v_eps != CE u_eps"
        assert (v64 == g64_matvec(C64, u64)).all(), "census: v64 != C64 u64"
        assert np.concatenate([u4, v4]).tobytes() in PBY, \
            "census: (u4,v4) not in joint image"
        wsum = int(WARR[h, u4, u64].sum() + WARR[ve, v4, v64].sum())
        assert wsum == r["w"], "census: weight formula"
    print(f"I3 PASS: all {len(census)} census rows satisfy the full component "
          f"relations ({time.time()-t00:.1f}s)")

    # ---- outer scan (I5): B0 over 256 x 4096 cells
    t0 = time.time()
    U4 = PAIRS[:, :8]
    V4 = PAIRS[:, 8:]
    survivors = []                            # (h, c, B0)
    b0_hist = Counter()
    for h in range(256):
        hb = np.array([(h >> i) & 1 for i in range(8)], dtype=np.int64)
        vb = (CE @ hb) % 2
        b0 = np.zeros(4096, dtype=np.int64)
        for i in range(8):
            b0 += WMIN[hb[i]][U4[:, i]]
            b0 += WMIN[vb[i]][V4[:, i]]
        for c in np.flatnonzero(b0 <= BUDGET):
            survivors.append((h, int(c), int(b0[c])))
            b0_hist[int(b0[c])] += 1
    # cross-check against the session-1 pre-filter (x16 kernel redundancy)
    expect16 = {0: 16, 6: 384, 8: 384, 10: 2688, 12: 13056, 14: 54528,
                16: 256368, 18: 865024}
    assert {k: 16 * v for k, v in sorted(b0_hist.items())} == expect16, b0_hist
    print(f"I5 PASS: outer scan: {len(survivors)} survivor cells "
          f"(B0 histogram x16 == session-1 measurement) "
          f"({time.time()-t0:.1f}s)")

    # ---- inner enumeration (I6)
    t0 = time.time()
    sols = {}                                 # (h, c, u64 bytes) -> weight
    strat_used = Counter()
    cand_total = 0
    UNIT = np.eye(8, dtype=np.int64)

    def est_count(codes, budget):
        poly = np.zeros(budget + 1, dtype=np.float64)
        poly[0] = 1.0
        for cd in codes:
            poly = np.convolve(poly, LEVELCNT[cd][:budget + 1])[:budget + 1]
        return float(poly.sum())

    def cartesian(arrays):
        n = 1
        for a in arrays:
            n *= len(a)
        out = np.empty((n, len(arrays)), dtype=np.int16)
        rep = n
        for j, a in enumerate(arrays):
            rep //= len(a)
            out[:, j] = np.tile(np.repeat(a, rep), n // (rep * len(a)))
        return out

    def enumerate_basis(coords, codes, budget, solve_mat, wu, wv, h, c):
        """All light solutions whose basis-coordinate excess sums <= budget.
        coords: 8 coordinate indices in 0..15 (u: 0-7, v: 8-15);
        solve_mat: 8x8 GF64, u64 = solve_mat @ (coord values); None = identity."""
        nonlocal cand_total
        levels = []
        for cd in codes:
            levels.append([(exc, vals) for exc, vals in
                           sorted(LEVELS[cd].items()) if exc <= budget])
        found = 0

        def rec(pos, rem, chosen):
            nonlocal found, cand_total
            if pos == 8:
                arrays = [vals for _, vals in chosen]
                cand = cartesian(arrays)
                cand_total += len(cand)
                if solve_mat is None:
                    w = cand
                else:
                    w = g64_batch(solve_mat, cand)
                v64 = g64_batch(C64, w)
                wt = np.zeros(len(cand), dtype=np.int64)
                for i in range(8):
                    wt += wu[i][w[:, i]]
                    wt += wv[i][v64[:, i]]
                keep = np.flatnonzero((wt <= BUDGET) & (wt > 0))
                for k in keep:
                    key = (h, c, w[k].astype(np.int8).tobytes())
                    if key not in sols:
                        sols[key] = int(wt[k])
                        found += 1
                return
            for exc, vals in levels[pos]:
                if exc <= rem:
                    rec(pos + 1, rem - exc, chosen + [(exc, vals)])

        rec(0, budget, [])
        return found

    def greedy_basis(costs, rows):
        """Min-cost information set of the [16,8] graph code (matroid greedy)."""
        order = sorted(range(16), key=lambda i: (costs[i], i))
        ech, chosen = [], []
        for ci in order:
            r = rows[ci].astype(np.int64).copy()
            for pc, er in ech:
                if r[pc]:
                    r ^= g64_scale(er, int(r[pc]))
            nz = np.flatnonzero(r)
            if nz.size:
                pc = int(nz[0])
                ech.append((pc, g64_scale(r, g64inv(int(r[pc])))))
                chosen.append(ci)
                if len(chosen) == 8:
                    return chosen
        raise RuntimeError("graph code rank < 8?!")

    graph_rows = [UNIT[i] for i in range(8)] + [C64[i] for i in range(8)]
    n_cells = 0
    for h, c, b0 in survivors:
        hb = np.array([(h >> i) & 1 for i in range(8)], dtype=np.int64)
        vb = (CE @ hb) % 2
        u4 = U4[c]
        v4 = V4[c]
        cu = [(int(hb[i]), int(u4[i])) for i in range(8)]
        cv = [(int(vb[i]), int(v4[i])) for i in range(8)]
        E = BUDGET - b0
        wu = [WARR[cd] for cd in cu]
        wv = [WARR[cd] for cd in cv]
        costs = [WMIN[cd] for cd in cu] + [WMIN[cd] for cd in cv]
        mixed = greedy_basis(costs, graph_rows)
        codes_m = [cu[i] if i < 8 else cv[i - 8] for i in mixed]
        E4 = 2 * (E // 4)
        ests = {
            "mixed": est_count(codes_m, E),
            "U": est_count(cu, E),
            "V": est_count(cv, E),
            "split": est_count(cu, E4) + est_count(cv, E4),
        }
        best = min(ests, key=ests.get)
        strat_used[best] += 1
        if best == "U":
            enumerate_basis(list(range(8)), cu, E, None, wu, wv, h, c)
        elif best == "V":
            enumerate_basis(list(range(8, 16)), cv, E, C64INV, wu, wv, h, c)
        elif best == "split":
            enumerate_basis(list(range(8)), cu, E4, None, wu, wv, h, c)
            enumerate_basis(list(range(8, 16)), cv, E4, C64INV, wu, wv, h, c)
        else:
            M8 = np.array([graph_rows[i] for i in mixed], dtype=np.int64)
            enumerate_basis(mixed, codes_m, E, g64_mat_inv(M8), wu, wv, h, c)
        n_cells += 1
        if n_cells % 5000 == 0:
            print(f"   [{n_cells}/{len(survivors)}] {len(sols)} light "
                  f"boundary vectors so far; {cand_total:,} candidates "
                  f"({time.time()-t0:.0f}s)", flush=True)
        if args.max_cells and n_cells >= args.max_cells:
            print(f"   DEBUG STOP after {n_cells} cells")
            break
    print(f"I6 PASS: inner enumeration complete: {len(sols)} light boundary "
          f"vectors from {cand_total:,} candidates; strategies {dict(strat_used)} "
          f"({time.time()-t0:.0f}s)")

    # ---- reconstruction (I7) + canonicalization
    t0 = time.time()
    PRE = {}
    for mask in range(512):
        u9 = np.array([(mask >> j) & 1 for j in range(9)], dtype=np.uint8)
        PRE[fib.fiber_triple(u9)] = mask
    AMAP = np.zeros((9, 2), dtype=np.int64)
    for j in range(9):
        a0 = (2 * j) % 9
        for s in range(2):
            AMAP[j, s] = a0 if a0 % 2 == s else a0 + 9
    PERMS = []
    for ta in range(LX):
        for ty in range(LY):
            p = np.zeros(2 * N2, dtype=np.int64)
            for q in range(2 * N2):
                blk, qq = divmod(q, N2)
                a, yy = divmod(qq, LY)
                p[q] = blk * N2 + ((a + ta) % LX) * LY + (yy + ty) % LY
            PERMS.append(p)

    def canon(supp):
        return min(tuple(sorted(p[supp].tolist())) for p in PERMS)

    def rebuild(h, c, u64):
        hb = [(h >> i) & 1 for i in range(8)]
        vb = (CE @ np.array(hb, dtype=np.int64)) % 2
        u4, v4 = U4[c], V4[c]
        v64 = g64_matvec(C64, u64)
        b = np.zeros(2 * N2, dtype=np.uint8)
        for i, (ss, yy) in enumerate(SITES):
            for blk, tri in ((0, (hb[i], int(u4[i]), int(u64[i]))),
                             (1, (int(vb[i]), int(v4[i]), int(v64[i])))):
                mask = PRE[tri]
                for j in range(9):
                    if (mask >> j) & 1:
                        b[blk * N2 + AMAP[j, ss] * LY + yy] = 1
        return b

    classes = {}                              # canon tuple -> weight
    for (h, c, ub), wt in sols.items():
        u64 = np.frombuffer(ub, dtype=np.int8).astype(np.int64)
        b = rebuild(h, c, u64)
        assert int(b.sum()) == wt, "reconstruction weight mismatch"
        key = canon(np.flatnonzero(b))
        if key not in classes:
            classes[key] = wt
    print(f"I7 PASS: {len(sols)} vectors -> {len(classes)} G-translation "
          f"classes ({time.time()-t0:.0f}s)")

    # ---- cross-validation against the SAT census
    census_canon = {}
    for r in census:
        census_canon[canon(np.array(r["b_support"]))] = r["w"]
    assert len(census_canon) == 1655, "census rows not translation-distinct?!"
    mine_by_w = Counter(classes.values())
    cens_by_w = Counter(census_canon.values())
    print(f"\nclass counts by weight:  mine {dict(sorted(mine_by_w.items()))}")
    print(f"                       census {dict(sorted(cens_by_w.items()))}")
    if not args.max_cells:
        low_mine = {k for k, w in classes.items() if w <= 16}
        low_cens = {k for k, w in census_canon.items() if w <= 16}
        assert low_mine == low_cens, (
            f"CLOSED-BAND MISMATCH: mine {len(low_mine)} vs census "
            f"{len(low_cens)}; missing {len(low_cens - low_mine)}, "
            f"extra {len(low_mine - low_cens)}")
        print(f"XVAL PASS: closed bands (w <= 16) match the SAT census "
              f"EXACTLY ({len(low_cens)} classes)")
        b18_mine = {k for k, w in classes.items() if w == 18}
        b18_cens = {k for k, w in census_canon.items() if w == 18}
        assert b18_cens <= b18_mine, \
            f"census band-18 classes missing from mine: {len(b18_cens - b18_mine)}"
        new18 = sorted(b18_mine - b18_cens)
        print(f"BAND 18 DEFINITIVE: {len(b18_mine)} classes "
              f"(SAT census had {len(b18_cens)}; {len(new18)} new)")

        with (OUT / "v7_complete_classes.jsonl").open("w") as fh:
            for key in sorted(classes, key=lambda k: (classes[k], k)):
                fh.write(json.dumps(
                    {"w": classes[key], "b_support": list(key),
                     "in_sat_census": key in census_canon}) + "\n")
        with (OUT / "v7_band18_new.jsonl").open("w") as fh:
            for key in new18:
                fh.write(json.dumps({"w": 18, "b_support": list(key)}) + "\n")
        (OUT / "v7_summary.json").write_text(json.dumps({
            "budget": BUDGET, "survivor_cells": len(survivors),
            "b0_hist": dict(sorted(b0_hist.items())),
            "light_vectors": len(sols), "candidates": cand_total,
            "strategies": dict(strat_used),
            "classes_by_w": {str(k): v for k, v in sorted(mine_by_w.items())},
            "band18_total": len(b18_mine), "band18_new": len(new18),
            "gf64_generator": GEN64,
            "secs": round(time.time() - t00, 1)}, indent=1))
        print(f"\nwrote {OUT / 'v7_complete_classes.jsonl'} "
              f"(+ v7_band18_new.jsonl, v7_summary.json)")
    print(f"TOTAL {time.time()-t00:.0f}s")


if __name__ == "__main__":
    main()
