"""A32 Part 0-4: the two-level tower slice calculus for Bravyi [[360,12,<=24]].

Frame:  C (30,6) --y-deck--> BY (30,3) --x-deck--> GB (15,3), lab convention
(H_X = [M_A|M_B] rows = X-stabilizers, H_Z = [M_B^T|M_A^T], X-logicals =
ker H_Z \\ rowspace H_X).  Both deck lifts are twisted (C: y^3 term; BY:
x^25, x^26 terms), which the calculus never notices: all transports are
chain-level ring-hom facts (sigma |-> 1 under the fold, (1+sigma) commutes
with lifted polynomials).

Parts (all falsify-first, hard asserts):
  0  frame construction + structural asserts (chain maps, twist-invariance,
     commuting folds via BX, embedding sections)
  1  slice/carry identities on random chains and cycles, both decks;
     parity lemma verified EXHAUSTIVELY via kernel bases (basis-even =>
     all cycles even)
  2  homology bookkeeping: sig machinery for C/BY/GB, the H1 maps
     p_y*, p_x*, tau_x*, the subspaces R_y / K_x / W, exactness
     (ker tau_x* = im p_x*, im tau_x* = ker p_x*), the reachability
     identity  p_x*^{-1}(W) = R_y, and translation-orbit structure
     (independent re-derivation of A24 SS2.2/SS2.6)
  3  banked (M)@24 stabilizer census (8,461 classes, read-only from the
     main checkout): x-deck decomposition b -> (beta, m2), transport
     assert (beta is a GB stabilizer), compression measurement
     (distinct GB translation-orbits of beta per band)
  4  A24 banked reachable/unreachable band censuses + the ISD weight-24
     witness bank, decomposed in tower coordinates (the falsification
     anchor for the sector trisection)

Outputs: data/a32/tower_validation.json  (summary numbers used in the note)
Read-only inputs: MAIN checkout data/a19 + data/a24 (never written).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.checks import circulant  # noqa: E402

from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a32"
DATA.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260810)


# ----------------------------------------------------------------- GF2 utils
def kernel_basis(M: np.ndarray) -> list[np.ndarray]:
    """Basis of ker M over GF(2), M uint8."""
    M = M.copy() % 2
    rows, n = M.shape
    piv, r = [], 0
    for c in range(n):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        M[[r, r + nz[0]]] = M[[r + nz[0], r]]
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
    Mr = M[:r]
    pivset = set(piv)
    out = []
    for c in [c for c in range(n) if c not in pivset]:
        v = np.zeros(n, dtype=np.uint8)
        v[c] = 1
        for i, pc in enumerate(piv):
            if Mr[i, c]:
                v[pc] ^= 1
        out.append(v)
    return out


def gf2_rank(vectors: list[int]) -> int:
    b, _ = rref_ints(list(vectors))
    return len(b)


def in_span(x: int, basis: list[int], piv: list[int]) -> bool:
    return reduce_int(x, basis, piv) == 0


# ------------------------------------------------------------------- codes
class BBCode:
    """One BB code in the lab convention, with homology helpers."""

    def __init__(self, name: str, lm: tuple[int, int], As: str, Bs: str):
        self.name = name
        self.G = AbelianGroup(lm)
        self.A = Poly.from_string(As, self.G)
        self.B = Poly.from_string(Bs, self.G)
        ng = self.G.cardinality
        self.ng, self.n = ng, 2 * ng
        MA = circulant(self.A).astype(np.uint8) % 2
        MB = circulant(self.B).astype(np.uint8) % 2
        self.HX = np.concatenate([MA, MB], axis=1) % 2          # ng x n
        self.HZ = np.concatenate([MB.T, MA.T], axis=1) % 2      # ng x n
        assert not ((self.HX @ self.HZ.T) % 2).any(), "CSS fails"
        self.rsHX_b, self.rsHX_p = rref_ints([v2i(r) for r in self.HX])
        self.rsHZ_b, self.rsHZ_p = rref_ints([v2i(r) for r in self.HZ])
        self.kerHZ = kernel_basis(self.HZ)   # X-side cycles
        self.kerHX = kernel_basis(self.HX)   # Z-side cycles
        self.k = len(self.kerHZ) - len(self.rsHX_b)
        # Z-logical reps: ker HX mod rowspace HZ
        bb, pp = list(self.rsHZ_b), list(self.rsHZ_p)
        zreps = []
        for kv in self.kerHX:
            x = reduce_int(v2i(kv), bb, pp)
            if x:
                bb.append(x)
                pp.append((x & -x).bit_length() - 1)
                zreps.append(kv)
        assert len(zreps) == self.k
        self.zreps = np.array(zreps, dtype=np.uint8)
        # X-logical reps: ker HZ mod rowspace HX
        bb, pp = list(self.rsHX_b), list(self.rsHX_p)
        xreps = []
        for kv in self.kerHZ:
            x = reduce_int(v2i(kv), bb, pp)
            if x:
                bb.append(x)
                pp.append((x & -x).bit_length() - 1)
                xreps.append(kv)
        assert len(xreps) == self.k
        self.xreps = np.array(xreps, dtype=np.uint8)
        # pairing nondegeneracy
        P = (self.xreps @ self.zreps.T) % 2
        assert gf2_rank([v2i(r) for r in P]) == self.k, "degenerate pairing"

    def sig(self, v: np.ndarray) -> np.ndarray:
        return (self.zreps @ v) % 2

    def is_cycle(self, v: np.ndarray) -> bool:
        return not ((self.HZ @ v) % 2).any()

    def is_stab(self, v: np.ndarray) -> bool:
        return in_span(v2i(v), self.rsHX_b, self.rsHX_p)

    def random_cycle(self) -> np.ndarray:
        c = np.zeros(self.n, dtype=np.uint8)
        for kv in self.kerHZ:
            if RNG.integers(2):
                c ^= kv
        return c


def fold_matrix(Gc: AbelianGroup, Gb: AbelianGroup, fold) -> np.ndarray:
    """Block-diagonal fold on 1-chains: cover coords -> base coords."""
    ngc, ngb = Gc.cardinality, Gb.cardinality
    P1 = np.zeros((ngb, ngc), dtype=np.uint8)
    for i, e in enumerate(Gc):
        P1[Gb.index(fold(e)), i] = 1
    Z = np.zeros_like(P1)
    return np.block([[P1, Z], [Z, P1]]) % 2


def fold_matrix0(Gc: AbelianGroup, Gb: AbelianGroup, fold) -> np.ndarray:
    """Fold on 0-cells (checks indexed by group elements)."""
    P0 = np.zeros((Gb.cardinality, Gc.cardinality), dtype=np.uint8)
    for i, e in enumerate(Gc):
        P0[Gb.index(fold(e)), i] = 1
    return P0


def embed_matrices(Gb: AbelianGroup, Gc: AbelianGroup, emb) -> list[np.ndarray]:
    """Two sheet embeddings on 1-chains: base coords -> cover coords."""
    out = []
    for s in (0, 1):
        E1 = np.zeros((Gc.cardinality, Gb.cardinality), dtype=np.uint8)
        for i, e in enumerate(Gb):
            E1[Gc.index(emb(e, s)), i] = 1
        Z = np.zeros_like(E1)
        out.append(np.block([[E1, Z], [Z, E1]]) % 2)
    return out


class Deck:
    """A free Z2 deck cover-base pair with slice/carry machinery."""

    def __init__(self, cover: BBCode, base: BBCode, fold, emb):
        self.cover, self.base = cover, base
        self.P = fold_matrix(cover.G, base.G, fold)       # 1-chains fold
        self.P0 = fold_matrix0(cover.G, base.G, fold)     # 0-cells fold
        self.EMB = embed_matrices(base.G, cover.G, emb)
        # sections: P o EMB_s = id
        for s in (0, 1):
            assert (self.P @ self.EMB[s] % 2 == np.eye(base.n,
                    dtype=np.uint8)).all(), "fold o emb != id"
        # chain-map transport: HZ_base o P = P0 o HZ_cover
        lhs = (base.HZ @ self.P) % 2
        rhs = (self.P0 @ cover.HZ) % 2
        assert (lhs == rhs).all(), "fold is not a chain map"
        # stabilizer transport: fold of each cover HX row is the base HX row
        for gidx in range(0, cover.ng, max(1, cover.ng // 24)):
            fb = (self.P @ cover.HX[gidx]) % 2
            gbar = base.G.index(fold(cover.G.from_index(gidx)))
            assert (fb == base.HX[gbar]).all(), "stab transport fails"
        # E-system for lifts of a base chain: E v0 = RHS beta
        self.E = (cover.HZ @ ((self.EMB[0] + self.EMB[1]) % 2)) % 2
        self.RHS = (cover.HZ @ self.EMB[1]) % 2
        # twist-invariance: tau(base stabilizer) is a cover cycle
        TAU = (self.EMB[0] + self.EMB[1]) % 2
        chk = (cover.HZ @ ((TAU @ base.HX.T) % 2)) % 2
        assert not chk.any(), "im S_base not <= ker E (twist argument fails)"
        self.TAU = TAU

    def sheets(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(v0, v1) with v = EMB0 v0 + EMB1 v1."""
        v0 = np.zeros(self.base.n, dtype=np.uint8)
        v1 = np.zeros(self.base.n, dtype=np.uint8)
        # invert via section: coordinates of v on each embedded copy
        # EMB[s] has one 1 per column; recover by transpose
        v0 = (self.EMB[0].T @ v) % 2
        v1 = (self.EMB[1].T @ v) % 2
        # exactness of the sheet decomposition
        rec = ((self.EMB[0] @ v0) + (self.EMB[1] @ v1)) % 2
        assert (rec == v).all(), "sheet decomposition failed"
        return v0, v1

    def slice_data(self, v: np.ndarray):
        """(shadow, overflow m, sheet0) with |v| = |shadow| + 2m."""
        v0, v1 = self.sheets(v)
        shadow = (v0 + v1) % 2
        m = int((v0 & v1).sum())
        assert int(v.sum()) == int(shadow.sum()) + 2 * m, "slice identity"
        assert (shadow == (self.P @ v) % 2).all(), "shadow != fold"
        return shadow, m, v0

    def lift(self, v0: np.ndarray, shadow: np.ndarray) -> np.ndarray:
        return ((self.EMB[0] @ v0) + (self.EMB[1] @ ((v0 + shadow) % 2))) % 2


def main():
    t0 = time.monotonic()
    out: dict = {}

    # ------------------------------------------------------------- Part 0
    C = BBCode("C", (30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26")
    BY = BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    GB = BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BX = BBCode("BX", (15, 6), "x^9 + y + y^2", "y^3 + x^10 + x^11")
    assert (C.k, BY.k, GB.k, BX.k) == (12, 8, 8, 8), \
        f"k mismatch {(C.k, BY.k, GB.k, BX.k)}"
    print(f"[{time.monotonic()-t0:5.1f}s] codes built: "
          f"k = {C.k}/{BY.k}/{GB.k}/{BX.k}")

    deck_y = Deck(C, BY, lambda e: (e[0], e[1] % 3),
                  lambda e, s: (e[0], e[1] + 3 * s))
    deck_x = Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                  lambda e, s: (e[0] + 15 * s, e[1]))
    # the other route C -> BX -> GB, for the overflow square
    deck_xc = Deck(C, BX, lambda e: (e[0] % 15, e[1]),
                   lambda e, s: (e[0] + 15 * s, e[1]))
    deck_yb = Deck(BX, GB, lambda e: (e[0], e[1] % 3),
                   lambda e, s: (e[0], e[1] + 3 * s))
    # commuting folds
    lhs = (deck_x.P @ deck_y.P) % 2
    rhs = (deck_yb.P @ deck_xc.P) % 2
    assert (lhs == rhs).all(), "folds do not commute"
    print(f"[{time.monotonic()-t0:5.1f}s] Part 0 OK: chain maps, stab "
          f"transport, twist-invariance, sections, commuting folds")
    out["part0"] = "all structural asserts pass"

    # ------------------------------------------------------------- Part 1
    # slice + carry identities on random cycles; parity exhaustive via bases
    for code in (C, BY, GB, BX):
        odd = [int(kv.sum()) % 2 for kv in code.kerHZ]
        assert not any(odd), f"odd-weight kernel basis vector in {code.name}"
    print(f"[{time.monotonic()-t0:5.1f}s] parity lemma: all kernel bases "
          f"even (C/BY/GB/BX) => every cycle has even weight  [EXHAUSTIVE]")

    n_carry = 0
    for _ in range(200):
        v = C.random_cycle()
        b, m1, v0 = deck_y.slice_data(v)
        assert BY.is_cycle(b), "shadow of cycle not a cycle (y)"
        # carry equation: E v0 = RHS b
        assert (((deck_y.E @ v0) + (deck_y.RHS @ b)) % 2 == 0).all(), \
            "carry equation fails (y)"
        beta, m2, b0 = deck_x.slice_data(b)
        assert GB.is_cycle(beta), "shadow of cycle not a cycle (x)"
        assert (((deck_x.E @ b0) + (deck_x.RHS @ beta)) % 2 == 0).all(), \
            "carry equation fails (x)"
        # two-level identity
        assert int(v.sum()) == int(beta.sum()) + 2 * (m1 + m2), \
            "two-level slice fails"
        # overflow square: other order through BX
        bx, mx, _ = deck_xc.slice_data(v)
        beta2, my2, _ = deck_yb.slice_data(bx)
        assert (beta2 == beta).all(), "composite shadows differ"
        assert m1 + m2 == mx + my2, "overflow square fails"
        n_carry += 1
    # converse direction: solve E v0 = RHS beta for random base cycles,
    # reconstruct, check cycle
    n_conv = 0
    for _ in range(60):
        b = BY.random_cycle()
        rhs = (deck_y.RHS @ b) % 2
        # solve E v0 = rhs by rref of [E | rhs]
        aug = np.concatenate([deck_y.E, rhs[:, None]], axis=1)
        R, piv = _rref_np(aug)
        if any(p == deck_y.E.shape[1] for p in piv):
            continue  # inconsistent: no lift (fine)
        v0 = np.zeros(deck_y.E.shape[1], dtype=np.uint8)
        for i, p in enumerate(piv):
            v0[p] = R[i, -1]
        v = deck_y.lift(v0, b)
        assert C.is_cycle(v), "reconstructed lift is not a cycle"
        n_conv += 1
    print(f"[{time.monotonic()-t0:5.1f}s] Part 1 OK: {n_carry} two-level "
          f"slice/carry/square checks, {n_conv} converse lifts")
    out["part1"] = {"random_cycles": n_carry, "converse_lifts": n_conv}

    # ------------------------------------------------------------- Part 2
    # H1 maps in sig coordinates
    def h1_map(deck: Deck, tau: bool = False) -> np.ndarray:
        """Matrix of p_* (or tau_*) : H1(cover) -> H1(base) (or reverse)."""
        src = deck.cover if not tau else deck.base
        dst = deck.base if not tau else deck.cover
        S = np.array([src.sig(r) for r in src.xreps], dtype=np.uint8)  # k x k
        op = (deck.P if not tau else deck.TAU)
        D = np.array([dst.sig((op @ r) % 2) for r in src.xreps],
                     dtype=np.uint8)
        # M @ S^T = D^T  =>  M = D^T (S^T)^{-1}
        Sinv = _gf2_inv(S.T)
        return (D.T @ Sinv) % 2

    My = h1_map(deck_y)            # H1(C) -> H1(BY), 8x12
    Mx = h1_map(deck_x)            # H1(BY) -> H1(GB), 8x8
    Tx = h1_map(deck_x, tau=True)  # H1(GB) -> H1(BY), 8x8
    rank_My = gf2_rank([v2i(c) for c in My.T])
    rank_Mx = gf2_rank([v2i(c) for c in Mx.T])
    rank_Tx = gf2_rank([v2i(c) for c in Tx.T])
    assert (rank_My, rank_Mx, rank_Tx) == (6, 4, 4), \
        f"H1 map ranks {(rank_My, rank_Mx, rank_Tx)} != (6,4,4)"
    # subspaces of F2^8 (BY-sig space and GB-sig space)
    Ry = _colspace(My)             # im p_y*  (BY sigs), dim 6
    Kx = _kernel_ints(Mx)          # ker p_x* (BY sigs), dim 4
    imTx = _colspace(Tx)           # im tau_x* (BY sigs)
    imMx = _colspace(Mx)           # im p_x*  (GB sigs), dim 4
    kerTx = _kernel_ints(Tx)       # ker tau_x* (GB sigs)
    assert _span_eq(Kx, imTx), "im tau_x* != ker p_x*"
    assert _span_eq(imMx, kerTx), "ker tau_x* != im p_x*"
    assert _span_leq(Kx, Ry), "K_x not <= R_y"
    W = [_apply(Mx, s) for s in Ry]
    Wb, Wp = rref_ints([x for x in W if x])
    assert len(Wb) == 2, f"dim W = {len(Wb)} != 2"
    # preimage of W under Mx == Ry
    pre = _preimage(Mx, Wb, Wp)
    Ryb, Ryp = rref_ints(list(Ry))
    assert _span_eq(pre, Ryb), "p_x*^{-1}(W) != R_y"
    assert len(Ryb) == 6
    print(f"[{time.monotonic()-t0:5.1f}s] Part 2 OK: rank p_y* = 6, "
          f"rank p_x* = 4 = rank tau_x*; im tau = ker p, ker tau = im p; "
          f"K_x <= R_y; dim W = 2; p_x*^{{-1}}(W) = R_y   [A24 SS2.6 "
          f"reproduced independently]")

    # translation action on H1(BY): orbit structure of R_y \ 0 and K_x \ 0
    trans_mats = _translation_action(BY)
    ry_set = _span_points(Ryb, Ryp)
    kx_b, kx_p = rref_ints(list(Kx))
    kx_set = _span_points(kx_b, kx_p)
    orb_all = _orbits(ry_set - {0}, trans_mats)
    orb_kx = _orbits(kx_set - {0}, trans_mats)
    sizes = sorted(len(o) for o in orb_all)
    out["part2"] = {
        "rank_py": 6, "rank_px": 4, "dimW": 2,
        "reachable_orbits": len(orb_all), "orbit_sizes": sizes,
        "Kx_orbits": len(orb_kx),
        "Kx_orbit_sizes": sorted(len(o) for o in orb_kx),
    }
    assert len(orb_all) == 11 and sizes == [3, 3, 3, 3, 3, 3, 9, 9, 9, 9, 9], \
        f"orbit structure {sizes} != A24's 6x3 + 5x9"
    # W \ 0 orbit structure under GB translations
    trans_gb = _translation_action(GB)
    w_set = _span_points(Wb, Wp)
    orb_w = _orbits(w_set - {0}, trans_gb)
    out["part2"]["W_orbits"] = [len(o) for o in orb_w]
    print(f"[{time.monotonic()-t0:5.1f}s]   reachable classes: 63 in 11 "
          f"orbits {sizes}; K_x\\0: 15 classes in "
          f"{len(orb_kx)} orbits {out['part2']['Kx_orbit_sizes']}; "
          f"W\\0 orbits {out['part2']['W_orbits']}")

    # ------------------------------------------------------------- Part 3
    m24 = []
    for line in (MAIN / "a19" / "m24_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            m24.append(r)
    assert len(m24) == 8461
    perms_gb = _translation_perms(GB)
    hist: dict[tuple[int, int], int] = {}
    orb_reps: dict[int, dict] = {}
    per_band_orbits: dict[int, set] = {}
    beta_zero = 0
    for e in m24:
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        assert BY.is_stab(b), "banked b not a BY stabilizer"
        beta, m2, _ = deck_x.slice_data(b)
        wbeta = int(beta.sum())
        assert GB.is_stab(beta), "beta not a GB stabilizer"
        hist[(e["w"], wbeta)] = hist.get((e["w"], wbeta), 0) + 1
        if wbeta == 0:
            beta_zero += 1
        key = _canon(beta, perms_gb)
        per_band_orbits.setdefault(e["band"], set()).add(key)
        if key not in orb_reps:
            orb_reps[key] = {"wbeta": wbeta, "count": 0}
        orb_reps[key]["count"] += 1
    n_orbits = len(orb_reps)
    wbeta_hist: dict[int, int] = {}
    for v in orb_reps.values():
        wbeta_hist[v["wbeta"]] = wbeta_hist.get(v["wbeta"], 0) + 1
    print(f"[{time.monotonic()-t0:5.1f}s] Part 3 OK: 8,461 banked BY "
          f"stabilizers decompose (transport + slice asserts all pass)")
    print(f"    distinct GB shadow orbits: {n_orbits}  (compression "
          f"{8461/n_orbits:.1f}x);  beta=0 records: {beta_zero}")
    print(f"    shadow-weight orbit histogram: "
          f"{dict(sorted(wbeta_hist.items()))}")
    print(f"    per-band distinct shadow orbits: "
          f"{ {b: len(s) for b, s in sorted(per_band_orbits.items())} }")
    out["part3"] = {
        "records": 8461, "distinct_beta_orbits": n_orbits,
        "beta_zero_records": beta_zero,
        "beta_weight_orbit_hist": {str(k): v for k, v in
                                   sorted(wbeta_hist.items())},
        "per_band_orbits": {str(b): len(s)
                            for b, s in sorted(per_band_orbits.items())},
        "joint_hist": {f"{w},{wb}": c
                       for (w, wb), c in sorted(hist.items())},
    }

    # ------------------------------------------------------------- Part 4
    # (a) A24 reachable band-16 census: sector trisection check
    b16 = [json.loads(x) for x in
           (MAIN / "a24" / "cell_census_reach_band16_mchecks.jsonl").open()]
    sectors = []
    for e in b16:
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        assert BY.is_cycle(b) and not BY.is_stab(b), "band16 b not logical"
        sb = v2i(b_sig := BY.sig(b))
        assert in_span(sb, Ryb, Ryp), "band16 class not reachable"
        beta, m2, b0 = deck_x.slice_data(b)
        wbeta = int(beta.sum())
        if wbeta == 0:
            gam = b0
            assert (deck_x.TAU @ gam % 2 == b).all()
            sec = ("B", int(gam.sum()), m2)
            assert GB.is_cycle(gam) and not GB.is_stab(gam)
            # class of gamma outside ker tau_x* = im p_x*
            sg = v2i(GB.sig(gam))
            imb, imp = rref_ints(list(imMx))
            assert not in_span(sg, imb, imp), \
                "diagonal band16 gamma class in im p_x* (would be stab lift)"
        elif GB.is_stab(beta):
            sec = ("C", wbeta, m2)
        else:
            sg = v2i(GB.sig(beta))
            assert in_span(sg, Wb, Wp) and sg != 0, "beta class not in W"
            sec = ("A", wbeta, m2)
        sectors.append({"wb": e["wb"], "m_exact_top": e["m_exact"],
                        "sector": sec[0], "wbeta_or_wgamma": sec[1],
                        "m2": sec[2]})
    cnt = {}
    for s in sectors:
        cnt[s["sector"]] = cnt.get(s["sector"], 0) + 1
    assert cnt == {"A": 3, "B": 2, "C": 1}, f"band16 sectors {cnt}"
    a_rows = [s for s in sectors if s["sector"] == "A"]
    assert all(s["wbeta_or_wgamma"] == 14 and s["m2"] == 1 for s in a_rows)
    b_rows = [s for s in sectors if s["sector"] == "B"]
    assert all(s["wbeta_or_wgamma"] == 8 and s["m2"] == 8 for s in b_rows)
    c_rows = [s for s in sectors if s["sector"] == "C"]
    assert all(s["wbeta_or_wgamma"] == 16 and s["m2"] == 0 for s in c_rows)
    print(f"[{time.monotonic()-t0:5.1f}s] Part 4a OK: band-16 reachable "
          f"census = 3 A (W-lift, |beta|=14, m2=1) + 2 B (diagonal, "
          f"|gamma|=8) + 1 C (flat stab-lift, |beta|=16)  [A24 SS2.5 "
          f"family split REPRODUCED from the calculus]")
    out["part4a"] = sectors

    # (b) unreachable light bands: all shadows must be non-W GB logicals
    unreach = []
    for fn in ("band12_census.jsonl", "cell_census_band14.jsonl"):
        for x in (MAIN / "a24" / fn).open():
            r = json.loads(x)
            if "b_support" in r:
                unreach.append(r)
    n_nonW = 0
    seen = set()
    for e in unreach:
        key = tuple(e["b_support"])
        if key in seen:
            continue
        seen.add(key)
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        assert BY.is_cycle(b) and not BY.is_stab(b)
        sb = v2i(BY.sig(b))
        assert not in_span(sb, Ryb, Ryp), "light band class reachable?!"
        beta = (deck_x.P @ b) % 2
        sg = v2i(GB.sig(beta))
        assert GB.is_cycle(beta) and not GB.is_stab(beta) and beta.sum() > 0
        assert sg != 0 and not in_span(sg, Wb, Wp), \
            "unreachable light b has W or trivial GB shadow"
        n_nonW += 1
    print(f"[{time.monotonic()-t0:5.1f}s] Part 4b OK: {n_nonW} distinct "
          f"light logicals (bands 12-14): every GB shadow is a nontrivial "
          f"non-W logical  [reachability filter power confirmed]")
    out["part4b"] = {"checked": n_nonW}

    # (c) ISD weight-24 witnesses in tower coordinates
    z = np.load(MAIN / "a19" / "isd_class_minima.npz")
    wts, vecs = z["weights"], z["vectors"]
    idx24 = np.nonzero(wts == 24)[0]
    prof: dict[str, int] = {}
    for i in idx24:
        v = vecs[i].astype(np.uint8)
        assert C.is_cycle(v) and not C.is_stab(v), "ISD vector not logical"
        b, m1, _ = deck_y.slice_data(v)
        wb = int(b.sum())
        if wb == 0 or BY.is_stab(b):
            key = f"dangerous(|b|={wb})"
        else:
            beta, m2, _ = deck_x.slice_data(b)
            wbeta = int(beta.sum())
            if wbeta == 0:
                sec = "B"
            elif GB.is_stab(beta):
                sec = "C"
            else:
                sg = v2i(GB.sig(beta))
                assert in_span(sg, Wb, Wp)
                sec = "A"
            key = f"{sec}(|b|={wb},m1={m1},|beta|={wbeta},m2={m2})"
        prof[key] = prof.get(key, 0) + 1
    print(f"[{time.monotonic()-t0:5.1f}s] Part 4c: 100 weight-24 ISD "
          f"witnesses in tower coordinates:")
    for k, v in sorted(prof.items(), key=lambda kv: -kv[1]):
        print(f"      {v:3d}  {k}")
    out["part4c"] = prof

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "tower_validation.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'tower_validation.json'}")


# ------------------------------------------------------- linear-algebra aux
def _rref_np(M: np.ndarray) -> tuple[np.ndarray, list[int]]:
    M = M.copy() % 2
    rows, n = M.shape
    piv, r = [], 0
    for c in range(n):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        M[[r, r + nz[0]]] = M[[r + nz[0], r]]
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
    return M[:r], piv


def _gf2_inv(M: np.ndarray) -> np.ndarray:
    k = M.shape[0]
    aug = np.concatenate([M % 2, np.eye(k, dtype=np.uint8)], axis=1)
    R, piv = _rref_np(aug)
    assert piv == list(range(k)), "matrix not invertible"
    return R[:, k:]


def _colspace(M: np.ndarray) -> list[int]:
    b, _ = rref_ints([v2i(c) for c in (M % 2).T])
    return b


def _kernel_ints(M: np.ndarray) -> list[int]:
    return [v2i(v) for v in kernel_basis(M % 2)]


def _apply(M: np.ndarray, sig_int: int) -> int:
    v = i2v(sig_int, M.shape[1])
    return v2i((M @ v) % 2)


def _span_eq(a: list[int], b: list[int]) -> bool:
    ab, ap = rref_ints(list(a))
    bb, bp = rref_ints(list(b))
    return sorted(ab) == sorted(bb)


def _span_leq(a: list[int], b: list[int]) -> bool:
    bb, bp = rref_ints(list(b))
    return all(in_span(x, bb, bp) for x in a)


def _preimage(M: np.ndarray, Wb: list[int], Wp: list[int]) -> list[int]:
    """Basis of {s : M s in span(Wb)} over GF2 = ker(F o M) where the rows
    of F are the functionals vanishing on span(Wb)."""
    n = M.shape[0]
    Wmat = np.array([i2v(x, n) for x in Wb], dtype=np.uint8) \
        if Wb else np.zeros((0, n), dtype=np.uint8)
    F = np.array(kernel_basis(Wmat), dtype=np.uint8)
    QM = (F @ M) % 2
    return [v2i(v) for v in kernel_basis(QM)]


def _translation_action(code: BBCode) -> list[np.ndarray]:
    """Sig-space matrices of the translation action on H1."""
    S = np.array([code.sig(r) for r in code.xreps], dtype=np.uint8)
    Sinv = _gf2_inv(S.T)
    mats = []
    for t in code.G:
        perm = _perm_for(code, t)
        D = np.array([code.sig(r[perm]) for r in code.xreps],
                     dtype=np.uint8)
        mats.append((D.T @ Sinv) % 2)
    return mats


def _perm_for(code: BBCode, t) -> np.ndarray:
    """index array: (t . v)[i] = v[perm[i]] for the translation by t."""
    ng = code.ng
    perm = np.zeros(2 * ng, dtype=np.int64)
    for i, e in enumerate(code.G):
        j = code.G.index(code.G.sub(e, t))
        perm[i] = j
        perm[ng + i] = ng + j
    return perm


def _translation_perms(code: BBCode) -> list[np.ndarray]:
    return [_perm_for(code, t) for t in code.G]


def _canon(v: np.ndarray, perms: list[np.ndarray]) -> int:
    return min(v2i(v[p]) for p in perms)


def _orbits(points: set[int], mats: list[np.ndarray]) -> list[set[int]]:
    left = set(points)
    orbs = []
    while left:
        x = left.pop()
        orb = {x}
        stack = [x]
        while stack:
            y = stack.pop()
            for M in mats:
                z = _apply(M, y)
                if z not in orb:
                    orb.add(z)
                    stack.append(z)
        left -= orb
        orbs.append(orb)
    return orbs


def _span_points(basis: list[int], piv: list[int]) -> set[int]:
    pts = {0}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return pts


if __name__ == "__main__":
    main()
