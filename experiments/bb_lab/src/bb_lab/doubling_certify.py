"""The doubling-aware certification front-end (composed around Tandem).

Given a BB code that is secretly a doubling cover — e.g. the [[360,4,?]]
codes over Z30xZ6 / Z15xZ12 — this module detects the base, certifies the
doubling-template inputs at base scale, and assembles the cover distance
d = 2*d_base with a certificate bundle, in minutes where monolithic
solving stalls for hours.  Pipeline (A28/A29/A30 theory):

  detect      axis literal-lift structure (per-block translation
              normalisation) + (R) k-preservation [A12]
  d_base      exact base distance, both CSS sides: logical-coset
              enumeration (certified two-window BZ) + witnesses
  census      every light stabilizer class up to 2 d_base - 2
              (translation-complete counting certificate)  [A28 lane]
  safe floor  seam-coset minima >= 2 d_base over the kernel orbit
              classes [Prop A14.1 + A29 seam offsets + A30 coset-BZ]
  rung pass   the dangerous sector: (M) at 2 d_base per census class
              (slice identity + sector decomposition)      [A30 §5.5]
  witness     a diagonal lift tau(u) of a weight-d_base logical,
              verified to be a nontrivial cover cycle
  assemble    sector completeness (stabilizer shadows -> rung pass;
              logical shadows lie in im p1 = im delta2 -> safe cosets)
              => d = 2 d_base, X-side; Z by the BB transpose duality.

Claim tier: certified computational data (deterministic enumeration with
composable counting invariants) — NOT kernel-checked Lean.  Every verdict
carries its tier.  Tandem composes in two ways: as the independent
witness lane (pass the certified floor as `-init-lb`, which deletes the
solver's proof phase — the incumbent search only has to FIND the
weight-2d witness), and as the monolithic fallback whenever detection or
any stage fails, refutes, or exceeds scale.

Scope: abelian bivariate codes; literal-lift axis doublings with (R);
|A|, |B| odd; d_base <= DBASE_CAP (the census at W = 2d-2 must fit the
enumeration budget); k <= 8 per side.  Refusals are explicit, never
silent.  Session provenance: A30 (scripts/a30_coset_bz.py,
scripts/a30_rung_pass.py, notes/A30_coset_bz_doubling_certificates.md).
"""

from __future__ import annotations

import itertools
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

import numpy as np

from .checks import bb_check_matrices, circulant
from .cosetbz import (
    NOFF_MAX, build_kernel, coset_base, disjoint_info_sets, in_rowspace,
    kernel_of, pair_radii, rref, run_window, unpack3,
)
from .fibering import kernel_basis, seam_offsets
from .group import AbelianGroup
from .poly import Poly

DBASE_CAP = 12          # d_base beyond this is out of front-end scope
WCAP_DBASE = 15         # logical-coset search cap for finding d_base
NODES_CAP = 2.5e12      # per-sweep enumeration cap (~20-25 min at C-core
                        # rate); larger runs need sharding or the theory
                        # route, and are refused with that message
Progress = Callable[..., None]


def _sweep_nodes(kappa: int, W: int) -> int:
    import math as _m
    r1, r2 = pair_radii(W)
    return (sum(_m.comb(kappa, s) for s in range(1, r1 + 1))
            + sum(_m.comb(kappa, s) for s in range(1, r2 + 1)))


# ------------------------------------------------------------ bitint utils
def v2i(v: np.ndarray) -> int:
    x = 0
    for i in np.nonzero(v)[0]:
        x |= 1 << int(i)
    return x


def i2v(x: int, n: int) -> np.ndarray:
    return np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)


def rref_ints(rows: list[int]) -> tuple[list[int], list[int]]:
    basis: list[int] = []
    piv: list[int] = []
    for r in rows:
        cur = r
        for b, p in zip(basis, piv):
            if (cur >> p) & 1:
                cur ^= b
        if cur:
            p = (cur & -cur).bit_length() - 1
            for i in range(len(basis)):
                if (basis[i] >> p) & 1:
                    basis[i] ^= cur
            basis.append(cur)
            piv.append(p)
    return basis, piv


def reduce_int(x: int, basis: list[int], piv: list[int]) -> int:
    for b, p in zip(basis, piv):
        if (x >> p) & 1:
            x ^= b
    return x


# ---------------------------------------------------------------- detection
@dataclass
class DoublingCandidate:
    axis: int
    shift_A: int                      # per-block axis translations applied
    shift_B: int
    base_group: tuple[int, ...]
    base_A: str
    base_B: str
    k_cover: int
    k_base: int
    R_holds: bool
    note: str = ""


def _axis_window_shift(exps: list[int], order2: int) -> Optional[int]:
    """A translation a with (e - a) mod order2 < order2/2 for all e."""
    ell = order2 // 2
    for a in range(order2):
        if all((e - a) % order2 < ell for e in exps):
            return a
    return None


def detect(G: AbelianGroup, A: Poly, B: Poly) -> list[DoublingCandidate]:
    """Axis literal-lift doubling structures of (G, A, B), with (R).

    For each even-order axis, look for per-block axis translations after
    which every exponent lies in the lower half-window; the base is the
    halved group with the same (shifted) supports.  Translating A and B
    independently is a coordinate permutation of the code (u-block by
    T_a, v-block by T_b), so weights and homology transport exactly.
    """
    out: list[DoublingCandidate] = []
    k_cover = 2 * kernel_basis(A, B).shape[0]
    for axis, order in enumerate(G.orders):
        if order % 2 or order < 4:
            continue
        ell = order // 2
        sA = _axis_window_shift([g[axis] for g in A.support], order)
        sB = _axis_window_shift([g[axis] for g in B.support], order)
        if sA is None or sB is None:
            continue
        base_orders = list(G.orders)
        base_orders[axis] = ell
        Gb = AbelianGroup(tuple(base_orders))

        def shifted(P: Poly, a: int) -> list[tuple[int, ...]]:
            sup = []
            for g in P.support:
                e = list(g)
                e[axis] = (e[axis] - a) % order
                assert e[axis] < ell
                sup.append(tuple(e))
            return sup

        Ab = Poly.from_support(shifted(A, sA), Gb)
        Bb = Poly.from_support(shifted(B, sB), Gb)
        if Ab.weight() != A.weight() or Bb.weight() != B.weight():
            continue  # monomials collide mod ell: not a literal lift
        k_base = 2 * kernel_basis(Ab, Bb).shape[0]
        out.append(DoublingCandidate(
            axis=axis, shift_A=sA, shift_B=sB,
            base_group=tuple(base_orders),
            base_A=Ab.canonical_string(), base_B=Bb.canonical_string(),
            k_cover=k_cover, k_base=k_base,
            R_holds=(k_cover == k_base),
            note="" if k_cover == k_base else
                 f"(R) fails: k_cover {k_cover} != k_base {k_base}",
        ))
    return out


# --------------------------------------------------------------- base tools
def code_rows(A: Poly, B: Poly) -> np.ndarray:
    MA = circulant(A).astype(np.uint8)
    MB = circulant(B).astype(np.uint8)
    return np.concatenate([MA.T, MB.T], axis=1) % 2


@dataclass
class BaseTools:
    G: AbelianGroup
    A: Poly
    B: Poly
    binp: Path = field(init=False)
    n: int = field(init=False)

    def __post_init__(self) -> None:
        self.binp = build_kernel()
        ng = self.G.cardinality
        self.n = 2 * ng
        self.MA = circulant(self.A).astype(np.uint8) % 2
        self.MB = circulant(self.B).astype(np.uint8) % 2
        self.MS = code_rows(self.A, self.B)
        (self.I1, self.G1, self.I2, self.G2, self.kappa) = \
            disjoint_info_sets(self.MS)
        self.S_rref, self.S_piv = rref(self.MS)

    def logical_reps(self, side: str) -> tuple[np.ndarray, np.ndarray, list]:
        """(H, Mrows, reps) for one CSS side; reps independent mod rowspace."""
        if side == "X":  # ker D1 / im D2
            H = np.concatenate([self.MB, self.MA], axis=1) % 2
            Mrows = self.MS
        else:            # ker D2^T / im D1^T
            H = np.concatenate([self.MA.T, self.MB.T], axis=1) % 2
            Mrows = np.concatenate([self.MB, self.MA], axis=1) % 2
        Rr, piv = rref(Mrows)
        reps: list[np.ndarray] = []
        Rcomb = [r.copy() for r in Rr]
        pcomb = list(piv)
        for v in kernel_of(H):
            w = v.copy()
            for i, pc in enumerate(pcomb):
                if w[pc]:
                    w ^= Rcomb[i]
            if w.any():
                reps.append(v)
                Rcomb.append(w)
                pcomb.append(int(np.nonzero(w)[0][0]))
        return H, Mrows, reps


# ---------------------------------------------------- stage: exact d_base
def d_side_exact(
    bt: BaseTools, side: str, Wcap: int, deadline: float, threads: int,
    workdir: Path,
) -> dict:
    """Exact side distance if <= Wcap (certified), else floor > Wcap."""
    H, Mrows, reps = bt.logical_reps(side)
    k = len(reps)
    if k > 14:
        return {"side": side, "error": f"k = {k}: 2^k logical cosets is "
                "beyond the front-end (needs per-code class structure)"}
    combos = []
    for t in range(1, 1 << k):
        L = np.zeros(bt.n, dtype=np.uint8)
        for j in range(k):
            if (t >> j) & 1:
                L ^= reps[j]
        combos.append(L)
    Rr, piv = rref(Mrows)
    I1, G1, I2, G2, kap = disjoint_info_sets(Mrows)
    nodes = []
    # progressive ladder: small radii find small d instantly; each rung's
    # complete pair certifies exactness of anything it finds
    ladder = [w for w in (9, 11, 13, Wcap) if w <= Wcap]
    for Wtry in ladder:
        best: Optional[np.ndarray] = None
        r1, r2 = pair_radii(Wtry)
        for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
            bases = [coset_base(Gs, window, L) for L in combos]
            for b in bases:
                if int(b.sum()) <= Wtry and (best is None
                                             or b.sum() < best.sum()):
                    best = b.copy()
            # the C kernel takes <= NOFF_MAX base words per invocation;
            # chunking repeats the (cheap, low-radius) walk per chunk and
            # each chunk is its own complete certified pass
            for c0 in range(0, len(bases), NOFF_MAX):
                res = run_window(bt.binp,
                                 f"dbase_{side}_w{wi}_c{c0 // NOFF_MAX}",
                                 Gs, bases[c0:c0 + NOFF_MAX], r, Wtry,
                                 deadline, threads=threads, workdir=workdir)
                nodes.append(res["nodes"])
                for j, hx in res.pop("hit_rows"):
                    c = unpack3(hx, bt.n)
                    assert not ((H @ c) % 2).any()
                    assert not in_rowspace(Rr, piv, c)
                    if best is None or c.sum() < best.sum():
                        best = c
        if best is not None:
            return {"side": side, "k": k, "d": int(best.sum()),
                    "witness": best, "nodes": nodes}
    return {"side": side, "k": k, "floor_gt": Wcap, "nodes": nodes}


# --------------------------------------------------------- stage: census
def census(
    bt: BaseTools, W: int, deadline: float, threads: int, workdir: Path,
    progress: Optional[Progress] = None,
) -> dict:
    """All G-translation classes of nonzero stabilizers of weight <= W."""
    zero = np.zeros(bt.n, dtype=np.uint8)
    raw: set[int] = set()
    r1, r2 = pair_radii(W)
    stats = []
    for wi, (window, Gs, r) in enumerate(
            [(bt.I1, bt.G1, r1), (bt.I2, bt.G2, r2)]):
        bases = [coset_base(Gs, window, zero)]
        assert not bases[0].any()
        res = run_window(bt.binp, f"census_w{wi}", Gs, bases, r, W, deadline,
                         threads=threads, workdir=workdir)
        for _, hx in res.pop("hit_rows"):
            raw.add(int(hx, 16))
        stats.append({k: res[k] for k in ("nodes", "hits", "wall_s")})
        if progress:
            progress("census-window", window=wi, **stats[-1])
    # translation canonicalization, vectorised
    ng = bt.G.cardinality
    perms = []
    for gel in bt.G:
        p = np.empty(bt.n, dtype=np.int64)
        for i, e in enumerate(bt.G):
            te = bt.G.reduce(tuple(a + b for a, b in zip(e, gel)))
            j = bt.G.index(te)
            p[j] = i          # translated vector reads coordinate i into j
            p[ng + j] = ng + i
        perms.append(p)
    vecs = np.array([i2v(x, bt.n) for x in sorted(raw)], dtype=np.uint8)
    if len(vecs) == 0:
        return {"W": W, "n_classes": 0, "classes": [], "windows": stats,
                "raw_hits": 0}

    def lexkey(V: np.ndarray) -> np.ndarray:
        packed = np.packbits(V, axis=1)
        pad = (-packed.shape[1]) % 8
        if pad:
            packed = np.pad(packed, ((0, 0), (0, pad)))
        packed = np.ascontiguousarray(packed)
        return packed.view(">u8").reshape(len(V), -1)

    minkey = None
    for p in perms:
        key = lexkey(vecs[:, p])
        if minkey is None:
            minkey = key.copy()
        else:  # lexicographic column-by-column running minimum
            lt = np.zeros(len(vecs), dtype=bool)
            eq = np.ones(len(vecs), dtype=bool)
            for c in range(key.shape[1]):
                lt |= eq & (key[:, c] < minkey[:, c])
                eq &= key[:, c] == minkey[:, c]
            minkey[lt] = key[lt]
    view = minkey.view([("", minkey.dtype)] * minkey.shape[1]).ravel()
    _, first = np.unique(view, return_index=True)
    classes = []
    for idx in sorted(first):
        v = vecs[idx]
        classes.append({"weight": int(v.sum()), "vec": v})
    classes.sort(key=lambda c: c["weight"])
    hist: dict[int, int] = {}
    for c in classes:
        hist[c["weight"]] = hist.get(c["weight"], 0) + 1
    return {"W": W, "n_classes": len(classes), "classes": classes,
            "weight_histogram": hist, "windows": stats,
            "raw_hits": len(raw)}


# ------------------------------------------------------ stage: safe floor
def safe_floor(
    bt: BaseTools, axis: int, target: int, deadline: float, threads: int,
    workdir: Path,
) -> dict:
    offs = seam_offsets(bt.A, bt.B, axis)
    labeled = [(f"rep{i}", np.concatenate([su, sv]))
               for i, (z, su, sv) in enumerate(offs)]
    t0s = [t for _, t in labeled]
    pars = set()
    for t in t0s:
        assert not in_rowspace(bt.S_rref, bt.S_piv, t), \
            "seam offset is a boundary — class 0?"
        pars.add(int(t.sum()) % 2)
    W = target - 1
    if len(pars) == 1 and (W - pars.pop()) % 2 == 1:
        W -= 1
    r1, r2 = pair_radii(W)
    best: dict[int, tuple[int, np.ndarray]] = {}
    stats = []
    for wi, (window, Gs, r) in enumerate(
            [(bt.I1, bt.G1, r1), (bt.I2, bt.G2, r2)]):
        bases = [coset_base(Gs, window, t) for t in t0s]
        for j, b in enumerate(bases):
            w = int(b.sum())
            if w <= W and (j not in best or w < best[j][0]):
                best[j] = (w, b.copy())
        res = run_window(bt.binp, f"safe_w{wi}", Gs, bases, r, W, deadline,
                         threads=threads, workdir=workdir)
        for j, hx in res.pop("hit_rows"):
            c = unpack3(hx, bt.n)
            w = int(c.sum())
            assert in_rowspace(bt.S_rref, bt.S_piv, (c + t0s[j]) % 2)
            if j not in best or w < best[j][0]:
                best[j] = (w, c)
        stats.append({k: res[k] for k in ("nodes", "hits", "wall_s")})
    per_class = {}
    for j, (lab, _) in enumerate(labeled):
        if j in best:
            w, c = best[j]
            per_class[lab] = {"verdict": "REFUTED", "min_weight_found": w}
        else:
            per_class[lab] = {"verdict": "CERTIFIED", "floor": target}
    return {"target": target, "W": W, "r_pair": [r1, r2],
            "n_orbit_reps": len(labeled), "per_class": per_class,
            "windows": stats,
            "certified": all(p["verdict"] == "CERTIFIED"
                             for p in per_class.values()),
            "min_refuted": min((p["min_weight_found"]
                                for p in per_class.values()
                                if p["verdict"] == "REFUTED"), default=None)}


# ------------------------------------------------------- stage: rung pass
class RungEngine:
    """The A30 §5.5 dangerous-sector checker, per (base, axis)."""

    def __init__(self, bt: BaseTools, axis: int, d: int):
        self.bt, self.axis, self.d = bt, axis, d
        G = bt.G
        gc = list(G.orders)
        gc[axis] *= 2
        self.Gc = AbelianGroup(tuple(gc))
        self.Ac = Poly.from_support(bt.A.support, self.Gc)
        self.Bc = Poly.from_support(bt.B.support, self.Gc)
        ng, nc = G.cardinality, self.Gc.cardinality
        self.n = 2 * ng
        MAc = circulant(self.Ac).astype(np.uint8) % 2
        MBc = circulant(self.Bc).astype(np.uint8) % 2
        self.D1c = np.concatenate([MBc, MAc], axis=1) % 2
        ell = G.orders[axis]
        emb = [np.zeros((nc, ng), dtype=np.uint8) for _ in range(2)]
        for i, e in enumerate(G):
            for s in (0, 1):
                ec = list(e)
                ec[axis] += s * ell
                emb[s][self.Gc.index(tuple(ec)), i] = 1
        z0 = np.zeros_like(emb[0])
        self.EMB = [np.block([[emb[s], z0], [z0, emb[s]]]) for s in (0, 1)]
        self.E = (self.D1c @ ((self.EMB[0] + self.EMB[1]) % 2)) % 2
        self.RHS_OP = (self.D1c @ self.EMB[1]) % 2
        self.MSc = code_rows(self.Ac, self.Bc)
        self.S_basis, self.S_piv = rref_ints([v2i(r) for r in bt.MS])
        self.Sc_basis, self.Sc_piv = rref_ints([v2i(r) for r in self.MSc])
        self.E_rows = [v2i(self.E[i]) for i in range(self.E.shape[0])]
        self.E_cols = [v2i(self.E[:, j]) for j in range(self.n)]
        kerE = kernel_of(self.E)
        for idx in (0, 1, min(7, bt.MS.shape[0] - 1)):
            assert not ((self.E @ bt.MS[idx]) % 2).any(), "im S !<= ker E"
        sec = []
        aug_b = list(self.S_basis)
        aug_p = list(self.S_piv)
        for kv in kerE:
            x = reduce_int(v2i(kv), aug_b, aug_p)
            if x:
                aug_b.append(x)
                aug_p.append((x & -x).bit_length() - 1)
                sec.append(v2i(kv))
        self.sector_basis = sec

    def solve_E(self, rhs: np.ndarray):
        n = self.n
        aug = [self.E_rows[i] | (int(rhs[i]) << n)
               for i in range(len(self.E_rows))]
        basis: list[int] = []
        piv: list[int] = []
        for r0 in aug:
            cur = r0
            for b, p in zip(basis, piv):
                if (cur >> p) & 1:
                    cur ^= b
            low = cur & ((1 << n) - 1)
            if low:
                p = (low & -low).bit_length() - 1
                for i in range(len(basis)):
                    if (basis[i] >> p) & 1:
                        basis[i] ^= cur
                basis.append(cur)
                piv.append(p)
            elif cur:
                return None
        x = np.zeros(n, dtype=np.uint8)
        for b, p in zip(basis, piv):
            x[p] = (b >> n) & 1
        assert not ((self.E @ x + rhs) % 2).any()
        return x

    def chain_of(self, v0: np.ndarray, b: np.ndarray) -> np.ndarray:
        return (self.EMB[0] @ v0 + self.EMB[1] @ ((v0 + b) % 2)) % 2

    def is_cover_boundary(self, chain: np.ndarray) -> bool:
        return reduce_int(v2i(chain), self.Sc_basis, self.Sc_piv) == 0

    def rung(self, b_vec: np.ndarray, M: int, deadline: float,
             threads: int, workdir: Path) -> dict:
        bt = self.bt
        n = self.n
        wb = int(b_vec.sum())
        rhs = (self.RHS_OP @ b_vec) % 2
        v0p = self.solve_E(rhs)
        if v0p is None:
            return {"verdict": "PASS", "lane": "vacuous", "w_b": wb, "M": M}
        r = len(self.sector_basis)
        v0p_i = v2i(v0p)
        sectors = []
        for t in range(1 << r):
            ti = 0
            for j in range(r):
                if (t >> j) & 1:
                    ti ^= self.sector_basis[j]
            vt = v0p_i ^ ti
            triv = self.is_cover_boundary(
                self.chain_of(i2v(vt, n), b_vec))
            sectors.append({"v": vt, "trivial": triv})
        nontriv = [s for s in sectors if not s["trivial"]]
        if not nontriv:
            return {"verdict": "PASS", "lane": "all-trivial", "w_b": wb,
                    "M": M}
        bmask = v2i(b_vec)
        viols: list[dict] = []

        def classify(v0_int: int) -> None:
            ov = bin(v0_int & ~bmask).count("1")
            if ov > M - 1:
                return
            v0v = i2v(v0_int, n)
            assert not ((self.E @ v0v + rhs) % 2).any()
            ch = self.chain_of(v0v, b_vec)
            if self.is_cover_boundary(ch):
                return
            assert int(ch.sum()) == wb + 2 * ov
            viols.append({"overflow": ov, "weight": wb + 2 * ov,
                          "v0_hex": f"{v0_int:x}"})

        if M - 1 <= 4:
            lane = f"restricted<={M-1}"
            bcols = [self.E_cols[j] for j in np.nonzero(b_vec)[0]]
            bb, bp = rref_ints(bcols)
            rhs_res = reduce_int(v2i(rhs), bb, bp)
            offb = [j for j in range(n) if not (bmask >> j) & 1]
            red = {j: reduce_int(self.E_cols[j], bb, bp) for j in offb}
            by_val: dict[int, list[int]] = {}
            for j in offb:
                by_val.setdefault(red[j], []).append(j)
            hits_X: set[tuple[int, ...]] = set()
            if rhs_res == 0:
                hits_X.add(())
            if M - 1 >= 1:
                for j in by_val.get(rhs_res, []):
                    hits_X.add((j,))
            if M - 1 >= 2:
                for j1 in offb:
                    for j2 in by_val.get(rhs_res ^ red[j1], []):
                        if j2 > j1:
                            hits_X.add((j1, j2))
            if M - 1 >= 3:
                for j1, j2 in itertools.combinations(offb, 2):
                    for j3 in by_val.get(rhs_res ^ red[j1] ^ red[j2], []):
                        if j3 > j2:
                            hits_X.add((j1, j2, j3))
            if M - 1 >= 4:
                pair_sum: dict[int, list[tuple[int, int]]] = {}
                for j1, j2 in itertools.combinations(offb, 2):
                    pair_sum.setdefault(red[j1] ^ red[j2], []).append((j1, j2))
                for val, prs in pair_sum.items():
                    for j3, j4 in pair_sum.get(rhs_res ^ val, []):
                        for j1, j2 in prs:
                            if j2 < j3:
                                hits_X.add((j1, j2, j3, j4))
            for X in sorted(hits_X):
                cols = list(np.nonzero(b_vec)[0]) + list(X)
                sub = [self.E_cols[j] for j in cols]
                b3: list[int] = []
                p3: list[int] = []
                h3: list[int] = []
                deps: list[int] = []
                for ci, cval in enumerate(sub):
                    cur, h = cval, 1 << ci
                    for bb3, pp3, hh in zip(b3, p3, h3):
                        if (cur >> pp3) & 1:
                            cur ^= bb3
                            h ^= hh
                    if cur:
                        b3.append(cur)
                        p3.append((cur & -cur).bit_length() - 1)
                        h3.append(h)
                    else:
                        deps.append(h)
                cur, hsel = v2i(rhs), 0
                for bb3, pp3, hh in zip(b3, p3, h3):
                    if (cur >> pp3) & 1:
                        cur ^= bb3
                        hsel ^= hh
                if cur:
                    continue
                if len(deps) > 14:
                    return {"verdict": "ABORT", "lane": lane, "w_b": wb,
                            "M": M, "reason": f"kernel 2^{len(deps)}"}
                for kt in range(1 << len(deps)):
                    sel = hsel
                    for j in range(len(deps)):
                        if (kt >> j) & 1:
                            sel ^= deps[j]
                    v0_int = 0
                    for ci, j in enumerate(cols):
                        if (sel >> ci) & 1:
                            v0_int |= 1 << int(j)
                    classify(v0_int)
        else:
            lane = "bz"
            Wp = M - 1 + wb
            r1 = min(Wp // 2, 8)
            r2 = max(Wp - r1 - 1, 0)
            if r1 + r2 + 2 <= Wp:
                return {"verdict": "ABORT", "lane": lane, "w_b": wb, "M": M,
                        "reason": "pair incomplete under cap"}
            bases_v = [i2v(s["v"], n) for s in nontriv]
            for wi, (window, Gs, rr) in enumerate(
                    [(bt.I1, bt.G1, r1), (bt.I2, bt.G2, r2)]):
                bases = [coset_base(Gs, window, bv) for bv in bases_v]
                for bv in bases:
                    classify(v2i(bv))
                res = run_window(bt.binp, f"rung_w{wi}", Gs, bases, rr, Wp,
                                 deadline, threads=threads, workdir=workdir)
                for j, hx in res.pop("hit_rows"):
                    classify(v2i(unpack3(hx, n)))
        if viols:
            return {"verdict": "VIOLATION", "lane": lane, "w_b": wb, "M": M,
                    "violations": viols[:5], "n_viol": len(viols)}
        return {"verdict": "PASS", "lane": lane, "w_b": wb, "M": M}


def rung_pass(
    bt: BaseTools, axis: int, d: int, classes: list[dict],
    deadline: float, threads: int, workdir: Path,
    progress: Optional[Progress] = None,
) -> dict:
    eng = RungEngine(bt, axis, d)
    verdicts: dict[str, int] = {}
    lanes: dict[str, int] = {}
    bad = None
    t0 = time.monotonic()
    for i, c in enumerate(classes):
        M = d - c["weight"] // 2
        if M <= 0:
            r = {"verdict": "PASS", "lane": "heavy-rung"}
        else:
            r = eng.rung(c["vec"], M, deadline, threads, workdir)
        verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
        if r["verdict"] != "PASS" and bad is None:
            bad = {**{k: v for k, v in r.items() if k != "violations"},
                   "class_index": i,
                   "violations": r.get("violations")}
            if r["verdict"] == "VIOLATION":
                break
        if progress and (i + 1) % 500 == 0:
            progress("rungs", done=i + 1, total=len(classes))
    return {"n_classes": len(classes), "verdicts": verdicts, "lanes": lanes,
            "sector_dim": len(eng.sector_basis),
            "all_pass": verdicts.get("PASS", 0) == len(classes),
            "first_bad": bad, "wall_s": round(time.monotonic() - t0, 1),
            "engine": eng}


# -------------------------------------------------------- stage: witness
def witness_lift(eng: RungEngine, u: np.ndarray) -> Optional[np.ndarray]:
    """tau(u): the diagonal cover chain; None unless a nontrivial cycle."""
    ch = (eng.EMB[0] @ u + eng.EMB[1] @ u) % 2
    if ((eng.D1c @ ch) % 2).any():
        return None
    if eng.is_cover_boundary(ch):
        return None
    return ch


# ------------------------------------------------------------- the verdict
def certify(
    orders: tuple[int, ...] | list[int],
    A_text: str,
    B_text: str,
    *,
    budget_s: float = 2400.0,
    threads: int = 8,
    workdir: Optional[Path] = None,
    progress: Optional[Progress] = None,
) -> dict:
    """The composed front-end: detection -> certificates -> d(cover).

    Returns a verdict dict; every claim carries its tier.  Refuses
    explicitly (verdict['status']) when out of scope, and always includes
    the Tandem fallback/cross-check block.
    """
    def emit(stage: str, **kw: Any) -> None:
        if progress:
            progress(stage, **kw)

    t_start = time.monotonic()
    deadline = t_start + budget_s
    wd = workdir or (Path(LAB_DATA) / "certify_runs" /
                     time.strftime("%Y%m%d-%H%M%S"))
    wd.mkdir(parents=True, exist_ok=True)
    G = AbelianGroup(tuple(orders))
    A = Poly.from_string(A_text, G)
    B = Poly.from_string(B_text, G)
    n_cover = 2 * G.cardinality
    k_cover = 2 * kernel_basis(A, B).shape[0]
    out: dict[str, Any] = {
        "input": {"orders": list(orders), "A": A.canonical_string(),
                  "B": B.canonical_string(), "n": n_cover, "k": k_cover},
        "tier": ("certified computational data (counting-invariant "
                 "enumeration + the doubling theorem); NOT kernel-checked"),
        "stages": {},
    }
    if A.weight() % 2 == 0 or B.weight() % 2 == 0:
        out["status"] = "REFUSED"
        out["reason"] = "|A|, |B| must be odd (parity layer of the theory)"
        out["tandem"] = _tandem_block(None, None)
        return out

    emit("detect")
    cands = detect(G, A, B)
    out["stages"]["detect"] = [c.__dict__ for c in cands]
    usable = [c for c in cands if c.R_holds]
    if not usable:
        out["status"] = "FALLBACK"
        out["reason"] = ("no literal-lift axis doubling with (R) detected; "
                         "use the monolithic Tandem lane")
        out["tandem"] = _tandem_block(None, None)
        return out

    out["stages"]["candidate_log"] = []
    refutation = None
    done = False
    for ci, cand in enumerate(usable):
        res = _certify_candidate(cand, deadline, threads, wd, emit, out)
        out["stages"]["candidate_log"].append({
            "candidate": ci, "axis": cand.axis,
            "outcome": (res or {}).get("status",
                                       out["stages"].pop("abort_reason",
                                                         "aborted")),
        })
        if res and res["status"] in ("CERTIFIED", "FLOOR-ONLY"):
            out.update(res)
            done = True
            break
        if res and res["status"] == "DOUBLING-REFUTED":
            refutation = res  # a certified negative; other axes may differ
    if not done:
        if refutation:
            out.update(refutation)
        else:
            out["status"] = "FALLBACK"
            out["reason"] = "no candidate completed in budget"

    d_claim = out.get("distance", {}).get("value")
    out["tandem"] = _tandem_block(d_claim, out.get("distance", {}).get("floor"))
    out["wall_s"] = round(time.monotonic() - t_start, 1)
    out["workdir"] = str(wd)
    return out


LAB_DATA = Path(__file__).resolve().parents[2] / "data"


def _certify_candidate(
    cand: DoublingCandidate, deadline: float, threads: int, wd: Path,
    emit: Progress, out: dict,
) -> Optional[dict]:
    """Run the full pipeline for one detected base; None = try next."""
    Gb = AbelianGroup(cand.base_group)
    Ab = Poly.from_string(cand.base_A, Gb)
    Bb = Poly.from_string(cand.base_B, Gb)
    bt = BaseTools(Gb, Ab, Bb)
    st = out["stages"]
    st["base"] = {"group": list(cand.base_group), "axis": cand.axis,
                  "n": bt.n, "k": cand.k_base, "kappa": bt.kappa}

    emit("d_base")
    sides = {}
    for side in ("X", "Z"):
        sides[side] = d_side_exact(bt, side, WCAP_DBASE, deadline, threads,
                                   wd)
        if "error" in sides[side]:
            st["abort_reason"] = sides[side]["error"]
            return None
    st["d_base"] = {s: {k: v for k, v in r.items() if k != "witness"}
                    for s, r in sides.items()}
    if any("floor_gt" in r for r in sides.values()):
        st["abort_reason"] = f"d_base > {WCAP_DBASE}: outside front-end scope"
        return None
    dX, dZ = sides["X"]["d"], sides["Z"]["d"]
    if dX != dZ:
        st["abort_reason"] = f"side distances differ ({dX} vs {dZ})?!"
        return None
    d = dX
    if d > DBASE_CAP:
        st["abort_reason"] = f"d_base = {d} > cap {DBASE_CAP}"
        return None
    est = _sweep_nodes(bt.kappa, 2 * d - 2)
    if est > NODES_CAP:
        st["abort_reason"] = (
            f"census at W={2*d-2} needs ~{est:.1e} nodes > cap "
            f"{NODES_CAP:.1e}: shard the sweep across machines (the "
            "counting certificate composes) or take the fibering/"
            "trisection theory route")
        return None
    emit("d_base-done", d=d)

    emit("census", W=2 * d - 2)
    cen = census(bt, 2 * d - 2, deadline, threads, wd, progress=emit)
    st["census"] = {k: v for k, v in cen.items() if k != "classes"}
    emit("census-done", n_classes=cen["n_classes"])

    emit("safe-floor", target=2 * d)
    sf = safe_floor(bt, cand.axis, 2 * d, deadline, threads, wd)
    st["safe_floor"] = sf
    if not sf["certified"]:
        st["abort_reason"] = (
            f"safe floor REFUTED at {sf['min_refuted']} < {2*d}: the "
            "doubling fails for this candidate (certified negative)")
        # a genuine refutation is a verdict, not a fallback
        return {"status": "DOUBLING-REFUTED",
                "distance": {"value": None,
                             "floor": None,
                             "upper": None,
                             "statement": (
                                 f"safe-class coset of weight "
                                 f"{sf['min_refuted']} exists; "
                                 f"d < {2*d} on the safe sector")}}
    emit("safe-floor-done")

    emit("rungs", n=cen["n_classes"])
    rp = rung_pass(bt, cand.axis, d, cen["classes"], deadline, threads, wd,
                   progress=emit)
    eng = rp.pop("engine")
    st["rung_pass"] = rp
    if not rp["all_pass"]:
        bad = rp["first_bad"]
        if bad and bad.get("verdict") == "VIOLATION":
            return {"status": "DOUBLING-REFUTED",
                    "distance": {"value": None, "floor": None, "upper": None,
                                 "statement": (
                                     "dangerous-sector violation: cover "
                                     f"logical of weight {bad['violations'][0]['weight']}"
                                     f" over census class {bad['class_index']}")}}
        st["abort_reason"] = f"rung pass aborted: {bad}"
        return None
    emit("rungs-done")

    emit("witness")
    wit = None
    for side in ("X",):
        u = sides[side]["witness"]
        wit = witness_lift(eng, u)
        if wit is not None:
            break
    st["witness"] = {
        "established": wit is not None,
        "weight": int(wit.sum()) if wit is not None else None,
    }
    upper = 2 * d if wit is not None else None
    emit("witness-done", established=wit is not None)

    stmt = (
        f"d = {2*d}" if upper is not None else f"d >= {2*d}"
    ) + (
        " (X-side certificates; Z by the BB transpose duality). "
        "Inputs: base d exact both sides + census (translation-complete) "
        "+ safe-class coset floors + dangerous-sector rung pass + "
        + ("diagonal-lift witness." if upper else "NO lift witness found "
           "— upper side needs the Tandem witness lane.")
    )
    return {"status": "CERTIFIED" if upper else "FLOOR-ONLY",
            "distance": {"value": 2 * d if upper else None,
                         "floor": 2 * d, "upper": upper,
                         "d_base": d, "statement": stmt}}


def scrub_json(o: Any) -> Any:
    """Verdicts carry numpy arrays; make them JSON-serialisable."""
    if isinstance(o, dict):
        return {k: scrub_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [scrub_json(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    return o


def _tandem_block(d_value: Optional[int], floor: Optional[int]) -> dict:
    """How Tandem composes with (or replaces) the certificate."""
    if floor:
        return {
            "role": "cross-check / witness lane",
            "suggested_flags": {"-init-lb": floor, "-cost-step": 2,
                                "-verb": 1},
            "acknowledge": ["-init-lb", "-cost-step"],
            "note": (
                f"passing the certified floor as -init-lb={floor} deletes "
                "the solver's proof phase: the run only has to FIND the "
                f"weight-{floor} witness, then stops at the floor. "
                "-cost-step=2 is sound (coset weight parity). Expect "
                "minutes, vs hours-to-days for the monolithic proof."
            ),
        }
    return {
        "role": "monolithic fallback",
        "suggested_flags": {"-cost-step": 2, "-verb": 1},
        "acknowledge": ["-cost-step"],
        "note": ("no certificate produced — solve the cover directly "
                 "(warning: UNSAT-side cost grows steeply with n)"),
    }
