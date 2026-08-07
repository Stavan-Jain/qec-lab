"""A30 — dangerous-sector rung pass: (M) certified over the census classes.

The teaching-doc dangerous sector (bb-doubling-theorem.tex §"The dangerous
sector"): a nontrivial dangerous cover logical v has sheets (v0, v1) with
v1 = v0 + b, b = p(v) a base stabilizer, and by the slice identity

    wt(v) = wt(b) + 2 * overflow,   overflow := |supp v0 \ supp b|.

(M) at target 2d requires, for every stabilizer b:  wt(b) + 2 m(b) >= 2d,
i.e. m(b) >= M(b) := d - wt(b)/2, where m(b) = min overflow over nontrivial
dangerous logicals with shadow b.  Rungs: b = 0 -> LogicalFloor d (certified,
A30 base floors); wt(b) >= 2d -> slice identity; middle b -> the census
(A28, translation-complete at W = 2d-2) pins b to a class rep, m() is
invariant under cover translations, and THIS script certifies the per-class
bound directly:

  * Sheets (v0, v0+b) form a cover cycle  <=>  E v0 = rhs(b), where
    E = D1_cover (embed0 + embed1) is CLASS-INDEPENDENT (one RREF per cell)
    and rhs(b) = D1_cover embed1(b).  Solutions: v0* + ker E.
  * im S (base stabilizers) <= ker E (diagonal lifts of stabilizers are
    cover boundaries — asserted numerically), and the cover homology class
    of v is constant on v0 + im S.  Sectors = ker E / im S (dim <= k);
    per-sector triviality = one cover-boundary membership test.  Violations
    can only come from nontrivial sectors.
  * Heavy lane (small M): a violation has supp(v0) <= supp(b) u X, |X| <=
    M-1.  Reduce E's columns mod span(E|_{supp b}); a consistent X exists
    iff rhs_res equals a <=(M-1)-subset XOR of reduced columns — found by
    meet-in-the-middle over 180 columns, then reconstructed exactly and
    sector-classified.  Complete by construction.
  * Light lane (large M): multi-offset coset-BZ (a30_coset_bz machinery,
    same kappa=88 windows as the safe floors) over v_t + im S for the
    nontrivial sectors at W' = M-1+wt(b); every hit is overflow-filtered
    and sector-checked.  Complete by the two-window counting invariant.

Validation: f2a6:y (d=8) must pass 113/113 — DangerousFloorNZ 16 is a Lean
theorem there — plus a soundness control that hand-builds a genuine
dangerous logical (tau(u) + S~(z~)) and confirms the checker FINDS it at an
inflated target.

Usage:
  uv run python scripts/a30_rung_pass.py validate           # f2a6 113/113 + control
  uv run python scripts/a30_rung_pass.py scope              # 50-class samples
  uv run python scripts/a30_rung_pass.py full [--only CID]  # the 3 cells
Outputs: data/a30/rungs_*.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
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

from a30_coset_bz import (  # noqa: E402
    DATA, build_kernel, code_rows, coset_base, disjoint_info_sets,
    poly_pair, run_window, unpack3,
)

# ------------------------------------------------------------ bitint helpers
def v2i(v: np.ndarray) -> int:
    x = 0
    for i in np.nonzero(v)[0]:
        x |= 1 << int(i)
    return x


def i2v(x: int, n: int) -> np.ndarray:
    return np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)


def rref_ints(rows: list[int]) -> tuple[list[int], list[int]]:
    """RREF over GF(2) on bit-int rows; returns (rows, pivot bit positions)."""
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


# --------------------------------------------------------------- cell object
class RungCell:
    """Everything precomputable for one (code, axis) cover instance."""

    def __init__(self, cid: str, group: tuple[int, int], As: str, Bs: str,
                 axis: int, d_base: int):
        self.cid, self.axis, self.d = cid, axis, d_base
        self.G = AbelianGroup(group)
        self.A = Poly.from_string(As, self.G)
        self.B = Poly.from_string(Bs, self.G)
        gc = list(group)
        gc[axis] *= 2
        self.Gc = AbelianGroup(tuple(gc))
        self.Ac = Poly.from_string(As, self.Gc)
        self.Bc = Poly.from_string(Bs, self.Gc)
        ng, nc = self.G.cardinality, self.Gc.cardinality
        self.n = 2 * ng          # base pair length (180)
        self.nc = 2 * nc         # cover pair length (360)
        MA = circulant(self.A).astype(np.uint8) % 2
        MB = circulant(self.B).astype(np.uint8) % 2
        MAc = circulant(self.Ac).astype(np.uint8) % 2
        MBc = circulant(self.Bc).astype(np.uint8) % 2
        self.D1c = np.concatenate([MBc, MAc], axis=1) % 2  # cover check
        # sheet embeddings: base index -> cover index, per sheet
        ell = group[axis]
        emb = [np.zeros((nc, ng), dtype=np.uint8) for _ in range(2)]
        for i, e in enumerate(self.G):
            for s in (0, 1):
                ec = list(e)
                ec[axis] += s * ell
                emb[s][self.Gc.index(tuple(ec)), i] = 1
        z0 = np.zeros_like(emb[0])
        self.EMB = [np.block([[emb[s], z0], [z0, emb[s]]]) for s in (0, 1)]
        # E (class-independent) and the rhs operator
        self.E = (self.D1c @ ((self.EMB[0] + self.EMB[1]) % 2)) % 2
        self.RHS_OP = (self.D1c @ self.EMB[1]) % 2
        # base stabilizer rows / cover stabilizer rows
        self.MS = code_rows(self.A, self.B)          # 90 x 180, rank 88
        self.MSc = code_rows(self.Ac, self.Bc)       # 180 x 360, rank 176
        self.S_basis, self.S_piv = rref_ints(
            [v2i(r) for r in self.MS])
        self.Sc_basis, self.Sc_piv = rref_ints(
            [v2i(r) for r in self.MSc])
        # E as row ints (for per-class solving) and column ints (heavy lane)
        self.E_rows = [v2i(self.E[i]) for i in range(self.E.shape[0])]
        self.E_cols = [v2i(self.E[:, j]) for j in range(self.n)]
        self.kerE = self._kernel_basis()
        # assertions: im S <= ker E; b -> E b == 0 for stabilizers
        for idx in (0, 1, 7):
            r = self.MS[idx]
            assert not ((self.E @ r) % 2).any(), "im S not <= ker E"
        # sector basis: ker E mod im S
        sec = []
        aug_basis = list(self.S_basis)
        aug_piv = list(self.S_piv)
        for kv in self.kerE:
            x = reduce_int(kv, aug_basis, aug_piv)
            if x:
                p = (x & -x).bit_length() - 1
                aug_basis.append(x)
                aug_piv.append(p)
                sec.append(kv)
        self.sector_basis = sec  # dim r <= k
        # windows for the BZ lane (same construction as the safe floors)
        (self.I1, self.G1, self.I2, self.G2, self.kappa) = \
            disjoint_info_sets(self.MS)
        self.binp = build_kernel()

    # -- linear algebra helpers -------------------------------------------
    def _kernel_basis(self) -> list[int]:
        M = self.E.copy()
        n = M.shape[1]
        # numpy rref
        piv = []
        r = 0
        for c in range(n):
            rows = np.nonzero(M[r:, c])[0]
            if len(rows) == 0:
                continue
            M[[r, r + rows[0]]] = M[[r + rows[0], r]]
            for i in np.nonzero(M[:, c])[0]:
                if i != r:
                    M[i] ^= M[r]
            piv.append(c)
            r += 1
            if r == M.shape[0]:
                break
        Mr = M[:r]
        free = [c for c in range(n) if c not in set(piv)]
        out = []
        for c in free:
            v = np.zeros(n, dtype=np.uint8)
            v[c] = 1
            for i, pc in enumerate(piv):
                if Mr[i, c]:
                    v[pc] ^= 1
            assert not ((self.E @ v) % 2).any()
            out.append(v2i(v))
        return out

    def solve_E(self, rhs: np.ndarray):
        """Particular solution of E v0 = rhs, or None if inconsistent.

        Fresh augmented row-RREF per class (rhs bit rides at position n);
        free variables set to 0, so x[pivot] = the reduced row's rhs bit."""
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
                return None  # 0 = 1: inconsistent
        x = np.zeros(n, dtype=np.uint8)
        for b, p in zip(basis, piv):
            x[p] = (b >> n) & 1
        assert not ((self.E @ x + rhs) % 2).any()
        return x

    def chain_of(self, v0: np.ndarray, b: np.ndarray) -> np.ndarray:
        return (self.EMB[0] @ v0 + self.EMB[1] @ ((v0 + b) % 2)) % 2

    def is_cover_boundary(self, chain: np.ndarray) -> bool:
        return reduce_int(v2i(chain), self.Sc_basis, self.Sc_piv) == 0

    # -- per-class rung ----------------------------------------------------
    def rung(self, b_vec: np.ndarray, M: int, deadline: float,
             bz_cap_r: int = 8) -> dict:
        """Certify: no nontrivial dangerous logical with shadow b has
        overflow <= M-1.  Returns verdict dict."""
        n = self.n
        wb = int(b_vec.sum())
        rhs = (self.RHS_OP @ b_vec) % 2
        v0p = self.solve_E(rhs)
        if v0p is None:
            return {"verdict": "PASS", "lane": "vacuous", "w_b": wb, "M": M}
        # sectors + triviality
        r = len(self.sector_basis)
        sectors = []
        v0p_i = v2i(v0p)
        for t in range(1 << r):
            ti = 0
            for j in range(r):
                if (t >> j) & 1:
                    ti ^= self.sector_basis[j]
            vt = v0p_i ^ ti
            vt_vec = i2v(vt, n)
            triv = self.is_cover_boundary(self.chain_of(vt_vec, b_vec))
            sectors.append({"t": t, "v": vt, "trivial": triv})
        nontriv = [s for s in sectors if not s["trivial"]]
        if not nontriv:
            return {"verdict": "PASS", "lane": "all-trivial", "w_b": wb,
                    "M": M, "sectors": len(sectors)}
        bmask = v2i(b_vec)
        viols = []

        def classify_and_check(v0_int: int) -> None:
            ov = bin(v0_int & ~bmask).count("1")
            if ov > M - 1:
                return
            v0v = i2v(v0_int, n)
            assert not ((self.E @ v0v + rhs) % 2).any(), "hit not a solution"
            ch = self.chain_of(v0v, b_vec)
            if self.is_cover_boundary(ch):
                return
            wt = int(ch.sum())
            assert wt == wb + 2 * ov
            viols.append({"overflow": ov, "weight": wt,
                          "v0_hex": f"{v0_int:x}"})

        if M - 1 <= 4:
            lane = f"restricted<= {M-1}"
            # reduce E columns mod span(E cols on supp b)
            bcols = [self.E_cols[j] for j in np.nonzero(b_vec)[0]]
            bb, bp = rref_ints(bcols)
            rhs_res = reduce_int(v2i(rhs), bb, bp)
            offb = [j for j in range(n) if not (bmask >> j) & 1]
            red = {j: reduce_int(self.E_cols[j], bb, bp) for j in offb}
            hits_X: set[tuple[int, ...]] = set()
            if rhs_res == 0:
                hits_X.add(())
            by_val: dict[int, list[int]] = {}
            for j in offb:
                by_val.setdefault(red[j], []).append(j)
            if M - 1 >= 1:
                for j in by_val.get(rhs_res, []):
                    hits_X.add((j,))
            if M - 1 >= 2:
                for j1 in offb:
                    tgt = rhs_res ^ red[j1]
                    for j2 in by_val.get(tgt, []):
                        if j2 > j1:
                            hits_X.add((j1, j2))
            if M - 1 >= 3:
                for j1, j2 in itertools.combinations(offb, 2):
                    tgt = rhs_res ^ red[j1] ^ red[j2]
                    for j3 in by_val.get(tgt, []):
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
            # reconstruct each consistent support set exactly
            for X in sorted(hits_X):
                cols = list(np.nonzero(b_vec)[0]) + list(X)
                sub = [self.E_cols[j] for j in cols]
                # solve sum x_j col_j = rhs on these columns
                aug, ap = [], []
                sol_bits = {}
                target = v2i(rhs)
                basis2, piv2 = [], []
                hist = []
                for ci, cval in enumerate(sub):
                    cur, h = cval, 1 << ci
                    for bb2, pp2, hh in zip(basis2, piv2, hist):
                        if (cur >> pp2) & 1:
                            cur ^= bb2
                            h ^= hh
                    if cur:
                        basis2.append(cur)
                        piv2.append((cur & -cur).bit_length() - 1)
                        hist.append(h)
                cur, hsel = target, 0
                for bb2, pp2, hh in zip(basis2, piv2, hist):
                    if (cur >> pp2) & 1:
                        cur ^= bb2
                        hsel ^= hh
                if cur:
                    continue
                # kernel of the restricted system
                kfree = []
                span_h, span_p = [], []
                for ci, cval in enumerate(sub):
                    x = reduce_int(cval, span_h, span_p)
                    if x:
                        span_h.append(x)
                        span_p.append((x & -x).bit_length() - 1)
                # kernel via dependency detection
                deps = []
                b3, p3, h3 = [], [], []
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
                if len(deps) > 14:
                    return {"verdict": "ABORT", "lane": lane,
                            "reason": f"restricted kernel 2^{len(deps)}",
                            "w_b": wb, "M": M}
                for kt in range(1 << len(deps)):
                    sel = hsel
                    for j in range(len(deps)):
                        if (kt >> j) & 1:
                            sel ^= deps[j]
                    v0_int = 0
                    for ci, j in enumerate(cols):
                        if (sel >> ci) & 1:
                            v0_int |= 1 << int(j)
                    classify_and_check(v0_int)
        else:
            lane = "bz"
            Wp = M - 1 + wb
            r1 = min(Wp // 2, bz_cap_r)
            r2 = max(Wp - r1 - 1, 0)
            bases_v = [i2v(s["v"], n) for s in nontriv]
            nodes = []
            for wi, (window, Gs) in enumerate(
                    [(self.I1, self.G1), (self.I2, self.G2)]):
                rr = r1 if wi == 0 else r2
                bases = [coset_base(Gs, window, bv) for bv in bases_v]
                for bv in bases:
                    classify_and_check(v2i(bv))  # the S=empty element
                res = run_window(self.binp, f"{self.cid}_rung_w{wi}", Gs,
                                 bases, rr, Wp, deadline)
                nodes.append(res["nodes"])
                for j, hx in res.pop("hit_rows"):
                    classify_and_check(v2i(unpack3(hx, n)))
            if r1 + r2 + 2 <= Wp:
                return {"verdict": "ABORT", "lane": lane, "w_b": wb, "M": M,
                        "reason": "bz pair incomplete under cap"}
        if viols:
            return {"verdict": "VIOLATION", "lane": lane, "w_b": wb, "M": M,
                    "violations": viols[:5], "n_viol": len(viols)}
        return {"verdict": "PASS", "lane": lane, "w_b": wb, "M": M,
                "sectors_nontrivial": len(nontriv)}


# ----------------------------------------------------------------- census IO
def census_vec(cell: RungCell, entry: dict) -> np.ndarray:
    ng = cell.G.cardinality
    v = np.zeros(cell.n, dtype=np.uint8)
    for g in entry["u_support"]:
        v[cell.G.index(tuple(g))] = 1
    for g in entry["v_support"]:
        v[ng + cell.G.index(tuple(g))] = 1
    assert int(v.sum()) == entry["weight"]
    # must be a base stabilizer
    assert reduce_int(v2i(v), cell.S_basis, cell.S_piv) == 0, \
        "census rep not a stabilizer"
    return v


CELLS = {
    "f2a6:y": dict(group=(5, 15), A="1 + y + x",
                   B="x*y^6 + x*y^10 + x^2*y^12", axis=1, d=8,
                   census="census_f2a6.json"),
    "37a70e02:x": dict(group=(15, 6), A="1 + y + x",
                       B="y^4 + x + x^11*y^2", axis=0, d=10,
                       census="census_docket37.json"),
    "5e50a976:x": dict(group=(15, 6), A="1 + y + x",
                       B="y^4 + x^8*y^2 + x^13", axis=0, d=10,
                       census="census_docket5e.json"),
    "5e50a976:y": dict(group=(15, 6), A="1 + y + x",
                       B="y^4 + x^8*y^2 + x^13", axis=1, d=10,
                       census="census_docket5e.json"),
}


def run_cell(name: str, sample: int | None, budget_s: float,
             seed: int = 0) -> dict:
    spec = CELLS[name]
    t0 = time.monotonic()
    cell = RungCell(name.split(":")[0], spec["group"], spec["A"], spec["B"],
                    spec["axis"], spec["d"])
    cen = json.load(open(DATA.parent / "a28" / spec["census"]))
    classes = cen["classes"]
    idx = list(range(len(classes)))
    if sample:
        random.Random(seed).shuffle(idx)
        # stratified: keep relative weights by sorting sample by weight
        idx = sorted(idx[:sample])
    deadline = time.monotonic() + budget_s
    out = {"cell": name, "d": spec["d"], "target": 2 * spec["d"],
           "n_classes_total": len(classes), "n_run": len(idx),
           "sector_dim": len(cell.sector_basis),
           "kerE_dim": len(cell.kerE), "per_lane": {}, "results": []}
    lanes: dict[str, list[float]] = {}
    worst = None
    for i in idx:
        e = classes[i]
        b = census_vec(cell, e)
        M = spec["d"] - e["weight"] // 2
        tc0 = time.monotonic()
        if M <= 0:
            r = {"verdict": "PASS", "lane": "heavy-rung", "w_b": e["weight"],
                 "M": M}
        else:
            r = cell.rung(b, M, deadline)
        dt = time.monotonic() - tc0
        r["class_index"] = i
        r["secs"] = round(dt, 4)
        lanes.setdefault(r["lane"], []).append(dt)
        out["results"].append(r)
        if r["verdict"] != "PASS":
            worst = r
            if r["verdict"] == "VIOLATION":
                break
        if time.monotonic() > deadline:
            out["aborted"] = f"budget after {len(out['results'])} classes"
            break
    out["per_lane"] = {k: {"n": len(v), "tot_s": round(sum(v), 2),
                           "max_s": round(max(v), 4)}
                       for k, v in lanes.items()}
    out["verdicts"] = {}
    for r in out["results"]:
        out["verdicts"][r["verdict"]] = out["verdicts"].get(r["verdict"], 0) + 1
    out["wall_s"] = round(time.monotonic() - t0, 1)
    out["all_pass"] = (out["verdicts"].get("PASS", 0) == len(out["results"])
                       and "aborted" not in out)
    if worst:
        out["first_bad"] = worst
    return out


# ------------------------------------------------------------------ control
def soundness_control() -> dict:
    """Hand-build a nontrivial dangerous logical on f2a6:y and confirm the
    checker finds it at an inflated target."""
    spec = CELLS["f2a6:y"]
    cell = RungCell("f2a6ctl", spec["group"], spec["A"], spec["B"],
                    spec["axis"], spec["d"])
    n = cell.n
    # a base X-logical u of weight 8: BZ over the logical cosets
    from a30_coset_bz import base_floor  # local import to reuse machinery
    # cheap direct search: enumerate logical cosets at W'=8 via run_window
    MA = circulant(cell.A).astype(np.uint8) % 2
    MB = circulant(cell.B).astype(np.uint8) % 2
    H = np.concatenate([MB, MA], axis=1) % 2
    # kernel reps independent mod im S
    Hn = H.copy()
    piv = []
    rr = 0
    for c in range(n):
        rows = np.nonzero(Hn[rr:, c])[0]
        if len(rows) == 0:
            continue
        Hn[[rr, rr + rows[0]]] = Hn[[rr + rows[0], rr]]
        for i in np.nonzero(Hn[:, c])[0]:
            if i != rr:
                Hn[i] ^= Hn[rr]
        piv.append(c)
        rr += 1
    Hn = Hn[:rr]
    free = [c for c in range(n) if c not in set(piv)]
    reps = []
    aug_b, aug_p = list(cell.S_basis), list(cell.S_piv)
    for c in free:
        v = np.zeros(n, dtype=np.uint8)
        v[c] = 1
        for i, pc in enumerate(piv):
            if Hn[i, c]:
                v[pc] ^= 1
        x = reduce_int(v2i(v), aug_b, aug_p)
        if x:
            p = (x & -x).bit_length() - 1
            aug_b.append(x)
            aug_p.append(p)
            reps.append(v)
    # find a weight-d logical: multi-offset BZ at W'=d over all 2^k-1 combos
    combos = []
    for t in range(1, 1 << len(reps)):
        L = np.zeros(n, dtype=np.uint8)
        for j in range(len(reps)):
            if (t >> j) & 1:
                L ^= reps[j]
        combos.append(L)
    deadline = time.monotonic() + 300
    best = None
    for wi, (window, Gs) in enumerate([(cell.I1, cell.G1),
                                       (cell.I2, cell.G2)]):
        bases = [coset_base(Gs, window, L) for L in combos]
        res = run_window(cell.binp, "ctl_log_w%d" % wi, Gs, bases,
                         cell.d // 2, cell.d, deadline)
        for j, hx in res.pop("hit_rows"):
            u = unpack3(hx, n)
            if best is None or u.sum() < best.sum():
                best = u
        if best is not None:
            break
    assert best is not None and int(best.sum()) == cell.d, "no weight-d logical found"
    u = best
    # v = tau(u) + S~(z~): z~ = a single cover monomial
    z = np.zeros(cell.Gc.cardinality, dtype=np.uint8)
    z[3] = 1
    Sz = np.concatenate([
        (circulant(cell.Ac).astype(np.uint8) @ z) % 2,
        (circulant(cell.Bc).astype(np.uint8) @ z) % 2])
    vch = ((cell.EMB[0] @ u + cell.EMB[1] @ u) % 2 + Sz) % 2
    assert not ((cell.D1c @ vch) % 2).any(), "control chain not a cycle"
    assert not cell.is_cover_boundary(vch), "control chain trivial"
    # its shadow and overflow
    p_fold = np.zeros(n, dtype=np.uint8)
    # fold: p(v) via EMB^T
    p_fold = ((cell.EMB[0].T @ vch) + (cell.EMB[1].T @ vch)) % 2
    wb = int(p_fold.sum())
    v0 = (cell.EMB[0].T @ vch) % 2
    ov = int(((v0 == 1) & (p_fold == 0)).sum())
    assert int(vch.sum()) == wb + 2 * ov, "slice identity violated"
    # the checker must FIND a violation at target M_test = ov+1
    r = cell.rung(p_fold, ov + 1, time.monotonic() + 600)
    found = (r["verdict"] == "VIOLATION"
             and any(v["overflow"] <= ov for v in r["violations"]))
    return {"control": "tau(u)+S(z)", "wt_v": int(vch.sum()), "w_b": wb,
            "overflow": ov, "checker": r["verdict"],
            "found_at_or_below": found, "ok": bool(found)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["validate", "scope", "full"])
    ap.add_argument("--only", default=None)
    ap.add_argument("--budget", type=float, default=3600.0)
    args = ap.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)

    if args.stage == "validate":
        ctl = soundness_control()
        print("control:", json.dumps(ctl))
        res = run_cell("f2a6:y", None, args.budget)
        res["control"] = ctl
        (DATA / "rungs_f2a6_validate.json").write_text(
            json.dumps(res, indent=1))
        print(f"f2a6:y — {res['verdicts']} of {res['n_run']} "
              f"in {res['wall_s']}s lanes={json.dumps(res['per_lane'])}")
        ok = res["all_pass"] and ctl["ok"]
        print("VALIDATE:", "OK" if ok else "** FAIL **")
        sys.exit(0 if ok else 1)

    names = [n for n in CELLS if n != "f2a6:y"
             and (not args.only or n.startswith(args.only))]
    for name in names:
        res = run_cell(name, 50 if args.stage == "scope" else None,
                       args.budget)
        tag = name.replace(":", "_")
        suffix = "_scope" if args.stage == "scope" else ""
        (DATA / f"rungs_{tag}{suffix}.json").write_text(
            json.dumps(res, indent=1))
        print(f"{name} — {res['verdicts']} of {res['n_run']}/"
              f"{res['n_classes_total']} in {res['wall_s']}s "
              f"lanes={json.dumps(res['per_lane'])}"
              + (" ABORTED " + res["aborted"] if "aborted" in res else ""))


if __name__ == "__main__":
    main()
