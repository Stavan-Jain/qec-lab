"""Scoping: does the A30 rung framework fit Bravyi [[360,12,<=24]]?

Adapts a30_rung_pass to the y-deck BY (30,3) -> C (30,6) — a TWISTED lift
(C's B = y^3+x^25+x^26 vs BY's 1+x^25+x^26; y^3 = the deck sigma, which acts
trivially on sigma-invariant chains, so tau(S_BY z) = S~_C(tau z) still
holds and the sector machinery is unchanged).  Stages:

  m12   the 7 banked M12 classes at target 12  — must reproduce A19 §8
  m24   the 8,310 banked classes (bands <= 20) at target 24 — reproduces
        A19 §9's SAT (M)@24 floors as counting certificates + measures
        per-stratum cost (the estimate anchors); the 151 stalled band-22
        records run too, labeled PARTIAL (band 22 census is incomplete).

Worktree-scoped scratch; no registry claim.  Census data copied read-only
from the main checkout into data/a19_scope/.
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

from a30_coset_bz import (  # noqa: E402
    build_kernel, code_rows, coset_base, disjoint_info_sets, run_window,
    unpack3,
)
from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402

DATA = LAB / "data" / "a19_scope"


class BravyiRungCell:
    """BY (30,3) --y--> C (30,6), explicit (twisted) cover polynomials."""

    def __init__(self):
        self.d_base = 12  # d(BY) exact, A19 §1 (solver-certified)
        self.G = AbelianGroup((30, 3))
        self.A = Poly.from_string("x^9 + y + y^2", self.G)
        self.B = Poly.from_string("1 + x^25 + x^26", self.G)
        self.Gc = AbelianGroup((30, 6))
        self.Ac = Poly.from_string("x^9 + y + y^2", self.Gc)
        self.Bc = Poly.from_string("y^3 + x^25 + x^26", self.Gc)
        ng, nc = self.G.cardinality, self.Gc.cardinality  # 90, 180
        self.n, self.nc = 2 * ng, 2 * nc                  # 180, 360
        # LAB convention (bb_check_matrices): H_X = [M_A | M_B], stabilizers
        # = its ROW space; X-logicals = ker H_Z \ rowspace H_X with
        # H_Z = [M_B^T | M_A^T].  (The a19 census lives here — transpose of
        # the a30 convention.)
        MA = circulant(self.A).astype(np.uint8) % 2
        MB = circulant(self.B).astype(np.uint8) % 2
        MAc = circulant(self.Ac).astype(np.uint8) % 2
        MBc = circulant(self.Bc).astype(np.uint8) % 2
        self.D1c = np.concatenate([MBc.T, MAc.T], axis=1) % 2  # H_Z(C)
        emb = [np.zeros((nc, ng), dtype=np.uint8) for _ in range(2)]
        for i, e in enumerate(self.G):
            for s in (0, 1):
                ec = (e[0], e[1] + 3 * s)
                emb[s][self.Gc.index(ec), i] = 1
        z0 = np.zeros_like(emb[0])
        self.EMB = [np.block([[emb[s], z0], [z0, emb[s]]]) for s in (0, 1)]
        self.E = (self.D1c @ ((self.EMB[0] + self.EMB[1]) % 2)) % 2
        self.RHS_OP = (self.D1c @ self.EMB[1]) % 2
        self.MS = np.concatenate([MA, MB], axis=1) % 2    # H_X(BY) rows
        self.MSc = np.concatenate([MAc, MBc], axis=1) % 2  # H_X(C) rows
        self.S_basis, self.S_piv = rref_ints([v2i(r) for r in self.MS])
        self.Sc_basis, self.Sc_piv = rref_ints([v2i(r) for r in self.MSc])
        self.E_rows = [v2i(self.E[i]) for i in range(self.E.shape[0])]
        self.E_cols = [v2i(self.E[:, j]) for j in range(self.n)]
        # structural asserts: im S_BY <= ker E (twist-invariance of tau)
        for idx in (0, 5, 41):
            assert not ((self.E @ self.MS[idx]) % 2).any(), \
                "im S not <= ker E — twist argument fails"
        self.kerE = self._kernel_basis()
        sec, aug_b, aug_p = [], list(self.S_basis), list(self.S_piv)
        for kv in self.kerE:
            x = reduce_int(kv, aug_b, aug_p)
            if x:
                aug_b.append(x)
                aug_p.append((x & -x).bit_length() - 1)
                sec.append(kv)
        self.sector_basis = sec
        # precompute sector shift chains (b-independent):
        # chain(v_t, b) = chain(v0*, b) XOR emb01(t)
        M01 = (self.EMB[0] + self.EMB[1]) % 2
        self.sector_chain_shift = []
        for j, t in enumerate(self.sector_basis):
            self.sector_chain_shift.append(v2i((M01 @ i2v(t, self.n)) % 2))
        (self.I1, self.G1, self.I2, self.G2, self.kappa) = \
            disjoint_info_sets(self.MS)
        self.binp = build_kernel()

    def _kernel_basis(self):
        M = self.E.copy()
        n = M.shape[1]
        piv, r = [], 0
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
        Mr = M[:r]
        out = []
        for c in [c for c in range(n) if c not in set(piv)]:
            v = np.zeros(n, dtype=np.uint8)
            v[c] = 1
            for i, pc in enumerate(piv):
                if Mr[i, c]:
                    v[pc] ^= 1
            out.append(v2i(v))
        return out

    def solve_E(self, rhs):
        n = self.n
        basis, piv = [], []
        for i in range(len(self.E_rows)):
            cur = self.E_rows[i] | (int(rhs[i]) << n)
            for b, p in zip(basis, piv):
                if (cur >> p) & 1:
                    cur ^= b
            low = cur & ((1 << n) - 1)
            if low:
                p = (low & -low).bit_length() - 1
                for k in range(len(basis)):
                    if (basis[k] >> p) & 1:
                        basis[k] ^= cur
                basis.append(cur)
                piv.append(p)
            elif cur:
                return None
        x = np.zeros(n, dtype=np.uint8)
        for b, p in zip(basis, piv):
            x[p] = (b >> n) & 1
        assert not ((self.E @ x + rhs) % 2).any()
        return x

    def chain_int(self, v0_int: int, b_vec: np.ndarray) -> int:
        v0 = i2v(v0_int, self.n)
        ch = (self.EMB[0] @ v0 + self.EMB[1] @ ((v0 + b_vec) % 2)) % 2
        return v2i(ch)

    def rung(self, b_vec: np.ndarray, M: int, deadline: float) -> dict:
        n = self.n
        wb = int(b_vec.sum())
        rhs = (self.RHS_OP @ b_vec) % 2
        v0p = self.solve_E(rhs)
        if v0p is None:
            return {"verdict": "PASS", "lane": "vacuous", "w_b": wb, "M": M}
        v0p_i = v2i(v0p)
        base_chain = self.chain_int(v0p_i, b_vec)
        r = len(self.sector_basis)
        nontriv = []
        for t in range(1 << r):
            ch = base_chain
            ti = 0
            for j in range(r):
                if (t >> j) & 1:
                    ch ^= self.sector_chain_shift[j]
                    ti ^= self.sector_basis[j]
            if reduce_int(ch, self.Sc_basis, self.Sc_piv) != 0:
                nontriv.append(v0p_i ^ ti)
        if not nontriv:
            return {"verdict": "PASS", "lane": "all-trivial", "w_b": wb,
                    "M": M}
        bmask = v2i(b_vec)
        viols = []

        def check(v0_int: int):
            ov = bin(v0_int & ~bmask).count("1")
            if ov > M - 1:
                return
            v0v = i2v(v0_int, n)
            assert not ((self.E @ v0v + rhs) % 2).any()
            ch = self.chain_int(v0_int, b_vec)
            if reduce_int(ch, self.Sc_basis, self.Sc_piv) == 0:
                return
            wt = bin(ch).count("1")
            assert wt == wb + 2 * ov, "slice identity violated"
            viols.append({"overflow": ov, "weight": wt,
                          "v0_hex": f"{v0_int:x}"})

        if M - 1 <= 6:
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
            import itertools
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
            if M - 1 >= 5:
                # size-5: triple (j1<j2<j3) + pair (j4<j5), j4 > j3
                for j1, j2, j3 in itertools.combinations(offb, 3):
                    tgt = rhs_res ^ red[j1] ^ red[j2] ^ red[j3]
                    for j4, j5 in pair_sum.get(tgt, []):
                        if j4 > j3:
                            hits_X.add((j1, j2, j3, j4, j5))
            if M - 1 >= 6:
                # size-6: triple + triple (min of 2nd > max of 1st)
                tri_sum: dict[int, list[tuple[int, int, int]]] = {}
                for tri in itertools.combinations(offb, 3):
                    tri_sum.setdefault(
                        red[tri[0]] ^ red[tri[1]] ^ red[tri[2]], []
                    ).append(tri)
                for t1 in itertools.combinations(offb, 3):
                    tgt = rhs_res ^ red[t1[0]] ^ red[t1[1]] ^ red[t1[2]]
                    for t2 in tri_sum.get(tgt, []):
                        if t2[0] > t1[2]:
                            hits_X.add(t1 + t2)
            for X in sorted(hits_X):
                cols = list(np.nonzero(b_vec)[0]) + list(X)
                sub = [self.E_cols[j] for j in cols]
                b3, p3, h3 = [], [], []
                deps = []
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
                    check(v0_int)
        else:
            lane = "bz"
            Wp = M - 1 + wb
            r1 = Wp // 2
            r2 = max(Wp - r1 - 1, 0)
            if len(nontriv) > 256:
                return {"verdict": "ABORT", "lane": lane, "w_b": wb, "M": M,
                        "reason": f"{len(nontriv)} offsets > 256"}
            bases_v = [i2v(v, n) for v in nontriv]
            for wi, (window, Gs) in enumerate(
                    [(self.I1, self.G1), (self.I2, self.G2)]):
                rr = r1 if wi == 0 else r2
                bases = [coset_base(Gs, window, bv) for bv in bases_v]
                for bv in bases:
                    check(v2i(bv))
                res = run_window(self.binp, f"bravyi_rung_w{wi}", Gs, bases,
                                 rr, Wp, deadline)
                for j, hx in res.pop("hit_rows"):
                    check(v2i(unpack3(hx, n)))
        if viols:
            return {"verdict": "VIOLATION", "lane": lane, "w_b": wb, "M": M,
                    "violations": viols[:5], "n_viol": len(viols)}
        return {"verdict": "PASS", "lane": lane, "w_b": wb, "M": M,
                "sectors_nontrivial": len(nontriv)}


def load_classes(path):
    out = []
    for line in open(path):
        r = json.loads(line)
        if "b_support" in r:
            out.append(r)
    return out


def main():
    t0 = time.monotonic()
    cell = BravyiRungCell()
    print(f"cell built {time.monotonic()-t0:.1f}s: kerE={len(cell.kerE)} "
          f"sectors=2^{len(cell.sector_basis)} kappa={cell.kappa}")
    # convention check: every m12 b must be a BY stabilizer
    m12 = load_classes(DATA / "m12_census_classes.jsonl")
    ok = 0
    for e in m12:
        b = np.zeros(cell.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        if reduce_int(v2i(b), cell.S_basis, cell.S_piv) == 0:
            ok += 1
    print(f"m12 convention check: {ok}/{len(m12)} stabilizers")
    assert ok == len(m12), "index convention mismatch"

    # stage 1: reproduce M12 (target 12)
    res12 = []
    for e in m12:
        b = np.zeros(cell.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        M = (12 - e["w"] + 1) // 2
        tc = time.monotonic()
        r = cell.rung(b, M, time.monotonic() + 600)
        r["secs"] = round(time.monotonic() - tc, 3)
        res12.append(r)
    v12 = {}
    for r in res12:
        v12[r["verdict"]] = v12.get(r["verdict"], 0) + 1
    print(f"M12 reproduction: {v12} "
          f"({sum(r['secs'] for r in res12):.2f}s total)")

    # stage 2: (M)@24 over the banked census
    m24 = load_classes(DATA / "m24_census_classes.jsonl")
    lanes: dict[str, list[float]] = {}
    verd: dict[str, int] = {}
    results = []
    deadline = time.monotonic() + 3 * 3600
    for i, e in enumerate(m24):
        b = np.zeros(cell.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        M = (24 - e["w"] + 1) // 2
        tc = time.monotonic()
        r = cell.rung(b, M, deadline)
        dt = time.monotonic() - tc
        r["band"] = e["band"]
        r["partial_band22"] = e["band"] == 22
        lanes.setdefault(f"w{e['w']}", []).append(dt)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        if r["verdict"] != "PASS":
            r["class_index"] = i
            results.append(r)
            print("NON-PASS:", json.dumps(r)[:400])
            if r["verdict"] == "VIOLATION" and not r["partial_band22"]:
                break
        if i % 1000 == 0:
            print(f"  ...{i}/{len(m24)} {time.monotonic()-t0:.0f}s")
    print(f"(M)@24: {verd} of {len(m24)}")
    for k in sorted(lanes, key=lambda s: int(s[1:])):
        v = lanes[k]
        print(f"  {k}: n={len(v)} tot={sum(v):.1f}s max={max(v):.3f}s")
    out = {"kerE": len(cell.kerE), "sector_dim": len(cell.sector_basis),
           "m12": v12, "m24": verd,
           "lanes": {k: {"n": len(v), "tot_s": round(sum(v), 2)}
                     for k, v in lanes.items()},
           "non_pass": results, "wall_s": round(time.monotonic() - t0, 1)}
    (DATA / "scope_results.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA/'scope_results.json'}")


if __name__ == "__main__":
    main()
