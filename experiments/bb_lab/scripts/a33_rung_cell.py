"""A33: the per-shadow rung engine for the class-Y tower (generic same-axis
Z2-deck port of a30_rung_pass / scope_bravyi_rung).

Two rung species per deck cell (base --y--> cover):

  rung(b, M)       b a BASE STABILIZER (dangerous sector): PASS iff no
                   NONTRIVIAL cover logical v with p(v) = b has overflow
                   < M (so |v| = |b| + 2*overflow >= |b| + 2M).  Full
                   sector dispatch (ker E / im S = H1(base), 2^k sectors;
                   trivial sectors excluded).  Lanes: vacuous /
                   all-trivial / restricted MITM (M-1 <= 6) / coset-BZ.

  seam_rung(w, M)  w a coset element with NONZERO H1 class (safe sector):
                   every cover cycle over w is automatically a nontrivial
                   logical (stabilizer transport: p(stab) is a stab), so
                   the rung is a PURE FEASIBILITY statement — PASS iff no
                   cover cycle over w has overflow < M.  No sector
                   dispatch; every enumerated solution is a violation
                   (asserted non-stab as a convention tripwire).

Completeness of the restricted lane is the exact-off-support subset-sum
argument (A32 SS4): a solution v0 with off-support X0, |X0| <= M-1, forces
XOR_{j in X0} red[j] = rhs_res over the on-support-reduced columns; every
such X is found (sizes 0..6 via MITM joins) and the full solution kernel
over supp(b) u X enumerated and re-verified.

Engine validations (run by main()):
  V1  sector-scan linear-reduction trick == direct reduction (256/256)
  V2  planted dangerous control: v = tau2(u10) + HXc[g] is a genuine
      nontrivial dangerous logical over the weight-6 stab shadow; the
      rung at M = ov(v)+1 must FIND a violation (exercises the BZ lane),
      and every violation found must weigh >= 20 (else d < 20!)
  V3  banked-convention gate: census rows load as lab Y4-stabilizers
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

from a30_coset_bz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, run_window, unpack3,
)
from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a33_tower_cells import build_tower  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a33"


class YRungCell:
    """One deck cell base --y--> cover with the rung machinery."""

    def __init__(self, name: str, base: TS.BBCode, cover: TS.BBCode,
                 deck: TS.Deck):
        assert deck.base is base and deck.cover is cover
        self.name, self.base, self.cover, self.deck = name, base, cover, deck
        self.n, self.nc = base.n, cover.n
        self.E = deck.E                    # (cover.ng x base.n)
        self.RHS_OP = deck.RHS
        self.S_basis, self.S_piv = base.rsHX_b, base.rsHX_p
        self.Sc_basis, self.Sc_piv = cover.rsHX_b, cover.rsHX_p
        self.E_rows = [v2i(self.E[i]) for i in range(self.E.shape[0])]
        self.E_cols = [v2i(self.E[:, j]) for j in range(self.n)]
        # sectors = ker E / im S  (ker E = ker HZ(base): tau z cycle <=> z)
        kerE = [v2i(kv) for kv in base.kerHZ]
        sec, aug_b, aug_p = [], list(self.S_basis), list(self.S_piv)
        for kv in kerE:
            x = reduce_int(kv, aug_b, aug_p)
            if x:
                aug_b.append(x)
                aug_p.append((x & -x).bit_length() - 1)
                sec.append(kv)
        self.sector_basis = sec
        M01 = (deck.EMB[0] + deck.EMB[1]) % 2
        self.sector_chain_shift = [
            v2i((M01 @ i2v(t, self.n)) % 2) for t in self.sector_basis]
        # linear-reduction trick: reduce_int is linear; precompute the
        # reductions of the sector chain shifts (validated in main(), V1)
        self.red_shift = [reduce_int(s, self.Sc_basis, self.Sc_piv)
                          for s in self.sector_chain_shift]
        (self.I1, self.G1, self.I2, self.G2, self.kappa) = \
            disjoint_info_sets(base.HX)
        self.binp = build_kernel()

    # ------------------------------------------------------------ solving
    def solve_E(self, rhs: np.ndarray):
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
        ch = (self.deck.EMB[0] @ v0
              + self.deck.EMB[1] @ ((v0 + b_vec) % 2)) % 2
        return v2i(ch)

    # -------------------------------------------- restricted-lane MITM core
    def _hits_X(self, b_vec: np.ndarray, rhs: np.ndarray, cap: int):
        """All off-support sets X, |X| <= cap <= 6, with
        XOR red[j] = rhs_res (the exact-subset-sum lane)."""
        assert cap <= 6
        n = self.n
        bmask = v2i(b_vec)
        bcols = [self.E_cols[j] for j in np.nonzero(b_vec)[0]]
        bb, bp = rref_ints(bcols)
        rhs_res = reduce_int(v2i(rhs), bb, bp)
        offb = [j for j in range(n) if not (bmask >> j) & 1]
        red = {j: reduce_int(self.E_cols[j], bb, bp) for j in offb}
        by_val: dict[int, list[int]] = {}
        for j in offb:
            by_val.setdefault(red[j], []).append(j)
        hits: set[tuple[int, ...]] = set()
        if rhs_res == 0:
            hits.add(())
        if cap >= 1:
            for j in by_val.get(rhs_res, []):
                hits.add((j,))
        if cap >= 2:
            for j1 in offb:
                for j2 in by_val.get(rhs_res ^ red[j1], []):
                    if j2 > j1:
                        hits.add((j1, j2))
        if cap >= 3:
            for j1, j2 in itertools.combinations(offb, 2):
                for j3 in by_val.get(rhs_res ^ red[j1] ^ red[j2], []):
                    if j3 > j2:
                        hits.add((j1, j2, j3))
        pair_sum: dict[int, list[tuple[int, int]]] = {}
        if cap >= 4:
            for j1, j2 in itertools.combinations(offb, 2):
                pair_sum.setdefault(red[j1] ^ red[j2], []).append((j1, j2))
            for val, prs in pair_sum.items():
                for j3, j4 in pair_sum.get(rhs_res ^ val, []):
                    for j1, j2 in prs:
                        if j2 < j3:
                            hits.add((j1, j2, j3, j4))
        if cap >= 5:
            for j1, j2, j3 in itertools.combinations(offb, 3):
                tgt = rhs_res ^ red[j1] ^ red[j2] ^ red[j3]
                for j4, j5 in pair_sum.get(tgt, []):
                    if j4 > j3:
                        hits.add((j1, j2, j3, j4, j5))
        if cap >= 6:
            tri_sum: dict[int, list[tuple[int, int, int]]] = {}
            for tri in itertools.combinations(offb, 3):
                tri_sum.setdefault(
                    red[tri[0]] ^ red[tri[1]] ^ red[tri[2]], []).append(tri)
            for t1 in itertools.combinations(offb, 3):
                tgt = rhs_res ^ red[t1[0]] ^ red[t1[1]] ^ red[t1[2]]
                for t2 in tri_sum.get(tgt, []):
                    if t2[0] > t1[2]:
                        hits.add(t1 + t2)
        return hits

    def _expand_X(self, b_vec: np.ndarray, X: tuple[int, ...],
                  rhs_i: int, kernel_cap: int = 16):
        """All solutions v0 with supp(v0) <= supp(b) u X (may be none)."""
        cols = [int(j) for j in np.nonzero(b_vec)[0]] + list(X)
        b3: list[int] = []
        p3: list[int] = []
        h3: list[int] = []
        deps: list[int] = []
        for ci, j in enumerate(cols):
            cur, h = self.E_cols[j], 1 << ci
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
        cur, hsel = rhs_i, 0
        for bb3, pp3, hh in zip(b3, p3, h3):
            if (cur >> pp3) & 1:
                cur ^= bb3
                hsel ^= hh
        if cur:
            return
        assert len(deps) <= kernel_cap, f"kernel 2^{len(deps)} at X={X}"
        for kt in range(1 << len(deps)):
            sel = hsel
            for jj in range(len(deps)):
                if (kt >> jj) & 1:
                    sel ^= deps[jj]
            v0_int = 0
            for ci, j in enumerate(cols):
                if (sel >> ci) & 1:
                    v0_int |= 1 << int(j)
            yield v0_int

    # ------------------------------------------------------ dangerous rung
    def rung(self, b_vec: np.ndarray, M: int, deadline: float,
             validate_sectors: bool = False) -> dict:
        n = self.n
        wb = int(b_vec.sum())
        rhs = (self.RHS_OP @ b_vec) % 2
        v0p = self.solve_E(rhs)
        if v0p is None:
            return {"verdict": "PASS", "lane": "vacuous", "w_b": wb, "M": M}
        v0p_i = v2i(v0p)
        base_chain = self.chain_int(v0p_i, b_vec)
        red_base = reduce_int(base_chain, self.Sc_basis, self.Sc_piv)
        r = len(self.sector_basis)
        nontriv = []
        for t in range(1 << r):
            acc = red_base
            ti = 0
            for j in range(r):
                if (t >> j) & 1:
                    acc ^= self.red_shift[j]
                    ti ^= self.sector_basis[j]
            if validate_sectors:  # V1: the linear trick == direct reduce
                ch = base_chain
                for j in range(r):
                    if (t >> j) & 1:
                        ch ^= self.sector_chain_shift[j]
                assert (reduce_int(ch, self.Sc_basis, self.Sc_piv) != 0) \
                    == (acc != 0), "linear sector scan mismatch"
            if acc != 0:
                nontriv.append(v0p_i ^ ti)
        if not nontriv:
            return {"verdict": "PASS", "lane": "all-trivial", "w_b": wb,
                    "M": M}
        bmask = v2i(b_vec)
        rhs_i = v2i(rhs)
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
            for X in sorted(self._hits_X(b_vec, rhs, M - 1)):
                for v0_int in self._expand_X(b_vec, X, rhs_i):
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
                res = run_window(self.binp, f"a33_{self.name}_rung_w{wi}",
                                 Gs, bases, rr, Wp, deadline)
                for j, hx in res.pop("hit_rows"):
                    check(v2i(unpack3(hx, n)))
        if viols:
            return {"verdict": "VIOLATION", "lane": lane, "w_b": wb, "M": M,
                    "violations": viols[:5], "n_viol": len(viols)}
        return {"verdict": "PASS", "lane": lane, "w_b": wb, "M": M,
                "sectors_nontrivial": len(nontriv)}

    # ----------------------------------------------------------- seam rung
    def seam_rung(self, w_vec: np.ndarray, M: int) -> dict:
        """Feasibility rung over a nonzero-class coset element w."""
        assert M - 1 <= 6, "seam rung implemented for restricted lane only"
        n = self.n
        ww = int(w_vec.sum())
        rhs = (self.RHS_OP @ w_vec) % 2
        v0p = self.solve_E(rhs)
        if v0p is None:
            return {"verdict": "PASS", "lane": "vacuous", "w_w": ww, "M": M}
        wmask = v2i(w_vec)
        rhs_i = v2i(rhs)
        viols = []
        for X in sorted(self._hits_X(w_vec, rhs, M - 1)):
            for v0_int in self._expand_X(w_vec, X, rhs_i):
                ov = bin(v0_int & ~wmask).count("1")
                if ov > M - 1:
                    continue
                v0v = i2v(v0_int, n)
                assert not ((self.E @ v0v + rhs) % 2).any()
                ch = self.chain_int(v0_int, w_vec)
                # stabilizer transport: w non-stab => chain non-stab
                assert reduce_int(ch, self.Sc_basis, self.Sc_piv) != 0, \
                    "cycle over non-stab w reduced to a cover stab?!"
                wt = bin(ch).count("1")
                assert wt == ww + 2 * ov
                viols.append({"overflow": ov, "weight": wt,
                              "v0_hex": f"{v0_int:x}"})
        if viols:
            return {"verdict": "VIOLATION", "lane": f"restricted<={M-1}",
                    "w_w": ww, "M": M, "violations": viols[:5],
                    "n_viol": len(viols)}
        return {"verdict": "PASS", "lane": f"restricted<={M-1}", "w_w": ww,
                "M": M}


def main():
    t0 = time.monotonic()
    out: dict = {}
    Y2, Y4, Y8, deck_top, deck_bot = build_tower()
    cell_top = YRungCell("top", Y4, Y8, deck_top)
    cell_bot = YRungCell("bot", Y2, Y4, deck_bot)
    print(f"[{time.monotonic()-t0:5.1f}s] cells built: top sectors=2^"
          f"{len(cell_top.sector_basis)} kappa={cell_top.kappa}; bot "
          f"sectors=2^{len(cell_bot.sector_basis)} kappa={cell_bot.kappa}")
    out["top"] = {"sector_dim": len(cell_top.sector_basis),
                  "kappa": cell_top.kappa}
    out["bot"] = {"sector_dim": len(cell_bot.sector_basis),
                  "kappa": cell_bot.kappa}

    # V3: banked-convention gate
    census = []
    for line in (MAIN / "a20" / "m_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            census.append(r)
    assert len(census) == 1655
    for e in census[:25] + census[-25:]:
        b = np.zeros(Y4.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        assert Y4.is_stab(b), "banked census b NOT a lab Y4-stabilizer"
    print(f"[{time.monotonic()-t0:5.1f}s] V3 convention gate: banked census "
          f"rows are lab Y4-stabilizers (50 spot-checked; full load in "
          f"a33_validate_banked)")

    # V1: sector-scan linear trick, validated on the hexagon rung
    hexb = next(e for e in census if e["w"] == 6)
    bhex = np.zeros(Y4.n, dtype=np.uint8)
    bhex[hexb["b_support"]] = 1
    tV = time.monotonic()
    res = cell_top.rung(bhex, 7, time.monotonic() + 600,
                        validate_sectors=True)
    print(f"[{time.monotonic()-t0:5.1f}s] V1 sector linear trick 256/256 + "
          f"hexagon rung at M=7: {res['verdict']} (lane {res['lane']}, "
          f"{time.monotonic()-tV:.1f}s)")
    assert res["verdict"] == "PASS"
    out["hexagon_M7"] = res

    # V2: planted dangerous control (exercises the BZ lane)
    wit = None
    for line in (MAIN / "a20" / "y144_ladder.log").read_text().splitlines():
        try:
            wit = json.loads(line)["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    u10 = np.zeros(Y4.n, dtype=np.uint8)
    u10[wit] = 1
    assert Y4.is_cycle(u10) and not Y4.is_stab(u10) and u10.sum() == 10
    v = ((deck_top.TAU @ u10) + Y8.HX[0]) % 2
    assert Y8.is_cycle(v) and not Y8.is_stab(v)
    b6 = (deck_top.P @ v) % 2
    assert Y4.is_stab(b6) and int(b6.sum()) == 6
    _, ov, _ = deck_top.slice_data(v)
    wv = int(v.sum())
    assert wv == 6 + 2 * ov and wv >= 20, f"planted |v| = {wv} < 20 ?!"
    print(f"[{time.monotonic()-t0:5.1f}s] V2 planted: v = tau(u10)+HXc[0], "
          f"|v| = {wv}, shadow w6 stab, overflow {ov}; rung at M = {ov+1} "
          f"must FIND it ...")
    tV = time.monotonic()
    res = cell_top.rung(b6, ov + 1, time.monotonic() + 900)
    dt = time.monotonic() - tV
    assert res["verdict"] == "VIOLATION", f"planted control NOT found: {res}"
    found = min(x["overflow"] for x in res["violations"])
    assert found <= ov
    assert all(x["weight"] >= 20 for x in res["violations"]), \
        "sub-20 violation in control — d < 20?!"
    print(f"[{time.monotonic()-t0:5.1f}s] V2 planted control FOUND "
          f"(lane {res['lane']}, min overflow {found}, all weights >= 20, "
          f"{dt:.1f}s)  [BZ lane exercised]")
    out["planted"] = {"pv_weight": wv, "overflow": ov,
                      "found_min_overflow": found, "lane": res["lane"],
                      "n_viol": res["n_viol"], "secs": round(dt, 1)}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "rung_validation.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'rung_validation.json'}")


if __name__ == "__main__":
    main()
