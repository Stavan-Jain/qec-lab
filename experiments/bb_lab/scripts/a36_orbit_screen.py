"""A36 — construction screen: engineer literal-lift doubling presentations.

For a base code (orders, A, B, d) and a doubling axis, sweep the
seam-relevant presentation orbit (A14 §15 equivalence set v1):

  - doubled-axis translations of A and B INDEPENDENTLY (undoubled-axis
    translations are exact cover symmetries — they cannot move seams),
  - diagonal unit automorphisms (x, y) -> (x^u, y^v),
  - A <-> B swap,

and screen every cell for the safe floor (the discriminating template
condition — A11 S2, A35 refutation anatomy) in tiers whose kills are all
sound (each tier exhibits a genuine coset element below the floor 2d):

  T0   S0 raw seam minimum over ALL 2^kappa - 1 kernel classes
       (kappa = k/2); reject iff min < 2d.  Vectorised numpy.
  T0.5 lazy k-gate on S0 survivors: k(cover) == k(base)  <=> (R)  (A12).
  T1   exact per-class coset decision at the floor: ONE SAT call per
       class ("exists v in seam-coset with |v| <= 2d-1-parity?"),
       CaDiCaL; reject on any SAT (witness weight recorded).

Survivors are certify() candidates (run separately — see
a36_certify_runner.py; T2 budget rule lives there).

Conventions are pinned by an in-run assertion against
bb_lab.fibering.seam_offsets on the identity variant of every (code,
axis) — if the grid seam math ever disagrees with the front-end's, the
screen aborts rather than screening the wrong space.

Usage (from experiments/bb_lab/):
    uv run python scripts/a36_orbit_screen.py controls
    uv run python scripts/a36_orbit_screen.py sweep --point T1 --axis x \
        [--pres 0] [--t1-cap 40] [--out data/a36/T1_x.json]
    uv run python scripts/a36_orbit_screen.py cell --orders 7,9 \
        --A "..." --B "..." --axis x --d 10        # screen one cell
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter
from math import gcd
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB_ROOT / "src"))

from pysat.card import CardEnc, EncType          # noqa: E402
from pysat.formula import CNF, IDPool            # noqa: E402
from pysat.solvers import Cadical195             # noqa: E402

from bb_lab.checks import circulant              # noqa: E402
from bb_lab.fibering import (                    # noqa: E402
    kernel_basis, kernel_orbit_reps, seam_offsets,
)
from bb_lab.group import AbelianGroup            # noqa: E402
from bb_lab.linalg import nullspace_f2           # noqa: E402
from bb_lab.poly import Poly                     # noqa: E402
from bb_lab.sat_distance import _xor_chain       # noqa: E402

DATA_DIR = LAB_ROOT / "data" / "a36"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ points
# Stored presentations from the A35 docket (data/a35/docket.jsonl of the
# sibling session, mined 2026-08-11 from the corpus store + merit CSVs).
POINTS: dict[str, dict] = {
    "T1": {  # [[126,12,10]] -> [[252,12,20]] q=19.05
        "orders": (7, 9), "d": 10, "k": 12,
        "pres": [
            ("x + x*y^4 + x^2*y^2", "1 + x*y^6 + x^6*y^4"),
            ("x*y^3 + x^3*y^2 + x^3*y^4", "x^2*y^2 + x^3*y^4 + x^6*y^8"),
            ("x*y^5 + x^4*y^7 + x^6*y", "y + x + x*y^2"),
        ],
    },
    "P72": {  # pair72 base [[36,4,4]] — Lean-proven doubling control
        "orders": (3, 6), "d": 4, "k": 4,
        "pres": [("x^2 + y + y^3", "1 + x + y^2")],
    },
    "T1s": {  # [[126,12,10]] extra store rows (maxsat-tandem exact d),
              # mined 2026-08-11 from bb_instances.duckdb; same point as T1
        "orders": (7, 9), "d": 10, "k": 12,
        "pres": [
            ("y^4 + y^5 + x", "x*y^3 + x^3*y + x^5*y^4"),
            ("y^4 + y^5 + x", "x*y^2 + x^3*y^4 + x^5*y"),
            ("y^4 + y^5 + x", "y^8 + x + x^2*y^3"),
            ("y^4 + y^5 + x", "1 + x^3*y + x^5*y^3"),
            ("y^4 + y^5 + x", "x^2 + x^4*y^7 + x^6*y"),
            ("y^4 + y^5 + x", "y^7 + x^3*y^8 + x^5*y"),
        ],
    },
    "T2": {  # [[120,8,12]] -> [[240,8,24]] q=19.20
        "orders": (6, 10), "d": 12, "k": 8,
        "pres": [
            ("x*y^3 + x^3*y^2 + x^4*y", "1 + x^5*y^2 + x^5*y^3"),
            ("x*y^4 + x^3*y^5 + x^4*y^6", "x*y^7 + x*y^8 + x^2"),
            ("x*y^6 + x^2*y^8 + x^2*y^9", "y^6 + x^3*y^4 + x^4*y^5"),
        ],
    },
    "T2b": {  # [[120,8,12]] -> [[240,8,24]] q=19.20
        "orders": (5, 12), "d": 12, "k": 8,
        "pres": [
            ("1 + x*y + x^2", "y^7 + x^2*y^11 + x^3*y^8"),
            ("x + x^2*y^11 + x^3*y^3", "y^5 + x*y^2 + x^3*y^3"),
            ("x*y + x^3*y^2 + x^4*y^11", "y^4 + x*y^3 + x^4"),
        ],
    },
    "T3": {  # [[90,8,10]] -> [[180,8,20]] q=17.78
        "orders": (15, 3), "d": 10, "k": 8,
        "pres": [
            ("y + y^2 + x^3", "1 + x + x^11"),
            ("y + y^2 + x^3", "1 + x + x^5"),
            ("y + y^2 + x^3", "1 + x^4 + x^5"),
        ],
    },
    "T4": {  # [[98,6,12]] -> [[196,6,24]] q=17.63
        "orders": (7, 7), "d": 12, "k": 6,
        "pres": [
            ("1 + y + x", "x*y^2 + x^4*y^5 + x^5*y"),
            ("1 + y + x", "y + x*y^4 + x^4"),
            ("1 + y + x", "y^2 + x*y^5 + x^4*y"),
        ],
    },
    "T5": {  # [[112,6,12]] -> [[224,6,24]] q=15.43
        "orders": (7, 8), "d": 12, "k": 6,
        "pres": [
            ("y + x + x^3", "1 + x*y^2 + x^3*y"),
            ("y + x + x^3", "1 + x^2*y^7 + x^6*y^4"),
            ("y + x + x^3", "x + x^2*y^6 + x^4*y^5"),
        ],
    },
    "T6": {  # [[84,6,10]] -> [[168,6,20]] q=14.29
        "orders": (6, 7), "d": 10, "k": 6,
        "pres": [
            ("y + y^3 + x", "1 + x*y + x^4*y^3"),
            ("y + y^3 + x", "1 + x*y^3 + x^2*y"),
            ("y + y^3 + x", "1 + x*y^5 + x^5*y^4"),
        ],
    },
}


# ------------------------------------------------------------------ helpers
def parse_poly(s: str) -> list[tuple[int, int]]:
    support: list[tuple[int, int]] = []
    for term in s.replace("−", "-").split("+"):
        term = term.strip()
        ex = ey = 0
        if term != "1":
            for factor in term.split("*"):
                factor = factor.strip()
                var, _, exp = factor.partition("^")
                e = int(exp) if exp else 1
                if var == "x":
                    ex += e
                elif var == "y":
                    ey += e
                else:
                    raise ValueError(f"bad factor {factor!r} in {s!r}")
        support.append((ex, ey))
    return support


def poly_str(sup) -> str:
    terms = []
    for (e, f) in sorted(sup):
        t = []
        if e:
            t.append(f"x^{e}" if e > 1 else "x")
        if f:
            t.append(f"y^{f}" if f > 1 else "y")
        terms.append("*".join(t) if t else "1")
    return " + ".join(terms)


def units(n: int) -> list[int]:
    return [u for u in range(1, n) if gcd(u, n) == 1]


# --------------------------------------------------------------- the screen
class AxisProblem:
    """One (code, axis), canonicalised so the doubled axis is axis 0.

    All screen math runs in the canonical frame; `to_original` maps a
    canonical support back for cover-spec emission.
    """

    def __init__(self, orders: tuple[int, int], A_sup, B_sup, axis: int,
                 d: int):
        self.axis_orig = axis
        self.orders_orig = tuple(orders)
        if axis == 1:  # swap coordinates so we always double axis 0
            orders = (orders[1], orders[0])
            A_sup = [(f, e) for (e, f) in A_sup]
            B_sup = [(f, e) for (e, f) in B_sup]
        self.ell, self.m = orders
        self.d = d
        self.target = 2 * d
        self.A0 = [(e % self.ell, f % self.m) for (e, f) in A_sup]
        self.B0 = [(e % self.ell, f % self.m) for (e, f) in B_sup]
        self.G = AbelianGroup((self.ell, self.m))
        self.Gc = AbelianGroup((2 * self.ell, self.m))
        self._seam_convention_checked = False

    # -- variant machinery ------------------------------------------------
    def variant_supports(self, u: int, v: int, swap: bool):
        Au = [((u * e) % self.ell, (v * f) % self.m) for (e, f) in self.A0]
        Bu = [((u * e) % self.ell, (v * f) % self.m) for (e, f) in self.B0]
        return (Bu, Au) if swap else (Au, Bu)

    def kernel_grids(self, Asup, Bsup) -> np.ndarray:
        """All 2^kappa - 1 nonzero kernel elements as (N, ell, m) grids."""
        A = Poly.from_support([tuple(t) for t in Asup], self.G)
        B = Poly.from_support([tuple(t) for t in Bsup], self.G)
        K = kernel_basis(A, B)  # (kappa, n)
        kappa = K.shape[0]
        combos = []
        for mask in range(1, 1 << kappa):
            z = np.zeros(K.shape[1], dtype=np.uint8)
            for j in range(kappa):
                if (mask >> j) & 1:
                    z ^= K[j]
            combos.append(z)
        return np.array(combos, dtype=np.uint8).reshape(-1, self.ell, self.m)

    def seam_grids(self, Psup, zgrids: np.ndarray) -> np.ndarray:
        """Carry grids (N, ell, m) of P~ . zeta_i~ (this cell's polys)."""
        N = zgrids.shape[0]
        ell, m = self.ell, self.m
        lifts = np.zeros((N, 2 * ell, m), dtype=np.uint8)
        lifts[:, :ell, :] = zgrids  # canonical section: coords < ell
        F = np.zeros_like(lifts)
        for (e, f) in Psup:
            F ^= np.roll(np.roll(lifts, e, axis=1), f, axis=2)
        return F[:, ell:, :]

    def seam_weights(self, Psup, zgrids: np.ndarray) -> np.ndarray:
        """W[i, t] = |carry((x^t P)~ . zeta_i~)| for every translation t.

        The lift of x^t P reads exponents REDUCED mod ell — wrapped
        monomials pick up a deck factor delta = x^ell, which is exactly
        what moves the seams.  Computed by direct rolls per t.
        """
        N = zgrids.shape[0]
        ell, m = self.ell, self.m
        lifts = np.zeros((N, 2 * ell, m), dtype=np.uint8)
        lifts[:, :ell, :] = zgrids  # canonical section: coords < ell
        W = np.zeros((N, ell), dtype=np.int32)
        for t in range(ell):
            F = np.zeros_like(lifts)
            for (e, f) in Psup:
                er = (e + t) % ell  # literal lift of the translated poly
                F ^= np.roll(np.roll(lifts, er, axis=1), f, axis=2)
            W[:, t] = F[:, ell:, :].reshape(N, -1).sum(axis=1)
        return W

    def s0_plane(self, u: int, v: int, swap: bool):
        """S0[ta, tb] = min over kernel classes of |seam_A| + |seam_B|."""
        Asup, Bsup = self.variant_supports(u, v, swap)
        zg = self.kernel_grids(Asup, Bsup)
        WA = self.seam_weights(Asup, zg)          # (N, ell)
        WB = self.seam_weights(Bsup, zg)          # (N, ell)
        s0 = (WA[:, :, None] + WB[:, None, :]).min(axis=0)  # (ell, ell)
        return s0, (Asup, Bsup)

    # -- per-cell exact machinery ------------------------------------------
    def cell_polys(self, Asup, Bsup, ta: int, tb: int):
        A2 = sorted(((e + ta) % self.ell, f) for (e, f) in Asup)
        B2 = sorted(((e + tb) % self.ell, f) for (e, f) in Bsup)
        return A2, B2

    def check_seam_convention(self, Asup, Bsup) -> None:
        """Assert grid seams == bb_lab.fibering.seam_offsets (once)."""
        if self._seam_convention_checked:
            return
        A = Poly.from_support([tuple(t) for t in Asup], self.G)
        B = Poly.from_support([tuple(t) for t in Bsup], self.G)
        offs = seam_offsets(A, B, 0)
        zg = self.kernel_grids(Asup, Bsup)
        flat = zg.reshape(zg.shape[0], -1)
        SU = self.seam_grids(Asup, zg).reshape(zg.shape[0], -1)
        SV = self.seam_grids(Bsup, zg).reshape(zg.shape[0], -1)
        for (zeta, su, sv) in offs:
            hits = np.where((flat == zeta[None, :]).all(axis=1))[0]
            assert hits.size == 1, "orbit rep not found among kernel combos"
            i = int(hits[0])
            assert (SU[i] == su).all() and (SV[i] == sv).all(), (
                "seam VECTOR mismatch vs bb_lab.fibering.seam_offsets — "
                "grid convention broken, aborting")
        self._seam_convention_checked = True

    def k_gate(self, A2, B2) -> tuple[int, int]:
        A = Poly.from_support([tuple(t) for t in A2], self.G)
        B = Poly.from_support([tuple(t) for t in B2], self.G)
        Ac = Poly.from_support([tuple(t) for t in A2], self.Gc)
        Bc = Poly.from_support([tuple(t) for t in B2], self.Gc)
        kb = 2 * kernel_basis(A, B).shape[0]
        kc = 2 * kernel_basis(Ac, Bc).shape[0]
        return kb, kc

    def t1_cell(self, A2, B2, time_budget: float = 90.0,
                conf_budget: int = 50_000) -> dict:
        """Conflict-budgeted safe-floor KILL screen at the floor 2d.

        For every nonzero kernel class: one budgeted SAT call 'exists v
        in seam-coset with |v| < 2d (parity-adjusted)?'.
          SAT   -> the cell is killed (a genuine light coset element —
                   sound, certificate-grade negative);
          UNSAT -> the class is proven heavy;
          UNDET (conflict/time budget) -> unknown, recorded — the cell
                   still 'passes' the screen and certify()'s BZ
                   safe-floor stage is the real decider.
        The kill direction is the cheap one (find a light vector); the
        UNSAT direction can be arbitrarily hard at these sizes, which is
        exactly why this is a screen, not a certifier.
        """
        t0 = time.time()
        A = Poly.from_support([tuple(t) for t in A2], self.G)
        B = Poly.from_support([tuple(t) for t in B2], self.G)
        MA = circulant(A).astype(np.uint8) % 2
        MB = circulant(B).astype(np.uint8) % 2
        MS = np.concatenate([MA.T, MB.T], axis=1) % 2  # code_rows
        dual = nullspace_f2(MS)
        # G-orbit reps of the kernel classes only (Prop A14.1 transport —
        # exactly what the front-end's safe_floor probes); cell polys
        # already carry their translations, so seams are read at t = 0.
        # Per class, probe several TRANSLATED members: the instances are
        # isomorphic (translation = coordinate permutation), so one
        # SAT kills and one UNSAT decides the class; extra members are
        # pure solver diversification against UNDET false-passes.
        zreps = kernel_orbit_reps(A, B)
        zg = np.array(zreps, dtype=np.uint8).reshape(
            -1, self.ell, self.m)
        N = zg.shape[0]
        SU = self.seam_grids(A2, zg).reshape(N, -1)
        SV = self.seam_grids(B2, zg).reshape(N, -1)
        raw = np.concatenate([SU, SV], axis=1)
        order = np.argsort(raw.sum(axis=1))  # lightest raw seams first
        tries = 6 if N == 1 else (3 if N <= 3 else 2)
        shifts = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 1), (1, 2),
                  (3, 2), (2, 3)]
        n_unsat = n_undet = 0
        for i in order:
            # cheap freeze probe at cap d (the A35 anatomy: most cells
            # die AT d_base)
            rep = raw[i].astype(np.uint8)
            verdict, w = self._coset_light(rep, dual, self.d,
                                           conf_budget)
            if verdict == "SAT":
                return {"pass": False, "kill_weight": int(w),
                        "kill_class": int(i), "lane": "freeze-probe",
                        "classes_done": n_unsat + n_undet + 1,
                        "wall_s": round(time.time() - t0, 2)}
            # decisive probes at the floor over translated members
            decided = None
            for gx, gy in shifts[:tries]:
                zt = np.roll(np.roll(zg[i], gx, axis=0), gy, axis=1)
                zt3 = zt[None, :, :]
                su = self.seam_grids(A2, zt3).reshape(-1)
                sv = self.seam_grids(B2, zt3).reshape(-1)
                rep_t = np.concatenate([su, sv]).astype(np.uint8)
                verdict, w = self._coset_light(
                    rep_t, dual, self.target - 1, conf_budget)
                if verdict == "SAT":
                    return {"pass": False, "kill_weight": int(w),
                            "kill_class": int(i), "lane": "floor-probe",
                            "classes_done": n_unsat + n_undet + 1,
                            "wall_s": round(time.time() - t0, 2)}
                if verdict == "UNSAT":
                    decided = "UNSAT"
                    break
                if time.time() - t0 > time_budget:
                    break
            if decided == "UNSAT":
                n_unsat += 1
            else:
                n_undet += 1
            if time.time() - t0 > time_budget:
                n_undet += N - (n_unsat + n_undet)
                break
        return {"pass": True, "n_classes": N, "n_unsat": n_unsat,
                "n_undet": n_undet,
                "wall_s": round(time.time() - t0, 2)}

    @staticmethod
    def _coset_light(rep: np.ndarray, dual: np.ndarray, cap: int,
                     conf_budget: int) -> tuple[str, int | None]:
        """Budgeted: is there v in rep + rowspace with |v| <= cap?

        Returns ("SAT", weight_found) / ("UNSAT", None) /
        ("UNDET", None).  The weight reported is the first witness's
        (an upper bound on the coset min — never report it as a
        minimum; A17 §8's lesson).
        """
        n = rep.shape[0]
        parity = int(rep.sum() % 2)
        bound = cap if (cap % 2) == parity else cap - 1
        if bound < parity:
            return ("UNSAT", None)
        pool = IDPool()
        qv = [pool.id() for _ in range(n)]
        cnf = CNF()
        rhs = (dual @ rep) % 2
        for row, r in zip(dual, rhs):
            idxs = np.flatnonzero(row)
            out = _xor_chain((qv[i] for i in idxs), pool, cnf)
            if out is None:
                continue
            cnf.append([out] if r else [-out])
        if bound < n:
            cnf.extend(CardEnc.atmost(
                lits=qv, bound=bound, vpool=pool,
                encoding=EncType.seqcounter).clauses)
        with Cadical195(bootstrap_with=cnf.clauses) as solver:
            solver.conf_budget(conf_budget)
            res = solver.solve_limited()
            if res is None:
                return ("UNDET", None)
            if res is False:
                return ("UNSAT", None)
            model = solver.get_model()
            truth = {abs(li): li > 0 for li in model}
            w = sum(1 for vv in qv if truth.get(vv, False))
        return ("SAT", w)

    # -- emission ----------------------------------------------------------
    def cover_spec(self, A2, B2) -> dict:
        """The certify() input for this cell (original orientation)."""
        if self.axis_orig == 1:
            orders = (self.orders_orig[0], 2 * self.orders_orig[1])
            A_o = [(f, e) for (e, f) in A2]
            B_o = [(f, e) for (e, f) in B2]
        else:
            orders = (2 * self.orders_orig[0], self.orders_orig[1])
            A_o, B_o = A2, B2
        return {"cover_orders": list(orders),
                "cover_A": poly_str(sorted(A_o)),
                "cover_B": poly_str(sorted(B_o))}


# ------------------------------------------------------------------- sweep
def sweep(point: str, axis: str, pres_idx: int | None, t1_cap: int,
          out: Path | None) -> None:
    spec = POINTS[point]
    orders, d = spec["orders"], spec["d"]
    pres_list = (spec["pres"] if pres_idx is None
                 else [spec["pres"][pres_idx]])
    t_start = time.time()
    all_records = []
    for pi, (a_s, b_s) in enumerate(pres_list):
        if pres_idx is not None:
            pi = pres_idx
        ap = AxisProblem(orders, parse_poly(a_s), parse_poly(b_s),
                         0 if axis == "x" else 1, d)
        us, vs = units(ap.ell), units(ap.m)
        print(f"== {point} pres {pi} axis {axis}: orbit "
              f"{len(us) * len(vs) * 2} variants x {ap.ell}^2 cells, "
              f"floor {ap.target} ==", flush=True)
        seen: set = set()
        s0_hist: Counter = Counter()
        survivors = []
        for u, v, swap in itertools.product(us, vs, (False, True)):
            s0, (Asup, Bsup) = ap.s0_plane(u, v, swap)
            if u == 1 and v == 1 and not swap:
                ap.check_seam_convention(Asup, Bsup)
            for ta in range(ap.ell):
                for tb in range(ap.ell):
                    A2, B2 = ap.cell_polys(Asup, Bsup, ta, tb)
                    key = (tuple(A2), tuple(B2))
                    if key in seen:
                        continue
                    seen.add(key)
                    val = int(s0[ta, tb])
                    s0_hist[val] += 1
                    if val >= ap.target:
                        survivors.append((val, A2, B2))
        print(f"  {len(seen)} unique cells; S0 histogram "
              f"{dict(sorted(s0_hist.items()))}; "
              f"S0 survivors (>= {ap.target}): {len(survivors)} "
              f"({time.time() - t_start:.0f}s)", flush=True)

        # T0.5 lazy k-gate + T1, strongest S0 first
        survivors.sort(key=lambda r: -r[0])
        kgate_fail = 0
        t1_run = 0
        finalists = []
        for val, A2, B2 in survivors:
            if t1_run >= t1_cap:
                break
            kb, kc = ap.k_gate(A2, B2)
            if kb != kc or kb != spec["k"]:
                kgate_fail += 1
                all_records.append({
                    "point": point, "pres": pi, "axis": axis,
                    "A": poly_str(A2), "B": poly_str(B2), "s0": val,
                    "stage": "k-fail", "k": kb, "k_cover": kc})
                continue
            t1_run += 1
            r = ap.t1_cell(A2, B2)
            rec = {"point": point, "pres": pi, "axis": axis,
                   "A": poly_str(A2), "B": poly_str(B2), "s0": val,
                   "stage": "t1", **r, **ap.cover_spec(A2, B2)}
            all_records.append(rec)
            status = ("SF-PASS" if r["pass"] else
                      f"kill@{r.get('kill_weight')}")
            print(f"  T1 [{t1_run}/{min(t1_cap, len(survivors))}] "
                  f"s0={val} {status} ({r['wall_s']}s) "
                  f"A={poly_str(A2)} B={poly_str(B2)}", flush=True)
            if r["pass"]:
                finalists.append(rec | ap.cover_spec(A2, B2))
                print(f"  *** SF-PASS finalist: "
                      f"{json.dumps(ap.cover_spec(A2, B2))}", flush=True)
        print(f"  pres {pi}: k-gate fails {kgate_fail}, T1 run {t1_run}, "
              f"finalists {len(finalists)}", flush=True)

    if out is None:
        out = DATA_DIR / f"{point}_{axis}_screen.json"
    out.write_text(json.dumps({
        "point": point, "axis": axis, "records": all_records,
        "wall_s": round(time.time() - t_start, 1)}, indent=1))
    print(f"DONE -> {out} ({time.time() - t_start:.0f}s)", flush=True)


# ---------------------------------------------------------------- controls
def controls() -> None:
    """Known-positive + known-negative + A35-consistency checks."""
    print("== control 1 (known-POSITIVE): A30-certified [[180,4,10]] "
          "Z15xZ6 x-axis ==", flush=True)
    ap = AxisProblem((15, 6), parse_poly("1 + y + x"),
                     parse_poly("y^4 + x + x^11*y^2"), 0, 10)
    ap.check_seam_convention(ap.A0, ap.B0)
    print("  seam convention vs bb_lab.fibering: OK", flush=True)
    zg = ap.kernel_grids(ap.A0, ap.B0)
    WA = ap.seam_weights(ap.A0, zg)[:, 0]
    WB = ap.seam_weights(ap.B0, zg)[:, 0]
    s0 = int((WA + WB).min())
    print(f"  S0 (stored cell) = {s0}  (need >= 20 to pass T0)", flush=True)
    r = ap.t1_cell(ap.A0, ap.B0)
    print(f"  T1: {r}  (must be pass=True)", flush=True)
    assert r["pass"], "known-positive control FAILED T1"

    print("== control 2 (known-NEGATIVE): bb_108 stored-y "
          "(A17 §8 exact d_safe = 14) ==", flush=True)
    ap2 = AxisProblem((9, 6), parse_poly("x^3 + y + y^2"),
                      parse_poly("y^3 + x + x^2"), 1, 10)
    ap2.check_seam_convention(
        *ap2.variant_supports(1, 1, False))
    zg2 = ap2.kernel_grids(ap2.A0, ap2.B0)
    WA2 = ap2.seam_weights(ap2.A0, zg2)[:, 0]
    WB2 = ap2.seam_weights(ap2.B0, zg2)[:, 0]
    print(f"  S0 (stored cell) = {int((WA2 + WB2).min())}", flush=True)
    r2 = ap2.t1_cell(ap2.A0, ap2.B0)
    print(f"  T1: {r2}  (must be pass=False; first-found witness weight "
          "is an upper bound)", flush=True)
    assert not r2["pass"], "known-negative control PASSED T1?!"
    # pin the exact value 14 (A17 §8): SAT at 14, UNSAT at 12
    A = Poly.from_support([tuple(t) for t in ap2.A0], ap2.G)
    B = Poly.from_support([tuple(t) for t in ap2.B0], ap2.G)
    MA = circulant(A).astype(np.uint8) % 2
    MB = circulant(B).astype(np.uint8) % 2
    dual2 = nullspace_f2(np.concatenate([MA.T, MB.T], axis=1) % 2)
    zo2 = np.array(kernel_orbit_reps(A, B), dtype=np.uint8).reshape(
        -1, ap2.ell, ap2.m)
    SU2 = ap2.seam_grids(ap2.A0, zo2).reshape(zo2.shape[0], -1)
    SV2 = ap2.seam_grids(ap2.B0, zo2).reshape(zo2.shape[0], -1)
    rep_k = np.concatenate([SU2[r2["kill_class"]],
                            SV2[r2["kill_class"]]]).astype(np.uint8)
    v14, _w14 = AxisProblem._coset_light(rep_k, dual2, 14, 500_000)
    v12, _w12 = AxisProblem._coset_light(rep_k, dual2, 12, 500_000)
    print(f"  exact pin on the killed class: <=14 {v14}, <=12 {v12} "
          "(A17 §8 says d_safe = 14: want SAT / UNSAT)", flush=True)
    assert v14 == "SAT", "bb108-y killed class not SAT at 14?!"
    assert v12 in ("UNSAT", "UNDET"), "bb108-y class SAT at 12?!"

    print("== control 3 (A35 consistency): [[126,12,10]] Z7xZ9 stored "
          "p0, x-axis (certify: DOUBLING-REFUTED) ==", flush=True)
    p = POINTS["T1"]["pres"][0]
    ap3 = AxisProblem((7, 9), parse_poly(p[0]), parse_poly(p[1]), 0, 10)
    ap3.check_seam_convention(ap3.A0, ap3.B0)
    zg3 = ap3.kernel_grids(ap3.A0, ap3.B0)
    WA3 = ap3.seam_weights(ap3.A0, zg3)[:, 0]
    WB3 = ap3.seam_weights(ap3.B0, zg3)[:, 0]
    s0 = int((WA3 + WB3).min())
    print(f"  S0 (stored cell) = {s0}", flush=True)
    if s0 >= 20:
        r3 = ap3.t1_cell(ap3.A0, ap3.B0)
        print(f"  T1: {r3}  (must be pass=False)", flush=True)
        assert not r3["pass"], "A35-refuted cell PASSED the screen?!"
    else:
        print("  killed at T0 (consistent with the certify refutation)",
              flush=True)
    print("ALL CONTROLS GREEN", flush=True)


# ------------------------------------------------------------------- cell
def one_cell(orders_s: str, a_s: str, b_s: str, axis: str, d: int) -> None:
    orders = tuple(int(t) for t in orders_s.split(","))
    ap = AxisProblem(orders, parse_poly(a_s), parse_poly(b_s),
                     0 if axis == "x" else 1, d)
    ap.check_seam_convention(*ap.variant_supports(1, 1, False))
    zg = ap.kernel_grids(ap.A0, ap.B0)
    WA = ap.seam_weights(ap.A0, zg)[:, 0]
    WB = ap.seam_weights(ap.B0, zg)[:, 0]
    print(f"S0 = {int((WA + WB).min())} (floor {ap.target})", flush=True)
    kb, kc = ap.k_gate(ap.A0, ap.B0)
    print(f"k-gate: k = {kb}, k_cover = {kc} "
          f"({'OK' if kb == kc else '(R) FAILS'})", flush=True)
    r = ap.t1_cell(ap.A0, ap.B0)
    print(f"T1: {r}", flush=True)
    if r["pass"]:
        print(f"cover spec: {json.dumps(ap.cover_spec(ap.A0, ap.B0))}",
              flush=True)


def main() -> None:
    apar = argparse.ArgumentParser()
    apar.add_argument("cmd", choices=["controls", "sweep", "cell"])
    apar.add_argument("--point", type=str)
    apar.add_argument("--axis", type=str, default="x")
    apar.add_argument("--pres", type=int, default=None)
    apar.add_argument("--t1-cap", type=int, default=40)
    apar.add_argument("--out", type=str, default=None)
    apar.add_argument("--orders", type=str)
    apar.add_argument("--A", type=str)
    apar.add_argument("--B", type=str)
    apar.add_argument("--d", type=int, default=10)
    args = apar.parse_args()
    if args.cmd == "controls":
        controls()
    elif args.cmd == "sweep":
        sweep(args.point, args.axis, args.pres, args.t1_cap,
              Path(args.out) if args.out else None)
    else:
        one_cell(args.orders, args.A, args.B, args.axis, args.d)


if __name__ == "__main__":
    main()
