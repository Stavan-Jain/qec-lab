#!/usr/bin/env python3
"""A42 S1 — exact windowed compact-phase machinery for the period-p
cylinder, and the falsify-first probes of the "floor >= 2p for all p"
claim (L-P's premise).

System convention (X-sector of the lab codes, TowerCode-compatible):
  cycles:      P f + Q g = 0   with (P, Q) = (Bbar, Abar)
  boundaries:  (f, g) = (R t, S t)  with (R, S) = (Abar, Bbar)
  (P R + Q S = 2 Abar Bbar = 0.)
where Abar = 1 + y^-1 + x^-3 y, Bbar = 1 + x^-1 + x y^3 (antipodes of
the fixed pair).  The A40 atlas floors are lane-symmetric, so the
banked calibration values (p=3: 6, p=6: 12; p in {2,4,5,7,8}: none)
apply to this sector verbatim.

Exactness of window triviality: for t with rightmost active column c,
(Bbar t) has rightmost column c+1 with content y^3 * t_c (a unit
multiple), and with leftmost active column c', (Abar t) has leftmost
column c'-3 with content y * t_{c'} — so a boundary supported in the
window forces t into a fixed enlargement of the window (no
telescoping), and window-restricted boundary space == cylinder
boundary space restricted to the window.  Hence the class functionals
built from the window boundary space are EXACT cylinder-triviality
functionals for window-supported cycles.

Probes:
  A. Calibration: p=3 min nontrivial = 6; p=5,7: Z == B (Theorem A
     machine check); p=6 min = 12.
  B. The (12,12) witness-gap probe: d((12,12)) = 18 is certified;
     hunt weight-18 X-logicals by CMS and test for cyclic x-gaps >= 4
     (union support of both blocks).  Any gapped witness windowizes
     (Lemma K move) into a weight-18 nontrivial compact period-12
     phase — refuting "floor >= 2p at p = 12" outright.
  C. Direct cylinder SAT at (p, W, wmax): min nontrivial compact
     weight with x-extent < W (exact within the window).

Every found object is re-verified end-to-end (cycle equations,
weight, functional pairing, embedding-torus nontriviality).
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

from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402
import pycryptosat  # noqa: E402

from bb_lab.linalg import nullspace_f2  # noqa: E402
from bb_lab.tower import TowerCode  # noqa: E402

DATA = LAB / "data" / "a42"
DATA.mkdir(parents=True, exist_ok=True)

# supports as (dx, dy): X-sector system
P_SUPP = ((0, 0), (-1, 0), (1, 3))    # Bbar
Q_SUPP = ((0, 0), (0, -1), (-3, 1))   # Abar
R_SUPP = Q_SUPP                       # boundary: f = Abar t
S_SUPP = P_SUPP                       # boundary: g = Bbar t

A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def rank_f2(M: np.ndarray) -> int:
    M = (M % 2).astype(np.uint8).copy()
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


class CylWindow:
    """Window [0, W) x Z_p of the cylinder, X-sector system."""

    def __init__(self, p: int, W: int):
        self.p, self.W = p, W
        self.n = 2 * W * p

    def vid(self, blk: int, c: int, j: int) -> int:
        return blk * self.W * self.p + c * self.p + (j % self.p)

    def build_H(self) -> np.ndarray:
        p, W = self.p, self.W
        rows = []
        for c in range(-4, W + 5):
            for j in range(p):
                row = np.zeros(self.n, dtype=np.uint8)
                any_on = False
                for (dx, dy) in P_SUPP:
                    cc = c - dx
                    if 0 <= cc < W:
                        row[self.vid(0, cc, j - dy)] ^= 1
                        any_on = True
                for (dx, dy) in Q_SUPP:
                    cc = c - dx
                    if 0 <= cc < W:
                        row[self.vid(1, cc, j - dy)] ^= 1
                        any_on = True
                if any_on and row.any():
                    rows.append(row)
        return np.array(rows, dtype=np.uint8)

    def build_boundary_in_window(self) -> np.ndarray:
        """Basis (rows) of the cylinder-boundary space restricted to
        window-supported vectors.  t ranges over a slack-enlarged
        window; ambient columns cover all image columns; rows of the
        image matrix are then combined so that out-of-window support
        vanishes (exact by the unit-edge argument)."""
        p, W = self.p, self.W
        SL = 5
        tcols = range(-SL, W + SL)
        # ambient x-range of images: t in [-SL, W+SL) shifted by dx in
        # [-3, 1] -> [-SL-3, W+SL+1); index shift
        lo = -SL - 3
        hi = W + SL + 1
        WA = hi - lo
        namb = 2 * WA * p

        def amb(blk, c, j):
            return blk * WA * p + (c - lo) * p + (j % p)

        rows = []
        for tc in tcols:
            for tj in range(p):
                row = np.zeros(namb, dtype=np.uint8)
                for (dx, dy) in R_SUPP:
                    row[amb(0, tc + dx, tj + dy)] ^= 1
                for (dx, dy) in S_SUPP:
                    row[amb(1, tc + dx, tj + dy)] ^= 1
                rows.append(row)
        D = np.array(rows, dtype=np.uint8)
        # split ambient columns into window vars (ordered like vid) and out
        win_idx = np.zeros(self.n, dtype=np.int64)
        for blk in (0, 1):
            for c in range(W):
                for j in range(p):
                    win_idx[self.vid(blk, c, j)] = amb(blk, c, j)
        mask = np.zeros(namb, dtype=bool)
        mask[win_idx] = True
        D_win = D[:, win_idx]
        D_out = D[:, ~mask]
        K = nullspace_f2(D_out.T)      # left-kernel of D_out
        if K.size == 0:
            return np.zeros((0, self.n), dtype=np.uint8)
        B = (K @ D_win) % 2
        B = B[B.any(axis=1)]
        return B.astype(np.uint8)

    def class_functionals(self, H: np.ndarray, Bnd: np.ndarray
                          ) -> tuple[np.ndarray, int, int]:
        """Rows L s.t. window cycle v is cylinder-nontrivial iff
        L v != 0.  Returns (L, dimZ, dimB)."""
        Z = nullspace_f2(H)
        dimZ = Z.shape[0]
        dimB = rank_f2(Bnd) if Bnd.size else 0
        if Bnd.size:
            N = nullspace_f2(Bnd)      # all functionals vanishing on B
        else:
            N = np.eye(self.n, dtype=np.uint8)
        if dimZ == dimB:
            return np.zeros((0, self.n), dtype=np.uint8), dimZ, dimB
        # pick N-rows spanning the dual of Z/B: pivoted reduce of N@Z^T
        Pmat = (N @ Z.T) % 2
        picked = []
        basis = []
        for i in range(Pmat.shape[0]):
            v = Pmat[i].copy()
            for (bv, bp) in basis:
                if v[bp]:
                    v ^= bv
            nz = np.flatnonzero(v)
            if nz.size:
                basis.append((v, nz[0]))
                picked.append(i)
            if len(picked) == dimZ - dimB:
                break
        L = N[picked]
        return L.astype(np.uint8), dimZ, dimB

    # ---- SAT ----
    def sat_min_nontrivial(self, H, L, wmax: int, log=print,
                           max_models_per_w: int = 1):
        """Exact min weight of a cylinder-nontrivial window cycle,
        searched downward from wmax.  Returns (min_w or None, witness)."""
        best = None
        best_v = None
        w = wmax
        while w >= 1:
            v = self._solve_at(H, L, w, log=log)
            if v is None:
                break
            wt = int(v.sum())
            assert wt <= w
            best, best_v = wt, v
            log(f"    SAT at <= {w}: found weight {wt}; descending")
            w = wt - 1
        return best, best_v

    def _solve_at(self, H, L, weight, log=print):
        pool = IDPool()
        qv = [pool.id() for _ in range(self.n)]
        solver = pycryptosat.Solver()
        for row in H:
            idx = np.flatnonzero(row)
            if idx.size:
                solver.add_xor_clause([qv[i] for i in idx], False)
        a_outs = []
        for row in L:
            idx = np.flatnonzero(row)
            if idx.size == 0:
                continue
            a = pool.id()
            solver.add_xor_clause([qv[i] for i in idx] + [a], False)
            a_outs.append(a)
        assert a_outs, "no functionals — nothing nontrivial to find"
        solver.add_clause(a_outs)
        if weight < self.n:
            card = CardEnc.atmost(lits=qv, bound=weight, vpool=pool,
                                  encoding=EncType.seqcounter)
            for cl in card.clauses:
                solver.add_clause(cl)
        sat, model = solver.solve()
        if not sat:
            return None
        return np.array([1 if model[q] else 0 for q in qv], dtype=np.uint8)

    # ---- verification ----
    def verify(self, v, H, L):
        assert not ((H @ v) % 2).any(), "not a cycle"
        assert ((L @ v) % 2).any(), "functional pairing zero"
        return True

    def torus_nontrivial(self, v, Lx: int) -> bool:
        """Embed the window vector into the torus Z_Lx x Z_p (X-sector
        system with x wrapped) and test v not in the boundary span."""
        p, W = self.p, self.W
        assert Lx >= W
        n = 2 * Lx * p

        def tid(blk, c, j):
            return blk * Lx * p + (c % Lx) * p + (j % p)

        vt = np.zeros(n, dtype=np.uint8)
        for blk in (0, 1):
            for c in range(W):
                for j in range(p):
                    if v[self.vid(blk, c, j)]:
                        vt[tid(blk, c, j)] = 1
        # torus cycle check
        for c in range(Lx):
            for j in range(p):
                s = 0
                for (dx, dy) in P_SUPP:
                    s ^= vt[tid(0, c - dx, j - dy)]
                for (dx, dy) in Q_SUPP:
                    s ^= vt[tid(1, c - dx, j - dy)]
                assert s == 0, "torus cycle fails"
        rows = []
        for tc in range(Lx):
            for tj in range(p):
                row = np.zeros(n, dtype=np.uint8)
                for (dx, dy) in R_SUPP:
                    row[tid(0, tc + dx, tj + dy)] ^= 1
                for (dx, dy) in S_SUPP:
                    row[tid(1, tc + dx, tj + dy)] ^= 1
                rows.append(row)
        Bt = np.array(rows, dtype=np.uint8)
        rb = rank_f2(Bt)
        rbv = rank_f2(np.vstack([Bt, vt]))
        return rbv > rb


def probe_1212_witnesses(nmax: int = 12, log=print):
    """Hunt weight-18 X-logicals of (12,12); report cyclic x/y gap
    structure of the union support."""
    c = TowerCode("m1212", (12, 12), frozenset(A_L), frozenset(B_L))
    assert c.k == 12
    HZ = c.HZ
    zreps = c.zreps
    n = c.n
    pool_found = []
    pool = IDPool()
    qv = [pool.id() for _ in range(n)]
    solver = pycryptosat.Solver()
    for row in HZ:
        idx = np.flatnonzero(row)
        solver.add_xor_clause([qv[i] for i in idx], False)
    a_outs = []
    for row in zreps:
        idx = np.flatnonzero(row)
        if idx.size == 0:
            continue
        a = pool.id()
        solver.add_xor_clause([qv[i] for i in idx] + [a], False)
        a_outs.append(a)
    solver.add_clause(a_outs)
    card = CardEnc.atmost(lits=qv, bound=18, vpool=pool,
                          encoding=EncType.seqcounter)
    for cl in card.clauses:
        solver.add_clause(cl)

    def gaps(cols_set, size):
        cols = sorted(cols_set)
        if not cols:
            return []
        gs = []
        for i in range(len(cols)):
            nxt = cols[(i + 1) % len(cols)]
            gap = (nxt - cols[i]) % size
            if i == len(cols) - 1:
                gap = (cols[0] + size) - cols[i]
            gs.append(gap - 1)
        return gs

    results = []
    for it in range(nmax):
        sat, model = solver.solve()
        if not sat:
            log(f"  no further weight-18 witnesses after {it}")
            break
        v = np.array([1 if model[q] else 0 for q in qv], dtype=np.uint8)
        wt = int(v.sum())
        # x-support union over blocks: qubit id = blk*144 + i*12 + j
        xcols = set()
        ycols = set()
        for blk in (0, 1):
            for i in range(12):
                for j in range(12):
                    if v[blk * 144 + i * 12 + j]:
                        xcols.add(i)
                        ycols.add(j)
        gx = max(gaps(xcols, 12)) if xcols else -1
        gy = max(gaps(ycols, 12)) if ycols else -1
        assert not ((HZ @ v) % 2).any()
        assert ((zreps @ v) % 2).any()
        results.append({"weight": wt, "max_xgap": gx, "max_ygap": gy,
                        "support": [int(i) for i in np.flatnonzero(v)]})
        log(f"  witness {it}: weight {wt}, max cyclic x-gap {gx}, "
            f"y-gap {gy}")
        # block this solution (its support pattern)
        blocking = [-qv[i] if v[i] else qv[i] for i in range(n)]
        solver.add_clause(blocking)
        pool_found.append(v)
    return results, pool_found


def direct_probe(p: int, W: int, wmax: int):
    """Probe C: exact min nontrivial compact weight at period p with
    x-extent < W, searched from wmax downward.  Writes incremental
    JSON checkpoints."""
    t0 = time.time()
    tag = f"p{p}_W{W}"
    ckpt = DATA / f"s1_probe_{tag}.json"
    cw = CylWindow(p, W)
    print(f"probe p={p} W={W} wmax={wmax}: building spaces...", flush=True)
    H = cw.build_H()
    Bnd = cw.build_boundary_in_window()
    L, dimZ, dimB = cw.class_functionals(H, Bnd)
    print(f"  dimZ={dimZ} dimB={dimB} classes={dimZ-dimB} "
          f"functionals={L.shape[0]} ({time.time()-t0:.0f} s)", flush=True)
    state = {"p": p, "W": W, "wmax": wmax, "dimZ": dimZ, "dimB": dimB,
             "classes": dimZ - dimB, "steps": []}
    if dimZ == dimB:
        state["verdict"] = "barren (Z == B)"
        ckpt.write_text(json.dumps(state, indent=1))
        print("  barren: no nontrivial compact cycles at any weight",
              flush=True)
        return

    import resource
    RSS_CAP = 2 * 1024 ** 3  # bytes (ru_maxrss is bytes on macOS)

    def rss_ok():
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < RSS_CAP

    # upward search: every UNSAT step banks "floor > w"; first SAT is
    # exact (previous w-1 already UNSAT).
    best = None
    for w in range(4, wmax + 1):
        ts = time.time()
        v = cw._solve_at(H, L, w)
        dt = round(time.time() - ts, 1)
        if v is None:
            state["steps"].append({"w": w, "sat": False, "s": dt})
            state["floor_gt"] = w
            print(f"  w<={w}: UNSAT ({dt} s) — floor > {w}", flush=True)
        else:
            wt = int(v.sum())
            cw.verify(v, H, L)
            Lx = p * ((W + p) // p) if p % 3 == 0 else 3 * ((W + 11) // 3)
            tv = cw.torus_nontrivial(v, Lx=max(Lx, W + 4))
            state["steps"].append({"w": w, "sat": True, "weight": wt,
                                   "s": dt,
                                   "torus_nontrivial": bool(tv),
                                   "support": [int(i) for i in
                                               np.flatnonzero(v)]})
            best = wt
            print(f"  w<={w}: SAT weight {wt} ({dt} s) torus-nt={tv} "
                  f"— EXACT min (prior weights UNSAT)", flush=True)
            ckpt.write_text(json.dumps(state, indent=1))
            break
        ckpt.write_text(json.dumps(state, indent=1))
        if not rss_ok():
            state["aborted"] = "RSS cap"
            ckpt.write_text(json.dumps(state, indent=1))
            print("  RSS cap reached — clean abort with checkpoint",
                  flush=True)
            return
    state["min_nontrivial_extent_lt_W"] = best
    state["wall_s"] = round(time.time() - t0, 1)
    ckpt.write_text(json.dumps(state, indent=1))
    print(f"probe done: min (extent < {W}) = {best} "
          f"({state['wall_s']} s); wrote {ckpt}", flush=True)


def main():
    t0 = time.time()
    out = {}
    args = set(sys.argv[1:])
    if sys.argv[1:2] == ["probe"]:
        direct_probe(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
        return

    print("== A: calibrations ==", flush=True)
    cal = []
    for p, W, wmax, expect in ((3, 10, 8, 6), (5, 12, 11, None),
                               (6, 14, 13, 12), (7, 12, 15, None)):
        cw = CylWindow(p, W)
        H = cw.build_H()
        Bnd = cw.build_boundary_in_window()
        L, dimZ, dimB = cw.class_functionals(H, Bnd)
        print(f"  p={p} W={W}: dimZ={dimZ} dimB={dimB} "
              f"classes={dimZ-dimB}", flush=True)
        if expect is None:
            assert dimZ == dimB, (p, dimZ, dimB, "Theorem A violated")
            cal.append({"p": p, "W": W, "dimZ": dimZ, "dimB": dimB,
                        "verdict": "barren (Z == B)"})
            print(f"    barren period confirmed (Theorem A)", flush=True)
            continue
        w, v = cw.sat_min_nontrivial(H, L, wmax,
                                     log=lambda s: print(s, flush=True))
        assert w == expect, (p, w, expect)
        cw.verify(v, H, L)
        # embedding torus must carry the omega-line: need 3 | Lx
        # (switch-on condition ord(alpha) = 3 | Lx; at 3 not dividing Lx
        # the class is torus-invisible and the control is blind)
        Lx = 3 * ((W + 11) // 3)
        tv = cw.torus_nontrivial(v, Lx=Lx)
        print(f"    min nontrivial = {w} (expect {expect}) "
              f"torus({Lx},{p})-nontrivial={tv}", flush=True)
        assert tv
        cal.append({"p": p, "W": W, "dimZ": dimZ, "dimB": dimB,
                    "min_nontrivial": w, "expect": expect})
    out["calibration"] = cal

    print("\n== B: (12,12) weight-18 witness gap probe ==", flush=True)
    res1212, vecs = probe_1212_witnesses(
        nmax=12, log=lambda s: print(s, flush=True))
    out["w1212"] = res1212

    gapped = [r for r in res1212 if r["max_xgap"] >= 4]
    print(f"  witnesses with x-gap >= 4: {len(gapped)} / {len(res1212)}")

    if gapped and "nowindow" not in args:
        print("\n== B2: windowize the first gapped witness ==", flush=True)
        r = gapped[0]
        v = np.zeros(288, dtype=np.uint8)
        v[np.array(r["support"])] = 1
        # rotate so the gap sits at the right edge: find gap start
        xcols = sorted({(i) for blk in (0, 1) for i in range(12)
                        for j in range(12) if v[blk * 144 + i * 12 + j]})
        # find a rotation with columns in [0, 12-5]
        best_rot = None
        for rot in range(12):
            cols = sorted(((c - rot) % 12) for c in xcols)
            if cols[-1] <= 12 - 5:
                best_rot = rot
                break
        assert best_rot is not None, "gap >= 4 must allow such a rotation"
        W = 12
        cw = CylWindow(12, W)
        vwin = np.zeros(cw.n, dtype=np.uint8)
        for blk in (0, 1):
            for i in range(12):
                for j in range(12):
                    if v[blk * 144 + i * 12 + j]:
                        vwin[cw.vid(blk, (i - best_rot) % 12, j)] = 1
        H = cw.build_H()
        Bnd = cw.build_boundary_in_window()
        L, dimZ, dimB = cw.class_functionals(H, Bnd)
        assert not ((H @ vwin) % 2).any(), "windowized vector not a cycle"
        nz = bool(((L @ vwin) % 2).any())
        tv = cw.torus_nontrivial(vwin, Lx=24)
        wt = int(vwin.sum())
        print(f"  windowized weight {wt}: cylinder-nontrivial={nz} "
              f"(functionals), torus(24,12)-nontrivial={tv}", flush=True)
        out["windowized"] = {"weight": wt, "cyl_nontrivial": nz,
                             "torus24_nontrivial": tv,
                             "classes_at_p12": dimZ - dimB}
        if nz:
            print("  >>> weight-18 nontrivial compact period-12 phase "
                  "EXISTS: 'floor >= 2p for all p' is REFUTED at p=12",
                  flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s1_cylfloor.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s1_cylfloor.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
