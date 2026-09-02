#!/usr/bin/env python3
"""A42 S4 — the FIBER DECOMPOSITION of the period-3q cylinder (the
L-band mixed half in its sharpest form).

p = 3q with 3 coprime to q.  Z_{3q} = Z_3 x Z_q, and F2[Z_3] = F2 x F4
(semisimple), so

    R_p = F2[Z_{3q}] = F2[Z_q]  x  F4[Z_q]
                       (barren)    (omega)

with the fibre over a cell (block, column, j in Z_q) of the period-q
cylinder being the three cells {j, j+q, j+2q} of the 3q-cylinder.
Per fibre, the triple (v0, v1, v2) in F2^3 maps bijectively to
(s, mu) in F2 x F4 with s = v0+v1+v2 (the barren bit) and
mu = v0 + zeta v1 + zeta^2 v2 (the omega content).  Hamming weight per
fibre:  0 (000), 1 (singleton: s=1, mu!=0), 2 (pair: s=0, mu!=0),
3 (full: s=1, mu=0).  Hence for every column, with S = {fibres with
mu != 0} and s = {fibres with s = 1}:

    pure(lambda) = 2|S|,   pure'(z') = 3|s|,
    wt(v) = n1 + 2 n2 + 3 n3 = 2|S| + 3|s| - 4|S cap s|.

This script (A) verifies the fibre identity against the banked CRT
pure table exhaustively at p = 6, 12 and on all 2^24 columns at
p = 24 (vectorized), (B) profiles every atlas cycle at p = 3, 6 by
(n1, n2, n3), and (C) runs PROBE-tier SAT searches at p = 6, 12 for
the extreme mixed shapes: section cycles (all fibres singleton) and
cycles with a prescribed omega-support |S| — the hiding-mass
frontier (weight, |S|) of the nontrivial sector.
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

import a42_lib as AL  # noqa: E402
import importlib.util

from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402
import pycryptosat  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


OM = _load("a42_s3b_omega_racer")
CF = _load("a42_s1_cylfloor")
AT = _load("a40_s4_phase_atlas")
DATA = LAB / "data" / "a42"
OMEGA2 = 0b111


def fibre_counts(col: int, p: int):
    q = p // 3
    n = [0, 0, 0, 0]
    for j in range(q):
        k = ((col >> j) & 1) + ((col >> (j + q)) & 1) + \
            ((col >> (j + 2 * q)) & 1)
        n[k] += 1
    return n


def part_A(log):
    out = {}
    for p in (6, 12):
        P, dim, Lmod = OM.crt_pure_table(p)
        q = p // 3
        Bmod = (1 << q) | 1
        bad = 0
        for col in range(1 << p):
            lam = AL.pmod(col, Lmod)
            n = fibre_counts(col, p)
            nS = n[1] + n[2]
            ns = n[1] + n[3]
            zp = AL.pmod(col, Bmod)
            if P[lam] != 2 * nS or bin(zp).count("1") != ns:
                bad += 1
        assert bad == 0, (p, bad)
        log(f"A p={p}: fibre identity pure(lambda) = 2|S| and "
            f"|fold| = |s| on all {1 << p} columns")
        out[f"p{p}"] = "exact on all columns"
    # p = 24 vectorized: lambda = linear in the column bits
    p = 24
    P, dim, Lmod = OM.crt_pure_table(p)
    q = 8
    T = []
    for k in range(3):
        T.append(np.array([AL.pmod(b << (8 * k), Lmod)
                           for b in range(256)], dtype=np.int64))
    v = np.arange(1 << p, dtype=np.int64)
    lam = T[0][v & 255] ^ T[1][(v >> 8) & 255] ^ T[2][(v >> 16) & 255]
    nS = np.zeros(v.size, dtype=np.int64)
    for j in range(q):
        b0 = (v >> j) & 1
        b1 = (v >> (j + q)) & 1
        b2 = (v >> (j + 2 * q)) & 1
        nS += ((b0 ^ b1) | (b1 ^ b2))
    assert (P[lam] == 2 * nS).all()
    log(f"A p=24: fibre identity pure(lambda) = 2|S| on all 2^24 "
        f"columns (dim Lambda_3 = {dim})")
    out["p24"] = "exact on all 2^24 columns"
    return out


def part_B(log):
    out = {}
    for p, cap in ((3, 8), (6, 13)):
        rows, _ = AT.atlas("AB", p, cap, keep_pts=True)
        q = p // 3
        prof = {}
        for row in rows:
            fib = {}
            for (c, y, blk) in row["pts"]:
                key = (blk, c, y % q)
                fib[key] = fib.get(key, 0) + 1
            n = [0, 0, 0, 0]
            for k in fib.values():
                n[k] += 1
            assert n[1] + 2 * n[2] + 3 * n[3] == row["weight"]
            key = (int(row["nontrivial"]), row["weight"],
                   n[1], n[2], n[3])
            prof[key] = prof.get(key, 0) + 1
        tab = sorted(prof.items())
        log(f"B p={p}: (nontrivial, weight, n1, n2, n3) -> count")
        for k, cnt in tab:
            nS = k[2] + k[3]
            ns = k[2] + k[4]
            log(f"    nt={k[0]} w={k[1]:2d} (n1,n2,n3)=({k[2]},{k[3]},"
                f"{k[4]}) |S|={nS:2d} |s|={ns:2d} : {cnt}")
        out[f"p{p}"] = [{"nontrivial": k[0], "weight": k[1],
                         "n1": k[2], "n2": k[3], "n3": k[4],
                         "count": cnt} for k, cnt in tab]
    return out


class FibreProbe:
    """PROBE tier: CylWindow SAT with fibre-shape side constraints."""

    def __init__(self, p, W, time_limit=600.0):
        self.p, self.W, self.q = p, W, p // 3
        self.cw = CF.CylWindow(p, W)
        self.H = self.cw.build_H()
        Bnd = self.cw.build_boundary_in_window()
        self.L, self.dimZ, self.dimB = self.cw.class_functionals(
            self.H, Bnd)
        assert self.L.shape[0] > 0
        self.time_limit = time_limit

    def fibres(self):
        for blk in (0, 1):
            for c in range(self.W):
                for j in range(self.q):
                    yield [self.cw.vid(blk, c, j + i * self.q)
                           for i in range(3)]

    def solve(self, wmax, section=False, Smax=None, n3zero=False,
              nontrivial=True, nopair=False):
        pool = IDPool()
        qv = [pool.id() for _ in range(self.cw.n)]
        solver = pycryptosat.Solver(time_limit=self.time_limit)
        for row in self.H:
            idx = np.flatnonzero(row)
            if idx.size:
                solver.add_xor_clause([qv[i] for i in idx], False)
        if nontrivial:
            a_outs = []
            for row in self.L:
                idx = np.flatnonzero(row)
                a = pool.id()
                solver.add_xor_clause([qv[i] for i in idx] + [a], False)
                a_outs.append(a)
            solver.add_clause(a_outs)
        else:
            solver.add_clause(qv)      # nonzero
        if nopair:
            for (a, b, c) in self.fibres():
                A, B, C = qv[a], qv[b], qv[c]
                solver.add_clause([-A, -B, C])
                solver.add_clause([-A, -C, B])
                solver.add_clause([-B, -C, A])
        card = CardEnc.atmost(lits=qv, bound=wmax, vpool=pool,
                              encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)
        ncs = []
        for (a, b, c) in self.fibres():
            A, B, C = qv[a], qv[b], qv[c]
            if section or n3zero:
                solver.add_clause([-A, -B, -C])
            if section:
                solver.add_clause([-A, -B])
                solver.add_clause([-A, -C])
                solver.add_clause([-B, -C])
            if Smax is not None:
                nc = pool.id()
                # not nc -> fibre constant
                solver.add_clause([nc, -A, B])
                solver.add_clause([nc, A, -B])
                solver.add_clause([nc, -B, C])
                solver.add_clause([nc, B, -C])
                ncs.append(nc)
        if Smax is not None:
            card = CardEnc.atmost(lits=ncs, bound=Smax, vpool=pool,
                                  encoding=EncType.seqcounter)
            for cl in card.clauses:
                solver.add_clause(cl)
        t0 = time.time()
        sat, model = solver.solve()
        dt = round(time.time() - t0, 1)
        if sat is None:
            return "timeout", None, dt
        if not sat:
            return "unsat", None, dt
        v = np.array([1 if model[q] else 0 for q in qv], dtype=np.uint8)
        if nontrivial:
            self.cw.verify(v, self.H, self.L)
        else:
            assert not ((self.H @ v) % 2).any() and v.any()
        return "sat", v, dt

    def profile(self, v):
        n = [0, 0, 0, 0]
        for fb in self.fibres():
            k = int(v[fb[0]]) + int(v[fb[1]]) + int(v[fb[2]])
            n[k] += 1
        return n


def part_C(log, p, W, wmax, budget_s):
    """The hiding-mass frontier at period p (probe tier)."""
    t_all = time.time()
    fp = FibreProbe(p, W, time_limit=budget_s)
    q = p // 3
    log(f"C p={p} W={W}: dimZ={fp.dimZ} dimB={fp.dimB} "
        f"classes={fp.dimZ - fp.dimB}")
    rec = {"p": p, "W": W, "classes": fp.dimZ - fp.dimB, "tier": "probe",
           "section": [], "Smax": []}
    # (i) section cycles: min weight nontrivial with all fibres <= 1
    for w in range(2 * q, wmax + 1, 2):
        st, v, dt = fp.solve(w, section=True)
        rec["section"].append({"w": w, "status": st, "s": dt})
        log(f"  section, weight <= {w}: {st} ({dt} s)")
        if st == "sat":
            n = fp.profile(v)
            rec["section"][-1]["profile"] = n
            rec["section"][-1]["weight"] = int(v.sum())
            log(f"    FOUND weight {int(v.sum())} profile {n} "
                f"(all-singleton nontrivial cycle)")
            break
        if st == "timeout":
            break
    # (ii) |S| <= Smax at weight <= wmax: which omega-support sizes
    # carry a floor-weight nontrivial cycle?
    for Smax in range(3 * q, wmax + 1):
        st, v, dt = fp.solve(wmax, Smax=Smax)
        e = {"Smax": Smax, "w": wmax, "status": st, "s": dt}
        if st == "sat":
            n = fp.profile(v)
            e["profile"] = n
            e["weight"] = int(v.sum())
            e["S"] = n[1] + n[2]
            e["s"] = n[1] + n[3]
        rec["Smax"].append(e)
        log(f"  |S| <= {Smax}, weight <= {wmax}: {st} ({dt} s)"
            + (f" profile {e['profile']} |S|={e['S']} |s|={e['s']} "
               f"wt={e['weight']}" if st == "sat" else ""))
        if st == "sat" or st == "timeout":
            break
    rec["wall_s"] = round(time.time() - t_all, 1)
    return rec


def part_D(log, p, W, wmax, budget_s):
    """Section exclusion (any NONZERO section cycle, trivial or not)
    and the halving-tight stratum (no pair fibres) floor — probe
    tier."""
    fp = FibreProbe(p, W, time_limit=budget_s)
    q = p // 3
    rec = {"p": p, "W": W, "tier": "probe", "any_section": [],
           "nopair": []}
    log(f"D p={p} W={W}: classes={fp.dimZ - fp.dimB}")
    st, v, dt = fp.solve(fp.cw.n, section=True, nontrivial=False)
    rec["any_section"].append({"w": "any", "status": st, "s": dt})
    log(f"  ANY nonzero section cycle (no weight cap, trivial allowed):"
        f" {st} ({dt} s)")
    if st == "sat":
        n = fp.profile(v)
        rec["any_section"][-1]["profile"] = n
        rec["any_section"][-1]["weight"] = int(v.sum())
        nz = bool(((fp.L @ v) % 2).any())
        rec["any_section"][-1]["nontrivial"] = nz
        log(f"    found weight {int(v.sum())} profile {n} "
            f"nontrivial={nz}")
    for w in range(2 * q, wmax + 1, 2):
        st, v, dt = fp.solve(w, nopair=True)
        e = {"w": w, "status": st, "s": dt}
        if st == "sat":
            n = fp.profile(v)
            e["profile"] = n
            e["weight"] = int(v.sum())
        rec["nopair"].append(e)
        log(f"  no-pair stratum, weight <= {w}: {st} ({dt} s)"
            + (f" profile {e['profile']}" if st == "sat" else ""))
        if st in ("sat", "timeout"):
            break
    return rec


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "AB"
    logf = open(DATA / f"s4_fiber_{which}.log", "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    out = {}
    if "A" in which:
        out["A"] = part_A(log)
    if "B" in which:
        out["B"] = part_B(log)
    if which.startswith("C"):
        p = int(sys.argv[2])
        W = int(sys.argv[3])
        wmax = int(sys.argv[4])
        budget = float(sys.argv[5]) if len(sys.argv) > 5 else 600.0
        out["C"] = part_C(log, p, W, wmax, budget)
    if which.startswith("D"):
        p = int(sys.argv[2])
        W = int(sys.argv[3])
        wmax = int(sys.argv[4])
        budget = float(sys.argv[5]) if len(sys.argv) > 5 else 600.0
        out["D"] = part_D(log, p, W, wmax, budget)
    (DATA / f"s4_fiber_{which}.json").write_text(json.dumps(out, indent=1))
    log(f"wrote {DATA / f's4_fiber_{which}.json'}")


if __name__ == "__main__":
    main()
