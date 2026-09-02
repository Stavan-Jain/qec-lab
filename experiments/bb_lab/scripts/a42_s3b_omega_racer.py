#!/usr/bin/env python3
"""A42 S3b — the OMEGA-QUOTIENT racer: the pure (h = 0) half of
L-band, run directly on the omega-factor.

Object: at p = 3m*2^a, the pure-lift cycles (all non-omega content
zero) are in bijection with compact omega-syzygies sigma of the pair
(Abar_omega, Bbar_omega) over Lambda_a = F2[y]/((y^2+y+1)^{2^a}),
and the true weight of the pure lift is  sum_cols pure(lambda_col)
with pure(lambda) = wt of the unique all-else-zero CRT lift.  So the
S4/S1h automaton, rebuilt with columns IN Lambda_a (the y-rotation
becomes multiplication by ybar), per-column costs pure(lambda), and
the S2e full-depth branch registers acting on the columns
THEMSELVES, enumerates all pure cycles by weight — with a state
space of 5*dim(Lambda_a) bits (40 at a = 2, 20 at a = 1) that never
sees the barren direction.  This is the register-quotient racer in
its exact-cost form on the h = 0 sector:

  * completing level c in both branches certifies "no nontrivial
    PURE compact cycle of weight <= c" (per-branch registers with
    empty joint kernel on H — s2_registers, a = 1 and a = 2);
  * its per-level state counts ARE the omega-direction growth curve
    (the S3 quotient hypothesis: subexponential, by Theorem H's
    finite class inventory);
  * at p = 18 (a = 1, dim 4) the full racer's 90-bit state is
    infeasible but this one is 20-bit — the first quantitative
    instrument at the r = 3 member period.

Controls: p = 3 pure min nontrivial = 6 (the 2p object is the pure
3-slot lift); p = 6: every PURE atlas cycle (zero barren content in
both blocks) replays through the omega automaton column-by-column
(forced columns, cost = weight, register = the jet replay register
of the full cycle), and the racer's return spectrum at cap 13
matches the pure-subset spectrum of the atlas; p = 12: the realized
weight-24 pure object must appear as a nonzero-register return at
24 when cap >= 24.

Cost-table construction: CRT idempotent — c_k = (y-lift of t^k) *
e mod (y^p + 1) with e = 1 on the omega-factor, 0 on every other
factor; pure(lambda) = popcount(sum of c_k over lambda's bits);
asserted: c_k == t^k on the omega-factor, == 0 on the complement,
pure >= 2 for lambda != 0 (cost-1 rigidity's pure face), rotation
invariance pure(ybar * lambda) = pure(lambda), and (p = 18) the
m-scaling theorem pure_18 = 3 * pure_6 checked table-wide.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a42_lib as AL  # noqa: E402
import importlib.util


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AT = _load("a40_s4_phase_atlas")
JR = _load("a42_s2_jet_racer")
CR = _load("a42_s3_corridor_racer")

DATA = LAB / "data" / "a42"
RSS_CAP = int(2.4 * 1024 ** 3)
CHUNK = 4096
OMEGA2 = 0b111


def crt_pure_table(p: int):
    """pure(lambda) for lambda in Lambda_a at period p = 3m*2^a:
    the weight of the unique lift with every non-omega factor zero.
    Returns (P, dim, Lmod)."""
    a = AL.v2(p)
    Lmod = OMEGA2
    for _ in range(a):
        Lmod = AL.pmul(Lmod, Lmod)
    dim = AL.pdeg(Lmod)
    ypmod = (1 << p) | 1          # y^p + 1
    comp, rem = AL.pdivmod(ypmod, Lmod)   # the barren complement
    assert rem == 0
    # e = comp * u with u = comp^{-1} mod Lmod (unit since coprime)
    compL = AL.pmod(comp, Lmod)
    u = None
    for cand in range(1, 1 << dim):
        if AL.pmod(AL.pmul(compL, cand), Lmod) == 1:
            u = cand
            break
    assert u is not None, "complement not invertible mod Lmod"
    e = AL.pmod(AL.pmul(comp, u), ypmod)
    assert AL.pmod(e, Lmod) == 1 and AL.pmod(e, comp) == 0
    basis = []
    for k in range(dim):
        ck = AL.pmod(AL.pmul(e, 1 << k), ypmod)
        assert AL.pmod(ck, Lmod) == (1 << k)
        assert AL.pmod(ck, comp) == 0
        basis.append(ck)
    P = np.zeros(1 << dim, dtype=np.int64)
    for lam in range(1, 1 << dim):
        c = 0
        for k in range(dim):
            if (lam >> k) & 1:
                c ^= basis[k]
        P[lam] = bin(c).count("1")
    assert (P[1:] >= 2).all(), "pure >= 2 violated"
    return P, dim, Lmod


class OmegaRacer:
    """Pure-half racer over Lambda_a columns at period p, one
    variety branch per instance (registers as in the jet engine,
    acting on the columns themselves)."""

    def __init__(self, p: int, branch: int):
        assert p % 3 == 0
        self.p = p
        self.branch = branch
        self.a = AL.v2(p)
        self.P, self.dim, self.Lmod = crt_pure_table(p)
        nl = 1 << self.dim
        assert 5 * self.dim <= 64

        def red(z):
            return AL.pmod(z, self.Lmod)

        def mul(x_, y_):
            return red(AL.pmul(x_, y_))

        self.red, self.Lmul = red, mul
        MT = np.zeros((nl, nl), dtype=np.uint8)
        for x_ in range(nl):
            for y_ in range(nl):
                MT[x_, y_] = mul(x_, y_)
        self.MT = MT
        inv = {}
        for x_ in range(1, nl):
            for y_ in range(1, nl):
                if mul(x_, y_) == 1:
                    inv[x_] = y_
                    break
        ybar = red(0b10)
        ybi = inv[ybar]
        # branch points: roots of Pbar_omega = 1 + ybar^{-1} +
        # x^{-3} ybar with residue != 1 (verbatim the jet init)
        roots = []
        for x_ in range(1, nl):
            if x_ not in inv:
                continue
            xi3 = mul(inv[x_], mul(inv[x_], inv[x_]))
            val = 1 ^ ybi ^ mul(xi3, ybar)
            if val == 0:
                resid = AL.pmod(x_, OMEGA2)
                if resid != 1:
                    roots.append((resid, x_))
        roots.sort()
        assert len(roots) == 2, ("branch points", roots)
        self.xstar = roots[branch][1]
        self.xresid = roots[branch][0]
        y3 = mul(ybar, mul(ybar, ybar))
        qv = 1 ^ inv[self.xstar] ^ mul(self.xstar, y3)
        assert AL.pmod(qv, OMEGA2) == 0, "xstar not on the variety"
        self.XI = inv[self.xstar]
        self.MXI = MT[self.XI]
        # step multiplier rows: acc = f0*(1 + ybar^{-1}) + o0*ybar^3
        #                             + o1 + a;  new_f = ybar^{-1}*acc
        self.R_F0 = MT[1 ^ ybi]
        self.R_O0 = MT[y3]
        self.R_TOP = MT[ybi]
        # rotation co-multipliers ybar^s, s in 0..p-1
        self.yrot = np.zeros(p, dtype=np.uint8)
        acc = 1
        for s in range(p):
            self.yrot[s] = acc
            acc = mul(acc, ybar)
        assert mul(int(self.yrot[p - 1]), ybar) == 1, "ybar^p != 1"
        self.MROT = [MT[self.yrot[s]] for s in range(p)]
        self.maskc = np.uint64(nl - 1)

    # ------------- scalar reference step (for validation) ---------
    def step_ref(self, cols, aa):
        f0, f1, f2, o0, o1 = cols
        acc = int(self.R_F0[f0]) ^ int(self.R_O0[o0]) ^ o1 ^ aa
        new_f = int(self.R_TOP[acc])
        cost = int(self.P[aa]) + int(self.P[new_f])
        return (f1, f2, new_f, o1, aa), cost, new_f

    # ------------- vectorized ------------------------------------
    def unpack(self, states):
        d = self.dim
        return [((states >> np.uint64(d * k)) & self.maskc)
                .astype(np.int64) for k in range(5)]

    def expand(self, states: np.ndarray, regs: np.ndarray):
        d = self.dim
        f0, f1, f2, o0, o1 = self.unpack(states)
        base = (self.R_F0[f0].astype(np.int64)
                ^ self.R_O0[o0].astype(np.int64) ^ o1)
        nreg_base = self.MXI[regs]
        outs = []
        for aa in range(1 << d):
            acc = base ^ aa
            new_f = self.R_TOP[acc].astype(np.int64)
            cost = self.P[aa] + self.P[new_f]
            ns = (f1.astype(np.uint64)
                  | (f2.astype(np.uint64) << np.uint64(d))
                  | (new_f.astype(np.uint64) << np.uint64(2 * d))
                  | (o1.astype(np.uint64) << np.uint64(3 * d))
                  | np.uint64(aa << (4 * d)))
            nreg = nreg_base ^ np.uint8(aa)
            outs.append((ns, nreg, cost))
        return outs

    def canon(self, states: np.ndarray, regs: np.ndarray):
        d = self.dim
        bs = states.copy()
        br = regs.copy()
        shifts = [np.uint64(d * k) for k in range(5)]
        cols = self.unpack(states)
        for s in range(1, self.p):
            M = self.MROT[s]
            rs = np.zeros_like(states)
            for k in range(5):
                rs |= M[cols[k]].astype(np.uint64) << shifts[k]
            rr = M[regs]
            better = (rs < bs) | ((rs == bs) & (rr < br))
            bs = np.where(better, rs, bs)
            br = np.where(better, rr, br)
        return bs, br

    def validate(self, nsamples=600, rng=None):
        rng = rng or np.random.default_rng(5)
        d = self.dim
        nl = 1 << d
        # P-table rotation invariance
        ybar = self.red(0b10)
        for lam in range(nl):
            assert self.P[self.Lmul(ybar, lam)] == self.P[lam]
        # vector step vs scalar reference
        for _ in range(nsamples):
            cols = [int(rng.integers(0, nl)) for _ in range(5)]
            aa = int(rng.integers(0, nl))
            sv = np.array([sum(c << (d * k)
                               for k, c in enumerate(cols))],
                          dtype=np.uint64)
            rv = np.array([int(rng.integers(0, nl))], dtype=np.uint8)
            ns, nreg, co = self.expand(sv, rv)[aa]
            st2, cost2, _ = self.step_ref(cols, aa)
            got = [(int(ns[0]) >> (d * k)) & (nl - 1)
                   for k in range(5)]
            assert tuple(got) == st2
            assert int(co[0]) == cost2
            assert int(nreg[0]) == int(self.MXI[rv[0]]) ^ aa
        # canon commutation: rotate columns+reg by ybar^s, canon equal
        for _ in range(300):
            cols = [int(rng.integers(0, nl)) for _ in range(5)]
            r0 = int(rng.integers(0, nl))
            s = int(rng.integers(1, self.p))
            sv = np.array([sum(c << (d * k)
                               for k, c in enumerate(cols))],
                          dtype=np.uint64)
            rv = np.array([r0], dtype=np.uint8)
            c1s, c1r = self.canon(sv, rv)
            M = self.MROT[s]
            rcols = [int(M[c]) for c in cols]
            sv2 = np.array([sum(c << (d * k)
                                for k, c in enumerate(rcols))],
                           dtype=np.uint64)
            rv2 = np.array([int(M[r0])], dtype=np.uint8)
            c2s, c2r = self.canon(sv2, rv2)
            assert int(c1s[0]) == int(c2s[0]) \
                and int(c1r[0]) == int(c2r[0])
        return True

    def replay(self, l1, l2):
        """March the omega automaton over a pure cycle given by its
        omega column dicts l1 (input block) and l2 (forced block);
        assert forced columns; return (cost, final register)."""
        cs = list(l1) + list(l2) or [0]
        lo, hi = min(cs) - 6, max(cs) + 6
        cols = (0, 0, 0, 0, 0)
        reg = 0
        cost = 0
        for g in range(lo, hi + 1):
            aa = l1.get(g + 1, 0)
            cols, c2, newf = self.step_ref(cols, aa)
            assert newf == l2.get(g + 3, 0), (g, newf)
            reg = int(self.MXI[reg]) ^ aa
            cost += c2
        assert cols == (0, 0, 0, 0, 0)
        return cost, reg

    # ------------- the BFS ---------------------------------------
    def run(self, cap: int, log=print, ckpt_path=None,
            time_budget_s=None, growth_est=2.5):
        t0 = time.time()
        zstate = np.uint64(0)
        buckets: dict[int, list] = {
            0: [(np.array([0], dtype=np.uint64),
                 np.zeros(1, dtype=np.uint8))]}
        seen_s = np.array([0], dtype=np.uint64)
        seen_r = np.zeros(1, dtype=np.uint8)
        returns: dict[int, set] = {}
        state = {"p": self.p, "cap": cap, "branch": self.branch,
                 "xstar": int(self.xstar), "levels": [],
                 "engine": "omega-quotient pure-half racer",
                 "lane": "AB (theta'-dual lane by banked duality)"}
        aborted = False
        last_level_t = 0.0
        for c in range(0, cap + 1):
            if time_budget_s is not None and c > 0:
                if (time.time() - t0) + last_level_t * growth_est \
                        > time_budget_s:
                    state["aborted"] = f"time budget before level {c}"
                    log(f"  time budget — stop before level {c}")
                    aborted = True
                    break
            lt0 = time.time()
            level_novel = 0
            while buckets.get(c):
                arrs = buckets.pop(c)
                sts = np.concatenate([x[0] for x in arrs])
                rgs = np.concatenate([x[1] for x in arrs])
                sts, rgs = JR.JetRacer.dedup(sts, rgs)
                fresh = ~CR.member_pairs(sts, rgs, seen_s, seen_r)
                if c == 0 and level_novel == 0:
                    fresh[:] = True
                sts, rgs = sts[fresh], rgs[fresh]
                if sts.size == 0:
                    break
                level_novel += int(sts.size)
                seen_s = np.concatenate([seen_s, sts])
                seen_r = np.concatenate([seen_r, rgs])
                o = np.lexsort((seen_r, seen_s))
                seen_s, seen_r = seen_s[o], seen_r[o]
                for lo in range(0, sts.size, CHUNK):
                    cst = sts[lo:lo + CHUNK]
                    crg = rgs[lo:lo + CHUNK]
                    outs = self.expand(cst, crg)
                    per_cost: dict[int, list] = {}
                    for aa in range(1 << self.dim):
                        ns, nreg, co = outs[aa]
                        outs[aa] = None
                        if c == 0 and aa == 0:
                            keep = ns != zstate
                            ns, nreg, co = ns[keep], nreg[keep], \
                                co[keep]
                        for w in np.unique(co):
                            cw = c + int(w)
                            if cw > cap:
                                continue
                            sel = co == w
                            zsel = sel & (ns == zstate)
                            if zsel.any() and cw > 0:
                                for rv in np.unique(nreg[zsel]):
                                    returns.setdefault(cw, set()) \
                                        .add(int(rv))
                            per_cost.setdefault(cw, []).append(
                                (ns[sel], nreg[sel]))
                    for cw, lst in per_cost.items():
                        ms = np.concatenate([x[0] for x in lst])
                        mr = np.concatenate([x[1] for x in lst])
                        ms, mr = self.canon(ms, mr)
                        ms, mr = JR.JetRacer.dedup(ms, mr)
                        fr = ~CR.member_pairs(ms, mr, seen_s, seen_r)
                        if fr.any():
                            buckets.setdefault(cw, []).append(
                                (ms[fr], mr[fr]))
                    del outs, per_cost
                    rss = resource.getrusage(
                        resource.RUSAGE_SELF).ru_maxrss
                    if rss > RSS_CAP:
                        state["aborted"] = "RSS"
                        log(f"  RSS cap ({rss//2**20}MB) — clean "
                            "abort")
                        aborted = True
                        break
                if aborted:
                    break
                for cw in list(buckets.keys()):
                    arrs2 = buckets[cw]
                    tot = sum(x[0].size for x in arrs2)
                    if len(arrs2) > 16 or tot > 1 << 20:
                        ms = np.concatenate([x[0] for x in arrs2])
                        mr = np.concatenate([x[1] for x in arrs2])
                        ms, mr = JR.JetRacer.dedup(ms, mr)
                        fr = ~CR.member_pairs(ms, mr, seen_s, seen_r)
                        buckets[cw] = [(ms[fr], mr[fr])] if fr.any() \
                            else []
            if aborted:
                break
            last_level_t = time.time() - lt0
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            state["levels"].append(
                {"c": c, "novel": level_novel,
                 "seen": int(seen_s.size), "rss_mb": rss // 2 ** 20,
                 "t": round(time.time() - t0, 1)})
            nontriv = sorted(w for w, rv in returns.items()
                             if any(v != 0 for v in rv))
            log(f"  level {c}: novel {level_novel} seen "
                f"{seen_s.size} rss {rss//2**20}MB "
                f"t {time.time()-t0:.0f}s returns {sorted(returns)} "
                f"nonzero-reg@{nontriv}")
            state["returns"] = {str(k): sorted(v)
                                for k, v in returns.items()}
            state["nonzero_reg_costs"] = nontriv
            if ckpt_path:
                ckpt_path.write_text(json.dumps(state, indent=1))
        state["wall_s"] = round(time.time() - t0, 1)
        if ckpt_path:
            ckpt_path.write_text(json.dumps(state, indent=1))
        return state


def controls():
    # p = 3, both branches: pure min nontrivial = 6
    mins = []
    for br in (0, 1):
        r = OmegaRacer(3, br)
        r.validate()
        st = r.run(8, log=lambda s: None)
        nt = st.get("nonzero_reg_costs", [])
        mins.append(nt[0] if nt else None)
    m = min(x for x in mins if x is not None)
    print(f"p=3 pure racer: branch minima {mins}, union {m} "
          "(expect 6)", flush=True)
    assert m == 6

    # p = 6: pure atlas subset replay + racer spectrum
    p = 6
    rT = OmegaRacer(p, 0)
    rV = OmegaRacer(p, 1)
    rT.validate()
    rV.validate()
    jT = JR.JetRacer(p, 0)
    jV = JR.JetRacer(p, 1)
    au = AT.Automaton(AT.A_L, AT.B_L, p)
    rows, _ = AT.atlas("AB", p, 13, keep_pts=True)
    # the omega/barren split at p = 6: Lmod = (y^2+y+1)^2,
    # complement = (y+1)^2
    Lmod = rT.Lmod
    compmod, _rem = AL.pdivmod((1 << p) | 1, Lmod)
    assert _rem == 0
    npure = 0
    pure_weights = {}
    for row in rows:
        v1, v2 = SW_cols(row["pts"])
        purecols = all(AL.pmod(c, compmod) == 0
                       for c in list(v1.values()) + list(v2.values()))
        if not purecols:
            continue
        npure += 1
        l1 = {c: AL.pmod(col, Lmod) for c, col in v1.items()}
        l2 = {c: AL.pmod(col, Lmod) for c, col in v2.items()}
        cT, gT = rT.replay(l1, l2)
        cV, gV = rV.replay(l1, l2)
        assert cT == cV == row["weight"], (cT, cV, row["weight"])
        jgT = JR.replay_cycle_through(au, jT, v1, v2)
        jgV = JR.replay_cycle_through(au, jV, v1, v2)
        assert (gT, gV) == (jgT, jgV), "omega vs jet register drift"
        assert ((gT != 0) or (gV != 0)) == row["nontrivial"]
        pure_weights.setdefault(row["weight"], [0, 0])
        pure_weights[row["weight"]][row["nontrivial"]] += 1
    print(f"p=6: {npure} PURE atlas cycles replay through the "
          f"omega automaton (cost = weight, registers = jet "
          f"registers, class verdicts match); weights "
          f"{ {k: tuple(v) for k, v in sorted(pure_weights.items())} }",
          flush=True)
    stT = rT.run(13, log=lambda s: None)
    stV = rV.run(13, log=lambda s: None)
    ntU = sorted(set(stT.get("nonzero_reg_costs", []))
                 | set(stV.get("nonzero_reg_costs", [])))
    atlas_nt = sorted(w for w, v in pure_weights.items() if v[1])
    assert ntU and ntU[0] == 12, (ntU,)
    assert set(ntU) == set(atlas_nt), (ntU, atlas_nt)
    print(f"p=6 pure racer cap 13: nonzero-register returns {ntU} "
          f"== pure-atlas nontrivial weights {atlas_nt} "
          "(min 12 = 2p)", flush=True)
    return True


def SW_cols(pts):
    v1, v2 = {}, {}
    for (c, y, blk) in pts:
        d = v1 if blk == 0 else v2
        d[c] = d.get(c, 0) | (1 << y)
    return v1, v2


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "controls"

    if which == "controls":
        controls()
        return

    if which == "run":
        p = int(sys.argv[2])
        tag = sys.argv[3]
        cap = int(sys.argv[4])
        budget = float(sys.argv[5]) * 3600
        br = {"T": 0, "V": 1}[tag]
        r = OmegaRacer(p, br)
        r.validate()
        print(f"p={p} omega pure racer branch {tag} "
              f"(xstar={r.xstar}, dim {r.dim}): cap {cap}, "
              f"budget {budget/3600:.1f} h", flush=True)
        st = r.run(cap, log=lambda s: print(s, flush=True),
                   ckpt_path=DATA / f"s3b_omega{p}_{tag}_ckpt.json",
                   time_budget_s=budget)
        done = max((lv["c"] for lv in st["levels"]), default=-1)
        nt = st.get("nonzero_reg_costs", [])
        print(f"p={p} branch {tag}: completed level {done}; "
              f"nonzero-register costs {nt}", flush=True)
        (DATA / f"s3b_omega{p}_{tag}.json").write_text(
            json.dumps(st, indent=1))
        return

    if which == "mscale18":
        # m-scaling cross-check: pure_18 == 3 * pure_6 as tables
        P6, d6, L6 = crt_pure_table(6)
        P18, d18, L18 = crt_pure_table(18)
        assert d6 == d18 == 4 and L6 == L18
        assert (P18 == 3 * P6).all(), "m-scaling violated at p=18"
        print("m-scaling checked table-wide: pure_18 = 3 * pure_6 "
              f"on all {1 << d6} contents (buckets "
              f"{sorted(set(P6[1:].tolist()))} -> "
              f"{sorted(set(P18[1:].tolist()))})", flush=True)
        return


if __name__ == "__main__":
    main()
