#!/usr/bin/env python3
"""A42 S2e — the JET racer: the class-complete register repair.

s2_registers.py proved (mechanically, at a = 1 and a = 2) that the
S1h residue register is blind to all of ker(pi) on H — 15 of 63
classes at a = 1, and at a = 2 it (and any mod-pi^2 functional)
vanishes on EVERY class, because the Tor-representatives at depth a
are pi^{2^a - 2^{order}}-divisible.  The repair, also validated
there: FULL-DEPTH Lambda_a-valued registers, one per variety branch —
the register accumulates  R = sum_c  x*^c * ev(v1-column at c)  in
Lambda_a = F_2[y]/((y^2+y+1)^{2^a}), where x* is an exact Lambda_a
root of the v1-block's boundary multiplier Pbar_omega (computed by
brute force, residue-anchored at omega / omega^2 — the two
non-vacuous branches).  Invariance is structural (x* kills boundary
contents exactly); the +c convention is forced (the -c convention
evaluates at x*^{-1}, a root on the OTHER branch, and loses
invariance — found mechanically in s2_registers).  Per-branch
kernels on H are nonzero but the JOINT kernel is EMPTY (verified at
a = 1, 2): two runs — one per branch — certify all classes.

Engine: the S1h/S2b automaton core (validated against
a40_s4_phase_atlas.Automaton.step), states in uint64 (5p bits, no
register bits), registers in a parallel uint8 (Lambda_a has <= 8
F2-dims for a <= 2), pair dedup/membership by NATIVE u64/u8 lexsort
and sort-merge (the first flight used 9-byte void keys and spent
80% of its time in numpy's per-element VOID_compare — sampled,
~20x slower), insertion-time dedup (the S2b memory fix), Dial
buckets, y-rotation canonicalization with register co-rotation by
ybar^s, RSS-capped clean aborts, one branch per process (ru_maxrss
is a lifetime peak).

Certificate semantics: completing level c in branch b enumerates all
compact cycles of weight <= c with their branch-b register; "no
nonzero-register return <= L" in BOTH branches + joint-kernel
emptiness => no nontrivial compact cycle of weight <= L, ANY class.

Controls: p=3 (union min nontrivial 6), p=6 (union min 12, PLUS
ground truth: every one of the ~315 atlas cycles at W<=12 replayed
through both branch registers must satisfy [embed-nontrivial <=>
some branch register nonzero]), p=9 regression (no nonzero-register
return <= 16 — floor(9) re-certified by the new engine), and random
explicit boundary cycles at p = 12 replay to register 0 (engine-exact
invariance at a = 2).
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

_spec = importlib.util.spec_from_file_location(
    "a40_s4_phase_atlas", Path(__file__).parent / "a40_s4_phase_atlas.py")
AT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AT)

DATA = LAB / "data" / "a42"

RSS_CAP = int(2.2 * 1024 ** 3)
CHUNK = 1024
OMEGA2 = 0b111


class JetRacer:
    """One variety branch's full-depth register racer at period p."""

    def __init__(self, p: int, branch: int):
        assert p % 3 == 0, "registers valid only when 3 | p"
        assert 5 * p <= 64
        self.p = p
        self.branch = branch
        self.a = AL.v2(p)
        # Lambda_a
        m = OMEGA2
        for _ in range(self.a):
            m = AL.pmul(m, m)
        self.Lmod = m
        self.dim = AL.pdeg(m)
        assert self.dim <= 8
        nl = 1 << self.dim

        def red(z):
            return AL.pmod(z, self.Lmod)

        def mul(x_, y_):
            return red(AL.pmul(x_, y_))

        self.red, self.Lmul = red, mul
        # automaton (verbatim S1h constants)
        self.au = AT.Automaton(AT.A_L, AT.B_L, p)
        au = self.au
        assert (au.nf, au.no) == (3, 2) and au.forced_blk == 1
        assert au.adv_f == 3 and au.adv_o == 1 and au.top_j == 1
        self.mask = np.uint64((1 << p) - 1)
        self.pc = np.array([bin(i).count("1") for i in range(1 << p)],
                           dtype=np.int64)
        # multiplication table over Lambda_a (uint8 x uint8)
        MT = np.zeros((nl, nl), dtype=np.uint8)
        for x_ in range(nl):
            for y_ in range(nl):
                MT[x_, y_] = mul(x_, y_)
        self.MT = MT
        # inverses of units
        inv = {}
        for x_ in range(1, nl):
            for y_ in range(1, nl):
                if mul(x_, y_) == 1:
                    inv[x_] = y_
                    break
        # the v1-block boundary multiplier Pbar_omega(x) =
        # 1 + ybar^{-1} + x^{-3} ybar  (bar(A_L) = [(0,0),(0,-1),(-3,1)])
        ybar = red(0b10)
        ybi = inv[ybar]
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
        # non-vacuousness: the v2-multiplier Qbar_omega at xstar has
        # residue 0 (the point IS on the residue variety)
        # Qbar = 1 + x^{-1} + x ybar^3   (bar(B_L) = [(0,0),(1,0),(-1,3)])
        y3 = mul(ybar, mul(ybar, ybar))
        qv = 1 ^ inv[self.xstar] ^ mul(self.xstar, y3)
        assert AL.pmod(qv, OMEGA2) == 0, "xstar not on the residue variety"
        self.XI = inv[self.xstar]
        self.MXI = MT[self.XI]          # row lookup: R -> XI * R
        # column content table ev[col] = col mod Lmod
        self.ev = np.array([red(c) for c in range(1 << p)],
                           dtype=np.uint8)
        # y-rotation co-multipliers ybar^s
        self.yrot = np.zeros(p, dtype=np.uint8)
        acc = 1
        for s in range(p):
            self.yrot[s] = acc
            acc = mul(acc, ybar)
        self.MROT = [MT[self.yrot[s]] for s in range(p)]

    # ---------------- vectorized step ----------------
    def expand(self, states: np.ndarray, regs: np.ndarray):
        p = self.p
        m = self.mask
        f0 = states & m
        f1 = (states >> np.uint64(p)) & m
        f2 = (states >> np.uint64(2 * p)) & m
        o0 = (states >> np.uint64(3 * p)) & m
        o1 = (states >> np.uint64(4 * p)) & m

        def rot(x, s):
            s %= p
            if s == 0:
                return x
            return ((x << np.uint64(s)) | (x >> np.uint64(p - s))) & m

        base = f0 ^ rot(f0, -1 % p) ^ o1 ^ rot(o0, 3)
        nreg_base = self.MXI[regs]      # XI * R  (input-independent)
        outs = []
        for a in range(1 << p):
            av = np.uint64(a)
            acc = base ^ av
            new_f = rot(acc, -1 % p)
            cost = self.pc[a] + self.pc[new_f.astype(np.int64)]
            ns = f1 | (f2 << np.uint64(p)) | (new_f << np.uint64(2 * p)) \
                | (o1 << np.uint64(3 * p)) | (av << np.uint64(4 * p))
            nreg = nreg_base ^ self.ev[a]
            outs.append((ns, nreg, cost))
        return outs

    # ---------------- canonicalization ----------------
    def canon(self, states: np.ndarray, regs: np.ndarray):
        p = self.p
        m = self.mask
        bs = states.copy()
        br = regs.copy()
        shifts = [np.uint64(p * k) for k in range(5)]
        for s in range(1, p):
            rs = np.zeros_like(states)
            for k in range(5):
                c = (states >> shifts[k]) & m
                c = ((c << np.uint64(s)) | (c >> np.uint64(p - s))) & m
                rs |= c << shifts[k]
            rr = self.MROT[s][regs]
            better = (rs < bs) | ((rs == bs) & (rr < br))
            bs = np.where(better, rs, bs)
            br = np.where(better, rr, br)
        return bs, br

    # ---------------- pair keys (native lexsort; the first-flight
    # V9-void-dtype path spent 80% of its time in VOID_compare) ------
    @staticmethod
    def dedup(ms: np.ndarray, mr: np.ndarray):
        """Sort (lex by state, then reg) and drop duplicate pairs."""
        o = np.lexsort((mr, ms))
        ms, mr = ms[o], mr[o]
        if ms.size <= 1:
            return ms, mr
        keep = np.empty(ms.size, dtype=bool)
        keep[0] = True
        keep[1:] = (ms[1:] != ms[:-1]) | (mr[1:] != mr[:-1])
        return ms[keep], mr[keep]

    @staticmethod
    def member(qs, qr, ss, sr):
        """Membership mask of pair-queries (qs, qr) — assumed
        pair-deduped — in the pair-unique sorted store (ss, sr):
        merge, lexsort, mark equal-neighbors.  Since each side is
        duplicate-free, a query equals a neighbor iff that neighbor
        is its store copy."""
        if ss.size == 0:
            return np.zeros(qs.size, dtype=bool)
        ns = ss.size
        all_s = np.concatenate([ss, qs])
        all_r = np.concatenate([sr, qr])
        o = np.lexsort((all_r, all_s))
        eq = (all_s[o][1:] == all_s[o][:-1]) & \
            (all_r[o][1:] == all_r[o][:-1])
        inm = np.zeros(all_s.size, dtype=bool)
        inm[1:] |= eq
        inm[:-1] |= eq
        back = np.empty_like(o)
        back[o] = np.arange(o.size)
        return inm[back[ns:]]

    # ---------------- validation ----------------
    def validate(self, nsamples=500, rng=None):
        rng = rng or np.random.default_rng(11)
        p = self.p
        for _ in range(nsamples):
            cols = [int(rng.integers(0, 1 << p)) for _ in range(5)]
            a = int(rng.integers(0, 1 << p))
            st = tuple(cols)
            st2, cost2, _ = self.au.step(st, a)
            sv = np.uint64(sum(c << (p * k) for k, c in enumerate(cols)))
            outs = self.expand(np.array([sv], dtype=np.uint64),
                               np.zeros(1, dtype=np.uint8))
            ns, nreg, co = outs[a]
            got = [(int(ns[0]) >> (p * k)) & ((1 << p) - 1)
                   for k in range(5)]
            assert tuple(got) == st2
            assert int(co[0]) == cost2
            # register rule: XI*0 + ev[a]
            assert int(nreg[0]) == int(self.ev[a])
        # rotation commutation incl. register co-rotation
        for _ in range(300):
            cols = [int(rng.integers(0, 1 << p)) for _ in range(5)]
            r0 = int(rng.integers(0, 1 << self.dim))
            s = int(rng.integers(1, p))
            sv = np.array([np.uint64(sum(c << (p * k)
                                         for k, c in enumerate(cols)))],
                          dtype=np.uint64)
            rv = np.array([r0], dtype=np.uint8)
            c1s, c1r = self.canon(sv, rv)
            # rotate then canon must agree
            rcols = [AT.rot(c, s, p) for c in cols]
            sv2 = np.array([np.uint64(sum(c << (p * k)
                                          for k, c in enumerate(rcols)))],
                           dtype=np.uint64)
            rv2 = np.array([self.MROT[s][r0]], dtype=np.uint8)
            c2s, c2r = self.canon(sv2, rv2)
            assert int(c1s[0]) == int(c2s[0]) and int(c1r[0]) == int(c2r[0])
        return True

    def replay_register(self, inputs):
        """Scalar replay of the register over an input-column sequence
        (the engine rule): R <- XI*R + ev[a]."""
        r = 0
        for a in inputs:
            r = int(self.MT[self.XI, r]) ^ int(self.ev[a])
        return r

    # ---------------- the BFS ----------------
    def run(self, cap: int, log=print, ckpt_path=None,
            time_budget_s=None, growth_est=3.2):
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
                 "engine": "jet (full-depth Lambda register)",
                 "lane": "AB (theta'-dual lane by banked duality)"}
        aborted = False
        last_level_t = 0.0
        for c in range(0, cap + 1):
            if time_budget_s is not None and c > 0:
                if (time.time() - t0) + last_level_t * growth_est \
                        > time_budget_s:
                    state["aborted"] = f"time budget before level {c}"
                    log(f"  time budget — stop before level {c} "
                        "(completed levels are certified)")
                    aborted = True
                    break
            lt0 = time.time()
            level_novel = 0
            while buckets.get(c):
                arrs = buckets.pop(c)
                sts = np.concatenate([x[0] for x in arrs])
                rgs = np.concatenate([x[1] for x in arrs])
                sts, rgs = self.dedup(sts, rgs)
                fresh = ~self.member(sts, rgs, seen_s, seen_r)
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
                    for a in range(1 << self.p):
                        ns, nreg, co = outs[a]
                        outs[a] = None
                        if c == 0 and a == 0:
                            keep = ns != zstate
                            ns, nreg, co = ns[keep], nreg[keep], co[keep]
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
                        ms, mr = self.dedup(ms, mr)
                        fr = ~self.member(ms, mr, seen_s, seen_r)
                        if fr.any():
                            buckets.setdefault(cw, []).append(
                                (ms[fr], mr[fr]))
                    del outs, per_cost
                    rss = resource.getrusage(
                        resource.RUSAGE_SELF).ru_maxrss
                    if rss > RSS_CAP:
                        state["aborted"] = "RSS"
                        log(f"  RSS cap ({rss//2**20}MB) — clean abort "
                            "(completed levels certified)")
                        aborted = True
                        break
                if aborted:
                    break
                for cw in list(buckets.keys()):
                    arrs2 = buckets[cw]
                    tot = sum(x[0].size for x in arrs2)
                    if len(arrs2) > 64 or tot > 1 << 22:
                        ms = np.concatenate([x[0] for x in arrs2])
                        mr = np.concatenate([x[1] for x in arrs2])
                        ms, mr = self.dedup(ms, mr)
                        fr = ~self.member(ms, mr, seen_s, seen_r)
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
            log(f"  level {c}: novel {level_novel} seen {seen_s.size} "
                f"rss {rss//2**20}MB t {time.time()-t0:.0f}s "
                f"returns {sorted(returns)} nonzero-reg@{nontriv}")
            state["returns"] = {str(k): sorted(v)
                                for k, v in returns.items()}
            state["nonzero_reg_costs"] = nontriv
            if ckpt_path:
                ckpt_path.write_text(json.dumps(state, indent=1))
        state["wall_s"] = round(time.time() - t0, 1)
        if ckpt_path:
            ckpt_path.write_text(json.dumps(state, indent=1))
        return state


# ------------------- ground-truth controls -------------------------
def boundary_cycle_inputs(p: int, rng, nterms=3, span=5):
    """Random explicit boundary (v1, v2) = (Pbar t, Qbar t) as
    column dicts; returns the v1-input sequence over a covering
    range plus the v2 columns for cross-checking."""
    Pb = [(-i, -j) for (i, j) in AT.A_L]
    Qb = [(-i, -j) for (i, j) in AT.B_L]
    v1, v2 = {}, {}
    for _ in range(nterms):
        tc = int(rng.integers(0, span))
        tj = int(rng.integers(0, p))
        for (dx, dy) in Pb:
            c = tc + dx
            v1[c] = v1.get(c, 0) ^ (1 << ((tj + dy) % p))
        for (dx, dy) in Qb:
            c = tc + dx
            v2[c] = v2.get(c, 0) ^ (1 << ((tj + dy) % p))
    v1 = {c: b for c, b in v1.items() if b}
    v2 = {c: b for c, b in v2.items() if b}
    return v1, v2


def replay_cycle_through(au, r: JetRacer, v1, v2):
    """March the automaton over the cycle's own columns, accumulating
    the jet register; assert forced columns match; return final reg."""
    cs = list(v1) + list(v2)
    lo, hi = min(cs) - 6, max(cs) + 6
    st = au.zero()
    reg = 0
    for g in range(lo, hi + 1):
        a = v1.get(g + 1, 0)
        st2, _, newf = au.step(st, a)
        assert newf == v2.get(g + 3, 0)
        reg = int(r.MT[r.XI, reg]) ^ int(r.ev[a])
        st = st2
    assert st == au.zero()
    return reg


def controls():
    rng = np.random.default_rng(23)
    # p=6 ground truth: all atlas cycles, embed-classified
    p = 6
    rT = JetRacer(p, 0)
    rV = JetRacer(p, 1)
    rT.validate()
    rV.validate()
    print(f"p=6 jet racers validated (xstar T={rT.xstar} V={rV.xstar})",
          flush=True)
    au = AT.Automaton(AT.A_L, AT.B_L, p)
    rows, _ = AT.atlas("AB", p, 12, keep_pts=True)
    nT = nV = nboth = nnt = 0
    for row in rows:
        v1, v2 = {}, {}
        for (c, y, blk) in row["pts"]:
            d = v1 if blk == 0 else v2
            d[c] = d.get(c, 0) | (1 << y)
        gT = replay_cycle_through(au, rT, v1, v2)
        gV = replay_cycle_through(au, rV, v1, v2)
        seen_nz = (gT != 0) or (gV != 0)
        assert seen_nz == row["nontrivial"], \
            (row["weight"], row["nontrivial"], gT, gV)
        if row["nontrivial"]:
            nnt += 1
            nT += gT != 0
            nV += gV != 0
            nboth += (gT != 0) and (gV != 0)
    print(f"p=6 GROUND TRUTH: {len(rows)} atlas cycles, "
          f"{nnt} nontrivial — [nontrivial <=> some jet register "
          f"nonzero] holds on every cycle (T sees {nT}, V sees {nV}, "
          f"both {nboth})", flush=True)

    # explicit boundary invariance at p = 12 (a = 2, engine-exact)
    p = 12
    rT12 = JetRacer(p, 0)
    rV12 = JetRacer(p, 1)
    rT12.validate()
    rV12.validate()
    au12 = AT.Automaton(AT.A_L, AT.B_L, p)
    for _ in range(60):
        v1, v2 = boundary_cycle_inputs(p, rng)
        if not v1 and not v2:
            continue
        gT = replay_cycle_through(au12, rT12, v1, v2)
        gV = replay_cycle_through(au12, rV12, v1, v2)
        assert gT == 0 and gV == 0, "boundary with nonzero jet register"
    print("p=12 invariance: 60 random explicit boundaries replay to "
          "register 0 in both branches (a=2 engine-exact)", flush=True)
    return rT12, rV12


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    if which in ("all", "controls"):
        # p=3: union min nontrivial 6
        mins = []
        for br in (0, 1):
            r = JetRacer(3, br)
            r.validate()
            st = r.run(8, log=lambda s: None)
            nt = st.get("nonzero_reg_costs", [])
            mins.append(nt[0] if nt else None)
        print(f"p=3 cap 8: branch minima {mins}; union min = "
              f"{min(x for x in mins if x is not None)} (expect 6)",
              flush=True)
        assert min(x for x in mins if x is not None) == 6
        # p=6 racer-level: union min 12
        mins6 = []
        for br in (0, 1):
            r = JetRacer(6, br)
            st = r.run(13, log=lambda s: None)
            nt = st.get("nonzero_reg_costs", [])
            mins6.append(nt[0] if nt else None)
        print(f"p=6 cap 13: branch minima {mins6}; union min = "
              f"{min(x for x in mins6 if x is not None)} (expect 12)",
              flush=True)
        assert min(x for x in mins6 if x is not None) == 12
        # p=9 regression: floor(9) by the new engine
        for br in (0, 1):
            r = JetRacer(9, br)
            st = r.run(16, log=lambda s: print(s, flush=True))
            assert not st.get("nonzero_reg_costs"), st
            assert max(lv["c"] for lv in st["levels"]) == 16
        print("p=9 cap 16 both branches: no nonzero-register return — "
              "floor(9) = 18 RE-CERTIFIED by the jet engine (a=0)",
              flush=True)
        controls()

    if which == "branch":
        # one branch per PROCESS: ru_maxrss is a lifetime peak, so
        # sequential branches in one process would inherit the first
        # branch's peak and abort instantly (the s2_jet12 first-launch
        # failure mode).
        tag = sys.argv[2]
        cap = int(sys.argv[3])
        budget = float(sys.argv[4]) * 3600
        br = {"T": 0, "V": 1}[tag]
        r = JetRacer(12, br)
        r.validate()
        print(f"p=12 branch {tag} (xstar={r.xstar}): cap {cap}, "
              f"budget {budget/3600:.1f} h", flush=True)
        st = r.run(cap, log=lambda s: print(s, flush=True),
                   ckpt_path=DATA / f"s2_jet12_{tag}_ckpt.json",
                   time_budget_s=budget)
        done = max((lv["c"] for lv in st["levels"]), default=-1)
        nt = st.get("nonzero_reg_costs", [])
        print(f"p=12 branch {tag}: completed level {done}; "
              f"nonzero-register costs {nt}", flush=True)
        (DATA / f"s2_jet12_{tag}.json").write_text(json.dumps(st,
                                                              indent=1))
        return

    if which == "summary":
        out = {}
        for tag in ("T", "V"):
            out[tag] = json.loads(
                (DATA / f"s2_jet12_{tag}.json").read_text())
        done = min(max((lv["c"] for lv in out[t]["levels"]),
                       default=-1) for t in ("T", "V"))
        ntU = sorted(set(out["T"].get("nonzero_reg_costs", []))
                     | set(out["V"].get("nonzero_reg_costs", [])))
        ev = done if done % 2 == 0 else done - 1
        print(f"branches complete through joint level {done}; "
              f"union nonzero-register costs {ntU}", flush=True)
        if not ntU and done >= 0:
            print(f">>> no nonzero-register return <= {done} in either "
                  f"branch: joint-kernel emptiness (s2_registers) => "
                  f"no nontrivial compact cycle of weight <= {done}; "
                  f"parity => floor(12) >= {ev + 2} unrestricted, ALL "
                  "classes (certificate tier, lane AB + banked "
                  "duality)", flush=True)
        (DATA / "s2_jet12.json").write_text(json.dumps(out, indent=1))
        print(f"wrote {DATA/'s2_jet12.json'}", flush=True)


if __name__ == "__main__":
    main()
