#!/usr/bin/env python3
"""A42 S3e — the two-word omega racer: the pure half at p = 24
(a = 3, the r = 4 member period — the next m = 1 rung).

Same object as a42_s3b_omega_racer (pure-lift cycles enumerated by
weight on the omega-quotient, branch registers, rotation canon) but
with 5 x dim(Lambda_a)-bit states held as TWO words (lo: columns
0..3, hi: column 4) so dim 16 fits, fixed-constant multiplications
as lo/hi byte-split tables, and INPUT-MAJOR expansion (per state,
one 2^dim-gather across all inputs) — the right shape when inputs
outnumber states.  Registers are uint16.

Soundness inputs, all banked: the a = 3 branch roots and the EMPTY
joint register kernel (s3d_registers_a3), Theorem H's dim-6 class
inventory at a = 3, the parity lemma, and the pure cost table via
the CRT idempotent (buckets asserted = 2*2^popcount(nu) on
pi^nu-units — the CMSS digit shadow).

Controls: (i) vector step vs scalar reference; (ii) p = 12
REGRESSION: this engine must reproduce the s3b one-word racer's
level table and returns EXACTLY (both branches, cap 26); (iii) the
p = 24 positive control: the pure lift of pi^7 sigma_0 realizes
weight 48 = 2p, so a completed level-48 run must return nonzero
registers at exactly 48 and nothing below (floor claim: no
nontrivial pure compact cycle of weight < 48).
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


OM = _load("a42_s3b_omega_racer")

DATA = LAB / "data" / "a42"
RSS_CAP = int(2.4 * 1024 ** 3)
OMEGA2 = 0b111


def pinv_mod(x: int, mod: int):
    """Inverse of x in F2[y]/(mod) by extended Euclid; None if not a
    unit."""
    r0, r1 = mod, AL.pmod(x, mod)
    s0, s1 = 0, 1
    while r1:
        q, r = AL.pdivmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, s0 ^ AL.pmul(q, s1)
    if r0 != 1:
        return None
    return AL.pmod(s0, mod)


def lin_tables(dim, mulK):
    """lo/hi byte-split tables of the F2-linear map v -> mulK(v)."""
    lo = np.array([mulK(v) for v in range(256)], dtype=np.uint16)
    if dim > 8:
        hi = np.array([mulK(v << 8) for v in range(256)],
                      dtype=np.uint16)
    else:
        hi = np.zeros(256, dtype=np.uint16)
    return lo, hi


def apply_lin(tabs, v):
    lo, hi = tabs
    return lo[v & 255] ^ hi[v >> 8]


class OmegaRacer2W:
    """Two-word pure-half omega racer at period p = 3*2^a (m = 1
    only here), one branch per instance."""

    def __init__(self, p: int, branch: int):
        assert p % 3 == 0 and p == 3 * (1 << AL.v2(p)), "m = 1 only"
        self.p = p
        self.branch = branch
        self.a = AL.v2(p)
        self.P, self.dim, self.Lmod = OM.crt_pure_table(p)
        self.P = self.P.astype(np.int64)
        d = self.dim
        nl = 1 << d
        assert d <= 16

        def red(z):
            return AL.pmod(z, self.Lmod)

        def mul(x_, y_):
            return red(AL.pmul(x_, y_))

        self.red, self.Lmul = red, mul
        # S4 fix: the session-3 draft found inverses by an O(nl^2)
        # scan (4e9 pmul calls at dim 16 — never finishes at p = 24);
        # extended Euclid over F2[y] instead, asserted per unit.
        inv = {}
        for x_ in range(1, nl):
            u_ = pinv_mod(x_, self.Lmod)
            if u_ is not None:
                assert mul(x_, u_) == 1
                inv[x_] = u_
        ybar = red(0b10)
        ybi = inv[ybar]
        roots = []
        for x_ in range(1, nl):
            if x_ not in inv:
                continue
            xi3 = mul(inv[x_], mul(inv[x_], inv[x_]))
            if 1 ^ ybi ^ mul(xi3, ybar) == 0:
                resid = AL.pmod(x_, OMEGA2)
                if resid != 1:
                    roots.append((resid, x_))
        roots.sort()
        assert len(roots) == 2, roots
        self.xstar = roots[branch][1]
        y3 = mul(ybar, mul(ybar, ybar))
        assert AL.pmod(1 ^ inv[self.xstar]
                       ^ mul(self.xstar, y3), OMEGA2) == 0
        XI = inv[self.xstar]
        self.T_XI = lin_tables(d, lambda v: mul(XI, v))
        self.T_F0 = lin_tables(d, lambda v: mul(1 ^ ybi, v))
        self.T_O0 = lin_tables(d, lambda v: mul(y3, v))
        self.T_TOP = lin_tables(d, lambda v: mul(ybi, v))
        # full TOP table over all inputs (for the input-major gather)
        self.TOP_ALL = np.array([mul(ybi, v) for v in range(nl)],
                                dtype=np.uint16)
        # rotation constants ybar^s
        self.rot_tabs = []
        acc = 1
        for s in range(p):
            self.rot_tabs.append(lin_tables(
                d, (lambda k: (lambda v: mul(k, v)))(acc)))
            acc = mul(acc, ybar)
        assert mul(acc, 1) == 1, "ybar^p != 1"
        self.ALL_A = np.arange(nl, dtype=np.int64)
        self.COST_A = self.P[self.ALL_A]

    # state = (lo: u64 = cols 0..3, hi: u16 = col 4); reg u16
    def cols_of(self, lo, hi):
        d = self.dim
        mask = np.uint64((1 << d) - 1)
        c = [((lo >> np.uint64(d * k)) & mask).astype(np.int64)
             for k in range(4)]
        c.append(hi.astype(np.int64))
        return c

    def pack_cols(self, cols):
        d = self.dim
        lo = np.uint64(0) if np.isscalar(cols[0]) else \
            np.zeros(len(cols[0]), dtype=np.uint64)
        for k in range(4):
            lo = lo | (np.asarray(cols[k]).astype(np.uint64)
                       << np.uint64(d * k))
        hi = np.asarray(cols[4]).astype(np.uint16)
        return lo, hi

    def step_ref(self, cols, aa):
        f0, f1, f2, o0, o1 = cols
        acc = int(apply_lin(self.T_F0, np.array([f0]))[0]) \
            ^ int(apply_lin(self.T_O0, np.array([o0]))[0]) ^ o1 ^ aa
        new_f = int(self.TOP_ALL[acc])
        cost = int(self.P[aa]) + int(self.P[new_f])
        return (f1, f2, new_f, o1, aa), cost, new_f

    def expand_state(self, cols, reg):
        """Input-major: all 2^dim transitions of ONE state.
        Returns (lo_arr, hi_arr, reg_arr, cost_arr, newf_arr)."""
        d = self.dim
        f0, f1, f2, o0, o1 = [int(x) for x in cols]
        base = int(apply_lin(self.T_F0, np.array([f0]))[0]) \
            ^ int(apply_lin(self.T_O0, np.array([o0]))[0]) ^ o1
        new_f = self.TOP_ALL[self.ALL_A ^ base]
        cost = self.COST_A + self.P[new_f.astype(np.int64)]
        # build lo vectorized: cols (f1, f2, new_f, o1)
        lo = (np.uint64(f1)
              | (np.uint64(f2) << np.uint64(d))
              | (new_f.astype(np.uint64) << np.uint64(2 * d))
              | (np.uint64(o1) << np.uint64(3 * d)))
        hi = self.ALL_A.astype(np.uint16)
        nreg = (apply_lin(self.T_XI, np.array([reg],
                                              dtype=np.uint16))[0]
                ^ self.ALL_A).astype(np.uint16)
        return lo, hi, nreg, cost, new_f

    def canon(self, lo, hi, reg):
        d = self.dim
        cols = self.cols_of(lo, hi)
        best_lo, best_hi, best_r = lo.copy(), hi.copy(), reg.copy()
        for s in range(1, self.p):
            tabs = self.rot_tabs[s]
            rc = [apply_lin(tabs, c.astype(np.uint16))
                  for c in cols]
            rlo, rhi = self.pack_cols(rc)
            rr = apply_lin(tabs, reg)
            better = (rlo < best_lo) | \
                ((rlo == best_lo) & (rhi < best_hi)) | \
                ((rlo == best_lo) & (rhi == best_hi) & (rr < best_r))
            best_lo = np.where(better, rlo, best_lo)
            best_hi = np.where(better, rhi, best_hi)
            best_r = np.where(better, rr, best_r)
        return best_lo, best_hi, best_r

    @staticmethod
    def dedup3(lo, hi, rg):
        o = np.lexsort((rg, hi, lo))
        lo, hi, rg = lo[o], hi[o], rg[o]
        if lo.size <= 1:
            return lo, hi, rg
        keep = np.empty(lo.size, dtype=bool)
        keep[0] = True
        keep[1:] = (lo[1:] != lo[:-1]) | (hi[1:] != hi[:-1]) \
            | (rg[1:] != rg[:-1])
        return lo[keep], hi[keep], rg[keep]

    @staticmethod
    def member3(qlo, qhi, qrg, slo, shi, srg):
        """Membership of triples in the lex-sorted store: searchsorted
        on lo, then scan hi/rg within equal-lo runs."""
        if slo.size == 0:
            return np.zeros(qlo.size, dtype=bool)
        left = np.searchsorted(slo, qlo, side="left")
        right = np.searchsorted(slo, qlo, side="right")
        res = np.zeros(qlo.size, dtype=bool)
        width = right - left
        wmax = int(width.max()) if width.size else 0
        for dd in range(wmax):
            sel = width > dd
            idx = left[sel] + dd
            res[sel] |= (shi[idx] == qhi[sel]) & (srg[idx] == qrg[sel])
        return res

    def validate(self, nsamples=400, rng=None):
        rng = rng or np.random.default_rng(9)
        d = self.dim
        nl = 1 << d
        for _ in range(nsamples):
            cols = [int(rng.integers(0, nl)) for _ in range(5)]
            aa = int(rng.integers(0, nl))
            r0 = int(rng.integers(0, nl))
            lo, hi, nreg, cost, new_f = self.expand_state(
                cols, np.uint16(r0))
            st2, cost2, nf2 = self.step_ref(cols, aa)
            got_cols = self.cols_of(np.array([lo[aa]]),
                                    np.array([hi[aa]]))
            got = tuple(int(c[0]) for c in got_cols)
            assert got == st2, (got, st2)
            assert int(cost[aa]) == cost2
            # register rule: XI*r + a
            xr = int(apply_lin(self.T_XI,
                               np.array([r0], dtype=np.uint16))[0])
            assert int(nreg[aa]) == (xr ^ aa)
        # canon commutation on singletons
        for _ in range(150):
            cols = [int(rng.integers(0, nl)) for _ in range(5)]
            r0 = int(rng.integers(0, nl))
            s = int(rng.integers(1, self.p))
            lo, hi = self.pack_cols([np.array([c]) for c in cols])
            rg = np.array([r0], dtype=np.uint16)
            c1 = self.canon(lo, hi, rg)
            tabs = self.rot_tabs[s]
            rc = [int(apply_lin(tabs,
                                np.array([c], dtype=np.uint16))[0])
                  for c in cols]
            lo2, hi2 = self.pack_cols([np.array([c]) for c in rc])
            rg2 = np.array([int(apply_lin(
                tabs, np.array([r0], dtype=np.uint16))[0])],
                dtype=np.uint16)
            c2 = self.canon(lo2, hi2, rg2)
            assert all(int(x[0]) == int(y[0])
                       for x, y in zip(c1, c2))
        return True

    def run(self, cap: int, log=print, ckpt_path=None,
            time_budget_s=None, growth_est=2.5):
        t0 = time.time()
        buckets: dict[int, list] = {
            0: [(np.zeros(1, dtype=np.uint64),
                 np.zeros(1, dtype=np.uint16),
                 np.zeros(1, dtype=np.uint16))]}
        slo = np.zeros(1, dtype=np.uint64)
        shi = np.zeros(1, dtype=np.uint16)
        srg = np.zeros(1, dtype=np.uint16)
        returns: dict[int, set] = {}
        state = {"p": self.p, "cap": cap, "branch": self.branch,
                 "xstar": int(self.xstar), "levels": [],
                 "engine": "two-word omega pure racer",
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
                lo = np.concatenate([x[0] for x in arrs])
                hi = np.concatenate([x[1] for x in arrs])
                rg = np.concatenate([x[2] for x in arrs])
                lo, hi, rg = self.dedup3(lo, hi, rg)
                fresh = ~self.member3(lo, hi, rg, slo, shi, srg)
                if c == 0 and level_novel == 0:
                    fresh[:] = True
                lo, hi, rg = lo[fresh], hi[fresh], rg[fresh]
                if lo.size == 0:
                    break
                level_novel += int(lo.size)
                slo = np.concatenate([slo, lo])
                shi = np.concatenate([shi, hi])
                srg = np.concatenate([srg, rg])
                o = np.lexsort((srg, shi, slo))
                slo, shi, srg = slo[o], shi[o], srg[o]
                per_cost: dict[int, list] = {}
                for i in range(lo.size):
                    cols = [int(x[0]) for x in
                            self.cols_of(lo[i:i + 1], hi[i:i + 1])]
                    zlo, zhi, zrg, zco, znf = self.expand_state(
                        cols, rg[i])
                    if c == 0 and cols == [0] * 5:
                        keep = ~((zlo == 0) & (zhi == 0))
                        zlo, zhi, zrg, zco = zlo[keep], zhi[keep], \
                            zrg[keep], zco[keep]
                    cw_all = c + zco
                    ok = cw_all <= cap
                    zsel = ok & (zlo == 0) & (zhi == 0)
                    if zsel.any():
                        for cw, rv in zip(cw_all[zsel], zrg[zsel]):
                            if cw > 0:
                                returns.setdefault(int(cw), set()) \
                                    .add(int(rv))
                    for w in np.unique(cw_all[ok]):
                        sel = cw_all == w
                        per_cost.setdefault(int(w), []).append(
                            (zlo[sel], zhi[sel],
                             zrg[sel].astype(np.uint16)))
                    # S4: flush every 64 states (1024 x 2^16-wide
                    # expansions was ~1.3 GB at dim 16 — the level-13
                    # RSS abort of the first p = 24 flight)
                    if (i & 63) == 63 or i == lo.size - 1:
                        for cw, lst in per_cost.items():
                            mlo = np.concatenate([x[0] for x in lst])
                            mhi = np.concatenate([x[1] for x in lst])
                            mrg = np.concatenate([x[2] for x in lst])
                            mlo, mhi, mrg = self.canon(mlo, mhi, mrg)
                            mlo, mhi, mrg = self.dedup3(mlo, mhi,
                                                        mrg)
                            fr = ~self.member3(mlo, mhi, mrg,
                                               slo, shi, srg)
                            if fr.any():
                                buckets.setdefault(cw, []).append(
                                    (mlo[fr], mhi[fr], mrg[fr]))
                        per_cost = {}
                        rss = resource.getrusage(
                            resource.RUSAGE_SELF).ru_maxrss
                        if rss > RSS_CAP:
                            state["aborted"] = "RSS"
                            log("  RSS cap — clean abort")
                            aborted = True
                            break
                if aborted:
                    break
                for cw in list(buckets.keys()):
                    arrs2 = buckets[cw]
                    tot = sum(x[0].size for x in arrs2)
                    if len(arrs2) > 16 or tot > 1 << 20:
                        mlo = np.concatenate([x[0] for x in arrs2])
                        mhi = np.concatenate([x[1] for x in arrs2])
                        mrg = np.concatenate([x[2] for x in arrs2])
                        mlo, mhi, mrg = self.dedup3(mlo, mhi, mrg)
                        fr = ~self.member3(mlo, mhi, mrg,
                                           slo, shi, srg)
                        buckets[cw] = [(mlo[fr], mhi[fr],
                                        mrg[fr])] if fr.any() else []
            if aborted:
                break
            last_level_t = time.time() - lt0
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            state["levels"].append(
                {"c": c, "novel": level_novel,
                 "seen": int(slo.size), "rss_mb": rss // 2 ** 20,
                 "t": round(time.time() - t0, 1)})
            nontriv = sorted(w for w, rv in returns.items()
                             if any(v != 0 for v in rv))
            log(f"  level {c}: novel {level_novel} seen {slo.size} "
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


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "regress12"

    if which == "regress12":
        # the decisive control: reproduce s3b at p = 12 exactly
        for br, tag in ((0, "T"), (1, "V")):
            r = OmegaRacer2W(12, br)
            r.validate()
            st = r.run(26, log=lambda s: None)
            ref = json.loads(
                (DATA / f"s3b_omega12_{tag}.json").read_text())
            assert st["returns"] == ref["returns"], \
                (tag, st["returns"], ref["returns"])
            lv1 = [(x["c"], x["novel"]) for x in st["levels"]]
            lv2 = [(x["c"], x["novel"]) for x in ref["levels"]]
            assert lv1 == lv2, (tag, "level drift")
            print(f"p=12 branch {tag}: two-word engine reproduces "
                  "the s3b level table and returns EXACTLY "
                  f"(nonzero@{st['nonzero_reg_costs']})", flush=True)
        return

    if which == "run24":
        tag = sys.argv[2]
        cap = int(sys.argv[3])
        budget = float(sys.argv[4]) * 3600
        br = {"T": 0, "V": 1}[tag]
        r = OmegaRacer2W(24, br)
        r.validate()
        # bucket sanity (S4 correction: the session-3 draft asserted
        # pure(pi^nu) itself = 2*2^popcount(nu), but pi = 1+y+y^2 has
        # y-weight 3 (pure 6) — the digit-shadow law is the MINIMUM
        # over the valuation bucket, attained e.g. by (1+t)^nu with
        # t = y^{1-q}; assert the ideal minima instead)
        pi = r.red(OMEGA2)
        nl_ = 1 << r.dim
        ideals = []
        z = 1
        for nu in range(9):
            ideals.append({r.Lmul(z, u) for u in range(nl_)})
            z = r.Lmul(z, pi)
        for nu in range(8):
            stratum = ideals[nu] - ideals[nu + 1]   # exact valuation
            mn = min(int(r.P[w]) for w in stratum)
            assert mn == 2 * (1 << bin(nu).count("1")), (nu, mn)
        print(f"p=24 pure racer branch {tag} (xstar={r.xstar}): "
              f"cap {cap}, budget {budget/3600:.1f} h; digit-shadow "
              "buckets asserted", flush=True)
        st = r.run(cap, log=lambda s: print(s, flush=True),
                   ckpt_path=DATA / f"s3e_omega24_{tag}_ckpt.json",
                   time_budget_s=budget)
        done = max((lv["c"] for lv in st["levels"]), default=-1)
        nt = st.get("nonzero_reg_costs", [])
        print(f"p=24 branch {tag}: completed level {done}; "
              f"nonzero-register costs {nt}", flush=True)
        (DATA / f"s3e_omega24_{tag}.json").write_text(
            json.dumps(st, indent=1))
        return


if __name__ == "__main__":
    main()
