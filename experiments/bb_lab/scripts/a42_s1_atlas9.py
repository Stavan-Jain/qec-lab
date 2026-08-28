#!/usr/bin/env python3
"""A42 S1h — the certified numpy racer: compact-cycle floor at p = 9
(and controls) by class-augmented min-cost BFS, no SAT, no parents.

Engine = the S4 atlas automaton (verbatim constants, validated
against `a40_s4_phase_atlas.Automaton.step` on random samples), with:
  * states bit-packed into uint64 ((nf+no) columns x p bits, p <= 11),
  * Dial-bucket min-cost BFS, np.unique dedup, sorted `seen` array,
  * NO parent tracking — instead the state carries a 4-bit omega-class
    accumulator (two F4 evaluation registers, one per x-point of
    gcd(Abar_w, Bbar_w)), updated linearly each step with a fixed
    x0^{-1} rescale so it is march-position independent.  A
    zero-STATE return with nonzero accumulator IS a nontrivial
    compact cycle (the accumulator equals a unit multiple of the
    syzygy parameter's evaluation h(x0); trivial cycles evaluate to
    0 at the gcd's roots).  Soundness of the functional is the exact
    §2.1/§2.2 theory; its implementation is validated on controls.
  * y-rotation canonicalization (quotient by the global y-rotation
    symmetry; accumulators co-rotate by omega^s) — verified to
    commute mechanically,
  * per-level RSS budget (2 GB) with clean abort; each completed
    level c certifies "all compact cycles of weight <= c enumerated".

Controls: p=3 cap 8 (expect nontrivial min 6), p=5 cap 11 (barren:
all returns trivial), p=6 cap 13 (expect min 12).  Production:
p=9 cap 16 — combined with the parity lemma and the realized
weight-18 object, "no nontrivial return <= 16" closes
floor(9) = 18 = 2p at certificate tier.  Lane: AB; the second lane
is its theta'-image (banked duality, †assumption noted in the JSON).
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

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a40_s4_phase_atlas", Path(__file__).parent / "a40_s4_phase_atlas.py")
AT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AT)

DATA = LAB / "data" / "a42"

RSS_CAP = 2 * 1024 ** 3

# ---- F4 arithmetic on 2-bit values: 0,1,2=w,3=w^2 -------------------
F4_MUL = np.zeros((4, 4), dtype=np.uint8)
for i in range(4):
    for j in range(4):
        if i == 0 or j == 0:
            F4_MUL[i, j] = 0
        else:
            # exponents: 1 -> 0, 2 -> 1, 3 -> 2  (w^k)
            ei = {1: 0, 2: 1, 3: 2}[i]
            ej = {1: 0, 2: 1, 3: 2}[j]
            F4_MUL[i, j] = {0: 1, 1: 2, 2: 3}[(ei + ej) % 3]


def f4_add(a, b):
    return a ^ b


class Racer:
    def __init__(self, p: int):
        assert p * 5 + 4 <= 64, "state too wide"
        self.p = p
        # the omega-class registers are eigen-evaluations at the F4
        # y-point: valid exactly when 3 | p (the switch-on condition).
        # For barren p the classification is Theorem A (all trivial),
        # not a computation.
        self.classify = (p % 3 == 0)
        self.au = AT.Automaton(AT.A_L, AT.B_L, p)
        au = self.au
        assert (au.nf, au.no) == (3, 2), (au.nf, au.no)
        assert au.forced_blk == 1 and au.adv_f == 3 and au.adv_o == 1
        assert au.top_j == 1
        assert sorted(au.terms_f) == [(0, -1), (0, 0)]
        assert sorted(au.terms_o) == [(-1, 3), (0, 0), (1, 0)]
        self.mask = np.uint64((1 << p) - 1)
        # column popcount table (p-bit)
        self.pc = np.array([bin(i).count("1") for i in range(1 << p)],
                           dtype=np.int64)
        # omega-evals of p-bit columns: e(col) = sum col_y w^{y mod 3}
        ev = np.zeros(1 << p, dtype=np.uint8)
        wpow = [1, 2, 3]  # w^0, w^1, w^2 in 2-bit coding
        for col in range(1 << p):
            acc = 0
            for y in range(p):
                if (col >> y) & 1:
                    acc ^= wpow[y % 3]
            ev[col] = acc
        self.ev = ev
        # the two x-points and per-point functional data
        self._functional_constants()

    # ---------------- functional constants ----------------
    def _functional_constants(self):
        """The cycle equation is Qbar v1 + Pbar v2 = 0 with (P, Q) =
        (A, B).  omega-parts: v1 = (Pbar_w/g)h, v2 = (Qbar_w/g)h with
        g = gcd.  Registers accumulate, for each root x0 of g, the
        value  sum_c x0^{-c} * e(column emitted at step c, chosen
        block), which equals (unit) * h~(x0) — nonzero for some x0
        iff the cycle class is nontrivial.

        Pbar_w(x) (v2-multiplier) and Qbar_w(x) (v1-multiplier) at
        y -> omega:  A = 1 + y + x^3 y^-1 -> Abar = 1 + y^-1 + x^-3 y:
        Abar_w(x) = 1 + w^2 + x^-3 w = w + w x^-3 = w(1 + x^-3).
        B = 1 + x + x^-1 y^-3 -> Bbar = 1 + x^-1 + x y^3:
        Bbar_w(x) = 1 + x^-1 + x  (y^3 -> 1).
        Up to units: Abar_w ~ x^3 + 1, Bbar_w ~ x^2 + x + 1.
        gcd = x^2 + x + 1 = g; cofactors: co(v2-block) =
        Abar_w/g ~ (x+1), co(v1-block) = Bbar_w/g ~ unit.
        At both roots x0 in {w, w^2}: v1's cofactor is a nonzero
        constant, so the v1-block (the INPUT column block) alone
        detects h at both points.  (v2's cofactor (x+1) is also
        nonzero at both roots — either block works; we use v1 = the
        input columns, emitted at x-position c + adv_o.)"""
        # v1 columns are the INPUT a-columns iff forced_blk == 1
        assert self.au.forced_blk == 1
        # register update per step (march c -> c+1), for x0:
        #   R <- x0^{-1} * R  (+)  x0^{-(adv_o mod 3)} * e(a)
        # (global x0^{c}-rescale; only nonzero-ness matters).
        # x0 in {w, w^2} coded 2, 3; x0^{-1} = x0^2.
        self.x0 = [2, 3]
        self.x0inv = [F4_MUL[x, x] for x in self.x0]
        # prefactor x0^{-adv_o} with adv_o = 1: = x0^{-1} = x0inv
        self.pref = [int(self.x0inv[0]), int(self.x0inv[1])]

    # ---------------- packing ----------------
    # layout: [f0 f1 f2 o0 o1 | acc0(2b) acc1(2b)]
    def pack(self, cols, acc0=0, acc1=0):
        s = 0
        for k, c in enumerate(cols):
            s |= int(c) << (self.p * k)
        s |= (int(acc0) | (int(acc1) << 2)) << (self.p * 5)
        return np.uint64(s)

    def unpack_state(self, s):
        p = self.p
        s = int(s)
        cols = [(s >> (p * k)) & ((1 << p) - 1) for k in range(5)]
        accb = s >> (p * 5)
        return cols, accb & 3, (accb >> 2) & 3

    # ---------------- vectorized step ----------------
    def expand(self, states: np.ndarray):
        """All transitions from `states` (uint64): returns
        (new_states, costs) as concatenated arrays over the 2^p
        inputs."""
        p = self.p
        m = self.mask
        pshift = np.uint64(p)
        f0 = states & m
        f1 = (states >> np.uint64(p)) & m
        f2 = (states >> np.uint64(2 * p)) & m
        o0 = (states >> np.uint64(3 * p)) & m
        o1 = (states >> np.uint64(4 * p)) & m
        accb = states >> np.uint64(5 * p)
        a0 = (accb & np.uint64(3)).astype(np.uint8)
        a1 = ((accb >> np.uint64(2)) & np.uint64(3)).astype(np.uint8)

        def rot(x, s):
            s %= p
            if s == 0:
                return x
            return ((x << np.uint64(s)) | (x >> np.uint64(p - s))) & m

        # acc-base: terms_f: f-col idx0 with j in {0,-1}; terms_o with
        # oo = (o0, o1, a): (a_=-1 -> oo[0]=o0, j=3), (0 -> o1, 0),
        # (1 -> a, 0)
        base = f0 ^ rot(f0, -1 % p) ^ o1 ^ rot(o0, 3)
        # per-step register rescale (input-independent part):
        na0 = F4_MUL[self.x0inv[0], a0]
        na1 = F4_MUL[self.x0inv[1], a1]

        outs = []
        costs = []
        n = states.shape[0]
        for a in range(1 << p):
            av = np.uint64(a)
            acc = base ^ av
            new_f = rot(acc, -1 % p)  # rot by -top_j
            cost = self.pc[a] + self.pc[new_f.astype(np.int64)]
            # new packed state: (f1, f2, new_f, o1, a)
            ns = f1 | (f2 << pshift) | (new_f << np.uint64(2 * p)) \
                | (o1 << np.uint64(3 * p)) | (av << np.uint64(4 * p))
            if self.classify:
                ea = self.ev[a]
                r0 = na0 ^ (F4_MUL[self.pref[0], ea] if ea else 0)
                r1 = na1 ^ (F4_MUL[self.pref[1], ea] if ea else 0)
                ns |= (r0.astype(np.uint64)
                       | (r1.astype(np.uint64) << np.uint64(2))) \
                    << np.uint64(5 * p)
            outs.append(ns)
            costs.append(cost)
        return outs, costs

    # ---------------- canonicalization ----------------
    def canon(self, states: np.ndarray) -> np.ndarray:
        """Min over the p global y-rotations (columns rotate together;
        accumulators multiply by w^{s mod 3})."""
        p = self.p
        m = self.mask
        best = states.copy()
        cols_shift = [np.uint64(p * k) for k in range(5)]
        for s in range(1, p):
            rs = np.zeros_like(states)
            for k in range(5):
                c = (states >> cols_shift[k]) & m
                c = ((c << np.uint64(s)) | (c >> np.uint64(p - s))) & m
                rs |= c << cols_shift[k]
            accb = states >> np.uint64(5 * p)
            a0 = (accb & np.uint64(3)).astype(np.uint8)
            a1 = ((accb >> np.uint64(2)) & np.uint64(3)).astype(np.uint8)
            if self.classify:
                wp = {0: 1, 1: 2, 2: 3}[s % 3]
                a0 = F4_MUL[wp, a0]
                a1 = F4_MUL[wp, a1]
                rs |= (a0.astype(np.uint64)
                       | (a1.astype(np.uint64) << np.uint64(2))) \
                    << np.uint64(5 * p)
            best = np.minimum(best, rs)
        return best

    # ---------------- validation ----------------
    def validate(self, nsamples=3000, rng=None):
        rng = rng or np.random.default_rng(7)
        p = self.p
        for _ in range(nsamples):
            cols = [int(rng.integers(0, 1 << p)) for _ in range(5)]
            a = int(rng.integers(0, 1 << p))
            st = tuple(cols)
            st2, cost2, newf = self.au.step(st, a)
            sv = np.array([self.pack(cols)], dtype=np.uint64)
            outs, costs = self.expand(sv)
            ns = int(outs[a][0])
            got_cols = [(ns >> (p * k)) & ((1 << p) - 1)
                        for k in range(5)]
            assert tuple(got_cols) == st2, (st, a, got_cols, st2)
            assert int(costs[a][0]) == cost2
        # rotation commutation
        for _ in range(500):
            cols = [int(rng.integers(0, 1 << p)) for _ in range(5)]
            a = int(rng.integers(0, 1 << p))
            s = int(rng.integers(1, p))
            st2, _, _ = self.au.step(tuple(cols), a)
            rcols = tuple(AT.rot(c, s, p) for c in cols)
            ra = AT.rot(a, s, p)
            rst2, _, _ = self.au.step(rcols, ra)
            assert tuple(AT.rot(c, s, p) for c in st2) == rst2
        return True

    # ---------------- the BFS ----------------
    def run(self, cap: int, log=print, ckpt_path=None):
        t0 = time.time()
        zero = self.pack([0] * 5)
        buckets: dict[int, list] = {0: [np.array([zero],
                                                dtype=np.uint64)]}
        seen = np.array([zero], dtype=np.uint64)
        returns = {}
        state = {"p": self.p, "cap": cap, "levels": [],
                 "lane": "AB (theta'-dual lane by banked duality)"}
        smask = np.uint64((1 << (5 * self.p)) - 1)
        for c in range(0, cap + 1):
            level_novel = 0
            # zero-cost closure: keep draining bucket c until empty
            while buckets.get(c):
                arrs = buckets.pop(c)
                batch = np.unique(np.concatenate(arrs))
                novel = batch[~np.isin(batch, seen,
                                       assume_unique=False)]
                if c == 0 and level_novel == 0:
                    novel = batch
                if novel.size == 0:
                    break
                level_novel += int(novel.size)
                seen = np.union1d(seen, novel)
                outs, costs = self.expand(novel)
                for a in range(1 << self.p):
                    ns, co = outs[a], costs[a]
                    if c == 0 and a == 0:
                        keep = ns != zero
                        ns, co = ns[keep], co[keep]
                    for w in np.unique(co):
                        cw = c + int(w)
                        if cw > cap:
                            continue
                        sel = ns[co == w]
                        zs = sel[(sel & smask) == 0]
                        if zs.size and cw > 0:
                            accs = (zs >> np.uint64(5 * self.p)) \
                                .astype(int)
                            for acv in np.unique(accs):
                                returns.setdefault(cw, set()) \
                                    .add(int(acv))
                        sel = self.canon(sel)
                        buckets.setdefault(cw, []).append(sel)
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            state["levels"].append(
                {"c": c, "novel": level_novel,
                 "seen": int(seen.size), "rss_mb": rss // 2 ** 20,
                 "t": round(time.time() - t0, 1)})
            nontriv = sorted(w for w, accs in returns.items()
                             if any(acv != 0 for acv in accs))
            log(f"  level {c}: novel {level_novel} seen {seen.size} "
                f"rss {rss//2**20}MB t {time.time()-t0:.0f}s "
                f"returns {sorted(returns)} nontrivial@{nontriv}",)
            state["returns"] = {str(k): sorted(v)
                               for k, v in returns.items()}
            state["nontrivial_costs"] = nontriv
            if ckpt_path:
                ckpt_path.write_text(json.dumps(state, indent=1))
            if rss > RSS_CAP:
                state["aborted"] = "RSS"
                log("  RSS cap — clean abort (levels completed are "
                    "certified)")
                break
        state["wall_s"] = round(time.time() - t0, 1)
        if ckpt_path:
            ckpt_path.write_text(json.dumps(state, indent=1))
        return state


def run_generic(p: int, cap: int):
    r = Racer(p)
    r.validate(nsamples=800)
    print(f"p={p}: validated; run cap {cap}", flush=True)
    st = r.run(cap, log=lambda s: print(s, flush=True),
               ckpt_path=DATA / f"s1_atlas_p{p}_ckpt.json")
    nt = st.get("nontrivial_costs", [])
    done_lv = max((lv["c"] for lv in st["levels"]), default=-1)
    print(f"p={p}: completed level {done_lv}; nontrivial costs {nt}",
          flush=True)
    out_path = DATA / f"s1_atlas_p{p}.json"
    out_path.write_text(json.dumps(st, indent=1))
    print(f"wrote {out_path}", flush=True)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if which in ("all", "controls"):
        for p, cap, expect_min in ((3, 8, 6), (5, 11, None),
                                   (6, 13, 12)):
            r = Racer(p)
            r.validate()
            print(f"p={p}: numpy step + rotation validated", flush=True)
            st = r.run(cap, log=lambda s: print(s, flush=True))
            nt = st.get("nontrivial_costs", [])
            if r.classify:
                mn = nt[0] if nt else None
                print(f"p={p} cap={cap}: min nontrivial = {mn} "
                      f"(expect {expect_min})", flush=True)
                assert mn == expect_min, (p, mn, expect_min)
            else:
                assert st.get("returns"), (p, "no returns at all?")
                print(f"p={p} cap={cap}: barren period — returns "
                      f"{sorted(st['returns'])} are ALL trivial by "
                      f"Theorem A (registers off)", flush=True)
            out[f"p{p}"] = st
    if which == "generic":
        run_generic(int(sys.argv[2]), int(sys.argv[3]))
        return
    if which in ("all", "p9"):
        r = Racer(9)
        r.validate()
        print("p=9: validated; production run cap 16", flush=True)
        st = r.run(16, log=lambda s: print(s, flush=True),
                   ckpt_path=DATA / "s1_atlas9_ckpt.json")
        out["p9"] = st
        nt = st.get("nontrivial_costs", [])
        print(f"p=9 cap 16: nontrivial return costs: {nt}", flush=True)
        if not nt and "aborted" not in st:
            print(">>> no nontrivial compact cycle of weight <= 16: "
                  "with parity + the weight-18 object, "
                  "floor(9) = 18 = 2p (certificate tier, lane AB + "
                  "banked duality)", flush=True)
    (DATA / "s1_atlas9.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s1_atlas9.json'}", flush=True)


if __name__ == "__main__":
    main()
