#!/usr/bin/env python3
"""A42 S3a — the CORRIDOR jet racer: register-reachability pruning
via a backward closure-distance table (the sound remnant of the
falsified sandwich).

Goal: close the p = 12 corner ({18 transverse-only, 20, 22} x
outside the h-DP envelope) by completing jet levels 19..22, which
are unreachable by brute growth (~x3.2/level).  The instrument is a
PRUNE, not a join — the join trilemma (S2 §2.4.1) does not apply.

Soundness (the corridor lemma).  For any compact cycle v of weight w
and any cut gamma, the mechanically-asserted identities of
a42_s2_sandwich (354/354 atlas cycles, every cut) give

    P_F(gamma) + Q_B(gamma) = w + s(gamma),

where s(gamma) = |v1[gamma-1]| + |v1[gamma]| + |v2[gamma]| +
|v2[gamma+1]| + |v2[gamma+2]| = POPCOUNT OF THE FORWARD STATE at the
cut, and the backward (x -> x^-1) automaton visits, at backward cost
Q_B(gamma), the SAME five columns in reversed field order (the map
phi below).  Let G(t) = the backward racer's min arrival cost of the
backward state t (min over all closure data; complete through level
g_max).  Then every cycle of weight w <= W through forward state s
at forward cost c satisfies

    w = c + Q_B - s(gamma) >= c + G(phi(s)) - popcount(s),

so a state may be DROPPED whenever c + G(phi(s)) - popcount(s) > W;
if phi(s) is absent from a table complete to g_max, then
G > g_max and the state may be dropped whenever
W - c + popcount(s) <= g_max.  No state lying on ANY compact cycle
of weight <= W is ever dropped, so the pruned jet run's returns
<= W (with their branch registers) are IDENTICAL to the unpruned
run's — the certificate semantics of a42_s2_jet_racer are preserved
verbatim.  This is the conservative form of register-reachability
pruning: a state that cannot close at all within budget a fortiori
cannot reach a nonzero-register closure.

Engine deltas vs a42_s2_jet_racer (validated by the battery):
  * membership by np.searchsorted against the lex-sorted seen store
    (fiber-width loop; the stock member() re-merge-sorted the full
    12.9M-pair store once per 1K-state chunk — the dominant cost of
    the S2 level-18 flight),
  * the corridor prune above (optional: off => stock semantics),
  * per-level instrumentation: prune counts + distinct omega- and
    barren-window projections of the kept frontier (the S3 quotient
    hypothesis: growth is subexponential in the omega-direction).

Battery (all must pass before production):
  1. validate_banked;
  2. BackRacer.step vs the generic Automaton on the x^-1 pair at
     p in {3, 6, 9, 12} (random states x inputs), canon commutation;
  3. every p in {3, 6} atlas cycle: fwd + bwd replays (the S2a
     asserts), the phi field-reversal correspondence at every cut,
     the popcount form of the cut identity, and G(phi(state)) <= Q_B
     against a full-depth backward table;
  4. 60 random explicit p = 12 boundaries replay through the
     backward automaton (forced-column asserts, cost = weight);
  5. pruned-vs-unpruned RETURN EQUALITY: p = 6 (cap 13, W = 13) and
     p = 9 (cap 16, W in {16, 18}), both branches — the pruned run
     must reproduce the stock jet returns EXACTLY (weights and
     register values).
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


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AT = _load("a40_s4_phase_atlas")
JR = _load("a42_s2_jet_racer")
SW = _load("a42_s2_sandwich")

DATA = LAB / "data" / "a42"
RSS_CAP = int(2.4 * 1024 ** 3)
CHUNK = 1024
BIG = np.int64(1 << 30)


def xinv(supp):
    return [(-i, j) for (i, j) in supp]


def phi(states: np.ndarray, p: int) -> np.ndarray:
    """Field-order reversal of the 5 p-bit columns: the forward
    state at cut gamma (v2[g], v2[g+1], v2[g+2], v1[g-1], v1[g])
    maps to the backward state at the same cut
    (v1[g], v1[g-1], v2[g+2], v2[g+1], v2[g])."""
    m = np.uint64((1 << p) - 1)
    out = np.zeros_like(states)
    for k in range(5):
        f = (states >> np.uint64(p * k)) & m
        out |= f << np.uint64(p * (4 - k))
    return out


def canon5(states: np.ndarray, p: int) -> np.ndarray:
    """Min over the p simultaneous y-rotations of 5 p-bit fields
    (register-free)."""
    m = np.uint64((1 << p) - 1)
    best = states.copy()
    shifts = [np.uint64(p * k) for k in range(5)]
    for s in range(1, p):
        rs = np.zeros_like(states)
        for k in range(5):
            c = (states >> shifts[k]) & m
            c = ((c << np.uint64(s)) | (c >> np.uint64(p - s))) & m
            rs |= c << shifts[k]
        best = np.minimum(best, rs)
    return best


def member_pairs(qs, qr, ss, sr):
    """Membership of pair queries in the lex-sorted pair store via
    searchsorted on the primary key + fiber-width scan."""
    if ss.size == 0:
        return np.zeros(qs.size, dtype=bool)
    left = np.searchsorted(ss, qs, side="left")
    right = np.searchsorted(ss, qs, side="right")
    res = np.zeros(qs.size, dtype=bool)
    width = right - left
    wmax = int(width.max()) if width.size else 0
    for d in range(wmax):
        sel = width > d
        res[sel] |= sr[left[sel] + d] == qr[sel]
    return res


class BackRacer:
    """Register-free racer for the REVERSED pair (x -> x^-1): the
    backward march consumes v2-columns right-to-left and emits
    v1-columns.  Arrival costs = closure distances G for the
    corridor prune.  State layout [f0 f1 o0 o1 o2] =
    (w1[d-1], w1[d], w2[d-3], w2[d-2], w2[d-1]), w_i[dd] = v_i[-dd].
    """

    def __init__(self, p: int):
        assert 5 * p <= 64
        self.p = p
        self.au = AT.Automaton(xinv(AT.A_L), xinv(AT.B_L), p)
        au = self.au
        assert (au.nf, au.no, au.forced_blk) == (2, 3, 0)
        assert au.adv_f == 1 and au.adv_o == 0 and au.top_j == 3
        assert sorted(au.terms_f) == [(-1, 0), (0, 0)]
        assert sorted(au.terms_o) == [(-3, 1), (0, -1), (0, 0)]
        self.mask = np.uint64((1 << p) - 1)
        self.pc = np.array([bin(i).count("1") for i in range(1 << p)],
                           dtype=np.int64)

    def pack(self, cols):
        s = 0
        for k, c in enumerate(cols):
            s |= int(c) << (self.p * k)
        return np.uint64(s)

    def expand(self, states: np.ndarray):
        p = self.p
        m = self.mask

        def rot(x, s):
            s %= p
            if s == 0:
                return x
            return ((x << np.uint64(s)) | (x >> np.uint64(p - s))) & m

        f0 = states & m
        f1 = (states >> np.uint64(p)) & m
        o0 = (states >> np.uint64(2 * p)) & m
        o1 = (states >> np.uint64(3 * p)) & m
        o2 = (states >> np.uint64(4 * p)) & m
        base = f0 ^ f1 ^ rot(o0, 1)
        outs = []
        for a in range(1 << p):
            av = np.uint64(a)
            acc = base ^ av ^ rot(av, -1 % p)
            new_f = rot(acc, -3 % p)
            cost = self.pc[a] + self.pc[new_f.astype(np.int64)]
            ns = f1 | (new_f << np.uint64(p)) | (o1 << np.uint64(2 * p)) \
                | (o2 << np.uint64(3 * p)) | (av << np.uint64(4 * p))
            outs.append((ns, cost))
        return outs

    def validate(self, nsamples=500, rng=None):
        rng = rng or np.random.default_rng(31)
        p = self.p
        for _ in range(nsamples):
            cols = [int(rng.integers(0, 1 << p)) for _ in range(5)]
            a = int(rng.integers(0, 1 << p))
            st2, cost2, _ = self.au.step(tuple(cols), a)
            sv = np.array([self.pack(cols)], dtype=np.uint64)
            ns, co = self.expand(sv)[a]
            got = [(int(ns[0]) >> (p * k)) & ((1 << p) - 1)
                   for k in range(5)]
            assert tuple(got) == st2, (cols, a, got, st2)
            assert int(co[0]) == cost2
        # canon commutation
        for _ in range(200):
            cols = [int(rng.integers(0, 1 << p)) for _ in range(5)]
            s = int(rng.integers(1, p))
            sv = np.array([self.pack(cols)], dtype=np.uint64)
            rcols = [AT.rot(c, s, p) for c in cols]
            sv2 = np.array([self.pack(rcols)], dtype=np.uint64)
            assert int(canon5(sv, p)[0]) == int(canon5(sv2, p)[0])
        return True

    def run_table(self, gmax: int, log=print, ckpt_stem=None,
                  time_budget_s=None, growth_est=2.9):
        """Dial-bucket BFS to level gmax; returns (sorted canonical
        states u64, arrival costs u8, meta).  Completing level g
        certifies: every backward-reachable canonical state with
        min closure cost <= g is in the table at its exact cost."""
        t0 = time.time()
        zero = np.uint64(0)
        buckets: dict[int, list] = {0: [np.array([zero],
                                                 dtype=np.uint64)]}
        seen_s = np.array([zero], dtype=np.uint64)
        seen_c = np.zeros(1, dtype=np.uint8)
        meta = {"p": self.p, "gmax_req": gmax, "levels": [],
                "engine": "backward closure-distance table "
                          "(x^-1 pair, standard charging)"}
        aborted = False
        last_level_t = 0.0
        for c in range(0, gmax + 1):
            if time_budget_s is not None and c > 0:
                if (time.time() - t0) + last_level_t * growth_est \
                        > time_budget_s:
                    meta["aborted"] = f"time budget before level {c}"
                    log(f"  bwd time budget — stop before level {c}")
                    aborted = True
                    break
            lt0 = time.time()
            level_novel = 0
            while buckets.get(c):
                arrs = buckets.pop(c)
                batch = np.unique(np.concatenate(arrs))
                idx = np.searchsorted(seen_s, batch)
                idx[idx >= seen_s.size] = seen_s.size - 1
                fresh = seen_s[idx] != batch
                if c == 0 and level_novel == 0:
                    fresh[:] = True
                novel = batch[fresh]
                if novel.size == 0:
                    break
                level_novel += int(novel.size)
                seen_s = np.concatenate([seen_s, novel])
                seen_c = np.concatenate(
                    [seen_c, np.full(novel.size, c, dtype=np.uint8)])
                o = np.argsort(seen_s, kind="stable")
                seen_s, seen_c = seen_s[o], seen_c[o]
                for lo in range(0, novel.size, CHUNK):
                    chunk = novel[lo:lo + CHUNK]
                    outs = self.expand(chunk)
                    per_cost: dict[int, list] = {}
                    for a in range(1 << self.p):
                        ns, co = outs[a]
                        outs[a] = None
                        if c == 0 and a == 0:
                            keep = ns != zero
                            ns, co = ns[keep], co[keep]
                        for w in np.unique(co):
                            cw = c + int(w)
                            if cw > gmax:
                                continue
                            per_cost.setdefault(cw, []).append(
                                ns[co == w])
                    for cw, lst in per_cost.items():
                        merged = np.unique(canon5(
                            np.concatenate(lst), self.p))
                        idx = np.searchsorted(seen_s, merged)
                        idx[idx >= seen_s.size] = seen_s.size - 1
                        merged = merged[seen_s[idx] != merged]
                        if merged.size:
                            buckets.setdefault(cw, []).append(merged)
                    del outs, per_cost
                    rss = resource.getrusage(
                        resource.RUSAGE_SELF).ru_maxrss
                    if rss > RSS_CAP:
                        meta["aborted"] = "RSS"
                        log(f"  bwd RSS cap ({rss//2**20}MB) — clean "
                            "abort")
                        aborted = True
                        break
                if aborted:
                    break
                for cw in list(buckets.keys()):
                    arrs2 = buckets[cw]
                    tot = sum(x.size for x in arrs2)
                    if len(arrs2) > 16 or tot > 1 << 20:
                        merged = np.unique(np.concatenate(arrs2))
                        idx = np.searchsorted(seen_s, merged)
                        idx[idx >= seen_s.size] = seen_s.size - 1
                        merged = merged[seen_s[idx] != merged]
                        buckets[cw] = [merged] if merged.size else []
            if aborted:
                break
            last_level_t = time.time() - lt0
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            meta["levels"].append(
                {"g": c, "novel": level_novel,
                 "seen": int(seen_s.size), "rss_mb": rss // 2 ** 20,
                 "t": round(time.time() - t0, 1)})
            log(f"  bwd level {c}: novel {level_novel} "
                f"seen {seen_s.size} rss {rss//2**20}MB "
                f"t {time.time()-t0:.0f}s")
            meta["gmax_completed"] = c
            if ckpt_stem and c >= 10:
                np.savez_compressed(
                    str(ckpt_stem) + ".npz",
                    states=seen_s, costs=seen_c,
                    meta=json.dumps(meta))
        meta["wall_s"] = round(time.time() - t0, 1)
        if ckpt_stem:
            np.savez_compressed(str(ckpt_stem) + ".npz",
                                states=seen_s, costs=seen_c,
                                meta=json.dumps(meta))
        return seen_s, seen_c, meta


class CorridorJet(JR.JetRacer):
    """The S2 jet racer + searchsorted membership + the corridor
    prune.  With gtab=None the semantics are stock (used for the
    pruned-vs-unpruned battery)."""

    def set_corridor(self, gtab_states, gtab_costs, gmax_completed,
                     W):
        self.gs = gtab_states
        self.gc = gtab_costs.astype(np.int64)
        self.gmax = int(gmax_completed)
        self.W = int(W)

    def corridor_keep(self, ms: np.ndarray, cw: int):
        """Keep mask + stats for canonical forward states ms at
        forward cost cw."""
        pcs = np.bitwise_count(ms).astype(np.int64)
        assert pcs.max(initial=0) <= cw, "state bits exceed paid cost"
        need = self.W - cw + pcs          # available Q_B budget
        bq = canon5(phi(ms, self.p), self.p)
        idx = np.searchsorted(self.gs, bq)
        idx[idx >= self.gs.size] = self.gs.size - 1
        found = self.gs[idx] == bq
        G = np.where(found, self.gc[idx], BIG)
        keep = np.where(found, G <= need, need >= self.gmax + 1)
        stats = dict(gen=int(ms.size),
                     pruned_known=int((found & ~keep).sum()),
                     pruned_absent=int((~found & ~keep).sum()),
                     kept_absent=int((~found & keep).sum()))
        return keep, stats

    def run_corridor(self, cap: int, log=print, ckpt_path=None,
                     time_budget_s=None, growth_est=3.2,
                     expected_returns=None, dump_kept_from=None,
                     dump_stem=None):
        """The stock jet BFS with fast membership + the prune.
        expected_returns: {cost: sorted register list} to assert
        against (banked regression).  dump_kept_from: level from
        which the kept novel frontier is saved (npz) for the
        Stage-2 census."""
        t0 = time.time()
        prune_on = getattr(self, "gs", None) is not None
        zstate = np.uint64(0)
        buckets: dict[int, list] = {
            0: [(np.array([0], dtype=np.uint64),
                 np.zeros(1, dtype=np.uint8))]}
        seen_s = np.array([0], dtype=np.uint64)
        seen_r = np.zeros(1, dtype=np.uint8)
        returns: dict[int, set] = {}
        state = {"p": self.p, "cap": cap, "branch": self.branch,
                 "xstar": int(self.xstar), "levels": [],
                 "engine": "corridor jet (searchsorted member"
                           + (", corridor prune" if prune_on else "")
                           + ")",
                 "W": getattr(self, "W", None),
                 "gmax": getattr(self, "gmax", None),
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
            lstats = dict(gen=0, pruned_known=0, pruned_absent=0,
                          kept_absent=0)
            level_kept: list[np.ndarray] = []
            while buckets.get(c):
                arrs = buckets.pop(c)
                sts = np.concatenate([x[0] for x in arrs])
                rgs = np.concatenate([x[1] for x in arrs])
                sts, rgs = self.dedup(sts, rgs)
                fresh = ~member_pairs(sts, rgs, seen_s, seen_r)
                if c == 0 and level_novel == 0:
                    fresh[:] = True
                sts, rgs = sts[fresh], rgs[fresh]
                if sts.size == 0:
                    break
                level_novel += int(sts.size)
                level_kept.append(sts.copy())
                seen_s = np.concatenate([seen_s, sts])
                seen_r = np.concatenate([seen_r, rgs])
                o = np.lexsort((seen_r, seen_s))
                seen_s, seen_r = seen_s[o], seen_r[o]
                for lo in range(0, sts.size, CHUNK):
                    cst = sts[lo:lo + CHUNK]
                    crg = rgs[lo:lo + CHUNK]
                    if c == cap and c > 0:
                        # FINAL LEVEL: only zero-cost continuations
                        # can matter (positive cost exceeds cap), and
                        # the zero-cost step is deterministic
                        # (input 0, forced 0).  Expand a = 0 only and
                        # keep only cost-0 outputs.
                        ns, nreg, co = self.expand(cst, crg)[0]
                        sel = co == 0
                        ns, nreg = ns[sel], nreg[sel]
                        zsel = ns == zstate
                        if zsel.any():
                            for rv in np.unique(nreg[zsel]):
                                returns.setdefault(c, set()) \
                                    .add(int(rv))
                        per_cost = {c: [(ns, nreg)]} if ns.size \
                            else {}
                    else:
                        outs = self.expand(cst, crg)
                        per_cost = {}
                        for a in range(1 << self.p):
                            ns, nreg, co = outs[a]
                            outs[a] = None
                            if c == 0 and a == 0:
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
                        ms, mr = self.dedup(ms, mr)
                        if prune_on:
                            keep, kst = self.corridor_keep(ms, cw)
                            for k in lstats:
                                lstats[k] += kst[k]
                            if not keep.all():
                                ms, mr = ms[keep], mr[keep]
                        if ms.size == 0:
                            continue
                        fr = ~member_pairs(ms, mr, seen_s, seen_r)
                        if fr.any():
                            buckets.setdefault(cw, []).append(
                                (ms[fr], mr[fr]))
                    del per_cost
                    rss = resource.getrusage(
                        resource.RUSAGE_SELF).ru_maxrss
                    if rss > RSS_CAP:
                        state["aborted"] = "RSS"
                        log(f"  RSS cap ({rss//2**20}MB) — clean "
                            "abort (completed levels certified)")
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
                        fr = ~member_pairs(ms, mr, seen_s, seen_r)
                        buckets[cw] = [(ms[fr], mr[fr])] if fr.any() \
                            else []
            if aborted:
                break
            last_level_t = time.time() - lt0
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            rec = {"c": c, "novel": level_novel,
                   "seen": int(seen_s.size),
                   "rss_mb": rss // 2 ** 20,
                   "t": round(time.time() - t0, 1)}
            if prune_on:
                rec["prune"] = lstats
            # omega/barren window projections of the level's novel
            if level_kept:
                lv = np.concatenate(level_kept)
                rec["n_omega"] = int(np.unique(
                    self._proj(lv, self.ev)).size)
                if self.p == 3 * (1 << self.a):
                    rec["n_barren"] = int(np.unique(
                        self._proj(lv, self._evb())).size)
            state["levels"].append(rec)
            nontriv = sorted(w for w, rv in returns.items()
                             if any(v != 0 for v in rv))
            log(f"  level {c}: novel {level_novel} seen {seen_s.size} "
                f"rss {rss//2**20}MB t {time.time()-t0:.0f}s "
                f"returns {sorted(returns)} nonzero-reg@{nontriv}"
                + (f" prune {lstats}" if prune_on else ""))
            state["returns"] = {str(k): sorted(v)
                                for k, v in returns.items()}
            state["nonzero_reg_costs"] = nontriv
            if expected_returns is not None:
                for k, v in returns.items():
                    if k in expected_returns:
                        assert sorted(v) == expected_returns[k], \
                            ("RETURN MISMATCH vs banked", k,
                             sorted(v), expected_returns[k])
            if dump_stem and level_kept and dump_kept_from is not None \
                    and c >= dump_kept_from:
                np.savez_compressed(
                    f"{dump_stem}_L{c}.npz",
                    states=np.concatenate(level_kept))
            if ckpt_path:
                ckpt_path.write_text(json.dumps(state, indent=1))
        state["wall_s"] = round(time.time() - t0, 1)
        if ckpt_path:
            ckpt_path.write_text(json.dumps(state, indent=1))
        return state

    def _proj(self, states: np.ndarray, tab: np.ndarray):
        """Pack per-column projections tab[col] (uint8) of the 5
        fields into one uint64 fingerprint."""
        p = self.p
        m = self.mask
        out = np.zeros_like(states)
        for k in range(5):
            f = ((states >> np.uint64(p * k)) & m).astype(np.int64)
            out |= tab[f].astype(np.uint64) << np.uint64(8 * k)
        return out

    def _evb(self):
        """Barren-part projection table col -> col mod (y+1)^{2^a}
        (valid at p = 3*2^a)."""
        if not hasattr(self, "_evb_tab"):
            import a42_lib as AL
            m = 0b11
            for _ in range(self.a):
                m = AL.pmul(m, m)
            self._evb_tab = np.array(
                [AL.pmod(c, m) for c in range(1 << self.p)],
                dtype=np.uint8)
        return self._evb_tab


# ----------------------- the battery -------------------------------
def battery():
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)

    # 2. backward step validation
    for p in (3, 6, 9, 12):
        BackRacer(p).validate()
    print("BackRacer step/canon validated at p in {3,6,9,12} "
          "(vs generic Automaton on the x^-1 pair)", flush=True)

    # 3. atlas replay battery with phi + popcount identity + G <= Q_B
    for p, W, gmax in ((3, 8, 16), (6, 12, 18)):
        au = AT.Automaton(AT.A_L, AT.B_L, p)
        bau = AT.Automaton(xinv(AT.A_L), xinv(AT.B_L), p)
        br = BackRacer(p)
        gs, gc, meta = br.run_table(gmax, log=lambda s: None)
        assert meta["gmax_completed"] == gmax
        rows, _ = AT.atlas("AB", p, W, keep_pts=True)
        nasrt = 0
        for row in rows:
            v1, v2 = SW.cols_from_pts(row["pts"])
            w = row["weight"]
            SW.check_cycle_conv(v1, v2, p)
            cuts, cf = SW.fwd_replay_with_cuts(au, v1, v2, p)
            assert cf == w
            qb, cb = SW.bwd_replay_with_cuts(bau, v1, v2, p)
            assert cb == w
            for (g, pf) in cuts:
                if g not in qb:     # one degenerate boundary cut
                    assert pf in (0, w)
                    continue
                fwd_cols = [v2.get(g, 0), v2.get(g + 1, 0),
                            v2.get(g + 2, 0), v1.get(g - 1, 0),
                            v1.get(g, 0)]
                bwd_cols = list(reversed(fwd_cols))
                fpk = np.array([br.pack(fwd_cols)], dtype=np.uint64)
                bpk = np.array([br.pack(bwd_cols)], dtype=np.uint64)
                assert int(phi(fpk, p)[0]) == int(bpk[0])
                s_here = int(np.bitwise_count(fpk[0]))
                assert pf + qb[g] == w + s_here, (g, pf, qb[g], w,
                                                 s_here)
                bq = canon5(bpk, p)
                i = int(np.searchsorted(gs, bq[0]))
                if i < gs.size and gs[i] == bq[0]:
                    assert int(gc[i]) <= qb[g], \
                        ("G exceeds Q_B", p, g, int(gc[i]), qb[g])
                else:
                    # absence certifies G > gmax; consistent only
                    # if this cycle's own Q_B exceeds the depth
                    assert qb[g] > gmax, \
                        ("cut state missing though Q_B <= gmax",
                         p, g, qb[g], gmax)
                nasrt += 1
        print(f"p={p}: {len(rows)} atlas cycles — phi correspondence, "
              f"popcount cut identity, G<=Q_B asserted at {nasrt} "
              "cuts (bwd table complete to "
              f"{meta['gmax_completed']})", flush=True)

    # 4. p=12 explicit boundaries through the backward automaton
    rng = np.random.default_rng(23)
    bau12 = AT.Automaton(xinv(AT.A_L), xinv(AT.B_L), 12)
    ndone = 0
    for _ in range(60):
        v1, v2 = JR.boundary_cycle_inputs(12, rng)
        if not v1 and not v2:
            continue
        qb, cb = SW.bwd_replay_with_cuts(bau12, v1, v2, 12)
        wt = sum(bin(c).count("1") for c in v1.values()) + \
            sum(bin(c).count("1") for c in v2.values())
        assert cb == wt
        ndone += 1
    print(f"p=12: {ndone} explicit boundaries replay through the "
          "backward automaton (forced columns + cost = weight)",
          flush=True)

    # 5. pruned-vs-unpruned return equality
    for p, cap, Ws, gmax in ((6, 13, (13,), 18),
                             (9, 16, (16, 18), 18)):
        br = BackRacer(p)
        gs, gc, meta = br.run_table(gmax, log=lambda s: None)
        for branch in (0, 1):
            stock = JR.JetRacer(p, branch)
            stock.validate()
            st0 = stock.run(cap, log=lambda s: None)
            for W in Ws:
                cj = CorridorJet(p, branch)
                cj.set_corridor(gs, gc, meta["gmax_completed"], W)
                st1 = cj.run_corridor(cap, log=lambda s: None)
                assert st1["returns"] == st0["returns"], \
                    (p, branch, W, st1["returns"], st0["returns"])
                assert st1["nonzero_reg_costs"] == \
                    st0["nonzero_reg_costs"]
                pr = [lv.get("prune") for lv in st1["levels"]
                      if lv.get("prune")]
                tot_pruned = sum(x["pruned_known"] +
                                 x["pruned_absent"] for x in pr)
                print(f"p={p} branch {branch} W={W}: pruned run "
                      f"returns == stock ({st0['returns']}); "
                      f"{tot_pruned} states pruned", flush=True)
    print("BATTERY GREEN", flush=True)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "battery"

    if which == "battery":
        battery()
        return

    if which == "bwd12":
        gmax = int(sys.argv[2]) if len(sys.argv) > 2 else 17
        budget = float(sys.argv[3]) * 3600 if len(sys.argv) > 3 \
            else 1.5 * 3600
        br = BackRacer(12)
        br.validate()
        print(f"p=12 backward table: gmax {gmax}, budget "
              f"{budget/3600:.1f} h", flush=True)
        gs, gc, meta = br.run_table(
            gmax, log=lambda s: print(s, flush=True),
            ckpt_stem=DATA / "s3_bwd12_table",
            time_budget_s=budget)
        print(f"bwd12 done: completed level "
              f"{meta.get('gmax_completed')} states {gs.size} "
              f"({meta['wall_s']} s)", flush=True)
        (DATA / "s3_bwd12_meta.json").write_text(
            json.dumps(meta, indent=1))
        return

    if which == "corridor":
        tag = sys.argv[2]
        cap = int(sys.argv[3])
        W = int(sys.argv[4])
        budget = float(sys.argv[5]) * 3600
        brch = {"T": 0, "V": 1}[tag]
        z = np.load(DATA / "s3_bwd12_table.npz")
        meta = json.loads(z["meta"].item())
        gs, gc = z["states"], z["costs"]
        assert meta["p"] == 12
        cj = CorridorJet(12, brch)
        cj.validate()
        cj.set_corridor(gs, gc, meta["gmax_completed"], W)
        banked = json.loads((DATA / "s2_jet12.json").read_text())[tag]
        expected = {int(k): v for k, v in banked["returns"].items()}
        print(f"p=12 corridor branch {tag} (xstar={cj.xstar}): "
              f"cap {cap} W {W} gmax {meta['gmax_completed']} "
              f"(bwd states {gs.size}), budget {budget/3600:.1f} h; "
              f"asserting banked returns {expected}", flush=True)
        st = cj.run_corridor(
            cap, log=lambda s: print(s, flush=True),
            ckpt_path=DATA / f"s3_corridor12_{tag}_ckpt.json",
            time_budget_s=budget, expected_returns=expected,
            dump_kept_from=19,
            dump_stem=str(DATA / f"s3_corridor12_{tag}_frontier"))
        done = max((lv["c"] for lv in st["levels"]), default=-1)
        nt = st.get("nonzero_reg_costs", [])
        print(f"p=12 corridor branch {tag}: completed level {done}; "
              f"nonzero-register costs {nt}", flush=True)
        (DATA / f"s3_corridor12_{tag}.json").write_text(
            json.dumps(st, indent=1))
        return

    if which == "summary":
        out = {}
        for tag in ("T", "V"):
            out[tag] = json.loads(
                (DATA / f"s3_corridor12_{tag}.json").read_text())
        done = min(max((lv["c"] for lv in out[t]["levels"]),
                       default=-1) for t in ("T", "V"))
        ntU = sorted(set(out["T"].get("nonzero_reg_costs", []))
                     | set(out["V"].get("nonzero_reg_costs", [])))
        print(f"corridor branches complete through joint level "
              f"{done}; union nonzero-register costs {ntU}",
              flush=True)
        if not ntU and done >= 0:
            ev = done if done % 2 == 0 else done - 1
            print(f">>> no nonzero-register return <= {done}, either "
                  "branch, corridor-pruned enumeration (soundness: "
                  "cut-identity corridor lemma, battery green) + "
                  "joint-kernel emptiness (s2_registers) => no "
                  f"nontrivial compact cycle of weight <= {done}; "
                  f"parity => floor(12) >= {ev + 2} unrestricted, "
                  "ALL classes (certificate tier, lane AB + banked "
                  "duality)", flush=True)
        (DATA / "s3_corridor12.json").write_text(
            json.dumps(out, indent=1))
        print(f"wrote {DATA/'s3_corridor12.json'}", flush=True)


if __name__ == "__main__":
    main()
