#!/usr/bin/env python3
"""A40 S7 — Stage 1 readouts: validation, specimens, and the
interface-tax table.

validate:  (i) pinch inequality on TC63's own heavy blocks (the
           only sub-rate-2 species that crosses the heavy line);
           (ii) TC63 crossing-links dominated by the analytic link
           grant used for u >= 5 strata in the s7 assembly;
           (iii) pinch on every enumerated link specimen (via the
           specimens subcommand's verification loop).
specimens: parented link march at modest caps; top-D links per h
           replayed row-by-row and INDEPENDENTLY verified through
           CoverFragment (E_j everywhere incl. the heavy block,
           slab classes vs the phase schedule, window rule, weight,
           drift); reports per-slab profiles and the two-ghost
           composition verdict (can both light stretches coast at
           ghost grade <= 2/slab?).
tax:       tau_u(h, L) = S6grant(h, L) - J_u(h, L), where
           S6grant = max_{a+b=h-L-1} [D^fwd_u(1+a) + D^bwd_1(b)]
           is what the s6 assembly effectively granted the same
           slabs (runs at standalone maxima, heavy slabs free), and
           J_u = the measured joint link deficit.  Also the T2-form
           readout: T_link = max certified [D - (6/7)(h - L)].
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

DATA = LAB / "data" / "a40"

from a40_s6_frontier import (  # noqa: E402
    wt, lsb, norm, seeds_full,
)
from a40_s7_link import LinkMarch, PRE, HEAVY, POST  # noqa: E402


# ---------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------

def validate():
    sys.argv = [sys.argv[0]]
    import numpy as np  # noqa: F401
    from a40_s6_drift import lift_phase, load_survivor
    from a40_s5_lightcore import Phase

    out = {}
    print("validate 1: TC63 (l=24) heavy blocks — pinch + link "
          "grants", flush=True)
    pts = load_survivor("s5_dense_l24p6.json", 6, 3, 10)
    ph = Phase.from_quotient_pts(24, 6, 3, pts)
    fr, s = lift_phase(ph, n_periods=4)
    assert s == 3
    slabs = fr.slabs()
    anchors = fr.anchors()
    n = len(slabs)
    print(f"  slabs (4 periods): {slabs}", flush=True)
    # classify + find heavy blocks with light flanks inside the lift
    heavy = [i for i, w in enumerate(slabs) if w >= 8]
    pinch_checks = []
    for i in range(1, n - 1):
        if slabs[i] >= 8:
            # find the maximal heavy block containing i
            a = i
            while a > 0 and slabs[a - 1] >= 8:
                a -= 1
            b = i
            while b < n - 1 and slabs[b + 1] >= 8:
                b += 1
            if a == i and a > 0 and b < n - 1:
                L = b - a + 1
                if L <= 3:
                    ok = slabs[a - 1] + slabs[b + 1] >= 8
                    pinch_checks.append(
                        dict(block=[a, b], L=L, exit=slabs[a - 1],
                             entry=slabs[b + 1], ok=bool(ok)))
                    assert ok, (a, b, slabs)
    print(f"  pinch on {len(pinch_checks)} blocks: all pass "
          f"{pinch_checks[:3]}", flush=True)
    out["tc63_pinch"] = pinch_checks

    # crossing links: seed at a min light slab, cross one block,
    # end at the next min light slab.  Check the ANALYTIC link
    # grant D <= (1+a)(2-u/4) + 0*L + b*(7/4)  [pre slabs >= u,
    # heavy slabs <= 0 deficit, post slabs granted at the u>=1
    # envelope 1.75] dominates the true link deficits.
    grants_ok = []
    lights = [i for i, w in enumerate(slabs) if w <= 7]
    for i0 in lights:
        for i1 in lights:
            if i1 <= i0:
                continue
            seg = slabs[i0:i1 + 1]
            if not any(w >= 8 for w in seg):
                continue
            if i1 - i0 + 1 > 12:
                continue
            L = sum(1 for w in seg if w >= 8)
            u = min(w for w in seg if w <= 7)
            h = len(seg)
            D = 2 * h - sum(seg) / 4
            delta = anchors[i1] - anchors[i0]
            # analytic grant with the pre/post split at the block
            first_heavy = next(j for j, w in enumerate(seg)
                               if w >= 8)
            a = first_heavy - 1
            b = h - L - 1 - a
            grant = (1 + a) * (2 - u / 4) + b * 1.75
            grants_ok.append(bool(D <= grant + 1e-9))
            assert D <= grant + 1e-9, (i0, i1, seg, D, grant)
    print(f"  analytic link grant dominates {len(grants_ok)} TC63 "
          f"crossing segments: PASS", flush=True)
    out["tc63_links_checked"] = len(grants_ok)
    (DATA / "s7_validate.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s7_validate.json'}", flush=True)


# ---------------------------------------------------------------------
# specimens (parented link march + independent verification)
# ---------------------------------------------------------------------

class ParentedLinkMarch(LinkMarch):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.parents = {}
        self.seed_of = {}

    def run(self, seed_states, **kw):
        self._cur_seed = None
        return super().run(seed_states, **kw)

    # capture seeds: LinkMarch.run pushes seeds directly; we tag
    # roots by intercepting push during seeding via nodes_pushed==0
    # bookkeeping — simpler: override record for roots is messy, so
    # we re-implement seeding by monkey-running: parents of root
    # keys are (None, seed_w8).
    def push(self, g2, state_rows, slab_mask, slab_hi, phase2, L2,
             h, dlt_new, M, anch, xlo, xhi, hin2, heap, seen):
        n_before = self.nodes_pushed
        super().push(g2, state_rows, slab_mask, slab_hi, phase2, L2,
                     h, dlt_new, M, anch, xlo, xhi, hin2, heap, seen)
        if self.nodes_pushed > n_before:      # actually pushed
            from a40_s6_frontier import lsb as _lsb
            nrows, S = norm(state_rows)
            anch_c = _lsb(slab_mask) - S
            key = (nrows, anch_c, phase2, L2, h + 1, dlt_new)
            self.parents[key] = (self._pkey, S, M)

    def expand(self, g, dyn, phase, L, h, dlt, anch, xlo, xhi, hin,
               heap, seen):
        self._pkey = (dyn, anch, phase, L, h, dlt)
        super().expand(g, dyn, phase, L, h, dlt, anch, xlo, xhi,
                       hin, heap, seen)


def replay(m, key, seeds):
    """Walk parents back to a seed; rebuild rows in the seed frame.
    Root keys are seeded states (h=1) whose dyn matches a seed."""
    chain = []
    k = key
    while k in m.parents:
        chain.append(k)
        k = m.parents[k][0]
    chain.append(k)               # the root (h=1)
    chain.reverse()
    root_dyn = chain[0][0]
    seed_w8 = None
    for w8 in seeds:
        dyn = (w8[0], w8[1], w8[2], w8[3], w8[5], w8[6], w8[7])
        if dyn == root_dyn:
            seed_w8 = w8
            break
    assert seed_w8 is not None, "root seed not found"
    offs = [0]
    for i in range(1, len(chain)):
        _, S, M = m.parents[chain[i]]
        offs.append(offs[-1] + M - S)
    rows = {}
    for r in range(4):
        rows[r] = (frozenset(c for c in range(80)
                             if seed_w8[r] >> c & 1),
                   frozenset(c for c in range(80)
                             if seed_w8[4 + r] >> c & 1))
    for i in range(1, len(chain)):
        dyn = chain[i][0]
        off = offs[i]
        v1new, v2new = dyn[3], dyn[6]
        rows[3 + i] = (
            frozenset((c - off) for c in range(v1new.bit_length())
                      if v1new >> c & 1),
            frozenset((c - off) for c in range(v2new.bit_length())
                      if v2new >> c & 1))
    return rows, chain


def specimens(args):
    """Sequential per-seed parented runs; each seed's top-D POST
    specimens are replayed and verified immediately, then the march
    is dropped (memory-bounded).  The verified records are merged
    per h across seeds."""
    t0 = time.time()
    sys.argv = [sys.argv[0]]
    from a40_s6_drift import CoverFragment
    import resource

    def rss_mb():
        return resource.getrusage(
            resource.RUSAGE_SELF).ru_maxrss // (1024 * 1024)
    sds = [w8 for w8 in seeds_full(args.u)
           if sum(wt(r) for r in w8) == args.u]
    lo, hi = 0, len(sds)
    if args.seeds:
        lo, hi = map(int, args.seeds.split(":"))
    best = {}                     # h -> verified record dict
    infos = []
    for si in range(lo, hi):
        w8 = sds[si]
        m = ParentedLinkMarch(args.u, kmax=args.kmax,
                              whcap=args.whcap, gcap=args.gcap,
                              hcap=args.hcap, dcap=args.dcap)
        info = m.run([w8], log=False)
        infos.append(info["popped"])
        print(f"  seed {si}: {info['popped']} nodes, "
              f"{len(m.tab_link)} link buckets, rss {rss_mb()}MB",
              flush=True)
        stop_after = rss_mb() > 1900
        cand = {}
        for key in m.parents:
            dyn, anch, phase, L, h, dlt = key
            if phase != POST:
                continue
            g = m.tab_link_L.get((h, L, dlt))
            if g is None:
                continue
            D = 2 * h - g / 4
            if h not in cand or D > cand[h][0]:
                cand[h] = (D, key, g)
        for h, (D, key, g) in cand.items():
            if h in best and best[h]["D"] >= D:
                continue
            rows, chain = replay(m, key, [w8])
            lo, hi = min(rows), max(rows)
            fr = CoverFragment([rows[j] for j in range(lo, hi + 1)],
                               lo)
            ok = fr.admissible()
            sl = fr.slabs()
            heavy_idx = [i for i, w in enumerate(sl) if w >= 8]
            pinch_ok = True
            slip = None
            if heavy_idx:
                a, b = heavy_idx[0], heavy_idx[-1]
                if a > 0 and b < len(sl) - 1 and b - a + 1 <= 3:
                    pinch_ok = sl[a - 1] + sl[b + 1] >= 8
                anchors = fr.anchors()
                if 0 < a and b + 1 < len(anchors):
                    slip = anchors[b + 1] - anchors[a - 1]
                pre_sl, post_sl = sl[:a], sl[b + 1:]
            else:
                pre_sl, post_sl = sl, []

            def longest_ghost(ws):
                run = bst = 0
                for w_ in ws:
                    run = run + 1 if w_ <= 2 else 0
                    bst = max(bst, run)
                return bst
            dlt_key = key[5]
            assert ok, f"specimen h={h} FAILS E-check"
            assert sum(sl) == g, (sum(sl), g)
            assert fr.drift() == dlt_key, (fr.drift(), dlt_key)
            assert not fr.window_prune_events()
            assert pinch_ok, (h, sl, "pinch violated!")
            best[h] = dict(
                h=h, D=round(D, 3), g=g, delta=dlt_key,
                slabs=sl, weight=fr.weight(),
                block_weight=(sum(sl[heavy_idx[0]:heavy_idx[-1]
                                     + 1]) if heavy_idx else 0),
                slip_across_block=slip,
                ghost_pre=longest_ghost(pre_sl),
                ghost_post=longest_ghost(post_sl),
                seed=si)
        del m
        if stop_after:
            print(f"  RED: RSS near budget after seed {si}; "
                  f"stopping early (specimens are demonstrations, "
                  f"partial coverage is honest)", flush=True)
            break
    ver = [best[h] for h in sorted(best)][-args.n:]
    for rec in ver:
        print(f"  h={rec['h']}: D={rec['D']} g={rec['g']} "
              f"delta={rec['delta']} slabs={rec['slabs']} "
              f"slip={rec['slip_across_block']} ghost(pre,post)="
              f"({rec['ghost_pre']},{rec['ghost_post']}) "
              f"E-verified", flush=True)
    out = dict(params=vars(args), n_nodes=infos, specimens=ver,
               wall_s=round(time.time() - t0, 1))
    p = DATA / (f"s7_specimens_u{args.u}_g{args.gcap}"
                f"_d{args.dcap}.json")
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({out['wall_s']} s)", flush=True)


# ---------------------------------------------------------------------
# tax table
# ---------------------------------------------------------------------

def load_s6_stratum(u):
    """(fwd table, bwd table, gcap) for stratum u from s6 prod."""
    p = DATA / f"s6_frontier_u{u}prod.json"
    d = json.loads(p.read_text())
    f = {tuple(map(int, k.split(","))): g
         for k, g in d["fwd"]["table"].items()}
    b = {tuple(map(int, k.split(","))): g
         for k, g in d["bwd"]["table"].items()}
    return f, b, d["fwd"]["params"]["gcap"]


def best_D(tab, h, gcap):
    """max deficit at exactly h slabs; absent h -> cap/envelope."""
    cands = [2 * h - g / 4 for (hh, d), g in tab.items() if hh == h]
    cap = 2 * h - (gcap + 1) / 4
    env = 1.75 * h
    return max(cands) if cands else min(cap, env)


def tax(args):
    d = json.loads((DATA / args.link).read_text())
    u = d["params"]["u"]
    gcap = d["params"]["gcap"]
    tabL = {tuple(map(int, k.split(","))): g
            for k, g in d["tab_link_L"].items()}
    fw_u, _, gcap_u = load_s6_stratum(u)
    _, bw_1, gcap_1 = load_s6_stratum(1)
    rows = []
    for L in sorted({L for (_, L, _) in tabL}):
        hs = sorted({h for (h, LL, _) in tabL if LL == L})
        for h in hs:
            Ds = [2 * h - g / 4 for (hh, LL, dd), g in tabL.items()
                  if hh == h and LL == L]
            J = max(Ds)
            # s6-style grant for the same slabs
            grant = -1e9
            for a in range(0, h - L - 1):
                b = h - L - 1 - a
                if b < 1:
                    continue
                gA = best_D(fw_u, 1 + a, gcap_u)
                gB = best_D(bw_1, b, gcap_1)
                grant = max(grant, gA + gB)
            rows.append(dict(u=u, L=L, h=h, J=round(J, 3),
                             s6grant=round(grant, 3),
                             tax=round(grant - J, 3)))
    print(f"tax table (u={u}, from {args.link}):", flush=True)
    print(f"  {'h':>3} {'L':>2} {'J':>7} {'s6grant':>8} "
          f"{'tax':>6}", flush=True)
    for r in rows:
        print(f"  {r['h']:>3} {r['L']:>2} {r['J']:>7} "
              f"{r['s6grant']:>8} {r['tax']:>6}", flush=True)
    # T2-form: certified link transient over the species rate
    T_link = -1e9
    T_wit = None
    for (h, L, dd), g in tabL.items():
        t = (2 * h - g / 4) - (6 / 7) * (h - L)
        if t > T_link:
            T_link, T_wit = t, (h, L, dd, g)
    print(f"T_link = max certified [D - (6/7)(h-L)] = "
          f"{T_link:.2f} at (h,L,delta,g)={T_wit}   "
          f"[s6 two-sided transient was ~7.5 per boundary]",
          flush=True)
    out = dict(src=args.link, rows=rows,
               T_link=round(T_link, 3), T_link_witness=T_wit)
    p = DATA / f"s7_tax_u{u}_g{gcap}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate")
    s = sub.add_parser("specimens")
    s.add_argument("--u", type=int, default=1)
    s.add_argument("--gcap", type=int, default=26)
    s.add_argument("--kmax", type=int, default=1)
    s.add_argument("--whcap", type=int, default=14)
    s.add_argument("--hcap", type=int, default=19)
    s.add_argument("--dcap", type=int, default=30)
    s.add_argument("--seeds", type=str, default="")
    s.add_argument("--n", type=int, default=10)
    s = sub.add_parser("tax")
    s.add_argument("--link", type=str, required=True)
    args = ap.parse_args()
    if args.cmd == "validate":
        validate()
    elif args.cmd == "specimens":
        specimens(args)
    else:
        tax(args)


if __name__ == "__main__":
    main()
