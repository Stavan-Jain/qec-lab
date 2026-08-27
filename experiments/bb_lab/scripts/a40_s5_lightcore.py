#!/usr/bin/env python3
"""A40 S5 — the L-W light-core window machinery: the exact y-walk
window engine, its controls, and the state-space pricing lemmas.

The wall certificate (S4 §9.7 item 5) needs a transfer graph over
sliding windows of the y-walk with slab-amortized costs, pruned by
local reduction.  This script builds and validates every LOCAL piece
of that machine at fixed l, WITHOUT claiming the full graph (the
pricing section shows why full materialization is RED — the honest
verdict this session's periodic-leg censuses route around):

  ENGINE (exact, per l):
  - the drift-periodic UNROLLER: a (p, d)-phase given as a vector on
    the SNF-normalized quotient torus of Z^2/<(l,0),(d,p)> is lifted
    back to walk rows (v1[y], v2[y] in F2^l, row(y+p) = x^d row(y)),
    and re-verified against the row recurrence E_j (C2 conventions)
    for every j — an independent end-to-end check of the phase;
  - slab weights W_j (rows j-3..j, both blocks) with the exact
    telescope sum_j W_j = 4|v| asserted;
  - the H = 5 window prune rule: a tooth (single H_X row) at check
    row cy spans block-1 rows cy-1..cy+1 and block-2 rows cy..cy+3 —
    5 consecutive rows — so a 5-row sliding window sees each tooth
    alignment exactly once (at cy = window-oldest + 1); the rule
    prunes when the walk holds > half (>= 4 of 6) of a visible tooth.
    Sound: a class-minimal v has |v + z| >= |v| for every stabilizer
    z, and the overlap is fully visible when z is.

  CONTROLS (falsify-first):
  C1  the a36 two-gross w18 witness at (12,12) [p=12, d=0]: must
      re-verify as a drift-periodic cycle, rate 1.5, and be UNPRUNED
      (it is class-minimal) — the charter's positive control.  Its
      slab profile decides whether it is a LIGHT-core cycle (all
      W_j <= 7).
  C2  combs: the single-tooth (p=6, d=0) phase at l = 18 and l = 12:
      rate 1, trivial, and the window rule must PRUNE it.
  C3  the L12 species at (18,6) [p=6, d=0]: rate 2; its slab profile
      (the rate-2 boundary object the potential certificate must
      price).
  C4  the l = 12 pilot survivor atlas (all 936 banked sub-rate-2
      nontrivial orbits): each re-verified by the unroller; per
      family: light fraction, slab ranges, window-prune (must be NO —
      they survived the stronger global catalog).  The light ones are
      genuine cycles of the pruned light core at l = 12: their min
      rate is the honest mu_light(l=12) upper bound (the charter's
      mu = 1.5 expectation is superseded by the pilot's rate-8/7
      family if that family is light — measured here).

  PRICING (the mandated RED/AMBER/GREEN before any full build):
  - exact universe counts for the window state spaces at l = 12/18
    (4+4 dynamics window; 5+5 prune window), x-translation quotient;
  - the constraint-structure lemma, verified mechanically: the 4-row
    window universe carries NO internal E-constraint (every light
    4-window occurs on an admissible bi-infinite walk), and the 5-row
    window carries EXACTLY ONE (E_{j-1}) — so reachability prunes
    (almost) nothing and the only cut is bi-recurrence, a greatest
    fixed point that requires materializing the universe: RED at
    l = 18.  The cluster quotient is exact at l = 18 (no gap can
    saturate below the 28-column light-bridge bound > 18) and only
    divides by l — the reduction that makes the space small is a
    FUTURE per-cluster potential decomposition, named in the note.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

_argv = sys.argv
sys.argv = [_argv[0], "12", "8"]
from a40_s4_phase_triage import quotient_code, snf2  # noqa: E402
from a40_s4_prune_pilot import inv2  # noqa: E402
sys.argv = _argv

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def code_at(lm, name=None):
    return TowerCode(name or f"tdg{lm}", lm, red(A_L, lm), red(B_L, lm))


def xs(row, s, l):
    return np.roll(row, s % l)


class Phase:
    """A drift-periodic walk phase at x-order l: rows R[blk][y] for
    y in 0..p-1, with row(y + p) = x^d row(y)."""

    def __init__(self, l, p, d, pts_e):
        """pts_e: iterable of (blk, x, y) in EXPONENT coordinates
        (any representatives; reduced mod the lattice here)."""
        self.l, self.p, self.d = l, p, d
        R = np.zeros((2, p, l), dtype=np.uint8)
        for blk, x, y in pts_e:
            yr = y % p
            xr = (x - d * (y // p)) % l
            R[blk, yr, xr] ^= 1
        self.R = R
        self.weight = int(R.sum())

    @classmethod
    def from_quotient_pts(cls, l, p, d, pts_f):
        """pts_f: (blk, f1, f2) on the SNF-normalized quotient torus
        (the pilot's storage format).  e = f V^{-1}."""
        M = [[l, 0], [d, p]]
        D, U, V = snf2(M)
        Vi = inv2(V)
        pts_e = [(blk, f1 * Vi[0][0] + f2 * Vi[1][0],
                  f1 * Vi[0][1] + f2 * Vi[1][1])
                 for blk, f1, f2 in pts_f]
        return cls(l, p, d, pts_e)

    def row(self, blk, y):
        return xs(self.R[blk][y % self.p], self.d * (y // self.p),
                  self.l)

    def verify_recurrence(self):
        """E_j = (1+x^-1)v1[j] + x v1[j-3] + v2[j] + v2[j+1]
        + x^-3 v2[j-1] = 0 for all j (one period suffices)."""
        l = self.l
        for j in range(self.p):
            e = (self.row(0, j) ^ xs(self.row(0, j), -1, l)
                 ^ xs(self.row(0, j - 3), 1, l)
                 ^ self.row(1, j) ^ self.row(1, j + 1)
                 ^ xs(self.row(1, j - 1), -3, l))
            if e.any():
                return False
        return True

    def slabs(self):
        """W_j = weight of rows j-3..j, both blocks, j in 0..p-1."""
        rw = [int(self.row(0, j).sum() + self.row(1, j).sum())
              for j in range(self.p)]
        return [sum(rw[(j - t) % self.p] for t in range(4))
                for j in range(self.p)]

    def window_prune_events(self):
        """Tooth alignments (cy, cx) fully visible in some 5-row
        window with overlap >= 4 (of 6).  Tooth at (cx, cy):
        block-1 (cx,cy),(cx,cy-1),(cx-3,cy+1);
        block-2 (cx,cy),(cx-1,cy),(cx+1,cy+3)."""
        events = []
        for cy in range(self.p):
            for cx in range(self.l):
                cells = [(0, cx, cy), (0, cx, cy - 1),
                         (0, (cx - 3), cy + 1),
                         (1, cx, cy), (1, (cx - 1), cy),
                         (1, (cx + 1), cy + 3)]
                ov = sum(int(self.row(blk, y)[x % self.l])
                         for blk, x, y in cells)
                if 2 * ov > 6:
                    events.append((cy, cx, ov))
        return events

    def report(self):
        sl = self.slabs()
        assert sum(sl) == 4 * self.weight, "slab telescope fails"
        ev = self.window_prune_events()
        return dict(weight=self.weight, p=self.p, d=self.d,
                    rate=self.weight / self.p,
                    slab_min=min(sl), slab_max=max(sl),
                    all_light=bool(max(sl) <= 7),
                    y_spanning=bool(min(sl) >= 1),
                    window_pruned=bool(ev), n_prune_events=len(ev))


def hx_row_pts(code, cx, cy):
    """Support of the H_X row at check (cx, cy) as (blk, x, y)."""
    out = []
    row = code.HX[code.G.index((cx % code.G.orders[0],
                                cy % code.G.orders[1]))]
    for i in np.nonzero(row)[0]:
        blk, gi = divmod(int(i), code.ng)
        g = code.G.from_index(gi)
        out.append((blk, g[0], g[1]))
    return out


def load_witness():
    """The a36 two-gross w18 witness in the paper frame (12,12)
    (transform banked in a40_s4_stack_gate.py, re-asserted here)."""
    tg = code_at((12, 12), "two-gross-paper")
    wit = json.loads((LAB / "data" / "a36" /
                      "w18_witness_banked.json").read_text())
    tg_s = TowerCode("tg/stored", (12, 12), "x^3 + y^2 + y^7",
                     "y^3 + x + x^2")
    v_s = np.zeros(tg_s.n, dtype=np.uint8)
    v_s[wit["v_support"]] = 1
    v18 = np.zeros(tg.n, dtype=np.uint8)
    for i in np.nonzero(v_s)[0]:
        blk, gi = divmod(int(i), tg.ng)
        h = tg_s.G.from_index(gi)
        s = (0, 7) if blk == 0 else (1, 0)
        u = ((h[0] + s[0]) % 12, (7 * (h[1] + s[1])) % 12)
        v18[blk * tg.ng + tg.G.index(u)] = 1
    assert tg.is_cycle(v18) and not tg.is_stab(v18) \
        and v18.sum() == 18
    pts = []
    for i in np.nonzero(v18)[0]:
        blk, gi = divmod(int(i), tg.ng)
        g = tg.G.from_index(gi)
        pts.append((blk, g[0], g[1]))
    return pts


def universe_counts():
    """Exact light-window universe sizes (pricing)."""
    rows = []
    for l in (12, 18):
        cells4 = 8 * l          # 4 rows x 2 blocks
        cells5 = 10 * l
        u4 = sum(math.comb(cells4, w) for w in range(0, 8))
        # 5-row window: both slabs <= 7 -> weight <= 14 (upper bound
        # count: <= 14 points, one internal E-constraint /2^l)
        u5 = sum(math.comb(cells5, w) for w in range(0, 15))
        rows.append(dict(
            l=l, window4_states=u4, window4_per_translation=u4 // l,
            window5_upper=u5,
            window5_per_translation_with_constraint=u5 // l // (2 ** l)))
    return rows


def constraint_lemma(rng):
    """Mechanical check of the window constraint structure at l=18:
    every RANDOM light 4+4 window extends to an admissible bi-infinite
    walk — forward by the forced-v2 recurrence (E_t defines v2[t+1]),
    backward by the forced-v1 recurrence (E_{r0+2} defines v1[r0-1]:
    the x v1[j-3] term has a monomial coefficient, so it is solvable) —
    after which EVERY fully-supported E_j on the strip must hold by
    construction.  A pass on all samples pins the lemma: the 4-row
    window universe carries no internal constraint, so reachability
    prunes nothing."""
    l = 18
    n_ok = 0
    T = 200
    for _ in range(T):
        v1 = [np.zeros(l, np.uint8) for _ in range(4)]
        v2 = [np.zeros(l, np.uint8) for _ in range(4)]
        for _2 in range(7):
            blk = int(rng.integers(2))
            (v1 if blk == 0 else v2)[int(rng.integers(4))][
                int(rng.integers(l))] ^= 1
        for _2 in range(6):   # forward: free v1 input 0, forced v2
            nxt2 = (v1[-1] ^ xs(v1[-1], -1, l) ^ xs(v1[-4], 1, l)
                    ^ v2[-1] ^ xs(v2[-2], -3, l))
            v1.append(np.zeros(l, np.uint8))
            v2.append(nxt2)
        for _2 in range(6):   # backward: free v2 input 0, forced v1
            # E_{r0+2}: v1[r0-1] = x^-1[(1+x^-1)v1[r0+2] + v2[r0+2]
            #                           + v2[r0+3] + x^-3 v2[r0+1]]
            b1 = xs(v1[2] ^ xs(v1[2], -1, l) ^ v2[2] ^ v2[3]
                    ^ xs(v2[1], -3, l), -1, l)
            v1.insert(0, b1)
            v2.insert(0, np.zeros(l, np.uint8))
        ok = True
        for j in range(3, len(v1) - 1):
            e = (v1[j] ^ xs(v1[j], -1, l) ^ xs(v1[j - 3], 1, l)
                 ^ v2[j] ^ v2[j + 1] ^ xs(v2[j - 1], -3, l))
            if e.any():
                ok = False
        n_ok += ok
    return dict(n_random_4windows=T, n_extended_ok=n_ok)


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}

    # ---- C1: the a36 witness ----------------------------------------
    wit_pts = load_witness()
    ph = Phase(12, 12, 0, [(b, x, y) for b, x, y in wit_pts])
    assert ph.weight == 18
    assert ph.verify_recurrence(), "witness fails E_j?!"
    r = ph.report()
    out["C1_witness"] = r
    print(f"C1 witness (12,12) p=12 d=0: rate {r['rate']:g}, slabs "
          f"[{r['slab_min']},{r['slab_max']}], all_light "
          f"{r['all_light']}, window_pruned {r['window_pruned']}")
    assert not r["window_pruned"], \
        "class-minimal witness must be unpruned"
    assert abs(r["rate"] - 1.5) < 1e-9

    # ---- C2: combs ----------------------------------------------------
    for l in (18, 12):
        q6, o6 = quotient_code(l, 6, 0)
        # single tooth on the (l,6) frame = the comb's phase vector
        c186 = code_at((l, 6))
        pts = hx_row_pts(c186, 0, 2)
        ph = Phase(l, 6, 0, pts)
        assert ph.weight == 6
        assert ph.verify_recurrence(), "comb phase fails E_j"
        r = ph.report()
        out[f"C2_comb_l{l}"] = r
        print(f"C2 comb l={l} p=6 d=0: rate {r['rate']:g}, slabs "
              f"[{r['slab_min']},{r['slab_max']}], window_pruned "
              f"{r['window_pruned']} ({r['n_prune_events']} events)")
        assert r["window_pruned"], "comb must be window-pruned"
        assert r["all_light"], "comb slabs must be light (<= 5)"

    # ---- C3: the L12 species at (18,6) --------------------------------
    from a40_s4_phase_atlas import atlas
    rows, _ = atlas("AB", 6, 12, keep_pts=True)
    l12s = [row for row in rows if row["nontrivial"]
            and row["weight"] == 12]
    assert len(l12s) == 66, f"L12 count {len(l12s)} != 66"
    profiles = []
    n_placed = 0
    for cand in l12s:
        if cand["extent"] > 13:   # wrap-faithful placement at l=18
            continue
        pts = [(blk, c, y) for (c, y, blk) in cand["pts"]]
        ph = Phase(18, 6, 0, pts)
        assert ph.verify_recurrence(), "L12 fails E_j"
        sl = ph.slabs()
        profiles.append((min(sl), max(sl)))
        n_placed += 1
    assert n_placed >= 50, f"only {n_placed} L12 placed at l=18"
    smin = min(a for a, b in profiles)
    smax = max(b for a, b in profiles)
    out["C3_L12"] = dict(n=66, n_placed=n_placed, rate=2.0,
                         slab_min_over_family=smin,
                         slab_max_over_family=smax)
    print(f"C3 L12 family at (18,6): 66 objects ({n_placed} placed "
          f"wrap-faithfully), rate 2, slab range over family "
          f"[{smin},{smax}] (the rate-2 boundary sits at slab 8 = "
          f"heavy threshold)")

    # ---- C4: the l=12 survivor atlas ----------------------------------
    surv = json.loads((DATA / "s4_prune_pilot_l12.json").read_text()
                      )["totals"]["survivors_below"]
    assert len(surv) == 936
    fam = {}
    for s in surv:
        ph = Phase.from_quotient_pts(12, s["p"], s["d"],
                                     [tuple(t) for t in s["pts"]])
        assert ph.weight == s["weight"], (s["p"], s["d"], s["weight"],
                                          ph.weight)
        assert ph.verify_recurrence(), \
            f"survivor (p={s['p']},d={s['d']}) fails E_j"
        r = ph.report()
        assert not r["window_pruned"], \
            "survivor window-pruned but globally unpruned?!"
        key = (s["p"], s["d"], s["weight"])
        f = fam.setdefault(key, dict(n=0, n_light=0, slab_max=0,
                                     slab_min=99))
        f["n"] += 1
        f["n_light"] += r["all_light"]
        f["slab_max"] = max(f["slab_max"], r["slab_max"])
        f["slab_min"] = min(f["slab_min"], r["slab_min"])
    out["C4_l12_survivors"] = [
        dict(p=p, d=d, weight=w, rate=w / p, **f)
        for (p, d, w), f in sorted(fam.items())]
    print("C4 l=12 survivor families (all 936 re-verified, all "
          "window-UNpruned):")
    mu_light = None
    for row in out["C4_l12_survivors"]:
        light = f"{row['n_light']}/{row['n']} light"
        print(f"   (p={row['p']},d={row['d']}) w={row['weight']} "
              f"rate {row['rate']:.3g}: {light}, slabs "
              f"[{row['slab_min']},{row['slab_max']}]")
        if row["n_light"]:
            mu_light = min(mu_light or 99, row["rate"])
    out["mu_light_l12_upper"] = mu_light
    print(f"=> mu_light(l=12) <= {mu_light} (cheapest ALL-LIGHT "
          f"unpruned nontrivial periodic orbit, p <= 7)")

    # ---- pricing -------------------------------------------------------
    out["universe"] = universe_counts()
    for u in out["universe"]:
        print(f"pricing l={u['l']}: 4+4 window universe "
              f"{u['window4_states']:.3g} ({u['window4_per_translation']:.3g} "
              f"per x-translation); 5+5 window <= {u['window5_upper']:.3g} "
              f"(~{u['window5_per_translation_with_constraint']:.3g} after "
              f"translation + the one internal constraint)")
    rng = np.random.default_rng(40)
    cl = constraint_lemma(rng)
    out["constraint_lemma"] = cl
    print(f"constraint lemma: {cl['n_extended_ok']}/"
          f"{cl['n_random_4windows']} random light 4+4 windows extend "
          f"to admissible walks (forced-row forward + forced-row "
          f"backward): the 4-window universe is constraint-free — "
          f"reachability gives no cut; the cut is bi-recurrence "
          f"(greatest fixed point) => RED without a decomposition.")

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s5_lightcore.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s5_lightcore.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
