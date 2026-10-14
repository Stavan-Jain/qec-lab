#!/usr/bin/env python3
"""A40 S6 — Stage 1a: the sound aperiodic DRIFT definition (cover
lift), its additivity, the species verification, and the
connectivity lemma.

DRIFT (the definition this session banks).  Work on the universal
cover of the x-circle: walk rows are subsets of Z x {blk} (no x-wrap).
For a fragment on rows [t0, t1] define the per-slab anchor
    A_j := min occupied column of rows [j-3, j]  (both blocks),
defined for j in [t0+3, t1] whenever the slab is nonempty, and the
fragment drift
    delta := A_{t1} - A_{t0+3}   (an integer).
Additivity is definitional: glue fragments overlapping in 4 rows
(exit window of one = entry window of the next) and the A-telescope
gives delta(glued) = delta(1) + delta(2).  On drift-periodic loops
(row(y+p) = row(y) + s pointwise on the cover) the anchor difference
over one period equals s for EVERY anchor convention — the loop drift
is gauge-free; open-fragment drift is anchored to min-column.

COVER LIFT of a torus phase: a (l, p, d)-phase given on the torus
lifts to the cover with integer period-drift s == d (mod l); the lift
is admissible for exactly the right s (the wrap-free recurrence is
strictly stronger than the torus one), which DETERMINES the integer
drift mechanically.  Species targets (charter): W7 must measure
s = -2 per 7 rows, TC63 must measure s = +3 per 6 rows, at both
l = 18 and l = 24; the a36 witness and L12 measure s = 0.

CONNECTIVITY LEMMA (the teleport killer, consumed by Stage 2).  The
recurrence E_j evaluated at column c touches exactly the cells
  v1[j]{c, c+1}, v1[j-3]{c-1}, v2[j]{c}, v2[j+1]{c}, v2[j-1]{c+3}:
x-span 4, y-span 4 — so two support cells not jointly inside one
footprint (dx > 4 or dy > 4 for every containing footprint) never
share a constraint.  If the
  support of a cycle v splits as S1 u S2 with no E-footprint meeting
  both, then E(v) = E(S1) + E(S2) on disjoint constraint sets, so S1
  and S2 are separately cycles.  Hence a CLASS-MINIMAL nontrivial
  logical has (footprint-)CONNECTED support: any component is a
  cycle; a trivial component could be subtracted (same class, lighter
  — contradiction), and if all components are nontrivial any single
  one is a lighter nontrivial logical unless there is only one.
  Mechanical checks here: the footprint spans, the splitting on
  fragments and on an embedding torus.
"""
from __future__ import annotations

import json
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
from a40_s5_lightcore import Phase, code_at  # noqa: E402
sys.argv = _argv

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


# ---------------------------------------------------------------------
# the cover fragment
# ---------------------------------------------------------------------

class CoverFragment:
    """Walk rows on the universal cover: rows[j] = (frozenset s1,
    frozenset s2) of INTEGER columns, j in [t0, t1]."""

    def __init__(self, rows, t0):
        self.t0 = t0
        self.rows = rows              # list of (set, set), index j-t0
        self.t1 = t0 + len(rows) - 1

    def s(self, blk, j):
        if j < self.t0 or j > self.t1:
            return frozenset()
        return self.rows[j - self.t0][blk]

    def weight(self):
        return sum(len(a) + len(b) for a, b in self.rows)

    def e_residual(self, j):
        """Support of E_j on the cover (must be empty when E_j is
        fully inside the fragment)."""
        acc = set()

        def add(cols):
            for c in cols:
                if c in acc:
                    acc.remove(c)
                else:
                    acc.add(c)
        add(self.s(0, j))
        add(c - 1 for c in self.s(0, j))          # x^-1 v1[j]
        add(c + 1 for c in self.s(0, j - 3))      # x    v1[j-3]
        add(self.s(1, j))
        add(self.s(1, j + 1))
        add(c - 3 for c in self.s(1, j - 1))      # x^-3 v2[j-1]
        return acc

    def admissible(self):
        """E_j = 0 for every j fully supported: j in [t0+3, t1-1]."""
        return all(not self.e_residual(j)
                   for j in range(self.t0 + 3, self.t1))

    def slab(self, j):
        return sum(len(self.s(b, t)) for b in (0, 1)
                   for t in range(j - 3, j + 1))

    def slabs(self):
        return [self.slab(j) for j in range(self.t0 + 3, self.t1 + 1)]

    def anchor(self, j):
        cols = [c for b in (0, 1) for t in range(j - 3, j + 1)
                for c in self.s(b, t)]
        return min(cols) if cols else None

    def anchors(self):
        return [self.anchor(j) for j in range(self.t0 + 3, self.t1 + 1)]

    def drift(self):
        a = self.anchors()
        assert a and a[0] is not None and a[-1] is not None
        return a[-1] - a[0]

    def extent(self):
        cols = [c for a, b in self.rows for c in a | b]
        return (max(cols) - min(cols) + 1) if cols else 0

    def window_prune_events(self):
        """Fully-visible tooth alignments held > half.  Tooth at
        (cx, cy): blk1 (cx,cy),(cx,cy-1),(cx-3,cy+1); blk2
        (cx,cy),(cx-1,cy),(cx+1,cy+3).  Visible iff rows cy-1..cy+3
        inside [t0, t1]."""
        ev = []
        cols = sorted({c for a, b in self.rows for c in a | b})
        if not cols:
            return ev
        for cy in range(self.t0 + 1, self.t1 - 2):
            for cx in range(cols[0] - 4, cols[-1] + 5):
                cells = [(0, cx, cy), (0, cx, cy - 1),
                         (0, cx - 3, cy + 1),
                         (1, cx, cy), (1, cx - 1, cy),
                         (1, cx + 1, cy + 3)]
                ov = sum(1 for blk, x, y in cells
                         if x in self.s(blk, y))
                if 2 * ov > 6:
                    ev.append((cy, cx, ov))
        return ev

    def subfragment(self, a, b):
        return CoverFragment(self.rows[a - self.t0:b - self.t0 + 1], a)

    def translated(self, dx):
        return CoverFragment(
            [(frozenset(c + dx for c in s1), frozenset(c + dx for c in s2))
             for s1, s2 in self.rows], self.t0)


def lift_phase(ph: Phase, n_periods=3):
    """Lift a torus phase to the cover, determining the INTEGER
    period-drift s == d (mod l) mechanically: try both candidates
    s in {d - l, d} (and all torus rotations to unstick the wrap);
    accept lifts whose 3-period fragment is wrap-free admissible.
    Returns (fragment, s); asserts s unique across accepted lifts."""
    l, p, d = ph.l, ph.p, ph.d
    cands = sorted({d - l, d, d + l}, key=abs)
    accepted = {}
    for s in cands:
        for rot in range(l):
            rows = []
            ok_compact = True
            for y in range(n_periods * p + 4):
                q, r = divmod(y, p)
                r1 = np.roll(ph.R[0][r], rot)
                r2 = np.roll(ph.R[1][r], rot)
                # torus row (rotated); lift reps in [0, l) then + s*q
                s1 = frozenset(int(c) + s * q
                               for c in np.nonzero(r1)[0])
                s2 = frozenset(int(c) + s * q
                               for c in np.nonzero(r2)[0])
                rows.append((s1, s2))
            fr = CoverFragment(rows, 0)
            # no false accepts are possible: a straddling rep
            # assignment breaks the wrap-free recurrence, so any
            # admissible lift is a genuine cover fragment.
            if fr.admissible():
                accepted.setdefault(s, 0)
                accepted[s] += 1
    if not accepted:
        # no compact cover lift: the phase genuinely winds x (e.g.
        # uses both wraps) — a correct detection, not an error.
        return None, None
    assert len(accepted) == 1, (l, p, d, accepted)
    s = next(iter(accepted))
    # rebuild one canonical accepted lift for return
    for rot in range(l):
        rows = []
        for y in range(n_periods * p + 4):
            q, r = divmod(y, p)
            r1 = np.roll(ph.R[0][r], rot)
            r2 = np.roll(ph.R[1][r], rot)
            s1 = frozenset(int(c) + s * q for c in np.nonzero(r1)[0])
            s2 = frozenset(int(c) + s * q for c in np.nonzero(r2)[0])
            rows.append((s1, s2))
        fr = CoverFragment(rows, 0)
        if fr.admissible():
            return fr, s
    raise AssertionError("unreachable")


# ---------------------------------------------------------------------
# species loading
# ---------------------------------------------------------------------

def load_survivor(fname, p, d, weight):
    recs = json.loads((DATA / fname).read_text())["frames"]
    for f in recs:
        if f["p"] == p and f["d"] == d:
            for s in f.get("survivors", []):
                if s["weight"] == weight:
                    return [tuple(t) for t in s["pts"]]
    raise KeyError((fname, p, d, weight))


def load_pilot_survivor(l, p, d, weight, want_light=None):
    surv = json.loads((DATA / f"s4_prune_pilot_l{l}.json").read_text()
                      )["totals"]["survivors_below"]
    for s in surv:
        if s["p"] == p and s["d"] == d and s["weight"] == weight:
            ph = Phase.from_quotient_pts(l, p, d,
                                         [tuple(t) for t in s["pts"]])
            if want_light is not None:
                if (max(ph.slabs()) <= 7) != want_light:
                    continue
            return ph
    raise KeyError((l, p, d, weight, want_light))


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}

    # ---- footprint spans (mechanical, from the Laurent supports) ----
    # E_j at column c touches: v1[j] {c, c-1}(+1 via x^-1: cell c+1?)
    # computed as the shifts used in e_residual: cells contributing to
    # E_j(c): v1[j] at c and c+1 (since x^-1 v1[j] contributes v1[j]
    # (c+1) to E_j(c)); v1[j-3] at c-1; v2[j] at c; v2[j+1] at c;
    # v2[j-1] at c+3.  x-span = [c-1, c+3] width 5; y-span rows
    # j-3..j+1, 5 rows.  Two cells share a constraint only if their
    # dx <= 4 and dy <= 4.
    cells = [(0, 0, 0), (0, 1, 0), (0, -1, -3), (1, 0, 0), (1, 0, 1),
             (1, 3, -1)]
    xs_ = [c for _, c, _ in cells]
    ys_ = [y for _, _, y in cells]
    span = dict(x_span=max(xs_) - min(xs_), y_span=max(ys_) - min(ys_))
    assert span["x_span"] == 4 and span["y_span"] == 4
    out["footprint"] = span
    print(f"footprint spans: dx <= {span['x_span']}, dy <= "
          f"{span['y_span']} => (4,4)-separated supports split",
          flush=True)

    # ---- species lifts --------------------------------------------------
    species = {}
    for name, l, p, d, w, expect in [
            ("W7_l18", 18, 7, 16, 8, -2),
            ("W7_l24", 24, 7, 22, 8, -2),
            ("TC63_l18", 18, 6, 3, 10, 3),
            ("TC63_l24", 24, 6, 3, 10, 3)]:
        fname = ("s5_dense_p7.json" if (l, p) == (18, 7) else
                 "s5_dense_p6.json" if (l, p) == (18, 6) else
                 f"s5_dense_l24p{p}.json")
        pts = load_survivor(fname, p, d, w)
        ph = Phase.from_quotient_pts(l, p, d, pts)
        assert ph.weight == w and ph.verify_recurrence()
        fr, s = lift_phase(ph, n_periods=3)
        assert s == expect, (name, s, expect)
        # per-period anchor drift (gauge-free on loops):
        anchors = fr.anchors()
        per = [anchors[i + p] - anchors[i]
               for i in range(len(anchors) - p)]
        assert all(x == s for x in per), (name, per)
        sl = fr.subfragment(4, 4 + p + 3).slabs()   # one period's slabs
        rec = dict(l=l, p=p, d=d, weight=w, cover_drift_per_period=s,
                   anchor_drift_checks=len(per),
                   one_period_slabs=sl,
                   all_light=bool(max(sl) <= 7),
                   extent_3periods=fr.extent(),
                   window_pruned=bool(fr.window_prune_events()),
                   deficit_per_period=2 * p - w)
        species[name] = rec
        print(f"{name}: cover drift {s}/period (expected {expect}) "
              f"OK; slabs {sl} all_light={rec['all_light']}, "
              f"3-period extent {rec['extent_3periods']}, "
              f"window_pruned {rec['window_pruned']}, deficit/period "
              f"{rec['deficit_per_period']}", flush=True)
    out["species"] = species

    # controls: zero-drift objects
    zero = {}
    from a40_s5_lightcore import load_witness
    wit_pts = load_witness()
    ph_w = Phase(12, 12, 0, wit_pts)
    fr, s = lift_phase(ph_w, n_periods=2)
    assert s is None, "a36 witness must have NO compact cover lift " \
        "(it uses both wraps — gate W)"
    zero["a36_witness"] = dict(compact_cover_lift=False,
                               note="x-winding, both wraps (gate W)")
    q6, _ = quotient_code(18, 6, 0)
    # L12: from the atlas with pts
    from a40_s4_phase_atlas import atlas
    rows, _ = atlas("AB", 6, 12, keep_pts=True)
    cand = next(r for r in rows if r["nontrivial"] and r["weight"] == 12
                and r["extent"] <= 13)
    ph_l = Phase(18, 6, 0, [(blk, c, y) for (c, y, blk) in cand["pts"]])
    fr_l, s_l = lift_phase(ph_l, n_periods=3)
    assert s_l == 0
    zero["L12"] = dict(drift=s_l, extent=fr_l.extent())
    # the l=12 light species for control (c): (4,4) w6 all-light —
    # MEASURED: no compact cover lift (both-wraps x-winding: this IS
    # the witness species; the a36 witness = 3 periods of it).  Its
    # deficit enters the assembly through the WINDING corner, not the
    # cover frontier — the b = 0 admitted term of control (c).
    ph44 = load_pilot_survivor(12, 4, 4, 6, want_light=True)
    fr44, s44 = lift_phase(ph44, n_periods=3)
    assert s44 is None, "expected the witness species to be winding"
    zero["l12_p4d4_w6"] = dict(compact_cover_lift=False,
                               note="both-wraps winding (the witness "
                                    "species); deficit 2/period rides "
                                    "the winding corner at l=12",
                               deficit_per_period=2 * 4 - 6)
    # the l=12 W7 twin (7,10) w8
    ph710 = load_pilot_survivor(12, 7, 10, 8, want_light=True)
    fr710, s710 = lift_phase(ph710, n_periods=3)
    assert s710 == -2
    zero["l12_p7d10_w8"] = dict(drift_per_period=s710,
                                extent=fr710.extent())
    out["zero_and_l12"] = zero
    print(f"a36 witness: NO compact cover lift (x-winding, both "
          f"wraps — matches gate W); L12 drift "
          f"{zero['L12']['drift']}; l=12 (4,4)w6: winding too (the "
          f"witness species, deficit 2/period in the winding "
          f"corner); l=12 (7,10)w8 drift {s710} — the W7 twin",
          flush=True)

    # ---- additivity (mechanical) ---------------------------------------
    # glue overlap-4: fragment rows [0, T1] and [T1-3, T2] agreeing on
    # the shared 4 rows; delta telescopes.  Verify on species: k-period
    # prefixes.
    checks = []
    for name in ("W7_l18", "TC63_l24"):
        rec = species[name]
        l, p, d = rec["l"], rec["p"], rec["d"]
        fname = ("s5_dense_p7.json" if (l, p) == (18, 7) else
                 f"s5_dense_l24p{p}.json")
        ph = Phase.from_quotient_pts(
            l, p, d, load_survivor(fname, p, d, rec["weight"]))
        fr, s = lift_phase(ph, n_periods=3)
        # fragment A = rows [0, 4+p], B = rows [1+p, 4+2p]:
        # exit window of A (rows [1+p, 4+p]) = entry window of B.
        A = fr.subfragment(0, 4 + p)
        B = fr.subfragment(1 + p, 4 + 2 * p)
        AB = fr.subfragment(0, 4 + 2 * p)
        dA, dB, dAB = A.drift(), B.drift(), AB.drift()
        assert dA + dB == dAB, (name, dA, dB, dAB)
        checks.append(dict(name=name, dA=dA, dB=dB, dAB=dAB))
        # independence of gauge on loops: anchor drift per period
        # equals s for both the min-column anchor and a max-column
        # anchor variant:
        amax = [max(c for b_ in (0, 1) for t in range(j - 3, j + 1)
                    for c in fr.s(b_, t))
                for j in range(fr.t0 + 3, fr.t1 + 1)]
        assert all(amax[i + p] - amax[i] == s
                   for i in range(len(amax) - p))
    out["additivity_checks"] = checks
    print(f"additivity: {len(checks)} glue checks pass "
          f"(delta(A)+delta(B) = delta(AB)); loop drift gauge-free "
          f"(min- and max-anchor agree per period)", flush=True)

    # ---- splitting / connectivity (mechanical) --------------------------
    # far-apart union of two admissible fragments is admissible, and
    # the components are separately cycles on an embedding torus.
    frW, _ = lift_phase(Phase.from_quotient_pts(
        18, 7, 16, load_survivor("s5_dense_p7.json", 7, 16, 8)), 2)
    frT, _ = lift_phase(Phase.from_quotient_pts(
        18, 6, 3, load_survivor("s5_dense_p6.json", 6, 3, 10)), 2)
    off = frW.extent() + 40
    T = min(frW.t1, frT.t1)
    un = CoverFragment(
        [(frW.rows[i][0] | frozenset(c + off for c in frT.rows[i][0]),
          frW.rows[i][1] | frozenset(c + off for c in frT.rows[i][1]))
         for i in range(T + 1)], 0)
    assert un.admissible(), "far union must be admissible"
    # component recovery under (4,4)-adjacency:
    pts = [(b, c, j) for j in range(T + 1) for b in (0, 1)
           for c in un.s(b, j)]
    comp = {q: q for q in pts}

    def find(q):
        while comp[q] != q:
            comp[q] = comp[comp[q]]
            q = comp[q]
        return q
    for i_, q in enumerate(pts):
        for q2 in pts[i_ + 1:]:
            if abs(q[1] - q2[1]) <= 4 and abs(q[2] - q2[2]) <= 4:
                comp[find(q)] = find(q2)
    ncomp = len({find(q) for q in pts})
    assert ncomp == 2, f"expected 2 components, got {ncomp}"
    out["splitting"] = dict(union_admissible=True, n_components=2)
    print("splitting: far union admissible, exactly 2 components "
          "under (4,4)-adjacency — each separately admissible by "
          "construction; connectivity lemma's mechanical side PASS",
          flush=True)

    # torus version of the splitting: two far x-placements of the
    # SAME closed cycle on an embedding torus — the union is a cycle
    # and each component is separately a cycle.
    big = code_at((48, 6), "embed-48x6")
    v_one = np.zeros(big.n, dtype=np.uint8)
    v_two = np.zeros(big.n, dtype=np.uint8)
    for (c, y, blk) in cand["pts"]:      # the L12 object, extent <= 13
        v_one[blk * big.ng + big.G.index((c % 48, y % 6))] ^= 1
        v_two[blk * big.ng + big.G.index((c % 48, y % 6))] ^= 1
        v_two[blk * big.ng + big.G.index(((c + 24) % 48, y % 6))] ^= 1
    assert big.is_cycle(v_one) and big.is_cycle(v_two)
    assert big.is_cycle(v_one ^ v_two)   # the far component alone
    out["torus_splitting"] = True
    print("torus splitting: far two-placement union is a cycle and "
          "each component is separately a cycle", flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s6_drift.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s6_drift.json'} ({out['wall_s']} s)", flush=True)


if __name__ == "__main__":
    main()
