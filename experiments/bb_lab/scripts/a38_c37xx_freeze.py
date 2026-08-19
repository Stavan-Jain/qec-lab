"""A38 S2: the [[720,4]] freeze certification (A35 SS8.2, GREEN, unexecuted).

Tower c37xx (A35 docket row; every level's polynomials by literal descent):

  L0 (60,6) [[720,4,?]]  <--x30--  the re-double of c37x along x
  L1 (30,6) [[360,4,20]] = A30 doubled code 37a70e02:x — d = 20
           END-TO-END CERTIFICATE TIER (A30, banked
           data/a30/decide_37a70e02.json + rungs_37a70e02_x.json)
  L2 (15,6) [[180,4,10]] = 37a70e02 — d = 10 certificate tier (same bank)
  L3 (15,3) [[90,4,?]]   — d(L3) measured exactly here (free by-product)

  A = 1 + y + x, B = y^4 + x + x^11*y^2 at every level (mod the level's
  orders); rungs x-U-(R), x-U-(R), y-T-(R); pairs (0,1) = R2 (S cap K
  = 0), (1,2) = R3 (S1 = ker p2*) — both asserted numerically here.

Questions (A14 SS13 freeze-vs-double, first re-double tower):

  Q1 (W = 18): no nontrivial L0-logical of weight <= 18  =>  d >= 20
      ("d = 20 if frozen", floor half).
  Q2 (W = 22): the same stack pushed to 22 with rungs at M = (24-|b|)/2,
      BOTH DIRECTIONS ARMED: all PASS => d >= 24 (freeze REFUTED);
      any find is an explicit nontrivial logical of weight 20/22,
      re-verified in-line => min find = d exactly (with Q1's floor).

Architecture (the A32/A33/A36 stack, two fiber layers composed; ALL
BZ censuses at n <= 180, every cap-8 MITM at the demonstrated n = 90
scale):

  L3: coset-BZ censuses (stab <= W; all 15 nonzero classes <= WC).
  L2 (descent): stab reps <= W and S1-coset reps <= W from L3-stab
      fibers + tau3-families; all-class cycle reps <= WC (the
      kernel-shift windows) additionally from L3-coset fibers.
      FALSIFY-FIRST GATE: a direct 4-offset BZ census at L2 (<= 14)
      must equal the descent-derived slice exactly (orbit key sets).
  L1 (descent): stab reps <= W from L2-stab fibers + tau2-families;
      SEAM-coset reps <= W (Q2 only; at W = 18 the seam branch is DEAD
      outright: seam elements are nontrivial L1-logicals, weight
      >= d(L1) = 20 > 18 by the A30 certificate).
  Rungs L1 -> L0 (library RungCell + the kernel-shift lane): dangerous
      per stab orbit rep at M = (target-|b|)/2, seam per orbit rep at
      the same M; verdicts transport along G (covariance spot-checked).

The kernel-shift lane (the F2b epsilon-recursion in executable form):
the solution set of the carry system E v0 = RHS b is v0p + ker E, and
ker E = Z(base) EXACTLY (tau is injective on 0-chains), so every
candidate with overflow <= M-1 is v0p (+) z with
|z| <= |b| + (M-1) + ov(v0p).  Choosing v0p = the row-decomposition
lift (b = sum of base HX rows => v0p = the same sum of cover rows)
makes the bound small for light b, and when it lands BELOW
d(base) = 20 every shift z is a BASE STABILIZER inside the
already-built census — the rung consumes level-(r-1) censuses instead
of a cap-8 MITM at n = 360 (C(354,4) ~ 6e8 subsets, ~300x the
demonstrated cap-8 instance; the A35 cap gate is n-blind — measured
here).  The same trick replaces the cap-7/8 fibers at n = 180.
Cross-validated in-run against the direct restricted lane on every
cell where both run.

Discipline: no SAT anywhere; witness weights are upper bounds only;
every census carries exact node-count asserts; d(L1) = 20 and
d(L2) = 10 are consumed at A30 certificate tier and TRIPWIRED in-run
(any nontrivial L1-cycle <= 18 or L2-cycle < 10 found anywhere
contradicts the bank => hard assert).

Usage: python a38_c37xx_freeze.py 18                  (Q1, one shot)
       python a38_c37xx_freeze.py 22 --census-only    (Q2 phase 1)
       python a38_c37xx_freeze.py 22 --rungs-only     (Q2 phase 2)
       python a38_c37xx_freeze.py 22                  (Q2, one shot)

The two-phase Q2 split exists because the harness kills background
tasks at ~1 h wall: phase 1 runs the censuses + gates and CHECKPOINTS
the L1 obligations (orbit-rep supports, data/a38/c37xx/ckpt_W22_*);
phase 2 reloads them (re-verifying every vector: weight, stab/seam
membership, class) and runs the rungs + assembly.

Output: data/a38/c37xx/freeze_W{W}.json + rungs_W{W}.jsonl
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    AxisDeck, RungCell, TowerCode, _as_support, batch_keys, colspace,
    enumerate_lifts_deep, fold_support, h1_map, i2v, in_span,
    kernel_ints, perm_for, rep_for, rref_ints, span_eq, span_points,
    translation_perms, v2i,
)

DATA = LAB / "data" / "a38" / "c37xx"
DATA.mkdir(parents=True, exist_ok=True)
A30 = LAB / "data" / "a30"

D_L1 = 20     # A30 certificate tier ([[360,4,20]], decide_37a70e02.json)
D_L2 = 10     # A30 certificate tier ([[180,4,10]], same bank)

SPEC = ((60, 6), "1 + y + x", "y^4 + x + x^11*y^2")
FOLDS = [(0, 30), (0, 15), (1, 3)]


# --------------------------------------------------------------- helpers
class Collector:
    """Collect vectors, dedupe to translation-orbit reps at the end."""

    def __init__(self, n: int):
        self.n = n
        self.vecs: list[np.ndarray] = []
        self.seen: set[int] = set()

    def add(self, v: np.ndarray) -> None:
        vi = v2i(v)
        if vi not in self.seen:
            self.seen.add(vi)
            self.vecs.append(v)

    def reps(self, perms) -> list[np.ndarray]:
        if not self.vecs:
            return []
        keys = batch_keys(np.array(self.vecs, dtype=np.uint8), perms)
        out: dict[bytes, np.ndarray] = {}
        for i, k in enumerate(keys):
            out.setdefault(bytes(k), self.vecs[i])
        return list(out.values())


def min_stab_decompose(code: TowerCode, b: np.ndarray) -> list[int]:
    """Row indices I with b = sum_{i in I} HX[i]; tries |I| <= 3 by
    direct search (weight-guided), falls back to the tracked-RREF
    decomposition."""
    wb = int(b.sum())
    bi = v2i(b)
    rows_i = [v2i(code.HX[i]) for i in range(code.ng)]
    if wb == 6:
        for i, ri in enumerate(rows_i):
            if ri == bi:
                return [i]
    if wb <= 12:
        by_val = {}
        for i, ri in enumerate(rows_i):
            by_val.setdefault(ri ^ bi, []).append(i)
        for i, ri in enumerate(rows_i):
            for j in by_val.get(ri, []):
                if j > i:
                    return [i, j]
    if wb <= 18:
        by_val2 = {}
        for i, ri in enumerate(rows_i):
            by_val2.setdefault(ri, []).append(i)
        for i, j in combinations(range(code.ng), 2):
            tgt = rows_i[i] ^ rows_i[j] ^ bi
            for k in by_val2.get(tgt, []):
                if k > j:
                    return [i, j, k]
    # fallback: tracked RREF
    basis: list[int] = []
    piv: list[int] = []
    hist: list[int] = []
    for i in range(code.ng):
        cur, h = rows_i[i], 1 << i
        for bb, pp, hh in zip(basis, piv, hist):
            if (cur >> pp) & 1:
                cur ^= bb
                h ^= hh
        if cur:
            basis.append(cur)
            piv.append((cur & -cur).bit_length() - 1)
            hist.append(h)
    x, hsel = bi, 0
    for bb, pp, hh in zip(basis, piv, hist):
        if (x >> pp) & 1:
            x ^= bb
            hsel ^= hh
    assert x == 0, "b is not a stabilizer"
    return [i for i in range(code.ng) if (hsel >> i) & 1]


def row_lift_v0(deck: AxisDeck, b: np.ndarray) -> tuple[np.ndarray, int]:
    """(v0p, ov) for the row-decomposition lift of the base stabilizer
    b: v0p = sheet-0 part of the sum of sheet-0-lifted cover rows."""
    I = min_stab_decompose(deck.base, b)
    vc = np.zeros(deck.cover.n, dtype=np.uint8)
    for i in I:
        g = deck.base.G.from_index(i)
        gi = deck.cover.G.index(deck.cover.G.reduce(deck.emb(g, 0)))
        vc ^= deck.cover.HX[gi]
    assert ((deck.P @ vc) % 2 == b).all(), "row lift shadow mismatch"
    v0p = deck.sheets(vc)[0]
    ov = int((v0p.astype(bool) & ~b.astype(bool)).sum())
    return v0p, ov


class KernelShift:
    """Kernel-shift candidate enumeration for one deck: candidates over
    a shadow b are v0p (+) z with z in the base-cycle census window
    (all G(base)-translates of the given orbit reps; the census must be
    COMPLETE to `complete_to`)."""

    def __init__(self, deck: AxisDeck, reps: list[np.ndarray],
                 complete_to: int):
        self.deck = deck
        self.complete_to = complete_to
        perms = translation_perms(deck.base)
        self.by_w: dict[int, list[int]] = {}
        seen: set[int] = set()
        for rep in reps:
            w = int(rep.sum())
            bucket = self.by_w.setdefault(w, [])
            for p in perms:
                zi = v2i(rep[p])
                if zi not in seen:
                    seen.add(zi)
                    bucket.append(zi)

    def candidates(self, b: np.ndarray, v0p: np.ndarray, cap: int):
        """Yields all solution ints v0 = v0p ^ z with
        |v0 off supp(b)| <= cap.  Sound iff B = |b| + cap + ov(v0p)
        <= complete_to (asserted)."""
        bmask = v2i(b)
        v0p_i = v2i(v0p)
        ovp = bin(v0p_i & ~bmask).count("1")
        B = int(b.sum()) + cap + ovp
        assert B <= self.complete_to, \
            f"kernel-shift window {B} > census bound {self.complete_to}"
        if bin(v0p_i & ~bmask).count("1") <= cap:
            yield v0p_i
        for w, bucket in self.by_w.items():
            if w <= B:
                for zi in bucket:
                    v0 = v0p_i ^ zi
                    if bin(v0 & ~bmask).count("1") <= cap:
                        yield v0


def census_pass(binp, code: TowerCode, offsets, W: int, tag: str,
                deadline_s: float = 3600.0) -> dict[str, set[int]]:
    """One complete two-window coset-BZ pass; {label -> element ints}.
    Node counts exact (asserted inside run_window)."""
    I1, G1, I2, G2, kappa = cosetbz.disjoint_info_sets(code.HX)
    assert not (set(I1) & set(I2))
    r1, r2 = W // 2, max(W - W // 2 - 1, 0)
    hits: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = []
        for lab, tv in offsets:
            cb = cosetbz.coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            if 0 < wcb <= W:      # the empty-window coset-base edge
                hits[lab].add(v2i(cb))
            bases.append(cb)
        res = cosetbz.run_window(binp, f"{tag}_w{wi}", Gs, bases, r, W,
                                 time.monotonic() + deadline_s)
        for j, hx in res.pop("hit_rows"):
            v = cosetbz.unpack3(hx, code.n)
            if v.any():
                hits[offsets[j][0]].add(v2i(v))
    return hits


def whist(vs) -> dict[str, int]:
    h: dict[int, int] = {}
    for v in vs:
        w = int(v.sum())
        h[w] = h.get(w, 0) + 1
    return {str(k): v for k, v in sorted(h.items())}


def main() -> None:
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    assert W in (18, 22)
    census_only = "--census-only" in sys.argv
    rungs_only = "--rungs-only" in sys.argv
    assert not (census_only and rungs_only)
    TARGET = W + 2
    WC = 16 if W == 18 else 18       # all-class window-census bound
    t0 = time.monotonic()
    out: dict = {"W": W, "target": TARGET, "WC": WC,
                 "phase": ("census" if census_only else
                           "rungs" if rungs_only else "full")}

    def log(msg: str) -> None:
        print(f"[{time.monotonic()-t0:7.1f}s] {msg}", flush=True)

    # ------------------------------------------------ 0. tower + structure
    orders, As, Bs = SPEC
    G_top = AbelianGroup(orders)
    levels = [(orders, _as_support(As, G_top), _as_support(Bs, G_top))]
    for axis, newmod in FOLDS:
        plm, pA, pB = levels[-1]
        nlm = tuple(newmod if a == axis else plm[a]
                    for a in range(len(plm)))
        levels.append((nlm, fold_support(pA, axis, newmod),
                       fold_support(pB, axis, newmod)))
    L = [TowerCode(f"L{i}", *lv) for i, lv in enumerate(levels)]
    assert [c.k for c in L] == [4, 4, 4, 4]
    assert [c.n for c in L] == [720, 360, 180, 90]
    deck0 = AxisDeck(L[0], L[1], 0)
    deck1 = AxisDeck(L[1], L[2], 0)
    deck2 = AxisDeck(L[2], L[3], 1)
    bank = json.loads((A30 / "decide_37a70e02.json").read_text())
    assert bank["group"] == [15, 6] and bank["d_base"] == D_L2 \
        and bank["floor"] == D_L1
    G2b = AbelianGroup((15, 6))
    assert L[2].A.support == Poly.from_string(bank["A"], G2b).support
    assert L[2].B.support == Poly.from_string(bank["B"], G2b).support
    log(f"tower built: n = {[c.n for c in L]}, k = 4/4/4/4; L2 == banked "
        f"A30 37a70e02 (d(L2) = {D_L2}, d(L1) = {D_L1}, certificate tier)")
    for c in L:
        assert not any(int(kv.sum()) % 2 for kv in c.kerHZ)
    Mp0, Mp1, Mp2 = h1_map(deck0), h1_map(deck1), h1_map(deck2)
    SEAM = colspace(Mp0)
    Sbb, Sbp = rref_ints(list(SEAM))
    seam_set = span_points(Sbb) - {0}
    K1 = kernel_ints(Mp1)
    S1 = colspace(Mp1)
    K2 = kernel_ints(Mp2)
    S1bb, _ = rref_ints(list(S1))
    S1set = span_points(S1bb) - {0}
    assert len(Sbb) == 2 and len(S1bb) == 2 and len(seam_set) == 3 \
        and len(S1set) == 3
    dSK = len(Sbb) + len(rref_ints(list(K1))[0]) - \
        len(rref_ints(list(SEAM) + list(K1))[0])
    assert dSK == 0, "pair (0,1) not R2"
    assert span_eq(S1, K2), "pair (1,2) not R3"
    assert not deck0.twisted() and not deck1.twisted() and deck2.twisted()
    log("structure: pairs (0,1) = R2 (S cap K = 0), (1,2) = R3 "
        "(S1 = ker p2*); rungs x-U/x-U/y-T all (R)  [A35 row REPRODUCED]")
    out["structure"] = {"SEAM_dim": 2, "R2_pair01": True,
                        "R3_pair12": True}
    perms1 = translation_perms(L[1])
    perms2 = translation_perms(L[2])
    perms3 = translation_perms(L[3])
    binp = cosetbz.build_kernel()

    # --------------------------------------------------- 1. L3 censuses
    hits3 = census_pass(binp, L[3], [("S", np.zeros(L[3].n, np.uint8))],
                        W, "c37xx_L3stab")
    s3v = [i2v(h, L[3].n) for h in sorted(hits3["S"])]
    for v in s3v[:: max(1, len(s3v) // 40)]:
        assert L[3].is_stab(v)
    mu3 = min(int(v.sum()) for v in s3v)
    assert mu3 >= 6, f"mu(L3) = {mu3} < 6: cap formulas need redesign"
    c3 = Collector(L[3].n)
    for v in s3v:
        c3.add(v)
    stab3_reps = c3.reps(perms3)
    log(f"L3 stab census <= {W}: {len(s3v)} vectors "
        f"{whist(s3v)}, {len(stab3_reps)} orbit reps, mu3 = {mu3}")
    out["L3_stab"] = {"vectors": len(s3v), "weight_hist": whist(s3v),
                      "orbits": len(stab3_reps), "mu": mu3}

    cls_reps3 = [rep_for(L[3], c) for c in range(1, 16)]
    hits3c = census_pass(
        binp, L[3],
        [(f"C{c}", cls_reps3[c - 1]) for c in range(1, 16)],
        WC, "c37xx_L3coset")
    cos3: dict[int, list[np.ndarray]] = {}
    d_l3 = None
    for c in range(1, 16):
        vs = [i2v(h, L[3].n) for h in sorted(hits3c[f"C{c}"])]
        for v in vs[:: max(1, len(vs) // 10)]:
            assert L[3].is_cycle(v) and not L[3].is_stab(v)
            assert v2i(L[3].sig(v)) == c
        cos3[c] = vs
        if vs:
            mn = min(int(v.sum()) for v in vs)
            d_l3 = mn if d_l3 is None else min(d_l3, mn)
    n_cos3 = sum(len(v) for v in cos3.values())
    assert d_l3 is not None and d_l3 <= WC
    log(f"L3 coset censuses <= {WC}: {n_cos3} elements over 15 classes; "
        f"d(L3) = {d_l3} EXACT (census-complete, free by-product)")
    out["d_L3"] = d_l3

    # ------------------------------------------- 2. L2 censuses (descent)
    coll_stab2 = Collector(L[2].n)
    coll_s1cos2 = Collector(L[2].n)
    coll_all2 = Collector(L[2].n)
    tripwire2 = [0]

    def classify_l2(b2: np.ndarray) -> None:
        w = int(b2.sum())
        if w == 0:
            return
        is_st = in_span(v2i(b2), L[2].rsHX_b, L[2].rsHX_p)
        if is_st:
            if w <= W:
                coll_stab2.add(b2)
        else:
            assert w >= D_L2, \
                f"nontrivial L2-cycle of weight {w} < {D_L2}: " \
                f"CONTRADICTS the A30 certificate!"
            tripwire2[0] += 1
            if w <= W and v2i(L[2].sig(b2)) in S1set:
                coll_s1cos2.add(b2)
        if w <= WC:
            coll_all2.add(b2)

    # tau3-family
    n_tau3 = 0
    for src in ([v for v in s3v] +
                [v for c in range(1, 16) for v in cos3[c]]):
        if 2 * int(src.sum()) <= max(W, WC):
            b2 = (deck2.TAU @ src) % 2
            assert L[2].is_cycle(b2)
            classify_l2(b2)
            n_tau3 += 1
    # fibers over L3-stab orbit reps
    Wmax2 = max(W, WC)
    nf = nl = 0
    for beta in stab3_reps:
        cap = (Wmax2 - int(beta.sum())) // 2
        if cap < 0:
            continue
        lifts = enumerate_lifts_deep(deck2, beta, cap=min(cap, 8))
        nf += 1
        for v0c, m2 in sorted(lifts.items()):
            classify_l2(deck2.lift(i2v(v0c, L[3].n), beta))
            nl += 1
        if nf % 100 == 0:
            log(f"  ... L3-stab fibers {nf}/{len(stab3_reps)} "
                f"({nl} lifts)")
    # fibers over L3-coset orbit reps (the 12-class species <= WC)
    nfc = 0
    for c in range(1, 16):
        cc = Collector(L[3].n)
        for v in cos3[c]:
            cc.add(v)
        for beta in cc.reps(perms3):
            cap = (WC - int(beta.sum())) // 2
            if cap < 0:
                continue
            lifts = enumerate_lifts_deep(deck2, beta, cap=min(cap, 8))
            nfc += 1
            for v0c, m2 in sorted(lifts.items()):
                classify_l2(deck2.lift(i2v(v0c, L[3].n), beta))
    stab2_reps = coll_stab2.reps(perms2)
    s1cos2_reps = coll_s1cos2.reps(perms2)
    all2_reps = coll_all2.reps(perms2)
    mu2 = min(int(b.sum()) for b in stab2_reps)
    assert mu2 >= 6
    d2_seen = min((int(b.sum()) for b in all2_reps
                   if not in_span(v2i(b), L[2].rsHX_b, L[2].rsHX_p)),
                  default=None)
    assert d2_seen == D_L2, f"lightest nontrivial L2 = {d2_seen} != {D_L2}"
    log(f"L2 descent censuses: stab reps <= {W}: {len(stab2_reps)} "
        f"{whist(stab2_reps)} (mu2 = {mu2}); S1-coset reps <= {W}: "
        f"{len(s1cos2_reps)} {whist(s1cos2_reps)}; all-class reps "
        f"<= {WC}: {len(all2_reps)}; lightest nontrivial = {d2_seen} == "
        f"d(L2) CROSS-CHECK; {nf}+{nfc} fibers -> {nl}+ lifts; "
        f"tau3 sources {n_tau3}; tripwires {tripwire2[0]} (all >= {D_L2})")
    out["L2"] = {"stab_orbits": len(stab2_reps),
                 "stab_whist": whist(stab2_reps),
                 "s1_orbits": len(s1cos2_reps),
                 "s1_whist": whist(s1cos2_reps),
                 "allclass_orbits": len(all2_reps), "mu2": mu2}

    # ---- 2X. banked A30 cross-checks (free falsify-first bonus):
    # (i) the <= 18 slice of our descent-derived L2-stab orbit census
    # must equal the banked A30 37a70e02:x dangerous-cell census band
    # by band (A30 derived it by its own direct method);
    # (ii) the S1-coset census <= 18 must be EMPTY — that emptiness IS
    # A30's SeamCosetFloor-20 safe-floor certificate, re-derived by the
    # composed L3->L2 fiber machinery.
    a30r = json.loads((A30 / "rungs_37a70e02_x.json").read_text())
    wh2 = whist(stab2_reps)
    lane_of = {"6": "bz", "10": "restricted<= 4", "12": "restricted<= 3",
               "14": "restricted<= 2", "16": "restricted<= 1",
               "18": "restricted<= 0"}
    n18 = 0
    for wkey, lane in lane_of.items():
        got = wh2.get(wkey, 0)
        want = a30r["per_lane"][lane]["n"]
        assert got == want, \
            f"A30 census cross-check FAIL at w{wkey}: {got} != {want}"
        n18 += got
    assert a30r["n_classes_total"] == n18 == 2203
    s1_le18 = [b for b in s1cos2_reps if int(b.sum()) <= 18]
    assert not s1_le18, \
        f"S1-coset elements <= 18 exist ({len(s1_le18)}): CONTRADICTS " \
        f"the banked A30 safe-floor-20 certificate!"
    log("A30 cross-checks: <= 18 L2-stab bands == banked 37a70e02:x "
        "dangerous census (2,203 cells, band-by-band EXACT); S1-coset "
        "<= 18 EMPTY == A30's SeamCosetFloor-20 re-derived through the "
        "composed fibers")
    out["A30_crosscheck"] = {"dangerous_2203": True,
                            "safe_floor_20_rederived": True}

    # ---- 2G. composed-machinery gate: direct L2 BZ census <= 14
    WG = 14
    g_offsets = [("S", np.zeros(L[2].n, np.uint8))]
    for c in sorted(S1set):
        rv = rep_for(L[2], c)
        g_offsets.append((f"C{c}", rv))
    hits2d = census_pass(binp, L[2], g_offsets, WG, "c37xx_L2gate")
    dir_stab = [i2v(h, L[2].n) for h in sorted(hits2d["S"])]
    for v in dir_stab[:: max(1, len(dir_stab) // 30)]:
        assert L[2].is_stab(v)
    cg = Collector(L[2].n)
    for v in dir_stab:
        cg.add(v)
    kd = {bytes(k) for k in
          batch_keys(np.array([v for v in cg.vecs], dtype=np.uint8),
                     perms2)} if cg.vecs else set()
    m14 = [b for b in stab2_reps if int(b.sum()) <= WG]
    km = {bytes(k) for k in
          batch_keys(np.array(m14, dtype=np.uint8), perms2)} \
        if m14 else set()
    assert kd == km, \
        f"L2-stab gate FAIL: direct {len(kd)} orbits != descent {len(km)}"
    cg1 = Collector(L[2].n)
    n_dir_s1 = 0
    for c in sorted(S1set):
        for h in sorted(hits2d[f"C{c}"]):
            v = i2v(h, L[2].n)
            assert L[2].is_cycle(v) and not L[2].is_stab(v)
            cg1.add(v)
            n_dir_s1 += 1
    kd1 = {bytes(k) for k in
           batch_keys(np.array(cg1.vecs, dtype=np.uint8), perms2)} \
        if cg1.vecs else set()
    m14s = [b for b in s1cos2_reps if int(b.sum()) <= WG]
    km1 = {bytes(k) for k in
           batch_keys(np.array(m14s, dtype=np.uint8), perms2)} \
        if m14s else set()
    assert kd1 == km1, \
        f"L2 S1-coset gate FAIL: direct {len(kd1)} != descent {len(km1)}"
    log(f"L2 GATE: direct BZ <= {WG} == descent slice EXACTLY (stab "
        f"orbits {len(kd)}; S1-coset orbits {len(kd1)} from {n_dir_s1} "
        f"elements) — composed fiber machinery census-complete")
    out["L2_gate"] = {"W": WG, "stab_orbits": len(kd),
                      "s1_orbits": len(kd1), "equal": True}

    # --------------------------------------- 3. L1 obligations (descent)
    ks_fib1 = KernelShift(deck1, all2_reps, complete_to=WC)
    coll_stab1 = Collector(L[1].n)
    coll_seam1 = Collector(L[1].n)
    tripwire1 = [0]

    def classify_l1(b1: np.ndarray) -> None:
        w = int(b1.sum())
        if w == 0 or w > W:
            return
        if in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
            coll_stab1.add(b1)
            return
        assert w >= D_L1, \
            f"nontrivial L1-cycle of weight {w} < {D_L1}: CONTRADICTS " \
            f"the A30 [[360,4,20]] certificate!"
        tripwire1[0] += 1
        if v2i(L[1].sig(b1)) in seam_set:
            coll_seam1.add(b1)

    # tau2-family over the all-class L2 census (reps expanded in KS)
    n_tau2 = 0
    for w, bucket in ks_fib1.by_w.items():
        if 2 * w <= W:
            for zi in bucket:
                b1 = (deck1.TAU @ i2v(zi, L[2].n)) % 2
                assert L[1].is_cycle(b1)
                classify_l1(b1)
                n_tau2 += 1
    log(f"L1 tau2-family: {n_tau2} sources processed")
    # fibers over L2-stab reps: direct deep for cap <= 6, kernel-shift
    # for light b2 (window inside the <= WC all-class census)
    nfd = nfk = 0
    for b2 in sorted(stab2_reps, key=lambda b: -int(b.sum())):
        wb2 = int(b2.sum())
        cap = (W - wb2) // 2
        if cap < 0:
            continue
        if cap <= 6:
            lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
            for v0c, m2 in sorted(lifts.items()):
                classify_l1(deck1.lift(i2v(v0c, L[2].n), b2))
            nfd += 1
        else:
            v0p, ovp = row_lift_v0(deck1, b2)
            seen_v0: set[int] = set()
            rhs = (deck1.RHS @ b2) % 2
            for v0i in ks_fib1.candidates(b2, v0p, cap):
                canon = min(v0i, v0i ^ v2i(b2))
                if canon in seen_v0:
                    continue
                seen_v0.add(canon)
                v0 = i2v(v0i, L[2].n)
                assert not (((deck1.E @ v0) + rhs) % 2).any()
                classify_l1(deck1.lift(v0, b2))
            nfk += 1
        if (nfd + nfk) % 200 == 0:
            log(f"  ... L2-stab fibers {nfd+nfk}/{len(stab2_reps)}")
    stab1_reps = coll_stab1.reps(perms1)
    mu1 = min(int(b.sum()) for b in stab1_reps)
    assert mu1 >= 6
    log(f"L1 stab census <= {W} (descent): {len(stab1_reps)} orbit reps "
        f"{whist(stab1_reps)} (mu1 = {mu1}); fibers {nfd} direct + "
        f"{nfk} kernel-shift; tripwires {tripwire1[0]} (all >= {D_L1})")
    out["L1_stab"] = {"orbits": len(stab1_reps),
                      "whist": whist(stab1_reps), "mu1": mu1,
                      "fibers_direct": nfd, "fibers_kernel_shift": nfk}

    # ---- 3G. second-layer gate: the L1-stab census <= 12 re-derived
    # INDEPENDENTLY through the OTHER quotient (the y-fold
    # (30,6) -> (30,3), n = 180 <= 192): a different deck, different
    # fibers, different tau-family — orbit key sets must agree exactly.
    # (A direct L1 BZ census is impossible: n = 360 > the 3x64-bit C
    # kernel; this cross-derivation is the A36 SS4 descent-lane pattern.)
    WG1 = 12
    LY = TowerCode("LY", (30, 3), fold_support(levels[1][1], 1, 3),
                   fold_support(levels[1][2], 1, 3))
    assert (1 << LY.k) - 1 <= 256, f"k(LY) = {LY.k} too large for offsets"
    decky = AxisDeck(L[1], LY, 1)
    permsy = translation_perms(LY)
    hitsy = census_pass(binp, LY, [("S", np.zeros(LY.n, np.uint8))],
                        WG1, "c37xx_LYgate")
    syv = [i2v(h, LY.n) for h in sorted(hitsy["S"])]
    for v in syv[:: max(1, len(syv) // 20)]:
        assert LY.is_stab(v)
    cy = Collector(LY.n)
    for v in syv:
        cy.add(v)
    # y-route tau-family: tau_y(u), u any LY-cycle <= 6 (coset census)
    cly = [rep_for(LY, c) for c in range(1, 1 << LY.k)]
    hitsyc = census_pass(
        binp, LY,
        [(f"C{c}", cly[c - 1]) for c in range(1, 1 << LY.k)],
        WG1 // 2, "c37xx_LYcos")
    ysrc = [v for v in syv if 2 * int(v.sum()) <= WG1]
    for c in range(1, 1 << LY.k):
        for h in sorted(hitsyc[f"C{c}"]):
            v = i2v(h, LY.n)
            if 2 * int(v.sum()) <= WG1:
                ysrc.append(v)
    coll_y = Collector(L[1].n)
    for u in ysrc:
        b1 = (decky.TAU @ u) % 2
        assert L[1].is_cycle(b1)
        if 0 < int(b1.sum()) <= WG1 \
                and in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
            coll_y.add(b1)
    nyf = 0
    for beta in cy.reps(permsy):
        cap = (WG1 - int(beta.sum())) // 2
        if cap < 0:
            continue
        lifts = enumerate_lifts_deep(decky, beta, cap=cap)
        nyf += 1
        for v0c, m2 in sorted(lifts.items()):
            b1 = decky.lift(i2v(v0c, LY.n), beta)
            w = int(b1.sum())
            if 0 < w <= WG1 \
                    and in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
                coll_y.add(b1)
    yreps = coll_y.reps(perms1)
    kd1g = {bytes(k) for k in
            batch_keys(np.array(yreps, dtype=np.uint8), perms1)} \
        if yreps else set()
    m12 = [b for b in stab1_reps if int(b.sum()) <= WG1]
    km1g = {bytes(k) for k in
            batch_keys(np.array(m12, dtype=np.uint8), perms1)} \
        if m12 else set()
    assert kd1g == km1g, \
        f"L1-stab gate FAIL: y-route {len(kd1g)} != x-route {len(km1g)}"
    log(f"L1 GATE: the <= {WG1} L1-stab census re-derived through the "
        f"INDEPENDENT y-quotient (30,3) ({nyf} fibers) == the x-route "
        f"descent slice EXACTLY ({len(kd1g)} orbits) — second-layer "
        f"completeness cross-derived")
    out["L1_gate"] = {"W": WG1, "orbits": len(kd1g), "equal": True,
                      "route": "independent y-quotient (30,3)"}

    # seam species (Q2 only; DEAD at W = 18 by d(L1) = 20)
    seam1_reps: list[np.ndarray] = []
    if W >= D_L1:
        for b2 in s1cos2_reps:
            cap = (W - int(b2.sum())) // 2
            if cap < 0:
                continue
            assert cap <= 6, f"seam fiber cap {cap} > 6?!"
            lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
            for v0c, m2 in sorted(lifts.items()):
                classify_l1(deck1.lift(i2v(v0c, L[2].n), b2))
        seam1_reps = coll_seam1.reps(perms1)
        log(f"L1 SEAM-coset census <= {W} (descent, orbit reps): "
            f"{len(seam1_reps)} {whist(seam1_reps)} — all >= {D_L1} "
            f"(tripwired)")
        out["L1_seam"] = {"orbits": len(seam1_reps),
                          "whist": whist(seam1_reps)}
    else:
        log(f"L1 SEAM branch: DEAD outright at W = {W} — seam elements "
            f"are nontrivial L1-logicals, weight >= d(L1) = {D_L1} > "
            f"{W} (A30 certificate consumed; no census needed)")
        out["L1_seam"] = {"dead_by_certificate": True, "d_L1": D_L1}

    # ---- checkpoint the L1 obligations (the harness kills background
    # tasks at ~1 h wall; phase 2 = a38_c37xx_rungs.py reloads them,
    # RE-VERIFIES every vector, and runs the rungs + assembly)
    with (DATA / f"ckpt_W{W}_stab1.jsonl").open("w") as f:
        for b in stab1_reps:
            f.write(json.dumps({
                "w": int(b.sum()),
                "support": sorted(int(j) for j in np.nonzero(b)[0]),
            }) + "\n")
    with (DATA / f"ckpt_W{W}_seam1.jsonl").open("w") as f:
        for v in seam1_reps:
            f.write(json.dumps({
                "w": int(v.sum()),
                "support": sorted(int(j) for j in np.nonzero(v)[0]),
            }) + "\n")
    log(f"checkpoint written: {len(stab1_reps)} stab + "
        f"{len(seam1_reps)} seam orbit-rep supports")
    if census_only:
        out["wall_s"] = round(time.monotonic() - t0, 1)
        (DATA / f"freeze_W{W}_census.json").write_text(
            json.dumps(out, indent=1))
        log(f"census phase complete ({out['wall_s']}s) -> "
            f"{DATA / f'freeze_W{W}_census.json'}; run "
            f"a38_c37xx_rungs.py {W} for phase 2")
        return

    # --------------------------------------------------- 4. the rungs
    cell = RungCell("c37xx_top", L[1], L[0], deck0)
    assert len(cell.sector_basis) == 4
    stab1_win = [b for b in stab1_reps]      # stabs <= W, complete
    ks_top = KernelShift(deck0, stab1_win, complete_to=W)

    def kernel_shift_rung(b: np.ndarray, M: int) -> dict:
        """The F2b lane: candidates = row-lift (+) Z(L1)-census window;
        sound when the window bound stays below d(L1) = 20 (every shift
        is then a censused stabilizer)."""
        wb = int(b.sum())
        v0p, ovp = row_lift_v0(deck0, b)
        B = wb + (M - 1) + ovp
        assert B < D_L1, \
            f"kernel-shift window {B} >= d(L1): cell |b|={wb} M={M} " \
            f"needs the logical censuses — escalate"
        rhs = (deck0.RHS @ b) % 2
        rhs_i = v2i(rhs)
        viols = []
        seen_v0: set[int] = set()
        bmask = v2i(b)
        for v0i in ks_top.candidates(b, v0p, M - 1):
            canon = min(v0i, v0i ^ bmask)
            if canon in seen_v0:
                continue
            seen_v0.add(canon)
            v0 = i2v(v0i, L[1].n)
            assert not (((deck0.E @ v0) + rhs) % 2).any()
            ch = (deck0.EMB[0] @ v0
                  + deck0.EMB[1] @ ((v0 + b) % 2)) % 2
            if in_span(v2i(ch), L[0].rsHX_b, L[0].rsHX_p):
                continue
            ov = bin(v0i & ~bmask).count("1")
            wt = int(ch.sum())
            assert wt == wb + 2 * ov, "slice identity violated"
            viols.append({"overflow": ov, "weight": wt,
                          "v0_hex": f"{v0i:x}"})
        if viols:
            return {"verdict": "VIOLATION", "lane": "kernel-shift",
                    "w_b": wb, "M": M, "violations": viols[:5],
                    "n_viol": len(viols),
                    "min_weight": min(x["weight"] for x in viols)}
        return {"verdict": "PASS", "lane": "kernel-shift", "w_b": wb,
                "M": M}

    rung_rows = []
    verd: dict[str, int] = {}
    lanes: dict[str, int] = {}
    viol_finds: list[dict] = []
    first3: list[tuple[np.ndarray, int, str]] = []
    n_xval = xval_cheap = xval_deep = 0
    n_done = 0
    tR = time.monotonic()
    for b in sorted(stab1_reps, key=lambda b: -int(b.sum())):
        wb = int(b.sum())
        M = (TARGET - wb) // 2
        if M <= 0:
            continue
        # cost-based lane assignment: direct restricted for small caps
        # (O(offb^2) at most), kernel-shift for deep caps when its
        # window closes below d(L1), direct <= 6 as the fallback.
        if (M - 1) <= 4:
            r = cell.rung(b, M, time.monotonic() + 3600)
            if xval_cheap < 3:            # cross-validate lanes
                v0p, ovp = row_lift_v0(deck0, b)
                if wb + (M - 1) + ovp < D_L1:
                    rks = kernel_shift_rung(b, M)
                    assert rks["verdict"] == r["verdict"], \
                        f"LANE MISMATCH |b|={wb} M={M}"
                    n_xval += 1
                    xval_cheap += 1
        else:
            v0p, ovp = row_lift_v0(deck0, b)
            B = wb + (M - 1) + ovp
            if B < D_L1 and B <= W:
                r = kernel_shift_rung(b, M)
                if xval_deep < 3 and (M - 1) <= 6:
                    rd = cell.rung(b, M, time.monotonic() + 3600)
                    assert rd["verdict"] == r["verdict"], \
                        f"LANE MISMATCH |b|={wb} M={M}: ks " \
                        f"{r['verdict']} vs direct {rd['verdict']}"
                    if r["verdict"] == "VIOLATION":
                        assert rd["n_viol"] == r["n_viol"], (rd, r)
                    n_xval += 1
                    xval_deep += 1
            else:
                assert (M - 1) <= 6, \
                    f"cell |b|={wb} M={M}: no sound lane (B={B}, " \
                    f"direct cap {M-1} > 6) — ESCALATE"
                r = cell.rung(b, M, time.monotonic() + 3600)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
        rung_rows.append({"species": "dangerous", "w": wb, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
        if r["verdict"] == "VIOLATION":
            viol_finds.append(r)
        if len(first3) < 3:
            first3.append((b, M, r["verdict"]))
        n_done += 1
        if n_done % 20000 == 0:
            log(f"  ... dangerous rungs {n_done}/{len(stab1_reps)} "
                f"({time.monotonic()-tR:.0f}s)")
    dtR = time.monotonic() - tR
    log(f"dangerous rungs: {sum(verd.values())} at target {TARGET} "
        f"({dtR:.1f}s): verdicts {verd}, lanes {lanes}; lane "
        f"cross-validation {n_xval} cells ({xval_deep} at M-1 in 3..6), "
        f"all equal")
    out["dangerous_rungs"] = {"rungs": sum(verd.values()),
                              "verdicts": verd, "lanes": lanes,
                              "xval_cells": n_xval,
                              "wall_s": round(dtR, 1)}

    verd2: dict[str, int] = {}
    if W >= D_L1 and seam1_reps:
        for w_el in sorted(seam1_reps, key=lambda v: -int(v.sum())):
            ww = int(w_el.sum())
            M = (TARGET - ww) // 2
            if M <= 0:
                continue
            r = cell.seam_rung(w_el, M)
            verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
            rung_rows.append({"species": "seam", "w": ww, "M": M,
                              "verdict": r["verdict"],
                              "lane": r["lane"]})
            if r["verdict"] == "VIOLATION":
                viol_finds.append(r)
        log(f"seam rungs: {sum(verd2.values())} at target {TARGET}: "
            f"{verd2}")
        out["seam_rungs"] = {"rungs": sum(verd2.values()),
                             "verdicts": verd2}
    if len(rung_rows) <= 50000:
        with (DATA / f"rungs_W{W}.jsonl").open("w") as f:
            for rr in rung_rows:
                f.write(json.dumps(rr) + "\n")
    else:                      # aggregate: rows carry no info beyond
        agg: dict[tuple, int] = {}     # (species, w, M, verdict, lane)
        for rr in rung_rows:
            key = (rr["species"], rr["w"], rr["M"], rr["verdict"],
                   rr["lane"])
            agg[key] = agg.get(key, 0) + 1
        with (DATA / f"rungs_W{W}.jsonl").open("w") as f:
            for (sp, w_, M_, vd, ln), cnt in sorted(agg.items()):
                f.write(json.dumps({"species": sp, "w": w_, "M": M_,
                                    "verdict": vd, "lane": ln,
                                    "count": cnt}) + "\n")

    # covariance spot-checks
    g = (7, 1)
    perm_g = perm_for(L[1], g)
    for b, M, v0_verd in first3:
        bt = b[perm_g]
        v0p, ovp = row_lift_v0(deck0, bt)
        B = int(bt.sum()) + (M - 1) + ovp
        rt = (kernel_shift_rung(bt, M) if B < D_L1 and B <= W
              else cell.rung(bt, M, time.monotonic() + 600))
        assert rt["verdict"] == v0_verd, "covariance broken"
    log("covariance: 3 translated reps re-rung, verdicts transport")

    # ------------------------------------------------- 5. the assembly
    n_viol_cells = len(viol_finds)
    if n_viol_cells == 0:
        log(f"ASSEMBLY (W = {W}): branch b = 0 dead (|v| = 2|u| >= "
            f"2 d(L1) = {2 * D_L1} > {W}); dangerous branch closed "
            f"({sum(verd.values())} rungs PASS + G-transport); seam "
            f"branch "
            + ("closed (" + str(sum(verd2.values())) + " rungs PASS + "
               "G-transport)" if W >= D_L1 else
               f"dead (weights >= {D_L1} > {W} by the A30 certificate)")
            + f" => NO nontrivial X-logical of weight <= {W}: "
            f"d([[720,4]]) >= {TARGET} at certificate tier (consuming "
            f"d(L1) = 20, d(L2) = 10 at A30 certificate tier). Z side "
            f"by BB transpose duality.")
        out["verdict"] = {"floor": TARGET, "all_pass": True}
    else:
        wts = []
        for r in viol_finds:
            wts.extend(x["weight"] for x in r["violations"])
            if "min_weight" in r:
                wts.append(r["min_weight"])
        wmin = min(wts)
        log(f"ASSEMBLY (W = {W}): VIOLATIONS — explicit nontrivial "
            f"L0-logicals found, min weight {wmin} ({n_viol_cells} "
            f"cells; every candidate re-verified in-line: E-system, "
            f"non-stab, slice identity). d([[720,4]]) <= {wmin}.")
        out["verdict"] = {"upper": wmin, "all_pass": False,
                          "viol_cells": n_viol_cells}

    # X<->Z transpose duality spot-check at L0
    ng0 = L[0].ng
    iota = np.zeros(L[0].n, dtype=np.int64)
    for i, e in enumerate(L[0].G):
        j = L[0].G.index(L[0].G.neg(e))
        iota[i] = j
        iota[ng0 + i] = ng0 + j
    swap = np.concatenate([np.arange(ng0, 2 * ng0), np.arange(0, ng0)])
    hzb, hzp = rref_ints([v2i(r) for r in L[0].HZ])
    for kv in L[0].kerHZ[:10]:
        d = kv[iota][swap]
        assert not ((L[0].HX @ d) % 2).any()
    for row in L[0].HX[::120]:
        d = row[iota][swap]
        assert in_span(v2i(d), hzb, hzp)
    log("X<->Z duality spot-check OK => Z side follows")

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / f"freeze_W{W}.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA / f'freeze_W{W}.json'}")


if __name__ == "__main__":
    main()
