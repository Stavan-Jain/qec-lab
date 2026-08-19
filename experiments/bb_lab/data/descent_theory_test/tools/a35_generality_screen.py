"""A35: how far does the tower slice calculus (A32/A33) generalize?

The precondition screen: for each candidate tower (a BB code + a chain of
free Z2 folds), measure every structural hypothesis the calculus consumes
and project every cost gate, WITHOUT running the closures themselves.

Measured per rung (cover --Z2--> base):
  - Lemma 1 transports (the TS.Deck constructor's hard asserts: chain map,
    stabilizer transport, twist-invariance im S <= ker E, sections)
  - twist status (does any lifted polynomial term carry the deck element)
  - (R) via k-preservation (A12: (R) <=> k(cover) = k(base) <=> Bezout)
  - H1 rank lattice: rank p*, rank tau*, exactness im tau* = ker p* and
    ker tau* = im p*, deck-triviality sigma* = id
  - NEW invariant: the liftable-cycle codimension  codim_lift :=
    k(base) - dim { [b] in H1(base) : E v0 = RHS b is consistent }.
    (All stabilizers lift -- rowspace transport -- so the obstruction
    lives on H1; codim_lift = rank of the carry obstruction on classes.
    Fibers over a non-liftable class are empty at EVERY overflow.)
  - sampled restricted-fiber infeasibility rate at cap <= 4 (the
    production enumerate_lifts lane on row/row-sum stabilizer shadows) --
    the "why the tower wins" rate (A32 measured 93-97%, A33 90.6%)

Measured per adjacent rung pair (two-level lattice, H1(mid)):
  S = im p_top*, K = ker p_bot*: dims, dim(S cap K), K <= S?
  (A32 regime: K <= S -- reachability decided one rung below);
  S cap K = 0? (A33/(R) regime: descent trisection collapses to one
  branch); dim p_bot*(S) (the W/W2 space) and the preimage identity
  p_bot*^{-1}(p_bot*(S)) = S (the A24 SS2.6 shape).

Projected per tower (cost gates, exact binomial arithmetic):
  kappa_l = n_l/2 - k_l/2 per level; census node counts at the bottom
  level for W_eff = d_top - 2 (parity) via the two-window r-pair rule
  (complete to r1+r2+1, balanced split); the same formula at the TOP
  level = the enumeration baseline the tower replaces; win factor;
  max fiber cap = (W_eff - mu)/2 (mu = lightest sampled stabilizer,
  <= |A|+|B|).

Falsify-first: the screen must first REPRODUCE the banked A32/A33
structural numbers (hard asserts, section [0]) and the census node
anchors within tolerance before any new tower is screened.

Anti-instances (structural, no compute): odd-|G| codes have no Z2 deck
at all; non-abelian mitten codes admit decks only at central involutions
(A26); |A| or |B| even breaks Lemma 2 (demo section [2]).

Output: data/a35/screen.json + the printed table.
Run:    cd experiments/bb_lab && uv run python scripts/a35_generality_screen.py
"""

from __future__ import annotations

import gzip
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

# --- PATH ADAPTATION (descent_theory_test extraction, 2026-08-17) ---------
# This file is `git show 2f063b0:experiments/bb_lab/scripts/a35_generality_screen.py`
# extracted verbatim into data/descent_theory_test/tools/.  The ONLY edits are
# this path block and the DATA output dir below: the original computed LAB
# from its own location (scripts/), which is wrong from tools/.  Dependency
# modules (a30_rung_pass, a32_tower_slice, a32_subclosures, a33_tower_cells,
# src/bb_lab) are byte-identical between 2f063b0 and this worktree
# (verified: `git diff --stat 2f063b0 HEAD -- <deps>` is empty).
LAB = Path(__file__).resolve().parent.parent.parent.parent  # experiments/bb_lab
assert (LAB / "src" / "bb_lab").is_dir(), f"LAB path adaptation broken: {LAB}"
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402

from a30_rung_pass import rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_subclosures import enumerate_lifts  # noqa: E402
from a33_tower_cells import h1_map, translation_mat  # noqa: E402

DATA = LAB / "data" / "descent_theory_test" / "validation"  # ADAPTED (was data/a35)
DATA.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260811)

T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:6.1f}s] {msg}", flush=True)


# ------------------------------------------------------- polynomial descent
def poly_terms(s: str, lm: tuple[int, int]) -> frozenset[tuple[int, int]]:
    G = AbelianGroup(lm)
    return frozenset((int(a), int(b)) for a, b in
                     Poly.from_string(s, G).support)


def fold_terms(terms, axis: int, newmod: int) -> frozenset[tuple[int, int]]:
    """Descend a polynomial support along the fold (mod-2 term merging)."""
    counts: dict[tuple[int, int], int] = {}
    for t in terms:
        u = list(t)
        u[axis] %= newmod
        u = tuple(u)
        counts[u] = counts.get(u, 0) + 1
    return frozenset(t for t, c in counts.items() if c % 2 == 1)


def terms_str(terms) -> str:
    def mono(t):
        a, b = t
        fs = []
        if a:
            fs.append(f"x^{a}" if a > 1 else "x")
        if b:
            fs.append(f"y^{b}" if b > 1 else "y")
        return "*".join(fs) if fs else "1"
    return " + ".join(mono(t) for t in sorted(terms)) if terms else "0"


def is_twisted(terms, axis: int, newmod: int) -> bool:
    """Does any cover term carry the deck element (exponent >= newmod)?"""
    return any(t[axis] >= newmod for t in terms)


# ----------------------------------------------------------- span utilities
def span_points(basis: list[int]) -> set[int]:
    pts = {0}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return pts


def span_dim(vs: list[int]) -> int:
    b, _ = rref_ints([v for v in vs if v])
    return len(b)


def intersect_dim(a: list[int], b: list[int]) -> int:
    """dim(span a  cap  span b) -- via dim A + dim B - dim(A+B)."""
    da, db = span_dim(a), span_dim(b)
    return da + db - span_dim(list(a) + list(b))


# ------------------------------------------------------------- cost formula
def census_nodes(kappa: int, W: int, n_bases: int = 1) -> int:
    """Exact node count of the two-window coset-BZ census, balanced r-pair
    (r1 + r2 + 1 = W); per window sum_{w=1..r} C(kappa, w), per coset base.
    Calibrated against the A32/A33 anchors in section [0] below."""
    if W < 1:
        return 0
    r1 = (W - 1 + 1) // 2
    r2 = (W - 1) - r1
    tot = 0
    for r in (r1, r2):
        tot += sum(math.comb(kappa, w) for w in range(1, r + 1))
    return tot * n_bases


# ------------------------------------------------------------- rung screen
def make_deck(cover: TS.BBCode, base: TS.BBCode, axis: int,
              newmod: int) -> TS.Deck:
    if axis == 0:
        return TS.Deck(cover, base, lambda e: (e[0] % newmod, e[1]),
                       lambda e, s: (e[0] + newmod * s, e[1]))
    return TS.Deck(cover, base, lambda e: (e[0], e[1] % newmod),
                   lambda e, s: (e[0], e[1] + newmod * s))


def liftable_codim(deck: TS.Deck) -> int:
    """rank of the carry obstruction on H1(base): k_base - dim of the
    classes [b] whose carry system E v0 = RHS b is consistent."""
    # functionals vanishing on im E = left null space of E
    F = np.array(TS.kernel_basis(deck.E.T), dtype=np.uint8)
    if F.size == 0:
        return 0
    O = (F @ deck.RHS) % 2                      # obstruction functionals
    # rank of O on the X-logical reps of the base (O vanishes on B1:
    # stabilizer rows lift by transport -- asserted cheaply here)
    for ridx in (0, deck.base.ng // 2):
        assert not ((O @ deck.base.HX[ridx]) % 2).any(), \
            "obstruction does not vanish on a stabilizer row"
    M = np.array([(O @ r) % 2 for r in deck.base.xreps], dtype=np.uint8)
    return TS.gf2_rank([v2i(r) for r in M])


def fiber_sample(deck: TS.Deck, W: int, want: int = 30,
                 tries: int = 600) -> dict:
    """Restricted-fiber probe at the production regime: HEAVY stabilizer
    shadows (|beta| in [W-6, W], caps 0-3), sampled as localized row sums
    (nearby translates overlap, keeping the weight near the band).  Light
    shadows are useless here -- single rows lift flat (m2 = 0), so any
    probe over them reads 0% empty; the production emptiness (A32: 93-97%
    of sector fibers, A33: 90.6%) lives at the heavy bands where the cap
    is tight.  Also samples mu (lightest nonzero stabilizer seen)."""
    base = deck.base
    lx, ly = base.G.orders
    tried = empty = 0
    by_cap: dict[int, list[int]] = {}
    mu = int(base.HX[0].sum())
    for _ in range(tries):
        if tried >= want:
            break
        j = int(RNG.integers(2, 6))
        g0 = (int(RNG.integers(lx)), int(RNG.integers(ly)))
        beta = np.zeros(base.n, dtype=np.uint8)
        for _ in range(j):
            da, db = int(RNG.integers(-2, 3)), int(RNG.integers(-2, 3))
            gi = base.G.index(((g0[0] + da) % lx, (g0[1] + db) % ly))
            beta = (beta + base.HX[gi]) % 2
        wb = int(beta.sum())
        if wb:
            mu = min(mu, wb)
        if not (W - 6 <= wb <= W):
            continue
        cap = (W - wb) // 2
        lifts = enumerate_lifts(deck, beta, cap)
        tried += 1
        if not lifts:
            empty += 1
        by_cap.setdefault(cap, [0, 0])
        by_cap[cap][0] += 0 if lifts else 1
        by_cap[cap][1] += 1
    return {"tried": tried, "empty": empty,
            "empty_rate": round(empty / tried, 3) if tried else None,
            "by_cap": {str(c): f"{e}/{t}" for c, (e, t)
                       in sorted(by_cap.items())},
            "mu_sampled": mu}


def screen_rung(cover: TS.BBCode, base: TS.BBCode, axis: int, newmod: int,
                W_eff: int, do_fibers: bool = True) -> dict:
    deck = make_deck(cover, base, axis, newmod)   # Lemma 1 asserts inside
    r: dict = {"axis": "xy"[axis], "fold_to": newmod,
               "lemma1": "PASS (constructor asserts)"}
    r["twisted"] = bool(is_twisted(cover.A.support, axis, newmod)
                        or is_twisted(cover.B.support, axis, newmod))
    r["k_cover"], r["k_base"] = cover.k, base.k
    r["R_holds"] = bool(cover.k == base.k)        # A12: (R) <=> k preserved
    if cover.k > 0 and base.k > 0:
        Mp = h1_map(deck)
        Mt = h1_map(deck, tau=True)
        r["rank_p"] = TS.gf2_rank([v2i(c) for c in Mp.T])
        r["rank_tau"] = TS.gf2_rank([v2i(c) for c in Mt.T])
        imP = TS._colspace(Mp)
        imT = TS._colspace(Mt)
        kerP = TS._kernel_ints(Mp)
        kerT = TS._kernel_ints(Mt)
        r["exact_cover"] = bool(TS._span_eq(imT, kerP))   # im tau* = ker p*
        r["exact_base"] = bool(TS._span_eq(kerT, imP))    # ker tau* = im p*
        sig_t = (newmod, 0) if axis == 0 else (0, newmod)
        St = translation_mat(cover, sig_t)
        r["sigma_id"] = bool((St == np.eye(cover.k, dtype=np.uint8)).all())
        r["_imP"] = imP
        r["_kerP"] = kerP
    else:
        r["rank_p"] = r["rank_tau"] = 0
        r["exact_cover"] = r["exact_base"] = None
        r["sigma_id"] = None
        r["_imP"] = []
        r["_kerP"] = []
    r["codim_lift"] = liftable_codim(deck)
    if do_fibers:
        r["fibers"] = fiber_sample(deck, W_eff)
    return r


# ------------------------------------------------------------ tower screen
def screen_tower(spec: dict) -> dict:
    name = spec["name"]
    lm, As, Bs = spec["top"]
    d_top = spec.get("d_top")            # None = unknown
    W_eff = spec.get("W_eff") or ((d_top - 2) if d_top else None)
    log(f"--- {name}: top ({lm[0]},{lm[1]}) A={As!r} B={Bs!r}")

    # build the level chain by literal descent
    levels = [(lm, poly_terms(As, lm), poly_terms(Bs, lm))]
    for axis, newmod in spec["folds"]:
        plm, pA, pB = levels[-1]
        assert plm[axis] == 2 * newmod, \
            f"{name}: fold {axis}/{newmod} is not an index-2 fold of {plm}"
        nlm = (newmod, plm[1]) if axis == 0 else (plm[0], newmod)
        levels.append((nlm, fold_terms(pA, axis, newmod),
                       fold_terms(pB, axis, newmod)))

    codes: list[TS.BBCode] = []
    lvl_rows = []
    for i, (glm, tA, tB) in enumerate(levels):
        code = TS.BBCode(f"{name}/L{i}", glm, terms_str(tA), terms_str(tB))
        codes.append(code)
        odd_ok = all(int(kv.sum()) % 2 == 0 for kv in code.kerHZ)
        kappa = len(code.rsHX_b)
        assert kappa == code.ng - code.k // 2, "kappa formula fails"
        lvl_rows.append({"lm": list(glm), "n": code.n, "k": code.k,
                         "wA": len(tA), "wB": len(tB),
                         "parity_scope": bool(len(tA) % 2 and len(tB) % 2),
                         "cycles_all_even": bool(odd_ok), "kappa": kappa})
        log(f"    L{i} ({glm[0]},{glm[1]}): [[{code.n},{code.k}]] "
            f"|A|={len(tA)} |B|={len(tB)} kappa={kappa} "
            f"even-cycles={odd_ok}")

    W_use = W_eff if W_eff else (2 * (spec.get("d_mid") or 0) - 2 or 16)
    rungs = []
    for i, (axis, newmod) in enumerate(spec["folds"]):
        r = screen_rung(codes[i], codes[i + 1], axis, newmod, W_use,
                        do_fibers=spec.get("fibers", True))
        rungs.append(r)
        log(f"    rung {i} ({r['axis']}->{newmod}): twisted={r['twisted']} "
            f"R={r['R_holds']} (k {r['k_cover']}->{r['k_base']}) "
            f"rank p*={r['rank_p']} tau*={r['rank_tau']} "
            f"exact(cov/base)={r['exact_cover']}/{r['exact_base']} "
            f"sigma*=id:{r['sigma_id']} codim_lift={r['codim_lift']}"
            + (f" fibers: {r['fibers']['empty']}/{r['fibers']['tried']} "
               f"empty (rate {r['fibers']['empty_rate']})"
               if "fibers" in r else ""))

    # adjacent-pair two-level lattice (in H1 of the middle code)
    pairs = []
    for i in range(len(rungs) - 1):
        S = rungs[i]["_imP"]           # im p_top*  in H1(level i+1)
        K = rungs[i + 1]["_kerP"]      # ker p_bot* in H1(level i+1)
        mid, bot = codes[i + 1], codes[i + 2]
        if mid.k == 0 or not S:
            pairs.append({"i": i, "note": "H1(mid) trivial"})
            continue
        dS, dK = span_dim(S), span_dim(K)
        dSK = intersect_dim(S, K)
        Mb = h1_map(make_deck(mid, bot, *spec["folds"][i + 1])) \
            if bot.k > 0 else None
        if Mb is not None:
            Wimg = [TS._apply(Mb, s) for s in span_points(
                rref_ints(list(S))[0])]
            dW = span_dim(Wimg)
            Wb2, Wp2 = rref_ints([w for w in Wimg if w])
            pre = TS._preimage(Mb, Wb2, Wp2)
            reach = bool(TS._span_eq(pre, rref_ints(list(S))[0])) \
                if dSK == dK else False
        else:
            dW, reach = 0, None
        Kb, Kp = rref_ints(list(K)) if K else ([], [])
        K_in_S = bool(all(TS.in_span(x, *rref_ints(list(S))) for x in K)) \
            if K else True
        pairs.append({"i": i, "dim_S": dS, "dim_K": dK, "dim_SK": dSK,
                      "K_in_S": K_in_S, "one_branch": bool(dSK == 0),
                      "dim_W": dW, "reach_preimage_eq_S": reach})
        log(f"    pair ({i},{i+1}): dim S={dS} K={dK} S^K={dSK} "
            f"K<=S:{K_in_S} one-branch:{dSK == 0} dim W={dW} "
            f"preimage=S:{reach}")

    # cost gates: one row per certification question (W value)
    W_list = spec.get("W_list") or ([W_eff] if W_eff else [])
    mu = min((r.get("fibers", {}).get("mu_sampled") or 99)
             for r in rungs) if rungs else 6
    mu = min(mu, lvl_rows[-1]["wA"] + lvl_rows[-1]["wB"])
    costs = []
    for W in W_list:
        lvl_nodes = [census_nodes(lv["kappa"], W) for lv in lvl_rows]
        nb, nt = lvl_nodes[-1], lvl_nodes[0]
        cap_max = (W - mu) // 2
        verdict = ("GREEN" if nb <= 2e11 and cap_max <= 8 else
                   "AMBER" if nb <= 1e14 and cap_max <= 12 else "RED")
        costs.append({
            "W": W, "mu": mu, "cap_max": cap_max, "verdict": verdict,
            "log10_nodes_per_level":
                [round(math.log10(x), 1) if x else None for x in lvl_nodes],
            "win_factor_vs_top": round(nt / nb, 1) if nb else None,
        })
        log(f"    costs @W={W}: per-level log10 nodes = "
            f"{costs[-1]['log10_nodes_per_level']} "
            f"win={costs[-1]['win_factor_vs_top']:.1e}x "
            f"cap_max={cap_max} -> {verdict}")

    for r in rungs:                       # strip non-JSON internals
        r.pop("_imP", None)
        r.pop("_kerP", None)
    return {"name": name, "tag": spec.get("tag", ""),
            "d_top": d_top, "levels": lvl_rows, "rungs": rungs,
            "pairs": pairs, "costs": costs, "notes": spec.get("notes", "")}


# ------------------------------------------------------------------ docket
TOWERS = [
    dict(name="bravyi360", tag="VAL-A32",
         top=((30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26"),
         folds=[(1, 3), (0, 15)], d_top=24, d_mid=12,
         notes="A32 instance: mixed-axis, top rung (R)-fails, d=24 CLOSED"),
    dict(name="ibm288Y", tag="VAL-A33",
         top=((18, 8), "1 + x*y^4 + x^14*y", "1 + x*y^2 + x^2*y^7"),
         folds=[(1, 4), (1, 2)], d_top=20, d_mid=10,
         notes="A33 instance: same-axis, all-(R), d=20 CLOSED"),
    dict(name="gross_xx", tag="NEW",
         top=((12, 6), "x^3 + y + y^2", "y^3 + x + x^2"),
         folds=[(0, 6), (0, 3)], d_top=12, d_mid=6,
         notes="retrospective target: gross over bb72 over (3,6)"),
    dict(name="bb288_yxx", tag="NEW",
         top=((12, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2"),
         folds=[(1, 6), (0, 6), (0, 3)], d_top=18, d_mid=12,
         notes="published record [[288,12,18]]; y-quotient IS gross "
               "(y^7 = y mod 6) -> 4-level tower to n=36"),
    dict(name="c37x_360420", tag="NEW",
         top=((30, 6), "1 + y + x", "y^4 + x + x^11*y^2"),
         folds=[(0, 15), (1, 3)], d_top=20, d_mid=10,
         notes="A30 doubled code 37a70e02:x = [[360,4,20]]; mixed-axis"),
    dict(name="e5e50yy_360420", tag="NEW",
         top=((15, 12), "1 + y + x", "y^4 + x^8*y^2 + x^13"),
         folds=[(1, 6), (1, 3)], d_top=20, d_mid=10,
         notes="A30 doubled code 5e50a976:y = [[360,4,20]]; same-axis"),
    dict(name="c37xx_720", tag="FRONTIER",
         top=((60, 6), "1 + y + x", "y^4 + x + x^11*y^2"),
         folds=[(0, 30), (0, 15), (1, 3)], d_top=None, d_mid=20,
         W_eff=18, W_list=[18, 22, 30, 38],
         notes="rung-2 re-double [[720,4,?]]; freeze-vs-double open "
               "(A14 SS13 / A33 SS6 ranking #1); W=18 certifies d=20 if "
               "frozen, W=38 = full doubling budget"),
    dict(name="a8_336", tag="NEW",
         top=((12, 14), "1 + y + x^3*y^3", "1 + x + x^2*y^7"),
         folds=[(0, 6), (1, 7)], d_top=12, d_mid=6,
         notes="A8/A29 [[336,12,12]] over its [[168,12,6]] base, then y"),
    dict(name="bravyi756", tag="FRONTIER",
         top=((21, 18), "x^3 + y^10 + y^17", "y^5 + x^3 + x^19"),
         folds=[(1, 9)], d_top=34, fibers=False,
         notes="[[756,16,<=34]]: ONE deck only (v2(18)=1, 21 odd); "
               "bottom n=378"),
    dict(name="cover300", tag="DEGENERATE-1LVL",
         top=((5, 30), "1 + y + x", "x*y^6 + x*y^10 + x^2*y^12"),
         folds=[(1, 15)], d_top=16, d_mid=8,
         notes="A15's [[300,8,16]]: v2(30)=1 -> one-level = the A15/A30 "
               "architecture (already closed there)"),
    dict(name="pair72", tag="NEW-TINY",
         top=((6, 6), "x^2 + y + y^3", "1 + x + y^2"),
         folds=[(0, 3), (1, 3)], d_top=8, d_mid=4,
         notes="[[72,4,8]] over [[36,4,4]] over (3,3); smallest 2-level"),
]

NO_DECK = [  # |G| odd: no free Z2 deck exists at all (layer L1 fails)
    ("bb90", (15, 3), "[[90,8,10]] -- A19 survey: 'not a Z2-cover at all'"),
    ("bb98", (7, 7), "[[98,6,12]] (A16 host)"),
    ("f2a6_base", (5, 15), "[[150,8,8]] -- the A24 odd-|G| lane"),
]


def main() -> None:
    out: dict = {"towers": [], "no_deck": [], "anchors": {}}

    # ------------------------------------------------- [0] cost anchors
    log("[0] census node formula vs banked anchors")
    a32_w = census_nodes(41, 22, 3)      # A32 W-coset <=22: reported 1.9e10
    a32_s = census_nodes(41, 22, 1)      # A32 stab <=22:    reported 6.3e9
    a33_h = census_nodes(68, 18, 3)      # A33 H5 pass:      reported 6.62e10
    for tag, got, want in [("a32_wcoset22", a32_w, 1.9e10),
                           ("a32_stab22", a32_s, 6.3e9),
                           ("a33_h5", a33_h, 6.62e10)]:
        ratio = got / want
        log(f"    {tag}: formula {got:.2e} vs banked {want:.2e} "
            f"(x{ratio:.2f})")
        assert 0.3 < ratio < 3.5, f"node formula off at {tag}"
        out["anchors"][tag] = {"formula": float(got), "banked": want,
                               "ratio": round(ratio, 2)}

    # ------------------------------------------------- [0b] fiber anchor
    # The shallow row-sum probe below CANNOT see the production emptiness
    # (its shadows are rowspace-shallow and mostly lift).  Validate the
    # enumerator itself bit-level against the banked A32 sector-C fiber
    # layer, and read the true census-population rates from the bank.
    stab_f = LAB / "data" / "a32" / "gb_census_stab.jsonl"
    if stab_f.exists():
        log("[0b] fiber enumerator vs banked A32 sector-C fibers")
        GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
        BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
        deck_x = make_deck(BY, GB, 0, 15)
        reps: dict[int, np.ndarray] = {}
        for line in stab_f.open():
            r = json.loads(line)
            if r["w"] in (14, 16) and r["canon"] not in reps:
                v = np.zeros(GB.n, dtype=np.uint8)
                v[r["support"]] = 1
                reps[r["canon"]] = v
        assert len(reps) == 397, f"{len(reps)} orbit reps != 64+333"
        tot = empty = 0
        hist: dict[int, int] = {}
        for beta in reps.values():
            cap = (22 - int(beta.sum())) // 2
            lifts = enumerate_lifts(deck_x, beta, cap)
            tot += len(lifts)
            empty += 0 if lifts else 1
            for m2 in lifts.values():
                hist[m2] = hist.get(m2, 0) + 1
        assert tot == 4132, f"sector-C lift total {tot} != banked 4,132"
        assert hist == {0: 416, 1: 554, 2: 561, 3: 1639, 4: 962}, hist
        log(f"    397 orbit fibers (caps 4/3): {tot} lifts == banked "
            f"4,132 EXACTLY, m2-hist matches; {empty} empty fibers")
        out["anchors"]["sectorC_1416_refiber"] = {
            "fibers": 397, "lifts": tot, "empty": empty,
            "m2_hist": {str(k): v for k, v in sorted(hist.items())}}
        band: dict[int, tuple[int, int]] = {}
        for line in gzip.open(LAB / "data" / "a32"
                              / "sectorAC_C18to22.jsonl.gz"):
            r = json.loads(line)
            e, t = band.get(r["wbeta"], (0, 0))
            band[r["wbeta"]] = (e + (0 if r.get("m2_hist") else 1), t + 1)
        out["anchors"]["sectorC_heavy_empty_rates_banked"] = {
            str(w): {"empty": e, "fibers": t, "rate": round(e / t, 3)}
            for w, (e, t) in sorted(band.items())}
        log("    banked heavy-band empty rates (caps 2/1/0): " + ", ".join(
            f"|b|={w}: {e}/{t} ({e/t:.0%})"
            for w, (e, t) in sorted(band.items())))

    # ------------------------------------------------- [1] tower screens
    for spec in TOWERS:
        res = screen_tower(spec)
        out["towers"].append(res)

        if spec["name"] == "bravyi360":   # falsify-first: banked structure
            r0, r1 = res["rungs"]
            assert [lv["k"] for lv in res["levels"]] == [12, 8, 8]
            assert (not r0["R_holds"]) and r1["R_holds"]
            assert r0["rank_p"] == 6 and r1["rank_p"] == 4 \
                and r1["rank_tau"] == 4
            assert r1["exact_cover"] and r1["exact_base"]
            assert r0["sigma_id"] is False and r0["twisted"] \
                and r1["twisted"]
            p = res["pairs"][0]
            assert p["dim_S"] == 6 and p["dim_K"] == 4 \
                and p["dim_SK"] == 4 and p["K_in_S"] \
                and p["dim_W"] == 2 and p["reach_preimage_eq_S"]
            log("    bravyi360: ALL banked A32 structure REPRODUCED")
        if spec["name"] == "ibm288Y":
            r0, r1 = res["rungs"]
            assert [lv["k"] for lv in res["levels"]] == [8, 8, 8]
            assert r0["R_holds"] and r1["R_holds"]
            assert r0["rank_p"] == r0["rank_tau"] == 4 == r1["rank_p"] \
                == r1["rank_tau"]
            assert r0["exact_cover"] and r0["exact_base"] \
                and r1["exact_cover"] and r1["exact_base"]
            assert r0["sigma_id"] and r1["sigma_id"]
            assert r0["twisted"] and r1["twisted"]
            p = res["pairs"][0]
            assert p["dim_S"] == 4 and p["dim_SK"] == 0 \
                and p["one_branch"] and p["dim_W"] == 4
            log("    ibm288Y: ALL banked A33 structure REPRODUCED")
        if spec["name"] == "gross_xx":    # A19 deck-survey k verdicts
            assert [lv["k"] for lv in res["levels"]] == [12, 12, 8], \
                "A19: gross x-deck R-holds (bb72), bb72 decks jump 8->12"
            log("    gross_xx: A19 deck-survey k pattern REPRODUCED")
        if spec["name"] == "bb288_yxx":
            assert res["levels"][1]["lm"] == [12, 6]
            assert res["rungs"][0]["R_holds"], "A19: bb288 all-(R) decks"
            log("    bb288_yxx: y-quotient = (12,6) gross frame confirmed")

    # ------------------------------------------------- [1b] bb288 bottom
    # The bb288 tower bottoms out at (3,6) with kappa = 14: its ENTIRE
    # stabilizer group has 2^14 = 16,384 elements and its cycle space
    # 2^22 -- full enumeration replaces every census species.  This pins
    # the recursion base exactly (d of the (3,6) code) and tests the
    # "light logicals concentrate outside im p*" phenomenon once more.
    log("[1b] bb288/gross bottom layer (3,6): FULL enumeration")
    bot = TS.BBCode("btm36", (3, 6), "1 + y + y^2", "y^3 + x + x^2")
    bb72 = TS.BBCode("bb72", (6, 6), "x^3 + y + y^2", "y^3 + x + x^2")
    deck_b = make_deck(bb72, bot, 0, 3)
    imP_b, _ = rref_ints(list(TS._colspace(h1_map(deck_b))))
    stab_ints = list(span_points(bot.rsHX_b))
    swh: dict[int, int] = {}
    for x in stab_ints:
        w = bin(x).count("1")
        if 0 < w <= 16:
            swh[w] = swh.get(w, 0) + 1
    K = np.array(bot.kerHZ, dtype=np.uint8)           # 22 x 36
    dimk = K.shape[0]
    assert dimk == 22
    d_bot = None
    light_min: dict[int, int] = {}
    for lo in range(0, 1 << dimk, 1 << 18):
        idx = np.arange(lo, lo + (1 << 18), dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(dimk)) & 1).astype(np.uint8)
        V = (bits @ K) % 2
        wts = V.sum(axis=1)
        sigs = (V @ bot.zreps.T) % 2
        svals = sigs @ (1 << np.arange(bot.k, dtype=np.int64))
        mask = svals != 0
        for w, s in zip(wts[mask].tolist(), svals[mask].tolist()):
            if s not in light_min or w < light_min[s]:
                light_min[s] = w
    d_bot = min(light_min.values())
    attain = sorted(s for s, w in light_min.items() if w == d_bot)
    outside = [s for s in attain if not TS.in_span(s, imP_b,
               [(b & -b).bit_length() - 1 for b in imP_b])]
    log(f"    (3,6): d = {d_bot} EXACT (all 2^22 cycles enumerated); "
        f"stab weight hist <=16: {dict(sorted(swh.items()))}")
    log(f"    {len(attain)} min-weight classes, {len(outside)} outside "
        f"im p* (bb72 rung)  [light-outside-im-p* test]")
    out["bb288_bottom_demo"] = {
        "d_exact": int(d_bot), "stab_whist_le16": dict(sorted(swh.items())),
        "min_classes": len(attain), "min_classes_outside_imp": len(outside),
        "per_class_minima_complete": True}

    # ------------------------------------------------- [2] no-deck rows
    log("[2] structural anti-instances (odd |G|: L1 fails, no deck)")
    for name, lm, note in NO_DECK:
        assert (lm[0] * lm[1]) % 2 == 1
        out["no_deck"].append({"name": name, "lm": list(lm), "note": note})
        log(f"    {name} ({lm[0]},{lm[1]}): |G| = {lm[0]*lm[1]} odd -> "
            f"no Z2 deck. {note}")

    # ------------------------------------------------- [3] parity demo
    log("[3] parity-scope demo: |A| even breaks Lemma 2")
    demo = None
    for glm, As, Bs in [((6, 6), "1 + x", "1 + y"),
                        ((6, 6), "1 + x + y + x*y", "y^3 + x + x^2"),
                        ((4, 4), "1 + x", "1 + y")]:
        code = TS.BBCode("paritydemo", glm, As, Bs)
        odd = [kv for kv in code.kerHZ if int(kv.sum()) % 2 == 1]
        if odd and code.k > 0:
            demo = {"lm": list(glm), "A": As, "B": Bs, "k": code.k,
                    "odd_cycle_weight": int(odd[0].sum())}
            log(f"    ({glm[0]},{glm[1]}) A={As!r} B={Bs!r}: k={code.k}, "
                f"odd-weight cycle of weight {int(odd[0].sum())} EXISTS "
                f"-> Lemma 2 fails outside the odd-|A|,|B| scope")
            break
    out["parity_demo"] = demo

    out["wall_s"] = round(time.monotonic() - T0, 1)
    (DATA / "screen.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA / 'screen.json'}")


if __name__ == "__main__":
    main()
