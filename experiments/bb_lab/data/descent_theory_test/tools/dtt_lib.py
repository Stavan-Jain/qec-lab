"""descent_theory_test — shared library for Phase 0 (cohort) + Phase 1
(frozen predictions) of the falsification-shaped test of the A32/A33/A35
tower-descent theory.

Everything here is PREDICTION-side machinery.  Contamination rules
(hard-coded below, see also PROTOCOL.md):
  * NO exact-distance work on any Level-0 (target) code, ever.
  * Exact distance may be consumed only for
      (a) quotient levels (j >= 1) with n <= FULL_ENUM_MAX_N via trivial
          full enumeration,
      (b) values already recorded in the corpus DB (read-only),
      (c) cheap randomized L1 d_ub sampling (upper-bound witnesses only).
  * The corpus DB in the MAIN checkout is opened read_only=True.

The structural screen is the A35 screen (extracted at 2f063b0, validated
against the banked anchors in validation/screen_run.log) driven per-row:
this module imports its functions rather than reimplementing them.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent                      # data/descent_theory_test
LAB = OUT.parent.parent                # experiments/bb_lab
assert (LAB / "src" / "bb_lab").is_dir(), f"bad path bootstrap: {LAB}"
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
sys.path.insert(0, str(HERE))

import a32_tower_slice as TS                      # noqa: E402
from a30_rung_pass import rref_ints, v2i          # noqa: E402
from a35_generality_screen import (               # noqa: E402
    census_nodes, fold_terms, is_twisted, make_deck, poly_terms,
    screen_rung, span_dim, span_points, intersect_dim, terms_str,
)
from bb_lab.automorphism import automorphisms     # noqa: E402
from bb_lab.canonical import build_perm_table, canonical_bits, \
    _bits_to_support, _support_to_bits            # noqa: E402
from bb_lab.checks import bb_check_matrices       # noqa: E402
from bb_lab.group import AbelianGroup             # noqa: E402
from bb_lab.l1_sampling import l1_distance_ub     # noqa: E402
from bb_lab.poly import Poly                      # noqa: E402
from bb_lab.store import canonical_hash           # noqa: E402

MAIN_DB = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/"
               "bb_instances.duckdb")

# ---- frozen protocol constants (mirrored in PROTOCOL.md) -----------------
FULL_ENUM_MAX_N = 40          # quotient full-enumeration ceiling (rule a)
L1_SAMPLES_QUOTIENT = 20_000  # rule (c) sample count for quotient d_ub
L1_SAMPLES_TARGET = 30_000    # rule (c) for fresh target rows (parity strip)
L1_SEED = 20260817
EXACT_MU_MAX_KAPPA = 18       # exact stabilizer-span mu when 2^kappa small
MU_SAMPLE_TRIES = 400
GREEN_NODES = 2e11            # A35 G1 threshold (verbatim)
AMBER_NODES = 1e14
GREEN_CAP = 8                 # A35 G2 demonstrated envelope
AMBER_CAP = 12
NODES_PER_SEC = 1.1e9         # A33 lane anchor: 6.6e10 nodes ~= 61 s
RUNG_OVERHEAD_S = 20.0        # A33 (~22 s/rung) / A36 (~10 s/rung) anchor
RNG = np.random.default_rng(20260817)


# ------------------------------------------------------------------ helpers
def v2(x: int) -> int:
    c = 0
    while x % 2 == 0 and x > 0:
        x //= 2
        c += 1
    return c


def bbcode_k(lm: tuple[int, int], tA, tB) -> int:
    """k of the BB code without constructing TS.BBCode (which asserts on
    k = 0).  Same lab convention: k = dim ker HZ - rank HX."""
    from bb_lab.checks import circulant
    G = AbelianGroup(lm)
    A = Poly(support=frozenset(tA), group=G)
    B = Poly(support=frozenset(tB), group=G)
    MA = circulant(A).astype(np.uint8) % 2
    MB = circulant(B).astype(np.uint8) % 2
    HX = np.concatenate([MA, MB], axis=1) % 2
    HZ = np.concatenate([MB.T, MA.T], axis=1) % 2
    kerHZ = TS.kernel_basis(HZ)
    rs, _ = rref_ints([v2i(r) for r in HX])
    return len(kerHZ) - len(rs)


# ----------------------------------------------------------- descent ladder
def descent_route(ell: int, m: int) -> list[tuple[int, int]]:
    """The FROZEN canonical route: all y-axis Z2 folds first, then all
    x-axis Z2 folds (matches the banked ibm288Y y,y / bb288 y,x,x /
    bravyi360 y,x conventions).  Returns [(axis, newmod), ...]."""
    folds = []
    mm = m
    while mm % 2 == 0:
        mm //= 2
        folds.append((1, mm))
    ee = ell
    while ee % 2 == 0:
        ee //= 2
        folds.append((0, ee))
    return folds


def build_ladder(ell: int, m: int, As: str, Bs: str, name: str) -> dict:
    """Build the maximal descent ladder along the frozen route, truncating
    when a level's k hits 0 (TS.BBCode is undefined there) or when
    construction fails (degenerate axis).  Returns levels metadata, the
    constructed TS.BBCode chain, and the folds actually used."""
    lm = (ell, m)
    terms = [(lm, poly_terms(As, lm), poly_terms(Bs, lm))]
    route = descent_route(ell, m)
    for axis, newmod in route:
        plm, pA, pB = terms[-1]
        nlm = (newmod, plm[1]) if axis == 0 else (plm[0], newmod)
        terms.append((nlm, fold_terms(pA, axis, newmod),
                      fold_terms(pB, axis, newmod)))

    codes: list[TS.BBCode] = []
    lvl_rows: list[dict] = []
    truncated = None
    for i, (glm, tA, tB) in enumerate(terms):
        if not tA or not tB:
            truncated = {"at_level": i, "lm": list(glm), "k": None,
                         "reason": "polynomial vanishes under fold "
                                   "(all terms cancel mod 2)"}
            break
        k_i = bbcode_k(glm, tA, tB)
        if k_i <= 0:
            truncated = {"at_level": i, "lm": list(glm), "k": int(k_i),
                         "reason": "k=0 quotient (no logical content; "
                                   "TS machinery undefined)"}
            break
        try:
            code = TS.BBCode(f"{name}/L{i}", glm, terms_str(tA),
                             terms_str(tB))
        except Exception as e:                      # degenerate axis etc.
            truncated = {"at_level": i, "lm": list(glm), "k": int(k_i),
                         "reason": f"BBCode construction failed: {e}"}
            break
        codes.append(code)
        kappa = len(code.rsHX_b)
        avail = [ax for ax, o in zip("xy", glm) if o % 2 == 0]
        lvl_rows.append({
            "level": i, "lm": list(glm), "n": code.n, "k": code.k,
            "A": terms_str(tA), "B": terms_str(tB),
            "wA": len(tA), "wB": len(tB),
            "parity_scope": bool(len(tA) % 2 and len(tB) % 2),
            "cycles_all_even": bool(all(int(kv.sum()) % 2 == 0
                                        for kv in code.kerHZ)),
            "kappa": kappa, "available_decks": avail,
        })
    folds_used = route[:len(codes) - 1] if len(codes) > 1 else []
    return {"levels": lvl_rows, "codes": codes, "folds": folds_used,
            "route_full": route, "truncated": truncated}


# ------------------------------------------------------------- mu probing
def mu_probe(code: TS.BBCode) -> dict:
    """Lightest nonzero X-stabilizer weight at this level.  Exact by
    stabilizer-span enumeration when kappa <= EXACT_MU_MAX_KAPPA, else
    sampled (single rows + localized row sums; upper bound only)."""
    kappa = len(code.rsHX_b)
    row_w = int(code.HX[0].sum())
    if kappa <= EXACT_MU_MAX_KAPPA:
        best = row_w
        for x in span_points(code.rsHX_b):
            if x:
                w = int(x).bit_count()
                if w < best:
                    best = w
        return {"mu": int(best), "provenance": "exact-span",
                "kappa": kappa}
    lx, ly = code.G.orders
    best = row_w
    for _ in range(MU_SAMPLE_TRIES):
        j = int(RNG.integers(2, 6))
        g0 = (int(RNG.integers(lx)), int(RNG.integers(ly)))
        beta = np.zeros(code.n, dtype=np.uint8)
        for _ in range(j):
            da, db = int(RNG.integers(-2, 3)), int(RNG.integers(-2, 3))
            gi = code.G.index(((g0[0] + da) % lx, (g0[1] + db) % ly))
            beta = (beta + code.HX[gi]) % 2
        wb = int(beta.sum())
        if 0 < wb < best:
            best = wb
    return {"mu": int(best), "provenance": "sampled-ub", "kappa": kappa}


# ------------------------------------------------------- structural screen
def _preimage_safe(M: np.ndarray, Wb: list[int], Wp: list[int]) -> list[int]:
    """TS._preimage with its degenerate case handled (AMENDMENTS.md #1).

    When span(Wb) is the FULL codomain, the annihilator basis is empty and
    ``np.array([])`` loses 2-D shape inside TS._preimage, crashing the
    matmul; mathematically ker(0-map) is the whole domain.  Byte-for-byte
    the same computation as TS._preimage otherwise (same i2v/v2i/
    kernel_basis conventions from a32_tower_slice).
    """
    n = M.shape[0]
    Wmat = (np.array([TS.i2v(x, n) for x in Wb], dtype=np.uint8)
            if Wb else np.zeros((0, n), dtype=np.uint8))
    F = np.array(TS.kernel_basis(Wmat), dtype=np.uint8)
    if F.size == 0:
        return [TS.v2i(v) for v in np.eye(M.shape[1], dtype=np.uint8)]
    QM = (F @ M) % 2
    return [TS.v2i(v) for v in TS.kernel_basis(QM)]


def regime_label(pair: dict) -> str:
    """A35 SS3 pair-regime taxonomy from the measured S/K lattice."""
    dS, dK, dSK = pair["dim_S"], pair["dim_K"], pair["dim_SK"]
    if dK == dSK and dS == dSK:
        return "R3"                      # S = K
    if dSK == 0:
        return "R2"                      # one-branch descent
    if pair["K_in_S"] and dK < dS:
        return "R1"                      # K proper subset of S
    return "R4"                          # partial overlap


def screen_structure(ell: int, m: int, As: str, Bs: str,
                     name: str = "row") -> dict:
    """Per-row structural screen: ladder + per-rung A35 measurements +
    pair regimes + mu.  NO distance work happens here."""
    t0 = time.monotonic()
    ladder = build_ladder(ell, m, As, Bs, name)
    codes = ladder["codes"]
    out: dict = {
        "levels": ladder["levels"],
        "folds": [[ax, nm] for ax, nm in ladder["folds"]],
        "route_full": [[ax, nm] for ax, nm in ladder["route_full"]],
        "truncated": ladder["truncated"],
        "depth_available": len(ladder["route_full"]),
        "depth_used": len(codes) - 1,
    }
    rungs = []
    lemma1_failures = []
    for i, (axis, newmod) in enumerate(ladder["folds"]):
        try:
            r = screen_rung(codes[i], codes[i + 1], axis, newmod,
                            W_eff=0, do_fibers=False)
        except AssertionError as e:      # Lemma 1 is a theorem: loud
            lemma1_failures.append({"rung": i, "error": str(e)})
            rungs.append({"rung": i, "axis": "xy"[axis],
                          "fold_to": newmod, "LEMMA1_FAILED": str(e)})
            continue
        r["rung"] = i
        # rank-law bookkeeping (A35 SS3): predicted values + violation flag
        kc = codes[i].k
        kb = codes[i + 1].k
        r["rank_law_predicted"] = kc // 2
        r["rank_law_holds"] = bool(r["rank_p"] == kc // 2
                                   and r["rank_tau"] == kc // 2)
        r["exact_cover_predicted"] = True          # theorem (transfer LES)
        r["R_predicted_iff"] = "exact_base and sigma_id iff R_holds (A12)"
        r["codim_lift_predicted"] = kb - kc // 2
        r["codim_lift_law_holds"] = bool(r["codim_lift"] == kb - kc // 2)
        rungs.append(r)
    out["rungs_lemma1_failures"] = lemma1_failures

    pairs = []
    for i in range(len(rungs) - 1):
        ra, rb = rungs[i], rungs[i + 1]
        if "LEMMA1_FAILED" in ra or "LEMMA1_FAILED" in rb:
            continue
        S = ra["_imP"]
        K = rb["_kerP"]
        mid, bot = codes[i + 1], codes[i + 2]
        if mid.k == 0 or not S:
            pairs.append({"i": i, "note": "H1(mid) trivial"})
            continue
        dS, dK = span_dim(S), span_dim(K)
        dSK = intersect_dim(S, K)
        if bot.k > 0:
            from a33_tower_cells import h1_map
            Mb = h1_map(make_deck(mid, bot, *ladder["folds"][i + 1]))
            Wimg = [TS._apply(Mb, s)
                    for s in span_points(rref_ints(list(S))[0])]
            dW = span_dim(Wimg)
            Wb2, Wp2 = rref_ints([w for w in Wimg if w])
            pre = _preimage_safe(Mb, Wb2, Wp2)
            reach = bool(TS._span_eq(pre, rref_ints(list(S))[0])) \
                if dSK == dK else False
        else:
            dW, reach = 0, None
        K_in_S = bool(all(TS.in_span(x, *rref_ints(list(S))) for x in K)) \
            if K else True
        p = {"i": i, "dim_S": dS, "dim_K": dK, "dim_SK": dSK,
             "K_in_S": K_in_S, "one_branch": bool(dSK == 0),
             "dim_W": dW, "reach_preimage_eq_S": reach}
        p["regime"] = regime_label(p)
        pairs.append(p)
    out["pairs"] = pairs
    out["regimes"] = sorted({p["regime"] for p in pairs if "regime" in p})

    mus = []
    for c in codes:
        mus.append(mu_probe(c))
    out["mu_per_level"] = mus
    out["mu_min"] = min((m_["mu"] for m_ in mus), default=None)

    for r in rungs:                      # strip non-JSON internals
        r.pop("_imP", None)
        r.pop("_kerP", None)
    out["rungs"] = rungs
    out["screen_wall_s"] = round(time.monotonic() - t0, 2)
    return out


# ------------------------------------------------------------ corpus lookup
_PERM_CACHE: dict[tuple[int, int], tuple] = {}
_LOOKUP_CACHE: dict[str, dict | None] = {}
_CORPUS_LABELS: set[str] | None = None


def _perms(lm: tuple[int, int]):
    if lm not in _PERM_CACHE:
        G = AbelianGroup(lm)
        _PERM_CACHE[lm] = (G, build_perm_table(G, auts=automorphisms(G)))
    return _PERM_CACHE[lm]


def canonical_id(lm: tuple[int, int], tA, tB) -> tuple[str, str, str, str]:
    """(instance_id, group_label, canonical A string, canonical B string)
    under the corpus's Aut x translation x swap dedup."""
    G, perms = _perms(lm)
    A_bits = _support_to_bits(set(tA), G)
    B_bits = _support_to_bits(set(tB), G)
    cA, cB, _ = canonical_bits(A_bits, B_bits, perms)
    A_str = Poly(support=frozenset(_bits_to_support(cA, G)),
                 group=G).canonical_string()
    B_str = Poly(support=frozenset(_bits_to_support(cB, G)),
                 group=G).canonical_string()
    label = G.label()
    return canonical_hash(label, A_str, B_str), label, A_str, B_str


def corpus_labels(con) -> set[str]:
    global _CORPUS_LABELS
    if _CORPUS_LABELS is None:
        _CORPUS_LABELS = {r[0] for r in con.execute(
            "select distinct group_struct from bb_instances").fetchall()}
    return _CORPUS_LABELS


def corpus_lookup(con, lm: tuple[int, int], tA, tB) -> dict | None:
    """Find the corpus row for this (group, A, B) up to the corpus's
    canonical equivalence, trying both axis orientations."""
    orientations = [(lm, tA, tB)]
    if lm[0] != lm[1]:
        sw = (lm[1], lm[0])
        orientations.append((sw, frozenset((b, a) for a, b in tA),
                             frozenset((b, a) for a, b in tB)))
    labels = corpus_labels(con)
    for olm, oA, oB in orientations:
        if AbelianGroup(olm).label() not in labels:
            continue
        iid, label, cA, cB = canonical_id(olm, oA, oB)
        key = iid
        if key in _LOOKUP_CACHE:
            hit = _LOOKUP_CACHE[key]
        else:
            row = con.execute(
                "select instance_id, group_struct, n, k, d_lb, d_ub, "
                "d_exact, d_method from bb_instances where instance_id = ?",
                [iid]).fetchone()
            hit = None if row is None else {
                "instance_id": row[0], "group_struct": row[1],
                "n": row[2], "k": row[3], "d_lb": row[4], "d_ub": row[5],
                "d_exact": row[6], "d_method": row[7]}
            _LOOKUP_CACHE[key] = hit
        if hit is not None:
            return hit
    return None


# --------------------------------------------------- quotient distance info
_ENUM_CACHE: dict[str, dict] = {}


def full_enumerate_distance(code: TS.BBCode) -> dict:
    """EXACT d by enumerating the whole X-cycle space (allowed only for
    quotient levels with n <= FULL_ENUM_MAX_N).  Pattern = A35 screen
    section [1b]."""
    assert code.n <= FULL_ENUM_MAX_N, "full enumeration size guard"
    K = np.array(code.kerHZ, dtype=np.uint8)
    dimk = K.shape[0]
    step = 1 << min(18, dimk)
    d_min = None
    for lo in range(0, 1 << dimk, step):
        idx = np.arange(lo, lo + step, dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(dimk)) & 1).astype(np.uint8)
        V = (bits @ K) % 2
        wts = V.sum(axis=1)
        sigs = (V @ code.zreps.T) % 2
        nz = sigs.any(axis=1)
        if nz.any():
            w = int(wts[nz].min())
            d_min = w if d_min is None else min(d_min, w)
    return {"d": int(d_min), "method": "full-enumeration",
            "space_dim": int(dimk)}


def quotient_distance(con, lm, tA, tB, level_index: int) -> dict:
    """Distance info for a QUOTIENT level (contamination guard:
    level_index >= 1).  Priority: corpus -> full enumeration (n <= 40)
    -> L1 d_ub sampling.  Constructs the TS.BBCode lazily (only the
    full-enumeration branch needs it)."""
    assert level_index >= 1, \
        "contamination guard: no distance work on Level-0 targets"
    iid, label, cA, cB = canonical_id(lm, tA, tB)
    if iid in _ENUM_CACHE:
        return _ENUM_CACHE[iid]
    hit = corpus_lookup(con, lm, tA, tB)
    n_code = 2 * lm[0] * lm[1]
    if hit is not None and hit["d_exact"] is not None:
        res = {"d": int(hit["d_exact"]), "d_is_exact": True,
               "provenance": f"corpus-exact ({hit['d_method']})",
               "corpus_id": hit["instance_id"]}
    elif n_code <= FULL_ENUM_MAX_N:
        code = TS.BBCode(f"q{label}", lm, terms_str(tA), terms_str(tB))
        e = full_enumerate_distance(code)
        res = {"d": e["d"], "d_is_exact": True,
               "provenance": "full-enumeration (n<=40 quotient)",
               "corpus_id": hit["instance_id"] if hit else None}
    elif hit is not None and hit["d_ub"] is not None:
        res = {"d": int(hit["d_ub"]), "d_is_exact": False,
               "provenance": "corpus-d_ub (L1 sampling witness)",
               "corpus_id": hit["instance_id"]}
    else:
        G = AbelianGroup(lm)
        A = Poly(support=frozenset(tA), group=G)
        B = Poly(support=frozenset(tB), group=G)
        r = l1_distance_ub(bb_check_matrices(A, B),
                           n_samples=L1_SAMPLES_QUOTIENT, seed=L1_SEED)
        res = {"d": int(r.distance_ub), "d_is_exact": False,
               "provenance": f"fresh-L1-d_ub ({L1_SAMPLES_QUOTIENT} "
                             "samples, witness only)",
               "corpus_id": None}
    res["instance_id_canonical"] = iid
    res["group_label"] = label
    _ENUM_CACHE[iid] = res
    return res


# ------------------------------------------------------------ cost + gates
def cost_block(levels: list[dict], W: int, mu_min: int) -> dict:
    """A35 cost gates at window W, verbatim thresholds."""
    assert levels, "cost_block needs a non-empty ladder (guard upstream)"
    lvl_nodes = [census_nodes(lv["kappa"], W) for lv in levels]
    nb, nt = lvl_nodes[-1], lvl_nodes[0]
    cap_max = max(0, (W - mu_min)) // 2
    verdict = ("GREEN" if nb <= GREEN_NODES and cap_max <= GREEN_CAP else
               "AMBER" if nb <= AMBER_NODES and cap_max <= AMBER_CAP
               else "RED")
    dispatch = max(2 ** lv["k"] for lv in levels[1:]) if len(levels) > 1 \
        else None
    t_est = nb / NODES_PER_SEC + RUNG_OVERHEAD_S * (len(levels) - 1)
    return {
        "W": W, "mu_min": mu_min, "cap_max": cap_max, "verdict": verdict,
        "log10_nodes_per_level":
            [round(math.log10(x), 1) if x else None for x in lvl_nodes],
        "bottom_nodes": float(nb),
        "win_factor_vs_top": round(nt / nb, 1) if nb else None,
        "sector_dispatch_max": dispatch,
        "wall_estimate_s": round(t_est, 1),
        "wall_model": f"bottom_nodes/{NODES_PER_SEC:.1e} nodes-per-s "
                      f"(A33 anchor) + {RUNG_OVERHEAD_S:.0f} s/rung; "
                      "order-of-magnitude only",
    }


def g5_window(levels: list[dict], dinfo: dict[int, dict],
              d_ub_top: int | None, parity_ok: bool) -> dict:
    """G5 tau-branch ceiling: certifiable-d window
    [2, min(d_ub_top, min_j 2^j d(Lj))] with per-term provenance."""
    terms = []
    for j in range(1, len(levels)):
        di = dinfo.get(j)
        if di is None:
            continue
        terms.append({"level": j, "d": di["d"],
                      "d_is_exact": di["d_is_exact"],
                      "provenance": di["provenance"],
                      "term": (2 ** j) * di["d"]})
    ceiling = min((t["term"] for t in terms), default=None)
    any_bounded = any(not t["d_is_exact"] for t in terms)
    stalls = []
    for j in range(1, len(levels) - 1):
        a, b = dinfo.get(j), dinfo.get(j + 1)
        if a and b and a["d_is_exact"] and b["d_is_exact"] \
                and a["d"] > 2 * b["d"]:
            stalls.append({"rung": j, "d_cover": a["d"], "d_base": b["d"],
                           "note": f"chain violation: d(L{j}) = {a['d']} "
                                   f"> 2 d(L{j+1}) = {2*b['d']}; route "
                                   f"floor at L{j} caps at {2*b['d']}"})
    hi = None
    if ceiling is not None and d_ub_top is not None:
        hi = min(d_ub_top, ceiling)
    elif ceiling is not None:
        hi = ceiling
    elif d_ub_top is not None:
        hi = d_ub_top
    if hi is not None and parity_ok:
        hi -= hi % 2
    return {"terms": terms, "ceiling": ceiling,
            "ceiling_is_upper_estimate": any_bounded,
            "window_lo": 2, "window_hi": hi,
            "chain_stalls": stalls,
            "note": ("window_hi uses bounded-only d at some level; the "
                     "true ceiling can only be lower" if any_bounded
                     else "all chain terms exact")}
