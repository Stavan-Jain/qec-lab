"""descent_theory_test Phase 2 — the generic tower-certificate closure engine.

Runs ONE pre-registered question ("certify d = T via the frozen route at
priced W = T - parity_slack") for one cohort row, exactly along the frozen
route recorded in the prediction (PROTOCOL.md §4).  The method is the
A30/A32/A33 architecture, generalized:

  1. bottom-level coset-BZ census over ALL 2^k logical sectors at window W
     (bb_lab.cosetbz: counting-invariant two-window Brouwer-Zimmermann;
     exact node-count asserts ON) -> the complete set of nonzero cycles of
     weight <= W at the bottom level;
  2. rung-by-rung transport up the tower: every cover cycle z with
     |z| <= W has slice data (beta = P z, m2) with |z| = |beta| + 2 m2, so
       beta != 0:  z is a bounded-overflow lift of a censused base cycle
                   (a32_subclosures.enumerate_lifts, complete by the
                   exact-off-support subset-sum argument; deep MITM lane
                   for caps 5..8 from a32_deep_fibers);
       beta  = 0:  z = tau(gamma) with gamma in ker E, |gamma| = |z|/2
                   (checked at runtime: ker E == base cycle space; else
                   the row stops with an honest engine-limit record);
     translation-orbit compression at every level (fold of translations is
     onto, spot-checked covariance per rung), sheet flips are deck
     translations and stay inside orbits;
  3. tiny levels (cycle-space dim <= FULL_ENUM_DIM_MAX) are replaced by
     direct full enumeration of the cycle space (complete by construction;
     the bottom BZ census is still executed as the criterion-(ii) datum);
  4. at the top: any non-stabilizer cycle of weight <= W is a verified
     COUNTEREXAMPLE (support recorded); none -> certificate-tier floor
     d >= W+1 (+1 more when the top cycle space is even-weight only).

Structural rung measurements for criterion (i) (rank p*, rank tau*,
exactness, sigma*=id, codim_lift) are re-measured here with the same
extracted-screen functions the predictions used.

Trust tier of a success: "certificate (counting-invariant coset-BZ census
+ deck-transport lifting; not kernel-checked)".  Witness weights are never
floors; UNKNOWN is UNKNOWN.

Provenance of copied code (copies, not forks of behavior):
  - enumerate_lifts / enumerate_lifts_deep: byte-level math and asserts of
    scripts/a32_subclosures.py::enumerate_lifts and
    scripts/a32_deep_fibers.py::enumerate_lifts_deep, with the
    deck-constant tables (E columns, RHS) hoisted out of the per-fiber
    loop into DeckTables (pure performance hoist; every assert kept).
  - batch_keys: scripts/a32_sectorAC_full.py::batch_keys verbatim.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DTT = HERE.parent                       # data/descent_theory_test
LAB = DTT.parent.parent                 # experiments/bb_lab
assert (LAB / "src" / "bb_lab").is_dir(), f"bad path bootstrap: {LAB}"
for p in (LAB / "src", LAB / "scripts", DTT / "tools", HERE):
    sys.path.insert(0, str(p))

import a32_tower_slice as TS                            # noqa: E402
from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402
from a35_generality_screen import (                     # noqa: E402
    census_nodes, fold_terms, liftable_codim, make_deck, poly_terms,
    terms_str,
)
from a33_tower_cells import h1_map, translation_mat     # noqa: E402
from bb_lab.cosetbz import (                            # noqa: E402
    NOFF_MAX, build_kernel, coset_base, disjoint_info_sets, pair_radii,
    run_window, unpack3,
)

FULL_ENUM_DIM_MAX = 22        # cycle-space dim <= this -> direct enumeration
CAP_ENVELOPE = 8              # demonstrated fiber-lane envelope (A35 G2)
HIT_GUARD = 4_000_000         # max expected/realized census or level set size
DEPS_GUARD = 18               # per-fiber kernel-dimension guard (2^deps enum)


def log(msg: str) -> None:
    print(msg, flush=True)


# ------------------------------------------------------------- copied utils
def batch_keys(vecs: np.ndarray, perms: list[np.ndarray]) -> np.ndarray:
    """Lexicographic-min packed key over translation perms, per row.
    [verbatim from scripts/a32_sectorAC_full.py]"""
    N, n = vecs.shape
    nwords = -(-(-(-n // 8)) // 8)  # ceil(ceil(n/8)/8)
    cur = None
    for p in perms:
        t = np.packbits(vecs[:, p], axis=1)
        pad = np.zeros((N, nwords * 8 - t.shape[1]), dtype=np.uint8)
        t64 = np.ascontiguousarray(
            np.concatenate([t, pad], axis=1)).view(">u8").reshape(N, nwords)
        if cur is None:
            cur = t64.copy()
        else:
            better = np.zeros(N, dtype=bool)
            tied = np.ones(N, dtype=bool)
            for w in range(nwords):
                better |= tied & (t64[:, w] < cur[:, w])
                tied &= t64[:, w] == cur[:, w]
            cur[better] = t64[better]
    return cur


class DeckTables:
    """Deck-constant tables for the lift enumerators (perf hoist only)."""

    def __init__(self, deck: TS.Deck):
        self.deck = deck
        n = deck.base.n
        self.n = n
        self.E_cols = [v2i(deck.E[:, j]) for j in range(n)]
        self.RHS = deck.RHS


def enumerate_lifts(dt: DeckTables, beta: np.ndarray, cap: int,
                    kernel_cap: int = 16) -> dict[int, int]:
    """All v0 with E v0 = RHS beta and |v0 off supp(beta)| <= cap.
    [math + asserts of a32_subclosures.enumerate_lifts; E_cols/RHS hoisted]
    """
    assert cap <= 4, "MITM lane implemented to size 4"
    n = dt.n
    E_cols = dt.E_cols
    rhs = (dt.RHS @ beta) % 2
    rhs_i = v2i(rhs)
    bsupp = [int(j) for j in np.nonzero(beta)[0]]
    bmask = v2i(beta)
    bcols = [E_cols[j] for j in bsupp]
    bb, bp = rref_ints(bcols)
    rhs_res = reduce_int(rhs_i, bb, bp)
    offb = [j for j in range(n) if not (bmask >> j) & 1]
    red = {j: reduce_int(E_cols[j], bb, bp) for j in offb}
    by_val: dict[int, list[int]] = {}
    for j in offb:
        by_val.setdefault(red[j], []).append(j)
    hits_X: set[tuple[int, ...]] = set()
    if rhs_res == 0:
        hits_X.add(())
    if cap >= 1:
        for j in by_val.get(rhs_res, []):
            hits_X.add((j,))
    if cap >= 2:
        for j1 in offb:
            for j2 in by_val.get(rhs_res ^ red[j1], []):
                if j2 > j1:
                    hits_X.add((j1, j2))
    if cap >= 3:
        for j1, j2 in itertools.combinations(offb, 2):
            for j3 in by_val.get(rhs_res ^ red[j1] ^ red[j2], []):
                if j3 > j2:
                    hits_X.add((j1, j2, j3))
    if cap >= 4:
        pair_sum: dict[int, list[tuple[int, int]]] = {}
        for j1, j2 in itertools.combinations(offb, 2):
            pair_sum.setdefault(red[j1] ^ red[j2], []).append((j1, j2))
        for val, prs in pair_sum.items():
            for j3, j4 in pair_sum.get(rhs_res ^ val, []):
                for j1, j2 in prs:
                    if j2 < j3:
                        hits_X.add((j1, j2, j3, j4))
    return _kernel_expand(dt, [bsupp + list(X) for X in sorted(hits_X)],
                          rhs_i, bmask, cap, kernel_cap, n)


def enumerate_lifts_deep(dt: DeckTables, beta: np.ndarray, cap: int,
                         kernel_cap: int = 20) -> dict[int, int]:
    """Deep lane, caps 5..8 (ordered-split MITM).
    [math + asserts of a32_deep_fibers.enumerate_lifts_deep; tables hoisted]
    """
    assert cap <= 8
    n = dt.n
    E_cols = dt.E_cols
    rhs = (dt.RHS @ beta) % 2
    rhs_i = v2i(rhs)
    bsupp = [int(j) for j in np.nonzero(beta)[0]]
    bmask = v2i(beta)
    bcols = [E_cols[j] for j in bsupp]
    bb, bp = rref_ints(bcols)
    rhs_res = reduce_int(rhs_i, bb, bp)
    offb = [j for j in range(n) if not (bmask >> j) & 1]
    red = {j: reduce_int(E_cols[j], bb, bp) for j in offb}
    half = (cap + 1) // 2
    by_size: list[dict[int, list[int]]] = [dict() for _ in range(half + 1)]
    by_size[0][0] = [0]
    for s in range(1, half + 1):
        for comb in itertools.combinations(offb, s):
            m = 0
            r = 0
            for j in comb:
                m |= 1 << j
                r ^= red[j]
            by_size[s].setdefault(r, []).append(m)
    hits_X: set[int] = set()
    for s in range(cap + 1):
        a, b = (s + 1) // 2, s // 2
        for asum, amasks in by_size[a].items():
            bucket = by_size[b].get(rhs_res ^ asum)
            if not bucket:
                continue
            for amask in amasks:
                if a:
                    alsb = (amask & -amask).bit_length() - 1
                for bmask2 in bucket:
                    if b == 0:
                        if a == s:
                            hits_X.add(amask)
                        continue
                    if bmask2.bit_length() - 1 < alsb \
                            and not (amask & bmask2):
                        hits_X.add(amask | bmask2)
    colsets = [bsupp + [j for j in range(n) if (X >> j) & 1]
               for X in sorted(hits_X)]
    return _kernel_expand(dt, colsets, rhs_i, bmask, cap, kernel_cap, n)


def _kernel_expand(dt: DeckTables, colsets, rhs_i: int, bmask: int,
                   cap: int, kernel_cap: int, n: int) -> dict[int, int]:
    """Shared per-X kernel enumeration (identical logic in both lanes)."""
    E_cols = dt.E_cols
    out: dict[int, int] = {}
    for cols in colsets:
        b3: list[int] = []
        p3: list[int] = []
        h3: list[int] = []
        deps: list[int] = []
        for ci, j in enumerate(cols):
            cur, h = E_cols[j], 1 << ci
            for bb3, pp3, hh in zip(b3, p3, h3):
                if (cur >> pp3) & 1:
                    cur ^= bb3
                    h ^= hh
            if cur:
                b3.append(cur)
                p3.append((cur & -cur).bit_length() - 1)
                h3.append(h)
            else:
                deps.append(h)
        cur, hsel = rhs_i, 0
        for bb3, pp3, hh in zip(b3, p3, h3):
            if (cur >> pp3) & 1:
                cur ^= bb3
                hsel ^= hh
        if cur:
            continue
        if len(deps) > DEPS_GUARD:
            raise EngineLimit(f"fiber kernel 2^{len(deps)} > 2^{DEPS_GUARD}")
        assert len(deps) <= kernel_cap, f"kernel 2^{len(deps)}"
        for kt in range(1 << len(deps)):
            sel = hsel
            for jj in range(len(deps)):
                if (kt >> jj) & 1:
                    sel ^= deps[jj]
            v0_int = 0
            for ci, j in enumerate(cols):
                if (sel >> ci) & 1:
                    v0_int |= 1 << j
            m2 = bin(v0_int & ~bmask).count("1")
            if m2 > cap:
                continue
            v0 = i2v(v0_int, n)
            # verify solution exactly as the originals do
            rhs_v = i2v(rhs_i, dt.deck.E.shape[0])
            assert not (((dt.deck.E @ v0) + rhs_v) % 2).any(), \
                "not a solution"
            canon = min(v0_int, v0_int ^ bmask)
            prev = out.get(canon)
            if prev is None or m2 < prev:
                out[canon] = m2
    return out


class EngineLimit(Exception):
    """Honest stop: the demonstrated lanes cannot complete this question."""


class Blowup(Exception):
    """Honest stop: state size beyond the engine's memory envelope."""


# --------------------------------------------------------------- the ladder
def build_codes(ell: int, m: int, As: str, Bs: str,
                folds: list[list[int]]) -> list[TS.BBCode]:
    """Literal descent along the FROZEN route (folds from the prediction
    record).  Any construction failure is an engine error (recorded)."""
    lm = (ell, m)
    terms = [(lm, poly_terms(As, lm), poly_terms(Bs, lm))]
    for axis, newmod in folds:
        plm, pA, pB = terms[-1]
        assert plm[axis] == 2 * newmod, \
            f"frozen fold {axis}->{newmod} is not index-2 from {plm}"
        nlm = (newmod, plm[1]) if axis == 0 else (plm[0], newmod)
        terms.append((nlm, fold_terms(pA, axis, newmod),
                      fold_terms(pB, axis, newmod)))
    codes = []
    for i, (glm, tA, tB) in enumerate(terms):
        codes.append(TS.BBCode(f"L{i}", glm, terms_str(tA), terms_str(tB)))
    return codes


def measure_rung(cover: TS.BBCode, base: TS.BBCode, deck: TS.Deck,
                 axis: int, newmod: int) -> dict:
    """Criterion-(i) structural measurements (same functions as the frozen
    screen; fibers OFF)."""
    r: dict = {"axis": "xy"[axis], "fold_to": newmod,
               "k_cover": cover.k, "k_base": base.k}
    if cover.k > 0 and base.k > 0:
        Mp = h1_map(deck)
        Mt = h1_map(deck, tau=True)
        r["rank_p"] = TS.gf2_rank([v2i(c) for c in Mp.T])
        r["rank_tau"] = TS.gf2_rank([v2i(c) for c in Mt.T])
        imP = TS._colspace(Mp)
        imT = TS._colspace(Mt)
        kerP = TS._kernel_ints(Mp)
        kerT = TS._kernel_ints(Mt)
        r["exact_cover"] = bool(TS._span_eq(imT, kerP))
        r["exact_base"] = bool(TS._span_eq(kerT, imP))
        sig_t = (newmod, 0) if axis == 0 else (0, newmod)
        St = translation_mat(cover, sig_t)
        r["sigma_id"] = bool((St == np.eye(cover.k, dtype=np.uint8)).all())
    else:
        r["rank_p"] = r["rank_tau"] = 0
        r["exact_cover"] = r["exact_base"] = None
        r["sigma_id"] = None
    r["codim_lift"] = liftable_codim(deck)
    r["rank_law_holds"] = bool(r["rank_p"] == cover.k // 2
                               and r["rank_tau"] == cover.k // 2)
    r["exactness_at_cover_holds"] = (bool(r["exact_cover"])
                                     if r["exact_cover"] is not None
                                     else None)
    return r


# ------------------------------------------------------------ level engines
def expected_hits(n: int, dim_cycles: int, W: int) -> float:
    """Expected number of cycles of weight <= W (binomial heuristic)."""
    tot = sum(math.comb(n, w) for w in range(1, W + 1))
    return tot * 2.0 ** (dim_cycles - n)


def full_enum_level(code: TS.BBCode, W: int) -> np.ndarray:
    """All nonzero cycles of weight <= W by direct enumeration of the
    cycle space (dim <= FULL_ENUM_DIM_MAX).  Complete by construction."""
    K = np.array(code.kerHZ, dtype=np.uint8)
    dimk = K.shape[0]
    assert dimk <= FULL_ENUM_DIM_MAX
    step = 1 << min(16, dimk)
    outs = []
    for lo in range(0, 1 << dimk, step):
        idx = np.arange(lo, lo + step, dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(dimk)) & 1).astype(np.uint8)
        V = (bits @ K) % 2
        wts = V.sum(axis=1)
        m = (wts > 0) & (wts <= W)
        if m.any():
            outs.append(V[m])
    if not outs:
        return np.zeros((0, code.n), dtype=np.uint8)
    return np.concatenate(outs, axis=0)


def bottom_census(code: TS.BBCode, W: int, threads: int, workdir: Path,
                  deadline: float) -> tuple[np.ndarray, dict]:
    """Coset-BZ census over ALL 2^k sectors at window W.  Returns
    (vectors, stats).  Node counts asserted exact inside run_window."""
    t0 = time.monotonic()
    binp = build_kernel()
    k = code.k
    assert k <= 8, f"bottom k = {k} > 8 (2^k sectors exceed one walk)"
    # sector representatives: all xrep combinations (incl. 0 = stab sector)
    reps = []
    for s in range(1 << k):
        v = np.zeros(code.n, dtype=np.uint8)
        for i in range(k):
            if (s >> i) & 1:
                v ^= code.xreps[i]
        reps.append(v)
    est = expected_hits(code.n, len(code.kerHZ), W)
    if est > HIT_GUARD:
        raise Blowup(f"census expected hits {est:.2e} > {HIT_GUARD:.0e}")
    single_window = False
    try:
        I1, G1, I2, G2, kappa = disjoint_info_sets(code.HX)
        windows = [(I1, G1), (I2, G2)]
        r1, r2 = pair_radii(W)
        radii = [r1, r2]
    except RuntimeError:
        # degenerate column matroid (no two disjoint bases reachable):
        # fall back to ONE window at r = W — complete unconditionally,
        # since |c|_I <= |c| <= W for every coset element of weight <= W.
        # Recorded as a fallback census; excluded from the criterion-(ii)
        # formula sample (different r-pair semantics), see PROTOCOL (ii).
        from bb_lab.cosetbz import rref as _rref
        R1n, piv0 = _rref(code.HX)
        kappa = len(piv0)

        def _systematic(window):
            S, p = _rref(np.concatenate([R1n[:, window], R1n], axis=1))
            assert p[: len(window)] == list(range(len(window)))
            return S[:, len(window):]

        windows = [(piv0, _systematic(piv0))]
        radii = [min(W, kappa)]
        if sum(math.comb(kappa, s) for s in range(1, radii[0] + 1)) > 5e8:
            raise EngineLimit(
                "single-window fallback census beyond node envelope")
        single_window = True
    r_used = radii
    hit_ints: set[int] = set()
    stats = {"kappa": kappa, "W": W, "r_pair": r_used, "windows": [],
             "n_offsets": len(reps), "single_window_fallback":
             single_window}
    for wi, ((window, Gs), r) in enumerate(zip(windows, radii)):
        bases = []
        for tv in reps:
            cb = coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            if 0 < wcb <= W:
                hit_ints.add(v2i(cb))
            bases.append(cb)
        assert len(bases) <= NOFF_MAX
        if r == 0:
            # r = 0 window: only the coset-base elements themselves
            # (collected above); the C kernel's walk is skipped (its
            # RMAX = 0 edge over-walks single rows — engine edge case).
            stats["windows"].append({"window": wi, "r": 0, "nodes": 0,
                                     "expect": 0, "wall_s": 0.0,
                                     "note": "r=0: bases only, walk "
                                             "skipped"})
            continue
        res = run_window(binp, f"cen_w{wi}", Gs, bases, r, W, deadline,
                         threads=threads, workdir=workdir)
        for _j, hx in res.pop("hit_rows"):
            v = unpack3(hx, code.n)
            if v.any():
                hit_ints.add(v2i(v))
        if len(hit_ints) > HIT_GUARD:
            raise Blowup(f"census hits > {HIT_GUARD:.0e}")
        stats["windows"].append({"window": wi, "r": r,
                                 "nodes": res["nodes"],
                                 "expect": res["expect"],
                                 "wall_s": res["wall_s"]})
    stats["nodes_total"] = sum(w["nodes"] for w in stats["windows"])
    stats["formula_nodes"] = census_nodes(kappa, W)
    stats["node_ratio_vs_formula"] = (
        None if single_window else
        round(stats["nodes_total"] / stats["formula_nodes"], 6)
        if stats["formula_nodes"] else None)
    stats["hits"] = len(hit_ints)
    stats["wall_s"] = round(time.monotonic() - t0, 2)
    vecs = np.array([i2v(h, code.n) for h in sorted(hit_ints)],
                    dtype=np.uint8) if hit_ints \
        else np.zeros((0, code.n), dtype=np.uint8)
    # spot soundness: every censused element is a cycle
    for row in vecs[:: max(1, len(vecs) // 20)]:
        assert code.is_cycle(row), "census element is not a cycle"
    return vecs, stats


def orbit_reps(code: TS.BBCode, vecs: np.ndarray,
               perms: list[np.ndarray], chunk: int = 400_000) -> np.ndarray:
    """Translation-orbit representatives (chunked to bound peak memory)."""
    if len(vecs) == 0:
        return vecs
    keys_parts = [batch_keys(vecs[lo:lo + chunk], perms)
                  for lo in range(0, len(vecs), chunk)]
    keys = np.concatenate(keys_parts, axis=0)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return vecs[idx]


def lift_level(dt: DeckTables, cover: TS.BBCode, base: TS.BBCode,
               reps_below: np.ndarray, W: int, kerE_eq_cycles: bool,
               deadline: float, progress_tag: str = "") -> tuple[np.ndarray, dict]:
    """Build all cycles of weight <= W at the cover from the complete
    orbit-rep set one level below (fibers + tau branch)."""
    t0 = time.monotonic()
    stats: dict = {"fibers": 0, "fibers_empty": 0, "lifts": 0,
                   "tau_added": 0, "deep_fibers": 0}
    if not kerE_eq_cycles:
        raise EngineLimit("ker E != base cycle space "
                          "(tau branch beyond demonstrated lane)")
    out_rows: list[np.ndarray] = []
    # ---- fiber branch (beta != 0)
    wts = reps_below.sum(axis=1).astype(int) if len(reps_below) else \
        np.zeros(0, dtype=int)
    if len(reps_below):
        max_cap = int((W - wts.min()) // 2)
        if max_cap > CAP_ENVELOPE:
            raise EngineLimit(
                f"fiber cap {max_cap} > {CAP_ENVELOPE} "
                f"(lightest shadow |beta| = {int(wts.min())} at W = {W})")
    n_reps = len(reps_below)
    for ri in range(n_reps):
        if time.monotonic() > deadline:
            raise TimeoutError("budget exhausted during fiber stage")
        beta = reps_below[ri]
        wb = int(wts[ri])
        cap = (W - wb) // 2
        if cap < 0:
            continue
        if cap <= 4:
            lifts = enumerate_lifts(dt, beta, cap, kernel_cap=DEPS_GUARD)
        else:
            lifts = enumerate_lifts_deep(dt, beta, cap,
                                         kernel_cap=DEPS_GUARD + 2)
            stats["deep_fibers"] += 1
        stats["fibers"] += 1
        if not lifts:
            stats["fibers_empty"] += 1
        bmask = v2i(beta)
        for v0c, m2 in lifts.items():
            for v0i in (v0c, v0c ^ bmask):
                b = dt.deck.lift(i2v(v0i, base.n), beta)
                wl = int(b.sum())
                # slice identity: |b| = |beta| + 2 m2 (off-support part is
                # identical for both sheet choices)
                assert wl == wb + 2 * m2, "slice weight identity fails"
                if wl <= W:
                    assert cover.is_cycle(b), "lift is not a cover cycle"
                    out_rows.append(b)
                    stats["lifts"] += 1
        if stats["lifts"] > HIT_GUARD:
            raise Blowup(f"lift set > {HIT_GUARD:.0e}")
        if progress_tag and n_reps > 2000 and ri % 2000 == 0 and ri:
            log(f"      {progress_tag}: fiber {ri}/{n_reps} "
                f"({time.monotonic()-t0:.0f}s)")
    # ---- tau branch (beta = 0): gamma in cycle space, |gamma| <= W/2
    half = W // 2
    if len(reps_below):
        sel = reps_below[wts <= half]
        for g in sel:
            z = (dt.deck.TAU @ g) % 2
            assert int(z.sum()) == 2 * int(g.sum()), "tau weight identity"
            assert cover.is_cycle(z), "tau lift is not a cover cycle"
            out_rows.append(z)
            stats["tau_added"] += 1
    vecs = (np.array(out_rows, dtype=np.uint8) if out_rows
            else np.zeros((0, cover.n), dtype=np.uint8))
    stats["wall_s"] = round(time.monotonic() - t0, 2)
    return vecs, stats


def classify_level(code: TS.BBCode, reps: np.ndarray, max_check: int = 0
                   ) -> dict:
    """Weight histogram + (optionally bounded) stab/logical split."""
    out: dict = {}
    if len(reps) == 0:
        out["whist"] = {}
        return out
    wts = reps.sum(axis=1).astype(int)
    wh: dict[int, int] = {}
    for w in wts.tolist():
        wh[w] = wh.get(w, 0) + 1
    out["whist"] = {str(k): v for k, v in sorted(wh.items())}
    order = np.argsort(wts)
    n_log = n_stab = 0
    min_log_w = None
    min_stab_w = None
    checked = 0
    for i in order.tolist():
        if max_check and checked >= max_check:
            break
        v = reps[i]
        checked += 1
        if code.is_stab(v):
            n_stab += 1
            if min_stab_w is None:
                min_stab_w = int(wts[i])
        else:
            n_log += 1
            if min_log_w is None:
                min_log_w = int(wts[i])
    out.update({"checked": checked, "n_stab": n_stab, "n_logical": n_log,
                "min_stab_w": min_stab_w, "min_logical_w": min_log_w})
    return out


# ------------------------------------------------------------ the question
def close_question(spec: dict) -> dict:
    """Run one pre-registered question.  spec keys:
    ell, m, A, B, folds (frozen route), W, threads, workdir, budget_s,
    expected_levels (frozen k/kappa per level), planted_support (optional),
    tag."""
    t0 = time.monotonic()
    deadline = t0 + spec["budget_s"]
    W = int(spec["W"])
    workdir = Path(spec["workdir"])
    workdir.mkdir(parents=True, exist_ok=True)
    out: dict = {"tag": spec.get("tag"), "W": W, "stages": {},
                 "outcome": None}

    codes = build_codes(spec["ell"], spec["m"], spec["A"], spec["B"],
                        spec["folds"])
    J = len(codes) - 1
    out["levels_built"] = [
        {"level": i, "lm": list(c.G.orders), "n": c.n, "k": c.k,
         "kappa": len(c.rsHX_b)} for i, c in enumerate(codes)]
    # engine-frame consistency vs the frozen screen
    exp = spec.get("expected_levels")
    if exp:
        for got, want in zip(out["levels_built"], exp):
            assert got["k"] == want["k"] and got["kappa"] == want["kappa"], \
                f"engine frame mismatch at level {got['level']}: " \
                f"{got} vs frozen {want}"
    top_even = bool(all(int(kv.sum()) % 2 == 0 for kv in codes[0].kerHZ))
    out["top_cycles_all_even"] = top_even

    # decks + criterion-(i) measurements
    decks = []
    rungs = []
    for i, (axis, newmod) in enumerate(spec["folds"]):
        deck = make_deck(codes[i], codes[i + 1], axis, newmod)
        decks.append(deck)
        r = measure_rung(codes[i], codes[i + 1], deck, axis, newmod)
        r["rung"] = i
        # tau-branch scope: ker E vs base cycle space
        kerHZ_b = np.array(codes[i + 1].kerHZ, dtype=np.uint8)
        ok_contained = not ((deck.E @ kerHZ_b.T) % 2).any() \
            if len(kerHZ_b) else True
        rankE = TS.gf2_rank([v2i(row) for row in deck.E.T])
        dim_kerE = codes[i + 1].n - rankE
        r["kerE_eq_cycles"] = bool(ok_contained
                                   and dim_kerE == len(kerHZ_b))
        rungs.append(r)
    out["rungs"] = rungs

    # per-level plan: census species
    plan = []
    for i, c in enumerate(codes):
        dim = len(c.kerHZ)
        if i == J:
            species = "bz-census"
        elif dim <= FULL_ENUM_DIM_MAX:
            species = "full-enum"
        else:
            species = "fibers+tau"
        plan.append({"level": i, "dim_cycles": dim, "species": species,
                     "expected_hits": round(expected_hits(c.n, dim, W), 1)})
    out["plan"] = plan

    try:
        # pre-flight state-size guard over the whole plan (heuristic;
        # runtime guards re-check with realized counts)
        for p in plan:
            if p["expected_hits"] > HIT_GUARD:
                raise Blowup(
                    f"level {p['level']} expected cycle-set "
                    f"{p['expected_hits']:.2e} > {HIT_GUARD:.0e} "
                    f"(beyond the engine memory envelope at W = {W})")
        # ---------- bottom census (criterion-(ii) datum, always executed)
        bot = codes[J]
        if len(bot.kerHZ) <= FULL_ENUM_DIM_MAX:
            # run the BZ census anyway for (ii), then cross-check counts
            vecs, cstats = bottom_census(bot, W, spec["threads"], workdir,
                                         deadline)
            fe = full_enum_level(bot, W)
            assert len(fe) == len(vecs), \
                f"census {len(vecs)} != full-enum {len(fe)} at bottom"
            k1 = {bytes(x) for x in batch_keys(
                fe, TS._translation_perms(bot))} if len(fe) else set()
            k2 = {bytes(x) for x in batch_keys(
                vecs, TS._translation_perms(bot))} if len(vecs) else set()
            assert k1 == k2, "census/full-enum orbit sets differ"
            cstats["cross_check_full_enum"] = "EQUAL"
        else:
            vecs, cstats = bottom_census(bot, W, spec["threads"], workdir,
                                         deadline)
        out["stages"][f"L{J}_census"] = cstats
        perms = TS._translation_perms(bot)
        reps = orbit_reps(bot, vecs, perms)
        lv_cls = classify_level(bot, reps, max_check=400)
        lv_cls["set_size"] = int(len(vecs))
        lv_cls["orbit_reps"] = int(len(reps))
        out["stages"][f"L{J}_set"] = lv_cls
        del vecs

        # ---------- upward
        for i in range(J - 1, -1, -1):
            cover, base = codes[i], codes[i + 1]
            if time.monotonic() > deadline:
                raise TimeoutError("budget exhausted before rung "
                                   f"{i}")
            if plan[i]["species"] == "full-enum":
                fev = full_enum_level(cover, W)
                perms = TS._translation_perms(cover)
                reps = orbit_reps(cover, fev, perms)
                st = {"species": "full-enum", "set_size": int(len(fev))}
                del fev
            else:
                dt = DeckTables(decks[i])
                # covariance spot-check (translation-equivariance)
                if len(reps) > 1:
                    b0 = reps[int(len(reps) // 2)]
                    tsl = base.G.from_index(1)
                    perm = TS._perm_for(base, tsl)
                    l0 = enumerate_lifts(dt, b0, min(
                        2, (W - int(b0.sum())) // 2), kernel_cap=DEPS_GUARD)
                    l1 = enumerate_lifts(dt, b0[perm], min(
                        2, (W - int(b0.sum())) // 2), kernel_cap=DEPS_GUARD)
                    assert sorted(l0.values()) == sorted(l1.values()), \
                        "fiber m2-profile not translation-covariant"
                lv, st = lift_level(dt, cover, base, reps, W,
                                    rungs[i]["kerE_eq_cycles"], deadline,
                                    progress_tag=f"rung{i}")
                st["species"] = "fibers+tau"
                perms = TS._translation_perms(cover)
                reps = orbit_reps(cover, lv, perms)
                del lv
            cls = classify_level(
                cover, reps, max_check=(0 if i == 0 else 400))
            cls.update(st)
            cls["orbit_reps"] = int(len(reps))
            out["stages"][f"L{i}_set"] = cls
            log(f"    L{i}: {cls.get('species')} reps={len(reps)} "
                f"wall={st.get('wall_s', '-')}")

        # ---------- top verdict
        top = codes[0]
        cex = None
        if len(reps):
            wts = reps.sum(axis=1).astype(int)
            order = np.argsort(wts)
            for idx in order.tolist():
                v = reps[idx]
                if not top.is_stab(v):
                    assert top.is_cycle(v)
                    cex = {"weight": int(wts[idx]),
                           "support": [int(x) for x in np.nonzero(v)[0]],
                           "sig": [int(x) for x in top.sig(v)]}
                    break
        # planted-control check
        if spec.get("planted_support") is not None:
            pv = np.zeros(top.n, dtype=np.uint8)
            pv[spec["planted_support"]] = 1
            assert top.is_cycle(pv) and not top.is_stab(pv), \
                "planted vector is not a logical"
            if int(pv.sum()) <= W:
                pk = bytes(batch_keys(pv[None, :],
                                      TS._translation_perms(top))[0])
                keys = ({bytes(x) for x in batch_keys(
                    reps, TS._translation_perms(top))} if len(reps)
                    else set())
                out["planted_found"] = bool(pk in keys)
                assert out["planted_found"], \
                    "PLANTED LOGICAL NOT FOUND — engine completeness bug"
        if cex is None:
            floor = W + 1
            if top_even and floor % 2 == 1:
                floor += 1
            out["outcome"] = "CERTIFIED_FLOOR"
            out["floor"] = floor
            out["trust_tier"] = ("certificate (counting-invariant coset-BZ "
                                 "census + deck-transport lifting; not "
                                 "kernel-checked)")
        else:
            out["outcome"] = "COUNTEREXAMPLE"
            out["counterexample"] = cex
            out["note"] = (f"verified logical of weight {cex['weight']} "
                           f"<= W = {W}: the floor at this W is FALSE; "
                           f"d <= {cex['weight']}")
    except TimeoutError as e:
        out["outcome"] = "BUDGET_KILL"
        out["reason"] = str(e)
    except EngineLimit as e:
        out["outcome"] = "ENVELOPE_STOP"
        out["reason"] = str(e)
    except Blowup as e:
        out["outcome"] = "BLOWUP_STOP"
        out["reason"] = str(e)

    out["wall_s"] = round(time.monotonic() - t0, 2)
    return out
