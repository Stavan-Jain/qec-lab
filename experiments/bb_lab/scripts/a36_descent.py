"""A36 Part 7: the descent lane — the gross-level censuses re-derived
from bb72 data (the A32 census-completeness bonus, third tower).

Every gross chain b with [b] in SEAM u {0} has beta = p_x(b) in
Stab(B72) u {0}  (SEAM = ker p_x*, the R3 identity), so the two census
species the direct closure consumed both decompose over the x-rung:

  beta != 0   bounded-overflow lift fibers over the bb72 stabilizer
              census <= 16 (orbit reps; caps (16-|beta|)/2 <= 5; shallow
              MITM lane to 4, deep ordered-split lane at 5 with the
              == shallow gate), every lift auto-classified stab-or-SEAM
              (p_x*[b] = 0 forces it — asserted per lift)
  beta  = 0   b = tau_x(gamma), |b| = 2|gamma|: the FULL bb72 cycle
              census <= 8 (single-window BZ over the dim-42 cycle space;
              complete since a weight-<=8 cycle has window weight <= 8),
              with the base-exactness prediction asserted per gamma:
              [gamma] in im p_x* <=> tau(gamma) is a gross stabilizer.

Comparison: G-canonical key sets (72 gross translations, sheet flip
absorbed) must equal the direct pass's — stab AND seam species.

Free bonus from the cycle census: d(B72) = 6 census-complete
(previously SAT-established), min nontrivial cycle weight + count.

Output: data/a36/descent.json
"""

# Provenance: copied verbatim 2026-08-18 (A38 S1) from the unmerged
# branch claude/tower-slice-calculus-generalize-410ed1 (the A35/A36
# session). That branch stays the source of truth until it merges;
# library-grade ports live in bb_lab.tower, not in edits here.


from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a30_coset_bz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, rref, run_window, unpack3,
)
from a30_rung_pass import i2v, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_sectorAC_full import batch_keys  # noqa: E402
from a32_subclosures import enumerate_lifts  # noqa: E402
from a32_deep_fibers import enumerate_lifts_deep  # noqa: E402
from a33_tower_cells import h1_map  # noqa: E402
from a36_tower_cells import build_tower, seam_info  # noqa: E402

DATA = LAB / "data" / "a36"
W = 16


def main():
    t0 = time.monotonic()
    out: dict = {}
    G8, GR, B72, B36, deck_y, deck_x, deck_b = build_tower()
    My = h1_map(deck_y)
    sd = seam_info(GR, My)
    perms = TS._translation_perms(GR)
    perms72 = TS._translation_perms(B72)
    binp = build_kernel()

    # ------------------------------------------- direct key sets (target)
    def load_keys(fn, field):
        vecs = []
        for line in (DATA / fn).open():
            r = json.loads(line)
            v = np.zeros(GR.n, dtype=np.uint8)
            v[r[field]] = 1
            vecs.append(v)
        return {bytes(k) for k in
                batch_keys(np.array(vecs, dtype=np.uint8), perms)}

    k_stab_direct = load_keys("stab_census_orbits.jsonl", "b_support")
    k_seam_direct = load_keys("seam_census.jsonl", "w_support")
    print(f"[{time.monotonic()-t0:6.1f}s] direct targets: "
          f"{len(k_stab_direct)} stab keys, {len(k_seam_direct)} seam keys")

    # --------------------------------------- bb72 stabilizer census <= 16
    I1, G1, I2, G2, kappa = disjoint_info_sets(B72.HX)
    assert kappa == 30
    hits: set[int] = set()
    nodes = 0
    for wi, (window, Gs, r) in enumerate([(I1, G1, 8), (I2, G2, 7)]):
        cb = coset_base(Gs, window, np.zeros(B72.n, dtype=np.uint8))
        assert not cb.any()
        res = run_window(binp, f"a36_b72stab_w{wi}", Gs, [cb], r, W,
                         time.monotonic() + 600)
        nodes += res["nodes"]
        for _, hx in res.pop("hit_rows"):
            v = unpack3(hx, B72.n)
            if v.any():
                hits.add(v2i(v))
    exp = sum(math.comb(30, s) for s in range(1, 9)) + \
        sum(math.comb(30, s) for s in range(1, 8))
    assert nodes == exp
    bvecs = np.array([i2v(h, B72.n) for h in sorted(hits)], dtype=np.uint8)
    bws = bvecs.sum(axis=1)
    for v in bvecs[:: max(1, len(bvecs) // 40)]:
        assert B72.is_stab(v)
    whist72 = {int(w): int((bws == w).sum()) for w in sorted(set(bws))}
    keys72 = batch_keys(bvecs, perms72)
    orb72: dict[bytes, int] = {}
    for i, k in enumerate(keys72):
        orb72.setdefault(bytes(k), i)
    print(f"[{time.monotonic()-t0:6.1f}s] bb72 stab census <= 16: "
          f"{len(bvecs)} vectors {whist72}, {len(orb72)} orbits "
          f"[{nodes:.2e} nodes]")
    out["b72_stab_census"] = {"vectors": int(len(bvecs)),
                              "weight_hist": whist72,
                              "orbits": len(orb72), "nodes": nodes}

    # ------------------------------------------- bb72 cycle census <= 8
    K = np.array(B72.kerHZ, dtype=np.uint8)
    assert K.shape[0] == 42
    Rn, piv = rref(K)
    assert len(piv) == 42
    res = run_window(binp, "a36_b72cyc", Rn, [np.zeros(B72.n,
                     dtype=np.uint8)], 8, 8, time.monotonic() + 600)
    assert res["nodes"] == sum(math.comb(42, s) for s in range(1, 9))
    gammas = []
    for _, hx in res.pop("hit_rows"):
        g = unpack3(hx, B72.n)
        if g.any():
            gammas.append(g)
    gws = np.array([int(g.sum()) for g in gammas])
    ghist = {int(w): int((gws == w).sum()) for w in sorted(set(gws))}
    # bonus: d(B72) census-complete
    ntl = [(int(w), g) for w, g in zip(gws, gammas) if not B72.is_stab(g)]
    d72 = min(w for w, _ in ntl)
    n_d72 = sum(1 for w, _ in ntl if w == d72)
    assert d72 == 6, f"d(B72) = {d72} != 6"
    print(f"[{time.monotonic()-t0:6.1f}s] bb72 cycle census <= 8: "
          f"{len(gammas)} cycles {ghist}; d(B72) = 6 CENSUS-COMPLETE "
          f"({n_d72} weight-6 nontrivial logicals) [solver-free "
          f"re-derivation; {res['nodes']:.2e} nodes]")
    out["b72_cycle_census"] = {"cycles": len(gammas), "weight_hist": ghist,
                               "d_B72": d72, "n_min": n_d72,
                               "nodes": res["nodes"]}

    # --------------------------- beta = 0 family: b = tau(gamma), <= 16
    imPx, imPxp = None, None
    Mx = h1_map(deck_x)
    imPx, imPxp = (lambda b: (b, [(x & -x).bit_length() - 1 for x in b]))(
        TS._colspace(Mx))[0], None
    from a30_rung_pass import rref_ints  # noqa: E402
    imb, imp = rref_ints(list(TS._colspace(Mx)))
    tau_stab, tau_seam = [], []
    n_exact_checks = 0
    for g in gammas:
        b = (deck_x.TAU @ g) % 2
        assert int(b.sum()) == 2 * int(g.sum()), "tau did not double weight"
        is_st = GR.is_stab(b)
        in_im = TS.in_span(v2i(B72.sig(g)), imb, imp)
        assert is_st == in_im, \
            "base exactness violated: [gamma] in im p* <=> tau(g) stab"
        n_exact_checks += 1
        if is_st:
            tau_stab.append(b)
        else:
            assert TS.in_span(v2i(GR.sig(b)), sd["basis"], sd["piv"]), \
                "nonzero tau-lift class outside SEAM (cover exactness?!)"
            tau_seam.append(b)
    print(f"[{time.monotonic()-t0:6.1f}s] tau-family: {len(tau_stab)} "
          f"stab + {len(tau_seam)} seam lifts; base-exactness prediction "
          f"asserted {n_exact_checks}/{n_exact_checks}")

    # ------------------------- beta != 0 fibers over bb72 stab orbit reps
    # deep-lane gate: == shallow lane at cap <= 4 on the first 8 fibers
    n_gate = 0
    for key, i in list(orb72.items())[:8]:
        beta = bvecs[i]
        cap = min(4, (W - int(bws[i])) // 2)
        if cap < 0:
            continue
        a = enumerate_lifts(deck_x, beta, cap)
        b = enumerate_lifts_deep(deck_x, beta, cap)
        assert a == b, "deep lane != shallow lane at equal cap"
        n_gate += 1
    print(f"[{time.monotonic()-t0:6.1f}s] deep-lane gate: {n_gate}/8 "
          f"fibers identical shallow vs deep at cap <= 4")

    fib_stab, fib_seam = [], []
    n_fibers = n_lifts = n_empty = 0
    tF = time.monotonic()
    for key, i in orb72.items():
        beta = bvecs[i]
        wb = int(bws[i])
        cap = (W - wb) // 2
        if cap < 0:
            continue
        lifts = (enumerate_lifts(deck_x, beta, cap) if cap <= 4
                 else enumerate_lifts_deep(deck_x, beta, cap))
        n_fibers += 1
        if not lifts:
            n_empty += 1
        for v0_int, m2 in lifts.items():
            b = deck_x.lift(i2v(v0_int, B72.n), beta)
            assert int(b.sum()) == wb + 2 * m2 <= W
            n_lifts += 1
            if GR.is_stab(b):
                fib_stab.append(b)
            else:
                assert TS.in_span(v2i(GR.sig(b)), sd["basis"], sd["piv"]), \
                    "fiber lift class outside SEAM u 0 ?!"
                fib_seam.append(b)
    print(f"[{time.monotonic()-t0:6.1f}s] fibers: {n_fibers} over bb72 "
          f"stab orbit reps, {n_lifts} lifts ({len(fib_stab)} stab / "
          f"{len(fib_seam)} seam), {n_empty} empty "
          f"({n_empty/n_fibers:.0%}) [{time.monotonic()-tF:.1f}s]")
    out["fibers"] = {"fibers": n_fibers, "lifts": n_lifts,
                     "empty": n_empty,
                     "empty_rate": round(n_empty / n_fibers, 3)}

    # ------------------------------------------------- key-set equality
    def keyset(vlist):
        if not vlist:
            return set()
        return {bytes(k) for k in
                batch_keys(np.array(vlist, dtype=np.uint8), perms)}

    k_stab_desc = keyset(tau_stab) | keyset(fib_stab)
    k_seam_desc = keyset(tau_seam) | keyset(fib_seam)
    assert k_stab_desc == k_stab_direct, \
        (len(k_stab_desc), len(k_stab_direct),
         len(k_stab_desc ^ k_stab_direct))
    assert k_seam_desc == k_seam_direct, \
        (len(k_seam_desc), len(k_seam_direct),
         len(k_seam_desc ^ k_seam_direct))
    print(f"[{time.monotonic()-t0:6.1f}s] KEY-SET EQUALITY: stab "
          f"{len(k_stab_desc)} == direct {len(k_stab_direct)}; seam "
          f"{len(k_seam_desc)} == direct {len(k_seam_direct)}  "
          f"[census completeness re-derived from n = 72 data]")
    out["equality"] = {"stab_keys": len(k_stab_desc),
                       "seam_keys": len(k_seam_desc),
                       "stab_equal": True, "seam_equal": True}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "descent.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'descent.json'}")


if __name__ == "__main__":
    main()
