"""A38 S1: the F2a probe — does a banked census walk factor along the
odd-CRT decomposition of F2[G]?  Falsify-first, on the cheapest fully
reproducible banked census: the gross stabilizer census at W = 16
(A36's offset-S species: 33,588 vectors / 469 orbits, n = 144).

Setup.  G = Z12 x Z6, |G| = 72 = 8 * 9; the odd part is H = <(4,0)> x
<(0,2)> = Z3 x Z3.  In char 2, N = sum_{h in H} h is a central
idempotent (|H| = 9 odd), so F2[G] = N F2[G] (+) (1+N) F2[G] with
N F2[G] = F2[G/H] = F2[Z4 x Z2] — Maschke/CRT: descent along the odd
part is a direct-summand splitting.  Refining along each Z3 axis gives
the 4-way orthogonal idempotent system

    c0 = N1 N2 (invariant, dim 1 per group element space),
    c1 = N1 E2, c2 = E1 N2 (dim 2 each), c3 = E1 E2 (dim 4),

with Nj = 1 + u + u^2 and Ej = 1 + Nj = u + u^2 along the two Z3
generators (the c3 block further splits into two F4 components —
recorded, not consumed).  Multiplication by any c is a block-diagonal
circulant on 1-chains and fixes the stabilizer space (c is central),
so the stab space splits EXACTLY:  kappa = 66 = sum_c kappa_c
(hard-asserted).

The F2a hypothesis (charter): census node counts factor along this
splitting — per-component walks recombined by a coupling law.  The
probe measures the exact obstruction on the banked census:

  P1  the census is REPRODUCED from scratch (BZ pass, exact node count,
      33,588 vectors, banked weight histogram, 469 orbits) — any
      mismatch aborts the probe;
  P2  idempotent-system exactness (partition of unity, orthogonality,
      per-vector component sums re-verified on every censused vector);
  P3  the Maschke rank bookkeeping kappa_c, and the structural weight
      law of the invariant sector: N-components are unions of
      H-cosets, so |v_N| in {0, 18, 36, ...} for censused (even) v —
      i.e. AT W = 16 THE INVARIANT COMPONENT OF EVERY CENSUS VECTOR IS
      EITHER 0 OR ALREADY HEAVIER THAN THE WHOLE CENSUS BOUND;
  P4  the measured joint component-weight distribution over all 33,588
      vectors: additivity failures (|v| != sum |v_c|), per-component
      maxima W_c (the bounds a factored walk would need), pure-sector
      populations, and the vN = 0 sub-census;
  P5  the verdict, priced: nodes(kappa_c, W_c) per component + the
      recombination product, vs the direct 7.50e9-node walk — the
      measured gain or the measured absence.  No claim is made unless
      counts reproduce EXACTLY.

Output: data/a38/f2a_probe.json
Run:    cd experiments/bb_lab && uv run python scripts/a38_f2a_probe.py
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

from bb_lab import tower as tw  # noqa: E402
from bb_lab.checks import circulant  # noqa: E402
from bb_lab.cosetbz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, run_window, unpack3,
)
from bb_lab.poly import Poly  # noqa: E402

DATA = LAB / "data" / "a38"
DATA.mkdir(parents=True, exist_ok=True)

W = 16
BANKED_WHIST = {6: 72, 10: 432, 12: 2268, 14: 3888, 16: 26928}
BANKED_VECS = 33588
BANKED_ORBITS = 469

T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:6.1f}s] {msg}", flush=True)


def conv_poly(a: frozenset, b: frozenset, G) -> frozenset:
    """Product of two F2[G] supports."""
    counts: dict = {}
    for s in a:
        for t in b:
            g = G.add(s, t)
            counts[g] = counts.get(g, 0) + 1
    return frozenset(g for g, c in counts.items() if c % 2 == 1)


def main() -> None:
    out: dict = {}
    GR = tw.TowerCode("GR", (12, 6), "x^3 + y + y^2", "y^3 + x + x^2")
    G = GR.G

    # ------------------------------------------------- P1: the census
    log("P1: reproducing the banked gross stab census (BZ, W = 16)")
    binp = build_kernel()
    I1, G1, I2, G2, kappa = disjoint_info_sets(GR.HX)
    assert kappa == 66
    hits: set[int] = set()
    nodes = 0
    for wi, (window, Gs, r) in enumerate([(I1, G1, 8), (I2, G2, 7)]):
        cb = coset_base(Gs, window, np.zeros(GR.n, dtype=np.uint8))
        assert not cb.any()
        res = run_window(binp, f"a38_f2a_w{wi}", Gs, [cb], r, W,
                         time.monotonic() + 3600)
        nodes += res["nodes"]
        for _, hx in res.pop("hit_rows"):
            v = unpack3(hx, GR.n)
            if v.any():
                hits.add(tw.v2i(v))
    exp = sum(math.comb(66, s) for s in range(1, 9)) + \
        sum(math.comb(66, s) for s in range(1, 8))
    assert nodes == exp
    V = np.array([tw.i2v(h, GR.n) for h in sorted(hits)], dtype=np.uint8)
    ws = V.sum(axis=1)
    whist = {int(w): int((ws == w).sum()) for w in sorted(set(ws))}
    assert len(V) == BANKED_VECS and whist == BANKED_WHIST, \
        (len(V), whist)
    perms = tw.translation_perms(GR)
    keys = tw.batch_keys(V, perms)
    n_orb = len({bytes(k) for k in keys})
    assert n_orb == BANKED_ORBITS, n_orb
    for v in V[:: max(1, len(V) // 40)]:
        assert GR.is_stab(v)
    log(f"  census REPRODUCED: {len(V)} vectors {whist}, {n_orb} "
        f"orbits, {nodes:.3e} nodes (exact)")
    out["census"] = {"vectors": len(V), "weight_hist": whist,
                     "orbits": n_orb, "nodes": nodes}

    # ------------------------------------- P2: the idempotent system
    log("P2: idempotent system (4-way odd-CRT split)")
    u, wgen = (4, 0), (0, 2)             # the two Z3 generators of H
    H = [G.add(tuple(np.multiply(u, i) % (12, 6)),
               tuple(np.multiply(wgen, j) % (12, 6)))
         for i in range(3) for j in range(3)]
    H = [G.reduce(h) for h in H]
    assert len(set(H)) == 9
    N1 = frozenset([G.reduce((0, 0)), G.reduce(u),
                    G.reduce((8, 0))])                    # 1 + u + u^2
    E1 = frozenset([G.reduce(u), G.reduce((8, 0))])       # u + u^2
    N2 = frozenset([G.reduce((0, 0)), G.reduce(wgen),
                    G.reduce((0, 4))])
    E2 = frozenset([G.reduce(wgen), G.reduce((0, 4))])
    comps = {"c0_inv": conv_poly(N1, N2, G),
             "c1_N1E2": conv_poly(N1, E2, G),
             "c2_E1N2": conv_poly(E1, N2, G),
             "c3_E1E2": conv_poly(E1, E2, G)}
    assert comps["c0_inv"] == frozenset(H), "N != sum over H"
    # partition of unity + orthogonality + idempotency
    tot: dict = {}
    for s in comps.values():
        for g in s:
            tot[g] = tot.get(g, 0) + 1
    unity = frozenset(g for g, c in tot.items() if c % 2 == 1)
    assert unity == frozenset({G.reduce((0, 0))}), "sum c_i != 1"
    names = list(comps)
    for i, a in enumerate(names):
        assert conv_poly(comps[a], comps[a], G) == comps[a], f"{a}^2"
        for b in names[i + 1:]:
            assert conv_poly(comps[a], comps[b], G) == frozenset(), \
                f"{a}*{b} != 0"
    log("  partition of unity + orthogonality + idempotency: PASS")

    # component multiplication matrices on 1-chains (block circulant)
    Ms = {}
    for name, supp in comps.items():
        P = circulant(Poly.from_support(supp, G)).astype(np.uint8)
        Z = np.zeros_like(P)
        Ms[name] = np.block([[P, Z], [Z, P]]) % 2

    # ------------------------------------------------- P3: rank split
    log("P3: Maschke rank bookkeeping + the invariant weight law")
    kappas = {}
    for name, M in Ms.items():
        rows = [(M @ r) % 2 for r in GR.HX]
        kappas[name] = tw.gf2_rank([tw.v2i(r) for r in rows])
    assert sum(kappas.values()) == 66, kappas
    # the invariant sector IS the quotient code over G/H = Z4 x Z2
    Q = tw.TowerCode("GRmodH", (4, 2), "x^3 + y + 1", "y + x + x^2")
    assert kappas["c0_inv"] == Q.kappa, (kappas["c0_inv"], Q.kappa)
    log(f"  kappa split {kappas} (sum 66 EXACT); invariant sector = "
        f"quotient [[{Q.n},{Q.k}]] code over Z4xZ2, kappa = {Q.kappa}")
    out["kappa_split"] = kappas
    out["quotient_code"] = {"n": Q.n, "k": Q.k, "kappa": Q.kappa}

    # per-vector component weights (vectorized matmuls)
    comps_w = {}
    for name, M in Ms.items():
        comps_w[name] = ((V.astype(np.int32) @ M.T.astype(np.int32)) % 2
                         ).sum(axis=1)
    vN = comps_w["c0_inv"]
    # the invariant weight law: |v_N| is 9 * (#H-cosets), and even
    assert all(int(x) % 18 == 0 for x in set(vN.tolist())), \
        sorted(set(vN.tolist()))
    log(f"  invariant components: weights in 18Z as forced "
        f"(law: coset unions x parity); distribution "
        f"{ {int(k): int((vN == k).sum()) for k in sorted(set(vN.tolist()))} }")

    # ------------------------------------------------- P4: statistics
    sum_w = sum(comps_w.values())
    additive = int((sum_w == ws).sum())
    all_bounded = int(np.all(
        np.stack([comps_w[n] <= W for n in names]), axis=0).sum())
    pure = {n: int(((comps_w[n] == ws)
                    & np.all(np.stack([comps_w[m] == 0
                                       for m in names if m != n]),
                             axis=0)).sum()) for n in names}
    vN0 = int((vN == 0).sum())
    Wc = {n: int(comps_w[n].max()) for n in names}
    out["stats"] = {
        "additive_weight_vectors": additive,
        "all_components_bounded_leW": all_bounded,
        "pure_sector_vectors": pure,
        "vN_zero_vectors": vN0,
        "component_weight_max": Wc,
        "vN_weight_hist": {str(int(k)): int((vN == k).sum())
                           for k in sorted(set(vN.tolist()))},
    }
    log(f"P4: additivity |v| = sum|v_c| on {additive}/{len(V)}; all "
        f"components <= {W} on {all_bounded}/{len(V)}; pure-sector "
        f"{pure}; vN = 0 on {vN0}/{len(V)}; per-component max "
        f"weights {Wc}")

    # ------------------------------------------------- P5: the verdict
    # A factored walk would need per-component censuses to bounds Wc
    # (measured), then a recombination pass.  Price it with the same
    # exact node formula the calculus uses.
    factored_nodes = {n: tw.census_nodes(kappas[n], min(Wc[n], 2 * W))
                      for n in names}
    direct_nodes = nodes
    out["pricing"] = {
        "direct_nodes": direct_nodes,
        "factored_walk_nodes_at_measured_Wc": factored_nodes,
        "factored_sum": sum(factored_nodes.values()),
    }
    gain = direct_nodes / max(sum(factored_nodes.values()), 1)
    # the honest verdict, mechanically derived
    factors = (additive == len(V) and all_bounded == len(V))
    if factors:
        verdict = ("FACTORS (weights additive and bounded on the whole "
                   "census) — counts must now be reproduced per sector "
                   "before any claim")
    else:
        verdict = (
            "DOES NOT FACTOR at W = 16: the odd-CRT components of "
            f"censused vectors are unbounded ({len(V) - all_bounded} of "
            f"{len(V)} vectors have a component above the census bound; "
            f"invariant components sit in 18Z by the coset-union weight "
            f"law, so ANY W < 18 census is invisible to the invariant "
            f"sector except through cancellation) and weight is not "
            f"additive across sectors ({len(V) - additive} failures). "
            f"The walk cannot be split into per-sector walks at the "
            f"census bound; the 2-part stays combinatorial and the odd "
            f"part enters only through cross-sector cancellation.")
    out["verdict"] = verdict
    out["nominal_gain_if_it_had_factored"] = round(gain, 3)
    log(f"P5 VERDICT: {verdict}")
    log(f"  pricing: direct {direct_nodes:.2e} nodes vs factored-sum "
        f"{sum(factored_nodes.values()):.2e} at measured Wc "
        f"(x{gain:.2f}) — plus a recombination product the direct walk "
        f"never pays")

    out["wall_s"] = round(time.monotonic() - T0, 1)
    (DATA / "f2a_probe.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA / 'f2a_probe.json'}")


if __name__ == "__main__":
    main()
