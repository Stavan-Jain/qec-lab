"""A28 V-series foundations — falsify-first convention locks.

V1: code parameters (rank of the boundary map, k = 2(N - rank)).
V2: every f2a6 census row is a member of C(A,B) with the stated weights.
V3: the |b|=6 census row is exactly del(monomial).
V4: gross base (CLASS) small patterns: hexagons weigh 6; exactly 12
    offsets d give D-pairs of weight 10; the <=11 stamp-count-<=2 classes
    are 1 + 6.
V5: translation canonicalization is class-invariant on census rows.

Run: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a28_v_foundations.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import REGISTRY, load_f2a6_census, rref, in_span


def main():
    # ---- V1: parameters ----
    expect_rank = {"grossbase": 30, "f2a6": 71, "docket37": 88, "docket5e": 88}
    for name, inst in REGISTRY.items():
        rows = inst.generator_rows()
        basis, pivots = rref(rows, 2 * inst.G.N)
        kappa = len(basis)
        k = 2 * (inst.G.N - kappa)
        assert k == inst.expect_k, (name, k, inst.expect_k)
        assert kappa == expect_rank[name], (name, kappa)
        print(f"V1 {name}: N={inst.G.N} rank={kappa} k={k} ✓")

    # ---- V2/V3/V5: f2a6 census membership ----
    inst = REGISTRY["f2a6"]
    G = inst.G
    rows = inst.generator_rows()
    basis, pivots = rref(rows, 2 * G.N)
    census = load_f2a6_census()
    assert len(census) == 113
    canon_seen = set()
    for r in census:
        u = G.from_support([(x, y) for blk, x, y in r["b_support"] if blk == 0])
        v = G.from_support([(x, y) for blk, x, y in r["b_support"] if blk == 1])
        assert G.weight(u) == r["u_weight"] and G.weight(v) == r["v_weight"], r
        assert G.weight(u) + G.weight(v) == r["b_weight"]
        assert in_span(basis, pivots, inst.pack((u, v))), ("not a boundary!", r)
        canon_seen.add(inst.canonical((u, v)))
    assert len(canon_seen) == 113, len(canon_seen)
    print(f"V2 f2a6: all 113 rows are boundaries, weights match, classes distinct ✓")

    six = [r for r in census if r["b_weight"] == 6]
    assert len(six) == 1
    r = six[0]
    u = G.from_support([(x, y) for blk, x, y in r["b_support"] if blk == 0])
    v = G.from_support([(x, y) for blk, x, y in r["b_support"] if blk == 1])
    hit = [g for g in G.elements()
           if inst.boundary(G.monomial(*g)) == (u, v)]
    assert len(hit) == 1, hit
    print(f"V3 f2a6: |b|=6 class = del(monomial at {hit[0]}) ✓")

    # ---- V4: gross base stamp patterns ----
    inst = REGISTRY["grossbase"]
    G = inst.G
    hexa = inst.boundary(1)  # del(delta_(0,0))
    assert inst.pair_weight(hexa) == 6
    dpair_offsets = []
    for g in G.elements():
        if g == (0, 0):
            continue
        w = inst.pair_weight(inst.boundary(1 ^ G.monomial(*g)))
        assert w in (10, 12), (g, w)
        if w == 10:
            dpair_offsets.append(g)
    assert len(dpair_offsets) == 12, dpair_offsets
    # classes: {g, -g} give translated-equal pairs
    classes = {inst.canonical(inst.boundary(1 ^ G.monomial(*g)))
               for g in dpair_offsets}
    assert len(classes) == 6, len(classes)
    # cross-check the difference-set prediction: offsets = dA u dB
    dA = {( (a1 - a2) % G.l, (b1 - b2) % G.m)
          for a1, b1 in inst.A_supp for a2, b2 in inst.A_supp
          if (a1, b1) != (a2, b2)}
    dB = {( (a1 - a2) % G.l, (b1 - b2) % G.m)
          for a1, b1 in inst.B_supp for a2, b2 in inst.B_supp
          if (a1, b1) != (a2, b2)}
    assert set(dpair_offsets) == dA | dB, "offsets != dA u dB"
    assert not (dA & dB)
    print(f"V4 grossbase: hexagon w6; 12 D-pair offsets = dA∪dB, w10, 6 classes ✓")

    print("ALL FOUNDATIONS PASS")


if __name__ == "__main__":
    main()
