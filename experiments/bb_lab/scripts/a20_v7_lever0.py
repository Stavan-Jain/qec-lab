"""A20 V7 Lever 0: pre-sweep measurements for the Y4 boundary-census
completeness problem (continuation of the a20-fibering session).

Three cheap measurements that decide the completeness-sweep design:

  M1 (the A22-P1 analog).  For each of the 1,655 SAT census classes
      (data/a20/m_census_classes.jsonl), min |f| and min site-support(f)
      over the 16 joint-kernel translates (f ranges over the coset
      f0 + ker d2, dim ker = 4).  A22 found every light delta-datum is
      A~m with |supp m| <= 3 sites; if the analog held here, completeness
      could reduce to a small-support-f sweep.  We also need the CONVERSE
      bound for that route to be complete -- so this is diagnostic only.

  M2 (kernel structure + the A23 gating fact).  Verify at Y4:
      ker(A*) = ker(B*) = ker d2 (each dim 4) -- the fact that makes the
      boundary set the graph code {(A*f, B*f)} with a single shared
      16-element redundancy.  Also weight/site-support histograms of the
      15 nonzero kernel elements (predicted: delta4-pure, weight 6 per
      active site).

  M3 (the A23 "Lever 5" feasibility profile).  Distribution of
      min(|u|, |v|) over the census classes, plus the flat-stratum cost
      C(72, b)/16 the A23 strata sweep would need at each min-side depth
      b <= 9.  (Total <= 18 even ==> min side <= 9.)

Usage: cd experiments/bb_lab && uv run python scripts/a20_v7_lever0.py
Writes data/a20/v7_lever0.json.
"""
import itertools
import json
import sys
from collections import Counter
from math import comb
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))
from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.group import AbelianGroup       # noqa: E402
from bb_lab.linalg import nullspace_f2      # noqa: E402
from bb_lab.poly import Poly                # noqa: E402

OUT = LAB / "data" / "a20"
CENSUS = OUT / "m_census_classes.jsonl"

Y4 = {"frame": (18, 4), "A": "1 + x + x^14*y", "B": "1 + x*y^2 + x^2*y^3"}


def main():
    G = AbelianGroup(Y4["frame"])
    checks = bb_check_matrices(Poly.from_string(Y4["A"], G),
                               Poly.from_string(Y4["B"], G))
    HX = checks.H_X % 2                    # 72 x 144, b = HX^T f
    D2 = HX.T.copy() % 2                   # 144 x 72
    DA, DB = D2[:72], D2[72:]              # u = DA f, v = DB f

    # ---- M2: kernel equality (the A23 gating fact) + kernel geometry
    KA = nullspace_f2(DA)
    KB = nullspace_f2(DB)
    K = nullspace_f2(D2)
    dims = (KA.shape[0], KB.shape[0], K.shape[0])
    # equality of the three spaces: check both one-sided containments by rank
    from bb_lab.linalg import rank_f2
    eqAB = rank_f2(np.vstack([KA, KB]) % 2) == KA.shape[0]
    eqAK = rank_f2(np.vstack([KA, K]) % 2) == KA.shape[0]
    print(f"M2: dim ker A* = {dims[0]}, ker B* = {dims[1]}, ker d2 = {dims[2]}; "
          f"ker A* = ker B*: {eqAB}; = ker d2: {eqAK}")

    # 16 kernel elements
    kernel = []
    for bits in itertools.product((0, 1), repeat=K.shape[0]):
        kv = np.zeros(72, dtype=np.uint8)
        for c, row in zip(bits, K):
            if c:
                kv ^= row.astype(np.uint8)
        kernel.append(kv)
    site_of = np.array([(a % 2) * 4 + b for a in range(18) for b in range(4)])

    def ssupp(fv):
        nz = np.flatnonzero(fv)
        return len(set(site_of[nz].tolist()))

    kw = sorted(int(kv.sum()) for kv in kernel if kv.any())
    ks = sorted(ssupp(kv) for kv in kernel if kv.any())
    print(f"    kernel element weights: {Counter(kw)}; site-supports: {Counter(ks)}")
    w_per_site = sorted(set(round(w / s, 2) for w, s in zip(kw, ks)))
    print(f"    weight per active site: {w_per_site} (delta4-pure prediction: 6.0)")

    # ---- M1: min |f| / min site-support over the kernel coset, per class
    rows = [json.loads(l) for l in CENSUS.read_text().splitlines()]
    classes = [r for r in rows if "w" in r]
    fmin_by_band = {}
    smin_by_band = {}
    per_class = []
    for r in classes:
        f0 = np.zeros(72, dtype=np.uint8)
        f0[r["f_support"]] = 1
        b = np.zeros(144, dtype=np.uint8)
        b[r["b_support"]] = 1
        assert ((D2 @ f0) % 2 == b).all(), "census row violates b = HX^T f"
        best_w, best_s = 99, 99
        for kv in kernel:
            fv = f0 ^ kv
            best_w = min(best_w, int(fv.sum()))
            best_s = min(best_s, ssupp(fv))
        per_class.append({"w": r["w"], "fmin": best_w, "smin": best_s})
        fmin_by_band.setdefault(r["w"], Counter())[best_w] += 1
        smin_by_band.setdefault(r["w"], Counter())[best_s] += 1
    print("\nM1: min |f| over kernel coset, by census band:")
    for w in sorted(fmin_by_band):
        print(f"    w={w}: {dict(sorted(fmin_by_band[w].items()))}")
    print("M1: min site-support(f) (8 sites max), by census band:")
    for w in sorted(smin_by_band):
        print(f"    w={w}: {dict(sorted(smin_by_band[w].items()))}")
    fmax = max(c["fmin"] for c in per_class)
    smax = max(c["smin"] for c in per_class)
    print(f"    => max over classes: min-|f| <= {fmax}, min-site-support <= {smax}")
    print(f"    sweep cost IF converse bound existed: C(72,{fmax}) = "
          f"{comb(72, fmax):,} raw supports (no converse bound is known -> "
          f"diagnostic only)")

    # ---- M3: min-side profile + flat-stratum cost (A23 Lever-5 feasibility)
    minside = Counter()
    for r in classes:
        b = np.zeros(144, dtype=np.uint8)
        b[r["b_support"]] = 1
        minside[min(int(b[:72].sum()), int(b[72:].sum()))] += 1
    print("\nM3: min(|u|,|v|) distribution over census classes:",
          dict(sorted(minside.items())))
    tot = 0
    print("    flat v-side stratum costs (A23-style, weight-b vectors in a "
          "codim-4 subspace of F2^72):")
    for bdepth in range(10):
        cost = comb(72, bdepth) // 16
        tot += cost
        print(f"      b={bdepth}: ~{cost:,}")
    print(f"    total to b=9: ~{tot:,}  (A23's f2a6 sweep was ~2x10^8 raw; "
          f"this is ~{tot / 2e8:,.0f}x that -- infeasible flat)")

    (OUT / "v7_lever0.json").write_text(json.dumps({
        "kernel_dims": dims, "kernel_eq": bool(eqAB and eqAK),
        "kernel_weights": kw, "kernel_site_supports": ks,
        "fmin_hist": {str(w): dict(c) for w, c in fmin_by_band.items()},
        "smin_hist": {str(w): dict(c) for w, c in smin_by_band.items()},
        "minside_hist": dict(minside)}))
    print(f"\nwrote {OUT / 'v7_lever0.json'}")


if __name__ == "__main__":
    main()
