"""A32 Part 5: GB-level (n = 90) censuses feeding the tower reduction.

Three census species, all via the a30 coset-BZ C kernel (two disjoint
information windows, r-pair complete to the stated weight):

  (a) LOGICAL census, weight <= 10, ALL 255 nonzero classes
      -> sector-B input: |gamma| in {8, 10} cycles with class outside
         im p_x*;  also re-derives d(GB) = 8 census-complete.
  (b) W-COSET census, weight <= 16, the 3 nonzero W-classes
      -> sector-A input (|beta| in {14, 16}); independently re-verifies
         A24 SS2.6's "W-coset minima = 14 exact" (SAT there, BZ here).
  (c) STABILIZER census, weight <= 16
      -> sector-C / flat-22 cost anchor (class counts per band, growth
         rate for the <= 22 extrapolation).

Completeness: two disjoint info windows I1, I2 of the 41-dim stabilizer
rowspace; a weight-w vector with w >= r1+1 on I1 and w >= r2+1 on I2 has
w >= r1+r2+2, so r-pair (r1, r2) is complete to W = r1+r2+1.  Node counts
are exact binomial sums (asserted by the kernel harness).

Output: data/a32/gb_census_{logical,wcoset,stab}.jsonl + summary json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a30_coset_bz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, run_window, unpack3,
)
from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402

DATA = LAB / "data" / "a32"
DATA.mkdir(parents=True, exist_ok=True)


def census(binp, tag, Gsys, I1, G1, I2, G2, bases, W, rpair, deadline):
    """Complete weight-<=W census of the given cosets; dedup across windows.

    The C kernel enumerates window patterns |S| >= 1 only, so the
    empty-window coset-base element of each window must be checked here
    (the A33 §8 erratum: this edge missed 15+2 vectors pre-fix)."""
    hits: dict[int, int] = {}  # vec_int -> base index
    nodes = 0
    for wi, (window, Gs, r) in enumerate(
            [(I1, G1, rpair[0]), (I2, G2, rpair[1])]):
        cb = [coset_base(Gs, window, b) for b in bases]
        for j, cbv in enumerate(cb):
            w = int(cbv.sum())
            if 0 < w <= W:
                hits[v2i(cbv)] = j
        res = run_window(binp, f"{tag}_w{wi}", Gs, cb, r, W, deadline)
        nodes += res["nodes"] * len(bases)
        for j, hx in res.pop("hit_rows"):
            v = unpack3(hx, Gsys.shape[1])
            hits[v2i(v)] = j
    return hits, nodes


def main():
    t0 = time.monotonic()
    out: dict = {}
    GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck_x = TS.Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                     lambda e, s: (e[0] + 15 * s, e[1]))
    perms = TS._translation_perms(GB)
    binp = build_kernel()
    I1, G1, I2, G2, kappa = disjoint_info_sets(GB.HX)
    assert kappa == 41, f"kappa {kappa}"
    print(f"[{time.monotonic()-t0:5.1f}s] GB frame: kappa=41, "
          f"windows {len(I1)}/{len(I2)} disjoint={not set(I1) & set(I2)}")

    # class-rep machinery: rep with prescribed sig
    S = np.array([GB.sig(r) for r in GB.xreps], dtype=np.uint8)  # 8x8
    SinvT = TS._gf2_inv(S.T)

    def rep_for(sig_int: int) -> np.ndarray:
        tvec = i2v(sig_int, 8)
        coeff = (SinvT @ tvec) % 2
        v = np.zeros(GB.n, dtype=np.uint8)
        for i in range(8):
            if coeff[i]:
                v ^= GB.xreps[i]
        assert v2i(GB.sig(v)) == sig_int
        return v

    # H1 maps needed for filters
    def h1_map(deck, tau=False):
        src = deck.cover if not tau else deck.base
        dst = deck.base if not tau else deck.cover
        Sm = np.array([src.sig(r) for r in src.xreps], dtype=np.uint8)
        op = (deck.P if not tau else deck.TAU)
        D = np.array([dst.sig((op @ r) % 2) for r in src.xreps],
                     dtype=np.uint8)
        return (D.T @ TS._gf2_inv(Sm.T)) % 2

    C = TS.BBCode("C", (30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26")
    deck_y = TS.Deck(C, BY, lambda e: (e[0], e[1] % 3),
                     lambda e, s: (e[0], e[1] + 3 * s))
    My = h1_map(deck_y)
    Mx = h1_map(deck_x)
    Ry = TS._colspace(My)
    Ryb, Ryp = rref_ints(list(Ry))
    imMx = TS._colspace(Mx)
    imb, imp = rref_ints(list(imMx))
    W_ints = sorted({TS._apply(Mx, s)
                     for s in TS._span_points(Ryb, Ryp)} - {0})
    Wb, Wp = rref_ints(list(W_ints))
    assert len(Wb) == 2 and len(W_ints) == 3

    # ---------------------------------------------------- (a) logical <= 10
    deadline = time.monotonic() + 3600
    bases = [rep_for(s) for s in range(1, 256)]
    hits, nodes = census(binp, "gb_log10", GB.HX, I1, G1, I2, G2,
                         bases, 10, (5, 4), deadline)
    # classify hits
    log_rows = []
    for vint in hits:
        v = i2v(vint, GB.n)
        w = int(v.sum())
        assert GB.is_cycle(v) and not GB.is_stab(v)
        sg = v2i(GB.sig(v))
        log_rows.append({"w": w, "sig": sg,
                         "in_im_px": TS.in_span(sg, imb, imp),
                         "in_W": TS.in_span(sg, Wb, Wp) and sg != 0,
                         "support": sorted(int(j) for j in np.nonzero(v)[0]),
                         "canon": TS._canon(v, perms)})
    wh = {}
    for r in log_rows:
        wh[r["w"]] = wh.get(r["w"], 0) + 1
    orb8 = {r["canon"] for r in log_rows if r["w"] == 8}
    orb10 = {r["canon"] for r in log_rows if r["w"] == 10}
    orb10_notim = {r["canon"] for r in log_rows
                   if r["w"] == 10 and not r["in_im_px"]}
    orb8_notim = {r["canon"] for r in log_rows
                  if r["w"] == 8 and not r["in_im_px"]}
    print(f"[{time.monotonic()-t0:5.1f}s] (a) GB logical census <= 10: "
          f"{len(log_rows)} vectors {wh}; orbits w8={len(orb8)} "
          f"w10={len(orb10)}; outside im p_x*: w8={len(orb8_notim)} "
          f"w10={len(orb10_notim)}  [{nodes:.2e} nodes]")
    with (DATA / "gb_census_logical.jsonl").open("w") as f:
        for r in sorted(log_rows, key=lambda r: (r["w"], r["canon"])):
            f.write(json.dumps(r) + "\n")
    out["logical"] = {"vectors": len(log_rows), "weight_hist": wh,
                      "orbits_w8": len(orb8), "orbits_w10": len(orb10),
                      "orbits_w8_notim": len(orb8_notim),
                      "orbits_w10_notim": len(orb10_notim), "nodes": nodes}

    # ---------------------------------------------------- (b) W-cosets <= 16
    bases = [rep_for(s) for s in W_ints]
    hits, nodes = census(binp, "gb_w16", GB.HX, I1, G1, I2, G2,
                         bases, 16, (8, 7), deadline)
    w_rows = []
    for vint in hits:
        v = i2v(vint, GB.n)
        w = int(v.sum())
        assert GB.is_cycle(v) and not GB.is_stab(v)
        sg = v2i(GB.sig(v))
        assert TS.in_span(sg, Wb, Wp) and sg != 0
        w_rows.append({"w": w, "sig": sg,
                       "support": sorted(int(j) for j in np.nonzero(v)[0]),
                       "canon": TS._canon(v, perms)})
    wh = {}
    for r in w_rows:
        wh[r["w"]] = wh.get(r["w"], 0) + 1
    min_w = min(wh) if wh else None
    orb14 = {r["canon"] for r in w_rows if r["w"] == 14}
    orb16 = {r["canon"] for r in w_rows if r["w"] == 16}
    assert min_w == 14, f"W-coset minimum {min_w} != 14 (A24 SS2.6)"
    print(f"[{time.monotonic()-t0:5.1f}s] (b) W-coset census <= 16: "
          f"{len(w_rows)} vectors {wh}; MINIMA = 14 REPRODUCED "
          f"(A24 SS2.6, SAT -> BZ); orbits w14={len(orb14)} "
          f"w16={len(orb16)}  [{nodes:.2e} nodes]")
    with (DATA / "gb_census_wcoset.jsonl").open("w") as f:
        for r in sorted(w_rows, key=lambda r: (r["w"], r["canon"])):
            f.write(json.dumps(r) + "\n")
    out["wcoset"] = {"vectors": len(w_rows), "weight_hist": wh,
                     "orbits_w14": len(orb14), "orbits_w16": len(orb16),
                     "nodes": nodes}

    # ---------------------------------------------------- (c) stabs <= 16
    hits, nodes = census(binp, "gb_stab16", GB.HX, I1, G1, I2, G2,
                         [np.zeros(GB.n, dtype=np.uint8)], 16, (8, 7),
                         deadline)
    s_rows = []
    for vint in hits:
        if vint == 0:
            continue
        v = i2v(vint, GB.n)
        w = int(v.sum())
        assert GB.is_stab(v)
        s_rows.append({"w": w, "canon": TS._canon(v, perms),
                       "support": sorted(int(j) for j in np.nonzero(v)[0])})
    orb_by_w: dict[int, set] = {}
    for r in s_rows:
        orb_by_w.setdefault(r["w"], set()).add(r["canon"])
    orbh = {w: len(s) for w, s in sorted(orb_by_w.items())}
    print(f"[{time.monotonic()-t0:5.1f}s] (c) GB stabilizer census <= 16: "
          f"{len(s_rows)} vectors; orbit histogram {orbh} "
          f"[{nodes:.2e} nodes]")
    with (DATA / "gb_census_stab.jsonl").open("w") as f:
        for r in sorted(s_rows, key=lambda r: (r["w"], r["canon"])):
            f.write(json.dumps(r) + "\n")
    out["stab"] = {"vectors": len(s_rows), "orbit_hist": orbh,
                   "nodes": nodes}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "gb_census_summary.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'gb_census_summary.json'}")


if __name__ == "__main__":
    main()
