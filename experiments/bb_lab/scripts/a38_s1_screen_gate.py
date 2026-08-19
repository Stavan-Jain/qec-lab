"""A38 S1: the F5 library gate — bb_lab.tower must reproduce the A35
generality screen before anything new runs.

Three layers, all falsify-first:

  [G1] bb_lab.tower.validate_banked(): node-count anchors (x1.00 / x1.01
       / the x3.00 A33 shared-walk datum), bit-level sector-C fiber-layer
       reproduction (397 fibers / 4,132 lifts / m2-hist), banked A32/A33
       tower structure, A19 deck-survey k-verdicts, no-deck rows.
  [G2] the full 11-tower A35 docket re-screened THROUGH THE LIBRARY with
       the banked RNG stream (seed 20260811, same call order), then
       deep-compared field-by-field against the banked
       data/a35/screen_banked.json (from branch
       claude/tower-slice-calculus-generalize-410ed1) — including the
       sampled fiber sections, which only match if the rank-generic code
       consumes the RNG exactly like the frozen rank-2 screen did.
       Only wall_s is exempt.
  [G3] the screen's demo sections re-run via the library: the (3,6)
       bottom full enumeration (d = 4 exact, stab weight hist, 27
       min-weight classes / 9 inside im p*) and the |A|-even parity demo
       (weight-7 cycle), compared against the banked values; plus a
       rank-3 smoke test (a k = 4 pair on Z4xZ3xZ3 with a twisted
       k-preserving Z2 x-fold, full rung screen) — a TYPE-level check
       that the group core is rank-generic, not a claim about any
       trivariate code.

Output: data/a38/screen_lib.json + pass/fail per layer on stdout.
Run:    cd experiments/bb_lab && uv run python scripts/a38_s1_screen_gate.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab import tower as tw  # noqa: E402

DATA = LAB / "data" / "a38"
DATA.mkdir(parents=True, exist_ok=True)
BANKED = LAB / "data" / "a35" / "screen_banked.json"

T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:6.1f}s] {msg}", flush=True)


def deep_diff(a, b, path="") -> list[str]:
    """All leaf differences between two JSON-ish values (wall_s exempt)."""
    diffs: list[str] = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k == "wall_s":
                continue
            if k not in a:
                diffs.append(f"{path}.{k}: missing in lib")
            elif k not in b:
                diffs.append(f"{path}.{k}: missing in banked")
            else:
                diffs += deep_diff(a[k], b[k], f"{path}.{k}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diffs.append(f"{path}: len {len(a)} != {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            diffs += deep_diff(x, y, f"{path}[{i}]")
    else:
        if a != b:
            diffs.append(f"{path}: {a!r} != {b!r}")
    return diffs


def main() -> None:
    rng = np.random.default_rng(20260811)   # the banked screen's seed

    # ---------------------------------------------------------- [G1]
    log("[G1] validate_banked (anchors + sector-C refiber + structure)")
    val = tw.validate_banked(LAB / "data", rng=rng, log=log)
    log("[G1] PASS")

    # ---------------------------------------------------------- [G2]
    log("[G2] full 11-tower docket through the library (shared stream)")
    towers = [val[n] for n in ("bravyi360", "ibm288Y", "gross_xx",
                               "bb288_yxx")]
    for spec in tw.A35_DOCKET[4:]:
        towers.append(tw.screen_tower(spec, rng=rng, log=log))

    banked = json.loads(BANKED.read_text())
    banked_towers = banked["towers"]
    assert [t["name"] for t in banked_towers] == \
        [t["name"] for t in towers], "docket order/name mismatch"
    n_diffs = 0
    for lib_t, bank_t in zip(towers, banked_towers):
        d = deep_diff(lib_t, bank_t, bank_t["name"])
        for line in d:
            log(f"    DIFF {line}")
        n_diffs += len(d)
    # anchors: same three ratios (banked stores them under 'anchors')
    for tag in ("a32_wcoset22", "a32_stab22", "a33_h5"):
        assert val["anchors"][tag] == banked["anchors"][tag], tag
    # sector-C refiber numbers vs the banked screen's record
    bk = banked["anchors"]["sectorC_1416_refiber"]
    lv = val["sectorC_refiber"]
    assert (lv["fibers"], lv["lifts"], lv["m2_hist"]) == \
        (bk["fibers"], bk["lifts"], bk["m2_hist"]), "sector-C refiber"
    assert val["sectorC_heavy_empty_rates_banked"] == \
        banked["anchors"]["sectorC_heavy_empty_rates_banked"]
    assert n_diffs == 0, f"{n_diffs} field diffs vs banked screen"
    log("[G2] PASS: 11/11 towers field-identical to the banked screen "
        "(fibers included — RNG stream reproduced)")

    # ---------------------------------------------------------- [G3]
    log("[G3] demo sections + rank-3 smoke")
    # (3,6) bottom full enumeration (the A35 [1b] demo)
    bot = tw.TowerCode("btm36", (3, 6), "1 + y + y^2", "y^3 + x + x^2")
    bb72 = tw.TowerCode("bb72", (6, 6), "x^3 + y + y^2", "y^3 + x + x^2")
    deck_b = tw.AxisDeck(bb72, bot, 0)
    imP_b, _ = tw.rref_ints(list(tw.colspace(tw.h1_map(deck_b))))
    stab_ints = list(tw.span_points(bot.rsHX_b))
    assert len(stab_ints) == 1 << 14, "the (3,6) stab group is 2^14"
    swh: dict[int, int] = {}
    for x in stab_ints:
        w = bin(x).count("1")
        if 0 < w <= 16:
            swh[w] = swh.get(w, 0) + 1
    K = np.array(bot.kerHZ, dtype=np.uint8)
    dimk = K.shape[0]
    assert dimk == 22
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
    piv = [(b & -b).bit_length() - 1 for b in imP_b]
    outside = [s for s in attain if not tw.in_span(s, imP_b, piv)]
    demo = {"d_exact": int(d_bot),
            "stab_whist_le16": {str(k): v for k, v in sorted(swh.items())},
            "min_classes": len(attain),
            "min_classes_outside_imp": len(outside),
            "per_class_minima_complete": True}
    bk_demo = banked["bb288_bottom_demo"]
    bk_demo = {**bk_demo,
               "stab_whist_le16": {str(k): v for k, v in
                                   bk_demo["stab_whist_le16"].items()}}
    assert demo == bk_demo, (demo, bk_demo)
    log(f"    (3,6) demo: d = {d_bot} exact, stab hist + 27/9 class "
        f"split == banked")

    # parity demo: |A| even -> an odd-weight cycle exists
    pd = None
    for glm, As, Bs in [((6, 6), "1 + x", "1 + y"),
                        ((6, 6), "1 + x + y + x*y", "y^3 + x + x^2"),
                        ((4, 4), "1 + x", "1 + y")]:
        code = tw.TowerCode("paritydemo", glm, As, Bs)
        odd = [kv for kv in code.kerHZ if int(kv.sum()) % 2 == 1]
        if odd and code.k > 0:
            pd = {"lm": list(glm), "A": As, "B": Bs, "k": code.k,
                  "odd_cycle_weight": int(odd[0].sum())}
            break
    assert pd == banked["parity_demo"], (pd, banked["parity_demo"])
    log(f"    parity demo == banked: k={pd['k']}, odd cycle weight "
        f"{pd['odd_cycle_weight']}")

    # rank-3 smoke: the group core is rank-generic (TYPE-level check
    # only, no distance claims).  A k = 4 trinomial pair on Z4xZ3xZ3
    # whose x-fold is a TWISTED k-preserving Z2 deck, so the full rung
    # screen (H1 rank law, exactness, sigma*, codim_lift) runs at rank 3.
    c3 = tw.TowerCode("r3cover", (4, 3, 3),
                      [(0, 1, 1), (0, 2, 2), (1, 2, 1)],
                      [(0, 2, 1), (3, 0, 1), (3, 1, 1)])
    b3 = tw.TowerCode("r3base", (2, 3, 3),
                      tw.fold_support(c3.A.support, 0, 2),
                      tw.fold_support(c3.B.support, 0, 2))
    assert (c3.k, b3.k) == (4, 4), (c3.k, b3.k)
    inv3 = tw.tower_inventory((4, 3, 3))
    assert inv3["v2_per_axis"] == [2, 0, 0] and inv3["odd_part"] == 9
    r3 = tw.screen_rung(c3, b3, 0, W_eff=6,
                        rng=np.random.default_rng(0))
    assert r3["twisted"] and r3["R_holds"]
    assert r3["rank_p"] == r3["rank_tau"] == 2 == c3.k // 2, r3
    assert r3["exact_cover"], "im tau* != ker p* at rank 3"
    assert r3["codim_lift"] == b3.k - c3.k // 2, r3
    d3 = tw.AxisDeck(c3, b3, 0)
    rng3 = np.random.default_rng(0)
    for _ in range(20):
        v = c3.random_cycle(rng3)
        b, m, _ = d3.slice_data(v)
        assert b3.is_cycle(b)
    log(f"    rank-3 smoke: [[{c3.n},{c3.k}]] -> [[{b3.n},{b3.k}]] "
        f"twisted (R) rung; rank law + cover exactness + codim_lift + "
        f"20 slice checks PASS at rank 3 (type-level only — no claims)")

    # ------------------------------------------------------------- bank
    out = {"anchors": val["anchors"],
           "sectorC_refiber": val["sectorC_refiber"],
           "sectorC_heavy_empty_rates_banked":
               val["sectorC_heavy_empty_rates_banked"],
           "towers": towers, "no_deck": val["no_deck"],
           "bb288_bottom_demo": demo, "parity_demo": pd,
           "rank3_smoke": {"cover": [c3.n, c3.k], "base": [b3.n, b3.k]},
           "gate": "PASS (G1+G2+G3)",
           "wall_s": round(time.monotonic() - T0, 1)}
    (DATA / "screen_lib.json").write_text(json.dumps(out, indent=1))
    log(f"GATE PASS total {out['wall_s']}s -> {DATA / 'screen_lib.json'}")


if __name__ == "__main__":
    main()
