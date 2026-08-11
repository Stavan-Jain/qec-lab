"""A32 Part 7: full sector-A closure + shallow sector-C fibers.

Extends a32_subclosures.py to the remaining engine-reachable fibers:

  A18/A20/A22  W-coset census to weight 22 (r-pair (11,10), complete),
               orbit reps at |beta| in {18, 20, 22}, lift caps (2, 1, 0),
               top rungs at M = ((24 - |b|) + 1) // 2.
               With 6b/6c this CLOSES SECTOR A: every SF24-y violation
               candidate with [beta] in W \\ 0 is excluded.

  C14/C16      GB-stabilizer fibers at |beta| in {14, 16} (caps 4, 3):
               lifts with [b] != 0 are x-dangerous BY-logicals in K_x
               classes -> top rungs; lifts with [b] = 0 and |b| <= 20
               must appear in the banked A19 (M)@24 census (bands <= 20
               COMPLETE) -> completeness cross-check, no rung needed
               (floors banked); lifts with [b] = 0, |b| = 22 are
               flat-22-dangerous candidates -> flat-top rung at M = 1.

  C18/C20/C22  GB-stabilizer census to weight 22 (base 0), fibers at
               |beta| in {18, 20, 22} (caps 2, 1, 0), same trichotomy.

Residue deliberately left open (deep fibers): sector C at
|beta| in {6, 10, 12} (caps 8, 6, 5 — beyond the size-4 MITM lane);
cost-modelled in the note.

Output: data/a32/sectorAC_summary.json + per-lane jsonl.
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
    build_kernel, disjoint_info_sets,
)
from a30_rung_pass import i2v, rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_gb_census import census  # noqa: E402
from a32_subclosures import enumerate_lifts  # noqa: E402
from scope_bravyi_rung import BravyiRungCell  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a32"


# ----------------------------------------------------- batched orbit canon
def batch_keys(vecs: np.ndarray, perms: list[np.ndarray]) -> np.ndarray:
    """Lexicographic-min packed key over translation perms, per row."""
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


def orbit_reps(vecs: np.ndarray, perms: list[np.ndarray]
               ) -> tuple[np.ndarray, np.ndarray]:
    """(reps, keys): one representative row per translation orbit."""
    keys = batch_keys(vecs, perms)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return vecs[idx], keys[idx]


def main():
    t0 = time.monotonic()
    out: dict = {}
    GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    C = TS.BBCode("C", (30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26")
    deck_x = TS.Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                     lambda e, s: (e[0] + 15 * s, e[1]))
    deck_y = TS.Deck(C, BY, lambda e: (e[0], e[1] % 3),
                     lambda e, s: (e[0], e[1] + 3 * s))
    perms_gb = TS._translation_perms(GB)
    perms_by = TS._translation_perms(BY)
    cell = BravyiRungCell()
    binp = build_kernel()
    I1, G1, I2, G2, kappa = disjoint_info_sets(GB.HX)
    print(f"[{time.monotonic()-t0:5.1f}s] frames built (kappa={kappa})")

    # W-class reps
    S = np.array([GB.sig(r) for r in GB.xreps], dtype=np.uint8)
    SinvT = TS._gf2_inv(S.T)

    def rep_for(sig_int: int) -> np.ndarray:
        tvec = i2v(sig_int, 8)
        coeff = (SinvT @ tvec) % 2
        v = np.zeros(GB.n, dtype=np.uint8)
        for i in range(8):
            if coeff[i]:
                v ^= GB.xreps[i]
        return v

    def h1_map(deck, tau=False):
        src = deck.cover if not tau else deck.base
        dst = deck.base if not tau else deck.cover
        Sm = np.array([src.sig(r) for r in src.xreps], dtype=np.uint8)
        op = (deck.P if not tau else deck.TAU)
        D = np.array([dst.sig((op @ r) % 2) for r in src.xreps],
                     dtype=np.uint8)
        return (D.T @ TS._gf2_inv(Sm.T)) % 2

    My = h1_map(deck_y)
    Mx = h1_map(deck_x)
    Ry = TS._colspace(My)
    Ryb, Ryp = rref_ints(list(Ry))
    W_ints = sorted({TS._apply(Mx, s)
                     for s in TS._span_points(Ryb, Ryp)} - {0})

    # banked A19 census canonical keys (for the [b]=0 completeness check)
    banked_by_w: dict[int, set] = {}
    banked_vecs, banked_ws = [], []
    for line in (MAIN / "a19" / "m24_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            b = np.zeros(BY.n, dtype=np.uint8)
            b[r["b_support"]] = 1
            banked_vecs.append(b)
            banked_ws.append(r["w"])
    BV = np.array(banked_vecs, dtype=np.uint8)
    BK = batch_keys(BV, perms_by)
    for i, w in enumerate(banked_ws):
        banked_by_w.setdefault(w, set()).add(bytes(BK[i]))
    print(f"[{time.monotonic()-t0:5.1f}s] banked census canonicalized "
          f"({len(banked_vecs)} classes)")

    def rung_lane(tag, fibers, caps_by_w, sector_c=False):
        """Run lift fibers + rungs; returns summary dict."""
        rows = []
        stats = {"fibers": 0, "lifts": 0, "rungs": 0, "violations": 0,
                 "flat22_dangerous_rungs": 0, "banked_checked": 0,
                 "lift_m2_hist": {}, "verdicts": {}}
        for rep_i, beta in enumerate(fibers):
            wbeta = int(beta.sum())
            cap = caps_by_w[wbeta]
            lifts = enumerate_lifts(deck_x, beta, cap=cap)
            stats["fibers"] += 1
            vh: dict[str, int] = {}
            m2h: dict[int, int] = {}
            for v0c, m2 in sorted(lifts.items()):
                b = deck_x.lift(i2v(v0c, GB.n), beta)
                wb = int(b.sum())
                assert wb == wbeta + 2 * m2 and BY.is_cycle(b)
                m2h[m2] = m2h.get(m2, 0) + 1
                stats["lifts"] += 1
                k = str(m2)
                stats["lift_m2_hist"][k] = stats["lift_m2_hist"].get(k, 0) + 1
                if sector_c and BY.is_stab(b):
                    if wb <= 20:
                        # must be in the banked complete census
                        key = bytes(batch_keys(b[None, :], perms_by)[0])
                        assert key in banked_by_w.get(wb, set()), \
                            f"stab lift w{wb} missing from banked census!"
                        stats["banked_checked"] += 1
                        continue
                    # wb == 22: flat-22-dangerous candidate
                    res = cell.rung(b, 1, time.monotonic() + 600)
                    assert res["verdict"] != "ABORT", res
                    stats["flat22_dangerous_rungs"] += 1
                    stats["rungs"] += 1
                    vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
                    stats["verdicts"][res["verdict"]] = \
                        stats["verdicts"].get(res["verdict"], 0) + 1
                    if res["verdict"] == "VIOLATION":
                        stats["violations"] += 1
                        print(f"  !! VIOLATION ({tag} flat22):",
                              json.dumps(res)[:300])
                    continue
                if sector_c:
                    assert not BY.is_stab(b)
                M = (24 - wb + 1) // 2
                res = cell.rung(b, M, time.monotonic() + 600)
                assert res["verdict"] != "ABORT", res
                stats["rungs"] += 1
                vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
                stats["verdicts"][res["verdict"]] = \
                    stats["verdicts"].get(res["verdict"], 0) + 1
                if res["verdict"] == "VIOLATION":
                    stats["violations"] += 1
                    print(f"  !! VIOLATION ({tag}):", json.dumps(res)[:300])
            rows.append({"wbeta": wbeta, "m2_hist": m2h, "verdicts": vh})
            if rep_i % 500 == 0 and rep_i:
                print(f"    {tag}: {rep_i} fibers "
                      f"({time.monotonic()-t0:.0f}s)")
        with (DATA / f"sectorAC_{tag}.jsonl").open("w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
        return stats

    # --------------------------------------- sector A: census to 22
    deadline = time.monotonic() + 7200
    bases = [rep_for(s) for s in W_ints]
    hits, nodes = census(binp, "gb_w22", GB.HX, I1, G1, I2, G2,
                         bases, 22, (11, 10), deadline)
    print(f"[{time.monotonic()-t0:5.1f}s] W-coset census <= 22: "
          f"{len(hits)} vectors [{nodes:.2e} nodes]")
    hv = np.array([i2v(h, GB.n) for h in hits], dtype=np.uint8)
    ws = hv.sum(axis=1)
    orb_hist = {}
    fibersA = []
    for w in (14, 16, 18, 20, 22):
        sel = hv[ws == w]
        if len(sel) == 0:
            orb_hist[w] = 0
            continue
        reps, _ = orbit_reps(sel, perms_gb)
        orb_hist[w] = len(reps)
        if w >= 18:
            fibersA.extend(list(reps))
    assert orb_hist[14] == 6 and orb_hist[16] == 68, \
        f"W-census orbit mismatch vs Part 5: {orb_hist}"
    print(f"    orbit histogram {orb_hist}  (14/16 match Part 5)")
    out["wcensus22"] = {"vectors": len(hits), "orbits": orb_hist,
                        "nodes": nodes}
    capsA = {18: 2, 20: 1, 22: 0}
    statsA = rung_lane("A18to22", fibersA, capsA, sector_c=False)
    print(f"[{time.monotonic()-t0:5.1f}s] sector A 18-22: {statsA}")
    out["sectorA_18to22"] = statsA

    # --------------------------------------- sector C: stab census to 22
    hits, nodes = census(binp, "gb_stab22", GB.HX, I1, G1, I2, G2,
                         [np.zeros(GB.n, dtype=np.uint8)], 22, (11, 10),
                         deadline)
    print(f"[{time.monotonic()-t0:5.1f}s] GB stab census <= 22: "
          f"{len(hits)} vectors [{nodes:.2e} nodes]")
    hv = np.array([i2v(h, GB.n) for h in hits if h != 0], dtype=np.uint8)
    ws = hv.sum(axis=1)
    orb_hist_s = {}
    fibers_by_w: dict[int, list] = {}
    for w in sorted(set(int(x) for x in ws)):
        sel = hv[ws == w]
        reps, _ = orbit_reps(sel, perms_gb)
        orb_hist_s[w] = len(reps)
        fibers_by_w[w] = list(reps)
    assert orb_hist_s[6] == 1 and orb_hist_s[10] == 6 \
        and orb_hist_s[12] == 21 and orb_hist_s[14] == 64 \
        and orb_hist_s[16] == 333, f"stab orbits {orb_hist_s} vs Part 5"
    print(f"    stab orbit histogram {orb_hist_s}")
    out["stabcensus22"] = {"vectors": int(len(hv)), "orbits": orb_hist_s,
                           "nodes": nodes}

    capsC = {14: 4, 16: 3, 18: 2, 20: 1, 22: 0}
    statsC1 = rung_lane("C14to16",
                        fibers_by_w[14] + fibers_by_w[16], capsC,
                        sector_c=True)
    print(f"[{time.monotonic()-t0:5.1f}s] sector C 14-16: {statsC1}")
    out["sectorC_14to16"] = statsC1
    statsC2 = rung_lane("C18to22",
                        fibers_by_w[18] + fibers_by_w[20] + fibers_by_w[22],
                        capsC, sector_c=True)
    print(f"[{time.monotonic()-t0:5.1f}s] sector C 18-22: {statsC2}")
    out["sectorC_18to22"] = statsC2

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "sectorAC_summary.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'sectorAC_summary.json'}")


if __name__ == "__main__":
    main()
