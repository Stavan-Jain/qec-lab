"""A32 Part 8: the deep sector-C fibers |beta| in {6, 10, 12} — the last
open cells of the SF24-y + flat-22 program.

Caps m2 <= (22-|beta|)/2 = 8 / 6 / 5 — beyond the size-4 MITM lane of
`a32_subclosures.enumerate_lifts`.  This script extends the lane to
size <= 8 via ordered-split meet-in-the-middle (every sorted off-support
set X of size s splits uniquely into low floor(s/2) and high ceil(s/2)
parts; join size-a masks with size-b buckets demanding lsb(high) >
msb(low)), then dispatches each lift b by the trichotomy:

  [b] != 0 (x-dangerous BY-logical, class in K_x, reachable):
      |b| <= 16  -> must already be known: bands <= 14 have NO reachable
                    members (A24 band censuses) and band 16 has exactly
                    6, all shallow-decomposing (Part 4a) — a deep-fiber
                    hit here would be a FALSIFICATION event (assert);
      |b| in {18, 20, 22} -> top rung at M = 3 / 2 / 1 (restricted lanes).
  [b] = 0 (BY-stabilizer lift):
      |b| <= 20  -> class must appear in the banked A19 (M)@24 census
                    (bands <= 20 complete; floors banked) — asserted;
      |b| = 22   -> flat-22-dangerous candidate: flat-top rung at M = 1.

Validation gate (falsify-first): at caps <= 4 the deep enumerator must
reproduce `enumerate_lifts` EXACTLY (asserted on |beta| = 14 and 16
fibers before any deep fiber runs).

ALL PASS closes sector C completely, and with §4-§5 the entire y-safe
sector + flat-22 residue.

Output: data/a32/deep_fibers.jsonl + summary.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a30_rung_pass import i2v, rref_ints, reduce_int, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_subclosures import enumerate_lifts  # noqa: E402
from a32_sectorAC_full import batch_keys  # noqa: E402
from scope_bravyi_rung import BravyiRungCell  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a32"


def enumerate_lifts_deep(deck: TS.Deck, beta: np.ndarray, cap: int,
                         kernel_cap: int = 20) -> dict[int, int]:
    """All v0 with E v0 = RHS beta, |v0 off supp(beta)| <= cap (cap <= 8).

    Ordered-split MITM: X (sorted, size s) = low part (size s//2) +
    high part (size (s+1)//2), lsb(high) > msb(low).  Complete by the
    exact-off-support subset-sum argument (see enumerate_lifts).
    """
    assert cap <= 8
    n = deck.base.n
    E_cols = [v2i(deck.E[:, j]) for j in range(n)]
    rhs = (deck.RHS @ beta) % 2
    rhs_i = v2i(rhs)
    bsupp = [int(j) for j in np.nonzero(beta)[0]]
    bmask = v2i(beta)
    bcols = [E_cols[j] for j in bsupp]
    bb, bp = rref_ints(bcols)
    rhs_res = reduce_int(rhs_i, bb, bp)
    offb = [j for j in range(n) if not (bmask >> j) & 1]
    red = {j: reduce_int(E_cols[j], bb, bp) for j in offb}
    half = (cap + 1) // 2
    # subsets by (size, reduced sum) -> list of index-masks
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
                        if a == s:  # X = high part only
                            hits_X.add(amask)
                        continue
                    if bmask2.bit_length() - 1 < alsb \
                            and not (amask & bmask2):
                        hits_X.add(amask | bmask2)
    # per-X kernel enumeration (as in enumerate_lifts)
    out: dict[int, int] = {}
    for X in sorted(hits_X):
        cols = bsupp + [j for j in range(n) if (X >> j) & 1]
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
        assert len(deps) <= kernel_cap, f"kernel 2^{len(deps)} at X={X:x}"
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
            assert not (((deck.E @ v0) + rhs) % 2).any()
            canon = min(v0_int, v0_int ^ bmask)
            prev = out.get(canon)
            if prev is None or m2 < prev:
                out[canon] = m2
    return out


def main():
    t0 = time.monotonic()
    out: dict = {}
    GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck_x = TS.Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                     lambda e, s: (e[0] + 15 * s, e[1]))
    perms_gb = TS._translation_perms(GB)
    perms_by = TS._translation_perms(BY)
    cell = BravyiRungCell()
    print(f"[{time.monotonic()-t0:5.1f}s] frames built")

    # ------------------------------------------------ validation gate
    stab_recs = [json.loads(x) for x in
                 (DATA / "gb_census_stab.jsonl").open()]
    for wtest in (14, 16):
        reps = {}
        for r in stab_recs:
            if r["w"] == wtest:
                reps.setdefault(r["canon"], r)
        r = sorted(reps.values(), key=lambda r: r["canon"])[0]
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[r["support"]] = 1
        cap = (22 - wtest) // 2
        l1 = enumerate_lifts(deck_x, beta, cap=min(cap, 4))
        l2 = enumerate_lifts_deep(deck_x, beta, cap=min(cap, 4))
        assert l1 == l2, f"deep enumerator mismatch at w{wtest}: " \
            f"{len(l1)} vs {len(l2)}"
    print(f"[{time.monotonic()-t0:5.1f}s] validation gate: deep enumerator "
          f"== size-4 lane on |beta|=14,16 fibers  [EXACT]")

    # banked census canonical keys + A24 band-16 reachable keys
    banked_by_w: dict[int, set] = {}
    bvs, bws = [], []
    for line in (MAIN / "a19" / "m24_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            b = np.zeros(BY.n, dtype=np.uint8)
            b[r["b_support"]] = 1
            bvs.append(b)
            bws.append(r["w"])
    BK = batch_keys(np.array(bvs, dtype=np.uint8), perms_by)
    for i, w in enumerate(bws):
        banked_by_w.setdefault(w, set()).add(bytes(BK[i]))
    reach16 = set()
    for x in (MAIN / "a24" /
              "cell_census_reach_band16_mchecks.jsonl").open():
        e = json.loads(x)
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        reach16.add(bytes(batch_keys(b[None, :], perms_by)[0]))
    print(f"[{time.monotonic()-t0:5.1f}s] banked keys ready")

    # ------------------------------------------------ the 28 deep fibers
    fibers = []
    for r in stab_recs:
        if r["w"] in (6, 10, 12):
            fibers.append(r)
    reps = {}
    for r in fibers:
        reps.setdefault(r["canon"], r)
    assert len(reps) == 28, f"{len(reps)} deep fibers != 28"
    rows = []
    stats = {"fibers": 0, "lifts": 0, "rungs": 0, "violations": 0,
             "flat22_rungs": 0, "banked_checked": 0,
             "light_logical_hits": 0, "verdicts": {},
             "lift_hist": {}}
    for rep_i, r in enumerate(sorted(reps.values(),
                                     key=lambda r: (r["w"], r["canon"]))):
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[r["support"]] = 1
        wbeta = r["w"]
        cap = (22 - wbeta) // 2
        tF = time.monotonic()
        lifts = enumerate_lifts_deep(deck_x, beta, cap=cap)
        vh: dict[str, int] = {}
        wbh: dict[str, int] = {}
        for v0c, m2 in sorted(lifts.items()):
            b = deck_x.lift(i2v(v0c, GB.n), beta)
            wb = int(b.sum())
            assert wb == wbeta + 2 * m2 and BY.is_cycle(b)
            stats["lifts"] += 1
            k = f"{wb}"
            wbh[k] = wbh.get(k, 0) + 1
            stats["lift_hist"][k] = stats["lift_hist"].get(k, 0) + 1
            if BY.is_stab(b):
                if wb <= 20:
                    key = bytes(batch_keys(b[None, :], perms_by)[0])
                    assert key in banked_by_w.get(wb, set()), \
                        f"stab lift w{wb} missing from banked census!"
                    stats["banked_checked"] += 1
                    continue
                res = cell.rung(b, 1, time.monotonic() + 600)
                assert res["verdict"] != "ABORT", res
                stats["flat22_rungs"] += 1
                stats["rungs"] += 1
                vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
                stats["verdicts"][res["verdict"]] = \
                    stats["verdicts"].get(res["verdict"], 0) + 1
                if res["verdict"] == "VIOLATION":
                    stats["violations"] += 1
                    print("  !! VIOLATION (deep flat22):",
                          json.dumps(res)[:300])
                continue
            # logical lift
            if wb <= 16:
                stats["light_logical_hits"] += 1
                assert wb == 16, \
                    f"NEW reachable logical at w{wb} — FALSIFIES A24 bands!"
                key = bytes(batch_keys(b[None, :], perms_by)[0])
                assert key in reach16, \
                    "NEW band-16 reachable class — FALSIFIES A24 census!"
                continue  # closed by A24 (m >= 5) + Part 4a
            M = (24 - wb + 1) // 2
            res = cell.rung(b, M, time.monotonic() + 600)
            assert res["verdict"] != "ABORT", res
            stats["rungs"] += 1
            vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
            stats["verdicts"][res["verdict"]] = \
                stats["verdicts"].get(res["verdict"], 0) + 1
            if res["verdict"] == "VIOLATION":
                stats["violations"] += 1
                print("  !! VIOLATION (deep logical):",
                      json.dumps(res)[:300])
        stats["fibers"] += 1
        rows.append({"wbeta": wbeta, "cap": cap, "lift_whist": wbh,
                     "verdicts": vh,
                     "wall_s": round(time.monotonic() - tF, 1)})
        print(f"    deep fiber {rep_i+1}/28 |beta|={wbeta}: "
              f"lifts_by_wb {wbh} verdicts {vh} "
              f"({rows[-1]['wall_s']}s)")
    print(f"[{time.monotonic()-t0:5.1f}s] deep fibers: {stats}")
    out["deep"] = stats
    with (DATA / "deep_fibers.jsonl").open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "deep_fibers_summary.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'deep_fibers_summary.json'}")


if __name__ == "__main__":
    main()
