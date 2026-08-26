#!/usr/bin/env python3
"""A41 Q6 — census self-similarity probe (S1 gate; banked data ONLY).

Question (A41 charter Q6): along the banked towers, does the level-r
cycle census renormalize to the level-(r-1) census under any simple
weight reindexing w -> a*w + b (a in {1,2}, small even b)?  A positive
signal is the precondition for a *vertical* closed-form family being
descent-provable (the census-carrying recursion becoming inductive);
no signal on every known tower => park the front, numbers recorded.

No new enumeration: reads the banked A38 c37xx freeze census, the A36
bb288 direct-close census, and the A32 GB censuses.  Verdict per tower
pair + a cross-tower note.  Output: data/a41/q6_selfsim_probe.json.
"""

from __future__ import annotations

import json
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = DATA / "a41" / "q6_selfsim_probe.json"


def load(p):
    p = DATA / p
    if not p.exists():
        return None
    with open(p) as fh:
        return json.load(fh)


def jsonl_whist(p, weight_key="w", count_orbits=True):
    """Weight histogram from a banked census jsonl (orbit rows)."""
    p = DATA / p
    if not p.exists():
        return None
    hist: dict[int, int] = {}
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            w = row.get(weight_key)
            if w is None:
                for k in ("weight", "wt", "|w|"):
                    if k in row:
                        w = row[k]
                        break
            if w is None:
                return None
            hist[int(w)] = hist.get(int(w), 0) + 1
    return dict(sorted(hist.items()))


def compare(hi: dict, lo: dict, wmax_hi: int, wmax_lo: int):
    """Best affine reindex w_hi = a*w_lo + b over a in {1,2}, even b in
    [-6,6]; score = spread of count ratios over the overlapping bands
    (1.0 = exactly proportional).  Windows are truncated so only bands
    complete on BOTH sides enter (censuses are complete only to their
    own W)."""
    hi = {int(k): v for k, v in hi.items()}
    lo = {int(k): v for k, v in lo.items()}
    best = None
    for a in (1, 2):
        for b in range(-6, 7, 2):
            pairs = []
            for wl, cl in lo.items():
                wh = a * wl + b
                if wh in hi and wh <= wmax_hi and wl <= wmax_lo:
                    pairs.append((wl, wh, cl, hi[wh]))
            if len(pairs) < 3:
                continue
            ratios = [ch / cl for (_, _, cl, ch) in pairs]
            spread = max(ratios) / min(ratios)
            cand = {"a": a, "b": b, "bands": len(pairs),
                    "ratio_min": round(min(ratios), 3),
                    "ratio_max": round(max(ratios), 3),
                    "spread": round(spread, 3),
                    "pairs": [(wl, wh, cl, ch) for wl, wh, cl, ch in pairs]}
            if best is None or cand["spread"] < best["spread"]:
                best = cand
    return best


def main() -> None:
    report = {"towers": {}, "cross_tower": {}, "verdict": None}

    # ---- c37xx (A38 S2): L3 (15,3) n=90 / L2 (15,6) n=180 / L1 (30,6)
    # n=360, stab censuses to W=22.  L3 hist is VECTORS; L2/L1 are
    # ORBITS -- compare like with like only (L2 vs L1, orbit-grade).
    c = load("a38/c37xx/freeze_W22_census.json")
    if c:
        l3v = {int(k): v for k, v in c["L3_stab"]["weight_hist"].items()}
        l2o = {int(k): v for k, v in c["L2"]["stab_whist"].items()}
        l1o = {int(k): v for k, v in c["L1_stab"]["whist"].items()}
        report["towers"]["c37xx"] = {
            "levels": {"L3(n=90)_vectors": l3v, "L2(n=180)_orbits": l2o,
                       "L1(n=360)_orbits": l1o},
            "mu_sequence": [c["L3_stab"]["mu"], c["L2"]["mu2"],
                            c["L1_stab"]["mu1"]],
            "d_sequence_known": [c["d_L3"], 10, 20, ">=24"],
            "L1_vs_L2_best_reindex": compare(l1o, l2o, 22, 22),
            "L2_vs_L3_note": "unit mismatch (orbits vs vectors); "
                             "shape-only comparison recorded",
            "L2_vs_L3_vector_shape": compare(
                l3v, l3v, 22, 22) and None,  # placeholder, see below
        }
        del report["towers"]["c37xx"]["L2_vs_L3_vector_shape"]

    # ---- a36 bb288 tower: mid-level (6,12) [[144,12]] stab census to
    # W=16, both vector- and orbit-grade; seam census to 16.
    a36 = load("a36/direct_close_banked.json")
    if a36:
        report["towers"]["bb288_mid"] = {
            "stab_orbit_hist": a36["stab_census"]["orbit_weight_hist"],
            "stab_vector_hist": a36["stab_census"]["weight_hist"],
            "seam_hist": a36["seam_census"]["weight_hist"],
        }

    # ---- a32 GB level (n=90) censuses, if present in banked form.
    for name, path in [("a32_gb_stab", "a32/gb_census_stab.jsonl"),
                       ("a32_gb_logical", "a32/gb_census_logical.jsonl")]:
        h = jsonl_whist(path)
        if h:
            report["towers"][name] = {"whist": h}

    # ---- cross-tower orbit-count prefix comparison (different codes,
    # same low bands): bb288 mid vs c37xx L2 -- a regularity check, not
    # self-similarity; flagged for what it is.
    if c and a36:
        l2o = {int(k): v for k, v in c["L2"]["stab_whist"].items()}
        m = {int(k): v for k, v in
             a36["stab_census"]["orbit_weight_hist"].items()}
        shared = sorted(set(l2o) & set(m))
        report["cross_tower"]["bb288mid_vs_c37xxL2_orbit_prefix"] = {
            str(w): [m[w], l2o[w]] for w in shared if w <= 16}

    # ---- verdict: does ANY tower pair renormalize (spread <= 1.5 over
    # >= 4 bands)?
    signals = []
    for tname, t in report["towers"].items():
        cmp_ = t.get("L1_vs_L2_best_reindex")
        if cmp_ and cmp_["spread"] <= 1.5 and cmp_["bands"] >= 4:
            signals.append((tname, cmp_))
    report["verdict"] = {
        "signal": bool(signals),
        "detail": signals or "no affine weight-reindex makes any banked "
                             "level pair proportional (spread > 1.5 or "
                             "< 4 comparable bands)",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=1, default=str)

    # ---- console
    for tname, t in report["towers"].items():
        print(f"== {tname}")
        for k, v in t.items():
            if k == "L1_vs_L2_best_reindex" and v:
                print(f"  best reindex w_L1 = {v['a']}*w_L2 + {v['b']}: "
                      f"{v['bands']} bands, ratio spread {v['spread']} "
                      f"[{v['ratio_min']}, {v['ratio_max']}]")
                for wl, wh, cl, ch in v["pairs"]:
                    print(f"    L2 w={wl}: {cl}  ->  L1 w={wh}: {ch} "
                          f"(x{ch / cl:.2f})")
            elif k in ("mu_sequence", "d_sequence_known"):
                print(f"  {k}: {v}")
    if report["cross_tower"]:
        print("== cross-tower orbit prefix (bb288 mid vs c37xx L2, w<=16)")
        for w, (a, b) in sorted(report["cross_tower"]
                                ["bb288mid_vs_c37xxL2_orbit_prefix"].items(),
                                key=lambda kv: int(kv[0])):
            print(f"    w={w}: {a} vs {b}")
    print("VERDICT:", "SIGNAL" if report["verdict"]["signal"] else
          "NO SIGNAL", "-", report["verdict"]["detail"] if not
          report["verdict"]["signal"] else report["verdict"]["detail"])
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
