"""The stratified twist-floor sweep for bb_288 — dual deck, all λ.

Per (deck, λ ≠ 0): certify F2(λ) via banded enumeration + masked-K
certificates (see bb_lab.twist_floors):

  band (1,10)  complete+empty ⟹ no light contribution; nonempty ⟹
               per-w window floors via budgeted K(T(w)) certs;
               incomplete ⟹ the item falls back to the v1 floor.
  band (11,12) complete ⟹ stratum ≥ 12 + 2·min masked_K(w, 3);
               incomplete ⟹ 12 (|v| ≥ |w|).
  band (13,14) complete ⟹ 14 + 2·min masked_K(w, 2); else 14.
  band (15,16) complete ⟹ 16 + 2·min masked_K(w, 1); else 16.
  floor(λ) = min(light, s12, s14, s16, 18).

A masked value of 0 in band 12/14/16 would exhibit a cover LOGICAL of
weight ≤ 16 (λ ≠ 0 forces nontriviality) — the sweep re-verifies and
reports it loudly instead of writing a floor (it would be a d ≤ 16
witness for bb_288).

The dual-deck assembly (rank [Λ06|Λ60] = 12, verified: no class is
dangerous for both decks) gives g[c] = max over decks with
Λ_D·c ≠ 0 of floor_D(Λ_D·c) — the dangerous fiber is never consulted.

Usage:
  uv run python scripts/twist_floor_sweep.py --out DIR [--jobs 8]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.descent_sat import compute_descent  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.shard_distance import (  # noqa: E402
    _bits_to_int,
    _int_to_bits,
    compute_class_action,
)
from bb_lab.twist_floors import (  # noqa: E402
    compute_twist,
    enum_band,
    masked_K_value,
    syndrome_min_budgeted,
    window_floor,
)

_CTX: dict = {}


def _setup():
    G = ZmZn(12, 12)
    A = Poly.from_string("x^3 + y^2 + y^7", G)
    B = Poly.from_string("y^3 + x + x^2", G)
    checks = bb_check_matrices(A, B)
    action = compute_class_action(checks)
    for sigma in ((0, 6), (6, 0)):
        dd = compute_descent(checks, action, sigma, A=A, B=B)
        td = compute_twist(checks, dd)
        _CTX[str(sigma)] = (dd, td)
    _CTX["k"] = action.k
    return checks, action


def _run_item(args):
    key, lam = args
    dd, td = _CTX[key]
    out = {"deck": key, "lam": int(lam), "strata": {}}
    t_start = time.perf_counter()

    def band(wmin, wmax, limit, wall):
        ws, complete = enum_band(
            td, dd.Lb_Z, dd.kb, lam, wmin, wmax,
            limit=limit, wall_budget=wall,
        )
        out["strata"][f"{wmin}-{wmax}"] = {
            "n": len(ws), "complete": complete,
        }
        return ws, complete

    # light band
    light_bound = 10**9
    ws, complete = band(1, 10, 500, 300)
    if not complete:
        out["floor"] = None
        out["reason"] = "light band incomplete"
        return out
    for w in ws:
        wt = int(w.sum())
        kthr = (18 + wt + 1) // 2
        K, exact = syndrome_min_budgeted(
            td.Hb_Z, (td.T @ w) % 2, cap=kthr, confl_budget=300_000
        )
        light_bound = min(light_bound, window_floor(wt, K))

    floors = [light_bound]
    for (lo, hi, kneed, limit, wall) in (
        (11, 12, 3, 2000, 400),
        (13, 14, 2, 2000, 500),
        (15, 16, 1, 4000, 700),
    ):
        ws, complete = band(lo, hi, limit, wall)
        if not complete:
            floors.append(hi)          # |v| ≥ |w| ≥ lo ⟹ ≥ hi by parity
            continue
        if not ws:
            continue                   # empty stratum contributes nothing
        mk = min(masked_K_value(td, w, kneed) for w in ws)
        if mk == 0:
            # a b ⊆ supp(w) solves the twist: cost-|w| LOGICAL exists!
            out["floor"] = None
            out["reason"] = f"RED-ALERT masked 0 in band {lo}-{hi}"
            return out
        # every member has |w| = hi (odd weights are parity-excluded),
        # so the stratum bound is hi + 2·(min masked value, ≤ kneed).
        floors.append(hi + 2 * mk)
    out["floor"] = int(min(min(floors), 18))
    out["seconds"] = round(time.perf_counter() - t_start, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    checks, action = _setup()
    k = action.k
    items = []
    for key in ("(0, 6)", "(6, 0)"):
        dd, _ = _CTX[key]
        lams = sorted(
            {
                _bits_to_int((_int_to_bits(c, k) @ dd.Lam) % 2)
                for c in range(1, 1 << k)
            }
            - {0}
        )
        items += [(key, lam) for lam in lams]
    print(f"{len(items)} (deck, λ) items", flush=True)

    results = []
    with mp.get_context("fork").Pool(args.jobs) as pool:
        for res in pool.imap_unordered(_run_item, items):
            results.append(res)
            print(
                f"[{len(results):3}/{len(items)}] {res['deck']} "
                f"λ={res['lam']}: floor={res.get('floor')} "
                f"{res.get('reason', '')} "
                f"{res.get('seconds', '')}s "
                f"{ {b: (s['n'], s['complete']) for b, s in res['strata'].items()} }",
                flush=True,
            )
    (outdir / "twist_floors.json").write_text(json.dumps(results))

    from collections import Counter
    per_deck = {}
    for res in results:
        per_deck.setdefault(res["deck"], {})[res["lam"]] = res.get("floor")
    for key, m in per_deck.items():
        vals = Counter(v for v in m.values())
        print(f"{key}: floor histogram {dict(sorted((str(k_), v) for k_, v in vals.items()))}")

    # dual-deck assembly over the 4095 classes
    dd06, _ = _CTX["(0, 6)"]
    dd60, _ = _CTX["(6, 0)"]
    g = [0] * (1 << k)
    fallback = 0
    for c in range(1, 1 << k):
        cb = _int_to_bits(c, k)
        best = 0
        for key, dd in (("(0, 6)", dd06), ("(6, 0)", dd60)):
            lam = _bits_to_int((cb @ dd.Lam) % 2)
            if lam == 0:
                continue
            f = per_deck.get(key, {}).get(lam)
            if f is not None:
                best = max(best, f)
        if best == 0:
            fallback += 1
        g[c] = best
    from collections import Counter as C2
    print("assembled g histogram:", dict(sorted(C2(g[1:]).items())),
          f"fallback classes: {fallback}")
    (outdir / "g_table.json").write_text(json.dumps(g))


if __name__ == "__main__":
    main()
