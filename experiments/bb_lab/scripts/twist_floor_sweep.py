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
    masked_K_sat_budgeted,
    masked_K_value,
    syndrome_min_budgeted,
    window_floor,
)

_CTX: dict = {}


def _setup(args):
    """Default: the bb_288 dual-deck configuration. With --base-id
    and --axis: the literal lift of that corpus row along the axis,
    with every available cover axis deck installed (dual when both
    cover axes are even). Sets _CTX['bands'] from (d̄, target)."""
    if args.base_id:
        import duckdb

        con = duckdb.connect(
            "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/"
            "bb_instances.duckdb", read_only=True,
        )
        ell, m, d_base, a_s, b_s = con.execute(
            "SELECT ell, m, d_exact, A_poly, B_poly FROM bb_instances "
            "WHERE instance_id = ?", [args.base_id],
        ).fetchone()
        con.close()
        ell, m, d_base = int(ell), int(m), int(d_base)
        lift = (2 * ell, m) if args.axis == 0 else (ell, 2 * m)
        G = ZmZn(*lift)
        A = Poly.from_string(a_s, G)
        B = Poly.from_string(b_s, G)
        target = args.target or 2 * d_base
    else:
        G = ZmZn(12, 12)
        A = Poly.from_string("x^3 + y^2 + y^7", G)
        B = Poly.from_string("y^3 + x + x^2", G)
        d_base, target = 12, args.target or 18
    checks = bb_check_matrices(A, B)
    action = compute_class_action(checks)
    ellc, mc = G.orders
    decks = []
    if ellc % 2 == 0:
        decks.append((ellc // 2, 0))
    if mc % 2 == 0:
        decks.append((0, mc // 2))
    for sigma in decks:
        dd = compute_descent(checks, action, sigma, A=A, B=B)
        td = compute_twist(checks, dd)
        _CTX[str(sigma)] = (dd, td)
    _CTX["k"] = action.k
    # Bands from the base distance to the target: the light band
    # (1, d̄−2) must be empty for λ ≠ 0 (base-distance theorem — a
    # member is a RED-ALERT); then two-wide bands up to target−2 with
    # masked need = ceil((target − hi)/2) each.
    bands = []
    hi = d_base
    while hi <= target - 2:
        need = (target - hi + 1) // 2
        bands.append((hi - 1, hi, need))
        hi += 2
    _CTX["light_max"] = d_base - 2
    _CTX["bands"] = bands
    _CTX["target"] = target
    _CTX["group_label"] = G.label()
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

    target = _CTX["target"]
    # light band: empty by the base-distance theorem for λ ≠ 0; any
    # member is either a base-distance violation or a machinery bug.
    light_max = _CTX["light_max"]
    light_bound = 10**9
    ws, complete = band(1, light_max, 500, 300)
    if not complete:
        out["floor"] = None
        out["reason"] = "light band incomplete"
        return out
    for w in ws:
        wt = int(w.sum())
        kthr = (target + wt + 1) // 2
        K, exact = syndrome_min_budgeted(
            td.Hb_Z, (td.T @ w) % 2, cap=kthr, confl_budget=300_000
        )
        light_bound = min(light_bound, window_floor(wt, K))

    floors = [light_bound]
    for (lo, hi, kneed) in _CTX["bands"]:
        limit = 2000 if kneed >= 2 else 4000
        wall = 400 + 100 * (3 - min(kneed, 3))
        ws, complete = band(lo, hi, limit, wall)
        if not complete:
            floors.append(hi)          # |v| ≥ |w| ≥ lo ⟹ ≥ hi by parity
            continue
        if not ws:
            continue                   # empty stratum contributes nothing
        mks = []
        for w in ws:
            if kneed <= 3:
                mks.append(masked_K_value(td, w, kneed))
            else:
                j, b = masked_K_sat_budgeted(td, w, kneed)
                if b is not None and j == 0:
                    mks.append(0)
                else:
                    mks.append(j)
        mk = min(mks)
        if mk == 0:
            # a b ⊆ supp(w) solves the twist: cost-|w| LOGICAL exists!
            out["floor"] = None
            out["reason"] = f"RED-ALERT masked 0 in band {lo}-{hi}"
            return out
        # every member has |w| = hi (odd weights are parity-excluded),
        # so the stratum bound is hi + 2·(min masked value, ≤ kneed).
        floors.append(hi + 2 * mk)
    out["floor"] = int(min(min(floors), target))
    out["seconds"] = round(time.perf_counter() - t_start, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument(
        "--resume", default=None,
        help="prior sweep log: completed '[i/n] (deck) λ=..: floor=..' "
        "lines are reused and their items skipped",
    )
    ap.add_argument(
        "--base-id", default=None,
        help="corpus instance_id of a BASE code: sweep its literal "
        "lift along --axis instead of bb_288",
    )
    ap.add_argument("--axis", type=int, default=0, choices=(0, 1))
    ap.add_argument(
        "--target", type=int, default=None,
        help="certification target (default 2·d̄ for lifts, 18 for bb_288)",
    )
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    checks, action = _setup(args)
    k = action.k
    deck_keys = [kk for kk in _CTX if kk.startswith("(")]
    print(
        f"group {_CTX['group_label']} target {_CTX['target']} "
        f"decks {deck_keys} bands {_CTX['bands']}",
        flush=True,
    )
    items = []
    for key in deck_keys:
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
    if args.resume:
        import re

        pat = re.compile(
            r"^\[\s*\d+/\d+\] \((\d+), (\d+)\) λ=(\d+): floor=(\d+|None)"
        )
        for line in Path(args.resume).read_text().splitlines():
            mo = pat.match(line)
            if not mo:
                continue
            key = f"({mo.group(1)}, {mo.group(2)})"
            lam = int(mo.group(3))
            fl = None if mo.group(4) == "None" else int(mo.group(4))
            if fl is not None:
                results.append(
                    {"deck": key, "lam": lam, "floor": fl,
                     "strata": {}, "resumed": True}
                )
        done = {(r["deck"], r["lam"]) for r in results}
        items = [it for it in items if it not in done]
        print(
            f"resumed {len(results)} items from log; {len(items)} to run",
            flush=True,
        )
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

    # per-class assembly (max over decks with nonzero Λ-image)
    g = [0] * (1 << k)
    fallback = 0
    for c in range(1, 1 << k):
        cb = _int_to_bits(c, k)
        best = 0
        for key in deck_keys:
            dd, _ = _CTX[key]
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
