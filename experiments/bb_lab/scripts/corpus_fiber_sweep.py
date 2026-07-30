"""Corpus A/B: Tandem naive+step vs +self-certified fiber tables.

Honest solve-mode protocol — the solver is told NOTHING it did not
derive itself: no -init-lb, no seeds, no stored distances. Per (code,
deck): the caller computes a greedy/probe upper bound U (its own
witness), emits the budgeted certificate with caps tied to U (v1 base
floors + v2 moving cost floors + dangerous-fiber pin + invariant
floor), then races naive+step against naive+step+fiber. Stored d_exact
is used only to ASSERT correctness after the fact.

The self-certified collapse condition min(invFloor, g_min) ≥ d_found
is the doubling signature: codes meeting it are flagged
DOUBLING-LIKE (the per-code analytic structure alone determines the
proof phase). Large fiber speedups without full collapse mark partial
structure.

Usage: uv run python scripts/corpus_fiber_sweep.py --binary BIN \
    --out DIR [--jobs 2] [--limit 36]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import statistics
import time
from pathlib import Path

import duckdb
import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.descent_sat import axis_decks  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.linalg import nullspace_f2, quotient_complement_basis  # noqa: E402
from bb_lab.maxsat_distance import (  # noqa: E402
    emit_fiber_certificate,
    maxsat_distance,
    write_wcnf,
)
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.sat_distance import find_logical_z  # noqa: E402
from bb_lab.shard_distance import (  # noqa: E402
    _greedy_upper_bound,
    compute_class_action,
)

SAMPLE_SQL = """
WITH ranked AS (
  SELECT *, ROW_NUMBER() OVER (
      PARTITION BY d_exact, n ORDER BY instance_id
  ) AS rk
  FROM bb_instances
  WHERE d_exact IS NOT NULL AND d_exact BETWEEN 6 AND 14
    AND k > 0 AND k <= 12 AND n <= 210
    AND (ell % 2 = 0 OR m % 2 = 0)
)
SELECT instance_id, ell, m, n, k, d_exact, A_poly, B_poly FROM ranked
WHERE rk <= 2
ORDER BY n ASC, d_exact DESC   -- small first: results stream early
LIMIT {limit}
"""

_BIN = None
_WORK = None


def _all_even(checks) -> bool:
    if any(int(r.sum()) % 2 for r in checks.H_X):
        return False
    V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
    L_Z = find_logical_z(checks)
    if V.shape[0] != L_Z.shape[0]:
        return False
    return not any(int(v.sum()) % 2 for v in V)


def sweep_one(row):
    iid, ell, m, n, k, d_exact, a_s, b_s = row
    G = ZmZn(int(ell), int(m))
    A = Poly.from_string(a_s, G)
    B = Poly.from_string(b_s, G)
    checks = bb_check_matrices(A, B)
    out = {"id": iid, "n": int(n), "k": int(k), "d": int(d_exact)}
    try:
        action = compute_class_action(checks)
    except (ValueError, AssertionError) as e:
        out["verdict"] = f"SKIP({e})"
        return out
    # self-computed upper bound (no stored knowledge): greedy witness,
    # tightened by cheap conflict-budgeted monolith probes — the caps
    # for the certificate budgets track the real optimum much closer
    # than the raw greedy weight.
    from bb_lab.sat_distance import find_logical_z as _flz
    from bb_lab.shard_distance import _monolith_probe

    wit = _greedy_upper_bound(action.V, checks.H_X)
    L_Z0 = _flz(checks)
    while int(wit.sum()) > 2:
        v2 = _monolith_probe(
            checks.H_Z, L_Z0, int(wit.sum()) - 1, 60_000
        )
        if v2 is None:
            break
        wit = v2
    U = int(wit.sum())
    out["probe_ub"] = U
    step = ("-cost-step=2",) if _all_even(checks) else ()

    best = None
    wdir = _WORK / iid          # per-code dir: group labels collide
    wdir.mkdir(parents=True, exist_ok=True)
    for sigma in axis_decks(checks):
        sig = "".join(str(int(x)) for x in sigma)
        tag = f"{iid}_{sig}"
        try:
            t0 = time.perf_counter()
            qv, a_lits = write_wcnf(
                checks, wdir / f"naive_{tag}.wcnf", mode="naive"
            )
            flb = wdir / f"fiber_{tag}.flb"
            emit_fiber_certificate(
                checks, A, B, sigma, qv, a_lits, flb,
                floor_cap=max(4, U - 1), floor_budget=60_000,
                cost_floors=True, cost_cap=max(4, U - 1),
                cost_budget=120_000,
            )
            emit_s = time.perf_counter() - t0
        except (AssertionError, ValueError) as e:
            out.setdefault("deck_errors", []).append((sig, str(e)))
            continue
        hdr = flb.read_text().split("\n", 4)
        inv = int(hdr[0].split()[4])
        gvals = [int(x) for x in hdr[3].split()[2:]]
        g_min = min(gvals) if gvals else 0

        ts, t0s = [], []
        okd = True
        for _ in range(2):
            rf = maxsat_distance(
                checks, _BIN, mode="naive", work_dir=wdir,
                extra_args=step + (f"-fiber-lb={flb}",),
            )
            r0 = maxsat_distance(
                checks, _BIN, mode="naive", work_dir=wdir,
                extra_args=step,
            )
            okd &= rf.distance == d_exact == r0.distance
            ts.append(rf.solver_seconds)
            t0s.append(r0.solver_seconds)
        rec = {
            "sigma": sig, "emit_s": round(emit_s, 1),
            "invFloor": inv, "g_min": g_min,
            "fiber_s": round(statistics.median(ts), 3),
            "step_s": round(statistics.median(t0s), 3),
            "ratio": round(
                statistics.median(t0s) / max(statistics.median(ts), 1e-9), 2
            ),
            "d_ok": okd,
            "collapse": bool(min(inv, g_min) >= d_exact),
        }
        if best is None or rec["ratio"] > best["ratio"]:
            best = rec
        out.setdefault("decks", []).append(rec)
    if best is None:
        out["verdict"] = out.get("verdict", "NO-DECK-USABLE")
    else:
        out["verdict"] = (
            "DOUBLING-LIKE" if best["collapse"]
            else ("FAST" if best["ratio"] >= 1.5 else "NEUTRAL")
        )
        out["best"] = best
    return out


def _init(binary, work):
    global _BIN, _WORK
    _BIN = binary
    _WORK = Path(work)
    _WORK.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--limit", type=int, default=36)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(DB, read_only=True)
    rows = con.execute(SAMPLE_SQL.format(limit=args.limit)).fetchall()
    con.close()
    print(f"{len(rows)} corpus rows", flush=True)

    results = []
    with mp.get_context("fork").Pool(
        args.jobs, initializer=_init,
        initargs=(args.binary, str(outdir / "wcnf")),
    ) as pool:
        for res in pool.imap_unordered(sweep_one, rows):
            results.append(res)
            print(json.dumps(res), flush=True)
    (outdir / "corpus_fiber.json").write_text(json.dumps(results))

    print("\n=== summary (best deck per code) ===", flush=True)
    for r in sorted(
        results,
        key=lambda r: -(r.get("best", {}).get("ratio", 0)),
    ):
        b = r.get("best")
        if not b:
            print(f"[[{r['n']},{r['k']},{r['d']}]] {r['id'][:8]}: "
                  f"{r['verdict']}")
            continue
        print(
            f"[[{r['n']},{r['k']},{r['d']}]] {r['id'][:8]} σ={b['sigma']}: "
            f"{b['step_s']}s → {b['fiber_s']}s ({b['ratio']}×) "
            f"inv={b['invFloor']} g_min={b['g_min']} "
            f"{'d-OK' if b['d_ok'] else 'D-MISMATCH!'} {r['verdict']}"
        )


if __name__ == "__main__":
    main()
