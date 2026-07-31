"""Fresh-code hunt inside one (group, k) cell, or along the merit line.

The corpus is only a *sample* of each cell -- n=168 carries 800 rows
per group out of a raw pair space of order 10^9 -- so "no such code
among the stored rows" is never "no such code".  This script draws new
weight-3 (A, B) pairs, filters on k, canonical-dedups, and settles the
survivors exactly with Tandem.

Two selection modes:

  --k K          keep codes with k exactly K.  Use to chase a distance
                 record in a cell where one is known to exist -- e.g.
                 Z7xZ11 at k=6 hosts the published Wang-Mueller
                 [[154,6,16]].
  --merit QMIN   keep codes whose k could reach q = k d^2/n >= QMIN at
                 the distance cap, i.e. k >= QMIN*n/DCAP^2.  Use to
                 chase the figure of merit rather than the distance.

Measured (2026-07-29, --merit 12 --dcap 14, six groups, 360k draws):
3,191 k-qualifying codes, all solved, zero timeouts, zero reaching
q >= 12 -- the k-filter and a high distance pull against each other,
and the surviving distance distribution peaks at d=8.

The run is two phases rather than one interleaved loop. Drawing,
k-filtering and canonical dedup are pure Python and cheap; solving is
everything else (the 2026-07-29 larger-n hunt spent 44.5 h here, one
row of it 6,157 s). So phase 1 assembles the survivor list serially and
phase 2 puts it on `bb_lab.sweep`'s dynamic queue. Splitting them also
buys resume, which this script never had: survivors are keyed by their
canonical (A_poly, B_poly) pair, which the CSV already carries, so
rerunning the same command skips what is already settled.

Note the draw is seeded (`--seed`), so phase 1 reproduces its candidate
list exactly; resume therefore lines up even when `--samples` grows.

Usage (from experiments/bb_lab):
  uv run python scripts/cell_hunt.py --binary <tandem> \\
      --ell 7 --m 11 --k 6 --samples 40000 --out hits.csv
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

from bb_lab.automorphism import automorphisms
from bb_lab.canonical import build_perm_table, canonical_pair
from bb_lab.checks import bb_check_matrices
from bb_lab.codeparams import code_params
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.sweep import bb_distance_task, default_jobs, run_sweep

# Unchanged from the serial writer, so existing hunt CSVs stay
# appendable. `bb_distance_task` also returns cost_step; run_sweep only
# writes the columns named here, so dropping it costs nothing.
FIELDNAMES = ["n", "k", "d", "q", "seconds", "status", "A_poly", "B_poly"]


def hunt_task(payload: dict) -> dict:
    """Distance solve plus the figure of merit it feeds."""
    row = bb_distance_task(payload)
    d = row.get("d")
    row["q"] = round(row["k"] * d * d / row["n"], 3) if d else None
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--ell", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--k", type=int, default=0, help="keep k exactly this")
    ap.add_argument("--merit", type=float, default=0.0,
                    help="keep k >= MERIT*n/DCAP^2 instead of a fixed k")
    ap.add_argument("--dcap", type=float, default=14.0)
    ap.add_argument("--samples", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=0)
    # Raised from 900 s with the move to a queue: concurrent solves run
    # ~20% slower each (shared cache / memory bandwidth), so a cap tuned
    # on serial runs drops rows that would otherwise settle -- and in a
    # distance hunt those are the interesting ones.
    ap.add_argument("--per-code-timeout", type=float, default=1800.0)
    ap.add_argument("--jobs", type=int, default=None,
                    help="workers (default: performance-core count)")
    ap.add_argument("--deadline", type=float, default=0.0)
    ap.add_argument("--work-dir", default="cell_hunt_work")
    ap.add_argument("--out", default="cell_hunt_results.csv")
    args = ap.parse_args()
    if not args.k and not args.merit:
        ap.error("give --k or --merit")

    G = ZmZn(args.ell, args.m)
    n = 2 * args.ell * args.m
    k_min = math.ceil(args.merit * n / (args.dcap ** 2)) if args.merit else 0
    els = [(i, j) for i in range(args.ell) for j in range(args.m)]
    rng = random.Random(args.seed)
    auts = automorphisms(G)
    perms = build_perm_table(G, auts)

    crit = f"k == {args.k}" if args.k else f"k >= {k_min} (q>={args.merit})"
    print(f"Z{args.ell}xZ{args.m} n={n}, keeping {crit}, "
          f"{args.samples} draws", flush=True)

    # Phase 1 -- draw, k-filter, canonical-dedup. Pure Python and cheap
    # next to the solving; kept serial so the candidate list is exactly
    # reproducible from --seed.
    seen: set = set()
    items: list[dict] = []
    for _ in range(args.samples):
        A = Poly.from_support(rng.sample(els, 3), G)
        B = Poly.from_support(rng.sample(els, 3), G)
        k = code_params(bb_check_matrices(A, B)).k
        if (args.k and k != args.k) or (not args.k and k < k_min):
            continue
        key = canonical_pair(A.support, B.support, G,
                             auts=auts, perms=perms).key
        if key in seen:
            continue
        seen.add(key)
        a_str, b_str = A.canonical_string(), B.canonical_string()
        items.append({
            "ell": args.ell, "m": args.m,
            "A_poly": a_str, "B_poly": b_str,
            "binary": args.binary, "mode": "naive",
            "timeout": args.per_code_timeout,
            "passthrough": {"n": n, "k": k,
                            "A_poly": a_str, "B_poly": b_str},
        })
    print(f"phase 1: {len(items)} distinct k-qualifying codes drawn",
          flush=True)

    # Phase 2 -- settle them on the dynamic queue. Resume keys off the
    # canonical (A_poly, B_poly) pair the CSV already carries.
    def report(row: dict) -> str | None:
        d = row.get("d")
        if d and d >= 16:
            return (f"  *** [[{row['n']},{row['k']},{d}]] q={row['q']} "
                    f"A={row['A_poly']} B={row['B_poly']} "
                    f"({row['seconds']:.0f}s)")
        if row.get("status") not in (None, "ok"):
            return (f"  !! A={row['A_poly']} B={row['B_poly']} "
                    f"{row['status']}")
        return None

    print(f"phase 2: solving, jobs={args.jobs or default_jobs()}", flush=True)
    run_sweep(
        items, hunt_task,
        out=Path(args.out), fieldnames=FIELDNAMES,
        key_field=("A_poly", "B_poly"),
        key=lambda it: (it["A_poly"], it["B_poly"]),
        jobs=args.jobs, work_root=Path(args.work_dir), report=report,
        deadline=args.deadline or None,
    )


if __name__ == "__main__":
    main()
