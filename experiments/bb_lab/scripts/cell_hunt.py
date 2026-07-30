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

Usage (from experiments/bb_lab):
  uv run python scripts/cell_hunt.py --binary <tandem> \\
      --ell 7 --m 11 --k 6 --samples 40000 --out hits.csv --work-dir w
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import subprocess
import time
from pathlib import Path

from bb_lab.automorphism import automorphisms
from bb_lab.canonical import build_perm_table, canonical_pair
from bb_lab.checks import bb_check_matrices
from bb_lab.codeparams import code_params
from bb_lab.group import ZmZn
from bb_lab.linalg import nullspace_f2, quotient_complement_basis
from bb_lab.maxsat_distance import maxsat_distance
from bb_lab.poly import Poly


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
    ap.add_argument("--per-code-timeout", type=float, default=900.0)
    ap.add_argument("--deadline", type=float, default=0.0)
    ap.add_argument("--work-dir", default=".")
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

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    out = Path(args.out)
    new = not out.exists()
    t_start = time.perf_counter()
    seen: set = set()
    tested = 0

    crit = f"k == {args.k}" if args.k else f"k >= {k_min} (q>={args.merit})"
    print(f"Z{args.ell}xZ{args.m} n={n}, keeping {crit}, "
          f"{args.samples} draws", flush=True)

    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["n", "k", "d", "q", "seconds", "status",
                        "A_poly", "B_poly"])
        for _ in range(args.samples):
            if args.deadline and time.perf_counter() - t_start > args.deadline:
                print("deadline reached", flush=True)
                break
            A = Poly.from_support(rng.sample(els, 3), G)
            B = Poly.from_support(rng.sample(els, 3), G)
            ch = bb_check_matrices(A, B)
            k = code_params(ch).k
            if (args.k and k != args.k) or (not args.k and k < k_min):
                continue
            key = canonical_pair(A.support, B.support, G,
                                 auts=auts, perms=perms).key
            if key in seen:
                continue
            seen.add(key)
            tested += 1
            hx_even = not any(int(r.sum()) % 2 for r in ch.H_X)
            V = quotient_complement_basis(ch.H_X, nullspace_f2(ch.H_Z))
            step = hx_even and not any(int(v.sum()) % 2 for v in V)
            t0 = time.perf_counter()
            try:
                r = maxsat_distance(
                    ch, args.binary, mode="naive", work_dir=work,
                    timeout=args.per_code_timeout,
                    extra_args=("-cost-step=2",) if step else (),
                )
                d, secs, status = r.distance, r.solver_seconds, "ok"
            except subprocess.TimeoutExpired:
                d, secs, status = None, time.perf_counter() - t0, "timeout"
            except Exception as e:
                d, secs, status = None, time.perf_counter() - t0, \
                    f"error:{type(e).__name__}"
            q = round(k * d * d / n, 3) if d else None
            w.writerow([n, k, d, q, round(secs, 2), status,
                        A.canonical_string(), B.canonical_string()])
            f.flush()
            if d and d >= 16:
                print(f"  *** [[{n},{k},{d}]] q={q} A={A.canonical_string()} "
                      f"B={B.canonical_string()} ({secs:.0f}s)", flush=True)
    print(f"tested {tested} distinct codes", flush=True)


if __name__ == "__main__":
    main()
