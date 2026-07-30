"""Merit sweep: settle every corpus row that could rival the gross code.

The comparison quantity is the BB figure of merit

    q = k * d^2 / n            gross [[144,12,12]] has q = 12.

Every corpus row has weight-6 checks (|A| = |B| = 3 on all 58,021
rows), so (n, k, d) are the only free parameters.  Under a distance
cap d <= D, a row can only reach q >= QMIN if k * D^2 / n >= QMIN --
rows failing that are excluded by arithmetic and never queried, which
is what makes a whole-corpus sweep affordable.

Per code: naive WCNF + Tandem, with `-cost-step=2` passed only when
the coset weight-parity premise is verified to hold for that code
(same gate as ladder_sweep.py).

Measured on the 2026-07-29 run: 28,003 rows at qmin=8, D=14, zero
timeouts, ~0.32 s/row mean over 6 workers.  Rows are ordered by
descending best-case q so the informative ones land first.

Usage (from experiments/bb_lab), sharded across N workers:
  uv run python scripts/merit_sweep.py --binary <tandem> \\
      --shard 0 --nshards 6 --qmin 8 --out shard0.csv --work-dir w0
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from pathlib import Path

import duckdb

from bb_lab.checks import bb_check_matrices
from bb_lab.group import ZmZn
from bb_lab.linalg import nullspace_f2, quotient_complement_basis
from bb_lab.maxsat_distance import maxsat_distance
from bb_lab.poly import Poly

DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

SQL = """
SELECT instance_id, group_struct, ell, m, n, k, A_poly, B_poly, d_ub
FROM bb_instances
WHERE d_exact IS NULL AND k > 0 AND k * {dcap} * {dcap} / n >= {qmin}
ORDER BY k * {dcap} * {dcap} / n DESC, n ASC, instance_id
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--qmin", type=float, default=8.0,
                    help="keep rows whose best-case q at d=DCAP is >= QMIN")
    ap.add_argument("--dcap", type=float, default=14.0,
                    help="distance cap defining 'best case' in the screen")
    ap.add_argument("--per-code-timeout", type=float, default=600.0)
    ap.add_argument("--work-dir", default=".")
    ap.add_argument("--out", default="merit_sweep_results.csv")
    args = ap.parse_args()

    con = duckdb.connect(DB, read_only=True)
    rows = con.execute(SQL.format(qmin=args.qmin, dcap=args.dcap)).fetchall()
    con.close()
    rows = [r for i, r in enumerate(rows) if i % args.nshards == args.shard]

    out = Path(args.out)
    done: set[str] = set()
    if out.exists():
        with out.open() as f:
            done = {ln.split(",")[0] for ln in f.read().splitlines()[1:]}
    new = not out.exists()

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    print(f"shard {args.shard}/{args.nshards}: {len(rows)} rows "
          f"({len(done)} already done)", flush=True)

    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["instance_id", "group", "n", "k", "d", "q",
                        "d_ub_known", "seconds", "cost_step", "status",
                        "A_poly", "B_poly"])
        for iid, gs, ell, m, n, k, ap_, bp_, dub in rows:
            if iid in done:
                continue
            G = ZmZn(int(ell), int(m))
            ch = bb_check_matrices(Poly.from_string(ap_, G),
                                   Poly.from_string(bp_, G))
            # -cost-step=2 is sound only when every H_X row and every
            # class representative has even weight; verify per code.
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
            w.writerow([iid, gs, n, k, d, q, dub, round(secs, 2),
                        int(step), status, ap_, bp_])
            f.flush()
            if d and q and q >= 8:
                print(f"  ** [[{n},{k},{d}]] q={q}  {gs}  A={ap_} B={bp_}",
                      flush=True)


if __name__ == "__main__":
    main()
