"""Insert cell_hunt discoveries into the corpus DB as new instances.

Unlike the sweep/battery write-back, these are codes the corpus never
contained, so they need an INSERT rather than a `d_exact` UPDATE.

Rows are stored the way `bb-lab enumerate` stores them, so a hunt
discovery is indistinguishable from an enumerated row afterwards:

  - the (A, B) actually drawn is replaced by its CANONICAL ORBIT
    REPRESENTATIVE under Aut(G) x G, so two draws from one orbit
    collapse to a single row (`Poly.canonical_string()` alone is only a
    deterministic rendering of one polynomial -- it is NOT the orbit
    canonicalisation, and using it would duplicate rows);
  - instance_id = canonical_hash(group, A, B), code_id carries a
    `bb_hunt_` prefix so the provenance stays legible next to
    `bb_enum_` / `bb_samp_` rows;
  - structural columns (k, ranks, kernel dims, orbit size) are
    recomputed here rather than trusted from the CSV.

`upsert_instance` COALESCEs distance fields, so an id that collides
with an existing corpus row will NOT have its stored d_exact
clobbered; any disagreement is reported and, unless --force, aborts.

Usage (from experiments/bb_lab):
  uv run python scripts/cell_hunt_ingest.py --csv 'results/bign_*.csv' --dry-run
  uv run python scripts/cell_hunt_ingest.py --csv 'results/bign_*.csv'
"""

from __future__ import annotations

import argparse
import csv
import datetime
import glob
import os
import re
import shutil
from pathlib import Path

import duckdb

from bb_lab.automorphism import automorphisms
from bb_lab.canonical import build_perm_table, canonical_pair
from bb_lab.checks import bb_check_matrices, circulant
from bb_lab.codeparams import code_params
from bb_lab.group import ZmZn
from bb_lab.linalg import rank_f2
from bb_lab.poly import Poly
from bb_lab.store import (
    StoredInstance,
    canonical_hash,
    connect,
    upsert_instance,
)

DB = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb")

# cell_hunt buffers its header until the first row flushes, and sibling
# workers on one cell can share a file, so parse positionally.
FIELDS = ["n", "k", "d", "q", "seconds", "status", "A_poly", "B_poly"]
NAME_RE = re.compile(r"_(\d+)x(\d+)_k(\d+)\.csv$")


def rows_from(path: str):
    """(ell, m, A_poly, B_poly, d) for each solved row; group from filename."""
    mo = NAME_RE.search(os.path.basename(path))
    if not mo:
        return
    ell, m, _k = (int(x) for x in mo.groups())
    with open(path) as fh:
        for row in csv.reader(fh):
            if not row or row[0] == "n" or len(row) < len(FIELDS):
                continue
            r = dict(zip(FIELDS, row))
            if r["status"] == "ok" and r["d"]:
                yield ell, m, r["A_poly"], r["B_poly"], int(r["d"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True, help="comma-separated globs")
    ap.add_argument("--db", default=str(DB))
    ap.add_argument("--method", default="maxsat-tandem@mse23+step2")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="proceed even if an id collides with a differing row")
    args = ap.parse_args()

    files = [f for pat in args.csv.split(",") for f in sorted(glob.glob(pat))]
    print(f"{len(files)} result files")

    tables: dict[tuple[int, int], tuple] = {}
    staged: dict[str, tuple[StoredInstance, int]] = {}
    n_rows = 0
    for f in files:
        for ell, m, ap_, bp_, d in rows_from(f):
            n_rows += 1
            if (ell, m) not in tables:
                G = ZmZn(ell, m)
                auts = automorphisms(G)
                tables[(ell, m)] = (G, auts, build_perm_table(G, auts))
            G, auts, perms = tables[(ell, m)]
            A0 = Poly.from_string(ap_, G)
            B0 = Poly.from_string(bp_, G)
            canon = canonical_pair(A0.support, B0.support, G,
                                   auts=auts, perms=perms)
            A = Poly(support=frozenset(canon.A_support), group=G)
            B = Poly(support=frozenset(canon.B_support), group=G)
            As, Bs = A.canonical_string(), B.canonical_string()
            iid = canonical_hash(G.label(), As, Bs)
            if iid in staged:
                continue
            ch = bb_check_matrices(A, B)
            p = code_params(ch)
            staged[iid] = (
                StoredInstance(
                    instance_id=iid,
                    code_id=f"bb_hunt_{G.label()}_{iid[:8]}",
                    group_struct=G.label(),
                    ell=ell, m=m, n=p.n, k=p.k,
                    A_poly=As, B_poly=Bs,
                    A_weight=A.weight(), B_weight=B.weight(),
                    rank_HX=p.rank_HX, rank_HZ=p.rank_HZ,
                    dim_ker_A=G.cardinality - rank_f2(circulant(A)),
                    dim_ker_B=G.cardinality - rank_f2(circulant(B)),
                    orbit_size=canon.orbit_size,
                    d_exact=d, d_method=args.method,
                ),
                d,
            )

    print(f"{n_rows} solved rows -> {len(staged)} distinct canonical orbits")

    con = duckdb.connect(args.db, read_only=True)
    fresh = collide_ok = collide_bad = 0
    for iid, (inst, d) in staged.items():
        got = con.execute(
            "SELECT d_exact FROM bb_instances WHERE instance_id=?", [iid]
        ).fetchone()
        if got is None:
            fresh += 1
        elif got[0] is None or got[0] == d:
            collide_ok += 1
        else:
            collide_bad += 1
            print(f"  CONFLICT {iid}: stored d={got[0]}, hunt d={d}")
    con.close()
    print(f"  {fresh} new rows, {collide_ok} already present and consistent, "
          f"{collide_bad} CONFLICTS")
    if collide_bad and not args.force:
        raise SystemExit("conflicts present — not writing (use --force to override)")
    if args.dry_run:
        return

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = Path(args.db).with_suffix(f".duckdb.bak-{stamp}")
    shutil.copy2(args.db, bak)
    print(f"backup: {bak}")

    with connect(args.db) as con:
        for inst, _d in staged.values():
            upsert_instance(con, inst)
    print(f"inserted/updated {len(staged)} instances")

    con = duckdb.connect(args.db, read_only=True)
    print("corpus coverage:", con.execute(
        "SELECT COUNT(d_exact), COUNT(*) FROM bb_instances").fetchone())
    print("d>=16 census:")
    for d, c in con.execute(
        "SELECT d_exact, COUNT(*) FROM bb_instances WHERE d_exact>=16"
        " GROUP BY 1 ORDER BY 1"
    ).fetchall():
        print(f"   d={d}: {c} rows")
    con.close()


if __name__ == "__main__":
    main()
