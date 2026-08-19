"""Order-144 sampled canonical enumeration (a18_sample_enum analogue).

Draws uniform random weight-3 (A, B) pairs on each requested order-144
group, keeps k >= 2, canonicalizes with the fast exact routine
(sweeplib), skips anything already carrying d_exact in the MAIN corpus
(read-only) plus the two published [[288]] codes, and upserts keepers
into THIS worktree's sweep DuckDB with code_id prefix `bb_samp_`
(sampled provenance, a18 convention).

  uv run python data/order144_sweep/sampler.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import duckdb
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import sweeplib  # noqa: E402
from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.checks import bb_check_matrices, circulant  # noqa: E402
from bb_lab.codeparams import code_params  # noqa: E402
from bb_lab.linalg import rank_f2  # noqa: E402
from bb_lab.store import StoredInstance, connect, upsert_instance  # noqa: E402

MAIN_DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"
SWEEP_DB = HERE / "sweep.duckdb"

# (ell, m, target); priority order.
QUOTAS = [
    (12, 12, 34),
    (18, 8, 6),
    (24, 6, 6),
    (36, 4, 4),
    (16, 9, 4),
    (48, 3, 4),
]

SEED = 20260817

# Published codes to skip even if absent from the corpus (both have
# certificate-tier exact d in the program's records).
PUBLISHED_SKIP = [
    # bb288 [[288,12,18]] (Bravyi et al table)
    ((12, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2"),
    # IBM class-Y [[288,8,20]] (A33; cover polys from DOUBLING_PAIRS_CONSOLIDATED)
    ((18, 8), "1 + x*y^4 + x^14*y", "1 + x*y^2 + x^2*y^7"),
]


def corpus_ids_with_exact(labels: list[str]) -> tuple[set[str], set[str]]:
    con = duckdb.connect(MAIN_DB, read_only=True)
    ph = ",".join("?" for _ in labels)
    rows = con.execute(
        f"SELECT instance_id, d_exact IS NOT NULL FROM bb_instances "
        f"WHERE group_struct IN ({ph})",
        labels,
    ).fetchall()
    con.close()
    with_exact = {r[0] for r in rows if r[1]}
    all_ids = {r[0] for r in rows}
    return with_exact, all_ids


def published_skip_ids() -> dict[str, str]:
    out = {}
    for orders, A_str, B_str in PUBLISHED_SKIP:
        A_supp = sweeplib.support_from_string(orders, A_str)
        B_supp = sweeplib.support_from_string(orders, B_str)
        _, cA, cB, cid, _ = sweeplib.canonicalize(*orders, A_supp, B_supp)
        out[cid] = f"published [[{2 * orders[0] * orders[1]}]] {A_str} ; {B_str}"
    return out


def sample_group(con, ell: int, m: int, target: int, *, skip_exact: set[str],
                 skip_pub: dict[str, str], corpus_all: set[str]) -> int:
    G = AbelianGroup((ell, m))
    N = G.cardinality
    label = G.label()
    t0 = time.time()
    gt = sweeplib.group_tables(ell, m)
    print(f"[{label}] |G|={N} |Aut|={gt.PHI.shape[0]} (tables {time.time() - t0:.1f}s)", flush=True)

    seen_ids: set[str] = {
        r[0] for r in con.execute(
            "SELECT instance_id FROM bb_instances WHERE group_struct = ?", [label]
        ).fetchall()
    }
    rng = np.random.default_rng(SEED + 1009 * ell + m)
    inserted = tried = k_rejects = dup_skips = exact_skips = 0
    t_start = time.time()
    max_tries = target * 600

    while inserted < target and tried < max_tries and time.time() - t_start < 420:
        tried += 1
        A_idx = np.sort(rng.choice(N, size=3, replace=False))
        B_idx = np.sort(rng.choice(N, size=3, replace=False))
        A_supp = tuple((int(i) // m, int(i) % m) for i in A_idx)
        B_supp = tuple((int(i) // m, int(i) % m) for i in B_idx)
        A = Poly(support=frozenset(A_supp), group=G)
        B = Poly(support=frozenset(B_supp), group=G)
        checks = bb_check_matrices(A, B)
        params = code_params(checks)
        if params.k < 2:
            k_rejects += 1
            continue
        _, A_str, B_str, iid, orbit = sweeplib.canonicalize(ell, m, A_supp, B_supp)
        if iid in seen_ids:
            dup_skips += 1
            continue
        seen_ids.add(iid)
        if iid in skip_pub:
            exact_skips += 1
            print(f"  [{label}] SKIP published: {skip_pub[iid]}", flush=True)
            continue
        if iid in skip_exact:
            exact_skips += 1
            print(f"  [{label}] SKIP corpus-exact instance {iid[:12]}", flush=True)
            continue
        dim_ker_A = N - rank_f2(circulant(A))
        dim_ker_B = N - rank_f2(circulant(B))
        upsert_instance(con, StoredInstance(
            instance_id=iid,
            code_id=f"bb_samp_{label}_{iid[:8]}",
            group_struct=label, ell=ell, m=m,
            n=params.n, k=params.k,
            A_poly=A_str, B_poly=B_str, A_weight=3, B_weight=3,
            rank_HX=params.rank_HX, rank_HZ=params.rank_HZ,
            dim_ker_A=dim_ker_A, dim_ker_B=dim_ker_B,
            orbit_size=orbit,
        ))
        if iid in corpus_all:
            con.execute(
                "UPDATE bb_instances SET d_method = 'corpus-overlap-no-exact' "
                "WHERE instance_id = ?", [iid])
        inserted += 1

    dt = time.time() - t_start
    print(f"[{label}] DONE {inserted}/{target} inserted, {tried} tried "
          f"(k<2: {k_rejects}, dup: {dup_skips}, exact-skip: {exact_skips}) {dt:.0f}s",
          flush=True)
    return inserted


def main() -> None:
    labels = [AbelianGroup((l, m)).label() for l, m, _ in QUOTAS]
    skip_exact, corpus_all = corpus_ids_with_exact(labels)
    print(f"corpus rows in sweep groups: {len(corpus_all)} ({len(skip_exact)} with d_exact)")
    skip_pub = published_skip_ids()
    print(f"published-skip canonical ids: {list(skip_pub.values())}")
    total = 0
    with connect(SWEEP_DB) as con:
        for ell, m, target in QUOTAS:
            total += sample_group(con, ell, m, target, skip_exact=skip_exact,
                                  skip_pub=skip_pub, corpus_all=corpus_all)
    print(f"TOTAL inserted: {total}")


if __name__ == "__main__":
    main()
