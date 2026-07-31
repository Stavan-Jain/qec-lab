"""Benchmark: descent-strengthened shard SAT vs the strengthened baseline.

For each code, runs the shard_distance coarse baseline and
descent_shard_distance for every available axis deck σ, asserts all
distances agree with the published value, and prints per-sector
breakdowns. Metric of record is solver-CPU (sum of task seconds);
wall includes the probe overhead both drivers share.

Usage (from experiments/bb_lab):
    uv run python scripts/descent_bench.py
    uv run python scripts/descent_bench.py --codes gross --no-baseline
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import yaml

LAB_ROOT = Path(__file__).resolve().parent.parent

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.descent_sat import (  # noqa: E402
    axis_decks,
    descent_shard_distance,
)
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.shard_distance import shard_distance  # noqa: E402


def load_instances() -> dict[str, dict]:
    path = LAB_ROOT / "instances" / "bravyi_table.yaml"
    data = yaml.safe_load(path.read_text())
    return {row["code_id"]: row for row in data["instances"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--codes", nargs="*",
        default=["bb_72_12_6", "bb_108_8_10", "gross"],
    )
    ap.add_argument("--no-baseline", action="store_true")
    ap.add_argument(
        "--variants", default="plain",
        help="comma list of plain,floors,scaffold,full",
    )
    args = ap.parse_args()

    table = load_instances()
    rows = []
    for code_id in args.codes:
        row = table[code_id]
        G = ZmZn(row["group"]["ell"], row["group"]["m"])
        A = Poly.from_string(row["polynomials"]["A"], G)
        B = Poly.from_string(row["polynomials"]["B"], G)
        checks = bb_check_matrices(A, B)
        d_pub = row["parameters"]["d"]

        base_cpu = base_wall = None
        if not args.no_baseline:
            t0 = time.perf_counter()
            rb = shard_distance(checks, backend="cms", verbose=False)
            base_wall = time.perf_counter() - t0
            base_cpu = sum(s.seconds for s in rb.shard_stats)
            assert rb.distance == d_pub, (code_id, rb.distance, d_pub)

        for sigma in axis_decks(checks):
            for vname in args.variants.split(","):
                opts = {
                    "plain": {},
                    "floors": {"use_floors": True},
                    "scaffold": {"scaffold_rest": True},
                    "full": {"use_floors": True, "scaffold_rest": True},
                }[vname]
                if opts.get("scaffold_rest") and len(axis_decks(checks)) < 2:
                    continue
                print(f"=== {code_id}  σ={sigma}  [{vname}]", flush=True)
                t0 = time.perf_counter()
                rd = descent_shard_distance(
                    checks, A, B, sigma=sigma, verbose=True, **opts
                )
                wall = time.perf_counter() - t0
                assert rd.distance == d_pub, (code_id, rd.distance, d_pub)
                for s in sorted(
                    rd.sector_stats, key=lambda s: -s.seconds
                )[:4]:
                    print(
                        f"    {s.sector:<4} grp={s.group} reps={s.n_reps:<4}"
                        f" w≤{s.weight:<3} {s.status:<6} {s.seconds:8.2f}s"
                    )
                rows.append(
                    (code_id, sigma, vname, d_pub, rd.solver_cpu,
                     rd.floor_seconds, wall, base_cpu, base_wall,
                     rd.a_killed_reps, rd.num_orbits)
                )
                print(
                    f"    descent CPU {rd.solver_cpu:.2f}s"
                    f" (+{rd.floor_seconds:.2f}s floors) wall {wall:.1f}s"
                    + (
                        f"  |  baseline CPU {base_cpu:.2f}s  |  ratio "
                        f"{base_cpu / max(rd.solver_cpu + rd.floor_seconds, 1e-9):.2f}×"
                        if base_cpu is not None else ""
                    ),
                    flush=True,
                )

    print("\n=== summary (CPU = SAT task seconds + floor precompute) ===")
    print(
        f"{'code':<14} {'σ':<9} {'variant':<9} {'d':>3} {'descent':>9}"
        f" {'baseline':>9} {'speedup':>8}  a-killed/orbits"
    )
    for (code_id, sigma, vname, d, cpu, fsec, wall, bcpu, bwall,
         akill, norb) in rows:
        tot = cpu + fsec
        sp = f"{bcpu / max(tot, 1e-9):.2f}×" if bcpu else "-"
        print(
            f"{code_id:<14} {str(sigma):<9} {vname:<9} {d:>3} {tot:>8.2f}s"
            f" {bcpu if bcpu else 0:>8.2f}s {sp:>8}  {akill}/{norb}"
        )


if __name__ == "__main__":
    main()
