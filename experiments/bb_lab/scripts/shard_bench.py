"""Benchmark: shard-decomposition solver vs the monolithic SAT baseline.

Runs both exact-distance paths on instances from
`instances/bravyi_table.yaml` and prints a comparison table plus the
per-shard breakdown for the largest code run.

Usage (from experiments/bb_lab):
    uv run python scripts/shard_bench.py                    # 3 small codes
    uv run python scripts/shard_bench.py --codes gross --jobs 8
    uv run python scripts/shard_bench.py --no-monolith --codes bb_288_12_18 --jobs 8
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import yaml

LAB_ROOT = Path(__file__).resolve().parent.parent

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.sat_distance import x_distance  # noqa: E402
from bb_lab.shard_distance import shard_distance  # noqa: E402


def load_instances() -> dict[str, dict]:
    path = LAB_ROOT / "instances" / "bravyi_table.yaml"
    data = yaml.safe_load(path.read_text())
    return {row["code_id"]: row for row in data["instances"]}


def build_checks(row: dict):
    G = ZmZn(row["group"]["ell"], row["group"]["m"])
    A = Poly.from_string(row["polynomials"]["A"], G)
    B = Poly.from_string(row["polynomials"]["B"], G)
    return bb_check_matrices(A, B)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--codes", nargs="*",
        default=["bb_72_12_6", "bb_90_8_10", "bb_108_8_10"],
    )
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--orbits-per-shard", type=int, default=None,
                    help="granularity dial; default = coarsest (best)")
    ap.add_argument(
        "--backend", default="auto", choices=["auto", "cms", "cadical"]
    )
    ap.add_argument("--no-monolith", action="store_true",
                    help="skip the monolithic baseline (for slow codes)")
    ap.add_argument("--json", type=Path, default=None,
                    help="also dump results as JSON")
    args = ap.parse_args()

    table = load_instances()
    rows_out = []
    for code_id in args.codes:
        row = table[code_id]
        checks = build_checks(row)
        d_pub = row["parameters"]["d"]
        print(f"\n=== {row['display_name']} ({code_id}) ===", flush=True)

        mono_t = mono_d = None
        if not args.no_monolith:
            t0 = time.perf_counter()
            mono = x_distance(checks)
            mono_t = time.perf_counter() - t0
            mono_d = mono.distance
            print(f"  monolith: d = {mono_d}  in {mono_t:8.2f}s", flush=True)

        t0 = time.perf_counter()
        res = shard_distance(
            checks, backend=args.backend, jobs=args.jobs, verbose=True,
            orbits_per_shard=args.orbits_per_shard,
        )
        shard_t = time.perf_counter() - t0
        n_tasks = len(res.shard_stats)
        solve_cpu = sum(s.seconds for s in res.shard_stats)
        print(
            f"  shards:   d = {res.distance}  in {shard_t:8.2f}s wall "
            f"({solve_cpu:.2f}s solver-CPU, {n_tasks} shard tasks, "
            f"{res.num_orbits} orbits of {res.num_classes} classes, "
            f"{res.rounds} rounds, initial d_ub={res.initial_upper_bound})",
            flush=True,
        )
        slowest = sorted(
            res.shard_stats, key=lambda s: s.seconds, reverse=True
        )[:5]
        for s in slowest:
            print(
                f"    slowest: rep={s.rep:>6}  case={s.case}  "
                f"w≤{s.weight:<3}  {s.status:<5}  {s.seconds:8.2f}s"
            )

        ok = res.distance == d_pub and (mono_d in (None, d_pub))
        print(f"  published d = {d_pub}  →  {'OK' if ok else 'MISMATCH ✗'}")
        rows_out.append({
            "code_id": code_id, "n": row["parameters"]["n"],
            "k": row["parameters"]["k"], "d_published": d_pub,
            "d_shard": res.distance, "d_monolith": mono_d,
            "t_shard_wall": round(shard_t, 3),
            "t_shard_solver_cpu": round(solve_cpu, 3),
            "t_monolith": None if mono_t is None else round(mono_t, 3),
            "num_orbits": res.num_orbits, "rounds": res.rounds,
            "shard_tasks": n_tasks,
            "initial_d_ub": res.initial_upper_bound,
            "jobs": args.jobs, "backend": args.backend,
            "orbits_per_shard": args.orbits_per_shard,
        })

    print("\n=== summary ===")
    hdr = (f"{'code':<14} {'n':>4} {'d':>3} {'shard(s)':>9} "
           f"{'mono(s)':>9} {'speedup':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows_out:
        sp = ("-" if r["t_monolith"] is None
              else f"{r['t_monolith'] / r['t_shard_wall']:.1f}×")
        mono = "-" if r["t_monolith"] is None else f"{r['t_monolith']:.2f}"
        print(
            f"{r['code_id']:<14} {r['n']:>4} {r['d_shard']:>3} "
            f"{r['t_shard_wall']:>9.2f} {mono:>9} {sp:>8}"
        )

    if args.json:
        args.json.write_text(json.dumps(rows_out, indent=2))
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
