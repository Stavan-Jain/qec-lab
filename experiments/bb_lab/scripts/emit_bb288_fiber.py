"""Emit -fiber-lb certificates for bb_288, both axis decks, in parallel.

The floor table is the expensive part (one witness-jumping base-SAT
per distinct Λ-image class, base n = 144); this script chunks the
classes over a process pool, merges, and writes files byte-compatible
with `emit_fiber_certificate`. Also prints kb / distinct-class counts
and the invariant floors — the analytic pre-read on which deck should
prune harder.

Usage: uv run python scripts/emit_bb288_fiber.py --out DIR [--jobs 8]
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import yaml

LAB_ROOT = Path(__file__).resolve().parent.parent

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.descent_sat import (  # noqa: E402
    axis_decks,
    base_coset_floors_budgeted,
    compute_descent,
    invariant_floor,
)
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.maxsat_distance import write_wcnf  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.shard_distance import (  # noqa: E402
    _bits_to_int,
    _int_to_bits,
    compute_class_action,
)

_DD = {}


def _floor_chunk(args):
    key, classes, cap = args
    return base_coset_floors_budgeted(
        _DD[key], set(classes), cap=cap, confl_budget=100_000
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--cap", type=int, default=16)
    ap.add_argument("--code", default="bb_288_12_18")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    table = {
        r["code_id"]: r
        for r in yaml.safe_load(
            (LAB_ROOT / "instances" / "bravyi_table.yaml").read_text()
        )["instances"]
    }
    row = table[args.code]
    G = ZmZn(row["group"]["ell"], row["group"]["m"])
    A = Poly.from_string(row["polynomials"]["A"], G)
    B = Poly.from_string(row["polynomials"]["B"], G)
    checks = bb_check_matrices(A, B)
    action = compute_class_action(checks)
    k = action.k

    wcnf = out / f"naive_{checks.group.label()}.wcnf"
    qv, a_lits = write_wcnf(checks, wcnf, mode="naive")
    print(f"{args.code}: k={k}, wcnf={wcnf}", flush=True)

    for sigma in axis_decks(checks):
        t0 = time.perf_counter()
        dd = compute_descent(checks, action, sigma, A=A, B=B)
        lam_of = {
            c: _bits_to_int((_int_to_bits(c, k) @ dd.Lam) % 2)
            for c in range(1, 1 << k)
        }
        distinct = sorted(set(lam_of.values()))
        inv = invariant_floor(dd)
        print(
            f"σ={sigma}: kb={dd.kb} distinct-Λ={len(distinct)} "
            f"invFloor={inv} ({time.perf_counter()-t0:.1f}s struct)",
            flush=True,
        )
        key = str(sigma)
        _DD[key] = dd
        chunks = [
            (key, distinct[i:: args.jobs], args.cap)
            for i in range(args.jobs)
        ]
        t0 = time.perf_counter()
        with mp.get_context("fork").Pool(args.jobs) as pool:
            merged: dict[int, int] = {}
            for part in pool.map(_floor_chunk, chunks):
                merged.update(part)
        vals = sorted(set(merged.values()))
        print(
            f"σ={sigma}: floors done in {time.perf_counter()-t0:.1f}s, "
            f"values {vals}",
            flush=True,
        )
        sig = "".join(str(int(x)) for x in sigma)
        flb = out / f"fiber_{checks.group.label()}_{sig}.flb"
        with flb.open("w") as f:
            f.write(f"p fiberlb {k} {dd.S.shape[0]} {inv}\n")
            f.write("a " + " ".join(str(int(x)) for x in a_lits) + "\n")
            tbl = [0] * (1 << k)
            for c in range(1, 1 << k):
                tbl[c] = merged[lam_of[c]]
            f.write("f " + " ".join(map(str, tbl)) + "\n")
            for b in range(dd.S.shape[0]):
                x, y = (int(qv[i]) for i in np.flatnonzero(dd.S[b]))
                f.write(f"{x} {y}\n")
        print(f"wrote {flb}", flush=True)


if __name__ == "__main__":
    main()
