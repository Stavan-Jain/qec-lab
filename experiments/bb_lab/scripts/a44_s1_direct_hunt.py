#!/usr/bin/env python3
"""A44 S1b — direct SAT witness hunt on a member torus (existence lane).

usage: a44_s1_direct_hunt.py L M Wmax [time_s] [dense|any] [nsol]

Finds nontrivial X-logicals of weight <= Wmax (gap-dense on both axes
when `dense`), re-verifies each in the kernel, blocks translation orbits,
records Laurent forms, gap structure and period lattice.  Witness weights
are UPPER bounds; UNSAT / timeout are observations (A40 §16.5: the solver
cannot certify emptiness at these sizes).
Output data/a44/s1_direct_{L}x{M}_W{Wmax}.{json,log}.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a44_common import (  # noqa: E402
    LAB, DATA, Logger, member_code, hunt, pair_of, poly_str, picture, rss_guard, C,
)
from bb_lab.tower import validate_banked  # noqa: E402


def periods(code, cells):
    l, m = code.G.orders
    cs = frozenset(cells)
    return [(dx, dy) for dx in range(l) for dy in range(m) if (dx, dy) != (0, 0)
            and frozenset((b, (x + dx) % l, (y + dy) % m) for b, x, y in cs) == cs]


def main():
    L, M, W = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 600.0
    dense = (sys.argv[5] if len(sys.argv) > 5 else "dense") == "dense"
    nsol = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    tag = f"s1_direct_{L}x{M}_W{W}" + ("" if dense else "_any")
    log = Logger(DATA / f"{tag}.log")
    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    code = member_code(L, M)
    log(f"=== ({L},{M}) n={code.n} k={code.k}: nontrivial <= {W}, dense={dense}, "
        f"cap {tl} s, nsol {nsol}")
    sols, status, dt = hunt(code, W, nsol=nsol, time_limit=tl, tag=tag, dense=dense)
    out = dict(L=L, M=M, W=W, dense=dense, status=status, t=dt, witnesses=[])
    for s in sols:
        cells = [tuple(c) for c in s["cells"]]
        v = C.vec_of(code, cells)
        assert code.is_cycle(v) and not code.is_stab(v)
        f, g = pair_of(code, v)
        nx, gx = C.gap_structure([c[1] for c in cells], L)
        ny, gy = C.gap_structure([c[2] for c in cells], M)
        per = periods(code, cells)
        rec = dict(w=int(v.sum()), n_x=nx, n_y=ny, gap_x=gx, gap_y=gy, periods=per,
                   blocks=[len(f), len(g)], f=poly_str(f), g=poly_str(g), cells=sorted(cells))
        out["witnesses"].append(rec)
        log(f"  witness w={rec['w']} gaps=({gx},{gy}) n_x={nx} n_y={ny} blocks={rec['blocks']} "
            f"periods={per}\n    f = {rec['f']}\n    g = {rec['g']}")
        log(picture(code, cells))
    log(f"=== status={status} ({dt} s), {len(sols)} witnesses, RSS {rss_guard('end'):.2f} GB")
    (DATA / f"{tag}.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
