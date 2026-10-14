#!/usr/bin/env python3
"""A42 S4 — sheet diagrams of the floor-weight MIXED cycles (p = 3, 6
atlas): the equality cases of the hiding-mass inequality

    4|s cap S| - 3|s|  <=  2 (|S| - 3q)          (HM)

(equivalent to wt >= 6q given the pure floor |S| >= 3q).  For every
nontrivial atlas cycle we tabulate (|S|, |s|, |s cap S|, excess,
discount) and print one sheet diagram per profile: for each block and
column, the q fibres as triples (sheet0 sheet1 sheet2) = cells
(j, j+q, j+2q).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AT = _load("a40_s4_phase_atlas")
DATA = LAB / "data" / "a42"


def analyse(p, cap):
    q = p // 3
    rows, _ = AT.atlas("AB", p, cap, keep_pts=True)
    out = []
    seen_prof = {}
    for row in rows:
        if not row["nontrivial"]:
            continue
        fib = {}
        for (c, y, blk) in row["pts"]:
            key = (blk, c, y % q)
            fib.setdefault(key, set()).add((y // q) % 3)
        S = {k for k, v in fib.items() if len(v) in (1, 2)}
        s = {k for k, v in fib.items() if len(v) in (1, 3)}
        n = [0, 0, 0, 0]
        for v in fib.values():
            n[len(v)] += 1
        rec = dict(weight=row["weight"], profile=n[1:], S=len(S),
                   s=len(s), cap=len(S & s), excess=len(S) - 3 * q,
                   discount=4 * len(S & s) - 3 * len(s))
        rec["slack"] = 2 * rec["excess"] - rec["discount"]
        assert rec["slack"] == row["weight"] - 6 * q
        out.append(rec)
        prof = (row["weight"], tuple(n[1:]))
        if prof not in seen_prof:
            seen_prof[prof] = (row["pts"], fib)
    # summary
    print(f"p={p}: {len(out)} nontrivial cycles (cap {cap})")
    agg = {}
    for r in out:
        k = (r["weight"], tuple(r["profile"]), r["S"], r["s"], r["cap"],
             r["excess"], r["discount"], r["slack"])
        agg[k] = agg.get(k, 0) + 1
    print("  (w, (n1,n2,n3), |S|, |s|, |S&s|, excess, discount, "
          "slack=w-6q) : count")
    for k, v in sorted(agg.items()):
        print(f"    {k} : {v}")
    # diagrams
    for prof, (pts, fib) in sorted(seen_prof.items()):
        w, n = prof
        if n[0] == 0 and n[2] == 0:
            continue   # pure — the familiar 3-slot lift
        print(f"\n  diagram: weight {w} profile {n}  "
              f"(fibre triples = sheets (j, j+q, j+2q); '.' empty)")
        cols = sorted({c for (_, c, _) in fib})
        for blk in (0, 1):
            print(f"    block {blk}:")
            for c in range(cols[0], cols[-1] + 1):
                cells = []
                for j in range(q):
                    sh = fib.get((blk, c, j), set())
                    cells.append("".join("x" if i in sh else "."
                                         for i in range(3)))
                tag = " ".join(cells)
                print(f"      col {c:3d}: {tag}")
    return out


def main():
    import contextlib
    import io
    buf = io.StringIO()
    res = {}
    with contextlib.redirect_stdout(buf):
        for p, cap in ((3, 8), (6, 13)):
            res[f"p{p}"] = analyse(p, cap)
    text = buf.getvalue()
    print(text, end="")
    (DATA / "s4_sheets.log").write_text(text)
    (DATA / "s4_sheets.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
