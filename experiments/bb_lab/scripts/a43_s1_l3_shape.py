#!/usr/bin/env python3
"""A43 S1 — the SHAPE of doubly-spanning minimizers (falsify-first for L3).

L3 (crossing lemma, A43 note §2): a nontrivial X-logical with no cyclic
x-gap >= 4 and no cyclic y-gap >= 4 ("gap-dense", S11's D sector) weighs
>= l + m - 6.  Before any proof attempt this script measures, on every
banked doubly-spanning minimizer (S11's populations, re-derived here):

  * occupancy: #occupied columns n_x, #occupied rows n_y, max gaps;
  * the per-column / per-row count histograms;
  * the DEFECT  n_x + n_y - w  (the number of qubits that would have to
    "pay twice" if every column and every row must be paid for once);
  * the cut states: for every vertical / horizontal cut the left-partial
    syndrome of the straddling Z-checks (dirty iff nonzero) and its weight;
  * block split, footprint components;

and prints ASCII pictures of representatives.  Positive controls: the
792 x-windowed gross minimizers (must show 2 per row and an x-gap >= 4),
the L12 stack at (18,12).

Populations (all re-verified end-to-end as nontrivial cycles):
  bb72 (6,6): full census W<=6 (84 objects, all 2D);
  two-gross (12,12): a36 witness + (12;4,4) shear pullbacks (49, all 2D);
  gross (12,6): full census W<=12 (1,884: 1,092 2D + 792 W_x).
Runtime ~2 min (the gross census dominates).  No SAT.

Output: data/a43/s1_l3_shape.json (+ .log).  Read-only on a40 data.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import TowerCode, v2i, validate_banked  # noqa: E402
import a40_s11_compare as C  # noqa: E402

DATA = LAB / "data" / "a43"
DATA.mkdir(parents=True, exist_ok=True)
LOG = DATA / "s1_l3_shape.log"
SPAN = C.SPAN  # 4: Lemma K's gap threshold (S11's gap-dense = no gap >= 4)


def log(msg: str):
    print(msg, flush=True)
    with LOG.open("a") as f:
        f.write(msg + "\n")


# ------------------------------------------------------------ geometry
def zcheck_cells(code: TowerCode, g):
    """Cells (blk, x, y) read by the Z-check at g (READ convention)."""
    l, m = code.G.orders
    return [(blk, (g[0] + s[0]) % l, (g[1] + s[1]) % m) for blk, s in C.READ]


def cut_states_simple(code: TowerCode, cells: frozenset, axis: int):
    """Same as cut_states but written plainly: unwrap the check position
    c so that the cut j sits inside its window [c+lo, c+hi]."""
    l, m = code.G.orders
    order = code.G.orders[axis]
    other_order = code.G.orders[1 - axis]
    offs = sorted({s[axis] for _, s in C.READ})
    lo, hi = offs[0], offs[-1]
    res = []
    for j in range(order):
        weight = 0
        # unwrapped check coordinate c_u with c_u + lo <= j < c_u + hi
        for c_u in range(j - hi + 1, j - lo + 1):
            c = c_u % order
            for q in range(other_order):
                g = (c, q) if axis == 0 else (q, c)
                par = 0
                for blk, s in C.READ:
                    cell = ((g[0] + s[0]) % l, (g[1] + s[1]) % m)
                    if c_u + s[axis] <= j and (blk, cell[0], cell[1]) in cells:
                        par ^= 1
                weight += par
        res.append((weight > 0, weight))
    return res


def picture(code: TowerCode, cells: frozenset) -> str:
    """ASCII: rows = y (top = m-1), cols = x; 'a' block 0, 'b' block 1,
    '*' both."""
    l, m = code.G.orders
    grid = [["." for _ in range(l)] for _ in range(m)]
    for blk, x, y in cells:
        ch = "a" if blk == 0 else "b"
        cur = grid[m - 1 - y][x]
        grid[m - 1 - y][x] = ch if cur == "." else "*"
    rows = ["".join(r) for r in grid]
    return "\n".join(f"y={m-1-i:2d} {r}" for i, r in enumerate(rows)) + \
        "\n     " + "".join(str(x % 10) for x in range(l))


def measure(code: TowerCode, v: np.ndarray, kind: str, sector: str,
            with_cuts=True) -> dict:
    l, m = code.G.orders
    cells = C.cells_of(code, v)
    xs = [c[1] for c in cells]
    ys = [c[2] for c in cells]
    n_x, gx = C.gap_structure(xs, l)
    n_y, gy = C.gap_structure(ys, m)
    colc = Counter(xs)
    rowc = Counter(ys)
    w = int(v.sum())
    rec = dict(w=w, kind=kind, sector=sector, n_x=n_x, n_y=n_y,
               gap_x=gx, gap_y=gy, defect=n_x + n_y - w,
               col_hist=dict(sorted(Counter(colc.values()).items())),
               row_hist=dict(sorted(Counter(rowc.values()).items())),
               blocks=[sum(1 for c in cells if c[0] == 0),
                       sum(1 for c in cells if c[0] == 1)],
               ncomp=len(C.components(code, cells)),
               cells=[list(c) for c in sorted(cells)])
    if with_cuts:
        cx = cut_states_simple(code, cells, 0)
        cy = cut_states_simple(code, cells, 1)
        rec["dirty_x"] = sum(1 for d, _ in cx if d)
        rec["dirty_y"] = sum(1 for d, _ in cy if d)
        rec["cutw_x"] = [wt for _, wt in cx]
        rec["cutw_y"] = [wt for _, wt in cy]
    return rec


def summarize(rows: list[dict], label: str) -> dict:
    out = {}
    by = {}
    for r in rows:
        by.setdefault((r["kind"], r["sector"], r["w"]), []).append(r)
    for key, rs in sorted(by.items(), key=lambda kv: str(kv[0])):
        kind, sector, w = key
        s = dict(n=len(rs),
                 defect_hist=dict(sorted(Counter(r["defect"] for r in rs).items())),
                 nx_hist=dict(sorted(Counter(r["n_x"] for r in rs).items())),
                 ny_hist=dict(sorted(Counter(r["n_y"] for r in rs).items())),
                 gapx_hist=dict(sorted(Counter(r["gap_x"] for r in rs).items())),
                 gapy_hist=dict(sorted(Counter(r["gap_y"] for r in rs).items())),
                 col_hist_union=dict(sorted(sum((Counter(
                     {int(k): v for k, v in r["col_hist"].items()})
                     for r in rs), Counter()).items())),
                 row_hist_union=dict(sorted(sum((Counter(
                     {int(k): v for k, v in r["row_hist"].items()})
                     for r in rs), Counter()).items())),
                 blocks_hist=dict(sorted(Counter(tuple(r["blocks"]) for r in rs).items())),
                 ncomp_hist=dict(sorted(Counter(r["ncomp"] for r in rs).items())))
        if "dirty_x" in rs[0]:
            s["dirty_x_hist"] = dict(sorted(Counter(r["dirty_x"] for r in rs).items()))
            s["dirty_y_hist"] = dict(sorted(Counter(r["dirty_y"] for r in rs).items()))
            s["cutw_x_range"] = [min(min(r["cutw_x"]) for r in rs),
                                 max(max(r["cutw_x"]) for r in rs)]
            s["cutw_y_range"] = [min(min(r["cutw_y"]) for r in rs),
                                 max(max(r["cutw_y"]) for r in rs)]
        s["blocks_hist"] = {str(k): v for k, v in s["blocks_hist"].items()}
        out[f"{kind}|{sector}|w{w}"] = s
        log(f"[{label}] {kind}|{sector}|w{w}: n={s['n']} defect={s['defect_hist']} "
            f"n_x={s['nx_hist']} n_y={s['ny_hist']} gap_x={s['gapx_hist']} "
            f"gap_y={s['gapy_hist']} cols={s['col_hist_union']} rows={s['row_hist_union']} "
            f"blocks={s['blocks_hist']} comps={s['ncomp_hist']}"
            + (f" dirty_x={s['dirty_x_hist']} dirty_y={s['dirty_y_hist']} "
               f"cutw_x∈{s['cutw_x_range']} cutw_y∈{s['cutw_y_range']}"
               if "dirty_x_hist" in s else ""))
    return out


def kinds_for(code: TowerCode, objs):
    """Class kind per object via S11's windowed-class decomposition
    (x-only = y-wrapping classes W_x^cl, y-only = x-wrapping W_y^cl,
    mixed = rank-2); 'n/a' when the decomposition is not direct."""
    dec = C.decomposition(code)
    if not dec["direct_sum"]:
        return ["n/a"] * len(objs), dec
    kinds = []
    for c, v in objs:
        cx, cy = C.class_split(c, dec["Wx_basis"], dec["Wy_basis"])
        kinds.append("x-only" if cy == 0 else "y-only" if cx == 0 else "mixed")
    return kinds, dec


def main():
    LOG.write_text(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')} a43_s1_l3_shape\n")
    t0 = time.time()
    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    binp = cosetbz.build_kernel()
    out: dict = {"span": SPAN, "frames": {}, "pictures": {}}

    def frame(name, code, objs, with_cuts=True, n_pics=2):
        kinds, dec = kinds_for(code, objs)
        rows = []
        pics = {}
        for (c, v), kind in zip(objs, kinds):
            assert code.is_cycle(v) and not code.is_stab(v)
            sec = C.classify_object(code, v, dirs=[(0, 1), (1, 0)])["sector"]
            rec = measure(code, v, kind, sec, with_cuts)
            rec["class"] = int(c)
            rows.append(rec)
            key = f"{kind}|{sec}|w{rec['w']}"
            if len(pics.get(key, [])) < n_pics:
                pics.setdefault(key, []).append(picture(code, C.cells_of(code, v)))
        log(f"=== {name}: {len(rows)} objects, k={code.k}, "
            f"decomposition={ {k: v for k, v in dec.items() if not k.endswith('basis')} }")
        summ = summarize(rows, name)
        for key, ps in pics.items():
            for p in ps:
                log(f"--- {name} picture {key}\n{p}")
        out["frames"][name] = dict(lm=list(code.G.orders), k=code.k,
                                   n_objects=len(rows), summary=summ,
                                   objects=rows)
        out["pictures"][name] = pics
        return rows

    # (6,6): full census W<=6
    bb72 = C.member_code(6, 6)
    C.assert_conventions(C.member_code(12, 12))
    objs = C.census_all_classes(binp, bb72, 6, "a43_s1_bb72")
    frame("bb72(6,6)", bb72, objs, n_pics=3)
    C.rss_guard("bb72")

    # (12,12): a36 witness + shear pullbacks (S11's population)
    tg = C.member_code(12, 12)
    wit = C.a36_witness(tg)
    objs = [(v2i(tg.sig(wit)), wit)]
    pb, _ = C.shear_pullbacks(binp, 12, 4, 4, 6, tg)
    seen = {v2i(wit)}
    for c, v in pb:
        if v2i(v) not in seen:
            seen.add(v2i(v))
            objs.append((v2i(tg.sig(v)), v))
    frame("two-gross(12,12)", tg, objs, n_pics=2)
    C.rss_guard("tg")

    # gross: full census W<=12
    gross = C.member_code(12, 6)
    objs = C.census_all_classes(binp, gross, 12, "a43_s1_gross")
    frame("gross(12,6)", gross, objs, n_pics=2)
    C.rss_guard("gross")

    # positive control: the L12 stack at (18,12) (windowed, weight 24)
    t21 = C.member_code(18, 12)
    v = C.l12_stack(t21)
    frame("tdg432(18,12) L12 stack", t21, [(v2i(t21.sig(v)), v)], n_pics=1)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s1_l3_shape.json").write_text(json.dumps(out, indent=1))
    log(f"wrote {DATA / 's1_l3_shape.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
