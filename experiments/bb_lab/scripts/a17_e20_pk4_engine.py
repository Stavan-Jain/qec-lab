#!/usr/bin/env python3
"""A17 E20 — P-K4 driver: the (2,4) overlay enumeration.

Geometry (Theorem E18.2 + doc §8). A (2,4) match forces: u_R a K₄
in Cay(G, dB), |σ| = 8, δ = g₂ − g₁ ∈ dA with Sidon pair
(a_hi, a_lo), part1 = (A ∖ {a_hi}) + g₁ and part2 = (A ∖ {a_lo}) + g₂
each a transversal of the four B-translates (one cell per translate,
2 per translate total), all six collision cells distinct and off σ.

Gauges: a₀ = b₀ = 0 (set translations), g₁ = 0 (joint cycle
translation), (hi, lo) = (4, 3) by a-relabeling, translates ordered
by their part1 label (π₁(j) = j−1), (b₁¹, b₁²) = (b₀, b₁) by
b-relabeling, first-use canonicalization on b-labels {2,3,4}.
Eliminations: t_j = a_{j−1} − b_j¹, g₂ = a₄ − a₃ — everything lives
in the same 8 coordinates (a₁..a₄, b₁..b₄) as the P-33 driver.

Branching: π₂ (bijection onto {0,1,2,4}) × ordered witness pairs
(b_j¹, b_j²) per translate × 6 K₄-edge witnesses (ordered dB-pairs),
with relations
    R_j:  a_{π₂(j)} + g₂ − b_j² − t_j = 0,
    E_jk: t_j − t_k − (b_s − b_u) = 0,
interleaved with an incrementally grown battery (static Sidon/D2 +
t-distinctness + m_A + m_B forms as translates appear). Terminals
are exhaustively exact-verified INCLUDING the witness-dependent
collision-cell forms (distinct ×15, off-σ ×48, m_B(v) = 2 ×60).

Soundness identical to the P-33 driver: sound d ∈ {1,2} kills only,
every kill exact-confirmed, every terminal exact-re-verified.
Expected terminals include the sub-±B families (u_R inside a
±B-translate), which are PROVEN dead frame-free by the §8
fifth-translate / m = 4 endgames — the classifier separates them.

Usage:
  uv run python scripts/a17_e20_pk4_engine.py --mode selftest
  uv run python scripts/a17_e20_pk4_engine.py --mode pk4 \
      --out data/a17/e20_pk4_table.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a17_e20_overlay_engine import (P, Battery, ZLattice, avec,
                                    bvec, killable, static_battery,
                                    vadd, vsub)

G2 = vsub(avec(4), avec(3))                      # g₂ = a₄ − a₃
AEDGE = [(p, q) for p in range(5) for q in range(5) if p != q]
PART2 = (0, 1, 2, 4)


def tvec(j: int, b1: int) -> tuple:
    """t_j for translate j (1-indexed) with part1 witness b_j¹."""
    return vsub(avec(j - 1), bvec(b1))


def cell_rel(j: int, l2lab: int, b1: int, b2: int) -> tuple:
    """R_j: a_{π₂(j)} + g₂ − b_j² − (a_{j−1} − b_j¹) = 0."""
    return vsub(vsub(vadd(avec(l2lab), G2), bvec(b2)),
                tvec(j, b1))


def mA_forms(j: int, l2lab: int) -> tuple:
    forms, kinds = [], []
    y1 = avec(j - 1)
    y2 = vadd(avec(l2lab), G2)
    for p in range(5):
        forms.append(vsub(vsub(y1, G2), avec(p)))   # part1 ∉ A+g₂
        kinds.append("eq")
        forms.append(vsub(y2, avec(p)))             # part2 ∉ A+g₁
        kinds.append("eq")
    return forms, kinds


def mB_forms(j: int, k: int, pi2, bp) -> tuple:
    """Cells of translates j vs translate k (both directions)."""
    forms, kinds = [], []
    for (cj, ck) in ((j, k), (k, j)):
        tv = tvec(ck, bp[ck][0])
        for y in (avec(cj - 1), vadd(avec(pi2[cj]), G2)):
            for p in range(5):
                forms.append(vsub(vsub(y, tv), bvec(p)))
                kinds.append("eq")
    return forms, kinds


def collision_forms(bp: dict, ewits: dict) -> tuple:
    """v_jk = b_u + t_j for edge witness b_s − b_u = t_j − t_k:
    pairwise distinct, off σ, in no third translate."""
    forms, kinds = [], []
    cells = {}
    for (j, k), (s, u) in ewits.items():
        cells[(j, k)] = vadd(bvec(u), tvec(j, bp[j][0]))
    vc = list(cells.values())
    for a, b in combinations(vc, 2):
        forms.append(vsub(a, b)); kinds.append("eq")
    sigma = [avec(i) for i in range(4)]
    sigma += [vadd(avec(l), G2) for l in collision_forms.pi2v]
    for v in vc:
        for y in sigma:
            forms.append(vsub(v, y)); kinds.append("eq")
    for (j, k), v in cells.items():
        for m in range(1, 5):
            if m in (j, k):
                continue
            tv = tvec(m, bp[m][0])
            for p in range(5):
                forms.append(vsub(vsub(v, tv), bvec(p)))
                kinds.append("eq")
    return forms, kinds


def pk4_search(max_nodes: int, progress: bool = True,
               relax: bool = False) -> dict:
    t0 = time.time()
    stat = static_battery()
    stats = {"nodes": 0, "dead": 0, "terminals": 0, "capped": False}
    terminals = []

    def bump() -> bool:
        stats["nodes"] += 1
        if progress and stats["nodes"] % 1000000 == 0:
            print(f"  ..{stats['nodes']} nodes dead={stats['dead']}"
                  f" term={stats['terminals']}"
                  f" {round(time.time() - t0)}s", flush=True)
        if stats["nodes"] >= max_nodes:
            stats["capped"] = True
        return not stats["capped"]

    EDGES = {2: [(1, 2)], 3: [(1, 3), (2, 3)],
             4: [(1, 4), (2, 4), (3, 4)]}

    def terminal(lat, forms, kinds, pi2, bp, ewits):
        collision_forms.pi2v = ([] if relax else
                                [pi2[j] for j in (1, 2, 3, 4)])
        cf, ck = collision_forms(bp, ewits)
        full = Battery(forms + cf, kinds + ck)
        if full.exact_all(lat):
            stats["dead"] += 1
            return
        stats["terminals"] += 1
        terminals.append({
            "pi2": [pi2[j] for j in (1, 2, 3, 4)],
            "bpairs": [list(bp[j]) for j in (1, 2, 3, 4)],
            "ewits": {f"{j},{k}": list(w)
                      for (j, k), w in ewits.items()},
            "rank": lat.rank(),
            "lattice": [r[:] for r in lat.rows]})

    def edge_layer(eidx, elist, lat, bat, RF, dset, pi2, bp, ewits,
                   after):
        if eidx == len(elist):
            after(lat, bat, RF, dset)
            return
        j, k = elist[eidx]
        target = vsub(tvec(j, bp[j][0]), tvec(k, bp[k][0]))
        base_rows, base_pivs = lat.basis_modp()
        for (s, u) in AEDGE:
            if not bump():
                return
            rel = vsub(target, vsub(bvec(s), bvec(u)))
            l2 = lat.clone()
            grew = l2.add(rel)
            if not grew:
                ewits[(j, k)] = (s, u)
                edge_layer(eidx + 1, elist, l2, bat, RF, dset,
                           pi2, bp, ewits, after)
                del ewits[(j, k)]
                if stats["capped"]:
                    return
                continue
            dead = False
            for idx in dset:
                if killable(l2.forced_denom(bat.forms[idx]),
                            bat.kinds[idx]):
                    dead = True
                    break
            if dead:
                stats["dead"] += 1
                continue
            rr = np.array(rel, dtype=np.int64) % P
            for br, bpv in zip(base_rows, base_pivs):
                rr = (rr - rr[bpv] * br) % P
            if not rr.any():
                nstate = (RF, dset)
            else:
                j0 = int(np.nonzero(rr)[0][0])
                rr = (rr * pow(int(rr[j0]), P - 2, P)) % P
                RF2 = (RF - np.outer(RF[:, j0], rr)) % P
                newly = np.nonzero(RF.any(axis=1)
                                   & ~RF2.any(axis=1))[0]
                d2 = list(dset)
                for idx in newly:
                    d = l2.forced_denom(bat.forms[idx])
                    if killable(d, bat.kinds[idx]):
                        dead = True
                        break
                    if d is not None:
                        d2.append(int(idx))
                if dead:
                    stats["dead"] += 1
                    continue
                nstate = (RF2, d2)
            ewits[(j, k)] = (s, u)
            edge_layer(eidx + 1, elist, l2, bat, nstate[0],
                       nstate[1], pi2, bp, ewits, after)
            del ewits[(j, k)]
            if stats["capped"]:
                return

    def extend_battery(lat, bat, newf, newk):
        """Grow the battery; None if a new form is already forced."""
        forms = bat.forms + newf
        kinds = bat.kinds + newk
        b2 = Battery(forms, kinds)
        if b2.first_kill(lat) is not None:
            return None
        RF = b2.reduced(lat)
        dset = [int(i) for i in np.nonzero(~RF.any(axis=1))[0]]
        return b2, RF, dset

    def level_j(j, lat, bat, pi2, bp, used2, used_b, ewits):
        if j == 5:
            terminal(lat, bat.forms, bat.kinds, pi2, bp, ewits)
            return
        labs = [l for l in PART2 if l not in used2]
        if relax:
            labs = labs[:1]
        for lab in labs:
            for (c, d) in AEDGE:
                # first-use canonicalization on b-labels {2,3,4}:
                # each NEW high label must be the smallest unused
                seen = set(used_b)
                ok = True
                for x in (c, d):
                    if x >= 2 and x not in seen:
                        if x != min({2, 3, 4} - seen):
                            ok = False
                            break
                        seen.add(x)
                if not ok:
                    continue
                if not bump():
                    return
                l2 = lat.clone()
                if not relax:
                    l2.add(cell_rel(j, lab, c, d))
                    if stat.first_kill(l2) is not None:
                        stats["dead"] += 1
                        continue
                nf, nk = mA_forms(j, lab)
                if relax:
                    nf = [f for i, f in enumerate(nf) if i % 2 == 0]
                    nk = nk[:len(nf)]
                tf = []
                tk = []
                for k in range(1, j):
                    tf.append(vsub(tvec(j, c), tvec(k, bp[k][0])))
                    tk.append("diff")
                    mf, mk = mB_forms(j, k, {**pi2, j: lab},
                                      {**bp, j: (c, d)})
                    if relax:      # keep only part1-cell (y1) forms
                        keep = [i for i in range(len(mf))
                                if (i // 5) % 2 == 0]
                        mf = [mf[i] for i in keep]
                        mk = [mk[i] for i in keep]
                    tf += mf
                    tk += mk
                ext = extend_battery(l2, bat, nf + tf, nk + tk)
                if ext is None:
                    stats["dead"] += 1
                    continue
                b2, RF, dset = ext
                pi2[j] = lab
                bp[j] = (c, d)
                edge_layer(
                    0, EDGES.get(j, []), l2, b2, RF, dset, pi2, bp,
                    ewits,
                    lambda l_, bt_, R_, D_: level_j(
                        j + 1, l_, bt_, pi2, bp,
                        used2 | {lab}, used_b | {x for x in (c, d)
                                                 if x >= 2},
                        ewits))
                del pi2[j]
                del bp[j]
                if stats["capped"]:
                    return

    # level 1: π₂(1) with gauge-fixed (b₁¹, b₁²) = (0, 1)
    for lab1 in (PART2[:1] if relax else PART2):
        if not bump():
            break
        lat = ZLattice()
        if not relax:
            lat.add(cell_rel(1, lab1, 0, 1))
            if stat.first_kill(lat) is not None:
                stats["dead"] += 1
                continue
        nf, nk = mA_forms(1, lab1)
        if relax:
            nf = [f for i, f in enumerate(nf) if i % 2 == 0]
            nk = nk[:len(nf)]
        ext = extend_battery(lat, stat, nf, nk)
        if ext is None:
            stats["dead"] += 1
            continue
        b2, RF, dset = ext
        pi2 = {1: lab1}
        bp = {1: (0, 1)}
        level_j(2, lat, b2, pi2, bp, {lab1}, set(), {})
        if stats["capped"]:
            break

    stats["secs"] = round(time.time() - t0, 1)
    stats["terminal_patterns"] = terminals
    return stats


def classify_subB(terminals: list) -> dict:
    """Sub-±B test: is t_j − t_k lattice-forced equal to
    ±(b_{f(j)} − b_{f(k)}) for an injective f?"""
    out = {"sub_B": 0, "sub_negB": 0, "other": []}
    for t in terminals:
        lat = ZLattice()
        lat.rows = [r[:] for r in t["lattice"]]
        bp = {j + 1: tuple(t["bpairs"][j]) for j in range(4)}
        tv = [tvec(j, bp[j][0]) for j in (1, 2, 3, 4)]
        found = None
        for sgn in (+1, -1):
            for f in permutations(range(5), 4):
                ok = True
                for (x, y) in combinations(range(4), 2):
                    diff = vsub(tv[x], tv[y])
                    bd = vsub(bvec(f[x]), bvec(f[y]))
                    if sgn < 0:
                        bd = tuple(-z for z in bd)
                    if lat.forced_denom(vsub(diff, bd)) != 1:
                        ok = False
                        break
                if ok:
                    found = sgn
                    break
            if found:
                break
        if found == 1:
            out["sub_B"] += 1
        elif found == -1:
            out["sub_negB"] += 1
        else:
            out["other"].append(t)
    return out


def selftest() -> None:
    # R_j consistency: a random integer solution of R_j must make
    # the part2 cell equation hold identically
    import random
    rng = random.Random(7)
    for _ in range(20):
        x = [rng.randrange(-50, 50) for _ in range(8)]
        av = lambda i: 0 if i == 0 else x[i - 1]
        bv = lambda i: 0 if i == 0 else x[4 + i - 1]
        j, lab, c, d = 2, 4, 2, 3
        rel = cell_rel(j, lab, c, d)
        val = sum(cf * xv for cf, xv in zip(rel, x))
        g2 = av(4) - av(3)
        tj = av(j - 1) - bv(c)
        lhs = av(lab) + g2
        rhs = bv(d) + tj
        assert (lhs - rhs) == val, "R_j encodes the cell equation"
    print("selftest OK")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["selftest", "pk4", "pk4-relaxed"])
    ap.add_argument("--max-nodes", type=int, default=200_000_000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.mode == "selftest":
        selftest()
        return
    res = pk4_search(args.max_nodes,
                     relax=args.mode == 'pk4-relaxed')
    summary = {k: v for k, v in res.items()
               if k != "terminal_patterns"}
    print(json.dumps(summary, indent=1))
    cls = classify_subB(res["terminal_patterns"])
    print(f"sub_B={cls['sub_B']} sub_negB={cls['sub_negB']} "
          f"other={len(cls['other'])}")
    if args.out:
        res["classify"] = {"sub_B": cls["sub_B"],
                           "sub_negB": cls["sub_negB"],
                           "n_other": len(cls["other"])}
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, indent=1))
    ok = not res["capped"]
    print(f"\nP-K4 table {'COMPLETE' if ok else 'CAPPED'}: "
          f"terminals={res['terminals']}, nodes={res['nodes']}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
