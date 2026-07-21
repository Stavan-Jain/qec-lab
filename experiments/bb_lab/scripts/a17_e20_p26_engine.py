#!/usr/bin/env python3
"""A17 E20 — P-26 driver: the (2,6) overlay enumeration.

Geometry (doc §9 profile theorem). For |u_L| = 2, |u_R| = 6 the
viable profiles (a, |σ|, n₁, n₂, n₃, n₄) are exactly
    (10,10,10,10,0,0)  (11,8,8,11,0,0)  (12,10,9,9,1,0)
    (13,8,7,10,1,0)    (14,10,10,8,0,1) (15,8,6,9,2,0)
    (15,8,8,9,0,1)
with n₅ = n₆ = 0 (mass 30 = n₁+2n₂+3n₃+4n₄, a = n₂+3n₃+6n₄,
|σ| = n₁+n₃, C2 cap n₁+3n₃ ≤ 12, a = 15 − |σ|/2 + 2n₃ + 4n₄).
Since Sidon-B gives every translate pair ≤ 1 shared cell, the a
collision pairs PARTITION into edge-disjoint cliques of shared
cells: n₂ K₂ (m = 2, off σ), n₃ K₃ (m = 3, ON σ), n₄ K₄ (m = 4,
off σ); per translate d_j = 5 − c₂ⱼ − c₄ⱼ ≤ 2 (C2), its σ-cells
in DISTINCT A-parts.

Coordinates: 9 = (a₁..a₄, b₁..b₄, g₂); gauges a₀ = b₀ = 0, g₁ = 0.
|σ| = 8 stratum: δ = g₂ = a₄ − a₃ ∈ dA (overlap gauge (hi,lo) =
(4,3), parts P₁ = {0,1,2,3}, P₂ = {0,1,2,4}; the overlap cell a₄
may host one even collision cell — off-σ forms exclude only the 8
σ-names). |σ| = 10 stratum: g₂ is a free coordinate, δ ∉ dA (21
"eq" ∉-forms — d = 2 is NOT a kill for δ), parts P₁ = P₂ = {0..4}.

Branching per (profile, clique-partition canon): translates in a
connectivity order (each anchored when processed); per translate
its σ-slot names + b-labels, then clique-cell b-labels; each
relation either ELIMINATES r_j or adds an integer lattice row;
battery grown incrementally: static Sidon/D2 + stratum δ-forms +
m = 1 sharpness (25/cell) + off-σ (s/cell) + outside-translate +
non-edge (21/pair) + collision-cell distinctness. Sound kills
d ∈ {1,2} only, every kill exact-confirmed, every terminal
exact-re-verified — the P-33/P-K4 soundness architecture verbatim.

First-use canonicalization: b-labels globally; a-labels within the
symmetric set ({0,1,2} at s = 8; all five at s = 10); at s = 10
additionally part 1 before part 2 (spends the g₁ ↔ g₂ swap).
Translate labels are spent by partition canon (Aut(partition)
overcount accepted — soundness unaffected). Every abstract (2,6)
match lands in this frame: the per-profile tables are COMPLETE.

Validation: --mode p26-relaxed drops the whole A-side (collision
half-system, gauge r₁ = 0) — realizable near-cliques exist in the
census, so terminals MUST appear. --mode p26-inject maps censused
candidate 6-sets from live members into the abstract collision
half, asserts every relation concretely in the group, and requires
the engine not to kill them.

Usage:
  uv run python scripts/a17_e20_p26_engine.py --mode selftest
  uv run python scripts/a17_e20_p26_engine.py --mode partitions
  uv run python scripts/a17_e20_p26_engine.py --mode p26 \
      --profile a12 --out data/a17/e20_p26_a12.json
  uv run python scripts/a17_e20_p26_engine.py --mode p26-relaxed \
      --profile a10
  uv run python scripts/a17_e20_p26_engine.py --mode p26-inject \
      --profile a10 --members data/a17/members_7x9.jsonl --cap 300
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

NV9 = 9
ZERO9 = tuple([0] * NV9)
G2 = tuple([0] * 8 + [1])


def pad(f) -> tuple:
    return tuple(f) + (0,)


def a9(i: int) -> tuple:
    return pad(avec(i))


def b9(i: int) -> tuple:
    return pad(bvec(i))


def static9() -> Battery:
    s = static_battery()
    return Battery([pad(f) for f in s.forms], list(s.kinds))


# name (i, p) -> cell element a_p + g_i  (g_1 = 0, g_2 = G2 coord)
def name_vec(i: int, p: int) -> tuple:
    return a9(p) if i == 1 else vadd(a9(p), G2)


PROFILES = {
    "a10":  (10, 10, 10, 10, 0, 0),
    "a11":  (11, 8, 8, 11, 0, 0),
    "a12":  (12, 10, 9, 9, 1, 0),
    "a13":  (13, 8, 7, 10, 1, 0),
    "a14":  (14, 10, 10, 8, 0, 1),
    "a15t": (15, 8, 6, 9, 2, 0),
    "a15q": (15, 8, 8, 9, 0, 1),
}


def part_names(s: int) -> list:
    """σ-cell names (i, p) for the stratum."""
    if s == 8:
        return ([(1, p) for p in (0, 1, 2, 3)]
                + [(2, p) for p in (0, 1, 2, 4)])
    return [(i, p) for i in (1, 2) for p in range(5)]


def stratum_forms(s: int) -> tuple:
    """δ-forms: s = 10 demands δ = g₂ ∉ dA ∪ {0} (21 "eq" forms —
    2-torsion δ is not excluded, so d = 2 must not kill).
    s = 8 needs none (g₂ = a₄ − a₃ is a lattice row; Sidon-A covers
    cross-part coincidences)."""
    if s == 8:
        return [], []
    forms = [G2]
    for p in range(5):
        for q in range(5):
            if p != q:
                forms.append(vsub(G2, vsub(a9(p), a9(q))))
    return forms, ["eq"] * len(forms)


def enum_partitions(prof: tuple) -> list:
    """Clique partitions up to S₆: n₄ K₄ + n₃ K₃ + n₂ K₂,
    edge-disjoint, per-vertex c₂+c₄ ≥ 3 and c₂+c₃+c₄ ≤ 5.
    Returns relabeled structures with a valid processing order
    1..6 (each translate σ-anchored or clique-adjacent to its
    prefix)."""
    a, s, n1, n2, n3, n4 = prof
    edges_of = {4: lambda q: set(combinations(sorted(q), 2)),
                3: lambda t: set(combinations(sorted(t), 2)),
                2: lambda e: {tuple(sorted(e))}}
    seen = set()
    raw = []

    def canon_key(k4s, k3s, k2s):
        best = None
        for pi in permutations(range(6)):
            key = (tuple(sorted(tuple(sorted(pi[v] for v in q))
                                for q in k4s)),
                   tuple(sorted(tuple(sorted(pi[v] for v in t))
                                for t in k3s)),
                   tuple(sorted(tuple(sorted((pi[u], pi[w])))
                                for (u, w) in k2s)))
            if best is None or key < best:
                best = key
        return best

    k4opts = list(combinations(range(6), 4)) if n4 else [None]
    for q4 in k4opts:
        used4 = edges_of[4](q4) if q4 else set()
        tri_pool = [t for t in combinations(range(6), 3)
                    if not (edges_of[3](t) & used4)]
        k3opts = ([()] if n3 == 0 else
                  [(t,) for t in tri_pool] if n3 == 1 else
                  [(t1, t2) for t1, t2 in combinations(tri_pool, 2)
                   if not (edges_of[3](t1) & edges_of[3](t2))])
        for tris in k3opts:
            used = set(used4)
            for t in tris:
                used |= edges_of[3](t)
            rest = [e for e in combinations(range(6), 2)
                    if e not in used]
            for e2s in combinations(rest, n2):
                c2 = [0] * 6
                c3 = [0] * 6
                c4 = [0] * 6
                for (u, w) in e2s:
                    c2[u] += 1
                    c2[w] += 1
                for t in tris:
                    for v in t:
                        c3[v] += 1
                if q4:
                    for v in q4:
                        c4[v] += 1
                if any(c2[v] + c4[v] < 3 for v in range(6)):
                    continue
                if any(c2[v] + c3[v] + c4[v] > 5 for v in range(6)):
                    continue
                key = canon_key([q4] if q4 else [], list(tris),
                                list(e2s))
                if key in seen:
                    continue
                seen.add(key)
                raw.append((([q4] if q4 else []), list(tris),
                            list(e2s), c2, c3, c4))

    out = []
    for k4s, k3s, k2s, c2, c3, c4 in raw:
        d = [5 - c2[v] - c4[v] for v in range(6)]
        c1 = [d[v] - c3[v] for v in range(6)]
        if any(x < 0 for x in c1):
            continue
        # processing order: greedy, always anchorable
        order = []
        rem = set(range(6))
        cliques = ([("K4", frozenset(q)) for q in k4s]
                   + [("K3", frozenset(t)) for t in k3s]
                   + [("K2", frozenset(e)) for e in k2s])

        def adj_prefix(v, pref):
            return any(v in mem and (mem & pref)
                       for _, mem in cliques)

        while rem:
            cands = [v for v in rem
                     if d[v] > 0 or adj_prefix(v, set(order))]
            if not cands:
                break
            cands.sort(key=lambda v: (
                -sum(1 for _, mem in cliques
                     if v in mem and (mem & set(order))),
                -d[v], v))
            order.append(cands[0])
            rem.discard(cands[0])
        if rem:
            continue                     # unanchorable (impossible)
        relab = {v: i + 1 for i, v in enumerate(order)}
        out.append({
            "k4": [sorted(relab[v] for v in q) for q in k4s],
            "k3": [sorted(relab[v] for v in t) for t in k3s],
            "k2": [sorted(relab[v] for v in e) for e in k2s],
            "c1": {relab[v]: c1[v] for v in range(6)},
            "d":  {relab[v]: d[v] for v in range(6)},
        })
    return out


def partition_meta(part: dict, relax: bool = False) -> dict:
    """Cliques with ids, ordered step list, non-edges. Full mode
    processes translates 1..6 (enum order: σ-anchorable); relax
    mode has no σ-cells, so it uses a BFS order over clique
    adjacency (root = translate in most cliques, gauged r = 0)."""
    cliques = []
    for q in part["k4"]:
        cliques.append({"mem": tuple(q), "size": 4, "sigma": False})
    for t in part["k3"]:
        cliques.append({"mem": tuple(t), "size": 3,
                        "sigma": not relax})
    for e in part["k2"]:
        cliques.append({"mem": tuple(e), "size": 2, "sigma": False})
    edge_pairs = set()
    for c in cliques:
        for u, w in combinations(c["mem"], 2):
            edge_pairs.add((u, w))
    non_edges = [e for e in combinations(range(1, 7), 2)
                 if e not in edge_pairs]
    if relax:
        deg = {j: sum(1 for c in cliques if j in c["mem"])
               for j in range(1, 7)}
        order = [max(deg, key=lambda j: (deg[j], -j))]
        while len(order) < 6:
            nxt = [j for j in range(1, 7) if j not in order
                   and any(j in c["mem"] and set(c["mem"])
                           & set(order) for c in cliques)]
            assert nxt, "clique graph disconnected (impossible)"
            nxt.sort(key=lambda j: (-deg[j], j))
            order.append(nxt[0])
    else:
        order = list(range(1, 7))
    pos = {j: k for k, j in enumerate(order)}
    steps = []
    for j in order:
        js = []
        if not relax:
            for _ in range(part["c1"][j]):
                js.append(("sigma1", j, None))
            for ci, c in enumerate(cliques):
                if j in c["mem"] and c["sigma"]:
                    js.append(("k3", j, ci))
        for ci, c in enumerate(cliques):
            if j not in c["mem"] or c["sigma"]:
                continue
            earlier = [k for k in c["mem"] if pos[k] < pos[j]]
            js.append(("join" if earlier else "cell", j, ci))
        # joins before fresh cells (a fresh cell needs r_j known)
        js.sort(key=lambda st: {"sigma1": 0, "k3": 1, "join": 2,
                                "cell": 3}[st[0]])
        steps += js
    return {"cliques": cliques, "non_edges": non_edges,
            "steps": steps, "order": order, "pos": pos}


def p26_search(pname: str, max_nodes: int, relax: bool = False,
               progress: bool = True) -> dict:
    prof = PROFILES[pname]
    _, s, n1, n2, n3, n4 = prof
    t0 = time.time()
    stat = static9()
    NAMES = part_names(s)
    NAMEV = {nm: name_vec(*nm) for nm in NAMES}
    ASYM = set(range(5)) if s == 10 else {0, 1, 2}
    BFULL = set(range(5))
    stats = {"profile": pname, "relax": relax, "nodes": 0,
             "dead": 0, "terminals": 0, "capped": False,
             "partitions": 0}
    terminals = []

    def bump() -> bool:
        stats["nodes"] += 1
        if progress and stats["nodes"] % 1_000_000 == 0:
            print(f"  ..{stats['nodes']} nodes dead={stats['dead']}"
                  f" term={stats['terminals']}"
                  f" {round(time.time() - t0)}s", flush=True)
        if stats["nodes"] >= max_nodes:
            stats["capped"] = True
        return not stats["capped"]

    def canon_a(p, seen):
        return p in seen or (p in ASYM and p == min(ASYM - seen))

    def canon_b(u, seen):
        return u in seen or u == min(BFULL - seen)

    MODS = (P, 2, 3, 5)

    def cand_mask(karr, FB):
        """Vectorized killability-necessity mask: a form can be
        forced with d ∈ {1,2} only if its reduced row vanishes
        mod p, 3, 5 — and mod 2 as well for "eq" kinds (d = 1
        only). Unconditional (no rank caveat): gcd(d, m) = 1."""
        m = (~FB[P].any(axis=1) & ~FB[3].any(axis=1)
             & ~FB[5].any(axis=1))
        return m & (karr | ~FB[2].any(axis=1))

    def push(lat, S, inspan, rows, nf, nk):
        """Apply rows + battery extension; None iff DEAD (every
        kill exact-confirmed). S = (forms, kinds, karr, FB, bas):
        FB[m] = battery reduced by the running mod-m row basis
        bas[m]; both updated INCREMENTALLY (one outer product per
        Z-growing row / one small block per extension). inspan:
        idx -> last exact denom for exact-checked forms; on growth
        only denominators ≤ 8 are re-checked — later drops are
        caught by the terminal mask + exact verification, the
        soundness floor."""
        forms, kinds, karr, FB, bas = S
        l2 = lat
        newrows = []
        if rows:
            l2 = lat.clone()
            for r in rows:
                if any(r) and l2.add(list(r)):
                    newrows.append(r)
        if newrows:
            for idx, dv in inspan.items():
                if dv <= 8 and killable(
                        l2.forced_denom(forms[idx]), kinds[idx]):
                    return None
        if not newrows and not nf:
            return l2, S, inspan
        FB2 = dict(FB)
        bas2 = dict(bas)
        for r in newrows:
            for m in MODS:
                rr = np.array(r, dtype=np.int64) % m
                for br, bp in bas2[m]:
                    if rr[bp]:
                        rr = (rr - rr[bp] * br) % m
                nzi = np.nonzero(rr)[0]
                if len(nzi) == 0:
                    continue
                p0 = int(nzi[0])
                rr = (rr * pow(int(rr[p0]), -1, m)) % m
                FB2[m] = (FB2[m]
                          - np.outer(FB2[m][:, p0], rr)) % m
                bas2[m] = bas2[m] + ((rr, p0),)
        if nf:
            forms = forms + nf
            kinds = kinds + nk
            karr = np.concatenate(
                [karr, np.array([k == "diff" for k in nk])])
            add = np.array(nf, dtype=np.int64)
            for m in MODS:
                Aq = add % m
                for br, bp in bas2[m]:
                    Aq = (Aq - np.outer(Aq[:, bp], br)) % m
                FB2[m] = np.vstack([FB2[m], Aq])
        S2 = (forms, kinds, karr, FB2, bas2)
        cand = cand_mask(karr, FB2)
        if not cand.any():
            return l2, S2, inspan
        ins2 = None
        for idx in np.nonzero(cand)[0]:
            idx = int(idx)
            if idx in inspan:
                continue
            d = l2.forced_denom(forms[idx])
            if killable(d, kinds[idx]):
                return None
            if ins2 is None:
                ins2 = dict(inspan)
            ins2[idx] = d if d is not None else 1 << 30
        return l2, S2, (inspan if ins2 is None else ins2)

    def cell_outside_forms(v, mem, rex):
        """v lies in NO translate outside its membership set."""
        nf = []
        for k2, rx in rex.items():
            if k2 in mem:
                continue
            for t in range(5):
                nf.append(vsub(vsub(v, rx), b9(t)))
        return nf, ["eq"] * len(nf)

    def anchor_forms(j, rex, cells, non_edges):
        """Forms due when r_j first becomes known: non-edge
        exclusions vs anchored translates + membership exclusions
        vs every existing cell not containing j."""
        nf = []
        for (u_, w_) in non_edges:
            if j not in (u_, w_):
                continue
            k2 = w_ if u_ == j else u_
            if k2 not in rex or k2 == j:
                continue
            d0 = vsub(rex[j], rex[k2])
            nf.append(d0)
            for x in range(5):
                for y_ in range(5):
                    if x != y_:
                        nf.append(vsub(d0, vsub(b9(x), b9(y_))))
        for (v, mem, _isc) in cells:
            if j in mem:
                continue
            for t in range(5):
                nf.append(vsub(vsub(v, rex[j]), b9(t)))
        return nf, ["eq"] * len(nf)

    def run_partition(pidx, part):
        meta = partition_meta(part, relax=relax)
        cliques = meta["cliques"]
        steps = meta["steps"]
        sfrm, skin = ([], []) if relax else stratum_forms(s)
        forms0 = stat.forms + sfrm
        kinds0 = stat.kinds + skin
        karr0 = np.array([kk == "diff" for kk in kinds0])
        F0 = np.array(forms0, dtype=np.int64)
        S0 = (forms0, kinds0, karr0,
              {m: F0 % m for m in MODS}, {m: () for m in MODS})

        def record(lat, st):
            terminals.append({
                "profile": pname, "partition": pidx,
                "k4": part["k4"], "k3": part["k3"],
                "k2": part["k2"],
                "assign": [list(x) for x in st["assign"]],
                "rank": lat.rank(),
                "lattice": [r[:] for r in lat.rows]})

        def emit(k, lat, S, inspan, st, rows, nf, nk, rec):
            res = push(lat, S, inspan, rows, nf, nk)
            if res is None:
                stats["dead"] += 1
                return
            l2, S2, i2 = res
            walk(k + 1, l2, S2, i2, rec(st))

        def walk(k, lat, S, inspan, st):
            if stats["capped"]:
                return
            if k == len(steps):
                forms, kinds, karr, FB, _bas = S
                dead = False
                for idx in np.nonzero(cand_mask(karr, FB))[0]:
                    idx = int(idx)
                    if killable(lat.forced_denom(forms[idx]),
                                kinds[idx]):
                        dead = True
                        break
                if dead:
                    stats["dead"] += 1
                else:
                    stats["terminals"] += 1
                    record(lat, st)
                return
            kind, j, ci = steps[k]
            rex = st["rex"]
            if kind in ("sigma1", "k3"):
                mem = (frozenset({j}) if kind == "sigma1"
                       else frozenset(cliques[ci]["mem"]))
                named = kind == "k3" and ci in st["cellname"]
                if named:
                    nm0 = st["cellname"][ci]
                    if nm0[0] in st["pused"][j]:
                        stats["dead"] += 1
                        return
                    cand_names = [nm0]
                else:
                    pinned = {st["cellname"][x][0]
                              for x, c in enumerate(cliques)
                              if c.get("sigma") and j in c["mem"]
                              and x in st["cellname"]}
                    cand_names = [
                        nm for nm in NAMES
                        if nm not in st["used_names"]
                        and nm[0] not in st["pused"][j]
                        and nm[0] not in pinned
                        and canon_a(nm[1], st["aseen"])
                        and not (s == 10 and not st["used_names"]
                                 and nm[0] != 1)]
                for nm in cand_names:
                    y = NAMEV[nm]
                    for u in sorted(BFULL - st["ub"][j]):
                        if not canon_b(u, st["bseen"]):
                            continue
                        if not bump():
                            return
                        rows, nf, nk = [], [], []
                        if j in rex:
                            rex2 = rex
                            rows.append(vsub(vsub(y, b9(u)),
                                             rex[j]))
                        else:
                            rex2 = {**rex, j: vsub(y, b9(u))}
                            af, ak = anchor_forms(
                                j, rex2, st["cells"],
                                meta["non_edges"])
                            nf += af
                            nk += ak
                        cells2 = st["cells"]
                        if not named:
                            of, ok_ = cell_outside_forms(y, mem,
                                                         rex2)
                            nf += of
                            nk += ok_
                            cells2 = cells2 + ((y, mem, False),)

                        def rec(stt, nm=nm, u=u, rex2=rex2,
                                cells2=cells2, named=named,
                                mem=mem):
                            s2 = dict(stt)
                            s2["rex"] = rex2
                            s2["cells"] = cells2
                            s2["used_names"] = (
                                stt["used_names"]
                                if named else
                                stt["used_names"] | {nm})
                            s2["pused"] = {
                                **stt["pused"],
                                j: stt["pused"][j] | {nm[0]}}
                            s2["ub"] = {**stt["ub"],
                                        j: stt["ub"][j] | {u}}
                            s2["aseen"] = stt["aseen"] | {nm[1]}
                            s2["bseen"] = stt["bseen"] | {u}
                            if not named and kind == "k3":
                                s2["cellname"] = {
                                    **stt["cellname"], ci: nm}
                            s2["assign"] = stt["assign"] + (
                                (kind, j, ci if ci is not None
                                 else -1, nm[0], nm[1], u),)
                            return s2

                        emit(k, lat, S, inspan, st, rows, nf, nk,
                             rec)
                        if stats["capped"]:
                            return
                return
            if kind == "cell":
                assert j in rex, "fresh cell before anchor"
                mem = frozenset(cliques[ci]["mem"])
                for u in sorted(BFULL - st["ub"][j]):
                    if not canon_b(u, st["bseen"]):
                        continue
                    if not bump():
                        return
                    v = vadd(b9(u), rex[j])
                    nf, nk = [], []
                    if not relax:
                        for nm in NAMES:
                            nf.append(vsub(v, NAMEV[nm]))
                            nk.append("eq")
                    for (v2, _m2, isc) in st["cells"]:
                        if isc:
                            nf.append(vsub(v, v2))
                            nk.append("eq")
                    of, ok_ = cell_outside_forms(v, mem, rex)
                    nf += of
                    nk += ok_

                    def rec(stt, u=u, v=v, mem=mem):
                        s2 = dict(stt)
                        s2["cells"] = stt["cells"] + (
                            (v, mem, True),)
                        s2["cellv"] = {**stt["cellv"], ci: v}
                        s2["ub"] = {**stt["ub"],
                                    j: stt["ub"][j] | {u}}
                        s2["bseen"] = stt["bseen"] | {u}
                        s2["assign"] = stt["assign"] + (
                            ("cell", j, ci, -1, -1, u),)
                        return s2

                    emit(k, lat, S, inspan, st, [], nf, nk, rec)
                    if stats["capped"]:
                        return
                return
            # kind == "join"
            v = st["cellv"].get(ci)
            if v is None and cliques[ci]["sigma"]:
                v = NAMEV[st["cellname"][ci]]
            assert v is not None, "join before cell creation"
            for u in sorted(BFULL - st["ub"][j]):
                if not canon_b(u, st["bseen"]):
                    continue
                if not bump():
                    return
                rows, nf, nk = [], [], []
                if j in rex:
                    rex2 = rex
                    rows.append(vsub(vsub(v, b9(u)), rex[j]))
                else:
                    rex2 = {**rex, j: vsub(v, b9(u))}
                    af, ak = anchor_forms(j, rex2, st["cells"],
                                          meta["non_edges"])
                    nf += af
                    nk += ak

                def rec(stt, u=u, rex2=rex2):
                    s2 = dict(stt)
                    s2["rex"] = rex2
                    s2["ub"] = {**stt["ub"], j: stt["ub"][j] | {u}}
                    s2["bseen"] = stt["bseen"] | {u}
                    s2["assign"] = stt["assign"] + (
                        ("join", j, ci, -1, -1, u),)
                    return s2

                emit(k, lat, S, inspan, st, rows, nf, nk, rec)
                if stats["capped"]:
                    return
            return

        root = {"rex": ({meta["order"][0]: ZERO9} if relax
                        else {}),
                "cellv": {}, "cellname": {},
                "used_names": frozenset(),
                "pused": {j: frozenset() for j in range(1, 7)},
                "ub": {j: frozenset() for j in range(1, 7)},
                "aseen": (frozenset({3, 4}) if s == 8
                          else frozenset()),
                "bseen": frozenset(),
                "cells": (), "assign": ()}
        walk(0, ZLattice(), S0, {}, root)

    parts = enum_partitions(prof)
    stats["partitions"] = len(parts)
    for pidx, part in enumerate(parts):
        if stats["capped"]:
            break
        run_partition(pidx, part)
        if progress:
            print(f"  partition {pidx + 1}/{len(parts)} done "
                  f"nodes={stats['nodes']} dead={stats['dead']}"
                  f" term={stats['terminals']}", flush=True)
    stats["secs"] = round(time.time() - t0, 1)
    stats["terminal_patterns"] = terminals
    return stats


def p26_inject(pname: str, member_paths: list, cap: int) -> dict:
    """Map censused candidate 6-sets (collision half-systems) from
    live members into the abstract frame; assert every lattice row
    concretely in the group; the battery must NOT kill them."""
    prof = PROFILES[pname]
    a_t, s_t, n1, n2, n3, n4 = prof
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent
                           / "src"))
    from a17_e18_k4_census import diffs
    from bb_lab.group import AbelianGroup
    from bb_lab.poly import Poly

    def gmul(G, k, x):
        acc = tuple(0 for _ in range(G.rank))
        step = x if k > 0 else G.neg(x)
        for _ in range(abs(k)):
            acc = G.add(acc, step)
        return acc

    stat = static9()
    out = {"candidates": 0, "profile_hits": 0, "off_partition": 0,
           "mapped": 0, "validated": 0, "killed": []}
    rows_in = []
    for path in member_paths:
        with open(path) as f:
            rows_in += [json.loads(l) for l in f if '"A"' in l]
    for rmem in rows_in:
        if out["mapped"] >= cap:
            break
        G = AbelianGroup(tuple(rmem["frame"]))
        B = sorted(Poly.from_string(rmem["B"], G).support)
        dBs = diffs(frozenset(B), G)
        dB = sorted(dBs)
        zero = tuple(0 for _ in range(G.rank))
        for s4 in combinations(dB, 4):
            if out["mapped"] >= cap:
                break
            p4 = sum(1 for u, v in combinations(s4, 2)
                     if G.sub(v, u) in dBs)
            need = a_t - 4 - p4
            if not 0 <= need <= 5:
                continue
            base = (zero,) + s4
            bset = frozenset(base)
            for v6 in G:
                if v6 in bset:
                    continue
                deg = sum(1 for u in base if G.sub(v6, u) in dBs)
                if deg != need:
                    continue
                T = list(base) + [v6]
                out["candidates"] += 1
                cells = {}
                for ji, t in enumerate(T):
                    for b in B:
                        y = G.add(b, t)
                        cells.setdefault(y, []).append((ji, b))
                sig = sum(1 for v in cells.values() if len(v) & 1)
                if sig != s_t:
                    continue
                from collections import Counter
                hist = Counter(len(v) for v in cells.values())
                if (hist.get(2, 0), hist.get(3, 0),
                        hist.get(4, 0)) != (n2, n3, n4) or \
                        hist.get(5, 0) or hist.get(6, 0):
                    out["off_partition"] += 1
                    continue
                # per-translate C2 σ-cap (a MATCH consequence):
                # candidates violating it never enter the abstract
                # frame — the engine's partitions have d_j ≤ 2
                sd = [0] * 6
                for y, mm in cells.items():
                    if len(mm) & 1:
                        for (ji, _b) in mm:
                            sd[ji] += 1
                if max(sd) > 2:
                    out["c2_violating"] = out.get(
                        "c2_violating", 0) + 1
                    continue
                out["profile_hits"] += 1
                if out["mapped"] >= cap:
                    continue
                out["mapped"] += 1
                # BFS order over clique adjacency
                shared = [(y, mm) for y, mm in cells.items()
                          if len(mm) >= 2]
                adj = {j: set() for j in range(6)}
                for y, mm in shared:
                    for (j1, _), (j2, _) in combinations(mm, 2):
                        adj[j1].add(j2)
                        adj[j2].add(j1)
                degc = {j: sum(1 for _, mm in shared
                               if j in {x[0] for x in mm})
                        for j in range(6)}
                order = [max(degc, key=lambda j: (degc[j], -j))]
                while len(order) < 6:
                    nxt = sorted((j for j in range(6)
                                  if j not in order
                                  and adj[j] & set(order)),
                                 key=lambda j: (-degc[j], j))
                    assert nxt, "collision graph disconnected"
                    order.append(nxt[0])
                pos = {j: k for k, j in enumerate(order)}
                bL: dict = {}

                def blab(b):
                    if b not in bL:
                        bL[b] = len(bL)
                    return bL[b]

                rex = {order[0]: ZERO9}
                rows = []
                cellrec = []
                for y, mm in sorted(
                        shared,
                        key=lambda x: min(pos[t[0]]
                                          for t in x[1])):
                    mm = sorted(mm, key=lambda t: pos[t[0]])
                    j0, b0 = mm[0]
                    assert j0 in rex, "unanchored first member"
                    vab = vadd(b9(blab(b0)), rex[j0])
                    mem = frozenset(t[0] for t in mm)
                    cellrec.append((vab, mem))
                    for (j, b) in mm[1:]:
                        assert G.add(b, T[j]) == y
                        if j not in rex:
                            rex[j] = vsub(vab, b9(blab(b)))
                        else:
                            rows.append(vsub(
                                vadd(b9(blab(b)), rex[j]), vab))
                # concrete verification of every abstract row
                bInv = {v: k for k, v in bL.items()}
                for ell in range(5):
                    if ell not in bInv:
                        bInv[ell] = next(b for b in B
                                         if b not in bL)
                        bL[bInv[ell]] = ell
                xs = [zero] * 4 + [G.sub(bInv[ell], bInv[0])
                                   for ell in (1, 2, 3, 4)] + [zero]
                for row in rows:
                    assert all(c == 0 for c in row[:4]) \
                        and row[8] == 0, "non-b coefficient"
                    acc = zero
                    for cf, xv in zip(row, xs):
                        if cf:
                            acc = G.add(acc, gmul(G, cf, xv))
                    assert all(v == 0 for v in acc), \
                        f"row fails concretely: {row}"
                # battery: static + outside + non-edge + distinct
                nf, nk = [], []
                for (vab, mem) in cellrec:
                    for j in range(6):
                        if j in mem:
                            continue
                        for t in range(5):
                            nf.append(vsub(vsub(vab, rex[j]),
                                           b9(t)))
                            nk.append("eq")
                epairs = set()
                for _y, mm in shared:
                    ms = sorted({t[0] for t in mm})
                    for pr in combinations(ms, 2):
                        epairs.add(frozenset(pr))
                for j1, j2 in combinations(range(6), 2):
                    if frozenset((j1, j2)) in epairs:
                        continue
                    d0 = vsub(rex[j1], rex[j2])
                    nf.append(d0)
                    nk.append("eq")
                    for x in range(5):
                        for y_ in range(5):
                            if x != y_:
                                nf.append(vsub(
                                    d0, vsub(b9(x), b9(y_))))
                                nk.append("eq")
                for (v1, _m1), (v2, _m2) in combinations(
                        cellrec, 2):
                    nf.append(vsub(v1, v2))
                    nk.append("eq")
                lat = ZLattice()
                for row in rows:
                    lat.add(list(row))
                bat = Battery(stat.forms + nf, stat.kinds + nk)
                hit = bat.first_kill(lat)
                if hit is None:
                    out["validated"] += 1
                else:
                    out["killed"].append(
                        {"frame": rmem["frame"], "form": int(hit),
                         "T": [list(t) for t in T]})
    return out


def selftest() -> None:
    import random
    stat = static9()
    lat = ZLattice()
    lat.add(list(vsub(a9(1), b9(1))))
    assert stat.first_kill(lat) is not None, "a1=b1 D2 kill"
    lat2 = ZLattice()
    lat2.add(list(vsub(a9(1), G2)))
    assert stat.first_kill(lat2) is None, "g2 axis is free"
    sf, sk = stratum_forms(10)
    b = Battery(stat.forms + sf, stat.kinds + sk)
    lat3 = ZLattice()
    lat3.add(list(vsub(G2, vsub(a9(1), a9(2)))))
    assert b.first_kill(lat3) is not None, "delta in dA kills s=10"
    lat4 = ZLattice()
    lat4.add([2 * x for x in vsub(G2, vsub(a9(1), a9(2)))])
    assert b.first_kill(lat4) is None, "2-torsion delta allowed"
    rng = random.Random(11)
    for _ in range(50):
        x = [rng.randrange(-40, 40) for _ in range(9)]

        def val(vec):
            return sum(c * z for c, z in zip(vec, x))

        def av(i):
            return 0 if i == 0 else x[i - 1]

        def bv(i):
            return 0 if i == 0 else x[4 + i - 1]

        g2v = x[8]
        i = rng.choice([1, 2])
        p, u = rng.randrange(5), rng.randrange(5)
        y = name_vec(i, p)
        assert val(y) == av(p) + (0 if i == 1 else g2v)
        rex = vsub(y, b9(u))
        i2 = rng.choice([1, 2])
        p2, t = rng.randrange(5), rng.randrange(5)
        row = vsub(vsub(name_vec(i2, p2), b9(t)), rex)
        lhs = av(p2) + (0 if i2 == 1 else g2v) - bv(t)
        rhs = av(p) + (0 if i == 1 else g2v) - bv(u)
        assert val(row) == lhs - rhs, "cell row encoding"
        vcell = vadd(b9(t), rex)
        row2 = vsub(vadd(b9(u), vcell), vadd(b9(t), vcell))
        assert val(row2) == bv(u) - bv(t), "join row encoding"
    for name, prof in PROFILES.items():
        parts = enum_partitions(prof)
        assert parts, f"no partitions for {name}"
        meta = partition_meta(parts[0])
        per = {j: 0 for j in range(1, 7)}
        for _kind, j, _ci in meta["steps"]:
            per[j] += 1
        assert all(v == 5 for v in per.values()), (name, per)
        metar = partition_meta(parts[0], relax=True)
        perr = {j: 0 for j in range(1, 7)}
        for _kind, j, _ci in metar["steps"]:
            perr[j] += 1
        want = {j: sum(1 for c in metar["cliques"]
                       if j in c["mem"]) for j in perr}
        assert perr == want, (name, perr, want)
    print("selftest OK")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["selftest", "partitions", "p26",
                             "p26-relaxed", "p26-inject"])
    ap.add_argument("--profile", default=None)
    ap.add_argument("--members", default=None)
    ap.add_argument("--cap", type=int, default=300)
    ap.add_argument("--max-nodes", type=int, default=400_000_000)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.mode == "selftest":
        selftest()
        return
    if args.mode == "partitions":
        for name, prof in PROFILES.items():
            parts = enum_partitions(prof)
            print(f"{name}: {prof} -> {len(parts)} partitions")
        return
    prof = PROFILES[args.profile]
    if args.mode == "p26-inject":
        res = p26_inject(args.profile, args.members.split(","),
                         args.cap)
        print(json.dumps(res, indent=1))
        ok = (res["mapped"] > 0 and res["killed"] == []
              and res["validated"] == res["mapped"])
        print(f"\nINJECTION {'OK' if ok else 'FAILED'}: "
              f"{res['validated']}/{res['mapped']} mapped "
              f"candidates survive ({res['candidates']} censused)")
        sys.exit(0 if ok else 1)
    relax = args.mode == "p26-relaxed"
    res = p26_search(args.profile, args.max_nodes, relax=relax)
    summary = {k: v for k, v in res.items()
               if k != "terminal_patterns"}
    print(json.dumps(summary, indent=1))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, indent=1))
    if relax:
        ok = res["terminals"] > 0
        print(f"\nVALIDATION {'OK' if ok else 'FAILED'}: relaxed "
              f"run {'found' if ok else 'found NO'} terminals")
        sys.exit(0 if ok else 1)
    ok = not res["capped"]
    print(f"\nP-26[{args.profile}] table "
          f"{'COMPLETE' if ok else 'CAPPED'}: "
          f"terminals={res['terminals']}, nodes={res['nodes']}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
