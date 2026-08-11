"""A32/M0: Lean-certification scoping for the [[150,30,10]] mitten code.

Ground-truth measurements that pin the statement shape and sweep budget
of the planned QECLean formalization (pipeline/attempts/mitten_150_30_10/,
plan.md §M0). Two independent parts:

**census** — exact enumeration of ALL weight-≤W kernel vectors of both
CSS directions (default W = 9 = d−1) on the sha-pinned shipped matrices,
via CryptoMiniSat native-XOR + seqcounter cardinality + blocking clauses
(same encoding family as `bb_lab.sat_distance`; no logical-indicator
constraint — boundaries are wanted here, they ARE the expected answer).
Statement under test (M4 floor shape, ISD-supported at 6k information
sets/side):

    the only nonzero weight-≤9 vectors of ker H_X are the 60 rows of
    H_Z (all weight exactly 9), and mirror.

Census cardinality is invariant under qubit relabeling, so the shipped
byte order (certificate-pinned) is authoritative for the rebuilt
Table XIII labeling the Lean instance will use (A26 §5 established the
Tanner isomorphism).

**symmetry** — deterministic re-derivation, on the canonical Table XIII
rebuild (group_30_1.txt, Eq. (J1) via a26_mitten_descent helpers), of
the code-automorphism facts the M4 split map normalizes against:

  1. blockwise L/R translation symmetries: exhaustive over all 30
     elements × 2^9 assignments (finding: exactly the 5 central ones);
  2. Aut(G)-induced symmetries: all 24 automorphisms of C₅×S₃ tested
     for setwise preservation of a0,a1,b0,b1 (finding: identity only);
  3. natural-family X↔Z dualities: grid permutation × antipode on/off
     per sector (finding: none);
  4. optional (pynauty importable, e.g. `uv run --with pynauty`): TRUE
     Tanner Aut group order + a definitive (Hx,Hz)≅(Hz,Hx) iso verdict
     via the colored-graph certificates of a26_mitten_descent.

Findings of record land in pipeline/attempts/mitten_150_30_10/ and
data-a32/m150_scoping.json.

Run (from experiments/bb_lab):
  uv run python scripts/a32_m150_lean_scoping.py census   [--max-weight 9]
  uv run python scripts/a32_m150_lean_scoping.py symmetry
  uv run python scripts/a32_m150_lean_scoping.py all
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pycryptosat
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

LAB_ROOT = Path(__file__).resolve().parent.parent
INSTANCE = LAB_ROOT / "instances" / "mitten_150_30_10"
GROUPS = LAB_ROOT / "instances" / "mitten_groups"
OUT = LAB_ROOT / "data-a32"  # tracked (plain data/ is gitignored)

# Table XIII sets for [[150,30,10]] (GAP Elements(G) 0-based indices) —
# single source of truth is a26_mitten_descent.TABLE; re-declared here so
# the census half runs without importing the a26 module.
SETS = dict(a0=(0, 14, 23), a1=(0, 2, 11), b0=(7, 20, 24), b1=(0, 2, 29))


def _load_a26():
    spec = importlib.util.spec_from_file_location(
        "a26", Path(__file__).parent / "a26_mitten_descent.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.MITTEN["150,30,10"] == dict(gid=(30, 1), **SETS), (
        "SETS drifted from a26_mitten_descent.MITTEN"
    )
    return mod


# ------------------------------------------------------------- census


def enumerate_light_kernel(
    H: np.ndarray, max_weight: int, log=print
) -> tuple[list[np.ndarray], float, bool]:
    """All v ≠ 0 with H v = 0 and |v| ≤ max_weight, by blocking clauses.

    Returns (solutions, seconds, exhausted). `exhausted` is True iff the
    final call returned UNSAT (the enumeration is complete)."""
    n = H.shape[1]
    pool = IDPool()
    qubit_vars = [pool.id() for _ in range(n)]
    solver = pycryptosat.Solver()
    for row in H:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([qubit_vars[i] for i in idx], False)
    solver.add_clause(qubit_vars)  # v ≠ 0
    card = CardEnc.atmost(
        lits=qubit_vars, bound=max_weight, vpool=pool,
        encoding=EncType.seqcounter,
    )
    for cl in card.clauses:
        solver.add_clause(cl)

    sols: list[np.ndarray] = []
    t0 = time.perf_counter()
    while True:
        sat, model = solver.solve()
        if not sat:
            return sols, time.perf_counter() - t0, True
        v = np.array([1 if model[qv] else 0 for qv in qubit_vars],
                     dtype=np.uint8)
        assert not (H @ v % 2).any() and 1 <= int(v.sum()) <= max_weight
        sols.append(v)
        if len(sols) % 20 == 0:
            log(f"    ... {len(sols)} solutions, "
                f"{time.perf_counter() - t0:.1f}s")
        # block exactly this assignment on the qubit vars
        solver.add_clause(
            [-qubit_vars[i] if v[i] else qubit_vars[i] for i in range(n)]
        )


def census(max_weight: int) -> dict:
    Hx = np.load(INSTANCE / "Hx.npy").astype(np.uint8) % 2
    Hz = np.load(INSTANCE / "Hz.npy").astype(np.uint8) % 2
    report: dict = {"max_weight": max_weight, "sides": {}}
    for side, (Hker, Hrows) in {
        "kerHx_vs_rowsHz": (Hx, Hz),
        "kerHz_vs_rowsHx": (Hz, Hx),
    }.items():
        print(f"[census] {side}: enumerating wt ≤ {max_weight} kernel "
              f"vectors ...", flush=True)
        sols, secs, exhausted = enumerate_light_kernel(Hker, max_weight)
        rows = {tuple(r) for r in Hrows}
        sol_set = {tuple(s) for s in sols}
        extra = sorted(sol_set - rows)
        missing_rows = sorted(rows - sol_set)
        wts = sorted({int(s.sum()) for s in sols})
        print(f"[census] {side}: {len(sols)} solutions in {secs:.1f}s "
              f"(exhausted={exhausted}); weights {wts}; "
              f"non-row extras: {len(extra)}; "
              f"generator rows missed: {len(missing_rows)}", flush=True)
        report["sides"][side] = {
            "solutions": len(sols),
            "seconds": round(secs, 1),
            "exhausted": exhausted,
            "weights": wts,
            "non_row_extras": [
                [int(i) for i in np.flatnonzero(np.array(e))] for e in extra
            ],
            "generator_rows_missed": len(missing_rows),
        }
        assert exhausted, f"{side}: enumeration did not terminate UNSAT"
    verdict = all(
        s["non_row_extras"] == [] and s["generator_rows_missed"] == 0
        and s["solutions"] == 60
        for s in report["sides"].values()
    )
    report["rows_only"] = verdict
    print(f"[census] VERDICT: light kernel = exactly the 60 opposite rows "
          f"on both sides: {verdict}", flush=True)
    return report


# ----------------------------------------------------------- symmetry


def _perm_matrix_family(G):
    n = G.n
    ar = np.arange(n)
    Lp = lambda g: G.mt[g, ar]        # x ↦ g·x  (index map)
    Rp = lambda g: G.mt[ar, g]        # x ↦ x·g
    return Lp, Rp


def translation_symmetries(G, HX, HZ) -> dict:
    """Exhaustive blockwise-L/R translation automorphism search."""
    import itertools

    n = G.n
    Lp, Rp = _perm_matrix_family(G)

    def qcperm(maps, H, rowmaps):
        # build permuted H: rows by rowmaps (2 blocks), cols by maps (5)
        Pc = np.empty(H.shape[1], dtype=int)
        for c, p in enumerate(maps):
            for x in range(n):
                Pc[c * n + x] = c * n + int(p[x])
        Pr = np.empty(H.shape[0], dtype=int)
        for i, p in enumerate(rowmaps):
            for x in range(n):
                Pr[i * n + x] = i * n + int(p[x])
        out = np.zeros_like(H)
        out[Pr[:, None], Pc[None, :]] = H
        return out

    found = []
    for g in range(1, n):
        hit = None
        for qt in itertools.product("LR", repeat=5):
            qm = [Lp(g) if t == "L" else Rp(g) for t in qt]
            for xt in itertools.product("LR", repeat=2):
                xm = [Lp(g) if t == "L" else Rp(g) for t in xt]
                if not (qcperm(qm, HX, xm) == HX).all():
                    continue
                for zt in itertools.product("LR", repeat=2):
                    zm = [Lp(g) if t == "L" else Rp(g) for t in zt]
                    if (qcperm(qm, HZ, zm) == HZ).all():
                        hit = ("".join(qt), "".join(xt), "".join(zt))
                        break
                if hit:
                    break
            if hit:
                break
        if hit:
            found.append({"g": g, "order": int(G.orders[g]),
                          "assignment": hit})
    return {
        "count": len(found),
        "elements": found,
        "center": [int(c) for c in G.center()],
    }


def aut_g_symmetries(G) -> dict:
    """All automorphisms of G; which preserve the four sets setwise."""
    n = G.n
    ords = G.orders
    g5 = [g for g in range(n) if ords[g] == 5]
    g3 = [g for g in range(n) if ords[g] == 3]
    g2 = [g for g in range(n) if ords[g] == 2]
    gens = [g5[0], g3[0], g2[0]]
    auts = []
    for i5 in g5:
        for i3 in g3:
            for i2 in g2:
                img = {0: 0, g5[0]: i5, g3[0]: i3, g2[0]: i2}
                ok = True
                frontier = list(img.keys())
                while frontier and ok:
                    new = []
                    for a in list(img.keys()):
                        for b in gens:
                            c = G.mul(a, b)
                            ci = G.mul(img[a], img[b])
                            if c in img:
                                if img[c] != ci:
                                    ok = False
                                    break
                            else:
                                img[c] = ci
                                new.append(c)
                        if not ok:
                            break
                    frontier = new
                if ok and len(img) == n and len(set(img.values())) == n:
                    auts.append(dict(img))
    preserving = [
        s for s in auts
        if all(sorted(s[x] for x in SETS[k]) == sorted(SETS[k])
               for k in SETS)
    ]
    return {"aut_order": len(auts), "set_preserving": len(preserving)}


def duality_search(G, HX, HZ) -> dict:
    """Natural-family X↔Z candidates: grid slot permutation × antipode
    on/off per sector; verdict by row-multiset equality."""
    import itertools

    n = G.n
    inv = G.inv
    rowsZ = set(map(tuple, HZ))
    tried = 0
    for gridmap in itertools.permutations(range(4)):
        for ag in (0, 1):
            for a5 in (0, 1):
                P = np.zeros(5 * n, dtype=int)
                for old in range(4):
                    for x in range(n):
                        xx = int(inv[x]) if ag else x
                        P[gridmap[old] * n + xx] = old * n + x
                for x in range(n):
                    xx = int(inv[x]) if a5 else x
                    P[4 * n + xx] = 4 * n + x
                tried += 1
                if set(map(tuple, HX[:, P])) == rowsZ:
                    return {"found": True,
                            "map": {"grid": gridmap, "antipode_grid": ag,
                                    "antipode_d5": a5}}
    return {"found": False, "candidates_tried": tried}


def tanner_aut(HX, HZ) -> dict:
    """TRUE Tanner Aut group + duality, if pynauty is importable."""
    try:
        import pynauty  # noqa: F401
    except ImportError:
        return {"available": False,
                "note": "pynauty not importable — run via "
                        "`uv run --with pynauty` for the definitive census"}
    a26 = _load_a26()
    nq = HX.shape[1]
    nx, nz = HX.shape[0], HZ.shape[0]
    adj = {}
    for r in range(nx):
        adj[nq + r] = [int(c) for c in np.flatnonzero(HX[r])]
    for r in range(nz):
        adj[nq + nx + r] = [int(c) for c in np.flatnonzero(HZ[r])]
    g = pynauty.Graph(
        nq + nx + nz, directed=False, adjacency_dict=adj,
        vertex_coloring=[set(range(nq)),
                         set(range(nq, nq + nx)),
                         set(range(nq + nx, nq + nx + nz))],
    )
    gens, grpsize1, grpsize2, _, _ = pynauty.autgrp(g)
    same = a26.tanner_certificate(HX, HZ) == a26.tanner_certificate(HZ, HX)
    return {
        "available": True,
        "aut_group_size": f"{grpsize1}e{int(grpsize2)}",
        "n_generators": len(gens),
        "xz_duality": bool(same),
    }


def symmetry() -> dict:
    a26 = _load_a26()
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    HX, HZ = a26.mitten_code(G, **SETS)
    print("[symmetry] translation search (30 × 2^9 exhaustive) ...",
          flush=True)
    trans = translation_symmetries(G, HX, HZ)
    print(f"[symmetry] blockwise L/R translation automorphisms: "
          f"{trans['count']} (center C5 expected: 4 non-identity)",
          flush=True)
    autg = aut_g_symmetries(G)
    print(f"[symmetry] |Aut(G)| = {autg['aut_order']}, set-preserving: "
          f"{autg['set_preserving']} (identity only expected)", flush=True)
    dual = duality_search(G, HX, HZ)
    print(f"[symmetry] natural-family X↔Z duality: {dual}", flush=True)
    nauty = tanner_aut(HX, HZ)
    print(f"[symmetry] tanner (nauty): {nauty}", flush=True)
    return {"translations": trans, "aut_g": autg, "duality_natural": dual,
            "tanner": nauty}


# ----------------------------------------------------------- leandict


def leandict() -> dict:
    """Canonical parameterization GAP index → z^i·r^j·s^k (i<5, j<3, k<2)
    for the planned Lean carrier `ZMod 5 × DihedralGroup 3`.

    Picks z = the first central order-5 element, and (r, s) an S₃
    complement: r of order 3, s of order 2, with s·r·s = r⁻¹ and both
    commuting with z (automatic: z central). Validates that
    g ↦ (i, j, k) is a bijection and a homomorphism against the full GAP
    table under the abstract relations — the emitter (M2) re-validates
    the *mathlib orientation* (r j vs sr j and the DihedralGroup mul
    convention) by an in-Lean `decide` at generation time, so no mathlib
    convention is baked in here."""
    a26 = _load_a26()
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    n = G.n
    z = next(g for g in G.center() if G.orders[g] == 5)
    r = next(g for g in range(n) if int(G.orders[g]) == 3)
    ri = int(G.inv[r])
    s = next(
        g for g in range(n)
        if int(G.orders[g]) == 2 and G.mul(G.mul(g, r), g) == ri
    )
    # enumerate z^i r^j s^k
    table = {}
    zi = 0
    for i in range(5):
        w = zi
        for j in range(3):
            for k in range(2):
                el = w
                for _ in range(j):
                    el = G.mul(el, r)
                if k:
                    el = G.mul(el, s)
                assert el not in table, "parameterization not injective"
                table[el] = (i, j, k)
        zi = G.mul(zi, z)
    assert len(table) == n
    # homomorphism check: product of parameterized elements agrees with
    # the parameterization of the GAP product, under the abstract
    # relations of C5 × S3 (s r = r^{-1} s).
    def mul_param(p, q):
        (i1, j1, k1), (i2, j2, k2) = p, q
        if k1 == 0:
            return ((i1 + i2) % 5, (j1 + j2) % 3, k2)
        return ((i1 + i2) % 5, (j1 - j2) % 3, (k1 + k2) % 2)
    for a in range(n):
        for b in range(n):
            assert mul_param(table[a], table[b]) == table[G.mul(a, b)], (
                f"hom failure at {a},{b}")
    print(f"[leandict] z = {z}, r = {r}, s = {s}; bijective "
          f"parameterization g ↦ z^i·r^j·s^k validated as a hom on all "
          f"{n}×{n} products", flush=True)
    sets_param = {k: [table[g] for g in v] for k, v in SETS.items()}
    print(f"[leandict] sets in (i,j,k) coords: {sets_param}", flush=True)
    return {
        "z": int(z), "r": int(r), "s": int(s),
        "index_to_ijk": {str(g): list(table[g]) for g in range(n)},
        "sets_ijk": {k: [list(t) for t in v]
                     for k, v in sets_param.items()},
    }


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["census", "symmetry", "leandict",
                                     "all"])
    ap.add_argument("--max-weight", type=int, default=9)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "m150_scoping.json"
    report = json.loads(out_path.read_text()) if out_path.exists() else {}

    if args.mode in ("census", "all"):
        report["census"] = census(args.max_weight)
    if args.mode in ("symmetry", "all"):
        report["symmetry"] = symmetry()
    if args.mode in ("leandict", "all"):
        report["leandict"] = leandict()
    report["generated_by"] = "scripts/a32_m150_lean_scoping.py"
    out_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"[a32] report → {out_path}", flush=True)


if __name__ == "__main__":
    main()
