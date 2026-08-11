#!/usr/bin/env python3
"""A32 M4 scoping — triple-equation architecture for the [[150,30,10]] floors.

The M4 floors ("0 < |v| <= 9 kernel vector => generator row", both sides)
reduce to classifying light solutions of SINGLE triple equations over
F2[G], G = C5xS3, |G| = 30:

  X side (ker H_X = dual cycles), instance beta in {0,1}:
      a0*u + a1*w = t*b_beta          (left-mults on the pair, right on t)
  Z side (ker H_Z = cycles), instance alpha in {0,1}:
      u*b0~ + w*b1~ = a_alpha~ * t    (~ = antipode; right-mults on pair)

because the two check blocks of each side couple ONLY through the shared
block t = v4, and a1 / b_beta (resp. b0~/b1~, a1~) are invertible, so any
two of (u, w, t) determine the third (up to a 4-element Ann(a0) coset in
one direction).  A full kernel vector v with |v| <= 9 restricts to a
<=9-light solution of EACH triple; classified triples then join along t.

This script is the falsify-first stage: it
  F1  rebuilds the block operator matrices and re-verifies every
      structural fact the Lean files will lean on (ranks, determinacy,
      Ann bases, parity, row-fragment families, join distinctness);
  F2  SAT-enumerates ALL <=9-light solutions of each of the 4 triples
      (exhaustive, blocking clauses) and checks they are EXACTLY the 60
      expected row fragments — the exact statement of the Lean
      classification lemmas;
  F3  emits the per-split sweep plan (which pair to sweep, class counts,
      native-op estimates) that the Lean sweep files will implement.

Outputs: data-a32/m4_scoping.json + stdout report.
Run:  cd experiments/bb_lab && uv run python scripts/a32_m4_scoping.py
"""

from __future__ import annotations

import importlib.util
import json
import time
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LAB_ROOT = HERE.parent
GROUPS = LAB_ROOT / "instances" / "mitten_groups"
OUT = LAB_ROOT / "data-a32"

SETS = dict(a0=(0, 14, 23), a1=(0, 2, 11), b0=(7, 20, 24), b1=(0, 2, 29))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- F1


def block_matrices(G):
    """All four 30x30 operator blocks per side, verified against the
    closed-form H entries used by the emitted Lean `Data.lean`."""
    n = 30
    mul, inv = G.mul, G.inv

    def Lm(a):  # (a*v)(h) = sum_x [h x^-1 in a] v(x)
        M = np.zeros((n, n), dtype=np.uint8)
        for h in range(n):
            for x in range(n):
                M[h, x] = int(mul(h, int(inv[x])) in a)
        return M

    def Rm(b):  # (v*b)(h) = sum_x [x^-1 h in b] v(x)
        M = np.zeros((n, n), dtype=np.uint8)
        for h in range(n):
            for x in range(n):
                M[h, x] = int(mul(int(inv[x]), h) in b)
        return M

    def Lms(a):  # (a~*v)(y) = sum_x [x y^-1 in a] v(x)
        M = np.zeros((n, n), dtype=np.uint8)
        for y in range(n):
            for x in range(n):
                M[y, x] = int(mul(x, int(inv[y])) in a)
        return M

    def Rms(b):  # (v*b~)(y) = sum_x [y^-1 x in b] v(x)
        M = np.zeros((n, n), dtype=np.uint8)
        for y in range(n):
            for x in range(n):
                M[y, x] = int(mul(int(inv[y]), x) in b)
        return M

    A0, A1 = Lm(set(SETS["a0"])), Lm(set(SETS["a1"]))
    B0, B1 = Rm(set(SETS["b0"])), Rm(set(SETS["b1"]))
    A0s, A1s = Lms(set(SETS["a0"])), Lms(set(SETS["a1"]))
    B0s, B1s = Rms(set(SETS["b0"])), Rms(set(SETS["b1"]))
    return dict(A0=A0, A1=A1, B0=B0, B1=B1, A0s=A0s, A1s=A1s,
                B0s=B0s, B1s=B1s)


def rank2(M) -> int:
    A = M.copy().astype(np.uint8)
    r = 0
    for c in range(A.shape[1]):
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        p = r + nz[0]
        A[[r, p]] = A[[p, r]]
        for q in np.flatnonzero(A[:, c]):
            if q != r:
                A[q] ^= A[r]
        r += 1
        if r == A.shape[0]:
            break
    return r


def nullspace2(M) -> np.ndarray:
    """Basis of {v : M v = 0} over GF(2), rows = basis vectors."""
    m, n = M.shape
    A = np.concatenate([M.copy().astype(np.uint8), np.eye(n, dtype=np.uint8)]).T
    # rows of A = [column of M | e_col]; eliminate on the first m coords
    r = 0
    for c in range(m):
        nz = [i for i in range(r, n) if A[i, c]]
        if not nz:
            continue
        A[[r, nz[0]]] = A[[nz[0], r]]
        for q in range(n):
            if q != r and A[q, c]:
                A[q] ^= A[r]
        r += 1
    ker = A[r:, m:] % 2
    assert not (M @ ker.T % 2).any()
    return ker


def inv2(M) -> np.ndarray:
    n = M.shape[0]
    A = np.concatenate([M.copy().astype(np.uint8), np.eye(n, dtype=np.uint8)], axis=1)
    for c in range(n):
        nz = [i for i in range(c, n) if A[i, c]]
        assert nz, "singular"
        A[[c, nz[0]]] = A[[nz[0], c]]
        for q in range(n):
            if q != c and A[q, c]:
                A[q] ^= A[c]
    W = A[:, n:] % 2
    assert ((M @ W) % 2 == np.eye(n, dtype=np.uint8)).all()
    return W


def f1_structure(G, M, HX, HZ) -> dict:
    n = 30
    mul, inv = G.mul, G.inv
    # closed-form H agreement: grid entry (check (beta,h), qubit (m,x))
    for beta in range(2):
        for m, blk in ((beta, "A0"), (2 + beta, "A1")):
            assert (HX[30 * beta:30 * beta + 30, 30 * m:30 * m + 30]
                    == M[blk]).all()
        assert (HX[30 * beta:30 * beta + 30, 120:150]
                == M[f"B{beta}"]).all()
    for alpha in range(2):
        for m, blk in ((2 * alpha, "B0s"), (2 * alpha + 1, "B1s")):
            assert (HZ[30 * alpha:30 * alpha + 30, 30 * m:30 * m + 30]
                    == M[blk]).all()
        assert (HZ[30 * alpha:30 * alpha + 30, 120:150]
                == M[f"A{alpha}s"]).all()

    ranks = {k: rank2(v) for k, v in M.items()}
    assert ranks["A1"] == ranks["B0"] == ranks["B1"] == 30
    assert ranks["A1s"] == ranks["B0s"] == ranks["B1s"] == 30
    assert ranks["A0"] == ranks["A0s"] == 28

    annA0 = nullspace2(M["A0"])          # dim 2
    annA0s = nullspace2(M["A0s"])
    ann_w = sorted(int(v.sum()) for v in _span(annA0) if v.any())
    ann_ws = sorted(int(v.sum()) for v in _span(annA0s) if v.any())
    # left-null certificates (solvability filters for the coset sweeps)
    lnA0 = nullspace2(M["A0"].T)
    lnA0s = nullspace2(M["A0s"].T)

    # row-fragment families + join distinctness
    def lset(s, y):  # {s*y : s in set}
        return frozenset(mul(s, y) for s in s)

    def rset(y, s):  # {y*s}
        return frozenset(mul(y, s) for s in s)

    a0, a1 = set(SETS["a0"]), set(SETS["a1"])
    b0, b1 = set(SETS["b0"]), set(SETS["b1"])
    b0i = {int(inv[s]) for s in b0}
    b1i = {int(inv[s]) for s in b1}
    a0i = {int(inv[s]) for s in a0}
    a1i = {int(inv[s]) for s in a1}
    tX = [lset(a0, y) for y in range(n)] + [lset(a1, y) for y in range(n)]
    tZ = [rset(h, b0i) for h in range(n)] + [rset(h, b1i) for h in range(n)]
    assert len(set(tX)) == 60, "X-side t-values not distinct"
    assert len(set(tZ)) == 60, "Z-side t-values not distinct"

    return dict(ranks=ranks, ann_weights=ann_w, ann_weights_star=ann_ws,
                annA0=annA0, annA0s=annA0s, lnA0=lnA0, lnA0s=lnA0s,
                sets=dict(a0=a0, a1=a1, b0=b0, b1=b1,
                          a0i=a0i, a1i=a1i, b0i=b0i, b1i=b1i))


def _span(basis: np.ndarray):
    k = basis.shape[0]
    for mask in range(1 << k):
        v = np.zeros(basis.shape[1], dtype=np.uint8)
        for i in range(k):
            if mask >> i & 1:
                v ^= basis[i]
        yield v


# ---------------------------------------------------------------- F2


def triple_matrix(P, Q, R) -> np.ndarray:
    """30x90 system  P u + Q w + R t = 0  (columns ordered u | w | t)."""
    return np.concatenate([P, Q, R], axis=1).astype(np.uint8)


def families(G, F1, side: str, idx: int) -> list[tuple]:
    """Expected 60 light solutions (u, w, t) as support-frozensets."""
    mul, inv = G.mul, G.inv
    s = F1["sets"]
    out = []
    if side == "X":
        b = s[f"b{idx}"]
        for y in range(30):
            yb = frozenset(mul(y, t) for t in b)
            out.append((yb, frozenset(), frozenset(mul(t, y) for t in s["a0"])))
            out.append((frozenset(), yb, frozenset(mul(t, y) for t in s["a1"])))
    else:
        a_i = s[f"a{idx}i"]
        for h in range(30):
            ah = frozenset(mul(t, h) for t in a_i)
            out.append((ah, frozenset(), frozenset(mul(h, t) for t in s["b0i"])))
            out.append((frozenset(), ah, frozenset(mul(h, t) for t in s["b1i"])))
    return out


def f2_censuses(G, M, F1, enum) -> dict:
    """SAT-exhaust each triple at weight <= 9; compare to families."""
    systems = dict(
        X0=triple_matrix(M["A0"], M["A1"], M["B0"]),
        X1=triple_matrix(M["A0"], M["A1"], M["B1"]),
        Z0=triple_matrix(M["B0s"], M["B1s"], M["A0s"]),
        Z1=triple_matrix(M["B0s"], M["B1s"], M["A1s"]),
    )
    report = {}
    for name, H in systems.items():
        side, idx = name[0], int(name[1])
        t0 = time.perf_counter()
        sols, secs, exhausted = enum(H, 9, log=lambda s: None)
        assert exhausted
        got = set()
        splits = {}
        for v in sols:
            u, w, t = v[:30], v[30:60], v[60:]
            got.add((frozenset(np.flatnonzero(u).tolist()),
                     frozenset(np.flatnonzero(w).tolist()),
                     frozenset(np.flatnonzero(t).tolist())))
            key = (int(u.sum()), int(w.sum()), int(t.sum()))
            splits[str(key)] = splits.get(str(key), 0) + 1
        fam = set(families(G, F1, side, idx))
        extras = got - fam
        missing = fam - got
        report[name] = dict(
            n_solutions=len(sols), seconds=round(secs, 1),
            splits=splits, n_extras=len(extras), n_missing=len(missing),
        )
        print(f"  [{name}] {len(sols)} light solutions in {secs:.1f}s; "
              f"extras={len(extras)} missing={len(missing)} splits={splits}")
        assert not missing, f"{name}: family member not found by SAT?!"
        if extras:
            print(f"    !! STRAYS: {sorted(map(sorted, list(extras)[:5]))}")
    return report


# ---------------------------------------------------------------- F3


def f3_sweep_plan(F1) -> dict:
    """Per-split sweep-pair choice + cost model for the Lean leaves.

    Split (p, q, r) = (|u|, |w|, |t|), p+q+r even (in-triple parity),
    p+q+r <= 8 (odd 9 impossible), plus the two weight-6 classify splits
    (3,0,3)/(0,3,3).  Sweep-pair preference: unique-derivation pairs
    [(u,t)->w, (u,w)->t on the X side; (u,t)->w, (w,t)->u on Z] before
    the Ann(a0)-coset pair; among available, smallest product."""
    plans = {}
    for side in ("X", "Z"):
        rows = []
        total_classes = 0
        for s in range(0, 9, 2):
            for p in range(s + 1):
                for q in range(s + 1 - p):
                    r = s - p - q
                    opts = []
                    # (u,t) -> w unique (X: w = A1inv(A0 u + t B); Z: mirror)
                    opts.append((comb(30, p) * comb(30, r), "ut", 1))
                    if side == "X":
                        # (u,w) -> t unique (right-mult by b inverse)
                        opts.append((comb(30, p) * comb(30, q), "uw", 1))
                        # (w,t) -> u in 4-coset of Ann(a0), if solvable
                        opts.append((comb(30, q) * comb(30, r) * 4, "wt", 4))
                    else:
                        # (w,t) -> u unique (B0s invertible)
                        opts.append((comb(30, q) * comb(30, r), "wt", 1))
                        # (u,w) -> t in 4-coset (alpha=0) / unique (alpha=1)
                        opts.append((comb(30, p) * comb(30, q) * 4, "uw", 4))
                    cost, pair, mult = min(opts)
                    rows.append(dict(split=(p, q, r), pair=pair,
                                     classes=cost, coset=mult))
                    total_classes += cost
        plans[side] = dict(splits=rows, total_classes=total_classes)
        worst = max(rows, key=lambda d: d["classes"])
        print(f"  [{side}] {len(rows)} splits, total classes "
              f"{total_classes:,} (x2 instances); worst {worst['split']} "
              f"-> {worst['classes']:,} via {worst['pair']}")
    return plans


def main() -> None:
    a26 = _load("a26_mitten_descent")
    assert a26.MITTEN["150,30,10"] == dict(gid=(30, 1), **SETS)
    scoping = _load("a32_m150_lean_scoping")
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    HX, HZ = a26.mitten_code(G, **SETS)

    print("[F1] block matrices + structure")
    M = block_matrices(G)
    F1 = f1_structure(G, M, HX, HZ)
    print(f"  ranks {F1['ranks']}")
    print(f"  Ann(a0) nonzero weights {F1['ann_weights']}, "
          f"Ann(a0~) {F1['ann_weights_star']}")

    print("[F2] exhaustive light-triple censuses (SAT, <=9)")
    F2 = f2_censuses(G, M, F1, scoping.enumerate_light_kernel)

    print("[F3] sweep plan")
    F3 = f3_sweep_plan(F1)

    OUT.mkdir(exist_ok=True)
    out = OUT / "m4_scoping.json"
    out.write_text(json.dumps(dict(
        ranks=F1["ranks"], ann_weights=F1["ann_weights"],
        ann_weights_star=F1["ann_weights_star"],
        censuses=F2, sweep_plan=F3), indent=1, default=str))
    print(f"[done] wrote {out}")


if __name__ == "__main__":
    main()
