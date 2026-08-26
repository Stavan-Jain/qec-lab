#!/usr/bin/env python3
"""A40 S4 gate 2 — the phase frame triage: k of every sheared
y-period-p quotient lattice, small p.

A y-periodic "phase" of the transfer system with period p and x-shift
d per period is exactly an X-cycle of the BB code on the quotient
lattice Z^2 / <(l,0), (d,p)>.  Cheap phases (weight < 2p per period)
at frames with k = 0 do not exist (no cycles outside stabilizers is
not enough — k = 0 still allows trivial phases, but NONTRIVIAL cheap
phases need k > 0).  This triage computes, for l in {12, 18} and
p = 1..8, all shears d in Z_l:
  - the rectangular normal form of the quotient group (Smith),
  - the transported polynomial supports,
  - k of the resulting code,
and reports the k > 0 survivors that a census must then examine.
Sound reduction: the quotient group is abelian of order l*p; we find
an explicit isomorphism to Z_o1 x Z_o2 and transport the fixed Laurent
supports through it (verified: CSS condition + k invariance under the
choice by cross-checking a known frame).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def snf2(M):
    """Smith normal form of an integer 2x2 matrix: U M V = D with
    U, V unimodular; returns (D, U, V)."""
    M = [list(map(int, r)) for r in M]
    U = [[1, 0], [0, 1]]
    V = [[1, 0], [0, 1]]

    def rowop(i, j, q):  # row_i -= q*row_j
        for c in range(2):
            M[i][c] -= q * M[j][c]
            U[i][c] -= q * U[j][c]

    def colop(i, j, q):  # col_i -= q*col_j
        for r in range(2):
            M[r][i] -= q * M[r][j]
            V[r][i] -= q * V[r][j]

    def swap_rows():
        M[0], M[1] = M[1], M[0]
        U[0], U[1] = U[1], U[0]

    def swap_cols():
        for r in range(2):
            M[r][0], M[r][1] = M[r][1], M[r][0]
            V[r][0], V[r][1] = V[r][1], V[r][0]

    for _ in range(200):
        if M[1][0] == 0 and M[0][1] == 0:
            break
        if M[0][0] == 0:
            if M[1][0]:
                swap_rows()
            else:
                swap_cols()
            continue
        if M[1][0] % M[0][0] == 0 and M[0][1] % M[0][0] == 0:
            rowop(1, 0, M[1][0] // M[0][0])
            colop(1, 0, M[0][1] // M[0][0])
            if M[1][0] == 0 and M[0][1] == 0:
                break
            continue
        if abs(M[1][0]) and abs(M[1][0]) < abs(M[0][0]):
            swap_rows()
            continue
        if abs(M[0][1]) and abs(M[0][1]) < abs(M[0][0]):
            swap_cols()
            continue
        if M[1][0]:
            rowop(1, 0, M[1][0] // M[0][0])
        if M[0][1]:
            colop(1, 0, M[0][1] // M[0][0])
    if M[0][0] and M[1][1] % M[0][0]:
        # enforce divisibility d1 | d2
        colop(0, 1, -1)
        return snf2_from(M, U, V)
    return M, U, V


def snf2_from(M, U, V):
    D2, U2, V2 = snf2(M)
    U3 = [[sum(U2[i][t] * U[t][j] for t in range(2)) for j in range(2)]
          for i in range(2)]
    V3 = [[sum(V[i][t] * V2[t][j] for t in range(2)) for j in range(2)]
          for i in range(2)]
    return D2, U3, V3


def quotient_code(l, p, d, name=None):
    """BB code with the fixed Laurent pair on Z^2 / <(l,0),(d,p)>.

    Coordinates: exponent vector e in Z^2 maps to the class of e.  We
    find the SNF <U M V = D> of the relation matrix M = [[l,0],[d,p]]
    (rows = relations).  Then Z^2/L = Z^2 / rowspace(M); substituting
    f = e V (new basis), relations become rows of U M V = D: the group
    is Z_{D00} x Z_{D11} in the f-coordinates: f = e V.
    """
    M = [[l, 0], [d, p]]
    D, U, V = snf2(M)
    o1, o2 = abs(D[0][0]), abs(D[1][1])
    assert o1 * o2 == l * p, (l, p, d, D)

    def tr(supp):
        return frozenset(((e[0] * V[0][0] + e[1] * V[1][0]) % o1,
                          (e[0] * V[0][1] + e[1] * V[1][1]) % o2)
                         for e in supp)

    return TowerCode(name or f"shear(l={l},p={p},d={d})", (o1, o2),
                     tr(A_L), tr(B_L)), (o1, o2)


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")

    # cross-check the reduction on known frames: d = 0 must reproduce
    # the rectangular codes (k = 12 at (12,6),(18,12); k known at others)
    for (l, p), k_expect in [((12, 6), 12), ((18, 12), 12),
                             ((18, 6), 12), ((12, 12), 12)]:
        c, o = quotient_code(l, p, 0)
        assert c.k == k_expect, (l, p, c.k, k_expect)
    print("cross-check: d = 0 shears reproduce rectangular k "
          "(12,6)/(18,6)/(18,12)/(12,12) = 12")

    rows = []
    survivors = []
    for l in (12, 18):
        for p in range(1, 9):
            ks = {}
            for d in range(l):
                c, o = quotient_code(l, p, d)
                ks.setdefault(c.k, []).append(d)
            for k, ds in sorted(ks.items()):
                rows.append(dict(l=l, p=p, k=k, shears=ds))
                if k > 0:
                    survivors.append(dict(l=l, p=p, k=k, shears=ds))
            kmax = max(ks)
            print(f"l={l} p={p}: k by shear: "
                  + ", ".join(f"k={k} at d={ds}" for k, ds in
                              sorted(ks.items())))
    print("\nSURVIVORS (k > 0): frames a phase census must examine:")
    for s in survivors:
        print(f"  l={s['l']} p={s['p']} k={s['k']} shears {s['shears']}")

    out = dict(rows=rows, survivors=survivors,
               wall_s=round(time.time() - t0, 1))
    (DATA / "s4_phase_triage.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s4_phase_triage.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
