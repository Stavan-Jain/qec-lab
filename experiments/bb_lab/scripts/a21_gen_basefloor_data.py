#!/usr/bin/env python
"""A21 session 2: emit Z5Z15F2A6/BaseFloorData.lean (coset-sweep data).

Discharges the engineering leaf of `weight6_cycle_is_boundary`
(notes/A21_analytic_base_floor.md §4 stage A): the translation-reduced
coset sweep for the [[150,8,8]] base `bb_neigh_z5z15_f2a6f17e`.

Everything is computed in the REPO convention (QECLean BBChainComplex):
  conv P v (r) = sum_h P(h) v(r-h);  cycle: B*u_L + A*u_R = 0;
  boundary: (A*f, B*f);  cell index g=(x,y) -> x*15 + y  (0..74).

Emitted tables (all hard-asserted here before writing):
  * wABmask  — particular solution  A ⋆ wAB = B   (packed 75-bit mask)
  * wBAmask  — particular solution  B ⋆ wBA = A
  * zmask1-4 — echelon basis of the shared kernel Ann(A) = Ann(B),
               δ-normalized on the free cells fc1-4
  * fc1-4    — the free cells (RREF free columns of conv_A)
  * pivListA/pivListB — row-combination pivot certificates: 71 entries
               (pivot cell j, row-combination mask w) with
               y := conv P̃ w  (P̃ = reflected stencil; the transpose
               pairing) satisfying y(j) = 1 and y = 0 on every LATER
               pivot cell.  These are the RREF rows (y) with their
               row-operation combinations (w), pivots ordered by cell
               index; RREF zeroes y at *all* other pivots, so any order
               is valid.

NOTE (negative result, kept for the record): the plain no-row-op peel
certificate (KernelCert `pivB` shape — each check row hitting its pivot
and no later pivot) does NOT exist for the full-torus systems: the peel
closure of every 4-cell free set stalls at ≤ 10 of 75 cells (checked
20k random + structured seeds).  The A17 experience that such orders
"always exist in practice" is a *window-with-boundary* phenomenon; on
the closed torus there are no corner seeds.  Hence the row-combination
form.

Verified sweep facts (numpy, exact Lean quantifier semantics — the five
`native_decide` sweeps of BaseFloorSweep.lean must match these):
  S0 : all 15 nonzero kernel combos have weight 40 (never 6)
  S1 : ∀ e:  wt(wAB + ker e) ≠ 5                       [(1,5) dead]
  S2 : ∀ t ∀ e:  wt(wAB + wAB(·-t) + ker e) ≠ 4        [(2,4) dead]
  S3 : ∀ t1 t2 ∀ e:  wt(cand) = 3 → generator column   [(3,3) classified]
  S1': ∀ e:  wt(wBA + ker e) ≠ 5                       [(5,1) dead]
  S2': ∀ t ∀ e:  wt(wBA + wBA(·-t) + ker e) ≠ 4        [(4,2) dead]
(S2/S2'/S3 quantify over ALL t including degenerate values — the
degenerate slots are verified true here, so the Lean statements carry
no side hypotheses.)

The Lean-side certificate Bool (stepOK/pivOK of BaseFloorKernel.lean)
is simulated bit-for-bit (`simulate_pivOK`) so drift fails at emission.

Usage:
  cd experiments/bb_lab
  uv run python scripts/a21_gen_basefloor_data.py [--force]
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
QECLEAN_ROOT = Path(os.environ.get(
    "QECLEAN_ROOT", str(LAB_ROOT.parent.parent / "QECLean"))).resolve()
OUT_PATH = (QECLEAN_ROOT / "QEC" / "Stabilizer" / "Codes" /
            "BivariateBicycle" / "Z5Z15F2A6" / "BaseFloorData.lean")

L, M = 5, 15
N = L * M  # 75

A_SUPP = [(0, 0), (0, 1), (1, 0)]        # 1 + y + x
B_SUPP = [(1, 6), (1, 10), (2, 12)]      # xy^6 + xy^10 + x^2 y^12


def idx(gx: int, gy: int) -> int:
    return (gx % L) * M + (gy % M)


def cell(i: int) -> tuple[int, int]:
    return divmod(i, M)


def sub(c1: tuple[int, int], c2: tuple[int, int]) -> tuple[int, int]:
    return ((c1[0] - c2[0]) % L, (c1[1] - c2[1]) % M)


def add(c1: tuple[int, int], c2: tuple[int, int]) -> tuple[int, int]:
    return ((c1[0] + c2[0]) % L, (c1[1] + c2[1]) % M)


def conv_matrix(supp) -> np.ndarray:
    """C[g,h] = P(g-h): (C f)(g) = (P ⋆ f)(g), repo convention."""
    C = np.zeros((N, N), dtype=np.uint8)
    for g in range(N):
        gx, gy = cell(g)
        for (px, py) in supp:
            C[g, idx(gx - px, gy - py)] ^= 1
    return C


def rref(Min: np.ndarray):
    A = Min.copy().astype(np.uint8)
    rows, cols = A.shape
    E = np.eye(rows, dtype=np.uint8)
    piv: list[int] = []
    r = 0
    for c in range(cols):
        pr = next((rr for rr in range(r, rows) if A[rr, c]), None)
        if pr is None:
            continue
        if pr != r:
            A[[r, pr]] = A[[pr, r]]
            E[[r, pr]] = E[[pr, r]]
        for rr in range(rows):
            if rr != r and A[rr, c]:
                A[rr] ^= A[r]
                E[rr] ^= E[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return A, E, piv


class Solver:
    def __init__(self, Mmat: np.ndarray):
        self.M = Mmat
        self.R, self.E, self.piv = rref(Mmat)
        self.rank = len(self.piv)
        self.cols = Mmat.shape[1]
        free = [c for c in range(self.cols) if c not in set(self.piv)]
        self.free = free
        ker = []
        for fc in free:
            v = np.zeros(self.cols, dtype=np.uint8)
            v[fc] = 1
            for i, pc in enumerate(self.piv):
                if self.R[i, fc]:
                    v[pc] = 1
            ker.append(v)
        self.ker = (np.array(ker, dtype=np.uint8) if ker
                    else np.zeros((0, self.cols), dtype=np.uint8))

    def solve(self, c: np.ndarray):
        y = (self.E @ c.astype(np.int64)) % 2
        if y[self.rank:].any():
            return None
        x = np.zeros(self.cols, dtype=np.uint8)
        x[self.piv] = y[: self.rank].astype(np.uint8)
        return x


def row_comb_certificate(Cmat: np.ndarray, free: list[int]):
    """RREF certificate avoiding the columns in `free`.

    Returns [(pivot cell j, w vector)] with y := Cmat.T @ w satisfying
    y[j] = 1 and y[j'] = 0 for every other pivot cell j'."""
    nonfree = [c for c in range(N) if c not in set(free)]
    perm = nonfree + list(free)
    Rp, Ep, pivp = rref(Cmat[:, perm])
    assert pivp == list(range(len(nonfree))), \
        "free-avoiding RREF did not pivot on every non-free column"
    out = []
    piv_cells = [perm[i] for i in pivp]
    for i in range(len(pivp)):
        j = perm[i]
        w = Ep[i]
        y = (Cmat.T.astype(np.int64) @ w.astype(np.int64)) % 2
        assert y[j] == 1
        for j2 in piv_cells:
            if j2 != j:
                assert y[j2] == 0
        out.append((j, w.astype(np.uint8)))
    return out


def adj3(supp, w: np.ndarray, c: int) -> int:
    """Lean `adj3`: (conv P̃ w)(c) = Σ_{s ∈ supp} w(c + s)."""
    acc = 0
    for s in supp:
        acc ^= int(w[idx(*add(cell(c), s))])
    return acc


def simulate_pivOK(supp, piv: list[tuple[int, np.ndarray]]) -> bool:
    """Bit-for-bit mirror of BaseFloorKernel.stepOK/pivOK."""
    for i, (j, w) in enumerate(piv):
        later = [jj for (jj, _) in piv[i + 1:]]
        if adj3(supp, w, j) != 1:
            return False
        for j2 in later:
            if adj3(supp, w, j2) != 0:
                return False
    return True


def pack(v: np.ndarray) -> int:
    return int(sum(1 << i for i in range(N) if v[i]))


def wt(v: np.ndarray) -> int:
    return int(v.sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if not QECLEAN_ROOT.is_dir():
        sys.exit(f"QECLean checkout not found at {QECLEAN_ROOT} "
                 "(set QECLEAN_ROOT)")

    CA, CB = conv_matrix(A_SUPP), conv_matrix(B_SUPP)
    solA, solB = Solver(CA), Solver(CB)

    # 1. ranks and the shared kernel
    assert solA.rank == 71 and solB.rank == 71, (solA.rank, solB.rank)
    assert solA.ker.shape[0] == 4 and solB.ker.shape[0] == 4
    for z in solA.ker:
        assert not ((CB @ z.astype(np.int64)) % 2).any(), "ker A ⊄ ker B"
    for z in solB.ker:
        assert not ((CA @ z.astype(np.int64)) % 2).any(), "ker B ⊄ ker A"
    print("[1] rank 71 / dim ker 4 both sides; Ann(A) = Ann(B) as spaces")

    kerAll = np.zeros((16, N), dtype=np.uint8)
    for i in range(16):
        v = np.zeros(N, dtype=np.uint8)
        for b in range(4):
            if i >> b & 1:
                v ^= solA.ker[b]
        kerAll[i] = v
    wts = sorted(int(v.sum()) for v in kerAll)
    assert wts == [0] + [40] * 15, wts
    print("[1b] kernel weight table: 0 + 15x40")

    # 2. free cells = RREF free columns of conv_A; echelon basis
    fset = solA.free
    assert len(fset) == 4
    Z = solA.ker  # already δ-normalized on fset by construction
    for i in range(4):
        for j in range(4):
            assert Z[i, fset[j]] == (1 if i == j else 0), "echelon fail"
        assert not ((CA @ Z[i].astype(np.int64)) % 2).any()
        assert not ((CB @ Z[i].astype(np.int64)) % 2).any()
        assert wt(Z[i]) == 40
    print(f"[2] free cells {[cell(f) for f in fset]}; "
          "echelon basis z1..z4 δ-normalized, both-kernel, wt 40")

    # 3. row-combination pivot certificates (both systems, same frees)
    pivA = row_comb_certificate(CA, fset)
    pivB = row_comb_certificate(CB, fset)
    assert len(pivA) == 71 and len(pivB) == 71
    assert simulate_pivOK(A_SUPP, pivA), "Lean pivOK(A) simulation FAILED"
    assert simulate_pivOK(B_SUPP, pivB), "Lean pivOK(B) simulation FAILED"
    assert sorted(j for j, _ in pivA) == sorted(set(range(N)) - set(fset))
    assert sorted(j for j, _ in pivB) == sorted(set(range(N)) - set(fset))
    print("[3] row-combination certificates pass the Lean pivOK simulation")

    # 4. particular solutions
    bvec = np.zeros(N, dtype=np.uint8)
    for s in B_SUPP:
        bvec[idx(*s)] = 1
    avec = np.zeros(N, dtype=np.uint8)
    for s in A_SUPP:
        avec[idx(*s)] = 1
    wAB = solA.solve(bvec)
    wBA = solB.solve(avec)
    assert wAB is not None and wBA is not None
    assert (((CA @ wAB.astype(np.int64)) % 2).astype(np.uint8) == bvec).all()
    assert (((CB @ wBA.astype(np.int64)) % 2).astype(np.uint8) == avec).all()
    print(f"[4] A⋆wAB = B (|wAB| = {wt(wAB)}); B⋆wBA = A (|wBA| = {wt(wBA)})")

    # 5. sweep facts, exact Lean quantifier semantics
    SUB = np.zeros((N, N), dtype=np.int64)   # SUB[g, t] = idx(g - t)
    for g in range(N):
        for t in range(N):
            SUB[g, t] = idx(*sub(cell(g), cell(t)))
    TW_AB = wAB[SUB]     # TW_AB[:, t] = wAB(· - t) as columns
    TW_BA = wBA[SUB]

    # kernel combos, Lean coefficient order (e1·z1 + e2·z2 + e3·z3 + e4·z4)
    K16 = np.zeros((16, N), dtype=np.uint8)
    for e in range(16):
        v = np.zeros(N, dtype=np.uint8)
        for b in range(4):
            if e >> b & 1:
                v ^= Z[b]
        K16[e] = v
    assert sorted(int(K16[e].sum()) for e in range(16)) == [0] + [40] * 15

    # S1 / S1'
    w1 = [int((wAB ^ K16[e]).sum()) for e in range(16)]
    assert all(w != 5 for w in w1), w1
    w1p = [int((wBA ^ K16[e]).sum()) for e in range(16)]
    assert all(w != 5 for w in w1p), w1p
    print(f"[5] S1 profile {sorted(set(w1))}; S1' profile {sorted(set(w1p))}")

    # S2 / S2' over ALL t (including t = 0)
    for t in range(N):
        c2 = wAB ^ TW_AB[:, t]
        ws = {int((c2 ^ K16[e]).sum()) for e in range(16)}
        assert 4 not in ws, (t, ws)
        c2 = wBA ^ TW_BA[:, t]
        ws = {int((c2 ^ K16[e]).sum()) for e in range(16)}
        assert 4 not in ws, (t, ws)
    print("[5b] S2/S2' hold over all 75 t (no side hypotheses needed)")

    # S3 over ALL ordered (t1, t2) incl. degenerate
    fired = []
    for t1 in range(N):
        base3 = wAB ^ TW_AB[:, t1]
        for t2 in range(N):
            c3 = base3 ^ TW_AB[:, t2]
            for e in range(16):
                candv = c3 ^ K16[e]
                if candv.sum() != 3:
                    continue
                S = {0, t1, t2}
                ok = False
                for t in range(N):
                    tc = cell(t)
                    At = {idx(*add(s, tc)) for s in A_SUPP}
                    Bt = np.zeros(N, dtype=np.uint8)
                    for s in B_SUPP:
                        Bt[idx(*add(s, tc))] = 1
                    if At == S and (Bt == candv).all():
                        ok = True
                        break
                assert ok, f"S3 FALSIFIED at t1={t1}, t2={t2}, e={e}"
                fired.append((t1, t2, e))
    assert len(fired) == 6, fired   # 3 translates x 2 orders, 1 e each
    print(f"[5c] S3: {len(fired)} fired class/e slots, all generator columns")

    # 6. sanity: generator columns are weight-6 cycles
    for g in range(N):
        colL = CA[:, g]
        colR = CB[:, g]
        assert colL.sum() == 3 and colR.sum() == 3
        s = (CB @ colL.astype(np.int64) + CA @ colR.astype(np.int64)) % 2
        assert not s.any()
    print("[6] all 75 generator columns are weight-6 cycles")

    # ── emit ─────────────────────────────────────────────────────────
    if OUT_PATH.exists() and not args.force:
        sys.exit(f"REFUSING to overwrite {OUT_PATH} (use --force)")

    def lean_cell(i: int) -> str:
        x, y = cell(i)
        return f"({x}, {y})"

    def lean_piv_list(piv, name: str, doc: str) -> str:
        items = [f"({lean_cell(j)}, {pack(w)})" for (j, w) in piv]
        body = ",\n    ".join(items)
        return (f"/-- {doc} -/\ndef {name} : List (G150 × ℕ) :=\n"
                f"  [{body}]\n")

    docA = ("Row-combination pivot certificate for `conv a150`: entries\n"
            "`(pivot cell j, mask w)` with `conv Ã w` equal to `1` at `j`\n"
            "and `0` at every later pivot cell.")
    docB = "Row-combination pivot certificate for `conv b150`."

    lines: list[str] = []
    lines.append(f"""/-
GENERATED FILE — DO NOT HAND-EDIT.
Generator: qec-lab:experiments/bb_lab/scripts/a21_gen_basefloor_data.py
Data: computed in-script (repo-convention conv matrices for the
Z5Z15F2A6 polynomials); every table is numpy-hard-asserted at emission.
Regen: cd experiments/bb_lab && uv run python scripts/a21_gen_basefloor_data.py --force
-/

import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.Defs

/-!
# Z5Z15F2A6 base-floor coset data (A21 weight-6 layer)

Data tables for the translation-reduced coset sweep that discharges
`weight6_cycle_is_boundary`: particular solutions `wAB` (`A⋆wAB = B`)
and `wBA` (`B⋆wBA = A`), the echelon basis `z1..z4` of the shared
4-dimensional kernel `Ann(A) = Ann(B)` (δ-normalized on the free cells
`fc1..fc4`), and row-combination pivot certificates
`pivListA`/`pivListB` (pairs `(pivot cell, row-combination mask)`)
certifying both kernels — consumed by `BaseFloorKernel.lean`.

Masks pack a chain `v : G150 → ZMod 2` as `Σ 2^(cell0Idx g) · v g` with
`cell0Idx (x, y) = x·15 + y`.
-/

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Z5Z15F2A6

/-- Packed-mask bit position of a cell: `(x, y) ↦ x·15 + y`. -/
def cell0Idx (g : G150) : ℕ := g.1.val * 15 + g.2.val

/-- Chain of a packed mask. -/
def maskFun (m : ℕ) : G150 → ZMod 2 :=
  fun g => if m.testBit (cell0Idx g) then 1 else 0

/-- Particular solution `A ⋆ wAB = B` (packed; weight {wt(wAB)}). -/
def wABmask : ℕ := {pack(wAB)}

/-- `wAB` as a chain. -/
def wABf : G150 → ZMod 2 := maskFun wABmask

/-- Particular solution `B ⋆ wBA = A` (packed; weight {wt(wBA)}). -/
def wBAmask : ℕ := {pack(wBA)}

/-- `wBA` as a chain. -/
def wBAf : G150 → ZMod 2 := maskFun wBAmask
""")
    for i in range(4):
        lines.append(f"""/-- Kernel basis element `z{i + 1}` (packed; weight 40). -/
def zmask{i + 1} : ℕ := {pack(Z[i])}

/-- `z{i + 1}` as a chain. -/
def zf{i + 1} : G150 → ZMod 2 := maskFun zmask{i + 1}
""")
    fcs = [lean_cell(f) for f in fset]
    lines.append(f"""/-- Free cell 1 (kernel coordinate 1). -/
def fc1 : G150 := {fcs[0]}

/-- Free cell 2. -/
def fc2 : G150 := {fcs[1]}

/-- Free cell 3. -/
def fc3 : G150 := {fcs[2]}

/-- Free cell 4. -/
def fc4 : G150 := {fcs[3]}

/-- The 16-element kernel span in free-cell coordinates:
`kerElt e1 e2 e3 e4 = e1·z1 + e2·z2 + e3·z3 + e4·z4`
(δ-normalization gives `kerElt e1 e2 e3 e4 fcᵢ = eᵢ`). -/
def kerElt (e1 e2 e3 e4 : ZMod 2) : G150 → ZMod 2 :=
  fun g => e1 * zf1 g + e2 * zf2 g + e3 * zf3 g + e4 * zf4 g
""")
    lines.append(lean_piv_list(pivA, "pivListA", docA))
    lines.append("")
    lines.append(lean_piv_list(pivB, "pivListB", docB))
    lines.append("""
end Z5Z15F2A6
end BB
end Homological
end Stabilizer
end Quantum
""")

    OUT_PATH.write_text("\n".join(lines))
    print(f"\nwrote {OUT_PATH}")


if __name__ == "__main__":
    main()
