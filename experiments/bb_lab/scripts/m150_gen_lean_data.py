"""A32/M0-M2 emitter: Lean data for the [[150,30,10]] mitten certification.

v0 scope = the **budget rehearsal** (plan.md §M0.4): emit a self-contained
Lean file with one representative obligation of each planned kind, sized
exactly like the real M1/M2 obligations, so its measured `lake env lean`
wall time calibrates the build_budget.md allocations BEFORE the M4 case
machine is designed. Obligations emitted:

  R1 (dictionary, native_decide): the GAP index → carrier map is a group
     hom — all 900 products against the GAP table;
  R2 (M1/M2 shape, native_decide): closed-form H entries — all 60 H_Z
     rows ∈ ker H_X (this is H_X·H_Zᵀ = 0 restated) + all rows weight 9;
  R3 (certificate shape, native_decide): Gaussian left-inverse pivot
     certificate W·H_X[:,P] = I₆₀ from emitted packed rows — the
     KernelCert-style rank fact all M2/M4 certificates will take;
  R4 (kernel decide probe): gapElem injectivity — 900 carrier
     comparisons through the kernel evaluator, calibrating the
     decide-vs-native crossover at this group's size.

Carrier: `Multiplicative (ZMod 5) × DihedralGroup 3` (ZMod 5 alone is
additive — the wrapper makes C₅ multiplicative so the product is a
mathlib `Group`).

Falsify-first validation before any emission:
  - Python model of mathlib's DihedralGroup multiplication; the
    dictionary map is verified to be a bijective hom in Python (900
    products) — and the *mathlib orientation* (the sr-index sign) is
    additionally verified in-Lean by R1 itself: a sign error emits a
    false proposition and native_decide fails the build, it cannot
    silently pass;
  - the closed-form H_X/H_Z entry formulas are checked entry-for-entry
    (2 × 60 × 150) against the matrix builder of record
    (a26_mitten_descent.mitten_code);
  - the pivot certificate W·H_X[:,P] = I is verified in numpy GF(2).

Output goes OUTSIDE the QECLean tree (default: qec-lab scratch path or
--out). The eventual M2 production emitter grows out of this file; the
emitted banner names this script per GENERATORS.md.

Run (from experiments/bb_lab):
  uv run python scripts/m150_gen_lean_data.py rehearsal --out /tmp/dir
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
GROUPS = LAB_ROOT / "instances" / "mitten_groups"
SETS = dict(a0=(0, 14, 23), a1=(0, 2, 11), b0=(7, 20, 24), b1=(0, 2, 29))
BANNER = (
    "/-\nGENERATED FILE — DO NOT HAND-EDIT.\n"
    "Emitted by qec-lab:experiments/bb_lab/scripts/m150_gen_lean_data.py "
    "(mode: rehearsal)\nfrom instances/mitten_groups/group_30_1.txt + "
    "arXiv:2607.28795 Table XIII sets.\nRegen: uv run python "
    "scripts/m150_gen_lean_data.py rehearsal --out <dir> --force\n"
    "Budget-rehearsal probe (pipeline/attempts/mitten_150_30_10/"
    "plan.md §M0.4) — NOT a library file.\n-/"
)


def _load_a26():
    spec = importlib.util.spec_from_file_location(
        "a26", Path(__file__).parent / "a26_mitten_descent.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.MITTEN["150,30,10"] == dict(gid=(30, 1), **SETS)
    return mod


# ---------------------------------------------------- dictionary layer


def build_dictionary(G):
    """GAP index → (i, j, k) with g = z^i r^j s^k, plus the mathlib
    carrier coordinates (i, dihedral) under the orientation
    r^j ↦ DihedralGroup.r j, r^j s ↦ DihedralGroup.sr (-j)."""
    n = G.n
    z = next(g for g in G.center() if G.orders[g] == 5)
    r = next(g for g in range(n) if int(G.orders[g]) == 3)
    ri = int(G.inv[r])
    s = next(g for g in range(n)
             if int(G.orders[g]) == 2 and G.mul(G.mul(g, r), g) == ri)
    table = {}
    zi = 0
    for i in range(5):
        for j in range(3):
            for k in range(2):
                el = zi
                for _ in range(j):
                    el = G.mul(el, r)
                if k:
                    el = G.mul(el, s)
                assert el not in table
                table[el] = (i, j, k)
        zi = G.mul(zi, z)
    assert len(table) == n

    # python model of mathlib DihedralGroup 3 multiplication
    # element = ('r'|'sr', idx mod 3);  ZMod-5 part is additive.
    def dmul(p, q):
        (t1, i1), (t2, i2) = p, q
        if t1 == "r" and t2 == "r":
            return ("r", (i1 + i2) % 3)
        if t1 == "r" and t2 == "sr":
            return ("sr", (i2 - i1) % 3)
        if t1 == "sr" and t2 == "r":
            return ("sr", (i1 + i2) % 3)
        return ("r", (i2 - i1) % 3)

    def carrier(ijk):
        i, j, k = ijk
        return (i, ("r", j) if k == 0 else ("sr", (-j) % 3))

    def cmul(p, q):
        return ((p[0] + q[0]) % 5, dmul(p[1], q[1]))

    for a in range(n):
        for b in range(n):
            assert cmul(carrier(table[a]), carrier(table[b])) == carrier(
                table[G.mul(a, b)]
            ), f"carrier-hom failure at {a},{b} — sr-orientation wrong?"
    return table, carrier, (z, r, s)


def lean_elem(c) -> str:
    i, (t, j) = c
    d = f"DihedralGroup.r {j}" if t == "r" else f"DihedralGroup.sr {j}"
    return f"(Multiplicative.ofAdd ({i} : ZMod 5), {d})"


# ---------------------------------------------- closed-form H entries


def closed_form_H(G):
    """H_X/H_Z entries from the Eq. (J1) closed forms, in GAP-index
    space; validated entry-for-entry against a26.mitten_code."""
    n = G.n
    inv = G.inv

    def hx(beta, h, m, x):
        if m == 4:
            return int(G.mul(int(inv[x]), h) in SETS[f"b{beta}"])
        if m % 2 != beta:
            return 0
        return int(G.mul(h, int(inv[x])) in SETS[f"a{m // 2}"])

    def hz(alpha, h, m, x):
        if m == 4:
            return int(G.mul(x, int(inv[h])) in SETS[f"a{alpha}"])
        if m // 2 != alpha:
            return 0
        return int(G.mul(int(inv[h]), x) in SETS[f"b{m % 2}"])

    HX = np.zeros((60, 150), dtype=np.uint8)
    HZ = np.zeros((60, 150), dtype=np.uint8)
    for i in range(2):
        for h in range(n):
            for m in range(5):
                for x in range(n):
                    HX[i * n + h, m * n + x] = hx(i, h, m, x)
                    HZ[i * n + h, m * n + x] = hz(i, h, m, x)
    return HX, HZ


def pivot_certificate(H: np.ndarray):
    """Pivot columns P (|P| = 60) and W = (H[:,P])⁻¹ over GF(2), with
    W·H[:,P] = I verified."""
    m = H.shape[0]
    A = H.copy().astype(np.uint8)
    E = np.eye(m, dtype=np.uint8)
    piv = []
    r = 0
    for c in range(H.shape[1]):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        A[[r, p]] = A[[p, r]]
        E[[r, p]] = E[[p, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
                E[i] ^= E[r]
        piv.append(c)
        r += 1
        if r == m:
            break
    assert r == m, f"rank {r} < {m}"
    W = E  # row-reduces H to RREF; on pivot columns RREF = I ⟹ W·H[:,P]=I
    assert (W.astype(int) @ H[:, piv] % 2 == np.eye(m, dtype=int)).all()
    return piv, W


# -------------------------------------------------------------- emit


def fmt_list(items, per_line=6, indent="  "):
    lines = []
    for i in range(0, len(items), per_line):
        lines.append(indent + ", ".join(items[i : i + per_line])
                     + ("," if i + per_line < len(items) else ""))
    return "[\n" + "\n".join(lines) + "]"


def emit_rehearsal(out_dir: Path, force: bool) -> None:
    a26 = _load_a26()
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    table, carrier, (z, r, s) = build_dictionary(G)
    print(f"[gen] dictionary validated (z={z}, r={r}, s={s}; carrier hom "
          f"checked on all 900 products incl. mathlib orientation model)")

    HX, HZ = closed_form_H(G)
    HXref, HZref = a26.mitten_code(G, **SETS)
    assert (HX == HXref).all() and (HZ == HZref).all(), (
        "closed-form H entries disagree with a26.mitten_code")
    print("[gen] closed-form H_X/H_Z validated entry-for-entry vs "
          "a26.mitten_code (2×60×150)")

    piv, W = pivot_certificate(HX)
    print(f"[gen] pivot certificate: |P| = {len(piv)}, W·H_X[:,P] = I "
          f"verified in GF(2)")

    elems = [lean_elem(carrier(table[g])) for g in range(30)]
    mulflat = [int(G.mt[a, b]) for a in range(30) for b in range(30)]
    set_lines = {
        k: fmt_list([lean_elem(carrier(table[g])) for g in v], per_line=1)
        for k, v in SETS.items()
    }
    piv_pairs = [f"(({c // 30} : Fin 5), {c % 30})" for c in piv]  # (block, gapIdx)
    wrows = [int("".join(str(b) for b in row[::-1]), 2) for row in W]

    lean = f"""{BANNER}
import Mathlib.GroupTheory.SpecificGroups.Dihedral
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Fintype.Prod

namespace A32Rehearsal

/-- The mitten group carrier: C₅ (multiplicative) × S₃ (= Dih₃). -/
abbrev M150G : Type := Multiplicative (ZMod 5) × DihedralGroup 3

/-- GAP `Elements(SmallGroup(30,1))` order → carrier, via z^i·r^j·s^k. -/
def gapElems : List M150G := {fmt_list(elems, per_line=1)}

def gapElem (i : Nat) : M150G := gapElems.getD i 1

/-- GAP multiplication table, row-major (30×30). -/
def mulFlat : List Nat := {fmt_list([str(v) for v in mulflat], per_line=20)}

/-- The four Table XIII sets, as carrier-element lists. -/
def a0 : List M150G := {set_lines["a0"]}
def a1 : List M150G := {set_lines["a1"]}
def b0 : List M150G := {set_lines["b0"]}
def b1 : List M150G := {set_lines["b1"]}

/-! ## R1 — dictionary is a hom (native_decide, 900 products) -/

theorem dict_hom : ∀ a : Fin 30, ∀ b : Fin 30,
    gapElem (mulFlat.getD (30 * a.val + b.val) 0)
      = gapElem a.val * gapElem b.val := by
  native_decide

/-! ## R2 — closed-form H entries; rows ∈ kernel + weight (M1/M2 shape) -/

def asets : List (List M150G) := [a0, a1]
def bsets : List (List M150G) := [b0, b1]

/-- `H_X` entry at X-check `(β, h)`, qubit `(m, x)` (Eq. (J1)). -/
def HXval (β : Fin 2) (h : M150G) (m : Fin 5) (x : M150G) : ZMod 2 :=
  if m.val = 4 then (if x⁻¹ * h ∈ bsets.getD β.val [] then 1 else 0)
  else if m.val % 2 = β.val then
    (if h * x⁻¹ ∈ asets.getD (m.val / 2) [] then 1 else 0)
  else 0

/-- `H_Z` entry at Z-check `(α, h)`, qubit `(m, x)`. -/
def HZval (α : Fin 2) (h : M150G) (m : Fin 5) (x : M150G) : ZMod 2 :=
  if m.val = 4 then (if x * h⁻¹ ∈ asets.getD α.val [] then 1 else 0)
  else if m.val / 2 = α.val then
    (if h⁻¹ * x ∈ bsets.getD (m.val % 2) [] then 1 else 0)
  else 0

/-- Every `H_Z` row is a cycle: `H_X · H_Zᵀ = 0` restated
(60 × 60 pairings × 150-qubit sums). -/
theorem hzRows_in_kerHX : ∀ α : Fin 2, ∀ β : Fin 2, ∀ h : Fin 30, ∀ h' : Fin 30,
    (∑ m : Fin 5, ∑ x : Fin 30,
      HXval β (gapElem h'.val) m (gapElem x.val)
        * HZval α (gapElem h.val) m (gapElem x.val)) = 0 := by
  native_decide

/-- Every `H_Z` row has weight exactly 9. -/
theorem hzRows_weight9 : ∀ α : Fin 2, ∀ h : Fin 30,
    (Finset.univ.filter fun q : Fin 5 × Fin 30 =>
      HZval α (gapElem h.val) q.1 (gapElem q.2.val) ≠ 0).card = 9 := by
  native_decide

/-! ## R3 — Gaussian left-inverse pivot certificate (rank H_X = 60) -/

/-- Pivot qubits, as `(block, GAP index)` pairs. -/
def pivots : List (Fin 5 × Nat) := {fmt_list(piv_pairs, per_line=6)}

/-- Rows of `W = (H_X[:,P])⁻¹`, packed little-endian in 60-bit Nats. -/
def wRows : List Nat := {fmt_list([str(v) for v in wrows], per_line=4)}

def checkIdx (k : Nat) : Fin 2 × Nat :=
  ((if k < 30 then (0 : Fin 2) else 1), k % 30)

/-- `W · H_X[:,P] = I₆₀` — the full-rank certificate. -/
theorem pivot_cert : ∀ i : Fin 60, ∀ j : Fin 60,
    (∑ k : Fin 60,
      (if ((wRows.getD i.val 0).testBit k.val) then (1 : ZMod 2) else 0)
        * (HXval (checkIdx k.val).1 (gapElem (checkIdx k.val).2)
            (pivots.getD j.val (0, 0)).1
            (gapElem (pivots.getD j.val (0, 0)).2)))
      = if i = j then 1 else 0 := by
  native_decide

/-! ## R4 — kernel-decide probe: dictionary injectivity (900 pairs) -/

theorem gapElem_inj : ∀ a : Fin 30, ∀ b : Fin 30,
    gapElem a.val = gapElem b.val → a = b := by
  decide

end A32Rehearsal
"""
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "A32Rehearsal.lean"
    if target.exists() and not force:
        sys.exit(f"refusing to overwrite {target} (use --force)")
    target.write_text(lean)
    baseline = out_dir / "A32Baseline.lean"
    baseline.write_text(
        f"{BANNER}\nimport Mathlib.GroupTheory.SpecificGroups.Dihedral\n"
        "import Mathlib.Data.ZMod.Basic\nimport Mathlib.Data.Fintype.Prod\n\n"
        "namespace A32Baseline\ntheorem trivial_probe : 1 + 1 = 2 := rfl\n"
        "end A32Baseline\n"
    )
    print(f"[gen] wrote {target} and {baseline}")
    print("[gen] time them from the QECLean checkout root:")
    print("        time lake env lean <dir>/A32Baseline.lean")
    print("        time lake env lean <dir>/A32Rehearsal.lean")


# --------------------------------------------------- instance (M2) mode


def rref_rows(M: np.ndarray):
    """(pivot column list, W) with W·M[:,piv] = I over GF(2), asserting
    full row rank (same routine as `pivot_certificate`, kept separate so
    the M2 lane is self-auditing)."""
    m = M.shape[0]
    A = M.copy().astype(np.uint8)
    E = np.eye(m, dtype=np.uint8)
    piv, r = [], 0
    for c in range(M.shape[1]):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        A[[r, p]] = A[[p, r]]
        E[[r, p]] = E[[p, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
                E[i] ^= E[r]
        piv.append(c)
        r += 1
        if r == m:
            break
    assert r == m, f"rank {r} < {m}"
    assert (E.astype(int) @ M[:, piv] % 2 == np.eye(m, dtype=int)).all()
    return piv, E


def nullspace(M: np.ndarray) -> np.ndarray:
    """Basis of ker over GF(2), rows = basis vectors."""
    m, n = M.shape
    A = M.copy().astype(np.uint8)
    piv, r = [], 0
    for c in range(n):
        if r >= m:
            break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        A[[r, p]] = A[[p, r]]
        for i in range(m):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        piv.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv]
    basis = np.zeros((len(free), n), dtype=np.uint8)
    for k, c in enumerate(free):
        basis[k, c] = 1
        for i, pc in enumerate(piv):
            basis[k, pc] = A[i, c]
    assert not (M.astype(int) @ basis.T % 2).any()
    return basis


def symplectic_basis(HX, HZ):
    """(Lx, Lz) with rows Lx ⊂ ker H_Z, Lz ⊂ ker H_X, Lx·Lzᵀ = I₃₀,
    by greedy symplectic Gram–Schmidt over the kernel bases."""
    KZ = [v.copy() for v in nullspace(HX)]  # candidate Z-logicals
    KX = [v.copy() for v in nullspace(HZ)]  # candidate X-logicals
    Lx, Lz = [], []
    while len(Lz) < 30:
        found = None
        for iz, z in enumerate(KZ):
            for ix, x in enumerate(KX):
                if int(x @ z) % 2 == 1:
                    found = (iz, ix)
                    break
            if found:
                break
        assert found, "symplectic pairing exhausted early"
        iz, ix = found
        z, x = KZ.pop(iz), KX.pop(ix)
        KZ = [w ^ z if int(x @ w) % 2 else w for w in KZ]
        KX = [w ^ x if int(w @ z) % 2 else w for w in KX]
        Lz.append(z)
        Lx.append(x)
    Lx, Lz = np.array(Lx, dtype=np.uint8), np.array(Lz, dtype=np.uint8)
    assert not (HZ.astype(int) @ Lx.T % 2).any()
    assert not (HX.astype(int) @ Lz.T % 2).any()
    assert (Lx.astype(int) @ Lz.T % 2 == np.eye(30, dtype=int)).all()
    return Lx, Lz


def weight10_witness(HX, HZ, KX):
    """A weight-10 v ∈ ker H_X with v ∉ rowspace(H_Z), plus a pairing
    x ∈ ker H_Z with ⟨x, v⟩ = 1.  CMS at exact weight 10 (the census
    proved no lighter kernel word exists beyond the weight-9 rows, and
    the anticommutation clause excludes those)."""
    from bb_lab.sat_distance import _solve_at_weight_cms

    v, _ = _solve_at_weight_cms(HX, KX, 10)
    assert v is not None and int(v.sum()) == 10
    assert not (HX.astype(int) @ v % 2).any()
    # not a boundary: rank jumps when stacked on H_Z
    piv_hz, _ = rref_rows(HZ)
    stacked = np.vstack([HZ, v[None, :]])
    A = stacked.copy().astype(np.uint8)
    r = 0
    for c in range(A.shape[1]):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        A[[r, p]] = A[[p, r]]
        for i in range(A.shape[0]):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    assert r == 61, "witness lies in rowspace(H_Z)"
    pair = next(x for x in KX if int(x @ v) % 2 == 1)
    assert not (HZ.astype(int) @ pair % 2).any()
    return v, pair


def emit_instance(out_dir: Path, force: bool) -> None:
    """M2 data feed: canonical-labeling decoder/logical/witness data for
    `QEC/Stabilizer/Codes/Mitten/M150/`, all validated before emission."""
    a26 = _load_a26()
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    table, carrier, (z, r, s) = build_dictionary(G)
    HX, HZ = closed_form_H(G)
    HXref, HZref = a26.mitten_code(G, **SETS)
    assert (HX == HXref).all() and (HZ == HZref).all()

    pivX, WX = rref_rows(HX)
    pivZ, WZ = rref_rows(HZ)
    # decoder form: phi·(∂ basis matrix) = I with phi supported on pivots;
    # equivalently W·H[:,piv] = I plus square-inverse symmetry — validate
    # the exact identity the Lean file will check: for all i, p:
    # Σ_j W[i,j]·H[p, piv_j] = [i = p] (H[:,piv]·W = I, square case).
    assert (HX[:, pivX].astype(int) @ WX % 2 == np.eye(60, dtype=int)).all()
    assert (HZ[:, pivZ].astype(int) @ WZ % 2 == np.eye(60, dtype=int)).all()

    Lx, Lz = symplectic_basis(HX, HZ)
    KX = nullspace(HZ)
    wit, witPair = weight10_witness(HX, HZ, KX)
    print(f"[gen] instance data validated: pivots/W both sides, "
          f"symplectic 30-basis (Lx wts {sorted(set(Lx.sum(1)))[:3]}…, "
          f"Lz wts {sorted(set(Lz.sum(1)))[:3]}…), witness wt "
          f"{int(wit.sum())} + pairing wt {int(witPair.sum())}")

    elems = [lean_elem(carrier(table[g])) for g in range(30)]
    sets_lean = {k: [lean_elem(carrier(table[g])) for g in v]
                 for k, v in SETS.items()}
    packedW = lambda W: [
        int("".join(str(b) for b in row[::-1]), 2) for row in W]
    sup = lambda v: [int(i) for i in np.flatnonzero(v)]

    def fmt_sup_rows(rows) -> str:
        # nested list-of-lists, every line ≤ 100 chars (repo linter)
        out = []
        for k, row in enumerate(rows):
            body = fmt_list([str(i) for i in sup(row)], per_line=14,
                            indent="   ")
            out.append("  " + body + ("," if k + 1 < len(rows) else ""))
        return "[\n" + "\n".join(out) + "]"

    instance_banner = (
        "/-\nGENERATED FILE — DO NOT HAND-EDIT.\n"
        "Emitted by qec-lab:experiments/bb_lab/scripts/m150_gen_lean_data.py "
        "(mode: instance)\nfrom instances/mitten_groups/group_30_1.txt + "
        "arXiv:2607.28795 Table XIII sets;\nall facts validated in numpy "
        "before emission (dictionary hom, closed-form H vs\n"
        "a26_mitten_descent.mitten_code, pivot inverses, symplectic basis, "
        "witness).\nRegen: uv run python scripts/m150_gen_lean_data.py "
        "instance --out <M150 dir> --force\n"
        "Attempt state: qec-lab:pipeline/attempts/mitten_150_30_10/.\n-/"
    )
    data = f"""{instance_banner}
import QEC.Stabilizer.Framework.Homological.LiftedProduct
import Mathlib.GroupTheory.SpecificGroups.Dihedral

namespace Quantum
namespace Stabilizer
namespace Homological
namespace LP
namespace M150

/-- The `[[150,30,10]]` mitten group carrier: C₅ (multiplicative) × S₃. -/
abbrev M150G : Type := Multiplicative (ZMod 5) × DihedralGroup 3

/-- GAP `Elements(SmallGroup(30,1))` order → carrier (z^i·r^j·s^k
parameterization; z = idx 2, r = idx 3, s = idx 1). -/
def gapElems : List M150G := {fmt_list(elems, per_line=1)}

/-- Carrier element of a GAP index (junk-total via identity). -/
def elemOf (i : Nat) : M150G := gapElems.getD i 1

/-- Qubit cell of a canonical qubit index `30·m + g`. -/
def qubitOf (c : Nat) : Fin 5 × M150G :=
  (⟨(c / 30) % 5, Nat.mod_lt _ (by omega)⟩, elemOf (c % 30))

/-- Check cell of a canonical check index `30·i + g`. -/
def checkOf (k : Nat) : Fin 2 × M150G :=
  (⟨(k / 30) % 2, Nat.mod_lt _ (by omega)⟩, elemOf (k % 30))

/-- Table XIII sets (paper order a0, a1, b0, b1). -/
def a0 : List M150G := {fmt_list(sets_lean["a0"], per_line=1)}
def a1 : List M150G := {fmt_list(sets_lean["a1"], per_line=1)}
def b0 : List M150G := {fmt_list(sets_lean["b0"], per_line=1)}
def b1 : List M150G := {fmt_list(sets_lean["b1"], per_line=1)}

/-- Pivot qubit indices for `H_X` (canonical `30·m + g`). -/
def pivX : List Nat := {fmt_list([str(c) for c in pivX], per_line=15)}

/-- Pivot qubit indices for `H_Z`. -/
def pivZ : List Nat := {fmt_list([str(c) for c in pivZ], per_line=15)}

/-- Rows of `(H_X[:,pivX])⁻¹`, packed little-endian 60-bit Nats. -/
def wX : List Nat := {fmt_list([str(v) for v in packedW(WX)], per_line=4)}

/-- Rows of `(H_Z[:,pivZ])⁻¹`. -/
def wZ : List Nat := {fmt_list([str(v) for v in packedW(WZ)], per_line=4)}

/-- Supports of the 30 logical-Z chains (rows of `Lz`, ⊂ ker H_X). -/
def logZsup : List (List Nat) := {fmt_sup_rows(Lz)}

/-- Supports of the 30 logical-X chains (rows of `Lx`, ⊂ ker H_Z);
`Lx·Lzᵀ = I₃₀` (validated offline, re-checked in Lean). -/
def logXsup : List (List Nat) := {fmt_sup_rows(Lx)}

/-- Support of the weight-10 distance witness (∈ ker H_X, ∉ rowspace H_Z). -/
def witSup : List Nat := {fmt_list([str(i) for i in sup(wit)], per_line=14)}

/-- Support of its dual pairing (∈ ker H_Z, ⟨·, wit⟩ = 1). -/
def witPairSup : List Nat := {fmt_list([str(i) for i in sup(witPair)], per_line=14)}

end M150
end LP
end Homological
end Stabilizer
end Quantum
"""
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "Data.lean"
    if target.exists() and not force:
        sys.exit(f"refusing to overwrite {target} (use --force)")
    target.write_text(data)
    print(f"[gen] wrote {target}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["rehearsal", "instance"])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.mode == "rehearsal":
        emit_rehearsal(args.out, args.force)
    else:
        emit_instance(args.out, args.force)


if __name__ == "__main__":
    main()
