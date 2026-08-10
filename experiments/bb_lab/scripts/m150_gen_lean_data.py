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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["rehearsal"])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    emit_rehearsal(args.out, args.force)


if __name__ == "__main__":
    main()
