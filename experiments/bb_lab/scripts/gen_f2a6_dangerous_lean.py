"""Emit the Z5Z15F2A6 dangerous-dispatch data + certificate modules.

Stage 2 of the (M) kernel route (A17 line): from the exhaustive
light-boundary classification (data/a17/f2a6_light_classes.jsonl, 113
translation classes = 94 small + 19 near-kernel), emit into the QECLean
checkout (env QECLEAN_ROOT, default ../QECLean):

  QEC/Stabilizer/Codes/BivariateBicycle/Z5Z15F2A6/ClassData.lean
    packed-Nat data: class reps, small-class (shift, seam-good f0)
    bundles, near-kernel window datasets (min-poke f0, window mask,
    dim-1 cycle table + preimages, t>=2 extra-cell cycles), the global
    d1-column table, and the testBit accessors.
  QEC/Stabilizer/Codes/BivariateBicycle/Z5Z15F2A6/CertSweep.lean
    the batched kernel certificates (native_decide) tying every table
    entry to the semantic objects: translate/seam/window/column/parity.

Everything emitted is hard-asserted in numpy first (repo convention);
a build failure downstream means drift, not silent wrongness.

Cell index convention (matches the lab and the Lean accessors):
  2-chain cell g = (x, y)        -> x*15 + y            (0..74)
  1-chain cell (g, blk)          -> blk*75 + x*15 + y   (0..149)

Usage:
  cd experiments/bb_lab
  uv run python scripts/gen_f2a6_dangerous_lean.py [--force]
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB_ROOT / "src"))

from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly
from bb_lab.checks import circulant
from bb_lab.linalg import nullspace_f2, rank_f2, rref_f2

A_STR, B_STR = "1 + y + x", "x*y^6 + x*y^10 + x^2*y^12"
ELL, M = 5, 15
N_CLASSES = 113

QECLEAN_ROOT = Path(os.environ.get("QECLEAN_ROOT",
                                   LAB_ROOT.parent.parent / ".." / "QECLean")).resolve()
OUT_DIR = QECLEAN_ROOT / "QEC" / "Stabilizer" / "Codes" / "BivariateBicycle" / "Z5Z15F2A6"

BANNER = """/-
GENERATED FILE — DO NOT HAND-EDIT.
Generator: qec-lab:experiments/bb_lab/scripts/gen_f2a6_dangerous_lean.py
Data source: qec-lab:experiments/bb_lab/data/a17/f2a6_light_classes.jsonl
Regen: cd experiments/bb_lab && uv run python scripts/gen_f2a6_dangerous_lean.py --force
-/
"""

# ---------------------------------------------------------------- setup
Gb = AbelianGroup((ELL, M))
Gc = AbelianGroup((ELL, 2 * M))
nb, nc = Gb.cardinality, Gc.cardinality
elems_b, elems_c = list(Gb), list(Gc)
base_idx = {g: i for i, g in enumerate(Gb)}
cover_idx = {g: i for i, g in enumerate(Gc)}
Ab, Bb = Poly.from_string(A_STR, Gb), Poly.from_string(B_STR, Gb)
Ac = Poly.from_support(Ab.support, Gc)
Bc = Poly.from_support(Bb.support, Gc)
MAb, MBb = circulant(Ab).astype(np.uint8), circulant(Bb).astype(np.uint8)
MAc, MBc = circulant(Ac).astype(np.uint8), circulant(Bc).astype(np.uint8)
D2b = np.vstack([MAb, MBb]) % 2
D1b = np.hstack([MBb, MAb]) % 2
D2c = np.vstack([MAc, MBc]) % 2

kerb = nullspace_f2(D2b).astype(np.uint8)
assert kerb.shape[0] == 4
ker_elems = []
for mask in range(16):
    z = np.zeros(nb, dtype=np.uint8)
    for i in range(4):
        if (mask >> i) & 1:
            z ^= kerb[i]
    ker_elems.append(z)

LIFT_COL = np.zeros((nc, nb), dtype=np.uint8)
for i, (x, y) in enumerate(elems_b):
    LIFT_COL[cover_idx[(x, y)], i] = 1
D2C_LIFT = (D2c @ LIFT_COL) % 2


def d2b(f):
    return (D2b @ f) % 2


def solve_f2(Amat, b):
    aug = np.hstack([Amat.astype(np.uint8) % 2,
                     (b.astype(np.uint8) % 2)[:, None]])
    R, piv = rref_f2(aug)
    ncols = Amat.shape[1]
    if ncols in piv:
        return None
    x = np.zeros(ncols, dtype=np.uint8)
    for r, c in enumerate(piv):
        x[c] = R[r, ncols]
    return x


def seam_of(f):
    """sheet0 of the lifted stabilizer, base coordinates (150-vector)."""
    L = (D2C_LIFT @ f) % 2
    s = np.zeros(2 * nb, dtype=np.uint8)
    for blk in range(2):
        for h in elems_b:
            s[blk * nb + base_idx[h]] = L[blk * nc + cover_idx[(h[0], h[1])]]
    return s


def pokes(f):
    return (int(((D2C_LIFT @ f) % 2).sum()) - int(d2b(f).sum())) // 2


def seam_good_coset(f):
    for z in ker_elems:
        f0 = (f ^ z) % 2
        if pokes(f0) == 0:
            return f0
    return None


def min_poke_coset(f):
    best = None
    for z in ker_elems:
        f0 = (f ^ z) % 2
        p = pokes(f0)
        if best is None or p < best[0]:
            best = (p, f0)
    return best


def translate_chain(b, tx, ty):
    """Lean translate1 convention: (translate1 c v) p = v (p.1 + c, p.2),
    so the translated chain at cell (g, blk) reads the original at
    (g + c, blk)."""
    out = np.zeros_like(b)
    for blk in range(2):
        for i, (gx, gy) in enumerate(elems_b):
            out[blk * nb + i] = b[blk * nb
                                  + base_idx[((gx + tx) % ELL, (gy + ty) % M)]]
    return out


def pack_bits(v) -> int:
    n = 0
    for i in np.nonzero(v)[0]:
        n |= 1 << int(i)
    return n


# ---------------------------------------------------------------- load
IN = LAB_ROOT / "data" / "a17" / "f2a6_light_classes.jsonl"
recs = []
complete = False
with open(IN) as fh:
    for line in fh:
        r = json.loads(line)
        if "complete" in r:
            complete = r["complete"]
        else:
            recs.append(r)
assert complete and len(recs) == N_CLASSES, (complete, len(recs))
print(f"loaded {len(recs)} classes (complete enumeration)")

# ---------------------------------------------------------------- analyze
KIND = []          # 0 = small (single-shape), 1 = near-kernel (window)
BW = []            # |b| per class
REPS = []          # packed 150-bit rep chain
SSH = []           # small: shift cell-id (x*15+y); windows: 0
SF0 = []           # small: packed 75-bit seam-good preimage; windows: 0
WF0 = []           # window: packed 75-bit min-poke preimage; smalls: 0
WMASK = []         # window: packed 150-bit window mask; smalls: 0
WZ = []            # window: packed 150-bit dim-1 cycle; smalls: 0
WZPRE = []         # window: packed 75-bit preimage of the cycle; smalls: 0
extras = []        # flat: (class_idx, cell, packed z, packed pre)

n_small = n_win = 0
for ci, r in enumerate(recs):
    b = np.zeros(2 * nb, dtype=np.uint8)
    for blk, gx, gy in r["b_support"]:
        b[blk * nb + base_idx[(gx, gy)]] = 1
    w = int(b.sum())
    assert w == r["b_weight"] and 0 < w <= 14
    REPS.append(pack_bits(b))
    BW.append(w)
    t = (16 - w) // 2

    # small stratum: find a seam-good translate + preimage
    found = None
    for tx in range(ELL):
        for ty in range(M):
            tb = translate_chain(b, tx, ty) if (tx, ty) != (0, 0) else b
            f = solve_f2(D2b, tb)
            assert f is not None
            f0 = seam_good_coset(f)
            if f0 is not None:
                found = ((tx, ty), tb, f0)
                break
        if found:
            break
    if r["coset_min"] <= 4:
        # small class: MUST have a seam-good translate
        assert found is not None, f"class {ci}: small but no S-translate"
        (tx, ty), tb, f0 = found
        # cert: d2 f0 = translate1 (tx,ty) rep  (Lean convention above)
        assert np.array_equal(d2b(f0), tb)
        assert np.array_equal(tb, translate_chain(b, tx, ty)) or (tx, ty) == (0, 0)
        assert pokes(f0) == 0
        KIND.append(0)
        SSH.append(tx * 15 + ty)
        SF0.append(pack_bits(f0))
        WF0.append(0); WMASK.append(0); WZ.append(0); WZPRE.append(0)
        n_small += 1
        continue

    # near-kernel class: no seam-good coset anywhere (assert), window data
    assert found is None, f"class {ci}: near-kernel but has S-translate"
    KIND.append(1)
    SSH.append(0)
    SF0.append(0)
    f = solve_f2(D2b, b)
    p, f0 = min_poke_coset(f)
    W = (b | seam_of(f0)) % 2
    Widx = np.nonzero(W)[0]
    assert 20 <= len(Widx) <= 26, (ci, len(Widx))
    # cycle table of the base window: dim must be 1
    D1W = D1b[:, Widx]
    kerW = nullspace_f2(D1W)
    assert kerW.shape[0] == 1, (ci, kerW.shape)
    z = np.zeros(2 * nb, dtype=np.uint8)
    z[Widx] = kerW[0]
    fz = solve_f2(D2b, z)
    assert fz is not None, f"class {ci}: window cycle is not a boundary"
    assert np.array_equal(d2b(fz), z)
    WF0.append(pack_bits(f0))
    WMASK.append(pack_bits(W))
    WZ.append(pack_bits(z))
    WZPRE.append(pack_bits(fz))
    n_win += 1

    # extras (t >= 2): single extra cells whose checks stay inside C_W
    if t >= 2:
        CW = set()
        for j in Widx:
            CW.update(np.nonzero(D1b[:, j])[0].tolist())
        for e in range(2 * nb):
            if W[e]:
                continue
            cols_e = set(np.nonzero(D1b[:, e])[0].tolist())
            if not cols_e <= CW:
                continue
            Wi = W.copy()
            Wi[e] = 1
            Widx_e = np.nonzero(Wi)[0]
            kerWe = nullspace_f2(D1b[:, Widx_e])
            if kerWe.shape[0] == 1:
                continue  # no new cycle through e
            assert kerWe.shape[0] == 2, (ci, e, kerWe.shape)
            # pick the basis vector supported at e
            epos = int(np.where(Widx_e == e)[0][0])
            pick = None
            for row in kerWe:
                if row[epos]:
                    pick = row
                    break
            assert pick is not None
            ze = np.zeros(2 * nb, dtype=np.uint8)
            ze[Widx_e] = pick
            fze = solve_f2(D2b, ze)
            assert fze is not None, f"class {ci}: extra cycle not a boundary"
            assert np.array_equal(d2b(fze), ze)
            extras.append((ci, e, pack_bits(ze), pack_bits(fze)))
    # t = 3: pair extensions must add no cycles (engine filters + checks;
    # asserted here so the Stage-3 table stays complete)
    if t == 3:
        CW = set()
        for j in Widx:
            CW.update(np.nonzero(D1b[:, j])[0].tolist())
        survivors = []
        off = [e for e in range(2 * nb) if not W[e]]
        for e1, e2 in itertools.combinations(off, 2):
            c1 = set(np.nonzero(D1b[:, e1])[0].tolist())
            c2 = set(np.nonzero(D1b[:, e2])[0].tolist())
            if not (c1 ^ c2) <= CW:
                continue
            survivors.append((e1, e2))
            Wi = W.copy()
            Wi[e1] = 1
            Wi[e2] = 1
            Widx_p = np.nonzero(Wi)[0]
            kerWp = nullspace_f2(D1b[:, Widx_p])
            # no cycle may use e1 or e2 beyond the tabulated singles
            for row in kerWp:
                p1 = int(np.where(Widx_p == e1)[0][0])
                p2 = int(np.where(Widx_p == e2)[0][0])
                assert not (row[p1] or row[p2]), (ci, e1, e2)
        print(f"  class {ci} (t=3): {len(survivors)} pair survivors, "
              f"0 pair cycles (asserted)")

assert n_small == 94 and n_win == 19, (n_small, n_win)
print(f"strata: {n_small} small (S-translates certified), {n_win} window; "
      f"{len(extras)} extra-cell cycles")

# global column table (each d1 column has weight exactly 3)
COLS = []
for j in range(2 * nb):
    ids = np.nonzero(D1b[:, j])[0].tolist()
    assert len(ids) == 3, (j, ids)
    COLS.extend(int(x) for x in ids)

# parity basis: every generator image has even weight
for g in range(nb):
    f = np.zeros(nb, dtype=np.uint8)
    f[g] = 1
    assert int(d2b(f).sum()) % 2 == 0

# ---------------------------------------------------------------- emit
def arr(name: str, xs, doc: str) -> str:
    body = ",".join(str(int(x)) for x in xs)
    return f"/-- {doc} -/\ndef {name} : Array Nat :=\n  #[{body}]\n\n"


def wrap_lines(s: str, width: int = 98) -> str:
    out = []
    for line in s.split("\n"):
        while len(line) > width:
            cut = line.rfind(",", 0, width)
            if cut <= 0:
                break
            out.append(line[:cut + 1])
            line = "    " + line[cut + 1:]
        out.append(line)
    return "\n".join(out)


classdata = BANNER + """
import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.Defs

/-!
# Z5Z15F2A6 dangerous-dispatch data

The 113 light-boundary translation classes of the `[[150,8,8]]` base and
their rung data (A17 §6; small classes carry a seam-good translate for
the single-shape rung, near-kernel classes carry window datasets for the
generalized window rung).  Chains are packed `Nat` bitmasks under
`cellIdx (g, blk) = blk*75 + x*15 + y`; 2-chains use `x*15 + y`.
-/

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Z5Z15F2A6

/-- Cell index of a 1-chain cell. -/
def cellIdx (p : G150 × Fin 2) : Nat :=
  p.2.val * 75 + p.1.1.val * 15 + p.1.2.val

/-- Cell index of a 2-chain cell. -/
def cell2Idx (g : G150) : Nat := g.1.val * 15 + g.2.val

"""
classdata += arr("KIND", KIND, "Rung kind per class: 0 = single-shape, 1 = window.")
classdata += arr("BW", BW, "Class boundary weight `|b|`.")
classdata += arr("REPS", REPS, "Canonical class representatives (150-bit chain masks).")
classdata += arr("SSH", SSH,
                 "Small classes: seam-good translate `c` as `x*15 + y` (0 for windows).")
classdata += arr("SF0", SF0,
                 "Small classes: seam-good preimage of `translate1 c rep` (75-bit; 0 for windows).")
classdata += arr("WF0", WF0, "Window classes: min-poke preimage of the rep (75-bit; 0 for smalls).")
classdata += arr("WMASK", WMASK, "Window classes: window mask `supp b ∪ seam` (150-bit).")
classdata += arr("WZ", WZ, "Window classes: the dim-1 window cycle (150-bit).")
classdata += arr("WZPRE", WZPRE, "Window classes: `∂₂`-preimage of the window cycle (75-bit).")
classdata += arr("XCLS", [e[0] for e in extras], "Extra-cell cycles: owning class index.")
classdata += arr("XCELL", [e[1] for e in extras], "Extra-cell cycles: the extra 1-chain cell id.")
classdata += arr("XZ", [e[2] for e in extras], "Extra-cell cycles: the cycle (150-bit).")
classdata += arr("XZPRE", [e[3] for e in extras], "Extra-cell cycles: `∂₂`-preimage (75-bit).")
classdata += arr("COLS", COLS,
                 "Flat `∂₁`-column table: checks `COLS[3j..3j+2]` of 1-chain cell `j` (each column weight 3).")
classdata += """
/-- Class representative as a chain. -/
def repChain (i : Fin 113) : G150 × Fin 2 → ZMod 2 := fun p =>
  if (REPS.getD i.val 0).testBit (cellIdx p) then 1 else 0

/-- Small-class seam-good preimage as a 2-chain. -/
def sF0Chain (i : Fin 113) : G150 → ZMod 2 := fun g =>
  if (SF0.getD i.val 0).testBit (cell2Idx g) then 1 else 0

/-- Small-class translate element. -/
def sShiftEl (i : Fin 113) : G150 :=
  (((SSH.getD i.val 0) / 15 : ℕ), ((SSH.getD i.val 0) % 15 : ℕ))

/-- Window-class min-poke preimage as a 2-chain. -/
def winF0Chain (i : Fin 113) : G150 → ZMod 2 := fun g =>
  if (WF0.getD i.val 0).testBit (cell2Idx g) then 1 else 0

/-- Window membership. -/
def winMem (i : Fin 113) (p : G150 × Fin 2) : Bool :=
  (WMASK.getD i.val 0).testBit (cellIdx p)

/-- The window cycle as a chain. -/
def winZChain (i : Fin 113) : G150 × Fin 2 → ZMod 2 := fun p =>
  if (WZ.getD i.val 0).testBit (cellIdx p) then 1 else 0

/-- Preimage of the window cycle. -/
def winZPreChain (i : Fin 113) : G150 → ZMod 2 := fun g =>
  if (WZPRE.getD i.val 0).testBit (cell2Idx g) then 1 else 0

/-- Number of extra-cell cycle entries. -/
def nExtras : Nat := XCLS.size

/-- Owning class of an extra-cell cycle (total via mod; the cert sweep
pins the semantic content). -/
def xClsIdx (e : Fin nExtras) : Fin 113 :=
  ⟨XCLS.getD e.val 0 % 113, Nat.mod_lt _ (by norm_num)⟩

/-- Extra-cell cycle as a chain. -/
def xZChain (e : Fin nExtras) : G150 × Fin 2 → ZMod 2 := fun p =>
  if (XZ.getD e.val 0).testBit (cellIdx p) then 1 else 0

/-- Preimage of an extra-cell cycle. -/
def xZPreChain (e : Fin nExtras) : G150 → ZMod 2 := fun g =>
  if (XZPRE.getD e.val 0).testBit (cell2Idx g) then 1 else 0

/-- The `∂₁` column of a 1-chain cell, from the flat table. -/
def colFn (j : G150 × Fin 2) : G150 → ZMod 2 := fun g =>
  if cell2Idx g = COLS.getD (3 * cellIdx j) 99 ∨
     cell2Idx g = COLS.getD (3 * cellIdx j + 1) 99 ∨
     cell2Idx g = COLS.getD (3 * cellIdx j + 2) 99 then 1 else 0

/-- Class `t`-value: `(16 − |b|) / 2`. -/
def tOf (i : Fin 113) : Nat := (16 - BW.getD i.val 0) / 2

end Z5Z15F2A6
end BB
end Homological
end Stabilizer
end Quantum
"""

certsweep = BANNER + """
import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.ClassData

/-!
# Z5Z15F2A6 dispatch certificates

Batched kernel certificates tying every `ClassData` table entry to the
semantic objects of the instance bundle.  Each theorem is one
`native_decide` over the finite index range; the per-class content is
documented in `ClassData.lean`.
-/

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Z5Z15F2A6

open scoped BigOperators

set_option maxRecDepth 4096

/-- Every class rep has the tabulated weight, and it is light and even. -/
theorem rep_weight_certs : ∀ i : Fin 113,
    (Finset.univ.filter fun j : G150 × Fin 2 => repChain i j ≠ 0).card
      = BW.getD i.val 0 := by
  native_decide

/-- Arithmetic shape of the class weights: `1 ≤ t` and `|b| + 2t = 16`. -/
theorem class_t_certs : ∀ i : Fin 113,
    1 ≤ tOf i ∧ BW.getD i.val 0 + 2 * tOf i = 16 := by
  native_decide

/-- Small classes: the tabulated preimage hits the tabulated translate of
the rep. -/
theorem s_translate_certs : ∀ i : Fin 113, KIND.getD i.val 0 = 0 →
    bbBoundary2Fn a150 b150 (sF0Chain i)
      = translate1 (sShiftEl i) (repChain i) := by
  native_decide

/-- Small classes: the tabulated preimage is seam-good. -/
theorem s_seam_certs : ∀ i : Fin 113, KIND.getD i.val 0 = 0 →
    ∀ j : G150 × Fin 2,
      coverData.sheet0 (coverData.liftStab (sF0Chain i)) j ≠ 0 →
      bbBoundary2Fn a150 b150 (sF0Chain i) j ≠ 0 := by
  native_decide

/-- Window classes: the tabulated preimage hits the rep itself. -/
theorem win_f0_certs : ∀ i : Fin 113, KIND.getD i.val 0 = 1 →
    bbBoundary2Fn a150 b150 (winF0Chain i) = repChain i := by
  native_decide

/-- Window classes: the mask is exactly `supp b ∪ seam`. -/
theorem win_mem_certs : ∀ i : Fin 113, KIND.getD i.val 0 = 1 →
    ∀ j : G150 × Fin 2,
      winMem i j = true ↔
        (bbBoundary2Fn a150 b150 (winF0Chain i) j ≠ 0 ∨
         coverData.sheet0 (coverData.liftStab (winF0Chain i)) j ≠ 0) := by
  native_decide

/-- Window classes: the tabulated cycle is a boundary supported in the
window. -/
theorem win_z_certs : ∀ i : Fin 113, KIND.getD i.val 0 = 1 →
    bbBoundary2Fn a150 b150 (winZPreChain i) = winZChain i ∧
    ∀ j : G150 × Fin 2, winZChain i j ≠ 0 → winMem i j = true := by
  native_decide

/-- Extra-cell cycles: each is a boundary supported in its window plus its
extra cell. -/
theorem x_z_certs : ∀ e : Fin nExtras,
    bbBoundary2Fn a150 b150 (xZPreChain e) = xZChain e ∧
    ∀ j : G150 × Fin 2, xZChain e j ≠ 0 →
      (winMem (xClsIdx e) j = true ∨ cellIdx j = XCELL.getD e.val 0) := by
  native_decide

/-- The flat `∂₁`-column table matches the boundary map on the `δ`-basis. -/
theorem col_certs : ∀ j : G150 × Fin 2,
    bbBoundary1Fn a150 b150 (Pi.single j 1) = colFn j := by
  native_decide

/-- Parity basis: every generator boundary has even total (`ZMod 2` sum
zero). -/
theorem parity_basis_certs : ∀ g : G150,
    (∑ j : G150 × Fin 2, bbBoundary2Fn a150 b150 (Pi.single g 1) j) = 0 := by
  native_decide

end Z5Z15F2A6
end BB
end Homological
end Stabilizer
end Quantum
"""

def write(path: Path, content: str, force: bool):
    if path.exists() and not force:
        print(f"REFUSING to overwrite {path} (use --force)")
        sys.exit(1)
    path.write_text(wrap_lines(content))
    print(f"wrote {path}")




# ---------------------------------------------------------------- leaves
LEAF_HEADER = BANNER + """
import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.WindowEngine

/-!
# Z5Z15F2A6 window sweep leaves ({name})

The per-class sweep certificates in the exact hypothesis shapes of
`window_sound_t1/t2/t3` ({content}).  One `native_decide` per theorem;
files are parallel build leaves.
-/

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Z5Z15F2A6

set_option maxRecDepth 4096

"""

LEAF_FOOTER = """
end Z5Z15F2A6
end BB
end Homological
end Stabilizer
end Quantum
"""


def sweep0_thm(i: int) -> str:
    K = f"({i} : Fin 113)"
    return f"""theorem win_sweep0_{i} :
    ∀ lam : Fin (2 ^ (winCellList {K}).length),
      syndFold (winCellList {K}) lam.val = 0 →
      (tableEntries {K} 150).any
        (fun pr => lam.val == localMaskOf (winCellList {K}) pr.2)
        = true := by
  native_decide

"""


def sweepE_thm(i: int) -> str:
    K = f"({i} : Fin 113)"
    return f"""theorem win_sweepE_{i} :
    ∀ e : G150 × Fin 2, winMem {K} e = false →
      survivorB {K} (cellIdx e) = true →
      ∀ lam : Fin (2 ^ (winCellList {K} ++ [cellIdx e]).length),
        syndFold (winCellList {K} ++ [cellIdx e]) lam.val = 0 →
        (tableEntries {K} (cellIdx e)).any (fun pr =>
          lam.val == localMaskOf (winCellList {K} ++ [cellIdx e]) pr.2)
          = true := by
  native_decide

"""


def sweepP_thm(i: int) -> str:
    K = f"({i} : Fin 113)"
    return f"""theorem win_sweepP_{i} :
    ∀ e₁ e₂ : G150 × Fin 2, winMem {K} e₁ = false →
      winMem {K} e₂ = false → e₁ ≠ e₂ →
      pairSurvivorB {K} (cellIdx e₁) (cellIdx e₂) = true →
      ∀ lam : Fin (2 ^ ((winCellList {K} ++ [cellIdx e₁])
          ++ [cellIdx e₂]).length),
        syndFold ((winCellList {K} ++ [cellIdx e₁]) ++ [cellIdx e₂])
          lam.val = 0 →
        ((tableEntries {K} (cellIdx e₁))
            ++ tableEntries {K} (cellIdx e₂)).any
          (fun pr => lam.val == localMaskOf
            ((winCellList {K} ++ [cellIdx e₁]) ++ [cellIdx e₂]) pr.2)
          = true := by
  native_decide

"""


def emit_leaves(force: bool):
    win_classes = [(ci, BW[ci], (16 - BW[ci]) // 2, bin(WMASK[ci]).count("1"))
                   for ci in range(N_CLASSES) if KIND[ci] == 1]
    t1 = sorted([w for w in win_classes if w[2] == 1],
                key=lambda w: -(2 ** w[3]))
    t2plus = [w for w in win_classes if w[2] >= 2]
    # greedy balance the t1 classes into two leaves by 2^|W|
    a, b, ca, cb = [], [], 0, 0
    for w in t1:
        if ca <= cb:
            a.append(w)
            ca += 2 ** w[3]
        else:
            b.append(w)
            cb += 2 ** w[3]
    print(f"leaf split: A={len(a)} classes (2^sum~{ca:.3g}), "
          f"B={len(b)} (~{cb:.3g}), C={len(t2plus)} t>=2 classes")
    fa = LEAF_HEADER.format(name="1 of 3",
                            content=f"base sweeps, t = 1 classes "
                                    f"{sorted(w[0] for w in a)}")
    for w in sorted(a):
        fa += sweep0_thm(w[0])
    fb = LEAF_HEADER.format(name="2 of 3",
                            content=f"base sweeps, t = 1 classes "
                                    f"{sorted(w[0] for w in b)}")
    for w in sorted(b):
        fb += sweep0_thm(w[0])
    fc = LEAF_HEADER.format(name="3 of 3",
                            content="the t >= 2 classes: base + "
                                    "single-extra sweeps, and the t = 3 "
                                    "pair sweep")
    for w in sorted(t2plus):
        fc += sweep0_thm(w[0])
        fc += sweepE_thm(w[0])
        if w[2] == 3:
            fc += sweepP_thm(w[0])
    write(OUT_DIR / "SweepWin1.lean", fa + LEAF_FOOTER, force)
    write(OUT_DIR / "SweepWin2.lean", fb + LEAF_FOOTER, force)
    write(OUT_DIR / "SweepWin3.lean", fc + LEAF_FOOTER, force)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--leaves-only", action="store_true")
    args = ap.parse_args()
    assert OUT_DIR.is_dir(), f"QECLean instance dir not found: {OUT_DIR}"
    if not args.leaves_only:
        write(OUT_DIR / "ClassData.lean", classdata, args.force)
        write(OUT_DIR / "CertSweep.lean", certsweep, args.force)
    emit_leaves(args.force)


if __name__ == "__main__":
    main()
