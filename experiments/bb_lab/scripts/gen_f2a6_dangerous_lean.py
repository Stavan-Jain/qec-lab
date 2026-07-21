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
SWEEP_HEADER = BANNER + """
import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.KernelCert

/-!
# Z5Z15F2A6 window sweeps via pivot certificates

Every sweep obligation of `window_sound_t1/t2/t3` — base sweeps of the
19 window classes, extension sweeps of the surviving extra cells, and
the t = 3 pair sweeps — discharged by Gaussian-elimination pivot
certificates instead of `2^L` mask enumeration.

Per system, one `native_decide` checks (`cert1B`/`cert2B`) that the
tabulated pivot order triangularizes the window's syndrome map and that
the tabulated kernel generators are δ-normalized on the free positions;
`kernel_classify_dim1/dim2` then classifies every zero-syndrome mask
into the generator span, and membership of each span element in the
candidate table is a second tabulated check.  The public theorem
statements are byte-identical to the enumeration versions.

`survivorB`/`pairSurvivorB` pass for only finitely many extra cells per
t ≥ 2 class; one cheap `native_decide` certifies each survivor list, and
the public `win_sweepE_*`/`win_sweepP_*` theorems dispatch through
`mem_of_filter_eq` to the per-survivor flats.  No sweep quantifier ever
puts a ball behind a `Decidable` arrow instance.
-/

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Z5Z15F2A6

"""

LEAF_FOOTER = """
end Z5Z15F2A6
end BB
end Homological
end Stabilizer
end Quantum
"""


def _win_cells(ci: int) -> list[int]:
    return [j for j in range(150) if (WMASK[ci] >> j) & 1]


def _col_mask(c: int) -> int:
    m = 0
    for r in COLS[3 * c: 3 * c + 3]:
        m |= 1 << r
    return m


def _cw_mask(ci: int) -> int:
    cw = 0
    for j in _win_cells(ci):
        cw |= _col_mask(j)
    return cw


def _survivors(ci: int) -> list[int]:
    """Cells passing the Lean gate  !winMem && survivorB, ascending —
    must match `(List.range 150).filter ...` exactly."""
    cw = _cw_mask(ci)
    out = []
    for c in range(150):
        if (WMASK[ci] >> c) & 1:
            continue
        cm = _col_mask(c)
        if cm & cw == cm:
            out.append(c)
    return out


def _pair_survivors(ci: int) -> list[int]:
    """Ordered pair indices n = c1*150 + c2 passing the Lean pair gate,
    ascending — must match `(List.range 22500).filter ...` exactly."""
    cw = _cw_mask(ci)
    win = WMASK[ci]
    out = []
    for n in range(150 * 150):
        c1, c2 = n // 150, n % 150
        if c1 == c2 or (win >> c1) & 1 or (win >> c2) & 1:
            continue
        x = _col_mask(c1) ^ _col_mask(c2)
        if x & cw == x:
            out.append(n)
    return out


# ------------------------------------------------- pivot certificates
def _kernel_f2(cols: list[int]) -> list[int]:
    """Kernel basis of the syndrome map as position masks."""
    done: list[tuple[int, int]] = []
    basis = []
    for i, v0 in enumerate(cols):
        v, m = v0, 1 << i
        for pv, pm in done:
            if v & (pv & -pv):
                v ^= pv
                m ^= pm
        if v == 0:
            basis.append(m)
        else:
            done.append((v, m))
    return basis


def _pivot_cert(cells: list[int]):
    """Greedy elimination-order pivot certificate for the window system.

    Returns (dim, frees, gens, piv) with `gens` δ-normalized on `frees`
    and `piv` a list of (position, row) pairs in elimination order
    satisfying the Lean `pivB` semantics (verified by `_check_cert`)."""
    cols = [_col_mask(c) for c in cells]
    basis = _kernel_f2(cols)
    d = len(basis)
    assert d in (1, 2), (len(cells), d)
    supports = [[i for i in range(len(cells)) if (b >> i) & 1]
                for b in basis]
    for frees in itertools.product(*supports):
        if len(set(frees)) != d:
            continue
        B = basis[:]
        M = [[(B[i] >> f) & 1 for f in frees] for i in range(d)]
        if d == 1:
            if M[0][0] != 1:
                continue
            Bn = [B[0]]
        else:
            if (M[0][0] * M[1][1]) ^ (M[0][1] * M[1][0]) != 1:
                continue
            Bn = [(B[0] if M[1][1] else 0) ^ (B[1] if M[0][1] else 0),
                  (B[0] if M[1][0] else 0) ^ (B[1] if M[0][0] else 0)]
        if not all(((Bn[i] >> frees[j]) & 1) == (1 if i == j else 0)
                   for i in range(d) for j in range(d)):
            continue
        rem = set(range(len(cells))) - set(frees)
        piv = []
        ok = True
        while rem:
            found = None
            for j in sorted(rem):
                other = 0
                for s in rem:
                    if s != j:
                        other |= cols[s]
                private = cols[j] & ~other
                if private:
                    found = (j, (private & -private).bit_length() - 1)
                    break
            if not found:
                ok = False
                break
            piv.append(found)
            rem.discard(found[0])
        if ok:
            return d, list(frees), Bn, piv
    raise AssertionError(f"no pivot certificate for {len(cells)}-cell system")


def _check_cert(cells: list[int], dim: int, frees: list[int],
                gens: list[int], piv: list[tuple[int, int]]):
    """Replicate the Lean cert1B/cert2B semantics bit-for-bit."""
    L = len(cells)
    cols = [_col_mask(c) for c in cells]

    def synd(mask):
        s = 0
        for i in range(L):
            if (mask >> i) & 1:
                s ^= cols[i]
        return s

    # pivB
    for t, (j, r) in enumerate(piv):
        assert j < L
        assert (cols[j] >> r) & 1
        sel = 0
        pos = 0
        for j2, _ in piv[t + 1:]:
            sel |= cols[j2]
            pos |= 1 << j2
        assert not (sel >> r) & 1
        assert not (pos >> j) & 1
    pm = 0
    for j, _ in piv:
        pm |= 1 << j
    fm = 0
    for f in frees:
        fm |= 1 << f
    assert pm | fm == (1 << L) - 1 and pm & fm == 0
    for i, g in enumerate(gens):
        assert g < (1 << L) and synd(g) == 0
        for k, f in enumerate(frees):
            assert ((g >> f) & 1) == (1 if i == k else 0)
    assert dim == len(frees) == len(gens)


def _body0(i: int, v: str) -> str:
    K = f"({i} : Fin 113)"
    return (f"syndFold (winCellList {K}) {v} = 0 →\n"
            f"      (tableEntries {K} 150).any\n"
            f"        (fun pr => {v} == localMaskOf (winCellList {K}) pr.2)\n"
            f"        = true")


def _bodyE(i: int, v: str) -> str:
    K = f"({i} : Fin 113)"
    return (f"syndFold (winCellList {K} ++ [cellIdx e]) {v} = 0 →\n"
            f"        (tableEntries {K} (cellIdx e)).any (fun pr =>\n"
            f"          {v} == localMaskOf (winCellList {K} ++ [cellIdx e]) pr.2)\n"
            f"          = true")


def _bodyP(i: int, v: str) -> str:
    K = f"({i} : Fin 113)"
    return (f"syndFold ((winCellList {K} ++ [cellIdx e₁]) ++ [cellIdx e₂])\n"
            f"          {v} = 0 →\n"
            f"        ((tableEntries {K} (cellIdx e₁))\n"
            f"            ++ tableEntries {K} (cellIdx e₂)).any\n"
            f"          (fun pr => {v} == localMaskOf\n"
            f"            ((winCellList {K} ++ [cellIdx e₁]) ++ [cellIdx e₂]) pr.2)\n"
            f"          = true")


def _piv_lit(piv: list[tuple[int, int]]) -> str:
    return "[" + ", ".join(f"({j}, {r})" for j, r in piv) + "]"


def _cert_and_flat(name: str, cell_ids: list[int], cells: str, table: str,
                   len_expr: str) -> str:
    """One batched native_decide certificate + the flat public theorem
    (byte-identical statement to the enumeration version) for one window
    system."""
    dim, frees, gens, piv = _pivot_cert(cell_ids)
    _check_cert(cell_ids, dim, frees, gens, piv)
    L = len(cell_ids)
    if dim == 1:
        certcall = (f"cert1B {cells}\n"
                    f"        {_piv_lit(piv)}\n"
                    f"        {frees[0]} {gens[0]}")
        spans = ["0", str(gens[0])]
        destr = "obtain ⟨⟨⟨hlen, hkc⟩, htab0⟩, htab1⟩ := hcert"
        rc = "rcases kernel_classify_dim1 hkc lam.val hlt hs with h | h"
        cases = ("  · rw [h]\n    exact htab0\n"
                 "  · rw [h]\n    exact htab1\n")
    else:
        certcall = (f"cert2B {cells}\n"
                    f"        {_piv_lit(piv)}\n"
                    f"        {frees[0]} {frees[1]} {gens[0]} {gens[1]}")
        spans = ["0", str(gens[0]), str(gens[1]),
                 f"({gens[0]} ^^^ {gens[1]})"]
        destr = ("obtain ⟨⟨⟨⟨⟨hlen, hkc⟩, htab0⟩, htab1⟩, htab2⟩, htab3⟩"
                 " := hcert")
        rc = ("rcases kernel_classify_dim2 hkc lam.val hlt hs"
              " with h | h | h | h")
        cases = ("  · rw [h]\n    exact htab0\n"
                 "  · rw [h]\n    exact htab1\n"
                 "  · rw [h]\n    exact htab2\n"
                 "  · rw [h]\n    exact htab3\n")
    tabs = "\n      && ".join(
        f"(({table}).any (fun pr => {m} == localMaskOf {cells} pr.2))"
        for m in spans)
    return f"""private theorem {name}_cert :
    (({len_expr} == {L})
      && ({certcall})
      && {tabs}) = true := by
  native_decide

theorem {name} :
    ∀ lam : Fin (2 ^ {len_expr}),
      syndFold {cells} lam.val = 0 →
      ({table}).any
        (fun pr => lam.val == localMaskOf {cells} pr.2)
        = true := by
  have hcert := {name}_cert
  simp only [Bool.and_eq_true, beq_iff_eq] at hcert
  {destr}
  rw [hlen]
  intro lam hs
  have hlt : lam.val < 2 ^ {len_expr} := by
    rw [hlen]
    exact lam.isLt
  {rc}
{cases}
"""


def sweep0_thm(i: int) -> str:
    K = f"({i} : Fin 113)"
    return _cert_and_flat(f"win_sweep0_{i}", _win_cells(i),
                          f"(winCellList {K})", f"tableEntries {K} 150",
                          f"(winCellList {K}).length")


def sweepE_flat(i: int, c: int) -> str:
    K = f"({i} : Fin 113)"
    return _cert_and_flat(f"win_sweepE_{i}_c{c}", _win_cells(i) + [c],
                          f"(winCellList {K} ++ [{c}])",
                          f"tableEntries {K} {c}",
                          f"(winCellList {K} ++ [{c}]).length")


def sweepP_flat(i: int, c1: int, c2: int) -> str:
    K = f"({i} : Fin 113)"
    return _cert_and_flat(
        f"win_sweepP_{i}_p{c1}_{c2}", _win_cells(i) + [c1, c2],
        f"((winCellList {K} ++ [{c1}]) ++ [{c2}])",
        f"(tableEntries {K} {c1}) ++ tableEntries {K} {c2}",
        f"((winCellList {K} ++ [{c1}]) ++ [{c2}]).length")


def _lit(lst: list[int]) -> str:
    return "[" + ", ".join(str(x) for x in lst) + "]"


def surv_cert(i: int) -> str:
    return f"""private theorem win_surv_{i} :
    (List.range 150).filter (fun c =>
      !(winMem ({i} : Fin 113) (coordOfC1 c)) && survivorB ({i} : Fin 113) c)
      = {_lit(_survivors(i))} := by
  native_decide

"""


def pairs_cert(i: int) -> str:
    return f"""private theorem win_pairs_{i} :
    (List.range 22500).filter (fun n =>
      !(n / 150 == n % 150) &&
      (!(winMem ({i} : Fin 113) (coordOfC1 (n / 150))) &&
        (!(winMem ({i} : Fin 113) (coordOfC1 (n % 150))) &&
          pairSurvivorB ({i} : Fin 113) (n / 150) (n % 150))))
      = {_lit(_pair_survivors(i))} := by
  native_decide

"""


def sweepE_pub(i: int) -> str:
    K = f"({i} : Fin 113)"
    surv = _survivors(i)
    assert surv, f"class {i}: no survivors — vacuous E dispatch unexpected"
    if len(surv) == 1:
        c = surv[0]
        cases = (f"  have hc : cellIdx e = {c} := by simpa using hmem\n"
                 f"  rw [hc]\n"
                 f"  exact win_sweepE_{i}_c{c}\n")
    else:
        ors = " ∨ ".join(f"cellIdx e = {c}" for c in surv)
        pat = " | ".join(["hc"] * len(surv))
        cases = (f"  have hm' : {ors} := by simpa using hmem\n"
                 f"  rcases hm' with {pat}\n")
        for c in surv:
            cases += (f"  · rw [hc]\n"
                      f"    exact win_sweepE_{i}_c{c}\n")
    return f"""theorem win_sweepE_{i} :
    ∀ e : G150 × Fin 2, winMem {K} e = false →
      survivorB {K} (cellIdx e) = true →
      ∀ lam : Fin (2 ^ (winCellList {K} ++ [cellIdx e]).length),
        {_bodyE(i, "lam.val")} := by
  intro e he hs
  have hp : (fun c =>
      !(winMem {K} (coordOfC1 c)) && survivorB {K} c)
      (cellIdx e) = true := by
    show (!(winMem {K} (coordOfC1 (cellIdx e))) &&
      survivorB {K} (cellIdx e)) = true
    rw [coordOfC1_cellIdx, he, hs]
    rfl
  have hmem : cellIdx e ∈ {_lit(surv)} :=
    mem_of_filter_eq win_surv_{i} (List.mem_range.mpr (cellIdx_lt e)) hp
{cases}
"""


def sweepP_pub(i: int) -> str:
    K = f"({i} : Fin 113)"
    pairs = _pair_survivors(i)
    assert pairs, f"class {i}: no pair survivors — vacuous P dispatch"
    ors = " ∨ ".join(f"cellIdx e₁ * 150 + cellIdx e₂ = {n}" for n in pairs)
    pat = " | ".join(["hc"] * len(pairs))
    cases = (f"  have hm' : {ors} := by simpa using hmem\n"
             f"  rcases hm' with {pat}\n")
    for n in pairs:
        c1, c2 = n // 150, n % 150
        cases += (f"  · have hc₁ : cellIdx e₁ = {c1} := by omega\n"
                  f"    have hc₂ : cellIdx e₂ = {c2} := by omega\n"
                  f"    rw [hc₁, hc₂]\n"
                  f"    exact win_sweepP_{i}_p{c1}_{c2}\n")
    return f"""theorem win_sweepP_{i} :
    ∀ e₁ e₂ : G150 × Fin 2, winMem {K} e₁ = false →
      winMem {K} e₂ = false → e₁ ≠ e₂ →
      pairSurvivorB {K} (cellIdx e₁) (cellIdx e₂) = true →
      ∀ lam : Fin (2 ^ ((winCellList {K} ++ [cellIdx e₁])
          ++ [cellIdx e₂]).length),
        {_bodyP(i, "lam.val")} := by
  intro e₁ e₂ h₁ h₂ hne hp
  have hlt₁ := cellIdx_lt e₁
  have hlt₂ := cellIdx_lt e₂
  have hdiv : (cellIdx e₁ * 150 + cellIdx e₂) / 150 = cellIdx e₁ := by omega
  have hmod : (cellIdx e₁ * 150 + cellIdx e₂) % 150 = cellIdx e₂ := by omega
  have hcne : cellIdx e₁ ≠ cellIdx e₂ := fun h => hne (cellIdx_inj h)
  have hp' : (fun n =>
      !(n / 150 == n % 150) &&
      (!(winMem {K} (coordOfC1 (n / 150))) &&
        (!(winMem {K} (coordOfC1 (n % 150))) &&
          pairSurvivorB {K} (n / 150) (n % 150))))
      (cellIdx e₁ * 150 + cellIdx e₂) = true := by
    show (!((cellIdx e₁ * 150 + cellIdx e₂) / 150 ==
        (cellIdx e₁ * 150 + cellIdx e₂) % 150) &&
      (!(winMem {K}
          (coordOfC1 ((cellIdx e₁ * 150 + cellIdx e₂) / 150))) &&
        (!(winMem {K}
            (coordOfC1 ((cellIdx e₁ * 150 + cellIdx e₂) % 150))) &&
          pairSurvivorB {K}
            ((cellIdx e₁ * 150 + cellIdx e₂) / 150)
            ((cellIdx e₁ * 150 + cellIdx e₂) % 150)))) = true
    rw [hdiv, hmod, coordOfC1_cellIdx, coordOfC1_cellIdx, h₁, h₂, hp]
    simp [hcne]
  have hmem : cellIdx e₁ * 150 + cellIdx e₂ ∈ {_lit(pairs)} :=
    mem_of_filter_eq win_pairs_{i} (List.mem_range.mpr (by omega)) hp'
{cases}
"""


def emit_leaves(force: bool):
    win_classes = [(ci, BW[ci], (16 - BW[ci]) // 2, bin(WMASK[ci]).count("1"))
                   for ci in range(N_CLASSES) if KIND[ci] == 1]
    t1 = sorted([w for w in win_classes if w[2] == 1])
    t2plus = sorted([w for w in win_classes if w[2] >= 2])
    # sanity: extras (kernel-growing cells) must be among the survivors
    for ci, e, _, _ in extras:
        assert e in _survivors(ci), (ci, e)
    nsys = 0
    body = SWEEP_HEADER
    for w in t1:
        body += sweep0_thm(w[0])
        nsys += 1
    for w in t2plus:
        body += sweep0_thm(w[0])
        nsys += 1
        for c in _survivors(w[0]):
            body += sweepE_flat(w[0], c)
            nsys += 1
        if w[2] == 3:
            for n in _pair_survivors(w[0]):
                body += sweepP_flat(w[0], n // 150, n % 150)
                nsys += 1
    for w in t2plus:
        body += surv_cert(w[0])
        if w[2] == 3:
            body += pairs_cert(w[0])
    for w in t2plus:
        body += sweepE_pub(w[0])
        if w[2] == 3:
            body += sweepP_pub(w[0])
    print(f"pivot certificates emitted for {nsys} window systems")
    write(OUT_DIR / "SweepWin.lean", body + LEAF_FOOTER, force)


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
