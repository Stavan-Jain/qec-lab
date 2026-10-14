"""Generate the gross seam-coset Brouwer–Zimmermann certificate data.

Emits `QEC/Stabilizer/Codes/BivariateBicycle/Gross/SeamBZData.lean` in the QECLean
checkout: the data behind `GrossBZ.grossSeamBZ`, an instance of the generic
`XDoubleCoverData.SeamBZCertificate` (`Framework/Homological/BBSeamBZ.lean`), which
discharges both `MImBound` and `LightStabilizerClassification` for the gross
[[144,12,12]] code.

Conventions follow QECLean: base group Z6 x Z6, qubit ((a, b), j) -> bit 6a + b + 36j,
d2 f = (A*f | B*f) with A = x^3 + y + y^2, B = y^3 + x + x^2, and row (p, j) of d2
reads f(p - a) for a in supp A (j = 0) / supp B (j = 1). The gross cover reduces x
mod 12 -> 6 with deck x^6 and section "canonical representative"; seamC is the
sheet-1 part of the lifted stabilizer (QECLean `seamCSparse`).

What is emitted (every item is checked here first; nothing is written unless all
checks pass):
  * the ker d2 basis (supports + free cells) and a 30-step peeling order;
  * the 42 generator words of the coset code K = seamC(ker d2) + im d2 — the 36
    boundaries d2(delta_g), then the 6 seam words tagged above bit 72 by their
    ker-d2 coordinates;
  * two disjoint information sets of K ("windows") with tagged systematic rows,
    each reconstructing every generator;
  * the 252 light boundaries (weight <= 10: 36 hexagons + 216 D-pairs) with shapes.
The two walks the Lean kernel runs (every XOR of <= 5 rows of window 1 and <= 4
rows of window 2: 510 414 words) are simulated with the same predicate `okGross`.

Run from experiments/bb_lab/:
    uv run python scripts/gen_gross_bz.py            # validate only
    uv run python scripts/gen_gross_bz.py --emit     # write (refuses to clobber)
    uv run python scripts/gen_gross_bz.py --emit --force
The output lands in the QECLean checkout (env QECLEAN_ROOT; default: the QECLean
directory sibling to this qec-lab repo, resolved from the script's location).
"""

import argparse
import itertools
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
QECLEAN_ROOT = os.environ.get(
    "QECLEAN_ROOT", os.path.normpath(os.path.join(_here, "..", "..", "..", "..", "QECLean")))
TARGET = os.path.join("QEC", "Stabilizer", "Codes", "BivariateBicycle", "Gross",
                      "SeamBZData.lean")

L, M = 6, 6
N = 72
LA = [(3, 0), (0, 1), (0, 2)]
LB = [(0, 3), (1, 0), (2, 0)]
G = [(a, b) for a in range(L) for b in range(M)]
PAIRDIRS = [(0, 1), (0, 5), (3, 1), (3, 2), (3, 4), (3, 5),
            (1, 0), (1, 3), (2, 3), (4, 3), (5, 0), (5, 3)]  # QECLean `pairDirections`
T1, T2, W = 5, 4, 10  # walk depths; every K-word of weight <= W is reached


def add(g, h): return ((g[0] + h[0]) % L, (g[1] + h[1]) % M)
def sub(g, h): return ((g[0] - h[0]) % L, (g[1] - h[1]) % M)
def qbit(h, j): return 6 * h[0] + h[1] + 36 * j
def rowcells(p, j): return [sub(p, a) for a in (LA if j == 0 else LB)]


def d2(f):
    """Mask of d2 of the 2-chain with support f."""
    m = 0
    for p in G:
        for j, LS in ((0, LA), (1, LB)):
            if sum(1 for a in LS if sub(p, a) in f) % 2:
                m |= 1 << qbit(p, j)
    return m


def seamC(xi):
    """Mask of the seam-crossing map on the 2-chain with support xi."""
    m = 0
    for p in G:
        for j, LS in ((0, LA), (1, LB)):
            s = 0
            for a in LS:
                g1 = (p[0] + 6 - a[0]) % 12
                g2 = (p[1] - a[1]) % 6
                if g1 < 6 and (g1 % 6, g2) in xi:
                    s ^= 1
            if s:
                m |= 1 << qbit(p, j)
    return m


# The systematic ker d2 basis (free cells first-to-last).
KB = [
    [(0, 0), (0, 1), (0, 3), (0, 4), (1, 0), (1, 1), (1, 3), (1, 4), (3, 0), (3, 1), (3, 3),
     (3, 4), (4, 0), (4, 1), (4, 3), (4, 4)],
    [(0, 0), (0, 2), (0, 3), (0, 5), (1, 0), (1, 2), (1, 3), (1, 5), (3, 0), (3, 2), (3, 3),
     (3, 5), (4, 0), (4, 2), (4, 3), (4, 5)],
    [(0, 0), (0, 4), (1, 1), (1, 5), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1), (3, 2),
     (3, 5), (4, 0), (4, 1), (4, 2), (4, 3), (5, 0), (5, 2)],
    [(0, 0), (0, 3), (0, 4), (0, 5), (1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (2, 4),
     (2, 5), (3, 2), (3, 4), (4, 0), (4, 2), (5, 1), (5, 3)],
    [(0, 1), (0, 5), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 5), (3, 0),
     (3, 1), (3, 2), (3, 3), (4, 0), (4, 2), (5, 0), (5, 4)],
    [(0, 0), (0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 0), (2, 1), (2, 2), (2, 3), (3, 1),
     (3, 2), (3, 3), (3, 4), (4, 1), (4, 3), (5, 1), (5, 5)],
]
FREE = [(4, 4), (4, 5), (5, 2), (5, 3), (5, 4), (5, 5)]
TWO_N = 1 << N


def build():
    KBS = [set(k) for k in KB]
    for k in KBS:
        assert d2(k) == 0, "kernel basis vector not in ker d2"
    for i, k in enumerate(KBS):
        for j, F in enumerate(FREE):
            assert (F in k) == (i == j), "kernel basis not systematic on the free cells"

    # peeling order (greedy): each row read has exactly one unknown cell
    known, steps = list(FREE), []
    while len(known) < len(G):
        for p, j in ((p, j) for p in G for j in (0, 1)):
            cells = rowcells(p, j)
            unk = [c for c in cells if c not in known]
            if len(unk) == 1 and cells.count(unk[0]) == 1:
                steps.append((unk[0], (p, j)))
                known.append(unk[0])
                break
        else:
            raise AssertionError("peeling stuck")

    # cross-check the seam convention against a mask certified in Lean history
    assert seamC(KBS[0]) == 0x1b0000006db, "seamC convention drifted"

    def tag(xi): return sum(1 << j for j, F in enumerate(FREE) if F in xi)
    gens = [d2({g}) for g in G] + [seamC(KBS[i]) + TWO_N * tag(KBS[i]) for i in range(6)]

    def window(forbid):
        rows, piv = list(gens), []
        for i in range(36):
            for c in range(N):
                if c in forbid or c in piv:
                    continue
                r = next((r for r in range(i, len(rows)) if rows[r] >> c & 1), None)
                if r is None:
                    continue
                rows[i], rows[r] = rows[r], rows[i]
                for rr in range(len(rows)):
                    if rr != i and rows[rr] >> c & 1:
                        rows[rr] ^= rows[i]
                piv.append(c)
                break
            else:
                raise AssertionError("no information set")
        assert all(r == 0 for r in rows[36:]), "dim K != 36"
        return list(zip(piv, rows[:36]))

    w1 = window(set())
    w2 = window({p for p, _ in w1})
    assert not ({p for p, _ in w1} & {p for p, _ in w2})

    def recon(w, m):
        x = 0
        for p, r in w:
            if m >> p & 1:
                x ^= r
        return x
    for g in gens:
        assert recon(w1, g % TWO_N) == g and recon(w2, g % TWO_N) == g, "window fails"

    light = {}
    for g in G:
        light.setdefault(d2({g}), (g, None))
    for g in G:
        for d in PAIRDIRS:
            light.setdefault(d2({g, add(g, d)}), (g, d))
    masks = sorted(light)
    assert len(masks) == 252 and {bin(x).count("1") for x in masks} == {6, 10}

    # simulate both kernel walks with `okGross`
    lightset = set(masks)

    def ok(x):
        heavy = bin(x % TWO_N).count("1") >= 11
        return x == 0 or heavy or (x >> N == 0 and x in lightset)
    for w, t in ((w1, T1), (w2, T2)):
        rows = [r for _, r in w]
        for k in range(1, t + 1):
            for comb in itertools.combinations(rows, k):
                x = 0
                for r in comb:
                    x ^= r
                assert ok(x), "walk rejects a word"
    return steps, gens, w1, w2, masks, light


def lit_group(g): return f"({g[0]}, {g[1]})"


def wrap(items, indent="  ", width=96):
    lines, cur = [], indent
    for k, it in enumerate(items):
        piece = it + ("," if k < len(items) - 1 else "")
        if len(cur) + len(piece) + 1 > width and cur.strip():
            lines.append(cur.rstrip())
            cur = indent
        cur += piece + " "
    lines.append(cur.rstrip())
    return "\n".join(lines)


BANNER = """/-
GENERATED FILE — DO NOT HAND-EDIT (edits WILL be clobbered by regen).
Generator : qec-lab:experiments/bb_lab/scripts/gen_gross_bz.py
Data      : none — computed from the code, validated in Python before emission
Regen     : in a qec-lab checkout's experiments/bb_lab (QECLEAN_ROOT points here):
            uv run python scripts/gen_gross_bz.py --emit --force
-/"""


def emit(steps, gens, w1, w2, masks, light):
    out = [BANNER + """
/-
# Gross seam-coset BZ certificate data

Data for `GrossBZ.grossSeamBZ` (`SeamBZ.lean`), in the base-qubit encoding
`((a, b), j) ↦ 6·a + b + 36·j`; words carry their seam tag above bit 72.
-/

import QEC.Stabilizer.Codes.BivariateBicycle.Gross.Defs

namespace Quantum.Stabilizer.Homological.BB.GrossBZ
"""]
    out.append("/-- The supports of the `ker ∂₂` basis vectors, in free-cell order. -/")
    out.append("def kbSupports : List (List BaseGroup) := [")
    out.append(",\n".join("  [" + wrap([lit_group(c) for c in k], indent="   ")[3:] + "]"
                          for k in KB) + "]")
    out.append("")
    out.append("/-- The free cells of the `ker ∂₂` basis. -/")
    out.append("def kbFree : List BaseGroup := [" + ", ".join(lit_group(c) for c in FREE) + "]")
    out.append("")
    out.append("/-- The peeling order: `(cell, row)` pairs. -/")
    out.append("def kbPeel : List (BaseGroup × (BaseGroup × Fin 2)) := [")
    out.append(wrap(["(" + lit_group(h) + ", (" + lit_group(p) + f", {j}))"
                     for h, (p, j) in steps]) + "]")
    out.append("")
    out.append("/-- The generator words: the 36 boundaries `∂₂ δ_g` (`g` in row-major order),")
    out.append("then the 6 tagged seam words of the kernel basis. -/")
    out.append("def gensData : List ℕ := [")
    out.append(wrap([hex(g) for g in gens]) + "]")
    for name, w in (("window₁", w1), ("window₂", w2)):
        out.append("")
        out.append(f"/-- The pivot bits of `{name}`. -/")
        out.append(f"def {name}Pivots : List ℕ := [")
        out.append(wrap([str(pv) for pv, _ in w]) + "]")
        out.append("")
        out.append(f"/-- The tagged systematic rows of `{name}`. -/")
        out.append(f"def {name}Rows : List ℕ := [")
        out.append(wrap([hex(r) for _, r in w]) + "]")
    out.append("")
    out.append("/-- The light boundaries (weight `≤ 10`) with their shapes: a hexagon `∂₂ δ_g`")
    out.append("(`none`) or a D-pair `∂₂ (δ_g + δ_{g+d})` (`some d`). -/")
    out.append("def lightData : List (ℕ × BaseGroup × Option BaseGroup) := [")
    ent = []
    for m in masks:
        g, d = light[m]
        ent.append(f"({hex(m)}, {lit_group(g)}, "
                   + ("none" if d is None else f"some {lit_group(d)}") + ")")
    out.append(wrap(ent) + "]")
    out.append("")
    out.append("/-- The light boundary masks. -/")
    out.append("def lightMasks : List ℕ := [")
    out.append(wrap([hex(m) for m in masks]) + "]")
    out.append("")
    out.append("end Quantum.Stabilizer.Homological.BB.GrossBZ")
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--emit", action="store_true", help="write the Lean data file")
    ap.add_argument("--force", action="store_true", help="overwrite an existing target")
    ap.add_argument("--out", help="write here instead of the QECLean target")
    args = ap.parse_args()
    data = build()
    print("validated: peel 30 steps, 42 generators, 2 windows, 252 light boundaries, "
          "walks 443703 + 66711 words all accepted")
    if not args.emit:
        return
    path = args.out or os.path.join(QECLEAN_ROOT, TARGET)
    if not args.out and not os.path.isdir(QECLEAN_ROOT):
        sys.exit(f"QECLean checkout not found at {QECLEAN_ROOT} (set QECLEAN_ROOT)")
    if os.path.exists(path) and not args.force:
        sys.exit(f"refusing to overwrite {path} without --force")
    with open(path, "w") as fh:
        fh.write(emit(*data))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
