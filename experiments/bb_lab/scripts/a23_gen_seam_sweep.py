"""A23 generator: emit Z5Z15F2A6/SeamSweepData.lean (seam-coset floor data).

Emits, into the QECLean checkout (env QECLEAN_ROOT):

  * `e0mask` — the weight-40 block idempotent e0 (75-bit, cell0Idx packing);
  * the 15-entry seamC dictionary: for each nonzero kernel pattern
    (e1..e4), a shift `ge` and 2-chain `fe` with
    seamC(kerElt e) = (translate ge e0 | 0) + d2 fe   (native_decide'd
    in SeamReduction.lean);
  * the 64 w-space relations (120-bit masks, widx packing) + packed
    constants;
  * the 16,384-mask sweep certificate table (row-combination
    inconsistency masks / RREF pivot certificates + particulars +
    delta-normalized kernel extras), padded to a 32,768-entry
    mask-indexed table.

Every table is numpy-hard-asserted at emission, including an independent
re-verification pass that mirrors the Lean checker semantics
(`checkCertAt`) rather than the construction path.

Usage:  cd experiments/bb_lab && uv run python scripts/a23_gen_seam_sweep.py [--force]
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

QECLEAN_ROOT = Path(os.environ.get(
    "QECLEAN_ROOT", str(LAB_ROOT.parent.parent / "QECLean"))).resolve()
OUT_PATH = (QECLEAN_ROOT / "QEC" / "Stabilizer" / "Codes" /
            "BivariateBicycle" / "Z5Z15F2A6" / "SeamSweepData.lean")
BFD_PATH = OUT_PATH.parent / "BaseFloorData.lean"

from a23_seam_calibration import (  # noqa: E402
    A_SUPP, B_SUPP, conv_matrix, seam_maps,
)
from a23_wspace_sweep import (  # noqa: E402
    E0, MSET, SITES, fpt, nibble_type, sweep, xbar, CTAB,
)
from bb_lab.linalg import nullspace_f2, rref_f2  # noqa: E402

N = 75
MA = conv_matrix(A_SUPP)
MB = conv_matrix(B_SUPP)
D2 = np.vstack([MA, MB])


def pack(v: np.ndarray) -> int:
    # NB: `1 << i` with numpy-int64 `i` overflows at 64 — force Python ints
    return sum(1 << int(i) for i in np.flatnonzero(v))


# ---------------------------------------------------------------- z-basis
def z_basis() -> np.ndarray:
    """The A21 kernel basis, read from BaseFloorData.lean and re-verified:
    each `zᵢ` must lie in `ker ∂₂` and the four must be δ-normalized on
    the free cells `fc1..fc4 = (4,11)..(4,14)`."""
    txt = BFD_PATH.read_text()
    B = np.zeros((4, N), dtype=np.uint8)
    for i in range(4):
        m = re.search(rf"def zmask{i + 1} : ℕ := (\d+)", txt)
        assert m, f"zmask{i + 1} not found in BaseFloorData.lean"
        v = int(m.group(1))
        for j in range(N):
            B[i, j] = (v >> j) & 1
    fcs = [4 * 15 + 11 + i for i in range(4)]
    assert not ((D2 @ B.T) % 2).any(), "z-basis not in ker d2"
    for i in range(4):
        for j in range(4):
            assert B[i, fcs[j]] == (1 if i == j else 0), "not delta-normalized"
    return B


# ------------------------------------------------------------- dictionary
def build_dictionary(Z: np.ndarray) -> list[dict]:
    """Per nonzero pattern e in F2^4: (ge, fe) with
    seamC(kerElt e) = (translate ge e0 | 0) + d2 fe."""
    # base solutions f_i for the 4 kernel basis elements
    fbase = []
    for i in range(4):
        _, sC = seam_maps(Z[i])
        sR = sC[75:]
        aug = np.hstack([MB, sR[:, None]])
        R, piv = rref_f2(aug)
        assert N not in piv
        f = np.zeros(N, dtype=np.uint8)
        for r, p in enumerate(piv):
            f[p] = R[r, N]
        assert np.array_equal((MB @ f) % 2, sR)
        fbase.append(f)
    out = []
    for patt in range(1, 16):
        e = [(patt >> i) & 1 for i in range(4)]
        zeta = np.zeros(N, dtype=np.uint8)
        fe = np.zeros(N, dtype=np.uint8)
        for i in range(4):
            if e[i]:
                zeta ^= Z[i]
                fe ^= fbase[i]
        _, sC = seam_maps(zeta)
        sL, sR = sC[:75], sC[75:]
        assert np.array_equal((MB @ fe) % 2, sR)
        b = (sL + MA @ fe) % 2
        assert not ((D2 @ b) % 2).any() and b.any()
        # find ge with b(h) = e0(h + ge)  [translate ge e0]
        ge = None
        for ga in range(5):
            for gb in range(15):
                ok = all(
                    b[x * 15 + y]
                    == E0[((x + ga) % 5) * 15 + ((y + gb) % 15)]
                    for x in range(5) for y in range(15)
                )
                if ok:
                    ge = (ga, gb)
                    break
            if ge:
                break
        assert ge is not None, f"no shift for pattern {patt}"
        # full identity re-assert, exactly the Lean native_decide shape
        for x in range(5):
            for y in range(15):
                h = x * 15 + y
                lhs0 = sC[h]
                rhs0 = (E0[((x + ge[0]) % 5) * 15 + ((y + ge[1]) % 15)]
                        + (MA @ fe)[h]) % 2
                assert lhs0 == rhs0, (patt, x, y)
                assert sC[75 + h] == (MB @ fe)[h] % 2, (patt, x, y)
        out.append({"patt": patt, "ge": ge, "femask": pack(fe)})
    return out


# ------------------------------------------------- independent cert verify
def verify_certs(out: dict) -> None:
    """Mirror the Lean checkCertAt semantics (not the construction)."""
    lam = out["lam"]
    consts = out["consts"]
    lmask = [pack(lam[j]) for j in range(64)]

    def fold_xor(sel: int) -> int:
        acc = 0
        for j in range(64):
            if (sel >> j) & 1:
                acc ^= lmask[j]
        return acc

    def on_site(m: int, widx: int) -> bool:
        return bool((m >> (widx // 8)) & 1)

    def cost_of_mask(x: int) -> int:
        tot = 0
        for sk in range(15):
            nu = (x >> (8 * sk)) & 15
            nv = (x >> (8 * sk + 4)) & 15
            tot += CTAB[nibble_type(nu)][nibble_type(nv)]
        return tot

    n_incons = n_cons = 0
    for m, cert in out["certs"].items():
        if cert["kind"] == 0:
            sel = pack(cert["rowsel"])
            comb = fold_xor(sel)
            assert all(not on_site(m, i) for i in range(120)
                       if (comb >> i) & 1), m
            par = sum((consts[j]) for j in range(64) if (sel >> j) & 1) % 2
            assert par == 1, m
            n_incons += 1
            continue
        n_cons += 1
        part = pack(cert["part"])
        k = len(cert["extras"])
        assert k <= 4
        ex = [pack(e) for e in cert["extras"]] + [0] * (4 - k)
        fr = list(cert["frees"]) + [0] * (4 - k)
        # supports on-S
        assert all(on_site(m, i) for i in range(120) if (part >> i) & 1)
        for c in range(k):
            assert all(on_site(m, i) for i in range(120) if (ex[c] >> i) & 1)
            assert fr[c] < 120 and on_site(m, fr[c])
        # relations
        wpart = np.array([(part >> i) & 1 for i in range(120)], dtype=np.uint8)
        assert np.array_equal((lam @ wpart) % 2, consts)
        for c in range(k):
            we = np.array([(ex[c] >> i) & 1 for i in range(120)],
                          dtype=np.uint8)
            assert not ((lam @ we) % 2).any()
        # delta-normalization
        for c in range(k):
            for c2 in range(k):
                assert ((ex[c] >> fr[c2]) & 1) == (1 if c == c2 else 0)
        # pivots
        pividx = [i for (i, _) in cert["piv"]]
        for (idx, combo) in cert["piv"]:
            sel = pack(combo)
            comb = fold_xor(sel)
            assert (comb >> idx) & 1
            for i in range(120):
                if (comb >> i) & 1 and on_site(m, i):
                    assert i == idx or i in fr[:k], (m, idx, i)
            for c in range(k):
                assert ((ex[c] >> idx) & 1) == ((comb >> fr[c]) & 1)
            assert all(fr[c] != idx for c in range(k))
        # coverage
        for i in range(120):
            if on_site(m, i):
                assert i in pividx or i in fr[:k], (m, i)
        # reps: all 16 tau
        for t in range(16):
            rep = part
            for c in range(4):
                if (t >> c) & 1:
                    rep ^= ex[c]
            assert cost_of_mask(rep) >= 16, (m, t)
    print(f"independent verify: {n_incons} incons + {n_cons} cons certs OK")


# ------------------------------------------------------------------ emit
def fmt_nat_array(name: str, vals: list[int], chunk: int = 2048) -> str:
    """Emit an Array ℕ, chunked to keep single literals manageable."""
    if len(vals) <= chunk:
        body = ",".join(str(v) for v in vals)
        lines = wrap(body)
        return f"def {name} : Array ℕ :=\n  #[{lines}]\n"
    parts = []
    names = []
    for ci in range(0, len(vals), chunk):
        cn = f"{name}c{ci // chunk}"
        names.append(cn)
        body = ",".join(str(v) for v in vals[ci:ci + chunk])
        parts.append(f"def {cn} : Array ℕ :=\n  #[{wrap(body)}]\n")
    lines = [" ++ ".join(names[li:li + 5]) for li in range(0, len(names), 5)]
    app = ("\n    ++ ").join(lines)
    parts.append(f"def {name} : Array ℕ :=\n  {app}\n")
    return "\n".join(parts)


def wrap(s: str, width: int = 88, indent: str = "    ") -> str:
    outl = []
    cur = ""
    for tok in s.split(","):
        add = tok if not cur else "," + tok
        if len(cur) + len(add) > width:
            outl.append(cur + ",")
            cur = tok
        else:
            cur += add
    outl.append(cur)
    return ("\n" + indent).join(outl)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not QECLEAN_ROOT.is_dir():
        sys.exit(f"QECLean checkout not found at {QECLEAN_ROOT} "
                 "(set QECLEAN_ROOT)")
    if OUT_PATH.exists() and not args.force:
        sys.exit(f"REFUSING to overwrite {OUT_PATH} (use --force)")

    t0 = time.time()
    Z = z_basis()
    print("z-basis matches BaseFloorData: OK")
    dic = build_dictionary(Z)
    print(f"dictionary: {len(dic)} patterns OK")
    out = sweep()
    st = out["stats"]
    assert st["violations"] == 0 and st["min_cost"] == 16
    assert st["cons"] == 300 and st["incons"] == 16084
    verify_certs(out)

    lam = out["lam"]
    consts = out["consts"]
    relm = [pack(lam[j]) for j in range(64)]
    relconsts = int(sum(int(consts[j]) << j for j in range(64)))

    rowsel = [0] * 32768
    cidx = [0] * 32768
    cons_ms = sorted(m for m, c in out["certs"].items() if c["kind"] == 1)
    for m, cert in out["certs"].items():
        if cert["kind"] == 0:
            rowsel[m] = pack(cert["rowsel"])
            assert rowsel[m] != 0
    cpart, ck = [], []
    cex = [[], [], [], []]
    cfr = [[], [], [], []]
    cpoff, cpidx, cpsel = [0], [], []
    for ci, m in enumerate(cons_ms):
        cert = out["certs"][m]
        cidx[m] = ci + 1
        cpart.append(pack(cert["part"]))
        k = len(cert["extras"])
        ck.append(k)
        for c in range(4):
            cex[c].append(pack(cert["extras"][c]) if c < k else 0)
            cfr[c].append(cert["frees"][c] if c < k else 0)
        for (idx, combo) in cert["piv"]:
            cpidx.append(idx)
            cpsel.append(pack(combo))
        cpoff.append(len(cpidx))
    print(f"tables: {len(cons_ms)} cons, {len(cpidx)} pivot certs, "
          f"e0mask weight {int(E0.sum())}")

    parts: list[str] = []
    parts.append(f"""/-
GENERATED FILE — DO NOT HAND-EDIT.
Generator: qec-lab:experiments/bb_lab/scripts/a23_gen_seam_sweep.py
Data: recomputed in-script (repo-convention conv matrices, the A22-style
(ε,δ) fibering in w-space, the 16,384-mask site sweep); every table is
numpy-hard-asserted at emission, including an independent verification
pass mirroring `checkCertAt`.
Regen: cd experiments/bb_lab && uv run python scripts/a23_gen_seam_sweep.py --force
-/

import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.BaseFloorData

/-!
# Z5Z15F2A6 seam-coset floor data (A23)

Data for the analytic discharge of `coverData.SeamCosetFloor 16`:

* `e0mask` — the weight-40 block idempotent `e₀` of the `F₁₆` kernel
  block (`ker(A⋆) = ker(B⋆) = ker ∂₂` = its 15 translates + 0);
* `dictGx/dictGy/dictFe` — the seam dictionary: for each nonzero kernel
  pattern `e`, `seamC (kerElt e) = (translate ge e₀ | 0) + ∂₂ fe`;
* `RELM`/`RELCONSTS` — the 64 affine relations cutting out the
  realizable δ-data space `{{w(f)}}` inside `F₂^120` (`widx` packing
  `8·(3·x + (y mod 3)) + r`, `r < 4` the `A`-side δ-coords at site
  `(x, y mod 3)`, `r ≥ 4` the `B`-side at site `(x+1, y mod 3)`);
* `ROWSEL`/`CIDX` + `C*` — the per-sitemask certificates of the
  ≤7-active-site sweep: row-combination inconsistency masks (16,084
  masks) or RREF pivot certificates + particular + δ-normalized kernel
  extras (300 masks, all at exactly 7 sites), whose 2^k τ-repricings
  all have site-cost ≥ 16.

Consumed by `SeamSweep.lean` (checker + soundness) and
`SeamReduction.lean` (the `SeamCosetFloor 16` assembly).
-/

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Z5Z15F2A6

-- large array literals exceed the default elaborator recursion depth
set_option maxRecDepth 16384

/-- The block idempotent `e₀` (75-bit `cell0Idx` packing; weight 40). -/
def e0mask : ℕ := {pack(E0)}

/-- `e₀` as a chain. -/
def e0f : G150 → ZMod 2 := maskFun e0mask

""")
    gx = [0] + [d["ge"][0] for d in dic]
    gy = [0] + [d["ge"][1] for d in dic]
    fe = [0] + [d["femask"] for d in dic]
    parts.append("/-- Dictionary shifts, x-part (index = kernel pattern "
                 "`e1+2e2+4e3+8e4`). -/\n" + fmt_nat_array("dictGx", gx))
    parts.append("\n/-- Dictionary shifts, y-part. -/\n"
                 + fmt_nat_array("dictGy", gy))
    parts.append("\n/-- Dictionary 2-chains `fe` (75-bit masks). -/\n"
                 + fmt_nat_array("dictFe", fe))
    parts.append("\n/-- The 64 relation masks (120-bit, `widx` packing). -/\n"
                 + fmt_nat_array("RELM", relm))
    parts.append(f"""
/-- Relation constants (bit `j` = the affine constant of relation `j`). -/
def RELCONSTS : ℕ := {relconsts}
""")
    parts.append("\n/-- Inconsistency row-combination masks (64-bit),\n"
                 "indexed by sitemask; `0` = not an inconsistency cert. -/\n"
                 + fmt_nat_array("ROWSEL", rowsel))
    parts.append("\n/-- 1-based consistent-certificate index by sitemask;\n"
                 "`0` = none. -/\n" + fmt_nat_array("CIDX", cidx))
    parts.append("\n/-- Consistent certs: particular solutions "
                 "(120-bit `w` masks). -/\n" + fmt_nat_array("CPART", cpart))
    parts.append("\n/-- Consistent certs: kernel dimension `k ≤ 4`. -/\n"
                 + fmt_nat_array("CK", ck))
    for c in range(4):
        parts.append(f"\n/-- Kernel extra {c} (120-bit; `0` beyond `k`). "
                     f"-/\n" + fmt_nat_array(f"CEX{c}", cex[c]))
    for c in range(4):
        parts.append(f"\n/-- Free coordinate {c} (raw `widx`; `0` beyond "
                     f"`k`). -/\n" + fmt_nat_array(f"CFR{c}", cfr[c]))
    parts.append("\n/-- Pivot-certificate offsets (prefix sums per "
                 "consistent cert). -/\n" + fmt_nat_array("CPOFF", cpoff))
    parts.append("\n/-- Pivot coordinates (raw `widx`, flattened). -/\n"
                 + fmt_nat_array("CPIDX", cpidx))
    parts.append("\n/-- Pivot row-combination masks (64-bit, flattened). "
                 "-/\n" + fmt_nat_array("CPSEL", cpsel))
    parts.append("""
end Z5Z15F2A6
end BB
end Homological
end Stabilizer
end Quantum
""")
    OUT_PATH.write_text("".join(parts))
    print(f"wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.0f} KB) "
          f"in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
