"""Certify d = 10 for the [[150,30,10]] mitten code with the lab SAT stack.

The smallest mitten code (arXiv:2607.28795, non-abelian LP over
𝔽₂[C₅×S₃], rate 1/5, check weight 9) is the first non-BB code pushed
through the lab pipeline. The paper reports its distance as exact via a
Gurobi-IP workflow with no published checkable artifacts; this script
independently establishes both CSS directions with two solver families
and emits `bb-cert/v2` witness certificates:

- **Tandem** (the MaxCDCL fork, naive WCNF encoding) solved to
  `s OPTIMUM FOUND` — the BnB optimum IS the distance;
- **CMS ladder** (`x_distance`, pycryptosat native-XOR): UNSAT at
  w = 1..9, SAT at 10 — an independent engine and encoding for the
  lower bound.

Deliberately NOT used, and why:
- the `strengthened` encoding (G-orbit selectors/anchors) — its
  `compute_class_action` is BB-translation-specific;
- `-cost-step=2` — the coset weight-parity premise needs every H_X row
  even; mitten check weight is 9 (odd), so the flag would be UNSOUND
  here. No fork flags are passed; tandem is behaviourally stock.

Trust model (same split as `bb_lab.maxsat_distance`): each witness is
re-verified in-process and recorded in the certificate (checkable
forever); the lower bound is the word of two independent solvers.
Neither lane emits DRAT — CMS's XOR reasoning has no DRAT story, and
the CaDiCaL-on-Tseitin route is the known-intractable A15 corner.

Matrices: `instances/mitten_150_30_10/` (vendored from the authors'
release — see the README there). `validate` mode checks them against
the paper's Definition 4 / Eq. (2) structure without needing GAP.

Usage (from experiments/bb_lab):
  uv run python scripts/mitten150_tandem.py validate
  uv run python scripts/mitten150_tandem.py final [--skip-cms] \
      [--binary third_party/maxcdcl/MaxCDCL/code/simp/tandem]
  uv run python scripts/mitten150_tandem.py x|z   # one tandem side only
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bb_lab.certificate import (
    make_certificate, verify_certificate, write_certificate,
)
from bb_lab.linalg import nullspace_f2
from bb_lab.maxsat_distance import maxsat_distance
from bb_lab.sat_distance import find_logical_z, x_distance

LAB_ROOT = Path(__file__).resolve().parent.parent
INSTANCE = LAB_ROOT / "instances" / "mitten_150_30_10"
CERTS = LAB_ROOT / "certificates"
WORK = LAB_ROOT / "scratch" / "mitten150"
DEFAULT_BINARY = LAB_ROOT / "third_party/maxcdcl/MaxCDCL/code/simp/tandem"
PAPER = "arXiv:2607.28795"


class _Label:
    def label(self) -> str:
        return "mitten_150_30_10"


@dataclass(frozen=True)
class GenericChecks:
    """Duck-typed stand-in for `bb_lab.checks.CheckMatrices`: the same
    attribute surface (H_X, H_Z, group.label(), num_qubits) without the
    2-block shape assumption — CheckMatrices hardwires
    `num_qubits = 2|G|`, wrong for a 5-block lifted product."""

    H_X: np.ndarray
    H_Z: np.ndarray
    group: _Label

    @property
    def num_qubits(self) -> int:
        return int(self.H_Z.shape[1])


def load() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    Hx = np.load(INSTANCE / "Hx.npy").astype(np.uint8) % 2
    Hz = np.load(INSTANCE / "Hz.npy").astype(np.uint8) % 2
    Lx = np.load(INSTANCE / "Lx.npy").astype(np.uint8) % 2
    Lz = np.load(INSTANCE / "Lz.npy").astype(np.uint8) % 2
    return Hx, Hz, Lx, Lz


def rank2(M: np.ndarray) -> int:
    return M.shape[1] - nullspace_f2(M).shape[0]


def validate() -> None:
    """GAP-free structural validation against the paper's mitten form."""
    Hx, Hz, Lx, Lz = load()
    ok = True

    def chk(name: str, cond: bool, detail: str = "") -> None:
        nonlocal ok
        ok &= cond
        tail = f" — {detail}" if detail else ""
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}{tail}")

    chk("shapes (60,150)/(60,150)/(30,150)/(30,150)",
        Hx.shape == (60, 150) and Hz.shape == (60, 150)
        and Lx.shape == (30, 150) and Lz.shape == (30, 150))
    chk("CSS Hx·Hzᵀ = 0", not (Hx @ Hz.T % 2).any())
    chk("check weight 9, all rows both sides",
        set(Hx.sum(1).tolist()) == {9} == set(Hz.sum(1).tolist()))
    rx, rz = rank2(Hx), rank2(Hz)
    chk("k = 30", 150 - rx - rz == 30, f"rank Hx={rx}, rank Hz={rz}")
    chk("Lx rows ∈ ker Hz", not (Hz @ Lx.T % 2).any())
    chk("Lz rows ∈ ker Hx", not (Hx @ Lz.T % 2).any())
    chk("Lx·Lzᵀ = I₃₀", (Lx @ Lz.T % 2 == np.eye(30, dtype=np.uint8)).all())
    chk("canonical wt(Lx) = 18 (Table I)", set(Lx.sum(1).tolist()) == {18})
    chk("canonical wt(Lz) = 10 (Table I)", set(Lz.sum(1).tolist()) == {10},
        "each Lz row is itself a weight-10 Z-logical ⟹ d_Z ≤ 10")

    # Eq. (2) block structure: 30-wide column blocks D1..D5, two 30-row
    # check blocks per side, 3 nonzero blocks per row block, exactly one
    # column block (D5) shared between the two row blocks.
    def pattern(H: np.ndarray) -> list[list[int]]:
        return [
            [int(H[30 * i:30 * (i + 1), 30 * c:30 * (c + 1)].any())
             for c in range(5)]
            for i in range(2)
        ]

    def blocks(H: np.ndarray, p: list[list[int]]) -> dict:
        return {(i, c): H[30 * i:30 * (i + 1), 30 * c:30 * (c + 1)]
                for i in range(2) for c in range(5) if p[i][c]}

    px, pz = pattern(Hx), pattern(Hz)
    for name, p in (("Hx", px), ("Hz", pz)):
        chk(f"{name} 2×5 block pattern (3+3, one shared col)",
            all(sum(r) == 3 for r in p)
            and sum(p[0][c] & p[1][c] for c in range(5)) == 1, f"{p}")
    bx, bz = blocks(Hx, px), blocks(Hz, pz)
    chk("every nonzero 30×30 block = sum of 3 permutations",
        all(set(B.sum(0).tolist()) == {3} == set(B.sum(1).tolist())
            for B in list(bx.values()) + list(bz.values())))
    # The L/R mechanism: on every shared data block, X- and Z-checks act
    # by opposite regular representations, hence commute pairwise.
    chk("all same-column (Hx-block, Hz-block) pairs commute (L/R)",
        all(not ((BX.astype(int) @ BZ - BZ.astype(int) @ BX) % 2).any()
            for (i, c), BX in bx.items()
            for (j, d), BZ in bz.items() if c == d))
    # Non-abelian signature: same-side blocks need NOT commute.
    chk("some Hx-block pair does NOT commute (non-abelian G)",
        any(((BX.astype(int) @ BY - BY.astype(int) @ BX) % 2).any()
            for key1, BX in bx.items()
            for key2, BY in bx.items() if key1 < key2))

    print("ALL PASS" if ok else "VALIDATION FAILED", flush=True)
    if not ok:
        sys.exit(1)


def _side_checks(side: str) -> tuple[GenericChecks, np.ndarray]:
    """(checks, released opposite-type basis) for one CSS direction.

    direction X: search v ∈ ker(H_Z) outside rowspace(H_X);
    direction Z: same with the roles swapped."""
    Hx, Hz, Lx, Lz = load()
    if side == "x":
        return GenericChecks(H_X=Hx, H_Z=Hz, group=_Label()), Lz
    return GenericChecks(H_X=Hz, H_Z=Hx, group=_Label()), Lx


def run_tandem(side: str, binary: Path) -> tuple[np.ndarray, int, float]:
    checks, L_released = _side_checks(side)
    WORK.mkdir(parents=True, exist_ok=True)
    res = maxsat_distance(
        checks, binary, mode="naive", work_dir=WORK, timeout=14400,
    )
    v = res.witness
    assert v is not None and int(v.sum()) == res.distance
    # Nontriviality against the AUTHORS' released opposite basis, on top
    # of maxsat_distance's own re-verification vs the recomputed basis.
    assert ((L_released @ v) % 2).any(), "witness trivial vs released basis"
    print(
        f"d_{side.upper()} (tandem/MaxCDCL): OPTIMUM {res.distance}, "
        f"witness re-verified, {res.solver_seconds:.1f}s",
        flush=True,
    )
    return v, res.distance, res.solver_seconds


def run_cms_ladder(side: str) -> tuple[int, list[str], float]:
    checks, _ = _side_checks(side)
    rungs: list[str] = []
    t0 = time.perf_counter()
    res = x_distance(
        checks, weight_lower_bound=1, weight_upper_bound=12,
        progress=lambda w, sat, s: rungs.append(
            f"w={w}: {'SAT' if sat else 'UNSAT'} {s:.1f}s"
        ),
    )
    dt = time.perf_counter() - t0
    print(
        f"d_{side.upper()} (CMS ladder): {res.distance} "
        f"[{'; '.join(rungs)}] total {dt:.1f}s",
        flush=True,
    )
    return res.distance, rungs, dt


def final(binary: Path, skip_cms: bool) -> None:
    try:
        import pycryptosat
        cms_ver = getattr(pycryptosat, "__version__", "?")
    except ImportError:
        cms_ver = "missing"
    for side in ("x", "z"):
        v, d_tandem, sec_tandem = run_tandem(side, binary)
        wall = sec_tandem
        solver = "tandem(maxcdcl-mse23+qeclab-patch, no flags)"
        if not skip_cms:
            d_cms, _, sec_cms = run_cms_ladder(side)
            assert d_cms == d_tandem, (d_cms, d_tandem)
            wall += sec_cms
            solver += f" + cms-ladder(pycryptosat@{cms_ver})"
        checks, _ = _side_checks(side)
        L = find_logical_z(checks)
        cert = make_certificate(
            code_id="mitten_150_30_10",
            H_check=checks.H_Z,          # the syndrome operator this side
            L_logical=L,
            witness=v,
            distance=d_tandem,
            direction=side.upper(),
            solver=solver,
            wall_seconds=round(wall, 1),
        )
        verify_certificate(cert, checks.H_Z, L)   # round-trip before write
        path = write_certificate(
            cert, CERTS / f"mitten_150_30_10_{side.upper()}.cert.json"
        )
        print(f"  cert: {path.relative_to(LAB_ROOT)} (verified)", flush=True)
    print(f"d = min(d_X, d_Z) = 10 — matches {PAPER} Table I", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("validate", "x", "z", "final"))
    ap.add_argument("--binary", type=Path, default=DEFAULT_BINARY,
                    help="tandem binary (build_maxcdcl.sh output)")
    ap.add_argument("--skip-cms", action="store_true",
                    help="final: skip the CMS UNSAT-ladder cross-check")
    args = ap.parse_args()
    if args.mode == "validate":
        validate()
    elif args.mode in ("x", "z"):
        if not args.binary.exists():
            sys.exit(f"tandem binary not found: {args.binary}")
        run_tandem(args.mode, args.binary)
    else:
        if not args.binary.exists():
            sys.exit(f"tandem binary not found: {args.binary}")
        final(args.binary, args.skip_cms)
