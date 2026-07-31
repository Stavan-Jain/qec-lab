"""DRAT proof emission via the `cadical` CLI + `drat-trim` independent
verification.

This is the v1-track (option b) soundness check. The pipeline:

  (CNF, distance bound) → cadical CLI → DRAT proof + DIMACS CNF
                                              ↓
                                     drat-trim verifier → 's VERIFIED'

When this test passes we know:
  - the proofs we emit are well-formed DRAT
  - they actually refute the CNF we built
  - the SAT distance algorithm is sound on the UNSAT direction
    (and not just trusting the solver's verdict)

Together with the SAT witness for the SAT direction, this gives a
checkable distance certificate. Lean-side consumption (via a verified
LRAT checker like Mario Carneiro's `lrat-check`) is the next step:
just run `drat-trim -L <out.lrat>` to convert DRAT→LRAT and feed Lean.

Tests are skipped if the `cadical` binary or `drat-trim` is not on PATH.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

from bb_lab.certificate import make_certificate, verify_certificate, WitnessCertificate
from bb_lab.checks import bb_check_matrices
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.sat_distance import find_logical_z, x_distance


CADICAL = shutil.which("cadical")
DRAT_TRIM = shutil.which("drat-trim") or "/tmp/drat-trim/drat-trim"
HAS_DRAT_TRIM = Path(DRAT_TRIM).exists()


pytestmark = pytest.mark.skipif(
    CADICAL is None,
    reason="cadical CLI not on PATH; install via `brew install cadical`",
)


def _bb_72_checks():
    G = ZmZn(6, 6)
    A = Poly.from_string("x^3 + y + y^2", G)
    B = Poly.from_string("y^3 + x + x^2", G)
    return bb_check_matrices(A, B)


def test_drat_emitted_for_each_unsat_weight():
    checks = _bb_72_checks()
    with tempfile.TemporaryDirectory() as td:
        result = x_distance(checks, proof_dir=td, code_id="bb_72_12_6")
        assert result.distance == 6
        # 5 UNSAT calls (w=1..5), then 1 SAT (w=6) — so 5 proofs
        assert len(result.unsat_proof_paths) == 5
        for w, proof_path in enumerate(result.unsat_proof_paths, start=1):
            cnf_path = proof_path.with_suffix(".cnf")
            assert proof_path.exists(), f"DRAT for w={w} missing"
            assert cnf_path.exists(), f"CNF for w={w} missing"
            assert proof_path.stat().st_size > 0, f"DRAT for w={w} empty"
            # Sanity-check CNF header
            with cnf_path.open() as f:
                header = f.readline()
            assert header.startswith("p cnf "), (
                f"w={w}: bad DIMACS header {header!r}"
            )
            # DRAT proofs from cadical --no-binary end with the empty-clause
            # derivation `0`
            with proof_path.open() as f:
                lines = f.readlines()
            assert lines[-1].strip() == "0", (
                f"w={w}: DRAT proof must end with empty-clause line '0', "
                f"got {lines[-1]!r}"
            )


def test_certificate_records_drat_paths():
    """End-to-end: SAT distance run → certificate with DRAT refs → JSON
    round-trip preserves the proof references with stable hashes."""
    checks = _bb_72_checks()
    L_Z = find_logical_z(checks)
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        result = x_distance(checks, proof_dir=td_path, code_id="bb_72_12_6")
        cert = make_certificate(
            code_id="bb_72_12_6",
            H_check=checks.H_Z,
            L_logical=L_Z,
            witness=result.witness,
            distance=result.distance,
            direction="X",
            solver="cadical@3.0.0",
            wall_seconds=0.5,
            unsat_drat_paths=result.unsat_proof_paths,
            cert_dir=td_path,
        )
        assert len(cert.unsat_proofs) == 5
        weights = sorted(p.weight_bound for p in cert.unsat_proofs)
        assert weights == [1, 2, 3, 4, 5]
        # JSON round-trip preserves hashes + paths
        again = WitnessCertificate.from_json(cert.to_json())
        assert again.unsat_proofs == cert.unsat_proofs

        # The certificate's own verify_certificate still passes — the
        # witness-level checks are independent of the DRAT material.
        verify_certificate(cert, checks.H_Z, L_Z)


def test_committed_cnfs_are_reproducible():
    """Re-emitting bb_72's proof bundle must reproduce the committed
    CNFs byte for byte.

    `certificates/.gitignore` commits only bb_72's bundle and justifies
    discarding the larger ones with "regeneration reproduces the same
    files". That claim quietly stopped holding when the 2026-07-29
    `quotient_complement_basis` rewrite changed which logical basis
    `find_logical_z` returns — the basis rows are encoded directly into
    the CNF, so every committed CNF became unreproducible with nothing
    to notice. This is that missing check.

    Only the CNFs are compared: they are built by our own encoder and so
    are fully determined by the code, whereas DRAT bytes depend on the
    solver build and would make this test fail on a different cadical.
    """
    committed = Path(__file__).resolve().parent.parent / "certificates" / "bb_72_12_6"
    if not committed.exists():
        pytest.skip("committed bb_72 proof bundle missing")
    checks = _bb_72_checks()
    with tempfile.TemporaryDirectory() as td:
        result = x_distance(checks, proof_dir=td, code_id="bb_72_12_6")
        for proof_path in result.unsat_proof_paths:
            fresh = proof_path.with_suffix(".cnf")
            gold = committed / fresh.name
            assert gold.exists(), f"committed bundle is missing {fresh.name}"
            assert fresh.read_bytes() == gold.read_bytes(), (
                f"{fresh.name} no longer reproduces the committed bundle — "
                "the CNF encoding (or the logical basis it encodes) changed; "
                "regenerate with `bb-lab bravyi-check --quick --emit-proofs`"
            )


def test_drat_sizes_grow_with_weight():
    """Cheap sanity: weight-(d-1) UNSAT should produce a meaningfully
    larger proof than weight-1."""
    checks = _bb_72_checks()
    with tempfile.TemporaryDirectory() as td:
        result = x_distance(checks, proof_dir=td, code_id="bb_72_12_6")
        first = result.unsat_proof_paths[0]
        last = result.unsat_proof_paths[-1]
        assert first.stat().st_size < last.stat().st_size, (
            "expected hardest UNSAT (w=d-1) to produce the largest proof"
        )


@pytest.mark.skipif(
    not HAS_DRAT_TRIM,
    reason="drat-trim binary not found (build from "
           "github.com/marijnheule/drat-trim then put on PATH or at "
           "/tmp/drat-trim/drat-trim)",
)
def test_drat_proofs_verified_by_drat_trim():
    """**The soundness gate.** Every DRAT proof we emit is independently
    accepted by drat-trim as a valid refutation of its companion CNF.

    Without this, the UNSAT direction of the distance bound is "the
    solver said so" — with it, an independent verifier confirms the
    proof, and we are one tool-swap away (drat-trim -L → LRAT → Lean
    LRAT checker) from a Lean-kernel-checked distance bound.
    """
    checks = _bb_72_checks()
    with tempfile.TemporaryDirectory() as td:
        result = x_distance(checks, proof_dir=td, code_id="bb_72_12_6")
        for w, drat_path in enumerate(result.unsat_proof_paths, start=1):
            cnf_path = drat_path.with_suffix(".cnf")
            proc = subprocess.run(
                [DRAT_TRIM, str(cnf_path), str(drat_path)],
                capture_output=True,
                text=True,
                timeout=120,
            )
            # drat-trim writes 's VERIFIED' on a clean line for accepted
            # proofs, 's NOT VERIFIED' for rejected ones. (Exit code is
            # not consistently informative.)
            assert "s VERIFIED" in proc.stdout, (
                f"w={w}: drat-trim REJECTED the proof.\n"
                f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
            )
            assert "s NOT VERIFIED" not in proc.stdout, (
                f"w={w}: drat-trim says NOT VERIFIED:\n{proc.stdout}"
            )
