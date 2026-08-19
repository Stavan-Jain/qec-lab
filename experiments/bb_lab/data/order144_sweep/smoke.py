"""Smoke tests for the order-144 sweep lanes. Run from experiments/bb_lab:
  uv run python data/order144_sweep/smoke.py
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import sweeplib  # noqa: E402
from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.canonical import (  # noqa: E402
    build_perm_table, canonical_bits, _bits_to_support,
)
from bb_lab.automorphism import automorphisms  # noqa: E402
from bb_lab.store import canonical_hash  # noqa: E402

MAIN_DB = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb")


def t1_canonical_equivalence() -> None:
    print("== T1: fast canonicalization == library canonical_bits ==")
    for (L, M) in [(6, 6), (12, 3), (6, 4)]:
        G = AbelianGroup((L, M))
        N = L * M
        t0 = time.time()
        perms = build_perm_table(G, auts=automorphisms(G))
        t_build = time.time() - t0
        rng = np.random.default_rng(7)
        n_mismatch = 0
        for trial in range(60):
            A_idx = np.sort(rng.choice(N, 3, replace=False))
            B_idx = np.sort(rng.choice(N, 3, replace=False))
            A_bits = int(sum(1 << int(i) for i in A_idx))
            B_bits = int(sum(1 << int(i) for i in B_idx))
            cA, cB, orb = canonical_bits(A_bits, B_bits, perms)
            fA, fB, forb = sweeplib.canonical_pair_fast(L, M, A_idx, B_idx)
            lib_A = tuple(sorted((i for i in range(N) if cA >> i & 1), reverse=True))
            lib_B = tuple(sorted((i for i in range(N) if cB >> i & 1), reverse=True))
            if (fA, fB, forb) != (lib_A, lib_B, orb):
                n_mismatch += 1
                print(f"  MISMATCH {L}x{M}: {A_idx} {B_idx}: fast {fA},{fB},{forb} lib {lib_A},{lib_B},{orb}")
        print(f"  Z{L}xZ{M}: 60 trials, {n_mismatch} mismatches (lib table build {t_build:.1f}s)")
        assert n_mismatch == 0


def t2_corpus_roundtrip() -> None:
    print("== T2: Z12xZ12 corpus row canonical round-trip ==")
    import duckdb
    con = duckdb.connect(str(MAIN_DB), read_only=True)
    rows = con.execute(
        "SELECT instance_id, A_poly, B_poly FROM bb_instances "
        "WHERE group_struct='Z12xZ12' AND d_exact IS NOT NULL"
    ).fetchall()
    con.close()
    t0 = time.time()
    sweeplib.group_tables(12, 12)
    print(f"  Z12xZ12 tables built in {time.time() - t0:.1f}s "
          f"(|Aut| = {sweeplib.group_tables(12, 12).PHI.shape[0]})")
    n_bad = 0
    t0 = time.time()
    for iid, A_str, B_str in rows:
        A_supp = sweeplib.support_from_string((12, 12), A_str)
        B_supp = sweeplib.support_from_string((12, 12), B_str)
        label, cA, cB, cid, orb = sweeplib.canonicalize(12, 12, A_supp, B_supp)
        if cid != iid:
            n_bad += 1
            print(f"  MISMATCH {iid}: got {cid} ({cA} ; {cB})")
    dt = time.time() - t0
    print(f"  {len(rows)} corpus rows: {len(rows) - n_bad} round-trip OK "
          f"({dt / max(len(rows), 1):.3f}s/row)")
    assert n_bad == 0


def t3_l1_timing() -> None:
    print("== T3: L1 d_ub timing at n=288 ==")
    from bb_lab.checks import bb_check_matrices
    from bb_lab.l1_sampling import l1_distance_ub
    G = AbelianGroup((12, 12))
    A = Poly.from_string("x^3 + y^2 + y^7", G)
    B = Poly.from_string("y^3 + x + x^2", G)
    checks = bb_check_matrices(A, B)
    for ns in (20000, 40000):
        t0 = time.time()
        res = l1_distance_ub(checks, n_samples=ns, seed=1)
        print(f"  [[288,12,18]] n_samples={ns}: d_ub={res.distance_ub} in {time.time() - t0:.1f}s")


def t4_sat_smoke() -> None:
    print("== T4: SAT worker smoke (gross quotient bb72 at n=72) ==")
    # [[72,12,6]] = gross's Z6xZ6 base; known d=6.
    prog = HERE / "smoke_sat_progress.txt"
    prog.unlink(missing_ok=True)
    t0 = time.time()
    p = subprocess.run(
        [sys.executable, str(HERE / "sat_worker.py"), "--ell", "6", "--m", "6",
         "--A", "x^3 + y + y^2", "--B", "y^3 + x + x^2",
         "--progress", str(prog)],
        capture_output=True, text=True, timeout=120,
    )
    print(f"  rc={p.returncode} out={p.stdout.strip()!r} ({time.time() - t0:.1f}s)")
    print(f"  progress: {prog.read_text().strip().splitlines()}")
    assert "DISTANCE 6" in p.stdout


def t5_certify_smoke() -> None:
    print("== T5: certify lane smoke (gross) ==")
    wd = HERE / "smoke_certify"
    t0 = time.time()
    p = subprocess.run(
        [sys.executable, str(HERE / "certify_worker.py"), "--ell", "12", "--m", "6",
         "--A", "x^3 + y + y^2", "--B", "y^3 + x + x^2",
         "--budget", "90", "--workdir", str(wd)],
        capture_output=True, text=True, timeout=150,
    )
    dt = time.time() - t0
    line = [l for l in p.stdout.splitlines() if l.startswith("CERTIFY_JSON")]
    print(f"  rc={p.returncode} ({dt:.1f}s)")
    if line:
        import json
        out = json.loads(line[0][len("CERTIFY_JSON "):])
        print(f"  status={out.get('status')} distance={out.get('distance')} wall={out.get('wall_s')}")
    else:
        print(f"  NO JSON; stderr tail: {p.stderr[-500:]}")


def t6_quotient_ladder() -> None:
    print("== T6: quotient ladder on bb288 [[288,12,18]] ==")
    A_supp = sweeplib.support_from_string((12, 12), "x^3 + y^2 + y^7")
    B_supp = sweeplib.support_from_string((12, 12), "y^3 + x + x^2")
    chain = sweeplib.quotient_chain((12, 12))
    print(f"  chain: {chain}")
    cur_o, cur_A, cur_B = (12, 12), A_supp, B_supp
    for q in chain:
        cur_A, cur_B = sweeplib.quotient_code(cur_o, cur_A, cur_B, q)
        p, checks, _, _ = sweeplib.code_params_of(q, cur_A, cur_B)
        print(f"  {q}: n={p.n} k={p.k} A={sweeplib.poly_string(q, cur_A)!r} B={sweeplib.poly_string(q, cur_B)!r}")
        cur_o = q


if __name__ == "__main__":
    which = sys.argv[1:] or ["t1", "t2", "t3", "t4", "t5", "t6"]
    for t in which:
        globals()[{"t1": "t1_canonical_equivalence", "t2": "t2_corpus_roundtrip",
                   "t3": "t3_l1_timing", "t4": "t4_sat_smoke",
                   "t5": "t5_certify_smoke", "t6": "t6_quotient_ladder"}[t]]()
