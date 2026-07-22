"""A23 phase 4: complete census of weight-8 base cycles by MITM.

Cycle condition (repo convention): d1 c = B*c_L + A*c_R = 0, i.e.
    B*c_L = A*c_R.
Enumerate all c_R supports of size r, hash the signature A*c_R; stream all
c_L supports of size l = 8 - r and look up B*c_L.  Splits covered:
(4,4), (5,3), (3,5).  (Lopsided (6,2)/(7,1)/(8,0) and their mirrors are
size-201M+ streams; deferred unless needed -- noted as a completeness
caveat.)

Outputs per found cycle: split, support, boundary?/logical?, H1 class
(vs a fixed logical basis), pairing with the seam subspace, S-cleanliness.
Writes the census + the reflect-swapped detector pool to data/a23/.
"""

from __future__ import annotations

import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rank_f2, rref_f2  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_seam_calibration import (  # noqa: E402
    A_SUPP,
    B_SUPP,
    conv_matrix,
    seam_maps,
)
from a23_detector_packing import reflect_swap  # noqa: E402
from a23_targeted_packing import S_GENS, is_s_clean  # noqa: E402


def pack75(v: np.ndarray) -> int:
    return int.from_bytes(np.packbits(v).tobytes(), "big")


def main() -> None:
    t00 = time.time()
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    K = nullspace_f2(D2)
    zeta = K[0]
    _, c1 = seam_maps(zeta)

    # per-generator packed signatures
    colA = [pack75(MA[:, g]) for g in range(75)]
    colB = [pack75(MB[:, g]) for g in range(75)]

    found: list[tuple[tuple[int, ...], tuple[int, ...]]] = []

    def mitm(l: int, r: int) -> None:
        t0 = time.time()
        # hash the smaller side
        if l <= r:
            table: dict[int, list[tuple[int, ...]]] = {}
            for cl in combinations(range(75), l):
                sig = 0
                for g in cl:
                    sig ^= colB[g]
                table.setdefault(sig, []).append(cl)
            cnt = 0
            for cr in combinations(range(75), r):
                sig = 0
                for g in cr:
                    sig ^= colA[g]
                for cl in table.get(sig, ()):
                    found.append((cl, cr))
                    cnt += 1
        else:
            table = {}
            for cr in combinations(range(75), r):
                sig = 0
                for g in cr:
                    sig ^= colA[g]
                table.setdefault(sig, []).append(cr)
            cnt = 0
            for cl in combinations(range(75), l):
                sig = 0
                for g in cl:
                    sig ^= colB[g]
                for cr in table.get(sig, ()):
                    found.append((cl, cr))
                    cnt += 1
        print(f"  split ({l},{r}): {cnt} cycles  [{time.time()-t0:.0f}s]")

    print("MITM census of weight-8 cycles:")
    mitm(4, 4)
    mitm(5, 3)
    mitm(3, 5)

    # classify
    Rb, pivb = rref_f2(D2.T)
    bnd_basis = Rb[: len(pivb)]
    rk_b = len(pivb)

    # primal logical basis for class coordinates
    Zprim = nullspace_f2(D1)
    cur = bnd_basis.copy()
    log_basis = []
    for row in Zprim:
        if rank_f2(np.vstack([cur, row])) > rank_f2(cur):
            cur = np.vstack([cur, row])
            log_basis.append(row)
    log_basis = np.array(log_basis, dtype=np.uint8)

    # class coordinates via dual pairing: use dual logical basis as functionals
    Zdual = nullspace_f2(D2.T)
    Rd, pivd = rref_f2(D1)
    dual_bnd = Rd[: len(pivd)]
    cur = dual_bnd.copy()
    dual_log_basis = []
    for row in Zdual:
        if rank_f2(np.vstack([cur, row])) > rank_f2(cur):
            cur = np.vstack([cur, row])
            dual_log_basis.append(row)
    dual_log_basis = np.array(dual_log_basis, dtype=np.uint8)

    # NOTE: pairing with dual_log_basis separates classes only up to the
    # radical of the pairing on this basis; for detector purposes the
    # pairing profile is exactly what matters.
    logicals = []
    n_bnd = 0
    class_profiles: dict[tuple[int, ...], int] = {}
    for (cl, cr) in found:
        v = np.zeros(150, dtype=np.uint8)
        for g in cl:
            v[g] = 1
        for g in cr:
            v[75 + g] = 1
        if rank_f2(np.vstack([bnd_basis, v])) == rk_b:
            n_bnd += 1
            continue
        prof = tuple(int((v & z).sum() % 2) for z in dual_log_basis)
        class_profiles[prof] = class_profiles.get(prof, 0) + 1
        logicals.append((v, prof))
    print(f"total: {len(found)} weight-8 cycles = {n_bnd} boundaries + "
          f"{len(logicals)} logicals")
    print(f"distinct pairing profiles among logicals: {len(class_profiles)}")
    for prof, cnt in sorted(class_profiles.items(), key=lambda kv: -kv[1])[:20]:
        print(f"    profile {''.join(map(str, prof))}: {cnt}")

    # detector pool: reflect-swap of logicals, keep odd pairing with c1
    dets = []
    n_sclean = 0
    for v, _ in logicals:
        z = reflect_swap(v)
        assert not ((D2.T @ z) % 2).any()
        if int((z & c1).sum()) % 2 == 1:
            dets.append(z)
            if is_s_clean(z):
                n_sclean += 1
    print(f"weight-8 detectors (odd pairing with c1): {len(dets)} "
          f"({n_sclean} S-clean)")

    outp = Path(__file__).resolve().parents[1] / "data/a23/w8_census.json"
    outp.parent.mkdir(exist_ok=True)
    outp.write_text(json.dumps({
        "splits": ["(4,4)", "(5,3)", "(3,5)"],
        "n_cycles": len(found),
        "n_boundaries": n_bnd,
        "n_logicals": len(logicals),
        "logical_supports": [sorted(map(int, np.flatnonzero(v)))
                             for v, _ in logicals],
        "detector_supports_odd_c1": [sorted(map(int, np.flatnonzero(z)))
                                     for z in dets],
    }, indent=1))
    print(f"wrote {outp}  [{time.time()-t00:.0f}s total]")


if __name__ == "__main__":
    main()
