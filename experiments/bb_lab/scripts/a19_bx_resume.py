"""A19 step 1: targeted resume of the BX (15,6) floor@11 round.

Parses the partial log for already-UNSAT classes and queries only the
outstanding orbit reps at w <= 11. Combined with the banked 14 UNSATs and the
weight-12 ISD witness, all-35-UNSAT proves d(BX) = 12 exactly.
"""
import sys, re, time, importlib.util, json
import multiprocessing as mp
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import numpy as np
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices
from bb_lab.sat_distance import find_logical_z

PARTIAL_LOG = LAB / "data" / "a19" / "bx_floor11_partial.log"
OUT_LOG = LAB / "data" / "a19" / "bx_floor11_resume.log"


def load_acd():
    spec = importlib.util.spec_from_file_location(
        "a15_coset_distance", LAB / "scripts" / "a15_coset_distance.py")
    acd = importlib.util.module_from_spec(spec)
    sys.modules["a15_coset_distance"] = acd
    spec.loader.exec_module(acd)
    return acd


def main():
    acd = load_acd()
    G = ZmZn(15, 6)
    ch = bb_check_matrices(Poly.from_string("x^9 + y + y^2", G),
                           Poly.from_string("y^3 + x^10 + x^11", G))
    L_Z = find_logical_z(ch)
    S = acd.x_class_reps(ch)
    k = S.shape[0]
    mats, _ = acd.signature_action(G, ch, L_Z, S)
    reps = [sig for sig, _ in acd.orbit_reps(mats, list(G.orders), k)]

    done = {int(m, 16) for m in re.findall(
        r"class (0x[0-9a-f]+): UNSAT", PARTIAL_LOG.read_text())}
    todo = [sig for sig in reps if sig not in done]
    lines = [f"resume: {len(reps)} reps total, {len(done)} banked UNSAT, "
             f"{len(todo)} outstanding"]
    print(lines[-1], flush=True)

    sat_hits = []
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ctx.Pool(8) as pool:
        work = [(ch, L_Z, sig, 11) for sig in todo]
        for sig, v, secs in pool.imap_unordered(acd._query_worker, work):
            tag = "UNSAT" if v is None else f"SAT wt {int(v.sum())}"
            line = f"  [R w<=11] class {sig:#x}: {tag} ({secs:.1f}s)"
            lines.append(line)
            print(line, flush=True)
            if v is not None:
                sat_hits.append((sig, int(v.sum())))
                np.save(LAB / "data" / "a19" / f"bx_wit_{sig:#x}.npy", v)

    verdict = ("d(BX) >= 12 PROVEN (all 35 orbit reps UNSAT@11; with the "
               "weight-12 witness: d(BX) = 12 EXACT)" if not sat_hits
               else f"REFUTED at 12: witnesses {sat_hits}")
    summary = json.dumps({"banked": len(done), "resumed": len(todo),
                          "sat_hits": sat_hits, "verdict": verdict,
                          "secs": round(time.time() - t0, 1)})
    lines += [verdict, summary]
    print(verdict, flush=True)
    print(summary, flush=True)
    OUT_LOG.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
