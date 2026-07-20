"""Phase 0.3 gate: try to beat weight 24 at the (30,6) cover via per-coset SAT
at w <= 22 on the lightest classes from the ISD census. Budget-bounded.
Any SAT hit => d(cover) <= hit weight, rerouting the certification target.
"""
import sys, time, importlib.util, json
import multiprocessing as mp
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
SCRATCH = Path(__file__).resolve().parent.parent / "data" / "a19"
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import numpy as np
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices
from bb_lab.sat_distance import find_logical_z


def load_acd():
    spec = importlib.util.spec_from_file_location(
        "a15_coset_distance", LAB / "scripts" / "a15_coset_distance.py")
    acd = importlib.util.module_from_spec(spec)
    sys.modules["a15_coset_distance"] = acd
    spec.loader.exec_module(acd)
    return acd


def main():
    acd = load_acd()
    GC = ZmZn(30, 6)
    chC = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GC),
                            Poly.from_string("y^3 + x^25 + x^26", GC))
    L_ZC = find_logical_z(chC)
    SC = acd.x_class_reps(chC)
    k = SC.shape[0]

    data = np.load(SCRATCH / "isd_class_minima.npz")
    sigs, wts = data["sigs"], data["weights"]
    # orbit-dedupe the lightest classes
    mats, _ = acd.signature_action(GC, chC, L_ZC, SC)
    reps = acd.orbit_reps(mats, list(GC.orders), k)
    rep_set = {sig for sig, _ in reps}
    order = np.argsort(wts)
    targets = []
    for i in order:
        sig = int(sigs[i])
        if sig in rep_set:
            targets.append((sig, int(wts[i])))
        if len(targets) >= 10:
            break
    print(f"targets (orbit reps, lightest ISD): {[(hex(s), w) for s, w in targets]}",
          flush=True)

    ctx = mp.get_context("spawn")
    deadline = time.time() + 3 * 3600
    results = []
    with ctx.Pool(6) as pool:
        async_res = [(sig, w0, pool.apply_async(
            acd._query_worker, ((chC, L_ZC, sig, 22),)))
            for sig, w0 in targets]
        for sig, w0, ar in async_res:
            left = deadline - time.time()
            if left <= 0:
                results.append((sig, w0, "TIMEOUT"))
                print(f"  class {sig:#05x} (isd {w0}): TIMEOUT", flush=True)
                continue
            try:
                _, v, secs = ar.get(timeout=left)
                out = f"SAT wt {int(v.sum())}" if v is not None else "UNSAT@22"
                if v is not None:
                    np.save(SCRATCH / f"beat24_wit_{sig:#05x}.npy", v)
                results.append((sig, w0, out))
                print(f"  class {sig:#05x} (isd {w0}): {out} ({secs:.0f}s)",
                      flush=True)
            except mp.TimeoutError:
                results.append((sig, w0, "TIMEOUT"))
                print(f"  class {sig:#05x} (isd {w0}): TIMEOUT", flush=True)
        pool.terminate()
    print(json.dumps({"beat24": [[f"{s:#05x}", w, o] for s, w, o in results]}),
          flush=True)


if __name__ == "__main__":
    main()
