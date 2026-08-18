"""A38 S1: headless reproduction driver for the banked closures.

Runs one of the three banked closure suites (A32 bravyi360 d=24,
A33 ibm288Y d=20, A36 bb288 d=18) end-to-end via the frozen session
scripts, logging per-script wall time and exit codes.  The scripts
carry their own falsify-first hard asserts (banked structure, exact
node counts, banked-SAT agreement); a nonzero exit from any of them
IS the finding.

Usage:  uv run python scripts/a38_repro_driver.py {a32|a33|a36}
Log:    caller redirects stdout (kept under data/a38/repro/).
"""

from __future__ import annotations

import subprocess
import sys
import time

SUITES = {
    "a32": ["a32_tower_slice", "a32_gb_census", "a32_subclosures",
            "a32_sectorAC_full", "a32_deep_fibers", "a32_dby_floor"],
    "a33": ["a33_tower_cells", "a33_validate_banked", "a33_h5_close",
            "a33_h5_descent", "a33_solver_free"],
    # a36_direct_close re-verifies the banked witness end-to-end from its
    # support, then a36_witness re-runs the exhaustive find-side ladder on
    # the fresh censuses (overwriting the witness bank with what it finds).
    "a36": ["a36_tower_cells", "a36_direct_close", "a36_witness",
            "a36_descent"],
}


def main() -> None:
    suite = sys.argv[1]
    t0 = time.monotonic()
    for s in SUITES[suite]:
        print(f"##### RUN {s} at t={time.monotonic()-t0:.1f}s", flush=True)
        r = subprocess.run([sys.executable, f"scripts/{s}.py"])
        print(f"##### EXIT {s} rc={r.returncode} "
              f"at t={time.monotonic()-t0:.1f}s", flush=True)
        if r.returncode != 0:
            sys.exit(r.returncode)
    print(f"##### {suite.upper()} SUITE DONE {time.monotonic()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    main()
