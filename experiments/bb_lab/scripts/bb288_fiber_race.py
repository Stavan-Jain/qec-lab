"""Time-boxed bb_288 race: Tandem naive+step vs naive+step+fiber-lb.

Legs run in parallel (single-threaded solver each), one stdout reader
thread per leg stamping wall-clock seconds onto every line, SIGKILL at
the budget. The interesting observables for a run that (expectedly)
does not finish: the o-line descent trajectory, the fiber-conflict
counter, and any leg that proves OPTIMUM inside the box.

Usage:
  uv run python scripts/bb288_fiber_race.py --dir $S/bb288 \
      --binary $S/bin/tandem_v4 --budget 3600
"""

from __future__ import annotations

import argparse
import subprocess
import threading
import time
from pathlib import Path


def run_leg(name: str, argv: list[str], log: Path, budget: float):
    t0 = time.perf_counter()
    with log.open("w") as f:
        f.write(f"# {' '.join(argv)}\n")
        proc = subprocess.Popen(
            argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True,
        )

        def reader():
            for line in proc.stdout:
                f.write(f"[{time.perf_counter()-t0:9.2f}] {line}")
                f.flush()

        th = threading.Thread(target=reader, daemon=True)
        th.start()
        try:
            proc.wait(timeout=budget)
        except subprocess.TimeoutExpired:
            proc.kill()
            f.write(f"# KILLED at budget {budget}s\n")
        th.join(timeout=5)
    return name


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--budget", type=float, default=3600)
    args = ap.parse_args()
    d = Path(args.dir)
    wcnf = d / "naive_Z12xZ12.wcnf"
    assert wcnf.exists()

    legs = [("control", [args.binary, "-cost-step=2", str(wcnf)])]
    for sig, lbl in (("60", "fiber60"), ("06", "fiber06")):
        flb = d / f"fiber_Z12xZ12_{sig}.flb"
        if flb.exists():
            legs.append(
                (lbl, [args.binary, "-cost-step=2",
                       f"-fiber-lb={flb}", str(wcnf)])
            )
    print(f"racing {len(legs)} legs, budget {args.budget}s", flush=True)
    threads = []
    for name, argv in legs:
        th = threading.Thread(
            target=run_leg,
            args=(name, argv, d / f"race_{name}.log", args.budget),
        )
        th.start()
        threads.append((name, th))
    for name, th in threads:
        th.join()
        print(f"leg {name} done", flush=True)

    for name, _ in legs:
        log = d / f"race_{name}.log"
        lines = log.read_text().splitlines()
        o_lines = [l for l in lines if "] o " in l]
        fib = [l for l in lines if "fiber conflicts" in l]
        s_lines = [l for l in lines if "] s " in l]
        print(f"--- {name}: {len(o_lines)} incumbents")
        for l in o_lines[-4:]:
            print("   ", l)
        if fib:
            print("   ", fib[-1])
        for l in s_lines:
            print("   ", l)


if __name__ == "__main__":
    main()
