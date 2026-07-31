"""Task functions for `test_sweep.py`.

They live in an importable module rather than in the test file because
`ProcessPoolExecutor` pickles the task by qualified name and the spawned
child re-imports it; a function defined in the pytest `__main__` would
not resolve there.
"""

from __future__ import annotations

import os


def square_task(item: dict) -> dict:
    """Trivial deterministic task."""
    return {
        "key": item["key"],
        "value": item["x"] ** 2,
        "status": "ok",
        "pid": os.getpid(),
    }


def flaky_task(item: dict) -> dict:
    """Raises for one designated key, to exercise the error path."""
    if item["key"] == "boom":
        raise ValueError("intentional worker failure")
    return {
        "key": item["key"],
        "value": item["x"] ** 2,
        "status": "ok",
        "pid": os.getpid(),
    }


def fanout_task(item: dict) -> dict:
    """Rendezvous task: prove the pool really runs workers concurrently.

    Each task registers its pid in the shared work root, then waits until
    at least two distinct pids are registered. If the runner were serial
    the wait would time out and `saw` would stay at 1 — so this tests
    concurrency directly instead of inferring it from timing, which a
    microsecond-long task cannot do reliably.
    """
    import time

    from bb_lab.sweep import worker_work_dir

    shared = worker_work_dir().parent
    (shared / f"pid-{os.getpid()}").touch()

    deadline = time.monotonic() + 10.0
    saw = 1
    while time.monotonic() < deadline:
        saw = len(list(shared.glob("pid-*")))
        if saw >= 2:
            break
        time.sleep(0.01)

    return {
        "key": item["key"],
        "value": saw,
        "status": "ok",
        "pid": os.getpid(),
    }


def pair_task(item: dict) -> dict:
    """Keyed by an (a, b) pair rather than a single id column."""
    return {
        "a": item["a"],
        "b": item["b"],
        "value": item["x"] ** 2,
        "status": "ok",
    }


def slow_task(item: dict) -> dict:
    """Sleeps `item['sleep']` seconds — for exercising the deadline."""
    import time

    time.sleep(item.get("sleep", 0.2))
    return {
        "key": item["key"],
        "value": item["x"] ** 2,
        "status": "ok",
        "pid": os.getpid(),
    }


def workdir_task(item: dict) -> dict:
    """Reports the per-worker scratch directory the runner assigned."""
    from bb_lab.sweep import worker_work_dir

    wd = worker_work_dir()
    return {
        "key": item["key"],
        "value": 0,
        "status": "ok",
        "pid": os.getpid(),
        "workdir": str(wd),
    }
