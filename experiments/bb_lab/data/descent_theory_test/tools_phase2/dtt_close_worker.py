"""Phase-2 closure worker: one pre-registered question, one process.

Reads a JSON spec file (argv[1]), runs dtt_close_lib.close_question under
os.nice(10), and prints the result as a single line tagged CLOSE_JSON
(plus streamed progress lines).  The parent enforces the hard budget with
a real kill; the in-process deadline is a courtesy check so stages can
stop cleanly and report BUDGET_KILL themselves when possible.
"""

from __future__ import annotations

import json
import os
import resource
import sys
import time


def main() -> None:
    os.nice(10)
    spec = json.loads(open(sys.argv[1]).read())
    t0 = time.time()
    from dtt_close_lib import close_question
    try:
        out = close_question(spec)
    except AssertionError as e:
        out = {"outcome": "ASSERTION", "reason": str(e)[:500],
               "tag": spec.get("tag")}
    except Exception as e:  # noqa: BLE001 — recorded, never silent
        out = {"outcome": "ERROR",
               "reason": f"{type(e).__name__}: {e}"[:500],
               "tag": spec.get("tag")}
    ru = resource.getrusage(resource.RUSAGE_SELF)
    rc = resource.getrusage(resource.RUSAGE_CHILDREN)
    out["cpu_s"] = round(ru.ru_utime + ru.ru_stime
                         + rc.ru_utime + rc.ru_stime, 2)
    out["wall_total_s"] = round(time.time() - t0, 2)
    out["maxrss_mb"] = round((ru.ru_maxrss) / 1e6, 1)
    print("CLOSE_JSON " + json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
