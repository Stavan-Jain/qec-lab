"""Doubling-certificate lane worker (one code), run under a parent watchdog.

Runs bb_lab.doubling_certify.certify on the cover and prints the
verdict dict as JSON on stdout (last line, tagged CERTIFY_JSON).
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def main() -> None:
    os.nice(10)
    ap = argparse.ArgumentParser()
    ap.add_argument("--ell", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--A", required=True)
    ap.add_argument("--B", required=True)
    ap.add_argument("--budget", type=float, required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    from bb_lab.doubling_certify import certify, scrub_json

    out = certify(
        (args.ell, args.m),
        args.A,
        args.B,
        budget_s=args.budget,
        threads=args.threads,
        workdir=Path(args.workdir),
    )
    st = out.pop("stages", {}) or {}
    # keep a slim, load-bearing extract; full detail stays in workdir
    out["base"] = st.get("base")
    out["d_base"] = st.get("d_base")
    out["n_detect_candidates"] = len(st.get("detect", []) or [])
    out["candidate_log"] = st.get("candidate_log")
    print("CERTIFY_JSON " + json.dumps(scrub_json(out)), flush=True)


if __name__ == "__main__":
    main()
