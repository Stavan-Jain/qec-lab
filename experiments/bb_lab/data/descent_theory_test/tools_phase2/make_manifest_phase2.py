"""Hash all Phase-2/3/4 outputs into MANIFEST_phase2.sha256."""

from __future__ import annotations

import hashlib
import time
from pathlib import Path

DTT = Path(__file__).resolve().parent.parent

FILES = [
    "phase2_results.jsonl",
    "phase3_groundtruth.jsonl",
    "phase3_armB.jsonl",
    "phase3_gt_selection.json",
    "phase3_joins.jsonl",
    "PHASE4_SCORECARD.md",
    "tools_phase2/dtt_close_lib.py",
    "tools_phase2/dtt_close_worker.py",
    "tools_phase2/run_phase2.py",
    "tools_phase2/run_armB.py",
    "tools_phase2/run_groundtruth.py",
    "tools_phase2/controls.py",
    "tools_phase2/analyze.py",
    "tools_phase2/sat_worker.py",
    "tools_phase2/cross_stage_joins.py",
    "tools_phase2/watch_sweep.sh",
    "tools_phase2/analysis_final.json",
    "tools_phase2/make_manifest_phase2.py",
]


def main() -> None:
    lines = [f"# descent_theory_test Phase 2-4 outputs, hashed "
             f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}"]
    for rel in FILES:
        p = DTT / rel
        if not p.exists():
            lines.append(f"# MISSING: {rel}")
            continue
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        lines.append(f"{h}  {rel}")
    (DTT / "MANIFEST_phase2.sha256").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
