"""Write MANIFEST.sha256: SHA-256 of every file under
data/descent_theory_test/ (except the manifest itself), sorted by
relative path, with an ISO-8601 UTC timestamp header. Re-run after any
legitimate re-freeze; the manifest is the immutability anchor
(PROTOCOL.md §7)."""

from __future__ import annotations

import datetime as dt
import hashlib
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent
MANIFEST = OUT / "MANIFEST.sha256"


def main() -> None:
    lines = [
        "# descent_theory_test MANIFEST -- frozen "
        + dt.datetime.now(dt.UTC).isoformat(),
        "# scope: Phase 0 cohort + Phase 1 frozen predictions "
        "(batch 1); batch 2 gets its own MANIFEST_batch2.sha256",
    ]
    files = sorted(p for p in OUT.rglob("*")
                   if p.is_file() and p != MANIFEST
                   and "__pycache__" not in p.parts)
    for p in files:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        lines.append(f"{h}  {p.relative_to(OUT)}")
    MANIFEST.write_text("\n".join(lines) + "\n")
    print(f"{len(files)} files -> {MANIFEST}")
    print("manifest self-hash:",
          hashlib.sha256(MANIFEST.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
