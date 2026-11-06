#!/usr/bin/env python3
"""Regenerate the mitten instance and its kernel-checked certificates."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True,
                        help="QECLean checkout to receive the generated Lean files")
    parser.add_argument("--force", action="store_true",
                        help="allow replacing existing generated files")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    lab = here.parent.parent
    output = args.root.resolve() / "QEC/Stabilizer/Codes/Mitten/M150"
    if not args.force and (output / "FloorData.lean").exists():
        parser.error("generated data already exists; pass --force to regenerate")
    output.mkdir(parents=True, exist_ok=True)
    for mode in ("instance", "m4"):
        subprocess.run([sys.executable, str(here.parent / "m150_gen_lean_data.py"),
                        mode, "--out", str(output), "--force"], check=True)
    subprocess.run([sys.executable, str(here / "generate_joins.py"),
                    "--data", str(lab / "data-a32/m4_triples.json"),
                    "--out", str(output), "--force"], check=True)
    subprocess.run([sys.executable, str(here / "generate_sweeps.py"),
                    "--root", str(args.root.resolve()), "--all", "--force"], check=True)


if __name__ == "__main__":
    main()
