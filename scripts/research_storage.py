#!/usr/bin/env python3
"""List, verify, and restore losslessly archived local research outputs.

Archives are local to this checkout and are not required for a fresh clone.
Restore recreates original paths without overwriting existing files.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / ".research-archive"
HASH = re.compile(r"[0-9a-f]{64}")


def safe_path(name: str) -> Path:
    """Reject traversal and symlink ancestors before accessing a repo path."""
    parts = PurePosixPath(name).parts
    if (
        not parts
        or PurePosixPath(name).is_absolute()
        or any(part in {"..", ".git", ".research-archive"} for part in parts)
        or "\\" in name
        or "\x00" in name
        or str(PurePosixPath(name)) != name
    ):
        raise ValueError(f"Unsafe original path: {name!r}")
    path = ROOT.joinpath(*parts)
    for candidate in [path, *path.parents]:
        if candidate == ROOT:
            break
        if candidate.is_symlink():
            raise ValueError(f"Refusing symlink path: {candidate}")
    if not path.resolve().is_relative_to(ROOT):
        raise ValueError(f"Original path escapes checkout: {name}")
    return path


def load_manifest() -> dict:
    if ARCHIVE.is_symlink():
        raise ValueError("Refusing a symlink archive directory")
    manifest_path = ARCHIVE / "manifest.json"
    if manifest_path.is_symlink():
        raise ValueError("Refusing a symlink manifest")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("version") != 1 or not isinstance(manifest.get("files"), list):
        raise ValueError("Unsupported archive manifest")
    seen = set()
    for entry in manifest["files"]:
        safe_path(entry["path"])
        if entry["path"] in seen:
            raise ValueError(f"Duplicate manifest path: {entry['path']}")
        seen.add(entry["path"])
        if not HASH.fullmatch(entry["sha256"]):
            raise ValueError("Invalid SHA-256 in manifest")
        if not isinstance(entry["bytes"], int) or entry["bytes"] < 0:
            raise ValueError("Invalid original byte count")
        if entry["object"] != f"objects/{entry['sha256']}.gz":
            raise ValueError("Archive objects must be addressed by SHA-256")
        object_path = ARCHIVE / entry["object"]
        if (ARCHIVE / "objects").is_symlink() or object_path.is_symlink():
            raise ValueError("Refusing symlink archive object")
        if not object_path.resolve().is_relative_to(ARCHIVE.resolve()):
            raise ValueError("Archive object escapes archive directory")
    return manifest


def check_object(entry: dict, output=None) -> None:
    digest = hashlib.sha256()
    size = 0
    with gzip.open(ARCHIVE / entry["object"], "rb") as source:
        while chunk := source.read(1024 * 1024):
            size += len(chunk)
            if size > entry["bytes"]:
                raise ValueError(f"Archive is larger than declared: {entry['path']}")
            digest.update(chunk)
            if output is not None:
                output.write(chunk)
    if size != entry["bytes"] or digest.hexdigest() != entry["sha256"]:
        raise ValueError(f"Archive verification failed: {entry['path']}")


def restore(entry: dict) -> None:
    destination = safe_path(entry["path"])
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite {entry['path']}")
    check_object(entry)
    destination.parent.mkdir(parents=True, exist_ok=True)
    safe_path(entry["path"])
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=".restore-", delete=False) as output:
            temp_path = Path(output.name)
            check_object(entry, output)
            output.flush()
            os.fsync(output.fileno())
        os.chmod(temp_path, entry.get("mode", 0o644) & 0o777)
        if "mtime_ns" in entry:
            os.utime(temp_path, ns=(entry["mtime_ns"], entry["mtime_ns"]))
        # link() fails when destination already exists, unlike replace().
        os.link(temp_path, destination)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("list", help="List original paths and archive sizes")
    commands.add_parser("verify", help="Verify every unique compressed object")
    restore_parser = commands.add_parser("restore", help="Restore original paths without overwriting")
    restore_parser.add_argument("paths", nargs="*", help="Exact original path(s), as shown by list")
    restore_parser.add_argument("--all", action="store_true", help="Restore every archived path")
    args = parser.parse_args()
    try:
        manifest = load_manifest()
        entries = manifest["files"]
        if args.command == "list":
            for entry in entries:
                print(f"{entry['bytes']:>10}  {entry['category']:<18}  {entry['path']}")
            print(f"{len(entries)} original files, {len({entry['object'] for entry in entries})} compressed objects")
        elif args.command == "verify":
            verified = set()
            for entry in entries:
                if entry["object"] not in verified:
                    check_object(entry)
                    verified.add(entry["object"])
            print(f"Verified {len(verified)} objects covering {len(entries)} original files.")
        else:
            if args.all == bool(args.paths):
                parser.error("Use exact paths or --all, but not both")
            by_path = {entry["path"]: entry for entry in entries}
            requested = list(by_path) if args.all else list(dict.fromkeys(args.paths))
            unknown = set(requested) - by_path.keys()
            if unknown:
                raise ValueError(f"Unknown archived paths: {', '.join(sorted(unknown))}")
            selected = [by_path[path] for path in requested]
            for entry in selected:
                destination = safe_path(entry["path"])
                if destination.exists():
                    raise FileExistsError(f"Refusing to overwrite {entry['path']}")
                check_object(entry)
            for entry in selected:
                restore(entry)
                print(f"Restored {entry['path']}")
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"research_storage: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
