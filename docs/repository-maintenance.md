# Repository storage and recovery

Keep research sources, certificates, result summaries, and negative results.
Dependency environments, compiler output, and large intermediate enumerations
have different storage requirements. A directory being ignored by Git does
**not** mean its contents can be discarded: some ignored experiment directories
also contain deliberately tracked evidence and unique local results.

## Cleanup result

The October 2026 audit reduced this checkout from **2.17 GiB to about
773 MiB (65.3%)**, including the **87.7 MiB local archive** retained inside
it. The reduction came from five inactive environments, lossless compression
of 2,326 raw output files, one duplicate array, and generated scratch files.
All worktrees, Git history, tracked evidence, and existing research edits
were preserved. Exact category counts and measurements are recorded in
`.research-archive/cleanup-summary.json`.

## Local archive

The cleanup of this checkout stores bulky local research outputs losslessly in
`.research-archive/`. This directory is ignored by Git and remains **inside this
checkout**. It is not a remote backup and does not appear in a fresh clone.
Preserve it when moving or backing up this checkout.

Each original path is recorded in `.research-archive/manifest.json`, with its
SHA-256, byte count, timestamps, permissions, producer, and restoration note.
Compressed objects are named by the original content hash, so identical
contents share one object. Every archive was decompressed and checked against
its original hash before the original file was removed.

From the repository root:

```bash
python3 scripts/research_storage.py list
python3 scripts/research_storage.py verify
python3 scripts/research_storage.py restore experiments/bb_lab/data/a42/s5_ckpt_W24_L2s1p.jsonl
```

Pass several exact paths to restore several files, or use `restore --all` to
restore all archived originals. Restoration recreates the original locations
and refuses to overwrite existing files. It also rejects path traversal and
symlink destinations. Consumers continue to use their original filenames;
restore their inputs before running them. The helper uses only Python's
standard library.

The archive covers these local outputs:

| Outputs | Why preserve them | When to restore |
| --- | --- | --- |
| Large A42 JSONL checkpoints | Re-running exact enumeration is expensive | Before resuming `a42_s5_mixed24.py` or `a42_s5_dangerous24.py` when the required checkpoint is absent |
| Historical DuckDB backups | They are distinct historical snapshots | To inspect or recover a prior database state; the current database remains in place |
| A35 raw `.hits` and `.mat` files in the high-FOM worktree | They retain raw enumeration data | Before resuming a historical run that needs its original raw inputs or hit streams |

Small summaries, verdicts, ledgers, tracked evidence, and the unique lifted-product
sweep data remain accessible. The A42 `s3_bwd12_table22.npz` name is a relative
symlink to the byte-identical `s3_bwd12_table.npz`; the canonical array remains
available to its existing consumers.

## Rebuildable environments and output

Five inactive dependency environments in historical worktrees were removed.
Their `pyproject.toml` and `uv.lock` files remain in place. The archive manifest
also records each removed environment's configuration, installed packages, and
lockfile hashes. Recreate one from its lab directory:

```bash
uv sync --frozen
```

Use `uv sync --frozen --extra dev` when that project's development dependencies are
needed. The main BB lab environment was retained. Historical worktrees, their
branches, and their uncommitted work were retained.

Python bytecode, pytest caches, LaTeX scratch files, and `dashboard/dist/` are
generated. Re-running their producer recreates them. Build the dashboard with:

```bash
python3 dashboard/build.py
```

The cleanup retained the main lab's `data/build/` kernel cache and the pinned
EC Zoo cache. `src/bb_lab/cosetbz.py` can rebuild its C kernel on demand. Follow
`pipeline/cache/PIN.md` when restoring the EC Zoo snapshot.

## Future cleanup

1. Check for running processes and identify the checkout that owns each path.
2. Check tracked **and nonignored untracked** files in that checkout. Never
   assume the main checkout's ignore rules describe another worktree.
3. Preserve unique research outputs. For large inactive outputs, create a
   lossless archive and verify its decompressed byte count and SHA-256 before
   removing originals.
4. Record the producer, consumers, original path, and restoration instructions.
5. Remove only verified generated caches or inactive environments with their
   dependency definitions preserved. Do not rewrite Git history to save space.

Avoid hardlinking writable source files between worktrees. Retiring a worktree
requires a separate recovery snapshot of its local work and ignored research
outputs; deleting its directory is not ordinary cache cleanup.
