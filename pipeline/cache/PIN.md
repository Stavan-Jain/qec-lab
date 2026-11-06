# Pinned eczoo_data snapshot

The catalog is built from a pinned snapshot of
[github.com/errorcorrectionzoo/eczoo_data](https://github.com/errorcorrectionzoo/eczoo_data).

**Pinned SHA:** `6f4da6d3284db17edb4100adfd38b73b0766c6d4`

## Reproduce the existing catalog

From the repository root, recreate a missing cache at the exact recorded
revision (cloning the current default branch alone is not reproducible):

```bash
git clone --no-checkout --depth 1 https://github.com/errorcorrectionzoo/eczoo_data.git pipeline/cache/eczoo_data
git -C pipeline/cache/eczoo_data fetch --depth 1 origin 6f4da6d3284db17edb4100adfd38b73b0766c6d4
git -C pipeline/cache/eczoo_data checkout --detach 6f4da6d3284db17edb4100adfd38b73b0766c6d4
uv run --project experiments/bb_lab python scripts/ingest_zoo.py
```

If the cache already exists, omit the clone and first check that it has no
local changes. Keep the cache out of Git; retain this pin and the generated
catalog. This is distinct from intentionally refreshing to a newer revision.

## Refresh procedure

To update the catalog to a newer snapshot:

```bash
git -C pipeline/cache/eczoo_data status --short  # must be clean
git -C pipeline/cache/eczoo_data fetch --depth 1 origin HEAD
git -C pipeline/cache/eczoo_data checkout --detach FETCH_HEAD
git -C pipeline/cache/eczoo_data rev-parse HEAD  # update every SHA in this file
uv run --project experiments/bb_lab python scripts/ingest_zoo.py
```

After regeneration, re-score by spawning the `qec-prioritizer` agent.
