#!/usr/bin/env python3
"""A42 S3c — the corner census: what survives the corridor at
levels 19..22 (the Stage-2a near-miss data).

Reads the corridor frontier dumps s3_corridor12_{T,V}_frontier_L*.npz
(the kept novel states at each level >= 19) plus the backward
closure-distance table, and characterizes the survivors:
  - popcount (window weight) distribution,
  - closure-distance G distribution (found / absent),
  - omega-slot structure: how many of the 5 window columns carry
    nonzero omega-content, barren-only content, or joint content —
    the halving-lemma shadow (columns with joint content are where
    tab can undercut pure),
  - distinct omega- and barren-window projection counts.

These states are exactly the (weight <= 22)-feasible corridor at
the deepest levels; with the run's verdict (no nonzero-register
returns) they are all segments of TRIVIAL near-cycles — the census
shows how close trivial traffic gets to the floor and through which
column shapes, feeding the hiding-mass analysis (§2.11).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a42_lib as AL  # noqa: E402
import importlib.util


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CR = _load("a42_s3_corridor_racer")
DATA = LAB / "data" / "a42"
P = 12


def main():
    z = np.load(DATA / "s3_bwd12_table.npz")
    meta = json.loads(z["meta"].item())
    gs, gc = z["states"], z["costs"].astype(np.int64)
    gmax = meta["gmax_completed"]

    # column content tables
    Lmod = 0b111
    for _ in range(2):
        Lmod = AL.pmul(Lmod, Lmod)
    Bmod = 0b11
    for _ in range(2):
        Bmod = AL.pmul(Bmod, Bmod)
    evL = np.array([AL.pmod(c, Lmod) for c in range(1 << P)],
                   dtype=np.uint8)
    evB = np.array([AL.pmod(c, Bmod) for c in range(1 << P)],
                   dtype=np.uint8)

    out = {}
    m = np.uint64((1 << P) - 1)
    for tag in ("T", "V"):
        for lvl in range(19, 23):
            f = DATA / f"s3_corridor12_{tag}_frontier_L{lvl}.npz"
            if not f.exists():
                continue
            sts = np.load(f)["states"]
            pcs = np.bitwise_count(sts).astype(np.int64)
            bq = CR.canon5(CR.phi(sts, P), P)
            idx = np.searchsorted(gs, bq)
            idx[idx >= gs.size] = gs.size - 1
            found = gs[idx] == bq
            G = np.where(found, gc[idx], -1)
            # per-column contents
            nsl = np.zeros(sts.size, dtype=np.int64)   # omega-only
            njt = np.zeros(sts.size, dtype=np.int64)   # joint
            nbr = np.zeros(sts.size, dtype=np.int64)   # barren-only
            for k in range(5):
                col = ((sts >> np.uint64(P * k)) & m).astype(np.int64)
                lo = evL[col] != 0
                bo = evB[col] != 0
                nsl += lo & ~bo
                njt += lo & bo
                nbr += ~lo & bo
            rec = dict(
                n=int(sts.size),
                pcs_hist={int(k): int(v) for k, v in
                          zip(*np.unique(pcs, return_counts=True))},
                G_hist={int(k): int(v) for k, v in
                        zip(*np.unique(G[found],
                                       return_counts=True))},
                n_absent=int((~found).sum()),
                slots_omega_only={int(k): int(v) for k, v in
                                  zip(*np.unique(
                                      nsl, return_counts=True))},
                slots_joint={int(k): int(v) for k, v in
                             zip(*np.unique(njt,
                                            return_counts=True))},
                cols_barren_only={int(k): int(v) for k, v in
                                  zip(*np.unique(
                                      nbr, return_counts=True))},
                n_omega=int(np.unique(_proj(sts, evL)).size),
                n_barren=int(np.unique(_proj(sts, evB)).size),
            )
            out[f"{tag}_L{lvl}"] = rec
            med = int(np.median(pcs))
            print(f"{tag} L{lvl}: n={rec['n']}  known-closers "
                  f"{rec['n'] - rec['n_absent']}  heavy-window "
                  f"unknowns {rec['n_absent']}  pcs median {med} "
                  f"(range {min(rec['pcs_hist'])}.."
                  f"{max(rec['pcs_hist'])})  "
                  f"joint-cols {rec['slots_joint']}", flush=True)
    out["gmax"] = gmax
    (DATA / "s3c_corner_census.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s3c_corner_census.json'}", flush=True)


def _proj(states, tab):
    m = np.uint64((1 << P) - 1)
    outv = np.zeros_like(states)
    for k in range(5):
        f = ((states >> np.uint64(P * k)) & m).astype(np.int64)
        outv |= tab[f].astype(np.uint64) << np.uint64(8 * k)
    return outv


if __name__ == "__main__":
    main()
