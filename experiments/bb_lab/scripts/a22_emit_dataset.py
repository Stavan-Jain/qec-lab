"""A22: emit the canonical delta-side classification dataset.

Writes data/a22/alpha_classes_full.json: per alpha-class everything a Lean
emitter (or the writeup) needs — canonical alpha, active set, per-site
values/types/costs, optimal h, m* (the <=3-site preimage), the h-flip
children with their |b| weights, and the match into the 113-class file.

Usage: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a22_emit_dataset.py
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a22_common import (
    A22_DATA, ATIL, AT_HAT, ETA0, LAB_ROOT, PREIMAGE, SIDX, SITES, WTAB,
    apply_theta, canon_alpha, compute_alpha_classes, config_cost, active_set,
    conv_site_gf16, cost_and_h, dft_vec, ginv, gmul, idft_vec, jofab,
    site_add, site_sub, site_type,
)

# ---- reconstruct the file classes for matching
DATA = LAB_ROOT / "data" / "a17" / "f2a6_light_classes.jsonl"
classes = [json.loads(l) for l in DATA.read_text().splitlines() if l.strip()]
classes = [c for c in classes if "b_weight" in c]
assert len(classes) == 113


def canon_uv(u: np.ndarray, v: np.ndarray) -> bytes:
    best = None
    for di in range(5):
        for dj in range(15):
            uu = np.roll(np.roll(u, di, 0), dj, 1)
            vv = np.roll(np.roll(v, di, 0), dj, 1)
            cand = uu.tobytes() + vv.tobytes()
            if best is None or cand < best:
                best = cand
    return best


file_canon: dict[bytes, int] = {}
for ci, c in enumerate(classes):
    u = np.zeros((5, 15), dtype=np.uint8)
    v = np.zeros((5, 15), dtype=np.uint8)
    for blk, i, j in c["b_support"]:
        (u if blk == 0 else v)[i, j] = 1
    file_canon[canon_uv(u, v)] = ci


def fiber_glue(eps: np.ndarray, dlt: np.ndarray) -> np.ndarray:
    u = np.zeros((5, 15), dtype=np.uint8)
    for k, (i, b) in enumerate(SITES):
        p = PREIMAGE[(int(eps[k]), int(dlt[k]))]
        for a in range(5):
            if (p >> a) & 1:
                u[i, jofab(a, b)] = 1
    return u


def build_uv(alpha: np.ndarray, h: np.ndarray):
    beta = apply_theta(alpha)
    ue, ud = h, alpha
    ve = np.array([h[SIDX[site_sub(s, (1, 0))]] for s in SITES])
    vd = np.zeros(15, dtype=np.int64)
    for kg, g in enumerate(SITES):
        vd[SIDX[site_add(g, (1, 0))]] = beta[kg]
    return fiber_glue(ue, ud), fiber_glue(ve, vd)


KH = np.zeros((5, 3), dtype=np.int64)
KH[ETA0] = 1
KER = idft_vec(KH)

reps = compute_alpha_classes()
assert len(reps) == 94

records = []
matched: set[int] = set()

# the pure-h class
h1 = np.zeros(15, dtype=np.int64)
h1[0] = 1
u, v = build_uv(np.zeros(15, dtype=np.int64), h1)
ci = file_canon[canon_uv(u, v)]
matched.add(ci)
records.append({"kind": "pure-h", "class_ids": [ci], "b_weight": 10})

for rep in reps:
    alpha = np.array(rep, dtype=np.int64)
    beta = apply_theta(alpha)
    S = active_set(alpha)
    # m* preimage
    ah = dft_vec(alpha)
    mh = np.zeros((5, 3), dtype=np.int64)
    for p in range(5):
        for q in range(3):
            if (p, q) != ETA0:
                mh[p, q] = gmul(int(ah[p, q]), ginv(int(AT_HAT[p, q])))
    m0 = idft_vec(mh)
    best, bs = None, 99
    for cc in range(16):
        m = m0.copy()
        if cc:
            for k in range(15):
                m[k] ^= gmul(cc, int(KER[k]))
        s = int(np.count_nonzero(m))
        if s < bs:
            best, bs = m, s
    assert bs <= 3
    assert np.array_equal(conv_site_gf16(ATIL, best), alpha)

    hstar = np.zeros(15, dtype=np.int64)
    pens, cost = [], 0
    site_rows = []
    for k in range(15):
        ck, hk, pk = cost_and_h(int(alpha[k]), int(beta[k]))
        cost += ck
        hstar[k] = hk
        pens.append(pk)
        if alpha[k] or beta[k]:
            site_rows.append({
                "site": list(SITES[k]), "alpha": int(alpha[k]),
                "beta": int(beta[k]),
                "types": site_type(int(alpha[k])) + site_type(int(beta[k])),
                "cost": ck, "h": hk, "flip_pen": pk})
    budget = 14 - cost
    children = []
    flippable = [k for k in range(15) if pens[k] <= budget]

    def rec(idx, rem, flips):
        hh = hstar.copy()
        for k in flips:
            hh[k] ^= 1
        u, v = build_uv(alpha, hh)
        w = int(u.sum() + v.sum())
        if 0 < w <= 14:
            ci = file_canon.get(canon_uv(u, v))
            assert ci is not None, "derived class not in file!"
            children.append({"flips": [list(SITES[k]) for k in flips],
                             "b_weight": w, "class_id": ci})
            matched.add(ci)
        for pos in range(idx, len(flippable)):
            k = flippable[pos]
            if pens[k] <= rem:
                rec(pos + 1, rem - pens[k], flips + [k])

    rec(0, budget, [])
    records.append({
        "kind": "alpha", "alpha": [int(x) for x in rep],
        "active_sites": [list(SITES[k]) for k in S], "n_active": len(S),
        "cost": cost, "sites": site_rows,
        "m_star": {str(SITES[k]): int(best[k]) for k in range(15) if best[k]},
        "m_star_nsites": bs,
        "children": children,
    })

assert len(matched) == 113, f"only matched {len(matched)} of 113"
out = A22_DATA / "alpha_classes_full.json"
out.write_text(json.dumps(records, indent=1))
print(f"wrote {out} ({len(records)} records: 1 pure-h + 94 alpha classes; "
      f"all 113 file classes matched)")

# summary table
by = Counter()
for r in records:
    if r["kind"] == "alpha":
        kids = Counter(ch["b_weight"] for ch in r["children"])
        by[(r["n_active"], r["cost"], r["m_star_nsites"],
            tuple(sorted(kids.items())))] += 1
print("\n(|S|, cost, |m*|, children-weight-multiset) -> #classes")
for k, cnt in sorted(by.items()):
    print(f"   {k}: {cnt}")
