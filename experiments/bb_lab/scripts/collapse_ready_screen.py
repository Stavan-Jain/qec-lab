"""Screen corpus d ≥ 14 hosts for collapse-ready doubling lifts.

For each base [[n, k, d̄]] and each axis lift (2ℓ, m) / (ℓ, 2m), in
falsify-first order:

  gate 1  k preservation (A12 (R)): k̃ = k by rank — else K-JUMP.
  gate 2  dangerous quick-kill: light λ=0 pushforwards (single check
          rows and row pairs, weight ≤ ~12) — budgeted SAT probe for a
          nontrivial completion of cost < 2d̄; finding one FAILS the
          doubling (and, re-verified, gives a d_cover upper bound).
  gate 3  λ ≠ 0 strata: |w| ≥ d̄ by the base-distance theorem, so
          sample λ's band (d̄, d̄+2) — sizes (sparsity) and per-w
          budgeted masked-K (min |b ∖ w| over the twist coset, the
          quantity that must reach (2d̄ − |w|)/2 for the window to
          certify 2d̄).

Verdicts: K-JUMP / DANGEROUS-KILL(cost) / PROMISING(sparse, masked
climbing) / GRAY. PROMISING candidates are the inputs to a full
twist_floor_sweep at target 2d̄.

Usage: uv run python scripts/collapse_ready_screen.py --out DIR [--jobs 2]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

import duckdb
import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

import pycryptosat  # noqa: E402
from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.descent_sat import compute_descent  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.linalg import nullspace_f2, quotient_complement_basis  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.sat_distance import find_logical_z  # noqa: E402
from bb_lab.shard_distance import (  # noqa: E402
    _bits_to_int,
    _int_to_bits,
    compute_class_action,
)
from bb_lab.twist_floors import compute_twist, enum_band  # noqa: E402

BASES = [  # (instance_id, tag) — k=6 d=16 first, then strong k=4
    ("1e6e791ebce284c9", "wm154_6_16"),
    ("344b6f450bdd03e9", "l168_6_16"),
    ("0f3ad144649e191b", "l168_4_16a"),
    ("5e31c1656d628e71", "l168_4_16b"),
]


def masked_K_sat(checks_cover, td, w, need, budget=150_000,
                 require_nontrivial=False, L_Z=None):
    """Budgeted certified floor on min |b∖w| over the twist coset
    (optionally excluding stabilizer completions). Ascending j:
    UNSAT at j ⟹ ≥ j+1; SAT returns (exact j, witness b)."""
    nb = td.Hb_Z.shape[1]
    supp = set(np.flatnonzero(w).tolist())
    outside = [i for i in range(nb) if i not in supp]
    syn = (td.T @ w) % 2
    for j in range(0, need):
        pool = IDPool()
        bv = [pool.id() for _ in range(nb)]
        s = pycryptosat.Solver(confl_limit=budget)
        for r, row in enumerate(td.Hb_Z):
            idx = np.flatnonzero(row)
            if idx.size:
                s.add_xor_clause([bv[i] for i in idx], bool(syn[r]))
            elif syn[r]:
                return need, None      # inconsistent: coset empty
        if require_nontrivial:
            # v = p*(b) + L0(w) nontrivial: μ(v) = Pμ·b + μ(L0 w) ≠ 0
            # encoded on the cover pairing via L_Z of the cover.
            pass  # handled by caller re-verification instead
        card = CardEnc.atmost(
            lits=[bv[i] for i in outside], bound=j, vpool=pool,
            encoding=EncType.seqcounter,
        )
        for cl in card.clauses:
            s.add_clause(cl)
        sat, model = s.solve()
        if sat is True:
            b = np.array([1 if model[v] else 0 for v in bv], np.uint8)
            return j, b
        if sat is not False:
            return j, None             # budget: certified ≥ j
    return need, None


_CTX = {}


def screen_one(args):
    iid, tag, axis = args
    con = duckdb.connect(DB, read_only=True)
    ell, m, n, k, d, a_s, b_s = con.execute(
        "SELECT ell, m, n, k, d_exact, A_poly, B_poly FROM bb_instances "
        "WHERE instance_id = ?", [iid],
    ).fetchone()
    con.close()
    ell, m, k, d = int(ell), int(m), int(k), int(d)
    target = 2 * d
    lift = (2 * ell, m) if axis == 0 else (ell, 2 * m)
    out = {
        "base": tag, "axis": axis, "lift_group": lift,
        "target": target, "n_cover": 4 * ell * m,
    }
    t0 = time.perf_counter()
    G = ZmZn(*lift)
    A = Poly.from_string(a_s, G)
    B = Poly.from_string(b_s, G)
    cover = bb_check_matrices(A, B)

    # gate 1: k preservation
    kt = find_logical_z(cover).shape[0]
    out["k_cover"] = int(kt)
    if kt != k:
        out["verdict"] = f"K-JUMP({kt})"
        out["seconds"] = round(time.perf_counter() - t0, 1)
        return out

    action = compute_class_action(cover)
    sigma = (ell, 0) if axis == 0 else (0, m)
    dd = compute_descent(cover, action, sigma, A=A, B=B)
    td = compute_twist(cover, dd)
    L_Z = action.L_Z

    # gate 2: dangerous quick-kill — light boundaries as w
    Hb_X = dd.base_checks.H_X
    kills = []
    probes = 0
    for i in range(0, Hb_X.shape[0], 7):
        cands = [Hb_X[i] % 2]
        if i + 1 < Hb_X.shape[0]:
            cands.append((Hb_X[i] + Hb_X[i + 1]) % 2)
        for w in cands:
            w = w.astype(np.uint8)
            wt = int(w.sum())
            if wt == 0 or wt >= target:
                continue
            probes += 1
            need = (target - wt + 1) // 2
            j, b = masked_K_sat(cover, td, w, need)
            if b is not None:
                v = ((dd.P_lift @ b) + (td.L0 @ w)) % 2
                cost = int(v.sum())
                nontriv = bool(((L_Z @ v) % 2).any())
                cyc = not ((cover.H_Z @ v) % 2).any()
                if cyc and nontriv and cost < target:
                    kills.append((wt, cost))
    out["dangerous_probes"] = probes
    if kills:
        out["verdict"] = f"DANGEROUS-KILL(d_cover<= {min(c for _, c in kills)})"
        out["kills"] = kills[:5]
        out["seconds"] = round(time.perf_counter() - t0, 1)
        return out

    # gate 3: λ ≠ 0 strata sparsity + masked climb
    lam_img = sorted(
        {
            _bits_to_int((_int_to_bits(c, k) @ dd.Lam) % 2)
            for c in range(1, 1 << k)
        }
        - {0}
    )
    samples = lam_img[:: max(1, len(lam_img) // 4)][:4]
    strata = {}
    climbs = []
    for lam in samples:
        ws, complete = enum_band(
            td, dd.Lb_Z, dd.kb, lam, d - 1, d + 2,
            limit=60, wall_budget=180,
        )
        strata[lam] = (len(ws), complete)
        for w in ws[:3]:
            wt = int(w.sum())
            need = min((target - wt + 1) // 2, 4)
            j, b = masked_K_sat(cover, td, w, need, budget=100_000)
            climbs.append((wt, j, b is not None))
    out["strata_sample"] = {str(kk): vv for kk, vv in strata.items()}
    out["masked_climbs"] = climbs
    sparse = all(nn <= 40 for nn, _ in strata.values())
    healthy = all(
        (not found) or (wt + 2 * j >= target)
        for wt, j, found in climbs
    )
    out["verdict"] = (
        "PROMISING" if (sparse and healthy and climbs)
        else ("GRAY" if healthy else "WINDOW-KILL")
    )
    out["seconds"] = round(time.perf_counter() - t0, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    items = [
        (iid, tag, axis) for iid, tag in BASES for axis in (0, 1)
    ]
    print(f"{len(items)} lift candidates", flush=True)
    results = []
    with mp.get_context("fork").Pool(args.jobs) as pool:
        for res in pool.imap_unordered(screen_one, items):
            results.append(res)
            print(json.dumps(res), flush=True)
    (outdir / "collapse_screen.json").write_text(json.dumps(results))
    for r in sorted(results, key=lambda r: r.get("verdict", "")):
        print(
            f"{r['base']:<12} axis={r['axis']} → {r['lift_group']} "
            f"[[{r['n_cover']},{r.get('k_cover', '?')},≤{r['target']}]]: "
            f"{r['verdict']}"
        )


if __name__ == "__main__":
    main()
