"""Order-144 sweep batch driver.

Per-code protocol under a HARD 120 s budget (cache misses the code
triggers are charged to it; cache hits are free):

  1. quotient ladder (order 72/36/18 by exponent reduction) with a
     shared base cache (sweep DB `base_cache`; corpus consulted
     read-only first, both axis orientations),
  2. L1-sampling d_ub at n = 288 (verified witness),
  3. doubling-certificate lane (bb_lab.doubling_certify) when a
     literal-lift (R) candidate exists and budget remains,
  4. SAT exact ladder with watchdog; on timeout the contiguous UNSAT
     prefix is recorded as an honest floor (solver-proved refutations).

DISCIPLINE: witness weights are upper bounds only; every exact claim
carries method + trust tier; UNKNOWN stays UNKNOWN.

Results append to results.jsonl as codes finish (crash-safe).

  caffeinate -ims uv run python data/order144_sweep/driver.py [--limit N]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import duckdb
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import sweeplib  # noqa: E402
from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.l1_sampling import l1_distance_ub, verify_witness_in_nontrivial_coset  # noqa: E402
from bb_lab.doubling_certify import detect  # noqa: E402
from bb_lab.store import canonical_hash  # noqa: E402

MAIN_DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"
SWEEP_DB = HERE / "sweep.duckdb"
RESULTS = HERE / "results.jsonl"
BUDGET = 120.0
GROUP_ORDER = ["Z12xZ12", "Z18xZ8", "Z24xZ6", "Z36xZ4", "Z16xZ9", "Z48xZ3"]

SAT_METHOD = "sat-cms-ladder (pycryptosat via bb_lab.sat_distance)"
UB_METHOD = "l1-sampling@100k (verified witness)"
CERT_METHOD = "doubling-certify (bb_lab.doubling_certify A30 front-end)"


def free_mb() -> float:
    return shutil.disk_usage("/Users/stavanjain").free / 1e6


# ---------------------------------------------------------------- corpus map

QUOT_LABELS: set[str] = set()
for l0, m0 in [(12, 12), (18, 8), (24, 6), (36, 4), (16, 9), (48, 3)]:
    for q in sweeplib.quotient_chain((l0, m0)):
        QUOT_LABELS.add(AbelianGroup(q).label())
        QUOT_LABELS.add(AbelianGroup((q[1], q[0])).label())


def load_corpus_map() -> dict[str, dict]:
    con = duckdb.connect(MAIN_DB, read_only=True)
    ph = ",".join("?" for _ in QUOT_LABELS)
    rows = con.execute(
        f"SELECT instance_id, group_struct, code_id, n, k, d_ub, d_exact, d_method "
        f"FROM bb_instances WHERE group_struct IN ({ph})",
        list(QUOT_LABELS),
    ).fetchall()
    con.close()
    return {
        r[0]: {"group": r[1], "code_id": r[2], "n": r[3], "k": r[4],
               "d_ub": r[5], "d_exact": r[6], "d_method": r[7]}
        for r in rows
    }


CORPUS = load_corpus_map()


# ---------------------------------------------------------------- base cache


def init_cache(con) -> None:
    con.execute("""
        CREATE TABLE IF NOT EXISTS base_cache (
          cache_key TEXT PRIMARY KEY,
          label TEXT, ell INT, m INT, n INT, k INT,
          A_poly TEXT, B_poly TEXT,
          d_ub INT, d_exact INT, d_method TEXT, floor_w INT,
          source TEXT, wall_s DOUBLE, updated_at TIMESTAMP DEFAULT now()
        )""")


def cache_get(con, key: str):
    r = con.execute(
        "SELECT label, ell, m, n, k, A_poly, B_poly, d_ub, d_exact, d_method, "
        "floor_w, source, wall_s FROM base_cache WHERE cache_key = ?", [key]
    ).fetchone()
    if r is None:
        return None
    return dict(zip(["label", "ell", "m", "n", "k", "A_poly", "B_poly", "d_ub",
                     "d_exact", "d_method", "floor_w", "source", "wall_s"], r))


def cache_put(con, key: str, e: dict) -> None:
    con.execute(
        "INSERT OR REPLACE INTO base_cache (cache_key, label, ell, m, n, k, "
        "A_poly, B_poly, d_ub, d_exact, d_method, floor_w, source, wall_s, "
        "updated_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,now())",
        [key, e["label"], e["ell"], e["m"], e["n"], e["k"], e["A_poly"],
         e["B_poly"], e.get("d_ub"), e.get("d_exact"), e.get("d_method"),
         e.get("floor_w"), e.get("source"), e.get("wall_s")])


# ---------------------------------------------------------------- SAT lane


def run_sat(ell: int, m: int, A_str: str, B_str: str, timeout: float,
            wmax: int | None, tag: str) -> dict:
    """sat_worker under a hard watchdog. Returns
    {outcome: exact|timeout|error, d?, floor_w?, wall_s, rungs}."""
    prog = HERE / f"prog_{tag}_{os.getpid()}.txt"
    prog.unlink(missing_ok=True)
    cmd = [sys.executable, str(HERE / "sat_worker.py"), "--ell", str(ell),
           "--m", str(m), "--A", A_str, "--B", B_str, "--progress", str(prog)]
    if wmax is not None:
        cmd += ["--wmax", str(wmax)]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        wall = time.time() - t0
        out = {"wall_s": round(wall, 1)}
        if p.returncode == 0 and "DISTANCE" in p.stdout:
            out["outcome"] = "exact"
            out["d"] = int(p.stdout.split("DISTANCE")[1].split()[0])
        else:
            out["outcome"] = "error"
            out["error"] = (p.stderr.strip().splitlines() or ["?"])[-1][:200]
    except subprocess.TimeoutExpired:
        out = {"outcome": "timeout", "wall_s": round(time.time() - t0, 1)}
    rungs = []
    if prog.exists():
        for line in prog.read_text().splitlines():
            w, s, dt = line.split(",")
            rungs.append((int(w), s, float(dt)))
        prog.unlink()
    out["rungs"] = rungs
    # honest floor: contiguous UNSAT prefix starting at the ladder's wmin
    floor_w = 0
    for w, s, _ in rungs:
        if s == "UNSAT" and w == floor_w + 1:
            floor_w = w
        else:
            break
    if floor_w:
        out["floor_w"] = floor_w        # solver-proved: d >= floor_w + 1
    return out


# ---------------------------------------------------------------- ladder


def corpus_lookup(orders, A_supp, B_supp):
    """Look the quotient up in the main corpus, both axis orientations."""
    L, M = orders
    _, _, _, iid, _ = sweeplib.canonicalize(L, M, A_supp, B_supp)
    hit = CORPUS.get(iid)
    if hit:
        return iid, hit, "native"
    At = tuple((b, a) for a, b in A_supp)
    Bt = tuple((b, a) for a, b in B_supp)
    _, _, _, iid_t, _ = sweeplib.canonicalize(M, L, At, Bt)
    hit = CORPUS.get(iid_t)
    if hit:
        return iid_t, hit, "transposed"
    return iid, None, None


def ladder_phase(con, orders, A_supp, B_supp, rem) -> tuple[list[dict], float]:
    """Descend the canonical chain; fill the shared cache within budget."""
    t0 = time.time()
    entries = []
    cur_o, cur_A, cur_B = orders, A_supp, B_supp
    for q in sweeplib.quotient_chain(orders):
        cur_A, cur_B = sweeplib.quotient_code(cur_o, cur_A, cur_B, q)
        cur_o = q
        p, checks, Aq, Bq = sweeplib.code_params_of(q, cur_A, cur_B)
        label, cA, cB, key, _ = sweeplib.canonicalize(*q, cur_A, cur_B)
        ent = {"orders": list(q), "n": p.n, "k": p.k, "A": cA, "B": cB,
               "cache_key": key[:12], "wA": len(cur_A), "wB": len(cur_B)}
        if p.k < 2:
            ent["d_note"] = "k < 2: no distance" + (
                " (unit poly)" if 1 in (len(cur_A), len(cur_B)) else "")
            entries.append(ent)
            continue
        hit = cache_get(con, key)
        if hit is not None:
            ent.update({k2: hit[k2] for k2 in ("d_ub", "d_exact", "d_method",
                                               "floor_w", "source")})
            ent["cache"] = "hit"
            entries.append(ent)
            continue
        # corpus first (free)
        ciid, chit, orient = corpus_lookup(q, cur_A, cur_B)
        centry = {"label": label, "ell": q[0], "m": q[1], "n": p.n, "k": p.k,
                  "A_poly": cA, "B_poly": cB, "wall_s": 0.0}
        if chit is not None and chit["d_exact"] is not None:
            centry.update({"d_exact": chit["d_exact"], "d_ub": chit["d_ub"],
                           "d_method": chit["d_method"],
                           "source": f"corpus:{chit['code_id']}({orient})"})
            cache_put(con, key, centry)
            ent.update({k2: centry.get(k2) for k2 in ("d_ub", "d_exact",
                                                      "d_method", "source")})
            ent["cache"] = "corpus"
            entries.append(ent)
            continue
        # compute within remaining budget
        budget_left = rem()
        wall0 = time.time()
        if p.n <= 36:
            to = min(12.0, budget_left)
        elif p.n <= 72:
            to = min(25.0, budget_left)
        else:
            to = None  # n = 144: L1 first, SAT gated below
        if p.n <= 72 and to > 2:
            r = run_sat(q[0], q[1], cA, cB, to, None, f"lad{p.n}")
            if r["outcome"] == "exact":
                centry.update({"d_exact": r["d"], "d_method": SAT_METHOD,
                               "source": "sat"})
            else:
                centry.update({"floor_w": r.get("floor_w"),
                               "d_method": f"sat-{r['outcome']}@{to:.0f}s",
                               "source": "sat-partial"})
        elif p.n == 144:
            res = l1_distance_ub(checks, n_samples=60_000, seed=11)
            centry["d_ub"] = int(res.distance_ub)
            centry["d_method"] = "l1-sampling@60k (ub only)"
            centry["source"] = "l1"
            if centry["d_ub"] <= 8 and rem() >= 90:
                r = run_sat(q[0], q[1], cA, cB, min(45.0, rem() - 30),
                            centry["d_ub"], "lad144")
                if r["outcome"] == "exact":
                    centry.update({"d_exact": r["d"], "d_method": SAT_METHOD,
                                   "source": "sat"})
                elif r.get("floor_w"):
                    centry["floor_w"] = r["floor_w"]
                    centry["d_method"] += f" + sat-floor>={r['floor_w'] + 1}"
        else:
            centry["d_method"] = "skipped (budget)"
            centry["source"] = "none"
        centry["wall_s"] = round(time.time() - wall0, 1)
        cache_put(con, key, centry)
        ent.update({k2: centry.get(k2) for k2 in ("d_ub", "d_exact", "d_method",
                                                  "floor_w", "source")})
        ent["cache"] = "miss-computed"
        ent["wall_s"] = centry["wall_s"]
        entries.append(ent)
    return entries, time.time() - t0


# ---------------------------------------------------------------- certify


def certify_phase(row, rem) -> dict:
    orders = (row["ell"], row["m"])
    G = AbelianGroup(orders)
    A = Poly.from_string(row["A_poly"], G)
    B = Poly.from_string(row["B_poly"], G)
    cands = detect(G, A, B)
    info: dict = {"n_candidates": len(cands),
                  "n_R": sum(c.R_holds for c in cands)}
    usable = [c for c in cands if c.R_holds]
    if not usable:
        info["outcome"] = "skip: no (R) literal-lift candidate"
        return info
    # cheap scope pre-gate: base d_ub must plausibly be <= 15
    c0 = usable[0]
    Gb = AbelianGroup(c0.base_group)
    bchecks = bb_check_matrices(Poly.from_string(c0.base_A, Gb),
                                Poly.from_string(c0.base_B, Gb))
    bub = int(l1_distance_ub(bchecks, n_samples=40_000, seed=5).distance_ub)
    info["base_d_ub_l1"] = bub
    if bub > 15:
        info["outcome"] = f"skip: base d_ub {bub} > 15 (front-end scope)"
        return info
    if rem() < 40:
        info["outcome"] = "skip: budget"
        return info
    budget = min(rem() - 20, 90.0)
    wd = HERE / "certify_runs" / row["instance_id"][:12]
    t0 = time.time()
    try:
        p = subprocess.run(
            [sys.executable, str(HERE / "certify_worker.py"),
             "--ell", str(orders[0]), "--m", str(orders[1]),
             "--A", row["A_poly"], "--B", row["B_poly"],
             "--budget", str(budget), "--workdir", str(wd), "--threads", "4"],
            capture_output=True, text=True, timeout=budget + 25)
        line = [l for l in p.stdout.splitlines() if l.startswith("CERTIFY_JSON")]
        if line:
            v = json.loads(line[0][len("CERTIFY_JSON "):])
            info["verdict"] = {k2: v.get(k2) for k2 in
                               ("status", "reason", "distance", "base",
                                "d_base", "wall_s", "candidate_log")}
            info["outcome"] = v.get("status", "?")
        else:
            info["outcome"] = "error"
            info["error"] = (p.stderr.strip().splitlines() or ["?"])[-1][:300]
    except subprocess.TimeoutExpired:
        info["outcome"] = "timeout"
    info["wall_s"] = round(time.time() - t0, 1)
    return info


# ---------------------------------------------------------------- per code


def process_code(con, row: dict) -> dict:
    t_start = time.time()
    rem = lambda: BUDGET - (time.time() - t_start)
    orders = (row["ell"], row["m"])
    out: dict = {
        "code_id": row["code_id"], "instance_id": row["instance_id"],
        "group": row["group_struct"], "orders": list(orders),
        "n": row["n"], "k": row["k"], "A": row["A_poly"], "B": row["B_poly"],
        "orbit_size": row["orbit_size"],
        "provenance": "sampled (bb_samp_, non-exhaustive)",
        "corpus_overlap": row["instance_id"] in CORPUS,
        "lanes": {},
    }
    A_supp = sweeplib.support_from_string(orders, row["A_poly"])
    B_supp = sweeplib.support_from_string(orders, row["B_poly"])

    # 1. ladder
    ladder, t_lad = ladder_phase(con, orders, A_supp, B_supp, rem)
    out["ladder"] = ladder
    out["lanes"]["ladder_wall_s"] = round(t_lad, 1)

    # 2. L1 d_ub at n = 288
    t0 = time.time()
    G = AbelianGroup(orders)
    checks = bb_check_matrices(Poly.from_string(row["A_poly"], G),
                               Poly.from_string(row["B_poly"], G))
    res = l1_distance_ub(checks, n_samples=100_000, seed=3)
    wv = bool(verify_witness_in_nontrivial_coset(checks, res.witness))
    d_ub = int(res.distance_ub)
    out["d_ub"] = d_ub
    out["d_ub_method"] = UB_METHOD
    out["d_ub_witness_verified"] = wv
    out["lanes"]["l1_wall_s"] = round(time.time() - t0, 1)

    d_exact = None
    d_method = None
    trust = None
    floor = None
    floor_method = None

    # 3. certificate lane
    cert = certify_phase(row, rem)
    out["lanes"]["certify"] = cert
    v = cert.get("verdict") or {}
    if cert.get("outcome") == "CERTIFIED" and v.get("distance"):
        d_exact = int(v["distance"]["value"])
        d_method = CERT_METHOD + f" = 2*d_base({v['distance'].get('d_base')})"
        trust = "certificate (counting-invariant enumeration + doubling theorem; not kernel-checked)"
        if d_exact > d_ub:
            out["ANOMALY"] = f"certified d {d_exact} > witness-backed d_ub {d_ub}"
    elif cert.get("outcome") == "FLOOR-ONLY" and v.get("distance"):
        floor = int(v["distance"]["floor"])
        floor_method = CERT_METHOD + " safe-floor (certificate tier)"

    # 4. SAT lane
    if d_exact is None and rem() >= 8:
        sat = run_sat(orders[0], orders[1], row["A_poly"], row["B_poly"],
                      rem(), d_ub, "cover")
        out["lanes"]["sat"] = sat
        if sat["outcome"] == "exact":
            d_exact = sat["d"]
            d_method = SAT_METHOD
            trust = "solver-exact (CMS SAT ladder, witness + UNSAT rounds)"
        elif sat.get("floor_w"):
            fw = sat["floor_w"] + 1
            if floor is None or fw > floor:
                floor = fw
                floor_method = f"sat-unsat rounds complete through w = {sat['floor_w']} (solver-proved refutations)"
            if fw == d_ub and wv:
                d_exact = d_ub
                d_method = (f"sat-unsat floor d >= {fw} + verified L1 witness "
                            f"at weight {d_ub}")
                trust = "solver-exact (UNSAT side) + explicit verified witness (SAT side)"
    elif d_exact is None:
        out["lanes"]["sat"] = {"outcome": "skip: budget exhausted"}

    out["d_exact"] = d_exact
    out["d_method"] = d_method
    out["trust_tier"] = trust
    out["floor"] = floor
    out["floor_method"] = floor_method
    if d_exact is None:
        out["d_status"] = f"UNKNOWN in [{floor or 1}, {d_ub}]"
    out["wall_s"] = round(time.time() - t_start, 1)
    out["ts"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    con.execute(
        "UPDATE bb_instances SET d_ub = ?, d_exact = ?, d_method = ?, "
        "d_lb = ?, updated_at = now() WHERE instance_id = ?",
        [d_ub, d_exact,
         d_method or (f"bounded-only; floor {floor}" if floor else "bounded-only"),
         floor, row["instance_id"]])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--only-group", default=None)
    args = ap.parse_args()

    con = duckdb.connect(str(SWEEP_DB))
    init_cache(con)
    done_ids = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            try:
                done_ids.add(json.loads(line)["instance_id"])
            except Exception:
                pass
    order_case = "CASE group_struct " + " ".join(
        f"WHEN '{g}' THEN {i}" for i, g in enumerate(GROUP_ORDER)) + " ELSE 99 END"
    rows = con.execute(
        f"SELECT instance_id, code_id, group_struct, ell, m, n, k, A_poly, "
        f"B_poly, orbit_size FROM bb_instances "
        f"ORDER BY {order_case}, k DESC, instance_id").fetchall()
    cols = ["instance_id", "code_id", "group_struct", "ell", "m", "n", "k",
            "A_poly", "B_poly", "orbit_size"]
    todo = [dict(zip(cols, r)) for r in rows if r[0] not in done_ids]
    if args.only_group:
        todo = [r for r in todo if r["group_struct"] == args.only_group]
    if args.limit:
        todo = todo[: args.limit]
    print(f"{len(todo)} codes to process ({len(done_ids)} already done)", flush=True)

    t_batch = time.time()
    for i, row in enumerate(todo):
        if free_mb() < 250:
            print("DISK GUARD: < 250 MB free, stopping", flush=True)
            break
        r = process_code(con, row)
        with open(RESULTS, "a") as f:
            f.write(json.dumps(r) + "\n")
        tag = (f"d={r['d_exact']}" if r["d_exact"] else
               f"{r.get('d_status', '?')}")
        print(f"[{i + 1}/{len(todo)}] {row['group_struct']} k={row['k']} "
              f"{tag} ({r['wall_s']}s) [batch {time.time() - t_batch:.0f}s]",
              flush=True)
    con.close()
    print(f"BATCH DONE in {time.time() - t_batch:.0f}s", flush=True)


if __name__ == "__main__":
    main()
