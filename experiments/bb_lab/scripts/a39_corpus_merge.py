#!/usr/bin/env python3
"""a39_corpus_merge.py — merge the A39 descent-theory-test outcomes into
`data/bb_instances.duckdb`.

Sources (all under experiments/bb_lab/data/, produced by the A39 line):
  * descent_theory_test/phase2_results.jsonl   — Stage-A closures:
      COUNTEREXAMPLE records (weight-complete transport => exact d at
      certificate tier) and CERTIFIED_FLOOR records (d_lb).
  * descent_theory_test/phase3_joins.jsonl     — cross-stage sandwich
      joins (Stage-A floor + re-verified witness => exact d).
  * order144_sweep/results.jsonl               — the 58 sampled order-144
      codes (row INSERTs; their own solver-/certificate-lane outcomes).
  * descent_theory_test/cohort.jsonl(+_batch2) — metadata for rows not in
      the corpus (parity-strip samples) and instance_id -> stratum map.

Rules (see PHASE4_SCORECARD.md; trust discipline is preserved verbatim):
  * existing d_exact is NEVER overwritten. Equal values count as
    agreement; different values ABORT the whole merge (nothing written).
  * exact values land in d_exact + d_method + cert_path; d_lb/d_ub are
    left as-is on exact rows (tandem-campaign convention).
  * floor-only outcomes raise d_lb (never lower), with cert_path set.
  * solver-lane sweep values keep the sweep's own solver d_method string;
    certificate-lane values get 'descent-cert@a39*' methods so the trust
    tier is readable from the DB.
  * new codes are INSERTed with bb_samp_-style code ids (sampled
    provenance, as in A18).

Default is a dry run. `--apply` first copies the DB to
`bb_instances.a39bak.duckdb` alongside, then writes in one transaction.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent          # experiments/bb_lab/scripts
LAB = HERE.parent                                # experiments/bb_lab
DTT = LAB / "data" / "descent_theory_test"
SWEEP = LAB / "data" / "order144_sweep"
DB = LAB / "data" / "bb_instances.duckdb"

CERT_PATH_DTT = "experiments/bb_lab/data/descent_theory_test/PHASE4_SCORECARD.md"
CERT_PATH_SWEEP = "experiments/bb_lab/data/order144_sweep/REPORT.md"

M_EXACT = "descent-cert@a39"
M_JOIN = "descent-cert@a39+witness-join"
M_SWEEP_CERT = "descent-cert@a39-sweep"


def jload(path: Path):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def collect():
    """Per-instance merged outcome dict:
    {exact, method, floor, cert_path, sources, insert_meta?}"""
    out: dict[str, dict] = {}

    def slot(iid):
        return out.setdefault(iid, {"exact": None, "method": None,
                                    "floor": 0, "cert_path": CERT_PATH_DTT,
                                    "sources": []})

    def put_exact(iid, d, method, src):
        s = slot(iid)
        if s["exact"] is not None and s["exact"] != d:
            raise SystemExit(f"INTERNAL CONFLICT {iid}: {s['exact']} vs {d} ({src})")
        # first writer wins on method; sources are recorded for the report
        if s["exact"] is None:
            s["exact"], s["method"] = d, method
        s["sources"].append(src)

    # Stage-A closures
    for r in jload(DTT / "phase2_results.jsonl"):
        iid = r.get("instance_id")
        e = r.get("eval") or {}
        if not iid or r.get("kind") != "closure":
            continue
        if e.get("outcome") == "CERTIFIED_FLOOR" and e.get("floor"):
            s = slot(iid)
            s["floor"] = max(s["floor"], int(e["floor"]))
            s["sources"].append("phase2:floor")
        if e.get("outcome") == "COUNTEREXAMPLE":
            w = (r.get("counterexample") or {}).get("weight")
            if w:
                put_exact(iid, int(w), M_EXACT, "phase2:cex")

    # sandwich joins (floor leg + re-verified witness leg)
    for r in jload(DTT / "phase3_joins.jsonl"):
        if not r.get("witness_verified"):
            continue
        put_exact(r["instance_id"], int(r["d_exact"]), M_JOIN, "phase3:join")

    # order-144 sweep rows: outcomes + insert metadata
    sweep_meta = {}
    for r in jload(SWEEP / "results.jsonl"):
        iid = r["instance_id"]
        sweep_meta[iid] = r
        s = slot(iid)
        s["cert_path"] = CERT_PATH_SWEEP
        if r.get("floor"):
            s["floor"] = max(s["floor"], int(r["floor"]))
            s["sources"].append("sweep:floor")
        if r.get("d_exact"):
            tier = (r.get("trust_tier") or "").lower()
            meth = M_SWEEP_CERT if "cert" in tier else (r.get("d_method") or "sat")
            put_exact(iid, int(r["d_exact"]), meth, f"sweep:{tier or 'solver'}")

    # cohort metadata (for inserts of rows the corpus has never seen)
    cohort_meta = {}
    for f in ("cohort.jsonl", "cohort_batch2.jsonl"):
        for r in jload(DTT / f):
            cohort_meta[r["instance_id"]] = r

    # sanity: floors never exceed exacts from the same line
    for iid, s in out.items():
        if s["exact"] is not None and s["floor"] > s["exact"]:
            raise SystemExit(f"INTERNAL CONFLICT {iid}: floor {s['floor']} > exact {s['exact']}")
    return out, sweep_meta, cohort_meta


def poly_weight(p: str) -> int:
    return len([t for t in p.split("+") if t.strip()])


def insert_row(iid, s, sweep_meta, cohort_meta, now):
    """Build the INSERT tuple for a row absent from the corpus."""
    if iid in sweep_meta:
        m = sweep_meta[iid]
        ell, mm = m.get("orders") or (None, None)
        group, n, k = m["group"], m["n"], m["k"]
        A, B = m["A"], m["B"]
        code_id = m.get("code_id") or f"bb_samp_{group}_{iid[:8]}"
        orbit = m.get("orbit_size")
        d_ub = m.get("d_ub")
    elif iid in cohort_meta:
        m = cohort_meta[iid]
        ell, mm = m.get("ell"), m.get("m")
        group, n, k = m["group_struct"], m["n"], m["k"]
        A, B = m["A_poly"], m["B_poly"]
        code_id = m.get("code_id") or f"bb_samp_{group}_{iid[:8]}"
        orbit = None
        d_ub = m.get("d_ub")
    else:
        return None
    return (iid, code_id, group, ell, mm, n, k, A, B,
            poly_weight(A), poly_weight(B), orbit,
            s["floor"] or None, d_ub, s["exact"], s["method"] if s["exact"] else None,
            s["cert_path"], now, now)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write (default: dry run)")
    ap.add_argument("--db", default=str(DB))
    args = ap.parse_args()

    import duckdb

    merged, sweep_meta, cohort_meta = collect()
    con = duckdb.connect(args.db, read_only=True)
    now = datetime.now(timezone.utc)

    plan = {"agree": [], "update_exact": [], "update_floor": [], "insert": [],
            "skip_nothing_new": []}
    conflicts = []

    for iid, s in sorted(merged.items()):
        row = con.execute(
            "select d_lb, d_ub, d_exact, d_method from bb_instances where instance_id = ?",
            [iid]).fetchone()
        if row is None:
            if s["exact"] is None and not s["floor"]:
                plan["skip_nothing_new"].append(iid)
            else:
                t = insert_row(iid, s, sweep_meta, cohort_meta, now)
                if t is None:
                    conflicts.append((iid, "no metadata for insert"))
                else:
                    plan["insert"].append(t)
            continue
        db_lb, db_ub, db_ex, db_meth = row
        if s["exact"] is not None:
            if db_ex is not None:
                if db_ex == s["exact"]:
                    plan["agree"].append((iid, db_ex, db_meth))
                else:
                    conflicts.append((iid, f"db d_exact={db_ex} ({db_meth}) vs a39 {s['exact']}"))
                continue
            if db_ub is not None and s["exact"] > db_ub:
                conflicts.append((iid, f"a39 exact {s['exact']} > db d_ub {db_ub}"))
                continue
            plan["update_exact"].append((iid, s["exact"], s["method"], s["cert_path"]))
        elif s["floor"]:
            if db_ex is not None:
                if s["floor"] > db_ex:
                    conflicts.append((iid, f"a39 floor {s['floor']} > db d_exact {db_ex}"))
                else:
                    plan["agree"].append((iid, db_ex, db_meth))
                continue
            if db_ub is not None and s["floor"] > db_ub:
                conflicts.append((iid, f"a39 floor {s['floor']} > db d_ub {db_ub}"))
                continue
            if db_lb is None or s["floor"] > db_lb:
                plan["update_floor"].append((iid, s["floor"], s["cert_path"]))
            else:
                plan["skip_nothing_new"].append(iid)
        else:
            plan["skip_nothing_new"].append(iid)
    con.close()

    print(f"plan: {len(plan['update_exact'])} exact updates, "
          f"{len(plan['insert'])} inserts, {len(plan['update_floor'])} floor raises, "
          f"{len(plan['agree'])} agreements with existing values, "
          f"{len(plan['skip_nothing_new'])} no-ops")
    meths = {}
    for _, _, m, _ in plan["update_exact"]:
        meths[m] = meths.get(m, 0) + 1
    for t in plan["insert"]:
        if t[15]:
            meths[t[15]] = meths.get(t[15], 0) + 1
    print("methods:", meths)
    if conflicts:
        print("\nCONFLICTS — NOTHING WRITTEN:")
        for iid, why in conflicts:
            print(f"  {iid}: {why}")
        sys.exit(2)

    if not args.apply:
        print("\ndry run only; re-run with --apply to write.")
        return

    bak = Path(args.db).with_name("bb_instances.a39bak.duckdb")
    shutil.copy2(args.db, bak)
    print(f"backup written: {bak}")

    con = duckdb.connect(args.db)
    con.execute("begin")
    for iid, d, meth, cp in plan["update_exact"]:
        con.execute("update bb_instances set d_exact=?, d_method=?, cert_path=?, "
                    "updated_at=? where instance_id=?", [d, meth, cp, now, iid])
    for iid, f, cp in plan["update_floor"]:
        con.execute("update bb_instances set d_lb=?, cert_path=?, updated_at=? "
                    "where instance_id=?", [f, cp, now, iid])
    for t in plan["insert"]:
        con.execute("insert into bb_instances (instance_id, code_id, group_struct, "
                    "ell, m, n, k, A_poly, B_poly, A_weight, B_weight, orbit_size, "
                    "d_lb, d_ub, d_exact, d_method, cert_path, inserted_at, updated_at) "
                    "values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", list(t))
    con.execute("commit")
    chk = con.execute("select d_method, count(*) from bb_instances "
                      "where d_method like 'descent-cert@a39%' group by 1").fetchall()
    total = con.execute("select count(*) from bb_instances").fetchone()
    con.close()
    print("post-apply certificate-tier rows:", chk)
    print("post-apply total rows:", total)


if __name__ == "__main__":
    main()
