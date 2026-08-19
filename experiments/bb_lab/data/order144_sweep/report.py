"""Render REPORT.md from results.jsonl + sweep.duckdb."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import duckdb

HERE = Path(__file__).resolve().parent


def main() -> None:
    rows = [json.loads(l) for l in (HERE / "results.jsonl").read_text().splitlines()]
    con = duckdb.connect(str(HERE / "sweep.duckdb"), read_only=True)

    by_group = defaultdict(list)
    for r in rows:
        by_group[r["group"]].append(r)

    def tier(r) -> str:
        if r["d_exact"] is None:
            return "bounded-only"
        t = r.get("trust_tier") or ""
        if "certificate" in t:
            return "certificate"
        return "solver-exact"

    tiers = Counter(tier(r) for r in rows)
    d_hist = Counter(r["d_exact"] for r in rows if r["d_exact"] is not None)
    ub_hist = Counter(r["d_ub"] for r in rows if r["d_exact"] is None)
    total_wall = sum(r["wall_s"] for r in rows)

    L = []
    L.append("# Order-144 BB breadth sweep — n = 288 (2026-08-17)\n")
    L.append("Weight-3 x weight-3 BB codes over abelian groups of order 144, "
             "sampled enumeration (`bb_samp_` provenance, a18 conventions), "
             "canonical Aut x translation x swap dedup, k >= 2, HARD 120 s "
             "per-code compute cap. This is a data errand (no approach "
             "number); all artifacts live in this directory.\n")

    L.append("## Headline\n")
    L.append(f"- codes processed: **{len(rows)}** "
             f"({', '.join(f'{g}: {len(v)}' for g, v in by_group.items())})")
    L.append(f"- exact d: **{tiers['solver-exact'] + tiers['certificate']}** "
             f"(solver-exact {tiers['solver-exact']}, certificate "
             f"{tiers['certificate']}); bounded-only: {tiers['bounded-only']}")
    L.append(f"- exact-d histogram: "
             + ", ".join(f"d={d}: {c}" for d, c in sorted(d_hist.items())))
    if ub_hist:
        L.append(f"- bounded-only rows by d_ub: "
                 + ", ".join(f"ub={d}: {c}" for d, c in sorted(ub_hist.items())))
    L.append(f"- total per-code wall: {total_wall / 60:.1f} min "
             f"(mean {total_wall / len(rows):.0f} s/code)\n")

    L.append("## Counts by group x outcome\n")
    L.append("| group | codes | solver-exact | certificate | bounded-only | exact d values |")
    L.append("|---|---|---|---|---|---|")
    for g in ["Z12xZ12", "Z18xZ8", "Z24xZ6", "Z36xZ4", "Z16xZ9", "Z48xZ3"]:
        v = by_group.get(g, [])
        if not v:
            continue
        t = Counter(tier(r) for r in v)
        ds = sorted(r["d_exact"] for r in v if r["d_exact"] is not None)
        L.append(f"| {g} | {len(v)} | {t['solver-exact']} | {t['certificate']} "
                 f"| {t['bounded-only']} | {ds} |")
    L.append("")

    L.append("## Notable rows\n")
    def fmt(r):
        d = (f"d = {r['d_exact']} ({tier(r)})" if r["d_exact"] is not None
             else r.get("d_status", "UNKNOWN"))
        return (f"- `{r['code_id']}` [[{r['n']},{r['k']}]] {d}; "
                f"A = `{r['A']}`, B = `{r['B']}`; d_ub = {r['d_ub']}; "
                f"floor = {r.get('floor')}; {r['wall_s']} s")
    hi = [r for r in rows if (r["d_exact"] or 0) >= 12 or
          (r["d_exact"] is None and (r.get("floor") or 0) >= 10)]
    hi.sort(key=lambda r: -(r["d_exact"] or r.get("floor") or 0))
    L.append("### High-distance (exact d >= 12 or floor >= 10)\n")
    L.extend(fmt(r) for r in hi) if hi else L.append("(none)")
    L.append("")
    hi_ub = [r for r in rows if r["d_exact"] is None and r["d_ub"] >= 14]
    hi_ub.sort(key=lambda r: -r["d_ub"])
    L.append("### Open tail (bounded-only, d_ub >= 14) — the interesting frontier\n")
    L.extend(fmt(r) for r in hi_ub) if hi_ub else L.append("(none)")
    L.append("")
    cert = [r for r in rows if tier(r) == "certificate"]
    L.append("### Certified doubles\n")
    if cert:
        for r in cert:
            v = r["lanes"]["certify"].get("verdict", {})
            L.append(fmt(r))
            L.append(f"  - base: {v.get('base')}, d_base: "
                     f"{(v.get('distance') or {}).get('d_base')}")
    else:
        L.append("(none — see certify-lane outcome mix below)")
    L.append("")
    anom = [r for r in rows if "ANOMALY" in r]
    L.append("### Anomalies\n")
    L.extend(f"- `{r['code_id']}`: {r['ANOMALY']}" for r in anom) if anom else L.append("(none)")
    L.append("")

    L.append("## Certify-lane outcome mix\n")
    cc = Counter()
    for r in rows:
        c = r["lanes"].get("certify", {})
        o = c.get("outcome", "?")
        if o.startswith("skip: no (R)"):
            o = "no (R) candidate"
        elif o.startswith("skip: base d_ub"):
            o = "base out of scope (d_ub > 15)"
        elif "FALLBACK" in o:
            log = ((c.get("verdict") or {}).get("candidate_log") or [])
            first = str(log[0].get("outcome", "")) if log else ""
            if "k =" in first:
                o = "front-end k cap (k > 14)"
            elif "census" in first:
                o = "census too large at required W (front-end scope)"
            else:
                o = "FALLBACK (budget/other)"
        cc[o] += 1
    for o, n in cc.most_common():
        L.append(f"- {o}: {n}")
    L.append("")
    L.append("Notes: the 3 DOUBLING-REFUTED verdicts are certified negatives "
             "('safe-class coset of weight d_base exists', so the safe-floor "
             "doubling fails on that axis); two of those rows were "
             "independently closed by SAT at exactly d = 8, consistent with "
             "the refutation's 'd < 2*d_base on the safe sector'. The 1 "
             "error row is an infrastructure timeout of certify's inner "
             "`cosetbz` native subprocess (bb_samp_Z36xZ4_23206eb7); the row "
             "fell through to the SAT lane normally.\n")

    L.append("## Base cache (shared quotient ladder, order 72 / 36 / 18)\n")
    bc = con.execute(
        "SELECT label, count(*), sum(CASE WHEN d_exact IS NOT NULL THEN 1 ELSE 0 END), "
        "sum(CASE WHEN source LIKE 'corpus%' THEN 1 ELSE 0 END), "
        "round(sum(wall_s), 1) FROM base_cache GROUP BY label ORDER BY label").fetchall()
    L.append("| quotient group | entries | with d_exact | from corpus | compute wall s |")
    L.append("|---|---|---|---|---|")
    for r in bc:
        L.append("| " + " | ".join(str(x) for x in r) + " |")
    n_entries, tot_wall = con.execute(
        "SELECT count(*), round(sum(wall_s),1) FROM base_cache").fetchone()
    L.append(f"\nCache: {n_entries} distinct quotient codes, {tot_wall} s total "
             "compute (amortized across covers; corpus hits free). Ladder d "
             "values are context, NOT floors/bounds for the cover (no bound "
             "transport claimed).\n")

    L.append("## Method / trust-tier legend\n")
    L.append("- **solver-exact** — CMS SAT ladder in `bb_lab.sat_distance`: "
             "witness at d + solver-proved UNSAT at every w < d (solver "
             "trust, same tier as corpus CaDiCaL d_exact). Rows closed as "
             "'UNSAT floor + verified L1 witness' are the same two sides "
             "established by different engines.")
    L.append("- **certificate** — `bb_lab.doubling_certify` (A30 front-end): "
             "d = 2 * d_base by the doubling theorem + counting-invariant "
             "enumeration; certificate tier, not kernel-checked.")
    L.append("- **bounded-only** — d in [floor, d_ub]: floor only when the "
             "solver genuinely refuted all weights < floor (contiguous "
             "UNSAT prefix); d_ub from L1 sampling with a verified logical "
             "witness. Witness weights are never reported as floors or "
             "exact. UNKNOWN stays UNKNOWN.\n")

    L.append("## Smoke tests (pre-batch)\n")
    L.append("- fast canonicalization == `bb_lab.canonical.canonical_bits` on "
             "180 random pairs (Z6xZ6, Z12xZ3, Z6xZ4): 0 mismatches; all 13 "
             "corpus Z12xZ12 exact rows round-trip to their instance_ids.")
    L.append("- SAT lane: gross base bb72 [[72,12,6]] -> DISTANCE 6, clean "
             "UNSAT ladder, 0.5 s.")
    L.append("- certify lane: gross [[144,12,12]] re-derived CERTIFIED "
             "d = 12 (= 2 * 6) in 3.7 s.")
    L.append("- L1 at n = 288: 100k samples ~5 s; finds d_ub = 18 on "
             "[[288,12,18]] (its true d).\n")

    L.append("## Caveats\n")
    L.append("- Sampled (not exhaustive) enumeration: `bb_samp_` provenance; "
             "orbit-size-biased uniform sampling per a18 conventions.")
    L.append("- Per-code 120 s cap means the interesting tail (d_ub >= 12) "
             "is mostly bounded-only; floors reflect only solver-completed "
             "UNSAT rounds at n = 288 (w <= ~9 typically).")
    L.append("- L1 d_ub at n = 288 with 100k samples is loose for low-d "
             "codes (SAT closed those anyway); for the open tail treat d_ub "
             "as an upper bound, not an estimate.")
    L.append("- Certify-lane REFUSED/FALLBACK outcomes are scope refusals, "
             "not refutations, EXCEPT explicit DOUBLING-REFUTED verdicts "
             "(certified negatives for that axis).")
    L.append("- Quotient-ladder d values are informational context; no "
             "cover-distance bound is claimed from them.")
    L.append("- Machine was contended by other sessions' solvers "
             "(load ~5-6); wall times are upper estimates of uncontended cost.\n")

    L.append("## Reproduction\n")
    L.append("```bash")
    L.append("cd experiments/bb_lab")
    L.append("uv run python data/order144_sweep/smoke.py          # lane smokes")
    L.append("uv run python data/order144_sweep/sampler.py        # sampled enumeration -> sweep.duckdb")
    L.append("caffeinate -ims uv run python data/order144_sweep/driver.py   # batch (resumable; skips done rows)")
    L.append("uv run python data/order144_sweep/report.py         # this report")
    L.append("```")
    L.append("\nMain corpus DB was opened READ-ONLY for dedup + quotient "
             "lookups; nothing under the main checkout was written. Skipped "
             "as already-exact: 13 corpus Z12xZ12 rows (incl. [[288,12,18]]) "
             "+ published IBM [[288,8,20]] (Z18xZ8) by canonical id.\n")

    (HERE / "REPORT.md").write_text("\n".join(L))
    print(f"wrote REPORT.md ({len(rows)} rows)")


if __name__ == "__main__":
    main()
