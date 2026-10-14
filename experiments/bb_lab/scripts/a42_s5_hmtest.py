#!/usr/bin/env python3
"""A42 S5 Stage 1 — the coset-leader reformulation of HM, tested on data.

For a period-3q cycle v (3 coprime to q) with omega-support S (fibres
with mu != 0), fold support s, excess eps = |S| - 3q, and the barren
STABILIZER code  Bq := X-cycles = X-stabilizers of the period-q barren
torus/cylinder (H1 = 0 by Theorem A), the S5 reformulation reads

    HM  <=>  for every T <= S:  3 d_S(T) >= |T| - 2 eps,
             d_S(T) := min{|x| : x in 1_T + Bq, x cap S = empty}

(proof: s in Bq, T := s cap S, x := s \\ S).  Two SUFFICIENT local
conditions, checkable per object:
    (B1)  3 d(1_T, Bq) >= |T| - 2 eps      (unrestricted coset leader)
    (B2)  |syn'(1_T)| >= |T| - 2 eps       (barren Z-syndrome weight;
          d >= |syn|/3 since every cell lies in exactly 3 Z-checks)
This script evaluates HM (must hold), B1, B2 at T = s cap S (the
object's own hiding set) and at T = S (the halving-tight stratum) for
  (a) all nontrivial atlas cycles at p = 6 (weight <= 13; the 66 floor
      cycles), embedded on the (24,6) torus with barren torus (24,2);
  (b) the (18,12) weight-24 census objects (s5_mixed24_W24.json and
      s5_dangerous24.json if present), barren torus (18,4).
Coset leaders by exact SAT minimization (n <= 144 variables; solver
tier, but every value is an UPPER bound verified by the returned x and
a LOWER bound by UNSAT at weight-1 — exact).
Output: data/a42/s5_hmtest.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util  # noqa: E402

from bb_lab.tower import TowerCode, v2i  # noqa: E402
from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402
import pycryptosat  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AT = _load("a40_s4_phase_atlas")
DATA = LAB / "data" / "a42"


def coset_leader(HZ: np.ndarray, syn: np.ndarray, forbid: set[int],
                 kmax: int = 40) -> int:
    """min |x| with HZ x = syn and x_f = 0 for f in forbid (exact by
    upward weight iteration; UNSAT at k-1, SAT at k)."""
    n = HZ.shape[1]
    for k in range(0, kmax + 1):
        pool = IDPool()
        qv = [pool.id() for _ in range(n)]
        solver = pycryptosat.Solver()
        for i in range(HZ.shape[0]):
            idx = np.flatnonzero(HZ[i])
            if idx.size:
                solver.add_xor_clause([qv[j] for j in idx], bool(syn[i]))
            elif syn[i]:
                return None
        for f in forbid:
            solver.add_clause([-qv[f]])
        if k < n:
            card = CardEnc.atmost(lits=qv, bound=k, vpool=pool,
                                  encoding=EncType.seqcounter)
            for cl in card.clauses:
                solver.add_clause(cl)
        sat, model = solver.solve()
        if sat:
            x = np.array([1 if model[q] else 0 for q in qv], dtype=np.uint8)
            assert ((HZ @ x) % 2 == syn).all() and int(x.sum()) <= k
            assert not any(x[f] for f in forbid)
            return int(x.sum())
    return None


def analyse(code: TowerCode, bar: TowerCode, v: np.ndarray, q: int, tag):
    """Fibre data of v on the barren torus bar = (l, q)."""
    l = code.G.orders[0]
    cnt: dict[int, int] = {}
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), code.ng)
        x, y = code.G.from_index(gi)
        f = blk * bar.ng + bar.G.index((x % l, y % q))
        cnt[f] = cnt.get(f, 0) + 1
    S = {f for f, c in cnt.items() if c in (1, 2)}
    s = {f for f, c in cnt.items() if c in (1, 3)}
    T = S & s
    eps = len(S) - 3 * q
    w = int(v.sum())
    assert w == 2 * len(S) + 3 * len(s) - 4 * len(T)
    rec = {"tag": tag, "w": w, "S": len(S), "s": len(s), "T": len(T),
           "eps": eps, "n3": len(s - S), "n2": len(S - s)}
    HZ = bar.HZ
    for name, TT in (("T", T), ("S", S)):
        x = np.zeros(bar.n, dtype=np.uint8)
        x[sorted(TT)] = 1
        syn = (HZ @ x) % 2
        d_free = coset_leader(HZ, syn, set())
        d_S = coset_leader(HZ, syn, S)
        need = len(TT) - 2 * eps
        rec[name] = {"size": len(TT), "syn": int(syn.sum()),
                     "d_free": d_free, "d_S": d_S, "need": need,
                     "HM_ok": bool(3 * d_S >= need),
                     "B1_ok": bool(3 * d_free >= need),
                     "B2_ok": bool(int(syn.sum()) >= need),
                     "HM_slack": 3 * d_S - need}
    return rec


def run_p6(log):
    Psupp, Qsupp = AT.PAIRS["AB"]
    rows, _ = AT.atlas("AB", 6, 13, keep_pts=True)
    lstar = 24
    code = TowerCode("t246", (lstar, 6),
                     frozenset((i % lstar, j % 6) for i, j in Psupp),
                     frozenset((i % lstar, j % 6) for i, j in Qsupp))
    bar = TowerCode("t242", (lstar, 2),
                    frozenset((i % lstar, j % 2) for i, j in Psupp),
                    frozenset((i % lstar, j % 2) for i, j in Qsupp))
    assert bar.k == 0, f"barren torus k = {bar.k}"
    log(f"p=6 atlas: {len(rows)} rows; barren torus (24,2) k = 0")
    recs = []
    for r in rows:
        if not r["nontrivial"]:
            continue
        v = np.zeros(code.n, dtype=np.uint8)
        c0 = min(c for c, _, _ in r["pts"])
        for c, y, blk in r["pts"]:
            v[blk * code.ng + code.G.index(((c - c0) % lstar, y % 6))] ^= 1
        assert code.is_cycle(v) and not code.is_stab(v)
        recs.append(analyse(code, bar, v, 2, f"p6w{r['weight']}"))
    return recs


def run_census(log):
    import a40_s11_compare as C
    code = C.member_code(18, 12)
    bar = TowerCode("t184", (18, 4), C.red(C.A_L, (18, 4)),
                    C.red(C.B_L, (18, 4)))
    assert bar.k == 0, f"barren torus (18,4) k = {bar.k}"
    recs = []
    for fn, tagp in (("s5_mixed24_W24.json", "seam"),
                     ("s5_dangerous24.json", "dang")):
        p = DATA / fn
        if not p.exists():
            log(f"  {fn} absent — skipped")
            continue
        d = json.loads(p.read_text())
        objs = d["census"]["objects"]
        log(f"  {fn}: {len(objs)} objects")
        for i, o in enumerate(objs):
            v = np.zeros(code.n, dtype=np.uint8)
            v[o["support"]] = 1
            assert code.is_cycle(v) and not code.is_stab(v)
            recs.append(analyse(code, bar, v, 4,
                                f"{tagp}|{o['kind']}|{o['gap_sector']}"))
            if (i + 1) % 200 == 0:
                log(f"    ... {i + 1}/{len(objs)}")
        # the tau0 family (pure-x pullbacks)
        if "controls" in d:
            for j, t in enumerate(d["controls"].get("tau0_family", [])):
                pass
    return recs


def summarize(recs, log):
    out = {}
    for name in ("T", "S"):
        stats = {"n": 0, "HM_fail": 0, "B1_fail": 0, "B2_fail": 0,
                 "HM_tight": 0, "B1_tight": 0, "B2_tight": 0}
        by_tag: dict = {}
        for r in recs:
            x = r[name]
            stats["n"] += 1
            stats["HM_fail"] += not x["HM_ok"]
            stats["B1_fail"] += not x["B1_ok"]
            stats["B2_fail"] += not x["B2_ok"]
            stats["HM_tight"] += (3 * x["d_S"] == x["need"])
            stats["B1_tight"] += (3 * x["d_free"] == x["need"])
            stats["B2_tight"] += (x["syn"] == x["need"])
            key = (r["tag"].split("|")[0], r["eps"], r["n3"])
            bt = by_tag.setdefault(str(key), {"n": 0, "B1_fail": 0,
                                              "B2_fail": 0, "profiles": set()})
            bt["n"] += 1
            bt["B1_fail"] += not x["B1_ok"]
            bt["B2_fail"] += not x["B2_ok"]
            bt["profiles"].add((x["size"], x["syn"], x["d_free"], x["d_S"],
                                x["need"]))
        for k in by_tag:
            by_tag[k]["profiles"] = sorted(by_tag[k]["profiles"])
        log(f"T = {name}: {stats}")
        for k, bt in sorted(by_tag.items()):
            log(f"   {k}: n={bt['n']} B1_fail={bt['B1_fail']} "
                f"B2_fail={bt['B2_fail']} (|T|,syn,d_free,d_S,need)="
                f"{bt['profiles'][:8]}")
        out[name] = {"stats": stats, "by_tag": by_tag}
    return out


def main():
    logf = (DATA / "s5_hmtest.log").open("a")

    def log(s):
        line = f"[{time.strftime('%H:%M:%S')}] {s}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()
    mode = sys.argv[1] if len(sys.argv) > 1 else "p6"
    out = {}
    if mode in ("p6", "all"):
        recs = run_p6(log)
        log(f"p=6: {len(recs)} nontrivial atlas cycles analysed")
        out["p6"] = summarize(recs, log)
        out["p6"]["records"] = recs
    if mode in ("census", "all"):
        recs = run_census(log)
        log(f"(18,12): {len(recs)} census objects analysed")
        out["census"] = summarize(recs, log)
        out["census"]["records"] = recs
    prev = {}
    if (DATA / "s5_hmtest.json").exists():
        prev = json.loads((DATA / "s5_hmtest.json").read_text())
    prev.update(out)
    (DATA / "s5_hmtest.json").write_text(json.dumps(prev, indent=1))
    log("-> s5_hmtest.json")


if __name__ == "__main__":
    main()
