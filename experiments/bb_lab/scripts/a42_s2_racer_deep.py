#!/usr/bin/env python3
"""A42 S2b — the deep racer: session-1 engine + the insertion-time
dedup fix (the identified memory bug-let), pushed as deep as the
budget allows at p = 12.

Fix: session 1 inserted one bucket array per (chunk, input, cost)
slice — ~4096 near-duplicate arrays per chunk — driving the observed
4 GB RSS at level 11 (the canonical state count was only 19,755).
Here every chunk groups its outputs per target cost ACROSS inputs,
canonicalizes, dedups, filters against `seen`, and inserts ONE array.
RSS is checked per chunk (clean abort keeps completed levels
certified); a soft wall-clock budget stops between levels.

Completing level c certifies "every compact cycle of weight <= c
enumerated as a return, with its residue-register value".

** S2 REGISTER CAVEAT (s2_registers.json, note §2.5): the residue
register inherited from S1h is class-blind at 2-adic depth a >= 1
(kernel = ker(pi) on H at a = 1; identically ZERO on H at a = 2).
At p = 12 this run's "nontrivial@[]" is therefore VACUOUS as a
classifier — the run's value is the complete return-weight spectrum
and the engine/memory validation.  The class-complete instrument is
a42_s2_jet_racer.py (two-branch full-depth registers). **

Controls re-run under the new insertion path before production:
p=3 (min nontrivial 6), p=6 (first nontrivial exactly 12),
p=9 to level 16 (no nontrivial return — the banked certificate,
valid at a = 0 where the residue register is complete).
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_atlas9", Path(__file__).parent / "a42_s1_atlas9.py")
A9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(A9)

DATA = LAB / "data" / "a42"

RSS_CAP = int(2.2 * 1024 ** 3)
CHUNK = 2048


class DeepRacer(A9.Racer):

    def run(self, cap: int, log=print, ckpt_path=None,
            time_budget_s: float | None = None, growth_est=2.9):
        t0 = time.time()
        zero = self.pack([0] * 5)
        buckets: dict[int, list] = {0: [np.array([zero],
                                                 dtype=np.uint64)]}
        seen = np.array([zero], dtype=np.uint64)
        returns = {}
        state = {"p": self.p, "cap": cap, "levels": [],
                 "engine": "deep (insertion-dedup)",
                 "lane": "AB (theta'-dual lane by banked duality)"}
        smask = np.uint64((1 << (5 * self.p)) - 1)
        aborted = False
        last_level_t = 0.0
        for c in range(0, cap + 1):
            if time_budget_s is not None and c > 0:
                if (time.time() - t0) + last_level_t * growth_est \
                        > time_budget_s:
                    state["aborted"] = f"time budget before level {c}"
                    log(f"  time budget — stop before level {c} "
                        "(completed levels are certified)")
                    aborted = True
                    break
            lt0 = time.time()
            level_novel = 0
            while buckets.get(c):
                arrs = buckets.pop(c)
                batch = np.unique(np.concatenate(arrs))
                novel = batch[~np.isin(batch, seen)]
                if c == 0 and level_novel == 0:
                    novel = batch
                if novel.size == 0:
                    break
                level_novel += int(novel.size)
                seen = np.union1d(seen, novel)
                for lo in range(0, novel.size, CHUNK):
                    chunk = novel[lo:lo + CHUNK]
                    outs, costs = self.expand(chunk)
                    # group across inputs per target cost
                    per_cost: dict[int, list] = {}
                    for a in range(1 << self.p):
                        ns, co = outs[a], costs[a]
                        if c == 0 and a == 0:
                            keep = ns != zero
                            ns, co = ns[keep], co[keep]
                        for w in np.unique(co):
                            cw = c + int(w)
                            if cw > cap:
                                continue
                            sel = ns[co == w]
                            zs = sel[(sel & smask) == 0]
                            if zs.size and cw > 0:
                                accs = (zs >> np.uint64(5 * self.p)) \
                                    .astype(int)
                                for acv in np.unique(accs):
                                    returns.setdefault(cw, set()) \
                                        .add(int(acv))
                            per_cost.setdefault(cw, []).append(sel)
                    for cw, lst in per_cost.items():
                        merged = np.unique(self.canon(
                            np.concatenate(lst)))
                        merged = merged[~np.isin(merged, seen)]
                        if merged.size:
                            buckets.setdefault(cw, []).append(merged)
                    del outs, costs, per_cost
                    rss = resource.getrusage(
                        resource.RUSAGE_SELF).ru_maxrss
                    if rss > RSS_CAP:
                        state["aborted"] = "RSS"
                        log(f"  RSS cap ({rss//2**20}MB) — clean "
                            "abort (completed levels certified)")
                        aborted = True
                        break
                if aborted:
                    break
                # per-batch compaction of future buckets
                for cw in list(buckets.keys()):
                    arrs2 = buckets[cw]
                    tot = sum(x.size for x in arrs2)
                    if len(arrs2) > 64 or tot > 1 << 22:
                        merged = np.unique(np.concatenate(arrs2))
                        merged = merged[~np.isin(merged, seen)]
                        buckets[cw] = [merged] if merged.size else []
            if aborted and "RSS" in state.get("aborted", ""):
                break
            if aborted:
                break
            last_level_t = time.time() - lt0
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            state["levels"].append(
                {"c": c, "novel": level_novel,
                 "seen": int(seen.size), "rss_mb": rss // 2 ** 20,
                 "t": round(time.time() - t0, 1)})
            nontriv = sorted(w for w, accs in returns.items()
                             if any(acv != 0 for acv in accs))
            log(f"  level {c}: novel {level_novel} seen {seen.size} "
                f"rss {rss//2**20}MB t {time.time()-t0:.0f}s "
                f"returns {sorted(returns)} nontrivial@{nontriv}")
            state["returns"] = {str(k): sorted(v)
                                for k, v in returns.items()}
            state["nontrivial_costs"] = nontriv
            if ckpt_path:
                ckpt_path.write_text(json.dumps(state, indent=1))
        state["wall_s"] = round(time.time() - t0, 1)
        if ckpt_path:
            ckpt_path.write_text(json.dumps(state, indent=1))
        return state


def main():
    # controls under the new insertion path
    for p, cap, expect_min in ((3, 8, 6), (6, 13, 12)):
        r = DeepRacer(p)
        r.validate(nsamples=600)
        st = r.run(cap, log=lambda s: print(s, flush=True))
        nt = st.get("nontrivial_costs", [])
        mn = nt[0] if nt else None
        print(f"control p={p} cap={cap}: min nontrivial {mn} "
              f"(expect {expect_min})", flush=True)
        assert mn == expect_min, (p, mn, expect_min)
    r = DeepRacer(9)
    r.validate(nsamples=600)
    st = r.run(16, log=lambda s: print(s, flush=True))
    assert not st.get("nontrivial_costs"), st["nontrivial_costs"]
    assert max(lv["c"] for lv in st["levels"]) == 16
    print("control p=9 cap 16: no nontrivial return — matches the "
          "banked floor(9) certificate", flush=True)

    # production: p = 12, deep
    cap = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    budget = float(sys.argv[2]) * 3600 if len(sys.argv) > 2 else 3.0 * 3600
    r = DeepRacer(12)
    r.validate(nsamples=600)
    print(f"p=12 production: cap {cap}, time budget {budget/3600:.1f} h",
          flush=True)
    st = r.run(cap, log=lambda s: print(s, flush=True),
               ckpt_path=DATA / "s2_racer12_deep_ckpt.json",
               time_budget_s=budget)
    done = max((lv["c"] for lv in st["levels"]), default=-1)
    nt = st.get("nontrivial_costs", [])
    print(f"p=12: completed level {done}; nontrivial {nt}", flush=True)
    if not nt:
        print(f">>> return-weight spectrum complete through level "
              f"{done}; NO floor claim — the residue register is "
              "class-blind at a = 2 (see the docstring caveat); the "
              "jet racer owns the p = 12 certificate", flush=True)
    (DATA / "s2_racer12_deep.json").write_text(json.dumps(st, indent=1))
    print(f"wrote {DATA/'s2_racer12_deep.json'}", flush=True)


if __name__ == "__main__":
    main()
