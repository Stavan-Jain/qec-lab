#!/usr/bin/env python3
"""A43 S1 — taxonomy of the banked doubly-spanning minimizers by class kind
(S11's pure-x / pure-y / pure-d / mixed via the three cylinder images,
K = 4) x shape (n_x, n_y, max gaps, block split).  Reads
data/a43/s1_l3_shape.json (the S11 populations re-derived with cells kept).
Output: data/a43/s1_taxonomy.json (+ .log) with ASCII pictures of one
representative per (kind, shape) family."""
from __future__ import annotations
import json, sys, time
from collections import Counter
from pathlib import Path
LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src")); sys.path.insert(0, str(LAB / "scripts"))
from bb_lab.tower import rref_ints, in_span  # noqa: E402
import a40_s11_compare as C  # noqa: E402
from a43_s1_l3_shape import picture  # noqa: E402

DATA = LAB / "data" / "a43"; DATA.mkdir(exist_ok=True)
LOG = DATA / "s1_taxonomy.log"
def log(m):
    print(m, flush=True); LOG.open("a").write(m + "\n")

def main():
    LOG.write_text(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')} a43_s1_taxonomy\n")
    cl = json.loads((DATA / "s1_l3_shape.json").read_text())  # cells kept
    out = {}
    for name, lm in [("gross(12,6)", (12, 6)), ("two-gross(12,12)", (12, 12)),
                     ("bb72(6,6)", (6, 6))]:
        code = C.member_code(*lm)
        W = {}
        for ab, tag in [((0, 1), "x"), ((1, 0), "y"), ((1, 1), "d")]:
            b, _ = C.image_classes(code, ab, K=4)
            W[tag] = rref_ints(list(b))
        fam = {}
        pics = {}
        for o in cl["frames"][name]["objects"]:
            sig = o["class"]; kind = "mixed"
            for tag, (bb, pp) in W.items():
                if in_span(sig, bb, pp):
                    kind = "pure-" + tag
            cells = frozenset(tuple(c) for c in o["cells"])
            xs = [c[1] for c in cells]; ys = [c[2] for c in cells]
            n_x, gx = C.gap_structure(xs, lm[0]); n_y, gy = C.gap_structure(ys, lm[1])
            blocks = (sum(1 for c in cells if c[0] == 0), sum(1 for c in cells if c[0] == 1))
            colmax = max(Counter(xs).values()); rowmax = max(Counter(ys).values())
            key = f"{o['sector']}|{kind}|w{o['w']}|nx{n_x}|ny{n_y}|gx{gx}|gy{gy}|blk{blocks}|cmax{colmax}|rmax{rowmax}"
            fam[key] = fam.get(key, 0) + 1
            if key not in pics:
                pics[key] = picture(code, cells)
        out[name] = dict(families=dict(sorted(fam.items(), key=lambda kv: -kv[1])),
                         pictures=pics)
        log(f"=== {name}: {sum(fam.values())} objects, {len(fam)} (kind,shape) families")
        for k, v in sorted(fam.items(), key=lambda kv: -kv[1]):
            log(f"  {v:5d}  {k}")
    (DATA / "s1_taxonomy.json").write_text(json.dumps(out, indent=1))
    log("wrote " + str(DATA / "s1_taxonomy.json"))
if __name__ == "__main__":
    main()
