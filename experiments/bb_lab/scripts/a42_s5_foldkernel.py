#!/usr/bin/env python3
"""A42 S5 — fold kernels vs transfer images at the b = 1 members.

For the member T = (l, m) = (6r+6, 6r) and its two axis folds
  p_y : T -> (l, m/2)   (deck y^{m/2}),   p_x : T -> (l/2, m)  (deck x^{l/2}),
compare ker(p_*) on H1(T) with S11's transfer images W_x = im pi_{(0,1)*}
(x-windowed classes) and W_y = im pi_{(1,0)*} (y-windowed classes):
the (18,12) run of a42_s5_mixed24.py found ker p_y* = W_x exactly.
Also records sigma_* = id and rank p_* for each deck.
Output: data/a42/s5_foldkernel.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab.tower import (  # noqa: E402
    AxisDeck, TowerCode, h1_map, kernel_ints, rref_ints, span_eq,
    translation_mat,
)
import a40_s11_compare as C  # noqa: E402

DATA = LAB / "data" / "a42"


def main():
    out = {}
    for lm in [(12, 6), (18, 12), (24, 18), (12, 12), (6, 6)]:
        t0 = time.time()
        l, m = lm
        T = C.member_code(l, m)
        Wx, _ = C.image_classes(T, (0, 1), K=4)
        Wy, _ = C.image_classes(T, (1, 0), K=4)
        rec = {"k": T.k, "dim_Wx": len(Wx), "dim_Wy": len(Wy)}
        for axis, name, Wsame, Wother in ((1, "y", "Wx", "Wy"),
                                          (0, "x", "Wy", "Wx")):
            new = list(lm)
            if new[axis] % 2 or new[axis] // 2 < 3:
                rec[f"fold_{name}"] = "no index-2 fold"
                continue
            new[axis] //= 2
            base = TowerCode(f"b{new}", tuple(new), C.red(C.A_L, tuple(new)),
                             C.red(C.B_L, tuple(new)))
            try:
                deck = AxisDeck(T, base, axis)
            except AssertionError as e:
                rec[f"fold_{name}"] = f"deck build failed: {e}"
                continue
            Mp = h1_map(deck)
            ker = kernel_ints(Mp)
            sig_id = bool((translation_mat(T, deck.sigma)
                           == np.eye(T.k, dtype=np.uint8)).all())
            rank = T.k - len(rref_ints(ker)[0]) if ker else T.k
            r = {"base": list(new), "k_base": base.k,
                 "sigma_star_id": sig_id, "rank_p": int(rank),
                 "dim_ker": len(rref_ints(ker)[0]) if ker else 0,
                 "ker_eq_" + Wsame: bool(span_eq(ker, list(Wx if Wsame == "Wx" else Wy))),
                 "ker_eq_" + Wother: bool(span_eq(ker, list(Wy if Wother == "Wy" else Wx)))}
            rec[f"fold_{name}"] = r
        out[str(lm)] = rec
        print(lm, json.dumps(rec), f"{time.time() - t0:.1f}s", flush=True)
    (DATA / "s5_foldkernel.json").write_text(json.dumps(out, indent=1))
    print("->", DATA / "s5_foldkernel.json")


if __name__ == "__main__":
    main()
