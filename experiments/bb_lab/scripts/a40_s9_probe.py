#!/usr/bin/env python3
"""A40 S9 — DP sensitivity probe (offline, banked tables only; NO
claim tier).  Question: which enumeration deepenings push the r>=2
DP-free grant to <= 99 quarters (floor 12)?

Scenarios (hypothetical caps injected into the banked link recs):
  base           — banked tables as-is (must reproduce 101 q).
  u1gG           — u=1 link march deepened to gcap G at hcap 19
                   (absent-bucket caps tighten; certs kept banked —
                   optimistic-for-us only through the cap side, so
                   the REAL run must confirm certs stay below).
  u1gG+certJ     — same, plus modeled certs 8h - (3.5h + 3) at
                   h = 8..18 (the measured J(h) trend, pessimistic
                   direction: certs can only raise grants).
  +u2short       — plus a u=2 short-h table (h <= H2, gcap G2,
                   whcap 16 class) with NO certs (pure caps; the
                   real run's certs would raise it — bracketed by
                   +u2certs below).
  +u2certs       — u=2 short-h certs modeled as u=1 certs shifted
                   +1 g per slab (a u=2 link's slabs are all >= 2
                   in pre... rough model, for bracketing only).
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import a40_s7_assemble as A  # noqa: E402


def main():
    s6 = A.load_s6()
    links = A.load_links()
    base = copy.deepcopy(links[14])

    def run(tag, lk):
        fl, det = A.floor_r2(lk, s6)
        print(f"{tag:28s} free {det['dp_free_bc']:6.2f} "
              f"(={det['dp_free_bc']*4:.0f}q) clo "
              f"{det['dp_closure']:5.2f} floor {det['floor']} "
              f"argmax {det['free_bc_argmax']}")
        return det

    run("base (banked)", base)

    for G in (38, 40, 42, 44, 46):
        lk = copy.deepcopy(base)
        for h in range(1, 20):
            lk[1]["gcap_by_h"][h] = max(lk[1]["gcap_by_h"].get(h, 0),
                                        G)
        lk[1]["gcap"] = G
        lk[1]["complete_h"] = 19
        run(f"u1 g{G} (caps only)", lk)

        lk2 = copy.deepcopy(lk)
        for h in range(8, 19):
            gmin = 3.5 * h + 3
            if gmin <= G:
                # modeled cert at (h, L=1, d=0)
                k = (h, 1, 0)
                g = int(round(gmin))
                if k not in lk2[1]["tabL"] or lk2[1]["tabL"][k] > g:
                    lk2[1]["tabL"][k] = g
        run(f"u1 g{G} + certs J-trend", lk2)

    # u2 short-h scenarios on top of u1 g42
    for G in (42,):
        for H2, G2 in ((5, 28), (6, 30), (7, 32)):
            lk = copy.deepcopy(base)
            for h in range(1, 20):
                lk[1]["gcap_by_h"][h] = max(
                    lk[1]["gcap_by_h"].get(h, 0), G)
            lk[1]["gcap"] = G
            lk[1]["complete_h"] = 19
            for h in range(8, 19):
                gmin = 3.5 * h + 3
                if gmin <= G:
                    k = (h, 1, 0)
                    g = int(round(gmin))
                    if k not in lk[1]["tabL"] or \
                            lk[1]["tabL"][k] > g:
                        lk[1]["tabL"][k] = g
            lk[2] = dict(tabL={}, gcap_by_h={h: G2 for h in
                                             range(1, H2 + 1)},
                         kmax=2, whcap=16, dcap=30, srcs=["probe"],
                         trunc_dcap=0, gcap=G2, complete_h=H2,
                         src="probe")
            run(f"u1 g{G}+certs, u2 h<={H2} g{G2}", lk)
            # + modeled u2 certs (u1 certs shifted +1/slab)
            lk3 = copy.deepcopy(lk)
            for (h, L, d), g in base[1]["tabL"].items():
                if h <= H2:
                    k = (h, L, d)
                    g2 = g + h - 1
                    if g2 <= G2 and (k not in lk3[2]["tabL"] or
                                     lk3[2]["tabL"][k] > g2):
                        lk3[2]["tabL"][k] = g2
            run(f"  + modeled u2 certs", lk3)


if __name__ == "__main__":
    main()
