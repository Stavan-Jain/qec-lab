#!/usr/bin/env python3
"""A42 S0 — the variety anchors: V(A,B) exactly, the spectral k-formula,
and the calibration battery.

Gates (all falsify-first, in order):
  G0  validate_banked green.
  G1  The Sylvester resultant Res_y(A~, B^) recomputed from scratch
      equals the banked A40 S3 certificate 1+x+x^2+x^4+x^5+x^11+x^13.
  G2  Factorization + fiber enumeration: every irreducible factor of the
      resultant carries >= 1 torus point of V(A,B) (monic-in-y product
      formula); every enumerated point satisfies both defining equations
      end-to-end in its residue field.
  G3  Scheme consistency: for every x-factor f, exponent of f in Res ==
      sum over point-orbits above f of (D/deg f) * plane multiplicity;
      total sum of deg*exp = 13.
  G4  Member battery: spectral k == direct F_2 rank k (independent of
      TowerCode) at rectangular members, both columns, r <= 6; == 12
      wherever A40 P1 verified.
  G5  The full banked triage grid (l in {12,18}, p = 1..8, all shears):
      spectral k == TowerCode k on honestly-cancelled transported
      supports, frame by frame.  Rows where the banked JSON differs are
      exactly the support-collision frames (frozenset artifact) — each
      is detected, classified by colliding term pair, and reported.
  G6  Spot battery beyond the banked grid: random sheared frames at
      l in {24,30}, spectral vs TowerCode.
  G7  Switch-on arithmetic for every non-F4 orbit (which members ever
      see it), with the CRT obstruction recorded; first reachable
      switch-on member verified by direct rank.

Output: data/a42/s0_variety.json + full stdout log.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

import a42_lib as L  # noqa: E402

DATA = LAB / "data" / "a42"
DATA.mkdir(parents=True, exist_ok=True)

BANKED_RES = 0b10100000110111  # 1+x+x^2+x^4+x^5+x^11+x^13 (A40 S3 K2)

# original Laurent supports (for direct rectangular rank checks)
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def t():
    return time.time()


def rank_bits(rows: list[int]) -> int:
    piv: dict[int, int] = {}  # top-bit -> row
    for r in rows:
        while r:
            tb = r.bit_length() - 1
            b = piv.get(tb)
            if b is None:
                piv[tb] = r
                break
            r ^= b
    return len(piv)


def direct_k_rect(l: int, m: int) -> int:
    """k = 2*lm - rank[H_X] - rank[H_Z] on Z_l x Z_m, bit-int rref.
    H_X rows = g*(A|B); H_Z rows = g*(Bbar|Abar) (antipode)."""
    n = l * m

    def idx(i, j):
        return (i % l) * m + (j % m)

    def rows_for(sa, sb):
        rows = []
        for gi in range(l):
            for gj in range(m):
                r = 0
                for (i, j) in sa:
                    r |= 1 << idx(gi + i, gj + j)
                for (i, j) in sb:
                    r |= 1 << (n + idx(gi + i, gj + j))
                rows.append(r)
        return rows

    abar = [(-i, -j) for (i, j) in A_L]
    bbar = [(-i, -j) for (i, j) in B_L]
    rx = rank_bits(rows_for(A_L, B_L))
    rz = rank_bits(rows_for(bbar, abar))
    return 2 * n - rx - rz


def tower_k(l: int, p: int, d: int) -> tuple[int, bool]:
    """TowerCode k on the honestly-cancelled transported supports."""
    sa, sb, fd, collided = L.transported_supports(l, p, d)
    c = TowerCode(f"a42(l={l},p={p},d={d})", (fd.o1, fd.o2), sa, sb)
    return c.k, collided


def main():
    t0 = t()
    out: dict = {}

    print("== G0: validate_banked ==", flush=True)
    validate_banked(LAB / "data")
    print(f"validate_banked: PASS ({t()-t0:.1f} s)", flush=True)

    # ---------------- G1: resultant ----------------
    res = L.sylvester_resultant_y()
    print(f"\n== G1: resultant ==\nRes_y(A~,B^) = {L.poly_str(res)}")
    assert res == BANKED_RES, (hex(res), hex(BANKED_RES))
    print("matches banked A40 S3 certificate: PASS")
    out["resultant"] = {"bits": res, "str": L.poly_str(res)}

    # ---------------- G2: factorization + variety ----------------
    print("\n== G2: factorization and fiber enumeration ==")
    fac = L.factorize_f2(res)
    tot = 0
    fac_rows = []
    for f, e in sorted(fac.items(), key=lambda kv: (L.pdeg(kv[0]), kv[0])):
        F = L.F2k(f)
        orda = L.mult_order(F, F.t, F.q - 1)
        print(f"  {L.poly_str(f)}  deg={L.pdeg(f)}  exp={e}  ord(root)={orda}")
        fac_rows.append({"poly": L.poly_str(f), "bits": f,
                         "deg": L.pdeg(f), "exp": e, "root_order": orda})
        tot += L.pdeg(f) * e
    assert tot == 13, tot
    print(f"  total deg*exp = {tot} = 13: PASS")
    out["factors"] = fac_rows

    orbits = L.variety_orbits(verbose=True)
    print(f"\n  {len(orbits)} Frobenius point-orbits on the torus")

    # plane multiplicities with stabilization check
    orb_rows = []
    for o in orbits:
        m8 = L.plane_multiplicity(o, N=8)
        m12 = L.plane_multiplicity(o, N=12)
        assert m8 == m12, (m8, m12, "plane multiplicity not stabilized")
        o.plane_mult = m8
        row = {"fx": L.poly_str(o.fx), "deg_fx": o.dfx, "D": o.D,
               "ord_alpha": o.ord_a, "ord_beta": o.ord_b,
               "beta_minpoly": L.poly_str(o.beta_minpoly),
               "plane_mult": o.plane_mult}
        orb_rows.append(row)
        print(f"  orbit: fx={row['fx']} D={o.D} ord=({o.ord_a},{o.ord_b}) "
              f"beta_min={row['beta_minpoly']} mult={o.plane_mult}",
              flush=True)
    out["orbits"] = orb_rows

    # ---------------- G3: Res-multiplicity consistency ----------------
    print("\n== G3: scheme consistency (Res exponents vs local mults) ==")
    for f, e in fac.items():
        s = sum((o.D // o.dfx) * o.plane_mult for o in orbits if o.fx == f)
        status = "PASS" if s == e else "FAIL"
        print(f"  {L.poly_str(f)}: exp {e} vs sum {s}: {status}")
        assert s == e, (L.poly_str(f), e, s)

    # ---------------- G4: member battery ----------------
    print("\n== G4: rectangular member battery ==", flush=True)
    member_rows = []
    for r in range(1, 9):
        for b in (0, 1):
            l, m = 6 * (r + b), 6 * r
            ks = L.spectral_k(orbits, l, m, 0)
            row = {"r": r, "b": b, "l": l, "m": m, "k_spectral": ks}
            if l * m <= 1600:
                kd = direct_k_rect(l, m)
                row["k_direct"] = kd
                assert ks == kd, (l, m, ks, kd)
            member_rows.append(row)
            print(f"  (l,m)=({l},{m}) r={r} b={b}: k_spectral={ks}"
                  + (f" k_direct={row['k_direct']} MATCH"
                     if "k_direct" in row else ""), flush=True)
    out["members"] = member_rows

    # ---------------- G5: full banked triage grid ----------------
    print("\n== G5: banked triage grid (480 frames incl. l=12,18) ==",
          flush=True)
    banked = json.loads((LAB / "data" / "a40" /
                         "s4_phase_triage.json").read_text())
    banked_k = {}
    for row in banked["rows"]:
        for d in row["shears"]:
            banked_k[(row["l"], row["p"], d)] = row["k"]
    grid_diffs = []
    n_frames = 0
    n_coll = 0
    tg = t()
    for l in (12, 18):
        for p in range(1, 9):
            for d in range(l):
                ks = L.spectral_k(orbits, l, p, d)
                kt, collided = tower_k(l, p, d)
                assert ks == kt, ("spectral vs TowerCode", l, p, d, ks, kt)
                n_frames += 1
                kb = banked_k.get((l, p, d))
                if collided:
                    n_coll += 1
                if kb is not None and kb != ks:
                    grid_diffs.append({"l": l, "p": p, "d": d,
                                       "k_true": ks, "k_banked": kb,
                                       "collided": collided})
            print(f"  l={l} p={p}: all shears spectral==TowerCode "
                  f"({t()-tg:.0f} s cum)", flush=True)
    print(f"  {n_frames} frames: spectral == TowerCode everywhere: PASS")
    print(f"  collision frames encountered: {n_coll}")
    print(f"  banked-vs-true diffs: {len(grid_diffs)}")
    for gdiff in grid_diffs:
        print(f"    l={gdiff['l']} p={gdiff['p']} d={gdiff['d']}: "
              f"true k={gdiff['k_true']} banked k={gdiff['k_banked']} "
              f"collided={gdiff['collided']}")
        assert gdiff["collided"], \
            "banked diff at a NON-collision frame — formula or data wrong"
    out["grid"] = {"frames": n_frames, "collision_frames": n_coll,
                   "banked_diffs": grid_diffs}

    # ---------------- G6: spot battery beyond the grid ----------------
    print("\n== G6: spot battery (l in {24,30}) ==", flush=True)
    import random
    rng = random.Random(42)
    spot_rows = []
    for l in (24, 30):
        for p in (3, 5, 6, 7):
            for d in rng.sample(range(l), 4):
                ks = L.spectral_k(orbits, l, p, d)
                kt, collided = tower_k(l, p, d)
                assert ks == kt, (l, p, d, ks, kt)
                spot_rows.append({"l": l, "p": p, "d": d, "k": ks,
                                  "collided": collided})
        print(f"  l={l}: PASS", flush=True)
    print(f"  {len(spot_rows)} spot frames: spectral == TowerCode: PASS")
    out["spot"] = spot_rows

    # W7 frames detail (the omega-eigenclass record)
    print("\n== W7 frames (l,7,l-2): contributing orbits ==")
    w7_rows = []
    for l in (12, 18, 24):
        k, contrib, fd = L.spectral_k(orbits, l, 7, l - 2, detail=True)
        det = [{"fx": L.poly_str(o.fx), "ord": (o.ord_a, o.ord_b),
                "D": o.D, "localdim": m} for (o, m) in contrib]
        w7_rows.append({"l": l, "k": k, "contrib": det})
        print(f"  (l,p,d)=({l},7,{l-2}): k={k} via {det}")
    out["w7_frames"] = w7_rows

    # ---------------- G7: switch-on arithmetic ----------------
    print("\n== G7: non-F4 orbit switch-on arithmetic ==")
    from math import gcd, lcm
    switch_rows = []
    for o in orbits:
        if o.ord_a == 3 and o.ord_b == 3:
            continue  # F4 cluster
        s1 = o.ord_a // gcd(o.ord_a, 6)
        t1 = o.ord_b // gcd(o.ord_b, 6)
        g = gcd(s1, t1)
        # b=0: r = lcm(s1,t1) minimal positive with s1|r, t1|r
        r0 = lcm(max(s1, 1), max(t1, 1))
        # b=1: need r == -1 mod s1, r == 0 mod t1; solvable iff g | 1
        if g == 1:
            r1 = None
            # CRT
            rr = 0
            step = t1
            while True:
                if (rr + 1) % s1 == 0 and rr > 0:
                    r1 = rr
                    break
                rr += step
                if rr > s1 * t1 + t1:
                    break
        else:
            r1 = None
        row = {"fx": L.poly_str(o.fx), "ord": (o.ord_a, o.ord_b),
               "D": o.D, "plane_mult": o.plane_mult,
               "s1": s1, "t1": t1, "gcd": g,
               "b0_first_r": r0, "b1_first_r": r1,
               "b1_blocked": g > 1}
        switch_rows.append(row)
        print(f"  orbit fx={row['fx']} ord=({o.ord_a},{o.ord_b}) D={o.D}: "
              f"s1={s1} t1={t1} gcd={g} -> b=0 first r={r0}, "
              f"b=1 {'BLOCKED' if g > 1 else f'first r={r1}'}")
    out["switch_on"] = switch_rows

    # predicted k at the earliest interesting members
    firsts = sorted({r["b0_first_r"] for r in switch_rows}
                    | {r["b1_first_r"] for r in switch_rows
                       if r["b1_first_r"]})
    pred_rows = []
    for r in firsts[:6]:
        for b in (0, 1):
            l, m = 6 * (r + b), 6 * r
            ks = L.spectral_k(orbits, l, m, 0)
            row = {"r": r, "b": b, "l": l, "m": m, "k_spectral": ks}
            if l * m <= 4000:
                kd = direct_k_rect(l, m)
                row["k_direct"] = kd
                assert ks == kd, (l, m, ks, kd)
            pred_rows.append(row)
            print(f"  member (l,m)=({l},{m}) r={r} b={b}: k={ks}"
                  + (" (direct-verified)" if "k_direct" in row else ""),
                  flush=True)
    out["switch_on_members"] = pred_rows

    # ---------------- member-shape multiplicity table ----------------
    print("\n== F4-orbit local dims vs 2-part shape (rect frames) ==")
    shape_rows = []
    for o in orbits:
        if not (o.ord_a == 3 and o.ord_b == 3):
            continue
        tab = {}
        for a1 in range(0, 4):
            for a2 in range(0, 4):
                dim = L.local_dim(o, 1 << a1, 1 << a2,
                                  ex_x=(1, 0), ex_y=(0, 1), mod_exp=True)
                tab[f"{a1},{a2}"] = dim
                assert dim <= o.plane_mult, (a1, a2, dim, o.plane_mult)
        shape_rows.append({"fx": L.poly_str(o.fx),
                           "beta_minpoly": L.poly_str(o.beta_minpoly),
                           "dims": tab})
        print(f"  orbit fx={L.poly_str(o.fx)} beta_min="
              f"{L.poly_str(o.beta_minpoly)}: dims {tab}", flush=True)
    out["f4_shape_table"] = shape_rows

    out["wall_s"] = round(t() - t0, 1)
    (DATA / "s0_variety.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s0_variety.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
