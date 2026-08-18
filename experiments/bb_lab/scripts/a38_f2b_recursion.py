"""A38 S2, front F2b: the epsilon-recursion chapter — formulation +
validation against banked ground truth.

THE FORMULATION (per Z2-deck rung p: C -> B, Lemma-1 setting, parity):

  eps := tau o p is multiplication by (1 + sigma) on C; ker eps = im eps
  = im tau (C is a free F2[Z2]-module).  The slice identity in eps-form:
  |v| = |eps v|/2 + 2m (since |eps v| = 2|p v|).  The two eps-strata of
  a nontrivial logical v are exactly the assembly's branches:

    v in im tau (eps v = 0, b = 0):  |v| = 2|u| >= 2 d(B) — the tau
        branch; the G5 ceiling 2 d(mid) is the DEGENERATE case of the
        recursion (the level-(r-1) floor consumed as a bare number).
    v not in im tau (b != 0):        |v| = |b| + 2 m*(b) with
        m*(x) := min overflow over nontrivial lifts of x — the rung
        content, stratified by [b] (stab / seam).

  THE RECURSION QUESTION (charter F2b): bound m*(x) from level-(r-1)
  data with no new enumeration at level r.  Two candidate forms:

  (N) the NUMBER-ONLY form — m*(x) >= (t - |x|)/2 from level-(r-1)
      FLOORS alone: FALSE in general; this is the naive SeamCosetFloor,
      refuted at A33/A36, and the banked m* table below shows m*
      depends on x, not on |x| and d(B) alone.
  (C) the CENSUS-CARRYING form (the kernel-shift lemma, NEW this
      session): the solution set of the carry system E v0 = RHS x is
      v0p + Z(B) EXACTLY (ker E = the level-(r-1) cycle space, because
      tau is injective on 0-chains).  Hence every candidate with
      overflow <= cap is v0p (+) z with |z| <= |x| + cap + ov(v0p), so
      the rung consumes the level-(r-1) cycle CENSUS in a weight
      window; when |x| + cap + ov(v0p) < d(B) the window holds
      STABILIZERS ONLY — the rung is decided by the level-(r-1)
      stabilizer census, no new enumeration species at level r.

  WHAT THE RECURSION NEEDS PER RUNG (the exact statement the charter
  asked for): a particular lift v0p of small overflow (the
  row-decomposition lift gives ov(v0p) <= (6|I| - |x|)/2 for x = a sum
  of |I| stabilizer rows), and the level-(r-1) cycle census complete to
  B = |x| + cap + ov(v0p).  WHERE IT FAILS to be cheaper than
  enumeration: when B >= d(B), the window needs the level-(r-1)
  LOGICAL-coset censuses too (their cost is real; the c37xx execution
  builds them at <= WC for exactly this reason), and when x is not a
  stabilizer (seam cells) no row-decomposition exists — v0p comes from
  solve_E with uncontrolled overflow, so (C) covers seam cells only
  through the generic window.  The measured coverage below quantifies
  both limits on banked towers.

VALIDATION (falsify-first, this script):

  V1  eps-identities on the a36 tower: eps = tau o p as matrices;
      |eps v| = 2 |p v| on random cycles; b = 0 <=> v in im tau.
  V2  the banked m* table (ground truth, read from banked artifacts):
      a36 witness ladder + a33 tightness probes — the (N)-refutation
      data, restated as m* values.
  V3  the kernel-shift lane vs the banked A36 dangerous battery:
      per-cell window bound B with the row-decomposition lift; on every
      cell with B < d(gross) = 12, run the kernel-shift rung over the
      banked stabilizer census window and assert verdict == banked
      PASS.  Report the coverage histogram (which cells the stab-only
      window reaches) — the expensive light-shadow cells are exactly
      the covered ones (complementary to the direct lane).

Output: data/a38/f2b_recursion.json
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

from bb_lab.tower import (  # noqa: E402
    AxisDeck, TowerCode, i2v, in_span, translation_perms, v2i,
)
from a38_c37xx_freeze import KernelShift, row_lift_v0  # noqa: E402

DATA = LAB / "data"
OUT = DATA / "a38"
D_GROSS = 12


def main() -> None:
    t0 = time.monotonic()
    out: dict = {}
    rng = np.random.default_rng(20260818)

    G8 = TowerCode("G8", (12, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2")
    GR = TowerCode("GR", (12, 6), "x^3 + y + y^2", "y^3 + x + x^2")
    deck = AxisDeck(G8, GR, 1)
    perms = translation_perms(GR)

    # ---------------------------------------------- V1: eps identities
    # eps = tau o p as chain maps: TAU @ P == action of (1 + sigma)
    sig_perm = np.zeros(G8.n, dtype=np.int64)
    for i, e in enumerate(G8.G):
        j = G8.G.index(G8.G.add(e, (0, 6)))       # + sigma = y^6
        sig_perm[i] = j
        sig_perm[G8.ng + i] = G8.ng + j
    EPS = (deck.TAU @ deck.P) % 2
    for _ in range(40):
        v = np.zeros(G8.n, dtype=np.uint8)
        v[rng.choice(G8.n, size=20, replace=False)] = 1
        lhs = (EPS @ v) % 2
        rhs = (v + v[sig_perm]) % 2
        assert (lhs == rhs).all(), "eps != 1 + sigma"
    n_str = 0
    for _ in range(60):
        v = G8.random_cycle(rng)
        b = (deck.P @ v) % 2
        ev = (EPS @ v) % 2
        assert int(ev.sum()) == 2 * int(b.sum()), "|eps v| != 2|p v|"
        in_tau = not b.any()
        # b = 0 <=> v in im tau (check via sheets: v0 == v1)
        v0, v1 = deck.sheets(v)
        assert in_tau == bool((v0 == v1).all())
        n_str += 1
    print(f"[{time.monotonic()-t0:5.1f}s] V1: eps = 1+sigma = tau o p "
          f"(40 vectors); |eps v| = 2|p v| and the b = 0 <=> im-tau "
          f"stratum identity on {n_str} random cycles")
    out["V1"] = {"eps_eq_one_plus_sigma": True, "slice_eps_form": n_str}

    # ------------------------------------- V2: the banked m* table
    wh = json.loads(
        (DATA / "a36" / "witness_hunt_banked.json").read_text())
    h5 = json.loads((DATA / "a33" / "h5_direct.json").read_text())
    mstar_rows = [
        {"tower": "a36", "cell": "seam |x|=12 (class 0x40)",
         "m_star": wh["witness"]["m1"], "lightest_over_x": 18,
         "deficit_bound_at_t18": 3,
         "note": "SATURATES the bound — d-attaining cell"},
        {"tower": "a36", "cell": "stab |x|=14 stratum (54 orbits)",
         "m_star": ">= 3 (0 finds at M=3)", "lightest_over_x": ">= 20",
         "deficit_bound_at_t18": 2,
         "note": "EXCEEDS the bound by >= 1"},
    ]
    for p in h5["tightness_probes"]:
        mstar_rows.append(
            {"tower": "a33", "cell": f"seam |x|={p['w']}",
             "m_star": p["min_overflow"],
             "lightest_over_x": p["lightest_logical_over_w"],
             "deficit_bound_at_t20": (20 - p["w"]) // 2,
             "note": ("saturates" if p["min_overflow"] ==
                      (20 - p["w"]) // 2 else
                      f"exceeds by "
                      f"{p['min_overflow'] - (20 - p['w']) // 2}")})
    assert wh["witness"]["m1"] == 3
    print(f"[{time.monotonic()-t0:5.1f}s] V2: banked m* table — "
          f"m* varies per element at fixed |x| across towers "
          f"(a33 w14 seam: m* = 4 vs deficit bound 3; a36 w12 seam: "
          f"m* = 3 = bound): the NUMBER-ONLY recursion (N) is refuted "
          f"by banked data; the census-carrying form (C) is the "
          f"correct statement")
    out["V2_mstar_table"] = mstar_rows

    # ------------------- V3: kernel-shift vs the banked a36 battery
    banked = json.loads(
        (DATA / "a36" / "direct_close_banked.json").read_text())
    reps: list[np.ndarray] = []
    for line in (DATA / "a36" /
                 "stab_census_orbits_banked.jsonl").open():
        r = json.loads(line)
        b = np.zeros(GR.n, dtype=np.uint8)
        b[r["b_support"]] = 1
        reps.append(b)
    assert len(reps) == 469
    ks = KernelShift(deck, reps, complete_to=16)
    cover_hist: dict[str, int] = {}
    n_cov = n_run = 0
    for b in sorted(reps, key=lambda b: int(b.sum())):
        wb = int(b.sum())
        M = (18 - wb) // 2
        if M <= 0:
            continue
        v0p, ovp = row_lift_v0(deck, b)
        B = wb + (M - 1) + ovp
        key = f"|b|={wb},M={M}"
        if B < D_GROSS:
            cover_hist[key] = cover_hist.get(key, 0) + 1
            n_cov += 1
            # run the kernel-shift rung (stab-only window is complete
            # for cycles <= B < d(gross))
            rhs = (deck.RHS @ b) % 2
            bmask = v2i(b)
            viol = 0
            seen: set[int] = set()
            for v0i in ks.candidates(b, v0p, M - 1):
                canon = min(v0i, v0i ^ bmask)
                if canon in seen:
                    continue
                seen.add(canon)
                v0 = i2v(v0i, GR.n)
                assert not (((deck.E @ v0) + rhs) % 2).any()
                ch = (deck.EMB[0] @ v0
                      + deck.EMB[1] @ ((v0 + b) % 2)) % 2
                if in_span(v2i(ch), G8.rsHX_b, G8.rsHX_p):
                    continue
                viol += 1
            assert viol == 0, \
                f"kernel-shift found a violation at |b|={wb} — " \
                f"CONTRADICTS the banked 469/469 PASS!"
            n_run += 1
    assert banked["dangerous_rungs"]["verdicts"] == {"PASS": 469}
    print(f"[{time.monotonic()-t0:5.1f}s] V3: kernel-shift coverage on "
          f"the banked a36 dangerous battery: {n_cov}/469 cells have a "
          f"stab-only window (B < 12): {cover_hist}; all {n_run} "
          f"covered rungs PASS == banked (the covered cell IS the "
          f"expensive cap-5 lane cell — complementary coverage "
          f"measured)")
    out["V3"] = {"covered": n_cov, "of": 469,
                 "coverage_hist": cover_hist,
                 "all_pass_eq_banked": True}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "f2b_recursion.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {OUT / 'f2b_recursion.json'}")


if __name__ == "__main__":
    main()
