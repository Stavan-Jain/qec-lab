"""Regression tests for the doubling-certification front-end.

Fast known-answer anchors:
  * pair72 cover [[72,4,8]] (Z6xZ6): detection must recover the Z3xZ6
    base, and the pipeline must certify d = 8 end-to-end (the PR #53
    Lean instance).
  * f2a6:y rung pass through the LIBRARY path: 113/113 must PASS
    (the Lean dangerousFloorNZ instance, as in scripts/a30_rung_pass.py).
  * by90 x-cover (Z60xZ3): the A14 §13 rung freeze — the front-end must
    NOT certify 16; expected DOUBLING-REFUTED (negative control).
  * detection on the [[360,4,?]] covers recovers the docket bases.

The full [[360,4,?]] certifications (~20 min each) are NOT tests; run
scripts/bb_certify_doubling.py for those.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from bb_lab.doubling_certify import (
    BaseTools, certify, detect, rung_pass,
)
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly

LAB = Path(__file__).resolve().parent.parent


def test_detect_recovers_docket_bases():
    G = AbelianGroup((30, 6))
    A = Poly.from_string("1 + y + x", G)
    B = Poly.from_string("y^4 + x + x^11*y^2", G)
    cands = detect(G, A, B)
    lifts = [c for c in cands if c.R_holds]
    assert len(lifts) == 1
    c = lifts[0]
    assert c.axis == 0 and c.base_group == (15, 6)
    assert c.k_cover == c.k_base == 4

    Gy = AbelianGroup((15, 12))
    Ay = Poly.from_string("1 + y + x", Gy)
    By = Poly.from_string("y^4 + x^8*y^2 + x^13", Gy)
    cands = [c for c in detect(Gy, Ay, By) if c.R_holds]
    assert len(cands) == 1 and cands[0].axis == 1
    assert cands[0].base_group == (15, 6)


def test_pair72_cover_certifies_d8(tmp_path):
    verdict = certify((6, 6), "x^2 + y + y^3", "1 + x + y^2",
                      budget_s=600, threads=8, workdir=tmp_path)
    assert verdict["status"] == "CERTIFIED", verdict.get("stages")
    assert verdict["distance"]["value"] == 8
    assert verdict["distance"]["d_base"] == 4
    assert verdict["stages"]["safe_floor"]["certified"]
    assert verdict["stages"]["rung_pass"]["all_pass"]
    assert verdict["stages"]["witness"]["established"]
    assert verdict["tandem"]["suggested_flags"]["-init-lb"] == 8


def test_f2a6_rung_pass_113(tmp_path):
    G = AbelianGroup((5, 15))
    A = Poly.from_string("1 + y + x", G)
    B = Poly.from_string("x*y^6 + x*y^10 + x^2*y^12", G)
    bt = BaseTools(G, A, B)
    cen = json.load(open(LAB / "data" / "a28" / "census_f2a6.json"))
    ng = G.cardinality
    classes = []
    for e in cen["classes"]:
        v = np.zeros(2 * ng, dtype=np.uint8)
        for g in e["u_support"]:
            v[G.index(tuple(g))] = 1
        for g in e["v_support"]:
            v[ng + G.index(tuple(g))] = 1
        assert int(v.sum()) == e["weight"]
        classes.append({"weight": e["weight"], "vec": v})
    import time
    rp = rung_pass(bt, 1, 8, classes, time.monotonic() + 600, 8, tmp_path)
    rp.pop("engine")
    assert rp["all_pass"] and rp["n_classes"] == 113, rp


def test_by90_rung_cover_is_not_certified(tmp_path):
    # (30,3) with these polynomials IS the [[180,8,12]] rung-1 cover of
    # the Bravyi tower; its base is the [[90,8,8]] bottom over (15,3).
    # A14 §13: the rung freezes at 12 < 16 — the front-end must not
    # certify 16 (expected: safe floor refutation).
    verdict = certify((30, 3), "x^9 + y + y^2", "1 + x^25 + x^26",
                      budget_s=900, threads=8, workdir=tmp_path)
    assert verdict["status"] != "CERTIFIED", verdict["distance"]
    assert verdict["status"] == "DOUBLING-REFUTED", verdict.get(
        "stages", {}).get("candidate_log")


def test_refuses_even_weight(tmp_path):
    v = certify((6, 6), "1 + x", "1 + y", budget_s=60, workdir=tmp_path)
    assert v["status"] == "REFUSED"
