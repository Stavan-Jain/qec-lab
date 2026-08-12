"""A36 — code-level light-class census + exact whole-orbit scan.

The session's Mode-1 data (T1/T3/T4/T2/T2b sweeps + 12 exact T1.5
verdicts, all REFUTED) suggested the right abstraction: safe-floor
minima are CODE invariants — the min weight of a homology class of
H1(base) — and a presentation cell passes the safe floor iff its
transfer image im Delta (the span of its kappa basis seam classes, a
k/2-dim subspace of the k-dim H1) avoids every "light" class (class
min <= 2d - 2; all class minima are even by the augmentation parity,
A17).  So:

  1. ONE per-code census (complete BZ enumeration over all 2^k - 1
     logical cosets at weight cap W) finds every light class, labelled
     by the symplectic pairing phi(v) = (v . zrep_j)_j in F_2^k — a
     presentation-canonical, stabilizer-invariant label.
  2. Per orbit cell, the kappa basis seams' phi-labels span im Delta
     (64 labels for k = 12) — the cell passes iff span cap lights = 0.
     Pure linear algebra, exact, ~ms per cell: the WHOLE orbit is
     decided exactly, no SAT, no UNDET, no sampling.

Tiering: a W = 14 census is ~walk-cheap and kills most cells (the
observed refutation weights are 10-14); survivors go to the exact
T1.5/certify() ladder which decides the full 2d floor.  A cell passing
the census-W filter is NOT yet SF-certified unless W = 2d - 2.

Usage:
    uv run python scripts/a36_light_census.py census --orders 7,9 \
        --A "..." --B "..." --axis x --d 10 [--W 14]
    uv run python scripts/a36_light_census.py orbit --point T1 --pres 0 \
        --axis x [--W 14] [--full]   # exact scan of the whole v1 orbit
    uv run python scripts/a36_light_census.py controls
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB_ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a36_orbit_screen import (                    # noqa: E402
    POINTS, AxisProblem, parse_poly, poly_str, units,
)
from bb_lab.cosetbz import (                      # noqa: E402
    NOFF_MAX, coset_base, pair_radii, run_window,
)
from bb_lab.doubling_certify import BaseTools     # noqa: E402
from bb_lab.group import AbelianGroup             # noqa: E402
from bb_lab.poly import Poly                      # noqa: E402

DATA_DIR = LAB_ROOT / "data" / "a36"
DATA_DIR.mkdir(parents=True, exist_ok=True)


class CodeCensus:
    """Light-class census of one base code, in the canonical frame
    (doubled axis = axis 0 of `orders`)."""

    def __init__(self, orders: tuple[int, int], A_sup, B_sup, d: int):
        import tempfile
        self.G = AbelianGroup(orders)
        self.A = Poly.from_support([tuple(t) for t in A_sup], self.G)
        self.B = Poly.from_support([tuple(t) for t in B_sup], self.G)
        self.d = d
        self.bt = BaseTools(self.G, self.A, self.B)
        # unique workdir per instance: run_window writes tag-named .mat
        # files, and concurrent censuses in a shared dir would collide
        self.wd = Path(tempfile.mkdtemp(prefix="census_",
                                        dir=str(DATA_DIR)))
        _, _, xreps = self.bt.logical_reps("X")
        _, _, zreps = self.bt.logical_reps("Z")
        self.k = len(xreps)
        assert len(zreps) == self.k, "side k mismatch"
        self.xreps = np.array(xreps, dtype=np.uint8)
        self.zreps = np.array(zreps, dtype=np.uint8)
        # phi(v) = (v . zrep_j)_j; must be injective on H1: the pairing
        # matrix phi(xrep_i) stacked must be invertible over F2
        M = (self.xreps @ self.zreps.T) % 2
        from bb_lab.linalg import rank_f2
        assert rank_f2(M) == self.k, "symplectic pairing degenerate?!"

    def phi(self, v: np.ndarray) -> int:
        bits = (self.zreps @ v) % 2
        out = 0
        for j, b in enumerate(bits):
            if b:
                out |= 1 << j
        return out

    def light_census(self, W: int, deadline_s: float = 1200.0,
                     threads: int = 4) -> dict[int, int]:
        """{phi-label: min found weight <= W} over all 2^k - 1 classes.

        Complete two-window BZ per class (same machinery as the
        front-end's d_side_exact); a class absent from the result has
        class-min > W (and by parity >= W + 2)."""
        import time as _t
        t0 = _t.monotonic()
        deadline = t0 + deadline_s
        bt = self.bt
        k = self.k
        # batched combos + labels (the per-element python loop was the
        # k = 12 bottleneck: 4095 x 2 coset reductions)
        masks = ((np.arange(1, 1 << k, dtype=np.int64)[:, None]
                  >> np.arange(k)[None, :]) & 1).astype(np.uint8)
        combos = (masks @ self.xreps) % 2                    # (N, n)
        lab_bits = (combos @ self.zreps.T) % 2               # (N, k)
        labels = (lab_bits.astype(np.int64)
                  @ (1 << np.arange(k, dtype=np.int64)))
        r1, r2 = pair_radii(W)
        found: dict[int, int] = {}
        for wi, (window, Gs, r) in enumerate(
                [(bt.I1, bt.G1, r1), (bt.I2, bt.G2, r2)]):
            # batched coset_base: systematic window => the sequential
            # reduction is linear, c = t0 + t0[window] @ Gs
            Cb = (combos + (combos[:, window] @ Gs)) % 2
            assert not Cb[:, window].any()
            for spot in (0, min(7, len(combos) - 1)):
                assert (Cb[spot] == coset_base(
                    Gs, window, combos[spot])).all()
            ws = Cb.sum(axis=1)
            for j in np.nonzero(ws <= W)[0]:
                lab = int(labels[j])
                if lab not in found or int(ws[j]) < found[lab]:
                    found[lab] = int(ws[j])
            bases = list(Cb)
            for c0 in range(0, len(bases), NOFF_MAX):
                res = run_window(
                    bt.binp, f"census_w{wi}_c{c0 // NOFF_MAX}", Gs,
                    bases[c0:c0 + NOFF_MAX], r, W, deadline,
                    threads=threads, workdir=self.wd)
                for j, hx in res.pop("hit_rows"):
                    from bb_lab.cosetbz import unpack3
                    v = unpack3(hx, bt.n)
                    w = int(v.sum())
                    lab = int(labels[c0 + j])
                    assert w % 2 == 0, "odd 1-cycle?! parity broken"
                    if lab not in found or w < found[lab]:
                        found[lab] = w
        return found

    # ---- per-cell im Delta ---------------------------------------------
    def imdelta_span(self, ap: AxisProblem, A2, B2) -> list[int]:
        """phi-labels of the nonzero elements of im Delta for this cell.

        Basis seams: seam(zeta_i) for a basis of ker d2 (seam is linear
        in zeta).  Sanity: each seam must be a 1-cycle and nonzero in
        homology unless the whole class collapses (checked)."""
        from bb_lab.fibering import kernel_basis
        A = Poly.from_support([tuple(t) for t in A2], self.G)
        B = Poly.from_support([tuple(t) for t in B2], self.G)
        K = kernel_basis(A, B)
        zg = K.reshape(-1, ap.ell, ap.m)
        SU = ap.seam_grids(A2, zg).reshape(K.shape[0], -1)
        SV = ap.seam_grids(B2, zg).reshape(K.shape[0], -1)
        gens = []
        for i in range(K.shape[0]):
            t = np.concatenate([SU[i], SV[i]]).astype(np.uint8)
            gens.append(self.phi(t))
        span = {0}
        for g in gens:
            span |= {s ^ g for s in span}
        span.discard(0)
        return sorted(span)


def _grid2(v: np.ndarray, ell: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    nb = ell * m
    return v[:nb].reshape(ell, m), v[nb:].reshape(ell, m)


def _transport_to_identity(t: np.ndarray, ap: AxisProblem, u_inv: int,
                           v_inv: int, swap: bool, ta: int, tb: int,
                           sign: int) -> np.ndarray:
    """Map a CELL-space vector into the identity-variant space.

    Cell code = translate(ta, tb) . swap . unit(u, v) applied to the
    identity presentation.  Transport = unit(u_inv, v_inv) . swap .
    translate(sign*ta, sign*tb) on (u-block, v-block) coordinates; the
    translation sign is fixed EMPIRICALLY by verify_transport (the
    hit-multiset equality against per-cell ground truth), not by
    convention archaeology."""
    ell, m = ap.ell, ap.m
    gu, gv = _grid2(t, ell, m)
    gu = np.roll(gu, sign * ta, axis=0)
    gv = np.roll(gv, sign * tb, axis=0)
    if swap:
        gu, gv = gv, gu
    # unit map inverse: coordinate permutation g -> (u_inv*gx, v_inv*gy)
    out_u = np.zeros_like(gu)
    out_v = np.zeros_like(gv)
    for x in range(ell):
        for y in range(m):
            out_u[(u_inv * x) % ell, (v_inv * y) % m] = gu[x, y]
            out_v[(u_inv * x) % ell, (v_inv * y) % m] = gv[x, y]
    return np.concatenate([out_u.reshape(-1), out_v.reshape(-1)])


def _cell_seam_labels(cc: CodeCensus, ap: AxisProblem, zg: np.ndarray,
                      Asup, Bsup, u: int, v: int, swap: bool, ta: int,
                      tb: int, sign: int) -> list[int]:
    """phi-labels (identity census) of the kappa basis seams of cell
    (u, v, swap, ta, tb).  zg = kernel basis grids OF THE CELL's pair."""
    A2, B2 = ap.cell_polys(Asup, Bsup, ta, tb)
    kappa = zg.shape[0]
    SU = ap.seam_grids(A2, zg).reshape(kappa, -1)
    SV = ap.seam_grids(B2, zg).reshape(kappa, -1)
    uinv = pow(u, -1, ap.ell)
    vinv = pow(v, -1, ap.m)
    gens = []
    for i in range(kappa):
        t = np.concatenate([SU[i], SV[i]]).astype(np.uint8)
        t = _transport_to_identity(t, ap, uinv, vinv, swap, ta, tb, sign)
        gens.append(cc.phi(t))
    return gens


def _span(gens: list[int]) -> set[int]:
    span = {0}
    for g in gens:
        span |= {s ^ g for s in span}
    span.discard(0)
    return span


def verify_transport(point: str, axis: str, pres_i: int, W: int,
                     threads: int) -> int:
    """Fix the translation sign and validate the whole transport chain:
    for sampled cells, hit-weight multisets (span cap lights) computed
    via (a) identity census + transport vs (b) the cell's OWN census
    must agree.  Returns the validated sign (+1 or -1); raises if
    neither validates."""
    spec = POINTS[point]
    orders, d = spec["orders"], spec["d"]
    a_s, b_s = spec["pres"][pres_i]
    ap = AxisProblem(orders, parse_poly(a_s), parse_poly(b_s),
                     0 if axis == "x" else 1, d)
    A0, B0 = ap.variant_supports(1, 1, False)
    cc0 = CodeCensus((ap.ell, ap.m), A0, B0, d)
    lights0 = cc0.light_census(W, threads=threads)
    us, vs = units(ap.ell), units(ap.m)
    rng = np.random.default_rng(7)
    cells = []
    for _ in range(4):
        cells.append((int(rng.choice(us)), int(rng.choice(vs)),
                      bool(rng.integers(2)), int(rng.integers(ap.ell)),
                      int(rng.integers(ap.ell))))
    cells.append((1, 1, False, 1, 2))  # pure-translation probe
    ok_sign = None
    for sign in (-1, 1):
        all_ok = True
        for (u, v, swap, ta, tb) in cells:
            Asup, Bsup = ap.variant_supports(u, v, swap)
            from bb_lab.fibering import kernel_basis
            Ap = Poly.from_support([tuple(t) for t in Asup], cc0.G)
            Bp = Poly.from_support([tuple(t) for t in Bsup], cc0.G)
            zg = kernel_basis(Ap, Bp).reshape(-1, ap.ell, ap.m)
            gens = _cell_seam_labels(cc0, ap, zg, Asup, Bsup, u, v,
                                     swap, ta, tb, sign)
            hits_a = sorted(lights0[h]
                            for h in _span(gens) & set(lights0))
            # ground truth: the cell's own census
            A2, B2 = ap.cell_polys(Asup, Bsup, ta, tb)
            cc_cell = CodeCensus((ap.ell, ap.m), A2, B2, d)
            lights_c = cc_cell.light_census(W, threads=threads)
            gens_c = [cc_cell.phi(np.concatenate([su, sv]).astype(np.uint8))
                      for su, sv in zip(
                          ap.seam_grids(A2, zg).reshape(zg.shape[0], -1),
                          ap.seam_grids(B2, zg).reshape(zg.shape[0], -1))]
            hits_b = sorted(lights_c[h]
                            for h in _span(gens_c) & set(lights_c))
            if hits_a != hits_b:
                all_ok = False
                break
        if all_ok:
            ok_sign = sign
            break
    if ok_sign is None:
        raise AssertionError("transport chain validates under NEITHER "
                             "sign — convention bug, do not scan")
    print(f"transport VALIDATED (sign {ok_sign:+d}) on {len(cells)} "
          "sampled cells incl. swap/unit/translation", flush=True)
    return ok_sign


def orbit_scan(point: str, axis: str, pres_i: int, W: int,
               threads: int, full_orbit_out: Path | None) -> None:
    """Exact census-filter scan of the whole v1 orbit of one
    presentation: ONE light census (identity variant), then pure
    linear algebra per cell via the validated transport."""
    spec = POINTS[point]
    orders, d = spec["orders"], spec["d"]
    a_s, b_s = spec["pres"][pres_i]
    ap = AxisProblem(orders, parse_poly(a_s), parse_poly(b_s),
                     0 if axis == "x" else 1, d)
    t0 = time.time()
    sign = verify_transport(point, axis, pres_i, W, threads)
    A0, B0 = ap.variant_supports(1, 1, False)
    cc0 = CodeCensus((ap.ell, ap.m), A0, B0, d)
    lights = cc0.light_census(W, threads=threads)
    hist = Counter(lights.values())
    light_set = set(lights)
    print(f"{point} p{pres_i} {axis}: census W={W} -> {len(lights)} "
          f"light classes, weights {dict(sorted(hist.items()))} "
          f"({time.time() - t0:.0f}s)", flush=True)
    us, vs = units(ap.ell), units(ap.m)
    n_cells = n_pass = 0
    survivors = []
    seen: set = set()
    from bb_lab.fibering import kernel_basis
    for u, v, swap in itertools.product(us, vs, (False, True)):
        Asup, Bsup = ap.variant_supports(u, v, swap)
        Ap = Poly.from_support([tuple(t) for t in Asup], cc0.G)
        Bp = Poly.from_support([tuple(t) for t in Bsup], cc0.G)
        zg = kernel_basis(Ap, Bp).reshape(-1, ap.ell, ap.m)
        for ta in range(ap.ell):
            for tb in range(ap.ell):
                A2, B2 = ap.cell_polys(Asup, Bsup, ta, tb)
                key = (tuple(A2), tuple(B2))
                if key in seen:
                    continue
                seen.add(key)
                n_cells += 1
                gens = _cell_seam_labels(cc0, ap, zg, Asup, Bsup, u, v,
                                         swap, ta, tb, sign)
                if not (_span(gens) & light_set):
                    n_pass += 1
                    survivors.append({
                        "point": point, "pres": pres_i, "axis": axis,
                        "u": u, "v": v, "swap": swap, "ta": ta,
                        "tb": tb, "A": poly_str(A2), "B": poly_str(B2),
                        "census_W": W})
    print(f"{point} p{pres_i} {axis}: {n_cells} unique cells, "
          f"census-W={W} filter pass = {n_pass} "
          f"({time.time() - t0:.0f}s)", flush=True)
    if survivors and full_orbit_out:
        full_orbit_out.write_text(
            "\n".join(json.dumps(s) for s in survivors))
        print(f"  -> {len(survivors)} survivors -> {full_orbit_out}",
              flush=True)


def scan_code(orders: tuple[int, int], A_sup, B_sup, d: int, axis: str,
              W: int, sign: int, threads: int,
              max_survivors: int = 12) -> dict:
    """Census + exact whole-orbit scan for one code presentation, with
    a pre-validated transport sign.  Returns summary + top survivors."""
    ap = AxisProblem(orders, A_sup, B_sup, 0 if axis == "x" else 1, d)
    A0, B0 = ap.variant_supports(1, 1, False)
    t0 = time.time()
    cc0 = CodeCensus((ap.ell, ap.m), A0, B0, d)
    lights = cc0.light_census(W, threads=threads)
    light_set = set(lights)
    hist = Counter(lights.values())
    us, vs = units(ap.ell), units(ap.m)
    n_cells = n_pass = 0
    survivors = []
    seen: set = set()
    from bb_lab.fibering import kernel_basis
    for u, v, swap in itertools.product(us, vs, (False, True)):
        Asup, Bsup = ap.variant_supports(u, v, swap)
        Ap = Poly.from_support([tuple(t) for t in Asup], cc0.G)
        Bp = Poly.from_support([tuple(t) for t in Bsup], cc0.G)
        zg = kernel_basis(Ap, Bp).reshape(-1, ap.ell, ap.m)
        for ta in range(ap.ell):
            for tb in range(ap.ell):
                A2, B2 = ap.cell_polys(Asup, Bsup, ta, tb)
                key = (tuple(A2), tuple(B2))
                if key in seen:
                    continue
                seen.add(key)
                n_cells += 1
                gens = _cell_seam_labels(cc0, ap, zg, Asup, Bsup, u, v,
                                         swap, ta, tb, sign)
                span = _span(gens)
                if len(span) != (1 << len(gens)) - 1:
                    # degenerate im Delta (dependent/zero seam classes):
                    # the (R) doubling structure is broken here — record
                    # as non-pass, never as a survivor
                    continue
                if not (span & light_set):
                    n_pass += 1
                    if len(survivors) < max_survivors:
                        survivors.append(
                            {"A": poly_str(A2), "B": poly_str(B2),
                             **ap.cover_spec(A2, B2)})
    lbasis: list[int] = []
    for lab in lights:
        cur = lab
        for b in lbasis:
            cur = min(cur, cur ^ b)
        if cur:
            lbasis.append(cur)
    return {"n_lights": len(lights),
            "light_hist": dict(sorted(hist.items())),
            "light_span_dim": len(lbasis), "k": cc0.k,
            "n_cells": n_cells, "n_pass": n_pass,
            "survivors": survivors,
            "wall_s": round(time.time() - t0, 1)}


def triage(input_path: Path, axis: str, sign: int, threads: int,
           w_override: int | None, only_group: str | None) -> None:
    """Batch code-level triage over a JSONL of presentations."""
    out_path = DATA_DIR / f"triage_{axis}_results.jsonl"
    done: set = set()
    if out_path.exists():
        for line in out_path.read_text().splitlines():
            r = json.loads(line)
            done.add((r["group"], r["A"], r["B"], r["axis"]))
    rows = [json.loads(x)
            for x in input_path.read_text().splitlines()]
    if only_group:
        rows = [r for r in rows if r["group"] == only_group]
    print(f"triage {axis}: {len(rows)} rows ({len(done)} already done)",
          flush=True)
    with out_path.open("a") as fh:
        for i, r in enumerate(rows):
            if (r["group"], r["A"], r["B"], axis) in done:
                continue
            W = w_override or (2 * r["d"] - 2)
            try:
                res = scan_code(tuple(r["orders"]), parse_poly(r["A"]),
                                parse_poly(r["B"]), r["d"], axis, W,
                                sign, threads)
            except Exception as exc:
                res = {"error": str(exc)[:200]}
            rec = {**r, "axis": axis, "W": W, **res}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            print(f"[{i + 1}/{len(rows)}] {r['group']} "
                  f"A={r['A'][:24]}.. lights={res.get('n_lights')} "
                  f"hist={res.get('light_hist')} "
                  f"pass={res.get('n_pass')}/{res.get('n_cells')} "
                  f"({res.get('wall_s')}s)", flush=True)


def one_census(orders_s: str, a_s: str, b_s: str, axis: str, d: int,
               W: int, threads: int) -> None:
    orders = tuple(int(t) for t in orders_s.split(","))
    ap = AxisProblem(orders, parse_poly(a_s), parse_poly(b_s),
                     0 if axis == "x" else 1, d)
    Asup, Bsup = ap.variant_supports(1, 1, False)
    cc = CodeCensus((ap.ell, ap.m), Asup, Bsup, d)
    t0 = time.time()
    lights = cc.light_census(W, threads=threads)
    hist = Counter(lights.values())
    # span dimension of the light labels
    basis: list[int] = []
    for lab in lights:
        cur = lab
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
    print(f"Z{ap.ell}xZ{ap.m} k={cc.k} d={d}: census W={W} -> "
          f"{len(lights)} light classes, weight histogram "
          f"{dict(sorted(hist.items()))}, span dim {len(basis)} "
          f"({time.time() - t0:.0f}s)", flush=True)
    # cell verdict for THIS presentation cell
    span = cc.imdelta_span(ap, Asup, Bsup)
    hit = sorted(set(span) & set(lights.keys()))
    hitw = [lights[h] for h in hit]
    print(f"  this cell: |im Delta| = {len(span)}, light hits {len(hit)}"
          f"{' (weights ' + str(sorted(hitw)) + ')' if hit else ''} -> "
          f"{'FAIL' if hit else 'PASS'} at census W", flush=True)


def controls(threads: int) -> None:
    print("== census control 1: pair72 base [[36,4,4]] Z3xZ6 x "
          "(Lean-proven doubling cover; cell must PASS census W=6) ==",
          flush=True)
    one_census("3,6", "x^2 + y + y^3", "1 + x + y^2", "x", 4, 6,
               threads)
    print("== census control 2: T1 pres0 stored (certify-refuted; cell "
          "must FAIL, hit weight <= 14... stored cell freeze) ==",
          flush=True)
    p = POINTS["T1"]["pres"][0]
    one_census("7,9", p[0], p[1], "x", 10, 14, threads)
    print("== census control 3: the T1.5-refuted c2 cell (min 14) — "
          "must FAIL census W=14 with a hit at 14 ==", flush=True)
    one_census("7,9", "y^2 + x^6 + x^6*y^4", "y^6 + x^5*y^4 + x^6", "x",
               10, 14, threads)


def main() -> None:
    apar = argparse.ArgumentParser()
    apar.add_argument("cmd", choices=["census", "orbit", "controls",
                                      "verify", "triage"])
    apar.add_argument("--orders", type=str)
    apar.add_argument("--A", type=str)
    apar.add_argument("--B", type=str)
    apar.add_argument("--axis", type=str, default="x")
    apar.add_argument("--d", type=int, default=10)
    apar.add_argument("--W", type=int, default=14)
    apar.add_argument("--point", type=str)
    apar.add_argument("--pres", type=int, default=0)
    apar.add_argument("--threads", type=int, default=4)
    apar.add_argument("--out", type=str, default=None)
    apar.add_argument("--sign", type=int, default=None,
                      help="validated transport sign (from `verify`)")
    apar.add_argument("--input", type=str,
                      default=str(DATA_DIR / "triage_input.jsonl"))
    apar.add_argument("--only-group", type=str, default=None)
    apar.add_argument("--w-override", type=int, default=None)
    args = apar.parse_args()
    if args.cmd == "controls":
        controls(args.threads)
    elif args.cmd == "census":
        one_census(args.orders, args.A, args.B, args.axis, args.d,
                   args.W, args.threads)
    elif args.cmd == "verify":
        verify_transport(args.point, args.axis, args.pres, args.W,
                         args.threads)
    elif args.cmd == "triage":
        if args.sign is None:
            raise SystemExit("triage requires --sign (run `verify` "
                             "first on a cheap point, both axes)")
        triage(Path(args.input), args.axis, args.sign, args.threads,
               args.w_override, args.only_group)
    else:
        out = (Path(args.out) if args.out else
               DATA_DIR / f"{args.point}_p{args.pres}_{args.axis}"
                          f"_censusscan.jsonl")
        orbit_scan(args.point, args.axis, args.pres, args.W,
                   args.threads, out)


if __name__ == "__main__":
    main()
