#!/usr/bin/env python3
"""A42 S3d — the register kernels at 2-adic depth a = 3 (the p = 24
scale): does the two-branch jet engine stay class-complete?

Mirrors s2_registers.run_a2_full at a = 3 (Lambda_3 =
F2[y]/((y^2+y+1)^8), dim 16): Teichmueller zeta by fixpoint,
u = (1+pi)^{1/3} by brute search over 1 + m (the cube map is
bijective there), the two full-depth A~-roots xt = zeta^2 u (tangent
branch) and zeta u (transverse), invariance of the block-1 register
on trivials, and the exact kernels of R(xt0), R(xt1) on the 63
nonzero H-classes.  Joint kernel EMPTY => the p = 24 two-branch
racer (once a 2-word state engine exists) is class-complete — the
r = 4 member period's register soundness, banked ahead of need.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib.util


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SR = _load("a42_s2_registers")
SY = _load("a42_s1_syzygy")
DATA = LAB / "data" / "a42"
OMEGA2 = 0b111


def main(Wx: int = 10):
    print("== a = 3: full-depth Lambda_3 registers (p = 24 scale) ==",
          flush=True)
    ow = SY.OmegaWindow(3, Wx)
    L3 = ow.L
    red, mul = L3.red, L3.mul
    assert L3.dim == 16
    pi = red(OMEGA2)
    # Teichmueller zeta: delta = pi + delta^2 fixpoint
    delta = 0
    for _ in range(10):
        delta = red(pi ^ mul(delta, delta))
    zeta = red(0b10 ^ delta)
    assert red(mul(zeta, zeta) ^ zeta ^ 1) == 0, "zeta"
    # u^3 = 1 + pi with u = 1 + m  (brute; the cube map is bijective
    # on 1 + m)
    target = red(1 ^ pi)
    u = None
    for cand in range(1 << L3.dim):
        if mul(mul(cand, cand), cand) == target:
            u = cand
            break
    assert u is not None, "no cube root of 1 + pi"
    xt0 = mul(mul(zeta, zeta), u)      # tangent branch
    xt1 = mul(zeta, u)                 # transverse branch
    for xt in (xt0, xt1):
        assert red(mul(mul(xt, xt), xt) ^ 1 ^ pi) == 0, "A~ root"
    print(f"  constants: zeta={zeta:#x} u={u:#x} xt0={xt0:#x} "
          f"xt1={xt1:#x} (A~-roots verified to full depth)",
          flush=True)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    print(f"  class space: dim H = {ncls} (Theorem H predicts 6)",
          flush=True)
    assert ncls == 6
    reps = SR.class_reps(ow, Z, Lrows, 64)
    for t in T[:40]:
        assert SR.R_full(ow, t, 1, xt0) == 0
        assert SR.R_full(ow, t, 1, xt1) == 0
    print("  invariance on trivials: block-1 register exact at both "
          "roots (40 checks)", flush=True)
    ker0 = [m for m, vv in reps.items()
            if SR.R_full(ow, vv, 1, xt0) == 0]
    ker1 = [m for m, vv in reps.items()
            if SR.R_full(ow, vv, 1, xt1) == 0]
    kerj = sorted(set(ker0) & set(ker1))
    print(f"  ker R(xt0) [tangent]: n={len(ker0)} {sorted(ker0)}",
          flush=True)
    print(f"  ker R(xt1) [transverse]: n={len(ker1)} {sorted(ker1)}",
          flush=True)
    verdict = ("EMPTY — the two-run engine is class-complete at a=3"
               if not kerj else "NONEMPTY")
    print(f"  joint kernel: {kerj} ({verdict})", flush=True)
    out = dict(a=3, Wx=Wx, ncls=ncls, zeta=zeta, u=u,
               xt0=xt0, xt1=xt1,
               ker_xt0=sorted(ker0), ker_xt1=sorted(ker1),
               joint=kerj)
    (DATA / "s3d_registers_a3.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s3d_registers_a3.json'}", flush=True)


if __name__ == "__main__":
    main()
