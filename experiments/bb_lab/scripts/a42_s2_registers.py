#!/usr/bin/env python3
"""A42 S2d — exact kernels of the racer's class registers on H, and
the two-register repair.

The S1h racer detects nontriviality through R_A = the pair of RESIDUE
evaluations of the input-block omega-content at the two x-points of
gcd(Abar_w, Bbar_w) (y -> omega kills pi).  At 2-adic depth a >= 1
this functional provably kills pi-divisible-representable classes; the
question with certificate consequences is its EXACT kernel on
H = F4(transverse) + F4[pi]/pi^2(tangency), and whether the repair
register R_B = the TANGENT-BRANCH pi-JET (evaluation over
W = F2[y]/((y^2+y+1)^2) = F4[pi]/pi^2 at the deformed variety point
x^0 = zeta^2(1+pi), an exact root of A~_w = (x^3+1) + pi over W)
covers the blind spot: expected ker R_A = pi-socle-representable
classes, ker R_B = the transverse summand, ker R_A ^ ker R_B = 0
(two runs see everything).

Computed exactly at a = 1 (p = 6 scale) and a = 2 (p = 12 scale) on
the OmegaWindow class space: per-class values of R_A (both x-points,
both blocks), R_B, pi-structure (im pi / ker pi), and the h-DP
min-weight from s1_classprofile.json.  Functional invariance under
trivial syzygies is asserted mechanically, as are the Teichmueller
constants (zeta = y^2+1, zeta^3 = 1, A~_w(x^0) = 0 over W).
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

import a42_lib as AL  # noqa: E402
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_syzygy", Path(__file__).parent / "a42_s1_syzygy.py")
SY = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SY)

DATA = LAB / "data" / "a42"

OMEGA2 = 0b111
W_MOD = AL.pmul(OMEGA2, OMEGA2)          # (y^2+y+1)^2 = y^4+y^2+1


def wred(z):
    return AL.pmod(z, W_MOD)


def wmul(x, y):
    return wred(AL.pmul(x, y))


def f4red(z):
    return AL.pmod(z, OMEGA2)


def f4mul(x, y):
    return f4red(AL.pmul(x, y))


# Teichmueller constants in W (asserted below)
ZETA = 0b101            # y^2 + 1
XHAT = 0b1101           # zeta^2 (1+pi) = y^3+y^2+1
XHATI = 0b1011          # its inverse


def sanity_W():
    assert wred(AL.pmul(ZETA, ZETA) ^ ZETA ^ 1) == 0, "zeta^2+zeta+1"
    z3 = wmul(wmul(ZETA, ZETA), ZETA)
    assert z3 == 1, ("zeta^3", z3)
    pi = wred(OMEGA2)
    x3 = wmul(wmul(XHAT, XHAT), XHAT)
    assert x3 == wred(1 ^ pi), ("xhat^3 != 1+pi", x3)   # A~_w(xhat)=0
    assert wmul(XHAT, XHATI) == 1
    print("W constants: zeta=y^2+1, xhat=y^3+y^2+1 verified "
          "(zeta^3=1, xhat^3=1+pi, xhat*xhati=1)", flush=True)


def content(ow, v, blk, c):
    lam = 0
    for i in range(ow.L.dim):
        if v[ow.bit(blk, c, i)]:
            lam |= 1 << i
    return int(lam)


def R_A(ow, v, blk):
    """Residue-evaluation register pair: for x0 in {w, w^2} (F4),
    sum_c x0^{-c} * (content mod OMEGA2 evaluated at y -> w)."""
    outs = []
    for x0 in (0b10, 0b11):                    # w, w^2 as F4 polys
        x0i = f4mul(x0, x0)                    # x0^-1 = x0^2
        acc = 0
        pw = 1
        for c in range(ow.Wx):
            lam = content(ow, v, blk, c)
            r = f4red(lam)
            # evaluate r = r0 + r1 y at y = w  (w == the class of y)
            acc ^= f4mul(pw, r)
            pw = f4mul(pw, x0i)
        outs.append(acc)
    return tuple(outs)


def R_B(ow, v, blk):
    """Tangent-branch jet register: sum_c xhat^{+c} * (content mod
    (y^2+y+1)^2) over W.  The +c weighting evaluates the generating
    function AT the shared tangent point xhat (both A~_w and B^_w
    vanish there over W — the scheme tangency); the -c weighting
    would evaluate at xhat^{-1} = zeta(1+pi), which is an A~-root on
    the TRANSVERSE branch and kills only block-1 boundaries (found
    mechanically: the first R_B draft failed invariance on block 0)."""
    acc = 0
    pw = 1
    for c in range(ow.Wx):
        lam = wred(content(ow, v, blk, c))
        acc ^= wmul(pw, lam)
        pw = wmul(pw, XHAT)
    return int(acc)


def pi_map(ow):
    """Matrix of content-wise multiplication by pi on window bits."""
    n = ow.nbits
    P = np.zeros((n, n), dtype=np.uint8)
    for blk in (0, 1):
        for c in range(ow.Wx):
            for i in range(ow.L.dim):
                z = ow.L.mul(ow.L.pi, 1 << i)
                for oi in range(ow.L.dim):
                    if (z >> oi) & 1:
                        P[ow.bit(blk, c, oi), ow.bit(blk, c, i)] = 1
    return P


def run(a: int, Wx: int, profile_key: str):
    print(f"== a = {a} (p = {3 * (1 << a)} scale) ==", flush=True)
    ow = SY.OmegaWindow(a, Wx)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    assert ncls == 6, ncls
    P = pi_map(ow)
    # class-rep table: greedy basis of labels
    reps = {}
    basis = []
    for z in Z:
        lab = int("".join(map(str, (Lrows @ z) % 2)), 2)
        if lab == 0:
            continue
        # try to extend basis
        cur = lab
        vv = z.copy()
        for (bl, bv) in basis:
            top = 1 << (bl.bit_length() - 1)
            if cur & top:
                cur ^= bl
                vv = (vv + bv) % 2
        if cur:
            basis.append((cur, vv))
        if len(basis) == 6:
            break
    assert len(basis) == 6, len(basis)
    for m in range(1, 64):
        acc = np.zeros(ow.nbits, dtype=np.uint8)
        lab = 0
        for i, (bl, bv) in enumerate(basis):
            pass
        # combine greedily: reduce m over basis labels
        cur = 0
        vv = np.zeros(ow.nbits, dtype=np.uint8)
        want = m
        # solve for combo with label m by Gaussian elim over labels
        labs = [bl for bl, _ in basis]
        vecs = [bv for _, bv in basis]
        x = want
        chosen = []
        for bl, bv in sorted(zip(labs, vecs),
                             key=lambda t: -t[0].bit_length()):
            top = 1 << (bl.bit_length() - 1)
            if x & top:
                x ^= bl
                chosen.append(bv)
        assert x == 0, (m, "labels do not span")
        for bv in chosen:
            vv = (vv + bv) % 2
        got = int("".join(map(str, (Lrows @ vv) % 2)), 2)
        assert got == m, (m, got)
        reps[m] = vv

    # invariance of the functionals under trivial syzygies
    nchk = 0
    for t in T[:40]:
        for blk in (0, 1):
            assert R_A(ow, t, blk) == (0, 0), (blk, "R_A not invariant")
            assert R_B(ow, t, blk) == 0, (blk, "R_B not invariant")
        nchk += 1
    print(f"  invariance: R_A, R_B vanish on {nchk} trivial-syzygy "
          f"basis rows, both blocks", flush=True)

    prof = json.load(open(DATA / "s1_classprofile.json"))[profile_key]
    pw = prof["profile"]

    impi = set()
    kerpi = set()
    for m, vv in reps.items():
        pv = (P @ vv) % 2
        # pi * syzygy is a syzygy
        assert not ((H @ pv) % 2).any()
        lab = int("".join(map(str, (Lrows @ pv) % 2)), 2)
        if lab == 0:
            kerpi.add(m)
        else:
            impi.add(lab)

    rows = []
    kerA = {"0": set(), "1": set()}
    kerB = {"0": set(), "1": set()}
    for m, vv in reps.items():
        ra0 = R_A(ow, vv, 0)
        ra1 = R_A(ow, vv, 1)
        rb0 = R_B(ow, vv, 0)
        rb1 = R_B(ow, vv, 1)
        for blk, ra, rb in ((0, ra0, rb0), (1, ra1, rb1)):
            if ra == (0, 0):
                kerA[str(blk)].add(m)
            if rb == 0:
                kerB[str(blk)].add(m)
        rows.append(dict(label=m, w=pw[str(m)], impi=m in impi,
                         kerpi=m in kerpi, RA_blk0=list(ra0),
                         RA_blk1=list(ra1), RB_blk0=rb0, RB_blk1=rb1))

    out = dict(a=a, impi=sorted(impi), kerpi=sorted(kerpi))
    for blk in ("0", "1"):
        kA, kB = kerA[blk], kerB[blk]
        both = kA & kB
        wA = min((pw[str(m)] for m in kA), default=None)
        wB = min((pw[str(m)] for m in kB), default=None)
        print(f"  block {blk}: |ker R_A| = {len(kA)} {sorted(kA)}; "
              f"min h-DP weight over ker R_A = {wA}", flush=True)
        print(f"            |ker R_B| = {len(kB)} {sorted(kB)}; "
              f"min weight over ker R_B = {wB}", flush=True)
        print(f"            ker R_A ^ ker R_B = {sorted(both)} "
              f"({'EMPTY - two runs cover H' if not both else 'NONEMPTY'})",
              flush=True)
        out[f"kerA_blk{blk}"] = sorted(kA)
        out[f"kerB_blk{blk}"] = sorted(kB)
        out[f"joint_blk{blk}"] = sorted(both)
        out[f"minw_kerA_blk{blk}"] = wA
        out[f"minw_kerB_blk{blk}"] = wB
    print(f"  im(pi) labels: {sorted(impi)}  (n={len(impi)}); "
          f"ker(pi) labels n={len(kerpi)}", flush=True)
    out["rows"] = rows
    return out


def class_reps(ow, Z, Lrows, nlab):
    """label -> window rep, via a greedy label-basis."""
    basis = []
    for z in Z:
        lab = int("".join(map(str, (Lrows @ z) % 2)), 2)
        if lab == 0:
            continue
        cur, vv = lab, z.copy()
        for (bl, bv) in basis:
            top = 1 << (bl.bit_length() - 1)
            if cur & top:
                cur ^= bl
                vv = (vv + bv) % 2
        if cur:
            basis.append((cur, vv))
        if len(basis) == Lrows.shape[0]:
            break
    assert len(basis) == Lrows.shape[0]
    reps = {}
    for m in range(1, nlab):
        x = m
        vv = np.zeros(ow.nbits, dtype=np.uint8)
        for bl, bv in sorted(basis, key=lambda t: -t[0].bit_length()):
            top = 1 << (bl.bit_length() - 1)
            if x & top:
                x ^= bl
                vv = (vv + bv) % 2
        assert x == 0
        got = int("".join(map(str, (Lrows @ vv) % 2)), 2)
        assert got == m
        reps[m] = vv
    return reps


def run_a0(Wx: int = 8):
    """a = 0 (p = 9 scale): is the residue register COMPLETE on
    H = F4^2?  (The floor(9) certificate rests on it.)"""
    print("== a = 0 (p = 9 scale): R_A completeness ==", flush=True)
    ow = SY.OmegaWindow(0, Wx)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    assert ncls == 4, ncls
    reps = class_reps(ow, Z, Lrows, 16)
    for t in T[:30]:
        for blk in (0, 1):
            assert R_A(ow, t, blk) == (0, 0)
    ker = {"0": [], "1": []}
    for m, vv in reps.items():
        for blk in (0, 1):
            if R_A(ow, vv, blk) == (0, 0):
                ker[str(blk)].append(m)
    for blk in ("0", "1"):
        print(f"  block {blk}: ker R_A on the 15 nonzero classes = "
              f"{ker[blk]} "
              f"({'COMPLETE - certificate stands' if not ker[blk] else 'INCOMPLETE'})",
              flush=True)
    return {"kerA_blk0": ker["0"], "kerA_blk1": ker["1"]}


# ---------------- full-depth Lambda_2 registers (a = 2) -------------
def lam2_constants(L2):
    """Teichmueller zeta and the two full-depth A~-roots in Lambda_2."""
    red, mul = L2.red, L2.mul
    pi = red(OMEGA2)
    # zeta = y + delta, delta = pi + delta^2 (fixpoint)
    delta = 0
    for _ in range(6):
        delta = red(pi ^ mul(delta, delta))
    zeta = red(0b10 ^ delta)
    assert red(mul(zeta, zeta) ^ zeta ^ 1) == 0, "zeta"
    # u = (1+pi)^{1/3} = 1 + pi + pi^2 + pi^3
    pi2 = mul(pi, pi)
    pi3 = mul(pi2, pi)
    u = red(1 ^ pi ^ pi2 ^ pi3)
    u3 = mul(mul(u, u), u)
    assert u3 == red(1 ^ pi), ("u^3", u3)
    z2 = mul(zeta, zeta)
    xt0 = mul(z2, u)          # tangent-branch full root
    xt1 = mul(zeta, u)        # transverse-branch full root
    for xt in (xt0, xt1):
        x3 = mul(mul(xt, xt), xt)
        # A~_w(x) = (x^3 + 1) + pi = 0
        assert red(x3 ^ 1 ^ pi) == 0, ("A~ root", xt)
    return zeta, xt0, xt1


def R_full(ow, v, blk, point):
    """Full Lambda_a register: sum_c point^c * content(blk, c),
    computed in Lambda_a (the +c convention evaluates at `point`;
    invariant on block-1 reads because boundaries there are A~ t and
    point is an exact A~-root)."""
    L = ow.L
    acc = 0
    pw = 1
    for c in range(ow.Wx):
        lam = content(ow, v, blk, c)
        acc ^= L.mul(pw, lam)
        pw = L.mul(pw, point)
    return int(acc)


def run_a2_full(Wx: int = 10):
    print("== a = 2: full-depth Lambda_2 registers ==", flush=True)
    ow = SY.OmegaWindow(2, Wx)
    L2 = ow.L
    zeta, xt0, xt1 = lam2_constants(L2)
    print(f"  constants: zeta={zeta:#x} xt0={xt0:#x} xt1={xt1:#x} "
          "(A~-roots verified to full depth)", flush=True)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    assert ncls == 6
    reps = class_reps(ow, Z, Lrows, 64)
    # invariance: block 1 must be exact; block 0 expected to fail
    inv = {"0": True, "1": True}
    for t in T[:40]:
        for blk in (0, 1):
            for pt in (xt0, xt1):
                if R_full(ow, t, blk, pt) != 0:
                    inv[str(blk)] = False
    print(f"  invariance on trivials: block0={inv['0']} (expected "
          f"False), block1={inv['1']} (required True)", flush=True)
    assert inv["1"], "block-1 full register not invariant!"
    prof = json.load(open(DATA / "s1_classprofile.json"))["12"]
    pw = prof["profile"]
    ker0, ker1, kerj = [], [], []
    for m, vv in reps.items():
        r0 = R_full(ow, vv, 1, xt0)
        r1 = R_full(ow, vv, 1, xt1)
        if r0 == 0:
            ker0.append(m)
        if r1 == 0:
            ker1.append(m)
        if r0 == 0 and r1 == 0:
            kerj.append(m)
    w0 = min((pw[str(m)] for m in ker0), default=None)
    w1 = min((pw[str(m)] for m in ker1), default=None)
    print(f"  ker R(xt0) [tangent]: n={len(ker0)} {ker0} minw={w0}",
          flush=True)
    print(f"  ker R(xt1) [transverse]: n={len(ker1)} {ker1} minw={w1}",
          flush=True)
    print(f"  joint kernel: {kerj} "
          f"({'EMPTY - two-run engine is class-complete' if not kerj else 'NONEMPTY'})",
          flush=True)
    return dict(ker_xt0=ker0, ker_xt1=ker1, joint=kerj,
                minw_xt0=w0, minw_xt1=w1,
                zeta=zeta, xt0=xt0, xt1=xt1)


def run_a1_full(Wx: int = 10):
    """Same full-register analysis at a = 1 (p = 6 scale) for the
    control ladder: here Lambda_1 = W and the full register at the
    branch roots should already separate what R_A+R_B could not."""
    print("== a = 1: full-depth Lambda_1 registers ==", flush=True)
    ow = SY.OmegaWindow(1, Wx)
    L1 = ow.L
    red, mul = L1.red, L1.mul
    pi = red(OMEGA2)
    delta = 0
    for _ in range(5):
        delta = red(pi ^ mul(delta, delta))
    zeta = red(0b10 ^ delta)
    assert red(mul(zeta, zeta) ^ zeta ^ 1) == 0
    u = red(1 ^ pi)                       # (1+pi)^{1/3} = 1+pi here
    assert mul(mul(u, u), u) == red(1 ^ pi)
    xt0 = mul(mul(zeta, zeta), u)
    xt1 = mul(zeta, u)
    for xt in (xt0, xt1):
        assert red(mul(mul(xt, xt), xt) ^ 1 ^ pi) == 0
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    reps = class_reps(ow, Z, Lrows, 64)
    for t in T[:40]:
        assert R_full(ow, t, 1, xt0) == 0
        assert R_full(ow, t, 1, xt1) == 0
    prof = json.load(open(DATA / "s1_classprofile.json"))["6"]
    pw = prof["profile"]
    ker0 = [m for m, vv in reps.items() if R_full(ow, vv, 1, xt0) == 0]
    ker1 = [m for m, vv in reps.items() if R_full(ow, vv, 1, xt1) == 0]
    kerj = sorted(set(ker0) & set(ker1))
    print(f"  ker R(xt0): n={len(ker0)} minw="
          f"{min((pw[str(m)] for m in ker0), default=None)}", flush=True)
    print(f"  ker R(xt1): n={len(ker1)} minw="
          f"{min((pw[str(m)] for m in ker1), default=None)}", flush=True)
    print(f"  joint kernel: {kerj}", flush=True)
    return dict(ker_xt0=ker0, ker_xt1=ker1, joint=kerj)


def main():
    t0 = time.time()
    sanity_W()
    out = {}
    out["a0"] = run_a0()
    out["a1"] = run(1, 10, "6")
    out["a1_full"] = run_a1_full()
    out["a2"] = run(2, 10, "12")
    out["a2_full"] = run_a2_full()
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s2_registers.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_registers.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
