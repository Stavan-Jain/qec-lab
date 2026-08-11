"""A32 M4 bridge-layer verification (falsify-first, run BEFORE emitting).

Verifies, in numpy/python, every fact the M4 chain-plumbing layer will
`decide`/`native_decide` in Lean (QECLean `Codes/Mitten/M150/`):

  1. the op→table bridge basis facts (maskOf (op δ_i) = tbl[i]) for all
     five blocks of all four triple instances, straight from the
     d2term/cmTerm entry semantics (independent re-derivation, not the
     emitter's Lm/Rm shortcut);
  2. the composite-table identities behind the derived-form lemmas
     (comp-id, inverse-composition, mutual-inverse pairs);
  3. the coset machinery: correction columns land in Ann, Ann/AnnT are
     xor-closed, left-null filters vanish on the singular map's columns;
  4. parity feeders: all raw operator columns have odd weight;
  5. split coverage: every even (p,q,r) with p+q+r ≤ 8 has an assigned
     mode from the instance's mode set;
  6. no classified triple has t = 0 (the mixed zero/listed join case is
     vacuous);
  7. rowsXpk/rowsZpk are exactly the packed H_Z/H_X row chains in the
     side's join layout (recomputed from cmTerm/d2term semantics);
  8. the NEW tables to be emitted (tRawAX, tCinvX0/1, tRawCZ0/1, tAinvZ)
     against their defining properties.

Exits nonzero on any failure.  Run (from experiments/bb_lab):
  uv run python scripts/a32_m4_bridge_check.py
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("gen", HERE / "m150_gen_lean_data.py")
GEN = importlib.util.module_from_spec(spec)
spec.loader.exec_module(GEN)

_msk, _cols, _xf, _pop = GEN._msk, GEN._cols, GEN._xf, GEN._pop
_inv2, _pseudo_inv = GEN._inv2, GEN._pseudo_inv

QECLEAN = Path(
    "/Users/stavanjain/Code/QECLean/.claude/worktrees/a32-mitten150"
)
FLOORDATA = QECLEAN / "QEC/Stabilizer/Codes/Mitten/M150/FloorData.lean"

fails = []


def check(name, ok):
    tag = "ok " if ok else "FAIL"
    print(f"[{tag}] {name}")
    if not ok:
        fails.append(name)


# ---------------------------------------------------------------- setup
a26 = GEN._load_a26()
G = a26.Group.from_file(GEN.GROUPS / "group_30_1.txt")
mul, inv = G.mul, G.inv
n = 30
S = {k: set(v) for k, v in GEN.SETS.items()}

# entry semantics as in QECLean StabilizerCode.lean (index world):
#   d2term (β,h) (m,x): m=4 → [x⁻¹h ∈ b_β]; grid, m%2=β → [h·x⁻¹ ∈ a_{m//2}]
#   cmTerm (α,y) (m,x): m=4 → [x·y⁻¹ ∈ a_α]; grid, m//2=α → [y⁻¹x ∈ b_{m%2}]


def d2term(beta, h, m, x):
    if m == 4:
        return int(mul(int(inv[x]), h) in S[f"b{beta}"])
    if m % 2 == beta:
        return int(mul(h, int(inv[x])) in S[f"a{m // 2}"])
    return 0


def cmterm(alpha, y, m, x):
    if m == 4:
        return int(mul(x, int(inv[y])) in S[f"a{alpha}"])
    if m // 2 == alpha:
        return int(mul(int(inv[y]), x) in S[f"b{m % 2}"])
    return 0


# rebuild the emitter's instance tables (identical code path)
def Lm(a):
    return np.array([[int(mul(h, int(inv[x])) in a) for x in range(n)]
                     for h in range(n)], dtype=np.uint8)


def Rm(b):
    return np.array([[int(mul(int(inv[x]), h) in b) for x in range(n)]
                     for h in range(n)], dtype=np.uint8)


def Lms(a):
    return np.array([[int(mul(x, int(inv[y])) in a) for x in range(n)]
                     for y in range(n)], dtype=np.uint8)


def Rms(b):
    return np.array([[int(mul(int(inv[y]), x) in b) for x in range(n)]
                     for y in range(n)], dtype=np.uint8)


A0, A1 = Lm(S["a0"]), Lm(S["a1"])
Rb = {0: Rm(S["b0"]), 1: Rm(S["b1"])}
B0s, B1s = Rms(S["b0"]), Rms(S["b1"])
As = {0: Lms(S["a0"]), 1: Lms(S["a1"])}
A1i = _inv2(A1)
B1si = _inv2(B1s)
B0si = _inv2(B0s)
As1i = _inv2(As[1])
PU, lnU, annU = _pseudo_inv(A0)
PT, lnT, annT = _pseudo_inv(As[0])

idcols = [1 << i for i in range(n)]
inst = {}
for beta in (0, 1):
    Ri = _inv2(Rb[beta])
    inst[f"X{beta}"] = dict(
        tUW=_cols(A1i @ A0 % 2), tTW=_cols(A1i @ Rb[beta] % 2),
        tUT=_cols(Ri @ A0 % 2), tWT=_cols(Ri @ A1 % 2),
        tWC=_cols(A1), tTC=_cols(Rb[beta]), tP=_cols(PU),
        ln0=_msk(lnU[0]), ln1=_msk(lnU[1]),
        ann=[_msk(a) for a in annU],
        tUC=[], tWC2=[], tPT=[], lnT0=0, lnT1=0, annT=[],
        modes=(0, 1, 2),
    )
for alpha in (0, 1):
    inst[f"Z{alpha}"] = dict(
        tUW=_cols(B1si @ B0s % 2), tTW=_cols(B1si @ As[alpha] % 2),
        tUT=[], tWT=[],
        tWC=_cols(B0si @ B1s % 2), tTC=_cols(B0si @ As[alpha] % 2),
        tP=list(idcols), ln0=0, ln1=0, ann=[0],
        tUC=_cols(B0s), tWC2=_cols(B1s),
        tPT=_cols(PT) if alpha == 0 else _cols(As1i),
        lnT0=_msk(lnT[0]) if alpha == 0 else 0,
        lnT1=_msk(lnT[1]) if alpha == 0 else 0,
        annT=[_msk(a) for a in annT] if alpha == 0 else [0],
        modes=(0, 2, 3),
    )

# the NEW tables the bridge layer needs emitted
NEW = dict(
    tRawAX=_cols(A0),
    tCinvX0=_cols(_inv2(Rb[0])),
    tCinvX1=_cols(_inv2(Rb[1])),
    tRawCZ0=_cols(As[0]),
    tRawCZ1=_cols(As[1]),
    tAinvZ=_cols(B0si),
)

# ------------------------------------------- 0. FloorData.lean cross-parse
txt = FLOORDATA.read_text()


def parse_list(name):
    m = re.search(rf"^def {name} : List Nat := \[(.*?)\]$", txt,
                  re.S | re.M)
    if m is None:
        m2 = re.search(rf"^def {name} : List Nat := \[\]$", txt, re.M)
        if m2:
            return []
        m3 = re.search(rf"^def {name} : List Nat := ([^\[\n]+)$", txt, re.M)
        assert m3, f"cannot parse {name}"
        parts = [p.strip() for p in m3.group(1).split("++")]
        return sum((parse_list(p) for p in parts), [])
    return [int(v) for v in re.split(r"[,\s]+", m.group(1).strip()) if v]


def parse_splits(name):
    m = re.search(rf"def {name} : List \(Nat × Nat × Nat × Nat\) := \[(.*?)\]\n",
                  txt, re.S)
    assert m, f"cannot parse {name}"
    return [tuple(int(v) for v in tup.split(","))
            for tup in re.findall(r"\(([^)]*)\)", m.group(1))]


emitted_ok = True
for name, I in inst.items():
    for fld in ("tUW", "tTW", "tUT", "tWT", "tWC", "tTC", "tP",
                "tUC", "tWC2", "tPT"):
        if parse_list(f"{fld}{name}") != I[fld]:
            emitted_ok = False
            print(f"    mismatch: {fld}{name}")
    for fld in ("ann", "annT"):
        if parse_list(f"{fld}{name}") != I[fld]:
            emitted_ok = False
            print(f"    mismatch: {fld}{name}")
check("emitted FloorData tables match recomputation", emitted_ok)

cls = {name: parse_list(f"cls{name}") for name in inst}
splits = {name: parse_splits(f"splits{name}") for name in inst}
rowsXpk = parse_list("rowsXpk")
rowsZpk = parse_list("rowsZpk")
triples = json.loads((GEN.OUT / "m4_triples.json").read_text())
for name in inst:
    want = sorted((u | (w << 30) | (t << 60))
                  for (u, w, t) in map(tuple, triples[name]["sols"])
                  if (u, w, t) != (0, 0, 0))
    check(f"cls{name} matches census ({len(want)})", cls[name] == want)

# ------------------------------------------- 1. op→table bridge basis facts
# X instance β: op_m(δ_i)(y) = dualBfn(embed m δ_i)(β, y) = d2term (β,y) (m, i)
# tables: m=β → tRawAX; m=2+β → tWC; m=4 → tTC; else []
for beta in (0, 1):
    I = inst[f"X{beta}"]
    tb = {beta: NEW["tRawAX"], 2 + beta: I["tWC"], 4: I["tTC"]}
    ok = True
    for m in range(5):
        tbl = tb.get(m, [])
        for i in range(n):
            col = _msk([d2term(beta, y, m, i) for y in range(n)])
            want = tbl[i] if tbl else 0
            ok &= col == want
    check(f"X{beta} bridge basis facts (5 blocks × 30)", ok)

# Z instance α: op_m(δ_i)(y) = cmTerm (α,y) (m,i)
# tables: m=2α → tUC; m=2α+1 → tWC2; m=4 → tRawCZα; else []
for alpha in (0, 1):
    I = inst[f"Z{alpha}"]
    tb = {2 * alpha: I["tUC"], 2 * alpha + 1: I["tWC2"],
          4: NEW[f"tRawCZ{alpha}"]}
    ok = True
    for m in range(5):
        tbl = tb.get(m, [])
        for i in range(n):
            col = _msk([cmterm(alpha, y, m, i) for y in range(n)])
            want = tbl[i] if tbl else 0
            ok &= col == want
    check(f"Z{alpha} bridge basis facts (5 blocks × 30)", ok)

# ------------------------------------------- 2. composite identities
for beta in (0, 1):
    I = inst[f"X{beta}"]
    ci = NEW[f"tCinvX{beta}"]
    check(f"X{beta} D1 tCinv∘tTC = id",
          all(_xf(ci, I["tTC"][i]) == 1 << i for i in range(n)))
    check(f"X{beta} D2 tCinv∘tRawAX = tUT",
          all(_xf(ci, NEW["tRawAX"][i]) == I["tUT"][i] for i in range(n)))
    check(f"X{beta} D3 tCinv∘tWC = tWT",
          all(_xf(ci, I["tWC"][i]) == I["tWT"][i] for i in range(n)))
    check(f"X{beta} D4 tTW∘tWT = id",
          all(_xf(I["tTW"], I["tWT"][i]) == 1 << i for i in range(n)))
    check(f"X{beta} D5 tTW∘tUT = tUW",
          all(_xf(I["tTW"], I["tUT"][i]) == I["tUW"][i] for i in range(n)))
for alpha in (0, 1):
    I = inst[f"Z{alpha}"]
    ai = NEW["tAinvZ"]
    rc = NEW[f"tRawCZ{alpha}"]
    check(f"Z{alpha} D1 tAinvZ∘tUC = id",
          all(_xf(ai, I["tUC"][i]) == 1 << i for i in range(n)))
    check(f"Z{alpha} D2 tAinvZ∘tWC2 = tWC",
          all(_xf(ai, I["tWC2"][i]) == I["tWC"][i] for i in range(n)))
    check(f"Z{alpha} D3 tAinvZ∘tRawC = tTC",
          all(_xf(ai, rc[i]) == I["tTC"][i] for i in range(n)))
    check(f"Z{alpha} D4 tUW∘tWC = id",
          all(_xf(I["tUW"], I["tWC"][i]) == 1 << i for i in range(n)))
    check(f"Z{alpha} D5 tUW∘tTC = tTW",
          all(_xf(I["tUW"], I["tTC"][i]) == I["tTW"][i] for i in range(n)))
    check(f"Z{alpha} D6 tP = id columns",
          I["tP"] == idcols)

# ------------------------------------------- 3. coset machinery
for beta in (0, 1):
    I = inst[f"X{beta}"]
    ann = set(I["ann"])
    corr = [_xf(I["tP"], NEW["tRawAX"][i]) ^ (1 << i) for i in range(n)]
    check(f"X{beta} D6c corr columns ∈ ann", all(c in ann for c in corr))
    check(f"X{beta} D7 ann xor-closed ∋ 0",
          0 in ann and all(a ^ b in ann for a in ann for b in ann))
    check(f"X{beta} D8 ln filters vanish on tRawAX columns",
          all(_pop(I["ln0"] & NEW["tRawAX"][i]) % 2 == 0
              and _pop(I["ln1"] & NEW["tRawAX"][i]) % 2 == 0
              for i in range(n)))
for alpha in (0, 1):
    I = inst[f"Z{alpha}"]
    rc = NEW[f"tRawCZ{alpha}"]
    annT_ = set(I["annT"])
    corrT = [_xf(I["tPT"], rc[i]) ^ (1 << i) for i in range(n)]
    check(f"Z{alpha} D7c corrT columns ∈ annT",
          all(c in annT_ for c in corrT))
    check(f"Z{alpha} D8 annT xor-closed ∋ 0",
          0 in annT_ and all(a ^ b in annT_ for a in annT_ for b in annT_))
    check(f"Z{alpha} D9 lnT filters vanish on tRawC columns",
          all(_pop(I["lnT0"] & rc[i]) % 2 == 0
              and _pop(I["lnT1"] & rc[i]) % 2 == 0 for i in range(n)))

# ------------------------------------------- 4. parity feeders (odd columns)
for name, tbls in (("X0", ["tRawAX", "tWC", "tTC"]),
                   ("X1", ["tRawAX", "tWC", "tTC"]),
                   ("Z0", ["tUC", "tWC2", "tRawCZ0"]),
                   ("Z1", ["tUC", "tWC2", "tRawCZ1"])):
    I = inst[name]
    ok = True
    for t in tbls:
        tbl = NEW.get(t, I.get(t))
        ok &= all(_pop(c) % 2 == 1 for c in tbl)
    check(f"{name} raw operator columns all odd", ok)

# ------------------------------------------- 5. split coverage per instance
for name, I in inst.items():
    sset = set(splits[name])
    got_modes = {m for (_, _, _, m) in sset}
    ok = got_modes <= set(I["modes"])
    for p, q, r in product(range(9), repeat=3):
        if p + q + r <= 8 and (p + q + r) % 2 == 0:
            ok &= any((p, q, r, m) in sset for m in I["modes"])
    check(f"{name} coverage: even (p,q,r) ≤ 8 all assigned, modes ⊆ "
          f"{I['modes']}", ok)

# ------------------------------------------- 6. no classified triple has t=0
for name in inst:
    check(f"cls{name} all have t ≠ 0",
          all(pk >> 60 != 0 for pk in cls[name]))

# ------------------------------------------- 7. rows = packed row chains
m60 = (1 << 60) - 1


def pack5(pk0, pk1):
    return (pk0 & m60) | ((pk1 & m60) << 60) | ((pk0 >> 60) << 120)


def packTriple(u, w, t):
    return u | (w << 30) | (t << 60)


rowsX_want = []
for alpha in (0, 1):
    for y in range(n):
        blk = [
            _msk([cmterm(alpha, y, m, x) for x in range(n)])
            for m in range(5)
        ]
        rowsX_want.append(
            pack5(packTriple(blk[0], blk[2], blk[4]),
                  packTriple(blk[1], blk[3], blk[4])))
check("rowsXpk = packed H_Z row chains (blocks 0,2|1,3|4)",
      sorted(rowsX_want) == rowsXpk)

rowsZ_want = []
for beta in (0, 1):
    for h in range(n):
        blk = [
            _msk([d2term(beta, h, m, x) for x in range(n)])
            for m in range(5)
        ]
        rowsZ_want.append(
            pack5(packTriple(blk[0], blk[1], blk[4]),
                  packTriple(blk[2], blk[3], blk[4])))
check("rowsZpk = packed H_X row chains (blocks 0,1|2,3|4)",
      sorted(rowsZ_want) == rowsZpk)

# --------------------------------- 8. end-to-end sanity on kernel samples
rng = np.random.default_rng(32)
HX = np.zeros((60, 150), dtype=np.uint8)
HZ = np.zeros((60, 150), dtype=np.uint8)
for beta in (0, 1):
    for h in range(n):
        for m in range(5):
            for x in range(n):
                HX[beta * 30 + h, m * 30 + x] = d2term(beta, h, m, x)
for alpha in (0, 1):
    for y in range(n):
        for m in range(5):
            for x in range(n):
                HZ[alpha * 30 + y, m * 30 + x] = cmterm(alpha, y, m, x)
check("H_X·H_Zᵀ = 0", not ((HX @ HZ.T) % 2).any())


def nullbasis(M):
    A = M.copy().astype(np.uint8)
    nn = A.shape[1]
    piv, r = [], 0
    for c in range(nn):
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        A[[r, r + nz[0]]] = A[[r + nz[0], r]]
        for q in np.flatnonzero(A[:, c]):
            if q != r:
                A[q] ^= A[r]
        piv.append(c)
        r += 1
    free = [c for c in range(nn) if c not in piv]
    out = []
    for f in free:
        v = np.zeros(nn, dtype=np.uint8)
        v[f] = 1
        for k, pc in enumerate(piv):
            v[pc] = A[k, f]
        out.append(v)
    return np.array(out, dtype=np.uint8)


KX = nullbasis(HX)
ok = True
for _ in range(60):
    coeff = rng.integers(0, 2, KX.shape[0]).astype(np.uint8)
    v = (coeff @ KX) % 2
    blk = [_msk(v[m * 30:(m + 1) * 30]) for m in range(5)]
    for beta in (0, 1):
        I = inst[f"X{beta}"]
        mu, mw, mt = blk[beta], blk[2 + beta], blk[4]
        ok &= (_xf(NEW["tRawAX"], mu) ^ _xf(I["tWC"], mw)
               ^ _xf(I["tTC"], mt)) == 0
        ok &= mw == _xf(I["tUW"], mu) ^ _xf(I["tTW"], mt)
        ok &= mt == _xf(I["tUT"], mu) ^ _xf(I["tWT"], mw)
        c = _xf(I["tWC"], mw) ^ _xf(I["tTC"], mt)
        ok &= _pop(I["ln0"] & c) % 2 == 0 and _pop(I["ln1"] & c) % 2 == 0
        ok &= (mu ^ _xf(I["tP"], c)) in set(I["ann"])
check("X side: random ker-H_X vectors satisfy every derived form", ok)

KZ = nullbasis(HZ)
ok = True
for _ in range(60):
    coeff = rng.integers(0, 2, KZ.shape[0]).astype(np.uint8)
    v = (coeff @ KZ) % 2
    blk = [_msk(v[m * 30:(m + 1) * 30]) for m in range(5)]
    for alpha in (0, 1):
        I = inst[f"Z{alpha}"]
        rc = NEW[f"tRawCZ{alpha}"]
        mu, mw, mt = blk[2 * alpha], blk[2 * alpha + 1], blk[4]
        ok &= (_xf(I["tUC"], mu) ^ _xf(I["tWC2"], mw) ^ _xf(rc, mt)) == 0
        ok &= mu == _xf(I["tWC"], mw) ^ _xf(I["tTC"], mt)
        ok &= mw == _xf(I["tUW"], mu) ^ _xf(I["tTW"], mt)
        c = _xf(I["tUC"], mu) ^ _xf(I["tWC2"], mw)
        ok &= (_pop(I["lnT0"] & c) % 2 == 0
               and _pop(I["lnT1"] & c) % 2 == 0)
        ok &= (mt ^ _xf(I["tPT"], c)) in set(I["annT"])
check("Z side: random ker-H_Z vectors satisfy every derived form", ok)

# ---------------------------------------------------------------- verdict
print()
if fails:
    print(f"FAILURES ({len(fails)}):")
    for f in fails:
        print(f"  - {f}")
    sys.exit(1)
print("ALL BRIDGE-LAYER FACTS VERIFIED — safe to emit and formalize.")
