"""A19 (M)@24 census redesign, session 1: the A22 (eps,delta) CRT fibering
ported to BY = (30,3), A = x^9 + y + y^2, B = 1 + x^25 + x^26.

Structure: Z_30 = Z_5<z> x Z_6<v> via CRT (z = x^6, v = x^25; x^a maps to
z^(a mod 5) v^(a mod 6)), so G = Z_5(fiber) x H with base H = Z_6<v> x Z_3<y>
(18 sites). F2[z]/(z^5-1) ~= F2 (eps: z->1) x GF(16) (delta: z->zeta),
zeta = t^3 in GF(16) = F2[t]/(t^4+t+1) (t primitive of order 15, zeta order 5).

Per-site exact weight table (A22, reused verbatim — same Z_5 fiber):
    W(0,0)=0   W(1,0)=5   W(0,mu5)=4   W(1,mu5)=1   W(0,other)=2   W(1,other)=3
and |u| = sum over the 18 sites of W(u_eps(site), u_delta(site)).

Polynomials in fibered coordinates:
    A = x^9 + y + y^2  = z^4 v^3 + y + y^2   (9 mod 5 = 4, 9 mod 6 = 3)
    B = 1 + x^25 + x^26 = 1 + v + z v^2      (25 -> (0,1); 26 -> (1,2))
    eps-side:  A_e = v^3 + y + y^2,  B_e = 1 + v + v^2      over F2[H]
    delta-side: At = zeta^4 v^3 + y + y^2,  Bt = 1 + v + zeta v^2  over GF(16)[H]

Session-1 scope: (1) verify the weight formula numerically; (2) component
analysis of the four base-ring operators (ranks/kernels; is there an A22-style
free-coordinate substitution? — expectation: NO on the eps-side, H has a Z_2
radical, so the (1+s)-layer engine is needed); (3) push the 7 known M12 census
classes through the fibered lens (delta-site occupancy = the m-site-theorem
calibration for the <= 23 classification).
"""
import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.checks import bb_check_matrices
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly

CENSUS = LAB / "data" / "a19" / "m12_census_classes.jsonl"

# ---------------------------------------------------------------- GF(16)
# F2[t]/(t^4 + t + 1), elements as ints 0..15 (bit i = coeff of t^i).
EXP = [1]
for _ in range(14):
    e = EXP[-1] << 1
    if e & 0x10:
        e ^= 0x13          # t^4 = t + 1
    EXP.append(e)
LOG = {e: i for i, e in enumerate(EXP)}


def gmul(a, b):
    if a == 0 or b == 0:
        return 0
    return EXP[(LOG[a] + LOG[b]) % 15]


ZETA = EXP[3]              # order 5
MU5 = {EXP[(3 * i) % 15] for i in range(5)}   # {1, zeta, ..., zeta^4}

W_TABLE = {}
for eps in (0, 1):
    for d in range(16):
        if d == 0:
            W_TABLE[(eps, d)] = 5 if eps else 0
        elif d in MU5:
            W_TABLE[(eps, d)] = 1 if eps else 4
        else:
            W_TABLE[(eps, d)] = 3 if eps else 2

# ---------------------------------------------------------------- fibering
ELL, M = 30, 3
G = AbelianGroup((ELL, M))
N = ELL * M                # 90 per block
H_SITES = [(q, b) for q in range(6) for b in range(3)]   # (v-exp, y-exp)
SITE_IDX = {s: i for i, s in enumerate(H_SITES)}


def g_of(p, q, b):
    """Group element (a, b) with a = CRT(p mod 5, q mod 6)."""
    for a in range(30):
        if a % 5 == p and a % 6 == q:
            return (a, b)
    raise AssertionError


G_OF = {(p, q, b): g_of(p, q, b) for p in range(5) for q in range(6)
        for b in range(3)}
IDX = {g: i for i, g in enumerate(G)}


def fiber_split(u):
    """u: F2 vector over one block (len 90) -> (eps[18] F2, delta[18] GF16)."""
    eps = np.zeros(18, dtype=np.uint8)
    dlt = np.zeros(18, dtype=np.int64)
    for (q, b), s in SITE_IDX.items():
        e, d = 0, 0
        for p in range(5):
            if u[IDX[G_OF[(p, q, b)]]]:
                e ^= 1
                d ^= EXP[(3 * p) % 15]   # zeta^p
        eps[s], dlt[s] = e, d
    return eps, dlt


def fibered_weight(u):
    eps, dlt = fiber_split(u)
    return sum(W_TABLE[(int(eps[s]), int(dlt[s]))] for s in range(18))


# ---------------------------------------------------------------- checks
A = Poly.from_string("x^9 + y + y^2", G)
B = Poly.from_string("1 + x^25 + x^26", G)
ch = bb_check_matrices(A, B)
HX = ch.H_X % 2

rng = np.random.default_rng(19)
print("== (1) weight-formula verification ==")
ok = 0
for _ in range(2000):
    f = rng.integers(0, 2, N).astype(np.uint8)
    bvec = (HX.T @ np.concatenate([f, np.zeros(0, dtype=np.uint8)])
            if False else None)
    # b = H_X^T f over checks: rows of H_X indexed by checks (len N)
    b = (HX.T @ f) % 2 if HX.shape[0] == N else None
    assert b is not None
    wl, wr = int(b[:N].sum()), int(b[N:].sum())
    fl, fr = fibered_weight(b[:N]), fibered_weight(b[N:])
    assert wl == fl and wr == fr, (wl, fl, wr, fr)
    ok += 1
print(f"  {ok}/2000 random stabilizers: |.| == sum W(eps,delta) on both "
      f"blocks  PASS", flush=True)

# ---------------------------------------------------------------- operators
print("== (2) base-ring operator components ==")
HG = AbelianGroup((6, 3))


def conv_rank(supp_coeffs, field):
    """Rank of convolution-by-P on field[H]; supp_coeffs: {(q,b): coeff}."""
    n = 18
    Mm = [[0] * n for _ in range(n)]
    for j, (qj, bj) in enumerate(H_SITES):
        for (dq, db), c in supp_coeffs.items():
            i = SITE_IDX[((qj + dq) % 6, (bj + db) % 3)]
            if field == "F2":
                Mm[i][j] ^= c
            else:
                Mm[i][j] ^= 0  # placeholder; GF16 handled below
    if field == "F2":
        Mat = np.array(Mm, dtype=np.uint8)
        R = Mat.copy()
        r = 0
        for cidx in range(n):
            piv = next((i for i in range(r, n) if R[i, cidx]), None)
            if piv is None:
                continue
            R[[r, piv]] = R[[piv, r]]
            for i in range(n):
                if i != r and R[i, cidx]:
                    R[i] ^= R[r]
            r += 1
        return r
    raise ValueError


def gf16_conv_rank(supp_coeffs):
    n = 18
    Mat = [[0] * n for _ in range(n)]
    for j, (qj, bj) in enumerate(H_SITES):
        for (dq, db), c in supp_coeffs.items():
            i = SITE_IDX[((qj + dq) % 6, (bj + db) % 3)]
            Mat[i][j] ^= c
    r = 0
    for cidx in range(n):
        piv = next((i for i in range(r, n) if Mat[i][cidx]), None)
        if piv is None:
            continue
        Mat[r], Mat[piv] = Mat[piv], Mat[r]
        inv = EXP[(15 - LOG[Mat[r][cidx]]) % 15]
        Mat[r] = [gmul(x, inv) for x in Mat[r]]
        for i in range(n):
            if i != r and Mat[i][cidx]:
                c = Mat[i][cidx]
                Mat[i] = [x ^ gmul(c, y) for x, y in zip(Mat[i], Mat[r])]
        r += 1
    return r


A_e = {(3, 0): 1, (0, 1): 1, (0, 2): 1}          # v^3 + y + y^2
B_e = {(0, 0): 1, (1, 0): 1, (2, 0): 1}          # 1 + v + v^2
A_d = {(3, 0): EXP[12], (0, 1): 1, (0, 2): 1}    # zeta^4 v^3 + y + y^2
B_d = {(0, 0): 1, (1, 0): 1, (2, 0): ZETA}       # 1 + v + zeta v^2

for name, supp, fld in [("A_eps", A_e, "F2"), ("B_eps", B_e, "F2")]:
    r = conv_rank(supp, fld)
    print(f"  {name}: rank {r}/18, kernel dim {18 - r}"
          f"{'  (INVERTIBLE — free-substitution available)' if r == 18 else ''}",
          flush=True)
for name, supp in [("A_delta", A_d), ("B_delta", B_d)]:
    r = gf16_conv_rank(supp)
    print(f"  {name}: GF(16)-rank {r}/18, kernel dim {18 - r}"
          f"{'  (INVERTIBLE — free-substitution available)' if r == 18 else ''}",
          flush=True)

# ---------------------------------------------------------------- census map
print("== (3) M12 census classes through the fibered lens ==")
for line in CENSUS.read_text().splitlines():
    row = json.loads(line)
    if "w" not in row:
        continue
    b = np.zeros(2 * N, dtype=np.uint8)
    b[row["b_support"]] = 1
    out = []
    for blk, u in (("L", b[:N]), ("R", b[N:])):
        eps, dlt = fiber_split(u)
        ne = int(eps.sum())
        nd = int((dlt != 0).sum())
        wf = sum(W_TABLE[(int(eps[s]), int(dlt[s]))] for s in range(18))
        out.append(f"{blk}: wt {int(u.sum())}=Σ{wf}, eps-sites {ne}, "
                   f"delta-sites {nd}")
    print(f"  w={row['w']}: " + "; ".join(out), flush=True)
