# A42 S16: every residual D excursion is impossible

2026-09-08. Bounded Codex proof continuation. New outputs are this
note and `scripts/a42_s16_codex_D_excursion.py`. The sole change to
an older file corrects "twenty" to "eighteen" in S14's section 4
heading; its ledger is unchanged.

## 0. Result and scope

**Theorem.** Let q be any positive integer. In the compact canonical
syndrome problem below, suppose aggregation loss is zero, the only
positive-slack rows are one entry with ID in `{0,1,4,5}` and one exit
with ID in `{9,10}`, each of slack one, and every other row has slack
zero. Then the actual row walk cannot contain `D->C`.

In particular, all the D-excursion cases assigned to S16 are empty.
This includes all four proposed height-eight word shapes, for both
entry residues and both BE/ED or BF/FD choices. The proof is stronger
than the proposed reductions: it needs neither dyadic q nor bounds
on the ABC, g-only, or D-loop lengths. It remains valid if every
occupied residue is relaxed from a monomial to an arbitrary Laurent
polynomial in that residue, while the empty residues stay empty.

This is **not** the full excess-two sharp bound, the general sharp
floor, or an IBM distance theorem. The same eight entry/exit pairs
with **no DC** remain outside this proof. S14's other eliminations
and the parent's terminal-source elimination are separate results.

## 1. Actual equations and conventions

Work over

```
F2[x,x^-1,t]/(t^q+1),
P=x^3+t+t^2,     Q=1+(x+x^2)t^3,
N_q=1+t+...+t^(q-1).
```

The x coordinate is unbounded. For `P f+Q g=N_q`, write `X=x^3`
and `c=x+x^2`. The exact equation on row j, with indices modulo q,
is

```
X f_j+f_(j-1)+f_(j-2)+g_j+c g_(j-3)=1.             (1.1)
```

Let `Pi_0` retain exactly the powers of x divisible by three.
Its values lie in `F2[X,X^-1]`; it does **not** evaluate X at 1.
In particular, for any residue-zero polynomial U,

```
Pi_0(c U)=0.                                       (1.2)
```

S7 uses states `(a_previous,a_earlier,d_previous)`, with d in F4.
Here

```
A=000, B=100, C=011, D=110,
E=11omega, F=11omega^2,
G=00omega, H=00omega^2, O=001.
```

The complete zero-edge list, retaining the transient states, is

```
AB BC BE BF CA DC DD ED FD GH HG AO IK IL JI EN FM.
```

The additional names in this list are `I=010, J=101, K=10omega,
L=10omega^2, M=01omega, N=01omega^2`. The remaining state is
`Z=111`. The relevant physical residue supports are:

| edge | f residues | g residues |
|---|---|---|
| AB | 0 | empty |
| BC | empty | 1,2 |
| BE | 2 | 1 |
| BF | 1 | 2 |
| ED | 2 | empty |
| FD | 1 | empty |
| DD | 0 | empty |
| DC | empty | 0 |
| CA | empty | empty |
| AO | empty | 0 |
| GH | empty | 2 |
| HG | empty | 1 |
| AG, ID 0 | 0,1 | 1 |
| AH, ID 1 | 0,2 | 2 |
| OG, ID 4 | 0,2 | 1 |
| OH, ID 5 | 0,1 | 2 |
| GB, ID 9 | 2 | empty |
| HB, ID 10 | 1 | empty |

Zero aggregation loss is essential: these empty residues really
are empty, not hidden cancelling physical pairs. There is no
restriction on the Laurent exponents of the indicated sites.

## 2. Why a DC must have one of two immediate entry contexts

No allowed charged edge ends at D. Consequently the row immediately
before any DC is DD, ED, or FD and has g=0.

Suppose the first three rows after DC were zero-slack. The zero
successor of DC is CA. The next edge is AB or AO; AO has no zero
successor, so the three-row suffix must be

```
CA, AB, (BC or BE or BF).
```

Let F be the residue-zero f polynomial on AB. The residue-zero
equation on AB is `X F=1`: its two preceding f rows are empty,
and its delayed g row is the incoming-D row, which is empty.
The residue-zero equation on the following B edge is `F=1`:
that edge has no residue-zero site, and its delayed g is DC,
whose product with c has no residue-zero part. This is impossible.
Thus the first following charge occurs within three rows of DC.

The only charged edges in this theorem are an A/O entry and a G/H
exit. Following the complete zero-edge list from C now leaves only

```
DC, CA, (AG or AH);
DC, CA, AO, (OG or OH).                             (2.1)
```

An exit cannot be the first charge: G and H are not zero-reachable
from C. Taking AB after CA would require at least three subsequent
zero rows after DC before an eligible charge, just excluded.

After the entry, the next row must be either its GH/HG zero
successor or the charged exit GB/HB. All four possibilities have
**no residue-zero f or g site**. This statement includes h=0:
the next row is then the charged exit, not a presumed clean row.

This proves the needed local coverage without an infinite-walk
truncation, a dyadic-height condition, or a D-loop endpoint lemma.

## 3. One two-row certificate eliminates both contexts

Number the entry row zero and put `F=Pi_0(f_0)`. Let R denote the
next row, one of the g-only zero rows or the appropriate exit.

For an A entry the exact five-row context is

```
row       -3            -2  -1   0          1
edge      DD/ED/FD      DC  CA   AG/AH      R.
```

In equation (1.1) on row zero, the two preceding f rows are empty,
the current g has nonzero residue, and the delayed g is empty.
Thus its residue-zero equation is `X F=1`.
On row one, the current f and g have no residue-zero part, the
earlier f row is CA and is empty, and the delayed g is DC.
By (1.2), its residue-zero equation is `F=1`.

For an O entry the context is

```
row       -3  -2  -1   0          1
edge      DC  CA  AO   OG/OH      R.
```

Again both preceding f rows on row zero are empty. This time the
delayed g is DC, but (1.2) makes its residue-zero contribution zero.
The row-zero equation is still `X F=1`. On row one, the earlier f
row is AO and is empty, and the delayed g is CA and is empty.
Hence its residue-zero equation is again `F=1`.

Equivalently, write the two projected syndrome residuals as

```
E_0 = X F+1,
E_1 = F+1.
```

They have the polynomial left-kernel certificate

```
E_0 + X E_1 = 1+X != 0.                            (3.1)
```

The last inequality holds in the unbounded Laurent ring. Thus
they cannot both vanish. This proves the theorem. It also proves
the stated polynomial-support relaxation; no monomial exponent
comparison was used. Every delayed g term has been retained.

If an original cyclic word is too short to realize the displayed
distinct edge roles, that word is already graph-impossible. The
forced D-excursion words have length at least six for A entries
and seven for O entries, so the five-row contexts have no row
identifications. No hidden small-height endpoint case is omitted.

## 4. Audit of the proposed longer reduction

The parent's reductions are valid under the stated S14 normal form,
but none is needed for (3.1). Here is the exact audit ledger.

1. The return is `(BC,CA,AB)^n, BK, KD, DD^r, DC, CA, [AO]`,
   where `(BK,KD)` is `(BE,ED)` or `(BF,FD)`. If n>=1, the
   final AB has `X M=1` because the g predecessor of the preceding
   BC is AB or the charged exit, both g-empty. The KD row two
   rows later has `Pi_0=M=1`. Thus n=0.
2. For r>=1, the CA equation after DC gives `f_lastDD=1`.
   The D-word recurrence starts at its third entry because KD
   has g=0. S8 Lemma 5.1 consequently gives r<=7. It is valid
   here to use the whole D word, not an incorrectly extended
   interior recurrence across a nonempty delayed-g endpoint.
3. For h>=2, define rho as the exit f residue, eta as the BK/KD
   f residue, and `Z_rho=x^-rho`, `G_rho=x^(3-2rho)`.
   Then `c Z_rho=1+G_rho`. For r>=1 the first DD gives J=M;
   the KD equation and its residue-zero part give eta=rho,
   `g_h=Z_rho`, and `E+G_rho=(X+1)M`. The exit equation gives
   `E=X^-1 G_rho`, hence M=E. The BK rho-projection would be
   `X M+E=0`, contradicting `G_rho+X^-1 G_rho!=0`.
   When h=2 the delayed g in the exit equation is the charged
   entry g, which is indeed one monomial in the required residue.
4. For r=0, CA gives `J=G_eta`, `g_BK=Z_eta`; DC gives M=J
   and `g_DC=1`. KD then forces eta=rho and `E=X G_rho`,
   whereas the exit equation requires `X E=G_rho`. The resulting
   `X^2 G_rho=G_rho` is impossible. These are actual row
   equations; no fictitious first DD row is used for r=0.

The h>=2 reductions retain the incoming monomial at h=2 and are
sound for both choices of eta. The proposed reduction to h<=1
and r<=7 therefore leaves only q=8 in the dyadic case, with:

| source of entry | h | r | rows after the entry |
|---|---:|---:|---|
| A | 0 | 2 | exit, BK, KD, DD, DD, DC, CA |
| A | 1 | 1 | g-only, exit, BK, KD, DD, DC, CA |
| O | 0 | 1 | exit, BK, KD, DD, DC, CA, AO |
| O | 1 | 0 | g-only, exit, BK, KD, DC, CA, AO |

All sixteen instantiations of this table, from the two entry
residues and the two BK choices, have the same certificate (3.1).
The cyclic eight-row verifier below checks exactly those templates,
in addition to the universal local contexts. Neither a 24-by-10
elimination nor a monomial-realization search is necessary.

## 5. Exact residual and verification scope

After combining this theorem with the separately banked S14 and
terminal-source eliminations, only the **DC-free** return shapes
for entries 0/1/4/5 and exits 9/10 remain in this branch:

```
entry, (alternating GH/HG)^h, exit,
(BC,CA,AB)^n, BC, CA, [AO].                         (5.1)
```

The optional AO occurs exactly for an O entry. The exit is fixed
by the G/H state reached after h alternating edges. Under the
existing S14 bounds, h<=9, n<=2, and dyadic q is in `{4,8,16}`.
All row equations (1.1) and the actual Laurent residue supports
remain required. The precise unresolved sharp assertion in this
subbranch is `|supp(f,g) intersect S|<=3`, where S consists of
f-column zero and g-columns zero and one. No claim about (5.1)
is proved here beyond excluding every return containing DC.

The new checker uses only the Python standard library. It verifies
the fixed 16-state Fourier supports and slacks, the exact projected
two-row identities in all sixteen local contexts, and all sixteen
explicit cyclic height-eight instantiations. Its symbolic variables
retain block, physical row, residue, and x shift. The certificate
is the exact nonzero polynomial `1+x^3`; it is not a rank test at
sampled values of X or an exponent-box computation.

Run from the repository root:

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a42_s16_codex_D_excursion.py
```

The script sets a 55-second wall timer and 50-second CPU limit and
checks peak RSS against 768 MiB. No solver, graph-family census,
support search, or height sweep is present. The old S14 ledger is
left as its original partial result, except for the authorized
heading-count correction. No parent or sibling proof is edited.

The single diagnostic replay passed in 0.000186 seconds of measured
checker time with 16.5 MiB peak RSS. It checked no additional word
layouts or heights beyond the sixteen prescribed cyclic templates.
