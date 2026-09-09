# Canonical HM: linear recovery must pay a width–load tradeoff

2026-09-08. Original binary cylinder, not a folded x quotient. This
note tests a concrete sufficient route to canonical HM: recover all
inside boundary bits linearly from outside bits, using each outside
bit at most three times. It proves that **no fixed-width recovery
map can do this at unbounded height**. It does not refute HM or rule
out nonlinear discharge or a growing-width analytic recovery map.

## 1. Recovery would imply the exact weighted inequality

Let q be a power of two and put

```
R_q=F2[t]/(t^q+1),       u=t+t^2,
P=x^3+u,                Q=1+(x+x^2)t^3,
B_q={(Qh,Ph): h in R_q[x^+/-1]},
S={(first block,0,j),(second block,0,j),(second block,1,j): j mod q}.
```

All boundary generators have finite x support. Suppose T is a
finite set of qubits outside S, and an F2-linear map

```
L:F2^T -> F2^S
```

satisfies `L(s|T)=s|S` for every s in B_q. Write c for the largest
Hamming weight of a column of L. Then

```
|s intersect S| <= c |s intersect T| <= c |s outside S|.
```

In particular c<=3 would prove canonical HM. No assumption that L
commutes with row shifts is needed in the theorem below.

Count a *block-column* as one of the two physical blocks at a fixed
integer x. If T is contained in K such block-columns, the following
necessary bound holds for every such recovery:

```
K c >= 3 ceil(q/3)-3.                                  (1.1)
```

Thus c<=3 requires `K>=ceil(q/3)-1`. If both blocks are allowed in
an x interval of W columns, it requires
`W>=ceil((ceil(q/3)-1)/2)`. The same obstruction applies to any
fixed x window even if its coefficients depend arbitrarily on q
and are nonlocal around the periodic row direction.

This is a necessary resource bound, not a lower bound on code
distance. Global recovery with growing K is not excluded.

## 2. A period-three obstruction to a local recovery row

Consider the same binary polynomials on a cylinder of row period
three, retaining the unbounded x coordinate. Define

```
k=1+t,       e=(k,(1+x)k).
```

Since `t^3=1`,

```
P e_first+Q e_second
  = (1+t+t^2)k = 1+t^3 = 0.                            (2.1)
```

The support of e lies entirely in S. For every individual coordinate
of S, some row translate of e has a nonzero value there.
It is not a compact boundary: a nonzero compact h makes Qh have
distinct lowest and highest x columns separated by at least two,
since Q's extreme coefficients are units. It cannot equal the
nonzero monocolumn first block of e.

Over the *plane* Laurent ring `F2[x^+/-1,t^+/-1]`, P and Q are
coprime. The fixed Smith/determinant calculation in
`A40_PARITY_TOR_PROOF.md` §2 proves this, as well as regularity.
Consequently every finite-support cochain phi annihilating all
plane boundaries has the form

```
phi=(bar P a,bar Q a),                                 (2.2)
```

where bar inverts both variables. Indeed boundary orthogonality
is `bar Q phi_first+bar P phi_second=0`; coprimality in the
Laurent UFD gives (2.2). Such a cochain is a finite combination
of syndrome rows. Its reduction to any row period annihilates
every cycle, including (2.1).

It follows that no finite plane cochain can both annihilate all
boundaries and restrict to **one single bit of S**. Reduction to
period three would pair it to one with a suitable translate of e,
contradicting (2.2). This uses a nontrivial cycle on a different
period as a locality obstruction, not as a counterexample to HM
on a binary-height cylinder.

## 3. Every recovery row must wind around the row direction

For each output bit i in S, extend row i of L to a cylinder
cochain phi_i by adding its unit coefficient at inside bit i and
setting all other inside coefficients to zero. The recovery
identity means phi_i annihilates every compact binary boundary.

**Lemma.** The occupied rows of phi_i have no cyclic gap of three
or more empty rows. In particular `|phi_i|>=ceil(q/3)`.

**Proof.** If such a gap existed, choose the origin so that its
support lies in rows `0,...,q-g-1`, with g>=3 empty rows. Unwrap
the support to the plane. Boundary orthogonality on the cylinder
says the reduction of

```
bar Q phi_first+bar P phi_second
```

is zero. All its row exponents lie in `[-3,q-g-1]`, an interval
containing at most q integers. Projection modulo q is injective
on this interval, so the expression is already zero on the plane.
The unwrapped cochain is finite because T is finite. Its restriction
to S is still one bit. This contradicts §2. Thus every cyclic
interval of three rows meets its support. If n rows are occupied,
the maximum allowed gaps give q<=3n, proving the weight bound.
For q=1,2 the same bound follows from its required inside bit. QED.

There are 3q output rows. Each phi_i has exactly one inside bit,
so the number of ones in L is at least

```
3q(ceil(q/3)-1).
```

On the other hand it is at most `c|T|<=c K q`. Dividing by q
proves (1.1). No row equivariance, representative minimization,
or numerical rank extrapolation is used.

## 4. A fixed-width recovery exists, but its load grows

The obstruction is not failure of linear erasure recovery itself.
Here is an explicit valid recovery from ten outside block-columns.
For inside row masks `(i,j,k)` representing `(i,j+xk)`, the four
actual syndrome columns are

```
z0=u i+j,
z1=k+t^3 j,
z2=t^3(j+k),
z3=i+t^3 k.                                            (4.1)
```

Set `Phi=1+t+t^2=1+u`. Direct elimination gives

```
j=Phi^-1 [z0+u z1+u(1+t^-3)z2+u z3],
k=z1+t^3 j,
i=z3+t^3 z1+t^6 j.                                   (4.2)
```

For every `3 ∤ q`, Phi is invertible modulo `t^q+1`. At a
binary height this also follows from `Phi^q=1`. For a boundary,
take z0,...,z3 to be the four central syndrome columns generated
by its *outside* qubits. Its total syndrome is zero, so these
are exactly the inside syndrome (4.1). Equations (4.2) recover
the inside bits correctly.

Only outside first-block columns `-3,-2,-1,1,2,3` and outside
second-block columns `-2,-1,2,3` can contribute to these four
syndrome columns. Thus K=10 suffices for recovery, but not for
a universal c=3 budget. By (1.1), **no** recovery using only these
ten block-columns can have c<=3 once `ceil(q/3)>11`, including
every dyadic q>=64.

The specific map (4.2) already spends far too much on the first
outside column. A single first-block input at x=-3 has z0=1 and
z1=z2=z3=0. Its three recovered masks are

```
(t^6 Phi^-1, Phi^-1, t^3 Phi^-1),
```

of combined weight `3|Phi^-1|`. The inverse is explicit:

```
q=3v+1: coefficients at rows 0 or 2 mod 3 are one;
q=3v+2: coefficients at rows 1 or 2 mod 3 are one.
```

Substitution in `Phi Phi^-1=1` proves these formulas, including
the wrap rows. Therefore this column alone has weight

```
2q+1 when q=1 mod 3,       2q-1 when q=2 mod 3.          (4.3)
```

This does not prove all possible recovery maps have this larger
load; the universal obstruction is only (1.1). It explains why
an algebraically exact fixed-width inverse does not give HM's
sharp factor three.

## 5. Disposition and verification

The linear-recovery route is still viable only with a growing
spatial span or an approach not constrained by a uniform linear
column load. A small parameter-independent stencil is now ruled
out analytically, rather than merely failing an optimization run.
The result does not contradict S8's fixed **generator** strip
certificate: restricting boundary generators and restricting a
global recovery map's input window are different assertions.

`scripts/ibm_canonical_recovery_obstruction.py` verifies (2.1),
the three inverse identities (4.2), the explicit Phi inverse,
the actual boundary-recovery equations, and the resulting finite
decoder loads at prescribed small heights. It searches neither
cosets nor code families. The all-q width–load bound is the
unwrapping and counting proof in §§2–3, not the sampled loads.

The independent audit in `A42_S9_tail_elimination.md` §9 passes
the plane lifting, period-three pairing, endpoint count, general
non-equivariant decoder count, and the explicit inverse. The parent
identity replay passes in 0.0004 seconds with about 17 MB RSS.
