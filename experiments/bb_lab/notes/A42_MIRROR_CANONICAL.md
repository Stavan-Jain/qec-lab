# A42: a separate all-height canonical theorem in the mirror cylinder

2026-09-08. Original IBM pair; proof scope is the opposite cylinder
orientation. This does not assume that the original code is its own
axis mirror. It proves the counterpart of the canonical-support
results in `A42_S7_spatial_equality.md`, not the entire cylinder floor.

Independent proof review passed, recorded at the opening of
`A42_S7_spatial_equality.md` §10. The parent's fixed convention and
nine-context diagnostic also passes in 0.0003 seconds.

## 1. Actual mirror equations

For `A=1+y+x^3 y^-1`, `B=1+x+x^-1 y^-3`, the coordinate rotation
`(x,y)->(y^-1,x)` produces `(B,bar A)`. Exchange its Koszul entries
and clear monomials to use

```
P_m=x^3 y bar A=x^3(1+y)+y^2,
Q_m=x y^3 B=1+(x+x^2)y^3.
```

Entry exchange and these unit multiplications are translations and
permutations of the two physical chain blocks, with a matching
generator translation. They preserve Hamming weight and cylinder
compactness. They do not identify this pair with the original pair.

Fix `q>=1`, `3` not dividing `q`, and put

```
R_q=F2[t]/(t^q+1),    N=1+t+...+t^(q-1),
P=x^3(1+t)+t^2,      Q=1+(x+x^2)t^3.
```

The binary boundary code is `B_m={(Qh,Ph): h in R_q[x^±]}`. All chains
are compact in x. The omega-twisted entry, after Fourier transformation
of the three-sheet factor of period `3q`, is

```
P_w=x^3(1+zeta t)+zeta^2 t^2,     Q_w=Q,
zeta^2+zeta+1=0.
```

The quaternary cycle

```
mu=(zeta N,(1+x)N)
```

has support `S={(first block,x=0,all rows),
(second block,x=0 or 1,all rows)}`, of size `3q`. It is a cycle because
on N, P_w acts as `zeta^2(x^3+1)` and Q as `1+x+x^2`. It is nontrivial:
a boundary would require `Qh=zeta N`, impossible since nonzero compact
h produces distinct lowest and highest output columns separated by
at least two. The extreme coefficients of Q are units.

The binary indicator `1_S=(N,(1+x)N)` has syndrome

```
P N+Q(1+x)N=x^3 N.                              (1.1)
```

## 2. Exact coset minimum and uniqueness

**Theorem M.** For every such q, among all compact binary solutions

```
P f+Q g=x^3 N,
```

the minimum Hamming weight is q, and the unique minimizer is
`e_0=(x^3 N,0)`. Consequently `d(1_S,B_m)=q` with the same unique leader.

### The unchanged finite certificate

At both `x=1` and `x=zeta`, P evaluates to `Phi=1+t+t^2`; Q evaluates
respectively to 1 and `(1+t)Phi`. The right-hand side evaluates to N
in both cases. Thus the exact same aggregation, Fourier triples,
16-state graph, and two-valued potential of A42 S7 §2 apply.
Summing its 128 verified edge inequalities on a q-step closed walk
gives `|f|+|g|>=q`. If equality holds, aggregation loses no weight
and all graph edges have zero slack. In particular an empty aggregated
residue is physically empty, and a singleton lifts to exactly one site.

Retain the state names A,B,C,D,E,F and residue-support table of that
note. There are two recurrent components: the six-state component
and the separate pure-g two-cycle. The latter is impossible for an
actual error: with f=0, `Qg=x^3 N` has no compact solution by the
two-extreme-column argument.

### A different spatial obstruction, not an assumed mirror symmetry

The actual row equation is now

```
x^3 f_j+x^3 f_(j-1)+f_(j-2)+g_j+(x+x^2)g_(j-3)=x^3.   (2.1)
```

Suppose an equality walk uses `D->C` at row zero. Its local context
in the six-state recurrent component must be

```
U->D->C->A->B->V->W,
U in {D,E,F}, V in {C,E,F}, W=A if V=C and W=D otherwise.
```

The edges are rows `-1,0,1,2,3,4`. Equality in aggregation and the
finite support table imply

```
g_-1=0, f_0=0, f_1=g_1=0, g_2=0,
f_2 has only residue-zero x exponents,
f_3 and f_4 have no residue-zero x exponents, g_4=0.
```

At row two, (2.1) reduces to `x^3 f_2=x^3`, hence `f_2=1`.
At row four, extract only x exponents divisible by three, retaining
their actual integer exponents. The terms `x^3 f_4,x^3 f_3` have
no such part; g_4 and g_1 vanish. The result is `f_2=x^3`, a
contradiction on the unbounded cylinder.

This verifies all nine choices of U,V at once. It is different from
S7's row-two/row-three contradiction, and proves the mirror statement
independently rather than transferring its metric conclusion.

Deleting `D->C` leaves only the D self-loop and the A-B-C three-cycle
as recurrent components of the main graph. A q-step closed walk cannot
use the latter because `3` does not divide q. Thus every row is the
D self-loop: physically g=0. Multiplication by P is injective on
compact chains: the lowest nonzero x column of an input produces the
same lowest column of output, multiplied by the unit `t^2`. Since
`P(x^3 N)=x^3 N`, the first block is uniquely determined. This proves
Theorem M for all q in its scope.

## 3. Matching physical representatives and the remaining gap

For any binary boundary s, put

```
e=1_S+s,      E_err=|e|,
k=|S minus s|=|e intersect S|,      e_out=|s minus S|.
```

The inverse three-sheet Fourier lift with omega component mu has
the exact weight

```
|v|=3q+k+3 e_out=3q+3 E_err-2k.                  (3.1)
```

Its omega class is nonzero, so it is a physical logical on the
cylinder. The unique leader e_0 lies outside S. Augmentation of
the syndrome equation gives `E_err=q mod 2`; hence any other error
has `E_err>=q+2`. It follows that, for every `3` not dividing q:

| holes k | physical lower bound |
|---|---|
| 0 | `6q` |
| 1 | `6q+4` |
| 2 | `6q+2` |
| 3 | `6q` |
| arbitrary | `4q+2` |

The no-pair stratum has exact minimum `6q`, attained by the unique
binary boundary

```
s_*=((1+x^3)N,(1+x)N)=(Qh_*,Ph_*),
h_*=(1+x)N.
```

It contains S and has exactly q outside fibres. These are statements
about the specified omega component and arbitrary compact binary
boundary changes, not all cylinder classes. The full fixed-support
HM requirement remains

```
k <= (3/2)(E_err-q).
```

Neither Theorem M nor S7 proves it for arbitrary k. General omega
supports, the pure floor, periods with extra factors of three, and
the toroidal member lower bound remain separate obligations.

For the binary tower, `P^q=1` by Frobenius, so P is a unit and every
binary cycle is a boundary. This permits replacing the boundary coset
by all solutions of (1.1) at `q=2^a`. No such identification is needed
for the minimum-weight or uniqueness proofs themselves.

## 4. Finite verification scope

The companion `a42_mirror_canonical_check.py` reuses the established
16-state graph constructor, checks all nine extended spatial contexts,
and verifies finite polynomial identities for the cycle, syndrome,
leader, boundary, and binary-tower inverse. It performs no search for
minimum distance or lifted errors. The universal proof is §§1–3;
the checks verify its small tables and conventions.
