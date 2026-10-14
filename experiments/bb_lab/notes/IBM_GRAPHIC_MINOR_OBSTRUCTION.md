# A local obstruction to a weight-preserving graph interpretation

2026-09-08. Original IBM pair, PR #27 research branch. This is a
finite exact certificate plus a uniform embedding argument, not a
distance bound or a refutation of the distance conjecture.

## 1. Statement and significance

**Theorem G.** On every rectangular original-pair BB torus with both
periods at least 12, the column matroid of the X-syndrome matrix has
an F7 (Fano) minor. In particular, for either IBM column and every
`r>=2`, this syndrome matrix cannot be converted to an ordinary
binary graph-incidence matrix by changing its row basis, adding or
removing redundant rows, and permuting physical qubits.

The analogous statement for Z syndromes follows by the usual block
swap and coordinate inversion. This uses CSS duality, not the false
axis-self-mirror implication discussed in A42-S6.

Consequently the especially simple route "reinterpret the same
physical qubits as graph edges and compute shortest nontrivial cycles"
is unavailable with only these operations. The theorem does **not**
exclude local Clifford maps, auxiliary qubits, more general chain maps,
weighted network descriptions, or a nonlinear normal-form theorem.
Those transformations require their own physical-weight accounting.

## 2. The fixed certificate

Write a qubit as `(block,x,y)`. In the unbounded plane its syndrome
column has the following three positions:

```
block 0: (x,y), (x-1,y), (x+1,y+3)
block 1: (x,y), (x,y-1), (x-3,y+1).
```

These are the original lab X-syndrome polynomials
`beta=1+x^-1+xy^3` and `alpha=1+y^-1+x^-3y`.

Contract the columns on this set C of 35 qubits:

```
(0,5,3) (0,3,5) (0,1,2) (0,6,1) (1,2,3)
(1,5,3) (1,4,5) (0,0,4) (1,7,6) (0,3,1)
(0,6,4) (1,3,1) (1,3,3) (0,6,2) (1,6,3)
(1,2,4) (1,4,3) (1,6,6) (0,2,3) (1,4,7)
(0,2,4) (0,4,2) (1,7,7) (0,5,1) (1,3,4)
(0,1,3) (1,2,2) (1,6,2) (0,5,0) (1,7,5)
(1,5,1) (0,2,2) (0,4,0) (1,4,6) (1,5,5).
```

Retain the following seven columns, in the displayed order, and delete
every other column:

```
(0,0,3) (0,3,4) (1,6,5) (0,3,2)
(1,3,2) (0,0,5) (1,4,2).
```

Call them `e_1,...,e_7`. Let `b_j` be the three-bit binary expansion
of the integer `j`, for `1<=j<=7`. Direct binary elimination gives

```
rank(C)=35,
rank(C union {e_j : j in J}) - 35 = rank({b_j : j in J})
```

for **all 128 subsets** J of `{1,...,7}`. In particular the union
has rank 38. The rank function after contraction is, by definition,
the rank difference on the left. Therefore the remaining seven-element
matroid is exactly F7, not merely a rank-three set with an unverified
dependency pattern.

`scripts/ibm_graphic_minor_check.py` contains the fixed lists and checks
every rank equality using integer-bit Gaussian elimination. Its default
mode performs no search. The replay passed in 0.0006 seconds, with
approximately 18 MB peak RSS. These finite equalities are the machine
certificate tier; they are not presented as a hand-derived list of
all elimination steps.

The A44 agent independently checked the certificate on both the plane
and the `12x12` torus, including all 128 subset ranks and independence
of the 35 contracted columns; see `A44_S1_codex_family.md`, §8.
The optional discovery mode now exits unsuccessfully if no minor is
found, rather than letting an inconclusive search look like a pass.

## 3. Uniform embedding and the graph contradiction

All 42 qubits lie in the patch `0<=x,y<=7`. Their syndrome positions
lie in `-3<=x<=8`, `-1<=y<=10`. Each interval contains exactly twelve
integer positions. Reduction modulo periods `ell,m>=12` is injective
on these syndrome positions and on the qubits. Thus the selected
matrix is unchanged up to row relabeling on every such torus. Delete
all other torus columns first, then use the fixed contraction. This
proves the F7-minor assertion uniformly, without enumerating members.

For completeness, F7 is not graphic by an elementary counting argument.
It has rank three, no loops, no parallel elements, and seven elements.
A graph realizing it would therefore be simple. For a simple graph
of matroid rank three, the nontrivial components have ranks `s_i`
summing to three and have at most `s_i(s_i+1)/2` edges each. Their
total number of edges is at most `3*4/2=6`, a contradiction. Graph
edge deletion and contraction realize the corresponding matroid
operations; hence a graphic matroid cannot have this minor.

Invertible row operations preserve column dependencies. Adding or
removing redundant rows does the same, and permuting qubits only
relabels the columns. None can remove this obstruction. This proves
Theorem G with the precise transformation scope in §1.

## 4. Discovery record and what remains open

An initial fixed-basis projection on a six-by-six qubit patch found
no F7. The lexicographic basis on an eight-by-eight patch also found
none. These incomplete tests did not establish graphicness. A seeded
change of column order (`random.Random(0)`) produced a different basis
and found the minor in about 2.3 seconds. Removing contraction rows
whose residual seven-bit pattern was already in the F7 row span
reduced its contraction set from 109 to 35. The final all-subset rank
verification does not depend on the discovery heuristic. Optional
`--discover` reruns that bounded discovery; it is not needed to check
the proof. Aggregate measured exploration was under 15 seconds.

The surviving graph-based proof program must retain more than an
ordinary incidence description on the unchanged qubits. A promising
replacement is a physical-weight-aware network normal form with its
own completeness theorem. No such normal form or exact IBM distance
has been established by this obstruction.
