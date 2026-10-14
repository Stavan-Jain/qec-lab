# Uniform stability of the original canonical syndrome leader

2026-09-08. Analytic continuation of A42 S7–S8. Actual Laurent x
coordinates are retained. No period or support search is used.

## 1. Statement and significance

For any positive q with `3` not dividing q, work in
`(F2[t]/(t^q+1))[x^+/-1]` and put

    P=x^3+t+t^2, Q=1+(x+x^2)t^3, N=1+t+...+t^(q-1),
    D(f,g)=Pf+Qg, e_star=(x^-3 N,0).

Let S consist of the whole first-block column zero and whole
second-block columns zero and one. Suppose a finite-x chain e
satisfies `D e=N`, and define

    delta=|e|-q, k=|supp(e) intersect S|.

**Theorem 1.1.** For every such solution,

    k <= 17 delta,
    |e+e_star| <= 55 delta.                            (1.1)

The second inequality is a physical Hamming-stability theorem:
the number of sites differing from the unique canonical leader
is controlled by the weight excess, independently of the cylinder
height and unbounded x support. It does not assert that a chosen
generator has bounded width or that normalization is nonexpansive.

The desired canonical hiding-mass inequality would require the
much sharper `k <= (3/2) delta`. That estimate remains open.
The new result supplies a uniform finite-excess reduction rather
than deriving the sharp inequality from uniqueness alone.

## 2. Marked rows and exact excess

Use the 16-state, 128-edge graph and potential of
`A42_S7_spatial_equality.md` section 2. For each actual row j,
let `L_j` be its physical weight minus the weight after parity
aggregation of x modulo three, separately in the two blocks.
Every `L_j` is an even nonnegative integer. Let `s_j` be its
graph-edge slack. These are nonnegative integers, and summing
the graph potential around the row circle gives exactly

    delta=L+sigma, L=sum_j L_j, sigma=sum_j s_j.         (2.1)

Mark a row if `L_j>0` or `s_j>0`, and let B count marked rows.
Then

    B <= L/2+sigma <= delta.                           (2.2)

On an unmarked, or clean, row every occupied residue lifts to
exactly one actual site, and the row is one of S7's seventeen
zero-slack edges. Inspection of that complete table gives two
useful bounds: a clean row has at most **one** site in S, and
at most two physical sites in total. Its distance from the
one-site row of e_star is therefore at most three.

If B=0, then delta=0 and S7's all-height uniqueness theorem
gives e=e_star, proving (1.1). Otherwise the unmarked rows
split into at most B maximal cyclic clean intervals. We bound
each interval uniformly, without assuming that its endpoints
lie in recurrent components of the zero-slack graph.

## 3. Removing the troublesome transition locally

Use S7's state names A,B,C,D,E,F; the number B of marked rows
is distinguished by context from the graph state named B.
Consider a clean interval of n rows numbered zero through n-1.
The edge `D->C` cannot occur at an index j satisfying

    1 <= j <= n-4.                                    (3.1)

Indeed its clean predecessor is one of `D,E,F -> D`, all with
empty second block. The next edge is necessarily `C->A`.
The following edge must be `A->B`: the only other clean option
ends at `(0,0,1)`, which has no clean outgoing edge, contradicting
the required row j+3. At j+3 the choices are `B->C,E,F`.
Thus the full five-edge context of S7's spatial obstruction is
present, including its endpoint hypotheses.

Writing the actual row equation as

    x^3 f_i+f_(i-1)+f_(i-2)+g_i+(x+x^2)g_(i-3)=1,

the equation at j+2 forces `f_(j+2)=x^-3`; its residue-zero
part at j+3 instead forces `f_(j+2)=1`. These are distinct
Laurent monomials. This proves (3.1), without a cyclic-walk
or recurrent-state assumption on the clean interval.

Discard its first row and last three rows for bookkeeping;
if n<=4, discard the whole interval. The discarded rows
contribute at most four holes and twelve sites of Hamming
difference. The remaining path has no `D->C` edge.

## 4. Complete path inventory after that removal

The complete zero-slack graph without `D->C` has the following
possibilities. They are an exact finite graph inventory, not a
supposition that only small recurrent components exist.

### 4.1 The A-B-C segment

Its directed cycle is `A->B->C->A`. Any segment with at least
eight edges contains the six consecutive edges
`A->B->C->A->B->C->A`. The unaggregated row equations for
that word are contradictory: at its row three they give
`x^3 f_3=1`, while at row five they give `f_3=1`.
The needed delayed second-block row is inside the six-row word.
This is exactly the local argument in S8 section 5.3.

Consequently such a segment has at most seven edges. The
successive upper bounds on holes are `1,1,0`, and on distance
from the canonical row they are `2,3,1`. In any of the three
phases a segment of length at most seven therefore has at most
five holes and Hamming difference at most fifteen.

### 4.2 Exits, transients, and the D loop

The main path can leave that cycle via `B->E` or `B->F`,
followed either by an edge to D or by a terminal edge. These
at most two edges add at most one hole and Hamming difference
at most five. The direct terminal edge `A->(0,0,1)` has no
larger cost. A path starting partway through an exit obeys
the same bounds. Once at D, only `D->D` remains.

On a run of r D-loop rows, write `f_i=M_i` with each M_i a
monomial in `X=x^3`; all g_i vanish. For i>=3 the row equation
is

    X M_i+M_(i-1)+M_(i-2)=1.                           (4.1)

Thus S8's monomial-word lemma applies to the subword
`M_1,...,M_(r-1)`, not necessarily to the entire run. That
subword has at most three entries equal to one. Allowing the
first row separately gives at most four holes on the run.

For Hamming difference, if r<=8 the bound is `2r<=16`.
If r>=9 the subword has length at least eight; the same lemma
forces every one of its entries to be `X^-1`, exactly the
canonical site. Only the first row may differ, by at most two.
Hence sixteen is a uniform bound for a D-loop run of any length.
This explicitly retains the delayed-g endpoint issue.

The only other transient path in the graph is

    (1,0,1) -> (0,1,0) -> (1,0,2) or (1,0,3),

ending at a state with no outgoing clean edge. It has at most
one hole and Hamming difference at most four. There are no
unlisted paths joining two different D-loop runs or returning
from a terminal state.

### 4.3 The isolated two-cycle

The remaining recurrent component alternates `(0,0,omega)`
and `(0,0,omega^2)`. Here f is zero and g is one actual
monomial per row. S8 section 5.2 applies to this actual clean
interval: after the first three rows the equation is
`g_i+(x+x^2)g_(i-3)=1`. Along each row class modulo three
the only possible chain of monomials is a subword of

    x^-2 -> x^-1 -> x.

The interval has at most nine rows and at most three holes.
Its Hamming difference from the canonical first-block sites
is at most eighteen. This component has no clean edges to
the preceding components or transient paths.

### 4.4 Uniform interval totals

The main path bounds dominate the other cases. Before restoring
the four discarded rows, they are

    holes <= 5+1+4=10,
    Hamming difference <= 15+5+16=36.

After restoring those rows, every maximal clean interval obeys

    holes <=14, Hamming difference <=48.               (4.2)

The graph inventory and these small integer maxima are checked
by `scripts/ibm_canonical_stability_check.py`. Spatial recurrence
and clean-word exclusions, rather than the graph alone, supply
the length-independent bounds for its recurrent components.

## 5. Summation proves stability

Every marked row has at most three sites in S. Its aggregated
weight is at most six, so its physical weight is at most
`6+L_j`, and its Hamming difference from the canonical row is
at most `7+L_j`. Summing these costs and (4.2) over at most B
clean intervals gives

    k <=3B+14B <=17 delta,
    |e+e_star| <=7B+L+48B
                 <=55 sigma+(55/2+1)L
                 <=55 delta.

This proves Theorem 1.1 for all stated heights and all finite
x supports. There is no weight-minimality, connectedness,
support-lifting, or generator-width hypothesis.

## 6. Consequence for the fixed canonical mixed lift

For an actual canonical binary boundary dressing, the mixed
physical weight is `W=6q+3 delta-2k`. Besides (1.1) there is
the trivial bound `k<=|e|=q+delta`. Splitting at delta=q/16
therefore gives

    W >= 65q/16.                                      (6.1)

Indeed below that split use `W>=6q-31 delta`; above it use
`W>=4q+delta`. For a noncanonical dressing the earlier parity
and uniqueness bound `W>=4q+2` remains valid too. Physical W
is even, so the combined bound can be stated as

    W >= max(4q+2, 2 ceil(65q/32)).                     (6.2)

The canonical dressing attains 6q: explicitly,
`e_star+1_S=delta_boundary(x^-3(1+x)N)`. These are improved lower
bounds for one fixed-support cylinder stratum, **not** the
quantum distance of the full IBM member or even a proof of
full canonical HM. In particular no exact-distance goal is
completed by (6.1).

More useful for the sharp problem is the uniform structural
bound: at excess two, a possible counterexample to the desired
next-shell inequality has at most 110 differing sites, at
every height. This does not assert completeness of a spatial
window or authorize searching all subsets of that size. Sharp
local accounting and the actual syndrome equations still need
to reduce that bound to at most three holes.

## 7. Verification scope

The fixed diagnostic checks the seventeen clean edges, all
possible transitions at the local obstruction, the full graph
after deleting `D->C`, the three phases of the at-most-seven-edge
triangle segment, and the constant arithmetic. It does not
enumerate cylinder configurations or infer a theorem from
sampled heights. Resource limits are 55 seconds wall, 50 CPU
seconds, and 768 MiB. The analytic proof uses the already-banked
and explicitly scoped S7 and S8 monomial lemmas.

The parent diagnostic passed in 0.00016 seconds at 17.9 MB peak
RSS. An independent read-only audit explicitly passed sections
2–5, including all seventeen zero-slack edges, the endpoint
trimming, triangle bound, D-loop subword, and both constants.
The connected-defect and large-height lifting argument below
also passed independent review. No finite support search was
used in either the proof or its review.

## 8. Connected defects and a complete fixed-excess reduction

This section restricts to dyadic q. In this case P is a unit:
with `u=t+t^2`, one has `u^q=0` and

    P^-1=sum_(j=0)^(q-1) x^(-3(j+1)) u^j.

Thus every finite-x binary cycle is a boundary. This is S6
Proposition 3.1, not an additional acyclicity assumption.

### 8.1 Each nonzero connected defect costs at least two

Put `s=e+e_star`. Make a graph on its occupied physical qubits,
joining two when their columns in D share a syndrome coordinate.
Each connected component c is separately a cycle: every syndrome
coordinate touched by it receives all of its contributions from
that same component, and the complete syndrome is zero.

Consequently `e_star+c` is another solution with syndrome N.
S7's uniqueness and minimum theorem imply

    delta_c=|c|-2|c intersect e_star| > 0.

Every binary cycle has even weight, by augmentation of
`P c_first+Q c_second=0`, since both P and Q augment to one.
Hence each positive delta_c is at least two. Distinct components
have disjoint support, so

    delta=sum_c delta_c, k=sum_c |c intersect S|.        (8.1)

There are at most delta/2 nonzero components. In particular an
excess-two error differs from the leader in exactly one connected
defect. This conclusion follows from strict leader uniqueness,
not from the generally false assertion that every class-minimal
logical must be connected.

### 8.2 An offending component must meet the leader line

If the sharp canonical inequality fails, then
`sum_c (3 delta_c-2|c intersect S|)<0`. At least one component
has negative summand. It is itself a canonical counterexample
via `e_star+c`, because binary cycles are boundaries here.

Such a component must intersect e_star. Otherwise delta_c=|c|
and `3 delta_c-2|c intersect S|>=|c|>0`. Theorem 1.1 applied
to `e_star+c` also gives

    |c| <=55 delta_c.                                 (8.2)

Fix a positive excess cap Delta. Any violation with
`delta<=Delta` therefore has a connected offending component
of at most `55 Delta` sites, anchored on the line x=-3.

### 8.3 A proved spatial window and height cutoff

Every monomial of P and Q has each coordinate between zero
and three. Adjacent vertices of the syndrome-interaction graph
therefore differ in x by at most three, and their row difference
has a representative between minus three and three. A connected
n-site defect meeting x=-3 is consequently contained in

    -3-3(n-1) <= x <= -3+3(n-1).                       (8.3)

If `q>3n`, lift its row coordinates along a spanning tree.
Any two lifted vertices differ by at most `3(n-1)` along
their tree path, so the lifted support has row span at most
`3(n-1)`. Multiplication by P or Q increases that span by
at most three. The lifted syndrome therefore has row span
at most `3n<q`. Its reduction modulo the row period is zero;
reduction is injective on that row interval. Hence the lifted
defect is an actual finite plane cycle. This argument does
not assume that a class lift gives an injective support lift.

Translate this finite plane support vertically and project it
to any dyadic height strictly greater than 3n. Its weight,
leader overlap, S overlap, zero syndrome, delta_c and negative
charge are all preserved. It is again a boundary because P
is a unit at the target height.

Define `Q_Delta` as the smallest power of two strictly greater
than `165 Delta`. A counterexample with excess at most Delta
exists at some dyadic height **if and only if** one exists
at one of the dyadic heights at most Q_Delta, with its
offending component contained in the actual x interval

    [-165 Delta, 165 Delta-6].                        (8.4)

The forward direction uses the original height if it is at
most `165 Delta`; otherwise (8.2) and the preceding lift
permit transplantation to Q_Delta. The reverse direction
is immediate. The excess of the selected component may be
smaller than that of the original counterexample, but remains
positive and no larger than Delta.

At Delta=2 the exact finite list is

    q in {1,2,4,8,16,32,64,128,256,512},
    x in [-330,324], |c|<=110.                        (8.5)

This is a **proved completeness reduction for bounded excess**.
It is not a claim that enumerating that large finite set is
cheap; no such enumeration has been run or proposed. It can
instead constrain a hand classification of the next shell.
The full HM question allows unbounded Delta, so (8.4) is not
a single finite gate proving full HM and does not complete
the original exact-distance objective.

### 8.4 Canonical-row surgery improves the height cutoff

There is a stronger height reduction using the actual clean
interval normal form, independently of the spanning-tree lift.
In the middle path of each clean interval, a D-loop run of
length r>=9 has all but its first row equal to the canonical
row `(x^-3,0)`, by section 4.2. Shorten that run to eight rows
by retaining its first four and last four rows. Every removed
row is canonical. All rows within distance three of the new
join are canonical as well. The row equation has memory three,
and is identically satisfied on these canonical rows:

    x^3*x^-3+x^-3+x^-3=1.

Thus this surgery preserves every actual syndrome equation.
It removes one unit from |e| and from q per deleted row, adds
no holes, and preserves delta and k exactly. This is surgery
on a proved homogeneous part of the actual representative;
it is not the false rule of clipping an arbitrary generator.

After all such shortenings, every clean interval has at most

    4+7+2+8=21

rows: the trimmed endpoints, triangle segment, exit, and a
now bounded D-loop. The isolated g-only component and pure
transient paths have smaller bounds. The resulting height
therefore obeys

    q_prime <= B+21B <=22 delta.                      (8.6)

The intermediate height need not be dyadic or coprime to three;
no uniqueness or invertibility assertion is applied at this
intermediate stage. The actual local equations, delta and k
are all retained directly.

If a shortening occurred, at least one of the eight-row runs
remains, with its last seven rows canonical. Insert arbitrarily
many canonical rows in its interior, keeping three canonical
rows on either side of the insertion. The same local argument
preserves all equations and both counts. In particular one can
increase q_prime to any chosen dyadic height at least q_prime.
At that final height P is a unit, so the reconstructed solution
is in the canonical boundary-dressing coset.

Let `Q'_Delta` be the least power of two at least `22 Delta`.
If an original counterexample of excess at most Delta already
has q<=Q'_Delta, no row surgery is needed. Otherwise (8.6)
forces at least one shortening and permits padding to
Q'_Delta. Combining with the connected offending-component
selection and its x bound gives the stronger complete gate

    q dyadic with q<=Q'_Delta,
    x in [-165 Delta,165 Delta-6], |c|<=55 Delta.       (8.7)

At excess two the height list in (8.5) can therefore be
replaced by **q=1,2,4,8,16,32,64**. The x window and the
110-site bound remain the same. The surgery changes neither
the negative-charge predicate nor boundary membership at the
final dyadic height. As before, this completeness result is
not a feasibility claim for enumerating the gate and does
not reduce unbounded-excess HM to finitely many tests.

An independent read-only review passed this stronger surgery
argument, including both three-row splice margins, the
21-row clean-interval count, and final dyadic acyclicity.
