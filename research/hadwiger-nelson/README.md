# Hadwiger–Nelson: a machine-checkable attack

**How many colours does the plane need, so that no two points at distance exactly 1
share a colour?**

Posed around 1950. Still open. The answer, written χ(ℝ²), is known only to lie in
**{5, 6, 7}**.

| bound | value | who, when | how |
|---|---|---|---|
| lower | ≥ 4 | Nelson, 1950 | the 7-vertex Moser spindle |
| lower | ≥ 5 | Aubrey de Grey, 2018 | a 1581-vertex unit-distance graph with no 4-colouring |
| lower | ≥ 5 | Jaan Parts, 2020 | the same, down to 509 vertices |
| upper | ≤ 7 | Isbell, 1950 | a hexagonal tiling of diameter just under 1 |

Nothing has moved since.

## What this is

Lowering the upper bound to 6 would mean constructing a 6-colouring of the entire
plane — an infinite object, not something a search can produce. **Raising the lower
bound to 6 is different: it needs a single finite unit-distance graph with no proper
5-colouring.** That is a finite object. Finite objects can be searched for, and found
ones can be checked by anybody.

So that is the target here, and the whole package is built around making any claim
falsifiable by a stranger:

- coordinates are **exact** — elements of ℚ(√3, √11), never floating point, so
  "distance exactly 1" is a decidable predicate rather than a tolerance;
- colourability is decided by a **SAT solver**;
- results ship as **certificates**: a colouring re-checked against edges re-derived
  in O(n²) from the published coordinates, or a **DRAT proof verified by drat-trim**
  against a formula rebuilt from those same coordinates.

A certificate never asks you to trust this code. It asks you to run `drat-trim`.

## Status

**The machinery works and is validated. The open problem is not solved.** Concretely:

- ✅ χ(ℝ²) ≥ 4 reproduced from scratch and **machine-checked end to end** —
  `certificates/moser_spindle_no3coloring.json`, verified by drat-trim.
- ✅ The de Grey spindling argument is **automated**, and rediscovers the classic
  step unaided: given only a unit rhombus, it finds that 3 colours force the two far
  tips to match, spindles on that pair, and lands on the Moser spindle.
- ✅ A **three-copy pigeonhole variant** of the argument, which needs only a forced
  *disjunction* rather than a forced pair, implemented and validated.
- ✅ χ(ℝ²) ≥ 5 **reproduced and machine-verified.** de Grey's 1581-vertex graph
  was rebuilt from the published 39-point set S and its recipe; every count
  matches (39 → 397 → 1581, with exactly one coincident point). kissat 4.0.4
  returns UNSATISFIABLE for 4-colourability under two configurations, and
  **drat-trim verifies the proof**: `s VERIFIED`, 13 140 458 lemmas, 2 016 499 in
  core, 130 857 426 resolution steps, 522 s.
  `certificates/degrey_1581_no4coloring.json`.

  Exactly what is verified, stated precisely: *no proper 4-colouring of this
  graph assigns colours 0, 1, 2 to one pinned triangle.* One step remains outside
  the proof — a triangle's three vertices are pairwise adjacent, so any proper
  4-colouring gives them three distinct colours, and colours are interchangeable,
  so permuting them to read 0, 1, 2 loses nothing. Standard, and still not
  something drat-trim checked. The plain unbroken formula removes even that step
  and was still solving.

- ❌ χ(ℝ²) ≥ 6 — the actual goal. Not found.

- ✅ A **19-vertex, 33-edge graph with no 3-colouring**, drat-trim verified and
  **vertex-critical** — `certificates/genuine_pair_19_no3coloring.json`. Not a
  record; the Moser spindle reaches χ ≥ 4 on seven. What is new is that its
  forced pair is forced **only jointly** — neither leg is forced on its own, so
  there is no single target to spindle and the classical argument has nothing to
  grip. Built rather than found, over ℚ(√3,√11)(√v) with v = (−66 + 30√33)/256;
  that extension is forced, since v has a negative conjugate and every
  multiquadratic field is totally real.

## Pressure, and why every k = 5 search was dead on arrival

Write **pressure(p)** for the least number of colours N(p) can be squeezed into,
over all k-colourings. `hn/forced.py` computes it in at most k incremental
solves, and two theorems make it the right invariant.

**A core of size r at p requires pressure(p) ≥ k − r.** Take a colouring
attaining the minimum: at least k − pressure(p) colours are free for p while
c(T) offers at most r values, so if k − pressure(p) > r, recolour p to a free
colour outside c(T). Properness at p asks only that its colour avoid c(N(p)), so
the colouring stays proper and T was not a core.

**A k-vertex-critical graph has no core of any size, at any vertex.** Colour
G − p with k−1 colours and give p the kth: p is then the *only* vertex carrying
it, so c(p) lies outside c(T) for every T at once. de Grey's G is
5-vertex-critical, so it is completely inert at five colours — checked directly,
pivot 0 is separable from all 1520 of its non-neighbours *simultaneously*, pivot
1 from all 1572. The 880 pivots of searching in this package were answering a
question whose answer was fixed before any solver started.

Measured pressure at k = 5: **exactly 2 at all 1581 vertices of G**, the two
hubs of degree 60 included, in two seconds — and still 2 on every union tried,
up to 5533 vertices. So the smallest core available at five colours is **three**.

### The core condition is itself a pressure measurement

T is a core of p exactly when **min |c(N(p) ∪ T)| = k**. If some colouring left a
colour unused on both, recolouring p to it is proper and puts c(p) outside c(T);
conversely the colouring witnessing a failure leaves c(p) itself unused on both.
So cores can be *built forwards* instead of shrunk: ask whether a colouring
leaves the first colour free on N(p) ∪ T — by symmetry that is the whole question
— and either it proves T is a core or it hands back a colouring, and any vertex
carrying that colour elsewhere kills it when added. One solve per vertex added.

It reproduces the 19-vertex construction's pair in six solves, and it **fails** on
de Grey's G at k = 5, which is exactly what the criticality corollary demands.

### Where blocking stops

Every leg's conflict graph on the copies has maximum degree 2 — a point of a
circle is one apart from at most two points of that circle — so α ≥ m/3, and the
counting certificate Σ α(G_L) < m needs r·m/3 < m, hence **r ≤ 2**. The
19-vertex configuration sits exactly on that boundary: alphas [2, 3] against 6
copies, blocking by one.

So the two bounds meet. Pressure 2 puts every core at five colours at three or
more; counting blocks only up to two.

**A core of three can still be blocked, by misalignment.** Rotations alone never
manage it — searched exhaustively over every N ≤ 30 and every triple of angles on
the Nth roots of unity, an escape always exists, while the same search finds 298
blocking *pairs*. Reflections break it, because each leg picks up its **own**
angle in the cross-orbit conflicts: 27 344 configurations block. The smallest is
N = 3, three legs at d² = 1/3, six copies over ℚ(√3) — alphas [2,2,2] against 6
copies, exactly the counting floor, so counting says nothing while SAT and a
brute force over all 3⁶ choices agree that it blocks.

**And three is the only size that blocks at all.** Four, five and six legs were
searched the same way — forty thousand random (N, angles, shifts) draws apiece —
and not one blocks, while three blocks on the fifth draw. The counting slack says
why: with each conflict graph of maximum degree 2 the independent sets available
total r·m/3, which at r = 3 is exactly the copy count, the knife edge where
overlap can still decide the question, and at r = 4 is a third more room than
there are copies to place.

Put beside the pressure bound this **pins the size exactly**. A core of r needs
pressure ≥ k − r, so pressure 2 at five colours forces r ≥ 3; blocking forces
r ≤ 3. Any configuration that could give χ(ℝ²) ≥ 6 through a pivot and its
isometries has a core of **exactly three**, with its legs' radii and angles
matching one of the 27 344 blocking patterns. That is a far smaller target than
"a forced pair somewhere".

What is missing is the core itself, not the block.

### The machine that makes pressure, taken apart

Exactly one configuration in this package reaches pressure 3, and it is worth
naming rather than measuring. Delete everything from Sa that can go while the
pivot's neighbourhood still refuses to be squeezed into two colours, and
**47 vertices** survive — `certificates/pressure3_witness_47.json`:

- the pivot;
- its circle of 30, which splits into **five hexagons** — the 60° orbits, as the
  bipartiteness argument requires — each with exactly two alternating
  2-colourings, so five independent orientation bits;
- **sixteen** further points, at d² = 1/3, (7 ± √33)/6 and (3 ± √33)/6, each
  adjacent to exactly two circle points — and in all sixteen cases those two lie
  in **different** hexagons, so each one reads the relative orientation of a pair
  of them.

Squeeze the circle into two colours. Every one of the sixteen that sees two
differently-coloured circle points is barred from both, so it is confined to the
remaining **k − 2**. And the sixteen among themselves form a 16-vertex, 27-edge
unit-distance graph that is 3-chromatic and **not bipartite**: an odd cycle.

**The mechanism is list colouring.** Squeezing the circle into two colours hands
each auxiliary point a *list* — the colours its circle neighbours leave it. Two
differently-coloured neighbours give a list of k − 2; two of the same colour give
k − 1. The auxiliary graph never faces one palette; it faces a list assignment,
and the squeeze is refused exactly when that assignment admits no proper
colouring. Over all 32 orientations:

| | lists | refused |
|---|---|---|
| k = 4 | sizes 2 and 3 | **32 of 32** → pressure 3 |
| k = 5 | sizes 3 and 4 | **0 of 32** → pressure 2 |

The odd cycle is the special case where every list is the same pair, and it is
load-bearing where it applies — deleting three of the sixteen makes their graph
bipartite and the pressure drops from 3 to 2 at once. But it covers only half the
orientations. In the other sixteen the confined set *is* bipartite, and the
minimal refusal there is 22 vertices using **eleven** auxiliaries of which only
six are confined; the other five carry lists of size 3 and sit at d² = 1/3. A
frame with only "confined" and "free" in it cannot express that. Lists can.

### de Grey's forcing, rebuilt from one angle — and a 127-vertex gadget

Three hexagon offsets suffice where de Grey uses five. Put hexagons at
**0, θ/2 and θ** on the pivot's unit circle, with cos θ = 5/6, and take as
auxiliaries every **u + v** with u, v in different hexagons — which is every
auxiliary there is, since a point one away from u and v is the pivot reflected
across the chord uv, i.e. u + v. That is 127 points and 528 edges over
ℚ(√3, √11), and it measures **pressure 3 at four colours** — a core of one is
possible, the classical spindle regime — and pressure 2 at five.
`certificates/three_hexagon_pressure3.json`.

### The hexagons are generated by the Moser angle

Reading the offsets off the witness, the five hexagons sit at

> **0, ±θ/2, ±(60° − θ)** modulo 60°, where **cos θ = 5/6**

and cos θ = 5/6 is the Moser rotation itself: the chord it subtends on the unit
circle squares to 2 − 2cos θ = **1/3**, the split-prime quotient at 3 in
ℚ(√−11) that the spindle's angle turned out to be. Verified in the field rather
than in degrees — the cosines between hexagon 0 and the rest come out exactly
1, √33/6 twice and (5 + √33)/12 twice, with 2(√33/6)² − 1 = 5/6 the half-angle
identity and (5 + √33)/12 = cos 60° cos θ + sin 60° sin θ.

So the same angle runs through the whole construction: it is the rotation the
spindle uses, the radius its targets sit on, and the offset that spaces the
hexagons whose orientation bits carry the forcing.

Two structural facts constrain any redesign. Auxiliaries are **Minkowski sums** —
the point one away from circle points u and v is the pivot reflected across the
chord, which is just u + v — so the auxiliary set is H₁ + H₂ + …, and two of them
are one apart exactly when |Δu + Δv| = 1, an equation with a *discrete* solution
set in the hexagons' relative angle. Scanning that angle finds nothing, because
the good values are measure zero; solving it gives ten per pair, and the
half-angles only appear once there are three hexagons, since they come from
edges *between* different hexagon pairs. And the always-confined points — those
whose two circle neighbours are adjacent — are exactly the ones at √3 from the
pivot, where the chord-1 angle is 2 arcsin(1/(2√3)) = θ, which does not divide
360°, so their graph is a union of paths: **bipartite, never even an odd cycle.**
Orientation-dependence is forced, not a quirk.

**So the lift is an exact and standard question, and it is about choosability
rather than chromatic number:** at five colours the auxiliary graph must fail to
be list-colourable with lists of sizes 3 and 4. A 4-chromatic subgraph all of
whose vertices carry the *same* list of three is one sufficient way — which is
why the Moser spindle and the 19-vertex jointly forced construction above are the
pieces to try — but it is not the only way, and the non-uniform lists are where
de Grey's construction gets its extra reach.

Searched over that family: hexagons chosen from de Grey's own offsets, all
Minkowski-sum auxiliaries, every orientation. Three hexagons at 0, θ/2, θ
refuse **8 of 8** orientations at four colours — that is the 127-vertex gadget
above. At five colours, with up to five hexagons and 360 auxiliaries, **0 of 32**
orientations refuse. The lists of size 3 and 4 are simply too generous for a
depth-one design.

The scope of that negative is worth stating exactly: it covers auxiliaries
adjacent to two circle points and nothing deeper.

**The second level was built, and it does not help.** Every point one away from
two circle points is a Minkowski sum u + v, so the depth-one design is already
exhausted — no further confined point exists without enlarging the circle. Points
attached to *auxiliaries* instead have no list at all, only propagation, so only
a solver sees them: adding all 498 of them to the gadget, 625 vertices and 3324
edges, leaves the pressure at 2.

### What a core of three actually needs, and why the gadget cannot host one

Pressure 3 was never the target — the theorems pin the core size at **exactly
three**, and pressure 2 is the condition that makes three right, not an obstacle.
But the gadget is **4-colourable**, so at five colours the pivot can always take
a fifth colour nobody else has and no core exists at all, of any size. That is
the criticality corollary again, from the other side: **a core at k = 5 needs a
host that is 5-chromatic with the pivot non-critical.** The gadget supplies
pressure; it cannot supply that.

And one honest measurement about how loose the pressure bound is. The gadget is
4-chromatic — not 3-colourable — with pressure 3 at four colours, which *permits*
a core of one. Its smallest core over all 127 pivots is **seven**. Necessary is a
long way from sufficient.

### The wall, stated precisely

Two halves that have never been in the same graph:

| | pressure | can host a core at k = 5 |
|---|---|---|
| three-hexagon gadget | 3 at k = 4, 2 at k = 5 | **no** — it is 4-colourable, so the pivot always finds a fifth colour of its own |
| de Grey's G | flat 2 | **no** — it is 5-vertex-critical, so *every* vertex finds one |

Grafting was tried both ways and neither cures the other. Landing the gadget's
pivot on G's hub leaves that hub critical, because G is 5-vertex-critical and so
G − hub is 4-colourable, and the gadget minus its pivot is 4-colourable too;
the counterexample builder runs out, as it must. Placing a translated copy of G
so it *avoids* the pivot does make the pivot non-critical in principle, but over
204 placements the two never share more than a vertex or two, and near-disjoint
copies colour independently — their colours can be permuted against each other,
so no fixed target set is ever forced. Adding the whole second level, 625
vertices, leaves the gadget 4-colourable still.

So the requirement is sharp and unmet: **a 5-chromatic host, with the pivot
non-critical, substantially overlapping a pressure gadget.** Each piece exists
separately. Nothing here puts them in one graph.

### What pressure 3 actually needs

de Grey's Sa reaches **pressure 3 at four colours** — k−1, the classical spindle
regime — with cores down to 7. And it does so with **no forced pair anywhere**:
the forced-different graph on a pivot's circle is exactly its 30 unit-distance
edges, every degree 2, no odd cycle, no forced non-edge at all. The ambient graph
kills every 2-colouring of the circle *jointly*, without pinning any single pair.

That removes the object this package kept failing to find. Pressure 3 at five
colours needs no rainbow and no forced non-edge — it needs ambient rigidity, and
pressure measures exactly that, and never decreases when points are added. What
it has not yielded to is stacking: translated copies at the most frequent
difference vector (|d|² = 1/3, folding 224 vertices), a 60-degree rotation about
a centre folding 789 of 1581, and de Grey's own dihedral group applied to G — all
still flat 2, up to 5533 vertices.


## What has been ruled out so far

Recording failures is the point of a search log; these are real constraints on where
a proof can come from, not just wasted cycles.

**Plain balls are 4-colourable.** Every set of points reachable from the origin by at
most *s* unit steps drawn from {ω^a σ^m}, ω = e^{iπ/3} and σ = e^{i·arccos(5/6)}, was
found 4-colourable, up to **30 715 vertices and 315 222 edges** (radius 4). Since a
graph containing a non-4-colourable subgraph is itself non-4-colourable, this is
conclusive for those sets: they contain no 5-chromatic subgraph at all.

**No forced monochromatic pair at the origin.** The spindling argument needs a pair
(p, q) that *every* 4-colouring paints alike. Searching from the origin over all
spindle-able distances — d² ∈ {1/3, 5/9, 1, 7/3, 3, 13/3, 7, 71/9, 31/3} — found
none, in balls up to **82 357 vertices and 931 158 edges**. Forcing is monotone
under taking supergraphs, so each of these rules out every subgraph too.

**Why all of that failed — the field was wrong.** A ball of 82 357 vertices came
back 4-colourable. If de Grey's 1581-vertex graph were inside it, the ball could not
be. So it is not inside, and no amount of extra radius or compute was ever going to
put it there: the *generating set* was wrong, not too small.

Reading the construction in [de Grey 2018](https://arxiv.org/abs/1804.02385) says
why. His rotations are 2·arcsin(1/4) and 2·arcsin(1/8) — the spindles at distance 2
and 4 — plus π/2 ± arcsin(1/8). Their sines are √15/8, 3√7/32 and 3√7/8. So the
construction does **not** live in ℚ(√3, √11); it needs **√7 and √15 = √3·√5** too.

The filter that kept rotations inside ℚ(√3, √11) — the ones indexed by Eisenstein
norms — excludes precisely the two angles de Grey used. `rotation_joining(4)` and
`rotation_joining(16)` raise `ValueError` against the small field, which was correct
behaviour serving a false assumption. `hn/geometry.py` now carries `DEGREY_FIELD =
Field((3, 5, 7, 11))`, dimension 16, in which all of his rotations are exact and a
point at distance 2 rotated by 2·arcsin(1/4) moves by exactly 1. The general
multiquadratic `Field` meant this was a parameter change rather than a rewrite.

**No forced disjunction either.** The weaker three-copy hypothesis was then tested
across the six most central pivots of each ball and every spindle-able distance,
using the single-query separation test. In every case a 4-colouring existed that
separates the pivot from *all* targets at that distance at once — so no disjunction
over any subset is forced, and neither spindle applies. Balls up to 25 675 vertices
and 253 842 edges, searched this way.

Limitations worth naming rather than burying: the σ-exponent stays at |m| ≤ 1 because
|m| ≤ 2 overruns the vertex budget at these radii, pivots are drawn from the centre
outwards rather than exhaustively, and the target distances are capped at d² ≤ 40.

## Method

**The field.** Every coordinate lives in ℚ(√3, √11), represented as four rationals
over the basis {1, √3, √11, √33}. √3 builds the triangular lattice (ω = (1+i√3)/2);
√11 builds the Moser hinge (σ = (5+i√11)/6, the rotation by arccos(5/6) that carries
two points at distance √3 to distance exactly 1). Both have modulus 1, so any
product ω^a σ^m is a unit vector and sums of them are the reachable points.

Which rotations are available is a question the field answers: a spindle at squared
distance d needs sin = √(4d−1)/(2d) to be representable. Over the integers that
allows d ∈ {1, 3, 7, 19, 25, 37, 61, 91, …} — the centred hexagonal numbers
3k²+3k+1 always work, with sine (2k+1)√3/(2d) — and over the rationals more besides.

**Scale.** `Fraction` arithmetic is exact but caps out near 10⁴ points. Everything
in one search shares a denominator, so multiplying through turns generation into
int64 numpy: ~80 000 vertices and ~930 000 exactly-confirmed edges in well under a
minute. An overflow guard refuses the integer path rather than returning a wrong
answer, and `tests/test_fast_agrees.py` pins the fast path to the slow one.

**Spindling, automated.** If every k-colouring of G paints p and q alike, and ρ is
the rotation about p taking q to distance 1 from itself, then G ∪ ρ(G) has *no*
k-colouring: both copies force colour(p) onto q and ρ(q), which are adjacent. This
is why searching balls for a 5-chromatic subgraph finds nothing while spindling one
of those same balls can succeed.

**A third copy weakens what has to be forced.** Demanding an outright forced pair is
a lot to ask. Ask instead that every k-colouring tie p to *q₁ or q₂* — a disjunction,
much weaker, and therefore much likelier to be found. Three rotated copies then
suffice by pigeonhole: each copy forces one of the two, so two copies force the same
one, and their images of it share p's colour while being adjacent.

That last step needs the rotated images to be pairwise at distance 1, and the
geometry permits it at exactly one radius. Two such points fit on any circle of
radius ≥ 1/2; **three fit only on the circumcircle of a unit equilateral triangle**,
radius 1/√3, with the rotations at 120° and 240°. So d² = 1/3 is special, and cores
of size 3 or more are unusable — a fourth point pairwise at distance 1 does not exist
on a circle.

**One query per family.** Give each target q a selector asserting p and q take
disjoint colours and assume them all at once. SAT means some colouring separates p
from every target, killing every disjunction over that set in a single query rather
than O(|Q|²). UNSAT means the disjunction is forced, and the solver's UNSAT core,
shrunk by re-solving and single deletions, is the minimal forcing subset — whose size
says which spindle applies.

The same test run on a unit triangle shows exactly why k=4 is hard: with 3 colours
the triangle exhausts the palette and the centre is forced to repeat one of the
corners; with 4 colours a spare colour remains and nothing is forced at all.

## Running it

```bash
pip install python-sat numpy pytest
python -m hn.cli demo                                  # rebuild and certify chi >= 4
python -m hn.cli verify certificates/moser_spindle_no3coloring.json \
    --drat-trim /path/to/drat-trim
pytest tests/ -q

python scripts/search_forced.py        # forced pairs, two-copy spindle
python scripts/search_disjunction.py   # forced disjunctions, three-copy spindle
#   HN_K=4 searches for chi >= 5; HN_K=5 is the open problem
```

For real verification, build the checker from
[marijnheule/drat-trim](https://github.com/marijnheule/drat-trim) (`gcc -O2 -o
drat-trim drat-trim.c`). Without it, certificates still report the SAT result, marked
explicitly as unverified.

## Layout

| file | what |
|---|---|
| `hn/field.py` | exact arithmetic in ℚ(√d₁,…,√d_k) |
| `hn/geometry.py` | points, rotations, the exact distance-1 predicate |
| `hn/graph.py` | unit-distance graphs, k-core and component reductions |
| `hn/generate.py` | vertex sets: lattice balls, reachable sets, rotation families |
| `hn/fast.py` | the int64 path — generation, complete edge finding, overflow guard |
| `hn/coloring.py` | SAT encoding, k-colourability, UNSAT-core minimisation |
| `hn/spindle.py` | forced pairs and disjunctions; two- and three-copy spindles |
| `hn/certify.py` | certificate creation and independent verification |

## A gradient, and what it measured

The search was blind for most of a day. Asking a SAT solver "is this pair forced?"
returns yes or no, so a sweep is a run of independent coin flips with no sense of
getting warmer — which is why 29 930 separation queries taught us nothing beyond
"not that one either".

The solver was measuring the missing quantity all along. A pair separated with zero
conflicts is wide open; one costing thousands is nearly forced, since almost every
colouring ties those vertices together and the solver had to work to find an
exception. `SeparationDifficulty` reads it off.

**Controlled against size.** A larger formula can cost more conflicts for no reason
but its size. Adding a translate 100 units away doubles the vertex count and adds no
constraint between the copies: score 88 → 129, a factor of 1.5. A genuine tightening
of the same size reaches 4216, a factor of 48. The score tracks constraint.

**Calibrated at k=4**, where the answer is known. Nested balls of de Grey's graph:

| radius | vertices | score |
|---|---|---|
| 1.0 | 48 | 0 |
| 1.5 | 210 | 4 |
| 2.0 | 534 | 19 |
| 2.5 | 915 | **2932** |
| 3.0 | 1302 | **forced** |

It does not saturate. It sits flat, explodes by two orders of magnitude, and crosses
on the next step — at d² = 1/3, the radius the geometry and the gradient had each
singled out independently.

**Localisation was wrong.** Tightening the whole graph doubles the vertex count every
round, so the plan was to tighten only a ball around the hardest pivot — forcing being
a local question. Run side by side, the tight radius saturates while the global one
keeps climbing:

| round | radius 1.6 | radius 3.0 |
|---|---|---|
| 1 | 88 | 67 |
| 2 | 262 | 2110 |
| 3 | 2604 | 22 697 |
| 4 | 2882 | — |

The pair under test is local; what hardens it is not. So the compute wall that
localisation was meant to dodge is still there, and the dodge does not work. A blind
sweep could not have told these two runs apart at all.

## The upper bound, and why searching it is closed

Everyone repeats that lowering the upper bound needs an infinite object and so
cannot be searched. That is false for *periodic* colourings — Isbell's 7-colouring
is periodic, one hexagon repeated — and a periodic colouring is determined by a
finite fundamental domain.

So: tile by a lattice, cut the fundamental domain into cells, and forbid two cells
sharing a colour whenever a unit distance between them is *possible*, over every
lattice translation. Closed cells make that conservative, so SAT yields a genuine
colouring while UNSAT rules out only that lattice at that resolution.

It calibrates. A first attempt with cells of diameter 0.35 came back UNSAT even at
k=7, which the encoding must satisfy — too coarse, since cells of diameter d forbid
the whole band [1−d, 1+d] rather than the circle. Measuring the slack in the known
colouring fixed the resolution: with hexagons of circumradius 0.45 the same-colour
distances avoid [0.8953, 1.2256], so cells under 0.105 fit. On the index-7
sublattice with 28 cells a side, 784 cells and 78 792 constrained pairs, a
7-colouring is found in 0.7 s.

**And it cannot reach six, by a theorem rather than by a compute limit.** Cell-based
colourings are map-type colourings with polygonal regions, and those need at least
seven colours: Woodall (1973) and Townsend (1981) give six,
[arXiv:2502.01958](https://arxiv.org/abs/2502.01958) raises it to seven for maps
whose boundaries are not arcs of unit circles, with arbitrary polygons as a
corollary. Any discretisation into cells lands squarely inside that.

Which is worth knowing precisely because of the exception the theorem is careful to
state. Boundaries that *are* arcs of unit circles are not excluded. So a 6-colouring
of the plane, if one exists, cannot have polygonal colour classes — it needs curved
boundaries of radius exactly 1, or non-measurable ones. No grid will ever find it.

Recorded here so the next person does not spend the afternoon on it, as this did.

## What the literature already had

Read after the fact, which was the wrong order:

- **Forced monochromatic pairs are standard.** Parts calls them *mono-pairs* and
  builds 5-chromatic graphs from "cycles connecting two or more mono-pairs and one
  unit edge", closed by a pigeonhole argument. The core of the spindle lemma here
  is established work. What is not written down in that form is the conflict-graph
  and independent-transversal statement, which systematises an idea that already
  existed rather than introducing one.
- **The order bound is classical.** (ℤ/n)* has exponent 2 exactly when n | 24.
- **No consensus on the answer.** Polymath16 participants split between 5 and 6.
- **Scale.** Parts derived N = 6906 as a bound tied to 6-colourability, which says
  something about the size of the object being hunted.

## An unexplored field

One lead from that thread has not been followed. Philip Gibbs notes that 5/3 is
monochromatic with the origin in any homomorphic 5-colouring, since 5/3 = η² + (−η)²
in ℤ[ω₁, ω₃]. Spindling that point needs |5/3|² = 25/9, whose rotation
(34.9152°, cos 41/50) has sine 3√91/50 — so it lives in **ℚ(√7, √13)**, a field that
appears nowhere in the sources read here. Everything in the literature sits in
ℚ(√3, √11), with √5 and √7 added for de Grey's rotations.

The caveat that keeps this a lead rather than a result: Gibbs' observation concerns
*homomorphic* colourings, a restricted class, and forcing there does not imply
forcing in general.

## Slack: why the difficulty jumps where it does

A unit-distance graph in the plane has clique number **3**. Four points pairwise at
distance 1 would be a regular tetrahedron, which does not fit in two dimensions.

Define the **slack** of a colouring problem as s = k − 3, and let **f(k)** be the
number of vertices in the smallest unit-distance configuration containing a pair
forced monochromatic under k colours.

At **s = 0** a triangle exhausts the palette, so forcing is *local*: two triangles
sharing an edge already do it. f(3) = 4, which meets the trivial bound f(k) ≥ k + 1
(with n ≤ k, colour everything differently) **exactly**.

At **s ≥ 1** no local configuration exhausts anything — a spare colour always
remains — so forcing has to be assembled combinatorially across many vertices, and
f jumps. The same four vertices force nothing at k = 4, and neither does the
seven-vertex Moser spindle.

| k | slack | f(k) | settled |
|---|---|---|---|
| 3 | 0 | **4** (tight) | trivial |
| 4 | 1 | measured by `scripts/measure_fk.py` | 1961 |
| 5 | 2 | > 1581 here | 2018, 57 years later |
| 6 | 3 | — | open |

**Measured, and then corrected.** Minimising from a ball known to carry forcing
gives a 359-vertex configuration in which the pivot is forced monochromatic with
*something* at d² = 1/3. The first reading of that was that spindling it would give
a 5-chromatic unit-distance graph on ≤ 717 vertices, which is wrong: the forced
object is a **disjunction over 34 targets**, not a pair. Spindling a disjunction of
size r needs r+1 rotated images pairwise at distance 1, and a circle holds at most
three. So that configuration is not spindle-able at all, and no bound on the
5-chromatic record follows from it.

**What that correction reveals is more useful than the claim was.** Forcing at k = 4
is not absent — it is the wrong *shape*. It appears as a wide disjunction, and the
spindle needs a narrow one: core size 1, or size 2 at exactly d² = 1/3. That gap
between "forcing exists" and "forcing is usable" is the real obstruction, and it
explains why de Grey's and Parts' constructions are elaborate rather than a single
spindle of some forced pair.

**Why the core stays wide: symmetry.** On the measured configuration the forced
core is 34 of 36 targets, and all 36 form a closed orbit under the 60° rotation
about the pivot — every one maps to another. A colouring argument cannot single
out targets that a symmetry permutes, so symmetry sets a floor on how narrow a
core can get.

That rotation is not a full automorphism: 280 of 359 vertices land back in the
graph, not all of them. And that partial asymmetry is exactly what bought the two
exclusions taking 36 down to 34. **Asymmetry buys exclusions; symmetry blocks
them.** Getting to 2 needs 32 more.

Which indicts the tightening used here for most of a day. Unioning a graph with
copies of itself rotated about the *same* pivot **raises** symmetry, pushing the
core the wrong way — the effort score climbed while the quantity that decides
anything could not move. de Grey's construction is asymmetric on purpose, two
different rotations about an off-centre pivot, and this is what that is for.

**Asymmetric tightening narrows it.** Acting on that, tightening about pivots
*other* than the forcing one — breaking its symmetry instead of reinforcing it —
moves the number that symmetric tightening never touched:

| round | vertices | forced core |
|---|---|---|
| start | 359 | **34** |
| 1 | 634 | **22** |
| 2 | 987 | **13** |

Twelve exclusions, then nine, where a day of symmetric tightening held it at 34
while the effort score climbed through 88, 2110, 22 697.

**And it saturates at 11.** Carried on, the descent goes 34, 29, 24, 21, 17, 14,
12, 11 and then stops dead. The first suspicion was exhausted candidates — each
round took the eight highest-degree vertices as pivots, which after seven rounds
are the same eight — so the search was rerun with **32 candidates sampled across
the whole degree range**. Two further rounds returned 11 from all 32.

That refutes the exhaustion explanation and settles the number: **11 is the floor
of this construction, not an artefact of where it was looking.** The mechanism is
real — two thirds of the way down from 34, where symmetric tightening moved nothing
at all — and its reach stops well above the 2 a spindle needs.

**Why 11, structurally.** Two points on the circle of radius 1/√3 are at distance
1 exactly when they are 120° apart. The eleven surviving targets sit at

```
 38.55  98.55  158.55  218.55  278.55  338.55     six, step 60°
    69.6  129.6  189.6  249.6  309.6              five, step 60°
```

— two interleaved six-fold orbits, 31.05° apart. Within an orbit there are
120° pairs; **between the orbits there are none**, since 31.05 + 60k is never 120.
The two families place no constraint on each other, so an independent choice always
exists and no number of rotated copies makes the spindle lemma fire. The floor is
the orbit structure, not the copy count.

Finding that also turned up a restriction the lemma never needed. It says every
chosen image carries the pivot's colour, so *any* adjacent pair among them is a
contradiction — including one arising from two **different** targets. The
implementation compared only images of the same target. Lifting that adds 30
cross-target conflict edges against 33 same-target ones here, so it is a real
strengthening; it still does not block, for the orbit reason above.

**The mechanism is real; the metric for it is not.** The obvious measure — what
fraction of vertices the rotation about the forcing pivot still preserves — does
not predict which tightening works, and is anti-correlated in round 2: the winning
candidate (core 13) was the *most* symmetric of the five by that measure, and the
least symmetric stayed at 22. Particular positions win repeatedly; an aggregate
quantity does not explain why. Recorded as a mechanism without a metric rather
than dressed in an explanation that fits.

Scope: this is at k = 4, which is settled territory. It is a proof of mechanism,
not the prize. Applying it at k = 5 needs forcing to exist there first, and on de
Grey's graph 29 930 queries found none — there is no core to narrow until one
appears.

It also replaces the search objective. Chasing effort asks "how hard is this pair to
separate"; the quantity that actually has to move is **the size of the minimal
forced core**, which must reach 2. Here it sits at 34.

The gap f(k) − (k+1) is what the plane's missing K₄ costs, and the history of the
problem is that gap widening. The anchors are pinned in `tests/test_slack.py`; the
clique bound is elementary and f(3) = 4 is immediate, so what is worth anything
here is f(4) as an actual number, which decides whether a search reaching ten
thousand vertices is close to s = 2 or missing orders of magnitude.

Forcing and chromatic number are not the same thing — the Moser spindle is
4-chromatic while forcing nothing at k = 4 — and the framing is about f, not χ.

## Every generator in the package symmetrises

A minimal forced core has to reach 2. Asymmetric tightening floors it at 11, and
the reason turned out to be upstream of the tightening: **the search space was
symmetric by construction.** `hex_ball` is the six-fold triangular lattice by
definition. `walk_ball` steps along the sixth roots of unity. de Grey's own `Sa`
is, in his words, "all points obtained by rotating S about the origin by multiples
of 60 degrees and/or negating their y-coordinates" — a twelve-element group
applied deliberately.

That is not a flaw in those constructions; symmetry is what makes them tractable
to state and to verify. But it means every configuration examined here was a union
of orbits, and the 11 surviving core targets confirm it: they sit at 38.55 + 60k
and 69.6 + 60k, two interleaved six-fold orbits 31.05 degrees apart, with no
cross-orbit pair at 120 degrees. A core of 1 or 2 cannot be a union of six-fold
orbits. The object was excluded by the generator, not by the mathematics.

`scripts/asym_grow.py` drops the group: it accretes points at distance 1 from a
Moser-spindle seed and never symmetrises. The first version attached at random and
traded one structural defect for another — 699 vertices with 1070 edges, average
degree 3.06 against de Grey's 9.97, and a sparse graph forces nothing because it
colours with room to spare. Sampling a pool of candidate attachments and keeping
whichever lands adjacent to the most existing points fixes that: 499 vertices,
2453 edges, average degree 9.83, matched to de Grey's graph while staying
asymmetric.

Still no forcing at k = 5. But the comparison is finally a fair one — the earlier
negatives were measuring graphs too loose to force anything, not graphs too
symmetric.

## The ceiling was the field, and the number is 24

The narrowing programme needed a forced core of 2 and stalled at 11. It was
not the search that failed.

Targets at one distance from the pivot put every rotated image on a single
circle, and a point of a circle lies at unit distance from at most two of its
points, so the unit-distance graph on the distinct images has maximum degree 2
-- a disjoint union of paths and cycles. Bipartite pieces never trap, since one
side is an independent set meeting every copy, so blocking lives entirely in
the odd cycles, and an odd cycle `C_q` traps at most `(q+1)/2` targets.

Which odd cycles exist is a question about the field. A rotation of order `n`
needs `zeta_n` inside `K(i)`; for multiquadratic `K` that Galois group is an
elementary abelian 2-group, which forces `(Z/n)*` to be one too, and that
happens exactly when **n divides 24**. The odd divisors of `24 = 2^3 * 3` are 1
and 3. Over any multiquadratic field, on any circle, with any number of copies,
**the only odd cycle available is the triangle** -- and a triangle traps 2.

de Grey's field is `Q(sqrt3, sqrt5, sqrt7, sqrt11)`. Every configuration
examined here lives in a multiquadratic field. The ceiling arrived with the
arithmetic, before any searching started.

One measurement corrected the first version of this. Conflict-graph degrees
grow 4, 10, 22, 46, 94 as copies are added and look like an escape; that is
multiplicity. At squared distance 1/3, 384 images collapse onto 12 distinct
points, and among distinct points the degree is 2 at every copy count tried.
`tests/test_transversal.py` pins it.

### Niven closes it, elementarily

There is a second bound, and it is the one that actually bound every search
here. Trapping needs the images to close into a cycle, so the angle must be
commensurable with `2 pi`; and `cos t = 1 - 1/(2 d^2)` is rational exactly when
`d^2` is. **Niven's theorem** allows only `0, ±1/2, ±1` as rational cosines of
rational multiples of `pi`, so a rational squared distance admits a
finite-order rotation at just four values:

| `d^2` | angle | cycle | traps |
|------:|------:|------:|------:|
| 1     | 60°   | 6  | 0 |
| 1/2   | 90°   | 4  | 0 |
| 1/3   | 120°  | 3  | **2** |
| 1/4   | 180°  | 2  | 0 |

At every other rational `d^2` the images form a path -- bipartite -- and
nothing is trapped at all, over any field, with any number of copies.

Twelve call sites in this package skip a target whose squared distance is
irrational: four in `hn/spindle.py` and eight across `scripts/`. **Every search
run here was inside a space capped at 2 before it started**, and no choice of
field would have rescued it. The escape needs targets at an *irrational*
squared distance, `1 / (4 sin^2(pi t/n))` for an `n` with an odd divisor of at
least 5 -- which is also, independently, where the Galois bound above stops
applying.

### A hierarchy of magic circles

Closing the images into a cycle at all pins the radius. Adjacency joins images
`beta = 2 arcsin(1/2d)` apart, and that is `t` rotation steps only when
`cos beta = 1 - 1/(2d^2) = cos(2 pi t/n)`, so

    d = 1 / (2 sin(pi t / n)).

Each rotation order owns a radius, and targets anywhere else see a path and
cannot be trapped at all.

| order | radius | traps |
|------:|-------:|------:|
| 3  | 0.5774 = 1/sqrt3 | 2 |
| 5  | 0.8507 | 3 |
| 7  | 1.1524 | **4** |
| 9  | 1.4619 | 4 |
| 11 | 1.7747 | **5** |
| 13 | 2.0900 | 5 |
| 15 | 2.4049 | 4 |

`n = 3` is the circle where two points are adjacent exactly when they are 120
degrees apart -- the one every spindle argument in the literature is built on,
and the only one `n | 24` permits.

The capacities are computed, not read off the independence number. An earlier
version of this table claimed `(n+1)/2` -- 8 at order 15, 11 at order 21 -- and
that is wrong from length 9 upwards. The copies are the `n` rotations of the
cycle, so their image sets are the `n` shifts of the target set `T`, and a
colouring escapes when some independent set `S` of `C_n` meets every shift,
which happens exactly when `S - T = Z_n`. Blocking means no independent `S`
covers, and that is stricter than `S` merely being large. **Capacity does not
grow with the cycle: it peaks at 5 and comes back down**, because a longer
cycle also has larger independent sets and covering gets easier faster than
trapping does. Past length 15 the enumeration stops being affordable, and
`trapping_bound` raises rather than guessing.

`hn/cyclotomic.py` builds points where the arithmetic is different: `Z[zeta_n]`,
where `|z| = 1` forces `z conj(z) = 1` exactly, so Kronecker makes `z` a root of
unity and the unit steps are precisely the `2n` elements `±zeta_n^k`. `n = 3` is
the Eisenstein lattice everything here was built from: 6 steps, capacity 2.

It is not yet a construction. `Z[zeta_15]` gives 30 unit steps and keeps unit
triangles, so the plane's clique number survives, and a three-step walk is
3901 vertices at average degree 11.4 -- denser than de Grey's graph. It is also
**3-chromatic**, because `sqrt11` is not in `Q(zeta_15)` and so no Moser spindle
lives there. Richer arithmetic on its own buys nothing; what the theory asks for
is the compositum, a field holding both the spindles and an odd-order rotation.

## Ramification decides which magic circles exist

A unit step satisfies `u conj(u) = 1`, so at any prime `P` that complex
conjugation fixes, `2 v_P(u) = 0` and `v_P(u) = 0`. Every point reachable from
the origin by unit steps then has `v_P >= 0`, and so does every squared
distance inside one connected component.

The order-`n` magic radius squared is `1/((1 - zeta_n)(1 - zeta_n^-1))`. For
`n` a prime power that is the ramified prime above `p`, so the radius has
valuation `-2` and lies outside reach. Measured, and it is what every
disconnection all day was saying: the pivot's component held 1306 of 5017
vertices and **none** of the 11-gon, and 36481 points reached nothing one step
from a pentagon vertex.

The way out is for the prime to split with conjugation swapping the factors --
and that is exactly what the classical construction is. 3 splits in
`Q(sqrt(-11))` since `-11 = 1 mod 3` is a residue, `(1 + sqrt(-11))/2` has norm
3, and its quotient by its conjugate is `(-5 + sqrt(-11))/6`: **the Moser
rotation, up to sign**. `|1 - rho|^2 = 1/3` exactly, so `1 - rho` sits on the
classical magic circle two steps from the origin. The spindle's angle was never
chosen -- it is the split-prime quotient at 3.

`15` is not a prime power, so `1 - zeta_15` is a unit, the radius is an
algebraic integer, and nothing forbids it. Counting points on each circle in
one ball, at depth 3:

| order | prime power? | points found | capacity |
|------:|:------------:|-------------:|---------:|
| 3  | yes, ramified | 4  | 2 |
| 5  | yes, ramified | **0** | 3 |
| 15 | no, unit      | **30** | **4** |

### How close it gets, and where it stops

Over `Q(zeta_15)(sqrt(-11))` -- degree 16, the smallest field with `omega` for
the triangles, `zeta_5` for the rotations, `sqrt5` for the radius and
`sqrt(-11)` for the spindle -- the order-15 circle carries two full orbits of
the order-15 rotation, and at `k = 3` the pivot's colour **is** forced onto an
orbit. Off-centre balls take the minimal core from the full orbit of 15 down to
**5**. The capacity is 4.

Five is not four, and the shortfall is structural rather than a matter of
searching harder. The core comes out as five *consecutive* positions of the
15-cycle, and an arc of five never blocks: `S = {0, 5, 10}` is independent and
meets every 5-arc. Of the 1365 four-subsets of `C_15` exactly **45 block**, in
three classes up to rotation -- gaps `(1,2,10,2)`, `(2,3,4,6)`, `(2,6,4,3)` --
and none of the 45 is forced on either orbit.

### One order, a family of circles -- and the same answer on all of them

Adjacency on the circle need not be one step of the rotation. At `t` steps the
radius is `1/(2 sin(pi t/n))` and the cycle is `C_{n/gcd(n,t)}` taken in the
order `0, t, 2t, ...`, so the *same angular arc* lands on a different subset of
the cycle. At `n = 15` the family is:

| t | radius | cycle | capacity | points found |
|--:|-------:|:-----:|---------:|-------------:|
| 1 | 2.4049 | C₁₅ | 4 | 30 |
| 2 | 1.2293 | C₁₅ | 4 | 30 |
| 3 | 0.8507 | C₅  | 3 | **0** — the order-5 radius, ramified |
| 4 | 0.6728 | C₁₅ | 4 | 30 |
| 5 | 0.5774 | C₃  | 2 | 4 — this is `1/sqrt3`, the classical circle |
| 7 | 0.5028 | C₁₅ | 4 | 30 |

The theory falls out of the table twice over: `t = 5` reproduces the classical
`1/sqrt3` circle with its capacity of 2, and `t = 3` is the order-5 radius with
**no points at all**, exactly as ramification predicts.

The four capacity-4 circles have genuinely different radii and different point
sets, and **all four give the same minimal core: the angular arc `[0,1,2,3,4]`**.
Relabelling by the cycle turns that arc into `[0,1,2,8,9]`, `[0,1,4,8,12]` and
`[0,7,9,11,13]` — four different shapes, none of them forced down to four.

So the forcing is *angularly local*: the constraint reaches the pivot from its
own neighbourhood, which sees a bounded range of directions, and a blocking
subset needs its targets spread around the circle. That is a shape mismatch, not
a size one, and it does not look like something a bigger ball fixes. Bigger, in
fact, is worse: raising the cap to 120000 vertices left 20 points on the circle
and **no complete orbit**, because a truncated breadth-first frontier cuts
orbits in half.

At `k = 4` there is no forcing at all in this field: 55363 vertices at average
degree 12.6, separable. Forcing at four colours needs a construction of de
Grey's kind, not a ball.

## A whole avenue, closed with a reason

Every search here chases one binary fact, and until it turns up there is
nothing to show. There is a continuous quantity giving the same conclusion,
and every finite graph reports a value for it.

Let `S` be measurable, avoiding distance 1, of upper density `d`, and `G` a
finite unit-distance graph on `n` vertices. Average over rigid motions: the
expected size of `sigma(V) ∩ S` is `n d`, and that intersection is independent
in a copy of `G`, so it never exceeds `alpha(G)`. Hence

    m_1(R^2) <= alpha(G) / n

for every finite unit-distance graph, and five measurable classes covering the
plane force one of density at least `1/5`. So a single `G` with
`alpha/n < 1/5` would give `chi_m(R^2) >= 6`, from one independent-set
computation.

**It cannot happen.** Croft's 1967 construction is a measurable 1-avoiding set
of density about `0.2293`, and averaging runs the other way too: that forces
`alpha(G)/n >= 0.2293` for *every* finite unit-distance graph. Since
`0.2293 > 1/5`, no graph ever brings the ratio under `0.2`, and the most this
argument can yield is `1/0.2293 = 4.36` — `chi_m >= 5`, which is already known.

### Every fractional method is blind to this problem

The same set closes more than one avenue. Translate Croft's independent set of
density `m_1` over the isometry group: the copies cover every point equally, so
weighting them is a *fractional colouring* of total weight `1/m_1`. Hence

    chi_f(R^2) <= 1 / m_1 <= 1 / 0.2293 = 4.36 < 5.

The fractional chromatic number of the plane is **below five**. Every LP and
SDP relaxation is bounded by `chi_f`, so no relaxation can see `chi >= 5`, let
alone `chi >= 6`. That is why de Grey's 2018 result had to be an integral
combinatorial argument, and why the LP line had stalled for decades without it
being a shortage of computation. **The only route to six is integral.**

### And what that says about the object being searched for

Turned around, the same bound is a fact about the target rather than about the
method. `alpha(G)/n < 1/5` would give `chi(G) >= 6` outright, with no measure
theory at all: a 5-colouring splits `n` vertices into five independent sets and
one of them has at least `n/5`. Croft's bound says no unit-distance graph ever
qualifies.

So **a 6-chromatic unit-distance graph has independent sets of at least
`0.2293 n`**, and five of them cover at least `1.1465 n` -- more than the whole
graph, with room to spare. Whatever stops such a graph from 5-colouring, it is
never counting. It has to be structure, and the classes have to fail to fit
together for reasons no cardinality argument can see.

Recorded rather than dropped, because knowing why an avenue closes is worth as
much as a search that fails quietly inside it. The ceiling here is a
construction, not a shortage of computation, and `tests/test_density.py` pins
it. The Moser spindle's `2/7 = 0.2857` and the published `m_1 <= 0.2470`
(Ambrus, Csiszárik, Matolcsi, Varga, Zsámboki 2023, by Fourier methods rather
than from a graph) sit between the two.

## Honest odds

Polymath16 worked on this for years. The chance that this finds a 6-chromatic
unit-distance graph is small. What it does provide is a correct, fast, fully
certifying search whose negative results are recorded precisely enough to be worth
something on their own — and against which any future claim, from anyone, can be
checked in one command.
