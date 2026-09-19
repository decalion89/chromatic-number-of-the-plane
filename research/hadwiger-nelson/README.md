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

## Honest odds

Polymath16 worked on this for years. The chance that this finds a 6-chromatic
unit-distance graph is small. What it does provide is a correct, fast, fully
certifying search whose negative results are recorded precisely enough to be worth
something on their own — and against which any future claim, from anyone, can be
checked in one command.
