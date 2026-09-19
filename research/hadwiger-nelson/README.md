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
- ❌ χ(ℝ²) ≥ 5 **not yet reproduced independently.** See below — the negative
  results are recorded rather than hidden.
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

Two limitations worth naming: the search used only the origin as pivot, and the
σ-exponent was capped at |m| ≤ 1 because |m| ≤ 2 overran the vertex budget. Both are
the obvious next levers, not conclusions.

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
k-colouring: both copies force colour(p) onto q and ρ(q), which are adjacent. One
CNF with a selector per candidate pair answers "is this pair forced?" for all pairs
at once. This is why searching balls for a 5-chromatic subgraph finds nothing while
spindling one of those same balls can succeed.

## Running it

```bash
pip install python-sat numpy pytest
python -m hn.cli demo                                  # rebuild and certify chi >= 4
python -m hn.cli verify certificates/moser_spindle_no3coloring.json \
    --drat-trim /path/to/drat-trim
pytest tests/ -q

python scripts/search_forced.py                        # HN_K=4 (chi>=5) or HN_K=5 (open)
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
| `hn/spindle.py` | forced monochromatic pairs, and the spindle union |
| `hn/certify.py` | certificate creation and independent verification |

## Honest odds

Polymath16 worked on this for years. The chance that this finds a 6-chromatic
unit-distance graph is small. What it does provide is a correct, fast, fully
certifying search whose negative results are recorded precisely enough to be worth
something on their own — and against which any future claim, from anyone, can be
checked in one command.
