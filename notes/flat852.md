# A Moser-spindle-free 5-chromatic unit-distance graph on 852 vertices

J. K. Haugland builds a 5-chromatic unit-distance graph from the 84 unit
vectors of a 7-fold symmetric graph on 21 points ([arXiv 2608.04542](https://arxiv.org/abs/2608.04542)).
His graph has 2131 vertices and no Moser spindle. This note gives a graph built
from the same vectors and their mirror images, found by asking where his
construction is forced to leave its lattice:

**Theorem.** The graph `core852` of `data/flat852/core852a.json` is a
unit-distance graph on 852 points with 4 487 edges, and its chromatic number is
5. It contains no Moser spindle. Its points lie in `ℚ(ζ₂₁) ⊂ ℂ`, so their
coordinates lie in `ℚ(ζ₈₄)⁺ = ℚ(cos 2π/7, √3, √7)`, a field without `√5` and
without `√11`.

The graph is 5-colourable by a stored colouring, checked on every edge. That it
is not 4-colourable is a computer proof: SAT solving with a DRAT proof checked
by drat-trim, done twice with separate encodings (§4). Nobody outside the project
has refereed it.

As far as we know it is the smallest Moser-spindle-free 5-chromatic unit-distance
graph (§5). It is not known to be vertex-critical: no vertex-by-vertex
minimisation has been run.

## 1. The 126 directions

Let `z = ζ₂₁ = e^{2πi/21}`, `F = ℚ(z)`, of degree 12, and

    ω = ζ₃/(ζ₇ − ζ₇⁻¹) − ζ₃²/(ζ₇² − ζ₇⁻²)
      = (4 − 3z + z² + 6z³ + z⁴ − 5z⁵ + z⁶ + 2z⁷ − 7z⁸ + 4z⁹ + z¹⁰ − 2z¹¹)/7,

a root of `7x¹² − 13x⁶ + 7`, with `|ω| = 1`. Haugland's 84 unit vectors are
`D = μ₄₂ ∪ ω·μ₄₂`: the 42nd roots of unity, and their products with `ω`. They
span the lattice `L = ℤ[z] + ℤ[z]·ω` of rank 12, whose only unit vectors are the
84 of `D` (a complete enumeration, `data/flat852/background/`). The mirror image
of Haugland's graph has the vectors `D̄ = μ₄₂ ∪ ω̄·μ₄₂`. The graph of this note
uses `U = D ∪ D̄`: 126 unit vectors, closed under negation and complex
conjugation. Each vertex is stored as 12 integers `c₀ … c₁₁`, the point
`(c₀ + c₁z + … + c₁₁z¹¹)/7`.

## 2. Why the mirror: the directions of `D` alone are 4-colourable

**Proposition** (reduction at the places above 2). Every unit-distance graph
whose edge vectors all lie in `D` is 4-colourable.

*Proof sketch* (details and checks in `data/flat852/background/haugland-2adic-report.md`). The
reduction is the one of `local_colourings.md`; the same reduction at a place above 2 is
Theorem A of [MildlyMeticulous/hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction)
(`notes/literature.md`).

1. The prime 2 is unramified in `F`, and `L/2L ≅ 𝔽₆₄²`, the residues at the two
   places of `F` above 2, which complex conjugation swaps. A vector `u` of norm 1
   reduces to `(t, 1/t)` with `t ∈ 𝔽₆₄^*`.
2. The 84 vectors of `D` reduce to 42 residues, with `t` in only two of the three
   cosets of `μ₂₁` in `𝔽₆₄^*`: `μ₂₁` for `μ₄₂`, and `r·μ₂₁` for `ω·μ₄₂`.
3. There are linear maps `φ: 𝔽₆₄² → 𝔽₂²` that vanish on none of these 42
   residues; there are exactly 42 of them. For each such `φ`, `x ↦ φ(x mod 2)`
   colours every coset of `L` properly with the four elements of `𝔽₂²`. ∎

Haugland's six colourings of `L` (his Table 2) are among these linear ones, so
they extend to the whole lattice, a point his paper leaves open. The same holds
for `D̄` alone, by conjugation.

The vectors `ω̄·μ₄₂` reduce to the third coset. For `U`, the residues fill
`𝔽₆₄^*`, no linear map avoids all 63, and the residue graph
`Cay(𝔽₆₄², {(t, 1/t)})` has chromatic number at least 6 by Hoffman's bound. So
reduction at 2 no longer colours a graph that uses all three orbits of `U`. The
graph of this note uses them all: 2 252 edges along `μ₄₂`, 1 848 along `ω·μ₄₂`
and 387 along `ω̄·μ₄₂`.

The places above 3 and 7 give no 4-colouring either (the report, §§3 and 6). This
analysis told us where to look; it proves nothing about the graph, which the
certificate of §4 settles.

## 3. How the graph was found

1. **Colouring-guided growth.** The search started from Haugland's 21-point graph
   `H` and its mirror image (42 points, 84 edges). Each round took a 4-colouring
   from a SAT solver, and added the points `x + u` (`u ∈ U`) whose neighbours in
   the graph already saw all four colours.
2. **Refutation.** After about 30 rounds, at 1 023 points and 5 440 edges, kissat found
   no 4-colouring. That graph, `g1023`, is certified like `core852`.
3. **Core.** The clauses in the core of the drat-trim check of `g1023` name 850
   vertices. With the fixed triangle and one peel of the 4-core, this gives
   852 vertices and 4 487 edges: `core852`.

**Shape of the graph.**
- 758 vertices lie in the coset of `L` that contains `H`, 92 in a coset of
  `L̄`, and 2 elsewhere.
- Degrees are 4 to 34, 10.5 on average.
- The drawing fits in a disc of radius 2.53.

The code is in `data/flat852/scripts/`, and the account of the search in
`data/flat852/REPORT-growth-agent.md`.

## 4. The certificate

- **Exact edges.** Every edge joins two stored points whose difference is one of
  the 126 vectors of `U`, and `d·d̄ = 1` holds exactly in `F`. The points are
  distinct: the closest two are 0.0095 apart.
- **5-colourable.** `data/flat852/core852.5colouring.json` is proper.
- **Not 4-colourable**, twice.
  - `core852a.cnf` is the formula "the graph has a proper 4-colouring in which
    the triangle `(8, 23, 34)` gets the colours 0, 1, 2". It is sound because
    any 4-colouring can be renamed to colour that triangle so. kissat 4.0.4
    answered UNSATISFIABLE in 25 minutes, and drat-trim VERIFIED its DRAT proof
    of 954 MB in 34 minutes (`data/flat852/logs/`).
  - A second encoding, written by separate code with the variable of vertex `v`
    and colour `c` at `c·n + v + 1`: kissat UNSATISFIABLE in 19 minutes, and
    drat-trim VERIFIED the proof of 1.0 GB in 34 minutes
    (`data/flat852/independent-check/`).
- **No Moser spindle.** A spindle has two rhombi at a common vertex `A` whose
  far tips `B`, `B'` lie at distance `√3` from `A` and at distance 1 from each
  other. The rotation `(B' − A)/(B − A)` then has real part `5/6` and lies in `F`,
  so `√−11 ∈ F`. But the quadratic subfields of `F` are `ℚ(√−3)`, `ℚ(√−7)` and
  `ℚ(√21)`: `(ℤ/21)^* ≅ ℤ/2 × ℤ/6` has exactly three subgroups of index 2. So no
  graph with points in `F` contains a Moser spindle.

`python3 scripts/verify_flat852.py` checks the exact edges, the colouring, and
that the stored CNF is exactly this graph's formula, in under a second with no
library. With `--kissat PATH --drat-trim PATH` it also writes the formula again,
solves it and checks the proof, in about an hour. `tests/test_flat852.py` runs
the fast checks.

## 5. Comparison

A Moser-spindle-free 5-chromatic unit-distance graph is one that contains no
copy of the Moser spindle as a subgraph. The sizes we found:

| graph | vertices | source |
|---|---|---|
| Haugland's heptagon graph | 2 131 | [arXiv 2608.04542](https://arxiv.org/abs/2608.04542) (August 2026) |
| the earlier example cited there | 1 441 | cited in arXiv 2608.04542 |
| a graph posted on 9 September 2026 | 1 435 | [ruturajr-raval/hadwiger-nelson-spindle-free](https://github.com/ruturajr-raval/hadwiger-nelson-spindle-free/releases/tag/v0.1.0) |
| a graph posted on GitHub by zach7036 | 1 299 | found in our survey of GitHub, 28 September 2026 |
| `core852` | 852 | this note |

The smallest 5-chromatic unit-distance graph known, with spindles allowed, is
the 509-vertex graph of Parts. `core852` does not approach it.

We searched arXiv and GitHub on 28 and 29 September 2026. A smaller
spindle-free graph may exist that we did not find.

## 6. Open

- **Minimisation.** How far below 852 does vertex-by-vertex minimisation go?
- **Symmetry.** Is there a 7-fold symmetric version, grown from `H ∪ H̄` by
  orbits?
- **The field.** What is the chromatic number of the plane over
  `ℚ(ζ₈₄)⁺ = ℚ(cos 2π/7, √3, √7)`? This graph shows it is at least 5; no upper
  bound below 7 is known to us.
- **Six colours.** The same reduction asks, for six colours, whether the finite
  planes of `F` at its other small places need six colours, starting with the
  place above 3 (`G₂₇`).

## 7. Files

| file | contents |
|---|---|
| `data/flat852/core852a.json` | the graph: exact points, floats, edges with their direction, the list `U`, the fixed triangle, the checks |
| `data/flat852/core852.5colouring.json` | a proper 5-colouring |
| `data/flat852/core852a.cnf`, `logs/` | the formula of §4 and the kissat and drat-trim logs |
| `data/flat852/independent-check/` | the second encoding, its logs, and the check that wrote it |
| `data/flat852/g1023.*` | the 1 023-vertex graph that `core852` was cut from, with the same checks |
| `data/flat852/scripts/` | the code of the search: the field, the seeds, the growth, the minimisation attempts, the certification |
| `data/flat852/background/` | the report on Haugland's family at the places above 2, and the reproduction of his construction |
| `scripts/verify_flat852.py`, `tests/test_flat852.py` | the checks of §4 |
