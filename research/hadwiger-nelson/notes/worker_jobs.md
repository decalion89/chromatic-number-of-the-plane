# Worker jobs: the Exoo–Ismailescu route at a repulsive distance

A self-contained brief for helper sessions. Nothing here proves `χ(ℝ²) ≥ 6`.
It says what would, and how to check it.

## The reduction

Fix a distance `d`. Suppose:

- **(W) a witness.** A finite graph `W` in the plane, with edges of length 1 and
  `d`, that is not 5-colourable. That is, `χ(ℝ², {1, d}) ≥ 6`.
- **(H) a gadget.** A finite unit-distance graph `H` with two vertices `A, B`
  at distance `d` that get different colours in every 5-colouring of `H`.

Then put a congruent copy of `H` on every `d`-edge of `W`. The union is a
unit-distance graph, and it is not 5-colourable, so `χ(ℝ²) ≥ 6`.
- `W` and `H` may live in different number fields.
- The copies are placed by Euclidean motions.
- Extra coincidental edges in the union only help.

**What is known for (W).** `χ(ℝ², {1, d}) ≥ 6` holds for:

| `d` | who | size |
|---|---|---|
| `(1+√5)/2` | Huddleston; Parts | 31 vertices ([arXiv 2010.12656](https://arxiv.org/abs/2010.12656)) |
| `2` | Exoo–Ismailescu | 426 vertices ([arXiv 1909.13177](https://arxiv.org/abs/1909.13177)) |
| `√3`, `(√3+1)/√2` | Ágoston–Pálvölgyi | Polymath16 |

**Why no one finished.** Every known `d` is a *two-step* distance
`|u + v|`, with `u, v` unit vectors. Two points at such a distance share a unit
neighbour (the middle point), so 5-colourings tend to colour them **alike**.
- In our L16 seed, `P(same)` near distance 2 is 0.3–0.4.
- (H) asks for the opposite, and it is hard for exactly these `d`.

**The repulsive distance.** In the L16 seed (`data/L16_seed.json`, 6 080
points), over tabu 5-colourings:

| distance | pairs | `P(same)` |
|---|---|---|
| `2/√3` | 17 138 | **0.081**, the most repulsive |
| `0.91485` | 5 009 | 0.097 |
| `4/√3` | 2 667 | 0.100 |
| near 2 | | 0.3–0.4 (attractive) |

At `2/√3`, 5 202 pairs are split, and lie in one Kempe component, in all 12
colourings. The figure at distance 2 is 347.

In L16 `2/√3` is **not** a two-step distance: that would need
`cos θ = −1/3`, hence `√−2`. So these pairs have no common neighbour.

## Fields

`L16 = ℚ(√−3, √−7, √−11, √−247)`, with real coordinates in
`Field((3, 7, 11, 247))`.
- It has no non-split place of norm below 53, so no reduction mod a prime
  5-colours it (`notes/local_colourings.md`, `scripts/fieldscreen.py`).
- It contains both of our 5-chromatic families, glued along their carrier.

For (W), a field with `√−2` and `√−3` is richer. There the isosceles triangle
`(1, 1, 2/√3)` exists, with `e^{iθ} = (−1 + 2√−2)/3`.
- Small Minkowski sums of hexagonal stars and their `e^{±iθ}` images, up to
  931 points, are 5-colourable at once.
- The bare lattice `ℤ[ω]/√−3` is 5-colourable, periodically mod 6.
- The L16 seed's `{1, 2/√3}`-graph (54 612 edges) is 5-colourable: kissat
  takes 23 s, and tabu fails.

## Tools

- `sh scripts/worker_setup.sh` installs the python packages, builds kissat,
  drat-trim and tabu2, and prints `KISSAT=…`.
- `scripts/grow_lean.py IN OUT [R] [near_weight]` is the colouring-guided growth.
  - env `MODE=plain|apart|same`. **The mode names the goal, not the constraint.**
    - `MODE=apart` imposes `c(A) = c(B)` for the JSON's `A, B`. UNSAT means
      they are forced **apart**: the gadget (H). Use it at repulsive
      distances such as `2/√3`.
    - `MODE=same` imposes `c(A) ≠ c(B)`. UNSAT means they are forced
      **same**. That alone gives `χ(ℝ²) ≥ 6`, with no witness needed: rotate a
      copy about `A` so that `B` moves by exactly 1. Use it at attractive
      distances.
  - env `DIST2=4/3 DIST2_GEN=3` adds the second distance `2/√3`, and grows a
    `{1, 2/√3}`-graph.
  - env `KISSAT`, `KTIME`.
  - It finds every edge exactly and checkpoints to `OUT` every 10 iterations.
- `scripts/backbone.py graph.json [rounds]` looks for forced pairs.
  - env `DIST2` scans a `{1, d}`-graph. A forced SAME pair `(x, y)` there,
    rotated about `x` so that `y` moves by 1 or `d`, gives a witness (W)
    (Parts' route).
  - env `APART_D2=4/3` tests pairs at `2/√3` for being forced apart.
- `scripts/verify_pair.py graph.json udg|apart|same|two [--d2 4/3] --kissat K --drat-trim D`
  checks a claim; the kind matches the growth's `MODE`. It
  rebuilds everything exactly, pins nothing, and asks three pysat solvers,
  plus kissat with a DRAT proof checked by drat-trim.

**Rule: an UNSAT is a claim, never a result, until `verify_pair.py` agrees
with all solvers and DRAT.**
