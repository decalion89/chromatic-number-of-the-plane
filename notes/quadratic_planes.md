# Planes over real quadratic fields that need four colours

For a field `K ⊂ ℝ`, `χ(K²)` is the least number of colours for the points of
the plane with coordinates in `K` such that two points at distance 1 get
different colours. For real quadratic fields `K = ℚ(√d)` the known values were
2 and 3. No real quadratic field was known to need four colours.

**Theorem.**
1. `χ(ℚ(√d)²) = 4` for `d = 11, 23, 35, 59, 71, 119, 131, 191, 239`.
2. `4 ≤ χ(ℚ(√47)²) ≤ 5`.

Each lower bound is a finite graph: a triangle-free unit-distance graph with
coordinates in `ℚ(√d)`, which is not 3-colourable and is vertex-critical. The
graphs have 106 vertices (`d = 11`), 686 (`d = 23`), 650 (`d = 35`), 872
(`d = 47`), 462 (`d = 59`), 674 (`d = 71`), 469 (`d = 119`), 511 (`d = 131`),
143 (`d = 191`) and 394 (`d = 239`).

So `χ(ℚ(√d)²)` is now known for every squarefree `d < 83` except `d = 47`
(§1): it is 2 if `d ≡ 1, 2 (mod 4)`, 3 if `d ≡ 3 (mod 4)` and `d ≢ 2 (mod 3)`,
and 4 for `d = 11, 23, 35, 59, 71`.

That they are not 3-colourable is a computer proof: SAT solving with a DRAT
proof checked by drat-trim, done for two separate encodings of each graph (§4).
The rest of each certificate is checked in exact integer arithmetic. Nobody
outside the project has refereed it.

## 1. What was known

Let `d ≥ 2` be squarefree and `K = ℚ(√d)`.

| `d` | `χ(K²)` | source |
|---|---|---|
| `d ≡ 1, 2 (mod 4)` | `= 2` | Johnson (1987); Moorhouse (2010), Lemma 8.4 and Theorem 8.5 |
| `d ≡ 0, 1 (mod 3)` | `≤ 3` | Fischer (1990), Theorem 9; Moorhouse, Corollary 8.3 |
| `d ≡ 0, 1, 2, 4 (mod 7)` | `≤ 4` | Moorhouse, Corollary 8.3 |
| `d ≡ 3 (mod 8)` | `≤ 4` | Fischer (1990), Theorem 10; Corollary D of hn-2adic-obstruction (2026); `notes/local_colourings.md`, Proposition A |
| `d ≡ 3 (mod 4)` | `≥ 3` | Fischer (1990), Theorem 8; Moorhouse, Theorem 8.6, for prime `d` (an odd cycle from Pell's equation) |

Madore (2015) proved `χ(ℚ(√3)²) = χ(ℚ(√7)²) = 3`. Moorhouse's Theorem 8.1,
`χ(K²) ≤ 4` unless `d ≡ 47, 59, 83 (mod 84)`, combines the rows on 3 and 7.
The public repository
[hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction)
(July 2026, not refereed; `notes/literature.md`, item 1) proves the rows on
`d ≡ 1, 2 (mod 4)` and `d ≡ 3 (mod 8)` again by reduction at 2 (its
Corollary D), and finds odd cycles over `ℚ(√11)`.

So a real quadratic field can need four colours only if `d ≡ 11 (mod 12)`:
`d ≡ 3 (mod 4)` (else 2 colours suffice) and `d ≡ 2 (mod 3)` (else 3 suffice).
The first such values are 11, 23, 35, 47, 59, 71, 83, …. For these fields the
known lower bound was 3. For `ℚ(√47)` the bounds were `3 ≤ χ ≤ 5`, and our own
experiments of 25 September had left it open (`docs/research-log.md`).

We searched Moorhouse (2010), Fischer (1990, 1994), Madore (2015), Payne (2007),
Cohen (2007), hn-2adic-obstruction, the Polymath16 threads and wiki, and the
sources listed in `notes/literature.md`, on 1 October 2026. None of them gives a
real quadratic field that needs four colours. An earlier proof may exist that we
did not find.

**A different question: `ℚ(√−d)` inside `ℂ`.** Cohen (2007) wrote that "the
natural conjecture is that `ℚ[α]` is 3-colorable for all `α` quadratic over `ℚ`",
for the set `ℚ[α] ⊂ ℂ`. For `α = √−d` this set is `{x + y√−d}`, of dimension 2
over `ℚ`; the plane `ℚ(√d)²` has dimension 4. So the theorem above does not
contradict the conjecture, and the conjecture holds, by the usual reduction at
one prime:
- for real `α`, the set lies on a line;
- for `d = 1, 2` it lies in a plane of the first row of the table;
- otherwise `d` has an odd prime factor `p`, whose prime `𝔭` in `ℚ(√−d)` is
  ramified. Every unit `z` of `ℚ(√−d)` is a `𝔭`-adic unit with `z ≡ z̄ (mod 𝔭)`,
  so `z² ≡ z z̄ = 1` and `z ≡ ±1 (mod 𝔭)`. Colour each point by its residue
  modulo `𝔭`, relative to a representative of its coset of the `𝔭`-integral
  elements, followed by a proper 3-colouring of the cycle `Cay(𝔽_p, {±1})`.

## 2. The graphs

A point is written `[a, b, c, e]` and stands for

    ((a + b√d)/D, (c + e√d)/D),

with integers `a, b, c, e` and a fixed denominator `D`. Two points are at
distance 1 exactly when the differences satisfy

    Δa² + d Δb² + Δc² + d Δe² = D²   and   Δa Δb + Δc Δe = 0.

These are two identities between integers, so every edge is checked exactly.

| `d` | `D` | vertices | edges | degrees | directions used |
|---|---|---|---|---|---|
| 11 | 30 | 106 | 246 | 3–31 | 78 |
| 23 | 120 | 686 | 1 777 | 3–55 | 102 |
| 35 | 390 | 650 | 1 711 | 3–121 | 252 |
| 47 | 240 | 872 | 2 280 | 3–76 | 102 |
| 59 | 210 | 462 | 1 163 | 3–53 | 72 |
| 71 | 120 | 674 | 1 737 | 3–90 | 132 |
| 119 | 240 | 469 | 1 210 | 3–85 | 116 |
| 131 | 390 | 511 | 1 191 | 3–82 | 134 |
| 191 | 240 | 143 | 344 | 3–40 | 52 |
| 239 | 480 | 394 | 1 007 | 3–44 | 82 |

In each graph:
- the edges are all the pairs of its points at distance 1;
- there is no triangle, and the girth is 4;
- `χ = 4`: there is a stored proper 4-colouring, and no proper 3-colouring;
- the graph is vertex-critical: for every vertex `v` there is a stored proper
  3-colouring of the graph minus `v`.

There is no triangle in any of these planes. A unit triangle would need
`√3 ∈ ℚ(√d)`, which fails for every `d` here.

## 3. How they were found

**The directions.** For a denominator `D`, take every unit vector
`((a + b√d)/D, (c + e√d)/D)` with integer coordinates. They form a finite set
`U_D`, closed under negation, under the quarter turn and under complex
conjugation. We used `|U_D| = 108` in each case (140 for `d = 71`, 148 for
`d = 119`, 132 for `d = 239`, 196 for `d = 131`, 276 for `d = 35`).

**The growth.**
1. Start from all sums of at most two vectors of `U_D`.
2. 3-colour the graph (tabucol, else kissat).
3. Add the points `p + u` whose neighbours already see all three colours.
4. Repeat until kissat answers that the graph is not 3-colourable. This took
   1 to 30 rounds and at most 12 211 points (21 895 for `d = 131`, 39 265 for
   `d = 35`). For `d = 35, 131, 239` we used a variant:
   when no candidate is blocked, it 3-colours the graph again from scratch, and
   if still none is blocked it adds candidates whose neighbours see two colours.
5. Shrink: keep the vertices whose clauses lie in the drat-trim core of a
   checked proof, then delete vertices one at a time while the rest stays not
   3-colourable.

**Why some direction sets fail.** Three places act as gates. When every unit
vector of a set is integral at such a place, reduction at that place colours
every graph the set builds.
- **At 2**, when `d ≡ 7 (mod 8)`, the place above 2 splits in `K(i)` with
  residue field `𝔽₂`. Integral unit vectors reduce to the edge `(1, 1)` of the
  graph `H₂ = Cay(𝔽₂², {(1, 1)})`, a perfect matching. With the 44 unit vectors
  `s/s̄` for `s` with coefficients in `{−1, 0, 1}`, the graphs over `ℚ(√47)`
  were bipartite.
- **At 3**, when `d ≡ 2 (mod 3)`, the place above 3 has residue field `𝔽₉` and
  splits in `K(i)`. Integral graphs map to
  `H₉ = Cay(𝔽₉², {(x, x⁻¹) : x ∈ 𝔽₉*})`, which is 3-colourable.
- **At 5**, when `d ≡ 0, 1, 4 (mod 5)`, the places above 5 have residue field
  `𝔽₅` and split in `K(i)`. Integral graphs map to `H₅`, which is
  3-colourable. For `d = 35` and `D = 174`, every direction is integral there.

The graphs `H_q` are those of §12 of `notes/local_colourings.md`. Their
chromatic number is 3 for `q = 3, 5, 9`, and at least 5 for
`q = 13, 17, 25, 29, 37, 41, 49` (kissat, this note).

When `d ≡ 11 (mod 12)`, the vector `s/s̄` with `s = 1 + i√d` is integral at
neither gate. For example:
- for `d = 11`, `s = 1 + i√11` gives `(−5 + i√11)/6`;
- for `d = 47`, `s = 7 + i√47` gives `(1 + 7i√47)/48`. Its length is 1
  because `48² − 47·7² = 1`, the Pell solution of Moorhouse's Theorem 8.6.

The gates do not explain everything. On 25 September we tried the balls of
radius 5 for the 48 unit vectors of height at most 40 over `ℚ(√47)`, up to
2 248 121 points, and they were 3-colourable (`docs/research-log.md`). That set
contains `(−23 + i√47)/24`, from `s = 1 + i√47`, which is integral at neither
gate. With all 108 unit vectors of denominator dividing 240, which include
`(1 + 7i√47)/48`, the colouring-guided growth reached a graph that is not
3-colourable in 30 rounds, at 9 139 points.

## 4. The certificates

For each `d`, `data/quadratic_planes/q{d}.json` holds the points, the edges, the
fixed edge, the 4-colouring and the 3-colourings of the vertex-deleted graphs.
`q{d}.cnf` is the formula "the graph has a proper 3-colouring in which the fixed
edge gets the colours 0 and 1". This loses no generality: rename the colours.

**Not 3-colourable**, twice for each `d`:
- `q{d}.cnf`, with the variable of vertex `v` and colour `c` at `3v + c + 1`;
- a second formula written by separate code, with the variable at `cn + v + 1`.

In both cases kissat 4.0.4 answered UNSATISFIABLE and drat-trim VERIFIED its
DRAT proof (`q{d}.logs/`), each in a few seconds.

`python3 scripts/verify_quadratic_planes.py` checks the points, the exact unit
lengths, the 4-colourings, the vertex-critical 3-colourings, that each stored
formula is exactly its graph's formula, and the upper bounds of §5. It uses no
library and runs in about a second. With
`--kissat PATH --drat-trim PATH` it also writes the second formula, solves it
and checks the proof, in about ten seconds for the five fields.
`tests/test_quadratic_planes.py` runs the fast checks.

## 5. The upper bounds

- **`d = 11, 23, 35, 71, 119, 191, 239`.** `d` is a nonzero square modulo 7
  (`d = 11, 23, 71, 191, 239`) or `d ≡ 0 (mod 7)` (`d = 35 = 5 · 7`,
  `d = 119 = 7 · 17`), and
  `7 ≡ 3 (mod 4)`. Moorhouse's Lemma 8.2 reduces the plane modulo a prime of
  norm 7 into the unit-distance graph of `𝔽₇²`. That graph is 4-colourable, so `χ(K²) ≤ 4`.
  `data/quadratic_planes/finite_planes.json` holds the colouring, and the
  checker checks it.
- **`d = 59, 131`.** `d ≡ 3 (mod 8)`, so `χ(K²) ≤ 4` by Fischer's Theorem 10 (which
  also covers `d = 11` and `d = 35`). The
  proof in `notes/local_colourings.md` reduces modulo the place above 2, whose
  residue plane is `K₄`.
- **`d = 47`.** `47 ≡ 3 = 5² (mod 11)` and `11 ≡ 3 (mod 4)`, so reduction modulo
  a prime of norm 11 gives `χ(K²) ≤ χ(𝔽₁₁²) = 5`. Madore's Lemma 4.5 prints a
  5-colouring of `𝔽₁₁²`; `finite_planes.json` holds one, checked by the
  checker. Reduction at 121 does not help (`notes/local_colourings.md`; the
  level-2 plane at 11 is not 4-colourable).

## 6. Open

- **`ℚ(√47)`.** Is `χ = 4` or `5`? It is the smallest open case of Moorhouse's
  table that this note does not settle.
- **All of `d ≡ 11 (mod 12)`.** Does every such field need four colours? For
  `d = 83, 95, 107, 155, 179, 203`, our growth has so far found only
  3-colourable graphs, with 60 to 212 directions. That may be a limit of the
  search: over `ℚ(√35)` and `ℚ(√131)` it failed in the same way with 60 to 180
  directions, and succeeded with the 276 and 196 directions of denominator
  dividing 390. Two of the
  failed direction sets fail at a gate of §3: for `d = 35` and `D = 174` every
  direction is integral at the places above 5, and for `d = 83` and `D = 410`
  at the places above 3. We do not know why the others fail.
- **Smaller witnesses.** The graph over `ℚ(√11)` has 106 vertices, and the one
  over `ℚ(√191)` 143. How small can a 4-chromatic unit-distance graph over a
  real quadratic field be?

## 7. Files

| file | contents |
|---|---|
| `data/quadratic_planes/q{d}.json` | the graph over `ℚ(√d)`: points, edges, fixed edge, 4-colouring, 3-colourings of `G − v` |
| `data/quadratic_planes/q{d}.cnf`, `q{d}.logs/` | the formula of §4 and the kissat and drat-trim logs for both encodings |
| `data/quadratic_planes/finite_planes.json` | proper colourings of the unit-distance graphs of `𝔽₇²` (4 colours) and `𝔽₁₁²` (5 colours) |
| `scripts/verify_quadratic_planes.py`, `tests/test_quadratic_planes.py` | the checks of §4 and §5 |
| `data/quadratic_planes/scripts/` | the code of the search of §3 (growth, shrinking, certification), as it was run, with a README |

## References

- D. Cohen, *The ℂ unit distance graph*, University of Chicago REU paper (2007),
  [math.uchicago.edu/~may/VIGRE/VIGRE2007/REUPapers/FINALFULL/Cohen.pdf](https://www.math.uchicago.edu/~may/VIGRE/VIGRE2007/REUPapers/FINALFULL/Cohen.pdf).
- K. G. Fischer, *Additive K-colorable extensions of the rational plane*,
  Discrete Math. 82 (1990), 181–195.
- *A 2-adic obstruction to 5-chromatic unit-distance graphs*, public repository
  [MildlyMeticulous/hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction)
  (July 2026), not refereed.
- P. D. Johnson Jr., Congr. Numer. 60 (1987), 51–58, as summarised by Payne
  (2007).
- D. A. Madore, *The Hadwiger–Nelson problem over certain fields*,
  [arXiv:1509.07023](https://arxiv.org/abs/1509.07023) (2015).
- G. E. Moorhouse, *On the chromatic numbers of planes*, draft of 3 March 2010,
  [ericmoorhouse.org/pub/chromatic.pdf](https://www.ericmoorhouse.org/pub/chromatic.pdf).
- T. H. Payne, *Unit distance graphs with ambiguous chromatic number*,
  [arXiv:0707.1177](https://arxiv.org/abs/0707.1177) (2007).
