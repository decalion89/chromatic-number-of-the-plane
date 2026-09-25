# Colouring the plane over a number field through one prime

A self-contained account of the arithmetic found on 24 September 2026. It does
not prove `χ(ℝ²) ≥ 6`. It explains why the searches for six in the field of
`five_rho7` could never succeed, and it says which fields are still in play.
Proofs of the finite-field facts are computer proofs (SAT, exact arithmetic).
Each one names the script that reproduces it.

## 1. Setting

Identify the plane with `ℂ`. A finite unit-distance graph can be moved so that
its vertices lie in a number field `K ⊂ ℂ` that is stable under complex
conjugation, with real subfield `L = K ∩ ℝ`. All the project's graphs live in
CM fields `K = ℚ(√−d₁, …, √−dₙ)`.

The **unit vectors of `K`** form the torus
`T(K) = {u ∈ K : u ū = 1}`. The **unit-distance graph `Γ(K)`** has vertex set
`K`, and `z ~ w` iff `z − w ∈ T(K)`. Every finite unit-distance graph with
vertices in `K` is a subgraph of `Γ(K)`, and `χ(ℝ²) = sup_K χ(Γ(K))`.

## 2. Reduction at a place that does not split

Let `v` be a finite place of `L`, and suppose `v` does not split in `K`: a
single place `w` of `K` lies above `v`, and `K_w` is a field, quadratic over
`L_v`. The decomposition group of `w` is all of `Gal(K/L)`, so the non-trivial
automorphism of `K_w/L_v` is complex conjugation and `N_{K_w/L_v}(u) = u ū`.
Hence the local torus `T(L_v) = {u ∈ K_w : u ū = 1}` is compact, and every
unit vector of `K` is a `w`-adic unit.

**Proposition A.**
1. **Unramified case.** If `K_w/L_v` is unramified and `L_v` has residue field
   `𝔽_q`, let `N₁ ⊂ 𝔽_{q²}` be the `q + 1` elements of norm 1, and
   `G_q = Cay(𝔽_{q²}, N₁)` the *finite plane*. Then `χ(Γ(K)) ≤ χ(G_q)`.
2. **Ramified case.** If `K_w/L_v` is ramified, `χ(Γ(K)) ≤ 3`.

*Proof.*
- **A residue map.** Let `O_w` be the ring of integers of `K_w`, `π` a
  uniformiser and `ρ : O_w → O_w/πO_w` the residue map. Fix a representative
  `rep(z)` of each coset `z + O_w`, and put `g(z) = ρ(z − rep(z))`. In the
  unramified case it lies in `O_w/πO_w = 𝔽_{q²}`.
- **Unit steps become unit steps.** If `u` is a unit vector, `z + u` lies in
  the same coset as `z`, so `g(z + u) = g(z) + ρ(u)`. The residue `ρ(u)` lies
  in `N₁`, because the norm reduces to the norm.
- **Pull back a colouring.** Composing `g` with a proper colouring of `G_q`
  therefore colours `Γ(K)` properly.
- **Ramified case.** The norm-one residues are `±1`. The same argument lands in
  `Cay(𝔽_q, {±1})`, a union of cycles. ∎

**Arcs are the linear case.** Circular colourings `⌊k·frac(φ(z))⌋`, with `φ`
additive, are exactly the colourings of `G_q` that come from a single linear
form. General colourings of `G_q` are stronger. At `q = 11` no arc works, yet
`χ(G₁₁) = 5`.

## 3. Consequences

**Theorem 1 (the Moser field).** `χ(Γ(ℚ(√−3, √−11))) = 4`.

*Proof.*
- **2 does not split.** `33 ≡ 1 (mod 8)`, so `√33 ∈ ℚ₂`. The field embeds in
  `ℚ₂(ω)`, which is unramified over `ℚ₂`, so 2 does not split.
- **Six classes of unit vectors.** The norm-one units of `ℤ₂[ω]` are, mod 4,
  exactly the six sixth roots of unity.
- **One character separates them.** `φ(α + βω) = frac₂((α + 2β)/4)` takes only
  the values 1/4, 1/2 and 3/4 on those units.
- **Lower bound.** The Moser spindle lies in the field. ∎

Checked in `hn/adelic.py`, `scripts/moser2adic.py` and
`tests/test_moser_field.py`.

**Theorem 2 (the field of `five_rho7`).**
`χ(Γ(ℚ(√−3, √−11, √−247))) = 5`.

*Proof.*
- **11 does not split.** Above 11, `L_v = ℚ₁₁(√33)` (`√741 ∈ ℚ₁₁`), and
  `K_w = L_v(√−3)` is unramified over `L_v` because −3 is a non-residue mod 11.
  So 11 does not split.
- **Upper bound.** `χ(G₁₁) = 5`: SAT, and 4 colours are UNSAT.
- **Lower bound.** `five_rho7` lies in the field. ∎

The pulled-back colouring is checked exactly, at both places above 11, on every
graph the project grew in this field: up to 32 312 points and 282 909 edges, all
with 0 monochromatic edges (`scripts/reduce11.py`, `tests/test_reduce11.py`).
So no search in this field could reach six.

**The denominator principle.** If a unit vector of `K` is not integral at
some place of `K` above a prime `p`, then some place of `L` above `p` splits in
`K`. Integrality is meant in `K`, not coordinatewise: `ω = (−1/2, √3/2)` has 2
in the denominators of its coordinates, yet it is an algebraic integer. Each
rung of the known constructions therefore had to add a rotation that is not
integral at a place over a new prime.

## 4. Finite planes

For `q ≡ 5 (mod 6)`, the case that occurs next to `√−3`, `G_q` has this
structure:
- **Triangles.** Every edge lies in exactly two triangles, there is no `K₄`,
  and `G_q` is the union of `(q+1)/6` rotated triangular tori.
- **Eigenvalues.** They are `λ(n) = Σ_c (1 − η(c² − 4n)) e(c/q)`, with `η` the
  quadratic character.
- **Ramanujan.** The graph satisfies `|λ| ≤ 2√q`; this was checked for every
  `q < 400`.

**Proposition B.** For every prime `q ≥ 53`, `χ(G_q) ≥ 6`. The same bound holds
for the deeper quotients `Cay(O_w/π^r, T mod π^r)`.

*Proof.* Hoffman gives `χ_f ≥ 1 + (q+1)/|λ_min|`. This exceeds 5 for `q > 62`
by the Ramanujan bound, and was computed for `q = 53, 59`. The deeper levels
have the same ratio: a primitive level-`j` character sums to at most `2q^{j−1}`
in modulus. ∎

The Delsarte LP gives the same bound as Hoffman, and adding the triangle
inequalities changes almost nothing. Below 53, therefore, only combinatorial
proofs help.

| `q` | `χ(G_q)` | evidence |
|---|---|---|
| 2 | 4 | the 2-adic analysis above |
| 5 | 4 | SAT |
| 11 | 5 | SAT |
| 17 | 6? | tabu finds a 6-colouring at once and no 5-colouring. Independent sets of 57 points are found easily, never 58. Since `5·57 < 289`, `α = 57` would prove `χ ≥ 6`. The SAT proof is running |
| 23 | ≥ 7? | tabu fails at 6 |
| 29, 41 | ≥ 6? | tabu fails at 5. For `q = 41` the best independent set found has 210 points, against 336 needed |
| ≥ 53 | ≥ 6 | Proposition B |

## 5. Which fields can hold a 6-chromatic graph

A `k`-chromatic graph can live only in fields with no non-split place of
local chromatic number below `k`. So the rungs demand:
- **For five:** 2 and 5 split. The Moser field fails at 2, and
  `ℚ(√−3, √−7, √−15)` fails at 5.
- **For six:** 2, 5 and 11 split, together with every `q < 53` whose finite
  plane is 5-colourable.

`scripts/fieldscreen.py` lists the non-split places of any multiquadratic CM
field:

| field | first non-split places | status for six |
|---|---|---|
| `ℚ(√−3, √−11, √−247)` | 11, 29 | dead (Theorem 2) |
| `ℚ(√−3, √−11, √−23)` | 11, 17 | dead (`χ ≤ 5`) |
| `ℚ(√−3, √−7, √−11)` | 17, 41, 83, 101 | open if `χ(G₁₇), χ(G₄₁) ≥ 6` |
| de Grey's `ℚ(√−3, √−7, √−11, √−15)` | 41, 101, 131 | open if `χ(G₄₁) ≥ 6` |
| `ℚ(√−3, √−7, √2717)` | 59, 83, 89 | open |
| `ℚ(√−3, √−7, √−11, √−247)` | 83, 173 | open |

**A 5-chromatic graph in `ℚ(√−3, √−7, √−11)`.**
1. Take the carrier's forced pair at `d² = 64/9`.
2. Compose it to `d² = 16` with the rotation `cos = 1/8`, `sin = 3√7/8`.
3. Spindle the result with `(31 + 3√−7)/32`.

It needs no `√5` and no `√247`. The graph is
`data/five_tuned_16_1_3_7_11.json`: 4 081 points and 27 242 edges, not
4-colourable by CaDiCaL or by kissat.

The field has rotations with every small prime in a denominator: 2, 3, 5, 7 and
11. The Exoo–Ismailescu `λ`-closure of the graph (378 unit vectors, rank 8)
kills every cyclic periodic colouring tried. Growth is running there. It has
passed 23 000 points, and kissat now needs 10–25 minutes per hard step.

## 6. A question

Is `χ(Γ(K))` the minimum, over the non-split places, of the local chromatic
numbers, whenever that minimum is at most `χ(ℝ²)`?

It holds in every case computed:

| field | `χ(Γ(K))` | where the minimum is attained |
|---|---|---|
| `ℚ(i)` | 2 | at 2 |
| `ℚ(√−3)` | 3 | at 3 |
| Moser | 4 | at 2 |
| `ℚ(i, √3, √11)`, i.e. the plane `ℚ(√3, √11)²` | 4 | at 2 (§8) |
| `ℚ(i, √2, √3)`, i.e. the plane `ℚ(√2, √3)²` | 4 | at 2 (§10) |
| `ℚ(√−3, √−11, √−247)` | 5 | at 11 |
| `ℚ(√−3, √−11, √−23)` | 5 | at 11 |

If it held in general, `ℚ(√−3, √−7, √−11)` would be 6-chromatic once
`χ(G₁₇) = 6` and `χ(G₄₁) ≥ 6`, and then `χ(ℝ²) ≥ 6`. We do not claim this.

The method of reducing to finite fields goes back to G. E. Moorhouse, *On the
chromatic numbers of planes* (draft, 2010). For the finite planes, see
Le Anh Vinh, *On chromatic number of unit-quadrance graphs*, arXiv
math/0510092: `√q/2 ≲ χ(G_q) ≲ q/2`.

## 7. Fields with no local obstruction (24 September, evening)

Extending `F8 = ℚ(√−3, √−7, √−11)` by one more `√−d` can remove the non-split
places 17 and 41. `scripts/fieldscreen.py` finds no non-split place of norm
below 53, and none ramified, for `F8(√−d)` with squarefree `d < 400` equal to

> 1, 2, 42, 43, 59, 66, 83, 86, 87, 103, 115, 118, 127, 154, 155, 166, 174,
> 185, 195, 203, 206, 213, 223, 230, 237, 247, 251, 254, …

In these fields Proposition B leaves no local 5-colouring at any place and any
level. The case `d = 247` is `L16 = ℚ(√−3, √−7, √−11, √−247)`:
- **Non-split places.** Below 300 there are only 83 and 173.
- **Both 5-chromatic families.** `five_tuned_16` and `five_rho7` share their
  402-point Moser-field carrier. Their union, 6 080 points and 37 474 edges,
  is `data/L16_seed.json`.
- **No circular colouring, numerically.** The seed's 471 unit directions span
  rank 12. The best `min_u ‖φ(u)‖` found over characters `φ` is 0.024, while a
  circular 5-colouring needs 0.2. The same search finds 0.204 at once in the
  Moser module. The exact MILP timed out, so this is evidence, not proof.

**Two more cautions.**
- **`F8(i)` is not automatically better.** `U ∪ iU` gives the Cartesian
  product `Γ(F8) □ Γ(F8)`, of chromatic number 5, and unit vectors mixing `F8`
  and `iF8` never close a triangle. `L16` is a better extension:
  `F8 ∩ ℚ(√−3, √−11, √−247)` is the whole Moser field, so the two halves
  share triangles and spindles.
- **Level 2 does not simply equal level 1.** A section `σ(a) = τ(a) + q·s(a)`
  of the level-2 plane over `G_q` would have to satisfy a linear system over
  `𝔽_q`. That system is inconsistent for `q = 5, 11, 17`, so `G_q` does not
  embed in level 2 that way. At `q = 17`, tabu search on level 2 (83 521
  vertices, degree 306) never improves on the lift of the best level-1
  colouring.

## 8. The plane over `ℚ(√3, √11)`: `χ = 4` (25 September)

This section answers a question that was open in print. It is not a step
towards `χ(ℝ²) ≥ 6`.

**The question.**
- Moorhouse ([draft, 2010](https://www.ericmoorhouse.org/pub/chromatic.pdf))
  noted that `ℚ(√3, √11)` is the smallest field whose plane contains a Moser
  spindle, so `χ ≥ 4`. He wrote: "We have not determined the exact value."
- Madore ([arXiv 1509.07023](https://arxiv.org/abs/1509.07023), Prop. 4.6)
  proved `4 ≤ χ(ℚ(√3, √11)²) ≤ 5`, reducing at a place over 11 and using
  `χ(𝔽₁₁²) ≤ 5`.
- Exoo and Ismailescu ([arXiv 1805.00157](https://arxiv.org/abs/1805.00157),
  DCG 2020) asked whether a 5-chromatic unit-distance graph embeds in
  `ℚ[√3, √11]²`. So did Mixon and Ismailescu in
  [Polymath16, thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/).
- Voronov, in [Polymath16, thread 17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/)
  (18 July 2021): "it seems likely that `χ(Q(i, √3, √11)) = 4` ... But as far
  as I know, nobody has proved this yet."

**Theorem 3.** `χ(ℚ(√3, √11)²) = 4`. So no 5-chromatic unit-distance graph
has all its coordinates in `ℚ(√3, √11)`.

*Proof.* Let `L = ℚ(√3, √11)`. The plane `L²` is the field `K = L(i)`, and
unit vectors are the `u ∈ K` with `u ū = 1`. Apply Proposition A.1 at a place
`v` of `L` over 2.
- **The completion.** `33 ≡ 1 (mod 8)`, so `√33 ∈ ℚ₂`, and 2 splits in
  `ℚ(√33)`. Hence `L` has two places over 2, each with completion
  `L_v = ℚ₂(√3)`. This field is ramified over `ℚ₂` with residue field `𝔽₂`.
  In it `√3` is a unit and `√3 − 1` is a uniformiser, of norm −2.
- **Inert, not split.** `i ∈ L_v` would need `−1` or `−3` to be a square in
  `ℚ₂`. Neither is, since both are `≢ 1 (mod 8)`. So
  `K_w = ℚ₂(√3)(i) = ℚ₂(√3)(√−3)`. This is unramified over `L_v`, because
  `ℚ₂(√−3)` is the unramified quadratic extension of `ℚ₂`.
- **The local ring.** `O_w = ℤ₂[√3][ω]` has residue field `𝔽₄`. Here
  `N₁ = 𝔽₄^×`, so `G₂ = Cay(𝔽₄, 𝔽₄^×) = K₄`, and `χ(Γ(K)) ≤ 4`.
- **Lower bound.** The Moser spindle lies in `L²` (Moorhouse, Prop. 1.4). ∎

**Explicitly.** Write `z = x + iy` on the `ℤ₂`-basis `1, √3, ω, √3ω` of
`O_w`, using `√11 = (s/3)√3` with `s² = 33` and `i = (1 + 2ω)√3/3`. Then
`z = a + b√3 + cω + d√3ω` with:
- `a = x_a + y_b`, `b = x_b + y_a/3`, `c = 2y_b`, `d = 2y_a/3`;
- `x = x_a + x_b√3`, where `x_a = x₁ + x₃₃s` and `x_b = x₃ + x₁₁s/3`, and
  likewise for `y`.

Colour `z` by `((a + b) mod 2, (c + d) mod 2)`, where each letter means the
units digit of the 2-adic integer part. This is `hn.adelic.q311_colour`.

**Checks.** `tests/test_q311.py` verifies:
- the decomposition `2 = P₁²P₂²` with residue fields `𝔽₂`, independently with
  sympy's `prime_decomp`;
- 810 unit vectors are 2-adic units with nonzero residue, at both places. The
  rotations include `(3 + 4i)/5` and `(√33 + 4i)/7`, which lie outside the
  Moser field.
- So are 300 random unit vectors of the form `t/t̄` (Hilbert 90).
- The colouring is proper on random unit steps, on the Moser spindle, and on
  Exoo–Ismailescu's 214-point graph (1 004 edges).
- On 638 further points and 3 012 edges (`data/ei_rho7.json`) there are no
  monochromatic edges, at either place.

**Why it was missed.** Moorhouse (Lemma 8.2, Lemma 8.4) and Madore
(Prop. 3.2, 3.8) reduce the *coordinates* `(x, y)`. That needs `x² + y²` to be
anisotropic modulo `𝔪` or `𝔪²`. At a place over 2 with `√3` this fails, since
`1² + 1² ≡ 0 (mod 𝔪²)`.
- Unit vectors such as `(−1/2, √3/2) = ω` have non-integral coordinates.
- They are integral in `O_w`, which is strictly larger than `O_v[i]`.

Reducing `z = x + iy` in `O_w` — the Hermitian form of the argument — sees the
`ω` that the coordinates hide. In the coordinates `α = x + y/√3`,
`β = 2y/√3` of `z = α + βω`, the squared distance is `α² − αβ + β²`, which is
anisotropic modulo `𝔪`. Theorem 3 is Madore's Prop. 3.2 with this form in
place of `x² + y²`; his proof applies verbatim. At odd places the two versions
agree, since `O_w = O_v[i]` there.

**Who came closest.** The ingredients are not new.
- **The reduction.** In
  [Polymath16, thread 2](https://dustingmixon.wordpress.com/2018/04/22/polymath16-second-thread-what-does-it-take-to-be-5-chromatic/#comment-4013)
  (25 April 2018), David Speyer 4-coloured the *Moser ring*
  `R = O_K[1/3]`, `K = ℚ(√−3, √−11)`. He observed `R/2R ≅ 𝔽₄ × 𝔽₄`, and that
  either projection is a proper 4-colouring. In
  [thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/)
  (3 May 2018) he gave 8-periodic variants. Philip Gibbs then reported
  (9 May 2018) that a computer search by Tamás Hubai found that all
  4-colourings of the ring have period 8.
- **The cosets.** Extending a colouring from a ring to the whole field by
  cosets is Madore's Prop. 3.2 and Moorhouse's Lemma 4.2.
- **The integrality.** Unit vectors are integral at a place that does not
  split because the norm-one torus is compact there, a standard fact.
- Dúcz ([arXiv 2606.12325](https://arxiv.org/abs/2606.12325), 2026) gave
  geometric 4-colourings of the Moser lattice and ring.

What was missing is one check: the places of `L` over 2 are inert in
`L(i) = ℚ(i, √3, √11)`, a field of degree 8 that contains the Moser field.
Then every unit vector of the plane is a 2-adic unit, not only those of the
ring, and Speyer's colouring extends to the whole plane. For an expert in local
fields this is a short observation. Its interest is that it settles questions
asked in print.

**Parts' paradox, resolved.**
- Exoo–Ismailescu's `G₄₀` forces a pair at distance 8/3 alike in every
  4-colouring *with no monochromatic pair at distance `√(11/3)`*.
- In [Polymath16, thread 13](https://dustingmixon.wordpress.com/2019/07/08/polymath16-thirteenth-thread-bumping-the-deadline/)
  (July 2019), Parts chained such pairs to get alike pairs at every distance
  `8/9ⁿ`, which sum to 1. He called this a "funny proof" that the field needs
  five colours, and asked why the colour is lost in the limit.
- Pálvölgyi pointed out that colour need not pass to the limit.
- The colouring above satisfies the hypothesis and makes the point concrete.
  Every pair at `√(11/3)` is apart, because `√33/3` is a 2-adic unit. Every
  pair at distance `8/9ⁿ` is alike, because `(8/9ⁿ)u ∈ 8·O_w`. The colouring is
  continuous for the 2-adic topology, not the real one.

**Corollaries for real fields `F`.** Proposition A at a place `v | 2` gives:
- **Ramified.** If `v` ramifies in `F(i)`, then `χ(F²) = 2`. The norm-one
  residues satisfy `ρ(u)² = 1`, so `ρ(u) = 1`, and any `𝔽₂`-linear form with
  `λ(1) = 1` 2-colours the plane. This recovers `χ(ℚ²) = 2` (Woodall).
  It also recovers Moorhouse's Theorems 7.1, 8.4 and 8.5, and Madore's
  Prop. 3.9.
- **Inert, residue field `𝔽₂`.** Then `χ(F²) ≤ 4`. For quadratic fields this
  happens exactly when `d ≡ 3 (mod 8)`. That extends Moorhouse's Theorem 8.1,
  which excluded `d ≡ 47, 59, 83 (mod 84)`, to `d ≡ 59, 83, 131 (mod 168)`.
  - For example, `3 ≤ χ(ℚ(√59)²) ≤ 4`. Reduction at 11 gave only 5
    (`59 ≡ 2² (mod 11)`).
  - The classes left open are `d ≡ 47, 143, 167 (mod 168)`.
- **Inert, larger residue fields.** The finite plane at an inert place over 2
  with residue field `𝔽_{2^f}` is `Cay(𝔽_{4^f}, μ_{2^f+1})`. SAT gives
  chromatic number 4 for `f = 1, 2, 3`. For `f = 2` it is the Clebsch graph.
- **A necessary condition for five.** Suppose a real field `F` has
  `χ(F²) ≥ 5`. Then:
  - no place of `F` over 2 ramifies in `F(i)`;
  - no place over 2 with residue field `𝔽₂`, `𝔽₄` or `𝔽₈` is inert in `F(i)`;
  - `F` has no place with residue field `𝔽₃` or `𝔽₇`, since `χ(𝔽₃²) = 3`
    and `χ(𝔽₇²) = 4`.

  Heule's `ℚ(√3, √5, √11)` and Exoo–Ismailescu's `ℚ(√3, √11, √247)` pass at 2:
  `−1` is a square in every completion over 2. There `√5` or `√247` supplies
  `i`, since `5 ≡ −3` and `247 ≡ −1 (mod 8)`.

We found no proof of Theorem 3 in any paper, in Polymath16 threads 1–18, on the
main page of the Polymath16 wiki, or by web search. The wiki's "Algebraic
formulation" page could not be retrieved. The closest work is Speyer's
colouring of the Moser ring, above. An independent check (25 September) found
no mathematical error. Nothing here has been refereed.

## 9. Split places colour modules (25 September)

§5 screens whole fields, and there only non-split places matter. At a split
place `w` the local torus is not compact, so unit vectors of `K` can have any
valuation there.

A search, however, never uses the whole field. It works in the module `M`
spanned by a finite unit set and the starting points. Suppose every generator
of `M` is integral at `w`.
- **A homomorphism.** Reduction mod `w` is a homomorphism of groups
  `M → 𝔽²`, and it sends every unit vector to the circle `a² + b² = 1`.
- **Every edge is covered.** So every unit-distance graph with vertices in
  `M` maps to `Cay(A, A ∩ circle)`, where `A` is the image of `M`. This
  includes edges that the search has not found yet.
- **The split finite plane.** At a split place with residue field `𝔽_q`,
  the target is `H_q = Cay(𝔽_q², {(a, 1/a)})`.

So a module is dead for `k` colours as soon as one integral place has a
`k`-colourable finite plane.

**Chromatic numbers of `H_q`.**

| `q` | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 29 | 37 | 49 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `χ(H_q)` | 2 | 3 | 3 | 4 | 4 | **5** | ? | ? | ? | ? | ? |
| tabu at 5 colours (best conflicts) | | | | | | | 2 | 4 | 222 | 1 290 | 4 581 |

- **Where the entries come from.** The exact values are from SAT.
- **Hoffman.** The eigenvalues are Kloosterman sums, of modulus at most
  `2√q` (Weil). So `χ(H_q) ≥ 1 + (q − 1)/(2√q)`, which is above 5 once
  `q ≥ 67`.
- **Split places of degree 2.** For every `p ≥ 11`, `H_{p²}` has
  `χ ≥ 7`.

**The two growth modules pass.** `scripts/module_gate.py` checks integrality
at each place exactly, with `p`-adic arithmetic; a coordinate can carry `p` in
its denominator and still be integral at one place over `p`. It then reduces
the growth units, the edge unit vectors and the component representatives.
- **`L16` seed** (918 growth units):
  - every place over 3, 5 and 7 has a non-integral generator;
  - at 11, 13, 17, 19, 23, 29 and 31, and at 41–61, the image is
    `𝔽_{p²}²`, which is safe by Hoffman;
  - the only small split place of degree 1 is at 37, and tabu finds no
    5-colouring of its image `H₃₇`.
- **`F8` growth** (822 units): 3 and 5 are killed. At 7 no unit has 7 in a
  denominator; the place is ramified with residue field `𝔽₄₉`, and tabu finds
  no 5-colouring of `H₄₉`.

So neither search is doomed by a reduction mod a prime, as far as tabu can
tell. What the gate does rule out: adding a rotation that removes `ρ₇` or
Moser's `σ` from `L16`'s unit set would reopen the places over 7 or 3, and
`H₇` or `H₃` would colour the module with 4 or 3 colours.

## 10. The plane over `ℚ(√2, √3)`: `χ = 4` (25 September)

Voronov, in [Polymath16, thread 17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/#comment-29291)
(18 July 2021), conjectured `χ = 4` for two fields: `ℚ(i, √3, √11)`, which is
Theorem 3, and "case (2, 3)", that is `ℚ(i, √2, √3) = ℚ(ζ₂₄)`. On 30 July
([comment 29476](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/#comment-29476))
he suggested colouring `ℤ[ζ₂₄, 1/3]` through a homomorphism to `(ℤ/4)[ζ₂₄]`,
a ring with `2¹⁶` elements, and added: "Perhaps there is a simpler way." The
argument of Theorem 3 is that simpler way: reduction modulo the place over 2,
into `𝔽₄`.

**Theorem 4.** `χ(ℚ(√2, √3)²) = 4`.

*Upper bound.* Let `L = ℚ(√2, √3)` and `K = L(i)`. Work in `ℚ₂*` modulo
squares.
1. **One place over 2.** The classes of 2, 3 and 6 (`3 ≡ −5`, `6 ≡ −10`)
   generate a group of order 4. So `L_v = ℚ₂(√2, √3)` has degree 4, and `L`
   has a single place over 2.
2. **Totally ramified.** The group misses 5 (≡ −3), the class of the
   unramified quadratic extension. So `v` is totally ramified with residue
   field `𝔽₂`, and `v(c) = v₂(N_{L/ℚ}(c))`.
3. **Not split.** The group misses −1, so `i ∉ L_v` and `v` does not split in
   `K`.
4. **Residue field 𝔽₄.** `K_w` contains `√−3 = i√3`, so `K_w = L_v(ω)`. This
   is unramified over `L_v`, with residue field `𝔽₄` and
   `O_w = O_v + O_v ω`.
5. **The colouring.** A unit vector `u` has `u ū = 1`, so `|u|_w = 1`: it is a
   `w`-adic unit, with nonzero residue. The residue of `z − rep(z)` is a
   proper 4-colouring.

Concretely, `z = x + iy = a + bω` with `a = x + y/√3` and `b = 2y/√3`, since
`i = (2ω + 1)/√3`.

**Choosing `rep(z)`.** `O_w = ℤ₂[ζ₂₄]` has the `ℤ₂`-basis `1, ζ, …, ζ⁷`
(`ζ = ζ₂₄`, with `ρ(ζ) = ρ(ω)²` since `ζ⁸ = ω`), and 2-adic fractional parts of the coordinates
in this basis give a representative. The basis `1, √2, √3, √6` of `L` does
not do this job: it is not a `ℤ₂`-basis of `O_v`, because `(√2 + √6)/2` is
integral (its square is `2 + √3`). Fractional parts in that basis are not
constant on cosets of `O_w`, and do not give a colouring.
`hn.adelic.q23_colour` avoids the choice: it colours the component of a base
point `z₀` by the residues of `a(z) − a(z₀)` and `b(z) − b(z₀)`, which are
integral there.

*Lower bound.* A 4-chromatic graph over `ℚ(√2, √3)` is already implicit in
Voronov, Neopryatnaya and Dergachev
([arXiv 2106.11824](https://arxiv.org/abs/2106.11824)). Their second series of
5-chromatic graphs starts from the 4-chromatic 10-vertex graph `L₁₀,₂` and the
unit vectors `ζ₂₄` and `(√6 + i√3)/3`. Already their set `M₂ = M₁ + M₁`
(2 593 points, 11 448 unit edges) has no proper 3-colouring
(`tests/test_q23.py`).

A shorter explicit example is a chain of three unit rhombi from `O` along
`u₁ = (1, 0)`, `u₂ = (0, 1)` and `u₃ = (−2/3 − √2/6, −2/3 + √2/6)`.
- Each rhombus forces its two tips alike in any 3-colouring.
- `|u₁ + u₂ + u₃|² = 1/3`, so the last tip `√3(u₁ + u₂ + u₃)` is at distance
  1 from `O`.

This is the closing condition of the three-rhombus chain in the research log.

The graph (`data/chain23.json`) has 10 vertices and 16 edges.
- Three pysat solvers report no 3-colouring.
- drat-trim verifies kissat's proof:
  `certificates/chain23_no3coloring.json` and
  `certificates/chain23_drat_trim_verification.txt`.

**Checks.** `tests/test_q23.py` verifies:
- the square-class facts above, and `2 = P⁴` independently with sympy's
  `prime_decomp`;
- that 312 listed unit vectors, and 300 random ones of the form `t/t̄`
  (Hilbert 90), are units with nonzero residue;
- that the colouring is proper on a 2 089-point ball of 24th roots of unity
  and on a 3 134-point graph mixing six families of unit vectors;
- the lower bound, for the rhombus chain and for the set `M₂` of Voronov,
  Neopryatnaya and Dergachev.

**What else uses this.** The criterion is general:
- **Upper bound.** A real field `L` with a place over 2 that does not split
  in `L(i)`, and whose extension has residue field `𝔽₄`, has `χ(L²) ≤ 4`.
- **Lower bound.** `√3 ∈ L` gives triangles and rhombi. A three-rhombus
  chain closing in `L` then gives `χ(L²) = 4`.

The lower bound was known, as above. We found no proof of the upper bound in
the literature.
