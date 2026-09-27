# Colouring the plane over a number field through one prime

A self-contained account of the arithmetic found on 24 September 2026. It does
not prove `χ(ℝ²) ≥ 6`. It explains why the searches for six in the field of
`five_rho7` could never succeed, and it says which fields are still in play.
Proofs of the finite-field facts are computer proofs (SAT, exact arithmetic).
Each one names the script that reproduces it.

## 1. Setting

Identify the plane with `ℂ`. Every finite unit-distance graph has an embedding,
not necessarily congruent to a given one, with its vertices in a number field
`K ⊂ ℂ` that is stable under complex conjugation, with real subfield
`L = K ∩ ℝ`. Most of the project's graphs live
in CM fields `K = ℚ(√−d₁, …, √−dₙ)`.

The **unit vectors of `K`** form the torus
`T(K) = {u ∈ K : u ū = 1}`. The **unit-distance graph `Γ(K)`** has vertex set
`K`, and `z ~ w` iff `z − w ∈ T(K)`. Every finite unit-distance graph with
vertices in `K` is a subgraph of `Γ(K)`, and `χ(ℝ²) = sup_K χ(Γ(K))` by the
compactness theorem of de Bruijn and Erdős (*A colour problem for infinite graphs
and a problem in the theory of relations*, Indag. Math. 13 (1951) 371–373; cf.
[Moorhouse](https://www.ericmoorhouse.org/pub/chromatic.pdf), Thm 1.2).

## 2. Reduction at a place that does not split

Let `v` be a finite place of `L`, and suppose `v` does not split in `K`: a
single place `w` of `K` lies above `v`, and `K_w` is a field, quadratic over
`L_v`. The decomposition group of `w` is all of `Gal(K/L)`, so the non-trivial
automorphism of `K_w/L_v` is complex conjugation and `N_{K_w/L_v}(u) = u ū`.
Hence the local torus `T(L_v) = {u ∈ K_w : u ū = 1}` is compact, and every
unit vector of `K` is a `w`-adic unit.

**Proposition A (cf. [Madore](https://arxiv.org/abs/1509.07023), Prop. 3.2, ¶6.6).**
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

**Arcs are a linear case.** A circular colouring `⌊k·frac(φ(z))⌋` whose
additive `φ` factors through the residue map at `w` colours `G_q` through one
linear form `λ : 𝔽_{q²} → 𝔽_q`, each colour an arc of consecutive values of
`λ`. General colourings of `G_q` are stronger. At `q = 11` no arc works, yet
`χ(G₁₁) = 5`.

## 3. Consequences

**Theorem 1 (the Moser field).** `χ(Γ(ℚ(√−3, √−11))) = 4`.

*Proof.* Here `L = ℚ(√33)`, and `ω = e^{2πi/3}`. (`hn/adelic.py` writes
elements on `1, ζ₆` with `ζ₆ = (1 + √−3)/2` instead; the formula for `φ` below
works in both conventions.)
- **The places over 2 do not split in `K`.** `33 ≡ 1 (mod 8)`, so `√33 ∈ ℚ₂`:
  2 splits in `L`, and both places have completion `ℚ₂`. Above either of
  them `K` embeds in `ℚ₂(ω)`, which is unramified over `ℚ₂`, so neither place
  splits in `K`.
- **Six classes of unit vectors.** The norm-one units of `ℤ₂[ω]` are, mod 4,
  exactly the six sixth roots of unity.
- **One character separates them.** `φ(α + βω) = frac₂((α + 2β)/4)`, an
  additive map from `ℚ₂(ω)` to `ℚ₂/ℤ₂`, takes only the values 1/4, 1/2 and
  3/4 on those units. So `c(z) = ⌊4φ(z)⌋` is a proper 4-colouring of `Γ(K)`: a
  unit step adds `k/4` to `φ`, with `k ∈ {1, 2, 3}`, and so adds `k` to `c`
  modulo 4.
- **Lower bound.** The Moser spindle lies in the field. ∎

Checked in `hn/adelic.py`, `scripts/moser2adic.py` and
`tests/test_moser_field.py`.

Theorem 1 also follows from Fischer's theorem (Theorem 3 below). The field
lies in `ℚ(i, √3, √11)`, the plane over `ℚ(√3, √11)`, since
`a + b√−3 + c√−11 + d√33 = (a + d√33) + i(b√3 + c√11)`. The proof above is a
local one.

**Theorem 2 (the field of `five_rho7`).**
`χ(Γ(ℚ(√−3, √−11, √−247))) = 5`.

*Proof.*
- **The places over 11 do not split in `K`.** Above 11, `L_v = ℚ₁₁(√33)`
  (`√741 ∈ ℚ₁₁`), and `K_w = L_v(√−3)` is unramified over `L_v` because −3 is
  a non-residue mod 11. So no place of `L` over 11 splits in `K`.
- **Upper bound.** Proposition A.1 with `χ(G₁₁) = 5`: SAT, and 4 colours are
  UNSAT.
- **Lower bound.** `five_rho7` lies in the field; so does its 803-vertex
  subgraph `five_247_c`, whose lack of a 4-colouring is certified by a DRAT
  proof (`certificates/five_247_c_no4coloring.json`). ∎

The pulled-back colouring is checked exactly, at both places above 11, on every
graph the project grew in this field: up to 32 312 points and 282 909 edges, all
with 0 monochromatic edges (`scripts/reduce11.py`, `tests/test_reduce11.py`).
So no search in this field could reach six.

Both bounds of Theorem 2 also follow from earlier work. `K` lies in the plane
over `F = ℚ(√3, √11, √247)`. Since `3 ≡ 5²` and `247 ≡ 4² (mod 11)`, a place of
`F` over 11 has residue field `𝔽₁₁`, so Madore's Cor. 3.4 and Lemma 4.5
([arXiv 1509.07023](https://arxiv.org/abs/1509.07023)) give `χ(F²) ≤ 5`, by the
argument of his Prop. 4.6 for `ℚ(√3, √11)`. Exoo and Ismailescu's 5-chromatic
graph ([arXiv 1805.00157](https://arxiv.org/abs/1805.00157)) lies in `F²`, and a
quarter turn moves it into `K`: its building blocks have coordinates
`((a√3 + b√11)/36, (c + d√33)/36)`, and its rotations, such as
`(119 + 3√−247)/128`, lie in `K`. What the proof above adds is the local
form of the argument.

**The denominator principle.** If a unit vector of `K` is not integral at
some place of `K` above a prime `p`, then some place of `L` above `p` splits in
`K`. Integrality is meant in `K`, not coordinatewise: `ω = (−1/2, √3/2)` has 2
in the denominators of its coordinates, yet it is an algebraic integer. Each
rung of the known constructions therefore had to add a rotation that is not
integral at a place over a new prime.

## 4. Finite planes

For `q ≡ 5 (mod 6)` (next to `√−3`, every non-split place over a prime
`p > 3` has such a residue field), `G_q` has this structure:
- **Triangles.** Every edge lies in exactly two triangles, there is no `K₄`,
  and `G_q` is the union of `(q+1)/6` rotated triangular tori.
- **Eigenvalues.** They are `λ(n) = Σ_c (1 − η(c² − 4n)) e(c/q)`, with `η` the
  quadratic character.
- **Ramanujan.** Every non-trivial eigenvalue satisfies `|λ| ≤ 2√q`. Parametrise
  `N₁` by `ℙ¹(𝔽_q)`; then `λ` is an exponential sum `Σ ψ(f(t))` of a rational
  function with two simple, conjugate poles, and Weil's bound (A. Weil, *On some
  exponential sums*, Proc. Natl. Acad. Sci. USA 34 (1948) 204–207) gives `2√q`. The
  bound was also checked numerically for every `q < 400`. For the finite
  Euclidean planes these spectral facts go back to Medrano, Myers, Stark and
  Terras (*Finite analogues of Euclidean space*, J. Comput. Appl. Math. 68
  (1996) 221–238), who expressed the eigenvalues through Gauss and Kloosterman
  sums and bounded them with Weil's estimate.

**Proposition B.** For every prime `q ≥ 53`, `χ(G_q) ≥ 6`. The same bound holds
for the deeper quotients `Cay(O_w/π^r, T mod π^r)`, given the character-sum
estimate sketched in the proof.

*Proof.* Hoffman's bound (A. J. Hoffman, *On eigenvalues and colorings of graphs*,
in *Graph Theory and its Applications* (B. Harris, ed.), Academic Press, 1970,
79–91), in the ratio form of §14, gives `χ_f ≥ 1 + (q+1)/|λ_min|`. This exceeds 5 for `q > 62`
by the Ramanujan bound, and was computed for `q = 53, 59, 61` (§14; the first
version of this proof omitted 61). For the deeper levels we use that a
primitive level-`j` character sums to at most `2q^{j−1}` in modulus: a
stationary-phase argument, sketched in `notes/rigidity.md` §10 and checked
numerically for small `q` and `j`. ∎

The Delsarte LP (P. Delsarte, *An algebraic approach to the association schemes of
coding theory*, Philips Res. Rep. Suppl. 10 (1973)) gives the same bound as Hoffman, and adding the triangle
inequalities changes almost nothing. Below 53 the three-point bound of §14
goes further: it proves `χ(G_q) ≥ 6` for `q = 29, 37, 41, 43, 47`.

| `q` | `χ(G_q)` | evidence |
|---|---|---|
| 2 | 4 | the 2-adic analysis above |
| 3 | 3 | Hoffman's bound, `α ≤ 3` (`scripts/finite_hoffman.py --inert 3`), and a 3-colouring (SAT; Moorhouse) |
| 4 | 4 | `Cay(𝔽₁₆, μ₅)` is the Clebsch graph, whose independence number is 5 `< 16/3`; a 4-colouring (SAT, `tests/test_biquadratic_bounds.py`) |
| 5 | 4 | Hoffman's bound, `α ≤ 7 < 25/3` (`scripts/finite_hoffman.py --inert 5`), and a 4-colouring (SAT) |
| 7 | 4 | SAT (Moorhouse) |
| 11 | 5 | SAT |
| 13 | 5–6 | a 6-colouring (SAT); the three-point bound gives `α ≤ 42.64` (§14, certified), so `α ≤ 42 < 169/4` and there is no 4-colouring. Tabu search finds no 5-colouring, and CaDiCaL ran 25 minutes without an answer; independent sets of 36 points exist, above `169/5`, so bounding `α` cannot decide it |
| 17 | 5–6 | no 4-colouring and a 6-colouring (SAT); tabu finds no 5-colouring. Independent sets of 57 points are found easily, never 58. Since `5·57 < 289`, `α = 57` would prove `χ ≥ 6`. Two kissat runs found 57-point sets in the same orbit: *rosettes*, made of a point, its whole circle `N = 12` and 38 points on six other circles around it. No independent set of 58 points contains a point together with its whole circle (DRAT, `scripts/g17_alpha.py`, part A); the other case is open. SAT runs on the 5-colouring have not finished; the three-point bound of §14 gives only `α ≤ 63.33` (numerically, not certified) |
| 19 | 5 | SAT; triangle-free, with the linear 5-colouring `(a, b) ↦ c(a + b mod 19)` (§13) |
| 23 | 5–8, likely ≥ 7 | no 4-colouring (Hoffman's bound, `α ≤ 131 < 529/4`, `scripts/finite_hoffman.py --inert 23`); an interval 8-colouring (§12); tabu finds no 7-colouring, and independent sets of 87 points, against `529/6 ≈ 88.2`; the three-point bound of §14 gives `α ≤ 107.04` (numerically, not certified), and proving `χ ≥ 6` this way needs `α ≤ 105` |
| 29 | ≥ 6 | three-point bound: `α ≤ 163.25 < 841/5` (§14) |
| 31 | 5–8 | `G₃₁ = 𝔽₃₁²` (§12); no 4-colouring (Hoffman's bound, `α ≤ 240 < 961/4`, `scripts/finite_hoffman.py --inert 31`); the three-point bound gives only `α ≤ 200.89` (numerically, not certified), against `961/5 = 192.2` |
| 37 | ≥ 6 | three-point bound: `α ≤ 263.64 < 1369/5` (§14) |
| 41 | ≥ 6 | three-point bound: `α ≤ 300.73 < 1681/5` (§14); local search finds 213 |
| 43 | 6–8 | `G₄₃ = 𝔽₄₃²`; three-point bound: `α ≤ 347.79 < 1849/5` (§14); a linear 8-colouring (§12) |
| 47 | 6–10 | `G₄₇ = 𝔽₄₇²`; three-point bound: `α ≤ 371.42 < 2209/5` (§14); an interval 10-colouring (§12) |
| ≥ 53 | ≥ 6 | Proposition B |

## 5. Which fields can hold a 6-chromatic graph

A `k`-chromatic graph can live only in fields with no non-split place of
local chromatic number below `k`; the local chromatic number is at most 3 at a
ramified place, and at most `χ(G_q)` at an unramified one. Next to `√−3`, a
non-split place has residue field `𝔽_q` with `q ≡ 2 (mod 3)` or `q` a power
of 3 (`notes/rigidity.md` §9). So the rungs demand:
- **For five:** no ramified non-split place, and no non-split place with
  residue field `𝔽₂`, `𝔽₃`, `𝔽₅` or `𝔽₈` (local chromatic number at most 4).
  The Moser field fails at 2, and `ℚ(√−3, √−7, √−15)` fails at 5.
- **For six:** also none with residue field `𝔽₁₁`, nor any `𝔽_q`, `q < 53`,
  whose finite plane is 5-colourable at some level. By §14 no prime `q` from
  29 to 52 other than 31 has a 5-colourable plane at level 1; the deeper
  levels below 53 are open.

`scripts/fieldscreen.py` lists the non-split places of any multiquadratic CM
field:

| field | first non-split places | status for six |
|---|---|---|
| `ℚ(√−3, √−11, √−247)` | 11, 29 | excluded (Theorem 2) |
| `ℚ(√−3, √−11, √−23)` | 11, 17 | excluded (`χ ≤ 5`) |
| `ℚ(√−3, √−7, √−11)` | 17, 41, 83, 101 | no obstruction at level 1 if `χ(G₁₇) ≥ 6` (`χ(G₄₁) ≥ 6`, §14); deeper levels at 17 and 41 open |
| de Grey's `ℚ(√−3, √−7, √−11, √−15)` | 41, 101, 131 | no obstruction at level 1 (`χ(G₄₁) ≥ 6`, §14); deeper levels at 41 open |
| `ℚ(√−3, √−7, √2717)` | 59, 83, 89 | open |
| `ℚ(√−3, √−7, √−11, √−247)` | 83, 173 | open |

**A 5-chromatic graph in `ℚ(√−3, √−7, √−11)`.**
1. Take the carrier's forced pair at `d² = 64/9`.
2. Compose it to `d² = 16` with the rotation `cos = 1/8`, `sin = 3√7/8`.
3. Spindle the result with `(31 + 3√−7)/32`.

It needs no `√5` and no `√247`. The graph is
`data/five_tuned_16_1_3_7_11.json`: 4 081 points and 27 242 edges, not
4-colourable: CaDiCaL and kissat agree, and drat-trim checks kissat's DRAT
proof (`certificates/data_no4_checks.txt`).

The field has rotations with every small prime in a denominator: 2, 3, 5, 7 and
11. The Exoo–Ismailescu `λ`-closure of the graph (378 unit vectors, rank 8)
admits none of the cyclic periodic colourings tried. A growth run there reached
27 614 points, with kissat taking 8–30 minutes per hard step; the growth
scripts were later found to miss edges outside their own unit set (research
log, "Blind edges, the level-2 plane at 17, and fields with no local
obstruction").

## 6. A question

Is `χ(Γ(K))` the minimum, over the non-split places, of the local chromatic
numbers, whenever that minimum is at most 6?

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

If it held in general, `ℚ(√−3, √−7, √−11)` would be 6-chromatic as soon as its
local planes at 17 and 41 need six colours at every level (so far `χ(G₄₁) ≥ 6`
is proved at level 1, §14, and `χ(G₁₇)` is 5 or 6) and so do those at its other
non-split places, above 83, 101, …: at level 1 by Weil's bound, at the deeper
levels only given the character-sum estimate sketched in `notes/rigidity.md`
(Proposition 12). Then `χ(ℝ²) ≥ 6`. We
do not claim this.

Reduction to a finite field goes back to Woodall (1973), who proved
`χ(ℚ²) = 2` by reducing modulo 2, and to Fischer (1990), who coloured
`ℚ(√N)²` through homomorphisms to `ℤ/3` and `ℤ/4`. G. E. Moorhouse, *On the
chromatic numbers of planes* (draft, 2010), and Madore (2015) developed it for
number fields. For the planes `𝔽_q²` with the form `x² + y²` (which is `G_q`
when `q ≡ 3 (mod 4)`), see Le Anh Vinh, *On chromatic number of
unit-quadrance graphs (finite Euclidean graphs)*, arXiv math/0510092:
`√q/2 ≲ χ ≲ q/2`.

## 7. Fields with no local obstruction

Extending `F8 = ℚ(√−3, √−7, √−11)` by one more `√−d` can remove the non-split
places 17 and 41. `scripts/fieldscreen.py` finds no non-split place of norm
below 53, and none ramified, for `F8(√−d)` with squarefree `d < 400` equal to

> 1, 2, 21, 33, 42, 43, 59, 66, 77, 83, 86, 87, 103, 115, 118, 127, 154, 155,
> 166, 174, 185, 195, 203, 206, 213, 223, 230, 237, 247, 251, 254, 271, …

Values of `d` whose ratio is 21, 33 or 77 give the same field: 1, 21, 33 and 77
give `F8(i)`, and 2, 42, 66 and 154 give `F8(√−2)`. In these fields
Proposition B leaves no local 5-colouring at any place and any level, given
the estimate sketched in its proof. The case `d = 247` is `L16 = ℚ(√−3, √−7, √−11, √−247)`:
- **Non-split places.** Below 300 there are only 83 and 173.
- **Both 5-chromatic families.** `five_tuned_16` and `five_rho7` both contain
  the same 402-point carrier in the Moser field, and share 404 points in all.
  Their union, 6 080 = 4 081 + 2 403 − 404 points and 37 474 edges, is
  `data/L16_seed.json`.
- **No circular colouring, numerically.** The seed's 918 unit vectors
  (459 directions) span rank 12. The best `min_u ‖φ(u)‖` found over characters `φ` is 0.024, while a
  circular 5-colouring needs 0.2. The same search finds 0.204 at once in the
  Moser module. The exact MILP timed out, so this is evidence, not proof.

**Two more cautions.**
- **`F8(i)` is not automatically better.** `U ∪ iU` gives the Cartesian
  product `Γ(F8) □ Γ(F8)`, whose chromatic number is that of `Γ(F8)`, and no
  triangle has edges in both `U` and `iU`, since `F8` does not contain
  `e^{iπ/6}`. `L16` is a better extension:
  `F8 ∩ ℚ(√−3, √−11, √−247)` is the whole Moser field, so the two halves
  share triangles and spindles.
- **Level 2 does not simply equal level 1.** A section `σ(a) = τ(a) + q·s(a)`
  of the level-2 plane over `G_q` would have to satisfy a linear system over
  `𝔽_q`. That system is inconsistent for `q = 5, 11, 17`, so `G_q` does not
  embed in level 2 that way. At `q = 17`, tabu search on level 2 (83 521
  vertices, degree 306) never improves on the lift of the best level-1
  colouring.

## 8. The plane over `ℚ(√3, √11)`: `χ = 4`

Theorem 3 below is a theorem of K. G. Fischer (1994), which later work treated
as open. This section gives a short proof. It is not a step towards
`χ(ℝ²) ≥ 6`.

**History.**
- Fischer, *A planar geometric graph of chromatic number four*, Congr. Numer.
  104 (1994) 73–79 ([Zbl 0836.05030](https://zbmath.org/?q=an:0836.05030)), proved that
  `ℚ(√p, √q)²` has an additive 4-colouring, with values in `ℤ/4`, for
  squarefree coprime `p ≡ 3`, `q ≡ 11 (mod 16)` with `pq ≡ 1 (mod 32)`, and
  concluded
  `χ(ℚ(√3, √11)²) = 4`.
  We found this only after writing this section; we have read the zbMATH
  summary, not the paper. None of the works below cites it.
- Moorhouse ([draft, 2010](https://www.ericmoorhouse.org/pub/chromatic.pdf))
  noted that `ℚ(√3, √11)` is the smallest field whose plane contains a Moser
  spindle, so `χ ≥ 4`. He wrote: "We have not determined the exact value …"
- Madore ([arXiv 1509.07023](https://arxiv.org/abs/1509.07023), Prop. 4.6)
  proved `4 ≤ χ(ℚ(√3, √11)²) ≤ 5`, reducing at a place over 11 and using
  `χ(𝔽₁₁²) ≤ 5`.
- Exoo and Ismailescu ([arXiv 1805.00157](https://arxiv.org/abs/1805.00157),
  DCG 2020) asked whether a 5-chromatic unit-distance graph embeds in
  `ℚ[√3, √11]²`. So did Ismailescu, recalling a question of Mixon, in
  [Polymath16, thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/).
- Voronov, in [Polymath16, thread 17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/)
  (18 July 2021): "it seems likely that `χ(Q(i, √3, √11)) = 4` ... But as far
  as I know, nobody has proved this yet."

**Theorem 3 (Fischer, 1994).** `χ(ℚ(√3, √11)²) = 4`. So no 5-chromatic
unit-distance graph has all its coordinates in `ℚ(√3, √11)`.

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
  Exoo–Ismailescu's 214-point graph `H` (1 004 edges;
  [arXiv 1909.13177](https://arxiv.org/abs/1909.13177), Claim 2.2).
- On 638 further points and 3 012 edges (`data/ei_rho7.json`) there are no
  monochromatic edges, at either place.

**Formal proof.** `lean/Q311.lean` proves `χ = 4` in Lean 4 (L. de Moura and S. Ullrich, *The Lean 4 theorem prover and
programming language*, CADE-28, LNCS 12699 (2021) 625–635) with Mathlib (The mathlib Community, *The Lean mathematical
library*, CPP 2020, 367–381), with only Lean's standard axioms
([`lean/README.md`](../lean/README.md)). There the place over 2 is a valuation subring given by Chevalley's
extension theorem, and its residue field is shown to be `𝔽₂` with `π = √3 − 1` and a Hensel-type argument for
`(1 + √33)/2`, whichever of the two places it is.

**Coordinates versus the Hermitian form.** Moorhouse (Lemma 8.2, Lemma 8.4) and Madore
(Prop. 3.2, 3.8) reduce the *coordinates* `(x, y)`. That needs `x² + y²` to be
anisotropic modulo `𝔪` or `𝔪²`. At a place over 2 with `√3` this fails, since
`1² + 1² ≡ 0 (mod 𝔪²)`.
- Unit vectors such as `(−1/2, √3/2) = ω` have non-integral coordinates.
- They are integral in `O_w`, which is strictly larger than `O_v[i]`.

Reducing `z = x + iy` in `O_w` — the Hermitian form of the argument — sees the
`ω` that the coordinates hide. In the coordinates `α = x + y/√3`,
`β = 2y/√3` of `z = α + βω`, the squared distance is `α² − αβ + β²`, which is
anisotropic modulo `𝔪`. Theorem 3 is the case of this form of Madore's ¶6.6,
which states his Prop. 3.2 for any quadratic form. At odd places the two versions
agree, since `O_w = O_v[i]` there.

**Prior work.** The ingredients are not new.
- **The reduction.** In
  [Polymath16, thread 2](https://dustingmixon.wordpress.com/2018/04/22/polymath16-second-thread-what-does-it-take-to-be-5-chromatic/#comment-4013)
  (25 April 2018), David Speyer 4-coloured the *Moser ring*
  `R = O_K[1/3]`, `K = ℚ(√−3, √−11)`. He observed `R/2R ≅ 𝔽₄ × 𝔽₄`, and that
  either projection is a proper 4-colouring. In
  [thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/)
  (3 May 2018) he gave 8-periodic variants. Philip Gibbs then reported
  (9 May 2018) that an analysis and computer search by Tamás Hubai found
  that all 4-colourings of the ring have period 8.
- **The cosets.** Extending a colouring from a subgroup or a ring to the
  whole field by cosets is Fischer's (1990, Thm 1, through the component of
  the origin), Moorhouse's Lemma 4.2 and Madore's Prop. 3.2.
- **The integrality.** Unit vectors are integral at a place that does not
  split because the norm-one torus is compact there, a standard fact.
- Dúcz ([arXiv 2606.12325](https://arxiv.org/abs/2606.12325), 2026) gave
  geometric 4-colourings of the Moser lattice and ring.

What the proof uses beyond these is one check: the places of `L` over 2 are
inert in `L(i) = ℚ(i, √3, √11)`, a field of degree 8 that contains the Moser
field. Then every unit vector of the plane is a 2-adic unit, not only those of
the ring, and Speyer's colouring extends to the whole plane. For an expert in
local fields this is a short observation. Fischer had settled the question in
1994; the observation gives a short proof, and it also covers `ℚ(√2, √3)`
(§10), which Fischer's hypotheses exclude.

**Parts' question about limits.**
- Exoo–Ismailescu's `G₄₀` forces a pair at distance 8/3 alike in every
  4-colouring *with no monochromatic pair at distance `√(11/3)`*.
- In [Polymath16, thread 13](https://dustingmixon.wordpress.com/2019/07/08/polymath16-thirteenth-thread-bumping-the-deadline/)
  (25 July 2019), Parts chained such pairs to get alike pairs at every
  distance `8/9ⁿ`, which sum to 1, and called this a "funny proof" that the
  field needs five colours.
- Pálvölgyi replied (26 July) that colour need not pass to the limit, and
  Parts then asked (27 July) why it is lost.
- The colouring above satisfies the hypothesis and makes the point concrete.
  Every pair at `√(11/3)` is apart, because `√33/3` is a 2-adic unit. Every
  pair at distance `8/9ⁿ` is alike, because `(8/9ⁿ)u ∈ 8·O_w`. The colouring is
  continuous for the 2-adic topology, not the real one.

**Corollaries for real fields `F`.** Proposition A at a place `v | 2` gives:
- **Ramified.** If `v` ramifies in `F(i)`, then `χ(F²) = 2`. The norm-one
  residues satisfy `ρ(u)² = 1`, so `ρ(u) = 1`, and any `𝔽₂`-linear form with
  `λ(1) = 1` 2-colours the plane. This recovers `χ(ℚ²) = 2` (Woodall).
  It also recovers Moorhouse's Theorem 7.1, Lemma 8.4 and Theorem 8.5, and
  Madore's Prop. 3.9.
- **Inert, residue field `𝔽₂`.** Then `χ(F²) ≤ 4`. For quadratic fields this
  happens exactly when `d ≡ 3 (mod 8)`, a bound already proved by Fischer
  (Discrete Math. 82 (1990); see Payne, arXiv 0707.1177). Together with
  Moorhouse's Theorem 8.1, which excluded `d ≡ 47, 59, 83 (mod 84)`, it covers
  `d ≡ 59, 83, 131 (mod 168)`.
  - For example, `3 ≤ χ(ℚ(√59)²) ≤ 4`. Reduction at 11 gave only 5
    (`59 ≡ 2² (mod 11)`).
  - The classes left open are `d ≡ 47, 143, 167 (mod 168)`.
- **Inert, larger residue fields.** The finite plane at an inert place over 2
  with residue field `𝔽_{2^f}` is `Cay(𝔽_{4^f}, μ_{2^f+1})`. Its chromatic
  number is 4 for `f = 1, 2, 3`: it is `K₄` for `f = 1` and the Clebsch graph
  for `f = 2`; for `f = 3` it contains `K₄` (the subfield `𝔽₄`, since
  `μ₃ ⊂ μ₉`), and SAT finds a 4-colouring.
- **A necessary condition for five.** Suppose a real field `F` has
  `χ(F²) ≥ 5`. Then:
  - no place of `F` over 2 ramifies in `F(i)`;
  - no place over 2 with residue field `𝔽₂`, `𝔽₄` or `𝔽₈` is inert in `F(i)`;
  - `F` has no place with residue field `𝔽₃` or `𝔽₇`, since `χ(𝔽₃²) = 3`
    and `χ(𝔽₇²) = 4`.

  Heule's `ℚ(√3, √5, √11)` and Exoo–Ismailescu's `ℚ(√3, √11, √247)` pass at 2:
  `−1` is a square in every completion over 2. There `√5` or `√247` supplies
  `i`, since `5 ≡ −3` and `247 ≡ −1 (mod 8)`.

Theorem 3 is Fischer's (1994); the proof above is a short alternative. Before
finding Fischer's paper we had searched Polymath16 threads 1–18, the Polymath16
wiki and the web without finding a proof. That the result was overlooked is
shown by the later work listed under **History**, which treats the value as
open. The wiki page
[Algebraic formulation of Hadwiger–Nelson problem](https://web.archive.org/web/20210412075722/https://asone.ai/polymath/index.php?title=Algebraic_formulation_of_Hadwiger-Nelson_problem)
colours rings such as the Moser ring, not whole planes. The closest work is
Speyer's colouring of the Moser ring, above. A second review, carried out
independently with AI assistance on 25 September, found no mathematical error.
No mathematician has checked this yet, and nothing here has been refereed.

## 9. Split places colour modules

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
- **The split finite plane.** At a split place over an odd prime, with
  residue field `𝔽_q`, the target is `H_q = Cay(𝔽_q², {(a, 1/a)})`. Over 2
  the circle degenerates to the line `a + b = 1`; Proposition C (§12), which
  reduces `z` and `z̄` instead of the coordinates, gives `H_q` there too.

So a module is `k`-colourable as soon as one integral place has a
`k`-colourable finite plane.

**Chromatic numbers of `H_q`.**

| `q` | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 29 | 37 | 49 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `χ(H_q)` | 2 | 3 | 3 | 4 | 4 | **5** | 5–6 | ≥ 5 | 5–6 | 6–8 | ≥ 5 |
| tabu at 5 colours (best conflicts) | | | | | | | 2 | 4 | 222 | 1 290 | 4 581 |

- **Where the entries come from.** The exact values are from SAT. For
  `q ≡ 1 (mod 4)`, `H_q ≅ 𝔽_q²` (§12), so the entries for 17, 29 and 37 are
  those of `𝔽_q²` in §12 and §14. SAT finds no 4-colouring of `H₁₇` or `H₁₉`.
  `H₄₉` needs at least five colours because the 5-chromatic `five_tuned_16`,
  which has no 7 in any denominator, maps into it at a place over 7.
- **Hoffman.** The eigenvalues are Kloosterman sums, of modulus at most
  `2√q` (Weil, 1948, cited in §4). So `χ(H_q) ≥ 1 + (q − 1)/(2√q)` (Hoffman, 1970, cited in §4), which is above 5 once
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
  - the only small split place of degree 1 is at 37; its image is all of
    `H₃₇ ≅ 𝔽₃₇²`, which needs six colours (§14).
- **`F8` growth** (822 units): 3 and 5 are excluded. At 7 no unit has 7 in a
  denominator; the place is ramified with residue field `𝔽₄₉`, and tabu finds
  no 5-colouring of `H₄₉`.

So neither search is ruled out by a reduction mod a prime: at 37 this is
proved, and at 49 it holds as far as tabu can tell. What the gate does rule out is a unit set without Moser's `σ`: the
places of `L16` over 3, with residue field `𝔽₉`, would become integral, and
`H₉` 3-colours. Without `ρ₇` the places over 7 would become integral too, but
their residue field is `𝔽₄₉`, and `H₄₉` needs at least five colours: the
5-chromatic `five_tuned_16`, which has no 7 in any denominator, maps into it.

## 10. The plane over `ℚ(√2, √3)`: `χ = 4`

Voronov, in [Polymath16, thread 17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/#comment-29291)
(18 July 2021), wrote that `χ = 4` "seems likely" for two fields:
`ℚ(i, √3, √11)`, which is Theorem 3, and "case (2, 3)", that is
`ℚ(i, √2, √3) = ℚ(ζ₂₄)`. On 30 July
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

**Formal proof.** `lean/Q23.lean` proves `χ = 4` in Lean 4 with Mathlib (references in §8), with only Lean's standard axioms
([`lean/README.md`](../lean/README.md)). The place over 2 is a valuation subring given by Chevalley's extension
theorem; the uniformiser `π = (√2 + √6)/2 − 1`, a root of the 2-Eisenstein polynomial
`X⁴ + 4X³ + 2X² − 4X − 2`, shows that its residue field is `𝔽₂`.

**What else uses this.** The criterion is general:
- **Upper bound.** A real field `L` with a place over 2 that does not split
  in `L(i)`, and whose extension has residue field `𝔽₄`, has `χ(L²) ≤ 4`.
- **Lower bound.** `√3 ∈ L` gives triangles and rhombi. A three-rhombus
  chain closing in `L` then gives `χ(L²) ≥ 4`, so `χ(L²) = 4` when the upper
  bound holds too.

The lower bound was known, as above. We found no proof of the upper bound in
the literature: Fischer's 1994 hypotheses (squarefree, coprime `p ≡ 3`,
`q ≡ 11 (mod 16)`, `pq ≡ 1 (mod 32)`) exclude this field.

## 11. Two square roots: Voronov's question

In [Polymath16, thread 17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/#comment-29283)
(17 July 2021) Voronov asked: "can we get a chromatic number 5 if the vertices
of the graph lie in an extension formed by two roots of prime numbers? Is it
possible to do it without `√3` (and without triangle)?"

**With `√3`.** Let `q > 1` be squarefree and prime to 3, and `L = ℚ(√3, √q)`.

**Theorem 5.**
- If `q ≡ 1 (mod 3)`, then `χ(L²) = 3`.
- If `q ≡ 2 (mod 3)`, then `χ(L²) ≥ 4`, with equality when `q` is even or
  `q ≡ 1, 3 (mod 8)`. For prime `q` this means `q = 2` or `q ≡ 11, 17 (mod 24)`.
- If `q ≡ 5, 23 (mod 24)`, then `χ(L²) ≥ 4`, and the upper bound is open.
  The first cases are `q = 5, 23, 29, 47, 53, 71`.

*Proof.*
- **`q ≡ 1 (mod 3)`.** 3 splits in `ℚ(√q)`, so a prime of `L` above 3 has
  residue field `𝔽₃`. Madore's Cor. 3.4 gives `χ(L²) ≤ χ(𝔽₃²) = 3`, and
  triangles give `χ(L²) ≥ 3`.
- **Upper bound 4.** If `q` is even or `q ≡ 1, 3 (mod 8)`, a prime of `L`
  above 2 has residue field `𝔽₂`, and the criterion of §8 (Theorem 2 of the
  note `papers/planes-4-chromatic/planes-4-chromatic.pdf`) gives `χ(L²) ≤ 4`.
- **Lower bound 4.** If `q ≡ 2 (mod 3)`, then 3 splits in `ℚ(√−q)`, a subfield
  of `L(i)`. Let `h` be its class number and `α` a generator of the `h`-th
  power of a prime above 3. Then `N(α) = 3^h`, and 3 does not divide
  `Tr(α²)`: otherwise that prime would divide `ᾱ`. So
  `α/ᾱ + ᾱ/α = Tr(α²)/3^h` is a sum of two unit vectors with exact denominator
  `3^h`, and `1/3` is a sum of unit vectors. The lemma below gives
  `χ(L²) ≥ 4`. ∎

`tests/test_q3q.py` checks the generator `α` for every `q ≡ 2 (mod 3)` up to
113, and a 271-vertex chain of 90 rhombi for `q = 17`. PARI/GP recomputes the
decompositions of 2 and 3 for every prime `q < 75` (`scripts/decompositions.gp`,
checked by `tests/test_decompositions.py`).

**More square roots.** The proof of the first two parts works for every
multiquadratic field `L = ℚ(√3, √q₁, …, √q_k)`. Write `q'_j` for `q_j` with any
factor 3 removed.
- If every `q'_j ≡ 1 (mod 3)`, then 3 splits in every quadratic subfield prime
  to 3, a prime of `L` above 3 has residue field `𝔽₃`, and `χ(L²) = 3`. For
  example `χ(ℚ(√3, √7, √13, √19)²) = 3`.
- Otherwise some `q'_j ≡ 2 (mod 3)`, `ℚ(√−q'_j) ⊂ L(i)`, and `χ(L²) ≥ 4`.

PARI/GP confirms the residue degree above 3 for nine such fields, listed in
`scripts/decompositions.gp`.

So the answer to Voronov's first question is no for every `ℚ(√3, √q)` outside
the class `q ≡ 5, 23 (mod 24)`. Fischer's family meets these fields only at
`q ≡ 11 (mod 32)`, and for `q ≠ 11` it gives the upper bound alone.

**A lemma for the lower bound.** If `√3 ∈ L` and `1/3` is a sum of unit vectors
of `L(i)`, then `L²` has no proper 3-colouring. The sums of unit vectors form a
ring `C₀` that contains `ζ₆ = e^{iπ/3}`, so `(1 + ζ₆)/3 = Σ uₖ` with unit
vectors `uₖ`. The unit rhombus whose long diagonal is `√3 uₖ` forces its tips
alike in every 3-colouring. Chaining these rhombi joins `0` to
`√3(1 + ζ₆)/3 = e^{iπ/6}`, a unit vector.

The same chain works at any number of colours. Suppose a finite graph `G` in
`L²` has two vertices `A`, `B` that share a colour in every `k`-colouring, and
`w/(B − A)` is a sum of unit vectors `u_j` for some unit vector `w`. Chain the
copies of `G` turned by the `u_j`, each starting where the previous one's image
of `B` lies. They join a point to a point at unit distance `w`, so
`χ(L²) > k`. In `ℝ²` two copies always suffice: this is spindling. Inside a
fixed field the chain needs no spindle rotation, only this arithmetic
condition; the rhombus is the case `k = 3`, `B − A = √3`.

**The plane over `ℚ(√3, √5)`.** Here the generator is explicit.
`τ = (2 + i√5)/3` is a unit vector, so
`1/3 = τ + τ̄ − 1` and
`(1 + ζ₆)/3 = τ + τ̄ − 1 + ζ₆τ + ζ₆τ̄ − ζ₆`. The chain of the six rhombi is
`data/chain35.json`: 19 vertices, 31 edges, no proper 3-colouring
(`tests/test_q35.py`). Since 3 and 5 are squares modulo 11, 11 splits
completely, and reduction at a place over 11 gives `χ ≤ χ(𝔽₁₁²) = 5`. So
`4 ≤ χ(ℚ(√3, √5)²) ≤ 5`, and whether it is 5 is the smallest open case of
Voronov's first question.

Reducing modulo 121 instead of 11 does not help: `Cay((ℤ/121)², U₂)`, with the
132 unit vectors modulo 121, has no proper 4-colouring (kissat, CaDiCaL and
Glucose; `scripts/experiments/level_plane.py`).

**A caution about local evidence.** Minkowski balls and colouring-guided
growths over `ℚ(√3, √5)` of up to 83 000 points stayed 3-colourable. The chain
reaches distance `√7 ≈ 2.65` from its start, and those graphs stayed nearer the
origin, so they said nothing about `χ(L²)`.

## 12. Split places see the integral edges: five over `ℚ(√3, √5)` needs 2 and 3

§9 colours a module through a split place when all its generators are
integral there. The argument needs only the edge vectors.

**Proposition C.** Let `w` be a place of `K = L(i)` over a place `v` of `L`
that splits in `K`, with residue field `𝔽_q`, and let
`H_q = Cay(𝔽_q², {(t, 1/t) : t ∈ 𝔽_q^×})`. If every edge vector of a
unit-distance graph `G` over `L` is a unit at `w`, then `χ(G) ≤ χ(H_q)`.

*Proof.*
- **Units at `w` are units at `w̄`.** A unit vector has `u ū = 1`, so
  `v_w(u) = −v_{w̄}(u)`.
- **A residue map.** Let `R` be the ring of elements integral at `w` and `w̄`,
  and `ρ` the reduction at `w`. Fix a representative `rep(z)` of each coset
  `z + R`, and put `g(z) = (ρ(z − rep z), ρ(z̄ − conj(rep z)))`.
- **Unit steps become edges of `H_q`.** An edge `z, z + u` stays in one coset,
  and `g` moves by `(ρ(u), ρ(ū)) = (ρ(u), ρ(u)^{−1})`. ∎

The points need not be integral. §9 is the case where they are.

**More values of `χ(H_q)`** (SAT; `tests/test_split_places.py`):

| `q` | 4 | 8 | 9 | 16 |
|---|---|---|---|---|
| `χ(H_q)` | 4 | 4 | 3 | 4 |

`χ(H₁₃) = 5` is certified there by a linear colouring:
`(a, b) ↦ c(a + 2b mod 13)`, where `c` colours the circulant on
`{t + 2/t} = {2, 3, 5, 8, 10, 11}`. For `q ≡ 1 (mod 4)` the substitution
`(a, b) = (x + iy, x − iy)` turns the circle `x² + y² = 1` into the hyperbola
`ab = 1`, so `H_q` is the unit-distance graph of `𝔽_q²`. The value therefore
settles the entry "5 or 6" for `𝔽₁₃²` in Moorhouse's Table 6.1.

**`𝔽₁₇²` needs at most six colours.** Moorhouse's entry is "5, 6 or 7".
- **The lines.** The line `αx + βy = r` meets the circle `x² + y² = 1`
  exactly when `α² + β² − r²` is a square (the discriminant of the
  intersection is `β²(α² + β² − r²)`).
- **Three lines per colour.** For `(α, β) = (3, 6)`, the values
  `45 − r² ≡ 11, 10, 7` (`r = 0, 1, 2`) are not squares mod 17. So `3x + 6y`
  never takes the values `0, ±1, ±2` on a unit vector, and each block of three
  consecutive lines `3x + 6y ∈ {3j, 3j + 1, 3j + 2}` is independent.
- **The colouring.** `(x, y) ↦ ⌊((3x + 6y) mod 17)/3⌋` is a proper
  6-colouring (tested).

The entry becomes "5 or 6". This is Vinh's colouring by pairs of parallel lines
with longer blocks: if `N, N − 1, N − 4, …, N − (m − 1)²` are all non-squares,
then `χ(𝔽_q²) ≤ ⌈q/m⌉`. By the Weil bound (1948, cited in §4) such an `N` exists once
`2^m (m + 1) < √q`. So the colouring gives `χ(𝔽_q²) ≤ (2 + o(1)) q / log₂ q`,
against `q (1/2 + o(1))` for pairs of lines. General theorems give the same
order: Molloy's bound for triangle-free graphs (M. Molloy, *The list chromatic
number of graphs with small clique number*, J. Combin. Theory Ser. B 134 (2019)
264–284) when `q ≡ ±7 (mod 12)`, and
Alon–Krivelevich–Sudakov's for locally sparse graphs (N. Alon, M. Krivelevich,
B. Sudakov, *Coloring graphs with sparse neighborhoods*, J. Combin. Theory
Ser. B 77 (1999) 73–82) otherwise. Linear
colourings `(x, y) ↦ c(αx + βy)` cannot do better for `q = 17`: the nine
admissible circulants `{αx + βy : x² + y² = 1}` all have independence number
at most 3 (exhaustive search), and `17/3 > 5`, so each needs six colours.
Whether five colours suffice is open. Tabu search stops at two monochromatic
edges. For comparison, a random graph with 289 vertices and 2312 edges has on
average `5^289 (4/5)^2312 ≈ e^{−51}` proper 5-colourings.

Two observations on five colours, neither a proof:
- **The two bad edges share an isotropic line.** Tabu search and a hybrid
  evolutionary search always stop at exactly two monochromatic edges. In all
  28 such colourings sampled, the two edges have endpoints on a common
  isotropic line `x + 4y = c` or `x − 4y = c`. For two random edges this
  happens with probability about 0.4, so the pattern is structural. It never
  reached one monochromatic edge in 5·10⁷ moves.
- **Cube and conquer isolates the hard case.** Pin one edge and break the
  symmetry of the other three colours by value precedence (colour `c + 1`
  first appears after colour `c`; Y. C. Law and J. H. M. Lee, *Global
  constraints for integer and set value precedence*, CP 2004, LNCS 3258,
  362–376). `march_cu` (M. J. H. Heule, O. Kullmann, S. Wieringa, A. Biere,
  *Cube and conquer: guiding CDCL SAT solvers by lookaheads*, HVC 2011,
  LNCS 7261 (2012) 50–65) then splits the CNF into 8614
  cubes, and the cubes provably cover every case (`drat-trim` checks the
  proof). Almost all cubes are degenerate cases, where some colour first
  appears late. The first 7000 were each refuted within 13 seconds, 1706 of
  them also with DRAT proofs checked by `drat-trim`. The last cube says only
  that colour 2 appears among the first 23 vertices, and it is as hard as the
  whole problem.

`H_q` is the *hyperbola graph* `HG(𝔽_q)` of Bardestani and Mallahi-Karai
([arXiv 1507.05300](https://arxiv.org/abs/1507.05300)), who define `HG(F)` for
every field `F`. They show that `HG(ℝ)` and `HG(ℚ_p)` have infinite Borel
chromatic number, and that in characteristic other than 2 `HG(F)` sits inside
the graph of every isotropic quadratic form over `F`. Since `H₉ ≅ 𝔽₉²`,
`χ(H₉) = 3` is Moorhouse's value for `𝔽₉²`.

**Moorhouse's table, continued** (`tests/test_finite_planes.py`,
`tests/test_finite_planes_slow.py` for 29, 31, 41 and 43, and
`tests/test_threepoint_certificates.py` for the sixth colour):

| `q` | 19 | 23 | 29 | 31 | 37 | 41 | 43 | 47 | 53 | 59 | 61 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `χ(𝔽_q²)` | 5 | 5–8 | 5–6 | 5–8 | 6–8 | 6–7 | 6–8 | 6–10 | 5–11 | 6–10 | 5–11 |

- **Upper bounds.** Interval colourings with `m = 4, 3, 5, 4, 5, 6` lines per
  colour for `q = 19, …, 41`, and `m = 5, 5, 6, 6` for `q = 47, …, 61`. For
  `q = 43` the linear colouring `(x, y) ↦ c(x + y mod 43)` with a circulant
  8-colouring `c` does better than intervals (9). Otherwise no linear colouring
  `(x, y) ↦ c(αx + βy)` beats the intervals: for the other `q ≤ 47` the circulant
  `{αx + βy : x² + y² = 1}` has independence number `α_c` with
  `⌈q/α_c⌉` equal to the interval count.
- **Lower bounds.** There is no 4-colouring for `q = 23, 31` (Hoffman's bound,
  in interval arithmetic: `scripts/finite_hoffman.py --inert 23 31`), for
  `q = 29` (SAT, with an edge pinned, since there is no unit triangle), for
  `q = 37, 41, 43` (they need six colours, below), nor for `q ≥ 47` (Hoffman's
  bound). Six colours are needed for `q = 37, 41, 43, 47` (the three-point
  bound) and for `q = 59` (Hoffman's bound), §14.
- **Optimal cases.** The interval colourings are optimal for `q = 7, 13, 19`.
- **`q = 23`.** Tabu finds no 7-colouring (six monochromatic edges at best), and
  its largest independent sets have 87 points. Since `529/87 > 6`, an
  independence number of 87 would give `χ ≥ 7`.

**Over `ℚ(√3, √5)`.** Both 2 and 3 split in `K = ℚ(i, √3, √5)`.
- **At 2.** `−15 ≡ 1 (mod 8)`, so `i ∈ ℚ₂(√15)`. The local field is
  `ℚ₂(√3, √5)`, with residue field `𝔽₄`.
- **At 3.** `−5 ≡ 1 (mod 3)`, so `i ∈ ℚ₃(√5)`. The residue field is `𝔽₉`.

Since `χ(H₉) = 3` and `χ(H₄) = 4`:
- a unit-distance graph over `ℚ(√3, √5)` with no proper 3-colouring has an
  edge vector that is not a unit above 3. In `chain35.json` it is
  `τ = (2 + i√5)/3`.
- one with no proper 4-colouring also has an edge vector that is not a unit
  above 2, such as `(1 + i√15)/4` (valuations `±2`) or de Grey's
  `(7 + i√15)/8` (`±4`).

The bound from `χ(H₄) = 4` is attained: `chain35.json` has no 3-colouring,
all its edge vectors are units above 2, and `z ↦ z mod w₂ ∈ 𝔽₄` is a proper
4-colouring of it (tested).

**Why the growths over `ℚ(√3, √5)` stopped at four.** Every unit set used
before this section (`ζ^a τ^b`, and the four-rhombus units below) consisted of
units above 2. Proposition C 4-colours every graph they can build, so no
amount of growth could have produced five.

**The known 5-chromatic graphs agree.** The coordinate fields of de Grey's
graph `ℚ(√3, √5, √7, √11)`, of Voronov–Neopryatnaya–Dergachev's
`ℚ(√2, √3, √5)` and of Exoo–Ismailescu's `ℚ(√3, √11, √247)` all have local
field `ℚ₂(√3, √5)` or a field containing it with residue field `𝔽₄`, and there
2 splits in `K`. So for every place `w` above 2, each of these graphs has an
edge vector that is not a unit at `w`. Their spindle rotations are the natural
candidates: `(7 + i√15)/8`, and `(119 + 3i√247)/128` for Exoo–Ismailescu.

**A shorter chain.** `1/√3` is a sum of four unit vectors of `K`:

    c₁ = (−(√3 + √15) + i(√15 − √3))/6,   c₂ = ((√15 − √3) + i(√3 + √15))/6,
    c₃ = ((2√3 − √5) − i(2 + √15))/6,     c₄ = ((2√3 + √5) + i(2 − √15))/6.

So there is a chain of four unit rhombi, 13 vertices, with no proper
3-colouring. Three do not come out of a search over small heights: three unit
vectors summing to `1/√3` are an `L`-point on a curve of genus 1.

**Finer 2-adic levels.** With units that are not units above 2, reduce at
`w^A w̄^B` after scaling by a power of the uniformiser
(`scripts/experiments/split_gate_q35.py`; the unit sets come from
`scripts/experiments/units_q35.py`). A proper 4-colouring of the finite quotient
4-colours every graph with those edge vectors, and colourability passes to
finer levels.

| unit set (with `ζ`-turns and products) | levels tried at 2 | 4-colourable? |
|---|---|---|
| `τ` | `A = B = 1` | yes (the image is `K₄`) |
| `τ, (1 + i√15)/4` | `3, 3` / `4, 4` | no / **yes** |
| `o₁ = ((√5 − √3) − i(√3 + √5))/4` (valuation `±1`) | `2, 2` | yes |
| `o₁, (1 + i√15)/4` | `4, 4`; `5, 4`; `5, 5` (1 048 576 points) | no (kissat) |

With `τ` added to the last set (252 units), the reductions at 3 (level `3, 2`)
and at 5 (`1, 1`, `2, 1` and `2, 2`, 390 625 points) have no proper 4-colouring
either, and 11 never gives a 4-colouring. So none of the quotients tried at 2,
3, 5 and 11 4-colours that module; finer levels at 3 and 5, and the other
primes, were not tried. It is the first unit set over `ℚ(√3, √5)` that no
reduction we ran 4-colours. A colouring-guided search inside it is running.

## 13. Two square roots without `√3`: the local bounds

Voronov's second question asks for a 5-chromatic graph over two square roots
of primes without `√3`, so without triangles. Proposition A, applied at the
places of `L = ℚ(√a, √b)` that do not split in `K = L(i)`, answers it for most
pairs.
- **Ramified above 2.** The residues of the unit vectors are `±1 = 1`, the
  target is a perfect matching, and `χ(L²) ≤ 2`. Moorhouse's Lemma 8.4 is the
  quadratic case, `χ(ℚ(√d)²) = 2` for `d ≡ 1 (mod 4)`.
- **Unramified with residue field `𝔽_q`.** Then `χ(L²) ≤ χ(G_q)`, with
  `χ(G_q) = 4, 4, 3, 4, 5, 5` for `q = 2, 4, 3, 7, 11, 19` (the table in §4).

`G₁₉` is new here. It has no triangle, needs five colours, and has the linear
5-colouring `(a, b) ↦ c(a + b mod 19)`, since `{a + b : a² + b² = 1}` avoids 0.

**Pairs of primes below 60**, places below 200
(`scripts/experiments/classify_biquadratic.py`; spot checks in
`tests/test_biquadratic_bounds.py`):

| local bound | 2 | 3 | 4 | 5 | none from `q ≤ 19` |
|---|---|---|---|---|---|
| fields without `√3` (120) | 28 | 14 | 33 | 28 | 17 |
| fields with `√3` (16) | 0 | 6 | 5 | 4 | 1 |

- **Most fields are settled.** 75 of the 120 fields without `√3` have
  `χ ≤ 4`, so they cannot answer Voronov's second question. For 28 of them
  the plane is even bipartite, for example `ℚ(√2, √5)` and `ℚ(√5, √13)`.
- **Bound 5.** 28 fields remain, such as `ℚ(√5, √7)`, whose first non-split
  place lies above 19. There the local plane `G₁₉` is itself triangle-free and
  5-chromatic, so the local picture does not rule out a triangle-free
  5-chromatic graph over `ℚ(√5, √7)`. That plane does have odd cycles. With
  `τ = (2 + i√5)/3` and `σ = (1 + i√35)/6`, `1 + 2 Re σ = 2 Re τ`, so the unit
  steps `1, σ̄, −τ, σ, −τ̄` close up into a 5-cycle. Hence
  `3 ≤ χ(ℚ(√5, √7)²) ≤ 5`.
- **Where a fourth colour must come from.** The places above 3 split in `K`.
  Since `√7 ∈ ℚ₃`, the local field is `ℚ₃(√5)`, the unramified quadratic
  extension. It contains `i`, and its residue field is `𝔽₉`.
  Since `χ(H₉) = 3`, Proposition C 3-colours every graph over `ℚ(√5, √7)` whose
  edge vectors are units above 3. So a 4-chromatic graph needs an edge vector
  like `τ`, which is not a unit above 3, as over `ℚ(√3, √5)` (§12).
- **Is it 3?** No 3-colouring can be additive. Since `τ + τ̄ = 4/3`, the
  number `1/3 = τ + τ̄ − 1` is a sum of unit vectors, and a homomorphism to
  `ℤ/3` would send the unit vector `1 = 3 · (1/3)` to 0. Yet the experiments
  find no obstruction:
  - the ball `U + U + U` of radius 1.5 for the 130 unit vectors of height at
    most 60 (187 003 points) is 3-colourable;
  - a colouring-guided growth with three colours stops finding points that see
    all three colours after 19 steps (research log).

  A 3-colouring would answer the Question of §6 in the negative. Every
  non-split place lies above a prime `q ≥ 19`, where Hoffman's bound gives
  `χ(G_q) ≥ 1 + (q + 1)/(2√q) > 3`, and the place above 19 gives exactly 5.
- **No small place.** For the other 17, such as `ℚ(√2, √31)`, `ℚ(√2, √47)` and
  `ℚ(√7, √41)`, the first non-split places lie above 23, 31, 43 or 59. There
  `χ(G₄₃) ≥ 6` and `χ(G₅₉) ≥ 6` (§14, §4), while `χ(G₂₃)` and `χ(G₃₁)` lie
  between 5 and 8.

**With `√3` the table reproduces Theorem 5.** It gives bound 3 exactly for
`q ≡ 1 (mod 3)`, and bound 4 exactly for `q = 2` and `q ≡ 11, 17 (mod 24)`. It
adds one fact: the first non-split place of `ℚ(√3, √29)` lies above 23, where
`5 ≤ χ(G₂₃) ≤ 8` and tabu finds no 7-colouring (§4, §12). So the local arguments available give no
bound below `χ(G₂₃)`, although the field has triangles. Among the fields
`ℚ(√3, √q)` with `q` a prime below 60 it is the only one without a non-split
place at `q ≤ 19`, which makes it a candidate for §5's screen for six.

## 14. Six colours for most finite planes

§4 and §12 left the lower bound at 5 for most finite planes beyond `q = 19`.
Both kinds of plane are vertex-transitive, so `χ ≥ χ_f = q²/α`, and
`α < q²/5` gives `χ ≥ 6`. Local search finds independent sets of only
`0.13–0.17 q²` for `23 ≤ q ≤ 43`, which suggests six colours, but SAT does not
prove it: kissat ran 40 minutes on `𝔽₂₉²` with five colours without an answer.
An upper bound on `α` can, and does for seven planes. Here `𝔽_q²` is the plane
of Moorhouse's table (`x² + y²`; his table stops at `q = 17`), and `G_q` the
anisotropic plane of §4; they coincide for `q ≡ 3 (mod 4)`.

**The spectral bound.** The eigenvalues are `λ_ξ = Σ_u cos(2π ξ·u/q)` over the
`q ∓ 1` unit vectors, and Hoffman's ratio bound gives
`α ≤ n (−λ_min)/(d − λ_min)` (Hoffman, 1970, cited in §4; in this form, W. H. Haemers,
*Hoffman's ratio bound*, Linear Algebra Appl. 617 (2021) 215–219). This is Vinh's method
([arXiv math/0510092](https://arxiv.org/abs/math/0510092)), with the spectra
of Medrano–Myers–Stark–Terras (§4) and, for the anisotropic planes `G_q`, of
Bannai, Shimabukuro and Tanaka (*Finite Euclidean graphs and Ramanujan graphs*,
Discrete Math. 309 (2009) 6126–6134); the thresholds below follow from it.
- **Weil.** Every non-trivial eigenvalue has `|λ| ≤ 2√q` (Weil, 1948, cited in §4): a Kloosterman sum
  for `q ≡ 1 (mod 4)`, the sum of §4 for `q ≡ 3`. Hence
  `χ(𝔽_q²) ≥ 1 + (q ∓ 1)/(2√q)`, which exceeds 5 for every prime `q ≥ 67` and
  6 for every prime `q ≥ 103`.
- **Exact spectra**, in interval arithmetic (`scripts/finite_hoffman.py`):
  `χ(𝔽_q²) ≥ 6` for `q = 59`, and `χ ≥ 7` for `q = 71, 97, 101`. For the
  anisotropic planes the same computation gives `χ_f(G_q) > 5` for
  `q = 53, 59, 61`, the cases of Proposition B below Weil's threshold
  (`χ_f(G₆₁) ≥ 5.20`; the first version of its proof omitted 61).
- **Where it stops.** For `q ≤ 53` and `q = 61` the ratio bound for `𝔽_q²`
  stays above `q²/5` (`0.2111 q²` at `q = 47`, `0.2036 q²` at `q = 61`). The
  Delsarte linear programme (Delsarte, 1973, cited in §4) on the circles `x² + y² = r` gives exactly the same
  numbers, and conditional constraints from unit triangles, Moser spindles and
  rigid templates move it by under 2%. Linear-programming bounds with triangle
  constraints were used earlier for planar sets avoiding unit distance, by
  Keleti, Matolcsi, de Oliveira Filho and Ruzsa (*Better bounds for planar sets
  avoiding unit distances*, Discrete Comput. Geom. 55 (2016) 642–661) and by
  DeCorte, de Oliveira Filho and Vallentin (*Complete positivity and
  distance-avoiding sets*, Math. Program. 191 (2022) 487–558).

**The three-point bound.** Translate an independent set `S` by a random vector
and turn it about the origin by a random isometry; let `v` be the indicator of
the moved set and `x₀` its value at the origin. Let `g(w)` be the probability
that `w` lies in the moved set given that `0` does, and `z(a, b)` the
probability that `a` and `b` both do, again given that `0` does. Then:
- `g` depends only on the circle of `w`, `g(0) = 1`, and `|S| = Σ_w g(w)`;
- `z(a, b)` depends only on the triangle `{0, a, b}` up to isometry and
  relabelling, and vanishes when a side is a unit vector;
- `M₁ = [z(a, b)]` over all `a, b`, with `z(0, 0) = 1`, `z(0, b) = g(b)` and
  `z(a, a) = g(a)`, and `M₀ = [g(b − a) − z(a, b)]` over `a, b ≠ 0`, with
  diagonal `1 − g(a)`, are positive semidefinite: up to the factor `n/|S|` they
  are the averages of `x₀ vvᵀ` and `(1 − x₀) vvᵀ`;
- `0 ≤ z ≤ g` on each side, `1 − g(a) − g(b) + z(a, b) ≥ 0`, and the Delsarte
  inequalities hold.

Maximising `Σ_w g(w)` under these constraints is Schrijver's three-point bound
for codes (IEEE Trans. Inform. Theory 51, 2005), moved to the finite plane.
Both matrices commute with the rotation group, which is cyclic of order `N`,
the number of unit vectors (`q ∓ 1` for `𝔽_q²`, `q + 1` for `G_q`). In the
eigenbasis of a rotation, made real by a reflection, each splits into
`N/2 + 1` real blocks of size about `q` (`scripts/threepoint.py`). The
variables are the triangle classes: 4 429 for `𝔽₃₇²`, 8 625 for `𝔽₄₇²`.

Where the plain bound falls short, more constraints can be added: conditional
triangle inequalities, and the localizing matrices of `1 − x₀ − x_e`,
`1 − x₀ − x_e − x_f` and `2 − Σ_t x_t` for a unit edge, a unit triangle and a
unit 5-cycle. A target `K₀` enters through the corners of these matrices,
which are valid when `|S| ≥ K₀`; a bound below `K₀` then proves `α < K₀`.

**Results.** Each bound below is proved by a dual solution stored in
`data/threepoint/`, which `scripts/threepoint_verify.py` checks.
`scripts/threepoint_verify_indep.py`, written independently of the solver and of that checker, checks
every one again (`certificates/threepoint_indep_checks.txt`):

| plane | `n = q²` | bound on `α` | `n/5` | so |
|---|---|---|---|---|
| `G₂₉` | 841 | 163.25 | 168.2 | `χ ≥ 6` |
| `𝔽₃₇²` | 1 369 | 259.90 | 273.8 | `χ ≥ 6` |
| `G₃₇` | 1 369 | 263.64 | 273.8 | `χ ≥ 6` |
| `𝔽₄₁²` | 1 681 | 327.68 | 336.2 | `χ ≥ 6` |
| `G₄₁` | 1 681 | 300.73 | 336.2 | `χ ≥ 6` |
| `𝔽₄₃² = G₄₃` | 1 849 | 347.79 | 369.8 | `χ ≥ 6` |
| `𝔽₄₇² = G₄₇` | 2 209 | 371.42 | 441.8 | `χ ≥ 6` |

The same programme gives `α(G₁₃) ≤ 42.64` (`inert13.npz`). This is above
`169/4 = 42.25`, but `α` is an integer, so `α(G₁₃) ≤ 42` and `χ(G₁₃) ≥ 169/42 > 4`.
Each file's `value` is the solver's floating-point bound; the verifier's
rigorous bound differs from it by less than `10⁻⁴` and lies below every number
quoted here.
With Proposition B and the spectral bound:
- **`χ(G_q) ≥ 6` for every prime `q ≥ 29` except `q = 31`**;
- **`χ(𝔽_q²) ≥ 6` for `q = 37, 41, 43, 47, 59` and every prime `q ≥ 67`.**
  With the colourings of §12, `χ(𝔽₄₁²) ∈ {6, 7}` and `χ(𝔽₃₇²), χ(𝔽₄₃²) ∈
  {6, 7, 8}`.

**Checks.**
- **Blocks.** For `q ≤ 13` the blocks equal the compressions of the explicit
  `n × n` matrices by an orthonormal basis (`tests/test_threepoint.py`), so
  they are positive semidefinite exactly when the matrices are. The same holds
  for the localizing blocks.
- **Feasibility.** Independent sets for `q = 11, 13` satisfy every constraint,
  with objective `|S|` (the same test), and so did an independent set of 57
  points in `G₁₇` (checked once, not part of the tests).
- **Sharpness.** For `q = 7` and 11 the bound is 14.09 and 29.18, against
  `α = 14` and 28.
- **Certificates.** DSDP solves the primal. The dual is then found on the null
  spaces of the blocks at the optimum, where the dual matrices must live. The
  verifier rebuilds the blocks in interval arithmetic, clips each dual matrix
  to its positive part, shifts it by `10⁻⁹` of its largest eigenvalue, proves
  it positive definite by an exact rational `LDLᵀ`, and bounds the rest using
  `0 ≤ g, z ≤ 1`. `tests/test_threepoint_certificates.py` runs it on every
  stored certificate (marked `slow`: up to four minutes each).

**What is left.** The values below are the solver's, not certified.
- **`G₂₃ = 𝔽₂₃²`.** The plain bound is 108.55. With every constraint above and
  the target 106 it is 107.04, and proving `χ ≥ 6` this way needs `α ≤ 105`.
  The triangle's localizing matrix gives almost all of the gain (107.07 on its
  own).
- **`G₃₁ = 𝔽₃₁²`.** 200.89, against `961/5 = 192.2`.
- **`𝔽₂₉²`.** 184.77, against 168.2; `G₂₉` is done.
- **`G₁₇`.** 64.56, and 63.33 with the conditional triangles and the
  localizing matrices, against 57.8. So `χ(G₁₇)` stays open, and with it the
  place above 17 of `ℚ(√−3, √−7, √−11)` (§5).
- **`𝔽₅₃²` and `𝔽₆₁²`.** The programmes have 12 827 and 19 501 triangle
  variables, against 8 625 at `q = 47`, where DSDP already used 6 GB. They have
  not been run.
- **Colourings.** Using that the five classes form a colouring, not just that
  each is independent, adds a positive semidefinite block (the standard
  representation of the symmetric group on the other four colours) and the
  inequality "three colours on a triangle" `1 − g(a) − g(b) − g(b − a) +
  2z(a, b) ≥ 0`. They move the bound by about 0.01 at most (`G₁₇`: 64.556
  to 64.545; `𝔽₂₃²`: 108.547 to 108.542) (research log;
  `scripts/experiments/threepoint_colouring.py`).

So over the prime fields, `χ(G_q) ≥ 6` is open only for `q = 13, 17, 23, 31`;
for `q = 2, 3, 5, 7, 11, 19` five colours or fewer suffice (§4).
