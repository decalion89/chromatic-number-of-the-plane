# The ruler-and-compass plane is not 5-colourable

**Status (9 October 2026).** Formally verified in Lean 4 and refereed once.

- **Lean.** The kernel accepts the proof, and `#print axioms` lists only `propext`, `Classical.choice` and
  `Quot.sound`. The proof was rebuilt from OpenAI's unmodified sources and the patch, and `leanchecker` replayed
  every module in the kernel.
- **Referee.** An independent mathematical reading found no error.
- **No human referee yet.**

The proof is OpenAI's proof that χ(ℝ²) ≥ 6
([*The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf),
23 September 2026, formalized in [openai/math](https://github.com/openai/math)), run on a smaller field. What is new
here is the change of field, the observation that square roots are enough, and the formal check of the result.

## Statement

Identify the plane with ℂ and consider four countable subfields of ℂ:

| field | definition |
|---|---|
| `E_c`, the constructible numbers | the smallest subfield of ℂ closed under square roots: the points constructible with ruler and compass from 0 and 1 |
| `E_o`, the origami numbers | the smallest subfield of ℂ closed under square roots and cube roots: the points constructible by paper folding (Alperin, 2000) |
| `E_r`, the numbers expressible by radicals | the smallest subfield of ℂ closed under n-th roots for every n (`solvableByRad ℚ ℂ` in Mathlib) |
| `ℚ̄`, the algebraic numbers | the field OpenAI's proof works with |

They are nested: `E_c ⊂ E_o ⊂ E_r ⊂ ℚ̄`. The unit-distance graph on a subset of ℂ joins two points at distance
exactly 1.

**Theorem.** Every colouring of the constructible points of the plane with five colours gives the same colour to
two points at distance 1. So the unit-distance graph on `E_c` has chromatic number 6 or 7.

A colouring of a larger plane restricts to `E_c`, so the theorem also holds for `E_o`, `E_r`, `ℚ̄` and ℝ². Isbell's
hexagonal 7-colouring of ℝ² restricts to each of these planes, which gives the upper bound 7.

The Lean statement, from `ConstructiblePlane.lean` in
[`lean/external/openai-five/fields/`](../lean/external/openai-five/fields/README.md), writes the field out:

```lean
theorem sqrtClosure_plane_not_five_colourable :
    ¬ ∃ c : (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ) → Fin 5,
      ∀ x y : (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ),
        ‖(x : ℂ) - (y : ℂ)‖ = 1 → c x ≠ c y
```

The set in `sInf` contains ℂ itself, so the infimum is the intersection of all subfields closed under square
roots: the constructible numbers. The norm is the usual modulus on ℂ. We compiled the origami and radical cases
first. They need fewer changes, and they follow from the constructible case.

**Corollary.** Some finite unit-distance graph with constructible coordinates has no proper 5-colouring.

*Proof.* By the de Bruijn–Erdős theorem, a graph is 5-colourable if all its finite subgraphs are. ∎

This step is formally verified too: `ConstructibleFinite.lean` proves the Corollary from the Theorem with Mathlib's
compactness theorem for graph homomorphisms, `SimpleGraph.nonempty_hom_of_forall_finite_subgraph_hom`. The axioms are
again `propext`, `Classical.choice` and `Quot.sound`.

The coordinates of such a graph lie in a field obtained from ℚ by finitely many square roots. The proof gives no
such graph, and we know of none.

For comparison, two fields obtained from ℚ by square roots have planes of chromatic number 4: ℚ(√2, √3) (the main
result of this repository) and ℚ(√3, √11) (Fischer, 1994). Adjoining to ℚ(√2, √3) finitely many further square
roots, enough to contain the coordinates of a graph as in the Corollary, gives a field whose plane has chromatic
number at least 6. We do not know how many square roots are needed.

## What OpenAI's proof uses of the field

OpenAI's paper works with `E = ℚ̄`, its real part `F = E ∩ ℝ`, and the group `K` of elements of `E` of modulus 1
(Definition 2.1). The proof uses three properties of `E`:

1. `E` is a countable subfield of ℂ, algebraic over ℚ, closed under complex conjugation, containing `i`.
2. `F` is closed under square roots of its nonnegative elements. This is used in four places:
   - the decomposition `x = k + k′` with `k = x/2 + i√(1 − x²/4)` in Lemma 2.7;
   - the number `L = √(AB)` in Lemma 3.3;
   - the parametrization `t + i√(1 − t²)` in §3.5;
   - the fact that `|z| ∈ F` for `z ∈ E`.
3. `K` is divisible:
   - Lemma 2.4 uses that a finite image of `K` is trivial;
   - Lemma 3.2 uses that the power maps `u ↦ uⁿ` used in the multiple averages are onto.

For `E_c`, items 1 and 2 hold.

Item 3 holds only for the exponents `±2^k`. The map `u ↦ u²` is onto `K`, because the square roots of `a + bi` with
`a² + b² = 1` have coordinates `±√((1 + a)/2)` and `±√((1 − a)/2)`. The map `u ↦ u³` is not onto, and this is
the impossibility of trisecting angles. For example, `u = (3 + 4i)/5 = (2 + i)/(2 − i)` has valuation 1 at the
prime `2 + i` of `ℤ[i]`. So `u` is not a cube in `ℚ(i)`, and `z³ − u` is irreducible over `ℚ(i)`. Its roots
therefore have degree 3 over `ℚ(i)` and are not constructible.

Two changes make squares enough.

### Change 1: in Lemma 2.4 an element of order 2 is enough

The coset argument of Lemma 2.4 needs one rotation `u ≠ 1` that fixes a coset `d + C` of positive mass. Then
`(u − 1)d ∈ C`, and so `d ∈ C`. The formalization already takes `u = −1` (`coset_fixed_minusOne`). That `−1` fixes
the coset follows from a group-theoretic fact. Let `G` be a commutative group in which every element is a square,
let `g ∈ G` satisfy `g² = 1`, and let `f` map `G` onto a finite group `H`. Then:

- squaring maps `H` onto itself, so it is a bijection of the finite set `H`;
- `f(g)² = f(g²) = 1`, so `f(g) = 1`.

In Lean this is `finite_image_involution`, which replaces OpenAI's `finite_image_trivial`. It is used twice:

- in `fixed_of_finite_orbit` and `fixed_of_positive_fiber`, for the coset argument;
- in `wild_line_null`, the null-set statement about lines in Lemma 2.4. There it shows that `K` is infinite: if `K`
  were finite, the identity map of `K` would send `−1` to 1.

### Change 2: the triple average at exponents −1, 0, 1

Lemma 3.3 of the paper proves the inequality `m(A)m(B) ≤ m(|A − B|)`. It applies the multiple-average lemma
(Lemma 3.2) to the product of characters `ψ_{−L} · U_u ψ_{A−B} · U_{u²} ψ_L`, at exponents 0, 1, 2. Its projection
is then replaced by `g_{−L} · U_u g_{A−B} · U_{u²} g_L`.

In the proof of Lemma 3.2, at (3.9), a common integer is added to the exponents of each term so that none is zero.
The formalization uses one shift for all terms and assumes that every nonzero power map is onto. A shift of 0, 1, 2
gives three consecutive nonzero integers, and one of them is a nonzero multiple of 3. So a single shift with every
shifted exponent onto needs `u ↦ u^{3k}` to be onto for some `k ≠ 0`, that is, cube roots. A shift chosen separately
for each term would avoid cube roots: by 1 when the centred factor sits at 0 or 1, by −3 when it sits at 2. The
proof of (3.5) needs surjectivity only at the exponent of the centred factor and at differences of exponents.

We avoid the shift altogether, by two remarks:

- `U_u` is a unitary, multiplicative operator (it composes with a measure-preserving map), and it commutes with the
  conditional expectation `P`. So `‖P(a · U_u b · U_{u²} c)‖ = ‖P(U_{u⁻¹} a · b · U_u c)‖`.
- For `|u| = 1`, `|−L + (A − B)u + Lu²| = |−Lu⁻¹ + (A − B) + Lu|`.

The same bound therefore follows from an average at exponents −1, 0, 1. In that average the middle factor is fixed,
and nothing is shifted. Write each factor as `v = Pv + (v − Pv)` and expand the product into terms:

- **The term with no centred factor** is `Y`-measurable. It is the product of the projections.
- **A term whose only centred factor is the fixed one**, `b − Pb`, has projection
  `U_{u⁻¹}Pa · U_u Pc · P(b − Pb) = 0`, because the other two factors are `Y`-measurable.
- **A term with a centred factor at exponent ±1.** By conditional independence, its squared norm is
  `∫_R (b′ ⊗ b̄′) · U_{u⁻¹}(a′ ⊗ ā′) · U_u(c′ ⊗ c̄′) dμ_R`, on the relative product `R = X ×_Y X`. This pairs the
  fixed function `b′ ⊗ b̄′` with a product of two factors at exponents −1 and 1. Statement (3.5) of the paper,
  applied on `R`, says that the average of that product tends to 0 in `L²(R)`. Its proof uses three things:
  - `R` is ergodic, by Lemma 3.1;
  - the pair criterion on `R` holds at `d = ±1`;
  - the case `l = 1`, through the power maps `u ↦ u^{±2}`.

  So the average of the squared norm tends to 0.

The power maps used are `u ↦ u^{±1}` and `u ↦ u^{±2}`, so square roots suffice. The other two averages of Lemma 3.3
need no new idea:

- the pair at exponents 0 and 2 is the pair criterion (3.3) with `d = 2`, which uses only the square map;
- for the pair at exponents 1 and −1, the lemma `conditional_integer_pair_energy` now takes the shift as an argument,
  and is called with shift 0. Its power maps are again `u ↦ u^{±1}` and `u ↦ u^{±2}`.

In Lean:

- the new lemmas `centered_conditional_mean_zero_fixed` and `conditional_multiple_mean_zero_fixed` are in
  `Energy.lean`;
- the same file has the triple energy `conditional_integer_triple_energy` at exponents −1, 0, 1;
- the identity it is used with is in `RadialEnergy.lean`.

### Changes common to the three fields

Each field is defined in a new file, `Field.lean`. It proves the facts of items 1 and 2 that the rest of the proof
uses:

- `E` contains the square roots of its elements (`E_sqrt_mem`), and the roots of order `2^a` (`E_root_mem`);
- `E` is algebraic and closed under conjugation, and it contains `i`;
- the real part, the imaginary part and the modulus of an element of `E` are in `F`;
- `F` contains ℚ and is closed under square roots of nonnegative elements.

Countability of `E` is derived from algebraicity, in `Basic.lean`.

The other changed files replace OpenAI's appeals to the algebraic closure by these facts:

- the power maps on `K` are proved onto only for the exponents `2^a`, in `K_pow_surjective` and
  `rot_zsmul_surjective`;
- the good selections of §3.5, built from the rotations in general position of Lemma 3.4, are taken inside `E`, in
  `uniform_selection_in` and `generatedField_le_E`.

The patches are in [`lean/external/openai-five/fields/`](../lean/external/openai-five/fields/README.md).

## The limit of the method

Items 1 and 2 need `F` to be closed under square roots of nonnegative elements, that is, `F` must be a Euclidean field.
For instance, with `B = 1` the number `L = √(AB)` of Lemma 3.3 is `√A`, for every `A ≥ 0` in `F`. The real
constructible numbers form the smallest Euclidean subfield of ℝ, and every such field contains them. So `E_c` is the
smallest field to which the proof applies in this form.

Smaller planes need a different argument. Examples are the plane over the Pythagorean closure of ℚ, or over any
explicitly given field of finite degree.

## References

- OpenAI, *The Euclidean plane is not five-colorable*, OpenAI Math Release preprint, 23 September 2026, and its Lean
  formalization in [openai/math](https://github.com/openai/math) (Apache License 2.0), commit
  `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.
- R. C. Alperin, *A mathematical theory of origami constructions and numbers*, New York J. Math. 6 (2000), 119–133.
- D. A. Cox, *Galois Theory*, 2nd ed., Wiley, 2012, Chapter 10 (constructible numbers).
- N. G. de Bruijn, P. Erdős, *A colour problem for infinite graphs and a problem in the theory of relations*,
  Indag. Math. 13 (1951), 371–373 ([doi](https://doi.org/10.1016/S1385-7258(51)50053-7)).
- K. G. Fischer, *A planar geometric graph of chromatic number four*, Congr. Numer. 104 (1994), 73–79.
