# A local–global question for the chromatic number of planes over number fields

*Working note, 3 October 2026. Theorem A is proved here (its "if" half was already in
`notes/local_colourings.md`); the rest is a question, with the evidence we have and the
tests under way. A separate agent refereed this note (its findings are applied below);
nobody outside the project has.*

Let `F` be a number field and give `F²` the unit-distance relation
`(x − x′)² + (y − y′)² = 1`. This is algebraic, so the same relation makes sense on
`F_v²` for every completion `F_v` of `F`. A proper colouring of `F_v²` restricts to one
of `F²`, so

> `χ(F²) ≤ χ(F_v²)` for every place `v`, and `χ(F²) ≤ χ(ℝ²) ≤ 7` when `F` is real.

Every upper bound below 7 that we know for a plane over a number field comes this way,
from a colouring that is *locally constant* at one finite place: reduction modulo a
prime power (Moorhouse's Corollary 8.3; Madore's coset colourings; Fischer's colourings,
2-adic in his Theorem 10 and 3-adic in his Theorem 9; the 11-adic levels for `ℚ(√47)`,
`notes/quadratic_planes.md` §6). The bound 7 comes from the real place. This note asks
whether that is the whole story.

## 1. The question

For a place `v` of `F` let `χ_loc(v)` be:
- at a finite place where `−1` is a square in `F_v` (the form `x² + y²` is isotropic),
  and at a complex place: `∞`. No locally constant colouring with finitely many colours
  exists at such a finite place (`docs/research-log.md`, "The other places of ℚ(√47)").
- at a finite place where `x² + y²` is anisotropic: the least `k` such that `F_v²` has
  a locally constant proper `k`-colouring. With `w` the place of `F(i)` above `v`, the
  unit vectors are the norm-one elements of `O_w`, and this is the least chromatic
  number of a finite level `Cay(O_w/𝔪_w^r, T_r)`. (At 2-adic places `O_w` can be larger
  than `O_v[i]`: over `ℚ₂(√3)`, `(√3/2, 1/2)` is a unit vector.)
- at a real place: the Borel chromatic number `χ_B(ℝ²)`, which is 5, 6 or 7
  (Falconer's measurable bound 5; the hexagonal 7-colouring is Borel).

Then `χ(F²) ≤ min_v χ_loc(v)` for every real number field `F`.

> **Question LG.** Is `χ(F²) = min_v χ_loc(v)` for every real number field `F`?

The part we can test is, for each `k`: *if no finite place of `F` has a locally constant
`k`-colouring, does `F²` contain a finite unit-distance graph with no proper
`k`-colouring?* (`LG_k`). For `k ≤ 4` this does not involve the real place, since
`χ_B(ℝ²) ≥ 5`. The real place is the speculative part: with it, LG would settle that
`χ(ℝ²) = χ_B(ℝ²)` (§4), a major open question.

Colourings through several places at once colour a categorical product of levels, so
they beat every single place only if the product needs fewer colours than its factors.
Moorhouse (2010, §9) noted that this needs a counterexample to Hedetniemi's conjecture.
Such counterexamples exist for large `k` (Shitov, 2019), so for large `k` stating LG with
single places is a real restriction. For `k = 3` the product of two non-3-colourable
graphs is not 3-colourable (El-Zahar and Sauer, 1985), so single places suffice. For
`k = 4` the question is open in general; Zhu's fractional theorem
(`χ_f(G × H) = min(χ_f(G), χ_f(H))`) excludes products of places whose levels all have
Hoffman ratio below `1/4`.

## 2. Theorem A: two colours

> **Theorem A.** Let `F` be a number field. Then `χ(F²) = 2` if and only if some prime of
> `F` above 2 ramifies in `F(i)`.

So `LG_2` holds for every number field. The "if" half is the corollary of Proposition A
in `notes/local_colourings.md` §3 (a ramified prime over 2 gives a locally constant
2-colouring; Madore's Proposition 3.9 is the case `ℚ(√2)`, and he did not state the
general case). The converse, Step 3 below, is new as far as we know. Together they
contain Fischer's Theorem 8 (`ℚ(√N)²` is 2-colourable exactly when `N ≢ 3 (mod 4)`), the
"if" part of Johnson (1987), which we know only through Payne's account, and Moorhouse's
Theorems 7.1 (odd degree: some `v | 2` has odd degree over `ℚ₂`, so `F_v(i)/F_v` is
ramified), 8.5 and 8.6 and Lemma 8.4.

*Proof.* Let `L = F[t]/(t² + 1)`: the field `F(i)` if `i ∉ F`, and `F × F` if `i ∈ F`.
Identify `F²` with `L` by `(x, y) ↦ x + ty`; then `x² + y²` is the norm `N(z) = z z̄`, and
the unit vectors are `T = {z ∈ L : N(z) = 1}`. `T` is a group under multiplication, so
its additive span `A = ℤ[T]` is a subring of `L`, and `F²` is a disjoint union of
translates of `A`, each a copy of `Cay(A, T)`. Every `u ∈ T` is a unit of `A`
(`u⁻¹ = ū ∈ T`).

*Step 1: `χ(F²) = 2` if and only if there is a ring homomorphism `φ : A → 𝔽₂`.* If
`Cay(A, T)` is bipartite, let `φ(a)` be the parity of the length of any walk from 0 to
`a` (since `−1 ∈ T`, every element of `A` is a sum of unit vectors). It is additive, and
multiplicative because a product of an `m`-term and an `n`-term sum of unit vectors is an
`mn`-term sum of unit vectors. Conversely, a ring homomorphism sends each `u ∈ T`, a unit
of `A`, to `1 ∈ 𝔽₂^×`, so `φ(a + u) = φ(a) + 1` and `φ` is a 2-colouring.

*Step 2: a ramified prime gives a colouring.* Let `v | 2` ramify in `L`, with `w` above
it; `w` is the only place of `L` above `v`, so conjugation preserves `|·|_w`. The
extension `L_w/F_v` is ramified, so the residue fields agree and conjugation acts
trivially on them. A unit vector `u` of `F_v²` therefore lies in `O_w`
(`|u|_w² = |u ū|_w = 1`), and `1 = u ū ≡ u²` modulo `𝔪_w`; in characteristic 2 this
gives `u ≡ 1 (mod 𝔪_w)`. Pick an additive map `λ` from the residue field to `𝔽₂` with
`λ(1) = 1`. Unit steps preserve the classes of `L_w` modulo `O_w`; fix a representative
`r` in each class and colour `z` by `λ((z − r) mod 𝔪_w)`. Each unit step changes the
colour. This is a locally constant 2-colouring of `F_v² ≅ L_w`, and so of `F²`.

*Step 3: a colouring gives a ramified prime.* Let `φ : A → 𝔽₂` be a ring homomorphism.

If `i ∈ F`, then `T = {(s, 1/s) : s ∈ F^×}` in `L = F × F`. For `s ≠ 0, 1` the sum of
`(s, 1/s)` and `(1 − s, 1/(1 − s))` is `(1, 1/(s(1 − s)))`; subtracting two of these
(`s = 2, 3`) gives `(0, x) ∈ A` with `x = −1/3 ≠ 0`. The set `{x : (0, x) ∈ A}` is an
ideal of the second projection of `A`, which is `F`; so `A ⊇ 0 × F`, by symmetry
`A = F × F`, and there is no homomorphism to `𝔽₂` (2 is invertible in `F`). So `i ∉ F`
and `L` is a field.

By Zorn's lemma the pairs `(B, ψ)` of a subring `B ⊇ A` of `L` and a homomorphism `ψ`
extending `φ` into an algebraic closure of `𝔽₂` have a maximal element, and by
Chevalley's theorem (Atiyah–Macdonald, Theorem 5.21) its ring `V` is a valuation ring of
`L` whose maximal ideal meets `A` in `ker φ`. `V` contains `ℤ` and is integrally closed,
so `V ⊇ O_L`, and `V ≠ L` because `2 ∈ ker φ`; so `V = (O_L)_𝔓` for a prime `𝔓 ∋ 2`,
that is, `V` is the valuation ring of a place `w | 2` of `L`, and every `u ∈ T` satisfies
`u ≡ 1 (mod 𝔪_w)`. Let `v` be the place of `F` below `w`. The map `p ↦ (p + i)/(p − i)`
is a homeomorphism from `F_v` minus the roots of `p² + 1` onto `T(F_v) ∖ {1}` and maps
`F` onto `T ∖ {1}`, so `T` is dense in `T(F_v)`.
- If `v` splits in `L`, then `T(F_v) ≅ F_v^×` through `z ↦ z_w`, which is unbounded, so
  some element of `T` is not integral at `w`: impossible.
- If `v` is inert, conjugation acts on the residue field `𝔽_{q²}` of `w` as `x ↦ x^q`,
  and the unit vectors `y/ȳ` (`y ∈ O_w^×`) reduce onto the group of `(q + 1)`-th roots of
  unity (`y ↦ y^{1−q}`). This group has `q + 1 ≥ 3` elements and each residue class is
  open, so some element of `T` is not `≡ 1`: impossible.

So `v` ramifies in `L`. ∎

**Biquadratic fields.** For `F = ℚ(√a, √b)` all primes above 2 are conjugate, and none
ramifies in `F(i)` exactly when the subgroup `⟨a, b⟩` of `ℚ₂^×/ℚ₂^×²` contains `−1` (then
`i ∈ F_v`) or `−5` (then `F_v(i) = F_v(√5)` is unramified), that is, when one of the
three quadratic subfields `ℚ(√c)` has `c ≡ 3 (mod 4)`. So `χ(ℚ(√a, √b)²) = 2` unless
a quadratic subfield already needs three colours: here Theorem A gives nothing beyond
the subfields.

**Checks.** `data/quadratic_planes/scripts/two_colour_criterion.gp` (PARI/GP) computes the
criterion from the relative discriminant of `F(i)/F`. It agrees with Fischer's Theorem 8
for all 242 squarefree `d ≤ 400`, with Moorhouse's Theorem 7.1 on nine fields of degree
3, 5 and 7, and with the biquadratic statement on 1 522 fields. These are cases known
before. The referee also tested the new half on twelve quartic fields where it predicts
`χ ≥ 3` although no quadratic subfield forces it (eight with Galois group `S₄`, some of
them with both split and inert primes above 2, and fields of type `C₄` and `D₄` that
contain the 2-colourable `ℚ(√5)`, `ℚ(√2)` or `ℚ(√53)`): in each, integer relations
among unit vectors `z/z̄` gave a closed walk of odd length, checked exactly in PARI
(length 5 for the roots of `α⁴ − 2α³ − 2α − 3` and of `α⁴ − 6α² + 4α + 2`, length 7 for
`ℚ(i)`), and fields predicted 2-colourable gave only even relations. `data/quadratic_planes/scripts/odd_walks.gp` repeats this for three of these fields
(with fixed seeds) and three 2-colourable controls.

## 3. Three colours: the evidence

For real quadratic `F = ℚ(√d)` (squarefree `d`), the places with a locally constant
3-colouring are known:
- `d ≢ 3 (mod 4)`: one above 2 (Theorem A), so `χ = 2`;
- `d ≡ 3 (mod 4)`, `d ≢ 2 (mod 3)`: one above 3 (Fischer, Theorem 9; Moorhouse,
  Corollary 8.3), and `χ = 3`;
- `d ≡ 11 (mod 12)`: none. Above 2 the plane is isotropic (`d ≡ 7 (mod 8)`) or contains
  the 76-vertex graph of `ℚ(√11)` (`d ≡ 3 (mod 8)`, since then `ℚ₂(√d) = ℚ₂(√11)`); 3 is
  inert; at `ℚ₇` every level needs four colours (`χ(ℚ₇²) = 4`), at `ℚ₇(√−7)` the graph
  of `ℚ(√35)` needs four, and levels 1 to 3 of `ℚ₇(√7)` have no 3-colouring
  (`data/quadratic_planes/scripts/ramified_levels.py`, CaDiCaL; level 3 has 117 649
  points and 23 059 204 edges); every other anisotropic place has Hoffman ratio below
  `1/3`. (For `ℚ₇(√7)`, which occurs for `d = 203`, levels 4 and up are not checked.)

So `LG_3` predicts `χ(ℚ(√d)²) ≥ 4` for every `d ≡ 11 (mod 12)` (for `d = 203`, unless a
deep level at 7 is 3-colourable). This is proved for the 27 values in
`notes/quadratic_planes.md` (11, 23, 35, 47, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251,
263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959), with certified graphs,
and open for `d = 83, 107, 143, 167, 203, …` (the growth searches have not yet found
graphs there). No real quadratic field is known where the prediction fails.

Other exact values agree with `min_v χ_loc(v)`: `χ(ℚ²) = 2`; Madore's
`χ(ℚ(√3)²) = χ(ℚ(√7)²) = 3`; our `χ(ℚ(√2, √3)²) = 4` and Fischer's
`χ(ℚ(√3, √11)²) = 4`, whose upper bounds are reductions at one place. In higher
dimension the same principle matches `χ(ℚ³) = 2 = χ(ℚ₂³)` and
`χ(ℚ⁴) = 4 = χ(ℚ₂⁴)` (Benda–Perles, as cited by Madore; `papers/padic-planes`), where
every odd place is isotropic.

## 4. Four colours: what the question predicts

A real field has a locally constant 4-colouring only at an anisotropic place whose
levels have Hoffman ratio at least `1/4`, that is above 2, 3, 7, 11 or 19 (research log,
"The other places of ℚ(√47)"; `data/quadratic_planes/scripts/hoffman_padic.py`). The
margins are thin near the threshold: the level-1 ratio is 0.249077 at 23 and 0.2499973
at 31 (Hoffman bound 4.000044), which matters for `ℚ(√2, √31)`, where 31 ramifies. Call
`F` *admissible* if `√3 ∉ F` and every place of `F` above 2, 3, 7, 11 and 19 contains
`i`. An admissible field has no locally constant 4-colouring at any place and no unit
triangle, so:

> **Consequence.** If `LG_4` holds for `ℚ(√167)`, there is a triangle-free 5-chromatic
> unit-distance graph in the plane.

Soifer asked for such a graph in a 2019 problem paper; chapter 55 of his *New
Mathematical Coloring Book* (2024) is titled "Triangle-Free 5-Chromatic Unit Distance
Graphs" (we have read only its abstract). MathWorld cites de Grey, "A 5-chromatic,
triangle-free unit-distance graph in ℝ³ with 61 vertices", Geombinatorics 35 (2026). We
know of no planar example.

The admissible quadratic fields are `ℚ(√d)` with `d ≡ 7 (mod 8)` a non-residue modulo 3,
7, 11 and 19 and prime to them: `d = 167, 887, 1055, 1319, 1823, …`. `ℚ(√167)` is the first
of them: the smallest `d` for which every place is shown to have no locally constant
4-colouring. (For `d = 47` and `143` the places above 11, and for 47 also above 19, are
undecided.) Its places above 2, 3, 7, 11, 19 and 31 are isotropic, and the anisotropic
ones lie above 167 (ramified) and the split primes `p ≡ 3 (mod 4)`, which are 23, 43,
59, 67, …; all have Hoffman ratio below `1/4` at every level. Fischer (1990, Example 4)
had already shown that `ℚ(√167)²` has no additive `k`-colouring for `k < 11`. Admissible
biquadratic fields include `ℚ(√2, √31)`, `ℚ(√2, √47)`, `ℚ(√10, √38)`, `ℚ(√2, √55)` and
`ℚ(√11, √13)` (`data/quadratic_planes/scripts/admissible.py`).

With the real place included, the question also predicts `χ(ℝ²) = χ_B(ℝ²)`: by
Dirichlet's theorem there are primes `d ≡ 7 (mod 8)` that are non-residues modulo every
prime `p ≡ 3 (mod 4)` below 100 (the smallest is `d = 635 087`), and every finite place
of `ℚ(√d)` then has `χ_loc ≥ 7` (Hoffman's bound `1 + (q + 1)/(2√q)` is at least 6.12
for `q ≥ 103`, and Zhu's theorem covers colourings through several places), so LG would
give `χ(ℚ(√d)²) = χ_B(ℝ²) ≥ χ(ℝ²)`. Moorhouse (2010, introduction, p. 1) wrote that he
suspected `χ(ℝ²) = 7`.

## 5. Tests under way

- `ℚ(√167)` at three colours: colouring-guided growth (`grow3r.py`, `QD = 167`) with the
  gate-open denominators 2784, 3480 and six others. A graph with no 3-colouring would
  confirm `LG_3` in Moorhouse's class `d ≡ 167 (mod 168)`. So far the runs stop, or time
  out, with 3-colourable graphs of 10 000 to 43 000 points.
- `ℚ(√167)` at four colours: points at distance exactly `1/2` lie 4 steps from the
  origin for `D = 2784` and `D = 4176` (e.g. `m = (−1058, 0, 0, −70)/2784`), so a pair
  forced to one colour at distance `1/2` closes a 5-chromatic graph after a half-turn
  (`m ↦ −m`, `|2m| = 1`). The growth stalls at 4 colours with these sparse direction
  sets (about 34 000 points, average degree 4.5, against 10 for `ℚ(√47)` at `D = 240`):
  the unit vectors of `ℚ(√167)` have larger denominators. `χ(H₂₅) ≥ 5`
  (`hyperbola_plane.py`), so the places above 5 are no gate there at four colours.
- `ℚ(√2, √47)`: its unit vectors include those of `ℚ(√47)`, and `√2` makes the places
  above 11 and 19 isotropic. But with the products of `ℚ(√47)`'s 108 directions
  (`D = 240`) and the units of `ℚ(√2)` with denominator 3, a new direction meets the
  72 997-point graph over `ℚ(√47)` (share C's colouring) in at most two points, so none
  of them blocks a colouring: the new directions add little at this scale.
- `ℚ(√47)`: the 19-adic level 2 (130 321 points) is with kissat; a 4-colouring would
  give `χ(ℚ(√47)²) = 4`.

## References

- M. F. Atiyah and I. G. Macdonald, *Introduction to Commutative Algebra*, Theorem 5.21.
- M. Benda and M. Perles, Colorings of metric spaces, Geombinatorics 9 (2000) 113–126
  (not read; cited through Madore).
- A. D. N. J. de Grey, A 5-chromatic, triangle-free unit-distance graph in ℝ³ with 61
  vertices, Geombinatorics 35 (2026) (not read; cited through MathWorld,
  "Grötzsch Graph").
- M. El-Zahar and N. Sauer, The chromatic number of the product of two 4-chromatic
  graphs is 4, Combinatorica 5 (1985) 121–126.
- K. J. Falconer, The realization of distances in measurable subsets covering `ℝⁿ`,
  J. Combin. Theory Ser. A 31 (1981) 184–189.
- K. G. Fischer, Additive K-colorable extensions of the rational plane, Discrete Math.
  82 (1990) 181–195.
- P. D. Johnson Jr., Two-colorings of real quadratic extensions of `ℚ²` that forbid many
  distances, Congr. Numer. 60 (1987) 51–58 (not read; known through Payne,
  arXiv:0707.1177).
- D. A. Madore, The Hadwiger–Nelson problem over certain fields, arXiv:1509.07023 (2015).
- G. E. Moorhouse, On the chromatic numbers of planes, draft of 3 March 2010, §§7–9.
- Y. Shitov, Counterexamples to Hedetniemi's conjecture, Ann. of Math. 190 (2019)
  663–667.
- A. Soifer, *The New Mathematical Coloring Book*, 2nd ed., Springer, 2024, chapter 55
  (abstract only).
- X. Zhu, The fractional version of Hedetniemi's conjecture is true, European J.
  Combin. 32 (2011) 1168–1175.

Not yet read, and needed before claiming novelty: Benda–Perles (2000) and Johnson's
survey of their problems (Geombinatorics 9 (2000) 170–179).
