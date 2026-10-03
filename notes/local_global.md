# A local–global question for the chromatic number of planes over number fields

*Working note, 3 October 2026. Theorem A is proved here; the rest is a conjecture with
the evidence we have and the tests under way. Nobody outside the project has refereed it.*

Let `F` be a number field and give `F²` the unit-distance relation
`(x − x′)² + (y − y′)² = 1`. This is algebraic, so the same relation makes sense on
`F_v²` for every completion `F_v` of `F`. A proper colouring of `F_v²` restricts to one
of `F²`, so

> `χ(F²) ≤ χ(F_v²)` for every place `v`, and `χ(F²) ≤ χ(ℝ²) ≤ 7` when `F` is real.

Every upper bound we know for a plane over a number field comes this way, and in fact
from a colouring that is *locally constant* at one place: reduction modulo a prime
power (Moorhouse's Corollary 8.3, Madore's coset colourings, Fischer's additive
colourings, which are 2-adic, and the 11-adic levels for `ℚ(√47)`,
`notes/quadratic_planes.md` §6). This note asks whether that is the whole story.

## 1. The question

For a place `v` of `F` let `χ_loc(v)` be:
- at a finite place where `−1` is a square in `F_v` (the form `x² + y²` is isotropic):
  `∞`. No locally constant colouring with finitely many colours exists there
  (`docs/research-log.md`, "The other places of ℚ(√47)").
- at a finite place where `x² + y²` is anisotropic: the least `k` such that `F_v²` has
  a locally constant proper `k`-colouring, that is, the least chromatic number of a
  finite level `Cay((O_v/π^r)², T_r)`.
- at a real place: the Borel chromatic number `χ_B(ℝ²)`, which is 5, 6 or 7
  (Falconer's measurable bound 5; the hexagonal 7-colouring is Borel).

Then `χ(F²) ≤ min_v χ_loc(v)` for every real number field `F`.

> **Conjecture LG.** For every real number field `F`, `χ(F²) = min_v χ_loc(v)`.

The weaker form we can test is: *if no place of `F` has a locally constant `k`-colouring,
then `F²` contains a finite unit-distance graph with no proper `k`-colouring* (`LG_k`).
`LG_k` for `k ≤ 4` does not involve the real place, since `χ_B(ℝ²) ≥ 5`.

Colourings through several places at once colour a categorical product of levels.
Moorhouse noticed that they can beat a single place only through a counterexample to
Hedetniemi's conjecture (his §9). For `k ≤ 4` and places whose Hoffman ratio is below
`1/k` this is excluded by Zhu's fractional Hedetniemi theorem
(`χ_f(G × H) = min(χ_f(G), χ_f(H))`), as in the research log; so the conjecture is stated
with single places.

## 2. Theorem A: two colours

> **Theorem A.** Let `F` be a number field. Then `χ(F²) = 2` if and only if some prime of
> `F` above 2 ramifies in `F(i)`.

So `LG_2` holds for every number field: a prime `v | 2` ramified in `F(i)` gives a
locally constant 2-colouring of `F_v²` (below), and if there is none, `F²` has an odd
cycle. Theorem A contains the earlier cases: `ℚ(√d)` with `d ≢ 3 (mod 4)` (Johnson 1987;
Fischer 1990, Theorem 8; Moorhouse 2010, Lemma 8.4 and Theorem 8.5), and fields of odd
degree (Moorhouse, Theorem 7.1: some `v | 2` has odd degree over `ℚ₂`, so `F_v(i)/F_v`
is ramified).

*Proof.* Let `L = F[t]/(t² + 1)`: the field `F(i)` if `i ∉ F`, and `F × F` if `i ∈ F`.
Identify `F²` with `L` by `(x, y) ↦ x + ty`; then `x² + y²` is the norm `N(z) = z z̄`, and
the unit vectors are `T = {z ∈ L : N(z) = 1}`. `T` is a group under multiplication, so
its additive span `A = ℤ[T]` is a subring of `L`, and `F²` is a disjoint union of
translates of `A`, each a copy of `Cay(A, T)`.

*Step 1: `χ(F²) = 2` if and only if there is a ring homomorphism `φ : A → 𝔽₂`.* If
`Cay(A, T)` is bipartite, let `φ(a)` be the parity of the length of any walk from 0 to
`a` (since `−1 ∈ T`, every element of `A` is a sum of unit vectors). It is additive, and
multiplicative because a product of an `m`-term and an `n`-term sum of unit vectors is an
`mn`-term sum of unit vectors. Conversely `φ(a + u) = φ(a) + 1` for `u ∈ T`, so `φ` is a
2-colouring.

*Step 2: a ramified prime gives a colouring.* Let `v | 2` ramify in `L`, with `w` above
it. The extension `L_w/F_v` is ramified, so the residue fields agree and conjugation
acts trivially on them. A unit vector `u` of `F_v²` lies in `O_w` (conjugation preserves
`|·|_w` and `|u|_w² = |N(u)|_w = 1`), and `1 = u ū ≡ u²` modulo `𝔪_w`; in characteristic 2
this gives `u ≡ 1 (mod 𝔪_w)`. Pick an additive map `λ` from the residue field to `𝔽₂`
with `λ(1) = 1`. Unit steps preserve the classes of `L_w` modulo `O_w`; fix a
representative `r` in each class and colour `z` by `λ((z − r) mod 𝔪_w)`. Each unit step
changes the colour. This is a locally constant 2-colouring of `F_v² ≅ L_w`, and so of
`F²`.

*Step 3: a colouring gives a ramified prime.* Let `φ : A → 𝔽₂` be a ring homomorphism.

If `i ∈ F`, then `T = {(s, 1/s) : s ∈ F^×}` in `L = F × F`. For `s ≠ 0, 1` the sum of
`(s, 1/s)` and `(1 − s, 1/(1 − s))` is `(1, 1/(s(1 − s)))`; subtracting two of these
(`s = 2, 3`) gives `(0, x) ∈ A` with `x = −1/3 ≠ 0`. The set `{x : (0, x) ∈ A}` is an ideal of the second projection of `A`,
which is `F`; so `A ⊇ 0 × F`, by symmetry `A = F × F`, and there is no homomorphism to
`𝔽₂` (2 is invertible in `F`). So `i ∉ F` and `L` is a field.

By Chevalley's extension theorem (Atiyah–Macdonald, Theorem 5.21) there is a valuation
ring `V ⊇ A` of `L` whose maximal ideal meets `A` in `ker φ`. Since `2 ∈ ker φ`, `V ≠ L`,
so `V` is the valuation ring of a place `w | 2` of `L`, and every `u ∈ T` satisfies
`u ≡ 1 (mod 𝔪_w)`. Let `v` be the place of `F` below `w`. The circle is a rational curve
with a rational point (`p ↦ (p + i)/(p − i)`), so `T` is dense in `T(F_v)`.
- If `v` splits in `L`, then `T(F_v) ≅ F_v^×` through `z ↦ z_w`, which is unbounded, so
  some element of `T` is not integral at `w`: impossible.
- If `v` is inert, conjugation acts on the residue field `𝔽_{q²}` of `w` as `x ↦ x^q`,
  and the unit vectors `y/ȳ` (`y ∈ O_w^×`) reduce onto the group of `(q + 1)`-th roots of
  unity (`y ↦ y^{1−q}`). This group has `q + 1 ≥ 3` elements and each residue class is
  open, so some element of `T` is not `≡ 1`: impossible.

So `v` ramifies in `L`. ∎

**Checks.** `data/quadratic_planes/scripts/two_colour_criterion.gp` (PARI/GP) computes the
criterion from the relative discriminant of `F(i)/F`. It agrees with "`χ = 2` iff
`d ≢ 3 (mod 4)`" for all 242 squarefree `d ≤ 400`, and with Moorhouse's Theorem 7.1 on
nine fields of degree 3, 5 and 7. For biquadratic `ℚ(√a, √b)` it predicts `χ ≥ 3` only
when one of the three quadratic subfields already has `χ ≥ 3`.

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

So `LG_3` predicts `χ(ℚ(√d)²) ≥ 4` for every `d ≡ 11 (mod 12)`. This is proved for the
27 values in `notes/quadratic_planes.md` (11, 23, 35, 47, 59, 71, 95, 119, 131, 155, 179,
191, 239, 251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959), with
certified graphs, and open for `d = 83, 107, 143, 167, 203, …` (the growth searches
have not yet found graphs there). No real quadratic field is known where the prediction
fails.

Other exact values agree with `min_v χ_loc(v)`: `χ(ℚ²) = 2`; Madore's
`χ(ℚ(√3)²) = χ(ℚ(√7)²) = 3`; our `χ(ℚ(√2, √3)²) = 4` and Fischer's
`χ(ℚ(√3, √11)²) = 4`, whose upper bounds are reductions at one place. In higher
dimension the same principle matches `χ(ℚ³) = 2 = χ(ℚ₂³)` and
`χ(ℚ⁴) = 4 = χ(ℚ₂⁴)` (Benda–Perles; `papers/padic-planes`), where every odd place is
isotropic.

## 4. Four colours: what the conjecture predicts

A real field has a locally constant 4-colouring only at an anisotropic place whose
levels have Hoffman ratio at least `1/4`, that is (research log, "The other places of
ℚ(√47)") above 2, 3, 7, 11 or 19. Call `F` *admissible* if `√3 ∉ F` and every place of
`F` above 2, 3, 7, 11 and 19 contains `i`. An admissible field has no locally constant
4-colouring at any place and no unit triangle, so:

> **Consequence.** `LG_4` implies that there is a triangle-free 5-chromatic
> unit-distance graph in the plane.

That is an open problem (Soifer, *The New Mathematical Coloring Book*, 2024, chapter
"Triangle-Free 5-Chromatic Unit Distance Graphs"; de Grey found one in `ℝ³`). The
admissible quadratic fields are `ℚ(√d)` with `d ≡ 7 (mod 8)` a non-residue modulo 3, 7,
11 and 19 and prime to them: `d = 167, 887, 1055, 1319, 1823, …`. `ℚ(√167)` is the first
real quadratic field with no local 4-colouring at all: the places above 2, 3, 7, 11, 19
and 31 are isotropic, and the anisotropic ones lie above 167 (ramified) and the split
primes `p ≡ 3 (mod 4)`, which are 23, 43, 59, 67, …; all have Hoffman ratio below `1/4`
at every level. Fischer (1990, Example 4) had already shown that `ℚ(√167)²` has no additive
`k`-colouring for `k < 11`. Admissible biquadratic fields include `ℚ(√2, √31)`,
`ℚ(√2, √47)`, `ℚ(√10, √38)`, `ℚ(√2, √55)` and `ℚ(√11, √13)`
(`data/quadratic_planes/scripts/admissible.py`). By contrast `ℚ(√47)` is exposed at 11
and 19, where finite levels might still be 4-colourable (level 2 at 19 is being decided).

With the real place included, the conjecture also predicts `χ(ℝ²) = χ_B(ℝ²)`: by
Dirichlet's theorem there are primes `d ≡ 7 (mod 8)` that are non-residues modulo every
prime `p ≡ 3 (mod 4)` below 100, and every finite place of `ℚ(√d)` then has
`χ_loc ≥ 7` (Hoffman's bound `1 + (q + 1)/(2√q)`), so `LG` would give
`χ(ℚ(√d)²) = χ_B(ℝ²) ≥ χ(ℝ²)`. Moorhouse wrote that he suspected `χ(ℝ²) = 7`.

## 5. Tests under way

- `ℚ(√167)` at three colours: colouring-guided growth (`grow3r.py`, `QD = 167`) with the
  gate-open denominators 2784, 3480 and six others. A graph with no 3-colouring would
  confirm `LG_3` in Moorhouse's class `d ≡ 167 (mod 168)`.
- `ℚ(√167)` at four colours: points at distance exactly `1/2` lie 4 steps from the
  origin for `D = 2784` and `D = 4176` (e.g. `m = (−1058, 0, 0, −70)/2784`), so a pair
  forced to one colour at distance `1/2` closes a 5-chromatic graph after a half-turn
  (`m ↦ −m`, `|2m| = 1`). The growth stalls at 4 colours with these sparse direction
  sets (about 34 000 points, average degree 4.5, against 10 for `ℚ(√47)` at `D = 240`):
  the unit vectors of `ℚ(√167)` have larger denominators.
- `ℚ(√47)`: the 19-adic level 2 (130 321 points) is with kissat; a 4-colouring would
  give `χ(ℚ(√47)²) = 4`, as `LG` then predicts.

## References

- M. F. Atiyah and I. G. Macdonald, *Introduction to Commutative Algebra*, Theorem 5.21.
- K. J. Falconer, The realization of distances in measurable subsets covering `ℝⁿ`,
  J. Combin. Theory Ser. A 31 (1981) 184–189.
- K. G. Fischer, Additive K-colorable extensions of the rational plane, Discrete Math.
  82 (1990) 181–195.
- P. D. Johnson Jr., Two-colorings of real quadratic extensions of `ℚ²` that forbid many
  distances, Congr. Numer. 60 (1987) 51–58.
- D. A. Madore, The Hadwiger–Nelson problem over certain fields, arXiv:1509.07023 (2015).
- M. Benda and M. Perles, Colorings of metric spaces, Geombinatorics 9 (2000) 113–126.
- G. E. Moorhouse, On the chromatic numbers of planes, draft of 3 March 2010, §§7–9.
- X. Zhu, The fractional version of Hedetniemi's conjecture is true, European J.
  Combin. 32 (2011) 1168–1175.
- A. Soifer, *The New Mathematical Coloring Book*, Springer, 2024.
