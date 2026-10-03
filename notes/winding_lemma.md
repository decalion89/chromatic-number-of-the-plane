# Three colours and characters: the winding lemma

*Working note, 3 October 2026. Steps 1 and 2 of the proof below are the discrete winding number of Krebs and Sankar
(arXiv:2410.11028, J. Combin. Theory Ser. B, 2026: Definition 3.5, Propositions 3.3 and 3.6, Remark 3.7, and the
parity used in the proof of their Theorem 4.1). The averaging step and the characterisation it gives, Theorem W, are
new as far as we know (§5 lists what we checked). An internal referee (a separate AI agent, with its own checker and
tests) went through the note, and its findings are applied. Nobody outside the project has refereed it.*

## 1. The statement

Let `Γ` be an abelian group, written additively, and `S = −S` a subset of `Γ ∖ {0}`. The Cayley graph
`Cay(Γ, S)` has vertex set `Γ` and an edge `{g, g + s}` for every `g ∈ Γ` and `s ∈ S`. A *character* is a
homomorphism `ξ : Γ → ℝ/ℤ`. For `x ∈ ℝ/ℤ` write `‖x‖` for the distance to `0`.

> **Theorem W.** Let `k ≥ 1`. `Cay(Γ, S)` has a homomorphism to the odd cycle `C_{2k+1}` if and only if some
> character `ξ` satisfies `ξ(s) ∈ [k/(2k+1), (k+1)/(2k+1)]` (modulo 1) for every `s ∈ S`.
>
> In particular (`k = 1`, `C₃ = K₃`): **`Cay(Γ, S)` is 3-colourable if and only if some character maps every
> `s ∈ S` into the closed arc `[1/3, 2/3]`.**

A homomorphism to `C_{2k+1}` is the same as a circular `(2k+1)/k`-colouring. Put
`κ(S) = sup_ξ inf_{s∈S} ‖ξ(s)‖`; the supremum is attained, because the dual group is compact and
`ξ ↦ inf_s ‖ξ(s)‖` is upper semicontinuous. The theorem says:

> `χ_c(Cay(Γ, S)) ≤ 2 + 1/k` if and only if `κ(S) ≥ k/(2k+1)`; in particular `χ ≤ 3` if and only if `κ(S) ≥ 1/3`.

The "if" direction is the classical bound `χ_c ≤ 1/κ` (Zhu's multiplier method for distance graphs, in D. D.-F. Liu,
"From Rainbow to Lonely Runner: A Survey on Coloring Parameters of Distance Graphs", §2; Alon, Eur. J. Combin. 2013,
Theorem 3.2, for Cayley graphs of `ℤ/p`). The content of Theorem W is the converse at the thresholds `2 + 1/k`,
which include 3.

**Recurrence: Katznelson's question for three colours.** In the language of recurrence, the case `k = 1` reads:
`S` is a set of 3-chromatic recurrence (every 3-colouring of `Γ` has two points of one colour that differ by an
element of `S`) if and only if `S` meets the Bohr set `{g : ‖ξ(g)‖ < 1/3}` of every single character `ξ`. A set of
Bohr recurrence meets every Bohr neighbourhood of 0 (every set `{g : ‖ξ₁(g)‖ < ε, …, ‖ξ_m(g)‖ < ε}` with characters
`ξⱼ` and `ε > 0`; Griesmer, arXiv:2108.02190, Definition 1.2), in particular these. So:

> **Corollary (three colours of Katznelson's question).** In every abelian group, every set of Bohr recurrence is a
> set of 3-chromatic recurrence. For `ℕ`: for every 3-colouring `ℕ = A₁ ∪ A₂ ∪ A₃` there is an `α` such that
> `(A₁ − A₁) ∪ (A₂ − A₂) ∪ (A₃ − A₃) ⊇ {n ∈ ℕ : ‖nα‖ < 1/3}`; in particular this union is a Bohr₀ set. More
> generally, every 3-colourable Cayley graph of an abelian group has a 3-colouring whose colour classes are the preimages
> `ξ⁻¹[0, 1/3)`, `ξ⁻¹[1/3, 2/3)`, `ξ⁻¹[2/3, 1)` under a single character `ξ`.

For `ℕ`, let `S = ℤ ∖ ⋃ᵢ (Aᵢ − Aᵢ)`; then `S = −S` (each `Aᵢ − Aᵢ` is symmetric) and `0 ∉ S` (some `Aᵢ` is
nonempty). Every finite piece of `Cay(ℤ, S)` can be translated into `ℕ`, where the given
colouring is proper on it, so `Cay(ℤ, S)` is 3-colourable (de Bruijn–Erdős). Theorem W then gives `α` with
`‖sα‖ ≥ 1/3` for all `s ∈ S`, that is, `{n : ‖nα‖ < 1/3}` misses `S`.

For `Γ = ℤ` no invariant mean is needed. If `S` is finite with largest element `L` and `c` is a 3-colouring of
`Cay(ℤ, S)`, two windows `[i, i + L)` and `[j, j + L)` with `i < j` carry the same colours, and repeating `c` on
`[i, j)` with period `N = j − i` is again a proper colouring. This is a 3-colouring of `Cay(ℤ/N, S)`, and Theorem W for
the finite group `ℤ/N`, where `m` is the plain average, gives `α ∈ (1/N)ℤ/ℤ`. For infinite `S`, apply this to finite
`S_m ↑ S` and take a limit point of the `α_m` in the circle; the conditions `‖sα‖ ≥ 1/3` are closed. So the answer to
Question 3 below needs only Theorem W for finite cyclic groups. This proof is formalised in Lean 4
(`lean/Recurrence.lean`, theorems `GKR.question3` and `GKR.chromatic_recurrence`), on top of the formal proof of
Theorem W for finite groups (`lean/TheoremW.lean`); both use only Lean's three standard axioms.

The radius `1/3` cannot be improved. For irrational `α`, the colouring `n ↦ ⌊3{nα}⌋` of `ℕ` has
`⋃ᵢ (Aᵢ − Aᵢ) ∩ ℕ = U := {n ≥ 1 : ‖nα‖ < 1/3}`, and `U` contains no set `{n ≥ 1 : ‖nβ‖ < δ}` with `δ > 1/3`. Indeed,
let `H` be the closure of the multiples of `(α, β)` in `𝕋²`; the multiples with `n ≥ 1` are dense in `H`, so it
suffices to find `(x, y) ∈ H` with `‖x‖ > 1/3` and `‖y‖ < δ`. If `H = 𝕋²`, take `(1/2, 0)`. Otherwise
`H = {(x, y) : ax + by = 0}` for a generator `(a, b)` of the annihilator of `H` in `ℤ²`, and `b ≠ 0` because `α` is
irrational. If `|b| ≥ 2`, some `(1/2, y) ∈ H` has `‖y‖ ≤ 1/4`. If `b = ±1`, then `H = {(x, cx)}` for an integer `c`,
and `x = 1/3 + ε` with small `ε > 0` has `‖cx‖ ≤ 1/3 + |c|ε < δ`. (The referee's check of this argument and of the
corollary is recorded in the research log.)

This answers Question 3 of Glasscock, Koutsogiannis and Richter ("On Katznelson's question for skew product
systems", Bull. Amer. Math. Soc. 59 (2022), 569–606; arXiv:2106.11393), who proved the case of two colours (their
Theorem 4.9) and wrote that it was not known for three. As they remark, it follows that the difference set
`A − A` of a set `A` with `A ∪ (A − ℓ₁) ∪ (A − ℓ₂) ⊇ ℕ` is a Bohr₀ set. The same holds in every abelian group `Γ`,
directly from Theorem W: for every 3-colouring `Γ = A₁ ∪ A₂ ∪ A₃` some character `ξ` has
`⋃ᵢ (Aᵢ − Aᵢ) ⊇ {g : ‖ξ(g)‖ < 1/3}`. So if three translates of `A ⊆ ℤ^d` cover `ℤ^d`, then `A − A` contains
`{z : ‖z · λ‖ < 1/3}` for some `λ ∈ 𝕋^d`: the case of three translates of their Question 8. Katznelson's question asks whether every
set of Bohr recurrence is a set of chromatic recurrence for every number of colours. It is open: Griesmer
(arXiv:2108.02190) notes that the answer is not known for any countably infinite abelian group. Alweiss
(arXiv:2511.21680, §3) notes that a counterexample needs at least 3 colours; by the corollary it needs at least 4.
For two colours the statement is easy (a bipartite Cayley graph has a character with `ξ(S) = {1/2}`). We have not
found the case of three colours in the literature we read: the papers above, Host, Kra and Maass ("Variations on
topological recurrence", Monatsh. Math. 179 (2016)), and Liu, Wu, Yang and Zhang (arXiv:2603.05490). A separate
agent made an independent search with the same result. We have not been able to read Katznelson's paper
(Combinatorica 21 (2001), 211–219).

**Periodic colourings.** Let `Γ` be finitely generated and `S` finite. In coordinates `Γ ≅ ℤ^r × F` the characters
with `ξ(S) ⊆ [k/n, (k+1)/n]` form a finite union of polytopes with rational vertices in the torus, so if there is
one there is one of finite order `N`, and the colouring `c = 2⌊nξ⌋ mod n` of the proof below is invariant under the
subgroup `ker ξ ⊇ NΓ`, of finite index. So:

> **Corollary (periodic colourings).** If a Cayley graph of a finitely generated abelian group with a finite
> connection set (for example a Cayley graph of `ℤ^d`) has a homomorphism to `C_{2k+1}`, in particular if it is
> 3-colourable, then it has one that is invariant under a subgroup of finite index; and whether it has one is
> decidable (enumerate the vertices of the polytopes).

For `d = 1` and finite `S` (distance graphs) there are periodic optimal colourings for any number of colours, by a
pigeonhole argument on windows (not every optimal colouring is periodic). For `d ≥ 2` that argument fails, and subshifts of finite type on `ℤ²` can be aperiodic. Abrishami,
Esperet, Giocanti, Hamann, Knappe and Möller (arXiv:2411.01951, Problem 4.6) ask whether every locally finite graph with a periodic
proper colouring has a periodic proper colouring with `χ(G)` colours. They show that the answer is yes for graphs of
bounded pathwidth (for instance Cayley graphs of 2-ended groups) and no for some ∞-ended graphs. Cayley graphs of
`ℤ^d`, `d ≥ 2`, are 1-ended; the corollary gives the answer yes for Cayley graphs of finitely generated abelian
groups with finite connection sets and `χ ≤ 3` (for `χ = 2`, the parity of `Cay(Γ′, S)` is a character of order 2
of `Γ′ = ⟨S⟩`; extended to a character of finite order of `Γ`, it gives a periodic 2-colouring). Without a finite
connection set this fails: for irrational `α` and `S = {n : ‖nα‖ ≥ 1/3}`, `Cay(ℤ, S)` is 3-colourable, but `S` meets
every `Mℤ`, so it has no periodic proper colouring at all. Vallentin,
Weißbach and Zimmermann (arXiv:2407.03513) note that it is not known whether the chromatic number of a lattice (the
Cayley graph of `Λ ≅ ℤⁿ` on its strict Voronoi vectors) is computable, and ask whether there is always a periodic
colouring with `χ(Λ)` colours. By the corollary, whether `χ(Λ) ≤ 3` is decidable (given the strict Voronoi vectors in
coordinates of `Λ`), and if it is, a periodic
3-colouring exists.

**Many colours: periodicity fails.** Three colours are special: with many colours, colourings can encode tilings,
and tilings can be aperiodic. For a finite `F ⊂ Γ` let
`S_F = (F − F) ∖ {0}`. Every translate `x + F` is a clique of `Cay(Γ, S_F)`, so `χ ≥ |F|`, and the proper
`|F|`-colourings are the tilings by `F`: if `c` is one, each translate `x + F` meets each colour class `A` exactly
once, which says `A ⊕ (−F) = Γ`; conversely, if `F ⊕ A = Γ`, colouring `f + a` by `f` is proper. A colouring
invariant under a subgroup of finite index gives a periodic tiling, and conversely. Greenfeld and Tao (Ann. of
Math. 200 (2024), arXiv:2211.15847) found a finite `F` that tiles `ℤ² × G₀` (`G₀` a finite abelian group), and one
that tiles `ℤ^d` for some large `d`, but neither tiles periodically. So:

> **Proposition (aperiodic colourings).** There are a finitely generated abelian group `Γ` (`ℤ^d` for some `d`, or
> `ℤ² × G₀` with `G₀` finite) and a finite `S ⊂ Γ` such that `Cay(Γ, S)` is `k`-colourable, `k = χ(Cay(Γ, S))`,
> but no proper `k`-colouring is invariant under a subgroup of finite index.

For `ℤ^d` we may take `S` generating (replace `ℤ^d` by the subgroup generated by `F − F`, which is free abelian;
a periodic tiling of it would extend to one of `ℤ^d`). Every automorphism of a connected Cayley graph of `ℤ^d` with
finite connection set is affine (Ryabchenko; Morris, Morris and Verret, New York J. Math. 22 (2016)), so a
colouring whose colour-preserving automorphisms act with finitely many orbits is invariant under a subgroup of
finite index. And `x ↦ x mod N` is a periodic proper colouring once `S ∩ Nℤ^d = ∅`. So Problem 4.6 of Abrishami
et al. has a negative answer for Cayley graphs of `ℤ^d`, while the corollary above gives the answer yes when
`χ ≤ 3`. We do not know the least `k` for which this happens (it is at least 4, by the corollary), nor whether
every `k`-colourable Cayley graph of `ℤ²` has a periodic `k`-colouring: Bhattacharya (Amer. J. Math. 142 (2020))
proved that a single tile of `ℤ²` always tiles periodically, so the examples above need `d ≥ 3` or torsion. Nor do
we know whether `k`-colourability of abelian Cayley graphs is decidable for `k ≥ 4`: Greenfeld and Tao (J. Eur.
Math. Soc., 2025, arXiv:2309.09504) proved that tiling a periodic subset of `ℤ² × G₀` by one tile is undecidable,
but a colouring problem cannot single out a subset. The question of Vallentin, Weißbach and Zimmermann concerns the
Voronoi connection set and is not touched.

## 2. Proof

*If.* Put `n = 2k + 1` and `c(g) = 2⌊n ξ(g)⌋ mod n`, with `ξ(g)` read in `[0, 1)`. Along an edge, `n ξ` moves by an
amount in `[k, k + 1]` modulo `n`, so `⌊n ξ⌋` moves by `k` or `k + 1`, and `c` by `2k ≡ −1` or `2k + 2 ≡ 1` (mod
`n`): `c` is a homomorphism to `C_n`.

*Only if.* Let `c : Γ → ℤ/n` be a homomorphism to `C_n` (adjacent vertices get colours that differ by `±1`). Let
`Γ'` be the subgroup generated by `S`. We work in `Γ'`, where `Cay(Γ', S)` is connected (step 2 uses this); a
character of `Γ'` extends to `Γ` because `ℝ/ℤ` is divisible. Averaging over all of `Γ` would be wrong when
`Cay(Γ, S)` is disconnected, since different components can wind in opposite directions. For `g ∈ Γ'` and `s ∈ S`
let `σ(g, s) ∈ {1, −1}` be the lift of `c(g + s) − c(g)`; then `σ(g + s, −s) = −σ(g, s)`.

1. *Squares do not wind* (Krebs–Sankar, Proposition 3.6). For `s, t ∈ S`, both `σ(g, s) + σ(g + s, t)` and
   `σ(g, t) + σ(g + t, s)` lie in `{−2, 0, 2}` and are congruent to `c(g + s + t) − c(g)` modulo the odd number
   `n ≥ 3`. Two different elements of `{−2, 0, 2}` differ by 2 or 4, so the two sums are equal.
2. *Winding* (Krebs–Sankar, Definition 3.5, Proposition 3.3 and Remark 3.7). Let `Z` be the abelian group generated
   by symbols `[s]`, `s ∈ S`, with the relations `[−s] = −[s]`. It is free on one representative of each pair
   `{s, −s}`, except that `2[s] = 0` when `2s = 0`. Let `R` be the kernel of `Z → Γ'`, `[s] ↦ s`. The parity of the
   number of symbols, `Z → ℤ/2`, is well defined.

   For a closed walk from `g` with steps `s₁, …, s_N` (so `Σ s_j = 0`), the sum
   `Σ_j σ(g + s₁ + … + s_{j−1}, s_j)` is `≡ 0 (mod n)`. By (1) it does not change when two consecutive steps are
   swapped, nor when a step `s` followed by `−s` is deleted. Every walk can be brought to a normal form by such
   moves: sort the steps by swaps, then cancel the pairs `s, −s`. So the sum depends only on the base point and on
   the class `r = Σ_j [s_j] ∈ R`. It does not depend on the base point either: going from `g + s` to `g`, doing the
   walk and coming back gives the same sum, and `Cay(Γ', S)` is connected. Write the sum as `n · w_c(r)`. Then
   `w_c : R → ℤ` is a homomorphism. Since `n · w_c(r)` is a sum of `N` odd numbers and `n` is odd, `w_c(r)` has the
   parity of `r`.
3. *Averaging.* `Γ'` is abelian, hence amenable: there is a translation-invariant mean `m` on the bounded functions
   on `Γ'`. For finite `Γ'` this is the plain average; for finitely generated `Γ'`, a limit along Følner sets. Put
   `f(s) = m(g ↦ σ(g, s)) ∈ [−1, 1]`; then `f(−s) = −f(s)`. For a closed walk that represents `r`, the function
   `g ↦ Σ_j σ(g + s₁ + … + s_{j−1}, s_j)` is the constant `n · w_c(r)` by (2). Applying `m`, by linearity and
   translation invariance, gives `Σ_j f(s_j) = n · w_c(r)`.
4. *The character.* Put `ξ(s) = (f(s)/n + 1)/2 ∈ [k/n, (k + 1)/n]`. Then `ξ(−s) = 1 − ξ(s)`. If `2s = 0`, then
   `f(s) = 0` and `ξ(s) = 1/2`. So `ξ` is well defined modulo 1 on `Z`. For `r ∈ R` represented by a walk of length
   `N`, `Σ_j ξ(s_j) = (w_c(r) + N)/2`, which is an integer by the parity in (2). So `ξ` vanishes on `R` modulo 1.
   It is therefore a character of `Γ' = Z/R` with `ξ(S) ⊆ [k/n, (k + 1)/n]`; extend it to `Γ`. ∎

For `ℤ²` with the standard generators, this is the height-function picture of 3-colourings of the square lattice,
and the averages are the slope. Squares not winding is special to odd cycles, and the analogue for four colours
fails: `K₄ = Cay((ℤ/2)², (ℤ/2)² ∖ {0})` is 4-colourable, but no character of `(ℤ/2)²` maps its three generators into
`[1/4, 3/4]`.

## 3. Consequences

**Groups of small exponent.** If `Γ` has exponent `m`, characters take values in `(1/m)ℤ/ℤ`. For `m = 2` and `m = 4`
the only such value in `[1/3, 2/3]` is `1/2`, and a character with `ξ(S) = {1/2}` makes `⌊2ξ⌋` a proper
2-colouring. So Theorem W gives at once:

- (Payan, 1992) a cube-like graph (a Cayley graph of `(ℤ/2)ⁿ`) is never 3-chromatic;
- (Krebs and Sankar, the exponent-4 case of their Theorem 1.2) a Cayley graph of an abelian group of exponent 4 is
  never 3-chromatic.

Krebs and Sankar's Theorem 1.2 covers all abelian groups: some Cayley graph is 3-chromatic if and only if the
exponent is not 1, 2 or 4. Their proof of the exponent-4 case uses the same winding number: they show that the
discrete fundamental group is torsion, and their Theorem 4.1 says that a non-bipartite graph with torsion
fundamental group needs four colours. Cervantes and Krebs (arXiv:2303.06267) gave another proof of Payan's theorem.

**Distance graphs.** For `D ⊂ ℤ_{>0}` finite, `G(ℤ, D)` is 3-colourable if and only if `κ(D) ≥ 1/3`, where
`κ(D) = max_α min_{d∈D} ‖αd‖`. The set `{α : ‖αd‖ ≥ 1/3 for all d ∈ D}` is a finite union of closed intervals whose
endpoints are among the points `(j + 1/3)/d` and `(j + 2/3)/d`. So it is non-empty exactly when it contains one of
these points, and the condition is a finite check. `data/quadratic_planes/winding/distgraph_test.py` compares it
with 3-colourability of the segment `[0, 12 max D + 59]` (SAT), for all 1 747 sets `D ⊆ [1, 24]` with `|D| = 3` and
`gcd(D) = 1`. They agree, and the 74 sets that need four colours are exactly Zhu's list (`D = {1, 2, 3m}` or
`D = {a, b, a + b}` with `a ≢ b (mod 3)`).

`circ_test.py` checks the theorem on 3 888 random Cayley graphs of `ℤ/m × ℤ/n`. The referee checked it for `C₃`,
`C₅` and `C₇` on all 132 874 symmetric connection sets of 42 groups of order at most 27, and on 24 240 random sets
in 101 groups of order at most 100, with no mismatch.

Liu's survey asks whether `χ_c(G(ℤ, D)) < 1/κ(D)` can happen when `|D| = 3` (Problem 3). By Theorem W, if one of
the two numbers is at most 3, then both equal 2 or both lie in the same interval `(2 + 1/(k + 1), 2 + 1/k]`.

**Planes over number fields.** Let `F` be a number field with `i ∉ F`, `L = F(i)`, `T ⊂ L` the unit vectors
(`x + iy` with `x² + y² = 1`) and `A = ℤ[T]`. The plane `F²` is a disjoint union of translates of `Cay(A, T)`, so:

> `χ(F²) ≤ 3` if and only if some character of `A` maps every unit vector into `[1/3, 2/3]`.

For a finite set `U` of unit vectors, the group `ℤU` is free of rank at most `[L : ℚ]`. This is 4 for a real
quadratic `F`, and the rank is 4 for each certificate below. So the test for `U` is a question on a torus of that
dimension. If no character of `ℤU` maps `U` into `[1/3, 2/3]`, then `Cay(ℤU, U)` is not 3-colourable and, by the
de Bruijn–Erdős theorem, neither is some finite unit-distance graph over `F`. A locally constant 3-colouring at a
place (`notes/local_global.md`) colours a finite level, and Theorem W applied to that level gives such a character
(for instance `ξ(t) = ±1/3` at level 1 above 3).

## 4. New values and bounds: `ℚ(√83)`, `ℚ(√107)`, `ℚ(√143)`, `ℚ(√167)`, `ℚ(√203)`

`notes/local_global.md` (§3) lists `d = 83, 107, 143, 167, 203` as the values `d ≡ 11 (mod 12)` below 210 that the
growth searches could not settle. `ℚ(√83)` is the smallest (`docs/research-log.md`). With `D = 510`, for example,
colouring-guided growth stops at 6 008 points, because a 3-colouring of the 2-ball around the origin extends to all
215 478 candidate points. The obstruction does not lie in a small ball around the origin.

For each of the five fields we take the unit vectors `((a + b√d)/D, (c + e√d)/D)` for one denominator `D`, one per
pair `±u`, and the integer relations among them (an LLL-reduced basis). A character with `ξ(U) ⊆ [1/3, 2/3]` would
give `f ∈ [1/3, 2/3]^U` (`f_u` the lift of `ξ(u)` in `[1/3, 2/3]`) with `⟨r, f⟩ ∈ ℤ` for every basis relation `r`.
`certify_w2.py` rules this out by branch and bound over the integers `⟨r, f⟩`. Each leaf carries a rational Farkas
vector, which shows that its equations have no solution in the box.

`check_w.py` checks a certificate in integer and rational arithmetic only:

- the vectors have length 1 (`a² + d b² + c² + d e² = D²`, `ab + ce = 0`);
- each relation sums to zero;
- every branch covers every integer value in its range;
- every Farkas vector is valid.

Any set of valid relations gives necessary conditions, so infeasibility is a proof whether or not the relations form
a basis. They do form one: the referee checked that they generate all relations.

| `d` | `D` | unit vectors (up to sign) | relations | nodes | leaves | file (`data/quadratic_planes/winding/`) | `χ(ℚ(√d)²)` |
|---|---|---|---|---|---|---|---|
| 83 | 510 | 54 | 50 | 6 886 | 4 664 | `cert_83_510_full.json.gz` | 4 |
| 107 | 1170 | 90 | 86 | 12 785 | 9 618 | `cert_107_1170.json.gz` | 4 |
| 143 | 1740 | 54 | 50 | 3 861 | 2 710 | `cert_143_1740.json.gz` | 4 or 5 |
| 167 | 1560 | 54 | 50 | 8 552 | 6 094 | `cert_167_1560.json.gz` | ≥ 4 |
| 203 | 1530 | 78 | 74 | 1 072 | 732 | `cert_203_1530.json.gz` | 4 |

So `Cay(ℤU, U)` is not 3-colourable and `χ(ℚ(√d)²) ≥ 4` in all five cases.

- **`d = 83, 107, 203`.** Here `d ≡ 3 (mod 8)`, so `χ(ℚ(√d)²) ≤ 4` (Fischer 1990, Theorem 10, as cited in
  `notes/quadratic_planes.md`, §1; also Axenovich et al., Graphs Combin. 30 (2014), Theorem 2.2; the repository's
  2-adic proof applies, as `(d + 1)/4` is odd). The value is 4.
- **`d = 143`.** The place above 11 is ramified, with residue field `𝔽₁₁`, and `i` is not in the completion. It
  reduces the unit circle to that of the anisotropic plane over `𝔽₁₁`, which is 5-colourable
  (`finite_planes.json`). So `4 ≤ χ(ℚ(√143)²) ≤ 5`, as for `d = 47`.
- **`d = 167`.** This is the first admissible field: no place gives a locally constant 4-colouring
  (`notes/local_global.md`, §4). We have not looked for an upper bound below 7.

With `d = 83` and `107`, `χ(ℚ(√d)²)` is known for every squarefree `d < 143` except `d = 47`. For `d = 203` the
result also shows that no level of the place above 7 gives a 3-colouring; `notes/local_global.md` had left this
open beyond level 3.

For `d = 83`, a subset of 28 of the 54 vectors (found by greedy deletion) already suffices. Its certificate has
98 300 nodes and is checked too, but not stored. For `ℚ(√11)` with `D = 30`, a 16-vector subset suffices
(`cert_11_30.json.gz`, 1 115 nodes), in agreement with the known 76-vertex graph.

These lower bounds rest on Theorem W, whose averaging step is not constructive, and on the certificates. No explicit
finite graph is given for the five new fields. *Speculation:* such graphs may be large, which would explain why the
growth searches failed; step 3 says nothing about their size. A quantitative version of step 3 over Følner boxes,
or an explicit finite graph certified like those of the other 27 fields, would remove the dependence on the infinite
averaging.

## 5. Checks and literature

- **Controls.**
  - For 40 pairs `(d, D)` with `d ∈ {2, 3, 5, 6, 7, 15, 19, 31, 39, 43}` (fields known to be 2- or 3-colourable),
    the test is feasible, as it must be.
  - For the five fields, other even denominators with irrational unit vectors of rank 4 give feasible tests:
    `(107, 1450)`, `(143, 1752)`, `(167, 1820)`, `(203, 1590)`. Odd `D` always gives a feasible test, since
    `a + b + c + e` is odd for every unit vector, so `ξ = 1/2` on all of them works.
  - The infeasible denominators found so far are all multiples of 6. The floating-point test (`kapparel.py`) was
    infeasible for other denominators of the five fields as well (`d = 107`: 1590, 1716, 1794; `143`: 1020, 1248,
    1260, 1272, 1440; `167`: 840, 1320, 1392, 1440, 1680, 1740).
- **The checker.** `check_w.py` raises an error, rather than asserting, on every failed check. It accepts only
  integer coordinates, relations and branch values, and exact fractions in the Farkas vectors. So it cannot be fooled
  by `python -O` or by floating-point entries. The referee also checked every certificate with a checker of its own.
- **Scan.** `scan_all.py` runs the test, with certificates, over every squarefree `d ≡ 11 (mod 12)` in a range. The
  runs up to 2000 are in progress.
- **Literature.**
  - Closest: Krebs and Sankar (the same winding number, used for torsion fundamental groups; no characters and no
    averaging).
  - The "if" direction is classical: Zhu's multiplier method in Liu's survey; Alon 2013, Theorem 3.2;
    García-Marco, Knauer and Menara, arXiv:2607.26942, Lemma 2.2.
  - Consistent with Theorem W, and not implying it: Cervantes and Krebs, arXiv:2303.06262, 2303.06267 and
    2303.06272. The last characterises 3-colourability for small dimension and rank by two forbidden subgraphs,
    diamond lanyards and `C₁₃(1, 5)`; `C₁₃(1, 5)` indeed has no suitable character. Also Krebs and Leyva,
    arXiv:2511.03028.
  - On Bohr recurrence: Katznelson, Combinatorica 21 (2001) (not seen by us); Glasscock, Koutsogiannis and Richter,
    Bull. Amer. Math. Soc. 59 (2022) (their Question 3 is answered in §1); Griesmer, arXiv:2108.02190; Host, Kra and
    Maass, Monatsh. Math. 2016; Alweiss, arXiv:2511.21680; Liu, Wu, Yang and Zhang, arXiv:2603.05490. Perarnau and
    Serra, arXiv:2409.20160, §3.5, has `χ_f ≤ χ_c ≤ 1/κ`.
  - On periodic colourings and computability: Abrishami, Esperet, Giocanti, Hamann, Knappe and Möller,
    arXiv:2411.01951; Vallentin, Weißbach and Zimmermann, arXiv:2407.03513 (see §1). On tilings: Greenfeld and Tao,
    Ann. of Math. 200 (2024) and J. Eur. Math. Soc. (2025); Bhattacharya, Amer. J. Math. 142 (2020). On automorphisms
    of Cayley graphs of `ℤ^d`: Morris, Morris and Verret, New York J. Math. 22 (2016), after Ryabchenko. We found no
    earlier statement of the proposition on aperiodic colourings, which is a direct translation of Greenfeld and
    Tao's theorem; the latest version of Abrishami et al. (June 2025) does not mention tilings.
  - We did not find Theorem W in any of these, nor in web searches on characters, circular colourings and odd-cycle
    homomorphisms of Cayley graphs. The proof uses only classical tools, so an earlier occurrence is possible; we
    would be glad to learn of one.
- **Complexity.** By Theorem W, whether a finite abelian Cayley graph, given by `Γ` and `S`, maps to `C_{2k+1}` can
  be decided by checking the `|Γ|` characters, in `O(|Γ| |S|)` operations: polynomial in the size of the graph.

## Reproducing

```sh
cd data/quadratic_planes/winding
python3 check_w.py cert_83_510_full.json.gz     # prints VERIFIED (likewise the other cert_*.json.gz)
python3 kapparel.py 83 510                      # the floating-point test: INFEASIBLE
python3 certify_w2.py full_83_510.json out.json # rebuilds a certificate (about 5 min)
python3 distgraph_test.py                       # distance graphs against Zhu's list
```
