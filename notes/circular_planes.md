# Circular colourings of planes over number fields

*Working note, 4 October 2026. Theorem C is proved: the upper bound by an explicit colouring, the lower bound by exact
certificates, each accepted by two exact checkers that share no code, together with Theorem W⁺
(`notes/winding_lemma.md` §2, proved in Lean in `lean/TheoremWInf.lean`). An internal referee (a separate AI agent,
with its own programs) checked Theorem C for `ℚ(√11)` and `ℚ₇` and found it correct; its corrections to this note are
applied (research log, 4 October). Proposition C2 and Question 1 are new and not yet refereed. Nobody outside the
project has checked any of it.*

## 1. Circular colourings and characters

For integers `p ≥ 2q ≥ 2`, the circular clique `K_{p/q}` has vertices `ℤ/p`, with `i ~ j` when
`q ≤ (j − i mod p) ≤ p − q`, and the circular chromatic number `χ_c(G)` is the infimum of the `p/q` with
`G → K_{p/q}`; `χ(G) − 1 < χ_c(G) ≤ χ(G)`. For the plane, DeVos, Ebrahimi, Ghebleh, Goddyn, Mohar and Naserasr
proved `χ_c(ℝ²) ≥ 4` (SIAM J. Discrete Math. 21 (2007) 461–465); since de Grey's theorem gives a finite subgraph with
`χ = 5`, in fact `χ_c(ℝ²) > 4`. Junosza-Szaniawski (arXiv:1506.01886) gives the upper bound `4 + 4√3/3 ≈ 6.31`. The
trivial values are classical (`χ_c = 2` for a bipartite plane, `χ_c = 3` for a plane with `χ = 3` that contains a
unit triangle, such as `ℚ(√3)²`); we found no earlier non-integral value of `χ_c` for the plane over a number field
or a `p`-adic field (web searches, 4 October 2026; the referee's searches found none either).

For a field `F` (with `i ∉ F`), `F²` is the Cayley graph `Cay(L, T)` of §2 of `notes/three_colours_number_fields.md`.
By Theorem W⁺, for `p/q < 4` a Cayley graph `Cay(Γ, S)` of an abelian group with finite `S = −S` maps to `K_{p/q}` if
and only if some character `ξ` has `ξ(S) ⊆ [q/p, 1 − q/p]`. Since `K_{p/q}` is finite, `F² → K_{p/q}` if and only if
every finite subgraph does (compactness); so for `p/q < 4`, `F² → K_{p/q}` if and only if every finite set of unit
vectors has a character into `[q/p, 1 − q/p]`.

## 2. Upper bounds from residue planes

> **Proposition C1.** Let `v` be a place of `F` whose residue field is `𝔽_p`, `p ≡ 3 (mod 4)` (any ramification), and
> let `λ : 𝔽_{p²} → 𝔽_p` be `𝔽_p`-linear with `λ(μ_{p+1}) ⊆ {m, m + 1, …, p − m}`. Then `F²` maps to `K_{p/m}`, so
> `χ_c(F²) ≤ p/m`.

*Proof.* As `−1` is not a square in `𝔽_p`, `i ∉ F_v`, and `L_w = F_v(i)` is unramified over `F_v` with residue field
`𝔽_{p²}`. A unit vector `z` has `z z̄ = 1`, so `w(z) = 0`; conjugation reduces to the Frobenius `x ↦ x^p` of
`𝔽_{p²}/𝔽_p`, so the residue `ẑ` satisfies `ẑ^{p+1} = 1`. Unit steps preserve the classes of `L_w` modulo `O_w`; fix
a representative `ρ` of each class and colour `z` by `λ((z − ρ) mod 𝔪_w) ∈ ℤ/p`. A unit step changes the colour by an
element of `{m, …, p − m}`: a homomorphism `F_v² → K_{p/m}`, which restricts to `F²`. ∎

The best `m` (`residue_circular.py`, and `kappa1` below): `p = 3`: `m = 1` (`K₃`, three colours); `p = 7`: `m = 2`,
with `λ(a + bi) = 2a + 3b` taking the values 2, 3, 4, 5 on `μ₈ = {±1, ±i, ±2 ± 2i}` (so `K_{7/2}`); `p = 11`: `K_{11/2}`;
`p = 19`: `K_{19/4}`; `p = 23`: `K_{23/3}`; `p = 31`: `K_{31/4}`; `p = 43`: `K_{43/5}`. For residue degree 1, only
`p = 3` and `p = 7` give a bound below 4 among all primes `p ≡ 3 (mod 4)` below 400.

So every number field with a place of residue field `𝔽₇` has `χ_c(F²) ≤ 7/2`, hence `χ(F²) ≤ 4` (as `K_{7/2} → K₄`;
for real quadratic fields with `d ≡ 0, 1, 2, 4 (mod 7)` this is Moorhouse's Corollary 8.3). Examples: `ℚ(√11)` and
`ℚ(√23)` (7 splits), `ℚ(√35)` (7 ramifies), `ℚ(√2, √7)`, `ℚ(2cos(2π/7), √7)`, and `ℚ₇` itself.

## 3. An exact value

> **Theorem C.** `χ_c(ℚ(√11)²) = χ_c(ℚ(√35)²) = χ_c(ℚ₇²) = 7/2`, while the three chromatic numbers are 4.

*Proof.* The upper bounds are Proposition C1 with `p = 7` (`11 ≡ 2² (mod 7)`, so `√11 ∈ ℚ₇`; `7 | 35`). For the lower
bound let `U` be the 70 unit vectors `G_25 V` (one per pair `±u`) with `V = {1, u_1, ū_1, u_7, ū_7, u_19, ū_19}`,
`u_n = (n + i√d)²/(n² + d)`, `d ∈ {11, 35}`, and `G_25` the rational unit vectors with denominator dividing 25. The
certificates `data/number_fields/circular/cert_sqrt11_open_7_2_N25.json.gz` and `cert_sqrt35_open_7_2_N25.json.gz`
show that no character of `ℤU` maps every vector of `U` into the open interval `(2/7, 5/7)`: they branch on the
integer values of 66 relations among the vectors (17 744 nodes for `d = 11`, 14 199 for `d = 35`) and close each of
their 12 850, resp. 10 241, leaves by an exact Farkas vector. Two exact checkers accept both: `check_open.py`, and
`check_open_indep.py`, written by the referee from scratch, which also rebuilds `U = G_25 V` from `d` in its own
arithmetic and checks that the relations are nonzero. Both reject corrupted copies (a changed unit, relation, field,
interval, leaf or branch, or an appended zero relation; `tests/test_three_colours.py`). If `F²` mapped to `K_{p/q}`
with `p/q < 7/2`, Theorem W⁺ would give a character of `ℤU` with values in `[q/p, 1 − q/p] ⊂ (2/7, 5/7)`. So
`χ_c(F²) ≥ 7/2`, and `ℚ(√11)² ⊂ ℚ₇²` gives the same for `ℚ₇`. Finally `K_{7/2} → K₄` gives `χ ≤ 4`, and
`χ ≥ χ_c = 7/2 > 3` gives `χ ≥ 4` (also Theorem 1 of `notes/four_colours_11_mod_12.md`, and
`notes/quadratic_planes.md` §5 for `ℚ₇`). ∎

The lower bound is about the infinite graph `Cay(ℤU, U)`: Theorem W⁺ is applied to it directly. Every finite subgraph
has `χ_c ≤ 7/2` and their supremum is `7/2`, but we do not know a finite subgraph with `χ_c = 7/2` exactly. The
16 characters that come from the 7-adic colourings (two embeddings, eight maps `λ`) lie in the closed box `[2/7, 5/7]`
and end at tight leaves of the tree; so the closed-interval statement is false and the open interval is essential.
The certificates were found with a floating-point solver (`kappa_max.py`: SCIP, through OR-tools, maximises
`min_u ‖ξ(u)‖` and returns `2/7` to solver precision; `cert_open.py` builds the exact tree); nothing in the proof
depends on the solver.

## 4. Locally constant colourings give only 2, 3 and 7/2 below 4

Let `K` be a finite extension of `ℚ_p` with `i ∉ K` and `K(i)/K` unramified (for odd `p` always; for `p = 2` exactly
when (a) of Theorem B fails), `𝔽_q` its residue field, `w` the place of `L = K(i)`, and

    κ₁(K) = max over additive λ : 𝔽_{q²} → (1/p)ℤ/ℤ of min over z ∈ μ_{q+1} of ‖λ(z)‖.

> **Proposition C2.** If `K²` has a homomorphism to `K_{P/Q}`, `P/Q < 4`, that is constant on the cosets of some power
> of the maximal ideal of `O_w`, then `Q/P ≤ κ₁(K)`. Conversely `K²` maps to `K_{p/m}` when `m/p = κ₁(K) > 0`.

*Proof.* The converse is Proposition C1 with `𝔽_q` in place of `𝔽_p`. The colouring induces a homomorphism
`Cay(O_w/𝔪^k, T̄) → K_{P/Q}`, so Theorem W⁺ gives a character `ξ` of `O_w/𝔪^k` with `‖ξ(z)‖ ≥ Q/P > 1/4 ≥ 1/(2p)` on
`T(K)`. Suppose `ξ(𝔪) ≠ 0`, and let `J ≥ 1` be largest with `ξ(𝔪^J) ≠ 0`; `ξ` induces a nonzero additive
`λ_J : 𝔪^J/𝔪^{J+1} = 𝔽_{q²} → (1/p)ℤ/ℤ`. For `x ∈ O_w`, `t = (1 + π^J x)/σ(1 + π^J x) ∈ T(K)` (`π` a uniformiser of
`K`, `σ` the conjugation) has `t ≡ 1 + π^J(x − σ(x)) (mod 𝔪^{J+1})`, and the residues of `x − σ(x)` fill the
`𝔽_q`-line `D = {y − y^q}`. So `ξ(zt) = ξ(z) + λ_J(ẑd)`, `d ∈ D`; the set `λ_J(ẑD)` is a subgroup of `(1/p)ℤ/ℤ`, and
if it were the whole group, `ξ(T(K))` would contain a point with `‖·‖ ≤ 1/(2p)`. So `λ_J` vanishes on `ẑD` for all
`ẑ ∈ μ_{q+1}`, hence on `D·span_{𝔽_q}(μ_{q+1}) = 𝔽_{q²}`: a contradiction. So `ξ` factors through `𝔽_{q²}` and
`Q/P ≤ κ₁(K)`. ∎

Values (`data/number_fields/circular/kappa1.py`): `κ₁ = 0` for `p = 2` (the sum of `λ` over `μ_{q+1}`, `q + 1` odd,
is `λ(0) = 0`); for odd `p`, `κ₁` depends only on `(p, f)`, and `κ₁ = 1/3` for `(3, 1)`, `2/7` for `(7, 1)`, `< 1/4`
for `f = 1` and every prime `p ≡ 3 (mod 4)` with `11 ≤ p < 400`, and `κ₁ = 0` (every `λ` vanishes somewhere on
`μ_{q+1}`) for `(p, f) = (3, 3), (3, 5), (7, 3), (11, 3)`. Checks on finite models (`level2.py`): at level 2 the best
character of `O/49` for `ℚ₇` is the level-1 one (`2/7`; the best that is nonzero on `7O` gives `2/49 < 1/14`), and
for `ℚ₂₇` no character of `O/9` keeps `T` away from 0.

So a colouring that is locally constant at one place gives, below 4, only the bounds 2 (Theorem B, (a)), 3 (Theorem
B, (b)) and `7/2` (a place above 7 with residue degree 1), at least for the residue fields computed.

## 5. Questions

1. *Circular local–global.* Is `χ_c(F²) ∈ {2, 3, 7/2}` whenever `χ_c(F²) < 4`? More precisely: is `χ_c(F²) = 7/2`
   exactly when (a) and (b) of Theorem B fail and some place above 7 has residue degree 1, and `χ_c(F²) ≥ 4` when (a)
   and (b) fail and no such place exists? Theorem B (with Corollary B6) is the part of this at 2 and 3. It predicts
   `χ_c(ℚ(√23)²) = 7/2` and `χ_c(ℚ(√59)²) = 4` (7 is inert in `ℚ(√59)`, and `χ(ℚ(√59)²) = 4`).
   Along the rational rotations with denominator `5^k` there are characters with `min ‖ξ(γ)‖ = 3/10` for every `k`
   we tried (`chartypes.py`, `lift.py`, up to `k = 6`); they look 5-adic, and 5-adic characters cannot survive on a
   whole field plane (every place above 5 splits in `F(i)`), but a proof needs an analogue of Proposition 1 for the
   intervals `[r, 1 − r]` with `1/4 < r < 1/3`.
2. What is `χ_c(F²)` for fields with `χ ≥ 4` and no place of residue field `𝔽₇`, such as `ℚ(√47)` (where
   `4 ≤ χ ≤ 5`)? Question 1 predicts `χ_c ≥ 4`; a value below 4 would give `χ(ℚ(√47)²) = 4`.

Measurements (floating-point MIP, not proofs): for `U = G_25{1, u_n, ū_n : n = 1, 7, 19}`, `max_ξ min_u ‖ξ(u)‖` is
`2/7` for `d = 35`; for `d = 59` the solver stopped with the optimum between 0.264 and 0.288; for `d = 71` (7 splits)
it reported 0.2958, so this `U` alone does not give `7/2` there; for `d = 23` and `d = 47` it did not finish.
