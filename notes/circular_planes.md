# Circular colourings of planes over number fields

*Working note, 4 October 2026. Theorem C is proved: the upper bound by an explicit colouring, the lower bound by exact
certificates, each accepted by two exact checkers that share no code, together with Theorem W⁺
(`notes/winding_lemma.md` §2, proved in Lean in `lean/TheoremWInf.lean`). An internal referee (a separate AI agent,
with its own programs) checked Theorem C for `ℚ(√11)` and `ℚ₇` and found it correct; its corrections to this note are
applied (research log, 4 October). Propositions C2, C3 and Corollaries C4, C5 were checked in the second referee
round of `papers/three-colours/`. Theorem D (§5) was found by a separate agent and refereed by another, which found no
error. Question 1 is open. Nobody outside the project has checked any of it.*

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

> **Corollary C5.** `χ_c(K²) = 7/2` for every finite extension `K` of `ℚ₇` with residue degree 1, and `χ_c(F²) = 7/2`
> for every number field `F` that contains `√11` or `√35` and has a place above 7 of residue degree 1 (for example
> `ℚ(√7, √11)`, `ℚ(√2, √35)`, or `ℚ(√11, √m)` for every `m` that is a nonzero square modulo 7).

*Proof.* `K ⊇ ℚ₇ ⊇ ℚ(√11)`, and `F ⊇ ℚ(√11)` or `ℚ(√35)`, give `χ_c ≥ 7/2` by Theorem C; the residue field `𝔽₇` gives
`χ_c ≤ 7/2` by Proposition C1. ∎

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

> **Proposition C3 (residue fields).** `κ₁(K) > 1/4` only for `(p, f) = (3, 1)`, where `κ₁ = 1/3`, and
> `(p, f) = (7, 1)`, where `κ₁ = 2/7`. In every other case (with `i ∉ K` and `K(i)/K` unramified) `κ₁(K) < 1/4`, and
> `κ₁(K) = 0` when `p = 2` or `f ≥ 3`.

*Proof.* `p = 2`: the elements of `μ_{q+1}` add up to 0 and `q + 1` is odd, so `λ ≡ 1/2` on `μ_{q+1}` is impossible.
For odd `p` (so `p ≡ 3 (mod 4)`, `f` odd) let `S(b) = Σ_{z ∈ μ_{q+1}} ψ(Tr(bz))`, `ψ(x) = e^{2πix/p}`. Writing the
indicator of `μ_{q+1} = ker N` with the characters `χ` of `𝔽_q^×` and using Hasse–Davenport
(`G(χ∘N, ψ∘Tr) = −G(χ, ψ_q)²`) gives `S(b) = −Σ_{y ∈ 𝔽_q^×} ψ_q(y + N(b)/y)`, minus a Kloosterman sum; so
`|S(b)| ≤ 2√q` (Weil). For `b ≠ 0` the number of zeros of `Tr(bz)` on `μ_{q+1}` is at least
`(q + 1)/p − 2(p − 1)√q/p > 0` when `f ≥ 3`. For `f = 1`, `p = 4k + 3`, the number of `z` with `‖Tr(bz)/p‖ < 1/4` is
at least `(p² − 1)/(2p) − 2√p (ln p + 1)` (Dirichlet kernel, `Σ_a 1/|sin(πa/p)| ≤ p(ln p + 1)`), positive for
`p ≥ 1001`; for `11 ≤ p < 3000` a direct computation gives `κ₁ < 1/4` (largest `4/19`, `kappa1_f1.py`). ∎

The bound `|S(b)| ≤ 2√q` and the identity with Kloosterman sums are checked numerically in `torus_sums.py`
(`(p, f) = (3, 1), (7, 1), (11, 1), (19, 1), (3, 3), (7, 3), (3, 5)` and `p = 2`, `f = 1, 3, 5`). `kappa1.py` computes
`κ₁(p, f)` for any small `(p, f)` (for example `κ₁ = 0` at `(3, 3), (3, 5), (7, 3), (11, 3)`), and `level2.py` checks the
level lemma on finite models: at level 2 the best character of `O/49` for `ℚ₇` is the level-1 one (`2/7`; the best
that is nonzero on `7O` gives `2/49 < 1/14`), and for `ℚ₂₇` no character of `O/9` keeps `T` away from 0.

> **Corollary C4 (local circular values).** Let `K` be a finite extension of `ℚ_p`. If `K²` has a locally constant
> homomorphism to some `K_{P/Q}` with `P/Q < 4`, then `p = 2` and `K(i)/K` is ramified (and `K²` is bipartite), or
> `(p, f) = (3, 1)` (and `P/Q ≥ 3`), or `(p, f) = (7, 1)` (and `P/Q ≥ 7/2`).

So the restriction to `F²` of a locally constant colouring of one completion `F_v²` gives, below 4, only the bounds 2
(Theorem B, (a)), 3 (Theorem B, (b)) and `7/2` (a place above 7 with residue degree 1).

## 5. A gap above 3

*Found by a separate agent (programs in `data/number_fields/circular/probe/`), refereed by another with its own exact
programs (`probe/indep/`): no error in Theorem D or Proposition D1; its corrections are applied below.*

> **Theorem D.** Let `F` be a number field with `χ(F²) ≥ 4`. Then (a) `χ_c(F²) ≥ 56/17 ≈ 3.294`, and (b)
> (computer-assisted: exact computations by two independent programs) `χ_c(F²) ≥ 10⁷/3001611 > 3.3315`. Hence for
> every number field `χ_c(F²) ∈ {2, 3}` or `χ_c(F²) ≥ 3.3315`.

*Proof.* If `i ∈ F` there is nothing to prove (`χ = ∞`, no homomorphism to any `K_{p/q}`). Let `F²` map to `K_{p/q}`
with `p/q < 56/17` (for (b): `p/q < 10⁷/3001611`); we show `χ(F²) ≤ 3`. If `p/q ≤ 3`, then
`χ(F²) ≤ χ(K_{p/q}) ≤ 3`. Otherwise `r = q/p ∈ (17/56, 1/3)` (for (b): `r > r* = 0.3001611`). Run the proof of
Theorem B with the interval `[r, 1 − r]` in place of `[1/3, 2/3]`. For a finite `V` and `U = G_N V`, Theorem W⁺ gives a character with `ξ(U) ⊆ [r, 1 − r]`,
so `φ(v) ∈ S_N^(r)` for `v ∈ V`, where `S_N^(r)` is `S_N` with `[r, 1 − r]`. Proposition D1 (for (b): its
computer-assisted form, with `k ≥ 17`) says that `φ(v) = Nε_v + x_v` with `ε_v ∈ E` (the same types c and q) and
`|x_v| ≤ R = max(s√10/3, √2 ε)`, and `R < 0.22` (`s < 0.2`, `ε < 1/30`). The exact-relations lemma (Lemma 1 of
`papers/three-colours/`, = Lemma 5 of `notes/four_colours_11_mod_12.md`) needs only this uniform bound: choose
`N = 5^k > 6R·max_v(|D_v| + Σ_j |μ_{v,j}|)` in Lemma B1 (and `k ≥ 17` for (b)). Nothing else in Lemmas B1–B5 depends
on the interval: case (c) gives (a), case (q) gives (b) of Theorem B, and `χ(F²) ≤ 3`. Finally, if
`χ_c(F²) < 56/17`, the definition of the infimum gives a homomorphism to some `K_{p/q}` with `p/q < 56/17`. ∎

Notation for the probe: `s = 1/2 − r`, `ε = 1/3 − r`, `h = (1 + i)/2`, `Sq_s = [−s, s]²`, `c*_N = Nh`,
`Q_N = {(N/3)(a + bi) : a, b ∈ {1, 2}}`, `ρ = (3 + 4i)/5`, `Λ± = ((2 ± i)/5)ℤ[i]`,
`P_k^(s) = {x : conj(x)ρ^j ∈ Sq_s, |j| ≤ k}`. Then `c ∈ S_N^(r)` iff `conj(c)ρ^j ∈ h + Sq_s + ℤ[i]` for `|j| ≤ k`. For
`q ∈ Q_N` let `m_j(q) ∈ ℤ[i]` be the lattice point with `conj(q)ρ^j ∈ h + Sq_s + m_j(q)` (the *cell* of `q` at `ρ^j`), and
`Y_k(q) = {y : conj(q + y)ρ^j ∈ h + Sq_s + m_j(q) for all |j| ≤ k}`.

> **Proposition D1 (the probe for `17/56 < r ≤ 1/3`).** For every `k ≥ 1`, `N = 5^k`,
> `S_N^(r) = (c*_N + P_k^(s) + Nℤ[i]) ⊔ ⋃_{q ∈ Q_N} (q + Y_k(q) + Nℤ[i])`, with `|x| ≤ s√10/3` on `P_k^(s)` and
> `|y| ≤ √2 ε` on `Y_k(q)`.

*Proof.*
- *Lemma 1 (values at type points).* `conj(c*_N)ρ^j ∈ h + ℤ[i]`; for `q = (N/3)(a + bi)`,
  `conj(q)ρ^j ≡ (1/3)(−1)^k(−i)^j(a − bi) (mod ℤ[i])` (mod 3, `5 ≡ −1` and `(2 ± i)² ≡ ±i`); in particular all values
  are in `{1/3, 2/3}²`, `p_{±(k−2)} ≡ −p_{±k}`, and for `q′ ∈ Q_{N/5}` the unique `q ∈ Q_N` with `q ≡ q′ (mod N/5)`
  (`a ≡ 2a′ mod 3`) has the same values as `q′` for `|j| ≤ k − 1`.
- *Lemma 2 (radii).* `P_1^(s) = s·P_1^(1)`, the 12-gon with vertices `(±1, ±1/3)`, `(±1/3, ±1)`, `(±5/7, ±5/7)` (times
  `s`), so `|x| ≤ s√10/3`. For `ε ≤ 1/30` (in fact also at `ε = 1/29`, `1/25`), `Y_1(q) = εZ_1(q)` is a hexagon of
  radius `√2 ε` (for `q = (5/3)(1 + i)`: vertices `(−5/4, 0)`, `(−5/7, −5/7)`, `(0, −5/4)`, `(1, −1/2)`, `(1, 1)`,
  `(−1/2, 1)`, times `ε`), and `Y_k(q) ⊆ Y_1(q)`.
- *Lemma C (`k ≥ 2`, `s < 11/56`).* If `x ∈ P_{k−1}^(s)` and `ν ∈ Λ+ \ ℤ[i]`, then `conj(x)ρ^k + ν ∉ Sq_s + ℤ[i]`;
  likewise with `ρ^{−k}`, `Λ−`. Put `u = conj(x)ρ^k`, `w = conj(x)ρ^{k−2} = u(−7 − 24i)/25 ∈ Sq_s`. As
  `|u| ≤ s√10/3` and `s(1 + √10/3) < 3/5`, only `ν ≡ ±(2/5, 1/5)` or `±(−1/5, 2/5)` (with no further shift) can occur.
  For `ν = σ(2/5, 1/5)` and `τ = −σ`: `τu₁ ≥ 2/5 − s`, `τu₂ ≥ 1/5 − s`, so `τ(24u₁ + 7u₂)/25 ≥ (11 − 31s)/25 > s`,
  contradicting `w ∈ Sq_s`; for `ν = σ(−1/5, 2/5)` the same bound holds for `τ Re w`. (`(11 − 31s)/25 > s` iff
  `s < 11/56`.)
- *Lemma Q (`k ≥ 2`, `ε < 5/168`).* If `q ∈ Q_N`, `|y| < 1/6`, the values of `q + y` at `ρ^{±(k−2)}` stay in the
  cells of `q`, and `ν ∈ Λ+ \ ℤ[i]`, then `conj(q + y)ρ^k + ν ∉ h + Sq_s + ℤ[i]` (likewise with `ρ^{−k}`, `Λ−`). By
  `y ↦ iy`, which preserves `S_N^(r)`, `Q_N` mod `N` and `Λ±`, reduce to `p_k = h + (1, 1)/6`, so
  `p_{k−2} = h − (1, 1)/6`. Put `u = conj(y)ρ^k`, so `conj(y)ρ^{k−2} = u(−7 − 24i)/25`. The condition at `ρ^k` reads
  `1/6 + ν_e + u_e ∈ [−s, s] + ℤ` (`e = 1, 2`). As `|u_e| ≤ |y| < 1/6 < 11/30 − s`, a coordinate `ν_e ≡ 1/5` or
  `2/5` is impossible, so of the four non-zero classes `±(2/5, 1/5)`, `±(−1/5, 2/5)` of `Λ+/ℤ[i]` only
  `ν ≡ (−2/5, −1/5)` remains, and it needs `u₁ ≥ 1/15 − ε`. The cell conditions at `ρ^{k−2}` read
  `24u₁ + 7u₂ ≤ 25ε` and `7u₁ − 24u₂ ≤ 25ε`, which give `625u₁ ≤ 775ε`, i.e. `u₁ ≤ 31ε/25`. So
  `1/15 − ε ≤ 31ε/25`, i.e. `ε ≥ 5/168`. (With `ρ^{−k}` and `Λ−` only `ν ≡ (−1/5, −2/5)` remains, and the same
  computation bounds `u₂`.) The hypothesis `|y| < 1/6` cannot be weakened to `|y|_∞ ≤ 1/6`: for `q = (25/3)(1 + i)`,
  `k = 2`, `y = (1/6, 1/6)` and `ν = (−2/5, −1/5)` the condition holds at `ε = 0`. In the induction `|y| ≤ √2 ε < 0.043`.
- *Base case `k = 1`.* Write `c = c*_5 + x + λ`, `x ∈ Sq_s`, `λ ∈ ℤ[i]/5` (25 classes), `μ± = conj(λ)ρ^{±1} ∈ Λ±`. For
  `s < 11/56`, `|x| ≤ s√2`, so a shift `ν` can only be a minimal representative (`s(1 + √2) < 3/5`). Both `μ±`
  integral: `x ∈ P_1^(s)`. Exactly one integral (8 classes): this is Lemma C for `k = 1` (`u = conj(x)ρ^{±1}`,
  `w = conj(x)ρ^{∓1} ∈ Sq_s`, with `|u| ≤ s√2` in place of `s√10/3`), so `s ≥ 11/56`. Both non-integral (16 classes, 4 orbits under `x ↦ ix`): with the minimal representatives, three
  orbits are infeasible for `s < 1/4`, `2/7`, `3/8` (three-term certificates in `probe/base_caseC.txt`; with all
  representatives the third threshold is `1/4`), and the fourth orbit is the component of a point of `Q_5`
  (`λ = 1 + i` gives `q = (10/3)(1 + i)`). `probe/base_case.txt` lists the 25 classes at `r = 17/56`.
- *Induction.* For `k ≥ 2`, `c ∈ S_N^(r) ⊆ S_{N/5}^(r)`. If `c = c*_N + x + (N/5)λ` with `x ∈ P_{k−1}^(s)`, the
  condition at `ρ^{±k}` and Lemma C give `(2 ± i) | λ`, so `c ≡ c*_N + x (mod N)` and `x ∈ P_k^(s)`. If
  `c = q′ + y + (N/5)λ` with `q′ ∈ Q_{N/5}`, `y ∈ Y_{k−1}(q′)`, take `q ∈ Q_N` with `q ≡ q′ (mod N/5)` (Lemma 1: the
  same values, so the same cells, for `|j| ≤ k − 1`); Lemma Q gives `c ≡ q + y (mod N)`, and `y ∈ Y_k(q)`. The
  converse inclusion is Lemma 1. ∎

*Sharpness, extensions and limits.*
- `17/56` is best possible for a statement valid for every `k`: at `r = 17/56`, `S_N^(r)` has further isolated points
  for `N = 5, 25, 125, 625` (8, 20, 12 and 4 of them; coordinates with denominators 8 and 56; `probe/first_new.txt`).
  The 8 points at `N = 5` and the family `(375/56)(1 + i)` at `N = 25` have values with denominator 56; the other 16
  points at `N = 25` have values with denominators 40, 280, 1400 (minimum 17/56).
- The constant comes from the small levels, and Theorem D only needs large `k`. Both programs find that
  `S_{5^5}^(r)` is the 5 main components plus 8 isolated points at `r = 333/1106`, so only the main components for
  `r > 333/1106`; the thresholds of Lemmas C and Q using all constraints `|j| ≤ k − 1` are non-increasing in `k`
  (in the frame `conj(x)ρ^k` the level-`(k+1)` system contains the level-`k` one) and equal `3303/10981` and
  `13183/43924` at `k = 6`. So D1 holds for all `k ≥ 5` when `r > 333/1106`, and Theorem D(a) holds with
  `1106/333 ≈ 3.3213`. At `r = 0.3001611` both programs find only the main components at `k = 17, 18, 19` (and extra
  components at `r = 0.3001609` for these `k`), and the lemma thresholds are at most `0.30005568` and `0.30000928` for `k ≥ 12`; this
  is Theorem D(b). (Both programs give the same number of components at every level `k ≤ 19` at `r = 0.3001611`:
  `probe/levels3_r3001611_19.txt` and `probe/indep/out_deep.txt`.)
- For `r ≤ 3/10` a classification with types c and q only fails for every `k`: the 5-adic point
  `c_k = c*_N + ρ^k/5 − 5^{k−1}(2 − i)` lies in `S_N^(3/10)` and has `min‖·‖ = 3/10` exactly. Indeed
  `5^{k−1}(2 + i)ρ^j ∈ ℤ[i]` for `−k ≤ j < k`, while `(2 + i)^{2k+1} ≡ 2 + i (mod 5)` (it is `≡ 0` mod `2 + i` and
  `≡ −1` mod `2 − i`); so the value of `c_k` at `ρ^j` is `h + ρ^{j−k}/5` (modulus `1/5`) for `j < k`, and
  `h + 1/5 − (2 + i)/5 = h − (1 + i)/5` at `j = k`. It lies at distance at least `√5·5^{k−1} − 1/5` from `c*_N + Nℤ[i]`
  and `(7/6)5^{k−1} − 1/5` from `Q_N + Nℤ[i]`, far outside the main components (the referee checked this exactly for `k ≤ 10`,
  `probe/indep/out_padic.txt`). This is a limit of the present method, not of the problem: on a whole field plane the
  places above 5 split in `F(i)`, which should kill such types, but the proof would have to add 5-adic types to D1
  and localise at 5.
- The earlier grid search (grid 0.02) had reported clusters with `min‖ξ(γ)‖ = 3/10`; exact polygons show they are the
  components with `κ = 17/56`, `333/1106`, `166/553`, which the grid underestimated.
- *Conjecture.* Let `r0(5^k)` be the least `r` such that `S_{5^k}^(r′)` consists of the main components for every
  `r′ > r` (`r0(5^k) = 17/56` for `k ≤ 4`, `r0(5^5) = 333/1106`, `r0(5^6) = 3303/10981`,
  `r0(5^17) = 452452525229/1507366479208 ≈ 0.30016093`). Then `r0(5^k) → 3/10`, which would give `χ_c ∉ (3, 10/3)` by
  this method.

## 6. Questions

1. *Circular local–global.* Is `χ_c(F²) ∈ {2, 3, 7/2}` whenever `χ_c(F²) < 4`? More precisely: is `χ_c(F²) = 7/2`
   exactly when (a) and (b) of Theorem B fail and some place above 7 has residue degree 1, and `χ_c(F²) ≥ 4` when (a)
   and (b) fail and no such place exists? Theorem B (with Corollary B6) is the part of this at 2 and 3. For `ℚ(√23)`,
   where 7 splits and `χ_c ≤ 7/2` by Proposition C1, it predicts `χ_c = 7/2`; for `ℚ(√59)` (7 inert, `χ = 4`) it
   predicts `χ_c = 4`.
   Theorem D settles the part below `3.3315`. The present method stops at `10/3`: for every `k` the probe has 5-adic
   points with `min ‖ξ(γ)‖ = 3/10` (§5), which the places above 5 should kill on a whole field plane (they split in
   `F(i)`), but the classification of the probe would have to include them; the gap up to `7/2` needs that or another
   idea.
2. What is `χ_c(F²)` for fields with `χ ≥ 4` and no place of residue field `𝔽₇`, such as `ℚ(√47)` (where
   `4 ≤ χ ≤ 5`)? Question 1 predicts `χ_c ≥ 4`; a value below 4 would give `χ(ℚ(√47)²) = 4`.

Measurements (floating-point MIP, not proofs): for `U = G_25{1, u_n, ū_n : n = 1, 7, 19}`, `max_ξ min_u ‖ξ(u)‖` is
`2/7` for `d = 11`, `35` and `23` (SCIP on an LLL-reduced relation basis, optimal; the exact certificate for `d = 23` is
being built), `0.2638` for `d = 59` (an optimal character with irrational-looking values, so this `U` does not reach
`χ_c ≥ 4`), and `0.2958` for `d = 71` (7 splits; this `U` alone does not give `7/2`). Without the LLL step the solver
had failed on `d = 23`, `47` (no solution in 300 s) and reported a wrong value elsewhere; those failures were
numerical.
