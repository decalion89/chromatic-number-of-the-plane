# Circular colourings of planes over number fields

*Working note, 4 October 2026. Theorem C is proved: the upper bound by an explicit colouring, the lower bound by exact
certificates, each accepted by two exact checkers that share no code, together with Theorem W⁺
(`notes/winding_lemma.md` §2, proved in Lean in `lean/TheoremWInf.lean`). An internal referee (a separate AI agent,
with its own programs) checked Theorem C for `ℚ(√11)` and `ℚ₇` and found it correct; its corrections to this note are
applied (research log, 4 October). Propositions C2, C3 and Corollaries C4, C5 were checked in the second referee
round of `papers/three-colours/`. Theorem D (§5) was found by a separate agent and refereed by another, which found no
error. Theorems E and F (§6), which decide `χ_c(F²)` below 4, were found by a separate agent and refereed by two
others, one for each theorem; neither found an error. Each has one computer-assisted step, checked by three exact
programs written separately. Nobody outside the project has checked any of it.*

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
> let `m ≥ 1` and `λ : 𝔽_{p²} → 𝔽_p` be `𝔽_p`-linear with `λ(μ_{p+1}) ⊆ {m, m + 1, …, p − m}`. Then `F²` maps to `K_{p/m}`, so
> `χ_c(F²) ≤ p/m`.

*Proof.* As `−1` is not a square in `𝔽_p`, `i ∉ F_v`, and `L_w = F_v(i)` is unramified over `F_v` with residue field
`𝔽_{p²}`. A unit vector `z` has `z z̄ = 1`, so `w(z) = 0`; conjugation reduces to the Frobenius `x ↦ x^p` of
`𝔽_{p²}/𝔽_p`, so the residue `ẑ` satisfies `ẑ^{p+1} = 1`. Unit steps preserve the classes of `L_w` modulo `O_w`; fix
a representative `x_A` of each class `A` and colour each `x ∈ A` by `λ((x − x_A) mod 𝔪_w) ∈ ℤ/p`. A unit step changes the colour by an
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

> **Proposition C6.** `χ_c(ℚ(√59)²) ≥ 53/15 > 7/2`, while `χ(ℚ(√59)²) = 4`. So the planes over `ℚ(√11)` and `ℚ(√59)`
> have the same chromatic number 4 and different circular chromatic numbers.

*Proof.* Take the same 70 vectors `U = G_25 V` (`n = 1, 7, 19`) for `d = 59`. The certificate
`cert_sqrt59_open_53_15_N25.json.gz` (127 181 nodes, 105 761 leaves) shows that no character of `ℤU` maps `U` into the
open interval `(15/53, 38/53)`; `check_open.py` and `check_open_indep.py` (which rebuilds `U` from `d`) accept it. By
Theorem W⁺, `ℚ(√59)²` has no homomorphism to `K_{p/q}` with `p/q < 53/15`. Here 7 is inert, so Proposition C1 gives
nothing; `χ ≥ 4` because `59 ≡ 11 (mod 12)`. And `χ ≤ 4` from the place `v` above 2: as `59·3 ≡ 1 (mod 8)`,
`F_v = ℚ₂(√59) = ℚ₂(√3)`, and `L_w = F_v(i) = ℚ₂(√3, √−3)` is unramified over `F_v` with residue field `𝔽₄`. Unit
vectors are units of `O_w` whose residues lie in `μ₃ = 𝔽₄^×`, so colouring each class of `L_w` modulo `O_w` by the
residue modulo `𝔪_w`, as in the proof of Proposition C1, is a proper 4-colouring (Fischer 1990;
`notes/local_colourings.md`). ∎

Theorem F (§6) gives more: `χ_c(ℚ(√59)²) = 4`. Proposition C6 is an explicit form of the separation: 70 explicit unit
vectors already push `χ_c` above `7/2`.

The lower bound is about the infinite graph `Cay(ℤU, U)`: Theorem W⁺ is applied to it directly. Every finite subgraph
has `χ_c ≤ 7/2` and their supremum is `7/2`; some finite subgraph has `χ_c = 7/2` exactly (Corollary F12, §6.8), and
§6.8 gives one with 155 vertices, found by computer and certified. The
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

So, below 4, the best upper bound on `χ_c(F²)` given by the restriction to `F²` of a locally constant colouring of
one completion `F_v²` is 2 (Theorem B, (a)), 3 (Theorem B, (b)) or `7/2` (a place above 7 with residue degree 1).

## 5. A gap above 3

*Found by a separate agent (programs in `data/number_fields/circular/probe/`), refereed by another with its own exact
programs (`probe/indep/`): no error in Theorem D or Proposition D1; its corrections are applied below. A third agent
refereed the version in Section 9 of `papers/three-colours/` (complete proofs) with a third implementation
(`probe/indep2/`): correct; its corrections are applied there. Theorem E (§6) supersedes the bounds of Theorem D;
part (a) keeps its proof by hand.*

> **Theorem D.** Let `F` be a number field with `χ(F²) ≥ 4`. Then (a) `χ_c(F²) ≥ 56/17 ≈ 3.294`, and (b)
> (computer-assisted: exact computations by three programs written separately) `χ_c(F²) ≥ 10⁷/3001611 > 3.3315`.
> Hence for every number field `χ_c(F²) ∈ {2, 3}` or `χ_c(F²) ≥ 56/17` (and, with (b), `χ_c(F²) > 3.3315`).

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
  is Theorem D(b). (The three programs give the same number of components at every level `k ≤ 19` at
  `r = 0.3001611`: `probe/levels3_r3001611_19.txt`, `probe/indep/out_deep.txt`, `probe/indep2/check5_levels_thetastar.txt`;
  and the same thresholds, `probe/lemma_2_14.txt` for `k ≤ 14`, `probe/indep/out_thresholds.txt` and
  `probe/indep2/check3_thresholds.txt` for `k ≤ 22`.)
- For `r ≤ 3/10` a classification with types c and q only fails for every `k`: the 5-adic point
  `c_k = c*_N + ρ^k/5 − 5^{k−1}(2 − i)` lies in `S_N^(3/10)` and has `min‖·‖ = 3/10` exactly. Indeed
  `5^{k−1}(2 + i)ρ^j ∈ ℤ[i]` for `−k ≤ j < k`, while `(2 + i)^{2k+1} ≡ 2 + i (mod 5)` (it is `≡ 0` mod `2 + i` and
  `≡ −1` mod `2 − i`); so the value of `c_k` at `ρ^j` is `h + ρ^{j−k}/5` (modulus `1/5`) for `j < k`, and
  `h + 1/5 − (2 + i)/5 = h − (1 + i)/5` at `j = k`. It lies at distance at least `√5·5^{k−1} − 1/5` from `c*_N + Nℤ[i]`
  and `(7/6)5^{k−1} − 1/5` from `Q_N + Nℤ[i]`, far outside the main components (the referee checked this exactly for `k ≤ 10`,
  `probe/indep/out_padic.txt`). This is a limit of the present method, not of the problem: on a whole field plane the
  places above 5 split in `F(i)`, which should kill such types, but the proof would have to add 5-adic types to D1
  and localise at 5. §6 avoids this: one rotation with another prime removes them.
- The earlier grid search (grid 0.02) had reported clusters with `min‖ξ(γ)‖ = 3/10`; exact polygons show they are the
  components with `κ = 17/56`, `333/1106`, `166/553`, which the grid underestimated.
- *Conjecture.* Let `r0(5^k)` be the least `r` such that `S_{5^k}^(r′)` consists of the main components for every
  `r′ > r` (`r0(5^k) = 17/56` for `k ≤ 4`, `r0(5^5) = 333/1106`, `r0(5^6) = 3303/10981`,
  `r0(5^17) = 452452525229/1507366479208 ≈ 0.30016093`). Then `r0(5^k) → 3/10`, which would give `χ_c ∉ (3, 10/3)` by
  this method.

## 6. Two primes: 7/2, and every value below 4

*Found by a separate agent (programs in `data/number_fields/circular/twoprime/`) and refereed by two others, one
for each theorem, each with two exact methods of its own (`twoprime/indep_E/`, `twoprime/indep_F/`). Neither found
an error; their corrections are applied below.*

> **Theorem E.** Let `F` be a number field. If `F²` has a homomorphism to `K_{p/q}` with `p/q < 7/2`, then
> `χ(F²) ≤ 3`. Hence `χ_c(F²) ≥ 7/2` whenever `χ(F²) ≥ 4`, so `χ_c(F²) ∈ {2, 3} ∪ [7/2, ∞]`; the value `7/2` is
> attained (Theorem C).

> **Theorem F.** If `F²` has a homomorphism to `K_{p/q}` with `p/q < 4`, then (a) or (b) of Theorem B holds, or
> **(7)** some place of `F` above 7 has residue degree 1. Consequently `χ_c(F²) = 2` if (a) holds; `3` if (b) holds
> and (a) fails; `7/2` if (7) holds and (a), (b) fail; and `χ_c(F²) ≥ 4` if (a), (b) and (7) all fail. So
> `χ_c(F²) ∈ {2, 3, 7/2} ∪ [4, ∞]`, and every value below 4 comes from a locally constant colouring of one completion
> (Corollary C4).

Theorem E has one computer-assisted step (Proposition E1), and Theorem F another (Proposition F1, which also gives
Theorem E again); both are exact finite computations, and the rest is by hand. Theorem E supersedes the bounds of Theorem D (whose part (a) keeps a proof by hand), and gives the lower bounds
of Theorem C again from `χ ≥ 4` alone. Theorem F answers affirmatively the question of earlier versions of this note:
`χ_c(F²) ∈ {2, 3, 7/2}` whenever `χ_c(F²) < 4`, with `7/2` exactly when (a) and (b) fail and (7) holds.

*The list from the first sentence of Theorem F.* Under (a), `χ = χ_c = 2`. Under (b) without (a), `χ = 3`
(Theorem B), so `χ_c = 3` (Corollary B6). Under (7) without (a) and (b), `χ ≥ 4` (Theorem B), so `χ_c ≥ 7/2`
(Theorem E), and `χ_c ≤ 7/2` by Proposition C1, as the residue field is `𝔽₇`. If (a), (b) and (7) all fail, no
`K_{p/q}` with `p/q < 4` receives a homomorphism, so `χ_c ≥ 4`. (If `i ∈ F`, all three fail, as every residue field
above 7 then contains `𝔽₄₉`, and `χ_c = χ = ∞`.)

*Examples.* For squarefree `d ≡ 11 (mod 12)`, (a) and (b) fail, so `χ_c(ℚ(√d)²) = 7/2` if `(d/7) ≠ −1` and
`χ_c(ℚ(√d)²) ≥ 4` if `(d/7) = −1`. If moreover `d ≡ 11 (mod 24)`, then `ℚ₂(√d) = ℚ₂(√3)` and the place above 2 gives
`χ ≤ 4` (as in the proof of Proposition C6), so `χ_c = 4`; only for `(d/7) = −1` and `d ≡ 23 (mod 24)` (`d = 47`,
`143`, `167`, `215`, `311`, …) do we not know whether `χ_c = 4` (for `d = 47`,
`4 ≤ χ_c ≤ 19/4`, §6.6). So `χ_c(ℚ(√23)²) = χ_c(ℚ(√71)²) = 7/2`,
`χ_c(ℚ(√59)²) = χ_c(ℚ(√83)²) = 4`, and `χ_c(ℚ(√47)²) ≥ 4`.

### 6.1 The window (1, 1)

Notation of §5 (`r`, `s = 1/2 − r`, `ε = 1/3 − r`, `h`, `Sq_s`, `ρ`, `c*_N`, `Q_N`, `P_k^(s)`, `Y_k(q)`), and
`σ = (3 + 2i)/(3 − 2i) = (5 + 12i)/13`. For `k, m ≥ 0` let

    G(k, m) = { iᵃρʲσˡ : a ∈ ℤ, |j| ≤ k, |l| ≤ m },
    S^(r)(k, m) = { c ∈ ℂ : Re(c̄γ) ∈ [r, 1 − r] + ℤ for every γ ∈ G(k, m) }.

These are rational rotations, so `G(k, m) ⊂ T` when `i ∉ F` (the automorphism of `L/F` restricts to complex
conjugation on `ℚ(i)`), and `S^(r)(k, 0)` is the probe `S_N^(r)` of §5. The *strip index* of `c` at `γ` is
`n_γ(c) = ⌊Re(c̄γ)⌋`; on `S^(r)(k, m)`, `Re(c̄γ) ∈ [n_γ + r, n_γ + 1 − r]`, and a point of `S^(r)` lies in `S^(r′)`
for `r′ < r` with the same strip indices. The *margins* `μ⁺_γ(c) = Re(c̄γ) − n_γ` and `μ⁻_γ(c) = n_γ + 1 − Re(c̄γ)`
are affine in `c` on each set of points with given strip indices. Those points of `S^(r)(k, m)` form a convex polygon,
and two different index vectors give polygons at positive distance (some functional takes values in two strips
separated by a gap `2r`); so the connected components of `S^(r)(k, m)` are the index vectors that occur. As
`65γ ∈ ℤ[i]` for `γ ∈ G(1, 1)`, `S^(r)(1, 1)` is invariant under `65ℤ[i]`; its *type points* are `65h` and
`(65/3)(a + bi)`, `a, b ∈ {1, 2}`. At a point of `S^(r)(1, 1)` and at a type point no value `Re(c̄γ)` is an integer,
so `n_{−γ} = −n_γ − 1` and `n_{iγ} = −n_{−iγ} − 1`: the strip indices at the 18 functionals below determine those at all
36 rotations.

> **Proposition E1 (computer-assisted).** Let `2/7 < r ≤ 1/3` and `c ∈ S^(r)(1, 1)`. Then for some type point `T`
> and some `m ∈ ℤ[i]`, `c` has the strip indices of `T + 65m` at every `γ ∈ G(1, 1)`.

*Proof.* Put `r₀ = 7/25 < 2/7`; then `S^(r)(1, 1) ⊆ S^(r₀)(1, 1)`. As `Re(c̄·iγ) = −Im(c̄γ)` and `[r, 1 − r] + ℤ` is
symmetric, only the 18 functionals `Re(c̄γ)`, `Im(c̄γ)` with `γ = ρʲσˡ`, `|j|, |l| ≤ 1`, matter.
- *Completeness* (`tree_r7_25.txt.gz`, checked by `verify_tree.py`). After a translation by `65ℤ[i]`, `c` lies in a
  cell `[a + r₀, a + 1 − r₀] × [b + r₀, b + 1 − r₀]` with `0 ≤ a, b ≤ 64`, which fixes the strip indices of `Re(c̄)` and
  `Im(c̄)` (those at `γ = 1` and `γ = −i`). For each of the other 16 functionals, in a fixed order, the candidate
  indices are the integers `n` for which `[n + r₀, n + 1 − r₀]` meets the range of the functional on the cell (the
  checker computes them from the four corners). The tree lists every candidate either as a branch to follow or as
  closed; a closed branch carries an exact Farkas certificate: three inequalities valid on that branch (sides of the
  cell, sides of the strips fixed so far) with rational multipliers `λ ≥ 0` whose combination reads `0 ≤ (a negative
  number)`. The open branches end in 13 leaves, so every `c ∈ S^(r₀)(1, 1)` has, up to `65ℤ[i]`, one of 13 explicit
  index vectors (4 225 cells, 15 620 branch nodes, 4 364 closed branches; the checker takes about a second).
- *The 13 vectors* (`prop7_certificates.txt`, checked by `verify_prop7.py`). Five are those of the type points (the
  checker recomputes the indices of `T + 65m`). Each of the other eight comes with three of its margins and rationals
  `λ₁, λ₂, λ₃ ≥ 0` with `Σλ_i = 1` and `Σλ_iμ_i ≡ 2/7` identically (the linear parts cancel). At a point of
  `S^(r)(1, 1)` every margin is at least `r > 2/7`, so no such point has one of these eight vectors. ∎

The eight discarded vectors are those of the 7-torsion points `65(a + bi)/7`,
`(a, b) ∈ {(1,2), (1,5), (2,1), (2,6), (5,1), (5,6), (6,2), (6,5)}`, whose least margin is exactly `2/7`: they are
on `G(1, 1)` the characters of the 7-adic colourings of Proposition C1 (`λ(a + bi) = 2a + 3b` gives `65(1 + 5i)/7`).
So `2/7` is optimal for this window, and for every finite set of rational rotations: their denominators are prime to
7, and if `D` is a common denominator and `d ≡ D⁻¹ (mod 7)`, the point `Dd(s + ti)/7` has the value `λ(γ mod 7)/7`
at each of them, for `λ(a + bi) = sa + tb` as in Proposition C1. (They are not characters of all of `ℚ(i)`: the
point `65(1 + 5i)/7` has margin `5/119` at `(−15 − 8i)/17`.) Cross-checks: the agent's
lifting (`probe2d.py`, in both orders of the primes) and direct enumeration (`indep2d.py`) find 13 components at
`r = 7/25` and 5 at `r = 2857143/10⁷`; the referee's clipping (`indep_E/ref_clip.py`) and enumeration of all
1 734 260 vertices of the line arrangement (`ref_vertex.py`) find 13 at `7/25` and at `2/7` (where the eight extra
ones shrink to the torsion points) and 5 at `2857143/10⁷`, `200001/700000`, `1429/5000`, `3/10` and `1/3`, and an
exact primal–dual computation (`ref_kappa.py`) gives the largest least margin `2/7` on each extra component, `1/3` on
type q and `1/2` on type c. The prime 13 can be replaced by 17 (`σ′ = (15 + 8i)/17`, modulus 85: only the 5 main
components just above `2/7`), but 5 cannot be dropped: with 13 and 17 the window keeps 12 extra components, with
largest least margins `3/10` and `2/5`. The mixed rotations matter: just above `2/7` the "cross"
`{1, ρ^{±1}, σ^{±1}}` alone keeps 180 extra components (228 at `7/25`), some with largest least margin `55/144`;
removing any one of the six rotation classes with `j = ±1` still leaves only the 5 main components, while removing
one with `j = 0` leaves 8 extra.

### 6.2 From the window to the probe

> **Lemma E2.** Let `2/7 < r ≤ 1/3`, `k ≥ 1`, `N = 5^k`, and let `C ∈ ℂ` satisfy `Re(C̄γ) ∈ [r, 1 − r] + ℤ` for every
> `γ ∈ G(k, 1)`. Then `C ∈ c*_N + P_k^(s) + Nℤ[i]`, or `C ∈ q + Y_k(q) + Nℤ[i]` for some `q ∈ Q_N`. In particular
> `C = Ne + x` with `e ∈ E` and `|x| ≤ R = max(s√10/3, √2 ε) < 0.23`.

*Proof.*
- *(a) Values and digits.* For `|j| ≤ k` the conditions at `ρʲ` and `−iρʲ` say that both coordinates of `C̄ρʲ` lie in
  `1/2 + [−s, s] + ℤ`, so `C̄ρʲ = h + g_j + n_j` with `g_j ∈ Sq_s` and `n_j ∈ ℤ[i]` (unique, as `s < 1/2`), and `n_j`
  is made of the strip indices at `ρʲ` and `−iρʲ`. For `−k < j ≤ k` the *digit* `ν_j = g_j − ρg_{j−1}` satisfies,
  as `C̄ρʲ = ρ·C̄ρ^{j−1}`,

      ν_j = (ρ − 1)h + ρn_{j−1} − n_j,

  so it depends only on strip indices. The same definitions apply to a point `c` of `S^(r)(1, 1)`: `g_t(c)`,
  `n_t(c)` for `|t| ≤ 1` and `ν_t(c)` for `t ∈ {0, 1}`, by the same formulas (for a type point `T` take `n_t(T)` from
  the floors, `g_t(T) = T̄ρᵗ − h − n_t(T)`).
- *(b) Type points.* Let `m ∈ ℤ[i]` and `|t| ≤ 1`. For `T = 65h + 65m`, `65ρᵗ` is a Gaussian integer of odd norm,
  so `T̄ρᵗ ≡ h (mod ℤ[i])`, `g_t(T) = 0` and `ν₀(T) = ν₁(T) = 0`. For `T = (65/3)(a + bi) + 65m`,
  `T̄ρᵗ ≡ (a − bi)u_t/3 (mod ℤ[i])` with `u_t ≡ 65ρᵗ ≡ −(−i)ᵗ (mod 3)` (`65 ≡ −1`, `ρ ≡ −i`); both coordinates lie in
  `{1/3, 2/3} + ℤ`, so `g_t(T) = η_t/6` with `η_t ∈ {±1 ± i}` and `η_{t+1} = −iη_t` (as `−ih ≡ h`), and
  `ν_t(T) = η_t/6 − ρη_{t−1}/6 = −η_{t−1}(1 + 3i)/10 ≠ 0`, `ν₁(T) = −iν₀(T)`.
- *(c) Windows.* For `−k < j < k` put `c = ρ^{−j}C`, so `c̄ = ρʲC̄`. For `γ ∈ G(1, 1)`, `ρʲγ ∈ G(k, 1)`, so
  `c ∈ S^(r)(1, 1)`, and by Proposition E1 `c` has the strip indices of some `T + 65m`. The strip indices of `c` at
  `ρᵗ`, `−iρᵗ` (`|t| ≤ 1`) are those of `C` at `ρ^{j+t}`, `−iρ^{j+t}`, so by (a)
  `(ν_j(C), ν_{j+1}(C)) = (ν₀(T + 65m), ν₁(T + 65m))`, which by (b) is `(0, 0)` or `(ν, −iν)` with
  `ν = −η(1 + 3i)/10`, `η ∈ {±1 ± i}`.
- *(d) Pure words.* The pairs `(ν_j, ν_{j+1})`, `−k < j < k`, cover `ν_{−k+1}, …, ν_k`, and consecutive pairs share a
  digit. So either every `ν_j` is 0, or every `ν_j` is non-zero and `ν_{j+1} = −iν_j` throughout.
- *(e) Identification.* If every digit is 0, then `g_j = ρʲg₀`; with `x̄ = g₀`, `x ∈ P_k^(s)` and
  `C̄ρʲ ≡ h + x̄ρʲ ≡ conj(c*_N + x)ρʲ (mod ℤ[i])` for `|j| ≤ k` (Lemma 1 of §5). Multiplying by units,
  `Re(w̄γ) ∈ ℤ` for `w = C − c*_N − x` and every `γ ∈ G(k, 0)`. The `ℤ`-span of `G(k, 0)` is `(1/N)ℤ[i]`: it is a
  `ℤ[i]`-module containing `ρ^{±k} = (2 ± i)^{2k}/N`, and `(2 + i)^{2k}`, `(2 − i)^{2k}` are coprime. With `1/N` and
  `i/N` this gives `w ∈ Nℤ[i]`. If every digit is non-zero, put `η_{j−1} = −10ν_j/(1 + 3i) ∈ {±1 ± i}` for
  `−k < j ≤ k` and `η_k = −iη_{k−1}`. Then `η_j = −iη_{j−1}` and `ν_j = η_j/6 − ρη_{j−1}/6`, so
  `g_j − η_j/6 = ρ(g_{j−1} − η_{j−1}/6)`, and with `ȳ = g₀ − η₀/6`, `g_j = η_j/6 + ȳρʲ` for `|j| ≤ k`. Let
  `q ∈ Q_N` be the point whose value at `ρ⁰` is `h + η₀/6` (Lemma 1 of §5: the four points of `Q_N` give the four
  values); its values `h + η_j(q)/6` also satisfy `η_{j+1}(q) = −iη_j(q)`, so `η_j(q) = η_j` for `|j| ≤ k`,
  `y ∈ Y_k(q)`, and `C̄ρʲ ≡ conj(q + y)ρʲ (mod ℤ[i])`; as before `C − q − y ∈ Nℤ[i]`.

The bounds `|x| ≤ s√10/3` and `|y| ≤ √2 ε` are Lemma 2 of §5 (the inclusion `Y_k(q) ⊆ εZ_1(q)` holds for every
`ε ≥ 0`; `papers/three-colours/`, Lemma `radii`). As `s < 3/14` and `ε < 1/21`, `R < 0.23`. ∎

### 6.3 Proof of Theorem E

Let `F²` map to `K_{p/q}` with `p/q < 7/2`. If `i ∈ F` there is no such map (`χ(F²) = ∞`). If `p/q ≤ 3`, then
`χ(F²) ≤ 3`. Otherwise `r = q/p ∈ (2/7, 1/3)` and `2q ≤ p < 4q`.

*Claim (the conclusion of Lemma B1).* For every finite `V ⊂ T` there is `β ∈ L` with `Tr(βv) ∈ E` for `v ∈ V`.
Enlarge `V` so that it contains a `ℚ(i)`-basis `b₁, …, b_n` of `L` inside `T`, write `D_v v = Σ_j μ_{v,j} b_j`
(`D_v ∈ ℤ ∖ {0}`, `μ_{v,j} ∈ ℤ[i]`), choose `N = 5^k > 6R·max_v(|D_v| + Σ_j |μ_{v,j}|)`, and let `U = G(k, 1)V`, a
finite symmetric set of unit vectors. Theorem W⁺ gives a character `ξ` of `ℤU` with `ξ(U) ⊆ [r, 1 − r] + ℤ`, and,
as in the proof of Lemma B1, a `ℚ(i)`-linear `φ : L → ℂ` with `ξ(γv) = Re(γφ(v))` for `v ∈ V`, `γ ∈ G(k, 1)`. So
`C_v = conj(φ(v))` satisfies the hypothesis of Lemma E2, and `φ(v) = Ne_v + x_v` with `e_v ∈ E` (`E` is closed under
conjugation) and `|x_v| ≤ R`. The exact-relations lemma (Lemma 1 of the paper), whose proof uses only `6E ⊆ ℤ[i]` and
the bound on `|x|`, gives `D_v e_v = Σ_j μ_{v,j} e_{b_j}`; with `β` defined by `Tr(βb_j) = e_{b_j}`,
`Tr(βv) = e_v ∈ E`.

The rest is the proof of Theorem B: Lemma B2 uses the colouring only through this claim, and the rest of the proof
only through (c) or (q). So (a) or (b) holds, and `χ(F²) ≤ 3`. If `χ(F²) ≥ 4` and `χ_c(F²) < 7/2`, the infimum gives
a homomorphism to some `K_{p/q}` with `p/q < 7/2`, which is excluded; with Corollary B6,
`χ_c(F²) ∈ {2, 3} ∪ [7/2, ∞]`. ∎

### 6.4 The places above 7

Let `𝔽₄₉ = ℤ[i]/7` and `μ₈` its elements of norm 1, the images of the rational rotations (`ρ ≡ 2 + 5i` has order 8
and `ρ² ≡ −i`). Put

    A′₇ = { e ∈ 𝔽₄₉ : Re(ē z) ∈ {2, 3, 4, 5} for every z ∈ μ₈ },
    E₇ = { e ∈ (1/7)ℤ[i] : 7e mod 7 ∈ A′₇ },   𝒜₇ = { e/7 : e ∈ ℤ₄₉, e mod 7 ∈ A′₇ } ⊂ ℚ₄₉.

Then:
- (A1) `A′₇ = {e : e ē = −1} = {2+3i, 2+4i, 3+2i, 3+5i, 4+2i, 4+5i, 5+3i, 5+4i}`, a single `μ₈`-orbit of 8
  elements, closed under conjugation, with sum 0; `0 ∉ A′₇`;
- (A2) `A′₇` contains no coset of a non-zero additive subgroup of `𝔽₄₉`;
- (A3) `Re` maps `A′₇` onto `{2, 3, 4, 5}`.

*Proof* (by hand; exact enumeration, `seven_local.py` and the referee's `indep_F/f49_and_digits.py`, agrees). For
`e = a + bi`, `z = x + yi` (`a, b, x, y ∈ 𝔽₇`), `Re(ē z) = ax + by` and `e ē = a² + b²`, which vanishes only for
`e = 0`, as `−1` is not a square modulo 7 (the non-zero squares are 1, 2, 4). `μ₈` is the conic `x² + y² = 1`. For
`e ≠ 0` the points of the line `ax + by = t` are `(t/eē)(a, b) + s(−b, a)`, `s ∈ 𝔽₇`, and such a point lies on the
conic iff `(eē)² s² = eē − t²`; so `t` is a value of `z ↦ Re(ē z)` on `μ₈` iff `eē − t²` is 0 or a square. Hence
`e ∈ A′₇` iff `e ≠ 0` and neither `eē` nor `eē − 1` is 0 or a square (`t = 0, ±1`), i.e. both lie in `{3, 5, 6}`:
`eē = 6 = −1`. The solutions of `a² + b² = −1` have `{a², b²} = {2, 4}`: the list. As the norm is multiplicative,
`A′₇ = e₀μ₈`; it is closed under conjugation, and its sum is `e₀ Σ_{z ∈ μ₈} z = 0` (as `μ₈ = −μ₈`). (A2): a coset of a non-zero additive subgroup contains a line `{e + sd : s ∈ 𝔽₇}`, `d ≠ 0`, on which
`(e + sd)·conj(e + sd) = eē + s(ed̄ + ēd) + s²dd̄` is a polynomial of degree 2 in `s` (`dd̄ ≠ 0`), so at most 2 of
its 7 points lie in `A′₇`. (A3) is read off the list. ∎

For `e ∈ E₇` the character `a ↦ Re(ē ã) mod 1` of `ℤ_(7)[i]`, where `ã ∈ ℤ[i]` is congruent to `a` modulo
`7ℤ_(7)[i]` (it factors through `𝔽₄₉`), takes values in `{2, …, 5}/7` at every rational rotation: these are the characters of the 7-adic colourings of Proposition C1.
`2E_c`, `3E_q` and `7E₇` lie in `ℤ[i]`, so `42E′ ⊆ ℤ[i]` for `E′ = E_c ∪ E_q ∪ E₇`, and `E′` is closed under
conjugation.

*Digit patterns* (by hand; `seven_patterns.py` and `f49_and_digits.py` agree). The values of these characters at
`ρʲ` are `h + ζ_j` with `ζ_j ∈ {±1/14, ±3/14}²` (both coordinates of an element of `A′₇` lie in `{2, …, 5}`), and
their digit words `ζ_j − ρζ_{j−1}` are the 8 shifts of one word `S` of period 8 with `S_{j+2} = −iS_j`, in which
non-zero digits `u(2 + i)/5` (`u ∈ μ₄`) alternate with zeros. Indeed `ρ̃² = −21 + 20i ≡ −i (mod 7)` and `−ih ≡ h`
(mod `ℤ[i]`), so `ζ_{j+2} = −iζ_j`; for `7e ≡ 2 + 3i`, `7ē ≡ 2 + 4i`, `(2 + 4i)ρ̃ ≡ 5 + 4i` and
`(2 + 4i)ρ̃² ≡ 4 + 5i`, so `ζ₀ = (−3 + i)/14`, `ζ₁ = (3 + i)/14`, `ζ₂ = (1 + 3i)/14`, `ζ₁ − ρζ₀ = (2 + i)/5` and
`ζ₂ − ρζ₁ = 0`; the other classes give the shifts, `ζ_t(e′) = ζ_{t+j}(e)` when `7e′ ≡ 7e·conj(ρ̃)^j (mod 7)`. Type c has
the zero word and type q a word of non-zero digits (also of the form `u(2 + i)/5`). So two consecutive digits
determine the family (c, q or 7) and the shift.

### 6.5 The window (2, 1)

Let `N₀ = 5²·13 = 325` and `W = G(2, 1)`. The type points modulo `N₀` are `N₀h` (family c),
`(N₀/3)(a + bi)`, `a, b ∈ {1, 2}` (family q), and the 8 points `N₀(a + bi)/7` with `a + bi ≡ 5e (mod 7)`,
`e ∈ A′₇` (family 7; `5 ≡ 325⁻¹ mod 7`). Their least margins are `1/2`, `1/3` and `2/7`.

> **Proposition F1 (computer-assisted).** Let `1/4 < r ≤ 1/3` and `c ∈ S^(r)(2, 1)`. Then for some type point `T`
> of family c, q or 7 and some `m ∈ ℤ[i]`, `c` has the strip indices of `T + 325m` at every `γ ∈ W`. If `r > 2/7`,
> then `T` has family c or q.

(Strip indices are floors, so they are defined at the 7-points also when `r > 2/7`; the digit formula uses only
them.)

*Proof.* As for Proposition E1, with `r₀ = 249/1000 < 1/4` and the 30 functionals of the 15 rotations `ρʲσˡ`,
`|j| ≤ 2`, `|l| ≤ 1`. The tree `tree_K2M1.txt.gz` over all `325²` cells (550 474 branch nodes, 199 920 closed
branches, 29 leaves) and `cert_K2M1.txt` are checked by `verify_window.py 2 1` (about 45 s): of the 29 index vectors,
5 are those of the c- and q-points, 8 those of the 7-points, each with all margins at least `2/7`, and each of the
other 16 has a certificate `λ ≥ 0`, `Σλ = 1`, `Σλ_iμ_i ≡ 1/4` on three margins, which excludes it for `r > 1/4`.
Each of the 8 vectors of family 7 has a certificate of the same kind, with three positive multipliers and
`Σλ_iμ_i = 2/7` identically (`cert_K2M1_seven.txt`, written
by `seven_certificates_K2M1.py` and checked by `verify_seven_K2M1.py`; added after the referees' reports, on a
suggestion of the referee of Theorem F), which excludes it for `r > 2/7`; these certificates use margins at
`ρ^{−2}σˡ`, outside `G(1, 1)`. ∎

The 16 excluded vectors are 2-adic: their largest least margin `1/4` is attained at 4- and 8-torsion points such as
`325(1 + i)/4` and `325(2 + i)/8`. So `1/4` is a limit of the method, consistent with `χ_c(ℚ(√59)²) = 4`. The
window (1, 1) does not suffice below `2/7` (at `r = 0.2501` it keeps extra components with largest least margins
`11/42` and `25/96`). Cross-checks: the lift computation (`kappa2d_seven.py`, `seven_r249_21.txt`), the run at
`r = 1/4 + 10⁻⁷` (13 vectors: c, q and 7), and the direct enumeration `indepN.py 249/1000 325`, which finds the 60
rotations of denominator dividing 325 and 29 components. The referee's vertex enumeration of the line arrangement (in
C with 128-bit integers, every point rechecked with fractions: `vertex_enum.c`, `analyze_components.py`) and its
lifting by rows (`lift_rows.py`) give the same 29 vectors at `r₀` and the same 13 at `1/4 + 10⁻⁷`
(`compare_methods.py`). The referee of the corollary on finite witnesses (§6.8) clipped `S^(2/7)(2, 1)` exactly,
cell by cell: 13 components (one 60-gon of family c, four 30-gons of family q, and the eight 7-points, each a single
point), with the stored index vectors.

> **Lemma F2.** Let `1/4 < r ≤ 1/3`, let `k` be a positive multiple of 6 (so `N = 5^k ≡ 1 mod 7`), and let `C ∈ ℂ` satisfy
> `Re(C̄γ) ∈ [r, 1 − r] + ℤ` for every `γ ∈ G(k, 1)`. Then modulo `Nℤ[i]`:
> - `C ≡ c*_N + x` with `x ∈ P_k^(s)`, or
> - `C ≡ q + x` with `q ∈ Q_N` and `x ∈ Y_k(q)`, or
> - `C ≡ Ne + x` with `e ∈ E₇` and `x ∈ X_k(e) = {x : ζ_j(e) + x̄ρʲ ∈ Sq_s for |j| ≤ k}`, where
>   `ζ_j(e) ∈ {±1/14, ±3/14}²` is given by `conj(Ne)ρʲ ∈ h + ζ_j(e) + ℤ[i]`.
>
> The third case does not occur if `r > 2/7`. In every case `|x| ≤ R′ = (1/4 + 3/14)√2 < 0.66`.

*Proof.* As for Lemma E2.
- *(b)* For `|j| ≤ k − 2`, `c = ρ^{−j}C ∈ S^(r)(2, 1)`, as `ρʲW ⊆ G(k, 1)`. By Proposition F1 the digits
  `ν_{j−1}, ν_j, ν_{j+1}, ν_{j+2}` of `C` are those of a type point `T` at `t = −1, 0, 1, 2`: `0000` for family c
  (`325ρᵗ` has odd norm); `(ν, −iν, −ν, iν)`, `ν ≠ 0`, for family q (as in Lemma E2(b), now with `325 ≡ 1 mod 3`); and
  four consecutive digits of a shift of `S` for family 7 (`325(a + bi) ≡ e (mod 7)` with `e ∈ A′₇`, so the values of
  `T` at `ρᵗ` are those of the character of `e/7`).
- *(c)* The windows `|j| ≤ k − 2` cover `ν_{−k+1}, …, ν_k`, and consecutive windows share three digits, at given
  positions in each window. Two consecutive digits at given positions determine the family and the class (for the
  7-adic classes, `ζ_t(e′) = ζ_{t+j}(e)` when `7e′ ≡ 7e·conj(ρ̃)^j (mod 7)`, `ρ̃ = 2 + 5i ≡ ρ`, as `ρ̃·conj(ρ̃) ≡ 1`). So
  if the window at `j` has family 7 and class `e_j`, the window at `j + 1` has family 7 and the class
  `e_{j+1} ≡ e_j·conj(ρ̃)`, and the digits of `C` are those of the class `e_0` at the same positions; likewise for
  the families c and q. The whole word is the zero word, a q-word, or the digits of one 7-adic class; by Proposition
  F1 the last case does not occur if `r > 2/7`.
- *(d)* The first two cases are as in Lemma E2(e). In the third, let `e ∈ E₇` be the class whose character has this
  shift. As `N ≡ 1 (mod 7)` and `Nρʲ ≡ ρʲ (mod 7ℤ_(7)[i])`, `conj(Ne)ρʲ = ē·Nρʲ` has, modulo `ℤ[i]`, the values of
  that character, so `g_j(Ne) = ζ_j(e)` has the digits of `C`, and `g_j − ζ_j(e) = ρ(g_{j−1} − ζ_{j−1}(e))`. With
  `x̄ = g₀ − ζ₀(e)`, `g_j = ζ_j(e) + x̄ρʲ`, so `x ∈ X_k(e)` and `C̄ρʲ ≡ conj(Ne + x)ρʲ`; as before
  `C − Ne − x ∈ Nℤ[i]`.

Bounds: `s√10/3`, `√2 ε`, and `|g₀| + |ζ₀(e)| ≤ (s + 3/14)√2`, all below `R′` as `s < 1/4`. ∎

### 6.6 Localisation at 2, 3 and 7

> **Lemma F3 (Lemma B1 with three types).** Let `i ∉ F`, and suppose `F²` maps to `K_{p/q}` with
> `1/4 < r = q/p ≤ 1/3`. Then for every finite `V ⊂ T` there is `β ∈ L` with `Tr(βv) ∈ E′` for every `v ∈ V`.

*Proof.* As the claim in §6.3, with `6 | k`, `N = 5^k > 42R′·max_v(|D_v| + Σ_j |μ_{v,j}|)`, Lemma F2 in place of
Lemma E2, and the exact-relations lemma with 42 in place of 6 (its proof uses only `42E′ ⊆ ℤ[i]`). ∎

> **Lemma F4.** Under the hypotheses of Lemma F3, (c) or (q) of Lemma B2 holds, or
> **(7′)** there is `β₇ ∈ L₇ = L ⊗ ℚ₇` with `Tr(β₇z) ∈ 𝒜₇` for every `z ∈ T₇ = ∏_{u | 7} T(F_u)`.

*Proof.* As Lemma B2, with a third prime. `K₇ = {β ∈ L₇ : Tr(βb_j) ∈ (1/7)ℤ₄₉ for all j}` is compact
(`Tr : L₇ → ℚ₇ ⊗ ℚ(i) = ℚ₄₉`). On `{c, q, 7}^T × K₂ × K₃ × K₇` the condition `P_v` reads:
`τ(v) = c`, `v_{1+i}(Tr β₂v) = −1`, `Tr β₃v ∈ O₃`, `Tr β₇v ∈ ℤ₄₉`; or `τ(v) = q`, `Tr β₂v ∈ O₂`, `Tr β₃v ∈ 𝒜`,
`Tr β₇v ∈ ℤ₄₉`; or `τ(v) = 7`, `Tr β₂v ∈ O₂`, `Tr β₃v ∈ O₃`, `Tr β₇v ∈ 𝒜₇`. These conditions are clopen, the `β` of
Lemma F3 satisfies them for `v ∈ V` (the elements of `E_c`, `E_q` and `E₇` satisfy them locally), and by
compactness some `(τ, β₂, β₃, β₇)` satisfies `P_v` for every `v ∈ T`. Let `a, a′` (at 2), `b, b′` (at 3) and
`c, c′` (at 7) be the local predicates; each pair is mutually exclusive (`c, c′` because `0 ∉ A′₇`). The clopen set
`Z = {(a ∧ b ∧ c) ∨ (a′ ∧ b′ ∧ c) ∨ (a′ ∧ b ∧ c′)} ⊆ T₂ × T₃ × T₇` contains the diagonal image of `T`, which is dense
(Lemma WA), so `Z` is everything. If `a(z₂⁰)` holds for some `z₂⁰`, then `b` and `c` hold everywhere, so `b′` and
`c′` hold nowhere and `a` holds everywhere: case (c). Otherwise `a′` holds everywhere; if `b′(z₃⁰)` holds for some
`z₃⁰`, then `c` holds everywhere, `c′` nowhere, and `b′` everywhere: case (q). Otherwise `b` and `c′` hold everywhere:
case (7′). ∎

Cases (c) and (q) give (a) and (b) as in the proof of Theorem B. For (7′) let `g_u(z) = Tr_{L_u/ℚ₄₉}(β_u z)` for
`u | 7`, so that `Tr(β₇y) = Σ_u g_u(y_u)`. At a place of even residue degree `i ∈ F_u` and `g_u = 0`, exactly as at 3
(the case "even f" in Section 5 of `papers/three-colours/`, with `7^m` in place of `3^m`). At a place of odd residue
degree `f`, `q = 7^f`, `L_u = F_u(i) = F_uℚ₄₉` is unramified of degree 2 over `F_u`, its non-trivial automorphism acts
on `ℚ₄₉` as the Frobenius, and `T(F_u) ⊂ O_{L_u}^×` contains the Teichmüller lifts of `μ_{q+1}`, among them `μ₈`
(`8 | q + 1`), and reduces onto `μ_{q+1}`; also `g_u(ζz) = ζg_u(z)` for `ζ ∈ μ₈ ⊂ ℚ₄₉`.

> **Lemma F5 (one place).** In case (7′) there are a place `u₀ | 7` of odd residue degree and `β ∈ L_{u₀}` with
> `Tr_{L_{u₀}/ℚ₄₉}(βz) ∈ 𝒜₇` for every `z ∈ T(F_{u₀})`.

*Proof.* As Lemma B3. Varying one component gives `g_u(z) − g_u(z′) ∈ 𝒜₇ − 𝒜₇ ⊆ (1/7)ℤ₄₉`; with `z′ = 1`, `z = i`
(and `i − 1` a unit at 7), `g_u(T(F_u)) ⊆ (1/7)ℤ₄₉`. Let `Z_u ⊆ 𝔽₄₉` be the set of residues of `7g_u(z)`,
`z ∈ T(F_u)`; it is `μ₈`-stable, and `Σ_u z_u ∈ A′₇` for every choice `z_u ∈ Z_u`. If `u ≠ u′` had non-zero `w ∈ Z_u`
and `w′ ∈ Z_{u′}`, fix the other components and put `c = ζ′w′ + Σ z_{u″}`, `ζ′ ∈ μ₈`: the 8 distinct elements
`c + ζw`, `ζ ∈ μ₈`, lie in `A′₇`, which has 8 elements, so they are all of `A′₇`, and summing gives
`8c = c = ΣA′₇ = 0` (`Σ_{ζ ∈ μ₈} ζ = 0`, characteristic 7). As this holds for every `ζ′ ∈ μ₈`, `w′ = 0`: a
contradiction. Since `0 ∉ A′₇`, exactly one place `u₀` has `Z_{u₀} ≠ {0}`, and it has odd residue degree. The other
`g_u` take values in `ℤ₄₉`, and `𝒜₇ + ℤ₄₉ = 𝒜₇`, so `β = β_{u₀}` works. ∎

> **Lemma F6 (level one).** Let `u | 7` have odd residue degree `f`, `q = 7^f`, and let `β ∈ L_u` with
> `Tr_{L_u/ℚ₄₉}(βz) ∈ 𝒜₇` for every `z ∈ T(F_u)`. Then `λ(x) = 7Tr(βx) mod 7` is defined on `O_{L_u}` and vanishes on
> its maximal ideal, and the induced `𝔽₄₉`-linear `λ₁ : 𝔽_{q²} → 𝔽₄₉` satisfies `λ₁(μ_{q+1}) ⊆ A′₇`.

*Proof.* The proof of Lemma B4 with `3, ℚ₉, 𝔽₉, A′` replaced by `7, ℚ₄₉, 𝔽₄₉, A′₇`. It uses that the prime is odd,
that the order of 7 modulo `q + 1` is `2f`, and (A2) in place of the fact that `A′` contains no coset of a non-zero
subgroup of `𝔽₉`. ∎

> **Lemma F7 (residue fields).** If `f` is odd and some `𝔽₄₉`-linear `λ₁ : 𝔽_{q²} → 𝔽₄₉`, `q = 7^f`, has
> `λ₁(μ_{q+1}) ⊆ A′₇`, then `f = 1`.

*Proof.* By (A3), `Re ∘ λ₁` is an additive map `𝔽_{q²} → 𝔽₇` with values in `{2, …, 5}` on `μ_{q+1}`; read in
`(1/7)ℤ/ℤ`, it stays at distance at least `2/7` from 0 there, so `κ₁ > 0` for the residue field `𝔽_q`. By
Proposition C3, `κ₁ = 0` for `p = 7` and `f ≥ 3`. ∎

For `f = 1`, multiplication by an element of `A′₇` works. The referee checked `f = 3` directly (no `𝔽₄₉`-linear map
`𝔽_{7⁶} → 𝔽₄₉` sends `μ₃₄₄` into `A′₇`, `indep_F/b5_f3_check.py`) and Lemma F6 on the finite models
`𝔽₄₉[π]/(πᵉ)`, `e ≤ 3` (`b4_model_check.py`).

*Proof of Theorem F.* Let `F²` map to `K_{p/q}` with `p/q < 4`. If `i ∈ F` this is impossible. If `p/q ≤ 3`, then
`χ(F²) ≤ 3` and (a) or (b) holds by Theorem B. Otherwise `r = q/p ∈ (1/4, 1/3)`, and by Lemmas F3 and F4 we are in
case (c), (q) or (7′). The first two give (a) or (b) as in the proof of Theorem B; in case (7′) Lemmas F5–F7 give a
place above 7 of residue degree 1, which is (7). ∎

> **Corollary F8 (local fields; local–global below 4).** Let `K` be a finite extension of `ℚ_p`. Then
> `χ_c(K²) = 2` if `p = 2` and `K(i)/K` is ramified; `3` if `p = 3` and `K` has residue degree 1; `7/2` if `p = 7`
> and `K` has residue degree 1; and `χ_c(K²) ≥ 4` in every other case. Consequently, for every number field `F`,
> `min(χ_c(F²), 4) = min(4, min_v χ_c(F_v²))`, the inner minimum over all places `v` of `F`: below 4 the circular
> chromatic number of `F²` is the least one of the planes over its completions.

*Proof.* The upper bounds in the first three cases are the local colourings (the 2-colouring of the two-colour paper,
the 3-colouring of Theorem B, Proposition C1), which colour `K²` itself; the lower bounds come from a bipartite plane
with an edge, from `ℚ(√7) ⊂ ℚ₃ ⊆ K` with `χ_c(ℚ(√7)²) = χ(ℚ(√7)²) = 3`, and from Corollary C5. Otherwise, if `i ∈ K`
then `χ = ∞`; if not, the Krasner construction of the corollary on local fields in `papers/three-colours/` (with a
product of quadratics generating `ℚ₄₉` added over `ℚ₇`, a cubic generating the unramified cubic extension of `ℚ₇` when
`p = 7` and `[K : ℚ₇]` is odd, and the auxiliary prime chosen outside `{2, 3, 7, p}`) gives a number field `F ⊂ K` for
which (a), (b) and (7) fail, so `χ_c(K²) ≥ χ_c(F²) ≥ 4` by Theorem F. For the last statement, `F² ⊆ F_v²` gives
`χ_c(F²) ≤ χ_c(F_v²)`, and if `χ_c(F²) < 4`, Theorem F gives a place where the completion has the same value. ∎

For example `χ_c(ℚ₁₉²)`, `χ_c(ℚ(√47)²)` and `χ_c(ℚ(√5, √7)²)` lie in `[4, 19/4]`: the lower bounds are Corollary F8 and
Theorem F, the upper bound Proposition C1 at 19 (`κ₁ = 4/19`; `47 ≡ 3² (mod 19)`, and 5 and 7 are squares modulo 19).
So if `χ(ℚ(√47)²) = 5`, which is open (`4 ≤ χ ≤ 5`: Corollary 2 of `papers/four-colours/`, and Moorhouse), its
circular chromatic number lies in `(4, 19/4]` and is not an integer; the same holds for `ℚ(√5, √7)`, where
`4 ≤ χ ≤ 5` (Theorem B and the place above 19).

Theorem E also follows from Proposition F1 (as the referee of Theorem F noted): for `r > 2/7` it leaves only the
families c and q, so Lemma F2 (with `6 | k`) gives the conclusion of Lemma E2, and the proof of Theorem E goes through.
The certificates that exclude the 8 vectors of family 7 there are not those of Proposition E1 (they use the rotations
`ρ^{−2}σˡ`). So Theorem F, whose proof uses Theorem E, rests on Proposition F1 alone, and Theorem E on either
proposition.

### 6.7 Why one prime is not enough

Context only; nothing above depends on it. Below `3/10` the one-prime probe has characters with no bounded type. At
`r = 3/10`, for every `t`, the sequence `g_j = ρ^{j−t}/5` (`j < t`), `g_t = −(1 + i)/5`, `g_j = iρ^{j−t}/5` (`j > t`)
is a character of `ℤ[1/5][i]` whose point at level `K` needs a type with denominator exactly `5^{K−t+1}`
(`glue_denominators.py`, exact for `K ≤ t + 9`); the points `c_k` of §5 are of this kind. Torsion characters of order
prime to 5 also beat `2/7` on `{iᵃρʲ}`: up to order 400, exactly those of orders 41 (`κ = 12/41`, for example
`c = (12 + 12i)/41`; `ρ` has order 5 modulo 41) and 76 (`κ = 11/38`), counted again at the multiples of these orders
(`torsion_exact.py`; the referee's `indep_E/ref_torsion.py` counts by exact order). This corrects the earlier file
`probe/torsion_260.txt` (research log, 4 October). Over `⟨i, ρ, σ⟩` none survives for orders up to 400 prime to 65
(the orbit of `12 + 12i` modulo 41 reaches 0; `torsion_twoprime.py`, `indep_E/ref_torsion.py`), and for `k ≤ 5` every
translate of the points `c_k` has least margin at most `41/250` on `G(k, 1)` (`indep_E/ref_ck.py`). By Lemma E2, a
torsion character `ξ` that keeps all of `⟨i, ρ, σ⟩` above `2/7` is of type c or q at every level, so `ξ(5^{−k})` and
`ξ(i5^{−k})` lie within `0.23·5^{−k}` of `1/2`, `1/3` or `2/3` for every `k`; as `ξ` has finite order, it has order 2
or 3.

### 6.8 Finite witnesses

By compactness `χ_c(F²)` is the supremum of `χ_c(H)` over the finite subgraphs `H` of `F²` (if every finite subgraph
maps to `K_{p/q}`, so does `F²`). Below 4 the supremum is attained (Corollary F12), and at 4 as well (§6.9). For a homomorphism `c` of a graph to `K_{p/q}`,
the *tight digraph* `D_c` has an arc `x → y` for every edge `xy` with `c(y) − c(x) ≡ q (mod p)`; a *tight cycle* is a
directed cycle of `D_c`. For a finite set `S` in an abelian group, `κ(S) = max_ξ min_{s∈S} ‖ξ(s)‖` over the characters
of `ℤS` (a maximum, as the dual group is compact); if `χ_c(Cay(ℤS, S)) < 4`, then `χ_c = 1/κ(S)` (Corollary 2 of the
winding paper).

> **Lemma F9 (tight cycles; the easy half of a characterisation of Guichard, J. Graph Theory 17 (1993); cf. Zhu's
> survey).** Let `H` be a finite graph with a homomorphism to `K_{p/q}`. If every homomorphism of `H` to `K_{p/q}`
> has a tight cycle, then `χ_c(H) = p/q`.

*Proof.* Suppose `χ_c(H) < p/q`: then `H` maps to some `K_{p′/q′}` with `p′/q′ < p/q`, by `c′`. Let
`g(x) = (p/p′)c′(x) ∈ [0, p)` (reading `c′(x)` in `{0, …, p′ − 1}`) and `c(x) = ⌊g(x)⌋ mod p`. For an edge `xy` the
lift `D ∈ [pq′/p′, p − pq′/p′]` of `g(y) − g(x)` modulo `p` has `q < D < p − q`, as `pq′/p′ > q`; the integer
`⌊g(x) + D⌋ − ⌊g(x)⌋` is `⌊D⌋` or `⌈D⌉`, so it lies in `[q, p − q]`, and it is congruent to `c(y) − c(x)`. So `c` maps
`H` to `K_{p/q}`, and `x → y` is tight exactly when this integer is `q`. Along a closed walk these integers have the
same sum as the lifts `D` (the floors telescope); on a tight cycle of length `m` the first sum is `mq` and the second
is larger. So `c` has no tight cycle. ∎

> **Lemma F10 (winding).** Let `Γ` be an abelian group, `S = −S ⊆ Γ` finite with `Γ = ℤS`, `p < 4q`, and `c` a
> homomorphism of `Cay(Γ, S)` to `K_{p/q}`. For `x ∈ Γ`, `s ∈ S` let `ℓ(x, s) ∈ [q, p − q]` be the integer congruent
> to `c(x + s) − c(x)` modulo `p`, and `a(s) = M_x ℓ(x, s)` for an invariant mean `M` on `Γ`. Then `a(−s) = p − a(s)`,
> and `s_1 + … + s_m ↦ (1/p) Σ a(s_k) mod 1` is a well-defined character `ξ` with `ξ(s) = a(s)/p ∈ [q/p, 1 − q/p]` on
> `S`; and if `u_1, …, u_r ∈ S` and positive integers `n_i` satisfy `Σ n_i u_i = 0` and `a(u_i) = q` for every `i`,
> every closed walk with `n_i` steps `u_i` is a closed walk in `D_c`, so `c` has a tight cycle.

*Proof.* This is the proof of Theorem W⁺ (`notes/winding_lemma.md`). For `s, t ∈ S` the integers
`ℓ(x, s) + ℓ(x + s, t)` and `ℓ(x, t) + ℓ(x + t, s)` are congruent modulo `p` and lie in `[2q, 2p − 2q]`, of length
`2p − 4q < p`, so they are equal. For a closed walk `W` with steps `s_1, …, s_m` from `x` (`x_k = x + s_1 + … + s_k`),
`Λ(W, x) = Σ_k ℓ(x_{k−1}, s_k)` is a multiple of `p`, and for `t ∈ S` the identity gives
`Λ(W, x + t) − Λ(W, x) = Σ_k (ℓ(x_k, t) − ℓ(x_{k−1}, t)) = 0`; so `Λ(W, x)` does not depend on `x`, and by
invariance it equals `Σ_k a(s_k)`. As `ℓ(x + s, −s) = p − ℓ(x, s)`, invariance gives `a(−s) = p − a(s)`, and `ξ` is
well defined. For the walk of the statement, `Λ(W, x) = q Σ n_i` for every `x`, and each of its `Σ n_i` terms is at
least `q`, so each is `q`. A closed walk in `D_c` contains a directed cycle. ∎

> **Lemma F11 (finite connection sets).** Let `Γ` be an abelian group, `S = −S ⊆ Γ` finite with `Γ = ℤS`, and
> `χ_c(Cay(Γ, S)) = p/q` with `2 < p/q < 4`. Then every homomorphism of `Cay(Γ, S)` to `K_{p/q}` has a tight cycle,
> and some finite subgraph of `Cay(Γ, S)` has circular chromatic number `p/q`.

*Proof.* Here `κ(S) = q/p`, so `Cay(Γ, S)` maps to `K_{p/q}`. Let `c` be such a homomorphism, `a`, `ξ` as in Lemma
F10, and `T = {s ∈ S : a(s) = q}`. Suppose no nonzero nonnegative integers `n_t` give `Σ_{t∈T} n_t t = 0`. Then the
images of `T` in `Γ ⊗ ℝ` have no nonzero nonnegative real relation (the cone of such relations is cut out by rational
equations, so a nonzero one would give a nonzero integral `(n_t)` with `Σ n_t t` torsion, killed by some `m ≥ 1`). By
Gordan's theorem there is a homomorphism `h : Γ → ℝ` with `h(t) > 0` on `T`. For small `ε > 0` the character
`ξ + εh mod 1` maps every element of the finite set `S` into the open interval `(q/p, 1 − q/p)` (nonempty, as
`p > 2q`): on `T` it is `q/p + εh(t)`, on `−T` it is `1 − q/p − εh(−t)`, and elsewhere `ξ` is already inside. Then
`κ(S) > q/p`: a contradiction. So such `n_t` exist, and Lemma F10 (on the `t` with `n_t > 0`) gives a tight cycle of
`c`. If every finite subgraph had a homomorphism to `K_{p/q}` without tight cycles, so would `Cay(Γ, S)` (for each
finite `H ⊂ Γ` the maps `Γ → ℤ/p` that are such a homomorphism on `H` form a nonempty closed subset of the compact
space `(ℤ/p)^Γ`; these have the finite intersection property, and a tight cycle is finite). So some finite subgraph
has a tight cycle for every homomorphism to `K_{p/q}`, and Lemma F9 applies. ∎

> **Corollary F12 (finite witnesses).** If `χ_c(F²) < 4`, some finite subgraph of `F²` has circular chromatic
> number `χ_c(F²)`. If `χ_c(F²) ∈ {3, 7/2}`, there is a finite `U ⊂ T` with `χ_c(Cay(ℤU, U)) = χ_c(F²)`, and the
> subgraph can be taken in `Cay(ℤU, U)`.

*Proof.* The value is 2, 3 or `7/2` (Theorem F), and `i ∉ F`; for 2 an edge is a witness. Fix a `ℚ(i)`-basis
`b_1, …, b_n ⊂ T` of `L` as in Lemma B2. It suffices to find a finite `U ⊂ T` with `χ_c(Cay(ℤU, U)) = χ_c(F²)`
(Lemma F11).

*`7/2`.* (7) holds and (a), (b) fail. If for every finite `V ⊇ {b_j}` there were `β ∈ L` with `Tr(βv) ∈ E` on `V`,
the proof of Lemma B2 and the proof of Theorem B would give (a) or (b); so some finite `V_0 ⊇ {b_j}` has no such `β`.
Let `k` (with `6 | k`) and `N` be as in Lemma F3 for `V_0` (they do not depend on `r`), and `U_0 = G(k, 1)V_0`. If
`κ(U_0) > 2/7`, a character with `ξ(U_0) ⊆ [r, 1 − r] + ℤ`, `2/7 < r ≤ min(κ(U_0), 1/3)`, gives through Lemma F2 (whose
third case does not occur for `r > 2/7`) and the proof of Lemma F3 a `β` with `Tr(βv) ∈ E` on `V_0`: impossible. So
`κ(U_0) ≤ 2/7`; Proposition C1 and Theorem W⁺ give `κ(U_0) ≥ 2/7`; so `χ_c(Cay(ℤU_0, U_0)) = 7/2`.

*`3`.* (b) holds and (a) fails. As in the proof of Corollary B6, if for every finite `V ⊇ {b_j}` there were `β` with
`Tr(βv) ∈ E_c` on `V`, then (c), hence (a), would hold; so some finite `V_0` has no such `β`. With `N` as in Lemma B1
for `V_0` and `U_0 = G_N V_0`: if `κ(U_0) > 1/3`, a character with `ξ(U_0) ⊆ [r, 1 − r] + ℤ`, `1/3 < r ≤ κ(U_0)`,
gives as in the proof of Corollary B6 points `φ(v) = Nε_v + x_v` with every `ε_v` of type c (type q has real part in
`{1/3, 2/3} + ℤ`, as `3 ∤ N`), hence such a `β`: impossible. So `κ(U_0) ≤ 1/3`, and `κ(U_0) ≥ 1/3` as `F²` is
3-colourable; so `χ_c(Cay(ℤU_0, U_0)) = 3`. ∎

The proof uses Proposition F1 only for `r > 2/7`, through Lemma F2, and not the certificates of value `2/7`
themselves. For `ℚ(√11)` and `ℚ(√35)` one can take for `U` the 140 vectors of Theorem C (its certificates give
`κ(U) ≤ 2/7`, the 7-adic character `κ(U) ≥ 2/7`): some finite subgraph of the Cayley graph of these 140 explicit unit
vectors has `χ_c = 7/2`. A finite extension `K` of `ℚ_p` with `χ_c(K²) < 4` (Corollary F8) has finite
witnesses too: an edge for the value 2, and for 3 and `7/2` those of Corollary F12 for `ℚ(√7)`, resp. `ℚ(√11)`, as
then `p = 3` and `K ⊇ ℚ₃ ⊃ ℚ(√7)`, where `χ_c(ℚ(√7)²) = 3` (Theorem B and Corollary B6), resp. `p = 7` and
`K ⊇ ℚ₇ ⊃ ℚ(√11)`. The proof gives no bound on the size of the subgraph; when `√3 ∈ F` and `χ_c(F²) = 3`, a
unit triangle is a witness, and over `ℚ(√7)`, which has none, there is one with nine vertices (below). The Moser spindle (7 vertices, over `ℚ(√3, √11)`) maps to `K_{7/2}` (the common vertex of
its two rhombi coloured 0, the far tips 6 and 1, the others 2, 4, 3, 5), and its independence number is 2, so
`χ_c ≥ χ_f ≥ 7/2` and `χ_c = 7/2`; but it has unit triangles, which `ℚ(√11)²` has not, and `χ_c(ℚ(√3, √11)²) ≥ 4` by
Theorem F ((a) and (b) fail as `χ = 4` there, and 7 is inert in `ℚ(√3)`). Measurements (SAT, colourings checked,
refutations not certified): the 76-vertex graph of `χ(ℚ(√11)²) = 4` has `χ_c = 16/5` (the decisive refutation, at
`67/21`, took 22 minutes; a referee's run there did not finish); the unit-distance graph induced on the union of its nine
images under `ρʲσˡ`, `|j|, |l| ≤ 1` (628 vertices, 1 506 edges; the nine copies alone have 1 494), maps to `K_{13/4}`
and not to `K_{16/5}`. Exactly (`finite_ball.py`): the ball of radius 2 in the Cayley graph of the 140 vectors
(9 941 vertices, 19 600 edges) is bipartite, and the unit-distance graph it induces (216 more edges) has `χ_c = 5/2`:
it has a 5-cycle and no triangle, and maps to `K_{5/2}`. **An explicit witness for `7/2`** (4 October; found by a
search with iterated Minkowski sums and lazy SAT, `data/number_fields/circular/finite_witness/`): with `A` the vertex
set of the 76-vertex graph, the unit-distance graph `H` on `A + A` (2 237 vertices, all 11 300 unit pairs) maps to
`K_{7/2}` by a 7-adic colouring, and every homomorphism of `H` to `K_{7/2}` has a tight cycle among 180 listed cycles
of lengths 14 to 42; so `χ_c(H) = 7/2` by Lemma F9. **Smaller witnesses** (4 October, evening; `grow.py`,
`minimise.py`): a colouring-guided growth with lazy SAT, started from `A`, asks for a `(7, 2)`-colouring in which every
listed directed cycle has a non-tight arc, lists the tight cycles of each answer, and, when the tight digraph of the
answer is acyclic, adds the points `x + u` (`x` a vertex, `u` a unit vector with the denominator of the starting
graph, here 30) whose neighbours leave no colour (or, if there are none, those whose neighbours leave one colour); it
stops after 52 rounds at 653 vertices. Deleting vertices while the formula stays unsatisfiable
leaves `H₁₁`: 170 vertices, 468 edges (all unit pairs), `χ_c(H₁₁) = 7/2`, and vertex-critical: for every vertex `v` a
`(7, 2)`-colouring of `H₁₁ − v` with an acyclic tight digraph is stored, so `χ_c(H₁₁ − v) < 7/2` (the other half of
Guichard's characterisation, explicitly: if `G` has `N` vertices and `pos` numbers them along a topological order of
the tight digraph of `c`, then `N·c + pos` is a homomorphism of `G` to `K_{7N/(2N+1)}`).
(The same growth with at most 100 new points per round stops at 573 vertices, and deletion leaves another
vertex-critical witness with 170 vertices; deletion from the union of the two, in two random orders, leaves ones with
157 and 161 vertices, and deletion from the union of these two leaves one with 155 vertices and 404 edges,
`witness_q11b.json.gz`, the smallest we know.) From the 71-vertex
graph over `ℚ(√455)` and the 96-vertex graph over `ℚ(√191)` the same procedure gives vertex-critical witnesses with 175
vertices (434 edges) and 293 vertices (803 edges) (denominators 780 and 240); over `ℚ(√911)` (327 vertices, denominator
1560) the growth stops at 873 vertices, and deletion with the cores of `kissat` refutations (extracted by `drat-trim`)
and then of CaDiCaL refutations leaves a vertex-critical witness with 324 vertices and 866 edges. So `χ_c = 7/2` for
these planes without Proposition F1: the witnesses give `≥ 7/2`, and the residue field `𝔽₇` gives `≤ 7/2`, as 7
ramifies in `ℚ(√455)` and splits in `ℚ(√191)` and `ℚ(√911)` (in agreement with Theorem F: (a) and (b) fail for
`d ≡ 11 mod 12`). Every lower bound is a SAT
computation, certified twice: the first checker (`check_witness.py`) and a second one sharing no code with it
(`verify_independent.py`) check the graph, the colouring and the cycles exactly and write the formula with different
variable numberings and arc indicators; `drat-trim` verifies the stored DRAT proof of the
first, and a new `kissat` proof of the second, which `cake_lpr` (a formally verified checker) also accepts in LRAT
form (`verification.txt`); `check_critical.py` checks the criticality certificates. A referee checked the witnesses
with 170, 175 and 293 vertices again with programs of its own (`finite_witness/indep_W/`: a third encoding, cycle
lists rebuilt from scratch, `cake_lpr` on the stored proofs, the explicit homomorphisms `N·c + pos`), found no error,
and reproduced `H₁₁` and the `ℚ(√455)` witness exactly; those with 155 and 324 vertices (over `ℚ(√11)` and `ℚ(√911)`),
found afterwards, were checked by the two programs above and, by us, with the referee's programs, with the same
results (`indep_W/results/witness_q11b.log`, `witness_q911.log`). As `√11 ∈ ℚ₇`, `H₁₁` is an explicit witness for `ℚ₇`
and for every finite extension `K` of `ℚ₇` with `χ_c(K²) < 4`. **A witness for 3 without
unit triangles** (4 October, night; `check_small.py`, `small_triangle_free.py`, `nine_vertices.py`,
`hexagon_search.py`): the same procedure with `K₃` in
place of `K_{7/2}` (an arc is tight when the colour increases by 1 mod 3), over `ℚ(√7)`, from 0, the 204 unit vectors
with denominator 160 and two points closing a 5-cycle, stops at 607 vertices, and deletion leaves nine points:
`P₀ = (−1, 0)`, `P₁ = (−(3 + √7)/8, −(5 + √7)/8)`, `P₂ = (1/4, −√7/4)`, `P₃ = (√7/4, 1/4)`, `P₄ = (−1/4, √7/4)`,
`P₅ = ((3 − √7)/8, (√7 − 5)/8)`, `P₆ = (1, 0)`, `P₇ = (0, 0)`, `S = (0, 1)`. Their unit pairs form the 8-cycle
`P₀⋯P₇`, the chords `P₀P₄`, `P₁P₅`, `P₂P₆` and the path `P₃SP₇`: `H₇` is the Wagner graph (the Möbius ladder on 8
vertices) with one chord subdivided. Each of its 84 proper 3-colourings has a tight 6-cycle, so `χ_c(H₇) = 3` by
Lemma F9. By hand, already for the subgraph `M = H₇ − P₁P₅`, the hexagon `v₀⋯v₅ = P₀P₄P₃P₂P₆P₇` with its long diagonals
`v_j v_{j+3}` subdivided by `m₀, m₁, m₂ = P₁, P₅, S`: put `δ(a, b) = ±1` with `c(b) − c(a) ≡ δ(a, b) (mod 3)` on each
arc; along a closed walk the sum of `δ` is divisible by 3, so it is ±3 on a 5-cycle and 0 or ±6 on a 6-cycle (±6:
tight in one direction). For the 5-cycles `Z_k = v_k v_{k+1} v_{k+2} v_{k+3} m_k`, as 1-chains `C_k = Z_{k+1} − Z_k` is
the 6-cycle `v_k m_k v_{k+3} v_{k+4} m_{k+1} v_{k+1}` and `Z₀ + Z₃` is the hexagon; `C_{k+3}` is `C_k` reversed, and
`Z₃ − Z₀ = C₀ + C₁ + C₂`. If none of `C₀, C₁, C₂` is tight in either direction, their sums are 0, `Z₀` and `Z₃` have
the same sum ±3, and the hexagon has sum ±6, so it is tight. (In `H₇`, `C₀ = P₀P₁P₂P₆P₅P₄` has the chord `P₁P₅` and is
never tight; the cycles stored with `H₇` are `C₁`, `C₂` and the hexagon, in both directions.) Equivalently `H₇` has no homomorphism to `K_{8/3}`, which is the Wagner graph itself (`χ_c` of a graph
with 9 vertices is a fraction with numerator at most 9, by Vince and by Bondy and Hell, and `8/3` is the largest one
below 3). `H₇` is
vertex-critical, and nine vertices are the fewest possible without unit triangles: every triangle-free graph with at
most 8 vertices maps to `K_{8/3}` (checked over the 4 682 270 triangle-free graphs on 8 labelled vertices). Of the
1 897 triangle-free graphs with 9 vertices (up to isomorphism) exactly three have `χ_c = 3` (`nine_vertices.py`, two
separate tests): `M`, `H₇ = M + m₀m₁` and `M + m₀m₁ + m₁m₂`. All three are unit-distance graphs over `ℚ(√7)` and over
`ℚ(√31)`: without unit triangles, nine points with the twelve unit pairs of `M` span `M` or `M` plus one or two
edges between midpoints (any other unit pair would close a triangle), and a direct search for equilateral hexagons
whose long diagonals are sums of two unit vectors (`hexagon_search.py`) finds all three over both fields. For
example, the hexagon `(√7/4, 1/4)`, `(√7/4 − 1, 1/4)`, `(−3/4, −√7/4)`, `(0, 0)`, `(−1, 0)`, `(−1/4, √7/4)` with
midpoints `(0, 1)`, `(√7/4 − 1, −3/4)`, `(−1/2 − √7/4, 1/4)` spans `M` itself (`witness_q7m.json.gz`), a witness with
nine vertices and twelve edges; the growth over `ℚ(√31)` gives nine points spanning `M + m₀m₁ + m₁m₂`
(`witness_q31.json.gz`, `q31_seed.json`). All of this is checked by enumeration, the tight cycles
also by certified SAT refutations; as `√7 ∈ ℚ₃`, `H₇` is an explicit witness for `ℚ₃` and its finite extensions `K`
with `χ_c(K²) < 4`. Lemmas F9–F11 were tested by two referees: on 291 small graphs and 766 tight relations; on 115
random finite sets `S` in `ℤ`, `ℤ²` and `ℤ × ℤ/m` with `2 < χ_c < 4` (at all 409 optimal points, vertices of the optimal set and
midpoints and centroids of them that stay optimal, the tight elements have a nonnegative relation;
`indep_FW2/test_kappa_relations.py 1 180`); and by SAT on finite pieces of 13 such Cayley graphs, all with finite witnesses
(for example `{0, …, 18}` for the distances `3, 4, 9, 12`, where `χ_c = 7/2`); `twoprime/indep_FW/`,
`twoprime/indep_FW2/`.

### 6.9 Finite witnesses at 4

*Refereed (6 October) by an independent checker, with two programs of its own on all abelian groups of order at
most 16 (`at_four/finite_groups/`): no error; its corrections are applied below.*

At `p = 4q` the proof of Lemma F10 breaks in one place only, and what breaks it is itself a tight cycle. For a
homomorphism `c` of `Cay(Γ, S)` to `K_4 = K_{4/1}`, with `ℓ(x, s) ∈ {1, 2, 3}` as in Lemma F10, a *tight square* is
a tight cycle `x → x + s → x + s + t → x + t → x` with `s, t ∈ S`, that is, `ℓ(x, s) = ℓ(x + s, t) = 1` and
`ℓ(x, t) = ℓ(x + t, s) = 3`. Its four vertices have the colours `c(x)`, `c(x) + 1`, `c(x) + 2`, `c(x) + 3`, so they
are distinct.

> **Lemma F13 (at 4: a tight square or a character).** Let `Γ` be an abelian group, `S = −S ⊆ Γ ∖ {0}` with
> `Γ = ℤS`, and `c` a homomorphism of `Cay(Γ, S)` to `K_4` without tight squares. With `a(s) = M_x ℓ(x, s)` as in
> Lemma F10, the conclusions of Lemma F10 hold with `p = 4`, `q = 1`: `ξ(s) = a(s)/4` defines a character of `Γ` with
> `ξ(S) ⊆ [1/4, 3/4]`, and if `u_1, …, u_r ∈ S` and positive integers `n_i` satisfy `Σ n_i u_i = 0` and `a(u_i) = 1`,
> then `c` has a tight cycle. Conversely, if a character `ξ` maps `S` into `[1/4, 3/4]`, then
> `c(x) = ⌊4ξ̃(x)⌋ mod 4` (`ξ̃(x) ∈ [0, 1)` the lift of `ξ(x)`) is a homomorphism to `K_4` without tight squares.

*Proof.* For `x ∈ Γ` and `s, t ∈ S` the integers `A = ℓ(x, s) + ℓ(x + s, t)` and `B = ℓ(x, t) + ℓ(x + t, s)` are
congruent modulo 4 and lie in `[2, 6]`. So `A = B` or `{A, B} = {2, 6}`; exchanging `s` and `t` exchanges `A` and
`B`, and `A = 2`, `B = 6` says exactly that `x → x + s → x + s + t → x + t → x` is a tight square (as
`ℓ(y + u, −u) = 4 − ℓ(y, u)`). So without tight squares `A = B` for all `x, s, t`, which is all that the proof of
Lemma F10 uses of `p < 4q`; the rest of that proof (independence of the base point, the invariant mean, the closed
walks) applies word for word. Conversely, for `s ∈ S` let `σ ∈ [1/4, 3/4]` be the lift of `ξ(s)`. As
`ξ̃(x + s) − ξ̃(x) − σ ∈ ℤ`, `ℓ(x, s) = ⌊4ξ̃(x) + 4σ⌋ − ⌊4ξ̃(x)⌋ ∈ {⌊4σ⌋, ⌊4σ⌋ + 1}` (it is 3 if `4σ = 3`), which
lies in `[1, 3]`; so `c` is a homomorphism to `K_4`, while a tight square would need `ℓ(x, s) = 1` and
`ℓ(x + t, s) = 3`. ∎

So at 4 Theorem W⁺ holds for homomorphisms without tight squares: `Cay(Γ, S)` has a homomorphism to `K_4` without
tight squares if and only if some character maps `S` into `[1/4, 3/4]`. Finiteness of `S` is not used, and the same
proof works for `K_{4q/q}` with any `q` (tight arcs are those of difference `q`). In the example
`K_4 = Cay((ℤ/2)², (ℤ/2)² ∖ {0})` of the winding paper every 4-colouring is a tight square. A tight 4-cycle need not
be a square: for `Γ = ℤ/8`, `S = {±1, ±2}` and `c(2k) = k`, `c(2k + 1) = k + 2`, the cycle `0 → 2 → 4 → 6 → 0` is
tight, there is no tight square, and `a/4` is the character `x ↦ 5x/8`, with `κ(S) = 1/4`.

> **Corollary F14 (finite connection sets at 4).** Let `Γ`, `S` be as in Lemma F13 with `S` finite, `Cay(Γ, S)`
> 4-colourable and `κ(S) ≤ 1/4`. Then every homomorphism of `Cay(Γ, S)` to `K_4` has a tight cycle, and some finite
> subgraph of `Cay(Γ, S)` has circular chromatic number 4. If `κ(S) < 1/4`, every such homomorphism has a tight square.

*Proof.* Let `c` be a homomorphism without tight cycles, `ξ`, `a` as in Lemma F13, and `T₁ = {s ∈ S : a(s) = 1}`
(disjoint from `−T₁`, as `a(−s) = 4 − a(s)`). If nonnegative integers `n_t`, not all 0, gave
`Σ_{t∈T₁} n_t t = 0`, Lemma F13 would give a tight cycle. Otherwise, as in the proof of Lemma F11 (the cone of
nonnegative real relations among the images of `T₁` in `Γ ⊗ ℝ` is rational, and a nonzero integral point `(n_t)` of
it has `Σ n_t t` torsion, killed by some `m ≥ 1`), Gordan's theorem gives a homomorphism `h : Γ → ℝ` with `h > 0` on
`T₁` (any `h` if `T₁` is empty). For small `ε > 0`, `ξ + εh` maps every `s ∈ S` into `(1/4, 3/4) + ℤ`: if `a(s) = 1`
it gives `1/4 + εh(s)`; if `a(s) = 3`, then `−s ∈ T₁` and it gives `3/4 − εh(−s)`; otherwise `ξ(s)` is already
inside. So `κ(S) > 1/4`, a contradiction. The rest is the compactness argument of Lemma F11, with Lemma F9 for
`p/q = 4`. If `κ(S) < 1/4`, Lemma F13 itself excludes homomorphisms without tight squares. ∎

Conversely, if a finite subgraph `H` of `F²` has `χ_c(H) = 4`, then `κ(U_H) ≤ 1/4` for the set `U_H` of its edge
vectors and their negatives: each component of `H` lies in a coset of `ℤU_H`, and if `κ(U_H) > 1/4`, then for a
rational `q/p ∈ (1/4, κ(U_H)]` a character with `ξ(U_H) ⊆ [q/p, 1 − q/p]` gives a homomorphism of `H` to `K_{p/q}`
(the *if* direction of Theorem W⁺ holds for every `p/q ≥ 2`), so `χ_c(H) ≤ p/q < 4`. So, when `χ(F²) ≤ 4`, some
finite subgraph of `F²` has `χ_c = 4` if and only if some finite `U ⊂ T` has `κ(U) ≤ 1/4`.

> **Proposition F15 (a finite set at 1/4).** Let `F` be a number field with `i ∉ F` for which (a), (b) and (7) fail.
> Then some finite `U ⊂ T` has `κ(U) ≤ 1/4`: no character of `ℤU` maps every element of `U` into the open interval
> `(1/4, 3/4)`.

*Proof.* As in the proof of Corollary F12 for `7/2`. If for every finite `V ⊇ {b_j}` there were `β ∈ L` with
`Tr(βv) ∈ E′` for every `v ∈ V`, the proof of Lemma F4, which uses Lemma F3 only through this conclusion (such a `β`
lies in `K₂ × K₃ × K₇`, and the types of the `Tr(βv)` give `τ`), and Lemmas F5–F7 would give (c), (q) or (7′), hence
(a), (b) or (7). So some finite `V_0 ⊇ {b_j}` has no such `β`. Let `k` be the least multiple of 6 with
`5^k > 42R′ max_{v ∈ V_0}(|D_v| + Σ_j |μ_{v,j}|)` (in the fixed basis; `R′` does not depend on `r`), `N = 5^k`, and
`U_0 = G(k, 1)V_0`, a finite symmetric subset of `T`. If `κ(U_0) > 1/4`, take `r` with
`1/4 < r ≤ min(κ(U_0), 1/3)` and a character `ξ` of `ℤU_0` with `ξ(U_0) ⊆ [r, 1 − r] + ℤ`; the proof of Lemma F3,
which uses the homomorphism only through such a character (Proposition F1 excludes its 16 extra vectors for each
`r > 1/4`, and Lemma F2 needs only `s = 1/2 − r < 1/4`), gives `β` with `Tr(βv) ∈ E′` on `V_0`: impossible. So
`κ(U_0) ≤ 1/4`. ∎

(For `i ∉ F` the converse holds as well: if (a), (b) or (7) holds, the corresponding colouring and Theorem W⁺ give
`κ(U) ≥ 2/7` for every finite `U ⊂ T`.)

> **Theorem F16 (finite witnesses up to 4).** For every number field `F` with `χ_c(F²) ≤ 4`, some finite
> unit-distance graph in `F²` has circular chromatic number `χ_c(F²)`. In particular the value 4 is attained by
> finite unit-distance graphs over `ℚ(√59)`, `ℚ(√83)` and `ℚ(√3, √11)`.

*Proof.* Below 4 this is Corollary F12. If `χ_c(F²) = 4`, every finite subgraph has `χ_c ≤ 4`, hence `χ ≤ 4`, and
`χ(F²) ≤ 4` by the theorem of de Bruijn and Erdős; so `i ∉ F`, and (a), (b) and (7) fail by Theorem F, as otherwise
`χ_c(F²) < 4`. Proposition F15 gives a finite `U ⊂ T` with `κ(U) ≤ 1/4`, `Cay(ℤU, U) ⊆ F²` is 4-colourable, and
Corollary F14 gives a finite subgraph `H` with `χ_c(H) = 4`; the unit-distance graph induced on its vertices also has
`χ_c = 4`, as `4 = χ_c(H) ≤ χ_c ≤ χ_c(F²) = 4`. For `ℚ(√59)` and `ℚ(√83)` see the examples after Theorem F; for
`ℚ(√3, √11)`, `χ = 4` and (7) fails, as 7 is inert in `ℚ(√3)`. ∎

The proof gives no bound on the size of the witness, and we have not found one: the growths of the research log (6
October), with the unit vectors of denominator `D = 210`, `1050` and `2730` over `ℚ(√59)`, stopped at their limits.
For `D = 210` this was bound to happen: the character `θ = (89/118, 1/2, 89/118, 1/2)` of `at_four/` keeps all 108
unit vectors with denominator 210 at distance at least `15/59` from `ℤ`, so every graph built from them maps to
`K_{59/15}` (checked by enumeration, and by the referee with another character).

## 7. Questions

1. *Above 4.* Theorem F decides `χ_c(F²)` below 4: there it is the least value of `χ_c(F_v²)` over the completions
   (Corollary F8), given by a locally constant colouring of one completion. Is that true above 4 as well, for `χ_c`
   or at least for `χ`? Theorem W⁺ is not available there (a homomorphism to `K_{p/q}` with `p/q ≥ 4` need not come
   from a character). The first case is `ℚ(√47)`: `4 ≤ χ_c(ℚ(√47)²) ≤ 19/4` and `4 ≤ χ ≤ 5`, and
   `χ_c(ℚ(√47)²) = 4` if and only if `χ(ℚ(√47)²) = 4` (Moorhouse's open case; for "only if", `χ_c = 4` gives
   `χ_c ≤ 4` on every finite subgraph, hence `χ ≤ 4` there, and `χ ≤ 4` by de Bruijn–Erdős). In particular: is
   there a number field with `4 < χ_c(F²) < 5`? Any `F` with `χ(F²) = 5` and a place of residue field `𝔽₁₉` would
   be one. The 5-chromatic unit-distance graphs we know lie over fields containing `√3`, or `√−3` and `√−7`, and 3
   and −7 are not squares modulo 19.
2. *Finite witnesses.* By Corollary F12 some finite unit-distance graph in `ℚ(√11)²` has `χ_c = 7/2`; the proof
   (compactness) gives no bound on its size, and the smallest we know has 155 vertices and is vertex-critical (§6.8). How small can it be? (For the value 3 in a plane without unit triangles nine vertices are needed, and they suffice over `ℚ(√7)` and
`ℚ(√31)`; §6.8. Over other such planes the least number is open; over `ℚ(√15)` the smallest witness we found has 13
vertices (`witness_q15.json.gz`), and none of the three nine-vertex graphs appeared with the denominators we tried.) When `χ_c(F²) = 4` the value is
   attained as well (Theorem F16, for example over `ℚ(√59)`), but the proof gives no bound on the size of a witness,
   and we know none: how small can one be over `ℚ(√59)`? By the converse in §6.9 its edge vectors `U` must satisfy
   `κ(U) ≤ 1/4`, which excludes those with denominator 210 (`at_four/theta59_210.json`).
3. *The one-prime probe.* Does `r0(5^k)` tend to `3/10` (the conjecture of §5)? Theorem E no longer needs it.

Measurements (floating-point MIP, not proofs): for `U = G_25{1, u_n, ū_n : n = 1, 7, 19}`, `max_ξ min_u ‖ξ(u)‖` is
`2/7` for `d = 11`, `35` and `23` (SCIP on an LLL-reduced relation basis, optimal; an exact certificate for `d = 23`
was started and stopped once Theorem F gave `χ_c(ℚ(√23)²) = 7/2`), `0.2638` for `d = 59` (an optimal character with
irrational-looking values, so this `U` does not reach `χ_c ≥ 4`; the exact certificate at `53/15` is Proposition C6,
and Theorem F gives `χ_c = 4`), and `0.2958` for `d = 71` (7 splits; this `U` alone does not give `7/2`, but Theorem F
does). Without the LLL step the solver had failed on `d = 23`, `47` (no solution in 300 s) and reported a wrong value
elsewhere; those failures were numerical.
