# Circular colourings of planes over number fields

*Working note, 4 October 2026. Theorem C is proved: the upper bound by an explicit colouring, the lower bound by an
exact certificate, checked by a separate exact checker, together with Theorem W⁺ (`notes/winding_lemma.md` §2,
proved in Lean in `lean/TheoremWInf.lean`). Nobody outside the project has checked it. The rest of the note is
questions and measurements.*

## 1. Circular colourings and characters

For integers `p ≥ 2q ≥ 2`, the circular clique `K_{p/q}` has vertices `ℤ/p`, with `i ~ j` when
`q ≤ (j − i mod p) ≤ p − q`, and the circular chromatic number `χ_c(G)` is the infimum of the `p/q` with
`G → K_{p/q}`; `χ(G) − 1 < χ_c(G) ≤ χ(G)`. For the plane, `χ_c(ℝ²) ≥ 4` (DeVos, Ebrahimi, Ghebleh, Goddyn, Mohar and
Naserasr, SIAM J. Discrete Math. 21 (2007)); an upper bound below 7 is in arXiv:1506.01886. We found no earlier
computation of `χ_c` for the plane over a number field or a `p`-adic field (one web search, 4 October 2026).

For a field `F` (with `i ∉ F`), `F²` is the Cayley graph `Cay(L, T)` of §2 of `notes/three_colours_number_fields.md`.
By Theorem W⁺, for `p/q < 4` a Cayley graph `Cay(Γ, S)` of an abelian group with finite `S = −S` maps to `K_{p/q}` if
and only if some character `ξ` has `ξ(S) ⊆ [q/p, 1 − q/p]`. Homomorphisms to a finite graph pass to the infinite graph
by compactness, so `χ_c(F²) ≤ p/q < 4` if and only if every finite set of unit vectors has such a character.

## 2. Upper bounds from residue planes

> **Proposition C1.** Let `v` be a place of `F` whose residue field is `𝔽_p`, `p ≡ 3 (mod 4)` (any ramification), and
> let `λ : 𝔽_{p²} → 𝔽_p` be `𝔽_p`-linear with `λ(μ_{p+1}) ⊆ {m, m + 1, …, p − m}`. Then `F²` maps to `K_{p/m}`, so
> `χ_c(F²) ≤ p/m`.

*Proof.* As `−1` is not a square in `𝔽_p`, `i ∉ F_v`, `L_w = F_v(i)` is unramified over `F_v` with residue field
`𝔽_{p²}`, and the unit vectors of `F_v²` are units of `O_w` whose residues lie in `μ_{p+1}`. Unit steps preserve the
classes of `L_w` modulo `O_w`; fix a representative `ρ` of each class and colour `z` by `λ((z − ρ) mod 𝔪_w) ∈ ℤ/p`. A
unit step changes the colour by an element of `{m, …, p − m}`: a homomorphism `F_v² → K_{p/m}`, which restricts to
`F²`. ∎

The best `m` at level 1 (`residue_circular.py`): `p = 3`: `m = 1` (`K₃`, three colours); `p = 7`: `m = 2`, with
`λ(a + bi) = 2a + 3b` taking the values 2, 3, 4, 5 on `μ₈` (so `K_{7/2}`); `p = 11`: `2/11`; `p = 19`: `4/19`;
`p = 23`: `3/23`; `p = 31`: `4/31`; `p = 43`: `5/43`. Only `p = 3` and `p = 7` give a bound below 4.

So every number field with a place of residue field `𝔽₇` has `χ_c(F²) ≤ 7/2` and `χ(F²) ≤ 4` (the second bound is
Moorhouse's Corollary 8.3, as `χ(𝔽₇-plane) = 4`). Examples: `ℚ(√11)` (7 splits), `ℚ(√35)` (7 ramifies), `ℚ(√2, √7)`,
`ℚ(2cos(2π/7), √7)`, and `ℚ₇` itself.

## 3. An exact value

> **Theorem C.** `χ_c(ℚ(√11)²) = χ_c(ℚ₇²) = 7/2`, while `χ(ℚ(√11)²) = χ(ℚ₇²) = 4`.

*Proof.* The upper bounds are Proposition C1 with `p = 7` (`11 ≡ 2² (mod 7)`, so `√11 ∈ ℚ₇`). For the lower bound let
`U` be the 70 unit vectors `G_25 V` (one per pair `±u`) with `V = {1, u_1, ū_1, u_7, ū_7, u_19, ū_19}`,
`u_n = (n + i√11)²/(n² + 11)`, and `G_25` the rational unit vectors with denominator dividing 25. The certificate
`data/number_fields/circular/cert_sqrt11_open_7_2_N25.json.gz` shows that no character of `ℤU` maps every vector of
`U` into the open interval `(2/7, 5/7)`: it branches on the integer values of 66 relations among the vectors
(17 744 nodes) and closes each of its 12 850 leaves by an exact Farkas vector; `check_open.py` checks it with integers
and fractions only (that every vector has length 1 in `ℚ(√11)`, the relations, that every branch covers all integers
strictly inside the range of its relation over the box, and every leaf), and rejects corrupted copies (a changed unit,
relation, field, interval, leaf or branch). If `F²` mapped to `K_{p/q}` with `p/q < 7/2`, Theorem W⁺ would give a
character of `ℤU` with values in `[q/p, 1 − q/p] ⊂ (2/7, 5/7)`. So `χ_c(ℚ(√11)²) ≥ 7/2`, and `ℚ(√11)² ⊂ ℚ₇²` gives
the same for `ℚ₇`. The chromatic numbers are Theorem 1 with Moorhouse's upper bound for `ℚ(√11)`, and
`notes/quadratic_planes.md` §5 (also in `lean/PadicPlanes.lean`) for `ℚ₇`. ∎

The certificate was found with a floating-point solver (`kappa_max.py`: SCIP, through OR-tools, maximises
`min_u ‖ξ(u)‖` and returns exactly `2/7`; `cert_open.py` builds the exact tree); nothing in the proof depends on the
solver.

## 4. What else is known

- **No value between 2 and 3.** By Corollary B6 of `notes/three_colours_number_fields.md`, `F²` maps to some
  `K_{p/q}` with `p/q < 3` only if it is bipartite; so `χ_c(F²) ∈ {2, 3}` or `χ_c(F²) > 3`, for every number field.
- **`ℚ(√35)`** (7 ramifies, `χ = 4`): the same 70-vector configuration gives `max_ξ min_u ‖ξ(u)‖ = 2/7` in the solver;
  an exact certificate is being built.
- **`ℚ(√23)`** (no place with residue field `𝔽₇`: 7 is inert): the solver did not finish in 300 s on the same
  configuration.

## 5. Questions

1. Is there a gap above 3: is `χ_c(F²) ≥ 7/2` whenever `χ(F²) ≥ 4`? Along the rational rotations with denominator
   `5^k` there are characters with `min ‖ξ(γ)‖ = 3/10` for every `k` we tried (`chartypes.py`, `lift.py`, up to
   `k = 6`); they look 5-adic, and 5-adic characters cannot survive on a whole field plane (every place above 5 is
   split in `F(i)`, which kills them as in the proof of Theorem B), but a proof needs an analogue of Proposition 1 for
   the intervals `[r, 1 − r]` with `2/7 < r < 1/3`.
2. What is `χ_c(F²)` for fields with `χ ≥ 4` and no place of residue field `𝔽₇`, such as `ℚ(√47)` (where
   `4 ≤ χ ≤ 5`)? A value below 4 there would give `χ(ℚ(√47)²) = 4`.
