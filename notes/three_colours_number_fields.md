# Three colours for the plane over every number field

*Working note, 4 October 2026. Theorem B extends Theorem 1 of `notes/four_colours_11_mod_12.md` from real quadratic
fields to every number field, with the same two tools: Theorem W (`notes/winding_lemma.md`) and the shape of the
characters along the rational rotations of denominator `5^k` (Proposition 1 and Lemma 5 there). Two internal referees
(separate AI agents, with their own programs) checked the proof and found no error; their corrections are applied
(research log, 4 October). Nobody outside the project has checked it. Proposition B9 (§8) is an elementary special
case with explicit vectors, and `χ(ℚ(√2, √7)²) = 4` (§7) has an exact certificate that does not depend on §§3–4.*

## 1. The statement

Let `F` be a number field and give `F²` the unit-distance relation `(x − x′)² + (y − y′)² = 1`.

> **Theorem B.** `χ(F²) ≤ 3` if and only if
> (a) some prime of `F` above 2 ramifies in `F(i)`, or
> (b) some prime of `F` above 3 has residue degree 1.
> Moreover `χ(F²) ≤ 2` if and only if (a) holds (Theorem A of `notes/local_global.md`).

So `χ(F²) ≥ 4` exactly when no prime above 2 ramifies in `F(i)` and every prime above 3 has residue degree at least
2. Whether this happens depends only on `F ⊗ ℚ₂` and `F ⊗ ℚ₃`; `data/number_fields/three_colours/classify.gp`
decides it from a defining polynomial (relative discriminant of `F(i)/F`, residue degrees above 3).

The "if" half is known in substance: (a) gives a locally constant 2-colouring at a place above 2 (Theorem A′ of
hn-2adic-obstruction; `notes/local_global.md`, Theorem A, Step 2), and (b) gives the 3-colouring of the residue plane
`𝔽₉ ⊃ μ₄` at a place above 3 (Madore, Cor. 3.4; Moorhouse, Corollary 8.3; Fischer 1990, Theorem 9, for quadratic
fields; the general argument is written out in `notes/quadratic_planes.md` §5). The new half is "only if".

**Consequences.**
- *Real quadratic fields.* For squarefree `d`, (a) holds iff `d ≢ 3 (mod 4)`, and (b) iff `d ≢ 2 (mod 3)`. So
  Theorem B gives Theorem 1 and the first sentence of Corollary 2 of `notes/four_colours_11_mod_12.md`
  (`χ(ℚ(√d)²) ≥ 4` iff `d ≡ 11 (mod 12)`); the upper bounds there come from other places.
- *`LG₃` for every number field* (`notes/local_global.md`, Question LG): `F²` contains a finite unit-distance graph
  with no proper 3-colouring unless some place of `F` has a locally constant 3-colouring (and the place can then be
  taken above 2 or 3). By de Bruijn–Erdős this is the same as: `χ(F²) ≤ 3` iff some place has a locally constant
  3-colouring.
- *Local fields* (Corollary B8). For a finite extension `K` of `ℚ_p`: `χ(K²) = 2` if `p = 2` and `K(i)/K` is
  ramified; `χ(K²) = 3` if `p = 3` and `K` has residue degree 1; otherwise `χ(K²) ≥ 4` (and `χ(K²) = ∞` when
  `i ∈ K`, by Davies's theorem). This extends Corollary 3 (`p ≥ 5`) and Madore's `χ(ℚ₂²) = 2`, `χ(ℚ₃²) = 3`. For
  example the plane over the unramified cubic extension `ℚ₂₇` of `ℚ₃` needs at least four colours, although `i ∉ ℚ₂₇`
  and `x² + y²` is anisotropic over it, as over `ℚ₃`.
- *Fields not covered before.* Two families were known to need four colours: fields containing `√d` with
  `d ≡ 11 (mod 12)` (Theorem 1), and fields containing `√3` and `√q` with `q ≡ 2 (mod 3)` (chains of rhombi,
  `notes/local_colourings.md` §11, Theorem 5). Theorem B adds, for example, `ℚ(√2, √7)`, `ℚ(√2, √31)` and
  `ℚ(√2, √55)` (no unit triangle, and every quadratic subfield has `χ ≤ 3`), and `ℚ(2cos(2π/7), √7)` of degree 6,
  where 3 has residue degree 3 and the obstruction is the one of Lemma B5. For `ℚ(√2, √7)` and the sextic field a
  place above 7 with residue field `𝔽₇` gives the upper bound 4 (§7), so `χ = 4` for both.

## 2. From colourings to a linear problem over `ℚ(i)`

If `i ∈ F`, then `χ(F²) = ∞` (Davies's theorem, through `(x, y) ↦ (x, iy)`; `notes/quadratic_planes.md` §5), (a)
and (b) fail (`F(i) = F`, and every prime above 3 has even residue degree, since `𝔽₉` lies in its residue field), and
Theorem B holds. From now on `i ∉ F`. Put `L = F(i)` and identify `F²` with `L` by `(x, y) ↦ x + iy`; then
`x² + y² = N(z) = z z̄`, the unit vectors form the group `T = {z ∈ L : z z̄ = 1}`, and the unit-distance graph is
`Cay(L, T)`. As a `ℚ(i)`-vector space `L` has dimension `n = [F : ℚ]`, and `Tr = Tr_{L/ℚ(i)}` gives a non-degenerate
pairing.

Let `G_N` (`N = 5^k`, `k ≥ 1`) be the rational unit vectors with denominator dividing `N`, `I = [1/3, 2/3]`,
`S_N = {c ∈ ℂ : Re(c̄γ) ∈ I + ℤ for all γ ∈ G_N}`, `r = √10/18`, and

    E = E_c ∪ E_q,   E_c = (1 + i)/2 + ℤ[i],   E_q = {(e₁ + e₂ i)/3 : e₁, e₂ ∈ {1, 2}} + ℤ[i]   (types c and q).

By Proposition 1 and Lemma 5 of `notes/four_colours_11_mod_12.md`, every `c ∈ S_N` is `Nε + x` with `ε ∈ E` and
`|x| ≤ r`, and if `σ_j = Nε_j + x_j ∈ S_N` and Gaussian integers `μ_j` satisfy `Σ μ_j σ_j = 0` and
`N > 6r Σ|μ_j|`, then `Σ μ_j ε_j = 0`. (Lemma 5 is stated there with rational integers; the proof is the same:
`Σ μ_j N ε_j` lies in the `ℤ[i]`-module `(N/6)ℤ[i]` and has modulus at most `r Σ|μ_j| < N/6`.) Only a uniform bound
on `|x|` matters below, so the cruder form of Proposition 1 proved in Lean (`|x|² < 1/18`) would do as well.
`S_N = conj(S_N)`, because `G_N` is closed under conjugation.

> **Lemma B1.** If `F²` is 3-colourable, then for every finite `V ⊂ T` there is `β ∈ L` with `Tr(βv) ∈ E` for
> every `v ∈ V`.

*Proof.* Enlarge `V` to contain a `ℚ(i)`-basis `b_1, …, b_n` of `L` contained in `T`. Such a basis exists: the
`ℚ(i)`-span `W` of `T` is a subring of `L` (as `T` is a group), finite-dimensional over `ℚ(i)` and contained in the
field `L`, hence itself a field; it contains `z = (p + i)/(p − i)` for every `p ∈ F`, hence `(z − 1)⁻¹` and
`p = i(z + 1)/(z − 1)`; so `W = L`. Write `D_v v = Σ_j μ_{v,j} b_j` with `D_v ∈ ℤ ∖ {0}`, `μ_{v,j} ∈ ℤ[i]`, and take
`N = 5^k` with `N > 6r max_v (|D_v| + Σ_j |μ_{v,j}|)` (so `k ≥ 1`). The set `U = G_N V` is finite, symmetric
(`−1 ∈ G_N`) and consists of unit vectors, and `Cay(ℤU, U)` is a subgraph of `Cay(L, T)`, so it is 3-colourable; by
Theorem W there is a character `ξ` of `ℤU` with `ξ(U) ⊆ I + ℤ`. Extend `ξ` to `L` (the circle is divisible) and
restrict it to a lattice `Λ = (1/M) ⊕_j ℤ[i] f_j ⊇ U`, where `f_1, …, f_n` is a `ℚ`-basis of `F` (so
`L = ⊕_j ℚ(i) f_j`). On `Λ ≅ ℤ^{2n}` every character lifts to a real-linear form, and every real-linear form on `ℂ^n`
is `Re` of a `ℂ`-linear one: there is a `ℚ(i)`-linear `φ : L → ℂ` with `ξ(w) = Re φ(w) (mod 1)` for `w ∈ Λ`. For
`v ∈ V` and `γ ∈ G_N`, `ξ(γv) = Re(γ φ(v))`, so `φ(v) ∈ conj(S_N) = S_N`. Write `φ(v) = Nε_v + x_v`
(Proposition 1). Applying `φ` to `D_v v − Σ_j μ_{v,j} b_j = 0` and Lemma 5 gives `D_v ε_v = Σ_j μ_{v,j} ε_{b_j}`. Let
`β ∈ L` be the element with `Tr(β b_j) = ε_{b_j}` for all `j`; by `ℚ(i)`-linearity `Tr(βv) = ε_v ∈ E` for every
`v ∈ V`. ∎

Only two properties of `E` matter below: for `ε ∈ E_c`, `v_{1+i}(ε) = −1` and `ε` is integral at 3; for `ε ∈ E_q`,
`ε` is integral at `1 + i` and `3ε ≡ ±1 ± i (mod 3)`, i.e. `3ε mod 3` lies in the set `A′ = {±1 ± i}` of
non-squares of `𝔽₉ = ℤ[i]/3`.

Lemma B1 can also be used on its own, by computer: if for some finite `V` no `β` exists, then `G_N V` is not
3-colourable for large `N`, and `χ(F²) ≥ 4`. The test is one-sided (a `β` for one `V` proves nothing); §7 uses it to
confirm the examples independently of §§3–4.

## 3. Localisation at 2 and 3

For a prime `p` write `L_p = L ⊗_ℚ ℚ_p = ∏_{w | p} L_w` and `T_p = {z ∈ L_p : z z̄ = 1} = ∏_{u | p} T(F_u)`, the
product over the places `u` of `F` above `p`. Since 2 ramifies and 3 is inert in `ℚ(i)`, `ℚ(i) ⊗ ℚ₂ = ℚ₂(i)` and
`ℚ(i) ⊗ ℚ₃ = ℚ₉`, and `Tr` extends to `Tr : L_2 → ℚ₂(i)` and `Tr : L_3 → ℚ₉`. Write `O₂ = ℤ₂[i]`, `O₃ = ℤ₉`, and
`𝒜 = {(±1 ± i)/3} + ℤ₉`.

*Weak approximation for `T`.* For every finite set `S` of places, `T` is dense in `∏_{u ∈ S} T(F_u)`. Indeed
`p ↦ (p + i)/(p − i)` is a homeomorphism from `F_u` minus the roots of `p² + 1` onto `T(F_u) ∖ {1}` (inverse
`z ↦ i(z + 1)/(z − 1)`; at a place where `i ∈ F_u`, so that `L_u ≅ F_u × F_u` and `T(F_u) = {(s, 1/s)}`, it is
`p ↦ (s, 1/s)` with `s = (p + ι)/(p − ι)`, `ι² = −1`), it maps `F` into `T`, `F` is dense in `∏_{u ∈ S} F_u`, and
`T(F_u)` has no isolated points.

> **Lemma B2.** If `F²` is 3-colourable, then one of the following holds:
> (c) there is `β₂ ∈ L₂` with `v_{1+i}(Tr(β₂ z)) = −1` for every `z ∈ T₂`;
> (q) there is `β₃ ∈ L₃` with `Tr(β₃ z) ∈ 𝒜` for every `z ∈ T₃`.

*Proof.* Fix a `ℚ(i)`-basis `b_1, …, b_n ⊂ T` of `L` and let `K₂ = {β ∈ L₂ : Tr(βb_j) ∈ (1 + i)^{−1}O₂ ∀ j}` and
`K₃ = {β ∈ L₃ : Tr(βb_j) ∈ (1/3)O₃ ∀ j}`; they are compact (`β ↦ (Tr(βb_j))_j` is a linear isomorphism onto
`ℚ₂(i)^n`, resp. `ℚ₉^n`). On the compact space `Ω = {c, q}^T × K₂ × K₃` consider, for `v ∈ T`, the condition
`P_v(τ, β₂, β₃)`:

    τ(v) = c, v_{1+i}(Tr(β₂v)) = −1 and Tr(β₃v) ∈ O₃;   or   τ(v) = q, Tr(β₂v) ∈ O₂ and Tr(β₃v) ∈ 𝒜.

Each `P_v` defines a clopen subset of `Ω`. For a finite `V ⊇ {b_j}`, Lemma B1 gives `β ∈ L`; its images in `L₂` and
`L₃` lie in `K₂` and `K₃` (because `E ⊂ (1 + i)^{−1}O₂` and `E ⊂ (1/3)O₃`), and with `τ(v)` the type of `Tr(βv)` the
triple satisfies `P_v` for every `v ∈ V` (by the two properties of `E` at the end of §2). So the closed sets
`{P_v for all v ∈ V}` have the finite intersection property, and some `(τ, β₂, β₃)` satisfies `P_v` for all `v ∈ T`.

Let `a(z)` be the statement `v_{1+i}(Tr(β₂z)) = −1` and `a′(z)` the statement `Tr(β₂z) ∈ O₂` (`z ∈ T₂`), and `b(y)`:
`Tr(β₃y) ∈ O₃`, `b′(y)`: `Tr(β₃y) ∈ 𝒜` (`y ∈ T₃`). These are locally constant, `a, a′` exclude each other, and so do
`b, b′`. The set `Z = {(z, y) : (a(z) ∧ b(y)) ∨ (a′(z) ∧ b′(y))}` is clopen in `T₂ × T₃` and contains the diagonal
image of `T`, which is dense (weak approximation); so `Z = T₂ × T₃`. If `a(z₀)` for some `z₀`, then `b(y)` for every
`y`, so no `y` satisfies `b′`, so `a(z)` for every `z`: this is (c). Otherwise `a′` holds everywhere, hence `b′`
holds everywhere: this is (q). ∎

**Case (c).** For `z ∈ L`, put `g(z) = (1 + i)Tr(β₂z)`; fix a representative `ρ` of each class of `ℚ₂(i)` modulo `O₂`
and colour `z` by `(g(z) − ρ) mod (1 + i)`, where `ρ` represents the class of `g(z)`. A unit step `v` adds `g(v)`, a
unit of `O₂`, which keeps the class modulo `O₂` and changes the residue modulo `1 + i`: a proper 2-colouring of
`L ≅ F²`. So `χ(F²) ≤ 2`, and by Theorem A of `notes/local_global.md` (its Step 3) some prime above 2 ramifies in
`F(i)`: (a) holds.

## 4. Case (q): the places above 3

Let `u` be a place of `F` above 3, `f = f_u` its residue degree, `q = 3^f`, `L_u = L ⊗_F F_u`, and
`g_u(z) = Tr_{L_u/ℚ₉}(β_u z)`, where `β₃ = (β_u)_u`. Then `Tr(β₃y) = Σ_u g_u(y_u)` for `y = (y_u) ∈ ∏_u T(F_u)`.

*Even `f`.* Then `𝔽₉` lies in the residue field, `i ∈ F_u`, `L_u ≅ F_u × F_u` and `T(F_u) = {(s, s^{−1}) : s ∈ F_u^×}`.
Write `g_u(s, t) = ℓ₁(s) + ℓ₂(t)` with `ℚ₃`-linear `ℓ₁, ℓ₂ : F_u → ℚ₉`. In (q), `g_u(y_u)` stays in a bounded set
while `y_u` varies alone (the sum lies in `𝒜` and the other terms are fixed). If `ℓ₁(s₀) ≠ 0`, then
`g_u(3^{−m}s₀, 3^m s₀^{−1}) = 3^{−m}ℓ₁(s₀) + 3^mℓ₂(s₀^{−1})` is unbounded as `m → ∞`; so `ℓ₁ = 0`, likewise
`ℓ₂ = 0`, and `g_u = 0`.

*Odd `f`.* Then `L_u = F_u(i) = F_u ℚ₉` is a field, unramified of degree 2 over `F_u`, with residue field `𝔽_{q²}`;
its non-trivial automorphism `σ` acts on `ℚ₉` as the Frobenius (`σ(i) = −i`); `T(F_u) ⊆ O_{L_u}^×` contains
`μ₄ ⊂ ℚ₉` and the Teichmüller lifts of `μ_{q+1} ⊂ 𝔽_{q²}^×` (`N(ζ) = ζ ζ^q = 1`), and reduces onto `μ_{q+1}`.
`g_u` is `ℚ₉`-linear, so `g_u(εz) = εg_u(z)` for `ε ∈ μ₄`.

> **Lemma B3 (one place).** In case (q) there is a single place `u₀ | 3`, with odd residue degree, and
> `β ∈ L_{u₀}` with `Tr_{L_{u₀}/ℚ₉}(βz) ∈ 𝒜` for every `z ∈ T(F_{u₀})`.

*Proof.* Only places with odd `f` contribute. For such `u`, and `z, z′ ∈ T(F_u)`, varying `y_u` alone in (q) gives
`g_u(z) − g_u(z′) ∈ 𝒜 − 𝒜 ⊆ (1/3)O₃`. With `z′ = 1` and `z = i`: `(i − 1)g_u(1) ∈ (1/3)O₃`, and `i − 1` is a unit
at 3, so `g_u(1) ∈ (1/3)O₃` and `g_u(T(F_u)) ⊆ (1/3)O₃`. Let `Z_u ⊆ 𝔽₉` be the set of residues `3g_u(z) mod 3`; it is
stable under `μ₄`, and (q) says `Σ_u z_u ∈ A′` for every choice of `z_u ∈ Z_u`. If two places had non-zero elements
`w ∈ Z_u`, `w′ ∈ Z_{u′}`, fix `z_{u″} ∈ Z_{u″}` at the other places and let `c = ε′w′ + Σ z_{u″}`; then the four
distinct elements `c + εw` (`ε ∈ μ₄`) lie in `A′`, which has four elements, so they are `A′`, and summing,
`4c = c = Σ A′ = 0` (in `𝔽₉`, `Σ_{ε ∈ μ₄} ε = 0`). This holds for every `ε′ ∈ μ₄`, so `w′ = iw′`, `w′ = 0`: a
contradiction. Since `0 ∉ A′`, exactly one place `u₀` has `Z_{u₀} ≠ {0}`; the other `g_u` take values in `O₃`, and
`β = β_{u₀}` works. ∎

> **Lemma B4 (level 1).** Let `u | 3` have odd residue degree `f`, and `β ∈ L_u` with `Tr_{L_u/ℚ₉}(βz) ∈ 𝒜` for all
> `z ∈ T(F_u)`. Then `λ(x) = 3Tr(βx) mod 3` is defined on `O_{L_u}`, vanishes on the maximal ideal `𝔪`, and the
> induced `𝔽₉`-linear map `λ₁ : 𝔽_{q²} → 𝔽₉` satisfies `λ₁(μ_{q+1}) ⊆ A′`.

*Proof.* The closed additive span `R` of `T(F_u)` is a closed subring of `O_{L_u}`. It contains `ℤ₃[μ_{q+1}]`, the
ring of integers `W` of the unramified extension of degree `2f` (the order of 3 modulo `q + 1` is `2f`). For a
uniformiser `π` of `F_u`, `z = (π + i)/(π − i) ∈ T(F_u)` and `σ(z) ∈ T(F_u)` give `2 − z − σ(z) = 4/(π² + 1) ∈ R`, a
unit of `O_{L_u}`; its inverse is a limit of its powers, so `(π² + 1)/4 ∈ R`, and
`π = (z − σ(z))(π² + 1)/(4i) ∈ R`. Hence `R ⊇ W[π] = O_{L_u}` (`π` is a uniformiser of `L_u`, as `L_u/F_u` is
unramified). By continuity `Tr(βO_{L_u}) ⊆ (1/3)O₃`, so `λ : O_{L_u} → 𝔽₉` is defined, `ℤ₉`-linear, and zero on
`3O_{L_u} = 𝔪^e`.

Suppose `λ(𝔪) ≠ 0` and let `j ≥ 1` be largest with `λ(𝔪^j) ≠ 0`; identify `𝔪^j/𝔪^{j+1}` with `𝔽_{q²}` by
`x ↦ (x/π^j) mod 𝔪`, so that `λ` induces a non-zero `𝔽₃`-linear map `λ_j : 𝔽_{q²} → 𝔽₉`. For `η ∈ O_{F_u}`,
`t = (1 + iπ^jη)/(1 − iπ^jη)` lies in `T(F_u)` (`σ(t) = 1/t`, as `σ(i) = −i` and `σ(π) = π`) and
`t ≡ 1 + 2iπ^jη (mod 𝔪^{j+1})`. For `z ∈ T(F_u)` with residue `ẑ ∈ μ_{q+1}`, `λ(zt) = λ(z) + λ_j(2i ẑ η̂)`, where `η̂`
is the residue of `η`. As `η̂` runs over `𝔽_q`, the second term runs over an `𝔽₃`-subspace `Σ_z` of `𝔽₉`, and
`λ(z) + Σ_z ⊆ A′`. `A′ = {a + bi : a, b ≠ 0}` contains no affine `𝔽₃`-line, so `Σ_z = 0`: `λ_j` vanishes on
`i 𝔽_q ẑ` for every `ẑ ∈ μ_{q+1}`. The `𝔽_q`-span of `μ_{q+1}` is `𝔽_{q²}` (`μ_{q+1} ⊄ 𝔽_q`), so `λ_j = 0`, a
contradiction. Hence `λ(𝔪) = 0`; `λ₁` is `𝔽₉`-linear, and `λ₁(μ_{q+1}) ⊆ A′` because `T(F_u)` reduces onto
`μ_{q+1}`. ∎

> **Lemma B5 (residue fields).** Let `f` be odd and `q = 3^f`. There is an `𝔽₉`-linear `λ₁ : 𝔽_{q²} → 𝔽₉` with
> `λ₁(μ_{q+1}) ⊆ A′` if and only if `f = 1`.

*Proof.* `A′` is the set of `y ∈ 𝔽₉` with `y⁴ = −1` (the non-squares). For `f = 1`, `μ_{q+1} = μ₄ ⊂ 𝔽₉` and
`λ₁(x) = (1 + i)x` works. Let `f ≥ 3`. Every `λ₁` is `x ↦ Tr_{𝔽_{q²}/𝔽₉}(bx) = Σ_{j<f} (bx)^{9^j}`. In
characteristic 3, `λ₁(z)⁴ = λ₁(z)³λ₁(z) = Σ_{j,k<f} b^{e_{jk}} z^{e_{jk}}` with `e_{jk} = 3·9^j + 9^k`, so if
`λ₁(μ_{q+1}) ⊆ A′`, the function `Φ(z) = 1 + Σ_{j,k} b^{e_{jk}} z^{e_{jk}}` vanishes on `μ_{q+1}`. On the cyclic
group `μ_{q+1}` (of order prime to 3) the characters `z ↦ z^e`, `e mod q + 1`, are linearly independent, so for every
class `e mod q + 1` the sum of the coefficients of the `z^{e_{jk}}` with `e_{jk} ≡ e` must vanish (or equal `−1` for
`e ≡ 0`). Modulo `q + 1`, `g = −3` satisfies `g^f = −3^f ≡ 1` (`f` odd), `9^j ≡ g^{2j}` and `3·9^j ≡ −g^{2j+1}`.
Take `e = e_{jj} = 4·9^j ≡ 4g^{2j}`. If `e_{j′k′} ≡ e_{jj}`, then multiplying by `g^{−2j}`, `4 ≡ −g^α + g^γ` with
`α ≡ 2j′ + 1 − 2j`, `γ ≡ 2k′ − 2j (mod f)`, `0 ≤ α, γ < f`; as integers `g^α = (−3)^α` and `g^γ = (−3)^γ` have
absolute value at most `3^{f−1}`, so `|4 + (−3)^α − (−3)^γ| ≤ 4 + 2·3^{f−1} < 3^f + 1` for `f ≥ 3`, and the
congruence is an equality `(−3)^γ − (−3)^α = 4`. It forces `γ = 0` (otherwise 3 divides the left side, or
`(−3)^γ = 5`) and then `(−3)^α = −3`, `α = 1`; so `j′ = j` and `k′ = j` (2 is invertible modulo `f`). Hence the class
of `4·9^j` contains only `e_{jj}`, it is not 0 (`q + 1 > 4` and `gcd(q + 1, 3) = 1`), and its coefficient is
`b^{4·9^j}`, which must vanish: `b = 0`. But then `Φ = 1 ≠ 0`. ∎

Exhaustive searches confirm the lemma for `f = 3` and `f = 5` (`residue_lemma.py`: no `b ∈ 𝔽₇₂₉`, resp. `𝔽₅₉₀₄₉`,
works), and the referees checked the uniqueness of the class of `4·9^j` for every odd `f` up to 79.

## 5. Proof of Theorem B, and two corollaries

*"Only if".* Let `F²` be 3-colourable, `i ∉ F`. By Lemma B2, (c) or (q) holds. In case (c), §3 gives (a). In case (q),
Lemmas B3 and B4 give a place `u₀ | 3` of odd residue degree `f` and an `𝔽₉`-linear `λ₁ : 𝔽_{q²} → 𝔽₉` with
`λ₁(μ_{q+1}) ⊆ A′`, and Lemma B5 gives `f = 1`: (b).

*"If".* (a) gives a 2-colouring (`notes/local_global.md`, Theorem A, Step 2). Under (b), let `u | 3` have residue
degree 1 and `w` the place of `L` above it. Then `i ∉ F_u` (−1 is not a square in `𝔽₃`), `L_w = F_u(i)` has residue
field `𝔽₉`, and unit vectors of `F_u²` are units of `O_w` whose residues lie in `μ₄`. Unit steps preserve the classes
of `L_w` modulo `O_w`; fix a representative `ρ` of each class and colour `z` by `ℓ((z − ρ) mod 𝔪_w)`, where
`ℓ(a + bi) = a + b` on `𝔽₉ = 𝔽₃ + 𝔽₃i`. Since `ℓ(±1) = ±1` and `ℓ(±i) = ±1`, every unit step changes the colour. This
colours `F_u² ⊇ F²` with three colours. Finally `χ(F²) ≤ 2` iff (a), by Theorem A. ∎

The proof uses Theorem W only through Lemma B1, for the finite sets `G_N V`; the rest is compactness, weak
approximation for the norm-one torus, the arithmetic of the places above 2 and 3, and, in case (c), Theorem A. The
rational rotations of denominator `5^k` act as a probe: Proposition 1 says that along them a character into
`[1/3, 2/3]` is, up to a bounded error, either the 2-adic half-character (type c) or a 3-adic level-1 character
(type q), and Lemma B2 turns this into the dichotomy between the places above 2 and those above 3.

> **Corollary B6 (no circular chromatic number between 2 and 3).** For every number field `F`, `F²` has a
> homomorphism to a circular clique `K_{p/q}` with `p/q < 3` only if `χ(F²) ≤ 2`. So `χ_c(F²) = 2`,
> `χ_c(F²) = 3` or `χ_c(F²) > 3`, and `F²` maps to an odd cycle only when it is bipartite.

*Proof.* Theorem W⁺ (`notes/winding_lemma.md` §2) and compactness give a character with `ξ(T) ⊆ [q/p, 1 − q/p]`, an
interval strictly inside `(1/3, 2/3)`. In the proof of Lemma B1, `φ(v)` then lies in the subset of `S_N` defined by
this interval, which contains no point of type q (at `γ = 1` their value `Re(c̄)` lies in `{1/3, 2/3} + ℤ`). So every `v` has type c, `τ ≡ c` in Lemma B2, case (c) holds, and `χ(F²) ≤ 2`. ∎

> **Corollary B7 (admissible fields).** If `i ∈ F_v` for every place `v` of `F` above 2 and above 3, then
> `χ(F²) ≥ 4`.

This covers the admissible fields of `notes/local_global.md` §4. Of those listed there, `ℚ(√167)`, `ℚ(√2, √47)`,
`ℚ(√11, √13)` and `ℚ(√10, √38)` (which contains `√95`) already contain some `√d` with `d ≡ 11 (mod 12)`; for
`ℚ(√2, √31)` and `ℚ(√2, √55)` the bound is new (also by Proposition B9).

## 6. Local fields

> **Corollary B8.** Let `K` be a finite extension of `ℚ_p`. Then `χ(K²) = 2` if `p = 2` and `K(i)/K` is ramified,
> `χ(K²) = 3` if `p = 3` and the residue degree of `K` is 1, and `χ(K²) ≥ 4` in every other case.

*Proof.* The colourings in the first two cases are those of §5 (they colour `K²` itself), and `χ(K²) ≥ 3` when
`p = 3` because `ℚ(√7) ⊂ ℚ₃ ⊆ K` (`7 ≡ 1 (mod 3)`) and `χ(ℚ(√7)²) = 3`. If `p ≥ 5`, `χ(K²) ≥ χ(ℚ_p²) ≥ 4`
(Corollary 3). If `i ∈ K`, `χ(K²) = ∞` (Davies). There remain `p = 2` with `i ∉ K` and `K(i)/K` unramified, and
`p = 3` with odd residue degree at least 3. In both we find a number field `F ⊂ K` for which (a) and (b) fail; then
`χ(K²) ≥ χ(F²) ≥ 4`.

Let `m = [K : ℚ_p]` and `p′` the other prime of `{2, 3}`. If `p = 2`, `m` is even: if `m` were odd, `ℚ₂^×/ℚ₂^{×2}`
would inject into `K^×/K^{×2}`, and as `−1` and `−5` are non-squares in `ℚ₂`, `i ∉ K` and `K(i) ≠ K(√5)`, the
unramified quadratic extension, so `K(i)/K` would be ramified. Choose monic polynomials over `ℚ_p`, `ℚ_{p′}` and `ℚ₅`
of a common even degree `n`:
- over `ℚ_p`: `h_K h′`, with `h_K` the minimal polynomial of a generator of `K` and `h′` a product of distinct
  irreducible quadratics whose roots generate `ℚ₂(i)` (`p = 2`) or `ℚ₉` (`p = 3`), coprime to `h_K`; if `p = 3` and
  `m` is odd, `h′` also contains one cubic factor generating `ℚ₂₇` (coprime to `h_K`, e.g. `x³ − x + 1` or
  `x³ − x − 1`), so `n = m + 3 + 2k`;
- over `ℚ_{p′}`: a product of distinct irreducible quadratics whose roots generate `ℚ₂(i)`, resp. `ℚ₉`;
- over `ℚ₅`: an irreducible polynomial.

By weak approximation choose a monic `h ∈ ℚ[x]` whose coefficients are close enough to all three at once that, by
Krasner's lemma (each root of `h` closer to a root of the target than the target's roots are to one another), `h`
factors over each of `ℚ_p`, `ℚ_{p′}`, `ℚ₅` with the same factor degrees and the same fields as the target. Then
`F = ℚ[x]/(h)` is a field (irreducible over `ℚ₅`), it embeds in `K` (a root of `h` close to the generator of `K`
generates `K`), its places above `p` are the place giving `K` and places containing `i` (or, for the cubic factor,
of residue degree 3), and its places above `p′` all contain `i`. So no place above 2 ramifies in `F(i)` and no place
above 3 has residue degree 1. ∎

For example `χ(ℚ₂₇²) ≥ 4`, and `χ(ℚ₂(√3)²) ≥ 4` (`ℚ₂(√3)(i) = ℚ₂(√3, √−3)` is unramified over `ℚ₂(√3)`), the latter
also from the 76-vertex graph of `ℚ(√11)` (`ℚ₂(√11) = ℚ₂(√3)`). The first referee ran the construction in PARI for
`K = ℚ₂(√3)` and `K = ℚ₂₇` (a quartic and a sextic field `F` with the predicted places) and observed that the
precision must really be chosen as above: with too little 2-adic precision the factor pattern over `ℚ₂` changed.

## 7. Checks

The scripts are in `data/number_fields/three_colours/` (README there); the tests in `tests/test_three_colours.py`.

- **Lemma B1 by computer** (§2, one-sided). `limit_test.py` decides, for a finite set `V`, whether some `β` has
  `Tr(βv) ∈ E` for all `v ∈ V`: the possible values `z = 6φ(v)` form the lattice `(column space) ∩ ℤ^{2m}` (PARI
  `matrixqz(A, −2)`), and the types are conditions on `z` modulo 2 and 3, so it enumerates `L/2L × L/3L`. `check.gp`
  repeats the computation in PARI/GP alone, with its own unit-vector check and another lattice routine (`matkerint`
  twice); the first referee wrote a third version (lattice by `matsolvemod`), which agreed on all 78 of its trials.
  - The quadratic cases of Theorem 1 come out infeasible (`d = 11, 23, 35, 47, 59` with the vectors of §§4–5 of
    `notes/four_colours_11_mod_12.md`), and `ℚ(√3)`, `ℚ(√7)` feasible (`known_cases.py`).
  - `ℚ(√2, √7)`: infeasible with 27 vectors (`V_sqrt2_sqrt7.json`) in all three implementations, and already with five
    of them (`V_sqrt2_sqrt7_min5.json`, after greedy deletion): `1`, `(7 + 4√2 i)/9`, `(23 + 10√2 i)/27` (in `ℚ(√2)²`,
    large at a place above 3) and `(1 + 3√7 i)/8`, `(9 + 5√7 i)/16` (in `ℚ(√7)²`, large at a place above 2). Their
    relations `54v₂ = 11 + 45v₁` and `24v₄ = 11 + 20v₃` force the vector 1 to have type c at 3 and type q at 2.
  - `ℚ(2cos(2π/7), √7)`: infeasible with 44 vectors (`V_c7_sqrt7.json`), and with 30 in the second referee's run; the
    obstruction is at the two places above 3 of residue degree 3 (Lemma B5).
  - `ℚ(√3, √5)`: infeasible with 21 vectors, as the known chain of rhombi requires.
  - Controls: `ℚ(√3, √7)` (3 has residue degree 1) stays feasible with all vectors of type q, and
    `ℚ(2cos(2π/7), √2)` (2 ramifies in `F(i)`) with all of type c, for every set tried.
  - The first referee tested 23 random quartic and sextic fields (78 sets `V`, built from 1, approximations of `√−1`
    at the places above 2 and 3 where `i ∈ F_u`, and random unit vectors): every field predicted to be 2- or
    3-colourable stayed feasible; infeasibility was found for all 3 quartic and 4 of the 6 sextic fields predicted to
    need four colours (including places above 3 of type `(e, f) = (1, 3)` and `(2, 3)`); for the other two sextic
    fields only the bipartite pattern survived the random search, which needs an odd relation it did not find.
  The obstructions appear only once `V` contains vectors that are large at the split places above 2 and 3, as the
  proof of Lemma B2 predicts.
- **Exact certificates for `ℚ(√2, √7)`.** The 50 vectors `G_25{1, u_2, ū_2, u_7, ū_7}` of §8 (one per pair `±u`),
  and also the 70 vectors `G_125 V` for the five vectors above, admit no character into `[1/3, 2/3]`:
  `cert_two_roots_2_7_N25.json` (1 579 nodes, 1 054 leaves) and `cert_sqrt2_sqrt7_N125.json.gz` (11 357 nodes,
  8 754 leaves; branch and bound over the integer values of the relations, exact Farkas vectors at the leaves,
  `cert_general.py`). `check_mq.py`, an exact checker written separately, accepts both; it also recomputes that every
  vector has length 1 in `ℚ(√2, √7)`, and it rejects corrupted copies (a changed unit, relation, field, leaf or
  branch). With Theorem W (proved in Lean for every abelian group), this proves `χ(ℚ(√2, √7)²) ≥ 4` without
  Proposition 1, Lemma 5 or anything in §§3–4.
- **Upper bounds.** At a place `v | 7` with residue field `𝔽₇` (so `i ∉ F_v`), the unit vectors reduce into
  `μ₈ ⊂ 𝔽₄₉`, and `a + bi ↦ 2a + 3b` (`= Tr_{𝔽₄₉/𝔽₇}((1 + 2i)(a + bi))`) takes only the values 2, 3, 4, 5 on `μ₈`.
  Colouring each class of `L_w` modulo `O_w` by this map gives a homomorphism of `F²` to the circular clique
  `K_{7/2}`, so `χ(F²) ≤ 4` and `χ_c(F²) ≤ 7/2`. This applies to `ℚ(√2, √7)` (7 ramifies in `ℚ(√7)`, `√2 ∈ ℚ₇`) and
  to `ℚ(2cos(2π/7), √7)` (7 is totally ramified): both have `χ = 4`.
- **The lemmas.** `residue_lemma.py` checks Lemma B5 exhaustively for `f = 1, 3, 5`. The referees checked, with their
  own programs: Lemma B5 for `f = 3` by a different model of `𝔽₇₂₉` (all `3^12` `𝔽₃`-linear maps), and for
  `f ≤ 9` one `b` per coset of `μ_{q+1}`; the uniqueness step for odd `f ≤ 79`; Lemma B4 on the finite models
  `O_{L_u}/3 = 𝔽_{q²}[π]/(π^e)` (`f = 1`, `e ≤ 7`: exactly four admissible maps, all of level 1; `f = 3`, `e ≤ 2`:
  none); Lemma B3 over all `μ₄`-stable subsets of `𝔽₉`.
- **The criterion.** `classify.gp` computes (a) and (b) with PARI; for squarefree `|d| < 200` it reproduces
  `(a) ⟺ d ≢ 3 (mod 4)` and `(b) ⟺ d ≢ 2 (mod 3)` (second referee), and it agrees with the known values
  `χ(ℚ(√2, √3)²) = χ(ℚ(√3, √11)²) = 4`, `χ(ℚ(√3)²) = χ(ℚ(√7)²) = 3`, and `χ = 2` for odd degree.

## 8. An explicit family: two square roots

For many fields the general argument can be replaced by two explicit relations, as in Theorem 1a. For `c ≥ 1` put

    u_c = (1 + i√c)²/(1 + c) = ((1 − c) + 2√c i)/(1 + c),     so that     (1 + c)(u_c + ū_c) = 2(1 − c)·1.

> **Proposition B9.** Let `F ⊂ ℝ` contain `√a` and `√b`, where `a ≡ 2 (mod 3)` and `b ≡ 7 (mod 8)` are positive
> integers (`a = b` is allowed). If `N = 5^k > (4√10/3)·max(a, b)`, no character of the group generated by
> `U = G_N {1, u_a, ū_a, u_b, ū_b}` maps every vector of `U` into `[1/3, 2/3]`. So `χ(F²) ≥ 4`.

*Proof.* As in §2 of `notes/four_colours_11_mod_12.md`, a character with `ξ(U) ⊆ I + ℤ` gives numbers
`w_v ∈ S_N` (`v ∈ {1, u_a, ū_a, u_b, ū_b}`) satisfying the two relations exactly, and Lemma 5 (the sums of the absolute
values of the coefficients are `4a` and `4b`, and `6r·4c = (4√10/3)c < N`) gives, for the types,

    (1 + a)(ε_{u_a} + ε_{ū_a}) = 2(1 − a) ε_1,       (1 + b)(ε_{u_b} + ε_{ū_b}) = 2(1 − b) ε_1.

*At 3.* `3 | 1 + a` and `3 ∤ 2(1 − a)`. Every `ε` lies in `(1/6)ℤ[i]`, so the left side of the first relation has
`v₃ ≥ 0`, hence `v₃(ε_1) ≥ 0`: `ε_1` is not of type q.
*At 2.* `v₂(1 + b) ≥ 3` and `v₂(2(1 − b)) = 2` (as `1 − b ≡ 2 (mod 8)`). In `ℤ[i]`, the left side of the second
relation has `v_{1+i} ≥ 2·3 − 1 = 5`, while `v_{1+i}(2(1 − b)ε_1) = 4 + v_{1+i}(ε_1)`, which is 3 if `ε_1` has type c.
So `ε_1` is not of type c either, a contradiction. ∎

For `a = b = d` this is Theorem 1a (`d ≡ 23 (mod 24)`). The proposition covers, for example, `ℚ(√2, √7)`,
`ℚ(√2, √31)`, `ℚ(√5, √7)`, `ℚ(√2, √55)` and `ℚ(√3, √5)` (`b = 15`), and any field containing one of them; it is the
case "both places split" of Theorem B (`√a` makes every place above 3 have even residue degree, `√b` puts `i` in every
completion above 2). With only real parts, as in the formal proof of Theorem 1 (`lean/FourColours.lean`), the
crude form of Proposition 1 needs `N > 24 max(a, b)`.

*Check.* For `ℚ(√2, √7)` and `N = 25` (smaller than the proposition needs), the 50 vectors `G_25{1, u_2, ū_2, u_7, ū_7}`
(one per pair `±u`; `u_2 = (−1 + 2√2 i)/3`, `u_7 = (−3 + √7 i)/4`) admit no character into `[1/3, 2/3]`:
`cert_two_roots_2_7_N25.json` (1 579 nodes, 1 054 leaves) is accepted by `check_mq.py`.

## 9. References

- J. Davies, *Chromatic number of spacetime*, arXiv:2308.16885.
- K. G. Fischer, Additive K-colorable extensions of the rational plane, Discrete Math. 82 (1990) 181–195.
- D. A. Madore, The Hadwiger–Nelson problem over certain fields, arXiv:1509.07023 (2015), Cor. 3.4.
- G. E. Moorhouse, On the chromatic numbers of planes, draft of 3 March 2010, Corollary 8.3.
- MildlyMeticulous (GitHub account), hn-2adic-obstruction, July 2026, Theorem A′.
- Within the project: Theorem W and Theorem W⁺ (`notes/winding_lemma.md`; `lean/TheoremWInf.lean`), Proposition 1,
  Lemma 5, Theorem 1 and Corollary 3 (`notes/four_colours_11_mod_12.md`; `lean/FourColours.lean`,
  `lean/PadicFour.lean`), Theorem A (`notes/local_global.md`), Theorem 5 of `notes/local_colourings.md` §11.
