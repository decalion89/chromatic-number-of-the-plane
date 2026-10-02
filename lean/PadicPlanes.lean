import Sqrt11

/-!
# The unit-distance graphs of the `p`-adic planes for `p = 7, 3, 2`

`QuadraticPlanes.sumSqGraph ℚ_[p]` is the graph on `ℚ_[p] × ℚ_[p]` in which `x ~ y` iff
`(x₁ - y₁)² + (x₂ - y₂)² = 1`. This file proves that its chromatic number is

* 4 for `p = 7` (`padicSeven_chromaticNumber`),
* 3 for `p = 3` (`padicThree_chromaticNumber`),
* 2 for `p = 2` (`padicTwo_chromaticNumber`).

*Upper bounds.* `O p = {x : ‖x‖ ≤ 1}` is a valuation subring of `ℚ_[p]` (it is `ℤ_[p]`), and the reduction
`red p : O p →+* ZMod p` (Mathlib's `PadicInt.toZMod`) has kernel the maximal ideal (`ker_red`).
* `p = 7`: this is what `QuadraticPlanes.sumSqGraph_colorable` needs (as `-1` is not a square modulo 7, unit
  vectors are integral, and a 4-colouring of the unit-distance graph of `𝔽₇²` colours the plane).
* `p = 3, 2`: `sumSqGraph_colorable_of` colours a point `z` by a colouring of `(ZMod p)²` at the residues of
  `z - rep z`, where `rep z` represents the coset `z + O²`, as soon as unit vectors are integral. They are
  (`mem_of_sum_sq_eq_one_of`), since `1 + t²` never vanishes in `ZMod 3`, nor in `ZMod 4` (for `p = 2`,
  through Mathlib's `PadicInt.toZModPow 2`). A unit vector `(a, b)` of `(ZMod p)²` has `a + b ≠ 0`, so
  `(x, y) ↦ x + y` is a proper colouring with `p` colours.

*Lower bounds.*
* `p = 7`: `11 ≡ 2² (mod 7)`, so Hensel's lemma gives `s ∈ ℤ_[7]` with `s² = 11` (`exists_sq_eq_eleven`).
  In any field of characteristic 0 with an element `s` such that `s² = 11`, the points `[a, b, c, e]` of
  `Sqrt11.P` go to `((a + bs)/30, (c + es)/30)`, and every listed edge of `Sqrt11.E` becomes a unit pair, by
  the integer identities that the kernel checked for `Sqrt11.checkEdges_E` (`adj_ptF`, which only uses
  `s² = 11`). This is a graph homomorphism, and the graph of `Sqrt11` is not 3-colourable
  (`Sqrt11.graph_not_colorable`, from an LRAT certificate): `not_colorable_three_of_sq_eq_eleven`.
* `p = 3`: `7 ≡ 1 (mod 3)` gives `t ∈ ℤ_[3]` with `t² = 7` (`exists_sq_eq_seven`), and
  `0, (1, 0), (1, 1), ((3 + t)/4, (3 - t)/4), (3/4, -t/4)` is a closed walk of length 5 with unit steps, so
  the graph is not 2-colourable.
* `p = 2`: `(0, 0) ~ (1, 0)`.
-/

namespace PadicPlanes

open LocalColouring QuadraticPlanes IsLocalRing Polynomial

instance fact_prime_seven : Fact (Nat.Prime 7) := ⟨by norm_num⟩

/-! ## The valuation subring `ℤ_[p]` of `ℚ_[p]` and its reduction -/

section Integers

variable (p : ℕ) [Fact p.Prime]

/-- The valuation subring `{x : ‖x‖ ≤ 1}` of `ℚ_[p]`, that is `ℤ_[p]`. -/
def O : ValuationSubring ℚ_[p] where
  toSubring := PadicInt.subring p
  mem_or_inv_mem' x := by
    by_cases h : ‖x‖ ≤ 1
    · exact Or.inl h
    · right
      show ‖x⁻¹‖ ≤ 1
      rw [norm_inv]
      exact inv_le_one_of_one_le₀ (le_of_lt (not_le.mp h))

lemma mem_O {x : ℚ_[p]} : x ∈ O p ↔ ‖x‖ ≤ 1 := Iff.rfl

/-- `O p` is `ℤ_[p]`. -/
def toPadicInt : O p ≃+* ℤ_[p] where
  toFun x := ⟨x.1, x.2⟩
  invFun x := ⟨x.1, x.2⟩
  left_inv _ := rfl
  right_inv _ := rfl
  map_mul' _ _ := rfl
  map_add' _ _ := rfl

lemma mem_maximalIdeal_iff {x : O p} :
    x ∈ maximalIdeal (O p) ↔ toPadicInt p x ∈ maximalIdeal ℤ_[p] := by
  rw [mem_maximalIdeal, mem_maximalIdeal, mem_nonunits_iff, mem_nonunits_iff, MulEquiv.isUnit_map]

/-- The reduction `O p →+* ZMod p`: Mathlib's `PadicInt.toZMod`. -/
noncomputable def red : O p →+* ZMod p := PadicInt.toZMod.comp (toPadicInt p).toRingHom

/-- The kernel of the reduction is the maximal ideal. -/
lemma ker_red : RingHom.ker (red p) = maximalIdeal (O p) := by
  ext x
  rw [RingHom.mem_ker, red, RingHom.comp_apply, ← RingHom.mem_ker, PadicInt.ker_toZMod,
    mem_maximalIdeal_iff]
  rfl

lemma red_sq_eq_zero {x : O p} (hx : x ∈ maximalIdeal (O p)) : red p (x ^ 2) = 0 := by
  rw [← ker_red] at hx
  rw [map_pow, RingHom.mem_ker.mp hx, zero_pow two_ne_zero]

end Integers

/-! ## Square roots, by Hensel's lemma -/

/-- `11 ≡ 2² (mod 7)`: a square root of `11` in `ℤ_[7]`. -/
lemma exists_sq_eq_eleven : ∃ s : ℤ_[7], s ^ 2 = 11 := by
  have h7 : ‖((-7 : ℤ) : ℤ_[7])‖ < 1 := by
    rw [PadicInt.norm_int_lt_one_iff_dvd]; norm_num
  have h4 : ‖((4 : ℤ) : ℤ_[7])‖ = 1 := by
    rw [PadicInt.norm_intCast_eq_one_iff]; norm_num
  have e1 : aeval (2 : ℤ_[7]) (X ^ 2 - C 11 : ℤ_[7][X]) = ((-7 : ℤ) : ℤ_[7]) := by
    simp; norm_num
  have e2 : aeval (2 : ℤ_[7]) (derivative (X ^ 2 - C 11 : ℤ_[7][X])) = ((4 : ℤ) : ℤ_[7]) := by
    simp; norm_num
  obtain ⟨z, hz, -⟩ := hensels_lemma (p := 7) (R := ℤ_[7]) (F := (X ^ 2 - C 11 : ℤ_[7][X])) (a := 2)
    (by rw [e1, e2, h4, one_pow]; exact h7)
  exact ⟨z, by simpa [sub_eq_zero] using hz⟩

/-- `7 ≡ 1² (mod 3)`: a square root of `7` in `ℤ_[3]`. -/
lemma exists_sq_eq_seven : ∃ t : ℤ_[3], t ^ 2 = 7 := by
  have h6 : ‖((-6 : ℤ) : ℤ_[3])‖ < 1 := by
    rw [PadicInt.norm_int_lt_one_iff_dvd]; norm_num
  have h2 : ‖((2 : ℤ) : ℤ_[3])‖ = 1 := by
    rw [PadicInt.norm_intCast_eq_one_iff]; norm_num
  have e1 : aeval (1 : ℤ_[3]) (X ^ 2 - C 7 : ℤ_[3][X]) = ((-6 : ℤ) : ℤ_[3]) := by
    simp; norm_num
  have e2 : aeval (1 : ℤ_[3]) (derivative (X ^ 2 - C 7 : ℤ_[3][X])) = ((2 : ℤ) : ℤ_[3]) := by
    simp; norm_num
  obtain ⟨z, hz, -⟩ := hensels_lemma (p := 3) (R := ℤ_[3]) (F := (X ^ 2 - C 7 : ℤ_[3][X])) (a := 1)
    (by rw [e1, e2, h2, one_pow]; exact h6)
  exact ⟨z, by simpa [sub_eq_zero] using hz⟩

/-! ## `p = 7` -/

/-- **Upper bound at 7.** `-1` is not a square modulo 7, so `QuadraticPlanes.sumSqGraph_colorable` applies to
`O 7` and the reduction to `ZMod 7`. -/
theorem colorable_four_seven : (sumSqGraph ℚ_[7]).Colorable 4 :=
  sumSqGraph_colorable (O 7) (red 7) (ker_red 7)

/-- The point `[a, b, c, e]` with denominator `D`, that is `((a + bs)/D, (c + es)/D)`, in a field with an
element `s`. -/
def ptF {F : Type*} [Field F] (s : F) (D : ℕ) (v : ℤ × ℤ × ℤ × ℤ) : F × F :=
  ((v.1 + v.2.1 * s) / D, (v.2.2.1 + v.2.2.2 * s) / D)

/-- The two integer identities of `QuadraticPlanes.unitPairB` give a unit distance in every field of
characteristic 0 with an element `s` such that `s² = d` (as `QuadraticPlanes.adj_pt` in `ℚ(√d)`). -/
lemma adj_ptF {F : Type*} [Field F] [CharZero F] {d D : ℕ} (hD : D ≠ 0) {s : F} (hs : s ^ 2 = d)
    {v w : ℤ × ℤ × ℤ × ℤ} (h : unitPairB d D v w = true) : (sumSqGraph F).Adj (ptF s D v) (ptF s D w) := by
  obtain ⟨h1, h2⟩ := unitPairB_spec h
  obtain ⟨a, b, c, e⟩ := v
  obtain ⟨a', b', c', e'⟩ := w
  dsimp only at h1 h2
  have h1' : ((a : F) - a') * ((a : F) - a') + d * (((b : F) - b') * ((b : F) - b'))
      + ((c : F) - c') * ((c : F) - c') + d * (((e : F) - e') * ((e : F) - e'))
      = (D : F) * D := by
    exact_mod_cast h1
  have h2' : ((a : F) - a') * ((b : F) - b') + ((c : F) - c') * ((e : F) - e') = 0 := by
    exact_mod_cast h2
  have hD' : (D : F) ≠ 0 := by exact_mod_cast hD
  show ((a + b * s) / D - (a' + b' * s) / D) ^ 2 + ((c + e * s) / D - (c' + e' * s) / D) ^ 2 = 1
  rw [div_sub_div_same, div_sub_div_same, div_pow, div_pow, ← add_div,
    div_eq_one_iff_eq (pow_ne_zero 2 hD')]
  linear_combination h1' + 2 * s * h2' + (((b : F) - b') ^ 2 + ((e : F) - e') ^ 2) * hs

/-- **Lower bound from `√11`.** If a field of characteristic 0 contains a square root of 11, the graph of
`x² + y²` over it is not 3-colourable: the 76-vertex graph of `Sqrt11` maps to it. -/
theorem not_colorable_three_of_sq_eq_eleven {F : Type*} [Field F] [CharZero F] {s : F} (hs : s ^ 2 = 11) :
    ¬ (sumSqGraph F).Colorable 3 := fun h =>
  Sqrt11.graph_not_colorable (h.of_hom (edgeGraph.hom (fun i => ptF s 30 (Sqrt11.P i))
    (fun e he => adj_ptF (by norm_num) (by rw [hs]; norm_num)
      (checkEdges_mem Sqrt11.checkEdges_E e he))))

/-- **Lower bound at 7.** -/
theorem not_colorable_three_seven : ¬ (sumSqGraph ℚ_[7]).Colorable 3 := by
  obtain ⟨s, hs⟩ := exists_sq_eq_eleven
  refine not_colorable_three_of_sq_eq_eleven (s := (s : ℚ_[7])) ?_
  rw [← PadicInt.coe_pow, hs]
  rfl

/-- **Theorem.** The unit-distance graph of the 7-adic plane has chromatic number 4. -/
theorem padicSeven_chromaticNumber : (QuadraticPlanes.sumSqGraph ℚ_[7]).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of colorable_four_seven not_colorable_three_seven

/-! ## Colouring through the residues, for `p = 3` and `p = 2` -/

section Residues

variable {F : Type*} [Field F] (O : ValuationSubring F)

/-- If `1 + t²` never vanishes in `R` and `g : O →+* R` kills the squares of the elements of the maximal
ideal, then `a² + b² = 1` forces `a ∈ O`. -/
lemma mem_of_sum_sq_eq_one_of {R : Type*} [CommRing R] (g : O →+* R)
    (hg : ∀ x ∈ maximalIdeal O, g (x ^ 2) = 0) (hR : ∀ t : R, 1 + t ^ 2 ≠ 0) {a b : F}
    (h : a ^ 2 + b ^ 2 = 1) : a ∈ O := by
  have key : ∀ t x : O, x ∈ maximalIdeal O → 1 + t ^ 2 ≠ x ^ 2 := fun t x hx he => by
    have h' := congrArg g he
    rw [map_add, map_one, map_pow, hg x hx] at h'
    exact hR _ h'
  by_contra ha
  have ha0 : a ≠ 0 := fun h0 => ha (h0 ▸ O.zero_mem)
  have hx : a⁻¹ ∈ O.nonunits := O.inv_mem_nonunits_iff.mpr (Or.inr ha)
  obtain ⟨hxO, hxm⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp hx
  by_cases ht : b / a ∈ O
  · -- `1 + (b/a)² = (a⁻¹)²`, with `a⁻¹ ∈ 𝔪`
    apply key ⟨b / a, ht⟩ ⟨a⁻¹, hxO⟩ hxm
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h
  · -- otherwise `a/b ∈ O`, and `1 + (a/b)² = (a/b · a⁻¹)²` with `a/b · a⁻¹ ∈ 𝔪`
    have hs : a / b ∈ O := by
      have := (O.mem_or_inv_mem (b / a)).resolve_left ht
      rwa [inv_div] at this
    have hb0 : b ≠ 0 := by
      rintro rfl
      exact ht (by simp)
    apply key ⟨a / b, hs⟩ (⟨a / b, hs⟩ * ⟨a⁻¹, hxO⟩) (Ideal.mul_mem_left _ _ hxm)
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h

variable {n : ℕ} (red : O →+* ZMod n)

open Classical in
/-- `red` extended by `0` outside `O`. -/
noncomputable def redExtN (x : F) : ZMod n := if h : x ∈ O then red ⟨x, h⟩ else 0

lemma redExtN_add {x y : F} (hx : x ∈ O) (hy : y ∈ O) :
    redExtN O red (x + y) = redExtN O red x + redExtN O red y := by
  simp only [redExtN, hx, hy, add_mem hx hy, dite_true]
  rw [← map_add]
  rfl

lemma redExtN_of_mem {x : F} (hx : x ∈ O) : redExtN O red x = red ⟨x, hx⟩ := by
  simp only [redExtN, hx, dite_true]

variable {ι : Type*} (c : ZMod n × ZMod n → ι)

/-- The colouring of `F²`: `c` of the residues of `z - rep z`. -/
noncomputable def colourN (z : F × F) : ι :=
  c (redExtN O red (z.1 - rep O z.1), redExtN O red (z.2 - rep O z.2))

lemma colourN_adj (hO : ∀ a b : F, a ^ 2 + b ^ 2 = 1 → a ∈ O)
    (hc : ∀ u w : ZMod n × ZMod n, (u.1 - w.1) ^ 2 + (u.2 - w.2) ^ 2 = 1 → c u ≠ c w)
    {p q : F × F} (hpq : (sumSqGraph F).Adj p q) : colourN O red c p ≠ colourN O red c q := by
  have h : (p.1 - q.1) ^ 2 + (p.2 - q.2) ^ 2 = 1 := hpq
  have ha : p.1 - q.1 ∈ O := hO _ _ h
  have hb : p.2 - q.2 ∈ O := hO _ (p.1 - q.1) (by linear_combination h)
  have hsq : red ⟨p.1 - q.1, ha⟩ ^ 2 + red ⟨p.2 - q.2, hb⟩ ^ 2 = 1 := by
    have e : (⟨p.1 - q.1, ha⟩ : O) ^ 2 + ⟨p.2 - q.2, hb⟩ ^ 2 = 1 := Subtype.ext (by simpa using h)
    have e' := congrArg red e
    rwa [map_add, map_pow, map_pow, map_one] at e'
  have e1 : p.1 - rep O p.1 = (p.1 - q.1) + (q.1 - rep O q.1) := by
    rw [rep_eq_of_sub_mem O ha]; ring
  have e2 : p.2 - rep O p.2 = (p.2 - q.2) + (q.2 - rep O q.2) := by
    rw [rep_eq_of_sub_mem O hb]; ring
  unfold colourN
  rw [e1, redExtN_add O red ha (sub_rep_mem O _), e2, redExtN_add O red hb (sub_rep_mem O _),
    redExtN_of_mem O red ha, redExtN_of_mem O red hb]
  apply hc
  simpa only [add_sub_cancel_right] using hsq

include red in
/-- **Colouring through the residues.** If unit vectors are integral for `O` and `c` is a proper colouring of
the graph of `x² + y²` on `(ZMod n)²`, then `sumSqGraph F` is `|ι|`-colourable. -/
theorem sumSqGraph_colorable_of [Fintype ι] (hO : ∀ a b : F, a ^ 2 + b ^ 2 = 1 → a ∈ O)
    (hc : ∀ u w : ZMod n × ZMod n, (u.1 - w.1) ^ 2 + (u.2 - w.2) ^ 2 = 1 → c u ≠ c w) :
    (sumSqGraph F).Colorable (Fintype.card ι) :=
  (SimpleGraph.Coloring.mk (colourN O red c) (fun h => colourN_adj O red c hO hc h)).colorable

end Residues

/-! ## `p = 3` -/

/-- Unit vectors of `ℚ_[3]²` are integral: `1 + t²` never vanishes in `ZMod 3`. -/
lemma mem_O_three {a b : ℚ_[3]} (h : a ^ 2 + b ^ 2 = 1) : a ∈ O 3 :=
  mem_of_sum_sq_eq_one_of (O 3) (red 3) (fun _ hx => red_sq_eq_zero 3 hx) (by decide) h

/-- **Upper bound at 3.** Colour by `x + y` modulo 3. -/
theorem colorable_three : (sumSqGraph ℚ_[3]).Colorable 3 := by
  simpa using sumSqGraph_colorable_of (O 3) (red 3) (fun u : ZMod 3 × ZMod 3 => u.1 + u.2)
    (fun _ _ h => mem_O_three h) (by decide)

/-- **Lower bound at 3.** With `t² = 7`, `0, (1, 0), (1, 1), ((3 + t)/4, (3 - t)/4), (3/4, -t/4)` is a closed
walk of length 5 with unit steps. -/
theorem not_colorable_two_three : ¬ (sumSqGraph ℚ_[3]).Colorable 2 := by
  obtain ⟨t₀, ht₀⟩ := exists_sq_eq_seven
  obtain ⟨t, ht⟩ : ∃ t : ℚ_[3], t ^ 2 = 7 := ⟨t₀, by rw [← PadicInt.coe_pow, ht₀]; rfl⟩
  rintro ⟨C⟩
  have h01 := C.valid (v := (0, 0)) (w := (1, 0)) (by
    show ((0 : ℚ_[3]) - 1) ^ 2 + (0 - 0) ^ 2 = 1
    ring)
  have h12 := C.valid (v := (1, 0)) (w := (1, 1)) (by
    show ((1 : ℚ_[3]) - 1) ^ 2 + (0 - 1) ^ 2 = 1
    ring)
  have h23 := C.valid (v := (1, 1)) (w := ((3 + t) / 4, (3 - t) / 4)) (by
    show ((1 : ℚ_[3]) - (3 + t) / 4) ^ 2 + (1 - (3 - t) / 4) ^ 2 = 1
    linear_combination (1 / 8 : ℚ_[3]) * ht)
  have h34 := C.valid (v := ((3 + t) / 4, (3 - t) / 4)) (w := (3 / 4, -t / 4)) (by
    show ((3 + t) / 4 - 3 / 4 : ℚ_[3]) ^ 2 + ((3 - t) / 4 - -t / 4) ^ 2 = 1
    linear_combination (1 / 16 : ℚ_[3]) * ht)
  have h40 := C.valid (v := (3 / 4, -t / 4)) (w := (0, 0)) (by
    show (3 / 4 - 0 : ℚ_[3]) ^ 2 + (-t / 4 - 0) ^ 2 = 1
    linear_combination (1 / 16 : ℚ_[3]) * ht)
  have odd : ∀ a b c d e : Fin 2, a ≠ b → b ≠ c → c ≠ d → d ≠ e → e ≠ a → False := by
    intro a b c d e
    omega
  exact odd _ _ _ _ _ h01 h12 h23 h34 h40

/-- **Theorem.** The unit-distance graph of the 3-adic plane has chromatic number 3. -/
theorem padicThree_chromaticNumber : (QuadraticPlanes.sumSqGraph ℚ_[3]).chromaticNumber = 3 := by
  rw [show (3 : ℕ∞) = ((2 : ℕ) : ℕ∞) + 1 by norm_num]
  exact SimpleGraph.chromaticNumber_eq_iff_colorable_not_colorable.mpr
    ⟨colorable_three, not_colorable_two_three⟩

/-! ## `p = 2` -/

/-- The reduction `O 2 →+* ZMod 4`: Mathlib's `PadicInt.toZModPow 2`. -/
noncomputable def red4 : O 2 →+* ZMod (2 ^ 2) := (PadicInt.toZModPow 2).comp (toPadicInt 2).toRingHom

/-- `red4` kills the squares of the elements of the maximal ideal `2 O`. -/
lemma red4_sq_eq_zero {x : O 2} (hx : x ∈ maximalIdeal (O 2)) : red4 (x ^ 2) = 0 := by
  rw [mem_maximalIdeal_iff, PadicInt.maximalIdeal_eq_span_p, Ideal.mem_span_singleton] at hx
  obtain ⟨y, hy⟩ := hx
  show PadicInt.toZModPow 2 (toPadicInt 2 (x ^ 2)) = 0
  rw [← RingHom.mem_ker, PadicInt.ker_toZModPow, Ideal.mem_span_singleton, map_pow, hy]
  exact ⟨y ^ 2, by ring⟩

/-- Unit vectors of `ℚ_[2]²` are integral: `1 + t²` never vanishes in `ZMod 4`. -/
lemma mem_O_two {a b : ℚ_[2]} (h : a ^ 2 + b ^ 2 = 1) : a ∈ O 2 :=
  mem_of_sum_sq_eq_one_of (O 2) red4 (fun _ hx => red4_sq_eq_zero hx) (by decide) h

/-- **Upper bound at 2.** Colour by `x + y` modulo 2. -/
theorem colorable_two : (sumSqGraph ℚ_[2]).Colorable 2 := by
  simpa using sumSqGraph_colorable_of (O 2) (red 2) (fun u : ZMod 2 × ZMod 2 => u.1 + u.2)
    (fun _ _ h => mem_O_two h) (by decide)

/-- **Lower bound at 2.** `(0, 0) ~ (1, 0)`. -/
theorem not_colorable_one_two : ¬ (sumSqGraph ℚ_[2]).Colorable 1 := by
  rintro ⟨C⟩
  exact C.valid (v := (0, 0)) (w := (1, 0)) (by
    show ((0 : ℚ_[2]) - 1) ^ 2 + (0 - 0) ^ 2 = 1
    ring) (Subsingleton.elim _ _)

/-- **Theorem.** The unit-distance graph of the 2-adic plane has chromatic number 2. -/
theorem padicTwo_chromaticNumber : (QuadraticPlanes.sumSqGraph ℚ_[2]).chromaticNumber = 2 := by
  rw [show (2 : ℕ∞) = ((1 : ℕ) : ℕ∞) + 1 by norm_num]
  exact SimpleGraph.chromaticNumber_eq_iff_colorable_not_colorable.mpr
    ⟨colorable_two, not_colorable_one_two⟩

end PadicPlanes
