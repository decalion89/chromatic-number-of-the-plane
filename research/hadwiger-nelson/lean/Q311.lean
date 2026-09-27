import LocalColouring

/-!
# The plane over `ℚ(√3, √11)` has chromatic number 4

**Theorem** (`Q311.chromaticNumber_eq_four`; K. G. Fischer, 1994). Let `L = ℚ(√3, √11) ⊆ ℝ`.
The unit-distance graph on `L²`, in which `p ~ q` iff `(p₁ - q₁)² + (p₂ - q₂)² = 1`, has
chromatic number 4.

*Upper bound.* We apply the criterion `LocalColouring.colorable_four_of_residue_two`: `√3 ∈ L`,
and every valuation subring `O` of `L` with `2` in its maximal ideal `𝔪` has residue field
`𝔽₂`. Chevalley's extension theorem provides such an `O`.

The field `L` has two places over 2. Each has completion `ℚ₂(√3)`, ramified over `ℚ₂` with
residue field `𝔽₂`. The proof follows this picture.
1. `π = √3 - 1` is a root of `X² + 2X - 2`. So `π ∈ 𝔪`, and `2 = π²(π + 3)`.
2. `β = (1 + √33)/2` satisfies `β(β - 1) = 8`, so `β ∈ 𝔪` or `β - 1 ∈ 𝔪`. Changing the sign of
   `√11` swaps `β` and `1 - β`, so we may assume `β ∈ 𝔪`. Then `β ∈ 8O`.
3. `β` is a 2-adic integer at `O` (Hensel's lemma): for every `k` some integer `b` has
   `β - b ∈ 2ᵏ⁺¹O`. Induct on `k`, using `β² = β + 8`.
4. `L` is spanned over `ℚ` by `1, π, β, πβ`. For `x ∈ O` write
   `x = (n₀ + n₁π + (n₂ + n₃π)β)/D` with `D = 2ʲD'` and `D'` odd. Replacing `β` by an integer
   `b` with `β - b ∈ 2ʲ⁺¹O` changes `x` by an element of `2O`. The result
   `y = (n₀ + n₂b + (n₁ + n₃b)π)/D` lies in `ℚ(√3)`.
5. In `ℚ(√3)` the residue is computed as for `ℚ(√2, √3)`. If `N₀ + N₁π ∈ 2O`, both `Nᵢ` are
   even; descent on the power of 2 in `D` gives `y ≡ N₀ (mod 𝔪)`. So `x ≡ 0` or `x ≡ 1`.

*Lower bound.* The Moser spindle (`certificates/moser_spindle_no3coloring.json`: 7 vertices,
11 edges) lies in `L²`. It is two unit rhombi with the common tip `0`. In a 3-colouring the
far tips `3` and `6` both take the colour of `0`, yet they are at distance 1.
-/

namespace Q311

open LocalColouring

/-- `L = ℚ(√3, √11)`, as a subfield of `ℝ`. -/
noncomputable def L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√3, √11}

lemma sqrt3_mem : √3 ∈ L := IntermediateField.subset_adjoin ℚ _ (by simp)
lemma sqrt11_mem : √11 ∈ L := IntermediateField.subset_adjoin ℚ _ (by simp)

/-- `√3` as an element of `L`. -/
noncomputable def r3 : L := ⟨√3, sqrt3_mem⟩
/-- `√11` as an element of `L`. -/
noncomputable def r11 : L := ⟨√11, sqrt11_mem⟩

@[simp] lemma coe_r3 : (r3 : ℝ) = √3 := rfl
@[simp] lemma coe_r11 : (r11 : ℝ) = √11 := rfl

lemma sqrt3_sq : √3 ^ 2 = (3 : ℝ) := Real.sq_sqrt (by norm_num)
lemma sqrt11_sq : √11 ^ 2 = (11 : ℝ) := Real.sq_sqrt (by norm_num)

lemma r3_sq : r3 ^ 2 = 3 := Subtype.ext (by push_cast; exact sqrt3_sq)
lemma r11_sq : r11 ^ 2 = 11 := Subtype.ext (by push_cast; exact sqrt11_sq)

/-- Every element of `L` is `(n₀ + n₁√3 + n₂t + n₃√3t)/d`, for either sign `t = ±√11`. -/
lemma exists_int_combo_L {t : L} (ht : t = r11 ∨ t = -r11) (x : L) :
    ∃ d : ℕ, 0 < d ∧ ∃ n₀ n₁ n₂ n₃ : ℤ,
      (d : L) * x = n₀ + n₁ * r3 + n₂ * t + n₃ * (r3 * t) := by
  obtain ⟨d, hd, n₀, n₁, n₂, n₃, h⟩ := exists_int_combo (m := 3) (n := 11)
    (by push_cast; exact sqrt3_sq) (by push_cast; exact sqrt11_sq) x.2
  have h' : (d : L) * x = n₀ + n₁ * r3 + n₂ * r11 + n₃ * (r3 * r11) :=
    Subtype.ext (by push_cast; exact h)
  rcases ht with rfl | rfl
  · exact ⟨d, hd, n₀, n₁, n₂, n₃, h'⟩
  · exact ⟨d, hd, n₀, n₁, -n₂, -n₃, by rw [h']; push_cast; ring⟩

/-! ## The elements `π` and `β` -/

/-- `π = √3 - 1`, a root of the 2-Eisenstein polynomial `X² + 2X - 2`. -/
noncomputable def piL : L := r3 - 1

lemma piL_sq : piL ^ 2 = -2 * piL + 2 := by
  unfold piL
  linear_combination r3_sq

lemma two_eq : (2 : L) = piL ^ 2 * (piL + 3) := by
  linear_combination -(piL + 1) * piL_sq

lemma piL_ne_zero : piL ≠ 0 := by
  intro h
  have := piL_sq
  rw [h] at this
  norm_num at this

/-- For `t = ±√11`, the element `β = (1 + √3 t)/2`; it satisfies `β² = β + 8`. -/
noncomputable def beta (t : L) : L := (1 + r3 * t) / 2

lemma beta_sq {t : L} (ht : t ^ 2 = 11) : beta t ^ 2 = beta t + 8 := by
  unfold beta
  linear_combination (t ^ 2 / 4) * r3_sq + (3 / 4) * ht

/-- In the basis `1, π, β, πβ`: every element of `L` is `(n₀ + n₁π + (n₂ + n₃π)β)/D`. -/
lemma exists_int_combo_beta {t : L} (ht : t = r11 ∨ t = -r11) (x : L) :
    ∃ D : ℕ, 0 < D ∧ ∃ n₀ n₁ n₂ n₃ : ℤ,
      (D : L) * x = n₀ + n₁ * piL + (n₂ + n₃ * piL) * beta t := by
  obtain ⟨d, hd, a, b, c, e, h⟩ := exists_int_combo_L ht x
  refine ⟨3 * d, by omega, 3 * a + 3 * b - c - 3 * e, 3 * b - c, 2 * c + 6 * e, 2 * c, ?_⟩
  unfold piL beta
  push_cast
  linear_combination 3 * h + (-(c : L) * t) * r3_sq

/-! ## The residue field is `𝔽₂` -/

section Residue

variable (O : ValuationSubring L)

lemma piL_mem : piL ∈ O :=
  mem_of_sq_eq O (neg_mem (ofNat_mem O 2)) (ofNat_mem O 2) piL_sq

lemma beta_mem {t : L} (ht : t ^ 2 = 11) : beta t ∈ O :=
  mem_of_sq_eq O O.one_mem (ofNat_mem O 8) (by rw [one_mul]; exact beta_sq ht)

variable (h2 : (2 : L) ∈ O.nonunits)
include h2

lemma piL_mem_nonunits : piL ∈ O.nonunits := by
  refine mem_nonunits_of_pow_mem O (piL_mem O) (n := 2) ?_
  rw [show piL ^ 2 = (1 - piL) * 2 by linear_combination piL_sq]
  exact mul_mem_nonunits_left O (sub_mem O.one_mem (piL_mem O)) h2

/-- Choosing the sign of `√11`, we may assume `β ∈ 𝔪`. -/
lemma exists_beta_mem_nonunits : ∃ t : L, (t = r11 ∨ t = -r11) ∧ beta t ∈ O.nonunits := by
  have hβ := beta_mem O r11_sq
  have h8 : beta r11 * (beta r11 - 1) ∈ O.nonunits := by
    rw [show beta r11 * (beta r11 - 1) = 4 * 2 by linear_combination beta_sq r11_sq]
    exact mul_mem_nonunits_left O (ofNat_mem O 4) h2
  rcases mem_nonunits_or_of_mul_mem O hβ (sub_mem hβ O.one_mem) h8 with h | h
  · exact ⟨r11, Or.inl rfl, h⟩
  · refine ⟨-r11, Or.inr rfl, ?_⟩
    rw [show beta (-r11) = -(beta r11 - 1) by unfold beta; ring]
    exact neg_mem h

/-- **Hensel's lemma for `β`.** If `β ∈ 𝔪`, then for every `k` there is an integer `b` with
`β - b ∈ 2ᵏ⁺¹O`. -/
lemma hensel {t : L} (ht : t ^ 2 = 11) (hβ : beta t ∈ O.nonunits) (k : ℕ) :
    ∃ b : ℤ, ∃ γ ∈ O, beta t - b = 2 ^ (k + 1) * γ := by
  induction k with
  | zero =>
    -- `β(β - 1) = 8` with `β - 1` a unit, so `β = 2 · 4(β - 1)⁻¹`
    have h1 : beta t - 1 ∉ O.nonunits := fun h => one_notMem_nonunits O (by
      have := sub_mem hβ h
      rwa [sub_sub_cancel] at this)
    have hne : beta t - 1 ≠ 0 := fun h => h1 (h ▸ zero_mem _)
    refine ⟨0, 4 * (beta t - 1)⁻¹, O.mul_mem _ _ (ofNat_mem O 4)
      (inv_mem_of_notMem_nonunits O h1), ?_⟩
    simp only [Int.cast_zero, sub_zero, zero_add, pow_one]
    field_simp
    linear_combination beta_sq ht
  | succ k ih =>
    obtain ⟨b, γ, hγ, hb⟩ := ih
    set G : O := ⟨γ, hγ⟩
    -- `N = b² - b - 8` is `2ᵏ⁺¹` times an element of `O`, hence divisible by `2ᵏ⁺¹`
    have hN : ((b ^ 2 - b - 8 : ℤ) : L) = 2 ^ (k + 1) * (-γ * (2 * b - 1 + 2 ^ (k + 1) * γ)) := by
      push_cast
      linear_combination beta_sq ht - (b + beta t + 2 ^ (k + 1) * γ - 1) * hb
    have hmem : -γ * (2 * b - 1 + 2 ^ (k + 1) * γ) ∈ O := by
      have := (-G * (2 * b - 1 + 2 ^ (k + 1) * G)).2
      push_cast at this
      exact this
    obtain ⟨c, hc⟩ := two_pow_dvd_of_eq O h2 hmem hN
    have hc' : (b : L) ^ 2 - b - 8 = 2 ^ (k + 1) * c := by exact_mod_cast hc
    refine ⟨b + 2 ^ (k + 1) * c, b * γ + 2 ^ k * γ ^ 2, ?_, ?_⟩
    · have := ((b : O) * G + 2 ^ k * G ^ 2).2
      push_cast at this
      exact this
    · push_cast
      linear_combination (b + beta t + 2 ^ (k + 1) * γ) * hb + hc' - beta_sq ht

/-- If `N₀ + N₁π ∈ 2O` with integers `Nᵢ`, then both are even. -/
lemma even_of_add_mul_piL {N₀ N₁ : ℤ} {b : L} (hb : b ∈ O) (h : (N₀ : L) + N₁ * piL = 2 * b) :
    Even N₀ ∧ Even N₁ := by
  have hπ := piL_mem_nonunits O h2
  have hπO := piL_mem O
  have hε : piL + 3 ∈ O := add_mem hπO (ofNat_mem O 3)
  have e0 : Even N₀ := even_of_intCast_mem_nonunits O h2 (by
    rw [show (N₀ : L) = piL * (piL * (piL + 3) * b - N₁) by linear_combination h + b * two_eq]
    exact mul_mem_nonunits_right O hπ
      (sub_mem (O.mul_mem _ _ (O.mul_mem _ _ hπO hε) hb) (intCast_mem O N₁)))
  obtain ⟨M₀, hM₀⟩ := e0
  have hM₀' : (N₀ : L) = M₀ + M₀ := by rw [hM₀]; push_cast; ring
  refine ⟨⟨M₀, hM₀⟩, even_of_intCast_mem_nonunits O h2 ?_⟩
  rw [show (N₁ : L) = piL * ((piL + 3) * (b - M₀)) by
    apply mul_left_cancel₀ piL_ne_zero
    linear_combination h - hM₀' + (b - M₀) * two_eq]
  exact mul_mem_nonunits_right O hπ (O.mul_mem _ _ hε (sub_mem hb (intCast_mem O M₀)))

/-- The residue computation in `ℚ(√3)`: if `y ∈ O` and `D y = N₀ + N₁π` with `D > 0`, then
`y ≡ 0` or `y ≡ 1`. -/
lemma residue_of_add_mul_piL {y : L} (hy : y ∈ O) : ∀ D : ℕ, 0 < D → ∀ N₀ N₁ : ℤ,
    (D : L) * y = N₀ + N₁ * piL → y ∈ O.nonunits ∨ y - 1 ∈ O.nonunits := by
  have hπ := piL_mem_nonunits O h2
  intro D
  induction D using Nat.strong_induction_on with
  | _ D ih =>
  intro hD N₀ N₁ h
  rcases Nat.even_or_odd D with ⟨D', rfl⟩ | ⟨D', rfl⟩
  · obtain ⟨⟨M₀, rfl⟩, ⟨M₁, rfl⟩⟩ := even_of_add_mul_piL O h2 (N₀ := N₀) (N₁ := N₁)
      (b := D' * y) (O.mul_mem _ _ (natCast_mem O D') hy) (by push_cast at h; linear_combination -h)
    refine ih D' (by omega) (by omega) M₀ M₁ ?_
    apply mul_left_cancel₀ (two_ne_zero (α := L))
    push_cast at h
    linear_combination h
  · have hyn : y - N₀ ∈ O.nonunits := by
      rw [show y - N₀ = piL * N₁ - ((D' : L) * y) * 2 by push_cast at h; linear_combination h]
      exact sub_mem (mul_mem_nonunits_right O hπ (intCast_mem O N₁))
        (mul_mem_nonunits_left O (O.mul_mem _ _ (natCast_mem O D') hy) h2)
    rcases Int.even_or_odd' N₀ with ⟨k, hk | hk⟩
    · left
      rw [show y = (y - N₀) + (k : L) * 2 by rw [hk]; push_cast; ring]
      exact add_mem hyn (mul_mem_nonunits_left O (intCast_mem O k) h2)
    · right
      rw [show y - 1 = (y - N₀) + (k : L) * 2 by rw [hk]; push_cast; ring]
      exact add_mem hyn (mul_mem_nonunits_left O (intCast_mem O k) h2)

/-- **The residue field is `𝔽₂`:** every `x ∈ O` is `≡ 0` or `≡ 1` modulo `𝔪`. -/
theorem mem_nonunits_or_sub_one_mem {x : L} (hx : x ∈ O) :
    x ∈ O.nonunits ∨ x - 1 ∈ O.nonunits := by
  obtain ⟨t, ht, hβ⟩ := exists_beta_mem_nonunits O h2
  have ht2 : t ^ 2 = 11 := by
    rcases ht with rfl | rfl
    · exact r11_sq
    · rw [neg_sq]; exact r11_sq
  obtain ⟨D, hD, n₀, n₁, n₂, n₃, h⟩ := exists_int_combo_beta ht x
  obtain ⟨j, D', hD', rfl⟩ := Nat.exists_eq_two_pow_mul_odd hD.ne'
  obtain ⟨b, γ, hγ, hb⟩ := hensel O h2 ht2 hβ j
  have hD'0 : (D' : L) ≠ 0 := by exact_mod_cast hD'.pos.ne'
  -- the odd part `D'` is a unit of `O`
  have hinv : (D' : L)⁻¹ ∈ O := inv_mem_of_notMem_nonunits O fun hm => by
    have := even_of_intCast_mem_nonunits O h2 (n := D') (by exact_mod_cast hm)
    exact Nat.not_even_iff_odd.mpr hD' ((Int.even_coe_nat D').mp this)
  -- `z = 2 (n₂ + n₃π) γ / D'` is in `2O`, and `x - z` lies in `ℚ(√3)`
  have hw : (n₂ + n₃ * piL) * γ * (D' : L)⁻¹ ∈ O :=
    O.mul_mem _ _ (O.mul_mem _ _ (add_mem (intCast_mem O n₂)
      (O.mul_mem _ _ (intCast_mem O n₃) (piL_mem O))) hγ) hinv
  have hz : 2 * ((n₂ + n₃ * piL) * γ * (D' : L)⁻¹) ∈ O.nonunits :=
    mul_mem_nonunits_right O h2 hw
  have hy : x - 2 * ((n₂ + n₃ * piL) * γ * (D' : L)⁻¹) ∈ O :=
    sub_mem hx (O.nonunits_subset hz)
  have hyD : ((2 ^ j * D' : ℕ) : L) * (x - 2 * ((n₂ + n₃ * piL) * γ * (D' : L)⁻¹))
      = ((n₀ + n₂ * b : ℤ) : L) + ((n₁ + n₃ * b : ℤ) : L) * piL := by
    have hDD : (D' : L) * (D' : L)⁻¹ = 1 := mul_inv_cancel₀ hD'0
    push_cast at h ⊢
    linear_combination h + (n₂ + n₃ * piL) * hb - 2 ^ (j + 1) * (n₂ + n₃ * piL) * γ * hDD
  rcases residue_of_add_mul_piL O h2 hy _ hD _ _ hyD with h' | h'
  · left
    have := add_mem h' hz
    rwa [sub_add_cancel] at this
  · right
    have := add_mem h' hz
    rwa [sub_right_comm, sub_add_cancel] at this

end Residue

/-- **Upper bound.** The unit-distance graph of `ℚ(√3, √11)²` is 4-colourable. -/
theorem colorable_four : (unitDistGraph L).Colorable 4 := by
  obtain ⟨O, h2⟩ := exists_valuationSubring_two_mem_nonunits L
  exact colorable_four_of_residue_two r3 r3_sq O h2
    fun x hx => mem_nonunits_or_sub_one_mem O h2 hx

/-! ## Lower bound: the Moser spindle -/

/-- `a + b√3 + c√11 + d√33` as an element of `L`. -/
noncomputable def mkL (a b c d : ℚ) : L := a + b * r3 + c * r11 + d * (r3 * r11)

@[simp] lemma coe_mkL (a b c d : ℚ) :
    (mkL a b c d : ℝ) = a + b * √3 + c * √11 + d * (√3 * √11) := by
  simp [mkL]

/-! The seven vertices of `certificates/moser_spindle_no3coloring.json`, on the basis
`1, √3, √11, √33`. -/

noncomputable def v0 : L × L := (mkL 0 0 0 0, mkL 0 0 0 0)
noncomputable def v1 : L × L := (mkL 1 0 0 0, mkL 0 0 0 0)
noncomputable def v2 : L × L := (mkL (1/2) 0 0 0, mkL 0 (1/2) 0 0)
noncomputable def v3 : L × L := (mkL (3/2) 0 0 0, mkL 0 (1/2) 0 0)
noncomputable def v4 : L × L := (mkL (5/6) 0 0 0, mkL 0 0 (1/6) 0)
noncomputable def v5 : L × L := (mkL (5/12) 0 0 (-1/12), mkL 0 (5/12) (1/12) 0)
noncomputable def v6 : L × L := (mkL (5/4) 0 0 (-1/12), mkL 0 (5/12) (1/4) 0)

lemma unitDist_iff (p q : L × L) :
    (unitDistGraph L).Adj p q ↔ ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1 := Iff.rfl

/-! The eleven unit distances. Each is a polynomial identity modulo `√3² = 3`, `√11² = 11`. -/

lemma adj01 : (unitDistGraph L).Adj v0 v1 := by
  rw [unitDist_iff]; simp only [v0, v1, coe_mkL]; push_cast
  ring
lemma adj02 : (unitDistGraph L).Adj v0 v2 := by
  rw [unitDist_iff]; simp only [v0, v2, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj04 : (unitDistGraph L).Adj v0 v4 := by
  rw [unitDist_iff]; simp only [v0, v4, coe_mkL]; push_cast
  linear_combination (1 / 36) * sqrt11_sq
lemma adj05 : (unitDistGraph L).Adj v0 v5 := by
  rw [unitDist_iff]; simp only [v0, v5, coe_mkL]; push_cast
  linear_combination (√11 ^ 2 / 144 + 25 / 144) * sqrt3_sq + (1 / 36) * sqrt11_sq
lemma adj12 : (unitDistGraph L).Adj v1 v2 := by
  rw [unitDist_iff]; simp only [v1, v2, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj13 : (unitDistGraph L).Adj v1 v3 := by
  rw [unitDist_iff]; simp only [v1, v3, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj23 : (unitDistGraph L).Adj v2 v3 := by
  rw [unitDist_iff]; simp only [v2, v3, coe_mkL]; push_cast
  ring
lemma adj36 : (unitDistGraph L).Adj v3 v6 := by
  rw [unitDist_iff]; simp only [v3, v6, coe_mkL]; push_cast
  linear_combination (√11 ^ 2 / 144 + 1 / 144) * sqrt3_sq + (1 / 12) * sqrt11_sq
lemma adj45 : (unitDistGraph L).Adj v4 v5 := by
  rw [unitDist_iff]; simp only [v4, v5, coe_mkL]; push_cast
  linear_combination (√11 ^ 2 / 144 + 25 / 144) * sqrt3_sq + (1 / 36) * sqrt11_sq
lemma adj46 : (unitDistGraph L).Adj v4 v6 := by
  rw [unitDist_iff]; simp only [v4, v6, coe_mkL]; push_cast
  linear_combination (√11 ^ 2 / 144 + 25 / 144) * sqrt3_sq + (1 / 36) * sqrt11_sq
lemma adj56 : (unitDistGraph L).Adj v5 v6 := by
  rw [unitDist_iff]; simp only [v5, v6, coe_mkL]; push_cast
  linear_combination (1 / 36) * sqrt11_sq

/-- The edges of the Moser spindle: the unit distances among the seven vertices. -/
def moserEdges : List (Fin 7 × Fin 7) :=
  [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (2, 3), (3, 6), (4, 5), (4, 6), (5, 6)]

/-- The Moser spindle: unit rhombi `0-{1,2}-3` and `0-{4,5}-6`, closed by the edge `3-6`. -/
abbrev moser : SimpleGraph (Fin 7) := edgeGraph moserEdges

/-- The Moser spindle is not 3-colourable. -/
theorem moser_not_colorable : ¬ moser.Colorable 3 := by
  rintro ⟨C⟩
  have v : ∀ i j, moser.Adj i j → C i ≠ C j := fun i j h => C.valid h
  have r1 : C 0 = C 3 := fin3_rhombus _ _ _ _ (v 0 1 (by decide)) (v 0 2 (by decide))
    (v 1 2 (by decide)) (v 3 1 (by decide)) (v 3 2 (by decide))
  have r2 : C 0 = C 6 := fin3_rhombus _ _ _ _ (v 0 4 (by decide)) (v 0 5 (by decide))
    (v 4 5 (by decide)) (v 6 4 (by decide)) (v 6 5 (by decide))
  exact v 3 6 (by decide) (r1.symm.trans r2)

/-- The positions of the seven vertices. -/
noncomputable def moserPt : Fin 7 → L × L := ![v0, v1, v2, v3, v4, v5, v6]

lemma moserPt_adj : ∀ e ∈ moserEdges, (unitDistGraph L).Adj (moserPt e.1) (moserPt e.2) := by
  simp only [moserEdges, List.forall_mem_cons, List.not_mem_nil, IsEmpty.forall_iff,
    implies_true, and_true]
  exact ⟨adj01, adj02, adj04, adj05, adj12, adj13, adj23, adj36, adj45, adj46, adj56⟩

/-- **Lower bound.** The unit-distance graph of `ℚ(√3, √11)²` is not 3-colourable. -/
theorem not_colorable_three : ¬ (unitDistGraph L).Colorable 3 := fun h =>
  moser_not_colorable (h.of_hom (edgeGraph.hom moserPt moserPt_adj))

/-- **Theorem** (Fischer). The unit-distance graph of `ℚ(√3, √11)²` has chromatic number 4. -/
theorem chromaticNumber_eq_four : (unitDistGraph L).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of colorable_four not_colorable_three

end Q311
