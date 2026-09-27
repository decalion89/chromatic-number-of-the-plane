import LocalColouring

/-!
# The plane over `ℚ(√2, √3)` has chromatic number 4

**Theorem** (`Q23.chromaticNumber_eq_four`). Let `L = ℚ(√2, √3) ⊆ ℝ`. The unit-distance graph
on `L²`, in which `p ~ q` iff `(p₁ - q₁)² + (p₂ - q₂)² = 1`, has chromatic number 4.

*Upper bound.* We apply the criterion `LocalColouring.colorable_four_of_residue_two`: `√3 ∈ L`,
and `L` has a valuation subring `O` with `2` in its maximal ideal and residue field `𝔽₂`.
1. Chevalley's extension theorem gives a valuation subring `O` with `2` in its maximal ideal.
2. Put `π = (√2 + √6)/2 - 1`. It is a root of the 2-Eisenstein polynomial
   `X⁴ + 4X³ + 2X² - 4X - 2`. So `π ∈ O`, `π` lies in the maximal ideal `𝔪`, and
   `2 = π⁴ε` with `ε = 43 - 2π - 35π² - 10π³`.
3. `L` is spanned over `ℚ` by `1, π, π², π³`. If `n₀ + n₁π + n₂π² + n₃π³ ∈ 2O` with integers
   `nᵢ`, then every `nᵢ` is even: peel off one power of `π` at a time.
4. Let `x ∈ O`, and write `x = (n₀ + n₁π + n₂π² + n₃π³)/d`. Descent on the power of 2 in `d`
   gives `x ≡ n₀ (mod 𝔪)`, so the residue field is `𝔽₂`.

The place is totally ramified over 2. The proof never uses `[L : ℚ] = 4` or a ring of integers.

*Lower bound.* The chain of three unit rhombi of `data/chain23.json` (10 vertices, 16 edges)
lies in `L²`. In a 3-colouring the two tips of each rhombus agree, so the first and last tips
agree, yet they are at distance 1.
-/

namespace Q23

open LocalColouring

/-- `L = ℚ(√2, √3)`, as a subfield of `ℝ`. -/
noncomputable def L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√2, √3}

lemma sqrt2_mem : √2 ∈ L := IntermediateField.subset_adjoin ℚ _ (by simp)
lemma sqrt3_mem : √3 ∈ L := IntermediateField.subset_adjoin ℚ _ (by simp)

/-- `√2` as an element of `L`. -/
noncomputable def s2 : L := ⟨√2, sqrt2_mem⟩
/-- `√3` as an element of `L`. -/
noncomputable def s3 : L := ⟨√3, sqrt3_mem⟩

@[simp] lemma coe_s2 : (s2 : ℝ) = √2 := rfl
@[simp] lemma coe_s3 : (s3 : ℝ) = √3 := rfl

lemma sqrt2_sq : √2 ^ 2 = (2 : ℝ) := Real.sq_sqrt (by norm_num)
lemma sqrt3_sq : √3 ^ 2 = (3 : ℝ) := Real.sq_sqrt (by norm_num)

lemma s2_sq : s2 ^ 2 = 2 := Subtype.ext (by push_cast; exact sqrt2_sq)
lemma s3_sq : s3 ^ 2 = 3 := Subtype.ext (by push_cast; exact sqrt3_sq)

/-- Every element of `L` is `(n₀ + n₁√2 + n₂√3 + n₃√6)/d`. -/
lemma exists_int_combo_L (x : L) : ∃ d : ℕ, 0 < d ∧ ∃ n₀ n₁ n₂ n₃ : ℤ,
    (d : L) * x = n₀ + n₁ * s2 + n₂ * s3 + n₃ * (s2 * s3) := by
  obtain ⟨d, hd, n₀, n₁, n₂, n₃, h⟩ := exists_int_combo (m := 2) (n := 3)
    (by push_cast; exact sqrt2_sq) (by push_cast; exact sqrt3_sq) x.2
  exact ⟨d, hd, n₀, n₁, n₂, n₃, Subtype.ext (by push_cast; exact h)⟩

/-! ## The uniformizer at the place over 2 -/

/-- The uniformizer `π = θ - 1`, where `θ = (√2 + √6)/2 = 2 cos(π/12)`. -/
noncomputable def piL : L := (s2 + s2 * s3) / 2 - 1

/-- `π` is a root of the 2-Eisenstein polynomial `X⁴ + 4X³ + 2X² - 4X - 2`. -/
lemma piL_eisenstein : piL ^ 4 + 4 * piL ^ 3 + 2 * piL ^ 2 - 4 * piL - 2 = 0 := by
  have h2 := s2_sq
  have h3 := s3_sq
  unfold piL
  linear_combination (s2^2*s3^4/16 + s2^2*s3^3/4 + 3*s2^2*s3^2/8 + s2^2*s3/4 + s2^2/16 + s3^4/8
    + s3^3/2 - s3^2/4 - 3*s3/2 - 7/8) * h2 + (s3^2/4 + s3 + 1/4) * h3

lemma s2_eq : s2 = piL ^ 3 + 3 * piL ^ 2 - 2 := by
  have h2 := s2_sq
  have h3 := s3_sq
  unfold piL
  linear_combination (-s2*s3^3/8 - 3*s2*s3^2/8 - 3*s2*s3/8 - s2/8) * h2
    + (-s2*s3/4 - 3*s2/4) * h3

lemma s3_eq : s3 = piL ^ 2 + 2 * piL - 1 := by
  have h2 := s2_sq
  have h3 := s3_sq
  unfold piL
  linear_combination (-s3^2/4 - s3/2 - 1/4) * h2 + (-1/2) * h3

lemma s6_eq : s2 * s3 = -piL ^ 3 - 3 * piL ^ 2 + 2 * piL + 4 := by
  have h2 := s2_sq
  have h3 := s3_sq
  unfold piL
  linear_combination (s2*s3^3/8 + 3*s2*s3^2/8 + 3*s2*s3/8 + s2/8) * h2 + (s2*s3/4 + 3*s2/4) * h3

/-- Every element of `L` is `(n₀ + n₁π + n₂π² + n₃π³)/d`. -/
lemma exists_int_combo_pi (x : L) : ∃ d : ℕ, 0 < d ∧ ∃ n₀ n₁ n₂ n₃ : ℤ,
    (d : L) * x = n₀ + n₁ * piL + n₂ * piL ^ 2 + n₃ * piL ^ 3 := by
  obtain ⟨d, hd, a, b, c, e, h⟩ := exists_int_combo_L x
  refine ⟨d, hd, a - 2 * b - c + 4 * e, 2 * c + 2 * e, 3 * b + c - 3 * e, b - e, ?_⟩
  rw [h]
  push_cast
  linear_combination (b : L) * s2_eq + (c : L) * s3_eq + (e : L) * s6_eq

/-- The unit `ε` with `2 = π⁴ε`. -/
noncomputable def epsL : L := 43 - 2 * piL - 35 * piL ^ 2 - 10 * piL ^ 3

lemma two_eq_pi_pow_four_mul : (2 : L) = piL ^ 4 * epsL := by
  unfold epsL
  linear_combination (10 * piL ^ 3 - 5 * piL ^ 2 + 2 * piL - 1) * piL_eisenstein

lemma piL_ne_zero : piL ≠ 0 := by
  intro h
  have := piL_eisenstein
  rw [h] at this
  norm_num at this

/-! ## The residue field is `𝔽₂` -/

section Residue

variable (O : ValuationSubring L)

/-- `π` is integral: it satisfies a monic equation. -/
lemma piL_mem : piL ∈ O := by
  by_contra hπ
  have hi : piL⁻¹ ∈ O.nonunits := O.inv_mem_nonunits_iff.mpr (Or.inr hπ)
  have hiO : piL⁻¹ ∈ O := O.nonunits_subset hi
  apply one_notMem_nonunits O
  have e : (1 : L) = piL⁻¹ * (-4 - 2 * piL⁻¹ + 4 * piL⁻¹ ^ 2 + 2 * piL⁻¹ ^ 3) := by
    have := piL_ne_zero
    field_simp
    linear_combination piL_eisenstein
  rw [e]
  refine mul_mem_nonunits_right O hi ?_
  have := (-4 - 2 * (⟨_, hiO⟩ : O) + 4 * (⟨_, hiO⟩ : O) ^ 2 + 2 * (⟨_, hiO⟩ : O) ^ 3).2
  push_cast at this
  exact this

lemma epsL_mem : epsL ∈ O := by
  have hπ := piL_mem O
  have := (43 - 2 * (⟨_, hπ⟩ : O) - 35 * (⟨_, hπ⟩ : O) ^ 2 - 10 * (⟨_, hπ⟩ : O) ^ 3).2
  push_cast at this
  exact this

variable (h2 : (2 : L) ∈ O.nonunits)
include h2

lemma piL_mem_nonunits : piL ∈ O.nonunits := by
  have hπ := piL_mem O
  have hu : piL ^ 4 = (1 + 2 * piL - piL ^ 2 - 2 * piL ^ 3) * 2 := by
    linear_combination piL_eisenstein
  refine mem_nonunits_of_pow_mem O hπ (n := 4) ?_
  rw [hu]
  refine mul_mem_nonunits_left O ?_ h2
  have := (1 + 2 * (⟨_, hπ⟩ : O) - (⟨_, hπ⟩ : O) ^ 2 - 2 * (⟨_, hπ⟩ : O) ^ 3).2
  push_cast at this
  exact this

/-- **Lemma A.** If `n₀ + n₁π + n₂π² + n₃π³ ∈ 2O`, then all `nᵢ` are even. -/
lemma even_coeffs {n₀ n₁ n₂ n₃ : ℤ} {b : L} (hb : b ∈ O)
    (h : (n₀ : L) + n₁ * piL + n₂ * piL ^ 2 + n₃ * piL ^ 3 = 2 * b) :
    Even n₀ ∧ Even n₁ ∧ Even n₂ ∧ Even n₃ := by
  have hπ := piL_mem_nonunits O h2
  have hε2 := two_eq_pi_pow_four_mul
  have hπ0 := piL_ne_zero
  -- an integer which is `π` times an element of `O` is even
  have key : ∀ (m : ℤ) (w : L), w ∈ O → (m : L) = piL * w → Even m := fun m w hw hm =>
    even_of_intCast_mem_nonunits O h2 (hm ▸ mul_mem_nonunits_right O hπ hw)
  set P : O := ⟨piL, piL_mem O⟩
  set E : O := ⟨epsL, epsL_mem O⟩
  set B : O := ⟨b, hb⟩
  have e0 : Even n₀ := by
    refine key n₀ (piL ^ 3 * epsL * b - n₁ - n₂ * piL - n₃ * piL ^ 2) ?_ ?_
    · have := (P ^ 3 * E * B - n₁ - n₂ * P - n₃ * P ^ 2).2
      push_cast at this
      exact this
    · linear_combination h + b * hε2
  obtain ⟨m₀, hm₀⟩ := e0
  have hm₀' : (n₀ : L) = m₀ + m₀ := by rw [hm₀]; push_cast; ring
  have e1 : Even n₁ := by
    refine key n₁ (piL ^ 2 * epsL * (b - m₀) - n₂ - n₃ * piL) ?_ ?_
    · have := (P ^ 2 * E * (B - m₀) - n₂ - n₃ * P).2
      push_cast at this
      exact this
    · apply mul_left_cancel₀ hπ0
      linear_combination h - hm₀' + (b - m₀) * hε2
  obtain ⟨m₁, hm₁⟩ := e1
  have hm₁' : (n₁ : L) = m₁ + m₁ := by rw [hm₁]; push_cast; ring
  have e2 : Even n₂ := by
    refine key n₂ (piL * epsL * (b - m₀ - m₁ * piL) - n₃) ?_ ?_
    · have := (P * E * (B - m₀ - m₁ * P) - n₃).2
      push_cast at this
      exact this
    · apply mul_left_cancel₀ (pow_ne_zero 2 hπ0)
      linear_combination h - hm₀' - piL * hm₁' + (b - m₀ - m₁ * piL) * hε2
  obtain ⟨m₂, hm₂⟩ := e2
  have hm₂' : (n₂ : L) = m₂ + m₂ := by rw [hm₂]; push_cast; ring
  have e3 : Even n₃ := by
    refine key n₃ (epsL * (b - m₀ - m₁ * piL - m₂ * piL ^ 2)) ?_ ?_
    · have := (E * (B - m₀ - m₁ * P - m₂ * P ^ 2)).2
      push_cast at this
      exact this
    · apply mul_left_cancel₀ (pow_ne_zero 3 hπ0)
      linear_combination h - hm₀' - piL * hm₁' - piL ^ 2 * hm₂'
        + (b - m₀ - m₁ * piL - m₂ * piL ^ 2) * hε2
  exact ⟨⟨m₀, hm₀⟩, ⟨m₁, hm₁⟩, ⟨m₂, hm₂⟩, e3⟩

/-- **Lemma B.** The residue field of `O` is `𝔽₂`: every `x ∈ O` is `≡ 0` or `≡ 1`. -/
lemma mem_nonunits_or_sub_one_mem {x : L} (hx : x ∈ O) :
    x ∈ O.nonunits ∨ x - 1 ∈ O.nonunits := by
  have hπO := piL_mem O
  have hπ := piL_mem_nonunits O h2
  obtain ⟨d, hd, n₀, n₁, n₂, n₃, h⟩ := exists_int_combo_pi x
  induction d using Nat.strong_induction_on generalizing n₀ n₁ n₂ n₃ with
  | _ d ih =>
  rcases Nat.even_or_odd d with ⟨d', rfl⟩ | ⟨d', rfl⟩
  · have hmem : (d' : L) * x ∈ O := O.mul_mem _ _ (natCast_mem O d') hx
    obtain ⟨⟨m₀, rfl⟩, ⟨m₁, rfl⟩, ⟨m₂, rfl⟩, ⟨m₃, rfl⟩⟩ :=
      even_coeffs O h2 hmem (by rw [← h]; push_cast; ring)
    refine ih d' (by omega) (by omega) m₀ m₁ m₂ m₃ ?_
    have e : (2 : L) * ((d' : L) * x) = 2 * (m₀ + m₁ * piL + m₂ * piL ^ 2 + m₃ * piL ^ 3) := by
      push_cast at h
      linear_combination h
    exact mul_left_cancel₀ two_ne_zero e
  · have hxn : x - n₀ ∈ O.nonunits := by
      have e : x - n₀ = piL * (n₁ + n₂ * piL + n₃ * piL ^ 2) - ((d' : L) * x) * 2 := by
        push_cast at h
        linear_combination h
      rw [e]
      refine sub_mem (mul_mem_nonunits_right O hπ ?_) (mul_mem_nonunits_left O ?_ h2)
      · have := ((n₁ : O) + n₂ * (⟨_, hπO⟩ : O) + n₃ * (⟨_, hπO⟩ : O) ^ 2).2
        push_cast at this
        exact this
      · exact O.mul_mem _ _ (natCast_mem O d') hx
    rcases Int.even_or_odd' n₀ with ⟨k, hk | hk⟩
    · left
      have e : x = (x - n₀) + (k : L) * 2 := by rw [hk]; push_cast; ring
      rw [e]
      exact add_mem hxn (mul_mem_nonunits_left O (intCast_mem O k) h2)
    · right
      have e : x - 1 = (x - n₀) + (k : L) * 2 := by rw [hk]; push_cast; ring
      rw [e]
      exact add_mem hxn (mul_mem_nonunits_left O (intCast_mem O k) h2)

end Residue

/-- **Upper bound.** The unit-distance graph of `ℚ(√2, √3)²` is 4-colourable. -/
theorem colorable_four : (unitDistGraph L).Colorable 4 := by
  obtain ⟨O, h2⟩ := exists_valuationSubring_two_mem_nonunits L
  exact colorable_four_of_residue_two s3 s3_sq O h2
    fun x hx => mem_nonunits_or_sub_one_mem O h2 hx

/-! ## Lower bound: a chain of three unit rhombi -/

/-- `a + b√2 + c√3 + d√6` as an element of `L`. -/
noncomputable def mkL (a b c d : ℚ) : L := a + b * s2 + c * s3 + d * (s2 * s3)

@[simp] lemma coe_mkL (a b c d : ℚ) :
    (mkL a b c d : ℝ) = a + b * √2 + c * √3 + d * (√2 * √3) := by
  simp [mkL]

/-! The ten vertices of `data/chain23.json`, on the basis `1, √2, √3, √6`. -/

noncomputable def p0 : L × L := (mkL 0 0 0 0, mkL 0 0 0 0)
noncomputable def p1 : L × L := (mkL 0 0 1 0, mkL 0 0 0 0)
noncomputable def p2 : L × L := (mkL 0 0 1 0, mkL 0 0 1 0)
noncomputable def p3 : L × L := (mkL 0 0 (1/3) (-1/6), mkL 0 0 (1/3) (1/6))
noncomputable def p4 : L × L := (mkL 0 0 (1/2) 0, mkL (1/2) 0 0 0)
noncomputable def p5 : L × L := (mkL 0 0 (1/2) 0, mkL (-1/2) 0 0 0)
noncomputable def p6 : L × L := (mkL (-1/2) 0 1 0, mkL 0 0 (1/2) 0)
noncomputable def p7 : L × L := (mkL (1/2) 0 1 0, mkL 0 0 (1/2) 0)
noncomputable def p8 : L × L := (mkL (1/3) (-1/12) (2/3) (-1/12), mkL (-1/3) (-1/12) (2/3) (1/12))
noncomputable def p9 : L × L := (mkL (-1/3) (1/12) (2/3) (-1/12), mkL (1/3) (1/12) (2/3) (1/12))

lemma unitDist_iff (p q : L × L) :
    (unitDistGraph L).Adj p q ↔ ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1 := Iff.rfl

/-! The sixteen unit distances. Each is a polynomial identity modulo `√2² = 2`, `√3² = 3`. -/

lemma adj03 : (unitDistGraph L).Adj p0 p3 := by
  rw [unitDist_iff]; simp only [p0, p3, coe_mkL]; push_cast
  linear_combination (√3 ^ 2 / 18) * sqrt2_sq + (1 / 3) * sqrt3_sq
lemma adj04 : (unitDistGraph L).Adj p0 p4 := by
  rw [unitDist_iff]; simp only [p0, p4, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj05 : (unitDistGraph L).Adj p0 p5 := by
  rw [unitDist_iff]; simp only [p0, p5, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj14 : (unitDistGraph L).Adj p1 p4 := by
  rw [unitDist_iff]; simp only [p1, p4, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj15 : (unitDistGraph L).Adj p1 p5 := by
  rw [unitDist_iff]; simp only [p1, p5, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj16 : (unitDistGraph L).Adj p1 p6 := by
  rw [unitDist_iff]; simp only [p1, p6, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj17 : (unitDistGraph L).Adj p1 p7 := by
  rw [unitDist_iff]; simp only [p1, p7, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj26 : (unitDistGraph L).Adj p2 p6 := by
  rw [unitDist_iff]; simp only [p2, p6, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj27 : (unitDistGraph L).Adj p2 p7 := by
  rw [unitDist_iff]; simp only [p2, p7, coe_mkL]; push_cast
  linear_combination (1 / 4) * sqrt3_sq
lemma adj28 : (unitDistGraph L).Adj p2 p8 := by
  rw [unitDist_iff]; simp only [p2, p8, coe_mkL]; push_cast
  linear_combination (√3 ^ 2 / 72 + 1 / 72) * sqrt2_sq + (1 / 4) * sqrt3_sq
lemma adj29 : (unitDistGraph L).Adj p2 p9 := by
  rw [unitDist_iff]; simp only [p2, p9, coe_mkL]; push_cast
  linear_combination (√3 ^ 2 / 72 + 1 / 72) * sqrt2_sq + (1 / 4) * sqrt3_sq
lemma adj38 : (unitDistGraph L).Adj p3 p8 := by
  rw [unitDist_iff]; simp only [p3, p8, coe_mkL]; push_cast
  linear_combination (√3 ^ 2 / 72 + 1 / 72) * sqrt2_sq + (1 / 4) * sqrt3_sq
lemma adj39 : (unitDistGraph L).Adj p3 p9 := by
  rw [unitDist_iff]; simp only [p3, p9, coe_mkL]; push_cast
  linear_combination (√3 ^ 2 / 72 + 1 / 72) * sqrt2_sq + (1 / 4) * sqrt3_sq
lemma adj45 : (unitDistGraph L).Adj p4 p5 := by
  rw [unitDist_iff]; simp only [p4, p5, coe_mkL]; push_cast
  ring
lemma adj67 : (unitDistGraph L).Adj p6 p7 := by
  rw [unitDist_iff]; simp only [p6, p7, coe_mkL]; push_cast
  ring
lemma adj89 : (unitDistGraph L).Adj p8 p9 := by
  rw [unitDist_iff]; simp only [p8, p9, coe_mkL]; push_cast
  linear_combination (1 / 18) * sqrt2_sq

/-- The edges of `data/chain23.json`. -/
def chainEdges : List (Fin 10 × Fin 10) :=
  [(0, 3), (0, 4), (0, 5), (1, 4), (1, 5), (1, 6), (1, 7), (2, 6), (2, 7), (2, 8), (2, 9),
   (3, 8), (3, 9), (4, 5), (6, 7), (8, 9)]

/-- The 10-vertex, 16-edge graph: unit rhombi `0-{4,5}-1`, `1-{6,7}-2`, `2-{8,9}-3`, closed by
the edge `0-3`. -/
abbrev chain23 : SimpleGraph (Fin 10) := edgeGraph chainEdges

/-- The chain of three rhombi is not 3-colourable. -/
theorem chain23_not_colorable : ¬ chain23.Colorable 3 := by
  rintro ⟨C⟩
  have v : ∀ i j, chain23.Adj i j → C i ≠ C j := fun i j h => C.valid h
  have r1 : C 0 = C 1 := fin3_rhombus _ _ _ _ (v 0 4 (by decide)) (v 0 5 (by decide))
    (v 4 5 (by decide)) (v 1 4 (by decide)) (v 1 5 (by decide))
  have r2 : C 1 = C 2 := fin3_rhombus _ _ _ _ (v 1 6 (by decide)) (v 1 7 (by decide))
    (v 6 7 (by decide)) (v 2 6 (by decide)) (v 2 7 (by decide))
  have r3 : C 2 = C 3 := fin3_rhombus _ _ _ _ (v 2 8 (by decide)) (v 2 9 (by decide))
    (v 8 9 (by decide)) (v 3 8 (by decide)) (v 3 9 (by decide))
  exact v 0 3 (by decide) (r1.trans (r2.trans r3))

/-- The positions of the ten vertices. -/
noncomputable def chainPt : Fin 10 → L × L := ![p0, p1, p2, p3, p4, p5, p6, p7, p8, p9]

lemma chainPt_adj : ∀ e ∈ chainEdges, (unitDistGraph L).Adj (chainPt e.1) (chainPt e.2) := by
  simp only [chainEdges, List.forall_mem_cons, List.not_mem_nil, IsEmpty.forall_iff,
    implies_true, and_true]
  exact ⟨adj03, adj04, adj05, adj14, adj15, adj16, adj17, adj26, adj27, adj28, adj29, adj38,
    adj39, adj45, adj67, adj89⟩

/-- **Lower bound.** The unit-distance graph of `ℚ(√2, √3)²` is not 3-colourable. -/
theorem not_colorable_three : ¬ (unitDistGraph L).Colorable 3 := fun h =>
  chain23_not_colorable (h.of_hom (edgeGraph.hom chainPt chainPt_adj))

/-- **Theorem.** The unit-distance graph of `ℚ(√2, √3)²` has chromatic number 4. -/
theorem chromaticNumber_eq_four : (unitDistGraph L).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of colorable_four not_colorable_three

end Q23
