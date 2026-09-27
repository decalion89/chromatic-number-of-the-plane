import Mathlib

/-!
# Colouring a plane over a real field through a place over 2

This file holds the part of the argument shared by `Q23` (the plane over `ℚ(√2, √3)`) and
`Q311` (the plane over `ℚ(√3, √11)`).

Let `K ⊆ ℝ` be a field containing `√3`. Its *unit-distance graph* has vertex set `K²`, with
`p ~ q` iff `(p₁ - q₁)² + (p₂ - q₂)² = 1`.

**Criterion** (`colorable_four_of_residue_two`). Suppose `K` has a valuation subring `O` with
`2` in its maximal ideal and residue field `𝔽₂`. Then the unit-distance graph of `K²` is
4-colourable.

*Proof.*
1. In the coordinates `a = x + y/√3`, `b = 2y/√3` the squared distance `x² + y²` becomes
   `a² - ab + b²` (`toEisHom`).
2. Over `𝔽₂` the form `a² + ab + b²` has no nontrivial zero. Hence if `a² - ab + b² = 1`, then
   `a, b ∈ O` and `(a, b)` does not reduce to `(0, 0)` (`core_lemma`). Otherwise divide by
   `s²`, where `s` is whichever of `a, b` has the larger absolute value for `O`; reducing gives
   a zero of `1 + t + t²` over `𝔽₂`, which has none.
3. Choose a representative `rep z` of each coset `z + O²` and colour `z` by the residue of
   `z - rep z` in `𝔽₂²` (`qGraph_colorable`). Two points at unit distance lie in the same coset
   and differ by a vector with nonzero residue, so they get different colours.

The file also provides the tools used on both fields:
* the existence of a valuation subring with `2` in its maximal ideal, from Chevalley's
  extension theorem in Mathlib;
* elementary facts about such subrings;
* every element of `ℚ(u, v)` with `u², v² ∈ ℤ` is `(n₀ + n₁u + n₂v + n₃uv)/d`;
* graphs on `Fin n` given by an edge list, for the lower bounds.
-/

namespace LocalColouring

open IsLocalRing

/-! ## The unit-distance graph -/

/-- The unit-distance graph on `K²`, for a subfield `K` of `ℝ`:
`p ~ q` iff `(p₁ - q₁)² + (p₂ - q₂)² = 1`, computed in `ℝ`. -/
def unitDistGraph (K : IntermediateField ℚ ℝ) : SimpleGraph (K × K) where
  Adj p q := ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1
  symm := ⟨fun p q h => by linear_combination h⟩
  loopless := ⟨fun p h => by simp at h⟩

/-- A graph is 4-chromatic once it is 4-colourable and not 3-colourable. -/
lemma chromaticNumber_eq_four_of {V : Type*} {G : SimpleGraph V} (h4 : G.Colorable 4)
    (h3 : ¬ G.Colorable 3) : G.chromaticNumber = 4 := by
  rw [show (4 : ℕ∞) = ((3 : ℕ) : ℕ∞) + 1 by norm_num]
  exact SimpleGraph.chromaticNumber_eq_iff_colorable_not_colorable.mpr ⟨h4, h3⟩

/-! ## The core lemma: unit vectors of `a² - ab + b²` are integral with nonzero residue -/

section Core

variable {F : Type*} [Field F] (O : ValuationSubring F) (red : O →+* ZMod 2)

lemma zmod2_one_sub_add_sq (r : ZMod 2) : 1 - r + r ^ 2 = 1 := by
  fin_cases r <;> decide

/-- Over `𝔽₂` the form `1 - t + t²` never vanishes, so it is not the square of an element
of the maximal ideal. -/
lemma one_sub_add_sq_ne (hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0) (t x : O)
    (hx : x ∈ maximalIdeal O) : 1 - t + t ^ 2 ≠ x ^ 2 := by
  intro h
  have h' := congrArg red h
  rw [map_add, map_sub, map_one, map_pow, map_pow, hred x hx, zmod2_one_sub_add_sq] at h'
  exact absurd h' (by decide)

/-- If `a² - ab + b² = 1`, then `a ∈ O`. -/
lemma mem_of_form_eq_one (hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0) {a b : F}
    (h : a ^ 2 - a * b + b ^ 2 = 1) : a ∈ O := by
  by_contra ha
  have ha0 : a ≠ 0 := fun h0 => ha (h0 ▸ O.zero_mem)
  have hx : a⁻¹ ∈ O.nonunits := O.inv_mem_nonunits_iff.mpr (Or.inr ha)
  obtain ⟨hxO, hxm⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp hx
  by_cases ht : b / a ∈ O
  · -- `t = b / a` is in `O` and `a⁻¹` is in `𝔪`
    apply one_sub_add_sq_ne O red hred ⟨b / a, ht⟩ ⟨a⁻¹, hxO⟩ hxm
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h
  · -- otherwise `s = a / b` is in `O`, and `b⁻¹ = s a⁻¹` is in `𝔪`
    have hs : a / b ∈ O := by
      have := (O.mem_or_inv_mem (b / a)).resolve_left ht
      rwa [inv_div] at this
    have hb0 : b ≠ 0 := by
      rintro rfl
      exact ht (by simp)
    apply one_sub_add_sq_ne O red hred ⟨a / b, hs⟩ (⟨a / b, hs⟩ * ⟨a⁻¹, hxO⟩)
      (Ideal.mul_mem_left _ _ hxm)
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h

/-- **Core lemma.** Let `O` be a valuation subring of `F` and `red : O →+* ZMod 2` have kernel
the maximal ideal. If `a² - ab + b² = 1`, then `a, b ∈ O` and `(red a, red b) ≠ (0, 0)`. -/
theorem core_lemma (hker : RingHom.ker red = maximalIdeal O) {a b : F}
    (h : a ^ 2 - a * b + b ^ 2 = 1) :
    ∃ (ha : a ∈ O) (hb : b ∈ O), (red ⟨a, ha⟩, red ⟨b, hb⟩) ≠ (0, 0) := by
  have hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0 := fun x hx => by
    rw [← hker] at hx
    exact hx
  have ha := mem_of_form_eq_one O red hred h
  have hb := mem_of_form_eq_one O red hred (a := b) (b := a) (by linear_combination h)
  refine ⟨ha, hb, ?_⟩
  intro h0
  simp only [Prod.mk.injEq] at h0
  have e : (⟨a, ha⟩ : O) ^ 2 - ⟨a, ha⟩ * ⟨b, hb⟩ + ⟨b, hb⟩ ^ 2 = 1 := Subtype.ext (by simpa using h)
  have e' := congrArg red e
  rw [map_add, map_sub, map_mul, map_pow, map_pow, h0.1, h0.2, map_one] at e'
  exact absurd e' (by decide)

end Core

/-! ## The colouring by residues of coset representatives -/

section Colouring

variable {F : Type*} [Field F]

/-- The graph on `F × F` of the norm form `a² - ab + b²`: `p ~ q` iff `Q(p - q) = 1`. -/
def qGraph (F : Type*) [Field F] : SimpleGraph (F × F) where
  Adj p q := (p.1 - q.1) ^ 2 - (p.1 - q.1) * (p.2 - q.2) + (p.2 - q.2) ^ 2 = 1
  symm := ⟨fun p q h => by linear_combination h⟩
  loopless := ⟨fun p h => by simp at h⟩

variable (O : ValuationSubring F) (red : O →+* ZMod 2)

open Classical in
/-- `red` extended by `0` outside `O`. -/
noncomputable def redExt (x : F) : ZMod 2 := if h : x ∈ O then red ⟨x, h⟩ else 0

lemma redExt_add {x y : F} (hx : x ∈ O) (hy : y ∈ O) :
    redExt O red (x + y) = redExt O red x + redExt O red y := by
  simp only [redExt, hx, hy, add_mem hx hy, dite_true]
  rw [← map_add]
  rfl

/-- A representative of the coset `x + O`. -/
noncomputable def rep (x : F) : F :=
  (QuotientAddGroup.mk x : F ⧸ O.toSubring.toAddSubgroup).out

lemma sub_rep_mem (x : F) : x - rep O x ∈ O := by
  obtain ⟨h, hh⟩ := QuotientAddGroup.mk_out_eq_add (s := O.toSubring.toAddSubgroup) x
  have : rep O x = x + h := hh
  rw [this, sub_add_cancel_left]
  exact O.neg_mem _ h.2

lemma rep_eq_of_sub_mem {x y : F} (h : x - y ∈ O) : rep O x = rep O y := by
  unfold rep
  congr 1
  rw [QuotientAddGroup.eq]
  have : -x + y = -(x - y) := by ring
  rw [this]
  exact O.neg_mem _ h

/-- The colouring: the residues of `z - rep z`, coordinatewise. -/
noncomputable def colour (z : F × F) : ZMod 2 × ZMod 2 :=
  (redExt O red (z.1 - rep O z.1), redExt O red (z.2 - rep O z.2))

lemma colour_adj (hker : RingHom.ker red = maximalIdeal O) {p q : F × F}
    (hpq : (qGraph F).Adj p q) : colour O red p ≠ colour O red q := by
  obtain ⟨ha, hb, hne⟩ := core_lemma O red hker hpq
  intro hc
  apply hne
  simp only [colour, Prod.mk.injEq] at hc
  have e1 : p.1 - rep O p.1 = (p.1 - q.1) + (q.1 - rep O q.1) := by
    rw [rep_eq_of_sub_mem O ha]; ring
  have e2 : p.2 - rep O p.2 = (p.2 - q.2) + (q.2 - rep O q.2) := by
    rw [rep_eq_of_sub_mem O hb]; ring
  rw [e1, redExt_add O red ha (sub_rep_mem O _)] at hc
  rw [e2, redExt_add O red hb (sub_rep_mem O _)] at hc
  have h1 : redExt O red (p.1 - q.1) = 0 := by linear_combination hc.1
  have h2 : redExt O red (p.2 - q.2) = 0 := by linear_combination hc.2
  simp only [redExt, ha, hb, dite_true] at h1 h2
  rw [h1, h2]

/-- The graph of `a² - ab + b²` is 4-colourable, with colours `𝔽₂²`. -/
theorem qGraph_colorable (hker : RingHom.ker red = maximalIdeal O) : (qGraph F).Colorable 4 := by
  have hc : (qGraph F).Colorable (Fintype.card (ZMod 2 × ZMod 2)) :=
    (SimpleGraph.Coloring.mk (colour O red) (fun h => colour_adj O red hker h)).colorable
  simpa using hc

end Colouring

/-! ## Eisenstein coordinates -/

/-- The coordinates `a = x + y/√3`, `b = 2y/√3`, where `r = √3` (note `1/√3 = √3/3`). -/
noncomputable def toEis {K : IntermediateField ℚ ℝ} (r : K) (p : K × K) : K × K :=
  (p.1 + p.2 * r / 3, 2 * p.2 * r / 3)

/-- In Eisenstein coordinates `x² + y² = a² - ab + b²`, so the change of coordinates maps the
unit-distance graph to the graph of the form `a² - ab + b²`. -/
noncomputable def toEisHom {K : IntermediateField ℚ ℝ} (r : K) (hr : r ^ 2 = 3) :
    unitDistGraph K →g qGraph K where
  toFun := toEis r
  map_rel' := by
    intro p q h
    have h' : (p.1 - q.1) ^ 2 + (p.2 - q.2) ^ 2 = (1 : K) := Subtype.ext (by push_cast; exact h)
    change ((p.1 + p.2 * r / 3) - (q.1 + q.2 * r / 3)) ^ 2
      - ((p.1 + p.2 * r / 3) - (q.1 + q.2 * r / 3)) * (2 * p.2 * r / 3 - 2 * q.2 * r / 3)
      + (2 * p.2 * r / 3 - 2 * q.2 * r / 3) ^ 2 = 1
    linear_combination h' + ((p.2 - q.2) ^ 2 / 3) * hr

/-! ## Valuation subrings with `2` in the maximal ideal -/

section Valuation

variable {F : Type*} [Field F] (O : ValuationSubring F)

/-- Every field of characteristic zero has a valuation subring with `2` in its maximal ideal
(Chevalley's extension theorem, applied to `ℤ` and the ideal `(2)`). -/
theorem exists_valuationSubring_two_mem_nonunits (F : Type*) [Field F] [CharZero F] :
    ∃ O : ValuationSubring F, (2 : F) ∈ O.nonunits := by
  have hne : Ideal.span {(2 : (⊥ : Subring F))} ≠ ⊤ := by
    rw [Ne, Ideal.span_singleton_eq_top]
    rintro ⟨u, hu⟩
    have hc2 : ((2 : (⊥ : Subring F)) : F) = 2 := rfl
    obtain ⟨n, hn⟩ := Subring.mem_bot.mp (u⁻¹ : (⊥ : Subring F)ˣ).1.2
    have h1 : ((u : (⊥ : Subring F)) : F) * ((u⁻¹ : (⊥ : Subring F)ˣ) : F) = 1 := by
      rw [← Subring.coe_mul, Units.mul_inv]; rfl
    rw [hu, ← hn, hc2] at h1
    have h2 : ((2 * n : ℤ) : F) = ((1 : ℤ) : F) := by push_cast; exact h1
    have h3 : (2 * n : ℤ) = 1 := Int.cast_injective h2
    omega
  obtain ⟨O, -, hO⟩ :=
    Ideal.image_subset_nonunits_valuationSubring (A := (⊥ : Subring F)) _ hne
  exact ⟨O, hO ⟨2, Ideal.subset_span rfl, rfl⟩⟩

lemma mul_mem_nonunits_left {x y : F} (hx : x ∈ O) (hy : y ∈ O.nonunits) :
    x * y ∈ O.nonunits := by
  obtain ⟨hy', hy''⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp hy
  exact O.coe_mem_nonunits_iff.mpr (Ideal.mul_mem_left _ (⟨x, hx⟩ : O) hy'')

lemma mul_mem_nonunits_right {x y : F} (hx : x ∈ O.nonunits) (hy : y ∈ O) :
    x * y ∈ O.nonunits := by
  rw [mul_comm]; exact mul_mem_nonunits_left O hy hx

lemma one_notMem_nonunits : (1 : F) ∉ O.nonunits := by
  rw [O.mem_nonunits_iff, map_one]; exact lt_irrefl 1

/-- The maximal ideal is prime. -/
lemma mem_nonunits_or_of_mul_mem {x y : F} (hx : x ∈ O) (hy : y ∈ O)
    (h : x * y ∈ O.nonunits) : x ∈ O.nonunits ∨ y ∈ O.nonunits := by
  obtain ⟨hxy, hm⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp h
  rcases (maximalIdeal.isMaximal O).isPrime.mem_or_mem (x := ⟨x, hx⟩) (y := ⟨y, hy⟩) hm with h | h
  · exact Or.inl (O.coe_mem_nonunits_iff.mpr h)
  · exact Or.inr (O.coe_mem_nonunits_iff.mpr h)

/-- The maximal ideal is radical. -/
lemma mem_nonunits_of_pow_mem {x : F} (hx : x ∈ O) {n : ℕ} (h : x ^ n ∈ O.nonunits) :
    x ∈ O.nonunits := by
  have hm : (⟨x, hx⟩ : O) ^ n ∈ maximalIdeal O := O.coe_mem_nonunits_iff.mp (by simpa using h)
  exact O.coe_mem_nonunits_iff.mpr ((maximalIdeal.isMaximal O).isPrime.mem_of_pow_mem n hm)

/-- Units of `O` have their inverse in `O`. -/
lemma inv_mem_of_notMem_nonunits {x : F} (hx : x ∉ O.nonunits) : x⁻¹ ∈ O := by
  rw [O.mem_nonunits_iff_or, not_or, not_not] at hx
  exact hx.2

/-- A root of a monic quadratic over `O` lies in `O`. -/
lemma mem_of_sq_eq {x a b : F} (ha : a ∈ O) (hb : b ∈ O) (h : x ^ 2 = a * x + b) : x ∈ O := by
  by_contra hx
  have hx0 : x ≠ 0 := fun h0 => hx (h0 ▸ O.zero_mem)
  have hi : x⁻¹ ∈ O.nonunits := O.inv_mem_nonunits_iff.mpr (Or.inr hx)
  have hiO : x⁻¹ ∈ O := O.nonunits_subset hi
  apply one_notMem_nonunits O
  have e : (1 : F) = x⁻¹ * (a + b * x⁻¹) := by
    field_simp
    linear_combination h
  rw [e]
  exact mul_mem_nonunits_right O hi (add_mem ha (O.mul_mem _ _ hb hiO))

variable (h2 : (2 : F) ∈ O.nonunits)
include h2

/-- An integer in the maximal ideal is even. -/
lemma even_of_intCast_mem_nonunits {n : ℤ} (hn : (n : F) ∈ O.nonunits) : Even n := by
  by_contra hne
  obtain ⟨k, rfl⟩ := Int.not_even_iff_odd.mp hne
  apply one_notMem_nonunits O
  have e : (1 : F) = ((2 * k + 1 : ℤ) : F) - (k : F) * 2 := by push_cast; ring
  rw [e]
  exact sub_mem hn (mul_mem_nonunits_left O (intCast_mem O k) h2)

/-- If an integer is `2ᵏ` times an element of `O`, then `2ᵏ` divides it. -/
lemma two_pow_dvd_of_eq [CharZero F] {k : ℕ} {n : ℤ} {γ : F} (hγ : γ ∈ O)
    (h : (n : F) = 2 ^ k * γ) :
    (2 : ℤ) ^ k ∣ n := by
  induction k generalizing n γ with
  | zero => simp
  | succ k ih =>
    have hn : Even n := even_of_intCast_mem_nonunits O h2 (by
      rw [h, pow_succ, mul_comm _ (2 : F), mul_assoc]
      exact mul_mem_nonunits_right O h2 (O.mul_mem _ _ (pow_mem (natCast_mem O 2) k) hγ))
    obtain ⟨m, rfl⟩ := hn
    have hm : (m : F) = 2 ^ k * γ := by
      apply mul_left_cancel₀ (two_ne_zero (α := F))
      push_cast at h
      linear_combination h
    obtain ⟨c, hc⟩ := ih hγ hm
    exact ⟨c, by rw [pow_succ]; linear_combination 2 * hc⟩

end Valuation

/-- A local ring with residue field `𝔽₂` (every element is `≡ 0` or `≡ 1`, and `2 ≡ 0`)
has a reduction map to `ZMod 2` with kernel the maximal ideal. -/
theorem exists_red_of_residue_two {R : Type*} [CommRing R] [IsLocalRing R]
    (h2 : (2 : R) ∈ maximalIdeal R)
    (h : ∀ x : R, x ∈ maximalIdeal R ∨ x - 1 ∈ maximalIdeal R) :
    ∃ red : R →+* ZMod 2, RingHom.ker red = maximalIdeal R := by
  classical
  have hm := (maximalIdeal.isMaximal R).isPrime
  let f : R → ZMod 2 := fun x => if x ∈ maximalIdeal R then 0 else 1
  have hf0 : ∀ z, z ∈ maximalIdeal R → f z = 0 := fun z hz => by simp only [f, hz, ite_true]
  have hf1 : ∀ z, z ∉ maximalIdeal R → f z = 1 := fun z hz => by simp only [f, hz, ite_false]
  have h1 : (1 : R) ∉ maximalIdeal R := (maximalIdeal R).ne_top_iff_one.mp hm.ne_top
  have hmul : ∀ x y, f (x * y) = f x * f y := by
    intro x y
    by_cases hx : x ∈ maximalIdeal R
    · rw [hf0 _ hx, hf0 _ (Ideal.mul_mem_right _ _ hx), zero_mul]
    by_cases hy : y ∈ maximalIdeal R
    · rw [hf0 _ hy, hf0 _ (Ideal.mul_mem_left _ _ hy), mul_zero]
    have hxy : x * y ∉ maximalIdeal R := fun hxy => (hm.mem_or_mem hxy).elim hx hy
    rw [hf1 _ hx, hf1 _ hy, hf1 _ hxy, mul_one]
  have hadd : ∀ x y, f (x + y) = f x + f y := by
    intro x y
    by_cases hx : x ∈ maximalIdeal R <;> by_cases hy : y ∈ maximalIdeal R
    · rw [hf0 _ hx, hf0 _ hy, hf0 _ (Ideal.add_mem _ hx hy), add_zero]
    · have hxy : x + y ∉ maximalIdeal R := fun hxy => hy (by
        have := Ideal.sub_mem _ hxy hx
        rwa [add_sub_cancel_left] at this)
      rw [hf0 _ hx, hf1 _ hy, hf1 _ hxy, zero_add]
    · have hxy : x + y ∉ maximalIdeal R := fun hxy => hx (by
        have := Ideal.sub_mem _ hxy hy
        rwa [add_sub_cancel_right] at this)
      rw [hf1 _ hx, hf0 _ hy, hf1 _ hxy, add_zero]
    · have hx1 := (h x).resolve_left hx
      have hy1 := (h y).resolve_left hy
      have hxy : x + y ∈ maximalIdeal R := by
        have := Ideal.add_mem _ (Ideal.add_mem _ hx1 hy1) h2
        convert this using 1
        ring
      rw [hf1 _ hx, hf1 _ hy, hf0 _ hxy]
      decide
  refine ⟨⟨⟨⟨f, hf1 1 h1⟩, hmul⟩, hf0 0 (maximalIdeal R).zero_mem, hadd⟩, ?_⟩
  ext x
  change f x = 0 ↔ x ∈ maximalIdeal R
  by_cases hx : x ∈ maximalIdeal R
  · simp only [hf0 _ hx, hx]
  · simp only [hf1 _ hx, hx, iff_false]
    decide

/-- **Criterion.** If `√3 ∈ K` and `K` has a valuation subring with `2` in its maximal ideal and
residue field `𝔽₂`, then the unit-distance graph of `K²` is 4-colourable. -/
theorem colorable_four_of_residue_two {K : IntermediateField ℚ ℝ} (r : K) (hr : r ^ 2 = 3)
    (O : ValuationSubring K) (h2 : (2 : K) ∈ O.nonunits)
    (hres : ∀ x ∈ O, x ∈ O.nonunits ∨ x - 1 ∈ O.nonunits) : (unitDistGraph K).Colorable 4 := by
  have hB : ∀ x : O, x ∈ maximalIdeal O ∨ x - 1 ∈ maximalIdeal O := by
    intro x
    rcases hres x x.2 with h | h
    · exact Or.inl (O.coe_mem_nonunits_iff.mp h)
    · exact Or.inr (O.coe_mem_nonunits_iff.mp (by simpa using h))
  obtain ⟨red, hker⟩ := exists_red_of_residue_two (O.coe_mem_nonunits_iff.mp h2) hB
  exact (qGraph_colorable O red hker).of_hom (toEisHom r hr)

/-! ## Elements of `ℚ(u, v)` -/

lemma isAlgebraic_of_sq_eq_intCast {u : ℝ} {m : ℤ} (hu : u ^ 2 = m) : IsAlgebraic ℚ u := by
  refine IsAlgebraic.of_pow (n := 2) two_pos ?_
  rw [hu, show (m : ℝ) = algebraMap ℚ ℝ m by simp]
  exact isAlgebraic_algebraMap _

/-- If `u² = m` and `v² = n` are integers, every element of `ℚ(u, v) ⊆ ℝ` is
`(n₀ + n₁u + n₂v + n₃uv) / d` with integers `nᵢ` and `d > 0`. -/
lemma exists_int_combo {u v : ℝ} {m n : ℤ} (hu : u ^ 2 = m) (hv : v ^ 2 = n) {y : ℝ}
    (hy : y ∈ IntermediateField.adjoin ℚ ({u, v} : Set ℝ)) :
    ∃ d : ℕ, 0 < d ∧ ∃ n₀ n₁ n₂ n₃ : ℤ,
      (d : ℝ) * y = n₀ + n₁ * u + n₂ * v + n₃ * (u * v) := by
  have halg : ∀ x ∈ ({u, v} : Set ℝ), IsAlgebraic ℚ x := by
    intro x hx
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
    rcases hx with rfl | rfl
    · exact isAlgebraic_of_sq_eq_intCast hu
    · exact isAlgebraic_of_sq_eq_intCast hv
  have hy' : y ∈ Algebra.adjoin ℚ ({u, v} : Set ℝ) := by
    rw [← IntermediateField.adjoin_toSubalgebra_of_isAlgebraic halg]
    exact hy
  clear hy halg
  induction hy' using Algebra.adjoin_induction with
  | mem x hx =>
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
    rcases hx with rfl | rfl
    · exact ⟨1, one_pos, 0, 1, 0, 0, by push_cast; ring⟩
    · exact ⟨1, one_pos, 0, 0, 1, 0, by push_cast; ring⟩
  | algebraMap r =>
    refine ⟨r.den, r.pos, r.num, 0, 0, 0, ?_⟩
    rw [eq_ratCast]
    simp only [Int.cast_zero, zero_mul, add_zero]
    exact_mod_cast Rat.den_mul_eq_num r
  | add x y _ _ hx hy =>
    obtain ⟨d₁, hd₁, a₁, b₁, c₁, e₁, h₁⟩ := hx
    obtain ⟨d₂, hd₂, a₂, b₂, c₂, e₂, h₂⟩ := hy
    refine ⟨d₁ * d₂, Nat.mul_pos hd₁ hd₂, d₂ * a₁ + d₁ * a₂, d₂ * b₁ + d₁ * b₂,
      d₂ * c₁ + d₁ * c₂, d₂ * e₁ + d₁ * e₂, ?_⟩
    push_cast
    linear_combination (d₂ : ℝ) * h₁ + (d₁ : ℝ) * h₂
  | mul x y _ _ hx hy =>
    obtain ⟨d₁, hd₁, a₁, b₁, c₁, e₁, h₁⟩ := hx
    obtain ⟨d₂, hd₂, a₂, b₂, c₂, e₂, h₂⟩ := hy
    refine ⟨d₁ * d₂, Nat.mul_pos hd₁ hd₂, a₁ * a₂ + m * b₁ * b₂ + n * c₁ * c₂ + m * n * e₁ * e₂,
      a₁ * b₂ + b₁ * a₂ + n * (c₁ * e₂ + e₁ * c₂), a₁ * c₂ + c₁ * a₂ + m * (b₁ * e₂ + e₁ * b₂),
      a₁ * e₂ + e₁ * a₂ + b₁ * c₂ + c₁ * b₂, ?_⟩
    push_cast
    linear_combination ((d₂ : ℝ) * y) * h₁
      + ((a₁ : ℝ) + b₁ * u + c₁ * v + e₁ * (u * v)) * h₂
      + ((b₁ : ℝ) * b₂ + e₁ * e₂ * v ^ 2 + v * (b₁ * e₂ + b₂ * e₁)) * hu
      + ((c₁ : ℝ) * c₂ + m * e₁ * e₂ + u * (c₁ * e₂ + c₂ * e₁)) * hv

/-! ## Finite graphs from edge lists -/

/-- In a 3-colouring, the two tips of a rhombus (two triangles sharing an edge) get the
same colour. -/
lemma fin3_rhombus (a b c d : Fin 3) (h₁ : a ≠ b) (h₂ : a ≠ c) (h₃ : b ≠ c) (h₄ : d ≠ b)
    (h₅ : d ≠ c) : a = d := by
  omega

/-- The graph on `Fin n` with the edges listed in `E`. -/
def edgeGraph {n : ℕ} (E : List (Fin n × Fin n)) : SimpleGraph (Fin n) :=
  SimpleGraph.fromRel fun i j => (i, j) ∈ E

instance {n : ℕ} (E : List (Fin n × Fin n)) : DecidableRel (edgeGraph E).Adj := fun i j =>
  decidable_of_iff (i ≠ j ∧ ((i, j) ∈ E ∨ (j, i) ∈ E)) (by simp [edgeGraph])

/-- A map sending each listed edge to an edge of `G` is a graph homomorphism. -/
def edgeGraph.hom {n : ℕ} {E : List (Fin n × Fin n)} {V : Type*} {G : SimpleGraph V}
    (f : Fin n → V) (hf : ∀ e ∈ E, G.Adj (f e.1) (f e.2)) : edgeGraph E →g G where
  toFun := f
  map_rel' := by
    intro i j h
    rw [edgeGraph, SimpleGraph.fromRel_adj] at h
    obtain ⟨-, h | h⟩ := h
    · exact hf _ h
    · exact (hf _ h).symm

end LocalColouring
