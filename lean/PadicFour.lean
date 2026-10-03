import QuadraticPlanes
import FourColours

/-!
# Four colours for the `p`-adic planes, `p ≥ 5`

Theorem 1 of `notes/four_colours_11_mod_12.md` says that the unit-distance graph of `ℚ(√d)²` is not
3-colourable for every `d ≡ 11 (mod 12)`. It is stated here as the hypothesis `Theorem1` (it is being
formalised separately); this file proves the transfer, Corollary 3 of the note:

* `not_colorable_three_of_sq`: if a field `F` of characteristic 0 contains a square root `s` of some
  `d ≡ 11 (mod 12)`, then `QuadraticPlanes.sumSqGraph F` is not 3-colourable.
  `d` is not a square (`d ≡ 3 (mod 4)`), so `X² - d` is irreducible over `ℚ` and is the minimal polynomial of
  `√d`; as `s` is a root of it in `F`, there is a `ℚ`-algebra map `φ : ℚ(√d) → F`
  (`IntermediateField.algHomAdjoinIntegralEquiv`). The map `(x, y) ↦ (φ x, φ y)` is a graph homomorphism from
  the unit-distance graph of `ℚ(√d)²` to `sumSqGraph F`, so a 3-colouring would pull back.
* `padic_not_colorable_three`: for every prime `p ≥ 5`, `sumSqGraph ℚ_[p]` is not 3-colourable. By the Chinese
  remainder theorem there is `d` with `d ≡ 11 (mod 12)` and `d ≡ 1 (mod p)`; Hensel's lemma at `a = 1` for
  `X² - d` over `ℤ_[p]` (`f(1) = 1 - d ≡ 0`, `f'(1) = 2` a unit as `p` is odd) gives `s ∈ ℤ_[p]` with `s² = d`.

Both theorems take `Theorem1` as a hypothesis. `FourColours.not_colorable_three` proves it, and
`padic_not_colorable_three_unconditional` and `not_colorable_three_of_sq_unconditional` are the unconditional
forms.
-/

namespace PadicFour

open Polynomial QuadraticPlanes

/-- Theorem 1 of `notes/four_colours_11_mod_12.md`, as a hypothesis. -/
def Theorem1 : Prop := ∀ d : ℕ, d % 12 = 11 →
  ¬ (LocalColouring.unitDistGraph (QuadraticPlanes.L d)).Colorable 3

/-- A number `≡ 3 (mod 4)` is not a square. -/
lemma not_isSquare {d : ℕ} (hd : d % 4 = 3) : ¬ IsSquare d := by
  rintro ⟨m, rfl⟩
  rw [Nat.mul_mod] at hd
  have : m % 4 < 4 := Nat.mod_lt m (by norm_num)
  interval_cases h : m % 4 <;> simp at hd

/-- `X² - d` has no rational root. -/
lemma not_sq_rat {d : ℕ} (hd : d % 4 = 3) (q : ℚ) : q ^ 2 ≠ d := by
  intro hq
  have hirr : Irrational √(d : ℝ) := irrational_sqrt_natCast_iff.mpr (not_isSquare hd)
  apply hirr
  refine ⟨|q|, ?_⟩
  push_cast
  rw [← Real.sqrt_sq_eq_abs]
  congr 1
  exact_mod_cast hq

/-- The minimal polynomial of `√d` over `ℚ`. -/
lemma minpoly_sqrt {d : ℕ} (hd : d % 4 = 3) : minpoly ℚ √(d : ℝ) = X ^ 2 - C (d : ℚ) := by
  have hm : (X ^ 2 - C (d : ℚ)).Monic := monic_X_pow_sub_C _ (by norm_num)
  have hdeg : (X ^ 2 - C (d : ℚ)).natDegree = 2 := natDegree_X_pow_sub_C
  refine (minpoly.eq_of_irreducible_of_monic ?_ ?_ hm).symm
  · rw [hm.irreducible_iff_roots_eq_zero_of_degree_le_three (by rw [hdeg]) (by rw [hdeg]; norm_num)]
    refine Multiset.eq_zero_of_forall_notMem fun q hq => ?_
    rw [mem_roots hm.ne_zero, IsRoot, eval_sub, eval_pow, eval_X, eval_C, sub_eq_zero] at hq
    exact not_sq_rat hd q hq
  · simp [Real.sq_sqrt (Nat.cast_nonneg d)]

lemma isIntegral_sqrt (d : ℕ) : IsIntegral ℚ √(d : ℝ) := by
  refine ⟨X ^ 2 - C (d : ℚ), monic_X_pow_sub_C _ two_ne_zero, ?_⟩
  simp [Real.sq_sqrt (Nat.cast_nonneg d)]

/-- A `ℚ`-algebra map `ℚ(√d) → F` when `F` contains a square root of `d ≡ 3 (mod 4)`. -/
lemma nonempty_algHom {F : Type*} [Field F] [CharZero F] {d : ℕ} (hd : d % 4 = 3) {s : F}
    (hs : s ^ 2 = d) : Nonempty (L d →ₐ[ℚ] F) := by
  refine ⟨(IntermediateField.algHomAdjoinIntegralEquiv ℚ (isIntegral_sqrt d)).symm ⟨s, ?_⟩⟩
  rw [minpoly_sqrt hd, mem_aroots]
  refine ⟨(monic_X_pow_sub_C _ (by norm_num)).ne_zero, ?_⟩
  simp [hs]

/-- A ring map `ℚ(√d) → F` maps the unit-distance graph of `ℚ(√d)²` to `sumSqGraph F`. -/
def hom {F : Type*} [Field F] {d : ℕ} (φ : L d →+* F) :
    LocalColouring.unitDistGraph (L d) →g sumSqGraph F where
  toFun x := (φ x.1, φ x.2)
  map_rel' {x y} h := by
    have h : ((x.1 : ℝ) - y.1) ^ 2 + ((x.2 : ℝ) - y.2) ^ 2 = 1 := h
    have h' : (x.1 - y.1) ^ 2 + (x.2 - y.2) ^ 2 = (1 : L d) := by
      apply Subtype.ext
      push_cast
      exact h
    show (φ x.1 - φ y.1) ^ 2 + (φ x.2 - φ y.2) ^ 2 = 1
    simpa using congrArg φ h'

/-- Corollary 3, first part: a field of characteristic 0 containing a square root of some `d ≡ 11 (mod 12)`. -/
theorem not_colorable_three_of_sq (h : Theorem1) {F : Type*} [Field F] [CharZero F]
    {d : ℕ} (hd : d % 12 = 11) {s : F} (hs : s ^ 2 = d) : ¬ (sumSqGraph F).Colorable 3 := by
  rintro ⟨C⟩
  obtain ⟨φ⟩ := nonempty_algHom (by omega) hs
  exact h d hd ⟨C.comp (hom φ.toRingHom)⟩

/-- A natural number `d ≡ 11 (mod 12)` with `d ≡ 1 (mod p)`, for a prime `p ≥ 5`. -/
lemma exists_d (p : ℕ) [hp : Fact p.Prime] (h5 : 5 ≤ p) : ∃ d : ℕ, d % 12 = 11 ∧ d % p = 1 := by
  have hco : Nat.Coprime 12 p := by
    refine Nat.Coprime.symm ((Nat.Prime.coprime_iff_not_dvd hp.out).mpr fun h => ?_)
    have := Nat.le_of_dvd (by norm_num) h
    have hp' := hp.out
    interval_cases p <;> first | omega | exact absurd hp' (by norm_num)
  obtain ⟨d, h1, h2⟩ := Nat.chineseRemainder hco 11 1
  refine ⟨d, h1, ?_⟩
  rw [Nat.ModEq, Nat.mod_eq_of_lt (by omega : 1 < p)] at h2
  exact h2

/-- Hensel's lemma: a square root of `d` in `ℤ_[p]` when `d ≡ 1 (mod p)` and `p` is odd. -/
lemma exists_sq_padicInt (p : ℕ) [hp : Fact p.Prime] (h3 : 3 ≤ p) {d : ℕ} (hd : d % p = 1) :
    ∃ s : ℤ_[p], s ^ 2 = d := by
  have hlt : ‖((1 - d : ℤ) : ℤ_[p])‖ < 1 := by
    rw [PadicInt.norm_int_lt_one_iff_dvd]
    have hd1 : 1 ≤ d := by
      rcases Nat.eq_zero_or_pos d with h | h
      · subst h; simp at hd
      · exact h
    have : p ∣ d - 1 :=
      (Nat.modEq_iff_dvd' hd1).mp (by rw [Nat.ModEq, hd, Nat.mod_eq_of_lt (by omega)])
    rw [show (1 - d : ℤ) = -((d - 1 : ℕ) : ℤ) by push_cast [hd1]; ring]
    exact (Int.natCast_dvd_natCast.mpr this).neg_right
  have h2 : ‖((2 : ℤ) : ℤ_[p])‖ = 1 := by
    rw [PadicInt.norm_intCast_eq_one_iff, Int.isCoprime_iff_gcd_eq_one]
    show Int.gcd ((2 : ℕ) : ℤ) (p : ℤ) = 1
    rw [Int.gcd_natCast_natCast]
    exact (Nat.coprime_primes Nat.prime_two hp.out).mpr (by omega)
  have e1 : aeval (1 : ℤ_[p]) (X ^ 2 - C (d : ℤ_[p]) : ℤ_[p][X]) = ((1 - d : ℤ) : ℤ_[p]) := by
    simp
  have e2 : aeval (1 : ℤ_[p]) (derivative (X ^ 2 - C (d : ℤ_[p]) : ℤ_[p][X])) = ((2 : ℤ) : ℤ_[p]) := by
    simp; norm_num
  obtain ⟨z, hz, -⟩ := hensels_lemma (p := p) (R := ℤ_[p]) (F := (X ^ 2 - C (d : ℤ_[p]) : ℤ_[p][X]))
    (a := 1) (by rw [e1, e2, h2, one_pow]; exact hlt)
  exact ⟨z, by simpa [sub_eq_zero] using hz⟩

/-- Corollary 3: for every prime `p ≥ 5`, the plane over `ℚ_[p]` needs four colours (given Theorem 1). -/
theorem padic_not_colorable_three (h : Theorem1) (p : ℕ) [Fact p.Prime] (hp : 5 ≤ p) :
    ¬ (sumSqGraph ℚ_[p]).Colorable 3 := by
  obtain ⟨d, hd12, hdp⟩ := exists_d p hp
  obtain ⟨s, hs⟩ := exists_sq_padicInt p (by omega) hdp
  refine not_colorable_three_of_sq h hd12 (s := (s : ℚ_[p])) ?_
  rw [← PadicInt.coe_pow, hs]
  simp

/-- Corollary 3, unconditionally: the plane over `ℚ_[p]` needs four colours for every prime `p ≥ 5`
(`FourColours.not_colorable_three` is Theorem 1). -/
theorem padic_not_colorable_three_unconditional (p : ℕ) [Fact p.Prime] (hp : 5 ≤ p) :
    ¬ (sumSqGraph ℚ_[p]).Colorable 3 :=
  padic_not_colorable_three FourColours.not_colorable_three p hp

/-- The same for every field of characteristic 0 containing a square root of some `d ≡ 11 (mod 12)`. -/
theorem not_colorable_three_of_sq_unconditional {F : Type*} [Field F] [CharZero F]
    {d : ℕ} (hd : d % 12 = 11) {s : F} (hs : s ^ 2 = d) : ¬ (sumSqGraph F).Colorable 3 :=
  not_colorable_three_of_sq FourColours.not_colorable_three hd hs

end PadicFour
