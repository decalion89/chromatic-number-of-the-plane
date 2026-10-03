import TheoremW
import Mathlib.Analysis.Normed.Group.AddCircle
import Mathlib.Algebra.Category.Grp.Injective
import Mathlib.Data.Fintype.Pigeonhole
import Mathlib.Topology.Compactness.Compact
import Mathlib.NumberTheory.DiophantineApproximation.Basic

/-!
# Chromatic recurrence for 3-colourings (Question 3 of Glasscock–Koutsogiannis–Richter)

We answer Question 3 of Glasscock, Koutsogiannis and Richter (Bull. Amer. Math. Soc. 59 (2022),
569–606) as a consequence of Theorem W (`TheoremW.theoremW_finite`):

`GKR.question3`: for every 3-colouring `c : ℕ → Fin 3` there is `α : ℝ` such that every
`n ≥ 1` with `‖n α‖ < 1/3` (the distance from `n α` to the nearest integer, written
`|n α - round (n α)|`) is a difference of two numbers of the same colour: `c m = c (m + n)` for
some `m`. In other words, the Bohr set `{n ≥ 1 : ‖n α‖ < 1/3}` consists of monochromatic
differences of `c`.

Equivalent form (`GKR.chromatic_recurrence`): every set `S ⊆ ℤ \ {0}` that contains, for every
`α`, some `s` with `‖s α‖ < 1/3` is chromatically recurrent for 3-colourings of `ℤ` (some
`s ∈ S` is a difference of two integers of the same colour).

Proof. Let `T` be the set of differences `t ≥ 1` that never occur inside a colour class.
1. (`finite_step`) Fix `N`. By pigeonhole two windows `[i, i + N)` and `[i + P, i + P + N)`
   (`P > 0`) carry the same colour pattern, so `c` is `P`-periodic on `[i, i + P + N)`
   (`periodic_of_window`). The induced colouring `x ↦ c (i + x.val)` of `ZMod P` is proper for
   the steps `±t`, `t ∈ T`, `t ≤ N`. Theorem W (`TheoremW.theoremW_finite`) gives a character `ξ` of
   the subgroup they generate with `ξ t ∈ [1/3, 2/3] (mod 1)`; since `ℝ/ℤ` is divisible (hence an
   injective `ℤ`-module) `ξ` extends to `F : ZMod P →+ ℝ/ℤ`, and `θ = F 1` satisfies
   `‖t θ‖ ≥ 1/3` for `t ∈ T`, `t ≤ N`.
2. (`exists_theta`) The sets `{θ ∈ ℝ/ℤ : ‖t θ‖ ≥ 1/3}`, `t ∈ T`, are closed and have the finite
   intersection property, so by compactness of `ℝ/ℤ` they have a common point `θ = α mod 1`.
-/

set_option autoImplicit false

namespace GKR

/-- On `[1/3, 2/3]` the distance to the nearest integer is at least `1/3`. -/
lemma third_le_abs_sub_round {x : ℝ} (h1 : 1/3 ≤ x) (h2 : x ≤ 2/3) :
    1/3 ≤ |x - round x| := by
  rcases le_or_gt (round x) 0 with h | h
  · have h' : (round x : ℝ) ≤ 0 := by exact_mod_cast h
    rw [abs_of_nonneg (by linarith)]
    linarith
  · have h'' : (1 : ℤ) ≤ round x := by omega
    have h' : (1 : ℝ) ≤ round x := by exact_mod_cast h''
    rw [abs_of_nonpos (by linarith)]
    linarith

/-- A colouring that agrees on the windows `[i, i + N)` and `[i + P, i + P + N)` is
`P`-periodic on `[i, i + P + N)`. -/
lemma periodic_of_window (c : ℕ → Fin 3) {i P N : ℕ} (hP : 0 < P)
    (hw : ∀ k < N, c (i + k) = c (i + P + k)) :
    ∀ m, m < P + N → c (i + m) = c (i + m % P) := by
  intro m
  induction m using Nat.strong_induction_on with
  | _ m ih =>
    intro hm
    rcases Nat.lt_or_ge m P with h | h
    · rw [Nat.mod_eq_of_lt h]
    · have h1 := hw (m - P) (by omega)
      rw [show i + P + (m - P) = i + m by omega] at h1
      rw [← h1, ih (m - P) (by omega) (by omega), ← Nat.mod_eq_sub_mod h]

/-- Step 1 (finite version, from Theorem W): for every `N` there is `θ ∈ ℝ/ℤ` with
`‖t • θ‖ ≥ 1/3` for every `t ∈ [1, N]` that never occurs as a difference inside a colour
class. -/
lemma finite_step (c : ℕ → Fin 3) (N : ℕ) :
    ∃ θ : UnitAddCircle, ∀ t : ℕ, 0 < t → t ≤ N → (∀ m, c m ≠ c (m + t)) →
      1/3 ≤ ‖t • θ‖ := by
  classical
  have : Fact ((0 : ℝ) < 1) := ⟨one_pos⟩
  -- pigeonhole on the colour patterns of the windows of length `N`
  obtain ⟨i, P, hP, hw⟩ : ∃ i P : ℕ, 0 < P ∧ ∀ k < N, c (i + k) = c (i + P + k) := by
    obtain ⟨a, b, hab, he⟩ :=
      Finite.exists_ne_map_eq_of_infinite (fun i : ℕ => fun k : Fin N => c (i + k))
    rcases Nat.lt_or_gt_of_ne hab with h | h
    · refine ⟨a, b - a, by omega, fun k hk => ?_⟩
      have := congr_fun he ⟨k, hk⟩
      simp only at this
      rw [this, show a + (b - a) + k = b + k by omega]
    · refine ⟨b, a - b, by omega, fun k hk => ?_⟩
      have := congr_fun he ⟨k, hk⟩
      simp only at this
      rw [← this, show b + (a - b) + k = a + k by omega]
  have hper := periodic_of_window c hP hw
  have : NeZero P := ⟨hP.ne'⟩
  -- the forbidden differences up to `N`, and the symmetric set of steps in `ZMod P`
  let T : Set ℕ := {t | 0 < t ∧ t ≤ N ∧ ∀ m, c m ≠ c (m + t)}
  let S : Set (ZMod P) := {x | ∃ t ∈ T, x = (t : ZMod P) ∨ x = -(t : ZMod P)}
  have hS : ∀ s ∈ S, -s ∈ S := by
    rintro s ⟨t, ht, rfl | rfl⟩
    · exact ⟨t, ht, Or.inr rfl⟩
    · exact ⟨t, ht, Or.inl (neg_neg _)⟩
  -- the periodic colouring of `ZMod P` (note `ZMod 3 = Fin 3` definitionally)
  let col : ZMod P → ZMod 3 := fun x => c (i + x.val)
  have hcol1 : ∀ g : ZMod P, ∀ t ∈ T, col (g + t) ≠ col g := by
    intro g t ht
    obtain ⟨-, htN, htc⟩ := ht
    have hv : (g + (t : ZMod P)).val = (g.val + t) % P := by
      rw [ZMod.val_add, ZMod.val_natCast, Nat.add_mod_mod]
    have hlt : g.val + t < P + N := by
      have := ZMod.val_lt g
      omega
    show c (i + (g + (t : ZMod P)).val) ≠ c (i + g.val)
    rw [hv, ← hper _ hlt, ← add_assoc]
    exact (htc _).symm
  have hcol : ∀ g, ∀ s ∈ S, col (g + s) ≠ col g := by
    rintro g s ⟨t, ht, rfl | rfl⟩
    · exact hcol1 g t ht
    · intro h
      apply hcol1 (g + -(t : ZMod P)) t ht
      rw [neg_add_cancel_right]
      exact h.symm
  obtain ⟨ξ, hξ⟩ := TheoremW.theoremW_finite S hS col hcol
  -- extend the character to all of `ZMod P` (`ℝ/ℤ` is divisible, hence injective)
  obtain ⟨F, hF⟩ := (Module.Baer.of_divisible UnitAddCircle).extension_property_addMonoidHom
    (AddSubgroup.closure S).subtype (AddSubgroup.subtype_injective _) ξ
  refine ⟨F 1, fun t ht0 htN htc => ?_⟩
  have hmem : (t : ZMod P) ∈ S := ⟨t, ⟨ht0, htN, htc⟩, Or.inl rfl⟩
  obtain ⟨x, hx1, hx2, hx⟩ := hξ _ hmem
  have key : t • F 1 = (x : UnitAddCircle) := by
    rw [hx, ← hF, ← map_nsmul, nsmul_eq_mul, mul_one]
    rfl
  rw [key, UnitAddCircle.norm_eq]
  exact third_le_abs_sub_round hx1 hx2

/-- Step 2 (compactness of `ℝ/ℤ`): one `θ` works for all forbidden differences. -/
lemma exists_theta (c : ℕ → Fin 3) :
    ∃ θ : UnitAddCircle, ∀ t : ℕ, 0 < t → (∀ m, c m ≠ c (m + t)) → 1/3 ≤ ‖t • θ‖ := by
  have : Fact ((0 : ℝ) < 1) := ⟨one_pos⟩
  let ι := {t : ℕ // 0 < t ∧ ∀ m, c m ≠ c (m + t)}
  let A : ι → Set UnitAddCircle := fun t => {θ | 1/3 ≤ ‖t.1 • θ‖}
  have hA : ∀ t, IsClosed (A t) := fun t =>
    isClosed_le continuous_const (continuous_norm.comp (continuous_nsmul t.1))
  obtain ⟨θ, -, hθ⟩ := isCompact_univ.inter_iInter_nonempty A hA fun u => by
    obtain ⟨θ, hθ⟩ := finite_step c (u.sup fun t => t.1)
    exact ⟨θ, Set.mem_univ _, Set.mem_iInter₂.2 fun t ht =>
      hθ t.1 t.2.1 (Finset.le_sup (f := fun t : ι => t.1) ht) t.2.2⟩
  exact ⟨θ, fun t ht0 htc => Set.mem_iInter.1 hθ ⟨t, ht0, htc⟩⟩

/-- **Question 3 of Glasscock–Koutsogiannis–Richter** (Bull. Amer. Math. Soc. 59 (2022),
569–606), answered for three colours: for every 3-colouring `c` of `ℕ` there is a real `α` such
that every `n ≥ 1` with `‖n α‖ < 1/3` — the distance from `n α` to the nearest integer, written
`|n α - round (n α)|` — is a difference of two numbers of the same colour: `c m = c (m + n)` for
some `m`. -/
theorem question3 (c : ℕ → Fin 3) :
    ∃ α : ℝ, ∀ n : ℕ, 0 < n → |(n : ℝ) * α - round ((n : ℝ) * α)| < 1/3 →
      ∃ m : ℕ, c m = c (m + n) := by
  obtain ⟨θ, hθ⟩ := exists_theta c
  obtain ⟨α, rfl⟩ := QuotientAddGroup.mk_surjective θ
  refine ⟨α, fun n hn hlt => ?_⟩
  by_contra h
  push Not at h
  have h1 := hθ n hn h
  rw [← AddCircle.coe_nsmul, nsmul_eq_mul, UnitAddCircle.norm_eq] at h1
  linarith

/-- **Chromatic recurrence** (equivalent form of `question3` for `ℤ`): if `S ⊆ ℤ \ {0}` contains,
for every `α : ℝ`, some `s` with `‖s α‖ < 1/3`, then for every 3-colouring `c` of `ℤ` some
element of `S` is a difference of two integers of the same colour. -/
theorem chromatic_recurrence (S : Set ℤ) (h0 : (0 : ℤ) ∉ S)
    (hS : ∀ α : ℝ, ∃ s ∈ S, |(s : ℝ) * α - round ((s : ℝ) * α)| < 1/3) (c : ℤ → Fin 3) :
    ∃ x : ℤ, ∃ s ∈ S, c x = c (x + s) := by
  by_contra hcon
  push Not at hcon
  obtain ⟨θ, hθ⟩ := exists_theta (fun m : ℕ => c m)
  obtain ⟨α, rfl⟩ := QuotientAddGroup.mk_surjective θ
  obtain ⟨s, hs, hlt⟩ := hS α
  have hs0 : s ≠ 0 := fun h => h0 (h ▸ hs)
  have h1 : ‖s • (α : UnitAddCircle)‖ < 1/3 := by
    rw [← AddCircle.coe_zsmul, zsmul_eq_mul, UnitAddCircle.norm_eq]
    exact hlt
  obtain ⟨n, rfl | rfl⟩ : ∃ n : ℕ, s = n ∨ s = -n := ⟨s.natAbs, Int.natAbs_eq s⟩
  · have hn : 0 < n := by omega
    have h2 := hθ n hn fun m (h : c (m : ℤ) = c ((m + n : ℕ) : ℤ)) =>
      hcon m n hs (by rw [h, Nat.cast_add])
    rw [natCast_zsmul] at h1
    linarith
  · have hn : 0 < n := by omega
    have h2 := hθ n hn fun m (h : c (m : ℤ) = c ((m + n : ℕ) : ℤ)) =>
      hcon (m + n : ℕ) (-n) hs (by
        rw [Nat.cast_add, add_neg_cancel_right, ← Nat.cast_add]
        exact h.symm)
    rw [neg_zsmul, natCast_zsmul, norm_neg] at h1
    linarith

/-- Sanity check: the hypotheses of `chromatic_recurrence` are satisfiable, e.g. by
`S = {1, 2, 3}` (Dirichlet's approximation theorem gives `‖k α‖ ≤ 1/4` for some `k ≤ 3`). -/
example : ∀ α : ℝ, ∃ s ∈ ({1, 2, 3} : Set ℤ), |(s : ℝ) * α - round ((s : ℝ) * α)| < 1/3 := by
  intro α
  obtain ⟨k, hk0, hk3, hk⟩ := Real.exists_nat_abs_mul_sub_round_le α (n := 3) (by norm_num)
  refine ⟨k, ?_, ?_⟩
  · simp only [Set.mem_insert_iff, Set.mem_singleton_iff]
    omega
  · push_cast
    exact lt_of_le_of_lt hk (by norm_num)

end GKR
