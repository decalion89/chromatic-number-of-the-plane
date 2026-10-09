import OAI.Geometry.PlaneColoring.Five
import PlaneSeven

/-!
# `χ(ℝ²) ∈ {6, 7}` in the formulation of formal-conjectures

The lower bound is OpenAI's theorem `OAI.EuclideanFiveColor.no_proper_five_coloring`
(openai/math, commit fd4aeeb2), stated for colourings of `ℂ`; it is moved to
`EuclideanSpace ℝ (Fin 2)` along the isometry `Complex.orthonormalBasisOneI.repr`.
The upper bound is `PlaneSeven.hadwigerNelson_atMostSeven`.
-/

namespace PlaneSeven

open SimpleGraph

/-- `6 ≤ χ(ℝ²)`, from OpenAI's theorem. -/
theorem hadwigerNelson_atLeastSix :
    6 ≤ SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) := by
  rw [show (6 : ℕ∞) = ((6 : ℕ) : ℕ∞) from rfl, le_chromaticNumber_iff_colorable]
  intro m hm
  by_contra h
  have hm6 : m < 6 := by
    rcases lt_or_ge m 6 with h' | h'
    · exact h'
    · exact absurd (by exact_mod_cast h') h
  obtain ⟨C⟩ := hm.mono (by omega : m ≤ 5)
  apply OAI.EuclideanFiveColor.no_proper_five_coloring
  refine ⟨fun z => C ⟨Complex.orthonormalBasisOneI.repr z, Set.mem_univ _⟩, fun x y hxy => C.valid ?_⟩
  change dist (Complex.orthonormalBasisOneI.repr x) (Complex.orthonormalBasisOneI.repr y) = 1
  rw [LinearIsometryEquiv.dist_map, dist_eq_norm, hxy]

/-- `Erdos508.HadwigerNelsonAtLeastFive` of formal-conjectures. -/
theorem hadwigerNelson_atLeastFive :
    5 ≤ SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) :=
  le_trans (by norm_num) hadwigerNelson_atLeastSix

/-- `Erdos508.HadwigerNelsonAtLeast4` of formal-conjectures. -/
theorem hadwigerNelson_atLeast4 :
    4 ≤ SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) :=
  le_trans (by norm_num) hadwigerNelson_atLeastSix

/-- The chromatic number of the plane is 6 or 7. -/
theorem hadwigerNelson_six_or_seven :
    SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) = 6 ∨
      SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) = 7 := by
  have h6 := hadwigerNelson_atLeastSix
  have h7 := hadwigerNelson_atMostSeven
  obtain ⟨n, hn⟩ : ∃ n : ℕ, SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) = n :=
    ENat.ne_top_iff_exists.mp (ne_top_of_le_ne_top (by simp) h7) |>.imp fun _ h => h.symm
  rw [hn] at h6 h7 ⊢
  have a : 6 ≤ n := by exact_mod_cast h6
  have b : n ≤ 7 := by exact_mod_cast h7
  interval_cases n <;> simp

end PlaneSeven

#print axioms PlaneSeven.hadwigerNelson_atLeastSix
#print axioms PlaneSeven.hadwigerNelson_six_or_seven
