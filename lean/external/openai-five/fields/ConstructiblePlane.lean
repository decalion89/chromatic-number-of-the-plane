import OAI.Geometry.PlaneColoring.Five

/-! Here `E` is the field of constructible numbers (the smallest subfield of `ℂ` closed under square roots: the
points constructible with ruler and compass). OpenAI's argument, with the field changed and its triple average
taken at the exponents `-1, 0, 1`, shows that the plane over that field is not 5-colourable. -/

namespace OAI
noncomputable section
open MeasureTheory
namespace EuclideanFiveColor

theorem constructible_plane_not_five_colourable :
    ¬ ∃ c : PlaneFiveColor.Spectral.E → Fin 5,
      ∀ x y : PlaneFiveColor.Spectral.E, ‖(x : ℂ) - (y : ℂ)‖ = 1 → c x ≠ c y := by
  rintro ⟨c, hc⟩
  let ω : PlaneFiveColor.LabelModel.Space 5 := ⟨c, hc⟩
  obtain ⟨H, hm, hp, ho⟩ := PlaneFiveColor.LabelModel.exists_label_field 5 ω
  obtain ⟨d, hd, hunit⟩ := sample_field (PlaneFiveColor.LabelModel.measure 5 ω) angularMeasure
    direction direction_continuous.measurable
    H hm hp (fun x t => ho x (AddCircle.toCircle t))
  exact no_measurable_weak_five_coloring hd ⟨hd.aemeasurable, hunit⟩

/-- The same statement with the field written out: the smallest subfield of `ℂ` closed under square roots. -/
theorem sqrtClosure_plane_not_five_colourable :
    ¬ ∃ c : (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ) → Fin 5,
      ∀ x y : (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ),
        ‖(x : ℂ) - (y : ℂ)‖ = 1 → c x ≠ c y :=
  constructible_plane_not_five_colourable

end EuclideanFiveColor
end
end OAI

#print axioms OAI.EuclideanFiveColor.constructible_plane_not_five_colourable
#print axioms OAI.EuclideanFiveColor.sqrtClosure_plane_not_five_colourable
