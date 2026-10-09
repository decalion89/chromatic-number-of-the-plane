import OAI.Geometry.PlaneColoring.Five

/-! Here `E` is the field of numbers expressible by radicals (`solvableByRad ℚ ℂ`): OpenAI's argument, with
the field changed, shows that the plane over that field is not 5-colourable. -/

namespace OAI
noncomputable section
open MeasureTheory
namespace EuclideanFiveColor

theorem radical_plane_not_five_colourable :
    ¬ ∃ c : PlaneFiveColor.Spectral.E → Fin 5,
      ∀ x y : PlaneFiveColor.Spectral.E, ‖(x : ℂ) - (y : ℂ)‖ = 1 → c x ≠ c y := by
  rintro ⟨c, hc⟩
  let ω : PlaneFiveColor.LabelModel.Space 5 := ⟨c, hc⟩
  obtain ⟨H, hm, hp, ho⟩ := PlaneFiveColor.LabelModel.exists_label_field 5 ω
  obtain ⟨d, hd, hunit⟩ := sample_field (PlaneFiveColor.LabelModel.measure 5 ω) angularMeasure
    direction direction_continuous.measurable
    H hm hp (fun x t => ho x (AddCircle.toCircle t))
  exact no_measurable_weak_five_coloring hd ⟨hd.aemeasurable, hunit⟩

/-- The same statement with the field written out: the points of the plane whose coordinates are expressible by
radicals cannot be coloured with 5 colours so that points at distance 1 get different colours. -/
theorem solvableByRad_plane_not_five_colourable :
    ¬ ∃ c : solvableByRad ℚ ℂ → Fin 5,
      ∀ x y : solvableByRad ℚ ℂ, ‖(x : ℂ) - (y : ℂ)‖ = 1 → c x ≠ c y :=
  radical_plane_not_five_colourable

end EuclideanFiveColor
end
end OAI

#print axioms OAI.EuclideanFiveColor.radical_plane_not_five_colourable
#print axioms OAI.EuclideanFiveColor.solvableByRad_plane_not_five_colourable
