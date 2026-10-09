import OAI.Geometry.PlaneColoring.Five

/-! Independent restatement (written here, not taken from the challenge file): no colouring of the plane `ℂ` with
five colours gives different colours to every two points at distance one. -/

theorem plane_not_five_colourable :
    ¬ ∃ c : ℂ → Fin 5, ∀ p q : ℂ, ‖p - q‖ = 1 → c p ≠ c q := by
  simpa [OAI.EuclideanFiveColor.ProperColoring] using OAI.EuclideanFiveColor.no_proper_five_coloring

/-- the same with any natural number `k ≤ 5` of colours -/
theorem plane_not_k_colourable (k : ℕ) (hk : k ≤ 5) :
    ¬ ∃ c : ℂ → Fin k, ∀ p q : ℂ, ‖p - q‖ = 1 → c p ≠ c q := by
  rintro ⟨c, hc⟩
  exact plane_not_five_colourable ⟨fun z => Fin.castLE hk (c z), fun p q h e => hc p q h (Fin.castLE_injective hk e)⟩

#print axioms plane_not_five_colourable
#print axioms plane_not_k_colourable
