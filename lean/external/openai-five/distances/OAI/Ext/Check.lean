import OAI.Ext.WeakColoringD

/-!
# Axiom check for the distance-set generalisation of OpenAI's Borel transfer

Restatements (written here) of the results of `OAI/Ext/BorelTransferD.lean` and
`OAI/Ext/WeakColoringD.lean`, followed by `#print axioms`. `angularMeasure` (normalised Haar measure
on `AddCircle (2π)`) and `direction` (`t ↦ exp (i t)`) are OpenAI's definitions from
`OAI/Geometry/PlaneColoring/Coloring.lean`; `DensityOne` is OpenAI's definition from
`OAI/Geometry/PlaneColoring/WeakColoring.lean`.
-/

open MeasureTheory OAI.EuclideanFiveColor

/-- Main theorem: a colouring of the plane avoiding a set `D` of positive real algebraic distances
yields one measurable colouring which, for each `r ∈ D`, separates almost every pair
`(x, x + r e^{it})` (almost every `x`, almost every angle `t`). -/
theorem borel_transfer_distances {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure,
      d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2) :=
  OAI.Ext.proper_to_borel_weak_distances hpos halg hc

/-- Density version: moreover no two density-one points of one colour class of `d` are at a
distance belonging to `D`. -/
theorem borel_transfer_distances_density {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧
      (∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2)) ∧
      (∀ i : Fin k, ∀ x y : ℂ, ‖x - y‖ ∈ D → DensityOne {z | d z = i} x → ¬ DensityOne {z | d z = i} y) :=
  OAI.Ext.proper_to_density_avoiding hpos halg hc

/-- Measurable classes: a measurable colouring which avoids every distance of `D` on a measurable
conull set. -/
theorem measurable_proper_on_conull {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∃ G : Set ℂ, MeasurableSet G ∧ volume Gᶜ = 0 ∧
      ∀ x ∈ G, ∀ y ∈ G, ‖x - y‖ ∈ D → d x ≠ d y :=
  OAI.Ext.proper_to_measurable_proper_on_conull hpos halg hc

/-- Every displacement vector whose length lies in `D` is separated almost everywhere. -/
theorem borel_transfer_all_vectors {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∀ v : ℂ, ‖v‖ ∈ D → ∀ᵐ x ∂(volume : Measure ℂ), d x ≠ d (x + v) :=
  OAI.Ext.proper_to_borel_all_vectors hpos halg hc

/-- `k` pairwise disjoint measurable sets avoiding every distance of `D`, covering almost all of
the plane. -/
theorem measurable_classes_avoiding {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ B : Fin k → Set ℂ, (∀ i, MeasurableSet (B i)) ∧ Pairwise (fun i j => Disjoint (B i) (B j)) ∧
      (∀ i, ∀ x ∈ B i, ∀ y ∈ B i, ‖x - y‖ ∉ D) ∧ volume (⋃ i, B i)ᶜ = 0 :=
  OAI.Ext.proper_to_measurable_classes hpos halg hc

/-- The case `D = {1}`: OpenAI's `proper_to_borel_weak` statement. -/
theorem borel_transfer_unit {k : ℕ} [NeZero k] (c : ℂ → Fin k) (hc : ProperColoring k c) :
    ∃ d : ℂ → Fin k, Measurable d ∧ WeakColoring k d :=
  OAI.Ext.proper_to_borel_weak_of_singleton hc

#print axioms borel_transfer_distances
#print axioms borel_transfer_distances_density
#print axioms measurable_proper_on_conull
#print axioms borel_transfer_all_vectors
#print axioms measurable_classes_avoiding
#print axioms borel_transfer_unit
#print axioms OAI.Ext.proper_to_borel_weak_distances
#print axioms OAI.Ext.proper_to_density_avoiding
#print axioms OAI.Ext.exists_measurable_classes_avoiding
#print axioms OAI.Ext.proper_to_measurable_proper_on_conull
#print axioms OAI.Ext.proper_to_borel_all_vectors
#print axioms OAI.Ext.proper_to_measurable_classes
#print axioms OAI.Ext.proper_to_borel_weak_of_singleton
#print axioms OAI.Ext.countable_of_algebraic
