import OAI.Ext.BorelTransferD
import OAI.Geometry.PlaneColoring.WeakColoring

/-!
# Density version, for a countable set of algebraic forbidden distances

OpenAI's file `OAI/Geometry/PlaneColoring/WeakColoring.lean` shows (`weak_density_proper`) that for
a measurable weak colouring (distance `1`), two points at distance `1` cannot both be density-one
points of the same colour class.

Provenance of the parts of this file.
* Copied from `WeakColoring.lean` and generalised from distance `1` to an arbitrary radius `r`:
  `weak_all_directions_radius` (from `weak_all_directions`) and `weak_density_proper_radius`
  (from `weak_density_proper`). OpenAI's generic lemmas `all_parameters_product_zero`,
  `truncatedVector*`, `densityOne_disjoint`, `densityOne_translate`, `direction_exists`,
  `ae_density_actual` are used unchanged.
* New: `exists_direction_radius`, `all_vectors_of_weak`, `density_classes_avoid`,
  `exists_measurable_classes_avoiding`, `proper_to_density_avoiding`, `proper_to_borel_all_vectors`,
  `proper_to_measurable_classes`, `proper_to_measurable_proper_on_conull` (combining the above with
  the main theorem `proper_to_borel_weak_distances` of `OAI/Ext/BorelTransferD.lean`), and the
  sanity check `proper_to_borel_weak_of_singleton` (the case `D = {1}` gives back OpenAI's
  `proper_to_borel_weak` statement).
-/

namespace OAI

noncomputable section
open MeasureTheory Filter Topology
namespace Ext
open EuclideanFiveColor
/-- every vector of length `r > 0` is `r * direction t` for some angle `t` (new) -/
lemma exists_direction_radius {r : ℝ} (hr : 0 < r) {v : ℂ} (hv : ‖v‖ = r) :
    ∃ t : Angles, (r:ℂ) * direction t = v := by
  have hr' : (r:ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hr.ne'
  have hu : ‖v / (r:ℂ)‖ = 1 := by
    rw [norm_div,Complex.norm_real,Real.norm_eq_abs,abs_of_pos hr,hv,div_self hr.ne']
  obtain ⟨t,ht⟩ := direction_exists hu
  exact ⟨t,by rw [ht]; exact mul_div_cancel₀ v hr'⟩
variable {k : ℕ} (c : ℂ → Fin k) (hc : Measurable c)
include hc
/-- Generalises OpenAI's `weak_all_directions` from distance `1` to distance `r`. -/
lemma weak_all_directions_radius (r : ℝ)
    (hweak : ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, c p.1 ≠ c (p.1 + (r:ℂ) * direction p.2)) :
    ∀ t : Angles, ∀ᵐ z ∂(volume : Measure ℂ), c z ≠ c (z+(r:ℂ) * direction t) := by
  have hdc : Continuous (fun t : Angles => (r:ℂ) * direction t) := continuous_const.mul direction_continuous
  have ho : ∀ᵐ t ∂angularMeasure, ∀ᵐ z ∂(volume : Measure ℂ), c z ≠ c (z+(r:ℂ) * direction t) := by
    apply (Measure.ae_ae_comm ?_).mp (Measure.ae_ae_of_ae_prod hweak)
    exact (measurableSet_eq_fun (hc.comp measurable_fst)
      (hc.comp (measurable_fst.add (hdc.measurable.comp measurable_snd)))).compl
  have hz (i : Fin k) (n : ℕ) : ∀ t : Angles, ∀ᵐ z ∂(volume : Measure ℂ),
      truncatedVector c hc i n z*truncatedVector c hc i n (z+(r:ℂ) * direction t)=0 := by
    apply all_parameters_product_zero angularMeasure (fun t => (r:ℂ) * direction t) hdc
    · filter_upwards [truncatedVector_ae c hc i n] with z hz
      rw [hz]
      exact Set.indicator_nonneg (fun _ _ => zero_le_one) _
    · filter_upwards [ho] with t ht
      filter_upwards [ht,truncatedVector_ae c hc i n,
        (measurePreserving_add_right volume ((r:ℂ) * direction t)).quasiMeasurePreserving.ae
          (truncatedVector_ae c hc i n)] with z hne h1 h2
      rw [h1,h2]
      by_cases h : z ∈ truncatedClass c i n
      · have hh : z+(r:ℂ) * direction t ∉ truncatedClass c i n := fun hh => hne (h.1.trans hh.1.symm)
        simp [Set.indicator_of_notMem hh]
      · simp [Set.indicator_of_notMem h]
  intro t
  have hi : ∀ᵐ z ∂(volume : Measure ℂ), ∀ i : Fin k, ∀ n : ℕ,
      ((truncatedClass c i n).indicator (fun _ => (1 : ℝ)) z)*
      ((truncatedClass c i n).indicator (fun _ => (1 : ℝ)) (z+(r:ℂ) * direction t))=0 := by
    apply ae_all_iff.mpr
    intro i
    apply ae_all_iff.mpr
    intro n
    filter_upwards [hz i n t,truncatedVector_ae c hc i n,
      (measurePreserving_add_right volume ((r:ℂ) * direction t)).quasiMeasurePreserving.ae
        (truncatedVector_ae c hc i n)] with z h0 h1 h2
    simpa only [h1,h2] using h0
  filter_upwards [hi] with z hz
  intro he
  obtain ⟨n,hn⟩ := exists_nat_gt (max ‖z‖ ‖z+(r:ℂ) * direction t‖)
  have h1 : z ∈ truncatedClass c (c z) n := ⟨rfl,by
    simpa only [Metric.mem_closedBall,dist_zero_right] using (lt_of_le_of_lt (le_max_left _ _) hn).le⟩
  have h2 : z+(r:ℂ) * direction t ∈ truncatedClass c (c z) n := ⟨he.symm,by
    simpa only [Metric.mem_closedBall,dist_zero_right] using (lt_of_le_of_lt (le_max_right _ _) hn).le⟩
  have h0 := hz (c z) n
  simp [Set.indicator_of_mem h1,Set.indicator_of_mem h2] at h0
/-- Generalises OpenAI's `weak_density_proper` from distance `1` to a distance `r > 0`. -/
lemma weak_density_proper_radius {r : ℝ} (hr : 0 < r)
    (hw : ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, c p.1 ≠ c (p.1 + (r:ℂ) * direction p.2))
    (i : Fin k) {x y : ℂ} (hxy : ‖x-y‖=r) (hx : DensityOne {z | c z=i} x) :
    ¬ DensityOne {z | c z=i} y := by
  let u := y-x
  have hr' : (r:ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hr.ne'
  have hu : ‖u / (r:ℂ)‖ = 1 := by
    have hn : ‖u‖ = r := by simpa only [u,norm_sub_rev] using hxy
    rw [norm_div,Complex.norm_real,Real.norm_eq_abs,abs_of_pos hr,hn,div_self hr.ne']
  obtain ⟨t,ht⟩ := direction_exists hu
  have hut : (r:ℂ) * direction t = u := by rw [ht]; exact mul_div_cancel₀ u hr'
  have hd : AEDisjoint volume {z | c z=i} ((fun z : ℂ => z+u) ⁻¹' {z | c z=i}) := by
    apply (measure_eq_zero_iff_ae_notMem).mpr
    filter_upwards [weak_all_directions_radius c hc r hw t] with z hz
    rintro ⟨h1,h2⟩
    exact hz (by simpa only [hut] using h1.trans h2.symm)
  have hh := densityOne_disjoint ((measurableSet_eq_fun hc measurable_const).preimage
    (measurable_id.add_const u)) hd hx
  simp only [id_eq] at hh
  rw [densityOne_translate] at hh
  simpa only [u,add_sub_cancel] using hh
omit hc
/-- For a measurable colouring `d` satisfying the conclusion of the main theorem, every
displacement `v` with `‖v‖ ∈ D` (not only almost every direction) is separated almost everywhere. -/
theorem all_vectors_of_weak {D : Set ℝ} (hpos : ∀ r ∈ D, 0 < r)
    (d : ℂ → Fin k) (hd : Measurable d)
    (hw : ∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2))
    (v : ℂ) (hv : ‖v‖ ∈ D) : ∀ᵐ x ∂(volume : Measure ℂ), d x ≠ d (x + v) := by
  obtain ⟨t,ht⟩ := exists_direction_radius (hpos _ hv) (rfl : ‖v‖ = ‖v‖)
  simpa only [ht] using weak_all_directions_radius d hd ‖v‖ (hw _ hv) t
/-- Density version: for a measurable colouring `d` satisfying the conclusion of the main theorem,
the density-one points of a colour class contain no two points at a distance belonging to `D`. -/
theorem density_classes_avoid {D : Set ℝ} (hpos : ∀ r ∈ D, 0 < r)
    (d : ℂ → Fin k) (hd : Measurable d)
    (hw : ∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2))
    (i : Fin k) {x y : ℂ} (hxy : ‖x - y‖ ∈ D) (hx : DensityOne {z | d z = i} x) :
    ¬ DensityOne {z | d z = i} y :=
  weak_density_proper_radius d hd (hpos _ hxy) (hw _ hxy) i rfl hx
/-- The colour classes of such a `d`, modified on null sets, are measurable sets avoiding every
distance of `D`, and their union is conull. -/
theorem exists_measurable_classes_avoiding {D : Set ℝ} (hpos : ∀ r ∈ D, 0 < r)
    (d : ℂ → Fin k) (hd : Measurable d)
    (hw : ∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2)) :
    ∃ B : Fin k → Set ℂ, (∀ i, MeasurableSet (B i)) ∧ (∀ i, B i ⊆ {z | d z = i}) ∧
      (∀ i, B i =ᵐ[volume] {z | d z = i}) ∧
      (∀ i, ∀ x ∈ B i, ∀ y ∈ B i, ‖x - y‖ ∉ D) ∧ volume (⋃ i, B i)ᶜ = 0 := by
  -- the points that are not density-one points of their own colour class form a null set `⊆ M`
  obtain ⟨M,hNM,hMm,hM⟩ := exists_measurable_superset_of_null (ae_iff.mp (ae_density_actual d hd))
  have key (i : Fin k) (x : ℂ) (hx : x ∈ {z | d z = i} \ M) : DensityOne {z | d z = i} x := by
    by_contra h
    have hdx : d x = i := hx.1
    refine hx.2 (hNM ?_)
    change ¬ DensityOne {z | d z = d x} x
    rw [hdx]
    exact h
  refine ⟨fun i => {z | d z = i} \ M,fun i => (measurableSet_eq_fun hd measurable_const).diff hMm,
    fun i => Set.sdiff_subset,fun i => sdiff_null_ae_eq_self hM,?_,?_⟩
  · intro i x hx y hy hxy
    exact density_classes_avoid hpos d hd hw i hxy (key i x hx) (key i y hy)
  · apply measure_mono_null _ hM
    intro x hx
    by_contra hxM
    exact hx (Set.mem_iUnion.mpr ⟨d x,Set.mem_sdiff_of_mem rfl hxM⟩)
/-- Density version of the main theorem: from a colouring avoiding all distances of `D`, a
measurable colouring `d` satisfying the conclusion of the main theorem, whose colour classes have no
two density-one points at a distance of `D`. -/
theorem proper_to_density_avoiding [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    {c : ℂ → Fin k} (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧
      (∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure, d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2)) ∧
      (∀ i : Fin k, ∀ x y : ℂ, ‖x - y‖ ∈ D → DensityOne {z | d z = i} x → ¬ DensityOne {z | d z = i} y) := by
  obtain ⟨d,hd,hw⟩ := proper_to_borel_weak_distances hpos halg hc
  exact ⟨d,hd,hw,fun i x y hxy hx => density_classes_avoid hpos d hd hw i hxy hx⟩
/-- From a colouring avoiding all distances of `D`, a measurable colouring `d` such that for every
vector `v` with `‖v‖ ∈ D`, `d x ≠ d (x + v)` for almost every `x`. -/
theorem proper_to_borel_all_vectors [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    {c : ℂ → Fin k} (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∀ v : ℂ, ‖v‖ ∈ D → ∀ᵐ x ∂(volume : Measure ℂ), d x ≠ d (x + v) := by
  obtain ⟨d,hd,hw⟩ := proper_to_borel_weak_distances hpos halg hc
  exact ⟨d,hd,fun v hv => all_vectors_of_weak hpos d hd hw v hv⟩
/-- From a colouring avoiding all distances of `D`: `k` pairwise disjoint measurable sets, none of
which contains two points at a distance belonging to `D`, whose union is conull. -/
theorem proper_to_measurable_classes [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    {c : ℂ → Fin k} (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ B : Fin k → Set ℂ, (∀ i, MeasurableSet (B i)) ∧ Pairwise (fun i j => Disjoint (B i) (B j)) ∧
      (∀ i, ∀ x ∈ B i, ∀ y ∈ B i, ‖x - y‖ ∉ D) ∧ volume (⋃ i, B i)ᶜ = 0 := by
  obtain ⟨d,hd,hw⟩ := proper_to_borel_weak_distances hpos halg hc
  obtain ⟨B,hBm,hBsub,-,hBavoid,hBnull⟩ := exists_measurable_classes_avoiding hpos d hd hw
  refine ⟨B,hBm,fun i j hij => Set.disjoint_left.mpr fun x hi hj => ?_,hBavoid,hBnull⟩
  have h1 : d x = i := hBsub i hi
  have h2 : d x = j := hBsub j hj
  exact hij (h1.symm.trans h2)
/-- From a colouring avoiding all distances of `D`, a measurable colouring which avoids all
distances of `D` on a measurable conull set. -/
theorem proper_to_measurable_proper_on_conull [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    {c : ℂ → Fin k} (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∃ G : Set ℂ, MeasurableSet G ∧ volume Gᶜ = 0 ∧
      ∀ x ∈ G, ∀ y ∈ G, ‖x - y‖ ∈ D → d x ≠ d y := by
  obtain ⟨d,hd,hw⟩ := proper_to_borel_weak_distances hpos halg hc
  obtain ⟨B,hBm,hBsub,-,hBavoid,hBnull⟩ := exists_measurable_classes_avoiding hpos d hd hw
  refine ⟨d,hd,⋃ i, B i,MeasurableSet.iUnion hBm,hBnull,?_⟩
  intro x hx y hy hxy he
  obtain ⟨i,hi⟩ := Set.mem_iUnion.mp hx
  obtain ⟨j,hj⟩ := Set.mem_iUnion.mp hy
  have hxi : d x = i := hBsub i hi
  have hyj : d y = j := hBsub j hj
  have hij : i = j := hxi.symm.trans (he.trans hyj)
  subst hij
  exact hBavoid i x hi y hj hxy
/-- Sanity check: the case `D = {1}` of the main theorem gives back OpenAI's
`EuclideanFiveColor.proper_to_borel_weak` (same hypothesis, same conclusion). -/
theorem proper_to_borel_weak_of_singleton [NeZero k] {c : ℂ → Fin k} (hc : ProperColoring k c) :
    ∃ d : ℂ → Fin k, Measurable d ∧ WeakColoring k d := by
  obtain ⟨d,hd,hw⟩ := proper_to_borel_weak_distances (D := {1})
    (fun r hr => by rw [Set.mem_singleton_iff.mp hr]; exact one_pos)
    (fun r hr => by rw [Set.mem_singleton_iff.mp hr]; exact isAlgebraic_one)
    (fun x y h => hc x y (Set.mem_singleton_iff.mp h))
  refine ⟨d,hd,hd.aemeasurable,?_⟩
  simpa only [Complex.ofReal_one,one_mul] using hw 1 (Set.mem_singleton 1)
end Ext
end

end OAI
