import OAI.Ext.LabelModelD
import OAI.Geometry.PlaneColoring.BorelTransfer

/-!
# Borel transfer for a countable set of algebraic forbidden distances

OpenAI's file `OAI/Geometry/PlaneColoring/BorelTransfer.lean` proves `proper_to_borel_weak`:
a colouring of the plane in which no two points at distance `1` share a colour can be replaced by a
measurable colouring `d` such that, for almost every point `x` and almost every direction `t`,
`d x ≠ d (x + direction t)`.

This file generalises that statement to an arbitrary set `D` of positive real algebraic numbers
(such a set is countable, `countable_of_algebraic` below): a colouring in which no two points at a
distance belonging to `D` share a colour can be replaced by ONE measurable colouring `d` such that,
for every `r ∈ D`, for almost every point `x` and almost every direction `t`,
`d x ≠ d (x + r * direction t)` (`Ext.proper_to_borel_weak_distances`).

Provenance of the parts of this file.
* Copied from `BorelTransfer.lean` (namespace `PlaneFiveColor.LabelModel`) and adapted to the label
  model `Ext.LabelModelD.Space D k` of `OAI/Ext/LabelModelD.lean` (extra argument `D`, otherwise
  unchanged): `colorValue`, `color_memLp`, `colorVector`, `color_ae`, `color_rotation_fixed`,
  `continuous_complete`, `projected`, `projected_range`, `projected_rotation_fixed`, `oneVector`,
  `one_ae`, `one_fixed`, `one_continuous`, `color_sum`, `projected_sum`, `projected_continuous`,
  `fieldVector`, `fieldVector_continuous`, `fieldVector_algebraic`, `fieldVector_sum`,
  `fieldVector_range`.
* Copied from `BorelTransfer.lean` and generalised from distance `1` to the distances in `D`:
  `color_correlation` (from `color_unit_correlation`: the correlation vanishes at every algebraic
  `z` with `‖z‖ ∈ D`, not only at `z ∈ K`), `projected_correlation` (from
  `projected_unit_correlation`; OpenAI's `residual_fourier` already holds for every nonzero algebraic
  `z`), `fieldVector_correlation` (from `fieldVector_unit_correlation`: radius `r ∈ D`; the density
  argument uses OpenAI's `circle_closed_induction` applied to the circle scaled by `r`, since
  `r * b` is algebraic for `b ∈ K`), `exists_label_field` (one field `H`, orthogonality at every
  radius of `D`), `sample_field_countable` (from `EuclideanFiveColor.sample_field`: a countable family
  of displacement maps sampled at a single label configuration), and the main theorem (from
  `proper_to_borel_weak`).
* New: `isAlgebraic_ofReal`, `realE`, `realE_coe`, `countable_of_algebraic`.
All other ingredients (`residual_fourier`, `projection_coefficient`, `projection_range`,
`projection_covariance_fixed`, `orbitExtension*`, `circle_closed_induction`, `dense_E`,
`exists_measurable_field`, `positive_orthogonal`, `positiveSelect*`, ...) are OpenAI's, used unchanged.
-/

namespace OAI

noncomputable section
open MeasureTheory Filter Topology
open scoped ENNReal NNReal InnerProductSpace
namespace Ext.LabelModelD
open PlaneFiveColor.Spectral EuclideanFiveColor EuclideanFiveColor.ContinuousVectors
variable (D : Set ℝ) (k : ℕ) (ω : Space D k)
def colorValue (i : Fin k) (η : Space D k) : ℂ := if η.val 0 = i then 1 else 0
lemma color_memLp (i : Fin k) : MemLp (colorValue D k i) 2 (measure D k ω) := by
  refine MemLp.of_bound ?_ 1 ?_
  · apply Measurable.aestronglyMeasurable
    apply Measurable.ite _ measurable_const measurable_const
    exact (measurableSet_singleton i).preimage (show Measurable (fun η : Space D k => η.val 0) from
      ((continuous_apply (0 : E)).comp continuous_subtype_val).measurable)
  · exact .of_forall fun η => by dsimp [colorValue]; split_ifs <;> simp
abbrev colorVector (i : Fin k) : Hilbert D k ω := (color_memLp D k ω i).toLp (colorValue D k i)
lemma color_ae (i : Fin k) : colorVector D k ω i =ᵐ[measure D k ω] colorValue D k i :=
  (color_memLp D k ω i).coeFn_toLp
lemma color_rotation_fixed (i : Fin k) (b : K) :
    (rotation D k ω b : Hilbert D k ω →L[ℂ] Hilbert D k ω) (colorVector D k ω i) = colorVector D k ω i := by
  rw [rotation_apply]
  apply Lp.ext
  filter_upwards [Lp.coeFn_compMeasurePreserving (colorVector D k ω i) (rotate_preserving D k ω b),
    (rotate_preserving D k ω b).quasiMeasurePreserving.ae (color_ae D k ω i),color_ae D k ω i] with η h1 h2 h3
  simp only [Function.comp_apply,h1,h2,h3,colorValue,rotate_apply,mul_zero]
/-- Generalises OpenAI's `color_unit_correlation`: the colour indicator is orthogonal to its
translate by any algebraic `z` whose length is a forbidden distance. -/
lemma color_correlation (i : Fin k) (z : E) (hz : ‖(z:ℂ)‖ ∈ D) :
    ⟪colorVector D k ω i,op (representation D k ω) z (colorVector D k ω i)⟫_ℂ = 0 := by
  rw [op_apply,L2.inner_def]
  apply integral_eq_zero_of_ae
  filter_upwards [color_ae D k ω i,Lp.coeFn_compMeasurePreserving (colorVector D k ω i) (trans_preserving D k ω z),
    (trans_preserving D k ω z).quasiMeasurePreserving.ae (color_ae D k ω i)] with η h1 h2 h3
  simp only [h1,h2,Function.comp_apply,h3,colorValue,trans_apply,zero_add]
  have hp : η.val 0 ≠ η.val z := η.property 0 z (by
    change ‖(0 : ℂ)-(z:ℂ)‖ ∈ D
    rwa [zero_sub,norm_neg])
  by_cases hi : η.val 0 = i
  · have hj : η.val z ≠ i := fun h => hp (hi.trans h.symm)
    simp [hi,hj]
  · simp [hi]
instance continuous_complete : CompleteSpace (space (A := E) (representation D k ω)) :=
  (isClosed_space (A := E) (representation D k ω)).completeSpace_coe
def projected (i : Fin k) : Hilbert D k ω :=
  (space (A := E) (representation D k ω)).starProjection (colorVector D k ω i)
lemma projected_range (i : Fin k) : ∀ᵐ η ∂measure D k ω,
    ((projected D k ω i) η).im = 0 ∧ 0 ≤ ((projected D k ω i) η).re ∧ ((projected D k ω i) η).re ≤ 1 := by
  apply projection_range (measure D k ω) (representation D k ω) (trans D k) (trans_preserving D k ω) (op_apply D k ω)
  filter_upwards [color_ae D k ω i] with η h
  rw [h]
  dsimp [colorValue]
  split_ifs <;> norm_num
lemma projected_rotation_fixed (i : Fin k) (b : K) :
    (rotation D k ω b : Hilbert D k ω →L[ℂ] Hilbert D k ω) (projected D k ω i) = projected D k ω i := by
  exact projection_covariance_fixed (representation D k ω) (rotation D k ω b)
    (Homeomorph.mulLeft₀ (kCoeff b) (kCoeff_ne_zero b)) (mul_zero _)
    (rotation_covariance D k ω b) _ (color_rotation_fixed D k ω i b)
/-- Generalises OpenAI's `projected_unit_correlation` from `z ∈ K` to every nonzero algebraic `z`
with `‖z‖ ∈ D` (OpenAI's `residual_fourier` holds for every nonzero algebraic `z`). -/
lemma projected_correlation (i : Fin k) (z : E) (hz0 : z ≠ 0) (hz : ‖(z:ℂ)‖ ∈ D) :
    ⟪projected D k ω i,op (representation D k ω) z (projected D k ω i)⟫_ℂ = 0 := by
  have hb := residual_fourier (representation D k ω) (rotation D k ω) (rotation_covariance D k ω)
    (colorVector D k ω i-projected D k ω i)
    ((space (A := E) (representation D k ω)).sub_starProjection_mem_orthogonal (colorVector D k ω i))
    (fun b => by rw [map_sub,color_rotation_fixed,projected_rotation_fixed]) z hz0
  exact (projection_coefficient (representation D k ω) (colorVector D k ω i) z hb).trans
    (color_correlation D k ω i z hz)
end Ext.LabelModelD
end

noncomputable section
open MeasureTheory Filter Topology
namespace Ext.LabelModelD
open PlaneFiveColor.Spectral EuclideanFiveColor EuclideanFiveColor.ContinuousVectors
variable (D : Set ℝ) (k : ℕ) (ω : Space D k)
def oneVector : Hilbert D k ω := (memLp_const (1 : ℂ)).toLp (fun _ : Space D k => (1 : ℂ))
lemma one_ae : oneVector D k ω =ᵐ[measure D k ω] fun _ => (1 : ℂ) := MemLp.coeFn_toLp _
lemma one_fixed (a : E) : op (representation D k ω) a (oneVector D k ω) = oneVector D k ω := by
  rw [op_apply]
  apply Lp.ext
  filter_upwards [Lp.coeFn_compMeasurePreserving (oneVector D k ω) (trans_preserving D k ω a),
    (trans_preserving D k ω a).quasiMeasurePreserving.ae (one_ae D k ω),one_ae D k ω] with η h1 h2 h3
  exact h1.trans (h2.trans h3.symm)
lemma one_continuous : oneVector D k ω ∈ space (A := E) (representation D k ω) := by
  change Tendsto (fun a : E => op (representation D k ω) a (oneVector D k ω)) (𝓝 0) (𝓝 (oneVector D k ω))
  simpa only [one_fixed] using
    (tendsto_const_nhds : Tendsto (fun _ : E => oneVector D k ω) (𝓝 0) (𝓝 (oneVector D k ω)))
lemma color_sum : ∑ i : Fin k, colorVector D k ω i = oneVector D k ω := by
  apply Lp.ext
  filter_upwards [Lp.coeFn_fun_finsetSum (Finset.univ) (colorVector D k ω),ae_all_iff.mpr (color_ae D k ω),
    one_ae D k ω] with η h1 h2 h3
  rw [h1,h3]
  simp_rw [h2,colorValue]
  simp
lemma projected_sum : ∑ i : Fin k, projected D k ω i = oneVector D k ω := by
  change ∑ i : Fin k, (space (A := E) (representation D k ω)).starProjection (colorVector D k ω i) = _
  rw [← map_sum,color_sum]
  exact (space (A := E) (representation D k ω)).starProjection_eq_self_iff.mpr (one_continuous D k ω)
end Ext.LabelModelD
end

/-! Real algebraic numbers as elements of `E` (new). -/
noncomputable section
namespace Ext
open PlaneFiveColor.Spectral
lemma isAlgebraic_ofReal {r : ℝ} (hr : IsAlgebraic ℚ r) : IsAlgebraic ℚ (r : ℂ) :=
  hr.algHom (Complex.ofRealHom.toRatAlgHom)
/-- a real algebraic number as an element of the field `E` of algebraic complex numbers -/
def realE (r : ℝ) (hr : IsAlgebraic ℚ r) : E := ⟨(r : ℂ),mem_algebraicClosure_iff.mpr (isAlgebraic_ofReal hr)⟩
@[simp] lemma realE_coe (r : ℝ) (hr : IsAlgebraic ℚ r) : ((realE r hr : E) : ℂ) = (r : ℂ) := rfl
/-- a set of real algebraic numbers is countable -/
lemma countable_of_algebraic {D : Set ℝ} (halg : ∀ r ∈ D, IsAlgebraic ℚ r) : D.Countable :=
  (Algebraic.countable ℚ ℝ).mono fun r hr => halg r hr
end Ext
end

noncomputable section
open MeasureTheory Filter Topology
open scoped InnerProductSpace
namespace Ext.LabelModelD
open PlaneFiveColor.Spectral EuclideanFiveColor EuclideanFiveColor.ContinuousVectors
variable (D : Set ℝ) (k : ℕ) (ω : Space D k)
lemma projected_continuous (i : Fin k) : projected D k ω i ∈ space (A := E) (representation D k ω) :=
  (space (A := E) (representation D k ω)).starProjection_apply_mem _
def fieldVector (i : Fin k) (x : ℂ) : Hilbert D k ω := orbitExtension (representation D k ω) (projected D k ω i) x
lemma fieldVector_continuous (i : Fin k) : Continuous (fieldVector D k ω i) :=
  orbitExtension_continuous _ (projected_continuous D k ω i)
lemma fieldVector_algebraic (i : Fin k) (a : E) :
    fieldVector D k ω i (a : ℂ) = op (representation D k ω) a (projected D k ω i) :=
  orbitExtension_eq _ (projected_continuous D k ω i) a
lemma fieldVector_sum (x : ℂ) : ∑ i : Fin k, fieldVector D k ω i x = oneVector D k ω := by
  refine dense_E.denseRange_val.induction_on x
    (isClosed_eq (continuous_finsetSum _ (fun i _ => fieldVector_continuous D k ω i)) continuous_const) ?_
  intro a
  simp only [fieldVector_algebraic,← map_sum,projected_sum,one_fixed]
lemma fieldVector_range (i : Fin k) (x : ℂ) : ∀ᵐ η ∂measure D k ω,
    (fieldVector D k ω i x η).im = 0 ∧ 0 ≤ (fieldVector D k ω i x η).re ∧ (fieldVector D k ω i x η).re ≤ 1 := by
  have he : clamp_lipschitz.compLp clamp_zero (fieldVector D k ω i x) = fieldVector D k ω i x := by
    refine dense_E.denseRange_val.induction_on x (isClosed_eq
      ((clamp_lipschitz.continuous_compLp clamp_zero).comp (fieldVector_continuous D k ω i))
      (fieldVector_continuous D k ω i)) ?_
    intro a
    rw [fieldVector_algebraic,op_apply]
    apply Lp.ext
    filter_upwards [clamp_lipschitz.coeFn_compLp clamp_zero
        (Lp.compMeasurePreserving (trans D k a) (trans_preserving D k ω a) (projected D k ω i)),
      Lp.coeFn_compMeasurePreserving (projected D k ω i) (trans_preserving D k ω a),
      (trans_preserving D k ω a).quasiMeasurePreserving.ae (projected_range D k ω i)] with η h1 h2 h3
    simp only [Function.comp_apply] at h1 h2
    rw [h1,h2]
    exact clamp_eq h3
  have ha := clamp_lipschitz.coeFn_compLp clamp_zero (fieldVector D k ω i x)
  rw [he] at ha
  filter_upwards [ha] with η hη
  exact hη.symm ▸ clamp_range _
/-- Generalises OpenAI's `fieldVector_unit_correlation` from the unit circle to the circle of any
radius `r ∈ D` (with `r` positive and algebraic). -/
lemma fieldVector_correlation (i : Fin k) {r : ℝ} (hr0 : 0 < r) (hra : IsAlgebraic ℚ r) (hrD : r ∈ D)
    (x : ℂ) (u : Circle) :
    ⟪fieldVector D k ω i x,fieldVector D k ω i (x+(r:ℂ)*(u:ℂ))⟫_ℂ = 0 := by
  rw [show ⟪fieldVector D k ω i x,fieldVector D k ω i (x+(r:ℂ)*(u:ℂ))⟫_ℂ =
    ⟪projected D k ω i,fieldVector D k ω i ((r:ℂ)*(u:ℂ))⟫_ℂ from
    orbitExtension_inner _ (projected_continuous D k ω i) (projected_continuous D k ω i) x ((r:ℂ)*(u:ℂ))]
  have hcl : IsClosed {w : ℂ | ⟪projected D k ω i,fieldVector D k ω i ((r:ℂ)*w)⟫_ℂ = 0} :=
    isClosed_eq (continuous_const.inner ((fieldVector_continuous D k ω i).comp (continuous_const_mul (r:ℂ))))
      continuous_const
  refine circle_closed_induction (p := fun w : ℂ => ⟪projected D k ω i,fieldVector D k ω i ((r:ℂ)*w)⟫_ℂ = 0)
    hcl (fun b => ?_) u
  have hrne : realE r hra ≠ 0 := by
    intro h
    have h' : (r : ℂ) = 0 := congrArg Subtype.val h
    exact hr0.ne' (Complex.ofReal_eq_zero.mp h')
  have hz0 : realE r hra * kCoeff b ≠ 0 := mul_ne_zero hrne (kCoeff_ne_zero b)
  have hzD : ‖((realE r hra * kCoeff b : E) : ℂ)‖ ∈ D := by
    change ‖(r:ℂ)*(b.val:ℂ)‖ ∈ D
    rw [norm_mul,Complex.norm_real,Real.norm_eq_abs,abs_of_pos hr0,Circle.norm_coe,mul_one]
    exact hrD
  change ⟪projected D k ω i,fieldVector D k ω i ((realE r hra * kCoeff b : E) : ℂ)⟫_ℂ = 0
  rw [fieldVector_algebraic]
  exact projected_correlation D k ω i _ hz0 hzD
end Ext.LabelModelD
end

noncomputable section
open MeasureTheory Filter Topology
namespace Ext.LabelModelD
open PlaneFiveColor.Spectral EuclideanFiveColor EuclideanFiveColor.ContinuousVectors
variable (D : Set ℝ) (k : ℕ) (ω : Space D k)
/-- Generalises OpenAI's `exists_label_field`: one measurable field, orthogonal along every circle
whose radius lies in `D`. -/
lemma exists_label_field (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r) :
    ∃ H : Fin k → ℂ × Space D k → ℝ,
    (∀ i, Measurable (H i)) ∧
    (∀ x : ℂ, ∀ᵐ η ∂measure D k ω, ∑ i, H i (x,η) = 1) ∧
    (∀ r ∈ D, ∀ x : ℂ, ∀ u : Circle, ∀ᵐ η ∂measure D k ω, ∀ i, H i (x,η)*H i (x+(r:ℂ)*(u:ℂ),η) = 0) := by
  have hex (i : Fin k) :=
    exists_measurable_field (measure D k ω) (fieldVector D k ω i) (fieldVector_continuous D k ω i)
  choose F hm hae using hex
  refine ⟨fun i p => (F i p).re,fun i => Complex.continuous_re.measurable.comp (hm i),?_,?_⟩
  · intro x
    have hs := Lp.coeFn_fun_finsetSum Finset.univ (fun i : Fin k => fieldVector D k ω i x)
    rw [fieldVector_sum] at hs
    filter_upwards [ae_all_iff.mpr (fun i => hae i x),hs,one_ae D k ω] with η he hs h1
    have hc : ∑ i : Fin k, F i (x,η) = 1 := by
      simp_rw [he]
      exact hs.symm.trans h1
    have hre := congrArg Complex.re hc
    simpa only [map_sum,Complex.re_sum,Complex.one_re] using hre
  · intro r hr x u
    apply ae_all_iff.mpr
    intro i
    have hz := positive_orthogonal (measure D k ω) (fieldVector D k ω i x) (fieldVector D k ω i (x+(r:ℂ)*(u:ℂ)))
      ((fieldVector_range D k ω i x).mono fun _ h => ⟨h.1,h.2.1⟩)
      ((fieldVector_range D k ω i (x+(r:ℂ)*(u:ℂ))).mono fun _ h => ⟨h.1,h.2.1⟩)
      (fieldVector_correlation D k ω i (hpos r hr) (halg r hr) hr x u)
    filter_upwards [hae i x,hae i (x+(r:ℂ)*(u:ℂ)),hz] with η h1 h2 h3
    simpa only [h1,h2] using h3
end Ext.LabelModelD
end

noncomputable section
open MeasureTheory Filter
namespace Ext
open EuclideanFiveColor
variable {Ω S : Type*} [MeasurableSpace Ω] [MeasurableSpace S]
  (μ : Measure Ω) [IsProbabilityMeasure μ] (ν : Measure S) [SFinite ν] {k : ℕ} [NeZero k]
/-- Generalises OpenAI's `EuclideanFiveColor.sample_field` from one displacement map `d` to a
countable family `d j` of displacement maps: one sample of the label configuration works for all
of them at once (a countable intersection of conull sets is conull). -/
lemma sample_field_countable {ι : Type*} [Countable ι] (d : ι → S → ℂ) (hd : ∀ j, Measurable (d j))
    (H : Fin k → ℂ × Ω → ℝ) (hm : ∀ i, Measurable (H i))
    (hp : ∀ x : ℂ, ∀ᵐ ω ∂μ, ∑ i, H i (x,ω) = 1)
    (ho : ∀ j, ∀ x : ℂ, ∀ s : S, ∀ᵐ ω ∂μ, ∀ i, H i (x,ω)*H i (x+d j s,ω) = 0) :
    ∃ c : ℂ → Fin k, Measurable c ∧ ∀ j, ∀ᵐ p ∂(volume : Measure ℂ).prod ν, c p.1 ≠ c (p.1+d j p.2) := by
  let c : ℂ × Ω → Fin k := positiveSelect H
  have hc : Measurable c := positiveSelect_measurable H hm
  have ha (j : ι) (p : ℂ × S) : ∀ᵐ ω ∂μ, c (p.1,ω) ≠ c (p.1+d j p.2,ω) := by
    filter_upwards [hp p.1,hp (p.1+d j p.2),ho j p.1 p.2] with ω h1 h2 h3
    intro he
    have hpos1 := positiveSelect_pos H (p.1,ω) h1
    have hpos2 := positiveSelect_pos H (p.1+d j p.2,ω) h2
    change 0 < H (c (p.1,ω)) (p.1,ω) at hpos1
    change 0 < H (c (p.1+d j p.2,ω)) (p.1+d j p.2,ω) at hpos2
    rw [← he] at hpos2
    exact (mul_pos hpos1 hpos2).ne' (h3 _)
  have hh (j : ι) : ∀ᵐ ω ∂μ, ∀ᵐ p ∂(volume : Measure ℂ).prod ν, c (p.1,ω) ≠ c (p.1+d j p.2,ω) := by
    have hm1 : Measurable (fun q : (ℂ × S) × Ω => c (q.1.1,q.2)) :=
      hc.comp ((measurable_fst.fst).prodMk measurable_snd)
    have hm2 : Measurable (fun q : (ℂ × S) × Ω => c (q.1.1+d j q.1.2,q.2)) :=
      hc.comp (((measurable_fst.fst).add ((hd j).comp measurable_fst.snd)).prodMk measurable_snd)
    exact (Measure.ae_ae_comm (measurableSet_eq_fun hm1 hm2).compl).mp (.of_forall (ha j))
  obtain ⟨ω,hω⟩ := (ae_all_iff.mpr hh).exists
  exact ⟨fun x => c (x,ω),hc.comp (measurable_id.prodMk measurable_const),hω⟩
end Ext
end

noncomputable section
open MeasureTheory
namespace Ext
open EuclideanFiveColor
/-- **Main theorem** (generalises OpenAI's `EuclideanFiveColor.proper_to_borel_weak` from the single
forbidden distance `1` to a set `D` of positive real algebraic distances). If no two points of the
plane at a distance belonging to `D` receive the same colour under `c`, then there is a measurable
colouring `d` with the same number of colours such that for every `r ∈ D`, for almost every point
`x` and almost every direction `t`, `d x ≠ d (x + r * direction t)`. -/
theorem proper_to_borel_weak_distances {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    {c : ℂ → Fin k} (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure,
      d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2) := by
  have : Countable D := (countable_of_algebraic halg).to_subtype
  let ω : LabelModelD.Space D k := ⟨fun x => c x,fun x y h => hc x y h⟩
  obtain ⟨H,hm,hp,ho⟩ := LabelModelD.exists_label_field D k ω hpos halg
  obtain ⟨d,hd,hw⟩ := sample_field_countable (LabelModelD.measure D k ω) angularMeasure
    (fun (r : D) (t : Angles) => (r.1:ℂ) * direction t)
    (fun _ => measurable_const.mul direction_continuous.measurable)
    H hm hp (fun r x t => ho r.1 r.2 x (AddCircle.toCircle t))
  exact ⟨d,hd,fun r hr => hw ⟨r,hr⟩⟩
end Ext
end

end OAI
