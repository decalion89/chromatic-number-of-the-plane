import OAI.Analysis.PlaneSpectrum.Representations

/-!
# The label model for a set `D` of forbidden distances

This file is a copy, with changes, of the part of OpenAI's file
`OAI/Analysis/PlaneSpectrum/Representations.lean` that lives in the namespace
`PlaneFiveColor.LabelModel` (the definitions `Proper`, `Space`, `trans`, `rotate`, the stationary
measure `measure`, the preservation lemmas, the Hilbert space `Hilbert`, the unitary representations
`translation`, `rotation`, `representation`, and `op_apply`, `rotation_covariance`).

The only change of substance is the notion of properness: OpenAI's
`Proper k ω := ∀ x y : E, ‖(x:ℂ)-(y:ℂ)‖ = 1 → ω x ≠ ω y` (one forbidden distance `1`) becomes
`Proper D k ω := ∀ x y : E, ‖(x:ℂ)-(y:ℂ)‖ ∈ D → ω x ≠ ω y` for an arbitrary set `D : Set ℝ`.
Translations and rotations by elements of `K` (algebraic points of the unit circle) preserve all
distances, so every construction goes through verbatim; every declaration takes the extra explicit
argument `D`. The generic results of OpenAI's library (`StationaryMean`, `Koopman`,
`ContinuousVectors`, ...) are used as they are, through the import.
-/

namespace OAI

noncomputable section
open MeasureTheory
namespace Ext.LabelModelD
open PlaneFiveColor.Spectral
variable (D : Set ℝ) (k : ℕ)
/-- (changed from OpenAI's `PlaneFiveColor.LabelModel.Proper`: the distance `1` is replaced by the set `D`) -/
def Proper (ω : E → Fin k) : Prop := ∀ x y : E, ‖(x:ℂ)-(y:ℂ)‖ ∈ D → ω x ≠ ω y
abbrev Space := {ω : E → Fin k // Proper D k ω}
lemma closed_proper : IsClosed {ω : E → Fin k | Proper D k ω} := by
  unfold Proper
  simp only [Set.ofPred_forall]
  refine isClosed_iInter fun (x:E) => isClosed_iInter fun (y:E) => isClosed_iInter fun _ => ?_
  have hc : Continuous (fun ω : E → Fin k => (ω x,ω y)) := (continuous_apply x).prodMk (continuous_apply y)
  exact (isClosed_discrete {p : Fin k × Fin k | p.1 ≠ p.2}).preimage hc
instance : CompactSpace (Space D k) := isCompact_iff_compactSpace.mp (closed_proper D k).isCompact
def trans (a : E) (ω : Space D k) : Space D k := ⟨fun x => ω.val (x+a),by
  intro x y h
  apply ω.property (x+a) (y+a)
  change ‖(x:ℂ)+(a:ℂ)-((y:ℂ)+(a:ℂ))‖ ∈ D
  simpa only [add_sub_add_right_eq_sub] using h⟩
def rotate (b : K) (ω : Space D k) : Space D k := ⟨fun x => ω.val (kCoeff b * x),by
  intro x y h
  apply ω.property (kCoeff b*x) (kCoeff b*y)
  change ‖(b.val:ℂ)*(x:ℂ)-(b.val:ℂ)*(y:ℂ)‖ ∈ D
  rw [← mul_sub,norm_mul,Circle.norm_coe,one_mul]
  exact h⟩
@[simp] lemma trans_apply (a : E) (ω : Space D k) (x : E) : (trans D k a ω).val x = ω.val (x+a) := rfl
@[simp] lemma rotate_apply (b : K) (ω : Space D k) (x : E) : (rotate D k b ω).val x = ω.val (kCoeff b*x) := rfl
lemma trans_continuous (a : E) : Continuous (trans D k a) := by
  apply Continuous.subtype_mk
  exact continuous_pi fun x => (continuous_apply (x+a)).comp continuous_subtype_val
lemma rotate_continuous (b : K) : Continuous (rotate D k b) := by
  apply Continuous.subtype_mk
  exact continuous_pi fun x => (continuous_apply (kCoeff b*x)).comp continuous_subtype_val
@[simp] lemma trans_zero (ω : Space D k) : trans D k 0 ω = ω := by apply Subtype.ext; funext x; simp
lemma trans_add (a b : E) (ω : Space D k) : trans D k a (trans D k b ω) = trans D k (a+b) ω := by
  apply Subtype.ext; funext x; simp [add_assoc]
@[simp] lemma rotate_one (ω : Space D k) : rotate D k 1 ω = ω := by
  apply Subtype.ext
  funext x
  have h1 : kCoeff 1 = 1 := rfl
  simp [h1]
lemma rotate_mul (a b : K) (ω : Space D k) : rotate D k a (rotate D k b ω) = rotate D k (a*b) ω := by
  apply Subtype.ext
  funext x
  change ω.val (kCoeff b * (kCoeff a * x)) = ω.val (kCoeff (a*b)*x)
  have he : kCoeff (a*b) = kCoeff a * kCoeff b := by rfl
  rw [he]; ring_nf
lemma covariance (a : E) (b : K) (ω : Space D k) :
    trans D k a (rotate D k b ω) = rotate D k b (trans D k (kCoeff b*a) ω) := by
  apply Subtype.ext
  funext x
  simp only [trans_apply,rotate_apply,mul_add]
def transHomeo (a : E) : Space D k ≃ₜ Space D k where
  toFun := trans D k a
  invFun := trans D k (-a)
  left_inv ω := by rw [trans_add,neg_add_cancel,trans_zero]
  right_inv ω := by rw [trans_add,add_neg_cancel,trans_zero]
  continuous_toFun := trans_continuous D k a
  continuous_invFun := trans_continuous D k (-a)
def rotateHomeo (b : K) : Space D k ≃ₜ Space D k where
  toFun := rotate D k b
  invFun := rotate D k b⁻¹
  left_inv ω := by rw [rotate_mul,inv_mul_cancel,rotate_one]
  right_inv ω := by rw [rotate_mul,mul_inv_cancel,rotate_one]
  continuous_toFun := rotate_continuous D k b
  continuous_invFun := rotate_continuous D k b⁻¹
def tr (a : Additive DiscreteE) : Space D k → Space D k := trans D k (Multiplicative.toAdd (Additive.toMul a))
def ro (b : Rot) : Space D k → Space D k := rotate D k (rotMul b)
def measure (ω : Space D k) : Measure (Space D k) :=
  EuclideanFiveColor.StationaryMean.measure (tr D k) (ro D k) ω
instance (ω : Space D k) : IsProbabilityMeasure (measure D k ω) :=
  inferInstanceAs (IsProbabilityMeasure (EuclideanFiveColor.StationaryMean.measure (tr D k) (ro D k) ω))
lemma trans_preserving (ω : Space D k) (a : E) : MeasurePreserving (trans D k a) (measure D k ω) (measure D k ω) := by
  apply EuclideanFiveColor.StationaryMean.preserving_of_average (tr D k) (ro D k) ω (transHomeo D k a)
  intro f
  apply EuclideanFiveColor.StationaryMean.average_translation (tr D k) (ro D k) ω
    (fun a b x => trans_add D k _ _ x)
    (fun b a => Additive.ofMul (Multiplicative.ofAdd (kCoeff (rotMul b)*Multiplicative.toAdd (Additive.toMul a))))
    (fun a b x => covariance D k _ _ x) (Additive.ofMul (Multiplicative.ofAdd a)) f (trans_continuous D k a)
lemma rotate_preserving (ω : Space D k) (b : K) : MeasurePreserving (rotate D k b) (measure D k ω) (measure D k ω) := by
  apply EuclideanFiveColor.StationaryMean.preserving_of_average (tr D k) (ro D k) ω (rotateHomeo D k b)
  intro f
  exact EuclideanFiveColor.StationaryMean.average_rotation (tr D k) (ro D k) ω
    (fun a b x => rotate_mul D k _ _ x) (mulRot b) f (rotate_continuous D k b)
end Ext.LabelModelD
end

noncomputable section
open MeasureTheory
namespace Ext.LabelModelD
open PlaneFiveColor.Spectral EuclideanFiveColor EuclideanFiveColor.ContinuousVectors
variable (D : Set ℝ) (k : ℕ) (ω : Space D k)
abbrev Hilbert := Lp ℂ 2 (measure D k ω)
def translation (a : E) : unitary (Hilbert D k ω →L[ℂ] Hilbert D k ω) :=
  Koopman.unitary (measure D k ω) (transHomeo D k a).toMeasurableEquiv (trans_preserving D k ω a)
def rotation (b : K) : unitary (Hilbert D k ω →L[ℂ] Hilbert D k ω) :=
  Koopman.unitary (measure D k ω) (rotateHomeo D k b).toMeasurableEquiv (rotate_preserving D k ω b)
lemma translation_apply (a : E) (v : Hilbert D k ω) :
    (translation D k ω a : Hilbert D k ω →L[ℂ] Hilbert D k ω) v =
      Lp.compMeasurePreserving (trans D k a) (trans_preserving D k ω a) v := rfl
lemma rotation_apply (b : K) (v : Hilbert D k ω) :
    (rotation D k ω b : Hilbert D k ω →L[ℂ] Hilbert D k ω) v =
      Lp.compMeasurePreserving (rotate D k b) (rotate_preserving D k ω b) v := rfl
def representation : DiscreteE →* unitary (Hilbert D k ω →L[ℂ] Hilbert D k ω) where
  toFun a := translation D k ω (Multiplicative.toAdd a)
  map_one' := by
    apply Subtype.ext
    apply ContinuousLinearMap.ext
    intro v
    change Lp.compMeasurePreserving (trans D k 0) (trans_preserving D k ω 0) v = v
    have he : trans D k 0 = id := funext (trans_zero D k)
    simp only [he,Lp.compMeasurePreserving_id_apply]
  map_mul' a b := by
    apply Subtype.ext
    apply ContinuousLinearMap.ext
    intro v
    change Lp.compMeasurePreserving (trans D k (Multiplicative.toAdd a+Multiplicative.toAdd b)) _ v =
      Lp.compMeasurePreserving (trans D k (Multiplicative.toAdd a)) _
        (Lp.compMeasurePreserving (trans D k (Multiplicative.toAdd b)) _ v)
    rw [← Lp.compMeasurePreserving_comp_apply]
    congr 2
    funext x
    exact (congrArg (fun c => trans D k c x) (add_comm _ _)).trans (trans_add D k _ _ x).symm
lemma op_apply (a : E) (v : Hilbert D k ω) : op (representation D k ω) a v =
    Lp.compMeasurePreserving (trans D k a) (trans_preserving D k ω a) v := rfl
lemma rotation_covariance (b : K) (a : E) (v : Hilbert D k ω) :
    (rotation D k ω b : Hilbert D k ω →L[ℂ] Hilbert D k ω) (op (representation D k ω) a v) =
      op (representation D k ω) (kCoeff b*a) ((rotation D k ω b : Hilbert D k ω →L[ℂ] Hilbert D k ω) v) := by
  simp only [rotation_apply,op_apply,← Lp.compMeasurePreserving_comp_apply]
  congr 2
  funext x
  exact covariance D k a b x
end Ext.LabelModelD
end

end OAI
