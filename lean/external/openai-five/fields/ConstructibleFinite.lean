import Check.ConstructiblePlane
import Mathlib.Combinatorics.SimpleGraph.Finsubgraph

/-! A finite form of `sqrtClosure_plane_not_five_colourable`: some finite set of constructible points has no proper
5-colouring. The step from the infinite statement is the de Bruijn–Erdős compactness theorem, in Mathlib as
`SimpleGraph.nonempty_hom_of_forall_finite_subgraph_hom` (homomorphisms to the complete graph on `Fin 5`). -/

namespace OAI
noncomputable section
namespace EuclideanFiveColor

theorem sqrtClosure_finite_set_not_five_colourable :
    ∃ S : Set (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ),
      S.Finite ∧ ¬ ∃ c : S → Fin 5, ∀ x y : S, ‖((x : _) : ℂ) - ((y : _) : ℂ)‖ = 1 → c x ≠ c y := by
  by_contra hcon
  push Not at hcon
  apply sqrtClosure_plane_not_five_colourable
  let G : SimpleGraph (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ) :=
    { Adj := fun x y => ‖(x : ℂ) - (y : ℂ)‖ = 1
      symm := ⟨fun x y h => by simpa only [norm_sub_rev] using h⟩
      loopless := ⟨fun x h => by simp at h⟩ }
  have h : ∀ G' : G.Subgraph, G'.verts.Finite → G'.coe →g (⊤ : SimpleGraph (Fin 5)) := by
    intro G' hfin
    refine ⟨fun v => Classical.choose (hcon G'.verts hfin) v, fun {v w} hvw => ?_⟩
    rw [SimpleGraph.top_adj]
    exact Classical.choose_spec (hcon G'.verts hfin) v w (G'.adj_sub hvw)
  obtain ⟨φ⟩ := SimpleGraph.nonempty_hom_of_forall_finite_subgraph_hom h
  refine ⟨fun x => φ x, fun x y hxy => ?_⟩
  have hadj : G.Adj x y := hxy
  have := φ.map_rel hadj
  rwa [SimpleGraph.top_adj] at this

end EuclideanFiveColor
end
end OAI

#print axioms OAI.EuclideanFiveColor.sqrtClosure_finite_set_not_five_colourable
