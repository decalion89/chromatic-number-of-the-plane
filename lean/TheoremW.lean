import Mathlib.Topology.Instances.AddCircle.Defs
import Mathlib.Data.ZMod.Basic
import Mathlib.Basic.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity

/-!
# Theorem W, case `k = 1`, for finite abelian groups

Let `G` be a finite abelian group, `S ⊆ G` a symmetric set and `c : G → ZMod 3` a proper
colouring of the Cayley graph `Cay(G, S)`, i.e. `c (g + s) ≠ c g` for all `g ∈ G`, `s ∈ S`
(this forces `0 ∉ S`). Then there is an additive character
`ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ)` with `ξ s ∈ [1/3, 2/3] (mod 1)` for `s ∈ S`.

Proof (winding numbers).
* `sig c g s = σ(g, s) ∈ {±1}` is the integer lift of `c (g + s) - c g ∈ ZMod 3`.
* `W c g L` is the winding sum of the walk from `g` with steps `L`.
* Squares do not wind (`sgn3_square`), hence `W c g` is invariant under permutations of the
  steps (`W_perm`); for closed walks this gives invariance under translations by `S`
  (`W_shift`) and so by `closure S` (`W_shift_closure`).
* For closed walks `W ≡ 0 (mod 3)` (telescoping in `ZMod 3`) and `W ≡ length (mod 2)`.
* Averaging over `H = closure S`: with `F s = ∑_{γ ∈ H} σ(γ, s)` and
  `ξ₀ s = (F s / (3 |H|) + 1) / 2 ∈ [1/3, 2/3]`, every closed walk has `∑ ξ₀ ∈ ℤ`.
* Hence `ξ₀` induces a homomorphism `H → ℝ/ℤ`.
-/

set_option autoImplicit false

namespace TheoremW

open Finset

/-! ### The `±1` lift of a nonzero difference in `ZMod 3` -/

/-- `sgn3 a b ∈ {1, -1}` is the integer lift of `b - a` (when `a ≠ b`). -/
def sgn3 (a b : ZMod 3) : ℤ := if b - a = 1 then 1 else -1

lemma sgn3_eq_one_or (a b : ZMod 3) : sgn3 a b = 1 ∨ sgn3 a b = -1 := by
  unfold sgn3; split_ifs <;> simp

lemma zmod3_cases : ∀ x : ZMod 3, x = 0 ∨ x = 1 ∨ x = -1 := by decide

lemma sgn3_cast {a b : ZMod 3} (h : b ≠ a) : ((sgn3 a b : ℤ) : ZMod 3) = b - a := by
  unfold sgn3
  split_ifs with h1
  · rw [h1, Int.cast_one]
  · rcases zmod3_cases (b - a) with h2 | h2 | h2
    · exact absurd (sub_eq_zero.1 h2) h
    · exact absurd h2 h1
    · rw [h2, Int.cast_neg, Int.cast_one]

lemma three_dvd_of_cast {x : ℤ} (h : ((x : ℤ) : ZMod 3) = 0) : (3 : ℤ) ∣ x := by
  have := (ZMod.intCast_zmod_eq_zero_iff_dvd x 3).1 h
  exact_mod_cast this

lemma sgn3_antisymm {a b : ZMod 3} (h : b ≠ a) : sgn3 b a = - sgn3 a b := by
  have h3 : (3 : ℤ) ∣ sgn3 b a + sgn3 a b := three_dvd_of_cast (by
    push_cast
    rw [sgn3_cast h, sgn3_cast h.symm]
    ring)
  rcases sgn3_eq_one_or a b with h4 | h4 <;> rcases sgn3_eq_one_or b a with h5 | h5 <;> omega

/-- Squares do not wind: both sides lie in `{-2, 0, 2}` and agree mod `3`. -/
lemma sgn3_square {a b c d : ZMod 3} (hab : b ≠ a) (hbd : d ≠ b) (hac : c ≠ a) (hcd : d ≠ c) :
    sgn3 a b + sgn3 b d = sgn3 a c + sgn3 c d := by
  have h3 : (3 : ℤ) ∣ sgn3 a b + sgn3 b d - sgn3 a c - sgn3 c d := three_dvd_of_cast (by
    push_cast
    rw [sgn3_cast hab, sgn3_cast hbd, sgn3_cast hac, sgn3_cast hcd]
    ring)
  rcases sgn3_eq_one_or a b with h1 | h1 <;> rcases sgn3_eq_one_or b d with h2 | h2 <;>
    rcases sgn3_eq_one_or a c with h4 | h4 <;> rcases sgn3_eq_one_or c d with h5 | h5 <;> omega

/-! ### Winding sums of walks -/

variable {G : Type*} [AddCommGroup G]

/-- `sig c g s = σ(g, s) ∈ {±1}`, the integer lift of `c (g + s) - c g`. -/
def sig (c : G → ZMod 3) (g s : G) : ℤ := sgn3 (c g) (c (g + s))

lemma sig_eq_one_or (c : G → ZMod 3) (g s : G) : sig c g s = 1 ∨ sig c g s = -1 :=
  sgn3_eq_one_or _ _

/-- `σ(g + s, -s) = -σ(g, s)`. -/
lemma sig_neg {S : Set G} {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {s : G} (hs : s ∈ S) (g : G) : sig c (g + s) (-s) = - sig c g s := by
  unfold sig
  rw [add_neg_cancel_right]
  exact sgn3_antisymm (hc g s hs)

/-- The winding sum `W(g, L) = ∑_j σ(g + s₁ + ⋯ + s_{j-1}, s_j)` of the walk from `g` with
steps `L = [s₁, …, s_N]`. -/
def W (c : G → ZMod 3) : G → List G → ℤ
  | _, [] => 0
  | g, s :: L => sig c g s + W c (g + s) L

@[simp] lemma W_nil (c : G → ZMod 3) (g : G) : W c g [] = 0 := rfl

@[simp] lemma W_cons (c : G → ZMod 3) (g s : G) (L : List G) :
    W c g (s :: L) = sig c g s + W c (g + s) L := rfl

lemma W_append (c : G → ZMod 3) (g : G) (L₁ L₂ : List G) :
    W c g (L₁ ++ L₂) = W c g L₁ + W c (g + L₁.sum) L₂ := by
  induction L₁ generalizing g with
  | nil => simp
  | cons x L ih => simp [ih, add_assoc]

/-- The winding sum is invariant under permutations of the steps (squares do not wind). -/
lemma W_perm {S : Set G} {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {L L' : List G} (hp : L.Perm L') : ∀ g, (∀ x ∈ L, x ∈ S) → W c g L = W c g L' := by
  induction hp with
  | nil => intro g _; rfl
  | cons x _ ih =>
    intro g hL
    simp only [W_cons]
    rw [ih (g + x) (fun y hy => hL y (List.mem_cons_of_mem x hy))]
  | swap x y l =>
    intro g hL
    have hx : x ∈ S := hL x (by simp)
    have hy : y ∈ S := hL y (by simp)
    have hyx : g + y + x = g + x + y := add_right_comm g y x
    have hbd : c (g + x + y) ≠ c (g + y) := by rw [← hyx]; exact hc (g + y) x hx
    have h := sgn3_square (hc g y hy) hbd (hc g x hx) (hc (g + x) y hy)
    simp only [W_cons, sig, hyx]
    omega
  | trans h₁ _ ih₁ ih₂ =>
    intro g hL
    rw [ih₁ g hL, ih₂ g (fun y hy => hL y (h₁.symm.subset hy))]

/-- Telescoping in `ZMod 3`. -/
lemma W_cast {S : Set G} {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) :
    ∀ (L : List G) (g : G), (∀ x ∈ L, x ∈ S) →
      ((W c g L : ℤ) : ZMod 3) = c (g + L.sum) - c g
  | [], g, _ => by simp
  | s :: L, g, hL => by
    have hs : s ∈ S := hL s (by simp)
    have ih := W_cast hc L (g + s) (fun x hx => hL x (List.mem_cons_of_mem s hx))
    rw [W_cons, Int.cast_add, ih, sig, sgn3_cast (hc g s hs), List.sum_cons,
      show g + (s + L.sum) = g + s + L.sum from (add_assoc _ _ _).symm]
    ring

lemma W_parity (c : G → ZMod 3) : ∀ (L : List G) (g : G), W c g L % 2 = (L.length : ℤ) % 2
  | [], g => by simp
  | s :: L, g => by
    have ih := W_parity c L (g + s)
    rw [W_cons, List.length_cons]
    push_cast
    rcases sig_eq_one_or c g s with h | h <;> omega

/-- A closed walk has winding sum `3 w` with `w ≡ N (mod 2)`. -/
lemma W_closed {S : Set G} {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {L : List G} (hL : ∀ x ∈ L, x ∈ S) (h0 : L.sum = 0) (g : G) :
    ∃ m : ℤ, W c g L = 3 * (2 * m - L.length) := by
  have h3 : (3 : ℤ) ∣ W c g L := by
    apply three_dvd_of_cast
    rw [W_cast hc L g hL, h0, add_zero, sub_self]
  have h2 := W_parity c L g
  exact ⟨(W c g L / 3 + L.length) / 2, by omega⟩

/-- Moving the base point of a closed walk by `s ∈ S` does not change the winding sum. -/
lemma W_shift {S : Set G} (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod 3}
    (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) {L : List G} (hL : ∀ x ∈ L, x ∈ S)
    (h0 : L.sum = 0) {s : G} (hs : s ∈ S) (g : G) : W c (g + s) L = W c g L := by
  have hmem : ∀ x ∈ L ++ [s, -s], x ∈ S := by
    intro x hx
    rw [List.mem_append] at hx
    rcases hx with hx | hx
    · exact hL x hx
    · simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
      rcases hx with rfl | rfl
      · exact hs
      · exact hS s hs
  have h := W_perm hc (List.perm_middle (l₁ := L) (a := s) (l₂ := [-s])) g hmem
  simp only [W_append, W_cons, W_nil, h0, add_zero, sig_neg hc hs] at h
  linarith

/-- The winding sum of a closed walk only depends on the coset of the base point. -/
lemma W_shift_closure {S : Set G} (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod 3}
    (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) {L : List G} (hL : ∀ x ∈ L, x ∈ S)
    (h0 : L.sum = 0) {γ : G} (hγ : γ ∈ AddSubgroup.closure S) :
    ∀ g, W c (g + γ) L = W c g L := by
  induction hγ using AddSubgroup.closure_induction with
  | mem x hx => exact W_shift hS hc hL h0 hx
  | zero => intro g; rw [add_zero]
  | add x y _ _ ihx ihy => intro g; rw [← add_assoc, ihy, ihx]
  | neg x _ ih => intro g; rw [← ih (g + -x), neg_add_cancel_right]

/-! ### Lists of generators -/

lemma list_sum_map_neg (L : List G) : (L.map (fun z => -z)).sum = -L.sum := by
  induction L with
  | nil => simp
  | cons x L ih => simp only [List.map_cons, List.sum_cons, ih, neg_add]

/-- Every element of the closure of a symmetric set `S` is the sum of a list of elements
of `S`. -/
lemma exists_list_of_mem_closure {S : Set G} (hS : ∀ s ∈ S, -s ∈ S) {γ : G}
    (hγ : γ ∈ AddSubgroup.closure S) : ∃ L : List G, (∀ x ∈ L, x ∈ S) ∧ L.sum = γ := by
  induction hγ using AddSubgroup.closure_induction with
  | mem x hx => exact ⟨[x], by simpa using hx, by simp⟩
  | zero => exact ⟨[], by simp, by simp⟩
  | add x y _ _ ihx ihy =>
    obtain ⟨L₁, h₁, e₁⟩ := ihx
    obtain ⟨L₂, h₂, e₂⟩ := ihy
    exact ⟨L₁ ++ L₂, fun z hz => (List.mem_append.1 hz).elim (h₁ z) (h₂ z), by
      rw [List.sum_append, e₁, e₂]⟩
  | neg x _ ih =>
    obtain ⟨L, h, e⟩ := ih
    refine ⟨L.map (fun z => -z), fun z hz => ?_, by rw [list_sum_map_neg, e]⟩
    obtain ⟨y, hy, rfl⟩ := List.mem_map.1 hz
    exact hS y (h y hy)

/-! ### Averaging over the finite subgroup `H = closure S` -/

section average

variable (H : AddSubgroup G) [Fintype H]

lemma sum_translate {M : Type*} [AddCommMonoid M] (φ : G → M) {p : G} (hp : p ∈ H) :
    ∑ γ : H, φ ((γ : G) + p) = ∑ γ : H, φ γ :=
  Fintype.sum_equiv (Equiv.addRight (⟨p, hp⟩ : H)) _ _ (fun γ => by simp)

/-- `F c s = ∑_{γ ∈ H} σ(γ, s)` (so `f(s) = F c s / |H|` is the average of `σ(·, s)`). -/
def F (c : G → ZMod 3) (s : G) : ℤ := ∑ γ : H, sig c γ s

lemma F_neg {S : Set G} (hSH : S ⊆ H) {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {s : G} (hs : s ∈ S) : F H c (-s) = - F H c s := by
  unfold F
  rw [← sum_translate H (fun h => sig c h (-s)) (hSH hs)]
  simp only [sig_neg hc hs, Finset.sum_neg_distrib]

lemma F_le (c : G → ZMod 3) (s : G) : F H c s ≤ Fintype.card H := by
  unfold F
  calc ∑ γ : H, sig c γ s ≤ ∑ _γ : H, (1 : ℤ) :=
        Finset.sum_le_sum (fun γ _ => by rcases sig_eq_one_or c γ s with h | h <;> omega)
    _ = Fintype.card H := by simp

lemma neg_card_le_F (c : G → ZMod 3) (s : G) : -(Fintype.card H : ℤ) ≤ F H c s := by
  unfold F
  calc -(Fintype.card H : ℤ) = ∑ _γ : H, (-1 : ℤ) := by simp
    _ ≤ ∑ γ : H, sig c γ s :=
        Finset.sum_le_sum (fun γ _ => by rcases sig_eq_one_or c γ s with h | h <;> omega)

/-- Summing the winding sums over all base points of a coset. -/
lemma sum_W {S : Set G} (hSH : S ⊆ H) (c : G → ZMod 3) :
    ∀ (L : List G) (p : G), p ∈ H → (∀ x ∈ L, x ∈ S) →
      ∑ γ : H, W c ((γ : G) + p) L = (L.map (F H c)).sum
  | [], p, _, _ => by simp
  | s :: L, p, hp, hL => by
    have hs : s ∈ S := hL s (by simp)
    have ih := sum_W hSH c L (p + s) (H.add_mem hp (hSH hs))
      (fun x hx => hL x (List.mem_cons_of_mem s hx))
    simp only [W_cons, add_assoc, Finset.sum_add_distrib, List.map_cons, List.sum_cons, ih]
    congr 1
    exact sum_translate H (fun h => sig c h s) hp

/-- For a closed walk, `∑_j F(s_j) = |H| · W(0, L)`. -/
lemma sum_F_closed {S : Set G} (hSH : S ⊆ H) (hHS : H ≤ AddSubgroup.closure S)
    (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {L : List G} (hL : ∀ x ∈ L, x ∈ S) (h0 : L.sum = 0) :
    (L.map (F H c)).sum = Fintype.card H * W c 0 L := by
  rw [← sum_W H hSH c L 0 H.zero_mem hL]
  calc ∑ γ : H, W c ((γ : G) + 0) L = ∑ _γ : H, W c 0 L :=
        Finset.sum_congr rfl (fun γ _ => by
          simpa using W_shift_closure hS hc hL h0 (hHS γ.2) 0)
    _ = Fintype.card H * W c 0 L := by simp

/-- `ξ₀(s) = (f(s)/3 + 1)/2` with `f(s) = F(s)/|H|`. -/
noncomputable def xi0 (c : G → ZMod 3) (s : G) : ℝ :=
  ((F H c s : ℝ) + 3 * Fintype.card H) / (6 * Fintype.card H)

lemma xi0_eq (c : G → ZMod 3) :
    xi0 H c = fun s => ((F H c s : ℝ) + 3 * Fintype.card H) / (6 * Fintype.card H) := rfl

lemma card_pos_real : (0 : ℝ) < Fintype.card H := by
  exact_mod_cast Fintype.card_pos_iff.2 ⟨0⟩

lemma xi0_bounds (c : G → ZMod 3) (s : G) : 1/3 ≤ xi0 H c s ∧ xi0 H c s ≤ 2/3 := by
  have hn := card_pos_real H
  have h1 : (F H c s : ℝ) ≤ Fintype.card H := by exact_mod_cast F_le H c s
  have h2 : -(Fintype.card H : ℝ) ≤ F H c s := by exact_mod_cast neg_card_le_F H c s
  unfold xi0
  constructor
  · rw [le_div_iff₀ (by positivity)]; linarith
  · rw [div_le_iff₀ (by positivity)]; linarith

lemma xi0_neg {S : Set G} (hSH : S ⊆ H) {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {s : G} (hs : s ∈ S) : xi0 H c (-s) = 1 - xi0 H c s := by
  have hn : (Fintype.card H : ℝ) ≠ 0 := (card_pos_real H).ne'
  unfold xi0
  rw [F_neg H hSH hc hs]
  push_cast
  field_simp
  ring

lemma list_sum_affine {α : Type*} (L : List α) (a : α → ℝ) (b d : ℝ) :
    (L.map (fun x => (a x + b) / d)).sum = ((L.map a).sum + L.length * b) / d := by
  induction L with
  | nil => simp
  | cons x L ih =>
    simp only [List.map_cons, List.sum_cons, List.length_cons, ih]
    push_cast
    ring

lemma cast_list_sum_map {α : Type*} (L : List α) (f : α → ℤ) :
    (((L.map f).sum : ℤ) : ℝ) = (L.map (fun x => (f x : ℝ))).sum := by
  induction L with
  | nil => simp
  | cons x L ih => simp [ih]

/-- The key integrality: along a closed walk, `∑ ξ₀ ∈ ℤ`. -/
lemma sum_xi0_closed {S : Set G} (hSH : S ⊆ H) (hHS : H ≤ AddSubgroup.closure S)
    (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {L : List G} (hL : ∀ x ∈ L, x ∈ S) (h0 : L.sum = 0) :
    ∃ m : ℤ, (L.map (xi0 H c)).sum = m := by
  obtain ⟨m, hm⟩ := W_closed hc hL h0 0
  refine ⟨m, ?_⟩
  have hn : (Fintype.card H : ℝ) ≠ 0 := (card_pos_real H).ne'
  have hF := sum_F_closed H hSH hHS hS hc hL h0
  have hF' : (L.map (fun s => (F H c s : ℝ))).sum
      = Fintype.card H * (3 * (2 * (m : ℝ) - L.length)) := by
    rw [← cast_list_sum_map, hF, hm]
    push_cast
    ring
  rw [xi0_eq, list_sum_affine L (fun s => (F H c s : ℝ)), hF']
  field_simp
  ring

lemma sum_xi0_neg {S : Set G} (hSH : S ⊆ H) {c : G → ZMod 3}
    (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) :
    ∀ (L : List G), (∀ x ∈ L, x ∈ S) →
      ((L.map (fun z => -z)).map (xi0 H c)).sum = L.length - (L.map (xi0 H c)).sum
  | [], _ => by simp
  | s :: L, hL => by
    have hs : s ∈ S := hL s (by simp)
    have ih := sum_xi0_neg hSH hc L (fun x hx => hL x (List.mem_cons_of_mem s hx))
    simp only [List.map_cons, List.sum_cons, List.length_cons, ih, xi0_neg H hSH hc hs]
    push_cast
    ring

/-- The value in `ℝ/ℤ` of a list of steps. -/
noncomputable def V (c : G → ZMod 3) (L : List G) : AddCircle (1 : ℝ) :=
  (((L.map (xi0 H c)).sum : ℝ) : AddCircle (1 : ℝ))

lemma V_append (c : G → ZMod 3) (L₁ L₂ : List G) :
    V H c (L₁ ++ L₂) = V H c L₁ + V H c L₂ := by
  simp only [V, List.map_append, List.sum_append, AddCircle.coe_add]

/-- Well-definedness: lists with the same sum have the same value in `ℝ/ℤ`. -/
lemma V_eq {S : Set G} (hSH : S ⊆ H) (hHS : H ≤ AddSubgroup.closure S)
    (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod 3} (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g)
    {L₁ L₂ : List G} (h₁ : ∀ x ∈ L₁, x ∈ S) (h₂ : ∀ x ∈ L₂, x ∈ S)
    (he : L₁.sum = L₂.sum) : V H c L₁ = V H c L₂ := by
  have hL : ∀ x ∈ L₁ ++ L₂.map (fun z => -z), x ∈ S := by
    intro x hx
    rcases List.mem_append.1 hx with hx | hx
    · exact h₁ x hx
    · obtain ⟨y, hy, rfl⟩ := List.mem_map.1 hx
      exact hS y (h₂ y hy)
  have h0 : (L₁ ++ L₂.map (fun z => -z)).sum = 0 := by
    rw [List.sum_append, list_sum_map_neg, he, add_neg_cancel]
  obtain ⟨m, hm⟩ := sum_xi0_closed H hSH hHS hS hc hL h0
  rw [List.map_append, List.sum_append, sum_xi0_neg H hSH hc L₂ h₂] at hm
  unfold V
  rw [← sub_eq_zero, ← AddCircle.coe_sub, AddCircle.coe_eq_zero_iff]
  refine ⟨m - L₂.length, ?_⟩
  rw [zsmul_eq_mul, mul_one]
  push_cast
  linarith

end average

/-! ### Main theorem -/

/-- Properness of the colouring forces `0 ∉ S`, so this hypothesis can be dropped. -/
lemma zero_not_mem_of_proper {S : Set G} {c : G → ZMod 3}
    (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) : (0 : G) ∉ S :=
  fun h0 => hc 0 0 h0 (by rw [add_zero])

/-- **Theorem W** (case `k = 1`, finite groups). Let `G` be a finite abelian group, `S ⊆ G`
symmetric and `c : G → ZMod 3` a proper colouring of `Cay(G, S)`. Then there is a character
`ξ : closure S →+ ℝ/ℤ` such that every `ξ s` (`s ∈ S`) is represented by a real number in
`[1/3, 2/3]`. (The hypothesis `0 ∉ S` is not needed: it follows from properness.) -/
theorem theoremW_finite {G : Type*} [AddCommGroup G] [Finite G]
    (S : Set G) (hS : ∀ s ∈ S, -s ∈ S)
    (c : G → ZMod 3) (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, 1/3 ≤ x ∧ x ≤ 2/3 ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  classical
  set H := AddSubgroup.closure S
  let _ : Fintype H := Fintype.ofFinite H
  have hSH : S ⊆ H := AddSubgroup.subset_closure
  have hHS : H ≤ AddSubgroup.closure S := le_rfl
  have hlist : ∀ γ : H, ∃ L : List G, (∀ x ∈ L, x ∈ S) ∧ L.sum = γ :=
    fun γ => exists_list_of_mem_closure hS γ.2
  choose Lf hLf hLs using hlist
  let ξf : H → AddCircle (1 : ℝ) := fun γ => V H c (Lf γ)
  have key : ∀ (L : List G), (∀ x ∈ L, x ∈ S) → ∀ γ : H, L.sum = γ → ξf γ = V H c L :=
    fun L hL γ he => V_eq H hSH hHS hS hc (hLf γ) hL ((hLs γ).trans he.symm)
  refine ⟨AddMonoidHom.mk' ξf (fun a b => ?_), fun s hs =>
    ⟨xi0 H c s, (xi0_bounds H c s).1, (xi0_bounds H c s).2, ?_⟩⟩
  · rw [key (Lf a ++ Lf b) (fun x hx => (List.mem_append.1 hx).elim (hLf a x) (hLf b x))
      (a + b) (by rw [List.sum_append, hLs a, hLs b, AddSubgroup.coe_add]), V_append]
  · show _ = ξf _
    rw [key [s] (by simpa using hs) _ (by simp)]
    simp [V]

/-- The statement with exactly the hypotheses of Theorem W (`S` a finite symmetric set with
`0 ∉ S`); the hypothesis `0 ∉ S` is redundant. -/
theorem theoremW_finite_finset {G : Type*} [AddCommGroup G] [Fintype G]
    (S : Finset G) (hS : ∀ s ∈ S, -s ∈ S) (_h0 : (0 : G) ∉ S)
    (c : G → ZMod 3) (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) :
    ∃ ξ : AddSubgroup.closure (S : Set G) →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, 1/3 ≤ x ∧ x ≤ 2/3 ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure (Finset.mem_coe.2 hs)⟩ :=
  theoremW_finite (S : Set G) (fun s hs => hS s hs) c (fun g s hs => hc g s hs)

/-- Sanity check: the hypotheses are satisfiable (`G = ZMod 3`, `S = {1, 2}`, `c = id`). -/
example : ∃ ξ : AddSubgroup.closure ({1, 2} : Set (ZMod 3)) →+ AddCircle (1 : ℝ),
    ∀ s (hs : s ∈ ({1, 2} : Set (ZMod 3))), ∃ x : ℝ, 1/3 ≤ x ∧ x ≤ 2/3 ∧
      (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  refine theoremW_finite ({1, 2} : Set (ZMod 3)) ?_ id ?_
  · intro s hs
    rcases hs with rfl | rfl
    · exact Or.inr (by decide : (-1 : ZMod 3) = 2)
    · exact Or.inl (by decide : (-2 : ZMod 3) = 1)
  · intro g s hs h
    have : s = 0 := by simpa using h
    rcases hs with rfl | rfl
    · exact absurd this (by decide)
    · exact absurd this (by decide)

end TheoremW
