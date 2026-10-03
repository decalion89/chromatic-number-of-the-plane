import Mathlib.Topology.Instances.AddCircle.Defs
import Mathlib.Data.ZMod.Basic
import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity

/-!
# Theorem W for circular cliques `K_{p/q}` with `p < 4q`, finite abelian groups

Let `G` be a finite abelian group, `S ⊆ G` a symmetric set and `c : G → ZMod p` a homomorphism
from `Cay(G, S)` to the circular clique `K_{p/q}`: along every edge `g → g + s` (`s ∈ S`) the
representative `δ(g, s) ∈ {0, …, p - 1}` of `c (g + s) - c g` lies in `[q, p - q]`. If
`p < 4q`, there is an additive character `ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ)` with
`ξ s ∈ [q/p, 1 - q/p] (mod 1)` for all `s ∈ S`:
* `theoremWplus_finite`: the statement for `p` odd and `2q < p < 4q`;
* `theoremWplus_general`: the same conclusion only assuming `0 < p < 4q` (oddness of `p` and
  `2q < p` are not used).

Proof (winding numbers, exactly as in `TheoremW.lean`, with the lift `δ` in place of `sgn3`; in
terms of the centred lift `σ = 2δ - p` of the informal proof, `δ = (σ + p) / 2`).
* `dl c g s = δ(g, s) ∈ [q, p - q]`, `δ(g, s) ≡ c (g + s) - c g (mod p)` and
  `δ(g + s, -s) = p - δ(g, s)` (`dl_bounds`, `dl_cast`, `dl_neg`).
* Squares do not wind (`dl_square`): `D = δ(g, s) + δ(g + s, t) - δ(g, t) - δ(g + t, s)` is
  `≡ 0 (mod p)` and `|D| ≤ 2 (p - 2q) < p` because `p < 4q`, so `D = 0`. Hence the winding
  sum `W c g L` of a walk is invariant under permutations of its steps (`W_perm`) and, for
  closed walks, under translations of the base point by `closure S` (`W_shift`,
  `W_shift_closure`).
* For closed walks `W ≡ 0 (mod p)` (telescoping in `ZMod p`, `W_closed`).
* Averaging over `H = closure S`: with `F s = ∑_{γ ∈ H} δ(γ, s)` and
  `ξ₀ s = F s / (p |H|) ∈ [q/p, 1 - q/p]`, every closed walk `s₁, …, s_N` has
  `∑ ξ₀ (s_j) = W / p ∈ ℤ`, and `ξ₀ (-s) = 1 - ξ₀ s`.
* Hence `ξ₀` induces a homomorphism `H → ℝ/ℤ` (`V_eq`, `AddMonoidHom.mk'`).
-/

set_option autoImplicit false

namespace TheoremWplus

open Finset

variable {G : Type*} [AddCommGroup G] {p q : ℕ}

/-- `c : G → ZMod p` is a homomorphism `Cay(G, S) → K_{p/q}`: along every edge `g → g + s`
(`s ∈ S`) the colour difference, read in `{0, …, p - 1}`, lies in `[q, p - q]`. -/
def IsCircHom (S : Set G) (p q : ℕ) (c : G → ZMod p) : Prop :=
  ∀ g, ∀ s ∈ S, q ≤ (c (g + s) - c g).val ∧ (c (g + s) - c g).val ≤ p - q

/-! ### The lift `δ(g, s) ∈ {0, …, p - 1}` of `c (g + s) - c g` -/

/-- `dl c g s = δ(g, s)`, the representative in `{0, …, p - 1}` of `c (g + s) - c g`. -/
def dl (c : G → ZMod p) (g s : G) : ℤ := ((c (g + s) - c g).val : ℤ)

lemma dl_cast [NeZero p] (c : G → ZMod p) (g s : G) :
    ((dl c g s : ℤ) : ZMod p) = c (g + s) - c g := by
  unfold dl
  rw [Int.cast_natCast, ZMod.natCast_zmod_val]

lemma dl_bounds {S : Set G} {c : G → ZMod p} (hc : IsCircHom S p q c) {g s : G} (hs : s ∈ S) :
    (q : ℤ) ≤ dl c g s ∧ dl c g s ≤ p - q := by
  obtain ⟨h1, h2⟩ := hc g s hs
  unfold dl
  constructor <;> omega

/-- `δ(g + s, -s) = p - δ(g, s)` (note `δ(g, s) ≠ 0` since `q ≥ 1`). -/
lemma dl_neg [NeZero p] {S : Set G} {c : G → ZMod p} (hc : IsCircHom S p q c) (h4 : p < 4 * q)
    {s : G} (hs : s ∈ S) (g : G) : dl c (g + s) (-s) = p - dl c g s := by
  have h1 := (hc g s hs).1
  have hlt := ZMod.val_lt (c (g + s) - c g)
  have hne : c (g + s) - c g ≠ 0 := by
    intro h0
    rw [h0, ZMod.val_zero] at h1
    omega
  unfold dl
  rw [add_neg_cancel_right, ← neg_sub (c (g + s)) (c g), ZMod.neg_val, ite_eq_right hne]
  omega

/-- Squares do not wind: the difference of the two sides is `≡ 0 (mod p)` and has absolute
value at most `2 (p - 2q) < p`. -/
lemma dl_square [NeZero p] {S : Set G} {c : G → ZMod p} (hc : IsCircHom S p q c)
    (h4 : p < 4 * q) {x y : G} (hx : x ∈ S) (hy : y ∈ S) (g : G) :
    dl c g x + dl c (g + x) y = dl c g y + dl c (g + y) x := by
  have hd : (p : ℤ) ∣ dl c g x + dl c (g + x) y - (dl c g y + dl c (g + y) x) := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    push_cast
    rw [dl_cast, dl_cast, dl_cast, dl_cast, add_right_comm g x y]
    ring
  have b1 := dl_bounds hc (g := g) hx
  have b2 := dl_bounds hc (g := g + x) hy
  have b3 := dl_bounds hc (g := g) hy
  have b4 := dl_bounds hc (g := g + y) hx
  have := Int.eq_zero_of_abs_lt_dvd hd (abs_lt.2 ⟨by omega, by omega⟩)
  omega

/-! ### Winding sums of walks -/

/-- The winding sum `W(g, L) = ∑_j δ(g + s₁ + ⋯ + s_{j-1}, s_j)` of the walk from `g` with
steps `L = [s₁, …, s_N]`. -/
def W (c : G → ZMod p) : G → List G → ℤ
  | _, [] => 0
  | g, s :: L => dl c g s + W c (g + s) L

@[simp] lemma W_nil (c : G → ZMod p) (g : G) : W c g [] = 0 := rfl

@[simp] lemma W_cons (c : G → ZMod p) (g s : G) (L : List G) :
    W c g (s :: L) = dl c g s + W c (g + s) L := rfl

lemma W_append (c : G → ZMod p) (g : G) (L₁ L₂ : List G) :
    W c g (L₁ ++ L₂) = W c g L₁ + W c (g + L₁.sum) L₂ := by
  induction L₁ generalizing g with
  | nil => simp
  | cons x L ih => simp [ih, add_assoc]

/-- The winding sum is invariant under permutations of the steps (squares do not wind). -/
lemma W_perm [NeZero p] {S : Set G} {c : G → ZMod p} (hc : IsCircHom S p q c) (h4 : p < 4 * q)
    {L L' : List G} (hperm : L.Perm L') : ∀ g, (∀ x ∈ L, x ∈ S) → W c g L = W c g L' := by
  induction hperm with
  | nil => intro g _; rfl
  | cons x _ ih =>
    intro g hL
    simp only [W_cons]
    rw [ih (g + x) (fun y hy => hL y (List.mem_cons_of_mem x hy))]
  | swap x y l =>
    intro g hL
    have hx : x ∈ S := hL x (by simp)
    have hy : y ∈ S := hL y (by simp)
    have h := dl_square hc h4 hy hx g
    simp only [W_cons, add_right_comm g y x]
    omega
  | trans h₁ _ ih₁ ih₂ =>
    intro g hL
    rw [ih₁ g hL, ih₂ g (fun y hy => hL y (h₁.symm.subset hy))]

/-- Telescoping in `ZMod p`. -/
lemma W_cast [NeZero p] (c : G → ZMod p) :
    ∀ (L : List G) (g : G), ((W c g L : ℤ) : ZMod p) = c (g + L.sum) - c g
  | [], g => by simp
  | s :: L, g => by
    have ih := W_cast c L (g + s)
    rw [W_cons, Int.cast_add, ih, dl_cast, List.sum_cons,
      show g + (s + L.sum) = g + s + L.sum from (add_assoc _ _ _).symm]
    ring

/-- A closed walk has winding sum divisible by `p`. -/
lemma W_closed [NeZero p] (c : G → ZMod p) {L : List G} (h0 : L.sum = 0) (g : G) :
    (p : ℤ) ∣ W c g L := by
  rw [← ZMod.intCast_zmod_eq_zero_iff_dvd, W_cast c L g, h0, add_zero, sub_self]

/-- Moving the base point of a closed walk by `s ∈ S` does not change the winding sum. -/
lemma W_shift [NeZero p] {S : Set G} (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod p}
    (hc : IsCircHom S p q c) (h4 : p < 4 * q) {L : List G} (hL : ∀ x ∈ L, x ∈ S)
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
  have h := W_perm hc h4 (List.perm_middle (l₁ := L) (a := s) (l₂ := [-s])) g hmem
  simp only [W_append, W_cons, W_nil, h0, add_zero, dl_neg hc h4 hs] at h
  linarith

/-- The winding sum of a closed walk only depends on the coset of the base point. -/
lemma W_shift_closure [NeZero p] {S : Set G} (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod p}
    (hc : IsCircHom S p q c) (h4 : p < 4 * q) {L : List G} (hL : ∀ x ∈ L, x ∈ S)
    (h0 : L.sum = 0) {γ : G} (hγ : γ ∈ AddSubgroup.closure S) :
    ∀ g, W c (g + γ) L = W c g L := by
  induction hγ using AddSubgroup.closure_induction with
  | mem x hx => exact W_shift hS hc h4 hL h0 hx
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

lemma sum_translate {M : Type*} [AddCommMonoid M] (φ : G → M) {a : G} (ha : a ∈ H) :
    ∑ γ : H, φ ((γ : G) + a) = ∑ γ : H, φ γ :=
  Fintype.sum_equiv (Equiv.addRight (⟨a, ha⟩ : H)) _ _ (fun γ => by simp)

/-- `F c s = ∑_{γ ∈ H} δ(γ, s)` (so `F c s / |H|` is the average of `δ(·, s)`). -/
def F (c : G → ZMod p) (s : G) : ℤ := ∑ γ : H, dl c γ s

lemma F_neg [NeZero p] {S : Set G} (hSH : S ⊆ H) {c : G → ZMod p} (hc : IsCircHom S p q c)
    (h4 : p < 4 * q) {s : G} (hs : s ∈ S) : F H c (-s) = Fintype.card H * p - F H c s := by
  unfold F
  rw [← sum_translate H (fun h => dl c h (-s)) (hSH hs)]
  simp only [dl_neg hc h4 hs, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul]

lemma F_bounds {S : Set G} {c : G → ZMod p} (hc : IsCircHom S p q c) {s : G} (hs : s ∈ S) :
    (Fintype.card H : ℤ) * q ≤ F H c s ∧ F H c s ≤ (Fintype.card H : ℤ) * (p - q) := by
  unfold F
  constructor
  · calc (Fintype.card H : ℤ) * q = ∑ _γ : H, (q : ℤ) := by simp
      _ ≤ ∑ γ : H, dl c γ s := Finset.sum_le_sum (fun γ _ => (dl_bounds hc hs).1)
  · calc ∑ γ : H, dl c γ s ≤ ∑ _γ : H, ((p : ℤ) - q) :=
          Finset.sum_le_sum (fun γ _ => (dl_bounds hc hs).2)
      _ = (Fintype.card H : ℤ) * (p - q) := by
          simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

/-- Summing the winding sums over all base points of a coset. -/
lemma sum_W {S : Set G} (hSH : S ⊆ H) (c : G → ZMod p) :
    ∀ (L : List G) (a : G), a ∈ H → (∀ x ∈ L, x ∈ S) →
      ∑ γ : H, W c ((γ : G) + a) L = (L.map (F H c)).sum
  | [], a, _, _ => by simp
  | s :: L, a, ha, hL => by
    have hs : s ∈ S := hL s (by simp)
    have ih := sum_W hSH c L (a + s) (H.add_mem ha (hSH hs))
      (fun x hx => hL x (List.mem_cons_of_mem s hx))
    simp only [W_cons, add_assoc, Finset.sum_add_distrib, List.map_cons, List.sum_cons, ih]
    congr 1
    exact sum_translate H (fun h => dl c h s) ha

/-- For a closed walk, `∑_j F(s_j) = |H| · W(0, L)`. -/
lemma sum_F_closed [NeZero p] {S : Set G} (hSH : S ⊆ H) (hHS : H ≤ AddSubgroup.closure S)
    (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod p} (hc : IsCircHom S p q c) (h4 : p < 4 * q)
    {L : List G} (hL : ∀ x ∈ L, x ∈ S) (h0 : L.sum = 0) :
    (L.map (F H c)).sum = Fintype.card H * W c 0 L := by
  rw [← sum_W H hSH c L 0 H.zero_mem hL]
  calc ∑ γ : H, W c ((γ : G) + 0) L = ∑ _γ : H, W c 0 L :=
        Finset.sum_congr rfl (fun γ _ => by
          simpa using W_shift_closure hS hc h4 hL h0 (hHS γ.2) 0)
    _ = Fintype.card H * W c 0 L := by simp

/-- `ξ₀(s) = F(s) / (p |H|)`, the average over `H` of `δ(·, s) / p`. -/
noncomputable def xi0 (c : G → ZMod p) (s : G) : ℝ := (F H c s : ℝ) / (p * Fintype.card H)

lemma xi0_eq (c : G → ZMod p) :
    xi0 H c = fun s => (F H c s : ℝ) / (p * Fintype.card H) := rfl

lemma card_pos_real : (0 : ℝ) < Fintype.card H := by
  exact_mod_cast Fintype.card_pos_iff.2 ⟨0⟩

lemma xi0_bounds [NeZero p] {S : Set G} {c : G → ZMod p} (hc : IsCircHom S p q c) {s : G}
    (hs : s ∈ S) : (q : ℝ) / p ≤ xi0 H c s ∧ xi0 H c s ≤ 1 - (q : ℝ) / p := by
  have hn := card_pos_real H
  have hp : (0 : ℝ) < p := by exact_mod_cast NeZero.pos p
  have hpn : (0 : ℝ) < p * Fintype.card H := mul_pos hp hn
  obtain ⟨h1, h2⟩ := F_bounds H hc hs
  have h1' : (Fintype.card H : ℝ) * q ≤ F H c s := by exact_mod_cast h1
  have h2' : (F H c s : ℝ) ≤ Fintype.card H * (p - q) := by exact_mod_cast h2
  unfold xi0
  constructor
  · rw [div_le_div_iff₀ hp hpn]; nlinarith
  · rw [one_sub_div hp.ne', div_le_div_iff₀ hpn hp]; nlinarith

lemma xi0_neg [NeZero p] {S : Set G} (hSH : S ⊆ H) {c : G → ZMod p} (hc : IsCircHom S p q c)
    (h4 : p < 4 * q) {s : G} (hs : s ∈ S) : xi0 H c (-s) = 1 - xi0 H c s := by
  have hn := card_pos_real H
  have hp : (0 : ℝ) < p := by exact_mod_cast NeZero.pos p
  unfold xi0
  rw [F_neg H hSH hc h4 hs, one_sub_div (mul_pos hp hn).ne']
  push_cast
  ring

lemma list_sum_div {α : Type*} (L : List α) (a : α → ℝ) (d : ℝ) :
    (L.map (fun x => a x / d)).sum = (L.map a).sum / d := by
  induction L with
  | nil => simp
  | cons x L ih => simp only [List.map_cons, List.sum_cons, ih, add_div]

lemma cast_list_sum_map {α : Type*} (L : List α) (f : α → ℤ) :
    (((L.map f).sum : ℤ) : ℝ) = (L.map (fun x => (f x : ℝ))).sum := by
  induction L with
  | nil => simp
  | cons x L ih => simp [ih]

/-- The key integrality: along a closed walk, `∑ ξ₀ = W(0, L) / p ∈ ℤ`. -/
lemma sum_xi0_closed [NeZero p] {S : Set G} (hSH : S ⊆ H) (hHS : H ≤ AddSubgroup.closure S)
    (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod p} (hc : IsCircHom S p q c) (h4 : p < 4 * q)
    {L : List G} (hL : ∀ x ∈ L, x ∈ S) (h0 : L.sum = 0) :
    ∃ m : ℤ, (L.map (xi0 H c)).sum = m := by
  obtain ⟨m, hm⟩ := W_closed c h0 0
  refine ⟨m, ?_⟩
  have hn := (card_pos_real H).ne'
  have hp : (p : ℝ) ≠ 0 := by exact_mod_cast NeZero.ne p
  have hF : (L.map (fun s => (F H c s : ℝ))).sum = Fintype.card H * (p * (m : ℝ)) := by
    rw [← cast_list_sum_map, sum_F_closed H hSH hHS hS hc h4 hL h0, hm]
    push_cast
    ring
  rw [xi0_eq, list_sum_div L (fun s => (F H c s : ℝ)), hF, div_eq_iff (mul_ne_zero hp hn)]
  ring

lemma sum_xi0_neg [NeZero p] {S : Set G} (hSH : S ⊆ H) {c : G → ZMod p}
    (hc : IsCircHom S p q c) (h4 : p < 4 * q) :
    ∀ (L : List G), (∀ x ∈ L, x ∈ S) →
      ((L.map (fun z => -z)).map (xi0 H c)).sum = L.length - (L.map (xi0 H c)).sum
  | [], _ => by simp
  | s :: L, hL => by
    have hs : s ∈ S := hL s (by simp)
    have ih := sum_xi0_neg hSH hc h4 L (fun x hx => hL x (List.mem_cons_of_mem s hx))
    simp only [List.map_cons, List.sum_cons, List.length_cons, ih, xi0_neg H hSH hc h4 hs]
    push_cast
    ring

/-- The value in `ℝ/ℤ` of a list of steps. -/
noncomputable def V (c : G → ZMod p) (L : List G) : AddCircle (1 : ℝ) :=
  (((L.map (xi0 H c)).sum : ℝ) : AddCircle (1 : ℝ))

lemma V_append (c : G → ZMod p) (L₁ L₂ : List G) :
    V H c (L₁ ++ L₂) = V H c L₁ + V H c L₂ := by
  simp only [V, List.map_append, List.sum_append, AddCircle.coe_add]

/-- Well-definedness: lists with the same sum have the same value in `ℝ/ℤ`. -/
lemma V_eq [NeZero p] {S : Set G} (hSH : S ⊆ H) (hHS : H ≤ AddSubgroup.closure S)
    (hS : ∀ s ∈ S, -s ∈ S) {c : G → ZMod p} (hc : IsCircHom S p q c) (h4 : p < 4 * q)
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
  obtain ⟨m, hm⟩ := sum_xi0_closed H hSH hHS hS hc h4 hL h0
  rw [List.map_append, List.sum_append, sum_xi0_neg H hSH hc h4 L₂ h₂] at hm
  unfold V
  rw [← sub_eq_zero, ← AddCircle.coe_sub, AddCircle.coe_eq_zero_iff]
  refine ⟨m - L₂.length, ?_⟩
  rw [zsmul_eq_mul, mul_one]
  push_cast
  linarith

end average

/-! ### Main theorems -/

/-- **Theorem W for `K_{p/q}`, general form.** Let `G` be a finite abelian group, `S ⊆ G`
symmetric, `0 < p < 4q`, and `c : G → ZMod p` a homomorphism `Cay(G, S) → K_{p/q}`. Then there
is a character `ξ : closure S →+ ℝ/ℤ` such that every `ξ s` (`s ∈ S`) is represented by a real
number in `[q/p, 1 - q/p]`. -/
theorem theoremWplus_general {G : Type*} [AddCommGroup G] [Finite G]
    (S : Set G) (hS : ∀ s ∈ S, -s ∈ S)
    (p q : ℕ) (hp0 : 0 < p) (h4 : p < 4 * q)
    (c : G → ZMod p)
    (hc : ∀ g, ∀ s ∈ S, q ≤ (c (g + s) - c g).val ∧ (c (g + s) - c g).val ≤ p - q) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, (q : ℝ) / p ≤ x ∧ x ≤ 1 - (q : ℝ) / p ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  classical
  have : NeZero p := ⟨hp0.ne'⟩
  set H := AddSubgroup.closure S
  let _ : Fintype H := Fintype.ofFinite H
  have hSH : S ⊆ H := AddSubgroup.subset_closure
  have hHS : H ≤ AddSubgroup.closure S := le_rfl
  have hlist : ∀ γ : H, ∃ L : List G, (∀ x ∈ L, x ∈ S) ∧ L.sum = γ :=
    fun γ => exists_list_of_mem_closure hS γ.2
  choose Lf hLf hLs using hlist
  let ξf : H → AddCircle (1 : ℝ) := fun γ => V H c (Lf γ)
  have key : ∀ (L : List G), (∀ x ∈ L, x ∈ S) → ∀ γ : H, L.sum = γ → ξf γ = V H c L :=
    fun L hL γ he => V_eq H hSH hHS hS hc h4 (hLf γ) hL ((hLs γ).trans he.symm)
  refine ⟨AddMonoidHom.mk' ξf (fun a b => ?_), fun s hs =>
    ⟨xi0 H c s, (xi0_bounds H hc hs).1, (xi0_bounds H hc hs).2, ?_⟩⟩
  · rw [key (Lf a ++ Lf b) (fun x hx => (List.mem_append.1 hx).elim (hLf a x) (hLf b x))
      (a + b) (by rw [List.sum_append, hLs a, hLs b, AddSubgroup.coe_add]), V_append]
  · show _ = ξf _
    rw [key [s] (by simpa using hs) _ (by simp)]
    simp [V]

-- `h2 : 2 * q < p` is part of the requested statement but is not needed by the proof.
set_option linter.unusedVariables false in
/-- **Theorem W for circular cliques `K_{p/q}`, `p` odd, `2q < p < 4q`** (finite groups).
Let `G` be a finite abelian group, `S ⊆ G` symmetric and `c : G → ZMod p` a homomorphism
`Cay(G, S) → K_{p/q}`. Then there is a character `ξ : closure S →+ ℝ/ℤ` such that every `ξ s`
(`s ∈ S`) is represented by a real number in `[q/p, 1 - q/p]`. The hypothesis `2q < p` is not
needed, and of `Odd p` only `0 < p` is used (see `theoremWplus_general`). -/
theorem theoremWplus_finite {G : Type*} [AddCommGroup G] [Finite G]
    (S : Set G) (hS : ∀ s ∈ S, -s ∈ S)
    (p q : ℕ) (hp : Odd p) (h2 : 2 * q < p) (h4 : p < 4 * q)
    (c : G → ZMod p)
    (hc : ∀ g, ∀ s ∈ S, q ≤ (c (g + s) - c g).val ∧ (c (g + s) - c g).val ≤ p - q) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, (q : ℝ) / p ≤ x ∧ x ≤ 1 - (q : ℝ) / p ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ :=
  theoremWplus_general S hS p q hp.pos h4 c hc

/-- Sanity check: the hypotheses are satisfiable (`G = ZMod 7`, `S = {2, 3, 4, 5}`, `p = 7`,
`q = 2`, `c = id`, i.e. `Cay(ℤ/7, {±2, ±3}) = K_{7/2}`). -/
example : ∃ ξ : AddSubgroup.closure ({2, 3, 4, 5} : Set (ZMod 7)) →+ AddCircle (1 : ℝ),
    ∀ s (hs : s ∈ ({2, 3, 4, 5} : Set (ZMod 7))), ∃ x : ℝ, 2 / 7 ≤ x ∧ x ≤ 5 / 7 ∧
      (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  have hS : ∀ s ∈ ({2, 3, 4, 5} : Set (ZMod 7)), -s ∈ ({2, 3, 4, 5} : Set (ZMod 7)) := by
    intro s hs
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hs ⊢
    rcases hs with rfl | rfl | rfl | rfl <;> decide
  have hc : ∀ g, ∀ s ∈ ({2, 3, 4, 5} : Set (ZMod 7)),
      2 ≤ ((id : ZMod 7 → ZMod 7) (g + s) - id g).val ∧
        ((id : ZMod 7 → ZMod 7) (g + s) - id g).val ≤ 7 - 2 := by
    intro g s hs
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hs
    simp only [id, add_sub_cancel_left]
    rcases hs with rfl | rfl | rfl | rfl <;> decide
  obtain ⟨ξ, hξ⟩ := theoremWplus_finite _ hS 7 2 ⟨3, rfl⟩ (by norm_num) (by norm_num) id hc
  refine ⟨ξ, fun s hs => ?_⟩
  obtain ⟨x, h1, h2, h3⟩ := hξ s hs
  refine ⟨x, ?_, ?_, h3⟩
  · norm_num at h1; linarith
  · norm_num at h2; linarith

/-- Sanity check: `p = 3`, `q = 1` is allowed (`2 < 3 < 4`), and `K_{3/1} = K_3`; this recovers
Theorem W of `TheoremW.lean` (`TheoremW.theoremW_finite`, same statement) for proper 3-colourings. -/
example {G : Type*} [AddCommGroup G] [Finite G]
    (S : Set G) (hS : ∀ s ∈ S, -s ∈ S)
    (c : G → ZMod 3) (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, 1/3 ≤ x ∧ x ≤ 2/3 ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  have hc' : ∀ g, ∀ s ∈ S, 1 ≤ (c (g + s) - c g).val ∧ (c (g + s) - c g).val ≤ 3 - 1 := by
    intro g s hs
    have h1 : (c (g + s) - c g).val ≠ 0 := by
      rw [Ne, ZMod.val_eq_zero, sub_eq_zero]
      exact hc g s hs
    have h2 := ZMod.val_lt (c (g + s) - c g)
    constructor <;> omega
  obtain ⟨ξ, hξ⟩ := theoremWplus_finite S hS 3 1 ⟨1, rfl⟩ (by norm_num) (by norm_num) c hc'
  refine ⟨ξ, fun s hs => ?_⟩
  obtain ⟨x, h1, h2, h3⟩ := hξ s hs
  refine ⟨x, ?_, ?_, h3⟩
  · norm_num at h1; linarith
  · norm_num at h2; linarith

/-! ### The converse: colourings from characters -/

/-- **Converse, for any lift.** Let `ξ : G →+ ℝ/ℤ` be a character of `G` and `r : G → ℝ` any
lift of `ξ`. If every `ξ s` (`s ∈ S`) is represented by a real number in `[q/p, 1 - q/p]`, then
`g ↦ ⌊p · r g⌋ mod p` is a homomorphism `Cay(G, S) → K_{p/q}`. -/
theorem converse_of_lift {G : Type*} [AddCommGroup G] (S : Set G) (p q : ℕ) (hp0 : 0 < p)
    (ξ : G →+ AddCircle (1 : ℝ)) (r : G → ℝ) (hr : ∀ g, (r g : AddCircle (1 : ℝ)) = ξ g)
    (hξ : ∀ s ∈ S, ∃ x : ℝ, (q : ℝ) / p ≤ x ∧ x ≤ 1 - (q : ℝ) / p ∧
      (x : AddCircle (1 : ℝ)) = ξ s) :
    ∀ g, ∀ s ∈ S,
      q ≤ ((⌊(p : ℝ) * r (g + s)⌋ : ZMod p) - (⌊(p : ℝ) * r g⌋ : ZMod p)).val ∧
      ((⌊(p : ℝ) * r (g + s)⌋ : ZMod p) - (⌊(p : ℝ) * r g⌋ : ZMod p)).val ≤ p - q := by
  intro g s hs
  have : NeZero p := ⟨hp0.ne'⟩
  obtain ⟨x, hx1, hx2, hx⟩ := hξ s hs
  have hp : (0 : ℝ) < p := by exact_mod_cast hp0
  -- `r (g + s) = r g + x + k` for an integer `k`
  have h0 : ((r (g + s) - r g - x : ℝ) : AddCircle (1 : ℝ)) = 0 := by
    rw [AddCircle.coe_sub, AddCircle.coe_sub, hr, hr, hx, map_add]
    abel
  obtain ⟨k, hk⟩ := (AddCircle.coe_eq_zero_iff (1 : ℝ)).1 h0
  rw [zsmul_eq_mul, mul_one] at hk
  -- `p x ∈ [q, p - q]`
  have hqp : (p : ℝ) * ((q : ℝ) / p) = q := by field_simp
  have hb1 : (q : ℝ) ≤ p * x := by nlinarith
  have hb2 : (p : ℝ) * x ≤ p - q := by nlinarith
  have hrs : (p : ℝ) * r (g + s) = (p : ℝ) * r g + p * x + ((p * k : ℤ) : ℝ) := by
    push_cast
    rw [hk]
    ring
  -- the colour difference is `D = ⌊p r g + p x⌋ - ⌊p r g⌋ ∈ [q, p - q]`
  have hdiff : ((⌊(p : ℝ) * r (g + s)⌋ : ZMod p) - (⌊(p : ℝ) * r g⌋ : ZMod p))
      = ((⌊(p : ℝ) * r g + p * x⌋ - ⌊(p : ℝ) * r g⌋ : ℤ) : ZMod p) := by
    rw [hrs, Int.floor_add_intCast]
    push_cast
    rw [ZMod.natCast_self]
    ring
  have hD1 : ⌊(p : ℝ) * r g⌋ + (q : ℤ) ≤ ⌊(p : ℝ) * r g + p * x⌋ := by
    rw [Int.le_floor]
    push_cast
    have := Int.floor_le ((p : ℝ) * r g)
    linarith
  have hD2 : ⌊(p : ℝ) * r g + p * x⌋ < ⌊(p : ℝ) * r g⌋ + ((p : ℤ) - q) + 1 := by
    rw [Int.floor_lt]
    push_cast
    have := Int.lt_floor_add_one ((p : ℝ) * r g)
    linarith
  rw [hdiff]
  have hv : (((⌊(p : ℝ) * r g + p * x⌋ - ⌊(p : ℝ) * r g⌋ : ℤ) : ZMod p).val : ℤ)
      = (⌊(p : ℝ) * r g + p * x⌋ - ⌊(p : ℝ) * r g⌋) % p := ZMod.val_intCast _
  have hlt := ZMod.val_lt (((⌊(p : ℝ) * r g + p * x⌋ - ⌊(p : ℝ) * r g⌋ : ℤ) : ZMod p))
  rcases Nat.eq_zero_or_pos q with hq | hq
  · subst hq
    constructor <;> omega
  · rw [Int.emod_eq_of_lt (by omega) (by omega)] at hv
    constructor <;> omega

/-- **Converse.** Let `ξ : G →+ ℝ/ℤ` be a character of `G` such that every `ξ s` (`s ∈ S`) is
represented by a real number in `[q/p, 1 - q/p]`. Then `c g = ⌊p · ρ(ξ g)⌋ mod p`, where
`ρ(ξ g) ∈ [0, 1)` is the representative of `ξ g`, is a homomorphism `Cay(G, S) → K_{p/q}`,
i.e. it satisfies the hypothesis `hc` of `theoremWplus_finite`. -/
theorem converse {G : Type*} [AddCommGroup G] (S : Set G) (p q : ℕ) (hp0 : 0 < p)
    (ξ : G →+ AddCircle (1 : ℝ))
    (hξ : ∀ s ∈ S, ∃ x : ℝ, (q : ℝ) / p ≤ x ∧ x ≤ 1 - (q : ℝ) / p ∧
      (x : AddCircle (1 : ℝ)) = ξ s)
    (c : G → ZMod p)
    (hcdef : ∀ g, c g = ⌊(p : ℝ) * (AddCircle.equivIco (1 : ℝ) 0 (ξ g) : ℝ)⌋) :
    ∀ g, ∀ s ∈ S, q ≤ (c (g + s) - c g).val ∧ (c (g + s) - c g).val ≤ p - q := by
  intro g s hs
  simp only [hcdef]
  exact converse_of_lift S p q hp0 ξ (fun g => (AddCircle.equivIco (1 : ℝ) 0 (ξ g) : ℝ))
    (fun _ => AddCircle.coe_equivIco) hξ g s hs

end TheoremWplus
