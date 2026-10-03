import TheoremW
import TheoremWplus
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Topology.Order.Compact
import Mathlib.Order.Filter.Cofinite

/-!
# Theorem W, case `k = 1`, for every abelian group and finite `S`

Let `G` be an abelian group (not necessarily finite), `S ⊆ G` a finite symmetric set and
`c : G → ZMod 3` a proper colouring of the Cayley graph `Cay(G, S)`. Then there is an additive
character `ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ)` with `ξ s ∈ [1/3, 2/3] (mod 1)` for
`s ∈ S` (`theoremW`). More generally, a homomorphism `Cay(G, S) → K_{p/q}` with `p < 4q` gives
a character with `ξ s ∈ [q/p, 1 - q/p] (mod 1)` (`theoremWplus`, the infinite analogue of
`TheoremWplus.theoremWplus_general`; `theoremW` is its case `p = 3`, `q = 1`, proved separately
from `TheoremW.sig` and recovered from `theoremWplus` in an example).

The proof reuses the winding sums of `TheoremW.lean` (`TheoremW.sig`, `TheoremW.W` and their
invariance lemmas, none of which use finiteness of `G`) and replaces the average over the finite
group `closure S` by Følner averaging and a compactness argument.

The averaging and compactness steps (`walk_avg`, `closed_avg`, `exists_F`) are stated for any
bounded real step function `d` and winding sum `Wd`, and are applied to `TheoremW.sig` (bounds
`[-1, 1]`) and to `TheoremWplus.dl` (bounds `[q, p - q]` on `S`).

* Averaging operators (`avg`): `avg n s f g = (1/(n+1)) ∑_{t ≤ n} f (g + t • s)` on functions
  `G → ℝ`. They do not increase a bound `|f| ≤ M` (`avg_bound`), commute with translations
  (`avg_shift`), and `|avg n s f (g + s) - avg n s f g| ≤ 2M/(n+1)` (`avg_step`, telescoping).
* `avgs n l` is the composition of the `avg n s` for `s` in a list `l` enumerating `S`
  (a box average). For every `s ∈ l` it is `2M/(n+1)`-almost invariant under translation by `s`
  (`avgs_lip`), and it fixes functions that are constant on `g + closure S` (`avgs_const`).
* For a closed walk `L` (`L.sum = 0`) the winding sum `W(·, L)` is constant on cosets of
  `closure S` (`TheoremW.W_shift_closure`), so with `a n s = avgs n l (σ(·, s)) 0`,
  `|W(0, L) - ∑_{s ∈ L} a n s| ≤ 2 |L|² / (n+1)` (`walk_avg`, `closed_avg`).
* Compactness: `a n s ∈ [-1, 1]`; we take limits along the ultrafilter `hyperfilter ℕ`
  (finer than `atTop`), simultaneously for all `s`, instead of extracting a subsequence.
  The limit `F` satisfies `∑_{s ∈ L} F s = W(0, L)` for every closed walk `L` (`exists_F`).
* `ξ₀ s = (F s + 3)/6 ∈ [1/3, 2/3]` has integer sums along closed walks, and so induces a
  character of `closure S` (`exists_char`, the construction of `TheoremW.theoremW_finite`
  generalised to any function with integer sums along closed walks).

Decisions: limits are taken along an ultrafilter rather than a convergent subsequence
(`IsCompact.tendsto_subseq`), since this gives all coordinates at once and needs no product
topology; the averaging operators use `n + 1` points so that no case `n = 0` arises.
-/

set_option autoImplicit false

namespace TheoremWInf

open Finset Filter Topology

variable {G : Type*} [AddCommGroup G]

/-! ### Averaging operators -/

/-- `avg n s f g = (1/(n+1)) ∑_{t ≤ n} f (g + t • s)`. -/
noncomputable def avg (n : ℕ) (s : G) (f : G → ℝ) (g : G) : ℝ :=
  (∑ t ∈ range (n + 1), f (g + t • s)) / (n + 1)

lemma avg_sub (n : ℕ) (s : G) (f f' : G → ℝ) (g : G) :
    avg n s (fun x => f x - f' x) g = avg n s f g - avg n s f' g := by
  simp only [avg, sum_sub_distrib, sub_div]

lemma avg_add (n : ℕ) (s : G) (f f' : G → ℝ) (g : G) :
    avg n s (fun x => f x + f' x) g = avg n s f g + avg n s f' g := by
  simp only [avg, sum_add_distrib, add_div]

lemma avg_shift (n : ℕ) (s : G) (f : G → ℝ) (h g : G) :
    avg n s (fun x => f (x + h)) g = avg n s f (g + h) := by
  unfold avg
  congr 1
  refine sum_congr rfl (fun t _ => ?_)
  rw [add_right_comm]

lemma avg_bound {n : ℕ} {s : G} {f : G → ℝ} {M : ℝ} (hf : ∀ x, |f x| ≤ M) (g : G) :
    |avg n s f g| ≤ M := by
  unfold avg
  have hn : (0 : ℝ) < n + 1 := by positivity
  rw [abs_div, abs_of_pos hn, div_le_iff₀ hn]
  calc |∑ t ∈ range (n + 1), f (g + t • s)| ≤ ∑ t ∈ range (n + 1), |f (g + t • s)| :=
        abs_sum_le_sum_abs _ _
    _ ≤ ∑ _t ∈ range (n + 1), M := sum_le_sum (fun t _ => hf _)
    _ = M * (n + 1) := by
        rw [sum_const, card_range, nsmul_eq_mul]
        push_cast
        ring

/-- Shifting by the averaging direction telescopes. -/
lemma avg_step {n : ℕ} {s : G} {f : G → ℝ} {M : ℝ} (hf : ∀ x, |f x| ≤ M) (g : G) :
    |avg n s f (g + s) - avg n s f g| ≤ 2 * M / (n + 1) := by
  have hn : (0 : ℝ) < n + 1 := by positivity
  have htel : ∑ t ∈ range (n + 1), f (g + s + t • s) - ∑ t ∈ range (n + 1), f (g + t • s)
      = f (g + (n + 1) • s) - f g := by
    rw [← sum_sub_distrib]
    have h := Finset.sum_range_sub (fun t : ℕ => f (g + t • s)) (n + 1)
    simp only [zero_smul, add_zero] at h
    rw [← h]
    refine sum_congr rfl (fun t _ => ?_)
    rw [succ_nsmul', ← add_assoc]
  unfold avg
  rw [← sub_div, htel, abs_div, abs_of_pos hn, div_le_div_iff_of_pos_right hn]
  have h1 := abs_le.1 (hf (g + (n + 1) • s))
  have h2 := abs_le.1 (hf g)
  rw [abs_le]
  constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]

/-- Averaging preserves almost-invariance under a translation. -/
lemma avg_lip {n : ℕ} {s h : G} {f : G → ℝ} {ε : ℝ} (hf : ∀ x, |f (x + h) - f x| ≤ ε)
    (g : G) : |avg n s f (g + h) - avg n s f g| ≤ ε := by
  rw [← avg_shift, ← avg_sub]
  exact avg_bound hf g

/-- The box average: the composition of the `avg n s` for `s ∈ l`. -/
noncomputable def avgs (n : ℕ) : List G → (G → ℝ) → G → ℝ
  | [], f => f
  | s :: l, f => avg n s (avgs n l f)

lemma avgs_bound {n : ℕ} {f : G → ℝ} {M : ℝ} (hf : ∀ x, |f x| ≤ M) :
    ∀ (l : List G) (g : G), |avgs n l f g| ≤ M
  | [], g => hf g
  | _ :: l, g => avg_bound (avgs_bound hf l) g

lemma avgs_add (n : ℕ) (f f' : G → ℝ) :
    ∀ l : List G, avgs n l (fun x => f x + f' x) = fun g => avgs n l f g + avgs n l f' g
  | [] => rfl
  | s :: l => by
    funext g
    simp only [avgs]
    rw [avgs_add n f f' l, avg_add]

lemma avgs_shift (n : ℕ) (f : G → ℝ) (h : G) :
    ∀ l : List G, avgs n l (fun x => f (x + h)) = fun g => avgs n l f (g + h)
  | [] => rfl
  | s :: l => by
    funext g
    simp only [avgs]
    rw [avgs_shift n f h l, avg_shift]

/-- The box average is almost invariant under translation by each direction of the box. -/
lemma avgs_lip {n : ℕ} {f : G → ℝ} {M : ℝ} (hf : ∀ x, |f x| ≤ M) {s : G} :
    ∀ (l : List G), s ∈ l → ∀ g, |avgs n l f (g + s) - avgs n l f g| ≤ 2 * M / (n + 1)
  | [], hs, _ => absurd hs (List.not_mem_nil)
  | t :: l, hs, g => by
    simp only [avgs]
    rcases List.mem_cons.1 hs with rfl | hs
    · exact avg_step (avgs_bound hf l) g
    · exact avg_lip (avgs_lip hf l hs) g

/-- The box average fixes functions that are constant on a coset of a subgroup containing the
directions. -/
lemma avgs_const {H : AddSubgroup G} {n : ℕ} {g : G} {f : G → ℝ} {C : ℝ}
    (hf : ∀ h ∈ H, f (g + h) = C) : ∀ (l : List G), (∀ s ∈ l, s ∈ H) →
      ∀ h ∈ H, avgs n l f (g + h) = C
  | [], _ => hf
  | t :: l, hl => fun h hh => by
    have ht : t ∈ H := hl t (by simp)
    have ih := avgs_const (n := n) hf l (fun s hs => hl s (List.mem_cons_of_mem t hs))
    have hn : (n : ℝ) + 1 ≠ 0 := by positivity
    have hk : ∀ k ∈ range (n + 1), avgs n l f (g + h + k • t) = C := fun k _ => by
      rw [add_assoc]
      exact ih _ (H.add_mem hh (H.nsmul_mem ht k))
    simp only [avgs, avg]
    rw [sum_congr rfl hk, sum_const, card_range, nsmul_eq_mul]
    push_cast
    field_simp

/-! ### Averaged winding sums

The walk part is stated for any real step function `d : G → G → ℝ` (the lift of a colour
difference along an edge) and any `Wd : G → List G → ℝ` satisfying the recursion of a winding
sum; it is used for `TheoremW.sig` (three colours) and `TheoremWplus.dl` (circular cliques). -/

lemma avg_mono {n : ℕ} {s : G} {f : G → ℝ} {a b : ℝ} (hf : ∀ x, a ≤ f x ∧ f x ≤ b) (g : G) :
    a ≤ avg n s f g ∧ avg n s f g ≤ b := by
  unfold avg
  have hn : (0 : ℝ) < n + 1 := by positivity
  have hc : ∀ e : ℝ, ∑ _t ∈ range (n + 1), e = e * (n + 1) := fun e => by
    rw [sum_const, card_range, nsmul_eq_mul]
    push_cast
    ring
  constructor
  · rw [le_div_iff₀ hn, ← hc]
    exact sum_le_sum (fun t _ => (hf _).1)
  · rw [div_le_iff₀ hn, ← hc]
    exact sum_le_sum (fun t _ => (hf _).2)

lemma avgs_mono {n : ℕ} {f : G → ℝ} {a b : ℝ} (hf : ∀ x, a ≤ f x ∧ f x ≤ b) :
    ∀ (l : List G) (g : G), a ≤ avgs n l f g ∧ avgs n l f g ≤ b
  | [], g => hf g
  | _ :: l, g => avg_mono (avgs_mono hf l) g

lemma list_sum_sub_le {α : Type*} (u v : α → ℝ) (ε : ℝ) :
    ∀ (L : List α), (∀ t ∈ L, |u t - v t| ≤ ε) →
      |(L.map u).sum - (L.map v).sum| ≤ L.length * ε
  | [], _ => by simp
  | x :: L, h => by
    have h1 := abs_le.1 (h x (by simp))
    have h2 := abs_le.1 (list_sum_sub_le u v ε L (fun t ht => h t (List.mem_cons_of_mem x ht)))
    simp only [List.map_cons, List.sum_cons, List.length_cons]
    push_cast
    rw [abs_le]
    constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]

section walks

variable (d : G → G → ℝ) (Wd : G → List G → ℝ) (hnil : ∀ x, Wd x [] = 0)
  (hcons : ∀ x s L, Wd x (s :: L) = d x s + Wd (x + s) L) {M : ℝ} (hM : ∀ x s, |d x s| ≤ M)
include hnil hcons hM

/-- The averaged winding sum of a walk is close to the sum of the averaged steps. -/
lemma walk_avg (n : ℕ) (l : List G) :
    ∀ (L : List G), (∀ x ∈ L, x ∈ l) → ∀ g,
      |avgs n l (fun x => Wd x L) g - (L.map (fun s => avgs n l (fun x => d x s) g)).sum|
        ≤ (L.length : ℝ) ^ 2 * (2 * M / (n + 1))
  | [], _, g => by
    have h := avgs_bound (n := n) (f := fun x => Wd x []) (M := 0)
      (fun x => by simp [hnil]) l g
    simpa using h
  | s :: L, hL, g => by
    have hs : s ∈ l := hL s (by simp)
    have ih := abs_le.1 (walk_avg n l L (fun x hx => hL x (List.mem_cons_of_mem s hx)) (g + s))
    have hstep := abs_le.1 (list_sum_sub_le
      (fun t => avgs n l (fun x => d x t) (g + s))
      (fun t => avgs n l (fun x => d x t) g) (2 * M / (n + 1)) L
      (fun t _ => avgs_lip (n := n) (fun x => hM x t) l hs g))
    have e1 : (fun x => Wd x (s :: L)) = fun x => d x s + Wd (x + s) L := by
      funext x
      exact hcons x s L
    have e2 := avgs_shift n (fun x => Wd x L) s l
    rw [e1, avgs_add, e2]
    simp only [List.map_cons, List.sum_cons, List.length_cons]
    have hM0 : 0 ≤ M := le_trans (abs_nonneg _) (hM 0 0)
    have hε : (0 : ℝ) ≤ 2 * M / (n + 1) := by positivity
    have hLl : (0 : ℝ) ≤ L.length := by positivity
    push_cast
    rw [abs_le]
    constructor <;> nlinarith [ih.1, ih.2, hstep.1, hstep.2]

/-- For a walk whose winding sum is constant on `closure S`,
`|Wd 0 L - ∑_{s ∈ L} a n s| ≤ 2 M |L|² / (n + 1)`. -/
lemma closed_avg {S : Set G} {l : List G} (hl : ∀ x, x ∈ l ↔ x ∈ S)
    {L : List G} (hL : ∀ x ∈ L, x ∈ S)
    (hconst : ∀ h ∈ AddSubgroup.closure S, Wd h L = Wd 0 L) (n : ℕ) :
    |Wd 0 L - (L.map (fun s => avgs n l (fun x => d x s) 0)).sum|
      ≤ (L.length : ℝ) ^ 2 * (2 * M / (n + 1)) := by
  have h1 := avgs_const (H := AddSubgroup.closure S) (n := n) (g := 0)
    (f := fun x => Wd x L) (C := Wd 0 L)
    (fun h hh => by rw [zero_add]; exact hconst h hh) l
    (fun s hs => AddSubgroup.subset_closure ((hl s).1 hs)) 0 (AddSubgroup.zero_mem _)
  rw [add_zero] at h1
  have h := walk_avg d Wd hnil hcons hM n l L (fun x hx => (hl x).2 (hL x hx)) 0
  rwa [h1] at h

/-! ### Compactness -/

/-- The limit of the averaged steps: `F s ∈ [a, b]` for `s ∈ S` (when the steps along `s` lie in
`[a, b]`) and `∑_{s ∈ L} F s = Wd 0 L` for every walk `L` with steps in `S` whose winding sum is
constant on `closure S`. -/
lemma exists_F {S : Set G} (hSfin : S.Finite) {a b : ℝ}
    (hab : ∀ x, ∀ s ∈ S, a ≤ d x s ∧ d x s ≤ b) :
    ∃ F : G → ℝ, (∀ s ∈ S, a ≤ F s ∧ F s ≤ b) ∧ ∀ L : List G, (∀ x ∈ L, x ∈ S) →
      (∀ h ∈ AddSubgroup.closure S, Wd h L = Wd 0 L) → (L.map F).sum = Wd 0 L := by
  classical
  set l := hSfin.toFinset.toList
  have hl : ∀ x, x ∈ l ↔ x ∈ S := fun x => by simp [l]
  set A : ℕ → G → ℝ := fun n s => avgs n l (fun x => d x s) 0 with hA_def
  have hA : ∀ n s, |A n s| ≤ M := fun n s => avgs_bound (fun x => hM x s) l 0
  set U : Ultrafilter ℕ := hyperfilter ℕ
  have hlim : ∀ s, ∃ y, Tendsto (fun n => A n s) U (𝓝 y) := by
    intro s
    obtain ⟨y, -, hle⟩ := (isCompact_Icc (a := -M) (b := M)).ultrafilter_le_nhds
      (U.map (fun n => A n s)) (by
        rw [Ultrafilter.coe_map, le_principal_iff, Filter.mem_map]
        exact univ_mem' (fun n => abs_le.1 (hA n s)))
    exact ⟨y, hle⟩
  choose F hF using hlim
  refine ⟨F, fun s hs => ⟨ge_of_tendsto' (hF s) (fun n => (avgs_mono (hab · s hs) l 0).1),
    le_of_tendsto' (hF s) (fun n => (avgs_mono (hab · s hs) l 0).2)⟩, fun L hL hconst => ?_⟩
  have hUle : (U : Filter ℕ) ≤ atTop :=
    hyperfilter_le_cofinite.trans Nat.cofinite_eq_atTop.le
  have t1 : Tendsto (fun n => (L.map (A n)).sum) U (𝓝 (L.map F).sum) :=
    tendsto_list_sum L (fun s _ => hF s)
  have e : Tendsto (fun n : ℕ => (L.length : ℝ) ^ 2 * (2 * M / ((n : ℝ) + 1))) atTop (𝓝 0) := by
    have := tendsto_one_div_add_atTop_nhds_zero_nat.const_mul ((L.length : ℝ) ^ 2 * (2 * M))
    rw [mul_zero] at this
    exact this.congr (fun n => by ring)
  have lo : Tendsto (fun n : ℕ => Wd 0 L - (L.length : ℝ) ^ 2 * (2 * M / ((n : ℝ) + 1)))
      U (𝓝 (Wd 0 L)) := by
    simpa using (tendsto_const_nhds.sub e).mono_left hUle
  have hi : Tendsto (fun n : ℕ => Wd 0 L + (L.length : ℝ) ^ 2 * (2 * M / ((n : ℝ) + 1)))
      U (𝓝 (Wd 0 L)) := by
    simpa using (tendsto_const_nhds.add e).mono_left hUle
  have key := fun n => abs_le.1 (closed_avg d Wd hnil hcons hM hl hL hconst n)
  have t2 : Tendsto (fun n => (L.map (A n)).sum) U (𝓝 (Wd 0 L)) :=
    tendsto_of_tendsto_of_tendsto_of_le_of_le lo hi (fun n => by linarith [(key n).2])
      (fun n => by linarith [(key n).1])
  exact tendsto_nhds_unique t1 t2

end walks

/-! ### The character -/

section character

variable {S : Set G} (x0 : G → ℝ)

/-- The value in `ℝ/ℤ` of a list of steps. -/
noncomputable def V (L : List G) : AddCircle (1 : ℝ) :=
  (((L.map x0).sum : ℝ) : AddCircle (1 : ℝ))

omit [AddCommGroup G] in
lemma V_append (L₁ L₂ : List G) : V x0 (L₁ ++ L₂) = V x0 L₁ + V x0 L₂ := by
  simp only [V, List.map_append, List.sum_append, AddCircle.coe_add]

lemma sum_x0_neg (hneg : ∀ s ∈ S, x0 (-s) = 1 - x0 s) :
    ∀ (L : List G), (∀ x ∈ L, x ∈ S) →
      ((L.map (fun z => -z)).map x0).sum = L.length - (L.map x0).sum
  | [], _ => by simp
  | s :: L, hL => by
    have hs : s ∈ S := hL s (by simp)
    have ih := sum_x0_neg hneg L (fun x hx => hL x (List.mem_cons_of_mem s hx))
    simp only [List.map_cons, List.sum_cons, List.length_cons, ih, hneg s hs]
    push_cast
    ring

/-- Well-definedness: lists with the same sum have the same value in `ℝ/ℤ`. -/
lemma V_eq (hS : ∀ s ∈ S, -s ∈ S) (hneg : ∀ s ∈ S, x0 (-s) = 1 - x0 s)
    (hcl : ∀ L : List G, (∀ x ∈ L, x ∈ S) → L.sum = 0 → ∃ m : ℤ, (L.map x0).sum = m)
    {L₁ L₂ : List G} (h₁ : ∀ x ∈ L₁, x ∈ S) (h₂ : ∀ x ∈ L₂, x ∈ S)
    (he : L₁.sum = L₂.sum) : V x0 L₁ = V x0 L₂ := by
  have hL : ∀ x ∈ L₁ ++ L₂.map (fun z => -z), x ∈ S := by
    intro x hx
    rcases List.mem_append.1 hx with hx | hx
    · exact h₁ x hx
    · obtain ⟨y, hy, rfl⟩ := List.mem_map.1 hx
      exact hS y (h₂ y hy)
  have h0 : (L₁ ++ L₂.map (fun z => -z)).sum = 0 := by
    rw [List.sum_append, TheoremW.list_sum_map_neg, he, add_neg_cancel]
  obtain ⟨m, hm⟩ := hcl _ hL h0
  rw [List.map_append, List.sum_append, sum_x0_neg x0 hneg L₂ h₂] at hm
  unfold V
  rw [← sub_eq_zero, ← AddCircle.coe_sub, AddCircle.coe_eq_zero_iff]
  refine ⟨m - L₂.length, ?_⟩
  rw [zsmul_eq_mul, mul_one]
  push_cast
  linarith

/-- A function with `x0 (-s) = 1 - x0 s` on `S` and integer sums along closed walks with steps in
`S` induces a character of `closure S`. -/
lemma exists_char (hS : ∀ s ∈ S, -s ∈ S) (hneg : ∀ s ∈ S, x0 (-s) = 1 - x0 s)
    (hcl : ∀ L : List G, (∀ x ∈ L, x ∈ S) → L.sum = 0 → ∃ m : ℤ, (L.map x0).sum = m) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ξ ⟨s, AddSubgroup.subset_closure hs⟩ = (x0 s : AddCircle (1 : ℝ)) := by
  classical
  set H := AddSubgroup.closure S
  have hlist : ∀ γ : H, ∃ L : List G, (∀ x ∈ L, x ∈ S) ∧ L.sum = γ :=
    fun γ => TheoremW.exists_list_of_mem_closure hS γ.2
  choose Lf hLf hLs using hlist
  let ξf : H → AddCircle (1 : ℝ) := fun γ => V x0 (Lf γ)
  have key : ∀ (L : List G), (∀ x ∈ L, x ∈ S) → ∀ γ : H, L.sum = γ → ξf γ = V x0 L :=
    fun L hL γ he => V_eq x0 hS hneg hcl (hLf γ) hL ((hLs γ).trans he.symm)
  refine ⟨AddMonoidHom.mk' ξf (fun a b => ?_), fun s hs => ?_⟩
  · rw [key (Lf a ++ Lf b) (fun x hx => (List.mem_append.1 hx).elim (hLf a x) (hLf b x))
      (a + b) (by rw [List.sum_append, hLs a, hLs b, AddSubgroup.coe_add]), V_append]
  · show ξf _ = _
    rw [key [s] (by simpa using hs) _ (by simp)]
    simp [V]

end character

/-! ### Main theorem -/

lemma abs_sig_le (c : G → ZMod 3) (s x : G) : |(TheoremW.sig c x s : ℝ)| ≤ 1 := by
  rcases TheoremW.sig_eq_one_or c x s with h | h <;> simp [h]

/-- **Theorem W** (case `k = 1`) for every abelian group `G` and finite symmetric `S ⊆ G`: a proper
colouring `c : G → ZMod 3` of `Cay(G, S)` gives a character `ξ : closure S →+ ℝ/ℤ` such that every
`ξ s` (`s ∈ S`) is represented by a real number in `[1/3, 2/3]`. -/
theorem theoremW {G : Type*} [AddCommGroup G]
    (S : Set G) (hSfin : S.Finite) (hS : ∀ s ∈ S, -s ∈ S)
    (c : G → ZMod 3) (hc : ∀ g, ∀ s ∈ S, c (g + s) ≠ c g) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, 1/3 ≤ x ∧ x ≤ 2/3 ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  obtain ⟨F, hFb, hF'⟩ := exists_F (fun x s => (TheoremW.sig c x s : ℝ))
    (fun x L => (TheoremW.W c x L : ℝ)) (fun x => by simp) (fun x s L => by simp)
    (M := 1) (fun x s => abs_sig_le c s x) hSfin (a := -1) (b := 1)
    (fun x s _ => abs_le.1 (abs_sig_le c s x))
  have hF : ∀ L : List G, (∀ x ∈ L, x ∈ S) → L.sum = 0 → (L.map F).sum = TheoremW.W c 0 L :=
    fun L hL h0 => hF' L hL (fun h hh => by
      have := TheoremW.W_shift_closure hS hc hL h0 hh 0
      rw [zero_add] at this
      rw [this])
  set x0 : G → ℝ := fun s => (F s + 3) / 6 with hx0
  have hFneg : ∀ s ∈ S, F (-s) = - F s := by
    intro s hs
    have h := hF [s, -s] (by simpa using ⟨hs, hS s hs⟩) (by simp)
    rw [TheoremW.W_cons, TheoremW.W_cons, TheoremW.W_nil, TheoremW.sig_neg hc hs 0] at h
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil] at h
    push_cast at h
    linarith
  have hneg : ∀ s ∈ S, x0 (-s) = 1 - x0 s := by
    intro s hs
    simp only [hx0, hFneg s hs]
    ring
  have hcl : ∀ L : List G, (∀ x ∈ L, x ∈ S) → L.sum = 0 → ∃ m : ℤ, (L.map x0).sum = m := by
    intro L hL h0
    obtain ⟨m, hm⟩ := TheoremW.W_closed hc hL h0 0
    refine ⟨m, ?_⟩
    have h := hF L hL h0
    rw [hx0, TheoremW.list_sum_affine L F 3 6, h, hm]
    push_cast
    ring
  obtain ⟨ξ, hξ⟩ := exists_char x0 hS hneg hcl
  refine ⟨ξ, fun s hs => ⟨x0 s, ?_, ?_, (hξ s hs).symm⟩⟩
  · have := (hFb s hs).1
    simp only [hx0]
    linarith
  · have := (hFb s hs).2
    simp only [hx0]
    linarith

/-- Sanity check: the hypotheses are satisfiable for an infinite group (`G = ℤ`, `S = {1, -1}`,
`c n = n mod 3`). -/
example : ∃ ξ : AddSubgroup.closure ({1, -1} : Set ℤ) →+ AddCircle (1 : ℝ),
    ∀ s (hs : s ∈ ({1, -1} : Set ℤ)), ∃ x : ℝ, 1/3 ≤ x ∧ x ≤ 2/3 ∧
      (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  refine theoremW ({1, -1} : Set ℤ) (Set.toFinite _) ?_ (fun n => (n : ZMod 3)) ?_
  · intro s hs
    rcases hs with rfl | rfl
    · exact Or.inr rfl
    · exact Or.inl rfl
  · intro g s hs h
    rcases hs with rfl | rfl
    · push_cast at h
      exact absurd (add_eq_left.1 h) (by decide)
    · push_cast at h
      exact absurd (add_eq_left.1 h) (by decide)

/-! ### Circular cliques -/

/-- **Theorem W for circular cliques `K_{p/q}`, `p < 4q`**, for every abelian group `G` and finite
symmetric `S ⊆ G`: a homomorphism `c : Cay(G, S) → K_{p/q}` gives a character
`ξ : closure S →+ ℝ/ℤ` such that every `ξ s` (`s ∈ S`) is represented by a real number in
`[q/p, 1 - q/p]`. This is `TheoremWplus.theoremWplus_general` without finiteness of `G`; the
case `p = 3`, `q = 1` is `theoremW` (see the example below). -/
theorem theoremWplus {G : Type*} [AddCommGroup G]
    (S : Set G) (hSfin : S.Finite) (hS : ∀ s ∈ S, -s ∈ S)
    (p q : ℕ) (hp0 : 0 < p) (h4 : p < 4 * q)
    (c : G → ZMod p)
    (hc : ∀ g, ∀ s ∈ S, q ≤ (c (g + s) - c g).val ∧ (c (g + s) - c g).val ≤ p - q) :
    ∃ ξ : AddSubgroup.closure S →+ AddCircle (1 : ℝ),
      ∀ s (hs : s ∈ S), ∃ x : ℝ, (q : ℝ) / p ≤ x ∧ x ≤ 1 - (q : ℝ) / p ∧
        (x : AddCircle (1 : ℝ)) = ξ ⟨s, AddSubgroup.subset_closure hs⟩ := by
  have : NeZero p := ⟨hp0.ne'⟩
  have hp : (0 : ℝ) < p := by exact_mod_cast hp0
  have hdM : ∀ x s, |(TheoremWplus.dl c x s : ℝ)| ≤ p := fun x s => by
    have h := ZMod.val_lt (c (x + s) - c x)
    unfold TheoremWplus.dl
    rw [abs_of_nonneg (by positivity)]
    push_cast
    exact_mod_cast h.le
  obtain ⟨F, hFb, hF'⟩ := exists_F (fun x s => (TheoremWplus.dl c x s : ℝ))
    (fun x L => (TheoremWplus.W c x L : ℝ)) (fun x => by simp) (fun x s L => by simp)
    hdM hSfin (a := q) (b := (p : ℝ) - q)
    (fun x s hs => by
      have h := TheoremWplus.dl_bounds hc (g := x) hs
      constructor <;> [exact_mod_cast h.1; exact_mod_cast h.2])
  have hF : ∀ L : List G, (∀ x ∈ L, x ∈ S) → L.sum = 0 →
      (L.map F).sum = TheoremWplus.W c 0 L :=
    fun L hL h0 => hF' L hL (fun h hh => by
      have := TheoremWplus.W_shift_closure hS hc h4 hL h0 hh 0
      rw [zero_add] at this
      rw [this])
  set x0 : G → ℝ := fun s => F s / p with hx0
  have hneg : ∀ s ∈ S, x0 (-s) = 1 - x0 s := by
    intro s hs
    have h := hF [s, -s] (by simpa using ⟨hs, hS s hs⟩) (by simp)
    rw [TheoremWplus.W_cons, TheoremWplus.W_cons, TheoremWplus.W_nil,
      TheoremWplus.dl_neg hc h4 hs 0] at h
    simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil] at h
    push_cast at h
    simp only [hx0]
    field_simp
    linarith
  have hcl : ∀ L : List G, (∀ x ∈ L, x ∈ S) → L.sum = 0 → ∃ m : ℤ, (L.map x0).sum = m := by
    intro L hL h0
    obtain ⟨m, hm⟩ := TheoremWplus.W_closed c h0 0
    refine ⟨m, ?_⟩
    rw [hx0, TheoremWplus.list_sum_div L F p, hF L hL h0, hm, div_eq_iff hp.ne']
    push_cast
    ring
  obtain ⟨ξ, hξ⟩ := exists_char x0 hS hneg hcl
  refine ⟨ξ, fun s hs => ⟨x0 s, ?_, ?_, (hξ s hs).symm⟩⟩
  · simp only [hx0]
    exact div_le_div_of_nonneg_right (hFb s hs).1 hp.le
  · simp only [hx0]
    rw [one_sub_div hp.ne']
    exact div_le_div_of_nonneg_right (hFb s hs).2 hp.le

/-- Sanity check: `p = 3`, `q = 1` (`K_{3/1} = K_3`) recovers `theoremW`. -/
example {G : Type*} [AddCommGroup G]
    (S : Set G) (hSfin : S.Finite) (hS : ∀ s ∈ S, -s ∈ S)
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
  obtain ⟨ξ, hξ⟩ := theoremWplus S hSfin hS 3 1 (by norm_num) (by norm_num) c hc'
  refine ⟨ξ, fun s hs => ?_⟩
  obtain ⟨x, h1, h2, h3⟩ := hξ s hs
  refine ⟨x, ?_, ?_, h3⟩
  · norm_num at h1; linarith
  · norm_num at h2; linarith

end TheoremWInf
