import LocalColouring

/-!
# Planes over real quadratic fields

Tools for the theorems `χ(ℚ(√d)²) = 4` of `notes/quadratic_planes.md`. `QuadraticPlanes.L d` is `ℚ(√d) ⊆ ℝ`.

*Upper bounds.* Each is a reduction at a place: colour a point by a colouring of a finite plane at the
residues of `z - rep z`, where `rep z` is a fixed representative of `z + O²`.
* `colorable_four` (Moorhouse, Lemma 8.2 at 7). If `d ≡ s² (mod 7)` with `7 ∤ s`, the unit-distance graph of
  `ℚ(√d)²` is 4-colourable.
  1. Chevalley's extension theorem gives a valuation subring `O` of `L d` with `7` in its maximal ideal `𝔪`.
  2. Its residue field is `𝔽₇` (`exists_int_sub_mem`): `√d ∈ O`, and one of `±√d - s` lies in `𝔪`, say `t - s`.
     Hensel's lemma gives integers `bₖ` with `t - bₖ ∈ 7ᵏ⁺¹O`, and writing `x = (n₀ + n₁t)/(7ʲD')` with `7 ∤ D'`
     shows that every `x ∈ O` is `≡` an integer modulo `𝔪`.
  3. If `x² + y² = 1`, then `x, y ∈ O` and `x̄² + ȳ² = 1` in `𝔽₇`, since `-1` is not a square modulo 7
     (`core_lemma_seven`). A proper 4-colouring of the unit-distance graph of `𝔽₇²`
     (`data/quadratic_planes/finite_planes.json`) colours `L²` (`sumSqGraph_colorable`).
* `colorable_four_ramified` (the same at 7 when `d = 7d'`, `7 ∤ d'`). Now `√d ∈ 𝔪`, and descent on the power of 7
  in the denominator shows again that the residue field is `𝔽₇` (`exists_int_sub_mem_ramified`).
* `colorable_four_two` (K. G. Fischer, 1990, Theorem 10: reduction at 2). If `d ≡ 3 (mod 8)`, take `O` with `2` in
  its maximal ideal. Its residue field is `𝔽₂`, since `π = √d - 1 ∈ 𝔪` and `N₀ + N₁π ∈ 2O` forces both `Nᵢ` even
  (`mem_nonunits_or_sub_one_mem_two`). In the coordinates `a = x + y/√d`, `b = 2y/√d` the squared distance is
  `a² - ab + κb²` with `κ = (d + 1)/4` odd (`toKHom`), and over `𝔽₂` the form `a² + ab + b²` has no nontrivial
  zero; so the unit vectors are integral with nonzero residue, and the residues in `𝔽₂²` colour `L²`
  (`kGraph_colorable`).

*Lower bounds.* A point `[a, b, c, e]` with denominator `D` stands for `((a + b√d)/D, (c + e√d)/D)`; two points
are at distance 1 when their differences satisfy two integer identities (`adj_pt`), which the kernel checks
for a list of edges (`checkEdges`). From a 3-colouring of the graph, `clause_facts` builds a valuation that
satisfies every clause of its colouring formula; each field's file refutes that formula with an LRAT proof
(Mathlib's `lrat_proof`).
-/

namespace QuadraticPlanes

open LocalColouring IsLocalRing

/-- `ℚ(√d)`, as a subfield of `ℝ`. -/
noncomputable def L (d : ℕ) : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√(d : ℝ)}

lemma sqrt_mem (d : ℕ) : √(d : ℝ) ∈ L d := IntermediateField.subset_adjoin ℚ _ (by simp)

/-- `√d` as an element of `L d`. -/
noncomputable def r (d : ℕ) : L d := ⟨√(d : ℝ), sqrt_mem d⟩

@[simp] lemma coe_r (d : ℕ) : (r d : ℝ) = √(d : ℝ) := rfl

lemma sqrt_sq (d : ℕ) : √(d : ℝ) ^ 2 = d := Real.sq_sqrt (Nat.cast_nonneg d)

lemma r_sq (d : ℕ) : r d ^ 2 = d := Subtype.ext (by push_cast; exact sqrt_sq d)

/-! ## Lower bounds: unit distances and edge lists -/

/-- The point `[a, b, c, e]` with denominator `D`, that is `((a + b√d)/D, (c + e√d)/D)`. -/
noncomputable def pt (d D : ℕ) (v : ℤ × ℤ × ℤ × ℤ) : L d × L d :=
  ((v.1 + v.2.1 * r d) / D, (v.2.2.1 + v.2.2.2 * r d) / D)

/-- The two integer identities for a unit distance between `[a, b, c, e]/D` and `[a', b', c', e']/D`, as a
Boolean (products written out, so that the kernel evaluates it quickly). -/
def unitPairB (d D : ℕ) (v w : ℤ × ℤ × ℤ × ℤ) : Bool :=
  let da := v.1 - w.1
  let db := v.2.1 - w.2.1
  let dc := v.2.2.1 - w.2.2.1
  let de := v.2.2.2 - w.2.2.2
  (da * da + d * (db * db) + dc * dc + d * (de * de) == (D : ℤ) * D) && (da * db + dc * de == 0)

lemma unitPairB_spec {d D : ℕ} {v w : ℤ × ℤ × ℤ × ℤ} (h : unitPairB d D v w = true) :
    (v.1 - w.1) * (v.1 - w.1) + d * ((v.2.1 - w.2.1) * (v.2.1 - w.2.1))
        + (v.2.2.1 - w.2.2.1) * (v.2.2.1 - w.2.2.1) + d * ((v.2.2.2 - w.2.2.2) * (v.2.2.2 - w.2.2.2))
        = (D : ℤ) * D ∧
      (v.1 - w.1) * (v.2.1 - w.2.1) + (v.2.2.1 - w.2.2.1) * (v.2.2.2 - w.2.2.2) = 0 := by
  simpa [unitPairB] using h

/-- A unit distance in `L²` is a unit distance in `ℝ²`. -/
lemma adj_of_sq_eq_one {K : IntermediateField ℚ ℝ} {p q : K × K}
    (h : (p.1 - q.1) ^ 2 + (p.2 - q.2) ^ 2 = 1) : (unitDistGraph K).Adj p q := by
  show ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1
  have := congrArg (fun x : K => (x : ℝ)) h
  push_cast at this
  exact this

/-- Since `(√d)² = d`, the two identities say that the distance is 1. -/
lemma adj_pt {d D : ℕ} (hD : D ≠ 0) {v w : ℤ × ℤ × ℤ × ℤ} (h : unitPairB d D v w = true) :
    (unitDistGraph (L d)).Adj (pt d D v) (pt d D w) := by
  obtain ⟨h1, h2⟩ := unitPairB_spec h
  obtain ⟨a, b, c, e⟩ := v
  obtain ⟨a', b', c', e'⟩ := w
  dsimp only at h1 h2
  have h1' : ((a : L d) - a') * ((a : L d) - a') + d * (((b : L d) - b') * ((b : L d) - b'))
      + ((c : L d) - c') * ((c : L d) - c') + d * (((e : L d) - e') * ((e : L d) - e'))
      = (D : L d) * D := by
    exact_mod_cast h1
  have h2' : ((a : L d) - a') * ((b : L d) - b') + ((c : L d) - c') * ((e : L d) - e') = 0 := by
    exact_mod_cast h2
  have hD' : (D : L d) ≠ 0 := by exact_mod_cast hD
  apply adj_of_sq_eq_one
  simp only [pt]
  rw [div_sub_div_same, div_sub_div_same, div_pow, div_pow, ← add_div, div_eq_one_iff_eq (pow_ne_zero 2 hD')]
  linear_combination h1' + 2 * r d * h2' + (((b : L d) - b') ^ 2 + ((e : L d) - e') ^ 2) * r_sq d

/-- Every listed edge satisfies the two identities. -/
def checkEdges {n : ℕ} (d D : ℕ) (P : Fin n → ℤ × ℤ × ℤ × ℤ) : List (Fin n × Fin n) → Bool
  | [] => true
  | e :: es => unitPairB d D (P e.1) (P e.2) && checkEdges d D P es

lemma checkEdges_mem {n : ℕ} {d D : ℕ} {P : Fin n → ℤ × ℤ × ℤ × ℤ} :
    ∀ {l : List (Fin n × Fin n)}, checkEdges d D P l = true → ∀ e ∈ l, unitPairB d D (P e.1) (P e.2) = true
  | [], _, _, he => by simp at he
  | _ :: _, h, e, he => by
    simp only [checkEdges, Bool.and_eq_true] at h
    rcases List.mem_cons.mp he with rfl | he
    · exact h.1
    · exact checkEdges_mem h.2 e he

/-- No listed edge is a loop. -/
def noLoops {n : ℕ} : List (Fin n × Fin n) → Bool
  | [] => true
  | e :: es => (e.1.val != e.2.val) && noLoops es

lemma noLoops_mem {n : ℕ} : ∀ {l : List (Fin n × Fin n)}, noLoops l = true → ∀ e ∈ l, e.1 ≠ e.2
  | [], _, _, he => by simp at he
  | _ :: _, h, e, he => by
    simp only [noLoops, Bool.and_eq_true, bne_iff_ne, ne_eq] at h
    rcases List.mem_cons.mp he with rfl | he
    · exact fun h' => h.1 (congrArg Fin.val h')
    · exact noLoops_mem h.2 e he

/-- The listed edges are unit distances in `L d`. -/
lemma pt_adj {n : ℕ} {d D : ℕ} (hD : D ≠ 0) {P : Fin n → ℤ × ℤ × ℤ × ℤ} {E : List (Fin n × Fin n)}
    (hE : checkEdges d D P E = true) : ∀ e ∈ E, (unitDistGraph (L d)).Adj (pt d D (P e.1)) (pt d D (P e.2)) :=
  fun e he => adj_pt hD (checkEdges_mem hE e he)

lemma adj_of_mem {n : ℕ} {E : List (Fin n × Fin n)} (hE : noLoops E = true) {i j : Fin n} (h : (i, j) ∈ E) :
    (edgeGraph E).Adj i j := by
  rw [edgeGraph, SimpleGraph.fromRel_adj]
  refine ⟨fun hij => ?_, Or.inl h⟩
  subst hij
  exact noLoops_mem hE _ h rfl

/-! ## Lower bounds: from a 3-colouring to a satisfying assignment -/

/-- Is `(a, b)` a listed edge? -/
def memE {n : ℕ} (a b : ℕ) : List (Fin n × Fin n) → Bool
  | [] => false
  | e :: es => (Nat.beq e.1.val a && Nat.beq e.2.val b) || memE a b es

lemma memE_spec {n : ℕ} {a b : ℕ} : ∀ {l : List (Fin n × Fin n)}, memE a b l = true →
    ∃ e ∈ l, e.1.val = a ∧ e.2.val = b
  | [], h => by simp [memE] at h
  | e :: _, h => by
    simp only [memE, Bool.or_eq_true, Bool.and_eq_true] at h
    rcases h with ⟨h1, h2⟩ | h
    · exact ⟨e, List.mem_cons_self .., Nat.eq_of_beq_eq_true h1, Nat.eq_of_beq_eq_true h2⟩
    · obtain ⟨e', he', h'⟩ := memE_spec h
      exact ⟨e', List.mem_cons_of_mem _ he', h'⟩

/-- The edge clauses `¬x ∨ ¬y` of a colouring formula: the same colour `x % 3 = y % 3` at the ends of an
edge. -/
def goodEdge {n : ℕ} (E : List (Fin n × Fin n)) (x y : ℕ) : Bool :=
  Nat.beq (x % 3) (y % 3) && (memE (x / 3) (y / 3) E || memE (y / 3) (x / 3) E)

/-- Renaming the colours so that `a ↦ 0` and `b ↦ 1`. -/
def rename (a b x : Fin 3) : Fin 3 := if x = a then 0 else if x = b then 1 else 2

lemma rename_ne : ∀ a b x y : Fin 3, a ≠ b → x ≠ y → rename a b x ≠ rename a b y := by decide

/-- A 3-colouring of `edgeGraph E`, with its colours renamed so that the edge `u₀v₀` gets 0 and 1, gives a
valuation (variable `x` says "vertex `x / 3` has colour `x % 3`") under which every clause of the colouring
formula holds: the vertex clauses `x ∨ x + 1 ∨ x + 2`, the edge clauses `¬x ∨ ¬y`, and the unit clauses
`3u₀` and `3v₀ + 1`. -/
theorem clause_facts {n : ℕ} {E : List (Fin n × Fin n)} (hE : noLoops E = true) {u₀ v₀ : Fin n}
    (h₀ : (u₀, v₀) ∈ E) (C : (edgeGraph E).Coloring (Fin 3)) :
    ∃ V : ℕ → Prop,
      (∀ x, (¬V x ∧ ¬V (x + 1) ∧ ¬V (x + 2)) → x % 3 = 0 → x / 3 < n → False) ∧
      (∀ x y, (V x ∧ V y) → goodEdge E x y = true → False) ∧
      V (3 * u₀.val) ∧ V (3 * v₀.val + 1) := by
  have hab : C u₀ ≠ C v₀ := C.valid (adj_of_mem hE h₀)
  let g : Fin n → Fin 3 := fun i => rename (C u₀) (C v₀) (C i)
  have hg0 : g u₀ = 0 := by simp [g, rename]
  have hg1 : g v₀ = 1 := by simp [g, rename, Ne.symm hab]
  have hgE : ∀ i j : Fin n, (i, j) ∈ E → g i ≠ g j := fun i j h =>
    rename_ne _ _ _ _ hab (C.valid (adj_of_mem hE h))
  let V : ℕ → Prop := fun x => if h : x / 3 < n then ((g ⟨x / 3, h⟩ : Fin 3) : ℕ) = x % 3 else False
  have hv : ∀ x (h : x / 3 < n), V x ↔ ((g ⟨x / 3, h⟩ : Fin 3) : ℕ) = x % 3 := fun x h => by
    simp only [V, h, dite_true]
  refine ⟨V, ?_, ?_, ?_, ?_⟩
  · -- every vertex clause holds
    intro x ⟨n0, n1, n2⟩ hx hlt
    rw [hv x hlt] at n0
    rw [hv (x + 1) (by omega)] at n1
    rw [hv (x + 2) (by omega)] at n2
    have h1 : (x + 1) / 3 = x / 3 := by omega
    have h2 : (x + 2) / 3 = x / 3 := by omega
    simp only [h1, h2] at n1 n2
    have : ((g ⟨x / 3, hlt⟩ : Fin 3) : ℕ) < 3 := (g ⟨x / 3, hlt⟩).2
    omega
  · -- every edge clause holds
    have hEdge1 : ∀ x y, V x → V y → Nat.beq (x % 3) (y % 3) = true → memE (x / 3) (y / 3) E = true →
        False := by
      intro x y vx vy hxy hE'
      have hxy' := Nat.eq_of_beq_eq_true hxy
      obtain ⟨⟨i, j⟩, he, hi, hj⟩ := memE_spec hE'
      simp only at hi hj
      have hx : x / 3 < n := hi ▸ i.2
      have hy : y / 3 < n := hj ▸ j.2
      rw [hv x hx] at vx
      rw [hv y hy] at vy
      have e1 : (⟨x / 3, hx⟩ : Fin n) = i := Fin.ext hi.symm
      have e2 : (⟨y / 3, hy⟩ : Fin n) = j := Fin.ext hj.symm
      rw [e1] at vx
      rw [e2] at vy
      exact hgE i j he (Fin.ext (by omega))
    intro x y ⟨vx, vy⟩ hc
    simp only [goodEdge, Bool.and_eq_true, Bool.or_eq_true] at hc
    obtain ⟨hxy, hE' | hE'⟩ := hc
    · exact hEdge1 x y vx vy hxy hE'
    · exact hEdge1 y x vy vx (by
        have := Nat.eq_of_beq_eq_true hxy
        rw [this]
        exact Nat.beq_refl _) hE'
  · -- the unit clauses
    have hu : (3 * u₀.val) / 3 < n := by omega
    rw [hv _ hu]
    have : (⟨(3 * u₀.val) / 3, hu⟩ : Fin n) = u₀ := Fin.ext (by simp only; omega)
    rw [this, hg0, Fin.val_zero]
    omega
  · have hu : (3 * v₀.val + 1) / 3 < n := by omega
    rw [hv _ hu]
    have : (⟨(3 * v₀.val + 1) / 3, hu⟩ : Fin n) = v₀ := Fin.ext (by simp only; omega)
    rw [this, hg1, Fin.val_one]
    omega

/-! ## Upper bound: reduction at a place over 7 -/

/-! ### The graph of `x² + y²` over a field with a place of residue field `𝔽₇` -/

section Seven

variable {F : Type*} [Field F] (O : ValuationSubring F) (red : O →+* ZMod 7)

/-- `-1` is not a square modulo 7. -/
lemma zmod7_one_add_sq_ne : ∀ t : ZMod 7, 1 + t ^ 2 ≠ 0 := by decide

lemma one_add_sq_ne (hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0) (t x : O)
    (hx : x ∈ maximalIdeal O) : 1 + t ^ 2 ≠ x ^ 2 := by
  intro h
  have h' := congrArg red h
  rw [map_add, map_one, map_pow, map_pow, hred x hx] at h'
  exact zmod7_one_add_sq_ne _ (by rw [h']; ring)

/-- If `a² + b² = 1`, then `a ∈ O`. -/
lemma mem_of_sum_sq_eq_one (hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0) {a b : F}
    (h : a ^ 2 + b ^ 2 = 1) : a ∈ O := by
  by_contra ha
  have ha0 : a ≠ 0 := fun h0 => ha (h0 ▸ O.zero_mem)
  have hx : a⁻¹ ∈ O.nonunits := O.inv_mem_nonunits_iff.mpr (Or.inr ha)
  obtain ⟨hxO, hxm⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp hx
  by_cases ht : b / a ∈ O
  · -- `1 + (b/a)² = (a⁻¹)²`, with `a⁻¹ ∈ 𝔪`
    apply one_add_sq_ne O red hred ⟨b / a, ht⟩ ⟨a⁻¹, hxO⟩ hxm
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h
  · -- otherwise `s = a/b ∈ O`, and `1 + s² = (s a⁻¹)²` with `s a⁻¹ = b⁻¹ ∈ 𝔪`
    have hs : a / b ∈ O := by
      have := (O.mem_or_inv_mem (b / a)).resolve_left ht
      rwa [inv_div] at this
    have hb0 : b ≠ 0 := by
      rintro rfl
      exact ht (by simp)
    apply one_add_sq_ne O red hred ⟨a / b, hs⟩ (⟨a / b, hs⟩ * ⟨a⁻¹, hxO⟩)
      (Ideal.mul_mem_left _ _ hxm)
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h

/-- **Core lemma at 7.** If `a² + b² = 1`, then `a, b ∈ O` and `(red a)² + (red b)² = 1`. -/
theorem core_lemma_seven (hker : RingHom.ker red = maximalIdeal O) {a b : F}
    (h : a ^ 2 + b ^ 2 = 1) :
    ∃ (ha : a ∈ O) (hb : b ∈ O), red ⟨a, ha⟩ ^ 2 + red ⟨b, hb⟩ ^ 2 = 1 := by
  have hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0 := fun x hx => by
    rw [← hker] at hx
    exact hx
  have ha := mem_of_sum_sq_eq_one O red hred h
  have hb := mem_of_sum_sq_eq_one O red hred (a := b) (b := a) (by linear_combination h)
  refine ⟨ha, hb, ?_⟩
  have e : (⟨a, ha⟩ : O) ^ 2 + ⟨b, hb⟩ ^ 2 = 1 := Subtype.ext (by simpa using h)
  have e' := congrArg red e
  rwa [map_add, map_pow, map_pow, map_one] at e'

open Classical in
/-- `red` extended by `0` outside `O`. -/
noncomputable def redExt7 (x : F) : ZMod 7 := if h : x ∈ O then red ⟨x, h⟩ else 0

lemma redExt7_add {x y : F} (hx : x ∈ O) (hy : y ∈ O) :
    redExt7 O red (x + y) = redExt7 O red x + redExt7 O red y := by
  simp only [redExt7, hx, hy, add_mem hx hy, dite_true]
  rw [← map_add]
  rfl

lemma redExt7_of_mem {x : F} (hx : x ∈ O) : redExt7 O red x = red ⟨x, hx⟩ := by
  simp only [redExt7, hx, dite_true]

/-- A proper 4-colouring of the unit-distance graph of `𝔽₇²`
(`data/quadratic_planes/finite_planes.json`, vertex `7x + y`). -/
def f7Colour : Fin 49 → Fin 4 := ![
    3, 2, 0, 1, 2, 1, 0, 0, 1, 2, 0, 3, 2, 1, 2, 3, 0, 3, 1, 0, 3, 3, 2, 1, 0,
    3, 1, 0, 1, 0, 3, 2, 1, 0, 2, 3, 2, 1, 0, 3, 2, 1, 1, 3, 2, 3, 0, 3, 2]

/-- The colour of a point of `𝔽₇²`. -/
def c7 (z : ZMod 7 × ZMod 7) : Fin 4 :=
  f7Colour ⟨7 * z.1.val + z.2.val, by
    have h1 := ZMod.val_lt z.1
    have h2 := ZMod.val_lt z.2
    omega⟩

/-- `c7` is a proper colouring: points of `𝔽₇²` at "distance" 1 get different colours. -/
lemma c7_proper : ∀ u w : ZMod 7 × ZMod 7, (u.1 - w.1) ^ 2 + (u.2 - w.2) ^ 2 = 1 → c7 u ≠ c7 w := by
  decide +kernel

/-- The colouring of `F²`: `c7` of the residues of `z - rep z`. -/
noncomputable def colour7 (z : F × F) : Fin 4 :=
  c7 (redExt7 O red (z.1 - rep O z.1), redExt7 O red (z.2 - rep O z.2))

/-- The graph on `F × F` of the form `x² + y²`: `p ~ q` iff `(p₁ - q₁)² + (p₂ - q₂)² = 1`. -/
def sumSqGraph (F : Type*) [Field F] : SimpleGraph (F × F) where
  Adj p q := (p.1 - q.1) ^ 2 + (p.2 - q.2) ^ 2 = 1
  symm := ⟨fun p q h => by linear_combination h⟩
  loopless := ⟨fun p h => by simp at h⟩

lemma colour7_adj (hker : RingHom.ker red = maximalIdeal O) {p q : F × F}
    (hpq : (sumSqGraph F).Adj p q) : colour7 O red p ≠ colour7 O red q := by
  obtain ⟨ha, hb, hsq⟩ := core_lemma_seven O red hker hpq
  have e1 : p.1 - rep O p.1 = (p.1 - q.1) + (q.1 - rep O q.1) := by
    rw [rep_eq_of_sub_mem O ha]; ring
  have e2 : p.2 - rep O p.2 = (p.2 - q.2) + (q.2 - rep O q.2) := by
    rw [rep_eq_of_sub_mem O hb]; ring
  unfold colour7
  rw [e1, redExt7_add O red ha (sub_rep_mem O _), e2, redExt7_add O red hb (sub_rep_mem O _),
    redExt7_of_mem O red ha, redExt7_of_mem O red hb]
  apply c7_proper
  simpa only [add_sub_cancel_right] using hsq

/-- The graph of `x² + y²` is 4-colourable when `F` has a place with residue field `𝔽₇`. -/
theorem sumSqGraph_colorable (hker : RingHom.ker red = maximalIdeal O) :
    (sumSqGraph F).Colorable 4 := by
  have hc : (sumSqGraph F).Colorable (Fintype.card (Fin 4)) :=
    (SimpleGraph.Coloring.mk (colour7 O red) (fun h => colour7_adj O red hker h)).colorable
  simpa using hc

end Seven

/-- Every field of characteristic zero has a valuation subring with `7` in its maximal ideal
(Chevalley's extension theorem, applied to `ℤ` and the ideal `(7)`). -/
theorem exists_valuationSubring_seven_mem_nonunits (F : Type*) [Field F] [CharZero F] :
    ∃ O : ValuationSubring F, (7 : F) ∈ O.nonunits := by
  have hne : Ideal.span {(7 : (⊥ : Subring F))} ≠ ⊤ := by
    rw [Ne, Ideal.span_singleton_eq_top]
    rintro ⟨u, hu⟩
    have hc7 : ((7 : (⊥ : Subring F)) : F) = 7 := rfl
    obtain ⟨n, hn⟩ := Subring.mem_bot.mp (u⁻¹ : (⊥ : Subring F)ˣ).1.2
    have h1 : ((u : (⊥ : Subring F)) : F) * ((u⁻¹ : (⊥ : Subring F)ˣ) : F) = 1 := by
      rw [← Subring.coe_mul, Units.mul_inv]; rfl
    rw [hu, ← hn, hc7] at h1
    have h2 : ((7 * n : ℤ) : F) = ((1 : ℤ) : F) := by push_cast; exact h1
    have h3 : (7 * n : ℤ) = 1 := Int.cast_injective h2
    omega
  obtain ⟨O, -, hO⟩ :=
    Ideal.image_subset_nonunits_valuationSubring (A := (⊥ : Subring F)) _ hne
  exact ⟨O, hO ⟨7, Ideal.subset_span rfl, rfl⟩⟩

/-- A local ring in which every element is `≡` an integer modulo `𝔪`, and an integer lies in `𝔪`
exactly when 7 divides it, has a reduction map to `ZMod 7` with kernel the maximal ideal. -/
theorem exists_red_of_residue_seven {R : Type*} [CommRing R] [IsLocalRing R]
    (hsurj : ∀ x : R, ∃ n : ℤ, x - n ∈ maximalIdeal R)
    (hint : ∀ n : ℤ, (n : R) ∈ maximalIdeal R ↔ (7 : ℤ) ∣ n) :
    ∃ red : R →+* ZMod 7, RingHom.ker red = maximalIdeal R := by
  let f : ℤ →+* ResidueField R := Int.castRingHom _
  have hf_apply : ∀ n : ℤ, f n = residue R n := fun n => by simp [f]
  have hf : Function.Surjective f := by
    intro y
    obtain ⟨x, rfl⟩ := residue_surjective y
    obtain ⟨n, hn⟩ := hsurj x
    refine ⟨n, ?_⟩
    rw [hf_apply, ← sub_eq_zero, ← map_sub, residue_eq_zero_iff, ← neg_mem_iff, neg_sub]
    exact hn
  have hker : RingHom.ker f = Ideal.span {((7 : ℕ) : ℤ)} := by
    ext n
    rw [RingHom.mem_ker, Ideal.mem_span_singleton, hf_apply, residue_eq_zero_iff, hint]
    norm_num
  let e1 : ℤ ⧸ RingHom.ker f ≃+* ResidueField R := RingHom.quotientKerEquivOfSurjective hf
  let e3 : ℤ ⧸ RingHom.ker f ≃+* ℤ ⧸ Ideal.span {((7 : ℕ) : ℤ)} := Ideal.quotEquivOfEq hker
  let e2 : ℤ ⧸ Ideal.span {((7 : ℕ) : ℤ)} ≃+* ZMod 7 := Int.quotientSpanNatEquivZMod 7
  refine ⟨((e2.toRingHom.comp e3.toRingHom).comp e1.symm.toRingHom).comp (residue R), ?_⟩
  ext x
  simp only [RingHom.mem_ker, RingHom.coe_comp, Function.comp_apply, RingEquiv.toRingHom_eq_coe,
    RingEquiv.coe_toRingHom, EmbeddingLike.map_eq_zero_iff, residue_eq_zero_iff]

/-- An integer prime to 7 has an inverse modulo 7. -/
lemma exists_inv_mod_seven {r : ℤ} (h : ¬ (7 : ℤ) ∣ r) : ∃ s k : ℤ, r * s = 1 + 7 * k := by
  have hdiv : 7 * (r / 7) + r % 7 = r := Int.mul_ediv_add_emod r 7
  have h0 : r % 7 ≠ 0 := fun h0 => h (Int.dvd_of_emod_eq_zero h0)
  have h1 : 0 ≤ r % 7 := Int.emod_nonneg _ (by norm_num)
  have h2 : r % 7 < 7 := Int.emod_lt_of_pos _ (by norm_num)
  generalize r % 7 = e at hdiv h0 h1 h2
  interval_cases e
  · exact absurd rfl h0
  · exact ⟨1, r / 7, by linear_combination -hdiv⟩
  · exact ⟨4, 4 * (r / 7) + 1, by linear_combination -4 * hdiv⟩
  · exact ⟨5, 5 * (r / 7) + 2, by linear_combination -5 * hdiv⟩
  · exact ⟨2, 2 * (r / 7) + 1, by linear_combination -2 * hdiv⟩
  · exact ⟨3, 3 * (r / 7) + 2, by linear_combination -3 * hdiv⟩
  · exact ⟨6, 6 * (r / 7) + 5, by linear_combination -6 * hdiv⟩

section SevenInt

variable {F : Type*} [Field F] (O : ValuationSubring F) (h7 : (7 : F) ∈ O.nonunits)
include h7

/-- An integer in `𝔪` is divisible by 7. -/
lemma seven_dvd_of_intCast_mem_nonunits {n : ℤ} (hn : (n : F) ∈ O.nonunits) : (7 : ℤ) ∣ n := by
  by_contra hnd
  obtain ⟨s, k, hs⟩ := exists_inv_mod_seven hnd
  apply one_notMem_nonunits O
  have hs' : (n : F) * s = 1 + 7 * k := by
    have := congrArg (fun z : ℤ => (z : F)) hs
    push_cast at this
    exact this
  rw [show (1 : F) = n * s - 7 * k by linear_combination -hs']
  exact sub_mem (mul_mem_nonunits_right O hn (intCast_mem O s))
    (mul_mem_nonunits_right O h7 (intCast_mem O k))

/-- If an integer is `7ᵏ` times an element of `O`, then `7ᵏ` divides it. -/
lemma seven_pow_dvd_of_eq [CharZero F] {k : ℕ} {n : ℤ} {γ : F} (hγ : γ ∈ O)
    (h : (n : F) = 7 ^ k * γ) : (7 : ℤ) ^ k ∣ n := by
  induction k generalizing n γ with
  | zero => simp
  | succ k ih =>
    have hn : (7 : ℤ) ∣ n := seven_dvd_of_intCast_mem_nonunits O h7 (by
      rw [h, pow_succ, mul_comm _ (7 : F), mul_assoc]
      exact mul_mem_nonunits_right O h7 (O.mul_mem _ _ (pow_mem (natCast_mem O 7) k) hγ))
    obtain ⟨m, rfl⟩ := hn
    have hm : (m : F) = 7 ^ k * γ := by
      apply mul_left_cancel₀ (show (7 : F) ≠ 0 by exact_mod_cast (by norm_num : (7 : ℤ) ≠ 0))
      push_cast at h
      linear_combination h
    obtain ⟨c, hc⟩ := ih hγ hm
    exact ⟨c, by rw [pow_succ]; linear_combination 7 * hc⟩

end SevenInt

/-! ### The residue field of a place of `ℚ(√d)` over 7 is `𝔽₇` when `d` is a nonzero square modulo 7 -/

/-- Every element of `L d` is `(n₀ + n₁t)/D`, for either sign `t = ±√d`. -/
lemma exists_int_combo_d {d : ℕ} {t : L d} (ht : t = r d ∨ t = -r d) (x : L d) :
    ∃ D : ℕ, 0 < D ∧ ∃ n₀ n₁ : ℤ, (D : L d) * x = n₀ + n₁ * t := by
  have hx : (x : ℝ) ∈ IntermediateField.adjoin ℚ ({√(d : ℝ), √(d : ℝ)} : Set ℝ) := by
    rw [Set.pair_eq_singleton]
    exact x.2
  obtain ⟨D, hD, n₀, n₁, n₂, n₃, h⟩ := exists_int_combo (m := d) (n := d)
    (by push_cast; exact sqrt_sq d) (by push_cast; exact sqrt_sq d) hx
  have h' : (D : L d) * x = (n₀ : L d) + n₃ * (r d * r d) + ((n₁ : L d) + n₂) * r d :=
    Subtype.ext (by push_cast [coe_r]; linear_combination h)
  rcases ht with rfl | rfl
  · exact ⟨D, hD, n₀ + d * n₃, n₁ + n₂, by
      rw [h']; push_cast; linear_combination (n₃ : L d) * r_sq d⟩
  · exact ⟨D, hD, n₀ + d * n₃, -(n₁ + n₂), by
      rw [h']; push_cast; linear_combination (n₃ : L d) * r_sq d⟩

section Residue

variable {d : ℕ} (O : ValuationSubring (L d))

lemma r_mem : r d ∈ O :=
  mem_of_sq_eq O O.zero_mem (natCast_mem O d) (by rw [zero_mul, zero_add]; exact r_sq d)

variable (h7 : (7 : L d) ∈ O.nonunits) {s : ℤ} (hs : ¬ (7 : ℤ) ∣ s) (hds : (7 : ℤ) ∣ (d : ℤ) - s ^ 2)
include h7 hs hds

omit hs in
/-- Choosing the sign of `√d`, we may assume `t - s ∈ 𝔪`. -/
lemma exists_t : ∃ t : L d, (t = r d ∨ t = -r d) ∧ t - s ∈ O.nonunits := by
  obtain ⟨c₀, hc₀⟩ := hds
  have hO := r_mem O
  have h : (r d - s) * (r d + s) ∈ O.nonunits := by
    have hc₀' : ((d : ℤ) : L d) - (s : L d) ^ 2 = 7 * c₀ := by exact_mod_cast hc₀
    rw [show (r d - s) * (r d + s) = 7 * c₀ by
      push_cast at hc₀'
      linear_combination r_sq d + hc₀']
    exact mul_mem_nonunits_right O h7 (intCast_mem O c₀)
  rcases mem_nonunits_or_of_mul_mem O (sub_mem hO (intCast_mem O s)) (add_mem hO (intCast_mem O s)) h
    with h | h
  · exact ⟨r d, Or.inl rfl, h⟩
  · refine ⟨-r d, Or.inr rfl, ?_⟩
    rw [show -r d - s = -(r d + s) by ring]
    exact neg_mem h

/-- **Hensel's lemma at 7.** If `t² = d` and `t - s ∈ 𝔪`, then for every `k` there is an integer
`b = s + 7m` with `t - b ∈ 7ᵏ⁺¹O`. -/
lemma hensel_seven {t : L d} (ht : t ^ 2 = d) (h2 : t - s ∈ O.nonunits) (k : ℕ) :
    ∃ m : ℤ, ∃ γ ∈ O, t - (s + 7 * m) = 7 ^ (k + 1) * γ := by
  obtain ⟨c₀, hc₀⟩ := hds
  have hc₀' : ((d : ℤ) : L d) - (s : L d) ^ 2 = 7 * c₀ := by exact_mod_cast hc₀
  push_cast at hc₀'
  -- `2s` is a unit modulo 7, with inverse `w`: `2sw = 1 + 7q`
  have h2s : ¬ (7 : ℤ) ∣ 2 * s := fun h => hs (by omega)
  obtain ⟨w, q, hw⟩ := exists_inv_mod_seven h2s
  have hw' : 2 * (s : L d) * w = 1 + 7 * q := by exact_mod_cast hw
  induction k with
  | zero =>
    -- `(t - s)(t + s) = 7c₀` with `t + s` a unit
    have hu : t + s ∉ O.nonunits := fun h => by
      have h2sm : ((2 * s : ℤ) : L d) ∈ O.nonunits := by
        rw [show ((2 * s : ℤ) : L d) = (t + s) - (t - s) by push_cast; ring]
        exact sub_mem h h2
      exact h2s (seven_dvd_of_intCast_mem_nonunits O h7 h2sm)
    have hne : t + s ≠ 0 := fun h => hu (h ▸ zero_mem _)
    refine ⟨0, c₀ * (t + s)⁻¹, O.mul_mem _ _ (intCast_mem O c₀) (inv_mem_of_notMem_nonunits O hu), ?_⟩
    simp only [Int.cast_zero, mul_zero, add_zero, zero_add, pow_one]
    field_simp
    linear_combination ht + hc₀'
  | succ k ih =>
    obtain ⟨m, γ, hγ, hb⟩ := ih
    -- `N = (s + 7m)² - d = 7ᵏ⁺¹ · (-γ (2(s + 7m) + 7ᵏ⁺¹γ))`
    have hN : (((s + 7 * m) ^ 2 - d : ℤ) : L d)
        = 7 ^ (k + 1) * (-γ * (2 * (s + 7 * m) + 7 ^ (k + 1) * γ)) := by
      push_cast
      linear_combination ht - ((s : L d) + 7 * m + 7 ^ (k + 1) * γ + t) * hb
    have hmem : -γ * (2 * (s + 7 * m) + 7 ^ (k + 1) * γ) ∈ O := by
      have := (-(⟨γ, hγ⟩ : O) * (2 * ((s : O) + 7 * (m : O)) + 7 ^ (k + 1) * (⟨γ, hγ⟩ : O))).2
      push_cast at this
      exact this
    obtain ⟨c, hc⟩ := seven_pow_dvd_of_eq O h7 hmem hN
    have hc' : ((s : L d) + 7 * m) ^ 2 - d = 7 ^ (k + 1) * c := by exact_mod_cast hc
    refine ⟨m - 7 ^ k * c * w, -(γ * q + 2 * γ * m * w + 7 ^ k * γ ^ 2 * w), ?_, ?_⟩
    · have := (-((⟨γ, hγ⟩ : O) * q + 2 * (⟨γ, hγ⟩ : O) * m * w + 7 ^ k * (⟨γ, hγ⟩ : O) ^ 2 * w)).2
      push_cast at this
      exact this
    · push_cast
      linear_combination (1 - (w : L d) * (t + (s + 7 * m) + 7 ^ (k + 1) * γ)) * hb - (w : L d) * hc'
        + (w : L d) * ht - 7 ^ (k + 1) * γ * hw'

/-- **The residue field is `𝔽₇`:** every `x ∈ O` is `≡` an integer modulo `𝔪`. -/
theorem exists_int_sub_mem {x : L d} (hx : x ∈ O) : ∃ N : ℤ, x - N ∈ O.nonunits := by
  obtain ⟨t, ht, ht2⟩ := exists_t O h7 hds
  have htsq : t ^ 2 = d := by
    rcases ht with rfl | rfl
    · exact r_sq d
    · rw [neg_sq]; exact r_sq d
  obtain ⟨D, hD, n₀, n₁, h⟩ := exists_int_combo_d ht x
  obtain ⟨j, D', hD', rfl⟩ := Nat.exists_eq_pow_mul_and_not_dvd hD.ne' 7 (by norm_num)
  obtain ⟨m, γ, hγ, hb⟩ := hensel_seven O h7 hs hds htsq ht2 j
  -- `7ʲ (D'x - 7n₁γ) = n₀ + n₁(s + 7m)`, an integer
  have hmem : (D' : L d) * x - 7 * n₁ * γ ∈ O :=
    sub_mem (O.mul_mem _ _ (natCast_mem O D') hx)
      (O.mul_mem _ _ (O.mul_mem _ _ (ofNat_mem O 7) (intCast_mem O n₁)) hγ)
  have hN : ((n₀ + n₁ * (s + 7 * m) : ℤ) : L d) = 7 ^ j * ((D' : L d) * x - 7 * n₁ * γ) := by
    push_cast at h ⊢
    linear_combination -h - (n₁ : L d) * hb
  obtain ⟨N', hN'⟩ := seven_pow_dvd_of_eq O h7 hmem hN
  have h7ne : (7 : L d) ≠ 0 := by exact_mod_cast (by norm_num : (7 : ℤ) ≠ 0)
  have hx' : (D' : L d) * x - 7 * n₁ * γ = N' := by
    apply mul_left_cancel₀ (pow_ne_zero j h7ne)
    rw [← hN, hN']
    push_cast
    ring
  -- `D'` is prime to 7, hence a unit of `O`, with an inverse `s'` modulo 7
  have hD'7 : ¬ (7 : ℤ) ∣ (D' : ℤ) := by exact_mod_cast hD'
  have hD'0 : (D' : L d) ≠ 0 := by
    intro h0
    apply hD'7
    have : (D' : ℕ) = 0 := by exact_mod_cast h0
    rw [this]
    simp
  have hinv : (D' : L d)⁻¹ ∈ O := inv_mem_of_notMem_nonunits O fun hm =>
    hD'7 (seven_dvd_of_intCast_mem_nonunits O h7 (by exact_mod_cast hm))
  obtain ⟨s', k', hs'⟩ := exists_inv_mod_seven hD'7
  have hs'' : (D' : L d) * s' = 1 + 7 * k' := by exact_mod_cast hs'
  refine ⟨N' * s', ?_⟩
  have e : x - ((N' * s' : ℤ) : L d) = 7 * ((n₁ * γ - N' * k') * (D' : L d)⁻¹) := by
    have hDD : (D' : L d) * (D' : L d)⁻¹ = 1 := mul_inv_cancel₀ hD'0
    push_cast
    linear_combination (D' : L d)⁻¹ * hx' - (N' : L d) * (D' : L d)⁻¹ * hs'' - (x - N' * s') * hDD
  rw [e]
  exact mul_mem_nonunits_right O h7 (O.mul_mem _ _ (sub_mem (O.mul_mem _ _ (intCast_mem O n₁) hγ)
    (O.mul_mem _ _ (intCast_mem O N') (intCast_mem O k'))) hinv)

end Residue

/-! ### The ramified case: `d = 7d'` with `7 ∤ d'` -/

section Ramified

variable {d : ℕ} (O : ValuationSubring (L d)) (h7 : (7 : L d) ∈ O.nonunits) {d' : ℕ} (hd : d = 7 * d')
  (hd' : ¬ 7 ∣ d')
include h7 hd hd'

omit h7 hd' in
lemma natCast_d : (d : L d) = 7 * d' := by
  have := congrArg (Nat.cast : ℕ → L d) hd
  push_cast at this
  exact this

omit hd' in
/-- `√d ∈ 𝔪`, since its square `7d'` is. -/
lemma r_mem_nonunits : r d ∈ O.nonunits :=
  mem_nonunits_of_pow_mem O (r_mem O) (n := 2) (by
    rw [r_sq d, natCast_d hd]
    exact mul_mem_nonunits_right O h7 (natCast_mem O d'))

/-- If `N₀ + N₁√d ∈ 7O` with integers `Nᵢ`, then 7 divides both. -/
lemma seven_dvd_of_add_mul_r {N₀ N₁ : ℤ} {b : L d} (hb : b ∈ O) (h : (N₀ : L d) + N₁ * r d = 7 * b) :
    (7 : ℤ) ∣ N₀ ∧ (7 : ℤ) ∣ N₁ := by
  have hr := r_mem_nonunits O h7 hd
  have hd'0 : (d' : L d) ≠ 0 := by
    have : d' ≠ 0 := fun h0 => hd' (h0 ▸ dvd_zero 7)
    exact_mod_cast this
  have hd'u : (d' : L d)⁻¹ ∈ O := inv_mem_of_notMem_nonunits O fun hm => hd' (by
    have := seven_dvd_of_intCast_mem_nonunits O h7 (n := d') (by exact_mod_cast hm)
    exact_mod_cast this)
  -- `7 = √d · √d / d'`
  have h7r : (7 : L d) = r d * (r d * (d' : L d)⁻¹) := by
    rw [← mul_assoc, ← sq, r_sq d, natCast_d hd, mul_assoc, mul_inv_cancel₀ hd'0, mul_one]
  have hrne : r d ≠ 0 := fun h0 => by
    have h1 : (7 : L d) * d' = 0 := by rw [← natCast_d hd, ← r_sq d, h0]; ring
    exact hd'0 ((mul_eq_zero.mp h1).resolve_left (by norm_num))
  have h0 : (7 : ℤ) ∣ N₀ := seven_dvd_of_intCast_mem_nonunits O h7 (by
    rw [show (N₀ : L d) = r d * (r d * (d' : L d)⁻¹ * b - N₁) by linear_combination h + b * h7r]
    exact mul_mem_nonunits_right O hr
      (sub_mem (O.mul_mem _ _ (O.mul_mem _ _ (r_mem O) hd'u) hb) (intCast_mem O N₁)))
  refine ⟨h0, ?_⟩
  obtain ⟨M₀, hM₀⟩ := h0
  have hM₀' : (N₀ : L d) = 7 * M₀ := by exact_mod_cast hM₀
  exact seven_dvd_of_intCast_mem_nonunits O h7 (by
    rw [show (N₁ : L d) = r d * ((d' : L d)⁻¹ * (b - M₀)) by
      apply mul_left_cancel₀ hrne
      linear_combination h - hM₀' + (b - M₀) * h7r]
    exact mul_mem_nonunits_right O hr (O.mul_mem _ _ hd'u (sub_mem hb (intCast_mem O M₀))))

/-- **The residue field is `𝔽₇`** in the ramified case: every `x ∈ O` is `≡` an integer modulo `𝔪`. -/
theorem exists_int_sub_mem_ramified {x : L d} (hx : x ∈ O) : ∃ N : ℤ, x - N ∈ O.nonunits := by
  have hr := r_mem_nonunits O h7 hd
  obtain ⟨D, hD, n₀, n₁, h⟩ := exists_int_combo_d (d := d) (t := r d) (Or.inl rfl) x
  -- descent on the power of 7 in `D`
  suffices key : ∀ D : ℕ, 0 < D → ∀ N₀ N₁ : ℤ, (D : L d) * x = N₀ + N₁ * r d →
      ∃ N : ℤ, x - N ∈ O.nonunits from key D hD n₀ n₁ h
  intro D
  induction D using Nat.strong_induction_on with
  | _ D ih =>
  intro hD N₀ N₁ h
  by_cases h7D : 7 ∣ D
  · obtain ⟨D'', rfl⟩ := h7D
    have h' : (N₀ : L d) + N₁ * r d = 7 * (D'' * x) := by push_cast at h; linear_combination -h
    obtain ⟨⟨M₀, rfl⟩, ⟨M₁, rfl⟩⟩ :=
      seven_dvd_of_add_mul_r O h7 hd hd' (O.mul_mem _ _ (natCast_mem O D'') hx) h'
    refine ih D'' (by omega) (by omega) M₀ M₁ ?_
    apply mul_left_cancel₀ (show (7 : L d) ≠ 0 by norm_num)
    push_cast at h
    linear_combination h
  · -- `D` is prime to 7, hence a unit, and `x ≡ N₀ / D`
    have hD7 : ¬ (7 : ℤ) ∣ (D : ℤ) := by exact_mod_cast h7D
    have hDne : (D : L d) ≠ 0 := by exact_mod_cast hD.ne'
    have hinv : (D : L d)⁻¹ ∈ O := inv_mem_of_notMem_nonunits O fun hm =>
      hD7 (seven_dvd_of_intCast_mem_nonunits O h7 (by exact_mod_cast hm))
    obtain ⟨s, k, hs⟩ := exists_inv_mod_seven hD7
    have hs' : (D : L d) * s = 1 + 7 * k := by exact_mod_cast hs
    refine ⟨N₀ * s, ?_⟩
    have e : x - ((N₀ * s : ℤ) : L d) = (N₁ * r d - 7 * (N₀ * k)) * (D : L d)⁻¹ := by
      have hDD : (D : L d) * (D : L d)⁻¹ = 1 := mul_inv_cancel₀ hDne
      push_cast
      linear_combination (D : L d)⁻¹ * h - (N₀ : L d) * (D : L d)⁻¹ * hs' - (x - N₀ * s) * hDD
    rw [e]
    exact mul_mem_nonunits_right O (sub_mem (mul_mem_nonunits_left O (intCast_mem O N₁) hr)
      (mul_mem_nonunits_right O h7 (O.mul_mem _ _ (intCast_mem O N₀) (intCast_mem O k)))) hinv

end Ramified

/-- `unitDistGraph K` maps to the graph of `x² + y²` over `K`. -/
def toSumSq (K : IntermediateField ℚ ℝ) : unitDistGraph K →g sumSqGraph K where
  toFun := id
  map_rel' := fun {p q} h => Subtype.ext (by push_cast; exact h)

/-- **Upper bound** (Moorhouse, Lemma 8.2 at 7). If `d ≡ s² (mod 7)` with `7 ∤ s`, the unit-distance
graph of `ℚ(√d)²` is 4-colourable. -/
theorem colorable_four {d : ℕ} {s : ℤ} (hs : ¬ (7 : ℤ) ∣ s) (hds : (7 : ℤ) ∣ (d : ℤ) - s ^ 2) :
    (unitDistGraph (L d)).Colorable 4 := by
  obtain ⟨O, h7⟩ := exists_valuationSubring_seven_mem_nonunits (L d)
  obtain ⟨red, hker⟩ := exists_red_of_residue_seven (R := O)
    (fun x => by
      obtain ⟨N, hN⟩ := exists_int_sub_mem O h7 hs hds x.2
      exact ⟨N, O.coe_mem_nonunits_iff.mp (by push_cast; exact hN)⟩)
    (fun n => ⟨fun h => seven_dvd_of_intCast_mem_nonunits O h7
        (by simpa using O.coe_mem_nonunits_iff.mpr h),
      fun ⟨c, hc⟩ => O.coe_mem_nonunits_iff.mp (by
        rw [hc]
        push_cast
        exact mul_mem_nonunits_right O h7 (intCast_mem O c))⟩)
  exact (sumSqGraph_colorable O red hker).of_hom (toSumSq (L d))

/-- **Upper bound** (Moorhouse, Lemma 8.2 at 7), the ramified case. If `d = 7d'` with `7 ∤ d'`, the
unit-distance graph of `ℚ(√d)²` is 4-colourable. -/
theorem colorable_four_ramified {d d' : ℕ} (hd : d = 7 * d') (hd' : ¬ 7 ∣ d') :
    (unitDistGraph (L d)).Colorable 4 := by
  obtain ⟨O, h7⟩ := exists_valuationSubring_seven_mem_nonunits (L d)
  obtain ⟨red, hker⟩ := exists_red_of_residue_seven (R := O)
    (fun x => by
      obtain ⟨N, hN⟩ := exists_int_sub_mem_ramified O h7 hd hd' x.2
      exact ⟨N, O.coe_mem_nonunits_iff.mp (by push_cast; exact hN)⟩)
    (fun n => ⟨fun h => seven_dvd_of_intCast_mem_nonunits O h7
        (by simpa using O.coe_mem_nonunits_iff.mpr h),
      fun ⟨c, hc⟩ => O.coe_mem_nonunits_iff.mp (by
        rw [hc]
        push_cast
        exact mul_mem_nonunits_right O h7 (intCast_mem O c))⟩)
  exact (sumSqGraph_colorable O red hker).of_hom (toSumSq (L d))

/-! ## Upper bound: reduction at the place over 2 when `d ≡ 3 (mod 8)` -/

/-! ### The graph of `a² - ab + κb²` over a field with a place of residue field `𝔽₂` -/

section TwoForm

variable {F : Type*} [Field F]

/-- The graph on `F × F` of the form `a² - ab + κb²`: `p ~ q` iff the form takes the value 1 at `p - q`. -/
def kGraph (κ : F) : SimpleGraph (F × F) where
  Adj p q := (p.1 - q.1) ^ 2 - (p.1 - q.1) * (p.2 - q.2) + κ * (p.2 - q.2) ^ 2 = 1
  symm := ⟨fun p q h => by linear_combination h⟩
  loopless := ⟨fun p h => by simp at h⟩

variable (O : ValuationSubring F) (red : O →+* ZMod 2) {κ : F} (hκ : κ ∈ O) (hκ1 : red ⟨κ, hκ⟩ = 1)
include hκ1

/-- Over `𝔽₂`, with `κ ≡ 1`, neither `1 - t + κt²` nor `t² - t + κ` vanishes; so neither is the square of an
element of the maximal ideal. -/
lemma kform_ne_sq (hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0) (t x : O) (hx : x ∈ maximalIdeal O) :
    1 - t + ⟨κ, hκ⟩ * t ^ 2 ≠ x ^ 2 ∧ t ^ 2 - t + ⟨κ, hκ⟩ ≠ x ^ 2 := by
  have key : ∀ u : ZMod 2, 1 - u + 1 * u ^ 2 ≠ 0 ^ 2 ∧ u ^ 2 - u + 1 ≠ 0 ^ 2 := by decide
  constructor
  · intro h
    have h' := congrArg red h
    rw [map_add, map_sub, map_mul, map_one, map_pow, map_pow, hred x hx, hκ1] at h'
    exact (key (red t)).1 h'
  · intro h
    have h' := congrArg red h
    rw [map_add, map_sub, map_pow, map_pow, hred x hx, hκ1] at h'
    exact (key (red t)).2 h'

/-- If `a² - ab + κb² = 1`, then `a, b ∈ O`. -/
lemma mem_of_kform_eq_one (hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0) {a b : F}
    (h : a ^ 2 - a * b + κ * b ^ 2 = 1) : a ∈ O ∧ b ∈ O := by
  -- if `a` dominates and `a ∉ O`, then `t = b/a ∈ O` and `1 - t + κt² = a⁻²` with `a⁻¹ ∈ 𝔪`
  have hA : a ≠ 0 → a⁻¹ ∈ O.nonunits → b / a ∈ O → False := by
    intro ha0 hx ht
    obtain ⟨hxO, hxm⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp hx
    apply (kform_ne_sq O red hκ hκ1 hred ⟨b / a, ht⟩ ⟨a⁻¹, hxO⟩ hxm).1
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h
  -- if `b` dominates and `b ∉ O`, then `s = a/b ∈ O` and `s² - s + κ = b⁻²` with `b⁻¹ ∈ 𝔪`
  have hB : b ≠ 0 → b⁻¹ ∈ O.nonunits → a / b ∈ O → False := by
    intro hb0 hx hs
    obtain ⟨hxO, hxm⟩ := O.mem_nonunits_iff_exists_mem_maximalIdeal.mp hx
    apply (kform_ne_sq O red hκ hκ1 hred ⟨a / b, hs⟩ ⟨b⁻¹, hxO⟩ hxm).2
    apply Subtype.ext
    push_cast
    field_simp
    linear_combination h
  have hout : ∀ x : F, x ∉ O → x ≠ 0 ∧ x⁻¹ ∈ O.nonunits := fun x hx =>
    ⟨fun h0 => hx (h0 ▸ O.zero_mem), O.inv_mem_nonunits_iff.mpr (Or.inr hx)⟩
  constructor
  · by_contra ha
    obtain ⟨ha0, hai⟩ := hout a ha
    by_cases ht : b / a ∈ O
    · exact hA ha0 hai ht
    · have hs : a / b ∈ O := by
        have := (O.mem_or_inv_mem (b / a)).resolve_left ht
        rwa [inv_div] at this
      have hb0 : b ≠ 0 := by
        rintro rfl
        exact ht (by simp)
      refine hB hb0 ?_ hs
      rw [show b⁻¹ = a / b * a⁻¹ by field_simp]
      exact mul_mem_nonunits_left O hs hai
  · by_contra hb
    obtain ⟨hb0, hbi⟩ := hout b hb
    by_cases hs : a / b ∈ O
    · exact hB hb0 hbi hs
    · have ht : b / a ∈ O := by
        have := (O.mem_or_inv_mem (a / b)).resolve_left hs
        rwa [inv_div] at this
      have ha0 : a ≠ 0 := by
        rintro rfl
        exact hs (by simp)
      refine hA ha0 ?_ ht
      rw [show a⁻¹ = b / a * b⁻¹ by field_simp]
      exact mul_mem_nonunits_left O ht hbi

/-- **Core lemma at 2.** Let `red : O →+* ZMod 2` have kernel the maximal ideal, and `κ ≡ 1`. If
`a² - ab + κb² = 1`, then `a, b ∈ O` and `(red a, red b) ≠ (0, 0)`. -/
theorem core_lemma_two (hker : RingHom.ker red = maximalIdeal O) {a b : F}
    (h : a ^ 2 - a * b + κ * b ^ 2 = 1) :
    ∃ (ha : a ∈ O) (hb : b ∈ O), (red ⟨a, ha⟩, red ⟨b, hb⟩) ≠ (0, 0) := by
  have hred : ∀ x : O, x ∈ maximalIdeal O → red x = 0 := fun x hx => by
    rw [← hker] at hx
    exact hx
  obtain ⟨ha, hb⟩ := mem_of_kform_eq_one O red hκ hκ1 hred h
  refine ⟨ha, hb, ?_⟩
  intro h0
  simp only [Prod.mk.injEq] at h0
  have e : (⟨a, ha⟩ : O) ^ 2 - ⟨a, ha⟩ * ⟨b, hb⟩ + ⟨κ, hκ⟩ * ⟨b, hb⟩ ^ 2 = 1 :=
    Subtype.ext (by simpa using h)
  have e' := congrArg red e
  rw [map_add, map_sub, map_mul, map_mul, map_pow, map_pow, h0.1, h0.2, hκ1, map_one] at e'
  exact absurd e' (by decide)

lemma colour_adj_two (hker : RingHom.ker red = maximalIdeal O) {p q : F × F}
    (hpq : (kGraph κ).Adj p q) : colour O red p ≠ colour O red q := by
  obtain ⟨ha, hb, hne⟩ := core_lemma_two O red hκ hκ1 hker hpq
  intro hc
  apply hne
  simp only [colour, Prod.mk.injEq] at hc
  have e1 : p.1 - rep O p.1 = (p.1 - q.1) + (q.1 - rep O q.1) := by
    rw [rep_eq_of_sub_mem O ha]; ring
  have e2 : p.2 - rep O p.2 = (p.2 - q.2) + (q.2 - rep O q.2) := by
    rw [rep_eq_of_sub_mem O hb]; ring
  rw [e1, redExt_add O red ha (sub_rep_mem O _)] at hc
  rw [e2, redExt_add O red hb (sub_rep_mem O _)] at hc
  have h1 : redExt O red (p.1 - q.1) = 0 := by linear_combination hc.1
  have h2 : redExt O red (p.2 - q.2) = 0 := by linear_combination hc.2
  simp only [redExt, ha, hb, dite_true] at h1 h2
  rw [h1, h2]

/-- The graph of `a² - ab + κb²` with `κ ≡ 1` is 4-colourable, with colours `𝔽₂²`. -/
theorem kGraph_colorable (hker : RingHom.ker red = maximalIdeal O) : (kGraph κ).Colorable 4 := by
  have hc : (kGraph κ).Colorable (Fintype.card (ZMod 2 × ZMod 2)) :=
    (SimpleGraph.Coloring.mk (colour O red) (fun h => colour_adj_two O red hκ hκ1 hker h)).colorable
  simpa using hc

end TwoForm

/-! ### The residue field of the place of `ℚ(√d)` over 2 is `𝔽₂` when `d ≡ 3 (mod 4)` -/

section Two

variable {d : ℕ} (O : ValuationSubring (L d)) (h2 : (2 : L d) ∈ O.nonunits) (hd : d % 4 = 3)
include h2 hd

/-- `π = √d - 1 ∈ 𝔪`: its square is `2((d + 1)/2 - √d)`. -/
lemma r_sub_one_mem_nonunits : r d - 1 ∈ O.nonunits := by
  refine mem_nonunits_of_pow_mem O (sub_mem (r_mem O) O.one_mem) (n := 2) ?_
  obtain ⟨c, hc⟩ : ∃ c : ℕ, d + 1 = 2 * c := ⟨(d + 1) / 2, by omega⟩
  have hc' : (d : L d) + 1 = 2 * c := by exact_mod_cast hc
  rw [show (r d - 1) ^ 2 = 2 * ((c : L d) - r d) by linear_combination r_sq d + hc']
  exact mul_mem_nonunits_right O h2 (sub_mem (natCast_mem O c) (r_mem O))

/-- If `N₀ + N₁π ∈ 2O` with integers `Nᵢ`, then both are even. -/
lemma even_of_add_mul_pi {N₀ N₁ : ℤ} {b : L d} (hb : b ∈ O) (h : (N₀ : L d) + N₁ * (r d - 1) = 2 * b) :
    Even N₀ ∧ Even N₁ := by
  have hπ := r_sub_one_mem_nonunits O h2 hd
  have e0 : Even N₀ := even_of_intCast_mem_nonunits O h2 (by
    rw [show (N₀ : L d) = 2 * b - N₁ * (r d - 1) by linear_combination h]
    exact sub_mem (mul_mem_nonunits_right O h2 hb) (mul_mem_nonunits_left O (intCast_mem O N₁) hπ))
  obtain ⟨M₀, hM₀⟩ := e0
  refine ⟨⟨M₀, hM₀⟩, ?_⟩
  -- `π(√d + 1) = d - 1 = 2(2k + 1)`, so `(2k + 1)N₁ = (√d + 1)(b - M₀) ∈ 𝔪`
  obtain ⟨k, hk⟩ : ∃ k : ℕ, d = 4 * k + 3 := ⟨d / 4, by omega⟩
  have hk' : (d : L d) = 4 * k + 3 := by exact_mod_cast hk
  have hM₀' : (N₀ : L d) = M₀ + M₀ := by exact_mod_cast hM₀
  have hm : ((N₁ * (2 * k + 1) : ℤ) : L d) ∈ O.nonunits := by
    rw [show ((N₁ * (2 * k + 1) : ℤ) : L d) = (b - M₀) * (r d + 1) by
      apply mul_left_cancel₀ (two_ne_zero (α := L d))
      push_cast
      linear_combination (r d + 1) * h - (r d + 1) * hM₀' - (N₁ : L d) * r_sq d - (N₁ : L d) * hk']
    refine mul_mem_nonunits_left O (sub_mem hb (intCast_mem O M₀)) ?_
    rw [show r d + 1 = (r d - 1) + 2 by ring]
    exact add_mem hπ h2
  rcases Int.even_mul.mp (even_of_intCast_mem_nonunits O h2 hm) with h1 | ⟨j, hj⟩
  · exact h1
  · omega

/-- The residue computation: if `y ∈ O` and `D y = N₀ + N₁π` with `D > 0`, then `y ≡ 0` or `y ≡ 1`. -/
lemma residue_of_add_mul_pi {y : L d} (hy : y ∈ O) : ∀ D : ℕ, 0 < D → ∀ N₀ N₁ : ℤ,
    (D : L d) * y = N₀ + N₁ * (r d - 1) → y ∈ O.nonunits ∨ y - 1 ∈ O.nonunits := by
  have hπ := r_sub_one_mem_nonunits O h2 hd
  intro D
  induction D using Nat.strong_induction_on with
  | _ D ih =>
  intro hD N₀ N₁ h
  rcases Nat.even_or_odd D with ⟨D', rfl⟩ | ⟨D', rfl⟩
  · obtain ⟨⟨M₀, rfl⟩, ⟨M₁, rfl⟩⟩ := even_of_add_mul_pi O h2 hd (N₀ := N₀) (N₁ := N₁)
      (b := D' * y) (O.mul_mem _ _ (natCast_mem O D') hy) (by push_cast at h; linear_combination -h)
    refine ih D' (by omega) (by omega) M₀ M₁ ?_
    apply mul_left_cancel₀ (two_ne_zero (α := L d))
    push_cast at h
    linear_combination h
  · have hyn : y - N₀ ∈ O.nonunits := by
      rw [show y - N₀ = (r d - 1) * N₁ - ((D' : L d) * y) * 2 by push_cast at h; linear_combination h]
      exact sub_mem (mul_mem_nonunits_right O hπ (intCast_mem O N₁))
        (mul_mem_nonunits_left O (O.mul_mem _ _ (natCast_mem O D') hy) h2)
    rcases Int.even_or_odd' N₀ with ⟨k, hk | hk⟩
    · left
      rw [show y = (y - N₀) + (k : L d) * 2 by rw [hk]; push_cast; ring]
      exact add_mem hyn (mul_mem_nonunits_left O (intCast_mem O k) h2)
    · right
      rw [show y - 1 = (y - N₀) + (k : L d) * 2 by rw [hk]; push_cast; ring]
      exact add_mem hyn (mul_mem_nonunits_left O (intCast_mem O k) h2)

/-- **The residue field is `𝔽₂`**: every `x ∈ O` is `≡ 0` or `≡ 1` modulo `𝔪`. -/
theorem mem_nonunits_or_sub_one_mem_two {x : L d} (hx : x ∈ O) :
    x ∈ O.nonunits ∨ x - 1 ∈ O.nonunits := by
  obtain ⟨D, hD, n₀, n₁, h⟩ := exists_int_combo_d (d := d) (t := r d) (Or.inl rfl) x
  exact residue_of_add_mul_pi O h2 hd hx D hD (n₀ + n₁) n₁ (by rw [h]; push_cast; ring)

end Two

/-- The coordinates `a = x + y/√d`, `b = 2y/√d` (note `1/√d = √d/d`). -/
noncomputable def toK (d : ℕ) (p : L d × L d) : L d × L d :=
  (p.1 + p.2 * r d / d, 2 * p.2 * r d / d)

/-- If `d + 1 = 4κ`, then in these coordinates `x² + y² = a² - ab + κb²`, so the change of coordinates maps
the unit-distance graph to the graph of that form. -/
noncomputable def toKHom {d κ : ℕ} (hdκ : d + 1 = 4 * κ) : unitDistGraph (L d) →g kGraph (κ : L d) where
  toFun := toK d
  map_rel' := by
    intro p q h
    have h' : (p.1 - q.1) ^ 2 + (p.2 - q.2) ^ 2 = (1 : L d) := Subtype.ext (by push_cast; exact h)
    have hκ' : (d : L d) + 1 = 4 * κ := by exact_mod_cast hdκ
    have hdd : (d : L d) * (d : L d)⁻¹ = 1 := mul_inv_cancel₀ (by
      have : d ≠ 0 := by omega
      exact_mod_cast this)
    change ((p.1 + p.2 * r d / d) - (q.1 + q.2 * r d / d)) ^ 2
      - ((p.1 + p.2 * r d / d) - (q.1 + q.2 * r d / d)) * (2 * p.2 * r d / d - 2 * q.2 * r d / d)
      + (κ : L d) * (2 * p.2 * r d / d - 2 * q.2 * r d / d) ^ 2 = 1
    linear_combination h' + ((p.2 - q.2) ^ 2 * (d : L d)⁻¹ ^ 2 * d) * r_sq d
      - ((p.2 - q.2) ^ 2 * r d ^ 2 * (d : L d)⁻¹ ^ 2) * hκ' + (p.2 - q.2) ^ 2 * ((d : L d)⁻¹ * d + 1) * hdd

/-- **Upper bound** (K. G. Fischer, 1990, Theorem 10: reduction at 2). If `d ≡ 3 (mod 8)`, the unit-distance
graph of `ℚ(√d)²` is 4-colourable. -/
theorem colorable_four_two {d : ℕ} (hd : d % 8 = 3) : (unitDistGraph (L d)).Colorable 4 := by
  obtain ⟨O, h2⟩ := exists_valuationSubring_two_mem_nonunits (L d)
  have hB : ∀ x : O, x ∈ maximalIdeal O ∨ x - 1 ∈ maximalIdeal O := by
    intro x
    rcases mem_nonunits_or_sub_one_mem_two O h2 (by omega) x.2 with h | h
    · exact Or.inl (O.coe_mem_nonunits_iff.mp h)
    · exact Or.inr (O.coe_mem_nonunits_iff.mp (by simpa using h))
  obtain ⟨red, hker⟩ := exists_red_of_residue_two (O.coe_mem_nonunits_iff.mp h2) hB
  obtain ⟨k, hk⟩ : ∃ k : ℕ, d + 1 = 4 * (2 * k + 1) := ⟨d / 8, by omega⟩
  have hκ : ((2 * k + 1 : ℕ) : L d) ∈ O := natCast_mem O _
  have hκ1 : red ⟨((2 * k + 1 : ℕ) : L d), hκ⟩ = 1 := by
    have e : (⟨((2 * k + 1 : ℕ) : L d), hκ⟩ : O) = ((2 * k + 1 : ℕ) : O) := Subtype.ext (by norm_cast)
    rw [e, map_natCast, ← ZMod.natCast_mod (2 * k + 1) 2, show (2 * k + 1) % 2 = 1 by omega, Nat.cast_one]
  exact (kGraph_colorable O red hκ hκ1 hker).of_hom (toKHom hk)

end QuadraticPlanes
