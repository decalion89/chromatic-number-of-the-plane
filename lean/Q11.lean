import LocalColouring

/-!
# The plane over `ℚ(√11)` has chromatic number 4

**Theorem** (`Q11.chromaticNumber_eq_four`). Let `L = ℚ(√11) ⊆ ℝ`. The unit-distance graph on `L²`, in which
`p ~ q` iff `(p₁ - q₁)² + (p₂ - q₂)² = 1`, has chromatic number 4.

`ℚ(√11)` is the first real quadratic field whose plane needs four colours: `notes/quadratic_planes.md`.

*Lower bound* (`not_colorable_three`). The graph of `data/quadratic_planes/q11.json` has 76 vertices
`[a, b, c, e]`, standing for `((a + b√11)/30, (c + e√11)/30)`, and 172 edges.
1. Each edge is a unit distance: the differences satisfy `Δa² + 11Δb² + Δc² + 11Δe² = 900` and
   `ΔaΔb + ΔcΔe = 0` (`adj_pt`); the kernel checks these integer identities for every edge (`checkEdges_E`).
2. The graph has no 3-colouring. The formula `data/quadratic_planes/q11.cnf` has the variable `3v + c + 1` for
   "vertex `v` has colour `c`", a clause per vertex, three per edge, and two unit clauses that give the edge
   `0-1` the colours 0 and 1. Mathlib's `lrat_proof` checks in the kernel that the LRAT proof
   `q11.lrat` refutes it, and states the result as a propositional theorem (`q11Refuted`). A
   3-colouring, with its colours renamed so that `0` has colour 0 and `1` colour 1, would make every clause
   true.

*Upper bound* (`colorable_four`; Moorhouse, Lemma 8.2 at 7). `11 ≡ 2² (mod 7)` and `-1` is not a square
modulo 7.
1. Chevalley's extension theorem gives a valuation subring `O` of `L` with `7` in its maximal ideal `𝔪`.
2. Its residue field is `𝔽₇` (`exists_int_sub_mem`): `√11 ∈ O`, and one of `±√11 - 2` lies in `𝔪`, say `t - 2`.
   Hensel's lemma gives integers `bₖ` with `t - bₖ ∈ 7ᵏ⁺¹O`, and writing `x = (n₀ + n₁t)/(7ʲD')` with `7 ∤ D'`
   shows that every `x ∈ O` is `≡` an integer modulo `𝔪`.
3. If `x² + y² = 1`, then `x, y ∈ O` and their residues satisfy `x̄² + ȳ² = 1` in `𝔽₇` (`core_lemma_seven`).
   Colour `z ∈ L²` by a proper 4-colouring of the unit-distance graph of `𝔽₇²`
   (`data/quadratic_planes/finite_planes.json`) at the residues of `z - rep z`, where `rep z` is a fixed
   representative of `z + O²` (`sumSqGraph_colorable`).
-/

namespace Q11

open LocalColouring IsLocalRing

/-- `L = ℚ(√11)`, as a subfield of `ℝ`. -/
noncomputable def L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√11}

lemma sqrt11_mem : √11 ∈ L := IntermediateField.subset_adjoin ℚ _ (by simp)

/-- `√11` as an element of `L`. -/
noncomputable def r11 : L := ⟨√11, sqrt11_mem⟩

@[simp] lemma coe_r11 : (r11 : ℝ) = √11 := rfl

lemma sqrt11_sq : √11 ^ 2 = (11 : ℝ) := Real.sq_sqrt (by norm_num)

lemma r11_sq : r11 ^ 2 = 11 := Subtype.ext (by push_cast; exact sqrt11_sq)

/-! ## Lower bound: a graph with 76 vertices -/

/-- The point `[a, b, c, e]`, that is `((a + b√11)/30, (c + e√11)/30)`. -/
noncomputable def pt (v : ℤ × ℤ × ℤ × ℤ) : L × L :=
  ((v.1 + v.2.1 * r11) / 30, (v.2.2.1 + v.2.2.2 * r11) / 30)

/-- The two integer identities for a unit distance between `[a, b, c, e]/30` and `[a', b', c', e']/30`,
as a Boolean (products written out, so that the kernel evaluates it quickly). -/
def unitPairB (v w : ℤ × ℤ × ℤ × ℤ) : Bool :=
  let da := v.1 - w.1
  let db := v.2.1 - w.2.1
  let dc := v.2.2.1 - w.2.2.1
  let de := v.2.2.2 - w.2.2.2
  (da * da + 11 * (db * db) + dc * dc + 11 * (de * de) == 900) && (da * db + dc * de == 0)

lemma unitPairB_spec {v w : ℤ × ℤ × ℤ × ℤ} (h : unitPairB v w = true) :
    (v.1 - w.1) * (v.1 - w.1) + 11 * ((v.2.1 - w.2.1) * (v.2.1 - w.2.1))
        + (v.2.2.1 - w.2.2.1) * (v.2.2.1 - w.2.2.1) + 11 * ((v.2.2.2 - w.2.2.2) * (v.2.2.2 - w.2.2.2)) = 900 ∧
      (v.1 - w.1) * (v.2.1 - w.2.1) + (v.2.2.1 - w.2.2.1) * (v.2.2.2 - w.2.2.2) = 0 := by
  simpa [unitPairB] using h

/-- A unit distance in `L²` is a unit distance in `ℝ²`. -/
lemma adj_of_sq_eq_one {p q : L × L} (h : (p.1 - q.1) ^ 2 + (p.2 - q.2) ^ 2 = 1) :
    (unitDistGraph L).Adj p q := by
  show ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1
  have := congrArg (fun x : L => (x : ℝ)) h
  push_cast at this
  exact this

/-- Since `r11² = 11`, the two identities say that the distance is 1. -/
lemma adj_pt {v w : ℤ × ℤ × ℤ × ℤ} (h : unitPairB v w = true) :
    (unitDistGraph L).Adj (pt v) (pt w) := by
  obtain ⟨h1, h2⟩ := unitPairB_spec h
  obtain ⟨a, b, c, e⟩ := v
  obtain ⟨a', b', c', e'⟩ := w
  dsimp only at h1 h2
  have h1' : ((a : L) - a') * ((a : L) - a') + 11 * (((b : L) - b') * ((b : L) - b'))
      + ((c : L) - c') * ((c : L) - c') + 11 * (((e : L) - e') * ((e : L) - e')) = 900 := by
    exact_mod_cast h1
  have h2' : ((a : L) - a') * ((b : L) - b') + ((c : L) - c') * ((e : L) - e') = 0 := by
    exact_mod_cast h2
  apply adj_of_sq_eq_one
  simp only [pt]
  linear_combination (1 / 900) * h1' + (2 * r11 / 900) * h2'
    + ((((b : L) - b') ^ 2 + ((e : L) - e') ^ 2) / 900) * r11_sq

/-- The 76 points of `data/quadratic_planes/q11.json`, as `[a, b, c, e]` with denominator 30. -/
def P : Fin 76 → ℤ × ℤ × ℤ × ℤ := ![
    (-19, 9, -3, -7), (0, 9, -3, 0), (18, 3, 9, -6), (-3, 0, 0, -9),
    (-19, 0, 0, -7), (-39, 3, 15, -3), (-9, 6, 18, 3), (-17, -8, -4, -11),
    (0, 0, 0, 0), (-21, -3, -9, -3), (-20, 3, 15, 4), (-37, -3, -9, -1),
    (9, -6, -18, -3), (2, -8, -4, -4), (-9, -3, 21, 3), (-18, -3, -9, 6),
    (0, -9, 3, 0), (-55, 0, 0, 5), (-39, 0, 0, 3), (-39, 6, -18, -3),
    (-36, 0, 0, 12), (-35, 5, -5, -5), (-33, 0, 0, -9), (-30, 0, 0, 0),
    (-30, 9, -3, 0), (-27, 0, 0, -9), (-26, -1, 13, -2), (-26, 1, -13, -2),
    (-25, 0, 0, 5), (-21, -4, -2, -3), (-21, 3, 9, -3), (-21, 4, 2, -3),
    (-20, 6, -18, 4), (-18, 3, 9, 6), (-17, 8, 4, -11), (-14, 0, 0, -2),
    (-12, -3, -9, -6), (-12, 3, 9, -6), (-11, -2, -22, 1), (-11, 0, 0, 7),
    (-11, 9, -3, 7), (-9, -6, -18, 3), (-9, -3, -15, -3), (-9, -3, -9, 3),
    (-9, 1, -7, 3), (-9, 3, -21, 3), (-9, 3, 9, 3), (-9, 3, 15, -3),
    (-9, 6, -18, -3), (-8, -2, -22, -8), (-2, -8, -4, 4), (-1, -1, 13, -7),
    (-1, 1, -13, -7), (0, -16, -8, 0), (0, -9, -27, 0), (0, -5, -25, 0),
    (0, 0, -30, 0), (0, 0, 30, 0), (0, 9, 27, 0), (5, 0, 0, 5),
    (6, 0, 0, -12), (8, -2, -22, 8), (9, -6, 12, -3), (9, -3, -15, 3),
    (9, -3, -9, -3), (9, 1, -7, -3), (9, 3, 9, -3), (9, 6, -18, 3),
    (9, 6, -12, -3), (9, 6, 18, -3), (9, 9, -9, -9), (11, -2, -22, -1),
    (11, 0, 0, -7), (18, -3, -9, -6), (25, 0, 0, -5), (30, 0, 0, 0)]

/-- The 172 edges of `data/quadratic_planes/q11.json`: all the unit distances among the points. -/
def E : List (Fin 76 × Fin 76) := [
    (0, 1), (0, 4), (0, 19), (1, 8), (1, 24), (1, 32), (1, 45), (1, 47), (1, 58), (2, 8),
    (2, 37), (2, 64), (2, 70), (3, 8), (3, 9), (3, 22), (3, 30), (4, 5), (4, 7), (4, 8),
    (4, 11), (4, 34), (4, 60), (4, 72), (5, 10), (5, 24), (5, 47), (6, 8), (6, 14), (6, 44),
    (6, 66), (7, 13), (7, 29), (7, 49), (8, 10), (8, 12), (8, 13), (8, 15), (8, 16), (8, 23),
    (8, 26), (8, 27), (8, 28), (8, 33), (8, 41), (8, 48), (8, 50), (8, 55), (8, 56), (8, 57),
    (8, 67), (8, 69), (8, 73), (8, 74), (8, 75), (9, 15), (9, 18), (9, 37), (9, 64), (10, 40),
    (11, 15), (11, 17), (11, 21), (11, 36), (11, 38), (12, 43), (12, 62), (13, 53), (13, 71), (14, 16),
    (14, 43), (14, 62), (15, 20), (15, 46), (15, 61), (16, 42), (16, 54), (16, 63), (17, 20), (17, 23),
    (17, 28), (18, 20), (18, 30), (18, 35), (19, 23), (19, 32), (19, 48), (20, 33), (20, 39), (21, 26),
    (21, 34), (21, 48), (22, 23), (22, 35), (23, 24), (23, 25), (23, 36), (23, 37), (23, 39), (24, 40),
    (25, 51), (25, 52), (26, 51), (27, 49), (27, 52), (28, 29), (28, 31), (28, 59), (29, 51), (30, 33),
    (30, 36), (30, 66), (31, 34), (31, 52), (32, 39), (33, 43), (34, 70), (35, 39), (35, 59), (35, 72),
    (36, 43), (36, 60), (36, 73), (37, 46), (37, 60), (38, 48), (38, 49), (38, 50), (38, 61), (38, 65),
    (39, 40), (39, 63), (40, 44), (40, 67), (41, 45), (41, 64), (42, 48), (42, 72), (43, 54), (44, 55),
    (44, 71), (45, 46), (45, 68), (46, 58), (46, 69), (47, 72), (48, 70), (49, 71), (49, 73), (50, 53),
    (51, 74), (52, 74), (53, 54), (54, 56), (54, 64), (55, 65), (56, 68), (57, 58), (57, 62), (58, 66),
    (59, 75), (60, 74), (61, 71), (62, 65), (62, 66), (63, 67), (64, 68), (65, 69), (66, 73), (67, 71),
    (68, 69), (72, 75)]

/-- Every listed edge satisfies the two identities. -/
def checkEdges : List (Fin 76 × Fin 76) → Bool
  | [] => true
  | e :: es => unitPairB (P e.1) (P e.2) && checkEdges es

/-- No listed edge is a loop. -/
def noLoops : List (Fin 76 × Fin 76) → Bool
  | [] => true
  | e :: es => (e.1.val != e.2.val) && noLoops es

lemma checkEdges_E : checkEdges E = true := by decide +kernel

lemma noLoops_E : noLoops E = true := by decide +kernel

lemma checkEdges_mem : ∀ {l : List (Fin 76 × Fin 76)}, checkEdges l = true →
    ∀ e ∈ l, unitPairB (P e.1) (P e.2) = true
  | [], _, _, he => by simp at he
  | _ :: _, h, e, he => by
    simp only [checkEdges, Bool.and_eq_true] at h
    rcases List.mem_cons.mp he with rfl | he
    · exact h.1
    · exact checkEdges_mem h.2 e he

lemma noLoops_mem : ∀ {l : List (Fin 76 × Fin 76)}, noLoops l = true → ∀ e ∈ l, e.1 ≠ e.2
  | [], _, _, he => by simp at he
  | _ :: _, h, e, he => by
    simp only [noLoops, Bool.and_eq_true, bne_iff_ne, ne_eq] at h
    rcases List.mem_cons.mp he with rfl | he
    · exact fun h' => h.1 (congrArg Fin.val h')
    · exact noLoops_mem h.2 e he

lemma pt_adj : ∀ e ∈ E, (unitDistGraph L).Adj (pt (P e.1)) (pt (P e.2)) :=
  fun e he => adj_pt (checkEdges_mem checkEdges_E e he)

lemma no_loop (i : Fin 76) : (i, i) ∉ E := fun h => noLoops_mem noLoops_E _ h rfl

lemma adj_of_mem {i j : Fin 76} (h : (i, j) ∈ E) : (edgeGraph E).Adj i j := by
  rw [edgeGraph, SimpleGraph.fromRel_adj]
  refine ⟨fun hij => ?_, Or.inl h⟩
  subst hij
  exact no_loop i h

-- `q11.cnf` is unsatisfiable: for all propositions `x₀, …, x_227`, some clause is false. The statement is
-- a disjunction, over the clauses, of the negations of the clauses.
lrat_proof q11Refuted
  (include_str "../data/quadratic_planes/q11.cnf")
  (include_str "../data/quadratic_planes/q11.lrat")

/-- Is `(a, b)` a listed edge? -/
def memEn (a b : ℕ) : List (Fin 76 × Fin 76) → Bool
  | [] => false
  | e :: es => (Nat.beq e.1.val a && Nat.beq e.2.val b) || memEn a b es

lemma memEn_spec {a b : ℕ} : ∀ {l : List (Fin 76 × Fin 76)}, memEn a b l = true →
    ∃ e ∈ l, e.1.val = a ∧ e.2.val = b
  | [], h => by simp [memEn] at h
  | e :: _, h => by
    simp only [memEn, Bool.or_eq_true, Bool.and_eq_true] at h
    rcases h with ⟨h1, h2⟩ | h
    · exact ⟨e, List.mem_cons_self .., Nat.eq_of_beq_eq_true h1, Nat.eq_of_beq_eq_true h2⟩
    · obtain ⟨e', he', h'⟩ := memEn_spec h
      exact ⟨e', List.mem_cons_of_mem _ he', h'⟩

/-- The edge clauses `¬x ∨ ¬y` of `q11.cnf`: the same colour `x % 3 = y % 3` at the ends of an edge. -/
def goodEdge (x y : ℕ) : Bool :=
  Nat.beq (x % 3) (y % 3) && (memEn (x / 3) (y / 3) E || memEn (y / 3) (x / 3) E)

/-- Renaming the colours so that `a ↦ 0` and `b ↦ 1`. -/
def rename (a b x : Fin 3) : Fin 3 := if x = a then 0 else if x = b then 1 else 2

lemma rename_ne : ∀ a b x y : Fin 3, a ≠ b → x ≠ y → rename a b x ≠ rename a b y := by decide

/-- The graph of `q11.json` is not 3-colourable. -/
theorem graph_not_colorable : ¬ (edgeGraph E).Colorable 3 := by
  rintro ⟨C⟩
  have hab : C 0 ≠ C 1 := C.valid (adj_of_mem (List.mem_cons_self ..))
  let g : Fin 76 → Fin 3 := fun i => rename (C 0) (C 1) (C i)
  have hg0 : g 0 = 0 := by simp [g, rename]
  have hg1 : g 1 = 1 := by simp [g, rename, Ne.symm hab]
  have hgE : ∀ i j : Fin 76, (i, j) ∈ E → g i ≠ g j := fun i j h =>
    rename_ne _ _ _ _ hab (C.valid (adj_of_mem h))
  -- the valuation: variable `x` (DIMACS `x + 1`) says "vertex `x / 3` has colour `x % 3`"
  let V : ℕ → Prop := fun x =>
    if h : x / 3 < 76 then ((g ⟨x / 3, h⟩ : Fin 3) : ℕ) = x % 3 else False
  have hv : ∀ x (h : x / 3 < 76), V x ↔ ((g ⟨x / 3, h⟩ : Fin 3) : ℕ) = x % 3 := fun x h => by
    simp only [V, h, dite_true]
  -- every vertex clause holds
  have hVert : ∀ x, (¬V x ∧ ¬V (x + 1) ∧ ¬V (x + 2)) → x % 3 = 0 → x / 3 < 76 → False := by
    intro x ⟨n0, n1, n2⟩ hx hlt
    rw [hv x hlt] at n0
    rw [hv (x + 1) (by omega)] at n1
    rw [hv (x + 2) (by omega)] at n2
    have h1 : (x + 1) / 3 = x / 3 := by omega
    have h2 : (x + 2) / 3 = x / 3 := by omega
    simp only [h1, h2] at n1 n2
    have : ((g ⟨x / 3, hlt⟩ : Fin 3) : ℕ) < 3 := (g ⟨x / 3, hlt⟩).2
    omega
  -- every edge clause holds
  have hEdge1 : ∀ x y, V x → V y → Nat.beq (x % 3) (y % 3) = true → memEn (x / 3) (y / 3) E = true →
      False := by
    intro x y vx vy hxy hE
    have hxy' := Nat.eq_of_beq_eq_true hxy
    obtain ⟨⟨i, j⟩, he, hi, hj⟩ := memEn_spec hE
    simp only at hi hj
    have hx : x / 3 < 76 := hi ▸ i.2
    have hy : y / 3 < 76 := hj ▸ j.2
    rw [hv x hx] at vx
    rw [hv y hy] at vy
    have e1 : (⟨x / 3, hx⟩ : Fin 76) = i := Fin.ext hi.symm
    have e2 : (⟨y / 3, hy⟩ : Fin 76) = j := Fin.ext hj.symm
    rw [e1] at vx
    rw [e2] at vy
    exact hgE i j he (Fin.ext (by omega))
  have hEdge : ∀ x y, (V x ∧ V y) → goodEdge x y = true → False := by
    intro x y ⟨vx, vy⟩ hc
    simp only [goodEdge, Bool.and_eq_true, Bool.or_eq_true] at hc
    obtain ⟨hxy, hE | hE⟩ := hc
    · exact hEdge1 x y vx vy hxy hE
    · exact hEdge1 y x vy vx (by
        have := Nat.eq_of_beq_eq_true hxy
        rw [this]
        exact Nat.beq_refl _) hE
  -- and the two unit clauses
  have hVu : V 0 := (hv 0 (by norm_num)).2 (by simpa using congrArg Fin.val hg0)
  have hVv : V 4 := (hv 4 (by norm_num)).2 (by simpa using congrArg Fin.val hg1)
  have H := q11Refuted
    (V 0) (V 1) (V 2) (V 3) (V 4) (V 5) (V 6) (V 7) (V 8) (V 9) (V 10) (V 11)
    (V 12) (V 13) (V 14) (V 15) (V 16) (V 17) (V 18) (V 19) (V 20) (V 21) (V 22) (V 23)
    (V 24) (V 25) (V 26) (V 27) (V 28) (V 29) (V 30) (V 31) (V 32) (V 33) (V 34) (V 35)
    (V 36) (V 37) (V 38) (V 39) (V 40) (V 41) (V 42) (V 43) (V 44) (V 45) (V 46) (V 47)
    (V 48) (V 49) (V 50) (V 51) (V 52) (V 53) (V 54) (V 55) (V 56) (V 57) (V 58) (V 59)
    (V 60) (V 61) (V 62) (V 63) (V 64) (V 65) (V 66) (V 67) (V 68) (V 69) (V 70) (V 71)
    (V 72) (V 73) (V 74) (V 75) (V 76) (V 77) (V 78) (V 79) (V 80) (V 81) (V 82) (V 83)
    (V 84) (V 85) (V 86) (V 87) (V 88) (V 89) (V 90) (V 91) (V 92) (V 93) (V 94) (V 95)
    (V 96) (V 97) (V 98) (V 99) (V 100) (V 101) (V 102) (V 103) (V 104) (V 105) (V 106) (V 107)
    (V 108) (V 109) (V 110) (V 111) (V 112) (V 113) (V 114) (V 115) (V 116) (V 117) (V 118) (V 119)
    (V 120) (V 121) (V 122) (V 123) (V 124) (V 125) (V 126) (V 127) (V 128) (V 129) (V 130) (V 131)
    (V 132) (V 133) (V 134) (V 135) (V 136) (V 137) (V 138) (V 139) (V 140) (V 141) (V 142) (V 143)
    (V 144) (V 145) (V 146) (V 147) (V 148) (V 149) (V 150) (V 151) (V 152) (V 153) (V 154) (V 155)
    (V 156) (V 157) (V 158) (V 159) (V 160) (V 161) (V 162) (V 163) (V 164) (V 165) (V 166) (V 167)
    (V 168) (V 169) (V 170) (V 171) (V 172) (V 173) (V 174) (V 175) (V 176) (V 177) (V 178) (V 179)
    (V 180) (V 181) (V 182) (V 183) (V 184) (V 185) (V 186) (V 187) (V 188) (V 189) (V 190) (V 191)
    (V 192) (V 193) (V 194) (V 195) (V 196) (V 197) (V 198) (V 199) (V 200) (V 201) (V 202) (V 203)
    (V 204) (V 205) (V 206) (V 207) (V 208) (V 209) (V 210) (V 211) (V 212) (V 213) (V 214) (V 215)
    (V 216) (V 217) (V 218) (V 219) (V 220) (V 221) (V 222) (V 223) (V 224) (V 225) (V 226) (V 227)
  casesm* _ ∨ _
  all_goals first
    | exact hVert _ ‹_› (by norm_num) (by norm_num)
    | exact hEdge _ _ ‹_› (by decide +kernel)
    | exact ‹¬V 0› hVu
    | exact ‹¬V 4› hVv

/-- **Lower bound.** The unit-distance graph of `ℚ(√11)²` is not 3-colourable. -/
theorem not_colorable_three : ¬ (unitDistGraph L).Colorable 3 := fun h =>
  graph_not_colorable (h.of_hom (edgeGraph.hom (fun i => pt (P i)) pt_adj))

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

/-! ### The residue field of a place of `ℚ(√11)` over 7 is `𝔽₇` -/

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

/-- Every element of `L` is `(n₀ + n₁t)/D`, for either sign `t = ±√11`. -/
lemma exists_int_combo_11 {t : L} (ht : t = r11 ∨ t = -r11) (x : L) :
    ∃ D : ℕ, 0 < D ∧ ∃ n₀ n₁ : ℤ, (D : L) * x = n₀ + n₁ * t := by
  have hx : (x : ℝ) ∈ IntermediateField.adjoin ℚ ({√11, √11} : Set ℝ) := by
    rw [Set.pair_eq_singleton]
    exact x.2
  obtain ⟨d, hd, n₀, n₁, n₂, n₃, h⟩ := exists_int_combo (m := 11) (n := 11)
    (by push_cast; exact sqrt11_sq) (by push_cast; exact sqrt11_sq) hx
  have h' : (d : L) * x = (n₀ : L) + n₃ * (r11 * r11) + ((n₁ : L) + n₂) * r11 :=
    Subtype.ext (by push_cast [coe_r11]; linear_combination h)
  rcases ht with rfl | rfl
  · exact ⟨d, hd, n₀ + 11 * n₃, n₁ + n₂, by
      rw [h']; push_cast; linear_combination (n₃ : L) * r11_sq⟩
  · exact ⟨d, hd, n₀ + 11 * n₃, -(n₁ + n₂), by
      rw [h']; push_cast; linear_combination (n₃ : L) * r11_sq⟩

section Residue

variable (O : ValuationSubring L)

lemma r11_mem : r11 ∈ O :=
  mem_of_sq_eq O O.zero_mem (ofNat_mem O 11) (by rw [zero_mul, zero_add]; exact r11_sq)

variable (h7 : (7 : L) ∈ O.nonunits)
include h7

/-- An integer in `𝔪` is divisible by 7. -/
lemma seven_dvd_of_intCast_mem_nonunits {n : ℤ} (hn : (n : L) ∈ O.nonunits) : (7 : ℤ) ∣ n := by
  by_contra hnd
  obtain ⟨s, k, hs⟩ := exists_inv_mod_seven hnd
  apply one_notMem_nonunits O
  have hs' : (n : L) * s = 1 + 7 * k := by exact_mod_cast hs
  rw [show (1 : L) = n * s - 7 * k by linear_combination -hs']
  exact sub_mem (mul_mem_nonunits_right O hn (intCast_mem O s))
    (mul_mem_nonunits_right O h7 (intCast_mem O k))

/-- If an integer is `7ᵏ` times an element of `O`, then `7ᵏ` divides it. -/
lemma seven_pow_dvd_of_eq {k : ℕ} {n : ℤ} {γ : L} (hγ : γ ∈ O) (h : (n : L) = 7 ^ k * γ) :
    (7 : ℤ) ^ k ∣ n := by
  induction k generalizing n γ with
  | zero => simp
  | succ k ih =>
    have hn : (7 : ℤ) ∣ n := seven_dvd_of_intCast_mem_nonunits O h7 (by
      rw [h, pow_succ, mul_comm _ (7 : L), mul_assoc]
      exact mul_mem_nonunits_right O h7 (O.mul_mem _ _ (pow_mem (natCast_mem O 7) k) hγ))
    obtain ⟨m, rfl⟩ := hn
    have hm : (m : L) = 7 ^ k * γ := by
      apply mul_left_cancel₀ (show (7 : L) ≠ 0 by norm_num)
      push_cast at h
      linear_combination h
    obtain ⟨c, hc⟩ := ih hγ hm
    exact ⟨c, by rw [pow_succ]; linear_combination 7 * hc⟩

/-- Choosing the sign of `√11`, we may assume `t - 2 ∈ 𝔪`. -/
lemma exists_t : ∃ t : L, (t = r11 ∨ t = -r11) ∧ t - 2 ∈ O.nonunits := by
  have hO := r11_mem O
  have h : (r11 - 2) * (r11 + 2) ∈ O.nonunits := by
    rw [show (r11 - 2) * (r11 + 2) = 7 by linear_combination r11_sq]
    exact h7
  rcases mem_nonunits_or_of_mul_mem O (sub_mem hO (ofNat_mem O 2)) (add_mem hO (ofNat_mem O 2)) h
    with h | h
  · exact ⟨r11, Or.inl rfl, h⟩
  · refine ⟨-r11, Or.inr rfl, ?_⟩
    rw [show -r11 - 2 = -(r11 + 2) by ring]
    exact neg_mem h

/-- **Hensel's lemma at 7.** If `t² = 11` and `t - 2 ∈ 𝔪`, then for every `k` there is an integer
`b = 2 + 7m` with `t - b ∈ 7ᵏ⁺¹O`. -/
lemma hensel_seven {t : L} (ht : t ^ 2 = 11) (h2 : t - 2 ∈ O.nonunits) (k : ℕ) :
    ∃ m : ℤ, ∃ γ ∈ O, t - (2 + 7 * m) = 7 ^ (k + 1) * γ := by
  induction k with
  | zero =>
    -- `(t - 2)(t + 2) = 7` with `t + 2` a unit
    have hu : t + 2 ∉ O.nonunits := fun h => one_notMem_nonunits O (by
      rw [show (1 : L) = 2 * (t + 2) - 2 * (t - 2) - 7 by ring]
      exact sub_mem (sub_mem (mul_mem_nonunits_left O (ofNat_mem O 2) h)
        (mul_mem_nonunits_left O (ofNat_mem O 2) h2)) h7)
    have hne : t + 2 ≠ 0 := fun h => hu (h ▸ zero_mem _)
    refine ⟨0, (t + 2)⁻¹, inv_mem_of_notMem_nonunits O hu, ?_⟩
    simp only [Int.cast_zero, mul_zero, add_zero, zero_add, pow_one]
    field_simp
    linear_combination ht
  | succ k ih =>
    obtain ⟨m, γ, hγ, hb⟩ := ih
    -- `N = (2 + 7m)² - 11 = 7ᵏ⁺¹ · (-γ (2(2 + 7m) + 7ᵏ⁺¹γ))`
    have hN : (((2 + 7 * m) ^ 2 - 11 : ℤ) : L)
        = 7 ^ (k + 1) * (-γ * (2 * (2 + 7 * m) + 7 ^ (k + 1) * γ)) := by
      push_cast
      linear_combination ht - (2 + 7 * (m : L) + 7 ^ (k + 1) * γ + t) * hb
    have hmem : -γ * (2 * (2 + 7 * m) + 7 ^ (k + 1) * γ) ∈ O := by
      have := (-(⟨γ, hγ⟩ : O) * (2 * (2 + 7 * (m : O)) + 7 ^ (k + 1) * (⟨γ, hγ⟩ : O))).2
      push_cast at this
      exact this
    obtain ⟨c, hc⟩ := seven_pow_dvd_of_eq O h7 hmem hN
    have hc' : (2 + 7 * (m : L)) ^ 2 - 11 = 7 ^ (k + 1) * c := by exact_mod_cast hc
    refine ⟨m - 7 ^ k * 2 * c, -(1 + 4 * m) * γ - 2 * 7 ^ k * γ ^ 2, ?_, ?_⟩
    · have := (-(1 + 4 * (m : O)) * (⟨γ, hγ⟩ : O) - 2 * 7 ^ k * (⟨γ, hγ⟩ : O) ^ 2).2
      push_cast at this
      exact this
    · push_cast
      linear_combination (1 - 2 * (2 + 7 * (m : L) + 7 ^ (k + 1) * γ + t)) * hb - 2 * hc' + 2 * ht

/-- **The residue field is `𝔽₇`:** every `x ∈ O` is `≡` an integer modulo `𝔪`. -/
theorem exists_int_sub_mem {x : L} (hx : x ∈ O) : ∃ N : ℤ, x - N ∈ O.nonunits := by
  obtain ⟨t, ht, ht2⟩ := exists_t O h7
  have htsq : t ^ 2 = 11 := by
    rcases ht with rfl | rfl
    · exact r11_sq
    · rw [neg_sq]; exact r11_sq
  obtain ⟨D, hD, n₀, n₁, h⟩ := exists_int_combo_11 ht x
  obtain ⟨j, D', hD', rfl⟩ := Nat.exists_eq_pow_mul_and_not_dvd hD.ne' 7 (by norm_num)
  obtain ⟨m, γ, hγ, hb⟩ := hensel_seven O h7 htsq ht2 j
  -- `7ʲ (D'x - 7n₁γ) = n₀ + n₁(2 + 7m)`, an integer
  have hmem : (D' : L) * x - 7 * n₁ * γ ∈ O :=
    sub_mem (O.mul_mem _ _ (natCast_mem O D') hx)
      (O.mul_mem _ _ (O.mul_mem _ _ (ofNat_mem O 7) (intCast_mem O n₁)) hγ)
  have hN : ((n₀ + n₁ * (2 + 7 * m) : ℤ) : L) = 7 ^ j * ((D' : L) * x - 7 * n₁ * γ) := by
    push_cast at h ⊢
    linear_combination -h - (n₁ : L) * hb
  obtain ⟨N', hN'⟩ := seven_pow_dvd_of_eq O h7 hmem hN
  have hx' : (D' : L) * x - 7 * n₁ * γ = N' := by
    apply mul_left_cancel₀ (pow_ne_zero j (show (7 : L) ≠ 0 by norm_num))
    rw [← hN, hN']
    push_cast
    ring
  -- `D'` is prime to 7, hence a unit of `O`, with an inverse `s` modulo 7
  have hD'7 : ¬ (7 : ℤ) ∣ (D' : ℤ) := by exact_mod_cast hD'
  have hD'0 : (D' : L) ≠ 0 := by
    intro h0
    apply hD'7
    have : (D' : ℕ) = 0 := by exact_mod_cast h0
    rw [this]
    simp
  have hinv : (D' : L)⁻¹ ∈ O := inv_mem_of_notMem_nonunits O fun hm =>
    hD'7 (seven_dvd_of_intCast_mem_nonunits O h7 (by exact_mod_cast hm))
  obtain ⟨s, k, hs⟩ := exists_inv_mod_seven hD'7
  have hs' : (D' : L) * s = 1 + 7 * k := by exact_mod_cast hs
  -- `x - N's = 7 (n₁γ - N'k) / D'`
  refine ⟨N' * s, ?_⟩
  have e : x - ((N' * s : ℤ) : L) = 7 * ((n₁ * γ - N' * k) * (D' : L)⁻¹) := by
    have hDD : (D' : L) * (D' : L)⁻¹ = 1 := mul_inv_cancel₀ hD'0
    push_cast
    linear_combination (D' : L)⁻¹ * hx' - (N' : L) * (D' : L)⁻¹ * hs' - (x - N' * s) * hDD
  rw [e]
  exact mul_mem_nonunits_right O h7 (O.mul_mem _ _ (sub_mem (O.mul_mem _ _ (intCast_mem O n₁) hγ)
    (O.mul_mem _ _ (intCast_mem O N') (intCast_mem O k))) hinv)

end Residue

/-- `unitDistGraph L` maps to the graph of `x² + y²` over `L`. -/
def toSumSq : unitDistGraph L →g sumSqGraph L where
  toFun := id
  map_rel' := fun {p q} h => Subtype.ext (by push_cast; exact h)

/-- **Upper bound** (Moorhouse, Lemma 8.2 at 7). The unit-distance graph of `ℚ(√11)²` is
4-colourable. -/
theorem colorable_four : (unitDistGraph L).Colorable 4 := by
  obtain ⟨O, h7⟩ := exists_valuationSubring_seven_mem_nonunits L
  obtain ⟨red, hker⟩ := exists_red_of_residue_seven (R := O)
    (fun x => by
      obtain ⟨N, hN⟩ := exists_int_sub_mem O h7 x.2
      exact ⟨N, O.coe_mem_nonunits_iff.mp (by push_cast; exact hN)⟩)
    (fun n => ⟨fun h => seven_dvd_of_intCast_mem_nonunits O h7
        (by simpa using O.coe_mem_nonunits_iff.mpr h),
      fun ⟨c, hc⟩ => O.coe_mem_nonunits_iff.mp (by
        rw [hc]
        push_cast
        exact mul_mem_nonunits_right O h7 (intCast_mem O c))⟩)
  exact (sumSqGraph_colorable O red hker).of_hom toSumSq

/-- **Theorem.** The unit-distance graph of `ℚ(√11)²` has chromatic number 4. -/
theorem chromaticNumber_eq_four : (unitDistGraph L).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of colorable_four not_colorable_three

end Q11
