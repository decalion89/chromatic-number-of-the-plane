import QuadraticPlanes

/-!
# The colouring formula, and lower bounds from its unsatisfiability

For a graph `edgeGraph E` on `Fin n` with a listed edge `u₀v₀`, `colourCNF E u₀ v₀` is the formula of
`data/quadratic_planes/q{d}.cnf` as a Lean object: clauses are lists of nonzero integers, as in the DIMACS
format, and variable `3v + c + 1` says "vertex `v` has colour `c`". It has, in this order, the vertex clauses
`3v + 1 ∨ 3v + 2 ∨ 3v + 3`, for each listed edge `ab` and colour `c` the clause `¬(3a + c + 1) ∨ ¬(3b + c + 1)`,
and the unit clauses `3u₀ + 1` and `3v₀ + 2`.

`not_colorable_of_unsatisfiable` proves that the graph is not 3-colourable if no assignment satisfies the
formula. The field files `Sqrt{d}.lean` whose LRAT proofs are too large for `lrat_proof` state their theorems
with this hypothesis, `Unsatisfiable (colourCNF E u₀ v₀)`; the files check, when they are built, that
`parseDimacs` reads exactly this formula from `q{d}.cnf`, and cake_lpr, a verified checker, checked an LRAT
proof that `q{d}.cnf` is unsatisfiable (`data/quadratic_planes/cake_lpr_checks.txt`).
-/

namespace QuadraticPlanes

open LocalColouring

/-- The literal `m + 1`, which says that variable `m + 1` is true. -/
def posLit (m : ℕ) : ℤ := ((m + 1 : ℕ) : ℤ)

/-- The literal `-(m + 1)`, which says that variable `m + 1` is false. -/
def negLit (m : ℕ) : ℤ := -((m + 1 : ℕ) : ℤ)

/-- The colouring formula of `edgeGraph E`, with the edge `u₀v₀` coloured `0, 1`. -/
def colourCNF {n : ℕ} (E : List (Fin n × Fin n)) (u₀ v₀ : Fin n) : List (List ℤ) :=
  (List.range n).map (fun k => [posLit (3 * k), posLit (3 * k + 1), posLit (3 * k + 2)]) ++
  E.flatMap (fun e => (List.range 3).map fun c => [negLit (3 * e.1.val + c), negLit (3 * e.2.val + c)]) ++
  [[posLit (3 * u₀.val)], [posLit (3 * v₀.val + 1)]]

/-- A literal under an assignment `σ` of the variables `1, 2, …`: `k > 0` holds if `σ k`, and `-k` if `¬σ k`. -/
def LitHolds (σ : ℕ → Prop) (l : ℤ) : Prop := if 0 < l then σ l.toNat else ¬σ l.natAbs

/-- No assignment satisfies every clause of `F`: for every assignment, some clause has no true literal. -/
def Unsatisfiable (F : List (List ℤ)) : Prop := ∀ σ : ℕ → Prop, ∃ c ∈ F, ∀ l ∈ c, ¬LitHolds σ l

lemma litHolds_posLit (σ : ℕ → Prop) (m : ℕ) : LitHolds σ (posLit m) ↔ σ (m + 1) := by
  unfold LitHolds posLit
  have h : (0 : ℤ) < ((m + 1 : ℕ) : ℤ) := by omega
  simp only [h, ↓reduceIte, Int.toNat_natCast]

lemma litHolds_negLit (σ : ℕ → Prop) (m : ℕ) : LitHolds σ (negLit m) ↔ ¬σ (m + 1) := by
  unfold LitHolds negLit
  have h : ¬(0 : ℤ) < -((m + 1 : ℕ) : ℤ) := by omega
  simp only [h, ↓reduceIte, Int.natAbs_neg, Int.natAbs_natCast]

lemma memE_of_mem {n : ℕ} {a b : Fin n} : ∀ {l : List (Fin n × Fin n)}, (a, b) ∈ l → memE a.val b.val l = true
  | [], h => by simp at h
  | _ :: _, h => by
    simp only [memE, Bool.or_eq_true, Bool.and_eq_true]
    rcases List.mem_cons.mp h with h | h
    · subst h
      exact Or.inl ⟨Nat.beq_refl _, Nat.beq_refl _⟩
    · exact Or.inr (memE_of_mem h)

lemma goodEdge_of_mem {n : ℕ} {E : List (Fin n × Fin n)} {e : Fin n × Fin n} (he : e ∈ E) {c : ℕ} (hc : c < 3) :
    goodEdge E (3 * e.1.val + c) (3 * e.2.val + c) = true := by
  have h1 : (3 * e.1.val + c) % 3 = c := by omega
  have h2 : (3 * e.2.val + c) % 3 = c := by omega
  have h3 : (3 * e.1.val + c) / 3 = e.1.val := by omega
  have h4 : (3 * e.2.val + c) / 3 = e.2.val := by omega
  have hm : memE e.1.val e.2.val E = true := memE_of_mem (by simpa using he)
  simp only [goodEdge, h1, h2, h3, h4, Nat.beq_refl, hm, Bool.true_and, Bool.true_or]

/-- **The lower bound from the formula.** If no assignment satisfies the colouring formula of `edgeGraph E`,
the graph is not 3-colourable. -/
theorem not_colorable_of_unsatisfiable {n : ℕ} {E : List (Fin n × Fin n)} (hE : noLoops E = true)
    {u₀ v₀ : Fin n} (h₀ : (u₀, v₀) ∈ E) (h : Unsatisfiable (colourCNF E u₀ v₀)) :
    ¬(edgeGraph E).Colorable 3 := by
  rintro ⟨C⟩
  obtain ⟨V, hVert, hEdge, hVu, hVv⟩ := clause_facts hE h₀ C
  -- variable `m + 1` is `V m`
  obtain ⟨cl, hcl, hF⟩ := h (fun k => V (k - 1))
  have hp : ∀ m, LitHolds (fun k => V (k - 1)) (posLit m) ↔ V m := fun m => by
    rw [litHolds_posLit]
    simp
  have hn : ∀ m, LitHolds (fun k => V (k - 1)) (negLit m) ↔ ¬V m := fun m => by
    rw [litHolds_negLit]
    simp
  simp only [colourCNF, List.mem_append, List.mem_map, List.mem_range, List.mem_flatMap, List.mem_cons,
    List.not_mem_nil, or_false] at hcl
  rcases hcl with (⟨k, hk, rfl⟩ | ⟨e, he, c, hc, rfl⟩) | rfl | rfl
  · -- a vertex clause
    have h0 := hF (posLit (3 * k)) (by simp)
    have h1 := hF (posLit (3 * k + 1)) (by simp)
    have h2 := hF (posLit (3 * k + 2)) (by simp)
    rw [hp] at h0 h1 h2
    exact hVert (3 * k) ⟨h0, h1, h2⟩ (by omega) (by omega)
  · -- an edge clause
    have h1 := hF (negLit (3 * e.1.val + c)) (by simp)
    have h2 := hF (negLit (3 * e.2.val + c)) (by simp)
    rw [hn, not_not] at h1 h2
    exact hEdge _ _ ⟨h1, h2⟩ (goodEdge_of_mem he hc)
  · -- the unit clause of `u₀`
    exact hF (posLit (3 * u₀.val)) (by simp) ((hp _).2 hVu)
  · -- the unit clause of `v₀`
    exact hF (posLit (3 * v₀.val + 1)) (by simp) ((hp _).2 hVv)

/-- The points of a graph as a balanced binary tree: `node k l r` holds the points of `l` (the first `k`) and then
those of `r`. The kernel finds a point in a number of steps that grows like the logarithm of their number, where a
list needs as many steps as the point's index. -/
inductive PtTree where
  | leaf (p : ℤ × ℤ × ℤ × ℤ)
  | node (k : ℕ) (l r : PtTree)

/-- The point with index `i` (counting from 0). -/
def PtTree.get : PtTree → ℕ → ℤ × ℤ × ℤ × ℤ
  | .leaf p, _ => p
  | .node k l r, i => if i < k then l.get i else r.get (i - k)

/-- The balanced tree of a list of points, the first half on the left; `f` bounds the depth (with `2 ^ f` at least
the length, every point is a leaf). -/
def PtTree.ofList : ℕ → List (ℤ × ℤ × ℤ × ℤ) → PtTree
  | 0, l => .leaf (l.headD (0, 0, 0, 0))
  | f + 1, l =>
    if l.length ≤ 1 then .leaf (l.headD (0, 0, 0, 0))
    else .node (l.length / 2) (ofList f (l.take (l.length / 2))) (ofList f (l.drop (l.length / 2)))

/-- The clauses of a DIMACS file with one clause per line (the header and comment lines are skipped). The
field files use it, when they are built, to check that `q{d}.cnf` is `colourCNF E u₀ v₀`. -/
def parseDimacs (s : String) : List (List ℤ) :=
  ((s.splitOn "\n").filter fun l => !(l.startsWith "c" || l.startsWith "p" || l.all Char.isWhitespace)).map fun l =>
    (((l.splitOn " ").filter (· ≠ "")).filterMap String.toInt?).filter (· ≠ 0)

end QuadraticPlanes
