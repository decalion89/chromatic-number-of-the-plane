import Mathlib

/-!
# The plane is 7-colourable: `χ(ℝ²) ≤ 7`

`SimpleGraph.UnitDistancePlaneGraph` is copied from google-deepmind/formal-conjectures
(`FormalConjecturesForMathlib/Combinatorics/SimpleGraph/UnitDistancePlaneGraph.lean`,
Apache License 2.0), and `PlaneSeven.hadwigerNelson_atMostSeven` is the statement of their
`Erdos508.HadwigerNelsonAtMostSeven`.

The colouring is a brick version of Isbell's hexagonal colouring (the square-tiling approach the
formal-conjectures docstring attributes to László Székely). Measured in units of `7/10`, row `j`
is the strip `j ≤ y < j + 1`, and brick `i` of row `j` is `i ≤ x - (11/20) j < i + 1`. Brick
`(i, j)` gets colour `i + 3 j` mod 7. Two points of one brick are less than `1` apart, since
`0.7² + 0.7² < 1`, and two different bricks of one colour are more than `1` apart.
-/

namespace SimpleGraph

/-- A unit distance graph in ℝ² (copied from formal-conjectures, Apache License 2.0). The only
changes are `Dist.dist` for `dist` and `_root_.dist_comm` for `dist_comm`: formal-conjectures imports only `Mathlib.Combinatorics.SimpleGraph.Basic`
and `Mathlib.Analysis.InnerProductSpace.PiL2`, while under `import Mathlib` the name `dist` inside
`namespace SimpleGraph` means the graph distance `SimpleGraph.dist` (and
likewise `dist_comm`). -/
def UnitDistancePlaneGraph (V : Set (EuclideanSpace ℝ (Fin 2))) : SimpleGraph V where
  Adj x y := Dist.dist x y = 1
  symm.symm x y := by simp [_root_.dist_comm]

end SimpleGraph

namespace PlaneSeven

open SimpleGraph

/-- The arithmetic core, in units of `7/10` (so distance `1` becomes `10/7`): points of bricks
`(i, j)` and `(i', j')` with `i + 3 j ≡ i' + 3 j' (mod 7)` are never at distance `10/7`. -/
theorem brick_key (X Y X' Y' : ℝ) (i j i' j' : ℤ)
    (hj : (j : ℝ) ≤ Y) (hj1 : Y < j + 1) (hj' : (j' : ℝ) ≤ Y') (hj'1 : Y' < j' + 1)
    (hi : (i : ℝ) ≤ X - j * (11 / 20)) (hi1 : X - j * (11 / 20) < i + 1)
    (hi' : (i' : ℝ) ≤ X' - j' * (11 / 20)) (hi'1 : X' - j' * (11 / 20) < i' + 1)
    (h7 : (7 : ℤ) ∣ (i' + 3 * j') - (i + 3 * j))
    (hd : (X - X') ^ 2 + (Y - Y') ^ 2 = 100 / 49) : False := by
  -- a horizontal gap `a` and a vertical gap `b` with `a² + b² > 100/49` contradict `hd`
  have gapX : ∀ a : ℝ, 0 ≤ a → a ^ 2 > 100 / 49 → (a < X - X' ∨ a < X' - X) → False := by
    intro a ha ha2 h
    rcases h with h | h <;> nlinarith [mul_self_lt_mul_self ha h, sq_nonneg (Y - Y')]
  have gapXY : ∀ a b : ℝ, 0 ≤ a → 0 ≤ b → a ^ 2 + b ^ 2 > 100 / 49 →
      ((a < X - X' ∧ b < Y - Y') ∨ (a < X' - X ∧ b < Y' - Y)) → False := by
    intro a b ha hb hab h
    rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · nlinarith [mul_self_lt_mul_self ha h1, mul_self_lt_mul_self hb h2]
    · nlinarith [mul_self_lt_mul_self ha h1, mul_self_lt_mul_self hb h2]
  have gapY : ∀ b : ℝ, 0 ≤ b → b ^ 2 > 100 / 49 → (b < Y - Y' ∨ b < Y' - Y) → False := by
    intro b hb hb2 h
    rcases h with h | h <;> nlinarith [mul_self_lt_mul_self hb h, sq_nonneg (X - X')]
  rcases (by omega : j' ≤ j - 3 ∨ j' = j - 2 ∨ j' = j - 1 ∨ j' = j ∨ j' = j + 1 ∨ j' = j + 2 ∨
      j + 3 ≤ j') with h | h | h | h | h | h | h
  · -- three or more rows down
    have hr : (j' : ℝ) ≤ j - 3 := by exact_mod_cast h
    exact gapY 2 (by norm_num) (by norm_num) (Or.inl (by linarith))
  · -- two rows down: `i' - i ≡ 6 (mod 7)`
    have hr : (j' : ℝ) = j - 2 := by exact_mod_cast h
    rcases (by omega : i' ≤ i - 1 ∨ i + 6 ≤ i') with hc | hc
    · have hcr : (i' : ℝ) ≤ i - 1 := by exact_mod_cast hc
      exact gapXY (11 / 10) 1 (by norm_num) (by norm_num) (by norm_num)
        (Or.inl ⟨by linarith, by linarith⟩)
    · have hcr : (i : ℝ) + 6 ≤ i' := by exact_mod_cast hc
      exact gapX (39 / 10) (by norm_num) (by norm_num) (Or.inr (by linarith))
  · -- one row down: `i' - i ≡ 3 (mod 7)`
    have hr : (j' : ℝ) = j - 1 := by exact_mod_cast h
    rcases (by omega : i' ≤ i - 4 ∨ i + 3 ≤ i') with hc | hc
    · have hcr : (i' : ℝ) ≤ i - 4 := by exact_mod_cast hc
      exact gapX (71 / 20) (by norm_num) (by norm_num) (Or.inl (by linarith))
    · have hcr : (i : ℝ) + 3 ≤ i' := by exact_mod_cast hc
      exact gapX (29 / 20) (by norm_num) (by norm_num) (Or.inr (by linarith))
  · -- the same row: `i' - i ≡ 0 (mod 7)`
    have hr : (j' : ℝ) = j := by exact_mod_cast h
    rcases (by omega : i' ≤ i - 7 ∨ i' = i ∨ i + 7 ≤ i') with hc | hc | hc
    · have hcr : (i' : ℝ) ≤ i - 7 := by exact_mod_cast hc
      exact gapX 6 (by norm_num) (by norm_num) (Or.inl (by linarith))
    · -- the same brick
      have hcr : (i' : ℝ) = i := by exact_mod_cast hc
      have hx1 : X - X' < 1 := by linarith
      have hx2 : X' - X < 1 := by linarith
      have hy1 : Y - Y' < 1 := by linarith
      have hy2 : Y' - Y < 1 := by linarith
      nlinarith [mul_pos (by linarith : (0 : ℝ) < 1 - (X - X')) (by linarith : (0 : ℝ) < 1 - (X' - X)),
        mul_pos (by linarith : (0 : ℝ) < 1 - (Y - Y')) (by linarith : (0 : ℝ) < 1 - (Y' - Y))]
    · have hcr : (i : ℝ) + 7 ≤ i' := by exact_mod_cast hc
      exact gapX 6 (by norm_num) (by norm_num) (Or.inr (by linarith))
  · -- one row up: `i' - i ≡ 4 (mod 7)`
    have hr : (j' : ℝ) = j + 1 := by exact_mod_cast h
    rcases (by omega : i' ≤ i - 3 ∨ i + 4 ≤ i') with hc | hc
    · have hcr : (i' : ℝ) ≤ i - 3 := by exact_mod_cast hc
      exact gapX (29 / 20) (by norm_num) (by norm_num) (Or.inl (by linarith))
    · have hcr : (i : ℝ) + 4 ≤ i' := by exact_mod_cast hc
      exact gapX (71 / 20) (by norm_num) (by norm_num) (Or.inr (by linarith))
  · -- two rows up: `i' - i ≡ 1 (mod 7)`
    have hr : (j' : ℝ) = j + 2 := by exact_mod_cast h
    rcases (by omega : i' ≤ i - 6 ∨ i + 1 ≤ i') with hc | hc
    · have hcr : (i' : ℝ) ≤ i - 6 := by exact_mod_cast hc
      exact gapX (39 / 10) (by norm_num) (by norm_num) (Or.inl (by linarith))
    · have hcr : (i : ℝ) + 1 ≤ i' := by exact_mod_cast hc
      exact gapXY (11 / 10) 1 (by norm_num) (by norm_num) (by norm_num)
        (Or.inr ⟨by linarith, by linarith⟩)
  · -- three or more rows up
    have hr : (j : ℝ) + 3 ≤ j' := by exact_mod_cast h
    exact gapY 2 (by norm_num) (by norm_num) (Or.inr (by linarith))

/-- the row of a point -/
noncomputable def row (p : EuclideanSpace ℝ (Fin 2)) : ℤ := ⌊p 1 * (10 / 7)⌋

/-- the brick of a point within its row -/
noncomputable def col (p : EuclideanSpace ℝ (Fin 2)) : ℤ :=
  ⌊p 0 * (10 / 7) - (row p : ℝ) * (11 / 20)⌋

/-- the colour of a point: `i + 3 j` mod 7 for brick `(i, j)` -/
noncomputable def colour (p : EuclideanSpace ℝ (Fin 2)) : ZMod 7 :=
  ((col p + 3 * row p : ℤ) : ZMod 7)

theorem colour_ne {p q : EuclideanSpace ℝ (Fin 2)} (h : dist p q = 1) : colour p ≠ colour q := by
  intro hc
  have h7 : (7 : ℤ) ∣ (col q + 3 * row q) - (col p + 3 * row p) := by
    have := (ZMod.intCast_eq_intCast_iff_dvd_sub _ _ 7).mp hc
    exact_mod_cast this
  have h2 : (p 0 - q 0) ^ 2 + (p 1 - q 1) ^ 2 = 1 := by
    have := EuclideanSpace.dist_sq_eq p q
    rw [h, Fin.sum_univ_two, Real.dist_eq, Real.dist_eq, sq_abs, sq_abs] at this
    linarith
  exact brick_key (p 0 * (10 / 7)) (p 1 * (10 / 7)) (q 0 * (10 / 7)) (q 1 * (10 / 7))
    (col p) (row p) (col q) (row q)
    (Int.floor_le _) (Int.lt_floor_add_one _) (Int.floor_le _) (Int.lt_floor_add_one _)
    (Int.floor_le _) (Int.lt_floor_add_one _) (Int.floor_le _) (Int.lt_floor_add_one _) h7
    (by linear_combination (100 / 49) * h2)

/-- the brick colouring of the unit-distance graph of the plane -/
noncomputable def brickColoring : (UnitDistancePlaneGraph Set.univ).Coloring (ZMod 7) :=
  Coloring.mk (fun p => colour p.1) fun {_ _} h => colour_ne h

/-- `Erdos508.HadwigerNelsonAtMostSeven` of formal-conjectures. -/
theorem hadwigerNelson_atMostSeven :
    SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) ≤ 7 := by
  have := brickColoring.colorable.chromaticNumber_le
  simpa [ZMod.card] using this

end PlaneSeven

