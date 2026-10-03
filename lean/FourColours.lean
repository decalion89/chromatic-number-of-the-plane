import Mathlib

namespace FourColours

/-! ## Gaussian-integer facts -/

/-- `(3 + 4i)^j = A j + B j i`. -/
def AB : ℕ → ℤ × ℤ
  | 0 => (1, 0)
  | j + 1 => (3 * (AB j).1 - 4 * (AB j).2, 4 * (AB j).1 + 3 * (AB j).2)

/-- The real part of `(3 + 4i)^j`. -/
def A (j : ℕ) : ℤ := (AB j).1
/-- The imaginary part of `(3 + 4i)^j`. -/
def B (j : ℕ) : ℤ := (AB j).2

@[simp] lemma A_zero : A 0 = 1 := rfl
@[simp] lemma B_zero : B 0 = 0 := rfl
lemma A_succ (j : ℕ) : A (j + 1) = 3 * A j - 4 * B j := rfl
lemma B_succ (j : ℕ) : B (j + 1) = 4 * A j + 3 * B j := rfl

lemma AB_sq (j : ℕ) : A j ^ 2 + B j ^ 2 = 25 ^ j := by
  induction j with
  | zero => simp
  | succ j ih => rw [A_succ, B_succ, pow_succ (25 : ℤ), ← ih]; ring

lemma AB_mod2 (j : ℕ) : A j % 2 = 1 ∧ B j % 2 = 0 := by
  induction j with
  | zero => simp
  | succ j ih => rw [A_succ, B_succ]; omega

lemma AB_mod5 (j : ℕ) (hj : 1 ≤ j) : A j % 5 = 3 ∧ B j % 5 = 4 := by
  induction j with
  | zero => omega
  | succ j ih =>
    rcases Nat.eq_zero_or_pos j with rfl | hj'
    · decide
    · have := ih hj'; rw [A_succ, B_succ]; omega

lemma AB_mod3 (j : ℕ) : (A j % 3 = 0 ∧ B j % 3 ≠ 0) ∨ (A j % 3 ≠ 0 ∧ B j % 3 = 0) := by
  induction j with
  | zero => simp
  | succ j ih => rw [A_succ, B_succ]; omega

/-! ## The sets `S_N` and the structure lemma -/

/-- `x ∈ [1/3, 2/3] + ℤ`. -/
def good (x : ℝ) : Prop := ∃ m : ℤ, 1 / 3 ≤ x - m ∧ x - m ≤ 2 / 3

lemma good_congr {x y : ℝ} (h : good x) (hxy : x = y) : good y := hxy ▸ h

/-- `w = p + qi` lies in `S_{5^k}`: the conditions for `γ = ρ^j, iρ^j, ρ̄^j, iρ̄^j`, `j ≤ k`. -/
def inS (k : ℕ) (p q : ℝ) : Prop :=
  ∀ j ≤ k, good ((A j * p + B j * q) / 5 ^ j) ∧ good ((A j * q - B j * p) / 5 ^ j) ∧
    good ((A j * p - B j * q) / 5 ^ j) ∧ good ((A j * q + B j * p) / 5 ^ j)

/-- The two types: `w = (5^k/6)(E₁ + E₂ i) + x` with either `E ≡ (3, 3) (mod 6)` and `|x|² < 1/18` (type c), or
`E₁, E₂` even and prime to 3 and `x = 0` (type q). -/
def Rep (k : ℕ) (p q : ℝ) : Prop :=
  ∃ E₁ E₂ : ℤ, ∃ x₁ x₂ : ℝ, p = 5 ^ k / 6 * E₁ + x₁ ∧ q = 5 ^ k / 6 * E₂ + x₂ ∧
    ((E₁ % 6 = 3 ∧ E₂ % 6 = 3 ∧ x₁ ^ 2 + x₂ ^ 2 < 1 / 18) ∨
      (E₁ % 2 = 0 ∧ E₁ % 3 ≠ 0 ∧ E₂ % 2 = 0 ∧ E₂ % 3 ≠ 0 ∧ x₁ = 0 ∧ x₂ = 0))

lemma good_half {x : ℝ} (h : good x) : ∃ m : ℤ, ∃ y : ℝ, x = m + 1 / 2 + y ∧ -1 / 6 ≤ y ∧ y ≤ 1 / 6 :=
  let ⟨m, h1, h2⟩ := h
  ⟨m, x - m - 1 / 2, by ring, by linarith, by linarith⟩

lemma rep_zero {p q : ℝ} (h : inS 0 p q) : Rep 0 p q := by
  obtain ⟨h1, h2, -, -⟩ := h 0 le_rfl
  simp only [A_zero, B_zero, Int.cast_one, Int.cast_zero, one_mul, zero_mul, add_zero, sub_zero,
    pow_zero, div_one] at h1 h2
  obtain ⟨m₁, y₁, rfl, hy₁, hy₁'⟩ := good_half h1
  obtain ⟨m₂, y₂, rfl, hy₂, hy₂'⟩ := good_half h2
  unfold Rep
  simp only [pow_zero]
  by_cases hc : y₁ ^ 2 = 1 / 36 ∧ y₂ ^ 2 = 1 / 36
  · obtain ⟨hc1, hc2⟩ := hc
    have e1 : y₁ = 1 / 6 ∨ y₁ = -1 / 6 := by
      have : (y₁ - 1 / 6) * (y₁ + 1 / 6) = 0 := by linear_combination hc1
      rcases mul_eq_zero.1 this with h | h
      · left; linarith
      · right; linarith
    have e2 : y₂ = 1 / 6 ∨ y₂ = -1 / 6 := by
      have : (y₂ - 1 / 6) * (y₂ + 1 / 6) = 0 := by linear_combination hc2
      rcases mul_eq_zero.1 this with h | h
      · left; linarith
      · right; linarith
    have key : ∀ (m : ℤ) (y : ℝ), (y = 1 / 6 ∨ y = -1 / 6) →
        ∃ E : ℤ, (m : ℝ) + 1 / 2 + y = 1 / 6 * E + 0 ∧ E % 2 = 0 ∧ E % 3 ≠ 0 := by
      intro m y hy
      rcases hy with rfl | rfl
      · exact ⟨6 * m + 4, by push_cast; ring, by omega, by omega⟩
      · exact ⟨6 * m + 2, by push_cast; ring, by omega, by omega⟩
    obtain ⟨E₁, h1, h1', h1''⟩ := key m₁ y₁ e1
    obtain ⟨E₂, h2, h2', h2''⟩ := key m₂ y₂ e2
    exact ⟨E₁, E₂, 0, 0, by rw [h1], by rw [h2],
      Or.inr ⟨h1', h1'', h2', h2'', rfl, rfl⟩⟩
  · refine ⟨6 * m₁ + 3, 6 * m₂ + 3, y₁, y₂, by push_cast; ring, by push_cast; ring,
      Or.inl ⟨by omega, by omega, ?_⟩⟩
    have b1 : y₁ ^ 2 ≤ 1 / 36 := by nlinarith
    have b2 : y₂ ^ 2 ≤ 1 / 36 := by nlinarith
    by_contra hlt
    exact hc ⟨by linarith, by linarith⟩

/-! ### The induction step -/

lemma sq_lb1 {a : ℤ} {s : ℝ} (h1 : -1 / 6 ≤ a / 5 + s) (h2 : a / 5 + s ≤ 1 / 6) (ha : a ≠ 0) :
    1 / 900 ≤ s ^ 2 := by
  rcases lt_or_gt_of_ne ha with ha | ha
  · have : (a : ℝ) ≤ -1 := by exact_mod_cast (show a ≤ -1 by omega)
    nlinarith
  · have : (1 : ℝ) ≤ a := by exact_mod_cast (show 1 ≤ a by omega)
    nlinarith

lemma sq_lb2 {a : ℤ} {s : ℝ} (h1 : -1 / 6 ≤ a / 5 + s) (h2 : a / 5 + s ≤ 1 / 6) (ha : 2 ≤ a ∨ a ≤ -2) :
    49 / 900 ≤ s ^ 2 := by
  rcases ha with ha | ha
  · have : (2 : ℝ) ≤ a := by exact_mod_cast ha
    nlinarith
  · have : (a : ℝ) ≤ -2 := by exact_mod_cast ha
    nlinarith

/-- Lemma 3(b), crude: a nonzero point of `Λ₊ ∪ Λ₋` modulo `ℤ[i]` is at squared distance `≥ 1/18` from `Sq₀`. -/
lemma five_dvd_of_small {a b : ℤ} {s t : ℝ} (hab : (a + 2 * b) % 5 = 0 ∨ (a - 2 * b) % 5 = 0)
    (h1 : -1 / 6 ≤ a / 5 + s) (h2 : a / 5 + s ≤ 1 / 6) (h3 : -1 / 6 ≤ b / 5 + t) (h4 : b / 5 + t ≤ 1 / 6)
    (hst : s ^ 2 + t ^ 2 < 1 / 18) : a % 5 = 0 ∧ b % 5 = 0 := by
  by_contra hne
  have hb0 : b ≠ 0 := by rintro rfl; apply hne; omega
  have ha0 : a ≠ 0 := by rintro rfl; apply hne; omega
  have hbig : 2 ≤ a ∨ a ≤ -2 ∨ 2 ≤ b ∨ b ≤ -2 := by
    by_contra h
    have ha : a = 1 ∨ a = -1 := by omega
    have hb : b = 1 ∨ b = -1 := by omega
    rcases ha with rfl | rfl <;> rcases hb with rfl | rfl <;> revert hab <;> decide
  have : ((2 ≤ a ∨ a ≤ -2) ∧ b ≠ 0) ∨ (a ≠ 0 ∧ (2 ≤ b ∨ b ≤ -2)) := by omega
  rcases this with ⟨ha, hb⟩ | ⟨ha, hb⟩
  · have := sq_lb2 h1 h2 ha
    have := sq_lb1 h3 h4 hb
    linarith
  · have := sq_lb1 h1 h2 ha
    have := sq_lb2 h3 h4 hb
    linarith

/-- A condition at level `k + 1` for type c: `g = 3h` with `h` odd. -/
lemma good_c {h X : ℤ} (hh : h % 2 = 1) {s y : ℝ} (hy : y = 3 * h / 6 + X / 5 + s) (hg : good y) :
    ∃ a : ℤ, (a - X) % 5 = 0 ∧ -1 / 6 ≤ a / 5 + s ∧ a / 5 + s ≤ 1 / 6 := by
  obtain ⟨m, h1, h2⟩ := hg
  obtain ⟨c, rfl⟩ : ∃ c, h = 2 * c + 1 := ⟨h / 2, by omega⟩
  refine ⟨X + 5 * (c - m), by omega, ?_, ?_⟩ <;> (subst hy; push_cast at h1 h2 ⊢; linarith)

/-- A condition at level `k + 1` for type q (no error term). -/
lemma good_q {g X : ℤ} {y : ℝ} (hy : y = g / 6 + X / 5) (hg : good y) :
    ∃ m : ℤ, 10 ≤ 5 * g + 6 * X - 30 * m ∧ 5 * g + 6 * X - 30 * m ≤ 20 := by
  obtain ⟨m, h1, h2⟩ := hg
  refine ⟨m, ?_, ?_⟩
  · have : (10 : ℝ) ≤ 5 * g + 6 * X - 30 * m := by subst hy; linarith
    exact_mod_cast this
  · have : 5 * g + 6 * X - 30 * m ≤ (20 : ℝ) := by subst hy; linarith
    exact_mod_cast this

lemma q_one {g X m : ℤ} (hg : g % 2 = 0 ∧ g % 3 ≠ 0) (h1 : 10 ≤ 5 * g + 6 * X - 30 * m)
    (h2 : 5 * g + 6 * X - 30 * m ≤ 20) : X % 5 = 0 ∨ X % 5 = 1 ∨ X % 5 = 4 := by
  obtain ⟨g₀, rfl⟩ : ∃ g₀, g = 2 * g₀ := ⟨g / 2, by omega⟩
  have h3 : g₀ % 3 = 1 ∨ g₀ % 3 = 2 := by omega
  rcases h3 with h3 | h3
  · obtain ⟨t, rfl⟩ : ∃ t, g₀ = 3 * t + 1 := ⟨g₀ / 3, by omega⟩
    omega
  · obtain ⟨t, rfl⟩ : ∃ t, g₀ = 3 * t + 2 := ⟨g₀ / 3, by omega⟩
    omega

lemma q_pair {g g' X Y m m' : ℤ} (hg : g % 2 = 0 ∧ g % 3 ≠ 0) (hg' : g' % 2 = 0 ∧ g' % 3 ≠ 0)
    (h1 : 10 ≤ 5 * g + 6 * X - 30 * m) (h2 : 5 * g + 6 * X - 30 * m ≤ 20)
    (h3 : 10 ≤ 5 * g' + 6 * Y - 30 * m') (h4 : 5 * g' + 6 * Y - 30 * m' ≤ 20)
    (hXY : (X + 2 * Y) % 5 = 0 ∨ (X - 2 * Y) % 5 = 0) : X % 5 = 0 ∧ Y % 5 = 0 := by
  have := q_one hg h1 h2
  have := q_one hg' h3 h4
  omega

lemma g_q {a b r t : ℤ} (hab3 : (a % 3 = 0 ∧ b % 3 ≠ 0) ∨ (a % 3 ≠ 0 ∧ b % 3 = 0))
    (hr : r = 2 ∨ r = 4) (ht : t = 2 ∨ t = 4) : (a * r + b * t) % 2 = 0 ∧ (a * r + b * t) % 3 ≠ 0 := by
  rcases hr with rfl | rfl <;> rcases ht with rfl | rfl <;> omega

lemma five_aux {X Y Z W : ℤ} (hX : X % 5 = 0) (hY : Y % 5 = 0) (h : X + Y = 5 * Z + W) : W % 5 = 0 := by
  omega

lemma five_aux' {X Y Z T : ℤ} (hX : X % 5 = 0) (hY : Y % 5 = 0) (h : X - Y = 5 * Z + 3 * T) :
    T % 5 = 0 := by
  omega

/-- The induction step, for `a + bi = (3 + 4i)^{k+1}` and `N = 5^k`. -/
lemma step_core {N : ℝ} (hN : 0 < N) {a b : ℤ} (hab2 : a % 2 = 1 ∧ b % 2 = 0) (hab5 : a % 5 = 3 ∧ b % 5 = 4)
    (hab3 : (a % 3 = 0 ∧ b % 3 ≠ 0) ∨ (a % 3 ≠ 0 ∧ b % 3 = 0)) (hsq : (a : ℝ) ^ 2 + (b : ℝ) ^ 2 = (5 * N) ^ 2)
    {E₁ E₂ : ℤ} {x₁ x₂ : ℝ}
    (htype : (E₁ % 6 = 3 ∧ E₂ % 6 = 3 ∧ x₁ ^ 2 + x₂ ^ 2 < 1 / 18) ∨
      (E₁ % 2 = 0 ∧ E₁ % 3 ≠ 0 ∧ E₂ % 2 = 0 ∧ E₂ % 3 ≠ 0 ∧ x₁ = 0 ∧ x₂ = 0))
    (c1 : good ((a * (N / 6 * E₁ + x₁) + b * (N / 6 * E₂ + x₂)) / (5 * N)))
    (c2 : good ((a * (N / 6 * E₂ + x₂) - b * (N / 6 * E₁ + x₁)) / (5 * N)))
    (c3 : good ((a * (N / 6 * E₁ + x₁) - b * (N / 6 * E₂ + x₂)) / (5 * N)))
    (c4 : good ((a * (N / 6 * E₂ + x₂) + b * (N / 6 * E₁ + x₁)) / (5 * N))) :
    ∃ E₁' E₂' : ℤ, N / 6 * E₁ = 5 * N / 6 * E₁' ∧ N / 6 * E₂ = 5 * N / 6 * E₂' ∧
      ((E₁' % 6 = 3 ∧ E₂' % 6 = 3 ∧ x₁ ^ 2 + x₂ ^ 2 < 1 / 18) ∨
        (E₁' % 2 = 0 ∧ E₁' % 3 ≠ 0 ∧ E₂' % 2 = 0 ∧ E₂' % 3 ≠ 0 ∧ x₁ = 0 ∧ x₂ = 0)) := by
  -- `E = 5r + 6T`
  obtain ⟨r₁, hr₁⟩ : ∃ r, r = (5 * E₁) % 6 := ⟨_, rfl⟩
  obtain ⟨r₂, hr₂⟩ : ∃ r, r = (5 * E₂) % 6 := ⟨_, rfl⟩
  obtain ⟨T₁, hT₁⟩ : ∃ T, E₁ = 5 * r₁ + 6 * T := ⟨(E₁ - 5 * r₁) / 6, by omega⟩
  obtain ⟨T₂, hT₂⟩ : ∃ T, E₂ = 5 * r₂ + 6 * T := ⟨(E₂ - 5 * r₂) / 6, by omega⟩
  -- the four integers `X₁, Y₁, X₂, Y₂` and their relations modulo 5
  obtain ⟨u, hu⟩ : ∃ u, a - 2 * b = 5 * u := ⟨(a - 2 * b) / 5, by omega⟩
  obtain ⟨v, hv⟩ : ∃ v, 2 * a + b = 5 * v := ⟨(2 * a + b) / 5, by omega⟩
  have m1 : (a * T₁ + b * T₂ + 2 * (a * T₂ - b * T₁)) % 5 = 0 := by
    rw [show a * T₁ + b * T₂ + 2 * (a * T₂ - b * T₁) = 5 * (u * T₁ + v * T₂) by
      linear_combination T₁ * hu + T₂ * hv]
    omega
  have m2 : (a * T₁ - b * T₂ - 2 * (a * T₂ + b * T₁)) % 5 = 0 := by
    rw [show a * T₁ - b * T₂ - 2 * (a * T₂ + b * T₁) = 5 * (u * T₁ - v * T₂) by
      linear_combination T₁ * hu - T₂ * hv]
    omega
  -- the expansions
  have hN6 : (5 * N) ≠ 0 := by positivity
  have e1 : (a * (N / 6 * E₁ + x₁) + b * (N / 6 * E₂ + x₂)) / (5 * N) =
      ((a * r₁ + b * r₂ : ℤ) : ℝ) / 6 + ((a * T₁ + b * T₂ : ℤ) : ℝ) / 5 + (a * x₁ + b * x₂) / (5 * N) := by
    rw [hT₁, hT₂]; push_cast; field_simp; ring
  have e2 : (a * (N / 6 * E₂ + x₂) - b * (N / 6 * E₁ + x₁)) / (5 * N) =
      ((a * r₂ - b * r₁ : ℤ) : ℝ) / 6 + ((a * T₂ - b * T₁ : ℤ) : ℝ) / 5 + (a * x₂ - b * x₁) / (5 * N) := by
    rw [hT₁, hT₂]; push_cast; field_simp; ring
  have e3 : (a * (N / 6 * E₁ + x₁) - b * (N / 6 * E₂ + x₂)) / (5 * N) =
      ((a * r₁ - b * r₂ : ℤ) : ℝ) / 6 + ((a * T₁ - b * T₂ : ℤ) : ℝ) / 5 + (a * x₁ - b * x₂) / (5 * N) := by
    rw [hT₁, hT₂]; push_cast; field_simp; ring
  have e4 : (a * (N / 6 * E₂ + x₂) + b * (N / 6 * E₁ + x₁)) / (5 * N) =
      ((a * r₂ + b * r₁ : ℤ) : ℝ) / 6 + ((a * T₂ + b * T₁ : ℤ) : ℝ) / 5 + (a * x₂ + b * x₁) / (5 * N) := by
    rw [hT₁, hT₂]; push_cast; field_simp; ring
  rw [e1] at c1
  rw [e2] at c2
  rw [e3] at c3
  rw [e4] at c4
  -- the divisibility `5 ∣ X₁, X₂`
  have hX : (a * T₁ + b * T₂) % 5 = 0 ∧ (a * T₁ - b * T₂) % 5 = 0 := by
    rcases htype with ⟨h1, h2, hx⟩ | ⟨h1, h1', h2, h2', rfl, rfl⟩
    · have hr1 : r₁ = 3 := by omega
      have hr2 : r₂ = 3 := by omega
      subst hr1 hr2
      obtain ⟨a₁, ha₁, ha₁', ha₁''⟩ := good_c (h := a + b) (X := a * T₁ + b * T₂)
        (s := (a * x₁ + b * x₂) / (5 * N)) (by omega) (by push_cast; ring) c1
      obtain ⟨a₂, ha₂, ha₂', ha₂''⟩ := good_c (h := a - b) (X := a * T₂ - b * T₁)
        (s := (a * x₂ - b * x₁) / (5 * N)) (by omega) (by push_cast; ring) c2
      obtain ⟨a₃, ha₃, ha₃', ha₃''⟩ := good_c (h := a - b) (X := a * T₁ - b * T₂)
        (s := (a * x₁ - b * x₂) / (5 * N)) (by omega) (by push_cast; ring) c3
      obtain ⟨a₄, ha₄, ha₄', ha₄''⟩ := good_c (h := a + b) (X := a * T₂ + b * T₁)
        (s := (a * x₂ + b * x₁) / (5 * N)) (by omega) (by push_cast; ring) c4
      have hp : (0 : ℝ) < (5 * N) ^ 2 := by positivity
      have hs1 : ((a * x₁ + b * x₂) / (5 * N)) ^ 2 + ((a * x₂ - b * x₁) / (5 * N)) ^ 2 < 1 / 18 := by
        rw [div_pow, div_pow, ← add_div, div_lt_iff₀ hp,
          show ((a : ℝ) * x₁ + b * x₂) ^ 2 + (a * x₂ - b * x₁) ^ 2 = (a ^ 2 + b ^ 2) * (x₁ ^ 2 + x₂ ^ 2) by ring,
          hsq]
        nlinarith
      have hs2 : ((a * x₁ - b * x₂) / (5 * N)) ^ 2 + ((a * x₂ + b * x₁) / (5 * N)) ^ 2 < 1 / 18 := by
        rw [div_pow, div_pow, ← add_div, div_lt_iff₀ hp,
          show ((a : ℝ) * x₁ - b * x₂) ^ 2 + (a * x₂ + b * x₁) ^ 2 = (a ^ 2 + b ^ 2) * (x₁ ^ 2 + x₂ ^ 2) by ring,
          hsq]
        nlinarith
      have k1 := five_dvd_of_small (Or.inl (by omega)) ha₁' ha₁'' ha₂' ha₂'' hs1
      have k2 := five_dvd_of_small (Or.inr (by omega)) ha₃' ha₃'' ha₄' ha₄'' hs2
      omega
    · have hr1 : r₁ = 2 ∨ r₁ = 4 := by omega
      have hr2 : r₂ = 2 ∨ r₂ = 4 := by omega
      simp only [mul_zero, add_zero, sub_zero, zero_div] at c1 c2 c3 c4
      obtain ⟨n₁, hn₁, hn₁'⟩ := good_q rfl c1
      obtain ⟨n₂, hn₂, hn₂'⟩ := good_q rfl c2
      obtain ⟨n₃, hn₃, hn₃'⟩ := good_q rfl c3
      obtain ⟨n₄, hn₄, hn₄'⟩ := good_q rfl c4
      have hab3' : ((a % 3 = 0 ∧ -b % 3 ≠ 0) ∨ (a % 3 ≠ 0 ∧ -b % 3 = 0)) := by omega
      have g1 := g_q hab3 hr1 hr2
      have g2 := g_q hab3' hr2 hr1
      have g3 := g_q hab3' hr1 hr2
      have g4 := g_q hab3 hr2 hr1
      simp only [neg_mul, ← sub_eq_add_neg] at g2 g3
      have k1 := q_pair g1 g2 hn₁ hn₁' hn₂ hn₂' (Or.inl m1)
      have k2 := q_pair g3 g4 hn₃ hn₃' hn₄ hn₄' (Or.inr m2)
      exact ⟨k1.1, k2.1⟩
  obtain ⟨α, hα⟩ : ∃ α, a = 5 * α + 3 := ⟨a / 5, by omega⟩
  obtain ⟨β, hβ⟩ : ∃ β, b = 5 * β + 4 := ⟨b / 5, by omega⟩
  have hT1 : T₁ % 5 = 0 := by
    have : (a * T₁ + b * T₂) + (a * T₁ - b * T₂) = 5 * (2 * α * T₁ + T₁) + T₁ := by
      linear_combination 2 * T₁ * hα
    exact five_aux hX.1 hX.2 this
  have hT2 : T₂ % 5 = 0 := by
    have : (a * T₁ + b * T₂) - (a * T₁ - b * T₂) = 5 * (2 * β * T₂ + T₂) + 3 * T₂ := by
      linear_combination 2 * T₂ * hβ
    exact five_aux' hX.1 hX.2 this
  obtain ⟨T₁', rfl⟩ : ∃ T, T₁ = 5 * T := ⟨T₁ / 5, by omega⟩
  obtain ⟨T₂', rfl⟩ : ∃ T, T₂ = 5 * T := ⟨T₂ / 5, by omega⟩
  refine ⟨r₁ + 6 * T₁', r₂ + 6 * T₂', ?_, ?_, ?_⟩
  · rw [hT₁]; push_cast; ring
  · rw [hT₂]; push_cast; ring
  · rcases htype with ⟨h1, h2, hx⟩ | ⟨h1, h1', h2, h2', hx1, hx2⟩
    · exact Or.inl ⟨by omega, by omega, hx⟩
    · exact Or.inr ⟨by omega, by omega, by omega, by omega, hx1, hx2⟩

lemma rep_succ {k : ℕ} {p q : ℝ} (hR : Rep k p q) (hS : inS (k + 1) p q) : Rep (k + 1) p q := by
  obtain ⟨E₁, E₂, x₁, x₂, rfl, rfl, htype⟩ := hR
  obtain ⟨c1, c2, c3, c4⟩ := hS (k + 1) le_rfl
  have h5N : (5 : ℝ) ^ (k + 1) = 5 * 5 ^ k := by rw [pow_succ]; ring
  rw [h5N] at c1 c2 c3 c4
  have hsq : ((A (k + 1) : ℤ) : ℝ) ^ 2 + ((B (k + 1) : ℤ) : ℝ) ^ 2 = (5 * 5 ^ k) ^ 2 := by
    have h' : ((A (k + 1) : ℝ)) ^ 2 + (B (k + 1) : ℝ) ^ 2 = 25 ^ (k + 1) := by exact_mod_cast AB_sq (k + 1)
    rw [h', show (25 : ℝ) = 5 ^ 2 by norm_num, ← pow_mul]
    ring
  obtain ⟨E₁', E₂', h1, h2, h3⟩ := step_core (by positivity) (AB_mod2 (k + 1)) (AB_mod5 (k + 1) (by omega))
    (AB_mod3 (k + 1)) hsq htype c1 c2 c3 c4
  refine ⟨E₁', E₂', x₁, x₂, ?_, ?_, h3⟩
  · rw [h1, h5N]
  · rw [h2, h5N]

theorem rep_of_inS (k : ℕ) {p q : ℝ} (h : inS k p q) : Rep k p q := by
  induction k with
  | zero => exact rep_zero h
  | succ k ih => exact rep_succ (ih fun j hj => h j (by omega)) h

/-! ## Exact relations and the arithmetic contradictions -/

lemma exact3 {N x₁ x₂ x₃ : ℝ} {m₁ m₂ m₃ E₁ E₂ E₃ : ℤ} (hx₁ : |x₁| ≤ 1 / 4) (hx₂ : |x₂| ≤ 1 / 4)
    (hx₃ : |x₃| ≤ 1 / 4) (hN : 2 * (|(m₁ : ℝ)| + |(m₂ : ℝ)| + |(m₃ : ℝ)|) < N)
    (h : (m₁ : ℝ) * (N / 6 * E₁ + x₁) + m₂ * (N / 6 * E₂ + x₂) + m₃ * (N / 6 * E₃ + x₃) = 0) :
    m₁ * E₁ + m₂ * E₂ + m₃ * E₃ = 0 := by
  by_contra hne
  have h1 : (1 : ℝ) ≤ |((m₁ * E₁ + m₂ * E₂ + m₃ * E₃ : ℤ) : ℝ)| := by
    rw [← Int.cast_abs]; exact_mod_cast Int.one_le_abs hne
  have h2 : N / 6 * ((m₁ * E₁ + m₂ * E₂ + m₃ * E₃ : ℤ) : ℝ) = -(m₁ * x₁ + m₂ * x₂ + m₃ * x₃) := by
    push_cast; linear_combination h
  have h3 : |(m₁ : ℝ) * x₁ + m₂ * x₂ + m₃ * x₃| ≤ (|(m₁ : ℝ)| + |(m₂ : ℝ)| + |(m₃ : ℝ)|) / 4 := by
    calc |(m₁ : ℝ) * x₁ + m₂ * x₂ + m₃ * x₃| ≤ |(m₁ : ℝ) * x₁| + |(m₂ : ℝ) * x₂| + |(m₃ : ℝ) * x₃| :=
          abs_add_three _ _ _
      _ = |(m₁ : ℝ)| * |x₁| + |(m₂ : ℝ)| * |x₂| + |(m₃ : ℝ)| * |x₃| := by simp only [abs_mul]
      _ ≤ _ := by
          have := abs_nonneg (m₁ : ℝ); have := abs_nonneg (m₂ : ℝ); have := abs_nonneg (m₃ : ℝ)
          nlinarith
  have hsum : 0 ≤ |(m₁ : ℝ)| + |(m₂ : ℝ)| + |(m₃ : ℝ)| := by positivity
  have hNpos : 0 < N := by linarith
  have h4 : |N / 6 * ((m₁ * E₁ + m₂ * E₂ + m₃ * E₃ : ℤ) : ℝ)| =
      N / 6 * |((m₁ * E₁ + m₂ * E₂ + m₃ * E₃ : ℤ) : ℝ)| := by
    rw [abs_mul, abs_of_pos (by positivity : 0 < N / 6)]
  rw [h2, abs_neg] at h4
  nlinarith

/-- The first coordinate of a point of `S_N`, and its type. -/
def typed (E : ℤ) : Prop := E % 6 = 3 ∨ (E % 2 = 0 ∧ E % 3 ≠ 0)

lemma even_iff_of_typed {E : ℤ} (h : typed E) : Even E ↔ ¬ (3 : ℤ) ∣ E := by
  rw [Int.even_iff]; unfold typed at h; omega

/-- The arithmetic contradiction, for `d ≡ 23 (mod 24)`. -/
lemma contra_23 {d a b E₀ E₁ E₂ : ℤ} (hd : d % 24 = 23) (ha : 4 * a = d + 1) (hb : 2 * b = d - 1)
    (h0 : typed E₀) (r1 : a * E₁ + a * E₂ + b * E₀ = 0) : False := by
  obtain ⟨a', rfl⟩ : ∃ a', a = 6 * a' := ⟨a / 6, by omega⟩
  have h2 : (2 : ℤ) ∣ b * E₀ := ⟨-(3 * a' * (E₁ + E₂)), by linear_combination r1⟩
  have h3 : (3 : ℤ) ∣ b * E₀ := ⟨-(2 * a' * (E₁ + E₂)), by linear_combination r1⟩
  rcases Int.prime_two.dvd_or_dvd h2 with h | h <;> rcases Int.prime_three.dvd_or_dvd h3 with h' | h' <;>
    unfold typed at h0 <;> omega

/-- The arithmetic contradiction, for `d ≡ 11 (mod 24)`. Here `d + 1 = t m` with `t = 3^s` and `3 ∤ m`, and
`n = 1 + 6t`. -/
lemma contra_11 {d a b t m E₀ E₁ E₂ E₃ : ℤ} (hd : d % 24 = 11) (ha : 4 * a = d + 1) (hb : 2 * b = d - 1)
    (ht : t % 2 = 1) (ht0 : t ≠ 0) (hm : m % 3 ≠ 0) (hdt : d + 1 = t * m)
    (h0 : typed E₀) (h1 : typed E₁) (h3 : typed E₃) (r1 : a * E₁ + a * E₂ + b * E₀ = 0)
    (r2 : ((1 + 6 * t) ^ 2 + d) * E₃ + (-(((1 + 6 * t) - 1) * ((1 + 6 * t) + d))) * E₀ +
      (-((1 + 6 * t) * (1 + d))) * E₁ = 0) : False := by
  -- step 1: `3 ∣ E₀`, so `E₀` is odd
  have h3a : (3 : ℤ) ∣ b * E₀ := ⟨-(a / 3) * (E₁ + E₂), by
    obtain ⟨a', rfl⟩ : ∃ a', a = 3 * a' := ⟨a / 3, by omega⟩
    rw [show 3 * a' / 3 = a' by omega]; linear_combination r1⟩
  have hE₀3 : (3 : ℤ) ∣ E₀ := by
    rcases Int.prime_three.dvd_or_dvd h3a with h | h
    · omega
    · exact h
  have hE₀ : ¬ Even E₀ := by rw [even_iff_of_typed h0]; tauto
  -- step 3 (at 3)
  have r3 : (12 + 36 * t + m) * E₃ = 3 * (2 * t * (6 + m)) * E₀ + ((1 + 6 * t) * m) * E₁ := by
    have : t * ((12 + 36 * t + m) * E₃ - 3 * (2 * t * (6 + m)) * E₀ - ((1 + 6 * t) * m) * E₁) = 0 := by
      have hd' : d = t * m - 1 := by linarith
      subst hd'
      linear_combination r2
    have := (mul_eq_zero.1 this).resolve_left ht0
    linear_combination this
  have k3 : (3 : ℤ) ∣ E₃ ↔ (3 : ℤ) ∣ E₁ := by
    have hP : ¬ (3 : ℤ) ∣ 12 + 36 * t + m := by omega
    have hR : ¬ (3 : ℤ) ∣ (1 + 6 * t) * m := by
      intro h
      rcases Int.prime_three.dvd_or_dvd h with h | h <;> omega
    constructor
    · intro h
      have : (3 : ℤ) ∣ ((1 + 6 * t) * m) * E₁ := by
        have h' : (3 : ℤ) ∣ (12 + 36 * t + m) * E₃ := dvd_mul_of_dvd_right h _
        rw [r3] at h'
        exact (dvd_add_right (dvd_mul_of_dvd_left (dvd_mul_right 3 _) _)).1 h'
      exact (Int.prime_three.dvd_or_dvd this).resolve_left hR
    · intro h
      have : (3 : ℤ) ∣ (12 + 36 * t + m) * E₃ := by
        rw [r3]
        exact dvd_add (dvd_mul_of_dvd_left (dvd_mul_right 3 _) _) (dvd_mul_of_dvd_right h _)
      exact (Int.prime_three.dvd_or_dvd this).resolve_left hP
  -- step 2 (at 2)
  obtain ⟨τ, rfl⟩ : ∃ τ, t = 2 * τ + 1 := ⟨t / 2, by omega⟩
  obtain ⟨δ, rfl⟩ : ∃ δ, d = 8 * δ + 3 := ⟨d / 8, by omega⟩
  have r4 : (36 * τ ^ 2 + 42 * τ + 2 * δ + 13) * E₃ =
      ((6 * τ + 3) * (6 * τ + 5 + 4 * δ)) * E₀ + ((12 * τ + 7) * (2 * δ + 1)) * E₁ := by
    apply mul_left_cancel₀ (show (4 : ℤ) ≠ 0 by norm_num)
    linear_combination r2
  have oP : ¬ Even (36 * τ ^ 2 + 42 * τ + 2 * δ + 13) := by
    rw [Int.not_even_iff_odd]; exact ⟨18 * τ ^ 2 + 21 * τ + δ + 6, by ring⟩
  have oQ : ¬ Even ((6 * τ + 3) * (6 * τ + 5 + 4 * δ)) := by
    rw [Int.even_mul]; rintro (h | h) <;> rw [Int.even_iff] at h <;> omega
  have oR : ¬ Even ((12 * τ + 7) * (2 * δ + 1)) := by
    rw [Int.even_mul]; rintro (h | h) <;> rw [Int.even_iff] at h <;> omega
  have e3 : Even ((36 * τ ^ 2 + 42 * τ + 2 * δ + 13) * E₃) ↔ Even E₃ := by
    rw [Int.even_mul]; exact ⟨fun h => h.resolve_left oP, Or.inr⟩
  have e0 : ¬ Even (((6 * τ + 3) * (6 * τ + 5 + 4 * δ)) * E₀) := by
    rw [Int.even_mul]; exact fun h => h.elim oQ hE₀
  have e1 : Even (((12 * τ + 7) * (2 * δ + 1)) * E₁) ↔ Even E₁ := by
    rw [Int.even_mul]; exact ⟨fun h => h.resolve_left oR, Or.inr⟩
  rw [r4, Int.even_add, e1] at e3
  rw [even_iff_of_typed h1, even_iff_of_typed h3] at e3
  by_cases hc : (3 : ℤ) ∣ E₁
  · exact (e3.1 ⟨fun h => absurd h e0, fun h => absurd hc h⟩) (k3.2 hc)
  · have := e3.2 (k3.1.mt hc)
    exact hc (by_contra fun h' => e0 (this.2 h'))

end FourColours
