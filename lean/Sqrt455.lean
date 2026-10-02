import QuadraticPlanes

/-!
# The plane over `ℚ(√455)` has chromatic number 4

**Theorem** (`Sqrt455.chromaticNumber_eq_four`). The unit-distance graph on `ℚ(√455)²` has chromatic number 4
(`notes/quadratic_planes.md`). This file is written from the data by `tools/field_lean.py`.

*Lower bound* (`not_colorable_three`). The graph of `data/quadratic_planes/q455.json` has 71 vertices `[a, b, c, e]`,
standing for `((a + b√455)/780, (c + e√455)/780)`, and 150 edges. The kernel checks that each edge has length 1
(`checkEdges_E`, `QuadraticPlanes.adj_pt`), and Mathlib's `lrat_proof` checks the LRAT proof
`data/quadratic_planes/q455.lrat` that the colouring formula `q455.cnf` (213 variables, 523 clauses) is
unsatisfiable. A 3-colouring would satisfy it (`QuadraticPlanes.clause_facts`).

*Upper bound* (`QuadraticPlanes.colorable_four_ramified`): `455 = 7 · 65`, so Moorhouse's reduction at 7 applies (the ramified case).
-/

namespace Sqrt455

open LocalColouring QuadraticPlanes

set_option maxHeartbeats 4000000 in
/-- The 71 points of `data/quadratic_planes/q455.json`, as `[a, b, c, e]` with denominator 780. -/
def points : List (ℤ × ℤ × ℤ × ℤ) := [
    (-1430, -6, -390, 22), (-1430, 0, -390, 0), (-1430, 0, 390, 0), (-1430, 6, 0, 0),
    (-1430, 6, 390, 22), (-1300, -4, -260, 20), (-1170, 0, 0, 18), (-1170, 6, 390, 18),
    (-1040, -2, -130, 16), (-1040, 2, 130, 16), (-715, -3, -195, 11), (-715, -3, 195, -11),
    (-715, -3, 585, 11), (-715, 3, -195, -11), (-715, 3, 195, 11), (-715, 3, 585, -11),
    (-585, -7, -455, 9), (-585, 7, 455, 9), (-455, -9, -585, 7), (-455, 9, 585, 7),
    (-260, 6, -390, -4), (-130, -10, -650, 2), (-130, -10, 650, -2), (-130, -4, 260, -2),
    (-130, 2, -130, -2), (-130, 10, -650, -2), (-130, 16, -1040, -2), (-130, 16, 1040, 2),
    (0, -6, -390, 0), (0, -6, 390, 0), (0, 0, 0, 0), (0, 6, -390, 0),
    (0, 6, 390, 0), (130, -16, -1040, -2), (130, -10, -650, -2), (130, -4, -260, -2),
    (130, -4, 260, 2), (130, 2, 130, -2), (455, 9, -585, 7), (585, -7, -455, -9),
    (585, -7, 455, 9), (715, -3, -195, -11), (715, 3, -195, 11), (715, 3, 195, -11),
    (1430, 6, 390, -22), (-845, 13, 845, 13), (-715, 9, 585, 11), (-715, 9, -195, 11),
    (-845, -13, -845, 13), (-455, 3, -195, -7), (455, 3, -195, 7), (585, -13, 845, 9),
    (-1885, 3, 195, 29), (-455, 3, 195, 7), (325, 7, -455, 5), (585, -1, 65, 9),
    (-715, 9, -585, -11), (715, 9, 585, -11), (845, -1, -65, -13), (-715, 9, 195, -11),
    (-585, 1, -65, -9), (-845, -1, 65, -13), (0, 12, 780, 0), (260, 0, 0, -4),
    (-260, 0, 0, -4), (455, -3, 195, 7), (715, -9, -585, -11), (845, -7, -455, -13),
    (-1040, 10, -650, -16), (0, 12, -780, 0), (0, 12, 0, 0)]

/-- The point with index `i`. -/
def P (i : Fin 71) : ℤ × ℤ × ℤ × ℤ := points.getD i.val (0, 0, 0, 0)

set_option maxHeartbeats 4000000 in
/-- The 150 edges of `data/quadratic_planes/q455.json`: all the unit distances among the points. -/
def E : List (Fin 71 × Fin 71) := [
    (0, 52), (0, 10), (0, 48), (1, 2), (1, 10), (1, 13), (2, 11), (2, 12), (2, 14), (2, 15),
    (3, 59), (3, 13), (3, 14), (3, 47), (4, 14), (4, 45), (4, 46), (5, 14), (5, 52), (5, 16),
    (5, 48), (6, 16), (6, 46), (6, 17), (6, 52), (6, 53), (7, 52), (7, 53), (7, 10), (7, 19),
    (8, 17), (8, 48), (8, 18), (9, 45), (9, 16), (9, 19), (10, 28), (10, 30), (10, 12), (10, 21),
    (11, 29), (11, 30), (11, 20), (11, 22), (12, 29), (13, 30), (13, 31), (13, 15), (13, 23), (13, 25),
    (14, 30), (14, 32), (15, 32), (16, 30), (16, 34), (16, 35), (17, 30), (17, 27), (18, 30), (18, 33),
    (19, 37), (19, 30), (20, 50), (20, 38), (20, 61), (21, 39), (21, 48), (22, 65), (22, 51), (22, 40),
    (23, 61), (23, 40), (23, 50), (23, 55), (24, 61), (24, 38), (24, 55), (24, 56), (25, 50), (25, 60),
    (26, 38), (26, 54), (26, 56), (27, 45), (27, 46), (28, 29), (28, 66), (28, 41), (28, 53), (29, 60),
    (29, 49), (29, 50), (29, 51), (30, 38), (30, 39), (30, 40), (30, 41), (30, 42), (30, 43), (31, 32),
    (31, 65), (31, 42), (31, 47), (31, 55), (31, 56), (32, 59), (32, 43), (32, 46), (32, 57), (33, 66),
    (33, 39), (34, 67), (34, 41), (35, 58), (35, 67), (35, 43), (35, 53), (36, 42), (36, 49), (36, 51),
    (36, 54), (37, 58), (37, 39), (37, 57), (43, 44), (44, 58), (44, 57), (46, 62), (46, 47), (47, 70),
    (49, 68), (49, 69), (50, 64), (50, 69), (53, 62), (53, 63), (54, 64), (56, 59), (56, 64), (56, 69),
    (57, 62), (57, 63), (59, 70), (60, 68), (61, 68), (62, 70), (63, 66), (63, 67), (64, 65), (69, 70)]

lemma checkEdges_E : checkEdges 455 780 P E = true := by decide +kernel

lemma noLoops_E : noLoops E = true := by decide +kernel

-- `q455.cnf` is unsatisfiable: for all propositions `x₀, …, x_212`, some clause is false. The statement is
-- a disjunction, over the clauses, of the negations of the clauses.
lrat_proof refuted
  (include_str "../data/quadratic_planes/q455.cnf")
  (include_str "../data/quadratic_planes/q455.lrat")

set_option maxHeartbeats 10000000 in
/-- The graph of `q455.json` is not 3-colourable. -/
theorem graph_not_colorable : ¬ (edgeGraph E).Colorable 3 := by
  rintro ⟨C⟩
  obtain ⟨V, hVert, hEdge, hVu, hVv⟩ := clause_facts noLoops_E (u₀ := 0) (v₀ := 52) (List.mem_cons_self ..) C
  have H := refuted
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
    (V 204) (V 205) (V 206) (V 207) (V 208) (V 209) (V 210) (V 211) (V 212)
  casesm* _ ∨ _
  all_goals first
    | exact hVert _ ‹_› (by norm_num) (by norm_num)
    | exact hEdge _ _ ‹_› (by decide +kernel)
    | exact ‹¬V 0› hVu
    | exact ‹¬V 157› hVv

/-- **Lower bound.** The unit-distance graph of `ℚ(√455)²` is not 3-colourable. -/
theorem not_colorable_three : ¬ (unitDistGraph (L 455)).Colorable 3 := fun h =>
  graph_not_colorable (h.of_hom (edgeGraph.hom (fun i => pt 455 780 (P i)) (pt_adj (by norm_num) checkEdges_E)))

/-- **Theorem.** The unit-distance graph of `ℚ(√455)²` has chromatic number 4. -/
theorem chromaticNumber_eq_four : (unitDistGraph (L 455)).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of (colorable_four_ramified (d' := 65) (by norm_num) (by norm_num)) not_colorable_three

end Sqrt455
