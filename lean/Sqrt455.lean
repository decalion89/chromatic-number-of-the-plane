import QuadraticPlanes

/-!
# The plane over `ℚ(√455)` has chromatic number 4

**Theorem** (`Sqrt455.chromaticNumber_eq_four`). The unit-distance graph on `ℚ(√455)²` has chromatic number 4
(`notes/quadratic_planes.md`). This file is written from the data by `tools/field_lean.py`.

*Lower bound* (`not_colorable_three`). The graph of `data/quadratic_planes/q455.json` has 74 vertices `[a, b, c, e]`,
standing for `((a + b√455)/780, (c + e√455)/780)`, and 161 edges. The kernel checks that each edge has length 1
(`checkEdges_E`, `QuadraticPlanes.adj_pt`), and Mathlib's `lrat_proof` checks the LRAT proof
`data/quadratic_planes/q455.lrat` that the colouring formula `q455.cnf` (222 variables, 559 clauses) is
unsatisfiable. A 3-colouring would satisfy it (`QuadraticPlanes.clause_facts`).

*Upper bound* (`QuadraticPlanes.colorable_four_ramified`): `455 = 7 · 65`, so Moorhouse's reduction at 7 applies (the ramified case).
-/

namespace Sqrt455

open LocalColouring QuadraticPlanes

set_option maxHeartbeats 4000000 in
/-- The 74 points of `data/quadratic_planes/q455.json`, as `[a, b, c, e]` with denominator 780. -/
def points : List (ℤ × ℤ × ℤ × ℤ) := [
    (-1430, 0, 390, 0), (-1430, 6, 0, 0), (-1300, -4, -260, 20), (-1170, 0, 0, -18),
    (-1170, 0, 0, 18), (-1170, 6, 390, 18), (-1040, -2, -130, 16), (-715, -3, -195, 11),
    (-715, -3, 195, -11), (-715, -3, 585, 11), (-715, 3, -195, -11), (-715, 3, 195, 11),
    (-715, 3, 585, -11), (-585, -7, -455, 9), (-585, -7, 455, -9), (-585, 7, -455, -9),
    (-585, 7, 455, 9), (-455, -9, -585, 7), (-455, 9, 585, 7), (-130, -10, -650, 2),
    (-130, -10, 650, -2), (-130, -4, 260, -2), (-130, 2, -130, -2), (-130, 16, -1040, -2),
    (0, -6, -390, 0), (0, -6, 390, 0), (0, 0, 0, 0), (0, 0, 780, 0),
    (0, 6, -390, 0), (0, 6, 390, 0), (130, -16, -1040, -2), (130, -10, -650, -2),
    (130, -4, 260, 2), (130, 4, -260, 2), (130, 10, 650, -2), (130, 16, 1040, -2),
    (260, 6, 390, -4), (390, -4, -260, -6), (585, -7, -455, -9), (585, -7, 455, 9),
    (585, 7, 455, -9), (715, -3, -195, -11), (715, -3, 195, 11), (715, 3, -195, 11),
    (715, 3, 195, -11), (1170, -6, 390, 18), (1170, 0, 0, -18), (1170, 0, 0, 18),
    (1430, 6, 390, -22), (-715, 9, 585, 11), (-715, 9, -195, 11), (-845, -13, -845, 13),
    (845, -13, -845, -13), (-455, 3, -195, -7), (325, -7, 455, 5), (455, 3, -195, 7),
    (585, -13, 845, 9), (-1885, 3, 195, 29), (845, 13, 845, -13), (-455, 3, 195, 7),
    (325, 7, -455, 5), (455, -3, -195, -7), (-715, 9, -585, -11), (715, 9, 585, -11),
    (-715, 9, 195, -11), (0, 12, 780, 0), (260, 0, 0, -4), (-455, -3, -195, 7),
    (-260, 0, 0, -4), (455, -3, 195, 7), (715, -9, -585, -11), (975, 3, 195, -15),
    (0, 12, -780, 0), (0, 12, 0, 0)]

/-- The point with index `i`. -/
def P (i : Fin 74) : ℤ × ℤ × ℤ × ℤ := points.getD i.val (0, 0, 0, 0)

set_option maxHeartbeats 4000000 in
/-- The 161 edges of `data/quadratic_planes/q455.json`: all the unit distances among the points. -/
def E : List (Fin 74 × Fin 74) := [
    (0, 8), (0, 9), (0, 11), (0, 12), (1, 50), (1, 10), (1, 11), (1, 64), (2, 51), (2, 57),
    (2, 11), (2, 13), (3, 15), (3, 53), (3, 62), (3, 14), (4, 16), (4, 59), (4, 49), (4, 57),
    (4, 67), (4, 13), (5, 18), (5, 59), (5, 7), (5, 57), (6, 51), (6, 16), (6, 17), (7, 24),
    (7, 26), (7, 19), (7, 9), (8, 25), (8, 20), (8, 26), (9, 27), (9, 25), (10, 26), (10, 28),
    (10, 21), (10, 12), (11, 26), (11, 29), (12, 27), (12, 29), (13, 26), (13, 31), (14, 26), (14, 32),
    (14, 22), (15, 26), (15, 33), (15, 23), (16, 26), (16, 34), (17, 26), (17, 30), (18, 26), (18, 35),
    (18, 36), (19, 51), (19, 38), (19, 61), (20, 69), (20, 56), (20, 39), (21, 55), (21, 39), (22, 54),
    (22, 62), (23, 62), (23, 60), (24, 70), (24, 59), (24, 25), (24, 41), (25, 55), (25, 56), (25, 53),
    (25, 42), (26, 38), (26, 39), (26, 40), (26, 27), (26, 41), (26, 42), (26, 43), (26, 44), (28, 69),
    (28, 50), (28, 43), (28, 62), (28, 29), (29, 44), (29, 49), (29, 63), (29, 67), (29, 64), (29, 61),
    (30, 70), (30, 52), (30, 38), (31, 41), (31, 52), (31, 67), (32, 43), (32, 56), (32, 53), (32, 60),
    (33, 54), (33, 42), (34, 44), (34, 58), (34, 59), (35, 40), (35, 63), (35, 58), (36, 58), (36, 59),
    (36, 41), (36, 71), (37, 52), (37, 40), (37, 71), (38, 46), (39, 47), (40, 46), (43, 45), (44, 48),
    (45, 69), (45, 56), (46, 70), (46, 63), (46, 61), (47, 69), (47, 55), (48, 58), (48, 63), (49, 50),
    (49, 65), (50, 73), (53, 72), (54, 68), (55, 72), (55, 68), (59, 65), (59, 66), (60, 68), (62, 72),
    (62, 64), (62, 68), (63, 65), (63, 66), (64, 73), (65, 73), (66, 70), (66, 67), (66, 71), (68, 69),
    (72, 73)]

lemma checkEdges_E : checkEdges 455 780 P E = true := by decide +kernel

lemma noLoops_E : noLoops E = true := by decide +kernel

-- `q455.cnf` is unsatisfiable: for all propositions `x₀, …, x_221`, some clause is false. The statement is
-- a disjunction, over the clauses, of the negations of the clauses.
lrat_proof refuted
  (include_str "../data/quadratic_planes/q455.cnf")
  (include_str "../data/quadratic_planes/q455.lrat")

set_option maxHeartbeats 10000000 in
/-- The graph of `q455.json` is not 3-colourable. -/
theorem graph_not_colorable : ¬ (edgeGraph E).Colorable 3 := by
  rintro ⟨C⟩
  obtain ⟨V, hVert, hEdge, hVu, hVv⟩ := clause_facts noLoops_E (u₀ := 0) (v₀ := 8) (List.mem_cons_self ..) C
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
    (V 204) (V 205) (V 206) (V 207) (V 208) (V 209) (V 210) (V 211) (V 212) (V 213) (V 214) (V 215)
    (V 216) (V 217) (V 218) (V 219) (V 220) (V 221)
  casesm* _ ∨ _
  all_goals first
    | exact hVert _ ‹_› (by norm_num) (by norm_num)
    | exact hEdge _ _ ‹_› (by decide +kernel)
    | exact ‹¬V 0› hVu
    | exact ‹¬V 25› hVv

/-- **Lower bound.** The unit-distance graph of `ℚ(√455)²` is not 3-colourable. -/
theorem not_colorable_three : ¬ (unitDistGraph (L 455)).Colorable 3 := fun h =>
  graph_not_colorable (h.of_hom (edgeGraph.hom (fun i => pt 455 780 (P i)) (pt_adj (by norm_num) checkEdges_E)))

/-- **Theorem.** The unit-distance graph of `ℚ(√455)²` has chromatic number 4. -/
theorem chromaticNumber_eq_four : (unitDistGraph (L 455)).chromaticNumber = 4 :=
  chromaticNumber_eq_four_of (colorable_four_ramified (d' := 65) (by norm_num) (by norm_num)) not_colorable_three

end Sqrt455
