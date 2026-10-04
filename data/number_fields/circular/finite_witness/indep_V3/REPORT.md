# Referee report: the nine-point witnesses for `χ_c = 3`

Written by a referee working separately, with programs of its own (this folder), on 4 October 2026 (night).
**Verdict:** every mathematical claim under review is correct; no error was found in the points, the graphs, the
values of `χ_c`, the proof by hand, the `K_{8/3}` argument, the criticality certificates, the bound for 8 vertices,
the classification for 9 vertices or the SAT records, and the reproduction over `ℚ(√7)` is exact. The findings
(Part A) concern reproducibility, wording and missed opportunities; how each was applied is listed at the end.

**Object:** `data/number_fields/circular/finite_witness/` at commit 17347f14 (`witness_q7.json.gz`,
`witness_q31.json.gz`, `q7_seed.json`, `check_small.py`, `small_triangle_free.py`, `nine_vertices.py`, the generalised
`grow.py`, `minimise.py`, `critical.py`, `README.md`, `verification.txt`), the value-3 paragraph of Section 10 of
`papers/three-colours/three-colours.tex`, the sentences of the introduction and of Question 3 on it, and §6.8 of
`notes/circular_planes.md`. The checkout was read only. Programs: `r1_points.py`, `r2_graphs.py`, `r3_cnf.py`,
`r4_repro_compare.py`, `rec_run.py`, `r5*.py`, `r6_fields.py`, `brute83.c`, `tf_small.c`, `tf_burnside.c`; outputs
in `results/`. The scripts that read `../repo_copy/` expect a copy of the repository's `finite_witness/` programs
and witness files there (as at 17347f14); `cycles31_classes.json`, read by `r5c` and `r5d`, is written by
`r5b_cycles31.py` (not stored, 221 kB).

## Part A: summary and findings

Every claim was confirmed by methods that share no code with the repository: brute force over all `2^28` labelled
graphs on 8 vertices; a Burnside orbit count for 9 vertices with no isomorphism testing; an exhaustive `8^8`
homomorphism search; a SAT encoding and a Kahn-algorithm tight-cycle test of its own; explicit homomorphisms
`N·c + pos` for the criticality certificates; a clause-by-clause check of the formulas.

1. **The `ℚ(√31)` growth cannot be reproduced from its description.** The README says "from 0, the unit vectors with
   denominator 80 and a 5-cycle; 367 vertices", but no seed is stored and the 5-cycle is not given. The 5-cycle
   `0, (−1, 0), (−2, 0), Q, R` of the `ℚ(√7)` seed does not exist over `ℚ(√31)` with denominator 80; there are
   47 040 directed 5-cycles through 0 with two extra points (2 972 classes up to symmetry). Two natural seeds (the two
   5-cycles of the stored witness through 0 and `(−1, 0)`) give 359 vertices, not 367; after deletion they give
   exactly the stored nine points, with the same coordinates. In 40 random seed classes the grown sizes ranged from
   168 to 559, never 367. *Fix:* store the seed with the exact command and output, as for `ℚ(√7)`, or drop "367".
2. **The `ℚ(√15)`/`ℚ(√39)` remark is unverifiable and depends on the seed.** "Witnesses with 22 and 25 vertices (not
   stored)" gives no seed, denominator or files. With denominator 80 and 12 random 5-cycle seeds each, the same
   programs gave 13–25 vertices over `ℚ(√15)` and 17–44 over `ℚ(√39)` (pipeline outputs, not certified). *Fix:* drop
   it, or say that the sizes depend on the seed and that the least size there is open.
3. **All three nine-vertex graphs are realisable; the text suggests otherwise** (a missed opportunity, not an error).
   Verified exactly, as induced unit-distance graphs, each passing all five checks of `check_small.py`:
   `M` itself (12 edges) over `ℚ(√31)`, `D = 80`: `(−88,0,64,8)`, `(−60,3,45,4)`, `(−60,8,120,4)`, `(−36,4,12,12)`,
   `(0,0,0,0)`, `(48,0,64,0)`, `(20,8,120,4)`, `(−36,4,92,12)`, `(−20,−2,30,4)`; `H₇` over `ℚ(√31)`, `D = 80`:
   `(−80,0,0,0)`, `(−68,−2,−34,4)`, `(−68,2,34,4)`, `(0,0,0,0)`, `(−12,−2,−34,−4)`, `(−136,0,0,8)`, `(−40,−6,−22,0)`,
   `(−40,6,22,0)`, `(−108,−8,−56,4)`; `M + m₀m₁ + m₁m₂` over `ℚ(√7)`, `D = 160`: `(−54,−40,−72,30)`, `(−2,−19,−111,58)`,
   `(0,0,0,0)`, `(70,−49,−165,18)`, `(91,−7,13,49)`, `(124,−9,−93,−12)`, `(52,21,−39,28)`, `(163,−37,−41,9)`,
   `(39,−28,52,21)` (point `[a, b, c, e]` is `((a + b√d)/D, (c + e√d)/D)`). The phrase "`M + m₀m₁ + m₁m₂`, which is
   the unit-distance graph of nine points of `ℚ(√31)²`" is true but reads as if that graph were special to `ℚ(√31)`.
   *Fix:* say that all three graphs occur; consider storing the edge-minimal `M`, to which the proof by hand applies
   directly.
4. **Citation for "numerator at most n".** Zhu's survey is acceptable as a secondary source, but the result is due to
   A. Vince (J. Graph Theory 12 (1988) 551–559) and J. A. Bondy and P. Hell (J. Graph Theory 14 (1990) 479–482);
   `check_small.py`'s docstring names them, the bibliography does not. *Fix:* add both.
5. **Question 3 wording.** "For the value 3 in a plane without unit triangles the answer is nine, over `ℚ(√7)`" can be
   read as a claim about every such plane. What is proved: at least nine in every plane without unit triangles, and
   exactly nine over `ℚ(√7)` (and `ℚ(√31)`); for `ℚ(√15)`, `ℚ(√39)` and others the least number is open. *Fix:*
   rephrase.
6. **Clarity of the proof by hand** (it is correct). The "six 6-cycles" `C_k = Z_{k+1} − Z_k` are three 6-cycles in
   both orientations (`C_{k+3} = −C_k` as 1-chains), and only `C₀, C₁, C₂` are needed, as `Z₃ − Z₀ = C₀ + C₁ + C₂`.
   The stored cycle list of `H₇` contains `C₁`, `C₂` and the hexagon but not `C₀ = P₀P₁P₂P₆P₅P₄`, which in `H₇` has
   the chord `P₁P₅` and is never tight; a sentence linking the proof to the listed cycles would help. "Already this
   subgraph `M`" would read better as "the subgraph `M = H₇ − P₁P₅`".
7. **Small wording points.** "Translated so that one of them is 0": the nine points as found already contain 0 (the
   image of `P₆`); the stored set is the found set translated by `(1, 0)`. The stored denominator of `H₇` is 160,
   although every coordinate has a denominator dividing 8. The note's passage heading cites `check_small.py` and
   `small_triangle_free.py` but not `nine_vertices.py`. `verification.txt` gives no size or core figures for
   `q31.cnf` (here: 47 variables, 115 clauses; 115 of 115 clauses and 46 of 144 lemmas in the core).
8. **`minimise.py` for `(7, 2)`.** The only behavioural change since 5cc8ae84 is that the first and last solver calls
   run without a conflict budget; the clause streams and outputs are otherwise identical. The README could say so.

## Part B: claim by claim

**Claim 1 (points, unit edges, induced, printed coordinates): confirmed.** `r1_points.py` uses exact `Fraction`s in
`ℚ(√d)` and none of the repository's code. `witness_q7` (`d = 7`, `D = 160`): 9 distinct points, 13 listed edges, 13
unit pairs, equal sets, so the graph is induced; no triangles. `witness_q31` (`d = 31`, `D = 80`): 9 distinct points,
14 edges, 14 unit pairs, equal. The paper's `P₀, …, P₇, S`, typed from the TeX source, equal the file's vertices
0–8; the README table and the note give the same points; the README's `ℚ(√31)` coordinates equal the file in order.

**Claim 2 (structure): confirmed** (`r2_graphs.py`). `H₇`'s edge set is exactly the 8-cycle `P₀⋯P₇`, the chords
`P₀P₄`, `P₁P₅`, `P₂P₆` and the path `P₃SP₇`; subdividing any one of the four chords of the Wagner graph gives a graph
isomorphic to `H₇`. `H₇ = M + m₀m₁` with `v = P₀, P₄, P₃, P₂, P₆, P₇` and `m = P₁, P₅, S`, as identical edge sets;
the `ℚ(√31)` graph is `M + m₀m₁ + m₁m₂` with `v = 0, 1, 8, 7, 2, 6` and `m = 3, 4, 5`. `M` plus any one midpoint edge
is isomorphic to `H₇`, plus any two to the `ℚ(√31)` graph. `M` is `K_{3,3}` with a perfect matching subdivided.
`|Aut| = 12, 4, 8` for `M`, `H₇` and the `ℚ(√31)` graph.

**Claim 3 (`χ_c = 3`, proof by hand, Lemma 20): confirmed.** Over all `3⁹` maps: 84 proper 3-colourings of `H₇` and
48 of the `ℚ(√31)` graph, each with a tight cycle (its own Kahn test) and a tight listed cycle. SAT finds no
homomorphism to any `K_{p/q}` with `2 ≤ p/q < 3` and `p ≤ 30`. Proof by hand: every `Z_k` is a 5-cycle of `M` and
every `C_k` a 6-cycle; `Z_{k+1} − Z_k = C_k` as 1-chains for `k = 0, …, 5` (indices of `v` modulo 6, of `m` modulo
3), `Z₆ = Z₀`, `Z₀ + Z₃` is the hexagon and `C_{k+3} = −C_k`. On all 126 proper 3-colourings of `M`: `δ(Z_k) = ±3`,
`δ(C_k) ∈ {0, ±6}`, the identities hold numerically, and whenever all `δ(C_k) = 0` the hexagon has `δ = ±6` (78
colourings through some `C_k`, 48 through the hexagon); all 126 have a tight cycle. The parity and divisibility
steps are correct. Lemma 20's proof was checked line by line (the floor construction, `q < D < p − q`, the
telescoping sums); it applies, so `χ_c(M) = χ_c(H₇) = 3`.

**Claim 4 (`K_{8/3}`): confirmed; citation secondary (finding 4).** `K_{8/3}` has 12 edges, `0, 3, 6, 1, 4, 7, 2, 5`
is an 8-cycle and the other edges join opposite vertices: it is the Wagner graph. `8/3` is the largest fraction below
3 with numerator at most 9. No homomorphism to `K_{8/3}` for `H₇`, `M` or the `ℚ(√31)` graph, by SAT, by the
exhaustive `8⁸` search `brute83.c` (count 0 for each) and by the Burnside enumeration below. Each `H₇ − v` maps to
`K_{8/3}` (16, 240, 480, 240, 16, 240, 480, 16 homomorphisms for `v = 1, …, 8`).

**Claim 5 (vertex-criticality): confirmed** for both files: every certificate is −1 exactly at `v`, proper on
`H − v`, with an acyclic tight digraph (Kahn); `N·c + pos` with `N = 8` is a homomorphism to `K_{24/9} = K_{8/3}` for
every `v`; SAT confirms directly that each `H − v` maps to `K_{8/3}`.

**Claim 6 (triangle-free graphs on at most 8 vertices): confirmed by a different scheme.** `tf_small.c` runs over all
`2^28 = 268 435 456` labelled graphs on 8 vertices and tests every triangle-free one, not only the maximal ones:
4 682 270 triangle-free, 15 247 maximal, 15 120 maximal and not bipartite; 0 without a homomorphism to `K_{8/3}`; 0
not 3-colourable; as a sanity check, 22 680 have no homomorphism to `K_{5/2}`, among them the 2 520 labelled Wagner
graphs. For `n = 3, …, 7` the triangle-free counts are 7, 41, 388, 5 789, 133 501. The reduction in
`small_triangle_free.py` (pad with isolated vertices, extend to a maximal triangle-free graph, restrict the
homomorphism) is valid; the repository's program gives the same output.

**Claim 7 (nine vertices: 1 897 classes, exactly three with `χ_c = 3`): confirmed with no isomorphism testing.**
`tf_burnside.c` enumerates, for one permutation of each of the 30 cycle types of `S₉`, the invariant triangle-free
graphs edge orbit by edge orbit. `Σ (class size × Fix) = 688 383 360 = 1 897 × 9!`; the same sum for "no
homomorphism to `K_{8/3}`" is `1 088 640 = 3 × 9!`; for `n = 1, …, 8` the counts are 1, 2, 3, 7, 14, 38, 107, 410,
none without a homomorphism. Identity term: 246 348 115 labelled triangle-free graphs; those without a homomorphism
number 30 240 with 12 edges, 90 720 with 13 and 45 360 with 14, that is `9!/12`, `9!/4` and `9!/8`: exactly the
labelled copies of `M`, `H₇` and `M + 2`. None fails to be 3-colourable. `nine_vertices.py` rerun: same output.

**Claim 8 (field facts): confirmed** (`r6_fields.py`). 21 and 93 are not squares, so `√3 ∉ ℚ(√7)` and
`√3 ∉ ℚ(√31)`, and a unit triangle forces `det² = 3/4`, so neither plane has one. Conditions (a) and (b) of the paper:
(a) fails and (b) holds for `d = 7` and 31; Theorem B and Corollary 1 give `χ = 3` and `χ_c = 3`. `√7` and `√31` lie
in `ℤ₃` (Hensel lift to `3^40`). By Corollary 6 (local fields), a finite extension `K` of `ℚ₃` has `χ_c(K²) < 4` only
if its residue degree is 1, and then `χ_c(K²) = 3`; `ℚ(√7)²` embeds in `K²` keeping unit and non-unit pairs, so `H₇`
is a witness there.

**Claim 9 (SAT records): confirmed.** `q7.cnf` (sha256 `258bcc7f…`): kissat 4.0.4 `s UNSATISFIABLE`; drat-trim
`s VERIFIED` (126 of 126 clauses and 49 of 158 lemmas in the core, as recorded); `cake_lpr` `s VERIFIED UNSAT`.
`q31.cnf` (sha256 `7c893b27…`): UNSAT, VERIFIED (115 of 115 clauses, 46 of 144 lemmas), VERIFIED UNSAT. The formula
of the graph as found (before translation and relabelling) also passes all three. `r3_cnf.py`: each clause set equals
exactly the expected one, with no duplicates; a model gives a proper colouring with a non-tight arc on every listed
cycle, and conversely. Without the cycle clauses the formula is satisfiable, and dropping any single listed cycle
makes it satisfiable.

**Claim 10 (reproduction and generalised programs): confirmed for `ℚ(√7)`; unclear for `ℚ(√31)` (finding 1).**
`grow.py q7_seed.json A G7 2000 200 3 1`: 2 growth rounds of 200 points, all by the rule "one colour left", 607
points, UNSAT with 7 cycles; `minimise.py`: 9 points, 13 edges, 6 cycles; `critical.py`: none non-critical. The found
points are the stored points translated by `(−1, 0)` under the map `{0:0, 1:1, 2:7, 3:2, 4:5, 5:6, 6:8, 7:4, 8:3}`;
edges, colouring, cycles and certificates correspond. Old (5cc8ae84) against new, with the clause-recording wrapper
`rec_run.py`: `grow.py` for `(7, 2)` on `q11.json` up to 401 points gives byte-identical logs (30 645 clauses, 35
solver calls); `minimise.py` identical clause streams, only the first and last calls changed from budgeted to
unbudgeted, identical outputs; `critical.py` identical logs.

**Claim 11 (other wording, novelty, seeds).** No novelty claims in these passages. The `ℚ(√31)` result is robust:
over 40 seed classes, 19 ended at nine vertices (14 edges 15 times, 13 edges 3 times, 12 edges once). Over `ℚ(√7)`,
25 random seeds gave nine-vertex witnesses with 14 edges (11 times) and 13 edges (7 times), never `M` itself.

## How the findings were applied

1. `q31_seed.json` is stored (0, the 156 unit vectors with denominator 80, and the two points that close the 5-cycle
   `0, (−1, 0), (−8/5, −4/5), (−8/5 − √31/16, 11/80), (−17/20 − √31/40, −17/40 + √31/20)`); `grow.py q31_seed.json A
   G31 2000 200 3 1` gives 367 vertices in 2 rounds, and `minimise.py` gives exactly the stored nine points, with the
   same coordinates, edges, colouring and cycles (a new test, `test_finite_witness_q31_found_by_growth`).
2. The remark is replaced: the growth also finds witnesses over `ℚ(√15)` and `ℚ(√39)`, of sizes that depend on the
   seed; the least size there is open.
3. All three graphs are realised over `ℚ(√7)` and over `ℚ(√31)`: `hexagon_search.py` counts the realisations of each,
   and `M` itself over `ℚ(√7)` is stored (`witness_q7m.json.gz`), checked by `check_small.py` and a certified SAT
   refutation. The paper, the note and the README say so.
4. Vince and Bondy–Hell are cited next to Zhu's survey.
5. Question 3 now says: at least nine in every plane without unit triangles, nine over `ℚ(√7)` and `ℚ(√31)`, open
   elsewhere (for example `ℚ(√15)`).
6. The proof is written with `C_k`, `C_{k+3}` reversed, only `C₀, C₁, C₂` needed, and the remark on `C₀` in `H₇`.
7. "Translated by `(1, 0)`"; the denominator remark; `nine_vertices.py` and `hexagon_search.py` in the note's heading;
   the `q31.cnf` figures in `verification.txt`.
8. The README says that `minimise.py` makes its first and last solver calls without a conflict budget.
