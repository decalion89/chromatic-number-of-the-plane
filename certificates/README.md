# Certificates: machine-checked claims

Each certificate is a JSON file that states one claim about a finite unit-distance graph and stores
the exact coordinates of its vertices, so that every edge can be recomputed. A claim that a graph
has no proper k-colouring is checked by a SAT solver on a formula rebuilt from the coordinates; a
claim that a colouring exists is checked from a stored colouring. The DRAT proofs are not stored,
since a solver regenerates them; the verification logs listed below are kept here.

Every claim can be rechecked without this repository's code, by the recipe below. The repository's
own checker is `python3 -m hn.cli verify <file> [--drat-trim PATH]`, which calls
`hn.certify.verify_certificate`; it checks six of the nine certificates completely, and the table
says what it does on the other three. It prints `VERIFIED` (exit status 0) when a colouring has been
checked on the exact edges or a DRAT proof by drat-trim, `REJECTED` (1) when the claim fails, and
the solver's `UNSAT` with the note that the proof is unverified (2) when drat-trim is not installed.

## Format

- `claim`: the statement; `k`: the number of colours; `n`, `m`: the numbers of vertices and of
  unit-distance pairs; `notes`: how the graph was built and checked.
- `field_generators` `[a, b, …]`: the coordinates lie in ℚ(√a, √b, …); `basis` lists the basis
  1, √a, √b, √ab, … on which they are written.
- `vertices`: a list of `{"x": …, "y": …}`; each coordinate is the list of its rational
  coefficients on `basis`, written as strings such as `"-109/48"`. In `chain23_no3coloring.json`
  the coefficients are `[numerator, denominator]` pairs instead, as in `data/`.
- In `genuine_pair_19_no3coloring.json` (`"kind": "real_quadratic_extension"`) the field is K(√v)
  with K = ℚ(√3, √11) and v = (−66 + 30√33)/256, given as `radicand` on `base_basis`; each
  coordinate is `{"a": …, "b": …}`, meaning a + b√v with a, b in K.
- `coloring`: a proper k-colouring (colours 0, …, k − 1 in vertex order), or `null`.

## The certificates

| file | claim | how it was checked | `verify_certificate` |
|---|---|---|---|
| `moser_spindle_no3coloring.json` | The Moser spindle (7 vertices, 11 edges, over ℚ(√3, √11)) has no proper 3-colouring, so χ(ℝ²) ≥ 4. | kissat 4.0.4 UNSAT on the 21-variable, 40-clause formula; its DRAT proof, checked by drat-trim (`moser_spindle_drat_trim_verification.txt`). | Handled: Glucose UNSAT, then drat-trim when it is installed. |
| `moser_spindle_4coloring.json` | The Moser spindle has a proper 4-colouring. | The stored colouring, checked on all 11 edges. | Handled: checks the colouring. |
| `degrey_1581_no4coloring.json` | De Grey's graph, rebuilt from the 39-point set S of arXiv:1804.02385 (1581 vertices, 7877 edges, over ℚ(√3, √5, √7, √11)), has no proper 4-colouring, so χ(ℝ²) ≥ 5. | kissat 4.0.4 UNSAT on the formula in which the triangle on vertices 0, 5 and 10 is pinned to colours 0, 1 and 2; DRAT proof checked by drat-trim (`degrey_1581_drat_trim_verification.txt`). | Not practical: it rebuilds the formula without the pinned triangle and gives it to Glucose, which had not finished after 25 minutes in a test (about 5 of them spent deriving the edges); a kissat run on the unpinned formula was stopped after 110 minutes without a verdict (research log, Status). |
| `genuine_pair_19_no3coloring.json` | A vertex-critical graph with 19 vertices and 33 edges over ℚ(√3, √11)(√v), v = (−66 + 30√33)/256, has no proper 3-colouring; it is the two-orbit block on a pair of targets that is forced only jointly, neither target being forced on its own. | SAT solver UNSAT on the 57-variable, 118-clause formula; kissat's DRAT proof, checked by drat-trim (`genuine_pair_19_drat_trim_verification.txt`). | Handled: Glucose UNSAT, then drat-trim when it is installed. |
| `two_orbit_409_no3coloring.json` | A graph with 409 vertices and 1062 edges over ℚ(√3, √11), the union of six copies of a 91-vertex graph under two orbits of the 120° rotation about a pivot, has no proper 3-colouring. | DRAT proof of the 3 595-clause formula checked by drat-trim (`two_orbit_409_drat_trim_verification.txt`); as that log notes, the pair it blocks is already forced by the classical rhombus argument, so the block is not needed for this conclusion. | Handled: Glucose UNSAT, then drat-trim when it is installed. |
| `pressure3_witness_47.json` | In the 47-vertex subgraph of de Grey's `Sa` formed by the pivot (vertex 0), its 30 unit neighbours and 16 further points, every proper 4-colouring uses at least three colours on the 30 neighbours (pressure 3 at the pivot). | SAT solver UNSAT, measured with `hn/forced.py` and tested in `tests/test_forced.py`; kissat's DRAT proof for the formula of step 5 below, checked by drat-trim (`pressure3_witness_47_drat_trim_verification.txt`). No colouring is stored. | Rejected: a file with no colouring is read as a claim that the graph has no proper k-colouring, and this graph is 4-colourable. |
| `three_hexagon_pressure3.json` | In a 127-vertex, 528-edge graph over ℚ(√3, √11), formed by the pivot (vertex 0), three hexagons on its unit circle at angles 0, θ/2 and θ with cos θ = 5/6, and every sum u + v of points u, v in different hexagons, every proper 4-colouring uses at least three colours on the pivot's 18 unit neighbours. | SAT solver UNSAT, measured with `hn/forced.py` and tested in `tests/test_mixed.py`, which rebuilds the same 127 points; kissat's DRAT proof for the formula of step 5 below, checked by drat-trim (`three_hexagon_pressure3_drat_trim_verification.txt`). The stored proper 4-colouring uses exactly three colours on the 18 neighbours, so the pressure is exactly 3. | Accepted, but it checks only that the stored colouring is proper, not the pressure claim. |
| `chain23_no3coloring.json` | A chain of three unit rhombi (10 vertices, 16 edges over ℚ(√2, √3), the graph of `data/chain23.json`) has no proper 3-colouring, so χ(ℚ(√2, √3)²) ≥ 4. | Three pysat solvers (CaDiCaL, Glucose, MiniSat) UNSAT; kissat's DRAT proof of the 30-variable, 58-clause formula, with no colour pinned, checked by drat-trim (`chain23_drat_trim_verification.txt`). | Handled: Glucose UNSAT, then drat-trim when it is installed. The file stores coefficients as `[numerator, denominator]` pairs, which the loader reads as well as decimal strings. |
| `five_247_c_no4coloring.json` | The 803-vertex graph `five_247_c` (`data/five_247_c.json`: vertex-critical, 4 065 edges, over ℚ(√3, √11, √247)) has no proper 4-colouring. It lies in the plane over ℚ(√−3, √−11, √−247) as well, so both planes need five colours. | kissat 4.0.4 UNSAT on the plain 3 212-variable, 17 063-clause formula, with no colour pinned; its DRAT proof, checked by drat-trim (`five_247_c_drat_trim_verification.txt`). CaDiCaL and Glucose also report UNSAT. | Handled: Glucose UNSAT (17 minutes on a loaded machine), then drat-trim when it is installed. |

The claim in `pressure3_witness_47.json` calls the subgraph minimal among those that keep the 30
unit neighbours of the pivot: deleting any one of the 16 further points lowers the pressure to 2.
It is not minimal for deletions on the circle: a SAT check shows that 16 of the 30 neighbours can
be deleted, leaving 31 vertices, with the pressure still 3.

## Verification logs

Most logs belong to a certificate above. Four record checks of files in `data/`: the two
multi-distance witnesses, the forced pair of Exoo and Ismailescu's graph H, and every graph there
said to have no proper 4-colouring (`data_no4_checks.txt`). Four record `α(G₁₃) ≤ 36` (`notes/g13.md`): its two cases, the cube tree of the
second, and the recheck by cake_lpr (`g13_*`). One records part A of
`scripts/g17_alpha.py`, a step towards α(G₁₇) ≤ 57 (`g17_part_a_checks.txt`).

| file | what it records |
|---|---|
| `chain23_drat_trim_verification.txt` | kissat (with `--no-binary`) and drat-trim on the 58-clause formula of `chain23_no3coloring.json`: `s VERIFIED`, with 18 of 37 lemmas in the core. |
| `data_no4_checks.txt` | `scripts/check_no4.py` on the 27 graphs of `data/` described as having no proper 4-colouring, from 803 to 39 313 vertices: each rebuilt from its coordinates, with one triangle pinned; kissat finds every formula unsatisfiable and drat-trim verifies all 27 DRAT proofs. Each line gives the SHA-256 of the formula. |
| `degrey_1581_drat_trim_verification.txt` | drat-trim on the 33 101-clause formula with the pinned triangle and a binary DRAT proof of 1 319 558 301 bytes: `s VERIFIED`, with 2 016 499 of 13 140 458 lemmas in the core, in 522 s. |
| `degrey_1581_cnf_head.txt` | The first 14 and the last 14 lines of that formula: the header, the first 13 vertex clauses, the last two edge clauses and the 12 unit clauses that pin the triangle. |
| `ei_H214_forced_pair_drat_trim_verification.txt` | kissat and drat-trim on the formula of `tests/test_denominator_five.py` for `data/ei_H214.json`, Exoo and Ismailescu's graph H with edges at distances 1 and 2 and its pair A, B at distance 5 required to differ, one triangle pinned: unsatisfiable, so A and B share a colour in every such 5-colouring; `s VERIFIED` with 436 214 of 570 197 lemmas in the core. |
| `five_247_c_drat_trim_verification.txt` | kissat and drat-trim on the 17 063-clause formula of `five_247_c_no4coloring.json`, rebuilt from the certificate: `s VERIFIED`, with 1 455 257 of 4 006 246 lemmas in the core, in 463 s; the file gives the SHA-256 of the formula and of the proof. |
| `g13_alpha_part_a_checks.txt` | `scripts/g13/certify.py`, case A of `α(G₁₃) ≤ 36`: the four formulas `E37_A_c` of `scripts/g13/enum_cert.py` (an independent dominating set of 37 points containing 0 and the whole circle N = c, c = 6, 7, 9, 11), with the SHA-256 of each, kissat UNSAT and drat-trim VERIFIED. |
| `g13_alpha_part_b_cubes.icnf` | The 4 822 leaves of the cube tree of case B (`scripts/g13/cuber2.py`, at most 14 decisions along `lex_order()`), 564 of them closed by unit propagation; `scripts/verify_g13.py` checks that they cover every assignment. |
| `g13_alpha_part_b_checks.txt.gz` | One line per leaf of case B: formula E37_B of `scripts/g13/enum_cert.py` plus the leaf, its SHA-256, kissat UNSAT within 120 s and drat-trim VERIFIED. |
| `g13_cake_lpr_checks.txt.gz` | `scripts/verify_g13.py` with cake_lpr, on another machine: every formula of both cases (the four of case A and the 4 822 leaves of case B) rebuilt, compared with the logs, refuted again by kissat and its proof checked by drat-trim and by cake_lpr; and the cover formula's proof checked by cake_lpr. |
| `g17_part_a_checks.txt` | `scripts/g17_alpha.py`, part A: for each of the seven circles N = c of G₁₇ with no unit distance inside, the formula saying that 39 vertices of the region adjacent to none of {0} + C_c are independent; kissat finds all seven unsatisfiable and drat-trim verifies all seven DRAT proofs, so no independent set of 58 points contains a point together with its whole circle. Each line gives the SHA-256 of the formula, which `tests/test_g17.py` recomputes. |
| `genuine_pair_19_drat_trim_verification.txt` | kissat (with `--no-binary`) and drat-trim on the 118-clause formula of `genuine_pair_19_no3coloring.json`, rebuilt from the certificate: `s VERIFIED`, with 50 of 91 lemmas in the core. |
| `moser_spindle_drat_trim_verification.txt` | kissat (with `--no-binary`) and drat-trim on the 40-clause formula of `moser_spindle_no3coloring.json`, rebuilt from the certificate: `s VERIFIED`, with 12 of 26 lemmas in the core. |
| `pressure3_witness_47_drat_trim_verification.txt` | kissat (with `--no-binary`) and drat-trim on the 583-clause formula of step 5 for `pressure3_witness_47.json`: `s VERIFIED`, with 301 of 895 lemmas in the core. |
| `threepoint_indep_checks.txt` | `scripts/threepoint_verify_indep.py` on the eight three-point certificates of `data/threepoint/`: for each, the orbits and constraints it rebuilds, the dual blocks it proves positive definite exactly, and the exact rational bound on α, each within 4·10⁻⁶ of the solver's value and below the number quoted in `notes/local_colourings.md` §14. |
| `three_hexagon_pressure3_drat_trim_verification.txt` | The same for `three_hexagon_pressure3.json`, 2 275 clauses: `s VERIFIED`, with 429 of 1 880 lemmas in the core. |
| `two_orbit_409_drat_trim_verification.txt` | The one-line verdict of `hn.cli verify` on `two_orbit_409_no3coloring.json` (drat-trim verified the proof of the 3 595-clause formula; drat-trim's own transcript was not kept), with the construction and the pigeonhole argument for why six copies suffice; the script it names is now `scripts/experiments/validate_k3.py`. |
| `W_lattice_16_21_28_61_drat_trim_verification.txt` | kissat (with `--no-binary`) and drat-trim on `data/W_lattice_16_21_28_61.json` (72 points, 553 edges at distances 1, 4/√3, √7, √(28/3) and √(61/3), recomputed from the lattice coordinates): the 5-colouring formula of 2 837 clauses is unsatisfiable, `s VERIFIED` with 415 477 of 575 631 lemmas in the core. |
| `W_moser_orbit_9_33_drat_trim_verification.txt` | `scripts/orbit_witness_test.py` with kissat and drat-trim on `data/W_moser_orbit_9_33.json` (187 points, 508 edges at distance 1 and 495 at the orbit d² = (9 ∓ √33)/6): the 5-colouring formula of 5 202 clauses is unsatisfiable, `s VERIFIED` with 52 186 of 79 301 lemmas in the core. |

## Checking a certificate without this repository's code

1. Recompute the edges: every pair of vertices at squared distance exactly 1, by exact arithmetic
   in the stated field. Their number must equal `m`.
2. For a colouring, check that the two ends of every edge receive different colours.
3. For a claim that no proper k-colouring exists, build the formula with the variable
   x(v, c) = 1 + k·v + c for vertex v = 0, …, n − 1 (in file order) and colour c = 0, …, k − 1:
   one clause x(v, 0) ∨ … ∨ x(v, k − 1) per vertex, and one clause ¬x(u, c) ∨ ¬x(v, c) per edge
   {u, v} and colour c. It has n·k variables and n + k·m clauses.
4. For `degrey_1581_no4coloring.json`, add the 12 unit clauses that pin the triangle on vertices
   0, 5 and 10 to colours 0, 1 and 2; in DIMACS they are `1`, `22`, `43`, `-2`, `-3`, `-4`, `-21`,
   `-23`, `-24`, `-41`, `-42`, `-44`. The formula then has 6 324 variables and
   33 101 = 1581 + 4·7877 + 12 clauses. Written with the vertex clauses first, the edge clauses in
   lexicographic order of (u, v) with u < v, and the unit clauses last, it begins and ends with the
   lines in `degrey_1581_cnf_head.txt`. One step lies outside the DRAT proof: in any proper
   colouring the three vertices of a triangle receive distinct colours, and a permutation of the
   colours makes them 0, 1 and 2.
5. For the two pressure certificates (k = 4), add the unit clauses ¬x(v, 2) and ¬x(v, 3) for every
   unit neighbour v of the pivot, vertex 0: 188 variables and 583 clauses for
   `pressure3_witness_47.json`, 508 variables and 2 275 clauses for `three_hexagon_pressure3.json`.
   The claim is that this formula is unsatisfiable; by the same symmetry, a 4-colouring that uses at
   most two colours on those neighbours can be permuted to use only colours 0 and 1.
6. Run a DRAT-producing solver such as kissat on the formula, and check its proof with drat-trim.

`scripts/worker_setup.sh` installs kissat and drat-trim.
