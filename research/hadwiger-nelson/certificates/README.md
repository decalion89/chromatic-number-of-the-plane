# Certificates: machine-checked claims

Each certificate is a JSON file that states one claim about a finite unit-distance graph and stores
the exact coordinates of its vertices, so that every edge can be recomputed. A claim that a graph
has no proper k-colouring is checked by a SAT solver on a formula rebuilt from the coordinates; a
claim that a colouring exists is checked from a stored colouring. The DRAT proofs are not stored,
since a solver regenerates them; the logs of three drat-trim checks are kept here.

Every claim can be rechecked without this repository's code, by the recipe below. The repository's
own checker is `python3 -m hn.cli verify <file> [--drat-trim PATH]`, which calls
`hn.certify.verify_certificate`; it checks four of the eight certificates completely, and the table
says what it does on the other four.

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
| `moser_spindle_no3coloring.json` | The Moser spindle (7 vertices, 11 edges, over ℚ(√3, √11)) has no proper 3-colouring, so χ(ℝ²) ≥ 4. | SAT solver UNSAT on the 21-variable, 40-clause formula; the research log records that drat-trim verified the DRAT proof, but no log is stored. | Handled: Glucose UNSAT, then drat-trim when it is installed. |
| `moser_spindle_4coloring.json` | The Moser spindle has a proper 4-colouring. | The stored colouring, checked on all 11 edges. | Handled: checks the colouring. |
| `degrey_1581_no4coloring.json` | De Grey's graph, rebuilt from the 39-point set S of arXiv:1804.02385 (1581 vertices, 7877 edges, over ℚ(√3, √5, √7, √11)), has no proper 4-colouring, so χ(ℝ²) ≥ 5. | kissat 4.0.4 UNSAT on the formula in which the triangle on vertices 0, 5 and 10 is pinned to colours 0, 1 and 2; DRAT proof checked by drat-trim (`degrey_1581_drat_trim_verification.txt`). | Not practical: it rebuilds the formula without the pinned triangle and gives it to Glucose; that unpinned formula has not been solved to completion in this project (research log). |
| `genuine_pair_19_no3coloring.json` | A vertex-critical graph with 19 vertices and 33 edges over ℚ(√3, √11)(√v), v = (−66 + 30√33)/256, has no proper 3-colouring; it is the two-orbit block on a pair of targets that is forced only jointly, neither target being forced on its own. | SAT solver UNSAT on the 57-variable, 118-clause formula; the research log records that drat-trim verified the DRAT proof, but no log is stored. | Handled: Glucose UNSAT, then drat-trim when it is installed. |
| `two_orbit_409_no3coloring.json` | A graph with 409 vertices and 1062 edges over ℚ(√3, √11), the union of six copies of a 91-vertex graph under two orbits of the 120° rotation about a pivot, has no proper 3-colouring. | DRAT proof of the 3 595-clause formula checked by drat-trim (`two_orbit_409_drat_trim_verification.txt`); as that log notes, the pair it blocks is already forced by the classical rhombus argument, so the block is not needed for this conclusion. | Handled: Glucose UNSAT, then drat-trim when it is installed. |
| `pressure3_witness_47.json` | In the 47-vertex subgraph of de Grey's `Sa` formed by the pivot (vertex 0), its 30 unit neighbours and 16 further points, every proper 4-colouring uses at least three colours on the 30 neighbours (pressure 3 at the pivot). | SAT solver UNSAT, measured with `hn/forced.py` and tested in `tests/test_forced.py`; no DRAT proof and no colouring are stored. | Rejected: a file with no colouring is read as a claim that the graph has no proper k-colouring, and this graph is 4-colourable. |
| `three_hexagon_pressure3.json` | In a 127-vertex, 528-edge graph over ℚ(√3, √11), formed by the pivot (vertex 0), three hexagons on its unit circle at angles 0, θ/2 and θ with cos θ = 5/6, and every sum u + v of points u, v in different hexagons, every proper 4-colouring uses at least three colours on the pivot's 18 unit neighbours. | SAT solver UNSAT, measured with `hn/forced.py` and tested in `tests/test_mixed.py`, which rebuilds the same 127 points; no DRAT proof is stored. The stored proper 4-colouring uses exactly three colours on the 18 neighbours, so the pressure is exactly 3. | Accepted, but it checks only that the stored colouring is proper, not the pressure claim. |
| `chain23_no3coloring.json` | A chain of three unit rhombi (10 vertices, 16 edges over ℚ(√2, √3), the graph of `data/chain23.json`) has no proper 3-colouring, so χ(ℚ(√2, √3)²) ≥ 4. | Three pysat solvers (CaDiCaL, Glucose, MiniSat) UNSAT; kissat's DRAT proof of the 30-variable, 58-clause formula, with no colour pinned, checked by drat-trim (`chain23_drat_trim_verification.txt`). | Fails with `TypeError`: it expects coefficients as strings, and this file stores `[numerator, denominator]` pairs. |

The claim in `pressure3_witness_47.json` calls the subgraph minimal. That holds with the 30 unit
neighbours of the pivot kept: deleting any one of the 16 further points lowers the pressure to 2.
It does not hold for deletions on the circle: a SAT check shows that 16 of the 30 neighbours can be
deleted, leaving 31 vertices, with the pressure still 3.

## Verification logs

| file | what it records |
|---|---|
| `chain23_drat_trim_verification.txt` | kissat (with `--no-binary`) and drat-trim on the 58-clause formula of `chain23_no3coloring.json`: `s VERIFIED`, with 18 of 37 lemmas in the core. |
| `degrey_1581_drat_trim_verification.txt` | drat-trim on the 33 101-clause formula with the pinned triangle and a binary DRAT proof of 1 319 558 301 bytes: `s VERIFIED`, with 2 016 499 of 13 140 458 lemmas in the core, in 522 s. |
| `degrey_1581_cnf_head.txt` | The first 14 and the last 14 lines of that formula: the header, the first 13 vertex clauses, the last two edge clauses and the 12 unit clauses that pin the triangle. |
| `two_orbit_409_drat_trim_verification.txt` | The output of `hn.cli verify` on `two_orbit_409_no3coloring.json` (drat-trim verified the proof of the 3 595-clause formula), with the construction and the pigeonhole argument for why six copies suffice; the script it names is now `scripts/experiments/validate_k3.py`. |

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
   unit neighbour v of the pivot, vertex 0. The claim is that this formula is unsatisfiable; by the
   same symmetry, a 4-colouring that uses at most two colours on those neighbours can be permuted
   to use only colours 0 and 1.
6. Run a DRAT-producing solver such as kissat on the formula, and check its proof with drat-trim.

`scripts/worker_setup.sh` installs kissat and drat-trim.
