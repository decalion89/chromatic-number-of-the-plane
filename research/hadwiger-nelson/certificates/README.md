# Certificates

Each certificate states one claim and holds everything needed to recheck it from scratch:
- the exact vertex coordinates;
- the edge count, re-derived by exact arithmetic;
- a colouring, or a pointer to a DRAT proof that `drat-trim` has verified.

None of them asks you to trust this repository's code.

| file | claim | how it is certified |
|---|---|---|
| `moser_spindle_no3coloring.json` | χ(ℝ²) ≥ 4: the Moser spindle (7 vertices, 11 edges) has no proper 3-colouring | SAT UNSAT, DRAT verified |
| `moser_spindle_4coloring.json` | the Moser spindle is 4-colourable | an explicit colouring |
| `degrey_1581_no4coloring.json` | χ(ℝ²) ≥ 5: de Grey's 1581-vertex graph, rebuilt from its 39-point seed, has no proper 4-colouring | kissat UNSAT with one triangle's colours fixed (a standard symmetry step); DRAT verified, see `degrey_1581_drat_trim_verification.txt` and `degrey_1581_cnf_head.txt` |
| `genuine_pair_19_no3coloring.json` | a vertex-critical 19-vertex, 33-edge graph with no 3-colouring, whose forced pair is forced only jointly | SAT UNSAT, DRAT verified |
| `two_orbit_409_no3coloring.json` | a 409-vertex graph with no 3-colouring, by a two-orbit block | DRAT verified, see `two_orbit_409_drat_trim_verification.txt` |
| `pressure3_witness_47.json` | 47 vertices of de Grey's `Sa` keep pressure 3 at the pivot: no 4-colouring squeezes its 30-point circle into two colours | SAT, measured with `hn/forced.py` |
| `three_hexagon_pressure3.json` | pressure 3 at the pivot at four colours, on 127 vertices | SAT, measured with `hn/forced.py` |
| `chain23_no3coloring.json` | χ(ℚ(√2, √3)²) ≥ 4: a 10-vertex chain of three unit rhombi has no proper 3-colouring | three pysat solvers UNSAT; DRAT verified, see `chain23_drat_trim_verification.txt` |

To recheck a DRAT certificate:
1. Rebuild the CNF from the coordinates: 3 or 4 variables per vertex, one at-least-one clause per
   vertex, and one conflict clause per edge and colour.
2. Run a DRAT-producing solver such as kissat.
3. Check its proof with `drat-trim`.

`scripts/worker_setup.sh` installs both tools.
