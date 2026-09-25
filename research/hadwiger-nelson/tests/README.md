# Tests

The suite recomputes every stored result from its exact coordinates: field arithmetic, graph
construction, colourings, forced pairs, the local colourings and the two theorems.

```sh
cd research/hadwiger-nelson
python3 -m pytest -q                       # everything: several hours (it re-solves de Grey's graph)
python3 -m pytest -q tests/test_q23.py     # one file
```

GitHub Actions runs the files marked **CI** on every push and pull request: about 330 tests in
under two minutes. The others take minutes to hours and run locally.

| file | CI | what it checks |
|---|:---:|---|
| `test_certify.py` | ✓ | Certificates must catch wrong answers, not just bless right ones. |
| `test_cyclotomic.py` | ✓ | Z[zeta_n]: the unit steps, the triangles, and what the field buys. |
| `test_degrey.py` |  | The de Grey reconstruction, checked at every stage. |
| `test_denominator_five.py` |  | The denominator-five facts, recomputed. |
| `test_dense_family.py` |  | The dense 5-chromatic graphs, recomputed from their stored coordinates. |
| `test_density.py` | ✓ | The independence-ratio certificate, and the Moser spindle's 2/7. |
| `test_disjunctive_spindle.py` | ✓ | The three- and six-fold disjunctive spindles, checked in exact arithmetic. |
| `test_fast_agrees.py` | ✓ | The vectorised integer path must agree with exact Fraction arithmetic. |
| `test_field.py` | ✓ | The field is the foundation: if arithmetic here is wrong, every edge is. |
| `test_five_247.py` |  | The 5-chromatic unit-distance graphs in Q(sqrt3, sqrt11, sqrt247). |
| `test_forced.py` | ✓ | Forced colour relations, and the ladder chi >= 6 reduces to. |
| `test_geometry.py` | ✓ | Rotations must be exact rotations, and the lattice must be the lattice. |
| `test_graph_coloring.py` | ✓ | Graph construction, reductions, and colouring, checked against known values. |
| `test_homcol.py` |  | Homomorphism colourings: the structural screen, and what it says. |
| `test_kappa.py` | ✓ | kappa = omega * rho7 = (-11 + 5 sqrt-3)/14 is congruent to 1 mod 5. |
| `test_mixed.py` | ✓ | Mixed-distance machinery: denesting, conflict rotations, escape counts. |
| `test_moser_field.py` | ✓ | The unit-distance graph on the Moser field Q(sqrt-3, sqrt-11) is 4-colourable (hn/adelic.py). |
| `test_mu.py` | ✓ | mu, the number that is the problem, checked on objects with known answers. |
| `test_multispindle.py` | ✓ | The general spindle lemma, and the arithmetic that bounds it. |
| `test_neighbourhoods_are_bipartite.py` | ✓ | Neighbourhoods in the plane are bipartite, so rigidity is never local. |
| `test_order_six_tuning.py` |  | Tuning the ANGLE instead of the distance closes the chain into a group. |
| `test_q23.py` | ✓ | The plane over Q(sqrt2, sqrt3) is 4-colourable and holds a 4-chromatic graph, so its chromatic number is 4. |
| `test_q311.py` | ✓ | The plane over Q(sqrt3, sqrt11) is 4-colourable, so its chromatic number is exactly 4 (hn/adelic.py). |
| `test_quadext.py` | ✓ | Q(zeta_15, sqrt(-7), sqrt(-11)): de Grey's rotations and zeta_15 together. |
| `test_reduce11.py` | ✓ | The unit-distance graph on Q(sqrt-3, sqrt-11, sqrt-247) is 5-colourable by reduction at 11. |
| `test_ring_geometry.py` | ✓ | The ring theorems, recomputed rather than quoted. |
| `test_slack.py` | ✓ | The slack framing: why the plane's difficulty jumps where it does. |
| `test_spindle.py` | ✓ | The spindling arguments, checked on cases whose answers are known. |
| `test_tight_four.py` | ✓ | The 24-point graph with circular chromatic number exactly four. |
| `test_transfer.py` | ✓ | Gluing copies along an interface: the transfer-matrix construction. |
| `test_transversal.py` | ✓ | The same-distance ceiling, and the measurement that corrected it. |
| `test_tuned_family.py` |  | The tuned family: 5-chromatic graphs whose field is a parameter. |
| `test_two_tunings.py` |  | Both tunings at once: the group and the chromatic number, together. |
