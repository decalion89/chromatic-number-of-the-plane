# Tests

The tests check the project's constructions and results in exact arithmetic, recomputing each graph
from its stored coordinates. They cover field arithmetic, graph construction, colourings, forced
pairs, the local colourings, and the two theorems χ(ℚ(√3, √11)²) = 4, first proved by K. G. Fischer
(1994), and χ(ℚ(√2, √3)²) = 4.

```sh
cd research/hadwiger-nelson
python3 -m pytest -q                       # the whole suite, including the slow tests
python3 -m pytest -q -m "not slow"         # without the tests marked slow
python3 -m pytest -q tests/test_q23.py     # one file
```

The tests marked `slow` include a 4-colourability solve of de Grey's 1581-vertex graph
(`test_degrey.py`, with a 30-minute solver limit).

GitHub Actions (`.github/workflows/tests.yml` at the repository root) runs the 28 files marked CI
on every push to `main`, on every pull request, and on manual dispatch: 346 tests, which took 80
seconds in a local run. The DRAT test in `test_certify.py` is skipped there, because the
workflow does not install drat-trim. The other eight files (269 tests) are run locally; without their
two tests marked `slow` they took 32 minutes in a local run, 20 of them in one test of
`test_two_tunings.py`.

| file | CI | what it checks |
|---|:---:|---|
| `test_certify.py` | ✓ | The certificate checker in `hn/certify.py` accepts a proper colouring and rejects a monochromatic unit pair, a colour outside the palette and a claimed lower bound on a colourable graph; its DRAT test runs only when drat-trim is on the `PATH`. |
| `test_cyclotomic.py` | ✓ | Unit steps in ℤ[ζₙ]: they are the roots of unity, ℚ(ζ₅) alone has no unit triangle while ℚ(ζ₁₅) has one, a Gauss sum gives the spindle's radical, and the Moser spindle lies in a cyclotomic field. |
| `test_degrey.py` |  | De Grey's graph rebuilt from his 39-point set S has 39, 397 and 1581 vertices at the three stages and 7877 edges, and (marked `slow`) no proper 4-colouring. |
| `test_denominator_five.py` |  | Exoo–Ismailescu's graph rebuilt from their 23 points, the arithmetic of their rotation λ, the fact that a denominator of 5 blocks coset colourings while the integral 803-vertex graph admits them, and a 61-point graph on the Eisenstein lattice with edges at squared lengths 3, 4 and 7 that is 6-chromatic. |
| `test_dense_family.py` |  | The dense 5-chromatic graphs `five_dense_2.json` and `five_dense_10.json`, recomputed from their coordinates: the tuning angle, the vertex and edge counts, their structure as two copies of one carrier, and the absence of a proper 4-colouring. |
| `test_density.py` | ✓ | Independence ratios of unit-distance graphs (2/7 for the Moser spindle), and the consequence of Croft's 1-avoiding set of density 0.2293 that an independence-ratio argument cannot give χ ≥ 6 and a fractional one cannot give χ ≥ 5. |
| `test_disjunctive_spindle.py` | ✓ | The three- and six-copy disjunctive spindles, which need only c(v) = c(q₁) or c(v) = c(q₂) instead of a forced pair, checked in exact arithmetic and by an exhaustive sweep over the choices of the copies. |
| `test_fast_agrees.py` | ✓ | The vectorised integer arithmetic of `hn/fast.py` agrees with exact rational arithmetic on points, unit vectors, walks and edge sets. |
| `test_field.py` | ✓ | Exact arithmetic in multiquadratic fields: admissible generators, products of radicals, inverses through the Galois conjugates, exact equality and hashing. |
| `test_five_247.py` |  | The 5-chromatic unit-distance graphs over ℚ(√3, √11, √247) and ℚ(√3, √11, √23) in `data/`: edge counts recomputed exactly, chromatic number 5, C₆-invariance of the symmetric graphs, and the radicals that the spindle angles need. |
| `test_forced.py` | ✓ | Forced colour relations on graphs with known answers: pairs forced alike or different, cores, and the pressure at a vertex p (the least number of colours a k-colouring uses on the unit neighbours of p), including `certificates/pressure3_witness_47.json`. |
| `test_geometry.py` | ✓ | Exact rotations: the 60° rotation, the Moser spindle's angle arccos(5/6), rotations about a pivot, and de Grey's rotations, which need ℚ(√3, √5, √7, √11). |
| `test_graph_coloring.py` | ✓ | Graph construction, reductions and colouring against known values: the Moser spindle is 4-chromatic, the triangular lattice is 3-chromatic, and the forced pair of a unit rhombus is found and spindled. |
| `test_homcol.py` |  | Homomorphism (coset) colourings of unit-distance graphs and the blocking screen built on them, on the project's graphs and fields (197 tests). |
| `test_kappa.py` | ✓ | κ = ωρ₇ = (−11 + 5√−3)/14 is congruent to 1 mod 5, so no homomorphism from the module of `five_rho7.json` to ℤ/5 distinguishes a unit vector u from κu; with a sample of exact Stiemke certificates. |
| `test_mixed.py` | ✓ | Graphs with several distances: denesting of square roots, conflict rotations, escape counts, the two-orbit block, the 19-vertex jointly forced construction and the three-hexagon gadget with pressure 3 at four colours. |
| `test_moser_field.py` | ✓ | The unit-distance graph on the field ℚ(√−3, √−11), viewed in ℂ, is 4-colourable by a 2-adic colouring (`hn/adelic.py`); it contains the Moser spindle, so its chromatic number is 4. |
| `test_mu.py` | ✓ | μₖ(G, p), the least number of colours that a proper k-colouring of G uses on the unit neighbours of p, on cases with known answers: the hexagon, the Moser spindle and de Grey's `Sa` (μ = 2 at five colours). |
| `test_multispindle.py` | ✓ | The spindle lemma with several rotated copies: the rotation orders realisable over a multiquadratic field (the divisors of 24), the radicals they need, and which sets of targets they can block. |
| `test_neighbourhoods_are_bipartite.py` | ✓ | The unit neighbours of any vertex induce a bipartite graph (paths and hexagons), since two points on a unit circle are at distance 1 exactly when they subtend 60° at its centre, so no vertex is pinned by its neighbourhood at four colours. |
| `test_order_six_tuning.py` |  | Tuning the angle of the composite isometry to 60° gives it order 6, so the union of six copies of the carrier is C₆-invariant and the orbit of the pivot is forced monochromatic. |
| `test_q23.py` | ✓ | χ(ℚ(√2, √3)²) = 4: the plane has a proper 4-colouring by 2-adic residues and contains a 10-vertex chain of three unit rhombi with no proper 3-colouring (details below). |
| `test_q311.py` | ✓ | χ(ℚ(√3, √11)²) = 4, a theorem of K. G. Fischer (1994), by a short proof: the plane has a proper 4-colouring by 2-adic residues and contains the Moser spindle, which has no proper 3-colouring (details below). |
| `test_q35.py` | ✓ | The plane over ℚ(√3, √5) has no proper 3-colouring: 1/3 is a sum of unit vectors, and a chain of six unit rhombi (`data/chain35.json`, 19 vertices, 31 exact edges) joins the origin to a unit vector. 11 splits completely, which gives χ ≤ 5. A chain of four rhombi (13 vertices) built from four unit vectors summing to 1/√3 has no proper 3-colouring either. |
| `test_q3q.py` | ✓ | The planes over ℚ(√3, √q): for every q ≡ 2 (mod 3) up to 113, a generator α of norm 3^h in ℚ(√−q) makes 1/3 a sum of unit vectors, so χ ≥ 4; for q = 17, a chain of 90 rhombi (271 vertices) has no proper 3-colouring. |
| `test_split_places.py` | ✓ | Proposition C of `notes/local_colourings.md` §12: the hyperbola graphs H_q = Cay(𝔽_q², {(t, 1/t)}) have χ = 2, 3, 4, 3, 4, 4, 3, 4, 5, 4 for q = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16 (an explicit 5-colouring for q = 13); over ℚ(√3, √5), reduction at the prime above 2 is a proper 4-colouring of `chain35.json`, whose edge vectors are units there, one of its edge vectors is not a unit above 3, and reduction above 3 followed by a 3-colouring of H₉ 3-colours a graph whose edge vectors are units above 3. |
| `test_quadext.py` | ✓ | Arithmetic in ℚ(ζ₁₅, √−7, √−11), built as a tower of quadratic extensions of ℚ(ζ₁₅), which contains de Grey's rotations and ζ₁₅. |
| `test_reduce11.py` | ✓ | The unit-distance graph on the field ℚ(√−3, √−11, √−247), viewed in ℂ, is 5-colourable: reduction at a place above 11 maps it to the Cayley graph of 𝔽₁₂₁ on its 12 elements of norm one, which is 5- but not 4-colourable, and this colours the 803-vertex graph `five_247_c.json` properly at both places above 11. |
| `test_ring_geometry.py` | ✓ | Six facts behind the ring constructions, recomputed: the 90° rotation about the centre of a unit square moves one diagonal onto the other at distance 1, the map σ satisfies \|σu\|² = \|u\|²/3 and \|u − σu\| = \|u\|, the two-ring configuration is 3-chromatic, within each component of a neighbourhood in `five_247_c.json` two neighbours lie on the same side of the bipartition exactly when the angle between them is an even multiple of 60°, 2 is not of the form a² + ab + b², and that graph has exactly two vertices of degree 4. |
| `test_slack.py` | ✓ | With slack s = k − 3: unit-distance graphs have clique number 3, a pair forced alike under k colours needs at least k + 1 vertices, a unit rhombus attains this at k = 3, and at k = 4 neither the rhombus nor the Moser spindle forces a pair. |
| `test_spindle.py` | ✓ | The spindle constructions (two-copy, triple and local), separation tests and cores, on cases whose answers are known. |
| `test_tight_four.py` | ✓ | The 24-point graph of `data/tight_four.json` is 4-chromatic and has circular chromatic number exactly 4. |
| `test_transfer.py` | ✓ | The transfer-matrix construction, which glues copies of a graph along congruent sets of points, on the Moser spindle. |
| `test_transversal.py` | ✓ | Blocking by rotated copies of targets at one distance: rotation orders over a multiquadratic field divide 24, the triangle is the only odd cycle available, the capacity bounds that follow, and the counting certificate. |
| `test_tuned_family.py` |  | The tuned family of 5-chromatic graphs in `data/five_tuned_*.json`, where a rotation with cos φ = d²/(2D²) − 1 places the composite forced pair at a chosen distance d and so determines the field. |
| `test_two_tunings.py` |  | Tuning the forced distance to 1/√3 and then the angle to 60°: the pivot's orbit lies on a circle of radius 1/√3, six of its fifteen pairs are at distance 1, and the union of the six copies has no proper 4-colouring. |

## The two theorems

`test_q311.py`, for χ(ℚ(√3, √11)²) = 4 (K. G. Fischer, Congr. Numer. 104 (1994) 73–79), checks a
short proof by 2-adic residues:
- the local facts: ℚ(√3, √11) has two places over 2, each with completion ℚ₂(√3) and inert in
  ℚ(√3, √11, i), whose completions have residue field 𝔽₄; sympy's `prime_decomp` confirms that
  2 = 𝔭₁²𝔭₂² with residue degree 1;
- 810 unit vectors generated by rotations, and 300 random unit vectors
  ((a² − b²)/(a² + b²), 2ab/(a² + b²)), a form that every unit vector takes by Hilbert 90, are
  2-adic units with nonzero residue in 𝔽₄ at both places;
- the colouring is proper on 300 random unit steps, on the Moser spindle (which has no proper
  3-colouring, the lower bound) and on Exoo–Ismailescu's 214-point graph H (`data/ei_H214.json`);
- on random samples, it gives pairs at distance 8/9ⁿ (n ≤ 3) the same colour and pairs at distance
  √(11/3) different colours.

`test_q23.py`, for χ(ℚ(√2, √3)²) = 4, a case that Fischer's hypotheses exclude, checks:
- the local facts: ℚ(√2, √3) has one place over 2, totally ramified with residue field 𝔽₂ and not
  split in ℚ(√2, √3, i), whose completion has residue field 𝔽₄; sympy's `prime_decomp` confirms
  that 2 = 𝔭⁴;
- 312 listed unit vectors, and 300 random unit vectors of the Hilbert 90 form above, have nonzero
  residue;
- the colouring is proper on the 2 089 sums of at most three 24th roots of unity and on a graph of
  sums of two listed unit vectors;
- the set M₂ = M₁ + M₁ of Voronov, Neopryatnaya and Dergachev (arXiv:2106.11824, Series 2), with
  2 593 points and 11 448 edges, has no proper 3-colouring and is 4-coloured by the residue
  colouring;
- the chain of three unit rhombi in `data/chain23.json` (10 vertices, 16 edges) has no proper
  3-colouring, the lower bound.
