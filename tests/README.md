# Tests

The tests check the project's constructions and results in exact arithmetic, recomputing each graph
from its stored coordinates. They cover field arithmetic, graph construction, colourings, forced
pairs, the local colourings, and the two theorems χ(ℚ(√3, √11)²) = 4, first proved by K. G. Fischer
(1994), and χ(ℚ(√2, √3)²) = 4.

```sh
python3 -m pytest -q                       # the whole suite, including the slow tests
python3 -m pytest -q -m "not slow"         # without the tests marked slow
python3 -m pytest -q tests/test_q23.py     # one file
```

The tests marked `slow` include a 4-colourability solve of de Grey's 1581-vertex graph
(`test_degrey.py`), which can take CaDiCaL hours. Every solver call in the tests has a four-hour
limit; `hn.coloring` enforces it for CaDiCaL, which ignores pysat's interrupt, by running the
solver in a separate process.

GitHub Actions (`.github/workflows/tests.yml`) runs the 41 files marked CI on every push to
`main`, on every pull request, and on manual dispatch: 525 tests. It skips the DRAT test in
`test_certify.py`, because the workflow does not install drat-trim, and leaves out the six slow tests
of `test_threepoint_indep.py`. The other eleven files (283 tests) are run locally; without their
15 tests marked `slow` they took 32 minutes in a local run, 20 of them in one test of
`test_two_tunings.py`. The eight slow tests of `test_threepoint_certificates.py` take about 15
minutes together.

| file | CI | what it checks |
|---|:---:|---|
| `test_biquadratic_bounds.py` | ✓ | Local upper bounds for the planes over ℚ(√a, √b) (`notes/local_colourings.md` §13): χ(G_q) = 3, 4, 5 for q = 3, 7, 11, χ(G₄) = 4 and χ(G₁₉) = 5 (a linear 5-colouring; no 4-colouring); a ramified place above 2 makes the plane over ℚ(√2, √5) bipartite; the bounds for ℚ(√3, √q), q < 60, agree with Theorem 5 (§11); the first non-split places of ℚ(√3, √29) and ℚ(√5, √7) lie above 23 and 19; a unit 5-cycle over ℚ(√5, √7) gives χ ≥ 3 there. |
| `test_certify.py` | ✓ | The certificate checker in `hn/certify.py` accepts a proper colouring and rejects a monochromatic unit pair, a colour outside the palette and a claimed lower bound on a colourable graph; its DRAT test runs only when drat-trim is on the `PATH`. |
| `test_cli.py` | ✓ | The command line `python3 -m hn.cli`: `verify` accepts both Moser spindle certificates, reports a solver's `UNSAT` without calling it verified when drat-trim is missing, and rejects a tampered colouring; `demo` writes a certificate that verifies, by default to `HN_OUT`. |
| `test_cyclotomic.py` | ✓ | Unit steps in ℤ[ζₙ]: they are the roots of unity, ℚ(ζ₅) alone has no unit triangle while ℚ(ζ₁₅) has one, a Gauss sum gives the spindle's radical, and the Moser spindle lies in a cyclotomic field. |
| `test_decompositions.py` | ✓ | The decompositions of 2 and 3 recomputed with PARI/GP (`scripts/decompositions.gp`): the places over 2 of ℚ(√3, √11) and ℚ(√2, √3) and their inertness in L(i); for ℚ(√3, √q), q < 75 prime, a place over 2 with residue field 𝔽₂ exactly when q = 2 or q ≡ 1, 3 (mod 8), and one over 3 with residue field 𝔽₃ exactly when q ≡ 1 (mod 3); nine fields with more square roots. Skipped when gp is not installed. |
| `test_degrey.py` |  | De Grey's graph rebuilt from his 39-point set S has 39, 397 and 1581 vertices at the three stages and 7877 edges, and (marked `slow`) no proper 4-colouring. |
| `test_denominator_five.py` |  | Exoo–Ismailescu's graph rebuilt from their 23 points, the arithmetic of their rotation λ, the fact that a denominator of 5 blocks coset colourings while the integral 803-vertex graph admits them, and a 61-point graph on the Eisenstein lattice with edges at squared lengths 3, 4 and 7 that is 6-chromatic. |
| `test_dense_family.py` |  | The dense 5-chromatic graphs `five_dense_2.json` and `five_dense_10.json`, recomputed from their coordinates: the tuning angle, the vertex and edge counts, their structure as two copies of one carrier, and the absence of a proper 4-colouring. |
| `test_density.py` | ✓ | Independence ratios of unit-distance graphs (2/7 for the Moser spindle), and the consequence of Croft's 1-avoiding set of density 0.2293 (H. T. Croft, *Incidence incidents*, Eureka 30 (1967) 22–26) that an independence-ratio argument cannot give χ ≥ 6 and a fractional one cannot give χ ≥ 5 (for the fractional case, R. Hochberg and P. O'Donnell, Geombinatorics 2(4) (1993) 83–84). |
| `test_disjunctive_spindle.py` | ✓ | The three- and six-copy disjunctive spindles, which need only c(v) = c(q₁) or c(v) = c(q₂) instead of a forced pair, checked in exact arithmetic and by an exhaustive sweep over the choices of the copies. |
| `test_fast_agrees.py` | ✓ | The vectorised integer arithmetic of `hn/fast.py` agrees with exact rational arithmetic on points, unit vectors, walks and edge sets. |
| `test_field.py` | ✓ | Exact arithmetic in multiquadratic fields: admissible generators, products of radicals, inverses through the Galois conjugates, exact equality and hashing. |
| `test_finite_planes.py` | ✓ | Bounds for Moorhouse's table of χ(𝔽_q²) (`notes/local_colourings.md` §12): interval colourings, with m consecutive parallel lines ax + by = r per colour when a² + b² − r² is a non-square for r < m, for q = 7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61 (optimal for 7, 13, 19); a linear 8-colouring of 𝔽₄₃²; no 4-colouring of 𝔽₂₃² or 𝔽₃₇²; linear colourings of 𝔽₁₇² need six colours (every admissible circulant has independence number at most 3). Hoffman's bound with every eigenvalue in interval arithmetic (§14): no 4-colouring of 𝔽₂₃² or 𝔽₃₁², χ(𝔽₅₉²) ≥ 6, χ ≥ 7 for q = 71, 97, 101, and χ(G_q) ≥ 6 for the anisotropic planes with q = 53, 59, 61, the cases of Proposition B below Weil's threshold. |
| `test_finite_planes_slow.py` |  | No 4-colouring of 𝔽₂₉², 𝔽₃₁², 𝔽₄₁² or 𝔽₄₃² (marked `slow`; 15 s to 8 minutes each). |
| `test_five_247.py` |  | The 5-chromatic unit-distance graphs over ℚ(√3, √11, √247) and ℚ(√3, √11, √23) in `data/`: edge counts recomputed exactly, chromatic number 5, C₆-invariance of the symmetric graphs, and the radicals that the spindle angles need. |
| `test_flat852.py` | ✓ | The 852-vertex spindle-free graph (`notes/flat852.md`), with `scripts/verify_flat852.py`: the arithmetic of ℚ(ζ₂₁) (Φ₂₁, z²¹ = 1, conjugation); the 126 directions are distinct unit vectors closed under negation and conjugation, and ω is a root of 7x¹² − 13x⁶ + 7; every edge of `core852` (which uses all 126) and of `g1023` is one of them, recomputed exactly; the 5-colourings are proper; the stored CNF files are exactly the graphs' formulas; both drat-trim logs say VERIFIED; a moved point and a bad colouring are caught. |
| `test_forced.py` | ✓ | Forced colour relations on graphs with known answers: pairs forced alike or different, cores, and the pressure at a vertex p (the least number of colours a k-colouring uses on the unit neighbours of p), including `certificates/pressure3_witness_47.json`. |
| `test_g13.py` | ✓ | `α(G₁₃) = 36` (`notes/g13.md`): the graph, its independent circles N = 1, 6, 7, 9, 11 and its 4 732 automorphisms, and the 15 known 36-point sets; the encoders of `scripts/g13/g13cnf.py` by brute force (the totalizer both ways on up to 7 literals, the counter, the lex-leader chains on random permutations and with the automorphisms of G₁₃, value precedence); the four formulas of case A have the SHA-256 of their log, the leaves of case B cover every assignment and each leaf formula has a drat-trim VERIFIED line; the resolution proof of the cover ends in the empty clause. The log of the second run has one line for each of the 4 826 formulas, with the same SHA-256, where cake_lpr checked its proof, and a cake_lpr VERIFIED line for the cover. |
| `test_g13_audit.py` | ✓ | The audit `scripts/g13/g13_audit.py` of the five formulas of `α(G₁₃) ≤ 36`, which reads their text alone: it rebuilds G₁₃, its circles and its 4 732 automorphisms from the definitions and agrees with `scripts/g13/g13.py`; it passes the real formulas, with the exact counts of edge, domination, no-rosette and chain clauses and the prefix of `lex_order()`; and it rejects copies with one wrong clause: a non-edge, a moved domination or circle clause, a changed totalizer clause or threshold, a chain whose map is not an automorphism or that skips a moved position, stray units, and others. |
| `test_g13_chi.py` | ✓ | `χ(G₁₃) = 6` (`notes/g13_chi.md`): the three case formulas `F34`, `F35`, `F36` (their clauses, the counts and the chains with the functions of the audit of `α`) on intended models and on colourings far from any normal form; `verify_plan_D.py`, the packer and `scripts/verify_g13_chi.py` on a small made-up run; and, when `certificates/g13_chi_certlogs.tar.gz` is present, the archive against its SHA256SUMS file and the logs of `E37_*` in it against the certificates of `α(G₁₃) = 36`. |
| `test_g17.py` | ✓ | The anisotropic plane G₁₇ of `scripts/g17_alpha.py`: 289 vertices of degree 18, and rotations, reflections and translations are automorphisms; its independent circles are N = 4, 5, 9, 11, 12, 14, 15, with regions of 90 or 108 vertices; the 57-point rosette is independent, contains 0 and the whole circle N = 12, and breaks exactly one clause of the part-B formula; the totalizer is exact on every assignment of up to 8 literals; the part-A formulas are, byte for byte, those whose DRAT proofs `certificates/g17_part_a_checks.txt` records as verified. |
| `test_g17_slow.py` |  | The seven part-A formulas are unsatisfiable, solved again by CaDiCaL (marked `slow`; about six minutes on a loaded machine). |
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
| `test_quadratic_planes.py` | ✓ | The planes over ℚ(√d) for d = 11, 23, 47, 59, 71 (`notes/quadratic_planes.md`), with `scripts/verify_quadratic_planes.py`: each d is ≡ 11 (mod 12); each graph's points are distinct, its edges have length exactly 1 and are all the pairs at distance 1, and it has no triangle; the 4-colouring and the 3-colourings of every vertex-deleted graph are proper; each stored CNF is exactly its graph's formula, and the logs of both encodings say UNSATISFIABLE and VERIFIED; the colourings of 𝔽₇² and 𝔽₁₁² behind the upper bounds; exact unit lengths; tampered data are rejected. |
| `test_quadext.py` | ✓ | Arithmetic in ℚ(ζ₁₅, √−7, √−11), built as a tower of quadratic extensions of ℚ(ζ₁₅), which contains de Grey's rotations and ζ₁₅. |
| `test_reduce11.py` | ✓ | The unit-distance graph on the field ℚ(√−3, √−11, √−247), viewed in ℂ, is 5-colourable: reduction at a place above 11 maps it to the Cayley graph of 𝔽₁₂₁ on its 12 elements of norm one, which is 5- but not 4-colourable, and this colours the 803-vertex graph `five_247_c.json` properly at both places above 11. |
| `test_ring_geometry.py` | ✓ | Six facts behind the ring constructions, recomputed: the 90° rotation about the centre of a unit square moves one diagonal onto the other at distance 1, the map σ satisfies \|σu\|² = \|u\|²/3 and \|u − σu\| = \|u\|, the two-ring configuration is 3-chromatic, within each component of a neighbourhood in `five_247_c.json` two neighbours lie on the same side of the bipartition exactly when the angle between them is an even multiple of 60°, 2 is not of the form a² + ab + b², and that graph has exactly two vertices of degree 4. |
| `test_slack.py` | ✓ | With slack s = k − 3: unit-distance graphs have clique number 3, a pair forced alike under k colours needs at least k + 1 vertices, a unit rhombus attains this at k = 3, and at k = 4 neither the rhombus nor the Moser spindle forces a pair. |
| `test_spindle.py` | ✓ | The spindle constructions (two-copy, triple and local), separation tests and cores, on cases whose answers are known. |
| `test_small_plane_colourings.py` | ✓ | The colourings behind the upper bounds for the small finite planes, stored in `data/small_plane_colourings.json`: G_q for q = 3, 4, 5, 7, 8, 13 and H_q for q ≤ 16, q ≠ 13, each rebuilt and checked on every edge. Also the lower bounds that need no solver: G₃ has an odd cycle, G₄ is the Clebsch graph (independence number 5), Hoffman's bound for G₅, and K₄ in G₈. |
| `test_split_places.py` | ✓ | Proposition C of `notes/local_colourings.md` §12: the hyperbola graphs H_q = Cay(𝔽_q², {(t, 1/t)}) have χ = 2, 3, 4, 3, 4, 4, 3, 4, 5, 4 for q = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16 (an explicit 5-colouring for q = 13); over ℚ(√3, √5), reduction at the prime above 2 is a proper 4-colouring of `chain35.json`, whose edge vectors are units there, one of its edge vectors is not a unit above 3, and reduction above 3 followed by a 3-colouring of H₉ 3-colours a graph whose edge vectors are units above 3. |
| `test_threepoint.py` | ✓ | The three-point bound of `scripts/threepoint.py` (§14), without a solver: for q ≤ 13 the rotation blocks and the localizing blocks are compressions of the explicit matrices by an orthonormal basis, so they are positive semidefinite exactly when the matrices are; real independent sets of 𝔽₁₁² and 𝔽₁₃² satisfy every constraint, with objective \|S\|; the stored certificates match `data/threepoint/SHA256SUMS`. |
| `test_threepoint_certificates.py` |  | Each certificate in `data/threepoint/` proves the lower bound on χ listed in the file (χ ≥ 6 for 𝔽₃₇², 𝔽₄₁², 𝔽₄₃², 𝔽₄₇², G₂₉, G₃₇, G₄₁; χ(G₁₃) ≥ 5), checked by `scripts/threepoint_verify.py` in interval and exact rational arithmetic (marked `slow`: up to four minutes each). |
| `test_threepoint_indep.py` | ✓ | The independent checker `scripts/threepoint_verify_indep.py`: on 𝔽₇², G₁₁ and G₁₃ the true z of an actual independent set satisfies every constraint it rebuilds; misreading a certificate's labelling (the sign of the basis shifts, or the frequency of each block) more than doubles the bound, so the check is not vacuous; `inert13.npz` and `inert29.npz` give α(G₁₃) ≤ 42 and α(G₂₉) ≤ 163, and the other six certificates their bounds (marked `slow`, and left out of CI: up to two minutes each). |
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
  2 = 𝔭₁²𝔭₂² with residue degree 1, and, through the subfield ℚ(√−3, √−11), that both places are
  inert in ℚ(√3, √11, i);
- 810 unit vectors generated by rotations, and 300 random unit vectors
  ((a² − b²)/(a² + b²), 2ab/(a² + b²)), a form that every unit vector takes by Hilbert 90, are
  2-adic units with nonzero residue in 𝔽₄ at both places;
- the colouring is proper on 300 random unit steps, on the Moser spindle (which has no proper
  3-colouring, the lower bound), on Exoo–Ismailescu's 214-point graph H (`data/ei_H214.json`) and
  on H with its images under ρ₇ and ρ₇⁻¹ (`data/ei_rho7.json`, 638 points and 3 012 edges);
- on random samples, it gives pairs at distance 8/9ⁿ (n ≤ 3) the same colour and pairs at distance
  √(11/3) different colours.

`test_q23.py`, for χ(ℚ(√2, √3)²) = 4, a case that Fischer's hypotheses exclude, checks:
- the local facts: ℚ(√2, √3) has one place over 2, totally ramified with residue field 𝔽₂ and not
  split in ℚ(√2, √3, i), whose completion has residue field 𝔽₄; sympy's `prime_decomp` confirms
  that 2 = 𝔭⁴ in ℚ(√2, √3), and 2 = 𝔓⁴ with residue field 𝔽₄ in ℚ(√2, √3, i) = ℚ(ζ₂₄);
- 312 listed unit vectors, and 300 random unit vectors of the Hilbert 90 form above, have nonzero
  residue;
- the colouring is proper on the 2 089 sums of at most three 24th roots of unity and on a graph of
  3 134 sums of two listed unit vectors;
- the set M₂ = M₁ + M₁ of Voronov, Neopryatnaya and Dergachev (arXiv:2106.11824, Series 2), with
  2 593 points and 11 448 edges, has no proper 3-colouring and is 4-coloured by the residue
  colouring;
- the chain of three unit rhombi in `data/chain23.json` (10 vertices, 16 edges) has no proper
  3-colouring, the lower bound.
