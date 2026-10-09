# Changelog

Notable changes to this repository, newest first. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers follow
[Semantic Versioning](https://semver.org/) for the library, the tools and the
data formats.

Mathematical corrections are listed in each version. The full record,
including every retracted claim, is the research log,
[`docs/research-log.md`](docs/research-log.md).

## [Unreleased]

### Changed
- **Six colours on the ruler-and-compass plane (9 October).** OpenAI's proof that χ(ℝ²) ≥ 6 runs over the
  constructible numbers, the smallest subfield of ℂ closed under square roots. Every 5-colouring of the points
  constructible with ruler and compass has two points at distance 1 with the same colour, so some finite unit-distance
  graph with constructible coordinates has no proper 5-colouring. Formally verified in Lean, off CI; the axioms are
  `propext`, `Classical.choice` and `Quot.sound`. The origami numbers (square and cube roots) and the numbers expressible
  by radicals were done first.
  - The changes to OpenAI's proof: a finite-image lemma that needs only square roots, and the triple average of its
    radial inequality taken at the exponents −1, 0, 1, with two new lemmas for the fixed middle factor.
  - Files: `notes/six_over_fields.md`; `lean/external/openai-five/fields/` (patches against openai/math `fd4aeeb2`, the
    statements, `build_field.sh`, a copy of OpenAI's Apache 2.0 licence); README (state of the art, results, the
    section on six, Spanish summary); `lean/README.md`; `notes/README.md`; research log.
- **`κ(U_1050) = 1/4` over `ℚ(√59)`, certified (9 October).** An exact certificate (1 997 203 nodes; both checkers
  accept it) shows that no character maps the 300 unit vectors with denominator 1050 into `(1/4, 3/4)`; with the
  character `θ` of `at_four/`, `κ(U_1050) = 1/4`, so a finite subgraph of their Cayley graph has `χ_c = 4` (Corollary
  F14). `at_four/cert59_1050_open_4.json.gz`, `check_theta.py` (now checks a certificate's units), note, paper (PDF
  rebuilt), README, two slow tests, research log.
- **`χ(ℝ²) ≥ 6` (OpenAI, September 2026), rebuilt in Lean (9 October).** OpenAI's proof that the plane is not
  5-colourable ([openai/math](https://github.com/openai/math), family 158) is non-constructive and formalized in Lean
  4.34.1 with the Mathlib commit of `lean/`. We compiled the 70 modules of its import closure against our Mathlib
  build, restated the theorem ourselves and checked its axioms (`propext`, `Classical.choice`, `Quot.sound`).
  `lean/external/openai-five/` (the procedure). README (state of the art, the section on six, which is now the search
  for an explicit graph, the Spanish summary), the introductions and bibliographies of the six papers (PDFs rebuilt,
  same page counts), `notes/local_global.md`, `notes/local_colourings.md`, `notes/worker_jobs.md`,
  `notes/rigidity.md`, research log.
- **A character with values in `ℤ/4` at a place above 2 (9 October).** Proposition F15′ (paper, Proposition 11): if
  `F` has a place `v ∣ 2` with `F_v ≅ ℚ₂(√3)`, some homomorphism `ℤT → ℤ/4` maps `T` into `{1, 2, 3}`; so `κ(U) ≥ 1/4`
  for every finite `U` over `ℚ(√59)`, `ℚ(√83)` and `ℚ(√3, √11)` (for quadratic fields the colouring is Fischer's,
  1990, Theorem 10(i)), the tight-square alternative of Corollary 8 never applies there, and a finite witness at 4
  contains a cycle of length at least 8 for each of the twelve such maps. Not so over `ℚ(√2, √3)` modulo 4. Refereed
  (own programs, no error; wording findings applied). `at_four/quarter_character.py`, tests, paper (Section 10), note
  §6.9, research log (including a corrected floating-point value for the 60 vectors of `H₄′`: `κ = 1/4`, not 0.178).
- **Explicit subgraphs of the Cayley graphs with `χ_c = 4` (9 October).** `finite_witness/witness_q3_11_cayley` (1 874
  points, the 7 887 pairs that differ by one of the 54 vectors of the construction, 3 389 cycles) and
  `finite_witness/witness_q2_3_cayley` (1 657 points, 6 199 pairs, 120 vectors, 5 264 cycles) have `χ_c = 4`: drat-trim
  and `cake_lpr` verify the refutations. The first is the finite subgraph of Corollary 8 for the 54 vectors, explicit.
  `check_witness4.py` checks the new `generators` key; `construction_cayley/` reproduces both files and proofs byte for
  byte; tests in `tests/test_at_four.py`. Refereed (`indep_W4c/`: own encoding, kissat, drat-trim, `cake_lpr`, 29
  mutations per file; no discrepancy). Paper (Section 10), note §6.9 and Question 2, README files, research log.
- **Audit of the base-point periods passages (6–7 October, night).** No mathematical error. Corrected: `H₄` is not a
  subgraph of the Cayley graph of the 27 vectors (198 of its edges are other unit vectors, as its clausal cores used all
  unit pairs; 39 edges of `H₄′` lie outside its 120 vectors); the remark replaces only the compactness step of Corollary
  8 (not the one behind Proposition 10) and needs `S` without elements of order 2; the abstract (tight squares; the
  computer-assisted step is also used for the attainment of 4); the `ℚ(√59)` search sentence; only `χ(ℚ(√2, √3)²) ≤ 4`
  is used, and `H₄′` reproves `χ ≥ 4`; definitions before use and notation (`ξ`, `K_{a/b}`, `r_m`, `𝒮`, `|W|`, `w`).
  Paper (31 pages), note §6.9 and Question 2, `finite_witness/README.md`, research log.
- **A second explicit graph with `χ_c = 4`, over `ℚ(√2, √3)` (6 October, night).**
  `data/number_fields/circular/finite_witness/witness_q2_3.json.gz`: 1 657 points (denominator 36), all 6 238 unit pairs
  as edges, `χ_c = χ = 4`; so `χ_c(ℚ(√2, √3)²) = 4` with an explicit graph. Certified by `check_witness4.py` (which now
  reads the field from the file) with drat-trim, and by a referee's own programs (`indep_W4b/`: two encodings, kissat,
  drat-trim, `cake_lpr`, 77 sanity checks). Found with base-point periods from 60 unit vectors `ζ₂₄^j s^l` of `ℚ(ζ₂₄)`,
  `ζ₂₄ = e^{iπ/12}`, `s = (1 + 2√−2)/3`; `construction_q2_3/` reproduces it byte for byte. The docstrings of the
  construction scripts now say Lemma F17 (they said Lemma P). Note, paper (abstract, introduction, Section 10, Question
  3), README files, research log, `tests/test_at_four.py`.
- **The construction by base-point periods always works (6 October, night).** A remark after Lemma F17 (note §6.9;
  after Lemma 24 in the paper): if `κ(U) ≤ 1/4`, some finite set of relations leaves the period system without a
  solution (one relation `ρ_p` for each of the finitely many `p` allowed by the basis, from an integer point of a
  rational cone), and the subgraph of `Cay(ℤU, U ∪ −U)` on the walks of its chains has a tight cycle in every proper
  4-colouring. This replaces the compactness step of Corollary F14 by a finite construction. Refereed (correct;
  wording fixes applied) and tested on 15 distance graphs on `ℤ` (`finite_witness/indep_W4/remark_check.py`, two
  logs). Note, paper (remark, introduction, Question 3; 31 pages), README files, research log.
- **An explicit graph with `χ_c = 4` over `ℚ(√3, √11)` (6 October, night).**
  `data/number_fields/circular/finite_witness/witness_q3_11.json.gz`: 1 874 points, all 8 085 unit pairs as edges,
  `χ_c = χ = 4`, so the value 4 of Theorem F16 is attained by an explicit graph. Certified by `check_witness4.py`
  with drat-trim and by a referee's own programs (`indep_W4/`: second encoding, kissat, drat-trim, `cake_lpr`).
  Found with base-point periods (Lemma F17 of the note, with proof, refereed; `construction_q3_11/` reproduces it). Note
  §6.9 and Question 2, paper (the lemma on base-point periods, abstract, introduction, Question 3), README files,
  research log, `tests/test_at_four.py`.
- **27 explicit unit vectors with `κ = 1/4` over `ℚ(√3, √11)` (6 October, evening).** `R₆₀ʲ R_Aᵏ R_Gˡ (1, 0)`,
  `j, k, l ∈ {−1, 0, 1}` (`R_A` the Moser angle, `R_G` cosine 11/14): a character of order 12 at exactly `1/4`
  and an exact certificate (12 073 nodes) that none maps them into `(1/4, 3/4)`, so a finite `χ_c = 4` witness lies
  in their Cayley graph (Corollary F14) without Theorem F; the witness is not known yet.
  `data/number_fields/circular/at_four/` (`q3_11.py`, `theta311.json`, `cert311_open_4.json.gz`,
  `cert_open_units.py`); `check_open.py` and `check_open_indep.py` accept this biquadratic plane and `P/Q = 4`.
  Note §6.9, paper, research log, `tests/test_at_four.py`.
- **Finite witnesses at four colours (6 October).** At `p/q = 4` the winding argument of Theorem W⁺ fails only on a
  *tight square* `x → x + s → x + s + t → x + t → x`, which is itself a tight cycle; so a 4-colouring without tight
  squares gives a character into `[1/4, 3/4]`, and conversely (Lemma F13). Hence, for a finite connection set with
  `κ(S) ≤ 1/4` and a 4-colourable Cayley graph, some finite subgraph has `χ_c = 4` (Corollary F14), and, with a finite
  form of Theorem F (Proposition F15), whenever `χ_c(F²) = 4` some finite unit-distance graph in `F²` has `χ_c = 4`
  (Theorem F16; for example over `ℚ(√59)`, `ℚ(√83)`, `ℚ(√3, √11)`). This answers the question of the note and of
  Question 3 of the paper, without a bound on the size. Refereed by an independent checker, with two programs of
  its own on every abelian group of order at most 16 (`data/number_fields/circular/at_four/indep/`): no error.
  `notes/circular_planes.md` §6.9, `papers/three-colours/` (subsection "Finite witnesses at four", abstract,
  introduction, Question 3), `papers/winding/` (Proposition `prop:four`). Also: the character
  `θ = (89/118, 1/2, 89/118, 1/2)` keeps the 108 unit vectors of `ℚ(√59)²` with denominator 210 at margin `15/59`, so
  the growth with `D = 210` could not reach `χ_c = 4`; `θ = (3/8, 5/8, 5/8, 5/8)` gives margin exactly `1/4` for
  denominator 1050 (`at_four/`, `check_theta.py`, `tests/test_at_four.py`).
- **The Theorem W scan up to 10 000 (6 October).** `data/quadratic_planes/winding/scan` now holds an exact Theorem W
  certificate for each of the 155 squarefree `d ≡ 11 (mod 12)` below 2000 (the smallest file found for each `d`; 50 MB,
  31 MB of it for `d = 443`), each accepted by both checkers (`scan/checks.txt`), with a README. Between 2000 and 10 000
  the search certified 593 of the 611 values; `scan/beyond2000.tsv` lists all 611, with the branch, git blob id and
  sha256 of the 591 stored certificates (both checkers accept all of them; they stay on the branches
  `claude/winding-share-z1` and `claude/winding-share-z2`, 327 MB) and the two that were checked once and not kept (over
  50 MB). The search tools `single_scan.py`, `union_search.py` and `union_cert.py` are added;
  `tests/test_winding_scan.py` checks the coverage, three certificates (all 155 when slow tests run) and the table.
  Theorem 1 already proves the bound for every `d ≡ 11 (mod 12)`; the scan is an independent check, field by field,
  through Theorem W alone.
- **A fourth field with an explicit `7/2` witness: `ℚ(√911)` (5 October).** The growth from the 327-vertex graph over
  `ℚ(√911)` stops at 873 vertices; deleting vertices with the cores of `kissat` and then CaDiCaL refutations, in blocks
  of adaptive size, leaves a vertex-critical unit-distance graph with 324 vertices, 866 edges and `χ_c = 7/2`
  (`witness_q911.json.gz`), certified like the others (both checkers, two encodings, `drat-trim` on both DRAT proofs,
  `cake_lpr` on the LRAT form of one, criticality certificates). The paper, the note, the READMEs and the test of the
  grown witnesses include it. The referee's programs (`indep_W/`), run by us on it, confirm it as they did
  `witness_q11b` (`indep_W/results/witness_q911.log`).
- **`witness_q11b` now has 155 vertices (5 October, night).** Deleting vertices, in a random order, from the union of
  the 157-vertex witness and a 161-vertex one (another order on the same 205-point union) leaves a vertex-critical
  unit-distance graph with 155 vertices, 404 edges and `χ_c = 7/2`. It replaces the 157-vertex graph under the same
  file names and is certified the same way: both checkers, two encodings, `drat-trim` on both DRAT proofs, `cake_lpr`
  on the LRAT form of one, and criticality certificates for every `H − v`. The paper, the note, the READMEs and the
  test of the grown witnesses now give 155. Its cycle list first repeated 563 of its 1 836 cycles; the repeats are
  removed (1 273 cycles) and both formulas refuted and checked again. The referee's programs (`indep_W/`), run by us
  on it, confirm everything: the points and all unit pairs, a third encoding, cycle lists rebuilt from scratch, every
  formula clause by clause, `kissat`, `drat-trim` and `cake_lpr` on all of them and on the stored proof, and the
  criticality certificates through explicit `(7N, 2N + 1)`-colourings (`indep_W/results/witness_q11b.log`).
- **A full reading of the three-colours paper (4 October, night).** A referee read the whole paper again, with
  priority on Section 10 and the Questions, and reran the cheap checks with programs of its own: no mathematical
  error. Its corrections are applied: the introduction now says that Lemma 19 also uses Proposition 6, in its case
  `f ≥ 3`, which is proved by hand (only the case `f = 1`, `p ≥ 11` rests on a computation, used only for Corollary 5);
  the sentence on the referee of the `7/2` witnesses no longer covers the 157-vertex witness, which only the two
  checkers checked; Lemma 21 needs `p < 4q`, so it is unavailable already at 4; the notation of Lemma 22 no longer
  clashes with `T`; the bound `19/4` now has an explicit map (`4a + 5b` on `μ₂₀ ⊂ 𝔽₃₆₁`); and small presentation
  fixes. The README of the witnesses records that the programs reproduce `witness_q191` exactly.

### Added
- **A 13-vertex witness for `χ_c = 3` over `ℚ(√15)` (5 October, night).** The growth for the value 3 from 70 seeds over
  `ℚ(√15)`, with 174 deletions, gives at best 13 vertices and 18 edges, always the same graph; it is stored
  (`witness_q15.json.gz`, `q15_seed.json`), checked by `check_small.py` on all `3^13` maps and by a refutation that
  `kissat`, `drat-trim` and `cake_lpr` certify, and reproduced by the programs (new tests). None of the three
  nine-vertex graphs appears over `ℚ(√15)` with 13 denominators tried, so the least size there is between 9 and 13;
  Question 3 of the paper and the note now say so.
- **A smaller witness for `7/2` over `ℚ(√11)`: 157 vertices (4 October, night).** Deleting vertices, in a random
  order, from the union of the 170-vertex witness and a second one (from the same growth with at most 100 new points
  per round) leaves a vertex-critical unit-distance graph with 157 vertices, 409 edges and `χ_c = 7/2`
  (`witness_q11b.json.gz`), certified like the others: both checkers, two encodings, `drat-trim` on both DRAT proofs,
  `cake_lpr` on the LRAT form of one, and criticality certificates for every `H − v`. The paper, the note and the
  READMEs now give 157 as the smallest known size; the test of the grown witnesses covers it.
- **The subdivided hexagon over `ℚ(√7)`, and a referee's report on the value-3 witnesses (4 October, night).** `M`, the
  hexagon with its three long diagonals subdivided (9 vertices, 12 edges), is the unit-distance graph of nine points
  of `ℚ(√7)²` (`witness_q7m.json.gz`), and so are `M` plus one or two edges between midpoints: all three nine-vertex
  graphs with `χ_c = 3` are unit-distance graphs over `ℚ(√7)`, and over `ℚ(√31)`. Found by `hexagon_search.py`
  (equilateral hexagons whose long diagonals are sums of two unit vectors; exact arithmetic; it counts the
  realisations of each graph); the stored witness is checked by enumeration (`check_small.py`) and by a SAT refutation
  certified by `kissat`, `drat-trim` and `cake_lpr`. `M` is a witness with the fewest vertices and, among those, the
  fewest edges. A referee checked the value-3 material with programs of its own (`finite_witness/indep_V3/`: all
  `2^28` labelled graphs on 8 vertices, a Burnside count on 9 vertices, an exhaustive `8^8` homomorphism search, the
  proof by hand on all 126 colourings of `M`, the formulas clause by clause): no error. Its findings are applied: the
  seed of the `ℚ(√31)` growth is stored (`q31_seed.json`; the growth reproduces the witness exactly), the remark on
  `ℚ(√15)` and `ℚ(√39)` is corrected, Vince and Bondy–Hell are cited for the numerator bound, the proof by hand names
  the three 6-cycles it uses, and Question 3 is stated more precisely. Tests: `test_finite_witness_q7m`,
  `test_finite_witness_q31_found_by_growth`, and (slow) `test_hexagon_search`.
- **The triangle-free graphs with nine vertices and `χ_c = 3`, and a referee's report on the `7/2` witnesses
  (4 October, night).** Of the 1 897 triangle-free graphs with 9 vertices (up to isomorphism) exactly three have
  `χ_c = 3` (`nine_vertices.py`, two separate tests that agree on every graph): the hexagon with its three long
  diagonals subdivided, and it with one or two edges between the midpoints. The second is the nine-point witness over
  `ℚ(√7)`; the third is the unit-distance graph of nine points of `ℚ(√31)²` (`witness_q31.json.gz`, found by the same
  growth; checked by enumeration and by a SAT refutation certified by `kissat`, `drat-trim` and `cake_lpr`). The proof
  by hand now covers the subdivided hexagon: for the 5-cycles `Z_k`, `Z_{k+1} − Z_k` is a 6-cycle and `Z₀ + Z₃` is the
  hexagon. A referee checked the three `7/2` witnesses with programs of its own (`finite_witness/indep_W/`: a third
  encoding, cycle lists rebuilt from scratch, `cake_lpr` on the stored proofs, explicit `(7N, 2N+1)`-colourings of
  every `H − v`, exact reproduction of two of them): no error. Its remarks are applied: the paper and the README get
  `χ_c = 7/2` over `ℚ(√191)` and `ℚ(√455)` from the witnesses and Proposition 3, without the computer-assisted
  Proposition 9; they describe the second growth rule and the denominators, state the perturbation `N·c + pos`
  explicitly, and document the reproduction; `check_witness.py`, `verify_independent.py` and the other checkers make
  every check explicitly (under `python -O` the old `assert`s were skipped), accept only an exact `s VERIFIED` from
  `drat-trim`, and use a temporary file. Tests: `test_finite_witness_q31`, `test_nine_vertices`,
  `test_witness_checkers_under_python_O`, and the hand-proof test for the hexagon.
- **A nine-point witness for `χ_c = 3` in a plane without unit triangles (4 October, night).** Over `ℚ(√7)`, where
  `χ_c = 3` and there is no unit triangle, nine points span a unit-distance graph `H₇` with `χ_c(H₇) = 3`: the Wagner
  graph (the Möbius ladder on 8 vertices, which is `K_{8/3}`) with one chord subdivided. Each of its 84 proper
  3-colourings has a tight 6-cycle, by a short proof by hand (a 1-chain identity `2Z = Q₀ + Q₁ + A + B − C` between
  two 4-cycles, a 5-cycle and three 6-cycles); equivalently it has no homomorphism to `K_{8/3}`. It is vertex-critical, and nine
  vertices are the fewest possible: every triangle-free graph with at most 8 vertices maps to `K_{8/3}` (enumeration
  of the 4 682 270 triangle-free graphs on 8 labelled vertices). As `√7 ∈ ℚ₃`, it is also a witness for `ℚ₃`. Found
  by the growth and deletion programs, now for any value `p/q` (`grow.py ... p q`; the results for `7/2` are
  unchanged), from 207 points; checked by enumeration (`check_small.py`, no solver) and by a SAT refutation that
  `kissat`, `drat-trim` and `cake_lpr` certify; `small_triangle_free.py`. The three-colours paper (Section 10 and the
  introduction) and `notes/circular_planes.md` §6.8 updated; tests in `tests/test_two_primes.py`.
- **Smaller witnesses for `χ_c = 7/2`, and the first over `ℚ(√455)` and `ℚ(√191)` (4 October, evening).** A colouring-guided growth with
  lazy SAT (`grow.py`: list the tight cycles of each `(7, 2)`-colouring the solver returns; when its tight digraph is
  acyclic, add the points whose neighbours leave no colour), followed by vertex deletion (`minimise.py`), gives a
  unit-distance graph over `ℚ(√11)` with 170 vertices and 468 edges and `χ_c = 7/2` (the previous one had 2 237
  vertices), and ones over `ℚ(√455)` and `ℚ(√191)` with 175 and 293 vertices. All three are vertex-critical: for every vertex `v`, a
  stored `(7, 2)`-colouring of `H − v` with an acyclic tight digraph shows `χ_c(H − v) < 7/2` (`check_critical.py`).
  Each lower bound is certified twice (two checkers sharing no code, now reading `d` from the witness file; two
  encodings; `drat-trim` on both DRAT proofs; `cake_lpr` on one in LRAT form).
  `data/number_fields/circular/finite_witness/`, tests in `tests/test_two_primes.py`; the three-colours paper
  (Section 10, introduction, Question 3), the winding paper and `notes/circular_planes.md` §6.8 and §7 updated.
- **An explicit finite witness for `χ_c(ℚ(√11)²) = 7/2` (4 October).** The unit-distance graph on `A + A`, where `A` is
  the vertex set of the 76-vertex graph of `χ(ℚ(√11)²) = 4` (2 237 vertices, 11 300 edges), has circular chromatic
  number exactly `7/2`: a 7-adic `(7, 2)`-colouring, and a SAT computation that every `(7, 2)`-colouring has a tight
  cycle among 180 listed ones, certified twice (two checkers sharing no code, two encodings, `drat-trim` on both DRAT
  proofs, `cake_lpr` on one in LRAT form). Corollary 7 of `papers/three-colours/` had shown only that a finite witness
  exists. Found by a separate search (iterated Minkowski sums, lazy SAT); `data/number_fields/circular/finite_witness/`,
  tests in `tests/test_two_primes.py`; paper (now 26 pages) and `notes/circular_planes.md` §6.8 updated.
- **Three colours for every number field (4 October).** `notes/three_colours_number_fields.md`, Theorem B:
  `χ(F²) ≤ 3` iff a prime above 2 ramifies in `F(i)` or a prime above 3 has residue degree 1. It contains Theorem 1,
  proves the local–global question at three colours for every number field, and decides the local fields. Checked by
  two separate agents (no error; two small gaps closed). New fields that need four colours, among them `ℚ(√2, √7)`,
  with an exact certificate (`data/number_fields/three_colours/`) and `χ = 4`; an elementary family (Proposition B9);
  tests in `tests/test_three_colours.py`.
- **A circular chromatic number (4 October).** `notes/circular_planes.md`: `χ_c(ℚ(√11)²) = χ_c(ℚ(√35)²) = χ_c(ℚ₇²) = 7/2`,
  by a 7-adic colouring and exact certificates (`data/number_fields/circular/`, two checkers) with Theorem W⁺; no plane
  over a number field has `2 < χ_c < 3`; below 4, the best upper bound from a locally constant colouring of one completion is 2, 3 or `7/2`
  (Propositions C2, C3 and Corollary C4, with Weil's bound for Kloosterman sums); `χ_c = 7/2` for infinitely many
  number fields (Corollary C5).
- **`χ_c(ℚ(√59)²) ≥ 53/15` (4 October).** `notes/circular_planes.md`, Proposition C6: an exact certificate
  (`data/number_fields/circular/cert_sqrt59_open_53_15_N25.json.gz`, 127 181 nodes, accepted by both checkers) shows
  that `χ_c(ℚ(√59)²) > 7/2`, although `χ(ℚ(√59)²) = 4`, like `ℚ(√11)` with `χ_c = 7/2`.
- **A gap above 3 (4 October).** `notes/circular_planes.md` §5, Theorem D: if `χ(F²) ≥ 4` then `χ_c(F²) ≥ 56/17`, and
  `χ_c(F²) > 3.3315` with a computer-assisted step (exact computations by three programs written separately,
  `data/number_fields/circular/probe/`, `probe/indep/` and `probe/indep2/`). Refereed twice by separate agents (no
  error; their corrections applied).
- **The circular chromatic number below 4 (4 October).** `notes/circular_planes.md` §6, Theorems E and F: for every
  number field `F`, `χ_c(F²)` is 2, 3, `7/2` or at least 4; it is `7/2` exactly when `χ(F²) ≥ 4` and some place of
  `F` above 7 has residue degree 1. So `χ_c(ℚ(√23)²) = 7/2`, `χ_c(ℚ(√59)²) = 4` and `χ_c(ℚ(√47)²) ≥ 4`. The proof adds
  the rotation `σ = (5 + 12i)/13` to the probe of Theorem D; one computer-assisted step for each theorem (exact
  certificates for the windows modulo 65 and 325 with checkers written separately,
  `data/number_fields/circular/twoprime/`). Found by a research agent and refereed by two others, each with two exact
  methods of its own (no error; their corrections applied). Tests in `tests/test_two_primes.py`. A third referee read
  Section 10 of the paper: no error in a proof; its corrections are applied, and eight new dual certificates
  (`cert_K2M1_seven.txt`, with a separate checker) make Theorem F depend on the window modulo 325 alone; the referee's programs are in `twoprime/indep_S10/`. Corollary 6
  of the paper (`notes/circular_planes.md`, Corollary F8): the same values for the local fields, so below 4 `χ_c(F²)`
  is the least `χ_c` of the planes over the completions of `F`; and `4 ≤ χ_c(ℚ(√47)²) ≤ 19/4`. Corollary 7 of the
  paper (note Corollary F12): below 4, `χ_c(F²)` is the circular chromatic number of a finite subgraph of `F²`. It
  rests on a general lemma (paper Lemma 22, note Lemma F11, and Section 7 of `papers/winding/`): for an abelian
  Cayley graph with a finite connection set and `2 < χ_c < 4`, every homomorphism to `K_{χ_c}` has a tight cycle
  (an optimal character has a nonnegative relation among the elements where it is tight, as otherwise it could be
  perturbed), so `χ_c` is attained by a finite subgraph; for the planes a finite set of unit vectors with
  `κ = 1/χ_c` exists by compactness. A first proof through the certificates of value `2/7` was refereed and then
  replaced by this one, refereed separately (no mathematical error in either; programs in `twoprime/indep_FW/` and
  `twoprime/indep_FW2/`). Corrected: the paper said that no finite unit-distance graph with `χ_c = 7/2`
  was known, but the Moser spindle has `χ_c = 7/2`.
- **Paper draft `papers/three-colours/` (4 October).** Theorem B, Proposition B9, the corollaries, Theorem C,
  the local circular values, Theorem D (Section 9, by hand) and Theorems E and F (Section 10); refereed (no
  mathematical error), fixes applied; Section 10 refereed separately as written (no error in a proof; fixes applied).
- **Literature sweep (4 October).** `notes/literature.md`, last section: what the systematic search of arXiv found for
  each result, and the relevant papers the project did not cite before.

### Changed
- **A full reading of `papers/four-colours/` and `papers/quadratic-planes/` (4 October, evening).** A separate agent
  checked every proof step of Theorem 1 and its corollaries by hand and with its own exact programs (correct, no gap),
  and every number and graph of the quadratic-planes draft against the data (all agree; CaDiCaL refutes all 27
  formulas). Applied: the remark that the bound of Theorem 4(1) is sharp holds for `k ≤ 5` and not for `k ≥ 6`
  (for `N = 5⁶` the 22 values `65 639 ≤ d ≤ 66 143` admit no character although `d ≥ 21N/5`; checked here with the
  referee's exact decision, now `family/fc_case23.py`, and the binding rotation `ρ⁶` recomputed); "26 fields" was 27;
  Moorhouse's Theorem 8.1 cited correctly; the checks section says which certificates use the `N` of the proof; the
  second certificate checker the paper mentioned is now in the repository (`check_w_indep.py`, written separately by
  a referee; it accepts all 26 stored certificates and rejects six kinds of corruption, also under `python -O`); a false
  explanatory claim in the quadratic-planes draft (on denominators for `d = 83`) replaced by what holds; a gap in
  its Hoffman remark closed (exact eigenvalues at 23 and 31); the quadratic-planes draft now records, in dated
  remarks, what the later drafts proved; the repository pointers say the files are on the development branch;
  references (Davies, de Bruijn–Erdős, Exoo–Ismailescu, Isbell, the title of MildlyMeticulous's repository) and
  wording. Both PDFs rebuilt (6 and 7 pages, no warnings).
- **A full reading of `papers/winding/` (4 October, afternoon).** A separate agent read the whole paper,
  reran every certificate and wrote its own checker and tests: no false theorem and no broken proof. Applied: Question
  14 no longer asks as open what `papers/four-colours/` proves (it now asks the four-colour local–global question);
  the raw wind is not invariant under backtracks, so the exponent-4 argument uses the normalised wind, invariant
  under backtracks and under every replacement of two consecutive steps by two with the same sum (the move of Krebs
  and Sankar's homotopy; a second reading found that exchanges alone do not suffice); the remark that a locally constant colouring at a place gives a character now proves integrality at the
  place; related literature added and checked (Wrochna, Matsushita, Gao–Jackson–Krohne–Seward, Heuberger,
  Gujgiczer–Naserasr–S–Taruni, Ryabchenko, Berger, Robinson, Youngs, Day); an explicit 5-colouring of the plane over
  `𝔽₁₁` in Theorem 12; bibliography (de Bruijn–Erdős pages, also in two other papers); the abstract's novelty claim
  hedged; the certificate for 28 of the 54 vectors for `d = 83` regenerated, checked and stored
  (`data/quadratic_planes/winding/cert_83_510_min28.json.gz`, in the slow tests); smaller points. The PDF is rebuilt
  (12 pages).
- **A full reading of `papers/three-colours/` (4 October, afternoon).** A separate agent read the whole paper
  and found no mathematical error in a proof; its corrections are applied. Lemma 13 (the arithmetic of `𝔽₄₉` at the
  places above 7: `A′₇ = {e : eē = −1}`, no coset of a subgroup inside it, the digit patterns) now has a proof by hand
  (note §6.4), so Propositions 8 and 9 are the only computer-assisted steps of Theorems E and F; the sentence on
  `ℚ(√47)` now says `4 ≤ χ_c ≤ 19/4`; the remark that Corollary 2 is contained in `papers/four-colours/` holds for real
  quadratic fields only (`ℚ(√−73)` is split above 2 and 3); `χ(ℚ(√2, √7)²) ≤ 4` is stated in Section 7; the radius-2
  ball of the 140 vectors of Theorem C is now described exactly (`twoprime/finite_ball.py`: the Cayley part is
  bipartite, the induced unit-distance graph has `χ_c = 5/2`); attributions (Isbell, Fischer 1994, the special cases
  of (a)), references (Exoo–Ismailescu, Heule, Parts, Soifer), notation and wording are corrected.
  A second agent then refereed the hand proof of Lemma 13 (correct, no gap; its programs in `twoprime/indep_L13/`,
  with a test) and the new sentences; its corrections (Proposition 6 rests on a computation of `κ₁` for
  `11 ≤ p < 1001`, now said in the introduction, and ten minor points) are applied.
- **Theorem D in the paper (4 October).** It now states only the bound `56/17`, proved by hand; the computer-assisted
  bound `3.3315` is a remark, as Theorem E supersedes it. `data/number_fields/circular/twoprime/seven_patterns.py`:
  three checks that had produced the last lines of its stored output, but were missing from the stored program, were
  written again.
- **Credit for earlier work (28 September).** The upper bound of `χ(ℚ(√2, √3)²) = 4`
  and of Fischer's `χ(ℚ(√3, √11)²) = 4` is also a case of Corollary B′ of
  [hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction), a public,
  unrefereed repository of July 2026, whose Theorem A is the same reduction of
  `z = x + iy` at a place over 2. We found it after our note. The front page,
  `notes/literature.md` and the note (version 6) now say so; what remains ours is
  the explicit case `ℚ(√2, √3)` with its lower bound, and the Lean proofs.

### Changed

- **The winding paper after a second referee** (3 October). Correct mathematics; attribution fixed: Steps 1–2 for
  `K_{p/q}` are the wind of Brewster, McGuinness, Moore and Noel and of Brewster and Moore (for `K₃`, Krebs and Sankar);
  the exponent-2/4 corollary is a uniform proof of known facts; the formal-proof remark says exactly what Lean checks;
  citations and journal data updated; the character on all of `G` added in Lean (`theoremWplus_extended`).

### Added

- **Theorem 1 and Corollary 3 in Lean** (3 October, night). `lean/FourColours.lean` proves
  `χ(ℚ(√d)²) ≥ 4` for every `d ≡ 11 (mod 12)` (`FourColours.not_colorable_three`), `lean/PadicFour.lean` proves
  `χ(ℚ_p²) ≥ 4` for every prime `p ≥ 5` and for every field of characteristic 0 containing a square root of some
  `d ≡ 11 (mod 12)`, and `lean/TheoremWInf.lean` proves Theorem W for every abelian group and finite `S` (box averages
  and a limit along an ultrafilter; also for `K_{p/q}`, `p < 4q`). Only Lean's three standard axioms. The Lean
  workflow's checks (build, axiom diff, kernel replay) passed on ec7b021a, run on a separate machine because the
  month's CI minutes were used up (`lean/VERIFY_ec7b021a.md`). The formal proof uses a cruder form of the structure lemma (described in
  §7 of the note and in the paper's checks section), with `data/quadratic_planes/winding/family/crude_check.py`
  for its numerical facts.

- **Four colours are needed for every `d ≡ 11 (mod 12)`** (3 October, evening). `χ(ℚ(√d)²) ≥ 4` for every
  `d ≡ 11 (mod 12)`, by a proof by hand from Theorem W (`notes/four_colours_11_mod_12.md`, draft
  `papers/four-colours/`). So a real quadratic plane needs four colours exactly when `d ≡ 11 (mod 12)`, and
  `χ(ℚ(√d)²) = 4` for all these `d` except possibly `d ≡ 47, 143, 167 (mod 168)`. The sets of unit vectors are
  explicit (Fischer's vector `(1 + i√d)²/(1 + d)`, its mirror image and, for `d ≡ 11 (mod 24)`, one more vector,
  rotated by the rational rotations with denominators dividing `5^k`); the bound on `k` is sharp for
  `d ≡ 23 (mod 24)`. Two internal referees checked the proofs (no error in the proof of Theorem 1; one wrong
  generality in the statement of Theorem 1b, for negative `n`, corrected). New: `data/quadratic_planes/winding/family/`
  (exact checks of the structure lemma, the configurations, the large-`N` criterion, 19 exact certificates) and
  `tests/test_winding_family.py` (in CI, except the slow certificates). Corollary: `χ(ℚ_p²) ≥ 4` for every prime
  `p ≥ 5` (choose `d ≡ 11 (mod 12)` with `d ≡ 1 (mod p)`, so `ℚ(√d) ⊂ ℚ_p`); checked by the second referee; the
  p-adic draft has a remark on it.

- **Liu's Problem 3 in Lean** (`lean/DistLiu.lean`, 3 October). For every set `D` of three positive integers,
  `G(ℤ, D)` maps to `K_{p/q}` if and only if some `α` has `‖dα‖ ≥ q/p` for all `d ∈ D`, that is,
  `χ_c(G(ℤ, D)) = 1/κ(D)` (`DistLiu.liu_problem3_iff_unconditional`). Proved from `lean/TheoremWplus.lean` through
  one periodic window, with Theorem W⁺ for every distance graph with `2q ≤ p < 4q` (`DistLiu.wplus_distance`), its
  converse, and a short proof of the lonely runner theorem for three real speeds (`DistLiu.lonely_runner_real`),
  also written out in the paper (Lemma 10). Standard axioms only; CI builds it, checks its axioms and replays it.
- **Theorem W⁺ and Liu's Problem 3** (`papers/winding`, `notes/winding_lemma.md`, 3 October). For `p/q < 4`, an
  abelian Cayley graph maps to `K_{p/q}` if and only if a character maps the connection set into `[q/p, 1 − q/p]`;
  so `χ_c = 1/κ` below 4, with 4 sharp. Consequence: `χ_c(G(ℤ, D)) = 1/κ(D)` and `χ(G(ℤ, D)) = ⌈1/κ(D)⌉` for every
  set `D` of three distances, a negative answer to Problem 3 of Liu's 2008 survey on distance graphs. Refereed
  inside the project (correct; no earlier publication found). The paper is rewritten around the theorem (9 pages);
  `tests/test_winding.py` compares `⌈1/κ⌉` with Zhu's formula for all three-distance sets up to 30.
- **Theorem W for circular cliques below 4, in Lean** (`lean/TheoremWplus.lean`, 3 October). For `p < 4q`, a
  homomorphism from a Cayley graph of a finite abelian group to the circular clique `K_{p/q}` gives a character `ξ`
  of the subgroup generated by `S` with `ξ(S) ⊆ [q/p, 1 − q/p] (mod 1)`, and conversely. The lift `δ ∈ [q, p − q]` of
  each edge makes squares rigid (two sides differ by a multiple of `p` of size at most `2(p − 2q) < p`), and `ξ(s)`
  is the average of `δ(·, s)` over the group divided by `p`. For `p = 3, q = 1` this is Theorem W. Tests compare it
  with SAT on random Cayley graphs and check that it fails at `p/q = 4` (`K₄`). The consequences (the circular
  chromatic number below 4) go into the note after the referee's report.
- **Katznelson's question for three colours, in Lean** (`lean/Recurrence.lean`, 3 October). `GKR.question3`: for
  every 3-colouring of `ℕ` there is a real `α` such that every `n ≥ 1` with `‖nα‖ < 1/3` is a difference of two
  numbers of the same colour (Question 3 of Glasscock, Koutsogiannis and Richter, Bull. Amer. Math. Soc. 59
  (2022)); `GKR.chromatic_recurrence`: the same for `ℤ` and any set that meets every `{s : ‖sα‖ < 1/3}`. Proved from
  `lean/TheoremW.lean` by periodic windows, Theorem W on `ZMod P` and compactness; only the three standard axioms.
- **Aperiodic colourings for many colours** (`notes/winding_lemma.md`, `papers/winding`, 3 October). The
  proper `|F|`-colourings of `Cay(Γ, (F − F) ∖ {0})` are the tilings by `F`, so Greenfeld and Tao's counterexample
  to the periodic tiling conjecture gives Cayley graphs of `ℤ^d` with no periodic colouring with `χ` colours: a
  negative answer to Problem 4.6 of Abrishami et al. for these graphs, against the positive answer for `χ ≤ 3` that
  Theorem W gives. A referee of the two corollaries found them correct and asked for fixes, now applied: "finite
  connection set" in the periodicity statements, a complete argument that the radius `1/3` is sharp, and smaller
  points (research log).
- **Theorem W in Lean** (`lean/TheoremW.lean`, 3 October). The case of three colours and finite abelian
  groups is formalised in Lean 4 with Mathlib: a proper 3-colouring of `Cay(G, S)` gives a character `ξ` of the
  subgroup generated by `S` with `ξ(S) ⊆ [1/3, 2/3] (mod 1)`. The proof follows the note (winding sums, squares
  do not wind, averaging over the group); no `sorry`, only Lean's three standard axioms, replayed by
  `leanchecker` in CI.
- **Theorem W: three colours and characters** (`notes/winding_lemma.md`, 3 October). A Cayley graph of an
  abelian group is 3-colourable if and only if some character maps every generator into `[1/3, 2/3]` (and maps to
  `C_{2k+1}` if and only if some character maps them into `[k/(2k+1), (k+1)/(2k+1)]`). Steps 1–2 of the proof are
  the discrete winding number of Krebs–Sankar; the averaging step is new. It gives Payan's theorem, the exponent-4
  case of Krebs–Sankar, for distance graphs `χ ≤ 3` iff `κ(D) ≥ 1/3`, and the three-colour case of Katznelson's
  question: every set of Bohr recurrence is a set of 3-chromatic recurrence, and for every 3-colouring of `ℕ` the
  union of the difference sets of the classes contains `{n : ‖nα‖ < 1/3}` for some `α` (Question 3 of Glasscock,
  Koutsogiannis and Richter, Bull. Amer. Math. Soc. 59 (2022)). Also: 3-colourable Cayley graphs of finitely
  generated abelian groups have periodic 3-colourings, and 3-colourability is decidable. Over `ℚ(√d)` it turns
  3-colourability into a question on a 4-torus, settled by exact branch-and-bound certificates
  (`data/quadratic_planes/winding/`: `certify_w2.py`, the independent checker `check_w.py`, six certificates):
  `χ(ℚ(√d)²) = 4` for `d = 83, 107, 203`, `4 ≤ χ(ℚ(√143)²) ≤ 5`, `χ(ℚ(√167)²) ≥ 4`, the fields the graph searches
  could not settle. `tests/test_winding.py` (in CI) checks the theorem on random Cayley graphs and distance graphs
  and four of the certificates (the other two are marked slow), and that `check_w.py` rejects mutated certificates.
  A separate agent refereed the note; its corrections are applied (attribution, an invariant mean in place of an
  ergodic measure, the `κ(D)` endpoint statement, the controls, a checker that raises errors instead of asserting
  and accepts only integers and exact fractions). Paper draft: `papers/winding/`.

- **Paper: *Two-colourable planes over number fields*** (`papers/two-colour-planes/`, 5 pages, draft of
  3 October): `χ(F²) = 2` if and only if some prime of `F` above 2 ramifies in `F(i)`, its corollaries, and the
  local–global question. The "if" direction is Theorem A′ of hn-2adic-obstruction, Lemma 2 is Fischer's
  Theorem 1(iii), and the multiquadratic case is that repository's Corollary B with Fischer's theorem; the
  converse (Proposition 4) is the new part as far as we know.

- **A local–global question** (`notes/local_global.md`, 3 October): is `χ(F²)` the least of
  the local chromatic numbers? Proved for two colours: for every number field `F`, `χ(F²) = 2` iff a prime of
  `F` above 2 ramifies in `F(i)` (the "if" half is Theorem A′ of hn-2adic-obstruction and our Proposition A;
  the converse is new as far as we know),
  which contains Fischer's and Moorhouse's cases (`two_colour_criterion.gp`, `odd_walks.gp`; refereed by a
  separate agent). At four colours the question would give a triangle-free 5-chromatic
  unit-distance graph in the plane, through `ℚ(√167)`, the smallest admissible real quadratic field: no
  local 4-colouring at any place and no unit triangle (`admissible.py`); `ramified_levels.py` and `hyperbola_plane.py` check the local cases used.

- **Real quadratic planes that need four colours** (`notes/quadratic_planes.md`,
  `data/quadratic_planes/`): `χ(ℚ(√d)²) = 4` for
  `d = 11, 23, 35, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959`,
  and `4 ≤ χ(ℚ(√47)²) ≤ 5`; so `χ(ℚ(√d)²)` is now known for every squarefree
  `d < 83` except 47. For real quadratic fields the values known before were 2
  and 3, and a field can need four colours only if `d ≡ 11 (mod 12)`. Each lower
  bound is a triangle-free, vertex-critical unit-distance graph (76, 393, 580,
  816, 406, 611, 1 404, 399, 356, 1 281, 259, 96, 338, 291, 394, 715, 331, 703,
  71, 835, 659, 712, 898, 538, 327, 252 and 513 vertices, for d = 11, 23, 35,
  47, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 443, 455,
  491, 599, 611, 791, 851, 911, 935, 959) with no 3-colouring: kissat with DRAT
  proofs checked by drat-trim, twice, with separate encodings. The upper bounds
  are known: Moorhouse's reduction at 7
  (`d = 11, 23, 35, 71, 95, 119, 155, 179, 191, 239, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959`),
  Fischer's Theorem 10 (`d = 59, 131, 251`) and the reduction at 11 (`d = 47`).
  The note also observes that Cohen's conjecture (2007) on the sets `ℚ(√−d) ⊂ ℂ`
  is a different question, which reduction at a ramified prime settles.
  `scripts/verify_quadratic_planes.py` checks everything again (with
  `--cake-lpr`, every proof is also checked by the verified checker cake_lpr; it
  accepted all of them; its log names each formula by its SHA-256 hash), and
  `tests/test_quadratic_planes.py` runs the fast checks. The graph over ℚ(√11)
  is drawn in the paper (Figure 1) and on the front page (Figure 3). As a
  corollary, a real number field that contains one of these `√d` and has a place
  with residue field `𝔽₇` has `χ = 4`; for example `ℚ(√11, 7^{1/m})`, of degree
  `2m`, so every even degree occurs. For `ℚ(√47)` the note narrows down where a
  4-colouring by reduction at one place could come from: not where `−1` is a
  square in the completion (no locally constant colouring with finitely many
  colours exists there), and not above 47 or any prime `p ≥ 23` (Hoffman's bound
  excludes every level of the local plane at once). That leaves the 11-adic
  plane, whose levels 1, 2 and 3 need five colours (an obstruction of 69 points
  at level 1 lifts to 244 points at level 2 and to 29 524 points at level 3,
  `data/quadratic_planes/padic11.json`), and the 19-adic plane, whose level 1
  needs five. The same graphs give `χ(ℚ₇²) = 4`; `χ(ℚ₂²) = 2` and `χ(ℚ₃²) = 3` follow from Madore's paper
  (arXiv 1509.07023; we had first called all three the first exact values), and `4 ≤ χ(ℚ₁₁²), χ(ℚ₁₉²) ≤ 5`
  (note §5, `data/quadratic_planes/scripts/padic_planes.py`); the three exact values are proved in Lean
  (`lean/PadicPlanes.lean`). For measurable colourings the `p`-adic planes need more and more colours:
  `χ_m(ℚ_p²) ≥ 1 + (p + 1)/(2√p)` for `p ≡ 3 (mod 4)`, so at least 5 from `p = 23`, 6 at `p = 59` and from 67, 7 at 71
  and from 103 (`padic_measurable.py`); for the real plane six measurable colours is open. This answers
  Question 1 of Bardestani and Mallahi-Karai (arXiv 1507.05300), the `p`-adic Hadwiger–Nelson problem, in the
  negative: the Borel chromatic number of `ℚ_p²` is not bounded over `p ≡ 3 (mod 4)`; the coset colouring also
  gives `χ_Bor(ℚ_p²) ≤ χ(𝔽_p²)`, which is at most `(p + 1)/2` for `p > 3` (Le Anh Vinh), where they had `O(p²)`. With Davies's theorem (arXiv 2308.16885),
  `χ(ℚ_p²) = ∞` for `p ≡ 1 (mod 4)`, so `χ(ℚ_p²)` is finite exactly for `p = 2` and `p ≡ 3 (mod 4)`; and by
  Chebotarev's theorem no finite set of number fields gives `χ(ℚ_p²) ≥ 4` for every `p ≡ 3 (mod 4)` (note §5). The 27
  certified fields reach every prime `p ≡ 3 (mod 4)` from 7 below `2 129 503 819`, the first one they miss
  (`padic_reach.c`, checked by `padic_reach.py`); below `3·10¹⁰` they miss seventeen, as the density `2⁻²⁵` from
  Chebotarev's theorem predicts (the 27 values of `d` have rank 25 modulo squares). A seven-page draft, *Colouring the p-adic plane*
  (`papers/padic-planes/`), proves these results; it also determines in every dimension `n` when `χ(ℚ_p^n)` is
  finite (`n = 1`, or `n = 2` and `p ≢ 1 (mod 4)`, or `p = 2` and `n ≤ 4`), with `χ(ℚ₂⁴) = 4` through the residue
  field `𝔽₄` of the 2-adic quaternions, and `χ(ℚ_p¹) = 2` but `χ_B(ℚ_p¹) = 3` for odd `p`.
- **`χ(ℚ(√d)²) = 4` in Lean 4 for ten fields**, `d = 11, 119, 131, 179, 191, 251,
  431, 455, 911, 935` (`lean/Sqrt{d}.lean`), with Mathlib. Each lower bound is the graph
  of `data/quadratic_planes/q{d}.json`: the kernel checks its unit distances and,
  through Mathlib's `lrat_proof`, the LRAT proof `data/quadratic_planes/q{d}.lrat`
  that its formula `q{d}.cnf` is unsatisfiable. The upper bounds are proved once,
  in `lean/QuadraticPlanes.lean`: Moorhouse's reduction at 7 when `d` is a nonzero
  square modulo 7 (Hensel's lemma for √d) or `d = 7d′` (`d = 119, 455`), and
  Fischer's reduction at 2 when `d ≡ 3 (mod 8)` (`d = 131, 251`), each from a
  valuation subring with the right residue field. The proofs depend only on Lean's
  three standard axioms; CI builds them, compares the axioms and replays them with
  `leanchecker`, and a test checks that `lean/tools/field_lean.py` writes the same
  files from the data.
- **The other graphs in Lean, given their formula's unsatisfiability.** For the other
  fourteen fields and the lower bound over `ℚ(√47)`, whose LRAT proofs are too large for
  `lrat_proof`, `lean/Sqrt{d}.lean` proves the theorem from the hypothesis that the
  graph's colouring formula is unsatisfiable. `lean/ColouringFormula.lean` defines that
  formula as a Lean object and proves that a 3-colouring would satisfy it; each file
  checks by evaluation (`#guard`), when it is built, that `q{d}.cnf` is exactly this
  formula, and
  `scripts/verify_quadratic_planes.py --cake-lpr` has the verified checker cake_lpr
  check an LRAT proof of `q{d}.cnf` itself. Typed coordinates and a balanced tree of
  points make every field file build in seconds to minutes.
- **A Moser-spindle-free 5-chromatic unit-distance graph with 852 vertices**
  (`notes/flat852.md`, `data/flat852/`). Its points lie in `ℚ(ζ₂₁)`, and its
  edges use the 84 unit vectors of J. K. Haugland's heptagon graph
  ([arXiv 2608.04542](https://arxiv.org/abs/2608.04542), 2 131 vertices) and
  their mirror images. Haugland's directions alone are 4-colourable by
  reduction at the places above 2, which is why the mirror is needed. No
  4-colouring: kissat with a DRAT proof checked by drat-trim, twice, with
  separate encodings. The field contains no `√−11`, so the graph has no Moser
  spindle. The smallest spindle-free example we found (the next has 1 299
  vertices); the record with spindles allowed is still Parts' 509. The graph is
  not known to be vertex-critical. `scripts/verify_flat852.py` checks it again,
  and `tests/test_flat852.py` runs the fast checks.
- **`χ(G₁₃) = 6`** (`notes/g13_chi.md`): the anisotropic plane over 𝔽₁₃ has no
  proper 5-colouring. A computer proof. The largest class of a 5-colouring would
  have 34, 35 or 36 points (`α(G₁₃) = 36`); in a colouring where it is as large as
  possible it is dominating, an automorphism makes it lex-leader, and renaming
  the other colours gives them value precedence. One formula for each size says
  this; cube and conquer refutes the three, 136 548 leaves in all,
  each by kissat with a DRAT proof checked by drat-trim. `scripts/g13/chi/` is the
  code of the run, `certificates/g13_chi_certlogs.tar.gz` holds its logs, and
  `scripts/verify_g13_chi.py` checks them again. `tests/test_g13_chi.py` reads
  the three formulas (the counts and the chains with the audit of `α`) and tests
  them on intended models, and tests the checkers on a small made-up run. The
  proofs have not been checked by cake_lpr.
- **`α(G₁₃) = 36`** (`notes/g13.md`): no 37 points of the anisotropic plane over
  𝔽₁₃ are independent, so its fractional chromatic number is 169/36. A computer
  proof: a case split on a point with its whole circle, then four formulas and
  a cube tree of 4 822 leaves, each refuted by kissat with a DRAT proof checked
  by drat-trim (`certificates/g13_alpha_*`); a second run had cake_lpr check
  every one of these proofs and the cover (`certificates/g13_cake_lpr_checks.txt.gz`).
  `scripts/g13/` writes the formulas, and `scripts/g13/g13_audit.py` checks
  from their text alone that every clause holds in the intended models;
  `scripts/verify_g13.py` checks everything again; `tests/test_g13.py` checks
  the encodings by brute force and every formula against the logs, and
  `tests/test_g13_audit.py` the audit. `χ(G₁₃)` stays 5 or 6.
- `data/small_plane_colourings.json` and `tests/test_small_plane_colourings.py`:
  the colourings behind the upper bounds for the small finite planes (`G_q` for
  `q = 3, 4, 5, 7, 8, 13`, `H_q` for `q ≤ 16`), found by SAT and until now
  recomputed or not stored, are stored and checked on every edge, with the
  lower bounds that need no solver.
- `scripts/threepoint_verify_indep.py`, an independent check of the eight
  three-point certificates. It imports neither the solver nor
  `threepoint_verify.py`: it rebuilds the programme from its definitions, takes
  only the dual multipliers from each file, proves the dual blocks positive
  definite exactly and bounds α by an exact rational. All eight bounds are
  confirmed, each within 4·10⁻⁶ of the solver's value
  (`certificates/threepoint_indep_checks.txt`). `tests/test_threepoint_indep.py`
  also checks it against actual independent sets, and shows that misreading a
  certificate destroys its bound.

- The Zenodo DOI of version 1.1.0, 10.5281/zenodo.22985036, in `CITATION.cff`
  and in the README's BibTeX entry.

### Changed
- **Layout.** The project moved from `research/hadwiger-nelson/` to the root of
  the repository, and the note to `papers/planes-4-chromatic/`, where each paper
  has its own folder. Commands now run from the root. Release 1.1.0, and
  version 5 of the note, which cites it, keep the old paths.
- **One front page.** The README merges the former front page and project page.
  Each result has a status (proved, formally verified, computer proof, known),
  the evidence behind it and where to read it; one table gives the command that
  checks each result.
- `python-flint` is now a requirement: the independent checker multiplies
  integer matrices with FLINT.

### Fixed
- **Citations**, after an audit against primary sources:
  - standard results are now credited where they are used: de Bruijn–Erdős,
    Hoffman (with Haemers for the ratio form), Weil, Delsarte, Stiemke, Croft,
    Falconer, Molloy, Alon–Krivelevich–Sudakov, value precedence, cube and
    conquer, Lean 4 and Mathlib;
  - the 103-vertex graph with edges at 1 and 2/√3 is Exoo's and Ismailescu's,
    not Ismailescu's alone;
  - Speyer's 2-adic colouring of the Moser ring is in Polymath16's thread 2,
    not thread 3;
  - χ_f(ℝ²) ≤ 4.36, derived in the research log as if new, is Hochberg and
    O'Donnell's (1993), and the log no longer says that the LP line had
    stalled;
  - the claim that the repulsion of a distance had not been measured before
    is withdrawn: it is an empirical counterpart of the probability `p_d` of
    Polymath16's probabilistic formulation;
  - the README's references gain journal data, and the works the repository
    cites.

  The research log lists every fix under "Citation audit (27 September)".

- **A wrong remark in the research log.** It said that the largest
  independent sets found in G₁₃ contain no point together with a whole circle;
  11 of the 15 known orbits of 36-point sets do. The log now says so under the
  remark, with the counts for G₁₁ and G₁₃ and the script that checks them
  (`scripts/experiments/largest_sets_whole_circles.py`). No result depended on
  the remark.

- **Statuses that said too little.** No 4-colouring of 𝔽₂₃² and 𝔽₃₁²,
  `χ(G₅) ≥ 4` and `χ(G₃) ≥ 3` follow from Hoffman's bound in interval
  arithmetic, and linear colourings of 𝔽₁₇² need six colours by an exhaustive
  count; the notes gave SAT for all of them. A conditional in
  `notes/local_colourings.md` §6 and `notes/rigidity.md` §11 now states the
  character-sum estimate it needs. `tests/test_finite_planes.py` checks the new
  statements.

## [1.1.0] - 2026-09-27

None of the results below has been refereed.

### Added
- **Formal proofs in Lean 4**, with Mathlib, of χ(ℚ(√2, √3)²) = 4 and of
  Fischer's theorem χ(ℚ(√3, √11)²) = 4 (`research/hadwiger-nelson/lean/`).
  Both depend only on Lean's three standard axioms. The workflow
  `.github/workflows/lean.yml` builds them, compares their axioms with
  `lean/axioms.expected` and replays them in Lean's kernel with `leanchecker`.
- Version 5 of the note, now in LaTeX (`docs/note/planes-4-chromatic.tex`):
  five pages, with a figure of the 10-vertex graph and a paragraph on the
  formal proofs, reviewed independently before release. Citations added or
  completed: L. and W. Moser (1961) for the spindle; Moorhouse's draft and its
  address; the Polymath16 threads, with their titles, comment numbers and
  addresses; issue and zbMATH numbers.
- `scripts/g17_alpha.py`, towards α(G₁₇) ≤ 57, which would give χ(G₁₇) ≥ 6:
  - the 57-point rosette found by kissat;
  - part A: no independent set of 58 points contains a point together with its
    whole circle, with the DRAT-checked log `certificates/g17_part_a_checks.txt`;
  - `tests/test_g17.py` (in CI) and `tests/test_g17_slow.py`.

  The remaining case, part B, is open (research log, "`α(G₁₇)`: the best sets
  are rosettes").
- The Zenodo DOIs: 10.5281/zenodo.22976636 for version 1.0.0 and
  10.5281/zenodo.22976635 for all versions. They are in `CITATION.cff` and in
  the README, as a badge and in the BibTeX entry.
- `CONTRIBUTING.md`: how to report an error or a result that does not
  reproduce, and what a pull request needs. A pull request template.

### Changed
- The note's HTML source, printed with Chromium, is replaced by the LaTeX
  source. In the note the prime is now written 𝔭 and the colouring κ, which
  were P and c, the names of points and of field elements.
- `hn.coloring`: a solve with no vertex subset, as in `is_k_colorable`, passes
  the vertex selectors as unit clauses instead of one assumption per vertex.
  CaDiCaL was 1.7 to 5.4 times faster on four 5-chromatic graphs of `data/` at
  k = 4. The core returned for such a solve is the whole vertex set, so
  `find_uncolorable_core` now names the subset for its first core.

## [1.0.0] - 2026-09-26

First versioned release. None of the results below has been refereed.

### Added
- **χ(ℚ(√2, √3)²) = 4**, Voronov's second case, and a short proof of
  K. G. Fischer's theorem χ(ℚ(√3, √11)²) = 4 (1994), in a three-page note
  (`research/hadwiger-nelson/docs/note/planes-4-chromatic.pdf`).
- Local proofs that χ(ℚ(√−3, √−11)) = 4 and χ(ℚ(√−3, √−11, √−247)) = 5 for
  the whole complex fields, by reduction at the primes 2 and 11. Both values
  also follow from earlier work: the first from Fischer's theorem, the second
  from Madore's reduction at 11 and Exoo–Ismailescu's graph.
- Necessary local conditions for a number field to hold a 6-chromatic
  unit-distance graph, and screens of multiquadratic fields.
- **Finite planes.** χ ≥ 6 for 𝔽₃₇², 𝔽₄₁², 𝔽₄₃², 𝔽₄₇² and the anisotropic
  planes G₂₉, G₃₇, G₄₁, from Schrijver's three-point bound with checked dual
  certificates. With Hoffman's bound, χ(G_q) ≥ 6 for every prime q ≥ 29 except
  31. Moorhouse's table of χ(𝔽_q²) is continued to q = 61, and χ(𝔽₁₃²) = 5.
- Reproductions, with DRAT proofs checked by `drat-trim`: the Moser spindle has
  no 3-colouring, and de Grey's 1581-vertex graph has no 4-colouring.
- A DRAT proof, checked by `drat-trim`, that the 803-vertex graph
  `five_247_c` over ℚ(√3, √11, √247) has no 4-colouring.
- Multi-distance graphs with no 5-colouring (187 and 72 points), for the
  Exoo–Ismailescu route to χ(ℝ²) ≥ 6.
- `scripts/check_no4.py`, which rechecks from the coordinates every graph in
  `data/` said to have no proper 4-colouring; its log
  (`certificates/data_no4_checks.txt`) has drat-trim verifying kissat's DRAT
  proof for all 27. drat-trim logs for the Moser spindle, the 19-vertex graph,
  the two pressure certificates, both multi-distance witnesses and the forced
  pair of Exoo–Ismailescu's graph H.
- `hn/`: exact number fields, unit-distance graphs, SAT colouring,
  certificates and local (adelic) colourings.
- 90 maintained scripts, among them `scripts/decompositions.gp` (PARI/GP) and
  `scripts/check_no4.py`, and 733 exploratory ones, kept as a record.
- Graphs and witnesses in exact coordinates; colourings, DRAT logs and
  three-point certificates with their checkers.
- 665 tests; GitHub Actions runs 383 of them on pushes to `main` and on pull
  requests.
- `requirements-lock.txt`: the exact environment in which the results were
  produced.

### Fixed
Corrections made before this release:
- K. G. Fischer had proved χ(ℚ(√3, √11)²) = 4 in 1994. We found his paper
  after the first version of the note had been sent to two mathematicians; the
  note and the project page now credit him.
- The two whole-field results for complex fields had been presented as new.
  They follow from earlier work, as stated above, and are now credited.
- The proof of Proposition B (χ(G_q) ≥ 6 for q ≥ 53) omitted q = 61, which
  Weil's bound does not reach. The statement holds (χ_f(G₆₁) ≥ 5.20); the
  proof and a test now cover it.
- The note said the colourings had been tested on graphs of up to about
  12 000 vertices, and that a second check had used PARI/GP; neither was
  recorded in the repository. The tested graphs have up to 3 134 vertices,
  and the PARI/GP check is now a script with a test.
- Solver time limits in `hn.coloring` did not stop CaDiCaL, the default
  solver, because pysat's interrupt does not reach it; a limited solve now runs
  it in a separate process. The slow tests' limits are four hours, and the
  documentation no longer promises a 30-minute solve of de Grey's graph.
- Descriptions of the data and certificates: `five_247_c` is not a subgraph of
  `five_247` (they share 317 points); `five_247_b` was built from 327 vertices
  of `Sa`, not 340; distances (9 ± √33)/6 in the 47-point pressure witness;
  404 shared points, not 402, in the seed of the L16 searches; the field of
  nine data files; the stale notes of de Grey's certificate; and the
  integrality step behind χ(G₁₃) ≥ 5.
- `python3 -m hn.cli verify` printed "VERIFIED" when drat-trim was missing and
  the proof had not been checked; it now prints the solver's UNSAT and exits
  with status 2. `demo` and `degrey` no longer write over the stored
  certificates.
- `scripts/six.py` shadowed the `six` package when `scripts/` came first on the
  import path, and two tests failed in a full run; it is now
  `scripts/degrey_forced_pair.py`.

[Unreleased]: https://github.com/decalion89/chromatic-number-of-the-plane/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/decalion89/chromatic-number-of-the-plane/releases/tag/v1.1.0
[1.0.0]: https://github.com/decalion89/chromatic-number-of-the-plane/releases/tag/v1.0.0
