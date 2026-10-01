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
- **Credit for earlier work (28 September).** The upper bound of `χ(ℚ(√2, √3)²) = 4`
  and of Fischer's `χ(ℚ(√3, √11)²) = 4` is also a case of Corollary B′ of
  [hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction), a public,
  unrefereed repository of July 2026, whose Theorem A is the same reduction of
  `z = x + iy` at a place over 2. We found it after our note. The front page,
  `notes/literature.md` and the note (version 6) now say so; what remains ours is
  the explicit case `ℚ(√2, √3)` with its lower bound, and the Lean proofs.

### Added

- **Real quadratic planes that need four colours** (`notes/quadratic_planes.md`,
  `data/quadratic_planes/`): `χ(ℚ(√d)²) = 4` for `d = 11, 23, 35, 59, 71, 119, 131,
  191, 239, 359`, and `4 ≤ χ(ℚ(√47)²) ≤ 5`; so `χ(ℚ(√d)²)` is now known for every squarefree
  `d < 83` except 47. For real quadratic fields the values known before were
  2 and 3, and a field can need four colours only if `d ≡ 11 (mod 12)`. Each
  lower bound is a triangle-free, vertex-critical unit-distance graph (106, 686,
  872, 462, 674, 469, 143, 394, 650, 511 and 741 vertices) with no 3-colouring: kissat with DRAT proofs checked
  by drat-trim, twice, with separate encodings. The upper bounds are known:
  Moorhouse's reduction at 7 (`d = 11, 23, 35, 71, 119, 191, 239, 359`), Fischer's Theorem 10
  (`d = 59, 131`) and the reduction at 11 (`d = 47`). The note also observes that
  Cohen's conjecture (2007) on the sets `ℚ(√−d) ⊂ ℂ` is a different question,
  which reduction at a ramified prime settles. `scripts/verify_quadratic_planes.py`
  checks everything again, and `tests/test_quadratic_planes.py` runs the fast
  checks.
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
