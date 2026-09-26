# Changelog

Notable changes to this repository, newest first. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers follow
[Semantic Versioning](https://semver.org/) for the library, the tools and the
data formats.

Mathematical corrections are listed in each version. The full record,
including every retracted claim, is the research log,
[`research/hadwiger-nelson/docs/research-log.md`](research/hadwiger-nelson/docs/research-log.md).

## [Unreleased]

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
- `hn/`: exact number fields, unit-distance graphs, SAT colouring,
  certificates and local (adelic) colourings.
- 89 maintained scripts, among them `scripts/decompositions.gp` (PARI/GP), and
  733 exploratory ones, kept as a record.
- Graphs and witnesses in exact coordinates; colourings, DRAT logs and
  three-point certificates with their checkers.
- 664 tests; GitHub Actions runs 382 of them on pushes to `main` and on pull
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

[Unreleased]: https://github.com/decalion89/chromatic-number-of-the-plane/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/decalion89/chromatic-number-of-the-plane/releases/tag/v1.0.0
