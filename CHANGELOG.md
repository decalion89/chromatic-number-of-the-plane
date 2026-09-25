# Changelog

Notable changes to this repository, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers follow
[Semantic Versioning](https://semver.org/) for the library, the tools and the
data formats.

Mathematical corrections are listed in each version. The full record,
including every retracted claim, is the research log,
[`research/hadwiger-nelson/docs/research-log.md`](research/hadwiger-nelson/docs/research-log.md).

## [1.0.0] – 2026-09-26

First versioned release. None of the results below has been refereed.

### Results
- **χ(ℚ(√2, √3)²) = 4**, Voronov's second case, and a short proof of
  K. G. Fischer's theorem χ(ℚ(√3, √11)²) = 4 (1994), in a three-page note
  (`research/hadwiger-nelson/docs/note/planes-4-chromatic.pdf`).
- χ(ℚ(√−3, √−11)) = 4 and χ(ℚ(√−3, √−11, √−247)) = 5 for the whole complex
  fields, by reduction at the primes 2 and 11.
- Necessary local conditions for a number field to hold a 6-chromatic
  unit-distance graph, and screens of multiquadratic fields.
- **Finite planes.** χ ≥ 6 for 𝔽₃₇², 𝔽₄₁², 𝔽₄₃², 𝔽₄₇² and the anisotropic
  planes G₂₉, G₃₇, G₄₁, from Schrijver's three-point bound with checked dual
  certificates. With Hoffman's bound, χ(G_q) ≥ 6 for every prime q ≥ 29 except
  31. Moorhouse's table of χ(𝔽_q²) is continued to q = 61, and χ(𝔽₁₃²) = 5.
- Reproductions, with DRAT proofs checked by `drat-trim`: the Moser spindle has
  no 3-colouring, and de Grey's 1581-vertex graph has no 4-colouring.
- Multi-distance graphs with no 5-colouring (187 and 72 points), for the
  Exoo–Ismailescu route to χ(ℝ²) ≥ 6.

### Contents
- `hn/`: exact number fields, unit-distance graphs, SAT colouring,
  certificates and local (adelic) colourings.
- 88 maintained scripts and 733 exploratory ones, kept as a record.
- Graphs and witnesses in exact coordinates; colourings, DRAT logs and
  three-point certificates with their checkers.
- 653 tests, 371 of them run by GitHub Actions on every push and pull request.
- `requirements-lock.txt`: the exact environment in which the results were
  produced.

### Corrections made before this release
- K. G. Fischer had proved χ(ℚ(√3, √11)²) = 4 in 1994. We found his paper
  after the first version of the note had been sent to two mathematicians; the
  note and the project page now credit him.
- The proof of Proposition B (χ(G_q) ≥ 6 for q ≥ 53) omitted q = 61, which
  Weil's bound does not reach. The statement holds (χ_f(G₆₁) ≥ 5.20); the
  proof and a test now cover it.

[1.0.0]: https://github.com/decalion89/darwin-50/releases/tag/v1.0.0
