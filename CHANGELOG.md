# Changelog

Notable changes to this repository, newest first. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Version numbers follow
[Semantic Versioning](https://semver.org/) for the library, the tools and the
data formats.

Mathematical corrections are listed in each version. The full record,
including every retracted claim, is the research log,
[`docs/research-log.md`](docs/research-log.md).

## [Unreleased]

### Added
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
