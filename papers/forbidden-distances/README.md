# Measurable density bounds hold for all colourings of the plane with algebraic forbidden distances

[`forbidden-distances.pdf`](forbidden-distances.pdf) is a fifteen-page note (version 7). OpenAI's proof that χ(ℝ²) ≥ 6
([*The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf),
September 2026) first transfers an arbitrary proper colouring of the plane to a weak measurable one. The note
observes that this transfer works, with the same proof, for every set `D` of algebraic forbidden distances and for
fractional colourings. So `χ(ℝ², D) ≥ 1/m_D(ℝ²)`, where `m_D(ℝ²)` is the largest density of a measurable planar
set with no two points at a distance in `D`, and for finite `D` the fractional chromatic number equals `1/m_D(ℝ²)`.
Known and new upper bounds for `m_D(ℝ²)` then give:

- **The fractional chromatic number of the plane is `1/m₁(ℝ²)`**, between 4.0485 and 4.36, both for finite
  unit-distance graphs and for the whole plane. Dúcz and Varga had already shown that it exceeds 4.0007
  (arXiv:2606.28157, June 2026), disproving a conjecture of Matolcsi, Ruzsa, Varga and Zsámboki; this is a second
  disproof, and it answers their question whether the two fractional chromatic numbers agree.
- **The least independence ratio of an `n`-point unit-distance graph tends to `m₁(ℝ²)`**: the first half of
  Conjecture 1 of Dúcz and Varga (that it is at least `m₁(ℝ²)` is classical).
- **Two forbidden distances can force twelve colours, with certified bounds**: `χ(ℝ², {1, 2}) ≥ 7`,
  `χ(ℝ², {1, √3}) ≥ 7`, `χ(ℝ², {1, (√6+√2)/2}) ≥ 7`, `χ(ℝ², {1, 11/6}) ≥ 8`, `χ(ℝ², {1, 40}) ≥ 9`,
  `χ(ℝ², {1, 100}) ≥ 10` and `χ(ℝ², {1, 1000}) ≥ 12`. The best bound known before for a pair was 6. Measurable
  versions of some of these were known (Oliveira Filho and Vallentin, 2010); the bounds for `{1, 2}` and `{1, √3}`
  use clique inequalities and seem new even there.
- **Three or four forbidden distances can force 48 colours**: `χ(ℝ², {1, 40, 1600}) ≥ 22`,
  `χ(ℝ², {1, 10, 100, 1000}) ≥ 26`, `χ(ℝ², {1, 100, 10⁴}) ≥ 27` and `χ(ℝ², {1, 30, 900, 27000}) ≥ 48`, from
  certified linear programs.
- **`k` forbidden distances can force more than `4.0485^k` colours**, in the plane and in every higher dimension. So
  Erdős's `L₂(r)` (Erdős Problem #706) grows exponentially, whereas Erdős expected polynomial growth. Davies proved
  the analogue in dimension at least 4 two days earlier by a different method (arXiv:2610.12301, 8 October 2026).
- **Two scales (version 5).** Using OpenAI's theorem a second time, as a statement about finite graphs, at the scale
  of `D₁`, and densities at the scale of `N·D₂`: for all large algebraic `N`, `χ(ℝ², D₁ ∪ N·D₂)` is at least the
  `c`-fold chromatic number of the plane with distances `D₂`, where `c = χ(ℝ², D₁)` — the lower half of a geometric
  version of Geller and Stahl's formula `χ(G[H]) = χ_{χ(H)}(G)` for lexicographic products. Hence
  **`χ(ℝ², {1, N}) ≥ 25` for every large algebraic `N`, although `χ_f(ℝ², {1, N}) < 19.01` for every `N`**: two
  distances can force about 1.3 times the fractional chromatic number. Iterating, three and four distances force 102
  and 413 colours, and `L₂(r) > 1.48·4.0486^r`. With explicit graphs the proof is effective: de Grey's graph and the
  Moser spindle give `χ(ℝ², {1, N}) ≥ 18` for every algebraic `N ≥ 10^25` (and 20 for `N ≥ 10^34` with the graph of
  Bellitto, Pêcher and Sédillot).
- **Uniform bounds for two distances (versions 6 and 7).** One certificate can serve a whole range of `d`: for every
  real algebraic `d`, `χ(ℝ², {1, d}) ≥ 7` for `d ≥ 5.82`, and the bound rises to 8, 9, 10, 11, 12, 13 and 14 from
  `d = 10.5`, `18.9`, `35.7`, `73.6`, `181`, `655` and `25 700`; also 7 colours for
  `1.72 ≤ d ≤ 2.04` (which contains `√3`, `(√6+√2)/2` and `2`) and 8 for `1.82 ≤ d ≤ 1.89` (Theorem 5). The same
  holds for `1/d`, and for measurable colourings for every real `d`. Each step is one function that works for every
  `d ≥ D` at once (Lemma 18: `J₀(x) ≥ −min(0.4028, √(2/(πx₀)))` for `x ≥ x₀`), with the configuration inequalities of
  unit equilateral triangles and of Moser spindles in the linear program; version 6, without them, had the thresholds
  6.52, 12.05, 24.3, 52.9, 140.2, 601 and 54 300 for 7 to 13 colours and could not give 14: without configurations
  the program's value tends to `1/12.13…`, the square of its one-distance value. With triangles and spindles it tends
  to about `1/13.27` in floating point, so these configurations do not give 15 colours.

  ![Certified lower bounds for χ(ℝ², {1, d})](uniform-bounds.png)

All the lower bounds rest on OpenAI's rigidity theorem (their Theorem 2.3), which OpenAI formalized in Lean. The
transfer of Theorem A (for proper colourings, not the fractional version) is formalized in Lean on top of their
library: [`lean/external/openai-five/distances/`](../../lean/external/openai-five/distances/README.md), standard
axioms only, replayed by `leanchecker`. No finite graph is given: as in OpenAI's proof, the graphs exist by
compactness.

Six AI-assisted referees read the note.
- The first read the transfer theorems: correct, with minor errors and gaps, all fixed. Among them were a false
  numerical inequality in the headline (`4.0486 < 1/0.2470`), a clash with the notation of Matolcsi et al. (which
  led to proving the equality for the fractional chromatic number of the whole plane as well), the equality of
  densities over discs and over squares, and the proof of Lemma 1 at radius `d`.
- The second read the corollaries and the certificates, with two checkers of its own: no mathematical error, and
  every certificate valid. It found that the priority claims were wrong: the draft did not know Dúcz–Varga or
  Davies. The novelty claims are corrected, Erdős's question is quoted from the original, the measurable
  precedents of Oliveira Filho and Vallentin are credited, and its other findings are applied.
- The third read version 3, with a checker of its own: no mathematical error. It found that Croft's constant was
  rounded up (`m₁ ≥ 0.229364`, not `0.229365`), that `f(n) ≥ m₁n` is classical and was not credited, and that the
  title said more than the theorem (now "with algebraic forbidden distances"); these and its minor findings are
  applied in version 4.
- The fourth read version 5 (the two-scale theorem, its strong form, the explicit pairs, the Lean paragraph and the
  48-colour certificate), with numerical checks and a certificate checker of its own: no mathematical error, nothing
  in the literature that already contains the theorem or any bound above 6 for a pair. It asked for wording and
  attribution fixes, all applied: the Lean paragraph names the declarations and says what is not formalized; the
  first 6-colour pair, `{1, (1+√5)/2}`, is credited to Owings, Tetiva and Huddleston; Geller–Stahl is cited; and
  the headline 25 is noted to need only a finite unit-distance graph with fractional chromatic number above 4.
- The fifth read version 6 (Theorem 5, Lemma 18, the certificates, the programs and the figure), with a checker of
  its own that shares no code or Bessel routine with ours (integer interval arithmetic, its own series for `J₀` and
  `J₁`, Taylor models): no mathematical error, all ten certificates valid (the smallest margin is `2.5·10⁻⁷`, for 13
  colours), and no earlier bound of 7 for a pair. It asked for the measurable statement to be argued, for one half of
  the limit in the remark to be proved (its construction is in the paper), and for an overclaim in the abstract
  ("each step by a single certificate", false for the step at 12.05) to be corrected; all applied.
- The sixth read version 7 (unit triangles and Moser spindles in the uniform certificates), with a checker of its own
  (no code or Bessel routine shared with ours): no mathematical error, all ten certificates valid (the smallest
  margin, `min h ≥ 5.9·10⁻⁷`, is for 14 colours), the eight tails rebuilt byte for byte. It corrected the value of a
  certificate in the remark (`z₀ + Σ α_c y_c`), asked for the conditions of the limiting construction to be stated,
  and found two overclaims (14 colours "follow" from a floating-point value; 13.27 called a limit of the program) and
  stale statements in Corollary 4 and the READMEs; all applied.

No mathematician has checked the note yet.

[`forbidden-distances.tex`](forbidden-distances.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run in this
folder, with TeX Live:

```sh
pdflatex forbidden-distances.tex && pdflatex forbidden-distances.tex
```

The files behind it:
- [`notes/forbidden_distances.md`](../../notes/forbidden_distances.md), the record: statements, what is new, what
  each referee found;
- [`data/forbidden_distances/`](../../data/forbidden_distances/README.md): the thirteen certificates, the programs
  that produce and check them (three independent checks), and `tests/test_forbidden_distances.py`;
- [`data/forbidden_distances/uniform/`](../../data/forbidden_distances/uniform/README.md): the ten certificates of
  Theorem 5, the two programs that check them, and the sources of the figure;
- [`lean/external/openai-five/distances/`](../../lean/external/openai-five/distances/README.md): the Lean proof of
  the transfer;
- the research log, entries of 10 October.
