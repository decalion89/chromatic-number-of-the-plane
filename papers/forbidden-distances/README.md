# Measurable density bounds hold for all colourings of the plane

[`forbidden-distances.pdf`](forbidden-distances.pdf) is a seven-page note. OpenAI's proof that χ(ℝ²) ≥ 6
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
  Conjecture 1 of Dúcz and Varga.
- **Two forbidden distances can force ten colours**: `χ(ℝ², {1, 2}) ≥ 7`, `χ(ℝ², {1, √3}) ≥ 7`,
  `χ(ℝ², {1, (√6+√2)/2}) ≥ 7`, `χ(ℝ², {1, 11/6}) ≥ 8`, `χ(ℝ², {1, 40}) ≥ 9` and `χ(ℝ², {1, 100}) ≥ 10`. The best
  bound known before for a pair was 6. Measurable versions of some of these were known (Oliveira Filho and
  Vallentin, 2010); the bounds for `{1, 2}` and `{1, √3}` use clique inequalities and seem new even there.
- **`k` forbidden distances can force more than `4.0485^k` colours**, in the plane and in space. So Erdős's
  `L₂(r)` (Erdős Problem #706) grows exponentially, `4.0485^r < L₂(r) ≤ 7^r`, whereas Erdős expected polynomial
  growth. Davies proved the analogue in dimension at least 4 two days earlier by a different method
  (arXiv:2610.12301, 8 October 2026).

All the results about arbitrary colourings rest on OpenAI's rigidity theorem (their Theorem 2.3). No finite graph
is given: as in OpenAI's proof, the graphs exist by compactness.

Two AI-assisted referees read the note.
- The first read the transfer theorems: correct, with minor errors and gaps, all fixed. Among them were a false
  numerical inequality in the headline (`4.0486 < 1/0.2470`), a clash with the notation of Matolcsi et al. (which
  led to proving the equality for the fractional chromatic number of the whole plane as well), the equality of
  densities over discs and over squares, and the proof of Lemma 1 at radius `d`.
- The second read the corollaries and the certificates, with two checkers of its own: no mathematical error, and
  every certificate valid. It found that the priority claims were wrong: the draft did not know Dúcz–Varga or
  Davies. The novelty claims are corrected, Erdős's question is quoted from the original, the measurable
  precedents of Oliveira Filho and Vallentin are credited, and its other findings are applied.

No mathematician has checked the note yet.

[`forbidden-distances.tex`](forbidden-distances.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run in this
folder, with TeX Live:

```sh
pdflatex forbidden-distances.tex && pdflatex forbidden-distances.tex
```

The files behind it:
- [`notes/forbidden_distances.md`](../../notes/forbidden_distances.md), the record: statements, what is new, what
  each referee found;
- [`data/forbidden_distances/`](../../data/forbidden_distances/README.md): the eight certificates, the programs that
  produce and check them (three independent checks), and `tests/test_forbidden_distances.py`;
- the research log, entry of 10 October.
