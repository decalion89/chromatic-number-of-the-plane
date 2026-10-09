# The ruler-and-compass plane is not 5-colourable

[`ruler-compass.pdf`](ruler-compass.pdf) is a five-page note. It shows that OpenAI's proof that χ(ℝ²) ≥ 6
([*The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf),
September 2026) works over the constructible numbers, the smallest subfield of ℂ closed under square roots. So the
points constructible with ruler and compass cannot be 5-coloured. By the de Bruijn–Erdős theorem, some finite
unit-distance graph with constructible vertices has no proper 5-colouring.

The proof is OpenAI's. The note gives the two changes that make square roots enough:
- a lemma about elements of order 2 in finite images (Lemma 3);
- the triple average of their Lemma 3.3, taken at the exponents −1, 0, 1 (Lemma 4).

It also describes the formal proof, which was checked in Lean 4, together with the finite form (Corollary 2).

An AI-assisted review of the paper found no mathematical error in the proof of Theorem 1. It found one wrong side
claim: that any shift of the exponents needs cube roots. That holds only for the single shift of the formalization;
the paper now says so. Its other findings, on wording and on how faithfully the paper reports what was checked, are
applied. No mathematician has checked the paper yet.

[`ruler-compass.tex`](ruler-compass.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run in this folder,
with TeX Live:

```sh
pdflatex ruler-compass.tex && pdflatex ruler-compass.tex && pdflatex ruler-compass.tex
```

The files behind it:
- [`notes/six_over_fields.md`](../../notes/six_over_fields.md), a fuller account of the same changes;
- [`lean/external/openai-five/fields/`](../../lean/external/openai-five/fields/README.md): the patches against
  openai/math, the Lean statements, the axioms printed, and a build script;
- the research log, entry of 9 October.
