# The ruler-and-compass plane is not 5-colourable

[`ruler-compass.pdf`](ruler-compass.pdf) is a six-page note. It shows that OpenAI's proof that χ(ℝ²) ≥ 6
([*The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf),
September 2026) works over the constructible numbers, the smallest subfield of ℂ closed under square roots. So the
points constructible with ruler and compass cannot be 5-coloured. By the de Bruijn–Erdős theorem, some finite
unit-distance graph with constructible vertices has no proper 5-colouring.

The proof is OpenAI's. The note gives the two changes that make square roots enough:
- a lemma about elements of order 2 in finite images (Lemma 3);
- the triple average of their Lemma 3.3, taken at the exponents −1, 0, 1 (Lemma 4).

It also describes the formal proof, which was checked in Lean 4, together with the finite form (Corollary 2).

Section 6 proves, on paper, that the plane over every Euclidean field is not 5-colourable either (Corollary 5),
because the real constructible numbers embed in every Euclidean field. For fields that contain a copy of the real
algebraic numbers, such as the real closed fields, this already follows from OpenAI's proof; for the others, such as
the real constructible numbers themselves, it rests on the ruler-and-compass theorem. A remark in Section 7 shows that
over the Pythagorean closure of ℚ (Hilbert's field) the analogue of OpenAI's rigidity theorem is false, so their first
step cannot be run there as it stands; whether that plane is 5-colourable stays open.

Three AI-assisted reviews read the paper.
- The first found no mathematical error in the proof of Theorem 1, and one wrong side claim: that any shift of the
  exponents needs cube roots. That holds only for the single shift of the formalization, and the paper now says so.
  Its other findings, on wording and on how faithfully the paper reports what was checked, are applied.
- The second and third read Section 6 and the remark on the Pythagorean closure. They found no mathematical error.
  One overclaim is corrected: the draft said that the corollary needs Theorem 1 for every Euclidean field that is not
  real closed. The other findings are applied.

No mathematician has checked the paper yet.

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
