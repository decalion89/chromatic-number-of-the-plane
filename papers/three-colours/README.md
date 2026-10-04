# Three colours for the plane over a number field

[`three-colours.pdf`](three-colours.pdf) is an eleven-page draft (4 October 2026) about the graph on `F²`, for a
number field `F`, in which two points are adjacent when `(x − x′)² + (y − y′)² = 1`:

- **Theorem B.** `χ(F²) ≤ 3` if and only if some prime of `F` above 2 ramifies in `F(i)` or some prime of `F` above 3
  has residue degree 1 (and `χ(F²) ≤ 2` exactly in the first case). The "if" part is classical (colourings through a
  residue field); the "only if" part is new: it extends to every number field the theorem that
  `χ(ℚ(√d)²) ≥ 4` for `d ≡ 11 (mod 12)` ([`papers/four-colours/`](../four-colours/)). So the local–global question
  is answered at three colours for every number field.
- **Consequences.** No plane over a number field has circular chromatic number strictly between 2 and 3
  (Corollary 1); every field in which `i` lies in all completions above 2 and 3 needs four colours (Corollary 2); for
  a finite extension `K` of `ℚ_p`, `χ(K²)` is 2, 3 or at least 4 according to explicit local conditions
  (Corollary 3; for example `χ(ℚ₂₇²) ≥ 4`); every real field containing `√a` and `√b` with `a ≡ 2 (mod 3)`,
  `b ≡ 7 (mod 8)` needs four colours, by two explicit relations (Proposition 2, also proved in Lean:
  `lean/TwoRoots.lean`). New fields that need four colours include `ℚ(√2, √7)`, `ℚ(√2, √31)`, `ℚ(√2, √55)` and
  `ℚ(2cos(2π/7), √7)`.
- **Theorem C.** `χ_c(ℚ(√11)²) = χ_c(ℚ(√35)²) = χ_c(ℚ₇²) = 7/2`, while the chromatic numbers are 4: as far as we
  know the first non-integral value of the circular chromatic number of such a plane. The lower bounds are computer
  proofs (exact certificates, accepted by two checkers that share no code, with Theorem W⁺ of
  [`papers/winding/`](../winding/)); the upper bound is a colouring through the residue field `𝔽₇`. Locally constant
  colourings of one completion give no upper bound on `χ_c` below 4 other than 2, 3 and `7/2` (Propositions 4 and 5,
  Corollary 5; Weil's bound for Kloosterman sums and a computation for `p < 3000`). Question 1 asks whether these are
  the only values below 4 for every number field.

The working notes are [`notes/three_colours_number_fields.md`](../../notes/three_colours_number_fields.md) and
[`notes/circular_planes.md`](../../notes/circular_planes.md); the programs and certificates are in
[`data/number_fields/`](../../data/number_fields/), and the tests in `tests/test_three_colours.py`. Theorem B was
refereed twice and the paper once (with a second look at the new Section 8) by separate AI agents with their own
programs; no mathematician has checked it yet.

[`three-colours.tex`](three-colours.tex) is the LaTeX source (`amsart`). To rebuild the PDF, run in this folder,
with TeX Live:

```sh
pdflatex three-colours.tex && pdflatex three-colours.tex
```
