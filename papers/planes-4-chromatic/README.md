# The planes over ℚ(√3, √11) and ℚ(√2, √3) are 4-chromatic

[`planes-4-chromatic.pdf`](planes-4-chromatic.pdf) is a five-page note proving χ(ℚ(√2, √3)²) = 4
and giving a short proof of K. G. Fischer's theorem χ(ℚ(√3, √11)²) = 4 (1994). Both theorems are also
proved in Lean 4 ([`lean/`](../../lean/)).
[`planes-4-chromatic.tex`](planes-4-chromatic.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run
in this folder, with TeX Live:

```sh
pdflatex planes-4-chromatic.tex && pdflatex planes-4-chromatic.tex && pdflatex planes-4-chromatic.tex
```

Figure 1 is drawn with TikZ from decimal coordinates; its exact coordinates are those of
[`data/chain23.json`](../../data/chain23.json), and `lean/Q23.lean` proves that its 16 edges have length 1.
Versions 1 to 4 were written in HTML and printed to PDF with headless Chromium; version 5 is the first in
LaTeX, and adds the formal proofs. Version 6 (28 September 2026) credits the repository
[hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction) (July 2026), found after
version 5: its Corollary B′ already gives the upper bounds, and its Theorem A is the same reduction at 2.

The fuller account, with the checks behind each step, is
[`notes/local_colourings.md`](../../notes/local_colourings.md) §8 and §10; the tests are
`tests/test_q311.py` and `tests/test_q23.py`. The note has not been refereed.
