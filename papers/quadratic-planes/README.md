# Real quadratic planes that need four colours

[`quadratic-planes.pdf`](quadratic-planes.pdf) is a seven-page draft proving χ(ℚ(√d)²) = 4 for twenty-four real
quadratic fields and 4 ≤ χ(ℚ(√47)²) ≤ 5. The lower bounds are computer proofs: unit-distance graphs over
ℚ(√d) with no 3-colouring, certified by DRAT proofs checked by drat-trim for two encodings. The upper bounds
are known reductions modulo a prime.

[`quadratic-planes.tex`](quadratic-planes.tex) is its LaTeX source (`amsart`). Table 1 is written from the
data by [`make_table.py`](make_table.py), and Figure 1 by [`make_figure.py`](make_figure.py). To rebuild the PDF, run in this folder, with TeX Live:

```sh
python3 make_table.py && python3 make_figure.py && pdflatex quadratic-planes.tex && pdflatex quadratic-planes.tex
```

The fuller account is [`notes/quadratic_planes.md`](../../notes/quadratic_planes.md); the data and the search
code are in [`data/quadratic_planes/`](../../data/quadratic_planes/), the checker is
`scripts/verify_quadratic_planes.py` and the tests are `tests/test_quadratic_planes.py`. The draft has not been
refereed.
