# Colouring the p-adic plane

[`padic-planes.pdf`](padic-planes.pdf) is a six-page draft about the graph on ℚ_p² in which two points are adjacent
when (x − x′)² + (y − y′)² = 1:

- For p ≡ 3 (mod 4) every measurable proper colouring needs at least 1 + (p + 1)/(2√p) colours, so the Borel
  chromatic number of ℚ_p² is unbounded. This answers Question 1 of Bardestani and Mallahi-Karai
  ([arXiv 1507.05300](https://arxiv.org/abs/1507.05300)), the "p-adic Hadwiger–Nelson problem", in the negative.
- χ(ℚ₇²) = 4, for ordinary and for Borel colourings; χ(ℚ₂²) = 2 and χ(ℚ₃²) = 3 follow from Madore's results
  (arXiv 1509.07023). The three ordinary values are also proved in Lean 4
  ([`lean/PadicPlanes.lean`](../../lean/PadicPlanes.lean)).
- With a theorem of Davies ([arXiv 2308.16885](https://arxiv.org/abs/2308.16885)), the graph of every isotropic form
  over a field of characteristic 0 has infinite chromatic number, so the dichotomy of Bardestani and Mallahi-Karai
  holds for ordinary colourings too; in particular χ(ℚ_p²) = ∞ for p ≡ 1 (mod 4), and χ(ℚ_p²) is finite exactly
  for p = 2 and p ≡ 3 (mod 4).
- χ(ℚ_p²) ≥ 4 for every prime p ≡ 3 (mod 4) with 7 ≤ p < 2 129 503 819, from the unit-distance graphs over 27 real
  quadratic fields of [`papers/quadratic-planes/`](../quadratic-planes/); by Chebotarev's theorem no finite set of
  number fields gives this for every p.

[`padic-planes.tex`](padic-planes.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run in this folder, with
TeX Live:

```sh
pdflatex padic-planes.tex && pdflatex padic-planes.tex
```

The computations are in [`data/quadratic_planes/scripts/`](../../data/quadratic_planes/scripts/):
`padic_planes.py` (the colourings of 𝔽₃², 𝔽₇², 𝔽₁₁², the 5-cycle over ℚ(√7), the table up to p = 83),
`padic_measurable.py` (Table 1, the least eigenvalues of 𝔽_p² in interval arithmetic) and `padic_reach.c` with
`padic_reach.py` (the first prime the 27 fields miss). The fuller account is
[`notes/quadratic_planes.md`](../../notes/quadratic_planes.md), §5; the tests are in
`tests/test_quadratic_planes.py`. The draft has not been refereed.
