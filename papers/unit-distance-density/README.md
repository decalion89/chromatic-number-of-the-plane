# A new upper bound for the density of planar sets avoiding unit distances

[`unit-distance-density.pdf`](unit-distance-density.pdf) is a five-page note. It proves

> every measurable planar set with no two points at distance 1 has upper density at most **0.24415**,

that is, `m₁(ℝ²) ≤ 0.24415`. The best bound before was `0.2470`, by Ambrus, Csiszárik, Matolcsi, Varga and
Zsámboki (*The density of planar sets avoiding unit distances*, Math. Program. 207 (2024) 303–327,
arXiv:2207.14179), who settled Erdős's conjecture `m₁(ℝ²) < 1/4`. The method is theirs: a linear program on the
patterns in which a random congruent copy of a finite set `X` meets the set. Their `X₂₃` has 23 points, found by a
beam search over sets of at most 24 points; the note adds 6 points with coordinates in the same field
`ℚ(√3, √11)`, found one at a time among the points at distance 1 from two points already present.

With the transfer theorem of OpenAI, as extended to fractional colourings in
[`papers/forbidden-distances/`](../forbidden-distances/README.md), the fractional chromatic number of the plane is
at least `1/m₁(ℝ²) ≥ 4.0959` (it was `4.0485`), and there are unit-distance graphs on `n` points with independence
number below `(0.24415 + ε)n` for every `ε > 0` and large `n`.

The proof is a rational certificate. Two programs written separately check it in exact arithmetic over
`ℚ(√3, √11)` and in interval arithmetic for the Bessel function `J₀`; both also accept the certificate for `X₂₃`
alone (`0.24698967`, the bound of Ambrus et al.) and reject corrupted copies. In floating point the value of the
program keeps decreasing as points are added (Table 2), so the method has not reached its limit.

The bound is weaker than the `0.2415` that Varga announced for the same group in December 2023 (blog post of the
Rényi Institute's AI group, *Solving a conjecture of Erdős*); we found no paper, preprint or certificate for that bound.
As far as we know, `0.24415` is the best upper bound with a published proof. The method is entirely that of Ambrus et
al.; what is new is the configuration and the certificate.

An internal referee (an AI agent working separately) proved the inclusion–exclusion bound (Proposition 3) in full,
wrote a third checker before reading ours, which accepts the three certificates, and ran 28 corrupted certificates
through the three checkers: every invalid one was rejected. Its findings are applied: `check_m1.py` no longer relies
on `assert` statements (which Python skips under `-O`), it cross-checks the number of independent sets and it is
self-contained; the announced `0.2415` is cited; the facts about `J₀` used by the checkers are given with proofs or
references. No mathematician has checked the note yet.

[`unit-distance-density.tex`](unit-distance-density.tex) is its LaTeX source (`amsart`), with the figure in
[`fig-config.tex`](fig-config.tex) (TikZ, written by `fig_config.py` from the configuration). To rebuild the PDF, run
in this folder, with TeX Live:

```sh
pdflatex unit-distance-density.tex && pdflatex unit-distance-density.tex
```

The files behind it are in [`data/forbidden_distances/m1/`](../../data/forbidden_distances/m1/README.md): the
configurations, the certificates, both checkers, the certificate generator and the search programs.
