# A new upper bound for the density of sets in three-space avoiding unit distances

[`unit-distance-density-3d.pdf`](unit-distance-density-3d.pdf) is a four-page note. It proves

> every measurable set in `ℝ³` with no two points at distance 1 has upper density at most **0.15158538**,

that is, `m₁(ℝ³) ≤ 0.15158538`. The best published bound was `0.1532996`, by DeCorte, Oliveira Filho and Vallentin
(*Complete positivity and distance-avoiding sets*, Math. Program. 191 (2022) 487–558, arXiv:1804.09099, Table 2);
before it, `0.165609` (Oliveira Filho and Vallentin, J. Eur. Math. Soc. 12 (2010)) and `3/16` (Székely and Wormald,
1989).

The method is the inclusion–exclusion linear program of Ambrus, Csiszárik, Matolcsi, Varga and Zsámboki (Math. Program.
207 (2024) 303–327, arXiv:2207.14179), who used it only in the plane: in `ℝ³` the Bessel function `J₀` becomes
`Ω(x) = sin x / x`, and we add the atom `δ²` that the spectral measure of a set of density `δ` has at the origin (as
Keleti, Matolcsi, Oliveira Filho and Ruzsa did in the plane). The finite set has 27 points `q/√2` with `q ∈ ℚ³`, grown
from a cuboctahedron and its centre, so all squared distances are rational.

Consequence: the measurable fractional chromatic number of `ℝ³` is at least `6.5969` (it was `6.5232`). The measurable
chromatic number bound `χ_m(ℝ³) ≥ 7` does not change (it would need `m₁(ℝ³) < 1/7`), and unlike the plane we know no
transfer theorem that would give a bound for all colourings.

The proof is a rational certificate. Two programs written separately check it, in exact rational arithmetic and in
interval arithmetic for `Ω`, and reject corrupted copies. An internal referee (an AI agent working separately) checked
the transfer of the argument to `ℝ³` against the sources and found no error; its two findings (continuity of the
averaged atom densities, and the attribution of the previous bounds) are in the text. No mathematician has checked the
note yet.

[`unit-distance-density-3d.tex`](unit-distance-density-3d.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run
in this folder, with TeX Live:

```sh
pdflatex unit-distance-density-3d.tex && pdflatex unit-distance-density-3d.tex
```

The files behind it are in [`data/forbidden_distances/m1r3/`](../../data/forbidden_distances/m1r3/README.md): the
configuration, the certificate, both checkers, the certificate generator and the search program.
