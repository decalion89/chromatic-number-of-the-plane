# Two-colourable planes over number fields

[`two-colour-planes.pdf`](two-colour-planes.pdf) is a five-page draft about the graph on `F²`, for a number field
`F`, in which two points are adjacent when `(x − x′)² + (y − y′)² = 1`:

- **Theorem 1.** `χ(F²) = 2` if and only if some prime of `F` above 2 ramifies in `F(i)`. The "if" direction, a
  2-colouring by reduction modulo that prime, is Theorem A′ of the public repository
  [hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction) (July 2026); Madore had the
  cases of `ℚ`, of fields with a prime above 2 unramified over `ℚ`, and of `ℚ(√2)`. The converse, Proposition 4, is
  the new part as far as we know: it finds an odd cycle whenever no prime above 2 ramifies. The unit vectors span a
  ring, a 2-colouring is a ring homomorphism to `𝔽₂` (Fischer 1990, Theorem 1(iii)), Chevalley's theorem puts it at
  a place above 2, and density of the circle in its completions shows that this place ramifies.
- It contains the known cases: real quadratic fields (Johnson; Fischer, Theorem 8), fields of odd degree
  (Moorhouse, Theorem 7.1) and multiquadratic fields (Corollary 5: `F²` is 2-colourable exactly when every
  quadratic subfield's plane is; the 2-colourable half is Corollary B of hn-2adic-obstruction). It decides fields
  that are not multiquadratic, such as quartic fields with Galois group `S₄` (Example 7). The same criterion holds
  for every nondegenerate binary form `x² + c y²`, with `F(√−c)` in place of `F(i)` (Remark 6).
- **A local–global question** (Question 8): is `χ(F²)` the least number of colours of a colouring of one
  completion `F_v²` that is locally constant (Borel at a real place)? Theorem 1 answers it when one side is 2. At
  three colours it predicts `χ(ℚ(√d)²) ≥ 4` exactly for `d ≡ 11 (mod 12)`, which 27 fields confirm. At four colours
  a positive answer for `ℚ(√167)` would give a triangle-free 5-chromatic unit-distance graph in the plane, which we
  have not found in the literature.

The working note behind it is [`notes/local_global.md`](../../notes/local_global.md). [`two-colour-planes.tex`](two-colour-planes.tex)
is the LaTeX source (`amsart`). To rebuild the PDF, run in this folder, with TeX Live:

```sh
pdflatex two-colour-planes.tex && pdflatex two-colour-planes.tex
```

The computations are in [`data/quadratic_planes/scripts/`](../../data/quadratic_planes/scripts/), with PARI/GP:

| script | what it checks |
|---|---|
| `two_colour_criterion.gp` | the criterion of Theorem 1 against Fischer's theorem (242 quadratic fields), Moorhouse's Theorem 7.1 (nine fields of odd degree) and Corollary 5 (1 522 biquadratic fields) |
| `odd_walks.gp` | exact closed walks of odd length in quartic fields where no quadratic subfield explains them, and only even relations in three 2-colourable controls |
| `admissible.py` | the admissible fields of Section 4 (`d = 167, 887, 1055, 1319, 1823`, and biquadratic ones) and the real-place example `d = 186 023` |

`tests/test_quadratic_planes.py` runs them (the PARI/GP tests skip when `gp` is not installed). Separate agents
refereed the proof in the working note and then two successive versions of the draft; the last pass found no
mathematical error. The revision fixed the credit for the "if" direction and for Lemma 2, which the first version
had attributed only to Madore and to our own notes. Nobody outside the project has refereed the draft.
