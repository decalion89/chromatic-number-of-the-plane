# Which real quadratic planes need four colours

[`four-colours.pdf`](four-colours.pdf) is a six-page draft (3 October 2026):

- **Theorem 1.** `χ(ℚ(√d)²) ≥ 4` for every positive integer `d ≡ 11 (mod 12)`.
- **Corollary 2.** For squarefree `d ≥ 2`, the plane over `ℚ(√d)` needs four colours exactly when
  `d ≡ 11 (mod 12)`; with the upper bounds of Fischer (1990) and Moorhouse (2010), `χ(ℚ(√d)²) = 4` for all these `d`
  except possibly `d ≡ 47, 143, 167 (mod 168)`, and `4 ≤ χ(ℚ(√47)²) ≤ 5`.
- **The proof** is explicit: one vector (`d ≡ 23 mod 24`) or two (`d ≡ 11 mod 24`) and their rotations by the
  rational rotations with denominators dividing `5^k` generate a group with no character mapping them into
  `[1/3, 2/3]`; Theorem W of [`papers/winding/`](../winding/) then excludes 3-colourings. The main step describes
  exactly the characters of `(1/5^k)ℤ[i]` that keep every such rotation in `[1/3, 2/3]`. For `d ≡ 23 (mod 24)` the
  bound on `k` is sharp for `k ≤ 5` and not for `k ≥ 6` (a third reading, 4 October).

[`four-colours.tex`](four-colours.tex) is its LaTeX source (`amsart`). To rebuild the PDF, run in this folder, with
TeX Live:

```sh
pdflatex four-colours.tex && pdflatex four-colours.tex
```

The fuller account is [`notes/four_colours_11_mod_12.md`](../../notes/four_colours_11_mod_12.md). The programs that
check the computations behind the proof, and exact certificates for small `d`, are in
[`data/quadratic_planes/winding/family/`](../../data/quadratic_planes/winding/family/); the tests are
`tests/test_winding_family.py`. Two internal referees (separate AI agents, with their own programs) checked the
proofs; the draft has not been refereed by anyone outside the project.
