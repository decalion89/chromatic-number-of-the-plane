# What the literature has, and what is new here (25 September 2026)

This note compares the project with the published record. Sources are:
- the papers below;
- all 18 Polymath16 threads;
- the Polymath16 wiki;
- the public repositories of the 2026 search efforts.

"New" means only that none of these sources has it. Nothing here has been
refereed.

## The strategy for six, and who had it

The main route is the Exoo–Ismailescu two-step reduction:
- a two-distance witness `W` with `χ(ℝ², {1, d}) ≥ 6`;
- a unit-distance gadget `H` that forces a pair at distance `d` apart.

The route is **not ours**.

| ingredient | who had it | where |
|---|---|---|
| The reduction `W + H ⇒ χ(ℝ²) ≥ 6` | Exoo–Ismailescu; Polymath16 ("virtual edges", "clamping") | [arXiv 1909.13177](https://arxiv.org/abs/1909.13177), [thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/) |
| `W` for `d = φ` (31 vertices) | Huddleston; Parts | [arXiv 2010.12656](https://arxiv.org/abs/2010.12656) |
| `W` for `d = 2` | Exoo–Ismailescu | [arXiv 1909.13177](https://arxiv.org/abs/1909.13177) |
| `W` for `d = √3` (33 vertices) and `(√6 + √2)/2` (117 vertices) | Ágoston–Pálvölgyi (Polymath16); De Neve et al. | [Monthly 2025](https://doi.org/10.1080/00029890.2025.2559554) |
| The range of `d` with `χ(ℝ², {1, d}) ≥ 6` | Ágoston | [arXiv 2601.07828](https://arxiv.org/abs/2601.07828) |
| Colouring-guided growth and minimisation, at 4 colours | Heule; Parts | [arXiv 1805.12181](https://arxiv.org/abs/1805.12181), [1907.00929](https://arxiv.org/abs/1907.00929), [2010.12665](https://arxiv.org/abs/2010.12665) |
| Unions of paths as a construction | Haugland (Moser-spindle-free, 2 131 vertices) | [arXiv 2608.04542](https://arxiv.org/abs/2608.04542) |
| Reduction of `F²` modulo a prime, for number fields `F` | Moorhouse (Lemma 8.2); Madore (Prop. 3.2, 3.8, 6.6) | [Moorhouse 2010](https://www.ericmoorhouse.org/pub/chromatic.pdf), [arXiv 1509.07023](https://arxiv.org/abs/1509.07023) |
| 2-adic 4-colourings of the Moser ring, with colours in `𝔽₄` | Speyer; Gibbs and Hubai (all have period 8); Dúcz (2026) | [Polymath16 thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/), [arXiv 2606.12325](https://arxiv.org/abs/2606.12325) |
| A lattice-like graph with no 5-colouring and bichromatic origin would give six | Frankl–Hubai–Pálvölgyi (Thm 27) | [arXiv 1912.02604](https://arxiv.org/abs/1912.02604) |
| A 6-chromatic unit-distance graph has at least 42 vertices | de Grey–Parts | [arXiv 2303.14714](https://arxiv.org/abs/2303.14714) |
| Large negative searches for six (2026) | Adler (`zeta42`); the Math Market bounty | [keithadler/zeta42](https://github.com/keithadler/zeta42), [math-market/chromatic-plane](https://github.com/math-market/chromatic-plane) |

## What we found nowhere else

1. **`χ(ℚ(√3, √11)²) = 4`.** This closes a question left open by:
   - Moorhouse (2010): "We have not determined the exact value";
   - Madore (2015), who proved `4 ≤ χ ≤ 5`;
   - Exoo–Ismailescu (2018), who asked whether a 5-chromatic unit-distance
     graph embeds in this plane;
   - Polymath16, thread 3 (2018), and Parts' "funny proof" in thread 13
     (2019);
   - Voronov, Polymath16 thread 17 (July 2021): "it seems likely that
     `χ(Q(i, √3, √11)) = 4` ... But as far as I know, nobody has proved this
     yet."

   The proof is the Hermitian form of the local reduction at an inert place
   over 2. See `notes/local_colourings.md` §8, `hn/adelic.py` and
   `tests/test_q311.py`.

   **The idea is not new; the theorem is.**
   - David Speyer used the same 2-adic reduction in
     [Polymath16, thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/)
     (May 2018) to 4-colour the Moser ring, with colours in `ℤ[ω]/2 = 𝔽₄`.
   - Gibbs and Hubai found that all such colourings have period 8.
   - Dúcz (2026) 4-coloured the Moser lattice and ring.

   What is new is the step to the whole plane: every unit vector of
   `ℚ(i, √3, √11)` is a 2-adic unit, because the place over 2 is inert, and
   denominators of 2 are handled by cosets.

   By-products:
   - `χ(ℚ(√d)²) ≤ 4` for `d ≡ 3 (mod 8)`, which extends Moorhouse's Theorem 8.1;
   - a necessary condition at 2, 3 and 7 for a real field to be 5-chromatic.
2. **Whole-field chromatic numbers of CM fields.**
   - `χ(ℚ(√−3, √−11)) = 4`, the Moser field. This is Speyer's colouring of
     the Moser ring, extended to the whole field;
   - `χ(ℚ(√−3, √−11, √−247)) = 5`, at the place over 11.

   Consequently no search in the second field can reach six. The finite
   planes `G_q` satisfy `χ(G_q) ≥ 6` for `q ≥ 53` (Hoffman), and
   `scripts/fieldscreen.py` lists fields with no small non-split place, such as
   `L16`.

   Dúcz ([arXiv 2606.12325](https://arxiv.org/abs/2606.12325)) 4-coloured the
   Moser *lattice* and *ring*. Our Theorem 1 covers the whole field, and
   Theorem 3 the 8-dimensional plane over `ℚ(√3, √11)`.
3. **Repulsion spectra at five colours, and the two-step explanation.** Every
   known witness distance is `|u + v|` for unit vectors `u, v`. Such pairs share
   a unit neighbour, and 5-colourings colour them alike 30–40% of the time.
   `2/√3` is the most repulsive distance measured, at 8%, and it is not a
   two-step distance in `L16`.
4. **A sparse 72-point lattice witness** for
   `{1, 4/√3, √7, √(28/3), √(61/3)}`, verified by three solvers and DRAT.
5. **Methods:**
   - rainbow growth at five colours;
   - `MODE` growth toward a forced pair;
   - skeleton growth around a witness;
   - Kempe backbones.

## What this means

The strategy for six is the known one. The new pieces are:
- the choice of distance;
- the arithmetic that says which fields can hold six;
- a first by-product theorem.

That theorem is modest, but it answers questions that were asked in print.
