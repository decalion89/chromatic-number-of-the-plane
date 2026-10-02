# What the literature has, and what is new here (25 September 2026, updated 28 September)

This note compares the project with the published record. Sources are:
- the papers below;
- all 18 Polymath16 threads;
- the Polymath16 wiki, including its page
  [Algebraic formulation of Hadwiger–Nelson problem](https://web.archive.org/web/20210412075722/https://asone.ai/polymath/index.php?title=Algebraic_formulation_of_Hadwiger-Nelson_problem),
  which colours rings such as the Moser ring (4 colours through `ℤ₄[ω]`, all
  4-colourings of period 8) but not whole planes;
- the public repositories of the 2026 search efforts.

"New" means only that none of these sources has it. Nothing here has been
refereed.

## The strategy for six, and who had it

The main route is the reduction of Exoo and Ismailescu, which Polymath16 calls
clamping onto a virtual edge:
- a two-distance witness `W` with `χ(ℝ², {1, d}) ≥ 6`;
- a unit-distance gadget `H` that forces a pair at distance `d` apart.

The route is **not ours**.

| ingredient | who had it | where |
|---|---|---|
| The reduction `W + H ⇒ χ(ℝ²) ≥ 6` | Exoo–Ismailescu, who proved five this way with `d = √(11/3)` (arXiv 1805.00157) and stated the step to six as Conjecture 8.1 of arXiv 1805.06055; Polymath16 ("virtual edges", "clamping"; de Grey's definitions in thread 3, comment 4172) | [arXiv 1805.00157](https://arxiv.org/abs/1805.00157), [arXiv 1805.06055](https://arxiv.org/abs/1805.06055), [thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/#comment-4172) |
| `W` for `d = φ` (31 vertices) | Huddleston; Parts | Owings–Tetiva–Huddleston, *Coloring the plane* (Problem 11236), Amer. Math. Monthly 115(2) (2008) 170–172, [JSTOR](https://www.jstor.org/stable/27642435); [arXiv 2010.12656](https://arxiv.org/abs/2010.12656) |
| `W` for `d = 2` | Exoo–Ismailescu | [arXiv 1909.13177](https://arxiv.org/abs/1909.13177) |
| `W` for `d = √3` (33 vertices) and `(√6 + √2)/2` (117 vertices) | Polymath16, thread 14, comment 24460 (21 October 2019); De Neve et al. | [comment 24460](https://dustingmixon.wordpress.com/2019/08/05/polymath16-fourteenth-thread-automated-graph-minimization/#comment-24460); Amer. Math. Monthly 132(10) (2025) 1007–1022, [doi](https://doi.org/10.1080/00029890.2025.2559554) |
| Which sets of distances `d` are the range of a two-distance graph (the semialgebraic sets with positive lower and upper bounds) | Ágoston | [arXiv 2601.07828](https://arxiv.org/abs/2601.07828) |
| Colouring-guided growth and minimisation, at 4 colours | Heule; Parts | [arXiv 1805.12181](https://arxiv.org/abs/1805.12181), [1907.00929](https://arxiv.org/abs/1907.00929), [2010.12665](https://arxiv.org/abs/2010.12665) |
| Unions of paths as a construction | Haugland (Moser-spindle-free, 2 131 vertices) | [arXiv 2608.04542](https://arxiv.org/abs/2608.04542) |
| Reduction of `F²` modulo a prime, for number fields `F`, and extension from a subgroup or a ring to the whole field by cosets | Woodall (`ℚ²` modulo 2); Fischer (1990, Thm 1: cosets of the component of the origin); Moorhouse (Lemma 4.2, Lemma 8.2); Madore (Cor. 3.4, Prop. 3.2, 3.8, ¶6.6) | Woodall: J. Combin. Theory Ser. A 14 (1973) 187–200; [Moorhouse 2010](https://www.ericmoorhouse.org/pub/chromatic.pdf), [arXiv 1509.07023](https://arxiv.org/abs/1509.07023) |
| `χ(ℚ(√N)²) = 2` for `N ≡ 1, 2 (mod 4)` | Johnson (1987) | Congr. Numer. 60, 51–58; see [Payne](https://arxiv.org/abs/0707.1177) |
| `χ(ℚ(√N)²) ≤ 3` for `N ≡ 0, 1 (mod 3)`, `≤ 4` for `N ≡ 3 (mod 8)` (Thms 9, 10; the second is a reduction at 2 into `ℤ/4`); colouring the plane through the component of the origin (Thm 1) | Fischer (1990) | Discrete Math. 82, 181–195 (read in full); see also [Payne](https://arxiv.org/abs/0707.1177) |
| No additive `k`-colouring of `ℚ(√N)²` for `k ≤ 6` when `N ≡ −1 (mod 24)`, for example `N = 47`, because `1/2` and `1/3` are sums of unit vectors | Fischer (1990), Thm 10(ii) | Discrete Math. 82, 181–195 |
| **`χ(ℚ(√3, √11)²) = 4`**; more generally an additive 4-colouring of `ℚ(√p, √q)²` for squarefree, relatively prime `p ≡ 3`, `q ≡ 11 (mod 16)` with `pq ≡ 1 (mod 32)` | Fischer (1994) | Congr. Numer. 104, 73–79; [Zbl 0836.05030](https://zbmath.org/?q=an:0836.05030) |
| The connected component of the origin in `ℚ(√N₁, …, √N_d)²` | Fischer (1990) | Congr. Numer. 72, 213–221; [Zbl 0733.05048](https://zbmath.org/?q=an:0733.05048) |
| A survey, as of 2000, of problems on colourings of `ℚⁿ` and its algebraic extensions, posed in or arising from Benda–Perles, *Colorings of metric spaces* | Johnson (2000); we have not seen it | Geombinatorics 9(4), 170–179; [Zbl 0974.05029](https://zbmath.org/?q=an:0974.05029); Benda–Perles: Geombinatorics 9(3), 113–126, [Zbl 0951.05037](https://zbmath.org/?q=an:0951.05037) |
| 2-adic 4-colourings of the Moser ring, with colours in `𝔽₄` | Speyer (thread 2, April 2018); Hubai's analysis and computer search, reported by Gibbs: all have period 8 (thread 3); Dúcz (2026) | [thread 2](https://dustingmixon.wordpress.com/2018/04/22/polymath16-second-thread-what-does-it-take-to-be-5-chromatic/#comment-4013), [thread 3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/), [arXiv 2606.12325](https://arxiv.org/abs/2606.12325) |
| A colouring of `ℤ[ζ₂₄, 1/3]` through `(ℤ/4)[ζ₂₄]`, proposed for the plane over `ℚ(√2, √3)` | Voronov (thread 17, 30 July 2021), who asked for "a simpler way" | [comment 29476](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/#comment-29476) |
| A 4-chromatic unit-distance graph over `ℚ(√2, √3)` | Voronov–Neopryatnaya–Dergachev (`L₁₀,₂` and its Minkowski sums) | [arXiv 2106.11824](https://arxiv.org/abs/2106.11824) |
| The question whether a field generated by two square roots of primes can carry a 5-chromatic unit-distance graph | Voronov (thread 17, 17 July 2021) | [comment 29283](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/#comment-29283) |
| A Galois automorphism that preserves unit distance maps one distance to another, so what holds for one distance holds for its conjugates (Tao used it to equate the densities `p_{d₁}` and `p_{d₂}` of two distances) | Tao (thread 7, comment 4893, 19 June 2018) | [comment 4893](https://dustingmixon.wordpress.com/2018/06/16/polymath16-seventh-thread-upper-bounds/#comment-4893) |
| A 103-vertex graph with edges at 1 and `2/√3` and no 4-colouring | Exoo and Ismailescu (Polymath16, thread 3, comment 4161) | [comment 4161](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/#comment-4161) |
| The eigenvalues of the finite Euclidean graphs, through Gauss and Kloosterman sums, and Weil's bound for them | Medrano–Myers–Stark–Terras | J. Comput. Appl. Math. 68 (1996) 221–238, [doi](https://doi.org/10.1016/0377-0427(95)00261-8) |
| The three-point semidefinite bound | Schrijver | IEEE Trans. Inform. Theory 51 (2005) 2859–2866 |
| A lattice-like graph with no 5-colouring and bichromatic origin would give six | Frankl–Hubai–Pálvölgyi (Thm 27) | [arXiv 1912.02604](https://arxiv.org/abs/1912.02604) |
| A 6-chromatic unit-distance graph has at least 42 vertices | de Grey–Parts | [arXiv 2303.14714](https://arxiv.org/abs/2303.14714) |
| Large negative searches for six (2026) | Adler (`zeta42`); the `math-market/chromatic-plane` public repository | [keithadler/zeta42](https://github.com/keithadler/zeta42), [math-market/chromatic-plane](https://github.com/math-market/chromatic-plane) |

## What we found nowhere else, and what turned out to be known

1. **`χ(ℚ(√2, √3)²) = 4`, and a short proof of Fischer's `χ(ℚ(√3, √11)²) = 4`.**

   **`χ(ℚ(√3, √11)²) = 4` is not new.** K. G. Fischer proved it in 1994
   (*A planar geometric graph of chromatic number four*, Congr. Numer. 104,
   73–79; Zbl 0836.05030). He proved that `ℚ(√p, √q)²` has an additive 4-colouring, with values in
   `ℤ/4`, for squarefree, relatively prime `p ≡ 3`, `q ≡ 11 (mod 16)` with
   `pq ≡ 1 (mod 32)`. zbMATH's summary adds that the plane contains finite
   graphs which require four colours; the Moser spindle is one. We found this only
   after the note had been sent out; we have read the zbMATH summary, not the
   paper. Johnson's 2000 survey (Geombinatorics 9(4), 170–179), on the problems of Benda and Perles (*Colorings of metric spaces*, Geombinatorics 9(3), 113–126), may record further
   results; we have not seen it. The result seems to have been overlooked by:
   - Moorhouse (2010): "We have not determined the exact value …";
   - Madore (2015), who proved `4 ≤ χ ≤ 5`;
   - Exoo–Ismailescu (2018), who asked whether a 5-chromatic unit-distance
     graph embeds in this plane;
   - Cranston–Rabern ([arXiv 1501.01647](https://arxiv.org/abs/1501.01647),
     Combinatorica 2017), who cite Fischer (1990) and ask for the fractional and
     the ordinary chromatic number of this plane;
   - Polymath16, thread 3 (2018), and Parts' "funny proof" in thread 13
     (2019);
   - Voronov, Polymath16 thread 17 (July 2021): "it seems likely that
     `χ(Q(i, √3, √11)) = 4` ... But as far as I know, nobody has proved this
     yet."

   We have since read Fischer's 1990 paper in full. It treats `ℚᵈ` and the
   quadratic fields `ℚ(√N)` only, and says nothing about `ℚ(√2, √3)`.

   Earlier, Fischer (*Additive K-colorable extensions of the rational plane*,
   Discrete Math. 82 (1990)) had shown `χ(ℚ(√N)²) ≤ 3` for `N ≡ 0, 1 (mod 3)`
   and `χ(ℚ(√N)²) ≤ 4` for `N ≡ 3 (mod 8)`, and Johnson (Congr. Numer. 60 (1987))
   `χ(ℚ(√N)²) = 2` for `N ≡ 1, 2 (mod 4)`, as summarised by Payne
   ([arXiv 0707.1177](https://arxiv.org/abs/0707.1177)).

   **What we had not found in the literature when we wrote the note, and what we found later.**
   On 28 September we found a public repository of July 2026, *A 2-adic obstruction to 5-chromatic
   unit-distance graphs* ([MildlyMeticulous/hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction); work dated 21–22 July, repository
   created 30 July; not refereed, no arXiv version). It contains most of the upper bounds below:
   - its Theorem A: if `K(i)` has a conjugation-stable valuation ring over 2 with residue field `k`, then
     `χ(K²) ≤ |k|`. This is our criterion in general form: the reduction of `z = x + iy` at a place over 2;
   - its Corollary B: for a real multiquadratic `K`, `χ(K²) ≤ 4` unless the squarefree kernel of some
     product of the radicands is `≡ 7 (mod 8)`. Its Corollary B′ covers the field generated by all `√d` and
     `√(2d)` with `d` squarefree and `d ≡ 1, 3 (mod 8)`, which contains `ℚ(√2, √3)` and `ℚ(√3, √11)`;
   - `χ(ℚ(√3, √11)²) = 4` (without Fischer 1994), and a dichotomy for spindles of the triangular lattice
     (its Theorem C).

   It does not state the case `ℚ(√2, √3)`. It also points to Axenovich, Choi, Lastrina, McKay, Smith and
   Stanton, *On the chromatic number of subsets of the Euclidean plane* (Graphs Combin. 2014), Thm 2.3, a
   statement for lattices with a `7 (mod 8)` hypothesis on the generators; we have not seen that paper.

   So, item by item:
   - `χ(ℚ(√2, √3)²) = 4`, Voronov's second case (`notes/local_colourings.md`
     §10). Fischer's hypotheses exclude it. Voronov had proposed a route
     through `(ℤ/4)[ζ₂₄]` and asked for a simpler one; the lower bound was
     already implicit in Voronov–Neopryatnaya–Dergachev. **The upper bound is a
     case of Corollary B′ of the repository above**, two months before our note.
     What is ours is the explicit statement and its Lean proof.
   - A partial answer to Voronov's question on two square roots of primes:
     `χ(ℚ(√3, √q)²)` is 3 for `q ≡ 1 (mod 3)` and 4 for `q = 2` or
     `q ≡ 11, 17 (mod 24)`, and at least 4 otherwise; and
     `4 ≤ χ(ℚ(√3, √5)²) ≤ 5` (`notes/local_colourings.md` §11). For
     `q ≡ 11 (mod 32)` the upper bound 4 is also a case of Fischer (1994).
     **Every upper bound 4 here is also a case of Corollary B of the repository
     above** (the products `3`, `q`, `3q` are never `≡ 7 (mod 8)` in these cases).
   - The criterion behind both proofs: if `√3 ∈ L` and some prime of `L` above 2
     has residue field `𝔽₂`, then `χ(L²) ≤ 4`. It is the case of this form of Madore's ¶6.6,
     his Prop. 3.2 for any quadratic form, in the coordinates `α = x + y/√3`,
     `β = 2y/√3`, where the squared distance is the form `α² − αβ + β²`,
     anisotropic modulo such a prime. For an expert in
     local fields this is a short step.
     **It is a case of Theorem A of the repository above.**

   Related work: Speyer used the same 2-adic reduction in
   [Polymath16, thread 2](https://dustingmixon.wordpress.com/2018/04/22/polymath16-second-thread-what-does-it-take-to-be-5-chromatic/#comment-4013)
   (25 April 2018) to 4-colour the Moser ring; Hubai's analysis and computer search, reported
   by Gibbs in thread 3, found that all its 4-colourings have period 8;
   Fischer's Thm 1 (1990), Moorhouse's Lemma 4.2 and Madore's Prop. 3.2 pass
   from a subgroup or a ring to the whole field by cosets; Dúcz (2026) gave
   geometric 4-colourings of the Moser lattice and ring.

   By-product: a necessary condition at 2, 3 and 7 for a real field to be
   5-chromatic.
2. **Whole-field chromatic numbers of two CM fields: known, with local proofs.**
   - `χ(ℚ(√−3, √−11)) = 4`, the Moser field. The field lies in the plane over
     `ℚ(√3, √11)`, so this follows from Fischer (1994). Our proof extends
     Speyer's colouring of the Moser ring to the whole field.
   - `χ(ℚ(√−3, √−11, √−247)) = 5`, at the place over 11. The upper bound
     follows from Madore's Cor. 3.4 and Lemma 4.5, by the argument of his
     Prop. 4.6; Exoo–Ismailescu's 5-chromatic graph lies in this field after
     a quarter turn (`notes/local_colourings.md` §3).

   What is ours is the local form of the argument and its consequence: no
   search in the second field can reach six. `scripts/fieldscreen.py` lists
   fields with no small non-split place, such as `L16`.

   Dúcz ([arXiv 2606.12325](https://arxiv.org/abs/2606.12325)) 4-coloured the
   Moser *lattice* and *ring*.
3. **Repulsion spectra at five colours, and the two-step explanation.** Every
   known witness distance is `|u + v|` for unit vectors `u, v`. Such pairs share
   a unit neighbour, and 5-colourings colour them alike 30–40% of the time.
   `2/√3` is the most repulsive distance measured, at 8%, and it is not a
   two-step distance in `L16`. Exoo and Ismailescu had already found a 103-vertex graph
   with edges at 1 and `2/√3` and no 4-colouring (Polymath16, thread 3, comment 4161).
   The quantity measured, how often two points at distance `d` get the same
   colour, is an empirical counterpart of the `p_d` of Polymath16's probabilistic
   formulation (Tao, thread 7, comment 4893; Ágoston, [arXiv 2112.07665](https://arxiv.org/abs/2112.07665)).
4. **A sparse 72-point lattice witness** for
   `{1, 4/√3, √7, √(28/3), √(61/3)}`, verified by four solvers and DRAT.
5. **Methods:**
   - rainbow growth at five colours;
   - `MODE` growth toward a forced pair;
   - skeleton growth around a witness;
   - Kempe backbones.
6. **Six colours for finite planes** (`notes/local_colourings.md` §14):
   `χ ≥ 6` for `𝔽₃₇²`, `𝔽₄₁²`, `𝔽₄₃²`, `𝔽₄₇²` and the anisotropic planes
   `G₂₉`, `G₃₇`, `G₄₁`, from Schrijver's three-point bound with checked dual
   certificates. `G₁₃` joins them, by a SAT proof rather than the three-point
   bound: `χ(G₁₃) = 6` (`notes/g13_chi.md`). Moorhouse's table stops at `q = 17`. We found no earlier bound
   of six for these planes, and no earlier use of the three-point bound for
   finite unit-distance graphs. For large `q` Hoffman's bound gives six, with
   the spectra of Medrano–Myers–Stark–Terras and, for the anisotropic planes, of
   Bannai–Shimabukuro–Tanaka (Discrete Math. 309 (2009) 6126–6134).
   Linear-programming bounds strengthened by triangle constraints were used
   earlier for planar sets avoiding unit distance: Keleti–Matolcsi–de Oliveira
   Filho–Ruzsa (Discrete Comput. Geom. 55 (2016) 642–661) and DeCorte–de Oliveira
   Filho–Vallentin (Math. Program. 191 (2022) 487–558).
7. **A Moser-spindle-free 5-chromatic unit-distance graph with 852 vertices**
   (`notes/flat852.md`), in `ℚ(ζ₂₁)`, from the unit vectors of Haugland's heptagon
   graph and their mirror images. The spindle-free 5-chromatic graphs we found
   before it have 2 131 vertices (Haugland, [arXiv 2608.04542](https://arxiv.org/abs/2608.04542)),
   1 441 (cited there), 1 435
   ([ruturajr-raval/hadwiger-nelson-spindle-free](https://github.com/ruturajr-raval/hadwiger-nelson-spindle-free/releases/tag/v0.1.0),
   9 September 2026) and 1 299 (posted on GitHub by zach7036). That Haugland's
   directions alone are 4-colourable is the reduction of Theorem A of
   hn-2adic-obstruction, applied to them; we found no earlier statement of it for
   this family. We searched arXiv and GitHub on 28 and 29 September 2026.
8. **Real quadratic planes that need four colours** (`notes/quadratic_planes.md`):
   `χ(ℚ(√d)²) = 4` for `d = 11, 23, 35, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959`, and `4 ≤ χ(ℚ(√47)²) ≤ 5`. For real
   quadratic fields the values in the literature are 2 (Johnson 1987; Moorhouse
   2010, §8) and 3 (Fischer 1990; Madore 2015); for `d ≡ 3 (mod 4)` the lower
   bound was 3 (Fischer 1990, Theorem 8). The upper bounds are known (Moorhouse's
   Corollary 8.3, Fischer's Theorem 10, reduction at 11); the lower bounds of four
   are ours. hn-2adic-obstruction (item 1) proves the upper bound 4 for
   `d ≡ 3 (mod 8)` again and finds odd cycles over `ℚ(√11)`; it gives no lower
   bound of four. D. Cohen (University of Chicago REU paper, 2007) conjectured
   that `ℚ[α] ⊂ ℂ` is 3-colourable for every quadratic `α`. That set has dimension
   2 over `ℚ`, not 4 like the plane `ℚ(√d)²`; reduction at a ramified prime
   proves the conjecture (`notes/quadratic_planes.md` §1). Moorhouse asked in a
   talk in 2010 what can be said about `χ(ℚ(√d)²)`, and about `χ(ℚ(√47)²)` in
   particular. We searched on 1 October 2026; the sources we read, and those we
   could not read, are listed in `notes/quadratic_planes.md` §1.

## What this means

The strategy for six is the known one. The new pieces are:
- the choice of distance;
- the arithmetic that rules fields out for six;
- the theorem χ(ℚ(√2, √3)²) = 4, and a short proof of Fischer's
  χ(ℚ(√3, √11)²) = 4 (1994), which later work had treated as open. The upper
  bounds, and the reduction of `z = x + iy` at 2 behind them, are also in a
  public repository of July 2026 that we found on 28 September (item 1);
- lower bounds of six for finite planes, from the three-point bound;
- real quadratic planes that need four colours (item 8).

All are modest. Voronov raised the case ℚ(√2, √3) in a Polymath16 comment.
The questions about ℚ(√3, √11) asked in print, by Exoo–Ismailescu and by
Cranston–Rabern, had been settled by Fischer.
