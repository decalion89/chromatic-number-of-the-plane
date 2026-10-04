# Three colours for the plane over a number field

[`three-colours.pdf`](three-colours.pdf) is a 26-page draft (4 October 2026) about the graph on `F²`, for a
number field `F`, in which two points are adjacent when `(x − x′)² + (y − y′)² = 1`:

- **Theorem B.** `χ(F²) ≤ 3` if and only if some prime of `F` above 2 ramifies in `F(i)` or some prime of `F` above 3
  has residue degree 1 (and `χ(F²) ≤ 2` exactly in the first case). The "if" part is known (colourings through a
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
  [`papers/winding/`](../winding/)); the upper bound is a colouring through the residue field `𝔽₇`. Below 4, the best
  upper bound on `χ_c` that a locally constant colouring of one completion gives is 2, 3 or `7/2` (Propositions 5
  and 6, Corollary 5; Weil's bound for Kloosterman sums and a computation for `p < 3000`). `χ_c(ℚ(√59)²) ≥ 53/15` with 70
  explicit unit vectors (Proposition 4).
- **Theorem D.** If `χ(F²) ≥ 4`, then `χ_c(F²) ≥ 56/17 ≈ 3.294`, by hand (Section 9): the probe of Theorem B keeps
  its shape for the intervals `[θ, 1 − θ]` with `θ > 17/56` (Proposition 7). Below `θ = 3/10` it has 5-adic
  characters, so this method stops at `10/3`; with exact computations by three programs it reaches `3.3315`
  ([`data/number_fields/circular/probe/`](../../data/number_fields/circular/probe/)).
- **Theorems E and F.** For every number field `F`, `χ_c(F²)` is 2, 3, `7/2` or at least 4: 2 or 3 as `χ(F²)`
  (Theorem B), `7/2` exactly when `χ(F²) ≥ 4` and some prime of `F` above 7 has residue degree 1, and at least 4
  otherwise; in particular `χ_c(F²) ≥ 7/2` whenever `χ(F²) ≥ 4` (Theorem E). So the circular chromatic number is
  decided by a single completion below 4, as the chromatic number is at three colours; for example
  `χ_c(ℚ(√23)²) = 7/2`, `χ_c(ℚ(√59)²) = 4` and `χ_c(ℚ(√47)²) ≥ 4`. Section 10 adds the rotation
  `σ = (5 + 12i)/13` to the probe. One computer-assisted step for each theorem, an exact finite computation along
  the rotations `iᵃρʲσˡ` modulo 65 and 325 (Propositions 8 and 9), shows that the probe has only the types c and q
  above `2/7`, and the types c, q and a 7-adic one above `1/4`; Lemmas 12–19 (by hand) and compactness at 2, 3 and 7
  give the theorems; Theorem F needs only Proposition 9, which also gives Theorem E again. Corollary 6: for every
  finite extension `K` of `ℚ_p`, `χ_c(K²)` is 2, 3, `7/2` or at least 4 by explicit local conditions, so below 4 the
  circular chromatic number of `F²` is the least one of the planes over the completions of `F`. Corollary 7: below 4
  the value is also the circular chromatic number of some finite unit-distance graph in `F²`, by Lemma 22 (for a
  finite connection set and `2 < χ_c < 4`, every homomorphism to `K_{χ_c}` has a tight cycle, so a finite subgraph
  attains `χ_c`) and a compactness argument that gives a finite set of unit vectors with `κ = 1/χ_c`; no bound on the
  size. The certificates are in [`data/number_fields/circular/twoprime/`](../../data/number_fields/circular/twoprime/).

The working notes are [`notes/three_colours_number_fields.md`](../../notes/three_colours_number_fields.md) and
[`notes/circular_planes.md`](../../notes/circular_planes.md); the programs and certificates are in
[`data/number_fields/`](../../data/number_fields/), and the tests in `tests/test_three_colours.py`,
`tests/test_gap_above_three.py` and `tests/test_two_primes.py`. Theorem B was refereed twice, the paper once (with a
second look at Section 8), Theorem D twice (the note's version and Section 9), and Theorems E and F once each, by
separate AI agents with their own programs; no mathematician has checked it yet. A third referee then read Section 10
as written in the paper, with programs of its own: no error in a proof; its corrections (two false sentences about
the 7-adic characters, a bound quoted for the wrong range, the dependence of Theorem F on Proposition 8, now removed
by eight new certificates, and smaller gaps) are applied. Corollary 6 was added after that reading. Corollary 7 was
refereed separately, in two versions (no mathematical error; the expository gaps are closed; the second, simpler
proof through Lemma 22 is the one in the paper); while preparing it we found that the
paper had wrongly said that no finite unit-distance graph with circular chromatic number `7/2` was known: the Moser
spindle is one (over `ℚ(√3, √11)`). A last reading of the whole paper by a separate agent found no mathematical
error; its corrections are applied: Lemma 13 (the arithmetic of `𝔽₄₉` at the places above 7), until then checked
only by computer, now has a proof by hand, so Propositions 8 and 9 are the only computer-assisted steps of Theorems
E and F; the sentence on `ℚ(√47)` now says `4 ≤ χ_c ≤ 19/4`; a remark after Corollary 2 holds for real quadratic
fields only (`ℚ(√−73)` is split above 2 and 3); and attributions (Isbell, Fischer 1994, the special cases of (a)),
references and wording are corrected. A second agent refereed the hand proof of Lemma 13 with a program of its own
(correct, no gap; `twoprime/indep_L13/`) and corrected one more sentence: Proposition 6 rests on a computation of
`κ₁` for the primes `11 ≤ p < 1001`. Section 10 now gives explicit finite witnesses for Corollary 7: over
`ℚ(√11)` the unit-distance graph on `A + A`, `A` the vertex set of the 76-vertex graph (2 237 vertices), and, from a
colouring-guided growth followed by vertex deletion, a vertex-critical one with 170 vertices; over `ℚ(√455)` and
`ℚ(√191)` ones with 175 and 293 vertices. Each has `χ_c = 7/2`, certified by two DRAT proofs for two encodings
(`data/number_fields/circular/finite_witness/`). For the value 3, over `ℚ(√7)`, whose plane has no unit triangle,
Section 10 gives nine points whose unit-distance graph, the Wagner graph with one chord subdivided, has `χ_c = 3`,
with a proof by hand; it is vertex-critical, no triangle-free graph with fewer vertices has `χ_c = 3`, and exactly
three triangle-free graphs with nine vertices have `χ_c = 3` (another is a unit-distance graph over `ℚ(√31)`). A
referee checked the `7/2` witnesses again with programs of its own (no error; its remarks applied: `χ_c = 7/2` over
`ℚ(√191)` and `ℚ(√455)` now follows from the witnesses and Proposition 3, without Proposition 9). The PDF has 27
pages.

[`three-colours.tex`](three-colours.tex) is the LaTeX source (`amsart`). To rebuild the PDF, run in this folder,
with TeX Live:

```sh
pdflatex three-colours.tex && pdflatex three-colours.tex
```
