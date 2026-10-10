# Measurable density bounds hold for all colourings of the plane with algebraic forbidden distances

**Status (10 October 2026).** Proved on paper, using OpenAI's rigidity theorem (their Theorem 2.3, formalized in Lean
by OpenAI); the transfer of Theorem A is formalized in Lean on top of their library
([`lean/external/openai-five/distances/`](../lean/external/openai-five/distances/README.md)); the density bounds are
computer proofs checked by several programs.

- **Referees.** Six AI-assisted referees inside the project. The first read the transfer theorems (correct, minor
  gaps, fixed). The second read the corollaries and the certificates with two checkers of its own (no mathematical
  error, every certificate valid; the priority claims were wrong and are corrected). The third read version 3 with a
  checker of its own (no mathematical error; Croft's constant was rounded up, one inequality was classical and
  uncredited, the title said more than the theorem; all fixed in version 4). The fourth read the two-scale theorem
  of version 5, its strong form, the explicit pairs, the Lean paragraph and the 48-colour certificate (no
  mathematical error, no prior occurrence in the literature; wording and attribution fixes, applied). The fifth read
  the uniform bounds of version 6 (Theorem 5, Lemma 18, the certificates and both programs) with a checker of its
  own sharing no code or Bessel routine with ours: no mathematical error, all ten certificates valid, no earlier bound
  of 7 for a pair; the measurable statement and one half of the limit in the remark are now argued, and an overclaim
  in the abstract is corrected. The sixth read version 7 (triangles and spindles in the uniform certificates) with a
  checker of its own: no mathematical error, all ten certificates valid (the smallest margin, `min h ≥ 5.9·10⁻⁷`, is
  for 14 colours), all eight tails rebuilt byte for byte; it corrected the value formula and an overclaim in the
  remark, Corollary 4 against the abstract, and statements in the README and these notes.
- **No human referee yet.**

The paper is [`papers/forbidden-distances/`](../papers/forbidden-distances/README.md); this note is the record
behind it: the statements, what is new and what is not, and what each referee found.

## The observation

OpenAI's proof that χ(ℝ²) ≥ 6
([*The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf),
23 September 2026, below [OAI]) has two parts. Their Theorem 1.3 transfers a proper `k`-colouring of the plane to a
*weak measurable* `k`-colouring (a measurable colouring whose classes contain only a null set of pairs at distance
1), and their Theorem 1.4 excludes weak measurable 5-colourings. The transfer rests on a rigidity theorem for
invariant measures on the character group of the algebraic plane ([OAI, Theorem 2.3]); its proof uses neither that
the forbidden distance is 1 nor that a colouring has finitely many colours rather than being fractional.

Write `F` for the real algebraic numbers, `χ(ℝ², D)` for the chromatic number of the plane with the distances in
`D` forbidden, and `m_D(ℝ²)` for the largest upper density of a measurable planar set with no two points at a
distance in `D` (`m₁ = m_{{1}}`).

- **Theorem A.** For every `D ⊆ F_{>0}`, a proper `k`-colouring of the algebraic plane with the distances `D`
  forbidden gives a weak measurable `k`-colouring of the plane avoiding `D`, and conversely. Hence
  `χ(ℝ², D) ≥ 1/m_D(ℝ²)`.
- **Theorem B.** For every `D ⊆ F_{>0}`, the finitary fractional chromatic number (the supremum over finite graphs
  drawn with edge lengths in `D`) is at least `1/m_D(ℝ²)`.
- **Theorem 1.** For finite `D ⊆ F_{>0}`: `χ(ℝ², D) ≥ χ_{f,0}(ℝ², D) = χ_f(ℝ², D) = 1/m_D(ℝ²)`. The upper bound
  `χ_f ≤ 1/m_D` is elementary: translates of a periodic set of density close to `m_D` (Lemma 13 of the paper).

Theorems A and B are [OAI, Section 4] with the unit circle replaced by the circles of radii `d ∈ D` (countably many,
since `D ⊆ F`) and, for B, with one indicator function in place of the `k` colour classes; the paper gives the
proofs.

## Consequences

**One distance.**
- `χ_{f,0}(ℝ²) = χ_f(ℝ²) = 1/m₁(ℝ²)`, so `4.0485 < χ_f(ℝ²) < 4.36` (Ambrus, Csiszárik, Matolcsi, Varga and
  Zsámboki: `m₁ ≤ 0.2470`, and their constants give `m₁ < 0.2469975`, so `χ_f > 4.0486`; Croft:
  `m₁ ≥ 0.2293647…`, so `m₁ ≥ 0.229364`; the value 0.229365 often quoted is rounded up).
- `f(n)`, the least independence number of the unit-distance graph of `n` points in the plane, satisfies
  `f(n) ≥ m₁ n` for every `n` (the classical averaging bound `m₁ ≤ α₁`) and `f(n)/n → m₁(ℝ²)`, which is new. The
  limit uses Theorem 2 of Matolcsi, Ruzsa, Varga and Zsámboki (the finitary fractional chromatic number equals the
  reciprocal of the finitary independence ratio) and Fekete's lemma.

**Two distances** (certified bounds of the next section, then Theorem A):
`χ(ℝ², {1, 2}) ≥ 7`, `χ(ℝ², {1, √3}) ≥ 7`, `χ(ℝ², {1, (√6+√2)/2}) ≥ 7`, `χ(ℝ², {1, 10}) ≥ 7`,
`χ(ℝ², {1, 11/6}) ≥ 8`, `χ(ℝ², {1, 46/25}) ≥ 8`, `χ(ℝ², {1, 40}) ≥ 9`, `χ(ℝ², {1, 100}) ≥ 10`,
`χ(ℝ², {1, 1000}) ≥ 12`; and `χ(ℝ², {1, d}) ≥ 17` for every large enough algebraic `d` (Bukh's theorem:
`m_{1,d} → m₁² < 1/16`), improved to 25 by the two-scale theorem below. Products of hexagonal colourings give
`χ(ℝ², {1, d}) ≤ 49`.

**Many distances.** For every `k` there are `k` positive integers `D = {1, t, …, t^{k−1}}` with
`χ(ℝ³, D) ≥ χ(ℝ², D) ≥ m₁(ℝ²)^{−k} > 4.0485^k` (Bukh's theorem: `m_D → m₁^k` as `t → ∞`), and
`sup{χ_f(ℝ², D) : D ⊆ F, |D| = k} = m₁(ℝ²)^{−k}` (Bukh also proves `m_D ≥ m₁^k` for every `k` distances). With the
product colouring, Erdős's `L₂(r)` satisfies `4.0485^r < L₂(r) ≤ 7^r` (and `L₂(r) > 1.48·4.0486^r` by the two-scale
theorem). Explicitly, from certified linear programs: `χ(ℝ², {1, 40, 1600}) ≥ 22`,
`χ(ℝ², {1, 10, 100, 1000}) ≥ 26`, `χ(ℝ², {1, 100, 10⁴}) ≥ 27` and `χ(ℝ², {1, 30, 900, 27000}) ≥ 48`. All of these hold
in every dimension `d ≥ 2`, by restriction to a plane.

## Two scales (version 5)

The density bounds can never exceed the fractional chromatic number, and for two distances
`χ_f(ℝ², {1, N}) = 1/m_{1,N} ≤ m₁^{−2} < 19.01` for every `N` (Bukh: `m_{1,N} ≥ m₁²`). Using OpenAI's theorem a
second time, as a statement about finite graphs, gives more.

- **Theorem C (two scales).** Let `D₁`, `D₂` be finite sets of positive real algebraic numbers, `c ≤ χ(ℝ², D₁)` and
  `W` a finite graph drawn in the plane with edge lengths in `D₂`. For all large algebraic `N`,
  `χ(ℝ², D₁ ∪ N·D₂) ≥ χ_c(W) ≥ c·χ_f(W)`, where `χ_c` is the `c`-fold chromatic number (Stahl's `c`-th chromatic
  number). Hence `χ(ℝ², D₁ ∪ N·D₂) ≥ χ_c(ℝ², D₂) ≥ c/m_{D₂}(ℝ²)` for all large `N`, with `c = χ(ℝ², D₁)`: the lower
  half of a geometric version of Geller and Stahl's formula `χ(G[H]) = χ_{χ(H)}(G)` for lexicographic products
  (J. Combin. Theory Ser. B 19 (1975)). Whether the upper half holds is open.
- **Corollary.** `χ(ℝ², {1, N}) ≥ 25` for every large algebraic `N` (`c = 6` by OpenAI, and a finite unit-distance
  graph `W` with `χ_f(W) > 4`; Dúcz and Varga's graph suffices), while `χ_f(ℝ², {1, N}) < 19.01`. Iterating with
  `χ_f(W) > 4.0486` (Theorem B and Ambrus et al.): integers `1 < d₂ < d₃ < d₄` with `χ ≥ 25, 102, 413` for
  `{1, d₂}`, `{1, d₂, d₃}`, `{1, d₂, d₃, d₄}`, and `r` integer distances with `χ ≥ 6·4.0486^{r−1}`, so
  `L₂(r) > 1.48·4.0486^r`. The bound `6/m₁` lies between 24.29 and 26.16 (it would be 26.16 if `m₁` is Croft's
  value, as Dúcz and Varga conjecture, giving 27).
- **Explicit pairs.** With de Grey's 1581-vertex graph (`c = 5`; all vertices within 2.2641 of one of them) and the
  Moser spindle (`χ₅ ≥ 18`), the proof gives `χ(ℝ², {1, N}) ≥ 18` for every algebraic `N ≥ 10^25` (threshold
  `6.75·10²⁴`); with the 607-vertex graph of Bellitto, Pêcher and Sédillot (`χ_f ≥ 1999983/512933`, so `χ₅ ≥ 20`)
  and the trivial edge bound, `χ(ℝ², {1, N}) ≥ 20` for `N ≥ 10^34`. These avoid `χ(ℝ²) ≥ 6` but still rest on
  Theorem A. The best certified LP bound for a pair is 12.

**Proof idea.** By Theorem A, a `k`-colouring gives measurable sets `A₁, …, A_k` avoiding `D₁ ∪ N·D₂` with conull
union. A finite graph `G` with edge lengths in `D₁` and `χ(G) ≥ c` (de Bruijn–Erdős) fits in a disc of radius `R₀`;
averaging over its translates shows that in every disc of radius `R = 2R₀` at least `c` sets `A_i` have density at
least `1/(4|V(G)|k)`. Let `S_i` be the set of centres where `A_i` has that density: every point lies in at least `c`
of them. Since `A_i` avoids `N·d`, its `R`-averaged density has pair density at distance `N·d` at most
`√(R/(πNd))` (Plancherel, with `|J₀(y)| ≤ √(2/(πy))` and `1 − (2J₁(s)/s)² ≤ min(1, s²/4)`; the true supremum of the
factor is 0.775). Placing `N·W` at random, the sets `{i : φ(w) ∈ S_i}` would be a `c`-fold colouring of `W` with
`k < χ_c(W)` colours unless some edge of `W` has both ends in some `S_i`, which happens with probability at most
`k|E(W)|ε⁻²√(R/(πN))`, below 1 for large `N`. For `c = 6`, `N₀` is not explicit, because OpenAI's proof gives no
finite graph.

## Uniform bounds for two distances (versions 6 and 7)

**Theorem 5 (version 7).** For every real algebraic `d ≥ 1`: `χ(ℝ², {1, d}) ≥ k` for `d ≥ d_k`, where
`d_7 = 5.82`, `d_8 = 10.5`, `d_9 = 18.9`, `d_10 = 35.7`, `d_11 = 73.6`, `d_12 = 181`,
`d_13 = 655`, `d_14 = 25 700`; moreover `χ ≥ 7` for `1.72 ≤ d ≤ 2.04` and `χ ≥ 8` for `1.82 ≤ d ≤ 1.89`. The
same holds with `1/d` in place of `d` (scale the plane), and, for measurable colourings, for every real `d`. Version 6
had, without configurations, `d_k = 6.52, 12.05, 24.3, 52.9, 140.2, 601, 54 300` for `k = 7, …, 13` and the ranges
`[1.75, 1.93]` and `[1.83, 1.855]`; version 7 adds the configuration inequalities of unit equilateral triangles (all
certificates) and of Moser spindles (the tails). Theorem 5 now contains every finite bound of Corollary 4 (the pairs
`{1, √3}`, `{1, (√6+√2)/2}`, `{1, 2}`, `{1, 11/6}`, `{1, 40}`, `{1, 1000}`) and raises `{1, 40}`, `{1, 100}` and
`{1, 1000}` to 10, 11 and 13 colours.

- **Why it is new.** For arbitrary colourings no pair was known to need 7 colours before this paper
  (Chybowska-Sokół, Junosza-Szaniawski and Węsek record it as open), and our earlier bounds above 6 were for isolated
  values of `d`; here a whole half-line of algebraic `d` needs 7 colours, and 14 from `25 700` on. For measurable
  colourings, Oliveira Filho and Vallentin had 8 colours for one pair and Bukh's theorem gives 14 or more for all `d`
  beyond an unspecified threshold, so there only the explicit thresholds and ranges are new.
- **The configurations (Proposition 17).** For a finite set `W` whose largest subset with no two points at a
  distance in `{1, d}` has `α(W)` points, and a point `p`, a measurable set of density `δ` avoiding `{1, d}` satisfies
  `Σ_{w∈W} f(|w − p|) ≤ α(W)δ` for its autocorrelation `f` (the inequality of Keleti, Matolcsi, Oliveira Filho and
  Ruzsa); in the dual this adds `y Σ_{w∈W} J₀(|w − p|t)`, `y ≥ 0`, to `g`, `α(W)y` to the numerator of the bound and
  `|W|y` to `s`. Two kinds of `W` are used, neither depending on `d`. *Unit triangles*: `W = {(0,0), (1,0),
  (1/2, √3/2)}` (`α = 1`) and `p_o` its centroid moved by a rational `o` along its axis, at distances `r₁(o)` (twice)
  and `r₃(o)`, `r₁(o)² = 1/3 + o/√3 + o²`, `r₃(o) = |1/√3 − o|`. *Moser spindles*:
  `M = {0, 1, ω, 1 + ω, ρ, ρω, ρ(1 + ω)}` in complex notation, `ω = e^{iπ/3}`, `ρ = (5 + i√11)/6` (`α = 2`: among any
  three of its points two are at distance 1), at a point `q` with rational coordinates. The linear program is offered
  `o = i/25`, `|i| ≤ 40`, and `q = (i/10, j/10)`, `i, j` odd from −7 to 21; the tails use 1 to 3 triangle positions
  and 1 to 3 spindle points.
- **Lemma 18 (one certificate for all `d ≥ D`).** With `μ = 0.4028 ≥ −min J₀` and `e(y) = min(μ, √(2/(πy)))`,
  `J₀(x) ≥ −e(y)` for every `x ≥ y > 0`. So if `z2 ≥ 0`, `s > 0` and `h(t) = S(t) − z2 e(Dt) ≥ 0` for `t > 0`, where
  `S(t) = z0 + z1 J₀(t) + Σ_c y_c (2J₀(r₁(o_c)t) + J₀(r₃(o_c)t)) + Σ_m u_m Σ_{w∈M} J₀(|w − q_m|t)`, then
  `g_d(t) = S(t) + z2 J₀(dt) ≥ h(t) ≥ 0` for every `d ≥ D`, and `m_{1,d}(ℝ²) ≤ (z0 + Σy + 2Σu)/s` with
  `s = z0 + z1 + z2 + 3Σy + 7Σu`. Checking `h ≥ 0` is a problem in one variable.
- **Interval certificates.** For `[1.72, 2.04]` (k = 7, 51 intervals) and `[1.82, 1.89]` (k = 8, 14 intervals)
  one function per interval of `d`, with triangles: on a box `[a, b] × [t0, t1]`, `dt ∈ [a t0, b t1]`, and `J₀` on an
  interval is enclosed by Taylor at the midpoint (`|J₀''| ≤ 1`), by `[−μ, 1]` and by the envelope.
- **How close to the method's limit.** In floating point, the program of Lemma 18 with triangles and spindles gives
  `k = 7, …, 14` colours from about `d = 5.797, 10.42, 18.81, 35.53, 73.36, 179.5, 652.6 and 25 590` on; with triangles alone, 7 to 13 from about 5.800, 10.67,
  20.04, 40.04, 90.24, 261.7 and 1544 on and never 14; without configurations, from about 6.503, 12.03, 23.87, 52.85,
  139.9, 599.2 and, for 13, from somewhere between 50 000 and 54 300. Without configurations its value tends to
  `μ₀²/(1 + μ₀)² = 1/12.13…` (`μ₀ = −min J₀ = 0.40275…`), the square of its value 0.287119 for one distance: for
  `d ≥ 20` it is above `1/12.2` (two constraints of the program, at `t = j_{1,1}/d` and near `t = j_{1,1}`), and the
  limit is attained by an explicit family of certificates (the fifth referee's). With configurations, the same
  construction with a one-distance certificate of value `v` gives bounds that tend to `v·μ₀/(1 + μ₀)`: with our
  triangles `v = 0.271823` (limit `1/12.81`; floating-point value about `1/12.80` at `D = 10⁸`), with triangles and
  spindles `v = 0.262476` (limit `1/13.27`; about `1/13.20` at `D = 10⁷`). So these configurations give 14 colours;
  in floating point they do not give 15, which would need `v < 0.2488` (Keleti, Matolcsi, Oliveira Filho and Ruzsa reach (numerically) 0.26305 with spindles
  alone and 0.258795 with further inequalities). For one `d` with triangles, the program gives 7 colours for
  `1 < d < 5.8` only for `d` between about 1.716 and 2.045, and again from about 5.55 on.
- **Programs.** `tail.py` and `uniform.py` (Arb, Taylor enclosures; they produce and re-verify) and
  `check_uniform.py` (separate code: point values and the Lipschitz bound `|J₁| ≤ 1/√2`, boxes in `t` and `d`, a
  floating-point counterexample search confirmed in Arb; radii recomputed from `o`, spindle distances from the points
  of `M`, built by complex multiplication). Both accept the ten certificates and reject eighteen corrupted ones (plus
  a slow nineteenth: a dip of `−10⁻⁵` that the floating-point scan misses). `tail.py k D --cliques --spindles`
  rebuilds the tail certificates byte for byte (the test checks `tail_k14.json`).
- **Figure.** `papers/forbidden-distances/uniform-bounds.pdf` (and `.png`): the certified staircase, the
  floating-point values of the program for one `d` with triangles and of the program of Lemma 18 with triangles and
  spindles, and the limit of the latter.

## Formal verification of the transfer

Theorem A in the form used by Theorem C — a proper `k`-colouring avoiding a set `D` of positive real algebraic
distances gives a measurable `k`-colouring whose classes weakly avoid `D`
(`OAI.Ext.proper_to_borel_weak_distances`), and `k` pairwise disjoint measurable sets avoiding `D` with conull union
(`OAI.Ext.proper_to_measurable_classes`) — is proved in Lean 4 on top of OpenAI's formalization
([`lean/external/openai-five/distances/`](../lean/external/openai-five/distances/README.md)): four new files (774
lines), only the standard axioms, built from source twice (37 modules) and replayed by `leanchecker`. Not
formalized: the converse of Theorem A, `χ ≥ 1/m_D`, Theorem B, the lemmas of the two-scale theorem and the density
bounds.

## What is new and what is not

Checked against the literature up to 10 October 2026 by the second referee and by us.

- **Dúcz and Varga**, *A unit-distance graph in the plane with independence ratio below 1/4*, arXiv:2606.28157
  (26 June 2026), already proved `χ_f(ℝ²) > 4.0007`, with a unit-distance graph `G₂₉` of geometric fractional
  chromatic number above 4.0007, and so disproved Conjecture 1 of Matolcsi, Ruzsa, Varga and Zsámboki
  (`χ_{f,0}(ℝ²) = 4`). Our draft v2 did not know this and claimed the disproof; the claim is withdrawn. What
  remains new: the exact value `1/m₁(ℝ²)` (so the bound `4.0485`), the equality `χ_{f,0} = χ_f` that Matolcsi et al.
  left open, and the first half of Dúcz and Varga's own Conjecture 1 (`f(n)/n = m₁(ℝ²) + o(1)`; their second half,
  `m₁ = 0.2293647…` from Croft's construction, stays open).
- **Davies**, *Measurable obstructions for unmeasurable colourings*, arXiv:2610.12301 (8 October 2026), proves by a
  different method that in `ℝ^d` with `d ≥ 4`, for finite `D` with rational squares, the chromatic number equals the
  measurable one and the independence ratio equals the measurable one; with Bukh's theorem he gets
  `sup_{|D|=k} χ_D(ℝ^d) ≥ m₁(ℝ^d)^{−k}` for `d ≥ 4`. He hopes to extend his method to `d = 2, 3`. Our results are for
  `d = 2` (and `d = 3` for the many-distance corollary), with arbitrary algebraic distances, through OpenAI's
  transfer.
- **Erdős's question.** P. Erdős, *On the combinatorial problems which I would most like to see solved*,
  Combinatorica 1 (1981), §I: "Does `L₂(r)` increase polynomially or exponentially? I expect that `L₂(r) < r^c` for
  some absolute constant `c`." It is [Erdős Problem #706](https://www.erdosproblems.com/706), not marked solved on
  10 October 2026; Croft, Falconer and Guy (Problem G11) ask the same. The best lower bound in the plane was of
  order `k√(log k)`. Exponential growth was known for measurable colourings (Bukh, 2008) and, since 8 October, for
  `d ≥ 4` (Davies). For the plane and for `ℝ³` it appears to be new, conditional on [OAI, Theorem 2.3].
- **Two distances.** The best lower bound for a pair `{1, d}` was 6: Exoo–Ismailescu for `d = 2`; De Neve, Vanden
  Kerchove, Colle, Tavernier and Pickavet (Amer. Math. Monthly 2025) for `√3` and `(√6+√2)/2`, cases claimed earlier
  in Polymath16. Chybowska-Sokół, Junosza-Szaniawski and Węsek (arXiv:2201.04499) record as open whether a pair can
  force 7. For *measurable* colourings, Oliveira Filho and Vallentin (2010, §4) had `m_{1, j_{1,2}/j_{1,1}} ≤ 0.141577`
  (so 8 colours) and, with Bukh's theorem, pairs needing 14 measurable colours (17 with the bound of Ambrus et al.).
  New: the bounds hold for arbitrary colourings, the distances are rational or algebraic, the certificates are
  rigorous, and the clique-strengthened bounds for `{1, 2}` and `{1, √3}` seem new even for measurable colourings.
- **Two scales.** The fourth referee searched for two-scale or lexicographic bounds of the form
  `χ(ℝⁿ, D₁ ∪ tD₂) ≥ χ(ℝⁿ, D₁)·(something)` and for any bound above 6 for a pair, for arbitrary or measurable
  colourings (Bukh 2008, Davies 2024 and 2026, Davies–McCarty–Pilipczuk, Naslund, Kula–Omar–Stockwell–Xie,
  Chybowska-Sokół et al., Parts, Exoo–Ismailescu, De Neve et al., Grytczuk et al., Bellitto et al.) and found none;
  Bukh's Theorem 1 is the measurable density version `m(t₁D₁ ∪ t₂D₂) → m(D₁)m(D₂)`.
- **Limits.** No finite graph is given. Dúcz and Varga's graph can in principle be made explicit but is
  astronomically large; an explicit `{1, 2}`-graph with no proper 6-colouring remains to be found. The distances
  must be algebraic: that is a limitation of the method, and probably of the statement, since Bukh conjectured that
  `χ(ℝ^d, D)` is bounded independently of `D` for algebraically independent `D`, which would make
  `χ(ℝ², D) ≥ 1/m_D` fail for such `D`.

## The certificates

Each bound is a dual solution of the linear program of Oliveira Filho and Vallentin (J. Eur. Math. Soc. 2010,
Theorem 1.1): if `g(t) = z₀ + Σ_k z_k J₀(d_k t) + Σ_c y_c Σ_{w∈W_c} J₀(|w − p_c| t) ≥ 0` for all `t ≥ 0`, with
`y_c ≥ 0`, then `m_D(ℝ²) ≤ (z₀ + Σ_c y_c)/(z₀ + Σ_k z_k + Σ_c y_c |W_c|)`. The clique terms come from point sets `W_c`
with all pairwise distances in `D`: the sets `{x ∈ A : x + w − p ∈ A}`, `w ∈ W_c`, are disjoint subsets of `A`
(Proposition 17 of the paper).

| `D` | clique sets | certified bound | `T` | gives |
|---|---|---|---|---|
| `{1, 2}` | 6 unit triangles, 2 triangles of side 2 | 0.159222 < 1/6 | 1297 | `χ ≥ 7` |
| `{1, √3}` | a unit triangle, a rhombus | 0.162866 < 1/6 | 757 | `χ ≥ 7` |
| `{1, (√6+√2)/2}` | none | 0.166018 < 1/6 | 107 | `χ ≥ 7` |
| `{1, 10}` | none | 0.149599 < 1/6 | 36 | `χ ≥ 7` |
| `{1, 11/6}` | none | 0.140995 < 1/7 | 171 | `χ ≥ 8` |
| `{1, 46/25}` | none | 0.139185 < 1/7 | 175 | `χ ≥ 8` |
| `{1, 40}` | none | 0.115511 < 1/8 | 36 | `χ ≥ 9` |
| `{1, 100}` | none | 0.103302 < 1/9 | 37 | `χ ≥ 10` |
| `{1, 1000}` | none | 0.089085 < 1/11 | 37 | `χ ≥ 12` |
| `{1, 40, 1600}` | none | 0.046389 < 1/21 | 36 | `χ ≥ 22` |
| `{1, 10, 100, 1000}` | none | 0.039641 < 1/25 | 37 | `χ ≥ 26` |
| `{1, 100, 10⁴}` | none | 0.037168 < 1/26 | 37 | `χ ≥ 27` |
| `{1, 30, 900, 27000}` | none | 0.021264 < 1/47 | 36 | `χ ≥ 48` |

Three programs check every certificate ([`data/forbidden_distances/`](../data/forbidden_distances/README.md)):
- `certify_lp2.py` (Arb balls; distances and radii enclosed from their exact expressions; on `[0, T]` point values
  and the Lipschitz bound from `|J₁(x)| ≤ 1/√2`, DLMF 10.14.1; on `[T, ∞)` the bound `|J₀(x)| ≤ √(2/(πx))`, from the
  Sonin–Pólya argument: `u = √x J₀(x)` solves `u″ + (1 + 1/(4x²))u = 0`, so `u² + u′²/(1 + 1/(4x²))` is nondecreasing
  with limit `2/π`). It reproduces each stored certificate byte for byte.
- `check2.py`, separate code in 60-digit floating point.
- `check3.py`, rigorous and independent of both: exact clique geometry in `ℚ(√3)`, Arb enclosures over whole
  intervals of `t`, adapted from the second referee's checker.

## History

- **v1 (10 October, morning).** The transfer for distance sets and for fractional colourings, the corollaries from
  published density bounds, and the first certificates (`certify_lp.py`, no cliques).
- **First referee** (transfer theorems): correct with minor gaps. It found a false inequality in the headline
  (`4.0486 < 1/0.2470`; `1/0.2470 = 4.04858…`), a clash with Matolcsi et al.'s notation (our `χ_f` was their
  `χ_{f,0}`), a vacuous remark about uncountable `D`, the missing equality of densities over discs and over squares,
  and a one-line proof of the density lemma at radius `d`. All fixed in v2, which also proves `χ_f = χ_{f,0}` for
  finite `D` (Lemma 9 and the translates).
- **v2.** Lemma 9; the clique-strengthened program and `certify_lp2.py`; `{1, 2}` and `{1, √3}` at 7.
  A floating-point run for `{1, 2}` once suggested a bound of 0.15; its `g` was nonnegative only on the sampled
  range of `t` (with `z₀ < 0` it is negative at infinity), and the certification rejected it. The float program now
  requires `z₀ ≥ 0` and samples `t` up to 400; the certified bound is 0.159222.
- **Second referee** (corollaries, certificates, literature): no mathematical error; all certificates valid, checked
  with two programs of its own; the priority claims wrong (Dúcz–Varga, Davies); Erdős's question to be quoted and
  stated explicitly; measurable precedents to be credited; `{1, 100}` gives 10; code and citation fixes (the bound
  `|J_α| ≤ 1/√2` holds for `α ≥ 1`, not for every `α > 0`; a clique check that used numerical identification; radii
  through decimal evaluation). All applied in v3.
- **v3 (10 October, afternoon).** The paper with the third checker.
- **Third referee** (version 3): no mathematical error; the new material (the limit `f(n)/n → m₁`, the supremum
  identity, the integer rounding, the Erdős bounds) checked line by line against the sources; all certificates
  accepted by its own checker, and broken ones rejected. Findings: Croft's density is `0.2293647…`, so `m₁ ≥ 0.229365`
  was rounded up (now `0.229364`); `f(n) ≥ m₁n` is classical and is now credited; the title now says "with algebraic
  forbidden distances"; a wrong bibliography entry (Parts); the paraphrase of Oliveira Filho and Vallentin §4;
  scripts that relied on `assert` and exited 0 when the bound failed (fixed: explicit checks, non-zero exit, and tests
  that broken certificates are rejected).
- **v4 (10 October, afternoon).** The third referee's findings, and certified bounds for three and four distances.
- **v5 (10 October, afternoon).** The two-scale theorem (first with fractional weights, then in the `c`-fold form),
  the explicit pairs, the Lean proof of the transfer, `{1, 1000}` at 12 and `{1, 30, 900, 27000}` at 48.
- **Fourth referee** (version 5, with numerical checks of every inequality of the smoothing lemma, a certificate
  checker of a different design, an exact recomputation of de Grey's graph and a literature search): no mathematical
  error; the true supremum in the smoothing lemma is 0.775, not 1, so the lemma has slack; nothing in the literature
  states a two-scale bound or any bound above 6 for a pair. Applied: the first 6-colour pair `{1, (1+√5)/2}` is now
  credited to Owings, Tetiva and Huddleston (Amer. Math. Monthly 2008; as reported by Chybowska-Sokół et al. and De
  Neve et al.), and Parts is described without saying which pair his small graph is for; the Lean paragraph names the
  declarations, versions and commit and lists what is not formalized; Geller–Stahl is cited; the window lemma no
  longer assumes disjoint sets; the threshold `6.76·10²⁴` is `6.75·10²⁴`; the Bellitto–Pêcher–Sédillot bound is
  cited as their weighting (their §4.4 heading has the inequality reversed); the headline 25 is noted to need only a
  finite unit-distance graph with `χ_f > 4` (Dúcz–Varga's), so it does not depend on Theorem B or on Ambrus et al.;
  Corollaries 5 and 8 are stated in every dimension. For `{1, (1+√5)/2}` our program with pentagon and triangle
  cliques gives only 0.1769 > 1/6.
- **v6 (10 October, evening).** Theorem 5 (uniform bounds for two distances), Lemma 18 (one certificate for all
  `d ≥ D`), the remark on the limit of the two-distance program, Figure 1, and the folder
  `data/forbidden_distances/uniform/` with ten certificates and two programs. The fifth referee found no mathematical
  error; applied: the measurable clause is argued, the remark states its hypotheses (`d ≥ 20`) and proves that the
  limit is attained, "each step by a single certificate" is restricted to `d ≥ 12.3` (the step at 12.05 uses interval
  certificates), the floating-point crossings are stated as "from about … on", the 13-colour crossing (between 50 000
  and 54 300) is given, and the measurable novelty is limited to the explicit thresholds.
- **v7 (10 October, night).** Unit triangles and Moser spindles in the uniform certificates (Proposition 17 with
  configurations), 14 colours from `d = 25 700` on. The sixth referee found no mathematical error; applied: the value
  of a certificate in the remark is `z₀ + Σ α_c y_c`, the conditions the construction needs on `Q` are stated, "14
  colours follow for all large `d`" now rests on a rigorous certificate (`S/S(0)` of `tail_k14.json`, value 0.267887,
  below the exact threshold 0.267894) instead of the floating-point 0.262476, Keleti et al.'s 0.26305 is marked
  numerical, Corollary 4 points to Theorem 5 for `{1, 40}`, `{1, 100}`, `{1, 1000}`, the figure caption no longer
  calls 13.27 a limit of the program, and the single-`d` window includes "from about 5.55 on".
