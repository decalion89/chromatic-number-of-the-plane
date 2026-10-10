# Measurable density bounds hold for all colourings of the plane

**Status (10 October 2026).** Proved on paper, conditional only on OpenAI's rigidity theorem (their Theorem 2.3,
formalized in Lean by OpenAI); the density bounds are computer proofs checked by three programs.

- **Referees.** Two AI-assisted referees inside the project. The first read the transfer theorems (correct, minor
  gaps, fixed); the second read the corollaries and the certificates with two checkers of its own (no mathematical
  error, every certificate valid; the priority claims were wrong and are corrected). A third reading of the final
  version is under way.
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
  `χ_f ≤ 1/m_D` is elementary: translates of a periodic set of density close to `m_D` (Lemma 9 of the paper).

Theorems A and B are [OAI, Section 4] with the unit circle replaced by the circles of radii `d ∈ D` (countably many,
since `D ⊆ F`) and, for B, with one indicator function in place of the `k` colour classes; the paper gives the
proofs.

## Consequences

**One distance.**
- `χ_{f,0}(ℝ²) = χ_f(ℝ²) = 1/m₁(ℝ²)`, so `4.0485 < χ_f(ℝ²) < 4.36` (Ambrus, Csiszárik, Matolcsi, Varga and
  Zsámboki: `m₁ ≤ 0.2470`, and their constants give `m₁ < 0.2469975`, so `χ_f > 4.0486`; Croft:
  `m₁ ≥ 0.229365`).
- `f(n)`, the least independence number of the unit-distance graph of `n` points in the plane, satisfies
  `f(n) ≥ m₁ n` for every `n` and `f(n)/n → m₁(ℝ²)`. The limit uses Theorem 2 of Matolcsi, Ruzsa, Varga and
  Zsámboki (the finitary fractional chromatic number equals the reciprocal of the finitary independence ratio) and
  Fekete's lemma.

**Two distances** (certified bounds of the next section, then Theorem A):
`χ(ℝ², {1, 2}) ≥ 7`, `χ(ℝ², {1, √3}) ≥ 7`, `χ(ℝ², {1, (√6+√2)/2}) ≥ 7`, `χ(ℝ², {1, 10}) ≥ 7`,
`χ(ℝ², {1, 11/6}) ≥ 8`, `χ(ℝ², {1, 46/25}) ≥ 8`, `χ(ℝ², {1, 40}) ≥ 9`, `χ(ℝ², {1, 100}) ≥ 10`; and
`χ(ℝ², {1, d}) ≥ 17` for every large enough algebraic `d` (Bukh's theorem: `m_{1,d} → m₁² < 1/16`; the threshold
is not computed). Products of hexagonal colourings give `χ(ℝ², {1, d}) ≤ 49`.

**Many distances.** For every `k` there are `k` positive integers `D = {1, t, …, t^{k−1}}` with
`χ(ℝ³, D) ≥ χ(ℝ², D) ≥ m₁(ℝ²)^{−k} > 4.0485^k` (Bukh's theorem: `m_D → m₁^k` as `t → ∞`), and
`sup{χ_f(ℝ², D) : D ⊆ F, |D| = k} = m₁(ℝ²)^{−k}` (Bukh also proves `m_D ≥ m₁^k` for every `k` distances). With the
product colouring, Erdős's `L₂(r)` satisfies `4.0485^r < L₂(r) ≤ 7^r`.

## What is new and what is not

Checked against the literature up to 10 October 2026 by the second referee and by us.

- **Dúcz and Varga**, *A unit-distance graph in the plane with independence ratio below 1/4*, arXiv:2606.28157
  (26 June 2026), already proved `χ_f(ℝ²) > 4.0007`, with a unit-distance graph `G₂₉` of geometric fractional
  chromatic number above 4.0007, and so disproved Conjecture 1 of Matolcsi, Ruzsa, Varga and Zsámboki
  (`χ_{f,0}(ℝ²) = 4`). Our draft v2 did not know this and claimed the disproof; the claim is withdrawn. What
  remains new: the exact value `1/m₁(ℝ²)` (so the bound `4.0485`), the equality `χ_{f,0} = χ_f` that Matolcsi et al.
  left open, and the first half of Dúcz and Varga's own Conjecture 1 (`f(n)/n = m₁(ℝ²) + o(1)`; their second half,
  `m₁ = 0.229365…` from Croft's construction, stays open).
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
(Proposition 10 of the paper).

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
- **v3 (10 October, afternoon).** The paper as it stands, with the third checker.
