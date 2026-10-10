# Uniform certificates for two forbidden distances

These certificates prove Theorem 5 of [`papers/forbidden-distances/`](../../../papers/forbidden-distances/README.md):
for every real algebraic `d ≥ 1`,

| colours | for every `d` with |
|---|---|
| `χ(ℝ², {1, d}) ≥ 7` | `1.72 ≤ d ≤ 2.04` or `d ≥ 5.82` |
| `≥ 8` | `1.82 ≤ d ≤ 1.89` or `d ≥ 10.5` |
| `≥ 9` | `d ≥ 18.9` |
| `≥ 10` | `d ≥ 35.7` |
| `≥ 11` | `d ≥ 73.6` |
| `≥ 12` | `d ≥ 181` |
| `≥ 13` | `d ≥ 655` |
| `≥ 14` | `d ≥ 25 700` |

and the same for `1/d`, since `χ(ℝ², {1, 1/d}) = χ(ℝ², {1, d})`. For measurable colourings the same holds for every
real `d`. For arbitrary colourings no pair of distances was known to need 7 colours before the paper (the bound was 6;
Chybowska-Sokół, Junosza-Szaniawski and Węsek record the question as open), and the pairs certified in
[`../certificates/`](../README.md) are isolated values of `d`. For measurable colourings, Oliveira Filho and Vallentin
had 8 colours for one pair, and Bukh's theorem gives 14 or more for all `d` beyond a threshold it does not specify; there
only the explicit thresholds and ranges are new.

Each certificate is a dual solution of the linear program of Oliveira Filho and Vallentin for the two distances
`1` and `d`, in exact fractions, with configuration inequalities (version 7 of the paper; version 6 had none, and the
thresholds `6.52, 12.05, 24.3, 52.9, 140.2, 601, 54 300` for 7 to 13 colours). Two configurations are used:

- **Unit equilateral triangles** (all certificates). For a rational `o` let `p_o` be the centroid of the triangle
  `W = {(0,0), (1,0), (1/2, √3/2)}` moved by `o` along its axis; its distances to the vertices are `r₁(o)` (twice) and
  `r₃(o)`, with `r₁(o)² = 1/3 + o/√3 + o²` and `r₃(o) = |1/√3 − o|`. A set avoiding the distance 1 contains at most
  one point of `W`.
- **Moser spindles** (tail certificates). `M = {0, 1, ω, 1 + ω, ρ, ρω, ρ(1 + ω)}` in complex notation, with
  `ω = e^{iπ/3}` and `ρ = (5 + i√11)/6`: two rhombi made of unit triangles, turned about `0` so that their far vertices
  `1 + ω` and `ρ(1 + ω)` are at distance 1. Among any three of its points two are at distance 1, so a set avoiding the
  distance 1 contains at most two of them. A spindle is used at a point `q` with rational coordinates, through the
  seven distances `|w − q|`, `w ∈ M`.

Put `S(t) = z0 + z1 J₀(t) + Σ_c y_c (2 J₀(r₁(o_c) t) + J₀(r₃(o_c) t)) + Σ_m u_m Σ_{w∈M} J₀(|w − q_m| t)` with
`y_c, u_m ≥ 0`. Then `g_d(t) = S(t) + z2 J₀(dt) ≥ 0` for all `t ≥ 0` gives
`m_{1,d}(ℝ²) ≤ (z0 + Σ y_c + 2 Σ u_m)/s` with `s = z0 + z1 + z2 + 3 Σ y_c + 7 Σ u_m` (Proposition 17 of the paper,
the configuration inequalities of Keleti, Matolcsi, Oliveira Filho and Ruzsa), and a bound below `1/(k−1)` gives
`χ(ℝ², {1, d}) ≥ k` by the transfer theorem. What is new is that one certificate covers many `d` at once.

![Certified lower bounds for χ(ℝ², {1, d})](../../../papers/forbidden-distances/uniform-bounds.png)

## The certificates

- **`certificates/tail_k7.json` … `tail_k14.json`: one certificate for all `d ≥ D`** (Lemma 18 of the paper). With
  `μ = 0.4028 ≥ −min J₀` and `e(y) = min(μ, √(2/(πy)))`, we have `J₀(x) ≥ −e(y)` for every `x ≥ y > 0`. So if
  `z2 ≥ 0` and `h(t) = S(t) − z2 e(Dt) ≥ 0` for all `t > 0`, then `g_d ≥ h ≥ 0` for every `d ≥ D`. Fields:
  `k`, `D`, `z0`, `z1`, `z2` (fractions), `cliques` (a list of `{oy, y}`: the position `o_c` and the weight `y_c`),
  `spindles` (a list of `{qx, qy, u}`: the point `q_m` and the weight `u_m`), `bound`, and `T`, `boxes` (the end of
  the box covering and the number of boxes of `tail.py`).

  | file | `D` | triangles | spindles | bound for all `d ≥ D` | proves for `d ≥ D` |
  |---|---|---|---|---|---|
  | `tail_k7.json` | 5.82 | 3 | 1 | 0.1665379 < 1/6 | `χ ≥ 7` |
  | `tail_k8.json` | 10.5 | 3 | 2 | 0.1426634 < 1/7 | `χ ≥ 8` |
  | `tail_k9.json` | 18.9 | 3 | 2 | 0.1249086 < 1/8 | `χ ≥ 9` |
  | `tail_k10.json` | 35.7 | 3 | 2 | 0.1110473 < 1/9 | `χ ≥ 10` |
  | `tail_k11.json` | 73.6 | 2 | 2 | 0.0999702 < 1/10 | `χ ≥ 11` |
  | `tail_k12.json` | 181 | 2 | 2 | 0.0908619 < 1/11 | `χ ≥ 12` |
  | `tail_k13.json` | 655 | 1 | 3 | 0.0833234 < 1/12 | `χ ≥ 13` |
  | `tail_k14.json` | 25 700 | 1 | 3 | 0.0769214 < 1/13 | `χ ≥ 14` |

- **`certificates/interval_k*_a_b.json`: certificates for intervals of `d`.** A list of contiguous intervals
  `[a, b]` covering the range of the file, each with `z0, z1, z2` and triangles such that `g_d ≥ 0` for every `d` in
  `[a, b]`: `interval_k7_1.72_2.04.json` (51 intervals, `χ ≥ 7`) and `interval_k8_1.82_1.89.json` (14,
  `χ ≥ 8`). They use triangles only.

## The programs

Run them in this folder.
- **`tail.py check certificates/tail_k*.json`** and **`uniform.py check certificates/interval_k*.json`** re-verify
  the stored certificates; **`tail.py k D --cliques --spindles`** and **`uniform.py k a b step --cliques --T=200`**
  produce them (the linear program is offered the triangle positions `o = i/25`, `−40 ≤ i ≤ 40`, and the spindle
  points `q = (i/10, j/10)` with `i, j` odd from −7 to 21). They choose the coefficients by a floating-point linear
  program (`lpcp.py` and `tail.py`'s own; nothing these compute is trusted), round them to fractions, raise `z0` by a
  third of the slack, and prove `h ≥ 0`, respectively `g_d ≥ 0` on `[a, b] × [0, ∞)`, in Arb ball arithmetic: they
  cover `[0, T]` by `t`-boxes, enclose `J₀` on an interval `[m − r, m + r]` by `J₀(m) ± (|J₁(m)| r + r²/2)` (since
  `|J₀''| ≤ 1`), by `[−μ, 1]` and by `±√(2/(πx₀))` beyond `x₀ > 0`, bisect until each box has a nonnegative lower
  bound, and use `|J₀(x)| ≤ √(2/(πx))` on `[T, ∞)`. `J₀ ≥ −0.4028` is proved by `check_j0min` in `uniform.py` (4000
  boxes on `[0, 4]`, the envelope beyond). The radii `r₁(o)`, `r₃(o)` and the spindle distances (from the exact
  coordinates of `M`) are enclosed by rational intervals from Arb. `tail.py 14 25700 --cliques --spindles`
  reproduces `tail_k14.json` byte for byte.
- **`check_uniform.py certificates/*.json`** is a second program, written separately. It recomputes each target
  `1/(k−1)` from `k`, each radius from `o` and each spindle distance from the points of `M`, which it builds by complex
  multiplication; it checks that the intervals of a file are contiguous and cover its range, looks for a
  counterexample in floating point and confirms it in Arb before rejecting, encloses `J₀` by point values and the
  Lipschitz bound `|J₁| ≤ 1/√2` instead of the Taylor bound, bisects boxes in both `t` and `d` for intervals, sweeps
  `h` with steps chosen from its value for tails, and checks `J₀ ≥ −0.4028` itself. It accepts all ten certificates.

Both kinds of program exit with a non-zero code when a check fails. They reject corrupted certificates (see the
tests): a lowered `z0`, a scaled `z2`, a negative or a missing triangle or spindle weight, a moved triangle or spindle,
a threshold `D` moved below the certified one, a raised `k`, a gap between intervals, a reversed interval, an interval
extended below its range, and a certificate whose `g_d` has minimum `−10⁻⁵` at a point that the floating-point scan
misses.

**`figure/`**: `uniform_fig.tex` draws the figure of the paper (`pdflatex uniform_fig.tex`) from `lpcurve.dat` and
`tailcurve.dat`, which `scan.py` computes from the floating-point values in `lowscan.tsv` (the linear program for one
`d`, with triangles, for `d` from 1.01 to 6.6) and `tailscan.tsv` (the program of Lemma 18, with triangles and
spindles, for `D` from 5 to 2·10⁵, and at 10⁶ and 10⁷). These values are not proofs. In floating point, the program
of Lemma 18 gives `k = 7, …, 14` colours from about `d = 5.797, 10.42, 18.81, 35.53, 73.36, 179.5, 652.6 and 25 590` on. With triangles alone it gives 7, …, 13
colours from about 5.800, 10.67, 20.04, 40.04, 90.24, 261.7 and 1544 on, and never 14; without triangles and spindles,
from about 6.503, 12.03, 23.87, 52.85, 139.9, 599.2 and, for 13, from between 50 000 and 54 300. For one `d` with
triangles, the program gives 7 colours for `d < 5.8` only between about 1.716 and 2.045. Without configurations the
value of the program of Lemma 18 tends to `μ₀²/(1 + μ₀)² = 1/12.13…` (`μ₀ = −min J₀`), the square of its value
0.287119 for one distance, as `D → ∞` (proved in the remark after Lemma 18). With configurations whose one-distance
program has value `v`, an explicit construction gives bounds that tend to `v·μ₀/(1 + μ₀)`: `v = 0.271823` with our
triangles (`1/12.81`; the floating-point value of the program is about `1/12.80` at `D = 10⁸`), and `v = 0.262476` with
triangles and spindles (`1/13.27`; the floating-point value is about `1/13.20` at `D = 10⁷`). So these configurations
give 14 colours; in floating point they do not give 15, which would need `v < 0.2488`.

**Tests.** From the root of the repository, `python3 -m pytest -q tests/test_forbidden_distances.py -k uniform` runs
both programs on every certificate and checks that both reject the corrupted ones (a few minutes); the corruption
that the floating-point scan misses is marked slow.
