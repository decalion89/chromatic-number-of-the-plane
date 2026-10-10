# Inclusion–exclusion certificates: 15, 16 and 17 colours for two distances

These certificates prove the steps `k = 15, 16, 17` of Theorem 5 of
[`papers/forbidden-distances/`](../../../../papers/forbidden-distances/README.md): for every real algebraic `d`,

| colours | for every `d` with | certificate | bound for `m_{1,d}(ℝ²)` |
|---|---|---|---|
| `χ(ℝ², {1, d}) ≥ 15` | `d ≥ 5·10⁵` | `ie_k15.json` (Proposition 19) | `0.0713614 < 1/14` |
| `≥ 16` | `d ≥ 7·10⁷` | `ie_k16.json` + `ie_B.json` (Proposition 21) | `0.0657685 < 1/15` |
| `≥ 17` | `d ≥ 8·10⁸` | `ie_k17.json` + `ie_B.json` (Proposition 21) | `0.0623330 < 1/16` |

and the same for `1/d` and, for measurable colourings, for every real `d`. Before, the explicit thresholds were 14 colours
from `d = 25 700` (the tails in [`../certificates/`](../README.md)) and 18 from `d = 10²⁵` (Theorem 8 of the paper,
with de Grey's graph); Bukh's theorem gives 17 for all `d` beyond a threshold it does not specify. Unlike the bound
with de Grey's graph, these are density bounds, so they also give `χ_f(ℝ², {1, d}) > 14.01, 15.20, 16.04` for these `d`.

## The configuration

`x23.py` holds the 23 points `X23` of Ambrus, Csiszárik, Matolcsi, Varga and Zsámboki
([arXiv:2207.14179](https://arxiv.org/abs/2207.14179), *The density of planar sets avoiding unit distances*, Math.
Program. 2024), copied from their verification notebook, with exact coordinates in `ℚ(√3, √11)` (`kfield.py`). Its
unit-distance graph has 47 edges and 13 552 independent sets (including `∅`); the 253 distances take 27 values. They
used `X23` to prove `m₁(ℝ²) ≤ 0.2470`: a random congruent copy of `X23` meets a set `A` avoiding the distance 1 in an
independent set `S`, and the probabilities of the patterns `S` satisfy linear equations ("inclusion–exclusion")
whose coefficients are values of the autocorrelation of `A`. `congruences.py` writes `cong23.json`: the independent
sets of size at least 2 split into 11 226 congruence classes (compared exactly), which give 2 302 constraints
`P(I ⊆ S) = P(J ⊆ S)`; they imply every congruence constraint between independent sets.

## Proposition 19: 15 colours

At the scale of 1 we use their program; at the scale of `d` only `f(d) = 0`, through Lemma 18
(`J₀(x) ≥ −e(y)` for `x ≥ y`, `e(y) = min(0.4028, √(2/(πy)))`), and `κ(0) = δ²`, through its tangent at `δ₀`. A
certificate gives numbers `ψ₀, ψ_T, λ_i, u_ij, w_c, η, τ` with

- `Φ(t) = ψ₀J₀(t) − Σλ_i − Σu_ij J₀(r_ij t) − η e(Dt) + 2δ₀τ ≥ 1` for `t > 0`,
- `Φ₀ = ψ₀ − Σλ − Σu + η + τ(2δ₀ − 1) ≥ 1`,
- `Q(S) = ψ_T + Σ_{i∈S} λ_i + Σ_{i<j∈S} u_ij + Σ_c w_c([I_c ⊆ S] − [J_c ⊆ S]) ≥ 0` for every independent `S`,

where the sums over `u_ij` run over the non-unit pairs, and then `m_{1,d}(ℝ²) ≤ ψ_T + τδ₀²` for every `d ≥ D`. In
floating point the program's value tends to about `0.0711 = 1/14.06` as `D → ∞` (about `0.24697 · 0.287`, the value of
the program of Ambrus et al. times the one-distance value `μ/(1 + μ)` of Lemma 18), so this method stops at 15 colours;
this is an observation, explained heuristically in the paper, not a theorem.

## Lemma 20 and Proposition 21: 16 and 17 colours

The new ingredient is the *hit set* `B_g = ⋃_i (A − g x_i)`, the set of points `y` such that `y + gX23` meets `A`.
Its density `q` is the probability that a random copy of `X23` meets `A`, and

- at the scale of 1, conditioning the inclusion–exclusion program on `S ≠ ∅` gives `δ ≲ 0.247·q`;
- at the scale of `d`, `B_g` almost avoids the distance `d`: since `1_{B_g}(y) ≤ Σ_i 1_A(y + g x_i)`, `f(d) = 0` and
  Lemma 18 bound its autocorrelation at `d` by `ε ≤ Σ_{t>0} κ(t) e(Dt)(529 − H(t))`, `H(t) = Σ_{i,j} J₀(r_ij t)`;
- so a second certificate, for sets that avoid the distance 1 up to `ε` (`ie_B.json`: `Φ_B ≥ 1` on `[0, ∞)` and
  `Q_B(R) ≥ 0` for **all** `2²³` subsets `R`, with a penalty `ϖ_ij` for each unit pair `{i, j}` inside `R`), gives
  `q ≤ ψ_T^B + L_B ε` with `L_B = max(0, ψ₀^B + Σ ϖ_ij)` (Lemma 20). In `ie_B.json`, `ψ_T^B = 0.2470097` and
  `L_B = 17.764` (rounded up); one penalty per unit pair, chosen by cutting planes over the subsets, gives a smaller
  `L_B` than a uniform penalty (29.87), which lowers the thresholds about threefold.

The A-side certificate has one more number `ζ ≥ 0` (with `ψ_T ≥ ζ`, and `Φ(t) − ζL_B e(Dt)(529 − H(t)) ≥ 1` in place
of `Φ(t) ≥ 1`), and then `m_{1,d}(ℝ²) ≤ ψ_T + τδ₀² − ζ(1 − ψ_T^B)` for every `d ≥ D`. As `D → ∞` the bound tends to
about `0.247² = 1/16.4`: `ie_B.json`, scaled by `1 + O(D^{-1/2})`, also serves on the A side with `ζ = ψ_T`, so
`lim sup m_{1,d} ≤ (ψ_T^B)²`. This is an effective form of the upper bound in Bukh's theorem `m_{1,d} → m₁(ℝ²)²`,
with `m₁(ℝ²)` replaced by the inclusion–exclusion bound `ψ_T^B` and with usable thresholds.

## Files and programs

| file | what it is |
|---|---|
| `x23.py`, `kfield.py` | the 23 points, exact arithmetic in `ℚ(√3, √11)` |
| `congruences.py` → `cong23.json` | the congruence constraints |
| `ie_lp.py` | the floating-point linear programs (A side, optionally with the hit-set row; B side with `--B`) |
| `ie_cert.py` | rounds a dual to a rational certificate (B side: exact integer check over all `2²³` subsets) |
| `check_ie.py` | **the rigorous checker** (exact arithmetic + Arb; uses `../uniform.py` for the `J₀` enclosures) |
| `check_ie_indep.py` | **a second checker, written separately**: its own interval arithmetic on integers and its own `J₀` series |
| `certificates/` | `ie_k15.json`, `ie_k16.json`, `ie_k17.json`, `ie_B.json` |

To re-verify: `python3 check_ie.py certificates/ie_k15.json 15`, `… ie_k16.json 16` and `… ie_k17.json 17` (about
half a minute for the first and five to ten minutes for the other two, which also check `ie_B.json` on all `2²³`
subsets; about 1 GB of memory), and the same with `check_ie_indep.py` (under a minute each; it needs only the
standard library and numpy). Every point index must lie in `1..23`: both checkers reject a certificate with an index
outside, which the seventh referee showed could otherwise fake a bound. `tests/test_inclusion_exclusion.py` runs both
checkers on the four certificates and on nine corrupted copies.
To rebuild: `python3 ie_lp.py --B --out=dualB.npz && python3 ie_cert.py dualB.npz ie_B.json --B`, then
`python3 ie_lp.py D d0 [--hit=ie_B.json] --out=dual.npz && python3 ie_cert.py dual.npz out.json [--hit=ie_B.json]`
(each linear program takes 10–25 minutes and about 1 GB).
