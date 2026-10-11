# The density of planar sets avoiding the distance 1

`m₁(ℝ²)` is the largest upper density of a measurable planar set with no two points at distance 1. Ambrus,
Csiszárik, Matolcsi, Varga and Zsámboki proved `m₁(ℝ²) ≤ 0.2470` (Math. Program. 207 (2024) 303–327,
arXiv:2207.14179) with an inclusion–exclusion linear program on the patterns in which a random congruent copy of a
23-point set `X₂₃` meets the set. This folder extends `X₂₃` by 8 points, all in the same field `ℚ(√3, √11)`,
and proves

> **`m₁(ℝ²) ≤ 0.24353`**, so `χ_f(ℝ²) ≥ 1/m₁(ℝ²) ≥ 4.1063`.

This improves the bound with a published proof. It does not improve the bound `0.2415` that Varga announced for the
same group in December 2023 (*Solving a conjecture of Erdős*, blog post of the Rényi Institute's AI group,
https://ai.renyi.hu/posts/solving-erdos-conjecture/), for which we found no paper, preprint or certificate.

The paper is [`papers/unit-distance-density/`](../../../papers/unit-distance-density/README.md).

| configuration | points | unit pairs | independent sets | relations | value (floating point) | certified bound |
|---|---|---|---|---|---|---|
| `X₂₃` (Ambrus et al.) | 23 | 47 | 13 552 | 2 122 | 0.246948 | 0.24698967 |
| first 25 points of `X₂₉` | 25 | 56 | 19 868 | 4 091 | 0.245583 | |
| `X₂₆` | 26 | 60 | 28 622 | 5 008 | 0.245169 | 0.24522384 |
| first 27 points of `X₂₉` | 27 | 66 | 33 898 | 6 044 | 0.244841 | |
| `X₂₉` | 29 | 76 | 55 610 | 10 012 | 0.244092 | 0.24414455 |
| `X₃₀` | 30 | 80 | 84 893 | 15 709 | 0.243743 | 0.24379682 |
| `X₃₁` | 31 | 84 | 118 421 | 23 863 | 0.243474 | **0.24352462** |

The value is that of the program with `κ` on the grid `t = 0, 0.1, …, 100` and `100.25, 100.5, …, 400` (the first
program of `m1cert2.py`, below); the certified bound is the `ψ_T` of the certificate, rounded up. `X₂₃` alone gives `0.24698967`, the bound of Ambrus et al.

## Files

- `configurations/x23.json`, `x26.json`, `x29.json`, `x30.json`, `x31.json`: the points, exactly. Each point is
  `[[a, b, c, d], [a', b', c', d']]` for `(a + b√3 + c√11 + d√33) + i(a' + b'√3 + c'√11 + d'√33)`, the rationals as
  strings. `x23.json` is Table 1 of Ambrus et al. in their order (the same points as
  [`../uniform/ie/x23.py`](../uniform/ie/x23.py)); each larger set starts with the smaller one.
- `certificates/m1_x23.json`, `m1_x26.json`, `m1_x29.json`, `m1_x30.json`, `m1_x31.json`: the certificates (format below).
- `check_m1.py CERT.json`: the rigorous check (exact arithmetic in `ℚ(√3, √11)`, Arb balls for `J₀`). It is
  self-contained (numpy and python-flint), cross-checks the number of independent sets, and does not use `assert`, so
  it also runs under `python3 -O`; run it as `python3 -I check_m1.py ...`.
- `check_m1_indep.py CERT.json` (with `m1_indep_core.py`): a second check, written separately with its own field
  arithmetic, its own integer interval arithmetic (power series and Hankel expansion of `J₀`) and its own sweep. It
  needs only Python and numpy.
- `m1cert2.py CONFIG.json OUT.json`: builds a certificate (floating point; nothing it computes is trusted).
- `ielp.py`: the program in floating point (HiGHS, or Clarabel with `--solver=clarabel`).
- `exactconf.py`: the exact unit-distance graph, the independent sets and the congruence relations of a configuration.
- `explore.py BASE.json OUT.jsonl WORKER NWORKERS`: the candidate points of Ambrus et al. (§7) for a configuration,
  each evaluated by the program; `extend.py BASE.json OUT.json re,im ...`: adds candidate points exactly.
- `logs/`: the output of both checkers on the five certificates.

## The certificate

A certificate is a JSON object with `side: "M1"`, the `points`, rationals `psi0` and `psiT`, `lam` (one value per
point), `u` (`[i, j, value]`, 1-based, only for pairs not at distance 1) and `wc` (`[I, J, value]`: two lists of
point indices, congruent in the given order). With `r_ij = |x_i − x_j|` it claims

- `Q(S) = psiT + Σ_{i∈S} lam_i + Σ_{i<j∈S} u_ij + Σ_wc value·([I ⊆ S] − [J ⊆ S]) ≥ 0` for every independent set `S` of
  the unit-distance graph of the points, the empty set included, and
- `Φ(t) = psi0·J₀(t) − Σ lam_i − Σ u_ij J₀(r_ij t) ≥ 1` for every `t ≥ 0`.

Then `m₁(ℝ²) ≤ psiT` (Proposition 3 of the paper, after Ambrus et al.): for a periodic set `A` of density `δ` avoiding
1, a random congruent copy of the points meets `A` in an independent set `S`; `E Q(S) ≥ 0` and the integral of `Φ`
against the spectral measure of `A` give `δ ≤ psiT`.

`check_m1.py` checks the unit-distance graph exactly, the ranges of all indices, that no `u_ij` sits on a unit pair
and that every relation is a congruence; `Q ≥ 0` on all independent sets in exact integers; and `Φ ≥ 1` by
bisection with Arb enclosures on `[0, T]` and `|J₀(x)| ≤ √(2/(πx))` beyond. Both checkers exit with a non-zero code
when a check fails.

## How the certificates were made

`m1cert2.py` solves the program in floating point with `κ` on `t = 0, 0.1, …, 100` and `100.25, 100.5, …, 400`
(without the relations between pairs, which repeat the pair equations), rounds the relation multipliers `w` to
multiples of `2⁻⁴⁰` and freezes them, then solves a small program in `(psi0, psiT, lam, u)` alone: minimize `psiT`
subject to `Q ≥ 0` on every independent set, `−Σ lam ≥ 1.03`, and `Φ ≥ 1 + 2·10⁻⁵` on a grid of step `0.005` up to
`t = 50` and `0.02` up to `1000` (rows added at the local minima of `Φ` until none is violated). It rounds the result
and raises `psiT` until `Q ≥ 0` holds exactly. Times on one core: about one minute for `X₂₃`, four for `X₂₆`,
seventeen for `X₂₉`, fifty-four for `X₃₀` and two hours for `X₃₁`. With `--reuse` it reads the first program's solution from `OUT.lp.npz` and redoes only the rest;
from the stored solutions it reproduces the certificates byte for byte.

## How the points were found

A candidate point is at distance 1 from at least two points of the configuration and at distance at least `0.1`
from all of them (the rule of Ambrus et al., §7). `X₂₃` has 208 candidates; 61 of them were evaluated one at a
time (`explore.py`). The six best of those at distance 1 from at least four points of `X₂₃` were added in the order
of their single gains: the new points of `X₂₆` are `x₂₄ = (3 − √33)/6 + i√3/3`, `x₂₅ = (11 − √33)/12 + i(√11 − √3)/12`
and `x₂₆ = 1/6 − i√11/6`; `X₂₉` adds `x₂₇ = x₂₄ + 1`, `x₂₈ = 4/3 + i(3√3 + √11)/6` and `x₂₉ = (2 − √33)/6 + i(2√3 + √11)/6`.
`X₃₀` adds the best of the others, `x₃₀ = x₆ + 1 = (17 − √33)/12 + i(5√3 + √11)/12` (the third largest single gain on
`X₂₃`), at distance 1 from `x₂`, `x₆`, `x₁₂` and `x₂₅`. `X₃₁` adds the next best of unit degree 3 on `X₂₃` (the sixth largest single gain), `x₃₁ = x₇ − (1 + i√11)/2 =
(9 − √33)/12 + i(5√3 − 3√11)/12`, at distance 1 from `x₃`, `x₁₇`, `x₂₁` and `x₂₄`.
