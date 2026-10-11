# The density of sets in three-space avoiding the distance 1

`m₁(ℝ³)` is the largest upper density of a measurable set in `ℝ³` with no two points at distance 1. The best published
bound was `m₁(ℝ³) ≤ 0.1532996` (DeCorte, Oliveira Filho and Vallentin, *Complete positivity and distance-avoiding
sets*, Math. Program. 191 (2022) 487–558, arXiv:1804.09099, Table 2). This folder proves

> **`m₁(ℝ³) ≤ 0.15158538`**, so the measurable fractional chromatic number of `ℝ³` is at least `6.5969`.

The method is the inclusion–exclusion linear program of Ambrus, Csiszárik, Matolcsi, Varga and Zsámboki (Math. Program.
207 (2024) 303–327, arXiv:2207.14179), which they used only in the plane, carried to `ℝ³` (the Bessel function `J₀`
becomes `Ω(x) = sin x / x`) and combined with the atom `κ({0}) = δ²` that the spectral measure of a set of density `δ`
has at the origin (as in Keleti, Matolcsi, Oliveira Filho and Ruzsa, Prop. 3.5). The configuration `X₂₇` has 27 points
`q/√2` with `q ∈ ℚ³`, grown from a cuboctahedron and its centre, so all squared distances are rational and everything
is exact.

The measurable chromatic number bound `χ_m(ℝ³) ≥ 7` is not changed (that would need `m₁(ℝ³) < 1/7`).

The paper is [`papers/unit-distance-density-3d/`](../../../papers/unit-distance-density-3d/README.md).

| points | independent sets | value (floating point) | certified bound |
|---|---|---|---|
| 1 (any single point) | 2 | 0.178467 | (Oliveira Filho–Vallentin's bound (8) for `n = 3`: 0.1785098…) |
| 13 (cuboctahedron and centre) | 109 | 0.178467 | |
| 19 | 1 459 | 0.156470 | |
| 23 | 8 430 | 0.152137 | 0.15218069 |
| 26 | 26 117 | 0.151549 | 0.15160217 |
| **27** (`X₂₇`) | **44 221** | **0.151456** | **0.15158538** |

The value is that of the program on the grid `t = 0, 0.1, …, 100` and `100.25, …, 400` used by `grow3.py`; the
certificates come from a finer program (step `0.02` up to `t = 60`), whose value for `X₂₇` is `0.15157938`.

## Files

- `configurations/x27.json`: the points. Each point is `[x, y, z]` (rationals as strings) for the point
  `(x, y, z)/√2` (`scale: "1/sqrt2"`). The first 13 are `q/√2` for the `q ∈ ℤ³` with `q₁ + q₂ + q₃` even and
  `|q|² ≤ 2`, in lexicographic order; each later point is at distance 1 from at least three earlier ones.
- `certificates/m1r3_x27.json`: the certificate (format below), and `m1r3_x27.json.lp.npz`, the solution of the first
  program, from which `m1cert3.py --reuse` rebuilds the certificate byte for byte.
- `check_m1r3.py CERT.json`: the rigorous check (exact rationals, Arb balls for `Ω`); self-contained (numpy and
  python-flint), no `assert`, run as `python3 -I check_m1r3.py ...`.
- `check_m1r3_indep.py CERT.json`: a second check, written separately without reading the first: its own enumeration
  and counts of the independent sets, exact evaluation of `Q`, and a third-order Taylor sweep with an explicit bound on
  `Φ'''` (Arb only for enclosures of `sin`, `cos` and `√`).
- `m1cert3.py CONFIG.json OUT.json`: builds a certificate (floating point; nothing it computes is trusted).
- `ielp3.py`: the program in floating point (HiGHS), with the tangent cuts `κ(0) ≥ 2aδ − a²`.
- `conf3.py`: exact configurations, the unit-distance graph, the independent sets and the congruence relations.
- `grow3.py SEED.npy OUT_PREFIX NMAX [MAXEVAL] [MINDEG] [--rational]`: the greedy growth (floating point).
- `referee_check.py CERT.json`: the internal referee's own program (brute force over all `2²⁷` subsets, exact `min Q`,
  the congruences, `Φ − m` on a dense grid: a floating-point check only).
- `logs/`: the output of both checkers and of the referee's program on the certificate, and the three growth runs
  (13 → 19, 19 → 24, 24 → 27 points).

## The certificate

A certificate is a JSON object with `side: "M1R3"`, `scale: "1/sqrt2"`, the `points`, rationals `psi0`, `psiT`, `m`,
`lam` (one value per point), `u` (`[i, j, value]`, 1-based, only for pairs not at distance 1) and `wc` (`[I, J, value]`:
two lists of point indices, congruent in the given order). With `r_ij = |x_i − x_j|` and `Ω(x) = sin x / x` it claims

- `Q(S) = psiT + Σ_{i∈S} lam_i + Σ_{i<j∈S} u_ij + Σ_wc value·([I ⊆ S] − [J ⊆ S]) ≥ 0` for every independent set `S` of
  the unit-distance graph, the empty set included;
- `Φ(t) = psi0·Ω(t) − Σ lam_i − Σ u_ij Ω(r_ij t) ≥ m` for every `t > 0`; and
- `c = Φ(0) − m > 0`.

Then `m₁(ℝ³)` is at most the positive root of `c x² + m x − psiT` (Proposition 2 of the paper). For a periodic set `A`
of density `δ` avoiding 1, choose a uniform translation and a Haar-random `g ∈ O(3)` (reflections included: 442 of the
1 078 relations of the certificate relate non-coplanar sets by an orientation-reversing isometry) and let `S` be the
set of points of `X` that land in `A`. Then `𝔼 Q(S) ≥ 0` gives `∫ Φ dκ ≤ psiT`, where `κ` is the radial spectral measure of
`A` (total mass `δ`, atom `δ²` at 0, `∫ Ω dκ = 0`), and `∫ Φ dκ ≥ m(δ − δ²) + Φ(0)δ² = mδ + cδ²`.

```
python3 -I check_m1r3.py certificates/m1r3_x27.json        # CERTIFIED: ... at most 0.15158538 = 1/6.59694 (13 s)
python3 -I check_m1r3_indep.py certificates/m1r3_x27.json  # ACCEPTED: upper density <= B <= 0.15158538 (1 s)
```

## Reproducing

```
python3 -c "import numpy as np, itertools; np.save('fcc13.npy', np.array([p for p in itertools.product(range(-3, 4),
  repeat=3) if sum(p) % 2 == 0 and sum(c * c for c in p) <= 2], float) / np.sqrt(2))"
python3 grow3.py fcc13.npy g 19 40                    # points 14-19 (as in logs/growth_1_seed13_to19.log)
python3 grow3.py g_19.npy r 24 40 3 --rational        # points 20-24
NPROC=2 python3 grow3.py r_24.npy r 27 48 3 --rational  # points 25-27
python3 conf3.py r_27.npy x27.json
python3 m1cert3.py x27.json m1r3_x27.json             # about 3 minutes; --reuse with the stored .lp.npz is exact
```

Floating-point ties may resolve differently on another machine; the configuration in `configurations/x27.json` is the
one that is certified.
