# Density certificates for planar sets avoiding two distances

Each file in `certificates/` is a dual solution of the linear program of Oliveira Filho and Vallentin
(J. Eur. Math. Soc. 2010, Theorem 1.1), with clique inequalities added for `{1, 2}` and `{1, √3}`. It proves an
upper bound `B` on `m_D(ℝ²)`, the largest density of a measurable planar set with no two points at a distance in
`D`. By the transfer theorem of [`papers/forbidden-distances/`](../../papers/forbidden-distances/README.md),
`χ(ℝ², D) ≥ 1/m_D(ℝ²)` for every set `D` of algebraic distances, so `B < 1/k` gives `χ(ℝ², D) ≥ k + 1`.

| file | `D` | clique sets | bound on `m_D(ℝ²)` | `T` | proves |
|---|---|---|---|---|---|
| `cl_2.cert.json` | `{1, 2}` | 6 unit triangles, 2 triangles of side 2 | 0.159222 < 1/6 | 1297 | `χ ≥ 7` |
| `cl_s3.cert.json` | `{1, √3}` | a unit triangle, a rhombus | 0.162866 < 1/6 | 757 | `χ ≥ 7` |
| `s_vs.cert.json` | `{1, (√6+√2)/2}` | none | 0.166018 < 1/6 | 107 | `χ ≥ 7` |
| `s_10.cert.json` | `{1, 10}` | none | 0.149599 < 1/6 | 36 | `χ ≥ 7` |
| `s_116.cert.json` | `{1, 11/6}` | none | 0.140995 < 1/7 | 171 | `χ ≥ 8` |
| `s_4625.cert.json` | `{1, 46/25}` | none | 0.139185 < 1/7 | 175 | `χ ≥ 8` |
| `s_40.cert.json` | `{1, 40}` | none | 0.115511 < 1/8 | 36 | `χ ≥ 9` |
| `s_100.cert.json` | `{1, 100}` | none | 0.103302 < 1/9 | 37 | `χ ≥ 10` |

The bounds are rounded up; the exact values are fractions, printed by the programs.

**Format.** A certificate is a JSON object with
- `D`: the distances, as expressions separated by `;`;
- `z0` and `z`: the coefficients, as exact fractions in strings;
- `cliques`: for each clique inequality, a shape of [`shapes.py`](shapes.py) (`W_c`), the offset `(ox, oy)` of the
  point `p_c` from the centroid of the shape, and the weight `y ≥ 0`;
- `T`: the end of the sweep.

The claim is that `g(t) = z0 + Σ_k z_k J₀(d_k t) + Σ_c y_c Σ_{w ∈ W_c} J₀(|w − p_c| t)` is nonnegative for all
`t ≥ 0`. Then `m_D(ℝ²) ≤ B = (z0 + Σ_c y_c) / (z0 + Σ_k z_k + Σ_c y_c |W_c|)` (Proposition 10 of the paper).

**Programs** (run them in this folder):
- `lp_simple_json.py "1;40" lp/s_40.json` and `lp4.py "1;2" tri111,tri222 lp/cl_2.json` solve the linear programs
  in floating point. They only choose a certificate; nothing they compute is trusted.
- `certify_lp2.py lp/NAME.json k [certificates/NAME.cert.json]` rounds the coefficients to fractions, checks the
  cliques exactly, proves `g ≥ 0` on `[0, T]` with Arb balls (point values and the Lipschitz bound
  `|J₁(x)| ≤ 1/√2`, DLMF 10.14.1) and on `[T, ∞)` with `|J₀(x)| ≤ √(2/(πx))`, and writes the certificate.
- `check2.py certificates/NAME.cert.json k` is separate code that repeats the check in 60-digit floating point.
- `check3.py certificates/NAME.cert.json k` is a rigorous check independent of both: it builds the cliques and
  their radii in exact arithmetic of `ℚ(√3)`, encloses `g` with Arb over whole intervals of `t`, bisecting until
  each enclosure is nonnegative, and checks the tail. It is adapted from the checker of the second referee.

`python3 check3.py certificates/cl_2.cert.json 6` takes about ten seconds; the other certificates take one second
or less. `certify_lp2.py` reproduces every stored certificate byte for byte.

**Tests.** From the root of the repository, `python3 -m pytest -q -m "not slow" tests/test_forbidden_distances.py`
runs `check3.py` on all eight certificates, `check2.py` on six of them and `certify_lp2.py` on one (under a minute;
CI runs these). Without `-m` the file also runs `check2.py` on the other two and rebuilds every certificate with
`certify_lp2.py` (about five minutes).

**Requirements.** `python-flint` (Arb), `sympy` and `mpmath`; the two linear-programming scripts also need `numpy`
and `scipy`.
