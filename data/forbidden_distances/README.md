# Density certificates for planar sets avoiding given distances

The folder [`uniform/`](uniform/README.md) holds the certificates of Theorem 5 of the paper, each valid for a whole range
of `d`: `χ(ℝ², {1, d}) ≥ 7` for every algebraic `d ≥ 5.82` and for `1.72 ≤ d ≤ 2.04`, up to `≥ 14` for every
algebraic `d ≥ 25 700`; its subfolder [`uniform/ie/`](uniform/ie/README.md) has the inclusion–exclusion certificates
for `≥ 15`, `≥ 16` and `≥ 17` from `d = 5·10⁵`, `7·10⁷` and `8·10⁸`.

The folder [`m1/`](m1/README.md) holds the certificates for one distance: `m₁(ℝ²) ≤ 0.24353`, from the
inclusion–exclusion program of Ambrus, Csiszárik, Matolcsi, Varga and Zsámboki for their 23 points and eight more
points of `ℚ(√3, √11)` ([`papers/unit-distance-density/`](../../papers/unit-distance-density/README.md)).

The folder [`m1r3/`](m1r3/README.md) does the same in three dimensions: `m₁(ℝ³) ≤ 0.15158538`, from the same program
carried to `ℝ³`, with the atom `δ²` at the origin, on 27 points `q/√2`, `q ∈ ℚ³`
([`papers/unit-distance-density-3d/`](../../papers/unit-distance-density-3d/README.md)).

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
| `s_1000.cert.json` | `{1, 1000}` | none | 0.089085 < 1/11 | 37 | `χ ≥ 12` |
| `m3_40.cert.json` | `{1, 40, 1600}` | none | 0.046389 < 1/21 | 36 | `χ ≥ 22` |
| `m4_10.cert.json` | `{1, 10, 100, 1000}` | none | 0.039641 < 1/25 | 37 | `χ ≥ 26` |
| `m3_100.cert.json` | `{1, 100, 10⁴}` | none | 0.037168 < 1/26 | 37 | `χ ≥ 27` |
| `m4_30.cert.json` | `{1, 30, 900, 27000}` | none | 0.021264 < 1/47 | 36 | `χ ≥ 48` |

The bounds are rounded up; the exact values are fractions, printed by the programs.

**Format.** A certificate is a JSON object with
- `D`: the distances, as expressions separated by `;`;
- `z0` and `z`: the coefficients, as exact fractions in strings;
- `cliques`: for each clique inequality, a shape of [`shapes.py`](shapes.py) (`W_c`), the offset `(ox, oy)` of the
  point `p_c` from the centroid of the shape, and the weight `y ≥ 0`;
- `T`: the end of the sweep.

The claim is that `g(t) = z0 + Σ_k z_k J₀(d_k t) + Σ_c y_c Σ_{w ∈ W_c} J₀(|w − p_c| t)` is nonnegative for all
`t ≥ 0`. Then `m_D(ℝ²) ≤ B = (z0 + Σ_c y_c) / (z0 + Σ_k z_k + Σ_c y_c |W_c|)` (Proposition 17 of the paper).

**Programs** (run them in this folder):
- `lp_simple_json.py "1;40" lp/s_40.json`, `lp_multi.py "1;100;10000" lp/m3_100.json` (a cutting-plane version
  for distances far apart) and `lp4.py "1;2" tri111,tri222 lp/cl_2.json` solve the linear programs in floating
  point. They only choose a certificate; nothing they compute is trusted.
- `certify_lp2.py lp/NAME.json k [certificates/NAME.cert.json]` rounds the coefficients to fractions, checks the
  cliques exactly, proves `g ≥ 0` on `[0, T]` with Arb balls (point values and the Lipschitz bound
  `|J₁(x)| ≤ 1/√2`, DLMF 10.14.1) and on `[T, ∞)` with `|J₀(x)| ≤ √(2/(πx))`, and writes the certificate.
- `check2.py certificates/NAME.cert.json k` is separate code that repeats the check in 60-digit floating point.
- `check3.py certificates/NAME.cert.json k` is a rigorous check independent of both: it builds the cliques and
  their radii in exact arithmetic of `ℚ(√3)`, encloses `g` with Arb over whole intervals of `t`, bisecting until
  each enclosure is nonnegative, and checks the tail. It is adapted from the checker of the second referee.

`python3 check3.py certificates/cl_2.cert.json 6` takes about ten seconds; the other certificates take one second
or less. `certify_lp2.py` reproduces every stored certificate byte for byte; for `s_1000` and the three- and
four-distance ones its sweep takes from twenty seconds to forty minutes (up to 26 million steps for `m4_30`), and
`check2.py`, whose Lipschitz bound is global, takes minutes for `s_1000` and hours for the others.

Every program exits with a non-zero code when a check fails, including when the bound is not below `1/k`; the
checks do not use `assert`, so they also run under `python3 -O`.

**Tests.** From the root of the repository, `python3 -m pytest -q -m "not slow" tests/test_forbidden_distances.py`
runs `check3.py` on all the certificates, `check2.py` on six of them and `certify_lp2.py` on one, and checks that
`check3.py` rejects three broken certificates and both checkers reject a claim that is too strong (about a minute;
CI runs these). Without `-m` the file also runs `check2.py` on the others and rebuilds every certificate with
`certify_lp2.py` (hours, because of `check2.py` on the three- and four-distance certificates).

**Requirements.** `python-flint` (Arb), `sympy` and `mpmath`; the two linear-programming scripts also need `numpy`
and `scipy`.
