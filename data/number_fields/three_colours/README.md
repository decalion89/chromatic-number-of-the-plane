# Three colours for the plane over a number field: computations

Companion to `notes/three_colours_number_fields.md` (Theorem B, Proposition B9). Everything is exact (rational
arithmetic in the field, lattices with PARI/GP); floating-point solvers only guide searches.

| file | what it does |
|---|---|
| `classify.gp`, `classify_examples.gp` | Theorem B's criterion with PARI: (a) a prime above 2 ramified in `F(i)` (relative discriminant), (b) a prime above 3 of residue degree 1; `gp -q classify_examples.gp` prints `[a, b]` for the fields of the note |
| `nf.py` | exact arithmetic in `F = ℚ[t₁, …, t_r]/(f₁(t₁), …, f_r(t_r))` and unit vectors `(p + i)²/(p² + 1)` |
| `limit_test.py` | the `N → ∞` test of Lemma B1 (one-sided): is there `β ∈ F(i)` with `Tr(βv) ∈ E` (types c and q) for every `v` of a finite set `V`? The values `6φ(v)` form the lattice `(column space) ∩ ℤ^{2m}` (PARI `matrixqz(A, −2)`); the types are conditions modulo 2 and 3, so `L/2L × L/3L` is enumerated. "INFEASIBLE" means that `G_N V` is not 3-colourable for large `N = 5^k`, so `χ(F²) ≥ 4`; "feasible" proves nothing |
| `check_saved.py V.json` | rebuilds a saved set, checks the unit vectors, runs the test |
| `check.gp`, `to_gp.py` | the same test in PARI/GP alone (own unit-vector check, lattice by `matkerint`): `python3 to_gp.py V.json v.gp && cat v.gp check.gp > run.gp && gp -q run.gp` |
| `known_cases.py` | Theorem 1's vectors: infeasible for `d = 11, 23, 35, 47, 59`; `ℚ(√3)`, `ℚ(√7)` feasible |
| `search.py` | random search for an infeasible set, with vectors that are large at the split places above 2 and 3 |
| `minimize_v.py` | greedy deletion inside an infeasible set |
| `residue_lemma.py` | Lemma B5 by exhaustive search for `f = 1, 3, 5` |
| `V_sqrt2_sqrt7.json`, `V_sqrt2_sqrt7_min5.json` | 27 and 5 unit vectors of `ℚ(√2, √7)²`: infeasible |
| `V_c7_sqrt7.json` | 44 unit vectors of `ℚ(2cos(2π/7), √7)²`: infeasible |
| `V_two_roots_2_7.json` | `{1, u_2, ū_2, u_7, ū_7}` of Proposition B9 for `ℚ(√2, √7)` |
| `gen_kappa.py V.json N [A]` | builds `G_N V` (one vector per pair `±u`) with its integer relations (PARI `matkerint`) and asks a MIP solver (SCIP) for a character into `[A, 1 − A]` (a guide only) |
| `cert_general.py V.json N OUT` | exact certificate (branch and bound over the integer values of the relations, exact Farkas vectors at the leaves) that no character of `ℤ(G_N V)` maps into `[1/3, 2/3]` (multiquadratic fields) |
| `check_mq.py CERT` | independent exact checker for those certificates (unit lengths in the multiquadratic field, relations, tree) |
| `cert_two_roots_2_7_N25.json` | `G_25{1, u_2, ū_2, u_7, ū_7}`, 50 vectors: no character into `[1/3, 2/3]` (1 579 nodes) |
| `cert_sqrt2_sqrt7_N125.json.gz` | `G_125 V` for the five vectors of `V_sqrt2_sqrt7_min5.json`, 70 vectors (11 357 nodes) |

With Theorem W (proved in Lean for every abelian group and finite `S`), each certificate proves `χ(ℚ(√2, √7)²) ≥ 4`.
