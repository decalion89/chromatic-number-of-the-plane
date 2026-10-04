# Circular colourings of field planes: computations

Companion to `notes/circular_planes.md` (Proposition C1, Theorems C, D, E and F). Uses `../three_colours/nf.py` and
`../three_colours/gen_kappa.py`.

| file | what it does |
|---|---|
| `residue_circular.py` | for `p ≡ 3 (mod 4)`, the best `𝔽_p`-linear `λ` on `𝔽_{p²}` keeping `μ_{p+1}` away from 0: `κ₁(p)` and the bound `χ_c ≤ 1/κ₁(p)` (Proposition C1) |
| `kappa_max.py d N n1,n2,…` | `max_ξ min_u ‖ξ(u)‖` over the characters of `ℤU`, `U = G_N{1, u_n, ū_n}` in `ℚ(√d)²` (MIP, SCIP; a guide) |
| `cert_open.py OUT d N n1,… P Q` | exact certificate that no character maps `U` into the open interval `(Q/P, 1 − Q/P)` |
| `check_open.py CERT` | independent exact checker (unit lengths in `ℚ(√d)`, relations, branches over the integers strictly inside each range, leaves with Farkas vectors valid for the open box) |
| `check_open_indep.py CERT [d,N,P,Q,n1,…]` | a second exact checker, written by the referee with no code shared with `check_open.py`; with the optional list it also rebuilds `U = G_N{1, u_n, ū_n}` from `d` and compares |
| `cert_sqrt11_open_7_2_N25.json.gz` | `ℚ(√11)`, `N = 25`, `V = {1, u_1, ū_1, u_7, ū_7, u_19, ū_19}`: 70 vectors, no character into `(2/7, 5/7)`; with Theorem W⁺, `χ_c(ℚ(√11)²) ≥ 7/2` |
| `cert_sqrt35_open_7_2_N25.json.gz` | the same for `ℚ(√35)` (14 199 nodes, 10 241 leaves): `χ_c(ℚ(√35)²) ≥ 7/2` |
| `cert_sqrt59_open_53_15_N25.json.gz` | the same 70 vectors for `ℚ(√59)`, interval `(15/53, 38/53)` (127 181 nodes, 105 761 leaves; `cert_open.py ... 59 25 1,7,19 53 15`, 73 min): `χ_c(ℚ(√59)²) ≥ 53/15 > 7/2` although `χ = 4` (Proposition C6). Check with `check_open.py` (9 min) and `check_open_indep.py CERT 59,25,53,15,1,7,19` (23 s) |
| `kappa1.py p,f …` | `κ₁(p, f)`: the best residue-field circular colouring at a place with residue field `𝔽_{p^f}` (Proposition C2) |
| `level2.py` | the level lemma of Proposition C2 on finite models: `O/49` for `ℚ₇`, `O/9` for `ℚ₂₇` |
| `kappa1_f1.py P` | `κ₁(p, 1)` for every prime `p ≡ 3 (mod 4)` below `P`, fast (Proposition C3 needs `P ≥ 1001`; we ran `P = 3000`) |
| `torus_sums.py` | the sums over the norm-one torus against `2√q` (Weil, through Kloosterman sums) |
| `chartypes.py k h thr`, `lift.py THR KMAX [types file]` | maps of the characters of `(1/N)ℤ[i]` that keep `G_N` in `[r, 1 − r]`, and their coherent lifts from `N` to `5N` (§5 of the note) |
| [`probe/`](probe/) | Theorem D (§5 of the note, Section 9 of `papers/three-colours/`): exact computation of the probe sets `S_N^(r)` for `r < 1/3`, the lemmas of Proposition D1, the computer-assisted bound `3.3315`; [`probe/indep/`](probe/indep/) is a referee's implementation that shares no code with it. See [`probe/README.md`](probe/README.md) |
| [`twoprime/`](twoprime/) | Theorems E and F (§6 of the note, Section 10 of the paper): the probe with the second rotation `σ = (5 + 12i)/13`, the certificates of the windows modulo 65 and 325 with their checkers, the facts at 7, and the two referees' programs. See [`twoprime/README.md`](twoprime/README.md) |
