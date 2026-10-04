# Circular colourings of field planes: computations

Companion to `notes/circular_planes.md` (Proposition C1, Theorem C). Uses `../three_colours/nf.py` and
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
| `kappa1.py p,f …` | `κ₁(p, f)`: the best residue-field circular colouring at a place with residue field `𝔽_{p^f}` (Proposition C2) |
| `level2.py` | the level lemma of Proposition C2 on finite models: `O/49` for `ℚ₇`, `O/9` for `ℚ₂₇` |
| `chartypes.py k h thr`, `lift.py THR KMAX [types file]` | maps of the characters of `(1/N)ℤ[i]` that keep `G_N` in `[r, 1 − r]`, and their coherent lifts from `N` to `5N` (§5 of the note) |
