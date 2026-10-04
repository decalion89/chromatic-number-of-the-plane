# The two-prime probe: Theorems E and F

Companion to `notes/circular_planes.md` §6 and to Section 10 of `papers/three-colours/`. All programs use exact
rational arithmetic (`fractions.Fraction`, or integers); linear-programming values are certified by exact primal and
dual solutions.

- **Theorem E.** If `χ(F²) ≥ 4` for a number field `F`, then `χ_c(F²) ≥ 7/2` (sharp: `χ_c(ℚ(√11)²) = 7/2`).
- **Theorem F.** If `F²` maps to `K_{p/q}` with `p/q < 4`, then some prime above 2 ramifies in `F(i)`, or some prime
  above 3 has residue degree 1, or some place above 7 has residue degree 1. So `χ_c(F²) ∈ {2, 3, 7/2}` or
  `χ_c(F²) ≥ 4`, decided by these local conditions.

The probe of Theorem D uses only the rotations `iᵃρʲ`, `ρ = (3 + 4i)/5`; below `3/10` it has 5-adic characters. Adding
one rotation with another prime, `σ = (3 + 2i)/(3 − 2i) = (5 + 12i)/13`, removes them. Each theorem has one finite
computer-assisted step, over a window `G(k, m) = {iᵃρʲσˡ : |j| ≤ k, |l| ≤ m}`:

- **Window (1, 1), modulus 65** (Proposition E1 of the note, Proposition 8 of the paper). For `2/7 < θ ≤ 1/3`, every
  point keeping `G(1, 1)` in `[θ, 1 − θ]` has the strip indices of one of the 5 type points (c or q). At `θ = 7/25`
  there are 13 components; the 8 others contain the 7-torsion points `65(a + bi)/7` and have largest least margin
  exactly `2/7`.
- **Window (2, 1), modulus 325** (Proposition F1 of the note, Proposition 9 of the paper). For `1/4 < θ ≤ 1/3`, the
  same with the type points of the families c, q and 7. At `θ = 249/1000` there are 29 components: 5 main, 8 of
  family 7, and 16 with largest least margin exactly `1/4` (2-adic: 4- and 8-torsion points).

Found by a research agent (the files in this folder); refereed by two separate agents, one for each theorem, each with
two exact methods of its own ([`indep_E/`](indep_E/) for Theorem E, [`indep_F/`](indep_F/) for Theorem F). Neither
found an error; their corrections are applied in the note and the paper. The agent's notes
`THEOREM_E_seven_halves.md` and `THEOREM_F_local_global_below_4.md` are kept as written before refereeing; they number
the results differently:

| agent's notes | `notes/circular_planes.md` | `papers/three-colours/` |
|---|---|---|
| Proposition 7, Lemma 12 | Proposition E1, Lemma E2 | Proposition 8, Lemma 12 |
| facts (F1)–(F3) and the digit patterns | (A1)–(A3), §6.4 | Lemma 13 |
| Proposition 7′, Lemma 12′ | Proposition F1, Lemma F2 | Proposition 9, Lemma 14 |
| the claim of §4, Lemmas B2′–B5′ | Lemmas F3–F7 | Lemmas 15–19 |

| file | what it does |
|---|---|
| `gen_tree.py` → `tree_r7_25.txt.gz`; `prop7_certificates.py` → `prop7_certificates.txt` | window (1, 1): a branch tree over all 65² cells whose closed branches carry exact Farkas certificates (4 225 cells, 15 620 branch nodes, 4 364 closed branches, 13 leaves), and the 13 index vectors with their type points or dual certificates |
| `verify_tree.py`, `verify_prop7.py` | checkers for the window (1, 1) that share no code with the generators (about a second together) |
| `certify_window.py` → `tree_K2M1.txt.gz`, `cert_K2M1.txt` | window (2, 1): the tree over all 325² cells (550 474 branch nodes, 199 920 closed branches, 29 leaves) and the type points and dual certificates (`python3 certify_window.py 2 1 249/1000`, about 2 minutes) |
| `verify_window.py` | checker for the window (2, 1) (`python3 verify_window.py 2 1`, 30 to 45 s); rejects corrupted copies |
| `seven_certificates_K2M1.py` → `cert_K2M1_seven.txt`; `verify_seven_K2M1.py` | the 8 index vectors of family 7 in the window (2, 1): dual certificates with value `2/7` (three margins each, among them margins at `ρ^{−2}σˡ`), so that for `θ > 2/7` the window (2, 1) leaves only the types c and q and gives Theorem E by itself; with a checker that shares no code with the generator (it rebuilds `325ρʲσˡ` as products of `2 ± i` and `3 ± 2i`). Added after the referees' reports: the referee of Section 10 noted that the stored certificates gave only the lower bound `2/7` for these vectors |
| `tight_walks.py` → `tight_walks.txt` | consistency checks for the corollary on finite witnesses (paper Corollary 7, note Corollary F12; its proof, through Lemma 22, does not need them): the eight certificates of value `2/7` have three positive multipliers and margins exactly `2/7` at their 7-points, and give positive integer relations among the tight rotations (`63, 52, 25` or `13, 14, 15`: closed walks of length 140 or 42); for type q, three tight rotations of `G_5` with multiplicities `3, 4, 5` in every class and for both signs; the Moser spindle has `χ = 4` and maps to `K_{7/2}` |
| `probe2d.py`, `lpexact.py`, `kappa2d.py`, `kappa2d_seven.py` | the two-prime probe by exact lifting (`5` and `13` steps in a given order), with labels c/q/7 and the certified largest least margin of every other component |
| `indep2d.py`, `indepN.py`, `hgeom_ref.py` → `indep_325_r249.txt` | direct enumeration over all cells (cross-check; the geometry is a referee's, from `../probe/indep/`) |
| `seven_local.py`, `seven_patterns.py` | the facts about `𝔽₄₉` (the set `A′₇` of 8 elements, one `μ₈`-orbit, no coset of a subgroup) and the digit patterns of the 7-adic characters |
| `subwin.py`, `kappa2d_r*.txt`, `seven_r*.txt` | sub-windows and other windows (context) |
| `gdyn.py`, `test_counts.py`, `glue_denominators.py` | the one-prime probe as digit words: at `3/10` its types have unbounded denominators (context) |
| `torsion_exact.py`, `torsion_twoprime.py` | torsion characters: with one prime, orders 41 (`κ = 12/41`) and 76 (`κ = 11/38`), and their multiples, beat `2/7`; with two primes none does (context) |

The stored outputs were produced with these arguments, and all of them reproduce exactly from the programs (the
two trees byte for byte): `torsion_exact.py 400`, `torsion_twoprime.py 300`, `glue_denominators.py`,
`kappa2d.py 2501/10000 5,13` (window (1, 1)), `… 5,13,13` (window (1, 2)), `… 5,13,5` (window (2, 1)),
`kappa2d.py 2500001/10000000 5,13,5`, `kappa2d_seven.py 249/1000 5,13,5`, `kappa2d_seven.py 13/50 5,13,5,5,5`,
`kappa2d_seven.py 1429/5000 5,13,5,5,5,5`, `prop7_certificates.py`, `gen_tree.py` and
`certify_window.py 2 1 249/1000`. The last three lines of `seven_patterns.txt` came from checks that were missing from
the stored copy of `seven_patterns.py`; they were written again when this folder was assembled, and the program now
reproduces the file.

A third referee read Section 10 of the paper with programs of its own (exact clipping of both windows, the lifting of
`S^θ(k, 1)` to the levels `k = 6` and `12`, and checks of the facts at 7); its corrections are applied in the paper and
the note.

Referees' programs:

| folder | what it checks |
|---|---|
| [`indep_E/`](indep_E/) | the window (1, 1) by exact clipping cell by cell (`ref_clip.py`) and by enumerating all vertices of the line arrangement (`ref_vertex.py`); the exact largest least margin of the 13 components (`ref_kappa.py`); the conclusion of Lemma E2 tested exactly up to `k = 7` (`ref_lemma12.py`); the 5-adic points `c_k` under `σ` (`ref_ck.py`); the 7-adic components (`ref_sevenadic.py`); torsion characters by exact order (`ref_torsion.py`) |
| [`indep_S10/`](indep_S10/) | the referee of Section 10 of the paper: its own exact clipping of both windows (`clip_window.py K M r0 cert`: the same 13 and 29 index vectors as the certificates), the exact lifting of `S^θ(k, 1)` to the levels `k = 6` and `12` for Lemma 14 (`lift_lemma14.py`; and at `k = 5`, components of family 7 not of the form of the lemma, so `6 \| k` is needed), the counterexample to the bound `0.23` near `θ = 1/4` (`lift_lemma14_witness.py`), the facts at 7 and the type points modulo 325 (`check_seven.py`, which reads `run/cert_K2M1.txt`, a copy of `../cert_K2M1.txt`), the restriction of the family-7 vectors to `G(1, 1)` (`check_remark1.py`), sub-windows (`subwindows.py`), the points `c_k` (`check_ck.py`), torsion (`check_torsion.py`), and the margin `5/119` of an archimedean point of `ℚ(i)` (`check_qi_char.py`). The lines that print `False` are expected: one digit does not determine the type, an earlier draft's reading of Lemma 13 (corrected), and above `2/7` only the 5 main vectors of the 13 remain |
| [`indep_FW/`](indep_FW/) | the referee of the corollary on finite witnesses (paper Corollary 7, note Corollary F11): the eight certificates of value `2/7` and the 16 of value `1/4` with conventions derived from scratch (`check_seven_certs.py`); type q (`check_type_q.py`); its own exact clipping of `S^θ(2, 1)` (`window21_components.py θ`, then `analyse_components.py` on its output: 13 components at `2/7`, the eight of family 7 single points; 29 at `249/1000`, 13 at `1/4 + 10⁻⁷`, 5 just above `2/7`); Lemmas F9 and F10 on small Cayley graphs (`test_lemmas_AB.py`, `test_lemmaB_72*.py`); the tight relations of the 7-adic colourings (`seven_adic_tight.py`); the Moser spindle (`moser.py`); and the 76-vertex graph over `ℚ(√11)` (`q11_circular.py`, `q11_partial.py`: `16/5` SAT, `19/6` and `35/11` UNSAT, not certified; the runs at `51/16` and `67/21` did not finish) |
| [`indep_FW2/`](indep_FW2/) | the referee of the general lemma behind the corollary (paper Lemma 22, note Lemma F11, winding paper Section 7): exact `κ(S)` for `ℤ^r × ℤ/m` by vertex enumeration (`kappa_exact.py`, checked by `sanity_kappa.py`); 115 random finite sets `S` with `2 < χ_c < 4`, whose 409 optimal characters all have a nonnegative relation among their tight elements (`test_kappa_relations.py` → `kappa_examples.json`); SAT on finite pieces of 13 Cayley graphs (`sat_tight.py`, `test_finite_pieces.py` → `finite_pieces_results.json`; for example `{0, …, 18}` for the distances `3, 4, 9, 12`, with `χ_c = 7/2`); random periodic homomorphisms (`test_lemma21_periodic.py`); the radii of Lemma 9 (`check_radii.py`: the supremum of `r_θ` on `(2/7, 1/3]` is `√10/14 < 0.23`); and the 140 vectors of Theorem C with the tight relations of their 7-adic characters (`thmC_vectors.py`, `thmC_tight7.py`) |
| [`indep_F/`](indep_F/) | the window (2, 1) by vertex enumeration (`vertex_enum.c`, 128-bit integers, with `analyze_components.py`, exact) and by lifting rows (`lift_rows.py`), compared by `compare_methods.py`; the facts at 7 and the digit patterns (`f49_and_digits.py`); Lemma F6 on finite models and Lemma F7 for `f = 3` (`b4_model_check.py`, `b5_f3_check.py`); known answers of Theorem D's probe (`an_kat_*.txt`) |

To recheck the computer-assisted steps (also `tests/test_two_primes.py`):

```sh
gunzip -k tree_r7_25.txt.gz tree_K2M1.txt.gz
python3 verify_tree.py && python3 verify_prop7.py      # window (1, 1)
python3 verify_window.py 2 1                           # window (2, 1)
python3 verify_seven_K2M1.py                           # the vectors of family 7 there, at 2/7
cd indep_E && python3 ref_vertex.py 7 25 && python3 ref_clip.py 7 25 clip.txt && python3 ref_kappa.py clip.txt
cd ../indep_F && cc -O2 -o vertex_enum vertex_enum.c   # then analyze_components.py on its output
cd ../indep_S10 && mkdir -p run && cp ../cert_K2M1.txt ../prop7_certificates.txt run/ && python3 check_seven.py
```
