# Referee report: the finite witnesses for `χ_c = 7/2` (`witness_q11`, `witness_q191`, `witness_q455`)

Written by a referee working separately, with programs of its own (this folder), on 4 October 2026. **Verdict:** all
claims about the three witnesses hold: each is an induced unit-distance graph in `ℚ(√d)²` with `χ_c(H) = 7/2`, and it
is vertex-critical. No error affecting the results was found. The findings (§8) concern wording, the description of
the search, and the robustness of the repository's checkers; all of them have since been applied (see `../README.md`
and `CHANGELOG.md`). The χ_c = 3 material added during the review was not part of it.

## 0. Scope and environment

- **Object:** `data/number_fields/circular/finite_witness/` at commit 5cc8ae84 (`README.md`, `witness_q11.*`,
  `witness_q191.*`, `witness_q455.*`, `*_critical.json.gz`), and in `papers/three-colours/three-colours.tex`
  Section 10 ("Finite witnesses"), Lemma 20, Corollary 7, Theorem F and the sentences of the introduction and of the
  Questions that cite the witnesses.
- **Concurrent edits:** commits c4811adc and 176d5594 landed during the review (a χ_c = 3 witness over `ℚ(√7)`;
  `grow.py`, `minimise.py`, `critical.py` generalised to any `p/q`). The data files of the three 7/2 witnesses, the
  three checkers and the Section 10 paragraph on them were unchanged.
- **Tools:** Python 3.11.15, python-sat 1.9.dev15, kissat 4.0.4, drat-trim, cake_lpr. Times are wall clock on a
  shared 4-core machine.

## 1. Points, unit edges, inducedness — confirmed

`referee_check.py` (written from scratch): `d` is a squarefree integer greater than 1, so distinct integer 4-tuples
are distinct points; points are integer 4-tuples and distinct; the edge list has no loops, duplicates or
out-of-range indices; all pairs `i < j` are compared, and with `(A, B, C, E) = p_j − p_i` the pair is at distance 1
if and only if `A² + dB² + C² + dE² = D²` and `AB + CE = 0` (integer arithmetic). `unit_vectors.py` enumerates every
unit vector with denominator `D` by brute force.

| witness | d, D | points | pairs | unit pairs = listed edges | listed but not unit / unit but not listed | unit vectors with denominator D / used as edges | degrees | connected |
|---|---|---|---|---|---|---|---|---|
| q11 | 11, 30 | 170 | 14 365 | 468 = 468 | 0 / 0 | 108 / 56 | 2–20 | yes |
| q191 | 191, 240 | 293 | 42 778 | 803 = 803 | 0 / 0 | 108 / 52 | 2–34 | yes |
| q455 | 455, 780 | 175 | 15 225 | 434 = 434 | 0 / 0 | 68 / 34 | 3–16 | yes |

Negative controls (`mutation_tests.py`): a changed coordinate, a removed edge, an added non-unit pair, a duplicated
point and a non-squarefree `d` are all rejected (`results/mutation_tests.log`).

## 2. The `(7, 2)`-colourings — confirmed

Colours are integers in 0..6, and `c(y) − c(x) mod 7 ∈ {2, 3, 4, 5}` on every edge of the three graphs (0
violations), so `χ_c(H) ≤ 7/2`. A changed colour is rejected.

## 3. Every `(7, 2)`-colouring has a tight cycle — confirmed

### 3.1 Listed cycles

Each listed cycle has at least 3 distinct vertices, consecutive vertices (cyclically) are adjacent, and every length
is a multiple of 7 (necessary for a tight cycle: `2m ≡ 0 mod 7`). No duplicates, even up to rotation.

| witness | cycles | lengths (count) | arcs | tight under the stored colouring | fixed_vertex (degree; max degree; first vertex of max degree?) |
|---|---|---|---|---|---|
| q11 | 879 | 7: 330, 14: 429, 21: 110, 28: 10 | 902 | 1 | 34 (20; 20; yes) |
| q191 | 489 | 7: 107, 14: 275, 21: 81, 28: 24, 35: 2 | 1335 | 1 | 56 (34; 34; yes) |
| q455 | 91 | 7: 3, 14: 68, 28: 20 | 616 | 1 | 34 (16; 16; yes) |

### 3.2 Fixing one colour is harmless

If `c` is a colouring with a non-tight arc on every listed cycle, `c' = c − c(v₀) mod 7` has the same colour
differences, hence is a `(7, 2)`-colouring with the same tight arcs, and `c'(v₀) = 0`. (Reflections `c ↦ −c` would
reverse the tight arcs; they are rightly not used.) Empirically, the q455 formula without the unit clause is refuted
too: kissat UNSAT in 33.4 s, drat-trim VERIFIED in 26.0 s, 169 337 of 376 427 lemmas in the core (21 996 of 43 107
with the unit clause).

### 3.3 The referee's encodings

- **E1** (`referee_check.py`): `x(v,k) = k·n + v + 1` ("v has colour k"), and `N(a)` for each arc `a = (u, w)` of a
  listed cycle with `N(a) → ¬(x(u,k) ∧ x(w,k+2))` for every `k`; at-least-one and at-most-one colour per vertex;
  `¬x(u,k) ∨ ¬x(w,k+δ)` for every edge, `k` and `δ ∈ {0, 1, 6}`; one clause `∨N(a)` per listed cycle; the unit clause
  `x(v₀, 0)`. (Written before reading the repository; its layout turned out close to `verify_independent.py`, hence
  E2.)
- **E2** (`referee_encode2.py`): two-sided tight indicators `(¬x(a,k) ∨ ¬x(b,k+2) ∨ T)`, `(¬T ∨ ¬x(a,k) ∨ x(b,k+2))`,
  cycle clauses `∨¬T`; then all variables randomly renamed and negated, and clauses and literals shuffled (seed
  20261004).
- **Rebuild from scratch** (`referee_cegar.py`): ignoring the stored cycles, CaDiCaL (pysat) is asked for a
  `(7, 2)`-colouring with vertex 0 coloured 0 and a non-tight arc on every cycle found so far; an acyclic tight
  digraph would be a counterexample (none occurred); otherwise, for each nontrivial strong component, shortest
  directed cycles through up to three of its vertices are added. The final formula has the E1 layout with the
  program's own cycles.
- **Semantic validation** (`validate_e1.py`): every clause of the E1 and rebuild formulas is satisfied by the intended
  assignment of any `(7, 2)`-colouring with an acyclic tight digraph; each cycle clause is a simple directed cycle of
  `H`. Result: 0 unexplained clauses in all six files.

### 3.4 Proof runs (`run_proof.sh`: `kissat --seed=7`, `drat-trim -L`, `cake_lpr`)

| formula | variables / clauses | sha256 prefix | kissat | drat-trim | lemmas in core | cake_lpr |
|---|---|---|---|---|---|---|
| q11 E1 | 2092 / 20762 | 1d64fad4 | UNSAT 27.4 s | VERIFIED 31.2 s | 126 807 / 234 983 | VERIFIED UNSAT 11.6 s |
| q191 E1 | 3386 / 33144 | 189de8ff | UNSAT 21.6 s | VERIFIED 31.2 s | 100 085 / 186 743 | VERIFIED UNSAT 9.0 s |
| q455 E1 | 1841 / 17368 | 510d06d6 | UNSAT 3.0 s | VERIFIED 3.7 s | 21 996 / 43 107 | VERIFIED UNSAT 1.1 s |
| q11 E2 | 2092 / 27076 | 6c4864f5 | UNSAT 27.2 s | VERIFIED 30.0 s | 120 718 / 228 006 | VERIFIED UNSAT 11.2 s |
| q191 E2 | 3386 / 42489 | 9c7e0fd2 | UNSAT 26.4 s | VERIFIED 38.1 s | 100 252 / 196 549 | VERIFIED UNSAT 8.3 s |
| q455 E2 | 1841 / 21680 | b6ebe3e4 | UNSAT 2.9 s | VERIFIED 6.0 s | 22 828 / 39 929 | VERIFIED UNSAT 1.4 s |
| q11 rebuild (279 rounds, 659 cycles, lengths 7–35) | 2100 / 20598 | 33551a4e | UNSAT 50.1 s | VERIFIED 52.3 s | 226 294 / 440 644 | VERIFIED UNSAT 28.3 s |
| q191 rebuild (257 rounds, 658 cycles, lengths 7–35) | 3495 / 34076 | 6b5afb34 | UNSAT 43.2 s | VERIFIED 51.2 s | 187 983 / 376 319 | VERIFIED UNSAT 25.1 s |
| q455 rebuild (64 rounds, 165 cycles, lengths 7/14/28) | 1892 / 17799 | fd5e6952 | UNSAT 13.4 s | VERIFIED 11.9 s | 78 842 / 158 182 | VERIFIED UNSAT 5.8 s |
| q455 E1 without unit clause | 1841 / 17367 | 97848c62 | UNSAT 33.4 s | VERIFIED 26.0 s | 169 337 / 376 427 | not run |

Full hashes: `results/*.ref.log`, `results/E2_gen.log`, `results/cegar_*.log`.

### 3.5 The stored certificates

`validate_stored.py` reads every clause of each stored formula under `x(v,k) = 7v + k + 1` and `t = [arc tight]`
(the arc of `t` read off its clauses). Only these clause types occur: at-least-one colour, edge conflict,
"tight ⇒ t", "some `t` of a listed cycle is false", and one `x(v₀, 0)`; each cycle clause is exactly the arc set of a
listed cycle, and every listed cycle is covered. So unsatisfiability of the stored formula implies the claim.
`run_stored.sh` then runs drat-trim (with LRAT output) and cake_lpr on the stored files.

| witness | stored CNF sha256 matches verification.txt | stored DRAT bytes / sha256 matches | clause types (ALO, EDGE, TDEF, CYC, UNIT; other) | drat-trim | cake_lpr |
|---|---|---|---|---|---|
| q11 | f36a7e6a…, yes | 28 824 442 / 6258c240…, yes | 170, 9828, 6314, 879, 1; 0 | VERIFIED 23.2 s | VERIFIED UNSAT 12.9 s |
| q191 | e6223861…, yes | 22 573 623 / a6f51b09…, yes | 293, 16863, 9345, 489, 1; 0 | VERIFIED 25.6 s | VERIFIED UNSAT 11.5 s |
| q455 | 0e06b849…, yes | 3 553 814 / 2c2738b1…, yes | 175, 9114, 4312, 91, 1; 0 | VERIFIED 4.2 s | VERIFIED UNSAT 1.6 s |

The repository had run cake_lpr only on its second encoding; these runs on the stored proofs are new.

## 4. Vertex-criticality — confirmed

`referee_critical.py` checks, for every vertex `v`, that `c_v(v) = −1`, that `c_v` is a `(7, 2)`-colouring of `H − v`,
finds a topological numbering `pos: V − v → {0, …, N − 1}` of its tight digraph (Kahn's algorithm, `N = n − 1`), and
checks on every edge that `f = N·c_v + pos (mod 7N)` is a `(7N, 2N + 1)`-colouring: with `δ = c_v(y) − c_v(x) mod 7`
and `e = pos(y) − pos(x)`, `|e| ≤ N − 1`, `f(y) − f(x) ≡ Nδ + e`; for `δ ∈ {3, 4}` this lies in `[2N + 1, 5N − 1]`;
for `δ = 2` the arc `x → y` is tight, so `e ≥ 1` and `Nδ + e ∈ [2N + 1, 3N − 1]`; for `δ = 5` the arc `y → x` is
tight, so `e ≤ −1` and `Nδ + e ∈ [4N + 1, 5N − 1]`. So `χ_c(H − v) ≤ 7N/(2N + 1) < 7/2`, explicitly. All 170, 293 and
175 certificates pass: `χ_c(H − v) ≤ 1183/339`, `2044/585` and `1218/349` (`results/critical.log`). Negative controls
are rejected. A proper induced subgraph `H'` misses some `v`, so `χ_c(H') ≤ χ_c(H − v) < 7/2`.

## 5. Lemma 20, the arithmetic, and the text

- **Lemma 20** re-derived: correct and correctly applied (if `χ_c(H) < p/q`, a homomorphism to `K_{p'/q'}` with
  `p'/q' < p/q` gives, after scaling and rounding, a `(p, q)`-colouring whose arcs have lifted differences strictly
  between `q` and `p − q` except possibly at `q`, and the sum along a closed walk shows that no cycle is tight).
- **Arithmetic** (`numtheory_and_seeds.py`, `results/numtheory.log`): `11 ≡ 191 ≡ 455 ≡ 11 (mod 12)`, so (a) and (b)
  fail; `(11/7) = (191/7) = 1`, so 7 splits in `ℚ(√11)` and `ℚ(√191)`, and 7 ramifies in `ℚ(√455)`; `√11, √191 ∈ ℚ₇`
  (Hensel lifts checked modulo `7¹⁰`), `√455 ∉ ℚ₇`. The seed graphs (76, 96, 71 vertices; `D = 30, 240, 780`) are not
  3-colourable, and 50, 75 and 50 of their points lie in the witnesses.
- **Section 10**: confirmed, apart from F2–F4, F9, F10 below; the statements about `ℚ₇` are correct (a field embedding
  `ℚ(√11) → ℚ₇ ⊆ K` preserves distance 1 in both directions, and by Corollary 6, `χ_c(K²) < 4` forces residue degree
  1 and the value 7/2).

## 6. Repository tests

`DRAT_TRIM=… python3 -m pytest -q tests/test_two_primes.py -k finite_witness`: 10 passed, 20 deselected (258 s),
including the slow proof checks.

## 7. Reproduction

- **q455**, with the programs at 5cc8ae84 and at c4811adc: `grow.py q455.json A G 12000 200` (71 rounds, 959
  vertices, 28 cycles), `minimise.py`, `critical.py`: points, edges, colouring, cycles and all 175 certificates
  identical to the stored files (`results/repro_q455_compare.log`). Growth rules used: "blocked" in 50 rounds (506
  points), "one colour left" in 3 rounds (382 points).
- **q11**, programs at 5cc8ae84: `grow.py q11.json A G11 3000 200` (52 rounds, 653 vertices), then the default
  `minimise.py` and `critical.py`: `witness_q11` and its 170 certificates exactly (`results/repro_q11_compare.log`).
  All 21 growth rounds used the "blocked" rule.

## 8. Findings (all applied since)

- **F1, F2** (README, paper): "Theorem F gives 7/2 as 7 ramifies/splits" omits that (a) and (b) must fail (they do,
  `d ≡ 11 mod 12`); better: the witness gives `≥ 7/2` and Proposition 3 gives `≤ 7/2`, without the computer-assisted
  Proposition 9.
- **F3** (README, paper): `grow.py` also adds, when no candidate is blocked, the candidates with exactly one colour
  left (382 of the 888 points added over `ℚ(√455)`).
- **F4** (paper): the denominator 30 is specific to `ℚ(√11)`; the other runs use 780 and 240.
- **F5** (checkers): `check_witness.py` and `verify_independent.py` relied on `assert`; under `python3 -O` a corrupted
  witness passed, and with an empty proof `check_witness.py` printed "s NOT VERIFIED" followed by the conclusion
  (`results/demo_O_before_fix.out`).
- **F6**: `check_witness.py` left `witness_check_tmp.cnf` in the current directory.
- **F7**: `grow.py`'s docstring named a nonexistent `witness_q11grow`.
- **F8** (README): reproduction commands were given only for q455 (q11 only for the growth), none for q191.
- **F9** (paper): the "two encodings" are the same direct encoding with different numbering and arc indicators.
- **F10** (paper): state the perturbation explicitly, `f = N·c + pos (mod 7N)`.
- **F11**: README wording ("written by the session that found …"); the JSON descriptions cite `witness_*.cnf` for the
  stored `witness_*.cnf.gz` and `witness_*.drat.xz`.
- **F12** (optional): `witness_q455` is a witness for the finite extensions of `ℚ₇` with residue degree 1 that contain
  `√455` (for example the ramified `ℚ₇(√455)`); Section 11 could point to `verification.txt`.

## 9. Files

The programs are in this folder; their outputs are in `results/`. The formulas and proofs they write (several hundred
MB) are not stored; the commands below regenerate them. With `FW=..` and `W = witness_q455` (likewise q11, q191):

    python3 referee_check.py $FW/$W.json.gz cnf
    python3 referee_encode2.py $FW/$W.json.gz cnf/${W}_E2.ref.cnf
    python3 referee_cegar.py $FW/$W.json.gz cnf/${W}_cegar.ref.cnf cnf/${W}_cegar_cycles.json 0
    bash run_proof.sh $W            # likewise ${W}_E2 and ${W}_cegar (KISSAT, DRAT_TRIM, CAKE_LPR from the environment)
    python3 validate_e1.py $FW/$W.json.gz cnf/$W.ref.cnf
    python3 validate_stored.py $FW/$W.json.gz $FW/$W.cnf.gz $FW/$W.drat.xz stored && bash run_stored.sh $W
    python3 referee_critical.py $FW/$W.json.gz $FW/${W}_critical.json.gz
    python3 mutation_tests.py $FW mut
    python3 numtheory_and_seeds.py ../../../../..
