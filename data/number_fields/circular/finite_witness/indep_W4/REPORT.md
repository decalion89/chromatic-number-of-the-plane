# Referee report: a finite unit-distance graph over Q(sqrt3, sqrt11) with circular chromatic number 4

**Verdict: ACCEPTED.** chi_c(H) = 4, and so chi(H) = 4.

Witness checked: `../witness_q3_11.json.gz` (at the time of the check `fin_j5.json.gz`, sha256
e592a0b13dfcecfc7fb09f43ae654b8661c3906cd783a33607a4267d0854addb).

## How to rerun

- `python3 ref_check4.py WITNESS.json.gz OUTDIR` runs tasks 1-5. The last line is `REFEREE: ACCEPTED` or
  `REFEREE: REJECTED <reason>`. Options: `--no-sat`, `--allow-edge-mismatch`, `--keep-proofs`, `--timeout`, and the
  paths of the three tools (`--kissat`, `--drat-trim`, `--cake-lpr`; by default the environment variables `KISSAT`,
  `DRAT_TRIM`, `CAKE_LPR`, else the names on the PATH).
- `python3 ref_tests4.py WITNESS.json.gz OUTDIR` runs the sanity tests.
- Both programs were written from scratch; no repository code is used.

## Numbers

- 1874 points, denominator 84, all distinct. The fixed vertex 953 is the origin.
- All 1 755 001 pairs were checked exactly: 8085 are at unit distance, identical to the declared `edges`.
- A floating-point pass found the same 8085 pairs; no other pair has squared distance within 1.76e-4 of 1.
- The colouring is proper, with class sizes 504/460/447/463.
- 4992 cycles (4916 of length 4, 52 of length 8, 24 of length 12). All are simple, all steps are edges, and all
  lengths are divisible by 4.
- CNF: 16 344 variables, 85 843 clauses (sha256 ffa6204cffe51aea3887cfdffa369461d86bb166e12e805e2188ff1e66c834fe).

## Tool output (final run, niceness 5, machine load about 5)

| Tool | Output | Exit | Wall | CPU | Proof file |
|---|---|---|---|---|---|
| kissat | `s UNSATISFIABLE` | 20 | 35.3 s | 23.6 s | 13.2 MB binary DRAT |
| drat-trim | `s VERIFIED` | 0 | 25.8 s | 16.8 s | 50.7 MB LRAT |
| cake_lpr | `s VERIFIED UNSAT` | 0 | 30.5 s | 20.3 s | |

drat-trim also reported 391 RAT lemmas in the core. An earlier identical run gave the same three lines. Both proof
files were deleted after the check.

## Task 1: points and edges

- All arithmetic is exact, on integer 4-tuples in the basis (1, r3, r11, r33), with r3^2 = 3, r11^2 = 11, r33^2 = 33,
  r3 r11 = r33, r3 r33 = 3 r11, r11 r33 = 11 r3.
- The basis is independent over Q: if r11 = a + b r3 with a, b rational, then ab = 0, so 11 or 11/3 would be a
  rational square. So two points are equal only if their integer 8-tuples are equal, and a pair is at distance 1
  exactly when dx^2 + dy^2 = (D^2, 0, 0, 0).
- Every edge was checked a second time with a general multiplication routine. The declared field and basis strings
  are checked too.
- Degrees range from 2 to 36; no vertex is isolated.
- The comparison with the declared edges is strict by default: any difference rejects the witness.

## Task 2: colouring

Proper on the recomputed edges.

## Task 3: cycles

Every cycle is simple, every step (including the closing one) is an edge, and every length is divisible by 4. They use
8848 distinct directed arcs and touch 1800 of the 1874 vertices; 304 of them are tight under the given colouring.

## Task 4: CNF

- Variables: colour variables x(v,k) = k n + v + 1 (1..7496); tight-arc variables t(a,b) (7497..16344), one per
  directed arc used by a cycle.
- Clauses:

| Group | Clause | Count |
|---|---|---|
| (A) | at least one colour per vertex | 1874 |
| (B) | at most one colour per vertex | 11244 |
| (C) | adjacent vertices differ, per edge and colour | 32340 |
| (D) | x(953,0) | 1 |
| (E) | x(a,k) and x(b,k+1 mod 4) imply t(a,b) | 35392 |
| (F) | one per cycle: OR of not t(a,b) over its arcs | 4992 |

- Why fixing colour(953) = 0 loses nothing: let c' = c − c(953) mod 4. Then c'(y) − c'(x) = c(y) − c(x) for every
  pair of vertices, so c' is proper exactly when c is, and c' has exactly the same tight arcs and tight cycles as c,
  with c'(953) = 0. Test (e) below confirms it for every shift s = 1, 2, 3.

## Task 5: solver and checkers

Output lines and timings are in the table above. The CNF hash was the same before and after the tools ran.

## Task 6: sanity tests (all passed)

- (a1) One coordinate of vertex 6 (on a cycle) changed: `REJECTED cycle 0 uses non-edge 0 -> 6`.
- (a2) The same with the edge comparison turned off: rejected for the same reason.
- (a3) One coordinate of vertex 10 (on no cycle) changed: `REJECTED declared edges differ ... non-unit 6`.
- (b) Edge 0-1 given the same colour on both ends: `REJECTED colouring not proper`.
- (c) Cycle [0,6,141,13] changed to [0,2,141,13]: `REJECTED cycle 0 uses non-edge 0 -> 2`.
- (d) Without the cycle clauses kissat finds `s SATISFIABLE` in 3.7 s; the model satisfies every clause and decodes
  to a proper 4-colouring with colour(953) = 0, which makes 332 listed cycles tight.
- (e) The assignment taken from the given colouring breaks exactly the 304 clauses of its tight cycles, and nothing
  else.

## Task 7: the mathematics

- UNSAT implies that every proper 4-colouring has a tight listed cycle: if some proper 4-colouring made no listed
  cycle tight, rotate it so that colour(953) = 0 and set x(v,k) = [c(v) = k], t(a,b) = [c(b) − c(a) = 1 mod 4]; this
  satisfies (A)-(F) (for (F): a cycle that is not tight has an arc with t false), a contradiction.
- Each listed cycle is a directed simple cycle of H, so a tight one is a tight cycle in the sense of Guichard and Zhu.
  Using only the listed cycles makes the claim stronger, not weaker; (E) only forces t to be true on tight arcs, which
  is all the argument needs.
- chi_c(H) <= 4 because of the proper 4-colouring of task 2.
- Guichard's direction, proved here so that the citation is not the only support: suppose chi_c(H) = r < 4. Take a
  circular r-colouring φ with 1 <= |φ(x) − φ(y)| <= r − 1 on every edge and set c = floor(φ), a proper 4-colouring.
  The function g = φ − c never decreases along a tight arc, and strictly increases (by 4 − r) along an arc that wraps
  from colour 3 to colour 0; every tight closed walk contains such a wrapping arc. So c has no tight cycle, a
  contradiction. Hence chi_c(H) = 4, and chi(H) = 4 as well.
- Gaps: none in the mathematics. Trust is needed in Python, in these programs and in the cake_lpr binary (not rebuilt
  here). kissat and drat-trim need not be trusted. The choice of real embedding of sqrt3 and sqrt11 does not matter:
  every embedding gives the same graph.
- Minor: H is not minimal (some vertices have degree 2, and 74 vertices are on no listed cycle). This does not affect
  correctness.

## Rerun from the repository copy

`results/` holds the outputs of both programs run from this folder on `../witness_q3_11.json.gz` (tool paths from the
environment; absolute paths replaced by `OUT/` and `../`): `REFEREE: ACCEPTED` (kissat `s UNSATISFIABLE`, drat-trim
`s VERIFIED`, cake_lpr `s VERIFIED UNSAT`) and `TESTS: ALL PASSED`. The CNF's first comment lines name the witness
file, so its sha256 is 1f28389f2aec69a069146c82b7c1a7870a5955fa72a380e21a8d3335efcdcde6 for `witness_q3_11.json.gz`
(and ffa6204c… for the working name `fin_j5.json.gz` above); the clauses are the same.

## The lemma on base-point periods (Lemma F17 of the note)

The same referee then checked the lemma that explains how the witness was found (its proof of `χ_c = 4` does not use
it), in the form: for a proper 4-colouring without tight cycles of a graph whose edges join `x` and `x + s`,
`s ∈ U ∪ −U`, the period `per(W) = Λ(W)/4 − N(W)` of closed walks is (a) an integer, (b) additive with
`per(W⁻¹) = −per(W)`, (c) invariant under backtracks and square swaps, and (d) strictly between `(P − 3N)/4` and
`(3P − N)/4`. Verdict: correct once three things are added, all of which the construction already satisfies.

- **A missing hypothesis:** `U ∩ (−U) = ∅`, so that every step lies in exactly one of `U`, `−U`. Without it (b) and
  the backtrack part of (c) fail: for `U = {u, −u}`, `H = {0, u}`, `c(0) = 0`, `c(u) = 1`, the walk `W = (0, u, 0)`
  has `W⁻¹ = W` but `per(W) = 1 − 2 = −1`, and deleting the backtrack gives 0. The 27 vectors satisfy it.
- **(a), (b):** correct, for every proper 4-colouring. **(c):** correct; the backtrack part holds for every proper
  colouring, the swap part needs "no tight cycle"; the proviso `s ≠ ±t` can be dropped (`t = s` changes nothing,
  `t = −s` exchanges two backtracks). **(d):** correct as stated; `L = P + N ≥ 1` is needed. **Proof:** correct;
  revisited vertices need no extra hypothesis.
- **Corollary:** false for walks of length 0 if a range is imposed on them (`0 < y < 0`); the empty walk gets
  `y = 0` and ranges are imposed only when `P + N ≥ 1`, as `per_build.py` does.
- **Use:** correct with the chain steps in both directions (`ρ′ = ρ ± r_j`), `0 ∈ H`, every relation joined to 0 by
  chain steps realised inside the graph, and `(r_j)` a ℤ-basis of the relation lattice (here 23 relations of rank
  `23 = 27 − 4`, Smith normal form all 1s). The character claim: a character of `ℤU` with representatives
  `f(u) ∈ (1/4, 3/4)` gives the integer homomorphism `p(ρ) = Σ n_u f(u)`, strictly inside every range with `ρ ≠ 0`;
  so infeasibility implies `κ(U) ≤ 1/4` (the converse does not follow).
- **Numerical check** (`lemmaP_check.py`, logs `results/lemmaP_check_run*.log`): 120 random small instances in `ℤ²`
  and `ℤ³` with every proper 4-colouring listed; on the colourings without tight cycles, 13 800 closed walks (with
  revisits) and about 316 000 moves (11 625 with `t = ±s`), no failure of (a)-(d); colourings with tight cycles break
  the swap part of (c) and (d), and a `U` meeting `−U` breaks (b), (c) and (d). For the use paragraph, 29 instances
  and 1 394 chain steps, each bubble-sorted (19 126 intermediate walks): `per` stayed constant through every sort,
  `y_ρ = Σ a_j p_j` and every range held, and for colourings `⌊4·frac⟨a, x⟩⌋` from characters,
  `per(W_ρ) = Σ n_u f(u)` exactly.

The note and the paper state the lemma with `U ∩ (−U) = ∅`, the ranges for `ρ ≠ 0`, the chain steps in both
directions with a ℤ-basis, and "infeasibility implies `κ(U) ≤ 1/4`".

## Files

`ref_check4.py`, `ref_tests4.py`, `lemmaP_check.py`; `results/main_run.log`, `results/tests_run.log`,
`results/summary.json`, `results/kissat.log`, `results/drat_trim.log`, `results/cake_lpr.log`,
`results/lemmaP_check_run1.log`, `results/lemmaP_check_run2.log`.
