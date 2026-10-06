# Referee report: a finite unit-distance graph over Q(sqrt2, sqrt3) with circular chromatic number 4

**Verdict: ACCEPTED.** I tested all 1,371,996 point pairs exactly. The 6238 declared edges are exactly the unit-distance
pairs. The colouring is proper. All 6062 cycles are simple, have length 4 or 8 and use only edges. My own CNF for
"proper 4-colouring, origin coloured 0, no listed cycle tight" is UNSAT, and cake_lpr verified the proof
(`s VERIFIED UNSAT`). With the proofs in section 7, this gives χ_c(H) = χ(H) = 4.

All programs used here were written from scratch. I read no other checker or construction script.

## Files

- `ref_check_z.py`: tasks 1–5. `--no-sat` skips the solver runs. `--encoding direct` selects the second (cross-check)
  encoding.
- `ref_tests_z.py`: task 6.
- `results/`: logs.
  - `check_nosat.log`, `check_full.log`, `check_full_direct.log`
  - `kissat_{arc,direct}.log`, `drat_trim_{arc,direct}.log`, `cake_lpr_{arc,direct}.log`
  - `tests.log`, `exit_codes.log` (all four runs exited 0), and JSON summaries.
- `work/ref_z.cnf` and `work/ref_z_direct.cnf`: the two formulas.
- `tests/*.json.gz`: the corrupted witnesses.

Input: sha256 `c84d17f4c97e488b1bc82df7d43334714117717e2423d8a76f1a8f0c2d359c05`, which matches. n = 1657, D = 36,
largest |coefficient| 76, all coordinates are JSON integers.

## 1. Points and edges

**Multiplication.** Put s = √2, t = √3 and √6 = st. The rules s² = 2 and t² = 3, with commutativity, give s·t = √6,
s·√6 = 2√3, t·√6 = 3√2 and √6·√6 = 6. So:

xy = (x0y0 + 2x1y1 + 3x2y2 + 6x3y3) + (x0y1 + x1y0 + 3x2y3 + 3x3y2)√2 + (x0y2 + x2y0 + 2x1y3 + 2x3y1)√3
+ (x0y3 + x3y0 + x1y2 + x2y1)√6

The program builds this table from s² = 2 and t² = 3 alone; the table is not typed in. It checks the table against
the real roots (error 8.9e-16). It also proves that its fast square formula equals the generic product: both are
polynomials of degree ≤ 2 in each variable, and they agree on the grid {−1, 0, 1, 2}⁴.

**Linear independence over Q.**
- (a) √2, √3 and √6 are irrational. If √k = p/q in lowest terms with k squarefree, then any prime ℓ dividing k
  divides both p and q, which is a contradiction.
- (b) √3 is not in Q(√2). If √3 = a + b√2, then 3 = a² + 2b² + 2ab√2, so ab = 0. If b = 0, √3 is rational. If a = 0,
  √6 = 2b is rational.
- (c) Suppose a + b√2 + c√3 + d√6 = (a + b√2) + (c + d√2)√3 = 0. If c + d√2 ≠ 0, then √3 is in Q(√2), against (b).
  So c = d = 0, and then a = b = 0.

So coefficients are unique. Two points are equal exactly when their 8-tuples are equal, and |p − q| = 1 exactly when
D²|p − q|² has coefficients (1296, 0, 0, 0). Because this is a coefficient test, every choice of signs for √2 and √3
gives the same graph.

**Results**

| Check | Result |
|---|---|
| Distinct points | 1657 of 1657 |
| Pairs tested exactly | 1,371,996, i.e. all pairs |
| Unit-distance pairs found | 6238 |
| Declared edges | 6238, valid, no repeats, no self-loops |
| Declared but not at distance 1 | 0 |
| At distance 1 but not declared | 0 |
| Pairs whose rational part is 1296 but whose irrational part is not 0 | 4734 (so all four coefficients must be tested, and they are) |
| Float pass (\|d − 1\| < 1e-9 counts as an edge) | 6238 edges, 0 disagreements with the exact pass |
| Largest \|d − 1\| on an edge (float) | 7.8e-16 |
| Smallest \|d − 1\| on a non-edge (the margin) | 9.10e-6, pair (1263, 1600), 7 pairs tie |
| Same pair with 60-digit decimals | \|d\|² = 0.9999817967177472949…; exactly, D²\|d\|² = 3312 + 144√2 + 144√3 − 1008√6 |
| Non-edges with \|d − 1\| < 1e-3 / < 1e-6 | 532 / 0 |
| Smallest distance between two points | 0.00609 |
| Fixed vertex | index 725, which is the origin |
| Degrees | 2 to 71 (the origin has degree 71), mean 7.53 |

## 2. Colouring

- Colour counts: 0: 331, 1: 365, 2: 489, 3: 472.
- 0 monochromatic edges among the 6238 recomputed edges, so the colouring is proper.

## 3. Cycles

- 6062 cycles: 5494 of length 4 and 568 of length 8.
- Every cycle is simple, has length ≥ 4 and divisible by 4, and every step is an edge, including the closing step.
- The cycles are distinct up to rotation, and each one's reverse is also listed. So this is 3031 cycles, each in both
  directions.
- 1650 of the 1657 vertices lie on some cycle.
- 389 cycles are tight under the given colouring, all of length 4.

## 4. Encoding

Variables:
- x(v,k) means "vertex v has colour k".
- y(u,v) is one variable per directed arc on a listed cycle (8824 arcs).

Clauses of the main ("arc") encoding:

| Clause | Count |
|---|---|
| x(725,0) | 1 |
| At least one colour per vertex | 1657 |
| At most one colour per vertex | 9942 |
| ¬x(u,k) ∨ ¬x(v,k) for each edge and colour | 24952 |
| ¬x(u,k) ∨ ¬x(v,k+1 mod 4) ∨ y(u,v) | 35296 |
| ¬y(a_0) ∨ … ∨ ¬y(a_{m−1}) for each cycle, closing arc included | 6062 |

That gives **15452 variables and 77910 clauses**, sha256
`a8e9d4bdee31892b03d2f6bbe4e865f38414434cd72be38768e3a65bf2483553`.

The cross-check ("direct") encoding has no auxiliary variables. For each cycle and each start colour k it adds
¬x(v_0,k) ∨ ¬x(v_1,k+1) ∨ … ∨ ¬x(v_{m−1},k+m−1). It has 6628 variables and 60800 clauses, sha256
`e572a24243a926b977c07e16c15430e08bf68e930ae384dae744e958267ee754`.

The program reads each file back and confirms it equals the clauses in memory.

**Fixing a colour loses no generality.** Let c be a proper 4-colouring and set c' = c − c(f) mod 4.
- c' is proper, because c'(u) = c'(v) exactly when c(u) = c(v).
- c'(f) = 0.
- Differences mod 4 are unchanged, so a listed cycle is tight for c' exactly when it is tight for c.

Only a rotation of the colours is needed; a general permutation would not preserve tightness.

## 5. Solver runs (final run)

| Step | Output | Exit | Time |
|---|---|---|---|
| kissat ref_z.cnf ref_z.drat | s UNSATISFIABLE | 20 | 109.0 s |
| drat-trim ref_z.cnf ref_z.drat -L ref_z.lrat | s VERIFIED | 0 | 54.0 s |
| cake_lpr ref_z.cnf ref_z.lrat | s VERIFIED UNSAT | 0 | 27.1 s |

- Proof sizes: DRAT 44,709,809 bytes, LRAT 88,741,487 bytes. Both were deleted afterwards.
- The formula's sha256 was identical before kissat, after kissat, after drat-trim and after cake_lpr.
- An earlier run gave the same results and proof sizes (98.6 s, 57.5 s and 40.9 s).

The direct encoding gave the same outcome: kissat `s UNSATISFIABLE` (exit 20, 34.3 s), drat-trim `s VERIFIED`
(exit 0, 19.4 s), cake_lpr `s VERIFIED UNSAT` (exit 0, 7.6 s), with the hash unchanged.

## 6. Sanity tests (77 checks, 0 failures)

Each corrupted file was checked by `ref_check_z.py --no-sat` in its own process. A test passes only if the set of
reasons reported is exactly the expected set.

| Corrupted file | Reasons reported |
|---|---|
| Unchanged witness | none, accepted |
| Vertex 0 (on 2 cycles) moved by 1/36 | edge_mismatch, cycle_nonedge |
| Vertex 188 (on no cycle) moved by 1/36 | edge_mismatch |
| Vertex 188 moved by about 7.1e-7 | edge_mismatch |
| Improper colour | colouring_improper |
| Cycle with an inner non-edge | cycle_nonedge (step 1) |
| Cycle with a non-edge closing step | cycle_nonedge (step 3 of 4, the closing step) |
| Closed walk [0,1,0,1] | cycle_not_simple |
| Triangle listed as a cycle | cycle_length |
| One declared edge removed | edge_mismatch |
| Closest non-edge declared as an edge | edge_mismatch |
| Duplicated point | points_not_distinct, edge_mismatch |
| Fixed vertex not the origin | fixed_not_origin |
| A coordinate written as a float | format |
| A colour set to 4 | colouring_malformed |

Encoding tests:
- **B1.** Without the cycle clauses the formula is satisfiable. Ten runs (one plain, nine on randomly renamed and
  shuffled copies) gave 10 distinct models. Each decodes to a proper 4-colouring with the origin coloured 0, and each
  has 316 to 413 tight listed cycles, as the UNSAT result predicts.
- **B2.** The given colouring, rotated so the origin has colour 0, satisfies every non-cycle clause. In both encodings
  it violates exactly the clauses of its 389 tight cycles.
- **B3.** With those 389 cycles dropped, the formula is satisfiable, and its model has no tight remaining cycle.

## 7. Mathematics

**(a) UNSAT implies a tight listed cycle.** Let c be a proper 4-colouring with c(f) = 0, and suppose no listed cycle is
tight. Set x(v,k) = [c(v) = k] and y(u,v) = [c(v) − c(u) ≡ 1]. Then every clause holds:
- The unit clause holds because c(f) = 0.
- The at-least-one and at-most-one clauses hold because each vertex has one colour.
- The edge clauses hold because c is proper.
- An arc clause holds because a tight arc makes its y true.
- A cycle clause holds because a non-tight cycle has an arc whose y is false.

So the formula would be satisfiable, which contradicts UNSAT.

The direct encoding works the same way. If one of its clauses were false, then c(v_t) = k + t for all t. Every inner
step would be +1, and the closing step is 1 − m ≡ 1 because 4 divides m, so the cycle would be tight.

**(b) Tight cycles force χ_c ≥ 4.** A circular r-colouring is f: V → [0, r) with 1 ≤ |f(u) − f(v)| ≤ r − 1 on every
edge, and χ_c = inf r. A (k,d)-colouring gives one with r = k/d via f = c/d, so the k/d definition is covered too.

Suppose f is a circular r-colouring with r < 4.
1. Let σ = 4/r > 1 and g = σf. Then g takes values in [0, 4), and σ ≤ |g(u) − g(v)| ≤ 4 − σ on every edge.
2. Floor map: c = ⌊g⌋. This is a proper 4-colouring, because equal floors mean a difference below 1, which is below σ.
3. By hypothesis c has a tight cycle v_0, …, v_{m−1}. Write g(v_t) = c(v_t) + e_t with 0 ≤ e_t < 1, and let
   δ_t = (g(v_{t+1}) − g(v_t)) mod 4.
   - Each step is an edge, so δ_t ≥ σ.
   - The step c(v_{t+1}) − c(v_t) is 1 or −3, so δ_t = 1 + e_{t+1} − e_t. This number lies in (0, 2).
4. Summing, the e-terms cancel and Σδ_t = m. But Σδ_t ≥ mσ > m, which is a contradiction.

So no r < 4 works, and χ_c ≥ 4. This proves directly the direction of the cited 1993 characterisation that is needed.

**(c) Conclusion.**
- The given colouring is proper, so χ ≤ 4.
- A proper k-colouring is a circular k-colouring, so χ_c ≤ χ ≤ 4.
- The verified UNSAT result with (a) and (b) gives χ_c ≥ 4.

So **χ_c(H) = χ(H) = 4**. No step fails.

**Trusted:**
- Python 3.11.15: integers, json, gzip, hashlib, and the machine.
- `ref_check_z.py`: field arithmetic (self-tested), the all-pairs loop, the edge, colouring and cycle checks, and the
  CNF writer.
- The cake_lpr binary, as a faithful build of the verified checker including its parsers.
- The proofs above.

kissat and drat-trim are checked by cake_lpr and need no trust. numpy is used only for the float cross-check.

**Doubts:** none of substance.

## Rerun from the repository copy

The copies here differ from the referee's files only in where they look for things: `ref_check_z.py` takes the three
tools from the environment variables `KISSAT`, `DRAT_TRIM` and `CAKE_LPR` (else from the PATH) and reads the witness
from the parent folder; `ref_tests_z.py` writes its corrupted files to a temporary folder (or `$REF_TESTDIR`). Run
from this folder,

    python3 ref_check_z.py --cnf /tmp/ref_z.cnf --results /tmp/res     # tasks 1-5; --no-sat for tasks 1-4
    python3 ref_tests_z.py                                             # task 6

gave `RESULT: PASS (tasks 1-5)` (kissat `s UNSATISFIABLE` in 102 s, drat-trim `s VERIFIED` in 53 s, cake_lpr
`s VERIFIED UNSAT` in 24 s, the formula hash `a8e9d4bd…` unchanged throughout; `results/repo_run.log`) and
`77 checks, 0 failures` (`results/repo_tests.log`). The referee's own runs are the other files in `results/`. The
two formulas (`work/`) and the corrupted witnesses (`tests/`) of the list above are not stored here; the programs
write them again.
