# Referee report: the two Cayley-subgraph witnesses at 4

This folder holds an independent referee's programs, logs and report for `../witness_q3_11_cayley.json.gz` and
`../witness_q2_3_cayley.json.gz`. The referee read only the data files and a description of the claims, and wrote
every program here from scratch. It did not open `check_witness4.py` or the construction scripts. It checked the files
of the parent folder, byte for byte the same; the hashes are in §1 and `results/data_sha256.txt`. The section
"Rerun from the repository copy" at the end was added afterwards.

## 0. Scope and independence

The referee read only the six data files and the task description. Its programs are:

- `biquad.py`: exact arithmetic in the two fields. `test_biquad.py` is its self-test.
- `check_geom.py`: claims 1–5 and the count of unit-distance pairs.
- `show_construction.py`: rebuilds the generators and prints each next to its match.
- `make_cnf.py`: the referee's own encoding.
- `compare_cnf.py`: compares that encoding with the stored formula.
- `encoding_test.py`: tests what the CNF files say.
- `model_check.py`: decodes solver models and checks the stored colouring.
- `sanity_tests.py`: mutation tests.
- `unit_numpy.py`: a second exact count of unit-distance pairs.
- `gen_usage.py`, `extra_units.py`, `unit_colour_check.py`: further statistics.

Every geometric decision uses exact integer or `Fraction` arithmetic. Floating point serves only as a cross-check, and
all cross-checks agree. Tools: kissat 4.0.4 (sha256 `cc4aa4ba…0f79`), drat-trim (`fcc6f026…d9d7`), `cake_lpr`
(`b25551c8…c8ca`).

## 1. Input hashes (sha256)

```
b24950821bd8ddd8951cb162a51108b11f0273163addde6235e2da7472ef64fd  witness_q3_11_cayley.json.gz
e5f201603e6ac470749b371a9b379ca3147a4dee5acb0f6ee3d1f42a88d0fd67  witness_q3_11_cayley.cnf.gz
d8f7b0c700626c257f1ea46554bd16a11ae9d27a6414fb4f6cecec0923505d07  witness_q3_11_cayley.drat (= xz -dc .drat.xz)
997b2a2c8130211d86542a3ed0cf6b102fb68190124f4f6c4b3538b4d3c63e2d  witness_q2_3_cayley.json.gz
cd3ada49739eee37cd5d6d7411239f323579ebd715eb26f6c5a47c83ee404ad5  witness_q2_3_cayley.cnf.gz
66a71e245a4ef937fbf136572c0a0a0860e40d3a44c32b153eccf6d480de02a0  witness_q2_3_cayley.drat (= xz -dc .drat.xz)
```

Both stored proofs are binary DRAT.

## 2. Method for claims 1–5

- **Representation.** Points and generators are integer 8-tuples over a common denominator `D`. The basis
  `{1, √A, √B, √AB}` is linearly independent over `ℚ`, so two tuples are equal exactly when the points are equal, and
  "length exactly 1" means `|d|² = (D², 0, 0, 0)` in the field.
- **Claim 1 (over `ℚ(√3, √11)`).** `R60 = (1/2, √3/2)`, `RA = (5/6, √11/6)` and `RG = (11/14, 5√3/14)` were built exactly
  and checked to have norm 1 (the conjugate is then the inverse). The 27 vectors `R60^j RA^k RG^l (1, 0)`, converted to
  integers over 84, were compared with the generators as sets of `±` pairs.
- **Claim 2 (over `ℚ(√2, √3)`).** `ζ = ((√6 + √2)/4, (√6 − √2)/4)` and `w = (1/3, 2√2/3)` were built exactly and checked:
  `|ζ| = |w| = 1`, `ζ² = (√3/2, 1/2)`, `ζ¹² = −1`, `ζ²⁴ = 1`, and `ζ` is a primitive 24th root of unity. The 120 vectors
  `ζ^j w^l` (`j ∈ ℤ/24`, `l ∈ {−2, …, 2}`) were compared with the generators.
- **Claim 3.** Cayley pairs were found by hashing `p + g` for every point `p` and signed generator `g`, and compared
  with the listed edges in both directions; an independent scan of all pairs tested whether each difference is `±` a
  generator. A breadth-first search from the fixed vertex (checked to be the origin) tested that every point is
  reached. Points are distinct, and there are no repeated edges, loops or indices out of range.
- **Claim 4.** The colouring has the right length, values in `{0, 1, 2, 3}`, and no edge with both ends alike.
- **Claim 5.** Each cycle is nonempty, has length divisible by 4 and distinct vertices, and every cyclic step,
  including the step back to the start, is an edge.

## 3. Geometry results (exact)

| | `ℚ(√3, √11)` | `ℚ(√2, √3)` |
|---|---|---|
| `D` | 84 | 36 |
| distinct points | 1 874 | 1 657 |
| edges | 7 887 | 6 199 |
| generators / distinct `±` pairs | 27 / 27 | 60 / 60 |
| vectors built / distinct `±` pairs | 27 / 27 | 120 / 60 (closed under negation) |
| mismatches with the construction | 0: claim 1 holds | 0: claim 2 holds |
| generators of length exactly 1 | 27 of 27 | 60 of 60 |
| Cayley pairs found | 7 887 | 6 199 |
| edges that are not Cayley pairs / Cayley pairs not listed | 0 / 0 | 0 / 0 |
| the all-pairs scan agrees | yes (7 887) | yes (6 199) |
| fixed vertex; is it the origin | 953; yes | 725; yes |
| points reached from it | 1 874 of 1 874 | 1 657 of 1 657 |
| least / largest degree | 2 / 36 | 2 / 71 |
| edges with both ends alike (claim 4) | 0 | 0 |
| sizes of the colour classes | 504, 460, 447, 463 | 331, 365, 489, 472 |
| cycles (of length 4 / of length 8) | 3 389 (3 370 / 19) | 5 264 (4 705 / 559) |
| defective cycles (claim 5) | 0 | 0 |
| cycles also listed in reverse | 3 006 | 4 986 |

Claims 1–5 hold for both files.

## 4. Claim 6: every proper 4-colouring has a tight listed cycle

**The encoding.** The variable `x(v, c) = c·n + v + 1` says "vertex `v` has colour `c`" (colour-major, deliberately
unlike the stored numbering). The clauses are:

1. each vertex has at least one colour;
2. each vertex has at most one colour (pairwise exclusions);
3. for each edge and each colour `c`, the two ends are not both `c`;
4. the unit clause `x(fixed, 0)`;
5. for each listed cycle and each shift `s ∈ {0, …, 3}`, the clause `OR_t ¬x(v_t, (s + t) mod 4)`: the cycle is not tight
   in its listed direction.

**Sizes and hashes.**

| CNF | variables | clauses | sha256 |
|---|---|---|---|
| `ref_q3_11.cnf` | 7 496 | 58 223 | `2528ba27e21f1769f9f19f409dbb7f02e6f48edc806ad027b43db0c0f0d8bf2b` |
| `ref_q2_3.cnf` | 6 628 | 57 452 | `7987ccdbd5587bfda27b3e2a3d29d39495453bade329d8d61bbf1962babceeeb` |
| `ref_q3_11_nocycles.cnf` | 7 496 | 44 667 | `a0999450d86820db3f5f2bf64d50bf8ac3ac8f9ad36c51a848f61a649c171ff1` |
| `ref_q2_3_nocycles.cnf` | 6 628 | 36 396 | `e3dbec73c0d82b00b105eb862a0e58b24cce838156c998fc1f77cf03ab9fcf6d` |
| the stored `q3_11` formula, decompressed | 7 496 | 58 223 | `b394ee93…1230` |
| the stored `q2_3` formula, decompressed | 6 628 | 57 452 | `802921ac…b85b` |

**Comparison with the stored formulas.** After renaming the variables to `4v + c + 1`, the referee's formula and the
stored one are the same set of clauses for both files, with no difference either way.

**What the formulas say.** On 386 colourings per file (the stored colouring, a solver model, all 24 permutations of
the colours, recolourings and random colourings), the number of clauses a colouring breaks equals the number of edges
with both ends alike, plus 1 if `c(fixed) ≠ 0`, plus the number of listed cycles that are tight in their listed
direction: no disagreement in any of the four formulas. An assignment giving a vertex no colour, or two, breaks one
more clause.

**Solvers and proofs on the referee's formulas.**

| step | `q3_11` | `q2_3` |
|---|---|---|
| kissat | `s UNSATISFIABLE`, exit 20, 12.3 s | `s UNSATISFIABLE`, exit 20, 22.3 s |
| `drat-trim -L` | `s VERIFIED` (core of 30 528 clauses, 118 989 lemmas, 584 RAT), 11.2 s | `s VERIFIED` (core of 33 849, 176 576 lemmas, 272 RAT), 28.0 s |
| `cake_lpr` | `s VERIFIED UNSAT`, 24.6 s | `s VERIFIED UNSAT`, 25.3 s |

drat-trim printed no warning. The hashes of the referee's proofs are in `results/proof_sha256.txt`.

**The stored proofs on the stored formulas.** drat-trim prints `s VERIFIED` for both (`q3_11`: 20.4 s, 822 RAT lemmas in
the core; `q2_3`: 16.2 s, 253). As an extra check, `cake_lpr` accepts the LRAT proofs that drat-trim produced from each
stored proof (`s VERIFIED UNSAT`).

**The checkers reject bad input.** `cake_lpr` on the formula without cycle clauses: "clause index unavailable: 44668",
the first cycle clause; on a proof with its last step removed: "empty clause not derived"; on a proof whose last step
lost a hint: "Checking failed". drat-trim on the formula without cycle clauses: `s NOT VERIFIED`. As `cake_lpr` exits
with code 0 even when it rejects, success was judged only by the line `s VERIFIED UNSAT`.

**Why fixing `c(fixed) = 0` loses nothing.** Subtracting `c(fixed)` from every colour modulo 4 keeps a colouring proper
and keeps every colour difference, so tight cycles stay tight. So the refutation covers every proper 4-colouring, and
claim 6 holds for both files.

## 5. Claim 7: `χ_c(H) = χ(H) = 4` (Guichard's argument)

Suppose that `χ_c(H) < 4`, so that `H` has a circular `r`-colouring `φ` with `r < 4`: adjacent vertices satisfy
`1 ≤ |φ(u) − φ(v)| ≤ r − 1` (distances on the circle of length `r`).

1. Let `c = ⌊φ⌋`. This is a proper 4-colouring.
2. Along a tight cycle of `c`, each step moves `φ` forward (modulo `r`) by some `δ_t` with `1 ≤ δ_t < 2`.
3. `φ` wraps around exactly when `c` goes from 3 to 0. If this happens `W` times, then `m = 4W`.
4. So `m ≤ Σ δ_t = W·r < 4W = m`, which is impossible: `c` has no tight cycle.

This contradicts claim 6, so `χ_c(H) ≥ 4`, and claim 4 gives `χ_c(H) ≤ χ(H) ≤ 4`. Hence `χ_c(H) = χ(H) = 4` for both
files. (Directly, `χ(H) ≠ 3`: a tight cycle uses all four colours, so a 3-colouring has none.)

## 6. Claim 8: sanity tests

**Mutations: 29 of 29 rejected for each file.** Moved points (several kinds, among them the origin, and a point moved
onto another); an added edge that is not a Cayley pair, among them one at distance exactly 1; a loop; removed Cayley
edges; improper colourings, a colour 4 and a short list of colours; cycles through a non-edge, a repeated cycle, an empty
cycle and a closed walk of length 6 with repeats; corrupted generators (wrong length, a wrong unit vector, one dropped);
an isolated point; a wrong index for the fixed vertex.

**Without the cycle clauses the formula is satisfiable** (`s SATISFIABLE`, exit 10, for both). The decoded models are
proper 4-colourings with `c(fixed) = 0`, with 202 and 321 tight listed cycles, as claim 6 requires; the stored colourings
have 227 and 332.

## 7. Unit-distance pairs that are not edges (exact)

| | `q3_11` | `q2_3` |
|---|---|---|
| pairs at distance exactly 1 | 8 085 | 6 238 |
| of which edges | 7 887 | 6 199 |
| unit-distance pairs that are not edges | 198 | 39 |
| distinct differences up to sign | 12 | 17 |

Two independent exact implementations agree on these counts, and so does the floating-point cross-check. The extra
differences are not of the form `±R60^j RA^k RG^l` with `|j|, |k|, |l| ≤ 4`, nor `±ζ^j w^l` with `|l| ≤ 6`; some are half
angles, for example `(√33/6, √3/6)` at `θ_A/2` and `(√6/3, √3/3)` at `arccos(1/3)/2`. The stored colouring is proper on
these pairs too, so the full unit-distance graph on each point set also has `χ = χ_c = 4`.

## 8. Discrepancies and remarks

No stated claim fails. To be aware of:

1. Not every generator is used: no edge uses 9 of the 27 generators (`q3_11`) or 24 of the 60 (`q2_3`), so only 18 and 36
   directions occur. This agrees with claim 3, but wording that implied that every direction is used would be wrong.
2. `H` is not the full unit-distance graph on its points; it misses 198 and 39 unit-distance pairs. This is harmless.
3. Many cycles are also listed in reverse; only the listed direction was encoded, as described.
4. The `.drat.xz` files decompress to the `.drat` files byte for byte.

**Verdicts.** `witness_q3_11_cayley`: accept. `witness_q2_3_cayley`: accept.

## Rerun from the repository copy

The referee's programs, run again from this folder on the files of the parent folder (`results/repo_rerun.log`):

- `test_biquad.py` passes.
- `check_geom.py ../witness_q3_11_cayley.json.gz` and `check_geom.py ../witness_q2_3_cayley.json.gz` print
  `OVERALL: PASS`.
- `make_cnf.py` writes the same two formulas: sha256 `2528ba27…bf2b` and `7987ccdb…ceeeb`, as in §4.
- kissat 4.0.4 refutes them (12 and 23 seconds), and `drat-trim -L` prints `s VERIFIED`.
- `cake_lpr` prints `s VERIFIED UNSAT` for both LRAT proofs.

Usage: `python3 check_geom.py WITNESS.json.gz [--no-unit] [--json OUT]`, `python3 make_cnf.py WITNESS.json.gz OUT.cnf
[--no-cycles] [--map OUT.map.json]`, `python3 compare_cnf.py REF.cnf STORED.cnf N`, `python3 sanity_tests.py
WITNESS.json.gz`; the other programs give their usage in their docstrings.
