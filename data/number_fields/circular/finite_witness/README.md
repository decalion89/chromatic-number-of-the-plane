# An explicit finite witness for `χ_c(ℚ(√11)²) = 7/2`

Corollary 7 of `papers/three-colours/` (Corollary F12 of `notes/circular_planes.md`) says that some finite
unit-distance graph over `ℚ(√11)` has circular chromatic number `χ_c(ℚ(√11)²) = 7/2`, but its proof, by compactness,
gives no bound on its size. Here is one.

**The graph.** Let `A` be the vertex set of the 76-vertex unit-distance graph of `data/quadratic_planes/q11.json`
(the graph that is not 3-colourable). `H` is the unit-distance graph on the sumset `A + A = {a + a′}`: 2 237 points
`((a + b√11)/30, (c + e√11)/30)` and all 11 300 pairs of them at distance exactly 1.

**The claim.** `χ_c(H) = 7/2`.

- *Upper bound.* `H` has a `(7, 2)`-colouring, that is, a homomorphism to `K_{7/2}`: the 7-adic colouring
  `λ(z mod w)`, `λ(x + iy) = 2x + 3y`, with `√11 ≡ 2` modulo the place `w` above 7 (stored in the witness file).
- *Lower bound.* Every `(7, 2)`-colouring of `H` has a *tight cycle* (a directed cycle along which each colour
  difference is `2 mod 7`), already among 180 listed cycles of lengths 14 to 42. This is a SAT computation, certified:
  the formula "a `(7, 2)`-colouring in which every listed cycle has a non-tight arc" is unsatisfiable. Then
  Lemma 20 of the paper (the easy half of a characterisation of Guichard) gives `χ_c(H) ≥ 7/2`.

So the graph attains the value of the whole plane. As `√11 ∈ ℚ₇`, it is also an explicit witness for
`χ_c(ℚ₇²) = 7/2` (Theorem C) and for every finite extension `K` of `ℚ₇` with `χ_c(K²) < 4`. For comparison (Section 10 of the paper): the 76-vertex graph
itself has `χ_c = 16/5` (SAT, not certified), the 628-vertex union of nine rotated copies has `χ_c ∈ (16/5, 13/4]`,
and the ball of radius 2 in the Cayley graph of the 140 vectors of Theorem C (9 941 vertices) has `χ_c = 5/2`.

**Verification** (`verification.txt`): two checks that share no code.
1. `check_witness.py` (written by the session that found the graph): exact integer checks of the points, of every
   edge and of the colouring and the cycles; rebuilds the stored CNF byte for byte; runs `drat-trim` on the stored
   DRAT proof: `s VERIFIED`.
2. `verify_independent.py` (written separately, from the file format only): rebuilds `A + A` from `q11.json`,
   recomputes all unit-distance pairs from the 108 unit vectors with denominator 30, checks the colouring and the
   cycles, and writes the formula with its own variable layout. `kissat` refutes it with a new DRAT proof, `drat-trim`
   verifies that proof and converts it to LRAT, and `cake_lpr` (the formally verified checker of CakeML) accepts the
   LRAT proof.

To recheck:

    xz -dk witness_q11sum.drat.xz
    python3 check_witness.py witness_q11sum.json.gz witness_q11sum.cnf.gz witness_q11sum.drat /path/to/drat-trim
    python3 verify_independent.py H.cnf   # then kissat H.cnf H.drat; drat-trim H.cnf H.drat -L H.lrat; cake_lpr H.cnf H.lrat

The witness was found by iterated Minkowski sums of the 76-vertex graph, tested by lazy SAT (tight cycles added to the
formula on demand); a search for a smaller witness inside `H` continues. The tests are in `tests/test_two_primes.py`
(`test_finite_witness_q11`; the proof check is marked slow).

| file | content |
|---|---|
| `witness_q11sum.json.gz` | denominator 30, the 2 237 points `[a, b, c, e]`, the 11 300 edges, the `(7, 2)`-colouring, the 180 cycles |
| `witness_q11sum.cnf.gz` | the formula of `check_witness.py` (sha256 of the uncompressed text in `verification.txt`) |
| `witness_q11sum.drat.xz` | its DRAT proof of unsatisfiability |
| `check_witness.py` | the first checker (see its docstring for the encoding and the logic) |
| `verify_independent.py` | the second checker and the second encoding |
| `verification.txt` | the outputs of both checks |
