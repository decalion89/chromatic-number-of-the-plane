# Explicit finite witnesses for `χ_c = 7/2` and `χ_c = 3`

Corollary 7 of `papers/three-colours/` (Corollary F12 of `notes/circular_planes.md`) says that below 4 the circular
chromatic number of the plane over a number field `F` is the circular chromatic number of a finite unit-distance graph
in `F²`; its proof, by compactness, gives no bound on the size. Here are explicit ones for `χ_c = 7/2`, over `ℚ(√11)`,
`ℚ(√191)` and `ℚ(√455)` (Theorem F gives `χ_c(ℚ(√d)²) = 7/2` for all three: some prime above 7 has residue degree 1),
and one with nine vertices for `χ_c = 3` over `ℚ(√7)`, whose plane has no unit triangle (last section).

**The claim.** Each graph `H` below is a unit-distance graph in `ℚ(√d)²`, induced (every pair of its points at
distance 1 is an edge), with `χ_c(H) = 7/2`.

- *Upper bound.* `H` has a `(7, 2)`-colouring, a homomorphism to `K_{7/2}` (stored in the file).
- *Lower bound.* Every `(7, 2)`-colouring of `H` has a *tight cycle* (a directed cycle along which each colour
  difference is `2 mod 7`), already among the listed cycles: the formula "a `(7, 2)`-colouring in which every listed
  cycle has a non-tight arc" is unsatisfiable. This is a SAT computation, certified (below). Lemma 20 of the paper (the
  easy half of a characterisation of Guichard) then gives `χ_c(H) ≥ 7/2`. For the grown witnesses the formula also
  fixes the colour of one vertex `v₀` (`fixed_vertex`, a vertex of largest degree) to 0, which loses nothing: rotating
  the colours keeps a `(7, 2)`-colouring and all its colour differences, hence its tight arcs; it makes the proofs
  about ten times shorter.

| witness | field | denominator | vertices | edges | listed cycles (lengths) | found by |
|---|---|---|---|---|---|---|
| `witness_q11` | `ℚ(√11)` | 30 | 170 | 468 | 879 (7, 14, 21, 28) | growth from the 76-vertex graph to 653 vertices, then vertex deletion |
| `witness_q191` | `ℚ(√191)` | 240 | 293 | 803 | 489 (7 to 35) | growth from the 96-vertex graph to 3 258 vertices, then vertex deletion |
| `witness_q455` | `ℚ(√455)` | 780 | 175 | 434 | 91 (7, 14, 28) | growth from the 71-vertex graph to 959 vertices, then vertex deletion |
| `witness_q11sum` | `ℚ(√11)` | 30 | 2 237 | 11 300 | 180 (14 to 42) | the sumset `A + A` of the 76-vertex graph (the first one found) |

The three grown witnesses are also *vertex-critical*: for every vertex `v`, `witness_*_critical.json.gz` gives a
`(7, 2)`-colouring of `H − v` whose tight digraph has no directed cycle, so `χ_c(H − v) < 7/2` (the other half of
Guichard's characterisation: perturb the colours along a topological order of the tight digraph), and no proper
induced subgraph of `H` is a witness; `check_critical.py` checks them. As `√11 ∈ ℚ₇`, `witness_q11` is also an explicit
witness for `ℚ₇` (Theorem C) and for every finite extension `K` of `ℚ₇` with `χ_c(K²) < 4`, and so is `witness_q191`
(`√191 ∈ ℚ₇`).

**How they were found.** `grow.py` (colouring-guided growth with lazy SAT): start from the vertex set `A` of the
4-chromatic graph `data/quadratic_planes/q<d>.json`; ask a SAT solver for a `(7, 2)`-colouring in which every listed
directed cycle has a non-tight arc; if the colouring has tight cycles, list them (a shortest one through a vertex of
each nontrivial strong component of the tight digraph) and ask again; if its tight digraph is acyclic, add the points
`p = x + u` (`x` a vertex, `u` a unit vector with the same denominator) whose neighbours leave `p` no colour, most
neighbours first, at most 200 a round. When the solver finds no such colouring, `H` is a witness. `minimise.py` then
deletes vertices one at a time, lowest degree first, while the formula stays unsatisfiable (one selector literal per
vertex in an incremental solver; a vertex outside the unsatisfiable core is dropped at once; the tight cycles it meets
are added to the list), and `critical.py` writes the criticality certificates. With python-sat 1.9 (CaDiCaL 1.5.3),
run in this folder,

    python3 grow.py ../../../quadratic_planes/q455.json A G 12000 200   # 71 rounds, 959 vertices: G.json
    python3 minimise.py G.json W.json                                    # 175 vertices
    python3 critical.py W.json C.json

reproduce the points, edges, colouring and cycles of `witness_q455` and its criticality certificates exactly, and
`q11.json` with `3000 200` gives the 653-vertex graph from which `witness_q11` was cut. (`fixed_vertex` is added
afterwards: a vertex of largest degree, the first one.) The three programs take any value `p/q` (`grow.py ... p q`;
the others read `p` and `q` from the file, `7/2` when absent), with the same results for `7/2`.

**Verification** (`verification.txt`): for every witness, two checks that share no code.
1. `check_witness.py` (written by the session that found `witness_q11sum`; it reads `d` from the file): exact
   integer checks of the points, of every edge and of all unit pairs, of the colouring and of the cycles; rebuilds
   the stored formula byte for byte; runs `drat-trim` on the stored DRAT proof: `s VERIFIED`.
2. `verify_independent.py` (written separately, from the file format only; it reads `d` and the denominator from the
   file): compares the points with `A`, recomputes all unit-distance pairs from the unit vectors with that
   denominator, checks the colouring and the cycles, and writes the formula with its own variable layout. `kissat`
   refutes it with a new DRAT proof, `drat-trim` verifies that proof and converts it to LRAT, and `cake_lpr` (the
   formally verified checker of CakeML) accepts the LRAT proof.

The stored proofs of the grown witnesses are the core proofs that `drat-trim -l` extracts from `kissat`'s proofs,
verified again by `check_witness.py`. To recheck one witness (`W` = `witness_q11`, `witness_q191`, `witness_q455`):

    xz -dk W.drat.xz
    python3 check_witness.py W.json.gz W.cnf.gz W.drat /path/to/drat-trim
    python3 check_critical.py W.json.gz W_critical.json.gz
    python3 verify_independent.py H.cnf W.json.gz   # then kissat H.cnf H.drat; drat-trim H.cnf H.drat -L H.lrat; cake_lpr H.cnf H.lrat

For `witness_q11sum` the commands are the same without `check_critical.py` (and `verify_independent.py H.cnf` alone
uses it by default). The tests are in `tests/test_two_primes.py` (`test_finite_witness_q11`,
`test_finite_witness_grown`, `test_finite_witness_critical_rejects`; the proof checks are marked slow).

| file | content |
|---|---|
| `witness_*.json.gz` | `d`, the denominator, the points `[a, b, c, e]` (the point `((a + b√d)/D, (c + e√d)/D)`), the edges, the `(7, 2)`-colouring, the cycles, and for the grown witnesses `fixed_vertex` |
| `witness_*.cnf.gz` | the formula of `check_witness.py` (sha256 of the uncompressed text in `verification.txt`) |
| `witness_*.drat.xz` | its DRAT proof of unsatisfiability |
| `witness_*_critical.json.gz` | for every vertex `v`, a `(7, 2)`-colouring of `H − v` with an acyclic tight digraph (`v` marked `−1`) |
| `check_witness.py` | the first checker (see its docstring for the encoding and the logic) |
| `verify_independent.py` | the second checker and the second encoding |
| `check_critical.py` | the checker of the criticality certificates |
| `grow.py`, `minimise.py`, `critical.py` | the programs that found the grown witnesses and their certificates |
| `witness_q7.json.gz` | the nine-point witness for `χ_c = 3` over `ℚ(√7)`: `d`, the denominator, `p = 3`, `q = 1`, the points, the edges, a 3-colouring, the cycles and the criticality certificates |
| `q7_seed.json` | the 207 points from which `grow.py` found it |
| `check_small.py` | the exhaustive checker for small witnesses (any `p/q`; no solver) |
| `small_triangle_free.py` | every triangle-free graph with at most 8 vertices maps to `K_{8/3}` |
| `verification.txt` | the outputs of both checks for every witness |

## The value 3 over `ℚ(√7)`: nine points

`ℚ(√7)²` has no unit triangle, as `√3 ∉ ℚ(√7)`, and `χ_c(ℚ(√7)²) = 3` (Theorem B and Corollary B6 of the paper).
`witness_q7.json.gz` is a unit-distance graph `H₇` on nine points of `ℚ(√7)²` (denominator 160), induced, with
`χ_c(H₇) = 3`; its vertices `0, …, 7, 8` are the points `P₀, …, P₇, S` of the paper:

| vertex | point | vertex | point |
|---|---|---|---|
| `P₀` | `(−1, 0)` | `P₅` | `((3 − √7)/8, (√7 − 5)/8)` |
| `P₁` | `(−(3 + √7)/8, −(5 + √7)/8)` | `P₆` | `(1, 0)` |
| `P₂` | `(1/4, −√7/4)` | `P₇` | `(0, 0)` |
| `P₃` | `(√7/4, 1/4)` | `S` | `(0, 1)` |
| `P₄` | `(−1/4, √7/4)` | | |

Its 13 edges, all the pairs at distance 1, are the 8-cycle `P₀P₁⋯P₇`, the chords `P₀P₄`, `P₁P₅`, `P₂P₆` and the path
`P₃SP₇`: `H₇` is the Wagner graph (the Möbius ladder on 8 vertices) with one of its four chords subdivided.

- *Upper bound.* A proper 3-colouring is stored.
- *Lower bound.* Each of the 84 proper 3-colourings has a tight cycle (a directed cycle along which the colour
  increases by 1 mod 3), already one of the six listed directed 6-cycles, so `χ_c(H₇) ≥ 3` by Lemma 20. Independently,
  `H₇` has no homomorphism to `K_{8/3}`, which is the Wagner graph itself; as the circular chromatic number of a graph
  with 9 vertices is a fraction with numerator at most 9 (Zhu's survey), and `8/3` is the largest such fraction below 3,
  this gives `χ_c(H₇) ≥ 3` again.
- *Vertex-critical.* For every vertex `v` a 3-colouring of `H₇ − v` with an acyclic tight digraph is stored.
- *Smallest possible.* Every triangle-free graph with at most 8 vertices maps to `K_{8/3}`, so a graph with `χ_c = 3`
  in a plane without unit triangles has at least 9 vertices.
- As `√7 ∈ ℚ₃`, `H₇` is also a witness for `ℚ₃` and for every finite extension `K` of `ℚ₃` with `χ_c(K²) < 4`.

The checks are finite enumerations (`verification.txt`):

    python3 check_small.py witness_q7.json.gz q7.cnf   # points, unit pairs, all 3^9 maps, K_{8/3}, criticality
    python3 small_triangle_free.py                     # the 4 682 270 triangle-free graphs on 8 labelled vertices

`check_small.py` also writes the formula "a 3-colouring in which every listed cycle has a non-tight arc" (the layout of
`check_witness.py`, with `(p, q) = (3, 1)`); kissat refutes it, and drat-trim and `cake_lpr` verify the proof.

**How it was found.** From `q7_seed.json` (0, the 204 unit vectors with denominator 160, and `(−2, 0)` and
`((−11 − √7)/8, (−5 − √7)/8)`, which close the 5-cycle `0, (−1, 0), (−2, 0), ((−11 − √7)/8, (−5 − √7)/8),
(−3/4, −√7/4)`),

    python3 grow.py q7_seed.json A G7 2000 200 3 1   # 2 growth rounds, 607 vertices, 7 cycles
    python3 minimise.py G7.json W7.json              # 9 vertices
    python3 critical.py W7.json C7.json

give the nine points translated by `(−1, 0)` and in another order, with the same edges, colouring, cycles and
certificates under that correspondence; the stored file uses the labels of the paper.
