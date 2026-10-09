# Explicit finite witnesses for `χ_c = 4`, `χ_c = 7/2` and `χ_c = 3`

Corollary 7 of `papers/three-colours/` (Corollary F12 of `notes/circular_planes.md`) says that below 4 the circular
chromatic number of the plane over a number field `F` is the circular chromatic number of a finite unit-distance graph
in `F²`; its proof, by compactness, gives no bound on the size. Here are explicit ones for `χ_c = 7/2`, over `ℚ(√11)`,
`ℚ(√191)`, `ℚ(√455)` and `ℚ(√911)`, and ones with nine vertices for `χ_c = 3` over `ℚ(√7)` and `ℚ(√31)`, whose planes
have no unit triangle (section on the value 3). For the four fields `χ_c(ℚ(√d)²) = 7/2`: the witnesses give `≥ 7/2`, and
Proposition 3 of the paper (a place with residue field `𝔽₇`) gives `≤ 7/2`, as 7 splits in `ℚ(√11)`, `ℚ(√191)` and
`ℚ(√911)` and ramifies in `ℚ(√455)`; this agrees with Theorem F, as (a) and (b) fail for `d ≡ 11 (mod 12)`, and does
not use its computer-assisted Proposition 9.

Theorem F16 of the note (Theorem G of the paper) adds the value 4: when `χ_c(F²) = 4`, some finite unit-distance
graph in `F²` has `χ_c = 4`, again with no bound on its size. The last three sections give explicit ones over
`ℚ(√3, √11)` and `ℚ(√2, √3)`, with 1 874 and 1 657 vertices, and their subgraphs in the Cayley graphs of the unit
vectors they were built from, which still have `χ_c = 4`.

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
| `witness_q11b` | `ℚ(√11)` | 30 | 155 | 404 | 1 273 (7, 14, 21, 28) | deletion from unions of `witness_q11` and other 7/2 witnesses (below) |
| `witness_q191` | `ℚ(√191)` | 240 | 293 | 803 | 489 (7 to 35) | growth from the 96-vertex graph to 3 258 vertices, then vertex deletion |
| `witness_q455` | `ℚ(√455)` | 780 | 175 | 434 | 91 (7, 14, 28) | growth from the 71-vertex graph to 959 vertices, then vertex deletion |
| `witness_q911` | `ℚ(√911)` | 1560 | 324 | 866 | 494 (7, 14, 21, 28, 35) | growth from the 327-vertex graph to 873 vertices, then deletion with solver cores (below) |
| `witness_q11sum` | `ℚ(√11)` | 30 | 2 237 | 11 300 | 180 (14 to 42) | the sumset `A + A` of the 76-vertex graph (the first one found) |

The five grown witnesses are also *vertex-critical*: for every vertex `v`, `witness_*_critical.json.gz` gives a
`(7, 2)`-colouring `c` of `H − v` whose tight digraph has no directed cycle, so `χ_c(H − v) < 7/2` (the other half
of Guichard's characterisation, explicitly: if `H − v` has `N` vertices and `pos` numbers them `0, …, N − 1` along a
topological order of the tight digraph, then `N·c + pos` is a homomorphism to `K_{7N/(2N+1)}`), and no proper induced
subgraph of `H` is a witness; `check_critical.py` checks them. As `√11 ∈ ℚ₇`, `witness_q11` is also an explicit
witness for `ℚ₇` (Theorem C) and for every finite extension `K` of `ℚ₇` with `χ_c(K²) < 4`, and so are `witness_q191`
and `witness_q911` (`√191`, `√911 ∈ ℚ₇`). `witness_q455` is one for every finite extension of `ℚ₇` with residue
degree 1 that contains `√455`, for example the ramified `ℚ₇(√455)` (Proposition 3 gives `χ_c ≤ 7/2` there).

**How they were found.** `grow.py` (colouring-guided growth with lazy SAT): start from the vertex set `A` of the
4-chromatic graph `data/quadratic_planes/q<d>.json`; ask a SAT solver for a `(7, 2)`-colouring in which every listed
directed cycle has a non-tight arc; if the colouring has tight cycles, list them (a shortest one through a vertex of
each nontrivial strong component of the tight digraph) and ask again; if its tight digraph is acyclic, add the points
`p = x + u` (`x` a vertex, `u` a unit vector with the same denominator) whose neighbours leave `p` no colour (or, if
there are none, those whose neighbours leave one colour), most neighbours first, at most 200 a round; over `ℚ(√455)`
the second rule supplied 382 of the 888 points added, over `ℚ(√11)` none. When the solver finds no such colouring, `H` is a witness. `minimise.py` then
deletes vertices one at a time, lowest degree first, while the formula stays unsatisfiable (one selector literal per
vertex in an incremental solver; a vertex outside the unsatisfiable core is dropped at once; the tight cycles it meets
are added to the list), and `critical.py` writes the criticality certificates. With python-sat 1.9 (CaDiCaL 1.5.3),
run in this folder,

    python3 grow.py ../../../quadratic_planes/q455.json A G 12000 200   # 71 rounds, 959 vertices: G.json
    python3 minimise.py G.json W.json                                    # 175 vertices
    python3 critical.py W.json C.json

reproduce the points, edges, colouring and cycles of `witness_q455` and its criticality certificates exactly; so do
`q11.json` with `3000 200` (52 rounds, 653 vertices) and then the same two commands for `witness_q11` (checked
again by the referee, below). For `witness_q191` the same commands with `q191.json` and `12000 200` (89 rounds, 3 258
vertices; `minimise.py` takes about 45 minutes) reproduce it exactly too, criticality certificates included. The files also carry `description` and `construction` (text) and
`fixed_vertex`, added afterwards (a vertex of largest degree, the first one); `critical.py` also writes
`not_critical` (empty). The three programs take any value `p/q` (`grow.py ... p q`;
the others read `p` and `q` from the file, `7/2` when absent), with the same results for `7/2`; the only change in
`minimise.py` for `7/2` is that its first and last solver calls now run without a conflict budget (the referee of the
value-3 witnesses compared the clause streams of both versions).

`witness_q11b` is smaller: the same growth from `q11.json` with at most 100 new points per round (`3000 100`) stopped
after 51 rounds at 573 vertices, and deletion left a second vertex-critical witness with 170 vertices (458 edges).
Deleting vertices from the union of the two (205 vertices) in random orders left vertex-critical witnesses with 157
and 161 vertices, and deleting vertices from the union of these two (170 vertices), again in a random order, left 155
vertices and 404 edges. These runs used the earlier working copies of the programs (deletion in a random order, with
one selector literal per vertex as in `minimise.py`) and are not reproduced here; the witness is certified like the
others (`verification.txt`), with `fixed_vertex` 108.

`witness_q911` comes from the same growth: `python3 grow.py ../../../quadratic_planes/q911.json A G 12000 200`
(63 rounds, 873 vertices, 54 cycles) reproduces that graph exactly. `kissat` refutes its formula in about a minute (with
the colour of a vertex of largest degree fixed), but the working copy of `minimise.py` needed half an hour for its first
refutation through python-sat and two hours to reach 747 vertices. The deletion therefore used `kissat` with DRAT
proofs: refute the current vertex set and keep the vertices whose clause "has a colour" lies in the clausal core of the
proof (`drat-trim -c`; from 873 vertices this alone left 730); then delete vertices, lowest degree first, one at a time
(to 707) and then in blocks of adaptive size up to 32 (a block is halved when the rest has a `(7, 2)`-colouring with an
acyclic tight digraph; the tight cycles of the colourings found are added to the list, as in `minimise.py`), keeping
the core after each refutation (to 489). From there the same deletion in blocks ran with CaDiCaL instead (python-sat,
one incremental solver for all calls, the failed-assumption core of each refutation in place of the `drat-trim` core,
at most 3 million conflicts per call), which was faster at that size (to 334). Near the end each refutation took
several minutes; splitting it into two cases made it much faster: with the colour of a vertex `r` of largest degree
fixed to 0, a neighbour `w` of `r` has colour 2, 3, 4 or 5, and as the cycle list is kept closed under reversal, the
reflection `c ↦ −c` exchanges the cases 2, 5 and 3, 4, so refuting the cases 2 and 3 suffices (from 334 vertices to
324 in 43 minutes, every vertex but `r` and `w` shown needed by a colouring of the rest). These working programs are
not stored; the witness is certified like the others, with the usual formula (`verification.txt`, `fixed_vertex`
87), and `critical.py` wrote its criticality certificates.

**Verification** (`verification.txt`): for every witness, two checks that share no code.
1. `check_witness.py` (written for `witness_q11sum`; it reads `d` from the file): exact
   integer checks of the points, of every edge and of all unit pairs, of the colouring and of the cycles; rebuilds
   the stored formula byte for byte; runs `drat-trim` on the stored DRAT proof: `s VERIFIED`.
2. `verify_independent.py` (written separately, from the file format only; it reads `d` and the denominator from the
   file): compares the points with `A`, recomputes all unit-distance pairs from the unit vectors with that
   denominator, checks the colouring and the cycles, and writes the formula with its own variable layout. `kissat`
   refutes it with a new DRAT proof, `drat-trim` verifies that proof and converts it to LRAT, and `cake_lpr` (the
   formally verified checker of CakeML) accepts the LRAT proof.

The stored proofs of the grown witnesses are the core proofs that `drat-trim -l` extracts from `kissat`'s proofs,
verified again by `check_witness.py`. To recheck one witness (`W` = `witness_q11`, `witness_q11b`, `witness_q191`,
`witness_q455`, `witness_q911`):

    xz -dk W.drat.xz
    python3 check_witness.py W.json.gz W.cnf.gz W.drat /path/to/drat-trim
    python3 check_critical.py W.json.gz W_critical.json.gz
    python3 verify_independent.py H.cnf W.json.gz   # then kissat H.cnf H.drat; drat-trim H.cnf H.drat -L H.lrat; cake_lpr H.cnf H.lrat

For `witness_q11sum` the commands are the same without `check_critical.py` (and `verify_independent.py H.cnf` alone
uses it by default). The tests are in `tests/test_two_primes.py` (`test_finite_witness_q11`,
`test_finite_witness_grown`, `test_finite_witness_critical_rejects`; the proof checks are marked slow).

**Referee** (`indep_W/`, with its report `REPORT.md`). A referee checked the three grown witnesses stored then
(`witness_q11`, `witness_q191`, `witness_q455`; `witness_q11b` and `witness_q911` came later: the two checkers above
check them, and we ran the referee's programs, unchanged, on them with the same results, `indep_W/results/`)
with programs of its own, written from the file format: all pairs of points (468, 803 and 434 unit pairs, none missing or extra), the
colourings and the cycles; a third encoding (two-sided arc indicators, variables renamed and negated at random); cycle
lists rebuilt from scratch by its own lazy SAT loop; a clause-by-clause validation of every formula, the stored ones
included; `kissat`, `drat-trim` and `cake_lpr` on all of them, and `cake_lpr` on the stored proofs; and the
criticality certificates through the explicit homomorphisms `N·c + pos` to `K_{7N/(2N+1)}`. It reproduced
`witness_q11` and `witness_q455` exactly with the programs here. It found no error; its remarks on the wording and on
the checkers are applied: `check_witness.py`, `verify_independent.py` and the other checkers now make every check
explicitly (an `assert` would be skipped under `python -O`), `check_witness.py` accepts only a line `s VERIFIED` from
`drat-trim` and writes the decompressed formula to a temporary file, which it removes.

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
| `witness_q31.json.gz` | the nine-point witness for `χ_c = 3` over `ℚ(√31)` (same fields as `witness_q7`) |
| `q31_seed.json` | the 159 points from which `grow.py` found it |
| `witness_q15.json.gz` | a 13-vertex witness for `χ_c = 3` over `ℚ(√15)` (same fields as `witness_q7`) |
| `q15_seed.json` | the 135 points from which `grow.py` found it |
| `witness_q7m.json.gz` | the subdivided hexagon `M` itself over `ℚ(√7)`: nine points, twelve edges (same fields) |
| `hexagon_search.py` | the search that found it: realisations of `M`, and of `M` plus one or two edges, over `ℚ(√d)` |
| `nine_vertices.py` | the three triangle-free graphs with 9 vertices and `χ_c = 3` (needs networkx) |
| `indep_W/` | a referee's programs and report for the three `7/2` witnesses (below) |
| `indep_V3/` | a referee's programs and report for the value-3 witnesses (below) |
| `verification.txt` | the outputs of both checks for every witness |
| `witness_q3_11.json.gz` | the witness for `χ_c = 4` over `ℚ(√3, √11)`: the denominator 84, the points (eight integers each; section on the value 4), the edges, a proper 4-colouring, the cycles and `fixed_vertex` (the origin) |
| `witness_q3_11.cnf.gz`, `witness_q3_11.drat.xz` | its formula (the one `check_witness4.py` builds) and a DRAT proof of unsatisfiability |
| `check_witness4.py` | the checker for the value 4 (see its docstring for the encoding and the logic) |
| `construction_q3_11/` | the programs that found it (section on the value 4) |
| `indep_W4/` | a referee's programs, results and report for it, for the lemma on base-point periods and for the finite construction (section on the value 4) |
| `witness_q2_3.json.gz` | a second witness for `χ_c = 4`, over `ℚ(√2, √3)`: 1 657 points, denominator 36, coordinates over `(1, √2, √3, √6)` (last section) |
| `witness_q2_3.cnf.gz`, `witness_q2_3.drat.xz` | its formula and a DRAT proof of unsatisfiability |
| `construction_q2_3/` | the programs that found it (last section) |
| `indep_W4b/` | a referee's programs, results and report for it (last section) |

## The value 3: nine points over `ℚ(√7)` and `ℚ(√31)`

`ℚ(√7)²` has no unit triangle, as `√3 ∉ ℚ(√7)`, and `χ_c(ℚ(√7)²) = 3` (Theorem B and Corollary 1 of the paper).
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
  increases by 1 mod 3), already one of the six listed directed cycles (the hexagon and two 6-cycles, each in both
  directions), so `χ_c(H₇) ≥ 3` by Lemma 20. The paper proves this by hand, for the subgraph `M = H₇ − P₁P₅`: the hexagon `v₀⋯v₅ = P₀P₄P₃P₂P₆P₇` with its three long
  diagonals subdivided by `m₀, m₁, m₂ = P₁, P₅, S`. With `δ = ±1` the colour step along an arc (mod 3), the sum of `δ`
  is ±3 on a 5-cycle and 0 or ±6 on a 6-cycle (±6: tight); the 5-cycles `Z_k = v_k v_{k+1} v_{k+2} v_{k+3} m_k` satisfy
  `C_k = Z_{k+1} − Z_k = v_k m_k v_{k+3} v_{k+4} m_{k+1} v_{k+1}` (a 6-cycle; `C_{k+3}` is `C_k` reversed),
  `Z₃ − Z₀ = C₀ + C₁ + C₂` and `Z₀ + Z₃ =` the hexagon, so if none of `C₀, C₁, C₂` is tight, `Z₀` and `Z₃` have the
  same sum and the hexagon has sum ±6. In `H₇`, `C₀ = P₀P₁P₂P₆P₅P₄` has the chord `P₁P₅` and is never tight; the
  listed cycles are `C₁`, `C₂` and the hexagon. Independently,
  `H₇` has no homomorphism to `K_{8/3}`, which is the Wagner graph itself; as the circular chromatic number of a graph
  with 9 vertices is a fraction with numerator at most 9 (Vince 1988; Bondy and Hell 1990; Zhu's survey), and `8/3` is
  the largest such fraction below 3,
  this gives `χ_c(H₇) ≥ 3` again.
- *Vertex-critical.* For every vertex `v` a 3-colouring of `H₇ − v` with an acyclic tight digraph is stored.
- *Smallest possible.* Every triangle-free graph with at most 8 vertices maps to `K_{8/3}`, so a graph with `χ_c = 3`
  in a plane without unit triangles has at least 9 vertices. Of the 1 897 triangle-free graphs with 9 vertices (up to
  isomorphism) exactly three have `χ_c = 3`: `M` (12 edges), `H₇ = M + m₀m₁` (13 edges) and `M + m₀m₁ + m₁m₂`
  (14 edges); `nine_vertices.py` finds them, with two separate tests of `χ_c = 3` (no homomorphism to `K_{8/3}`; every
  3-colouring has a tight cycle) that agree on every graph.
- As `√7 ∈ ℚ₃`, `H₇` is also a witness for `ℚ₃` and for every finite extension `K` of `ℚ₃` with `χ_c(K²) < 4`.

`witness_q7m.json.gz` is `M` itself, as an induced unit-distance graph on nine points of `ℚ(√7)²` (denominator 80):
the hexagon `v₀, …, v₅ = (√7/4, 1/4)`, `(√7/4 − 1, 1/4)`, `(−3/4, −√7/4)`, `(0, 0)`, `(−1, 0)`, `(−1/4, √7/4)` (vertices
0–5) and the midpoints `m₀, m₁, m₂ = (0, 1)`, `(√7/4 − 1, −3/4)`, `(−1/2 − √7/4, 1/4)` of the diagonals `v₀v₃`, `v₁v₄`,
`v₂v₅` (vertices 6–8), with exactly the twelve unit pairs of `M`. Its eight listed cycles are the hexagon in both
directions and the six directed 6-cycles `C₀, …, C₅` of the proof. It is a witness with the fewest vertices and,
among those, the fewest edges.

It was found by `hexagon_search.py`, a direct search for equilateral hexagons with unit vectors of one denominator
whose long diagonals are sums of two unit vectors (exact integer arithmetic). Without unit triangles, nine points
with the twelve unit pairs of `M` span `M` or `M` plus one or two edges between midpoints, since any other pair at
distance 1 would close a triangle (a pair `v_i v_{i+2}` with `v_{i+1}`, a pair `v_i v_{i+3}` with `m_i`, a pair
`v_i m_j` with `v_j` or `v_{j+3}`); the search counts the point sets of each kind. `python3 hexagon_search.py 7 80 400`
(about a minute) lists 400 realisations of `M` (the stored one is the 181st, translated) and counts 62 256 sets of
nine points of `ℚ(√7)²` with denominator 80 that realise the twelve edges of `M`: 46 512 span `M`, 14 688 span `M`
plus one edge and 1 056 span `M` plus two. With `31 80` the counts are 4 288, 1 728, 1 920 and 640. So all three
nine-vertex graphs with `χ_c = 3` are unit-distance graphs over `ℚ(√7)` and over `ℚ(√31)`; with denominator 80
there are none over `ℚ(√15)`, `ℚ(√23)` or `ℚ(√39)`.

`witness_q31.json.gz` is the third of these graphs, `M + m₀m₁ + m₁m₂`, as an induced unit-distance graph on nine points
of `ℚ(√31)²` (denominator 80; hexagon `v₀, …, v₅` = vertices 0, 1, 8, 7, 2, 6 and midpoints `m₀, m₁, m₂` = 3, 4, 5):
`(−8/5, −4/5)`, `(−1, 0)`, `((−60 − 3√31)/80, (45 − 4√31)/80)`, `(−3/5, −4/5)`, `(0, 0)`, `(√31/16, −15/16)`,
`((−30 + √31)/40, (−15 − 2√31)/40)`, `((−108 − 3√31)/80, (−19 − 4√31)/80)`, `(−1 + √31/16, −15/16)`. It is
vertex-critical too, and it was found by the same growth over `ℚ(√31)` and deletion, from `q31_seed.json` (0, the 156
unit vectors with denominator 80, and the two points `(−8/5, −4/5)` and `(−8/5 − √31/16, 11/80)` that close the
5-cycle `0, (−1, 0), (−8/5, −4/5), (−8/5 − √31/16, 11/80), (−17/20 − √31/40, −17/40 + √31/20)`):

    python3 grow.py q31_seed.json A G31 2000 200 3 1   # 2 growth rounds, 367 vertices, 7 cycles
    python3 minimise.py G31.json W31.json              # 9 vertices: the stored points, edges, colouring and cycles
    python3 critical.py W31.json C31.json              # the stored certificates

The growth also finds witnesses over `ℚ(√15)` and `ℚ(√39)`, with sizes that depend on the seed. Over `ℚ(√15)`, from
70 seeds (0, the unit vectors with denominator 80 or 160, and the points of one unit 5-cycle through 0; 63 of the
growths ended with a witness) and 174 deletions (the lowest-degree order for the 63, and three random orders for each
of the 37 from denominator 80), the smallest has 13 vertices and 18 edges, and the 30 of that size are all the same
graph up to isomorphism: a vertex `a` joined by paths of length 2 to five vertices `b, c, d, e, f`,
which are joined by the edges `bc`, `bd`, `de`, `ef` and the paths `b–f` and `c–e` of length 2 (every cycle of a
minimum cycle basis has length 5). It is stored as `witness_q15.json.gz` and found from `q15_seed.json` (0, the 132 unit
vectors with denominator 80, and the two points `((21 − √15)/40, (−7 − 3√15)/40)` and `((21 − 11√15)/40, (3 − 3√15)/40)`,
which close the 5-cycle `0, (−√15/8, −7/8), ((21 − √15)/40, (−7 − 3√15)/40), ((21 − 11√15)/40, (3 − 3√15)/40),
((3 − 4√15)/20, (4 + 3√15)/20)`):

    python3 grow.py q15_seed.json A G15 3000 200 3 1   # 29 rounds, 540 vertices, 23 cycles
    python3 minimise.py G15.json W15.json              # 13 vertices, 18 edges, 12 cycles: the stored ones
    python3 critical.py W15.json C15.json              # the stored certificates

None of the three nine-vertex graphs has a realisation over `ℚ(√15)` with coordinates of denominator 48, 68, 80, 104,
112, 120, 136, 160, 208, 221, 240, 272 or 320 (`hexagon_search.py`), although `ℚ(√15)²` has unit 5-cycles (for example
`2a + 2b + c = 0` with `a·b = −7/8`, which needs `√15`). So over `ℚ(√15)` the least size is between 9 and 13; over
`ℚ(√39)` the smallest found (30 seeds) has 17 vertices.

The checks are finite enumerations (`verification.txt`):

    python3 check_small.py witness_q7.json.gz q7.cnf    # points, unit pairs, all 3^9 maps, K_{8/3}, criticality
    python3 check_small.py witness_q31.json.gz q31.cnf  # the same for the witness over Q(sqrt31)
    python3 check_small.py witness_q7m.json.gz q7m.cnf  # and for M itself over Q(sqrt7)
    python3 check_small.py witness_q15.json.gz q15.cnf  # the 13-vertex witness over Q(sqrt15) (all 3^13 maps)
    python3 hexagon_search.py 7 80 400                  # the realisations of M and of M plus one or two edges
    python3 small_triangle_free.py                      # the 4 682 270 triangle-free graphs on 8 labelled vertices
    python3 nine_vertices.py                            # the triangle-free graphs on 9 vertices (needs networkx)

`check_small.py` also writes the formula "a 3-colouring in which every listed cycle has a non-tight arc" (the layout of
`check_witness.py`, with `(p, q) = (3, 1)`); kissat refutes it, and drat-trim and `cake_lpr` verify the proof.

**How it was found.** From `q7_seed.json` (0, the 204 unit vectors with denominator 160, and `(−2, 0)` and
`((−11 − √7)/8, (−5 − √7)/8)`, which close the 5-cycle `0, (−1, 0), (−2, 0), ((−11 − √7)/8, (−5 − √7)/8),
(−3/4, −√7/4)`),

    python3 grow.py q7_seed.json A G7 2000 200 3 1   # 2 growth rounds, 607 vertices, 7 cycles
    python3 minimise.py G7.json W7.json              # 9 vertices
    python3 critical.py W7.json C7.json

give the nine points translated by `(−1, 0)` and in another order, with the same edges, colouring, cycles and
certificates under that correspondence; the stored file uses the labels of the paper (and keeps the denominator 160
of the growth, although every coordinate has a denominator dividing 8).

**Referee** (`indep_V3/`, with its report `REPORT.md`). A referee checked the value-3 material with programs of its
own: exact arithmetic on the points; all `2^28` labelled graphs on 8 vertices; a Burnside count of the graphs on 9
vertices with no isomorphism test (1 897 classes, three with `χ_c = 3`); an exhaustive search for homomorphisms to
`K_{8/3}`; the proof by hand on all 126 colourings of `M`; the explicit homomorphisms `N·c + pos` of every `H − v`;
the formulas clause by clause; and the reproduction over `ℚ(√7)`. It found no error. Its findings were applied: the
seed of the `ℚ(√31)` growth is stored, the remark on `ℚ(√15)` and `ℚ(√39)` is corrected, all three graphs are shown to
be unit-distance graphs, Vince and Bondy–Hell are cited, the proof by hand names the cycles it uses, and Question 3
is worded more precisely.

## The value 4: a witness over `ℚ(√3, √11)`

`witness_q3_11.json.gz` is an induced unit-distance graph `H₄` in `ℚ(√3, √11)²` with 1 874 vertices and 8 085 edges, and
`χ_c(H₄) = χ(H₄) = 4`. Its points lie in `ℤU` for the 27 unit vectors `U` below, but 198 of its edges are unit vectors
other than those of `±U` (12 directions up to sign): the clausal cores of the construction were taken with all unit
pairs as edges. A point is `[a0, a1, a2, a3, b0, b1, b2, b3]`, meaning
`((a0 + a1√3 + a2√11 + a3√33)/84, (b0 + b1√3 + b2√11 + b3√33)/84)`; distance 1 is four integer equations in the basis
`(1, √3, √11, √33)`. As `χ(ℚ(√3, √11)²) = 4`, `χ_c(ℚ(√3, √11)²) = 4`, and `H₄` attains it; Theorem F16 says that some
finite graph does, without a bound.

- *Upper bound.* A proper 4-colouring (stored).
- *Lower bound.* Every proper 4-colouring of `H₄` has a tight cycle (colour difference `1 mod 4` along every arc)
  among the 4 992 listed ones (4 916 of length 4, 52 of length 8, 24 of length 12): the formula of `check_witness4.py`
  ("a proper 4-colouring, the origin coloured 0, in which no listed cycle is tight"; no auxiliary variables) is
  unsatisfiable, and Lemma 20 of the paper (Guichard) gives `χ_c(H₄) ≥ 4`.
- *Verification.* `xz -dk witness_q3_11.drat.xz` and
  `python3 check_witness4.py witness_q3_11.json.gz witness_q3_11.cnf.gz witness_q3_11.drat drat-trim` check the points,
  all 1 755 001 pairs, the colouring, the cycles and the formula (sha256 `febd6b3d…`), and drat-trim prints
  `s VERIFIED` (about 15 seconds). The referee's programs (`indep_W4/`, report `REPORT.md`, outputs in `results/`)
  share no code with it: they recompute every unit pair exactly and in floating point, write their own encoding
  (variables `x(v, k) = kn + v + 1` and one per tight arc), and kissat, `drat-trim -L` and `cake_lpr` print
  `s UNSATISFIABLE`, `s VERIFIED` and `s VERIFIED UNSAT`; their tests reject a moved point, an improper colouring and a
  cycle through a non-edge, and find the formula satisfiable without the cycle clauses.
- `H₄` has not been minimised (74 of its vertices lie on no listed cycle).

**How it was found** (`construction_q3_11/`; Lemma F17 of the note, base-point periods). Let `U` be the 27 unit
vectors `R₆₀ʲ R_Aᵏ R_Gˡ (1, 0)` of `../at_four/q3_11.py`, and `c` a proper 4-colouring without tight cycles of a finite
graph whose edges are differences in `±U`. For a closed walk `W` let `Λ(W)` be the sum of the representatives in
`{1, 2, 3}` of the colour differences along it, and `per(W) = Λ(W)/4 − N(W)`, with `N(W)` the number of steps in
`−U` (no vector of `U` is the negative of another). Then `per` is an integer, additive, unchanged by backtracks and by
exchanging two consecutive steps whose square lies in the graph (otherwise the square is a tight cycle), and
`(P − 3N)/4 < per(W) < (3P − N)/4` when `W` has `P` steps in `U` and `N` in `−U`, `P + N ≥ 1` (otherwise `W` or its
reverse is tight). Walks of relations `ρ = Σ n_u u = 0` joined by such moves to the walks of a ℤ-basis of the
relation lattice make `per` an integer homomorphism `p` on it with `p(ρ)` strictly inside these ranges for `ρ ≠ 0`; a
character of `ℤU` with representatives in `(1/4, 3/4)` would give one, so if none exists, `κ(U) ≤ 1/4`. With python-sat
1.9 (CaDiCaL 1.5.3), OR-Tools 9.15, PARI/GP 2.15, kissat 4.0.4 and drat-trim, run in a scratch folder (`KISSAT` and
`DRAT_TRIM` name the binaries for `iter_shrink.sh`; `C` is this folder),

    python3 $C/short_rel.py 12 rel_B12.json                  # the 4 641 relations with Σ n_u² ≤ 12
    python3 $C/per_cegar.py rel_B12.json 400 added_B12.json  # 41 more, of length 8; then no p is left
    python3 $C/per_ilp.py rel_B12.json added_B12.json        # core.json: 256 of them leave none
    python3 $C/short_rel.py 14 rel_B14.json                  # the pool of the chains
    python3 $C/per_build.py core.json rel_B14.json added_B12.json H_core1.json   # 13 230 points
    python3 $C/cycles_of.py H_core1.json Wc_core1.json       # 37 290 directed cycles
    cp $C/certify4.py $C/shrink_core.py $C/iter_shrink.sh .
    python3 certify4.py Wc_core1.json core1.cnf              # kissat refutes it in about a minute
    kissat --seed=1 core1.cnf core1.drat; drat-trim core1.cnf core1.drat -c core1.core
    python3 shrink_core.py Wc_core1.json core1.core W_s1.json; python3 certify4.py W_s1.json s1.cnf   # 2 688 points
    ./iter_shrink.sh s1 6; PFX=j ./iter_shrink.sh it6 5      # eleven more clausal cores: 1 874 points
    python3 $C/make_final4.py W_j5.json witness_q3_11

reproduce `witness_q3_11.json.gz` and its formula exactly (checked again from these copies). Each round of
`iter_shrink.sh` refutes the current formula with kissat and keeps the vertices that occur in the clausal core
(`drat-trim -c`). Each of the 41 added relations lies among the vectors with one value of `k`, so it is `R_A^k` times a
relation among the 18 vectors `ζ^j g^l` (`ζ = e^{iπ/3}`, `g = (11 + 5√−3)/14`, the rotation `R_G`); each has length 8
and a coefficient ±3 (absolute coefficients 3, 3, 1, 1 for 17 of them, such as `3ζ̄ + ζ + g − 3ζ̄g = 0`; 3, 2, 1, 1, 1
for 14; 3, 2, 2, 1 for 10), so its squared norm (20, 16 or 18) exceeds the bound 12 of the enumeration; with all 18 138
relations of squared norm at most 14 but without them a solution `p` exists. `core_min.py` reduces the infeasible system
to 37 relations of length 8, whose walks and chains give 4 134 points, but kissat did not refute that formula in ten
minutes. The integer solver (CP-SAT) only guided the construction; the proof of `χ_c(H₄) = 4` is the SAT refutation
above.

**The construction always works** (the remark at the end of §6.9 of the note, after Lemma 24 in the paper). When
`κ(U) ≤ 1/4` (over all characters of `ℤU`; `U ∩ −U = ∅`, so no element of `U ∪ −U` has order 2), some finite set of
relations leaves no `p`: the ranges of the basis relations allow finitely many `p`, and for each an integer point of a
rational cone gives a relation `ρ_p` outside whose range `p` falls. The subgraph of `Cay(ℤU, U ∪ −U)` induced by `0`,
the walks along the chains to the basis and to the `ρ_p`, and the corners of the swaps then has a tight cycle in every
proper 4-colouring, so this replaces the compactness step of Corollary F14 by a finite construction. The referee checked
the remark (`indep_W4/REPORT.md`, section on the remark) and tested it with `indep_W4/remark_check.py` on 15 sets of
integer distances: for each, the system on the relations found has no solution, and a SAT encoding that does not use the
lemma shows that every proper 4-colouring of the graph has a tight cycle; in 10 of them the chains of the basis alone
are not enough. `results/remark_check_all.log` is a rerun from this copy (`KISSAT` names the binary), and
`results/remark_check_referee.log` the referee's own run.

## A second witness at 4: over `ℚ(√2, √3)`

`witness_q2_3.json.gz` is an induced unit-distance graph `H₄′` in `ℚ(√2, √3)²` with 1 657 vertices and 6 238 edges, and
`χ_c(H₄′) = χ(H₄′) = 4`. A point is `[a0, a1, a2, a3, b0, b1, b2, b3]`, meaning
`((a0 + a1√2 + a2√3 + a3√6)/36, (b0 + b1√2 + b2√3 + b3√6)/36)`; `check_witness4.py` reads the field from the file. As
`χ(ℚ(√2, √3)²) ≤ 4` (a residue colouring at the prime above 2; `papers/planes-4-chromatic/`), `χ_c(ℚ(√2, √3)²) = 4`, and
`H₄′` attains it; it also reproves `χ(ℚ(√2, √3)²) ≥ 4`. The value also follows from Theorem F (no prime of `ℚ(√2, √3)`
above 7 has residue degree 1, as 3 is not a square modulo 7), but `H₄′` does not need it. 39 of its edges are unit
vectors other than the 120 below (17 directions up to sign).

- *Upper bound.* A proper 4-colouring (stored).
- *Lower bound.* Every proper 4-colouring of `H₄′` has a tight cycle among the 6 062 listed ones (5 494 of length 4,
  568 of length 8): the formula of `check_witness4.py` (sha256 `cd6382b7…`) is unsatisfiable. kissat 4.0.4
  (`--seed=1`) refutes it in about 40 seconds, and drat-trim verifies the proof `witness_q2_3.drat.xz` (20.6 MB
  uncompressed).
- *Verification.* `xz -dk witness_q2_3.drat.xz` and
  `python3 check_witness4.py witness_q2_3.json.gz witness_q2_3.cnf.gz witness_q2_3.drat drat-trim` (about 25
  seconds). The referee's programs (`indep_W4b/`, report `REPORT.md`, outputs in `results/`)
  were written from the file format alone and share no code with it: they test all 1 371 996 pairs exactly (and in
  floating point), write their own encoding (one variable per tight arc; 15 452 variables, 77 910 clauses) and a
  second one without auxiliary variables, and kissat, `drat-trim -L` and `cake_lpr` print `s UNSATISFIABLE`,
  `s VERIFIED` and `s VERIFIED UNSAT` for both; their 77 sanity checks (moved points, an improper colouring, cycles
  through non-edges, a missing or an extra edge, …) all pass.

**How it was found** (`construction_q2_3/`; the method of the previous section). The plane `ℚ(√2, √3)²` is the field
`ℚ(ζ₂₄)`, `ζ₂₄ = e^{iπ/12}` (the point `(x, y)` is `x + iy`), and its unit vectors are the elements of norm 1. Take the
unit vectors `ζ₂₄^j w^l` (`j` modulo 24, `|l| ≤ 2`) with `w = (1 + 2√−2)/3`, the rotation with cosine 1/3, one of each
pair `±u`: 60 vectors, with coordinates over the power basis `1, ζ₂₄, …, ζ₂₄⁷` and denominator 9, and a relation lattice
of rank 52 (`setup_z24.py`, with PARI/GP; `κ` of these vectors is about 1/6 by a floating-point integer program,
`kappa_z24.py`, which is not needed). Already the 1 552 relations `Σ n_u u = 0` with `Σ n_u² ≤ 8` leave the period
system of Lemma F17 without a solution, with no counterexample search; a minimal infeasible subset has 279 relations,
all of length 8. Their walks, tied to the basis relations by bubble sorts, span 2 220 points, and twelve rounds of
clausal cores reduce them to 1 657. With python-sat 1.9 (CaDiCaL 1.5.3), OR-Tools 9.15, PARI/GP 2.15, kissat 4.0.4 and
drat-trim, run in a scratch folder (`C` is this folder; `KISSAT` and `DRAT_TRIM` name the binaries),

    python3 $C/setup_z24.py s2:2 cfg_s2_2.json             # the 60 vectors and a reduced basis of the 52 relations
    export CFG=cfg_s2_2.json
    python3 $C/short_rel.py 8 rel8.json                      # the 1 552 relations with Σ n_u² ≤ 8
    python3 $C/per_ilp.py rel8.json                          # no p is left; core.json: 379 of them
    python3 $C/core_min.py core.json core_min.json 1         # 279
    python3 $C/short_rel.py 10 rel10.json                    # the pool of the chains
    python3 $C/per_build.py core_min.json rel10.json H_m.json   # 2 220 points
    python3 $C/cycles_of.py H_m.json Wc_m.json               # 7 854 directed cycles
    python3 $C/certify_z24.py Wc_m.json m.cnf; cp Wc_m.json W_m.json
    $C/iter_shrink_z.sh m 12                                 # twelve clausal cores: 1 657 points
    python3 $C/make_final_z.py W_z12.json witness_q2_3

reproduce `witness_q2_3.json.gz` byte for byte, and its formula (checked again from these copies). `make_final_z.py`
writes the points in plane coordinates over `(1, √2, √3, √6)` with denominator 36 (`to_xy.py`), recomputes every unit
pair exactly, and writes the formula of `check_witness4.py`; it is byte for byte the formula that `certify_z24.py`
built in the arithmetic of `ℚ(ζ₂₄)`, so the two computations of the unit pairs agree. `per_ilp.py`, `core_min.py`,
`cycles_of.py` and `shrink_core.py` are copies of those in `construction_q3_11/`; `short_rel.py` and `per_build.py`
read the vectors and the basis from `$CFG`.

## Inside the Cayley graphs: `witness_q3_11_cayley` and `witness_q2_3_cayley`

`H₄` and `H₄′` take all unit pairs among their points as edges, and some of these are not differences in `±U` for the
unit vectors `U` of their constructions: 198 of the 8 085 edges of `H₄` (12 directions up to sign) and 39 of the 6 238
edges of `H₄′` (17 directions). Their subgraphs `H₄°` and `H₄′°` on the same points, with only the pairs that differ by
an element of `U ∪ −U` as edges, still have `χ_c = 4`, with more listed cycles:

| witness | field | vertices | edges | generators `U` | listed cycles (lengths) | sha256 of the formula |
|---|---|---|---|---|---|---|
| `witness_q3_11_cayley` | `ℚ(√3, √11)` | 1 874 | 7 887 | the 27 vectors of `../at_four/q3_11.py` | 3 389 (4, 8) | `b394ee93…` |
| `witness_q2_3_cayley` | `ℚ(√2, √3)` | 1 657 | 6 199 | the 60 vectors `ζ₂₄^j w^l` | 5 264 (4, 8) | `802921ac…` |

So `H₄°` is a finite subgraph of the Cayley graph of the 54 vectors `±U` with `χ_c = 4`: the graph that Corollary 8 of
the paper (Corollary F14 of the note) gives for these vectors, here explicit. Every point is joined to the origin by
steps in `±U` inside the graph, so all points lie in `ℤU`. The files have one more key, `generators` (one vector of each
pair `±u`, in the coordinates of the points), and `check_witness4.py` then also checks (1b): the generators are unit
vectors, pairwise distinct up to sign; the edges are exactly the pairs of points that differ by an element of `U ∪ −U`;
and every point is joined to the fixed vertex (the origin) by edges. The formula is the same as for the other witnesses.

- *Verification.* `xz -dk witness_q3_11_cayley.drat.xz` and
  `python3 check_witness4.py witness_q3_11_cayley.json.gz witness_q3_11_cayley.cnf.gz witness_q3_11_cayley.drat drat-trim`
  (about 10 seconds; the same for `witness_q2_3_cayley`, about 20 seconds). kissat 4.0.4 refutes the two formulas in
  about 15 and 25 seconds, and `drat-trim -L` turns the stored proofs into LRAT proofs that `cake_lpr` accepts
  (`s VERIFIED UNSAT`). `tests/test_at_four.py` checks the generators against `../at_four/q3_11.py` and against the 60
  vectors rebuilt in its own arithmetic, and that an extra edge in another direction is rejected. A referee's programs
  (`indep_W4c/`, report `REPORT.md`, outputs in `results/`) share no code with these: they rebuild the 27 and the 60
  vectors in their own arithmetic, check the Cayley pairs, the colouring and the cycles exactly, write their own
  encoding, and kissat, `drat-trim -L` and `cake_lpr` print `s UNSATISFIABLE`, `s VERIFIED` and `s VERIFIED UNSAT` for
  both files; their 29 mutations of each file are all rejected. They also show that the stored colourings are proper on
  the full unit-distance graphs on the same points (with the 198 and 39 other unit pairs).
- *How they were found* (`construction_cayley/`). On the points of `H₄` (`H₄′`) with the Cayley edges only, the formula
  with the listed cycles is satisfiable. `cayley.py init` repeats: kissat finds a proper 4-colouring without a tight
  listed cycle, and every tight directed cycle that a depth-first search finds among the tight arcs of that colouring is
  added to the list; after three rounds (751 cycles added; for `H₄′` two rounds, 422) the formula is unsatisfiable.
  Three rounds of clausal cores (`shrink_cayley.sh`) then keep 3 389 (5 264) cycles, all of length 4 or 8. In a scratch
  folder, with `KISSAT` and `DRAT_TRIM` naming the binaries and `C` this folder's `construction_cayley/`,

    python3 $C/cayley.py init q3_11 cay_q3_11.json                 # 349 + 133 + 269 tight cycles, then unsatisfiable
    $C/shrink_cayley.sh cay_q3_11.json 3 witness_q3_11_cayley     # 5 743 -> 4 013 -> 3 617 -> 3 389 cycles
    python3 $C/cayley.py init q2_3 cay_q2_3.json                   # 165 + 257 tight cycles
    $C/shrink_cayley.sh cay_q2_3.json 3 witness_q2_3_cayley       # 6 484 -> 5 678 -> 5 420 -> 5 264 cycles

reproduce both files, their formulas and their DRAT proofs byte for byte (kissat 4.0.4; checked from these copies).
