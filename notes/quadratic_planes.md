# Planes over real quadratic fields that need four colours

For a field `K ⊂ ℝ`, `χ(K²)` is the least number of colours for the points of
the plane with coordinates in `K` such that two points at distance 1 get
different colours. For real quadratic fields `K = ℚ(√d)` the known values were
2 and 3. No real quadratic field was known to need four colours.

**Theorem.**
1. `χ(ℚ(√d)²) = 4` for `d = 11, 23, 35, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959`.
2. `4 ≤ χ(ℚ(√47)²) ≤ 5`.

Each lower bound is a finite graph: a triangle-free unit-distance graph with
coordinates in `ℚ(√d)`, which is not 3-colourable and is vertex-critical. The
graphs have 76 vertices (`d = 11`), 393 (`d = 23`), 580 (`d = 35`), 816
(`d = 47`), 406 (`d = 59`), 611 (`d = 71`), 1 404 (`d = 95`), 399 (`d = 119`),
356 (`d = 131`), 1 281 (`d = 155`), 259 (`d = 179`), 96 (`d = 191`), 338
(`d = 239`), 291 (`d = 251`), 394 (`d = 263`), 715 (`d = 359`), 331 (`d = 431`),
703 (`d = 443`), 71 (`d = 455`), 835 (`d = 491`), 659 (`d = 599`), 712
(`d = 611`), 898 (`d = 791`), 538 (`d = 851`), 327 (`d = 911`), 252 (`d = 935`)
and 513 (`d = 959`).

So `χ(ℚ(√d)²)` is now known for every squarefree `d < 83` except `d = 47`
(§1): it is 2 if `d ≡ 1, 2 (mod 4)`, 3 if `d ≡ 3 (mod 4)` and `d ≢ 2 (mod 3)`,
and 4 for `d = 11, 23, 35, 59, 71`.

That they are not 3-colourable is a computer proof: SAT solving with a DRAT
proof checked by drat-trim, done for two separate encodings of each graph (§4).
The rest of each certificate is checked in exact integer arithmetic. Nobody
outside the project has refereed it.

## 1. What was known

Let `d ≥ 2` be squarefree and `K = ℚ(√d)`.

| `d` | `χ(K²)` | source |
|---|---|---|
| `d ≡ 1, 2 (mod 4)` | `= 2` | Johnson (1987); Moorhouse (2010), Lemma 8.4 and Theorem 8.5 |
| `d ≡ 0, 1 (mod 3)` | `≤ 3` | Fischer (1990), Theorem 9; Moorhouse, Corollary 8.3 |
| `d ≡ 0, 1, 2, 4 (mod 7)` | `≤ 4` | Moorhouse, Corollary 8.3 |
| `d ≡ 3 (mod 8)` | `≤ 4` | Fischer (1990), Theorem 10; hn-2adic-obstruction (2026), `RESULTS.md` §0, item 3 (Corollary D of its `NOTES.md`); `notes/local_colourings.md`, Proposition A |
| `d ≡ 3 (mod 4)` | `≥ 3` | Fischer (1990), Theorem 8; Moorhouse, Theorem 8.6, for prime `d` (an odd cycle from Pell's equation); directly, the unit vector `z = (1 + i√d)/(1 − i√d)` gives `((d + 1)/4)(z + z̄) + ((d − 1)/2)·1 = 0`, a closed walk of odd length `d` |

Madore (2015) proved `χ(ℚ(√3)²) = χ(ℚ(√7)²) = 3`. Moorhouse's Theorem 8.1,
`χ(K²) ≤ 4` unless `d ≡ 47, 59, 83 (mod 84)`, combines the rows on 3 and 7.
The public repository
[hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction)
(July 2026, not refereed; `notes/literature.md`, item 1) proves the rows on
`d ≡ 1, 2 (mod 4)` and `d ≡ 3 (mod 8)` again by reduction at 2 (its
`RESULTS.md`, §0, item 3, and Corollary D of its `NOTES.md`), and finds odd
cycles over `ℚ(√11)`.

So a real quadratic field can need four colours only if `d ≡ 11 (mod 12)`:
`d ≡ 3 (mod 4)` (else 2 colours suffice) and `d ≡ 2 (mod 3)` (else 3 suffice).
The first such values are 11, 23, 35, 47, 59, 71, 83, …. For these fields the
known lower bound was 3. For `ℚ(√47)` the bounds were `3 ≤ χ ≤ 5`, and our own
experiments of 25 September had left it open (`docs/research-log.md`). Moorhouse
asked both questions in a talk in 2010 ("Colouring the plane", Designs, Codes &
Geometries 2010; slides at
[ericmoorhouse.org/slides/ebertfest.pdf](https://www.ericmoorhouse.org/slides/ebertfest.pdf)):
"What can we say about `χ(K²)` when `K = ℚ(√d)`?" and "What about `χ(ℚ(√47)²)`?"

**The literature search.** On 1 October 2026 we read, for statements about
`χ(ℚ(√d)²)`: Moorhouse's draft and his talk (2010); Fischer (1990, 1994); Madore
(2015); Payne (2009); Axenovich, Choi, Lastrina, McKay, Smith and Stanton
(Graphs Combin. 30 (2014), 71–81), whose Theorem 2.2 collects the upper bounds;
Currie and Eggleton (arXiv:1509.03667); Bardestani and Mallahi-Karai
(arXiv:1507.05300); Cohen (2007); hn-2adic-obstruction; the Polymath16 threads
and wiki; and the sources listed in `notes/literature.md`. We also searched
arXiv titles and abstracts, MathOverflow, Mathematics Stack Exchange and the
web. None of these gives a real quadratic field that needs four colours: the
largest lower bound they give for a real quadratic field is 3, and Madore (§2)
writes that practically the only useful graphs here are the triangle and the
Moser spindle. We could not read Johnson (1987), Fischer's paper on connected
components (Congr. Numer. 72, 1990), Chilakamarri's survey (1993), Benda and
Perles (2000) or Soifer's book (2024), and we could not query the citation lists
of zbMATH, Semantic Scholar or OpenAlex. An earlier proof may exist that we did
not find.

**A different question: `ℚ(√−d)` inside `ℂ`.** Cohen (2007) wrote that "the
natural conjecture is that `ℚ[α]` is 3-colorable for all `α` quadratic over `ℚ`",
for the set `ℚ[α] ⊂ ℂ`. For `α = √−d` this set is `{x + y√−d}`, of dimension 2
over `ℚ`; the plane `ℚ(√d)²` has dimension 4. So the theorem above does not
contradict the conjecture, and the conjecture holds, by the usual reduction at
one prime:
- for real `α`, the set lies on a line;
- for `d = 1, 2` it lies in a plane of the first row of the table;
- otherwise `d` has an odd prime factor `p`, whose prime `𝔭` in `ℚ(√−d)` is
  ramified. Every unit `z` of `ℚ(√−d)` is a `𝔭`-adic unit with `z ≡ z̄ (mod 𝔭)`,
  so `z² ≡ z z̄ = 1` and `z ≡ ±1 (mod 𝔭)`. Colour each point by its residue
  modulo `𝔭`, relative to a representative of its coset of the `𝔭`-integral
  elements, followed by a proper 3-colouring of the cycle `Cay(𝔽_p, {±1})`.

## 2. The graphs

A point is written `[a, b, c, e]` and stands for

    ((a + b√d)/D, (c + e√d)/D),

with integers `a, b, c, e` and a fixed denominator `D`. Two points are at
distance 1 exactly when the differences satisfy

    Δa² + d Δb² + Δc² + d Δe² = D²   and   Δa Δb + Δc Δe = 0.

These are two identities between integers, so every edge is checked exactly.

| `d` | `D` | vertices | edges | degrees | directions used |
|---|---|---|---|---|---|
| 11 | 30 | 76 | 172 | 3–26 | 46 |
| 23 | 156 | 393 | 1 011 | 3–47 | 76 |
| 35 | 390 | 580 | 1 501 | 3–107 | 218 |
| 47 | 240 | 816 | 2 134 | 3–72 | 94 |
| 59 | 210 | 406 | 993 | 3–48 | 66 |
| 71 | 120 | 611 | 1 557 | 3–85 | 128 |
| 95 | 480 | 1 404 | 3 780 | 3–99 | 126 |
| 119 | 240 | 399 | 1 019 | 3–62 | 106 |
| 131 | 390 | 356 | 804 | 3–61 | 126 |
| 155 | 510 | 1 281 | 3 526 | 3–103 | 124 |
| 179 | 390 | 259 | 622 | 3–39 | 90 |
| 191 | 240 | 96 | 212 | 3–21 | 50 |
| 239 | 480 | 338 | 840 | 3–41 | 78 |
| 251 | 390 | 291 | 715 | 3–56 | 98 |
| 263 | 1020 | 394 | 1 017 | 3–62 | 72 |
| 359 | 600 | 715 | 1 851 | 3–63 | 74 |
| 431 | 600 | 331 | 764 | 3–28 | 66 |
| 443 | 2652 | 703 | 1 765 | 3–83 | 116 |
| 455 | 780 | 71 | 150 | 3–14 | 34 |
| 491 | 2340 | 835 | 2 023 | 3–91 | 154 |
| 599 | 1020 | 659 | 1 686 | 3–102 | 160 |
| 611 | 1020 | 712 | 1 821 | 3–97 | 130 |
| 791 | 1020 | 898 | 2 237 | 3–85 | 120 |
| 851 | 2460 | 538 | 1 368 | 3–63 | 98 |
| 911 | 1560 | 327 | 777 | 3–47 | 102 |
| 935 | 1020 | 252 | 590 | 3–33 | 60 |
| 959 | 2460 | 513 | 1 233 | 3–57 | 94 |

In each graph:
- the edges are all the pairs of its points at distance 1;
- there is no triangle, and the girth is 4;
- `χ = 4`: there is a stored proper 4-colouring, and no proper 3-colouring;
- the graph is vertex-critical: for every vertex `v` there is a stored proper
  3-colouring of the graph minus `v`.

There is no triangle in any of these planes. A unit triangle would need
`√3 ∈ ℚ(√d)`, which fails for every `d` here.

## 3. How they were found

**The directions.** For a denominator `D`, take every unit vector
`((a + b√d)/D, (c + e√d)/D)` with integer coordinates. They form a finite set
`U_D`, closed under negation, under the quarter turn and under complex
conjugation. We used `|U_D| = 108` in each case (140 for `d = 71`, 148 for
`d = 119` and `d = 155`, 132 for `d = 95` and `d = 239`, 140 for `d = 359`, 116 for
`d = 251` and `d = 431`, 196 for `d = 131`, 212 for `d = 179`, 276 for `d = 35`,
68 for `d = 455`, 84 for `d = 935`, 156 for `d = 263` and `d = 959`, 252 for
`d = 599`, 132 for `d = 611`, 164 for `d = 791`, 276 for `d = 911`, 212 for `d = 443`, 180 for `d = 23` with
`D = 156`). The graph over `ℚ(√11)` came
from a different start (below, **Spindles**).

**The growth.**
1. Start from all sums of at most two vectors of `U_D`.
2. 3-colour the graph (tabucol, else kissat).
3. Add the points `p + u` whose neighbours already see all three colours.
4. Repeat until kissat answers that the graph is not 3-colourable. This took
   1 to 30 rounds and at most 12 211 points (18 041 for `d = 23`, 21 895 for
   `d = 131`, 24 785 for `d = 179`, 39 265 for `d = 35`; 96 rounds for
   `d = 431`; 38 rounds and 15 913 points for `d = 155`, 43 rounds and 12 477
   points for `d = 95`; for `d = 611`, `791` and `911`, 7, 46 and 8 rounds and
   10 580, 17 750 and 40 198 points; for `d = 443`, 15 rounds and 25 429 points; for `d = 491`, 27 rounds and
   22 135 points; for `d = 851`, 11 rounds and 10 992 points). For `d = 23, 35, 95, 131, 155, 179, 239,
   251, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959` we used a variant:
   when no candidate is blocked, it 3-colours the graph again from scratch, and
   if still none is blocked it adds candidates whose neighbours see two colours.
5. Shrink: keep the vertices whose clauses lie in the drat-trim core of a
   checked proof, then delete vertices one at a time while the rest stays not
   3-colourable. The published graphs come from a second pass with one
   incremental SAT solver (CaDiCaL through PySAT): a selector literal per
   vertex, each test a solve under assumptions, the core of failed assumptions
   of every refutation as the new vertex set, and up to 400 random deletion
   orders per grown graph, keeping the smallest critical graph (for `d = 95`,
   `d = 155` and `d = 251` so far one order). For `d = 191`, `239`, `455` and
   `935` a third pass did better (2 October). `min3fast2.py` keeps the vertex
   set of every round of drat-trim cores, and `min3multi.py` started from one
   of these (122, 449, 78 and 315 vertices) instead of from the grown graph.
   For `d = 191` and `455`, whose grown graphs have first cores of failed
   assumptions of 438 and 234 vertices, 394 of 1 000 orders gave 96 vertices,
   and 980 of 3 000 gave 71, where 3 000 orders from the grown graphs had given
   no fewer than 100 and 74. For `d = 239` and `935`, 29 of 729 orders gave
   338 vertices (355 before) and 14 of 1 000 gave 252 (257 before).

**Why some direction sets fail.** Three places act as gates. Let `𝔭` be a prime
of `K` that splits as `𝔓𝔓̄` in `K(i)`, with residue field `𝔽_q`. A unit vector
`z = x + iy` has `zz̄ = 1`, so `v_𝔓(z) = −v_𝔓̄(z)`; when every direction of a set
is integral at `𝔓` and `𝔓̄` (for instance when `p ∤ D`), reduction maps every
graph the set builds into `H_q = Cay(𝔽_q², {(x, x⁻¹)})`.
- **At 2**, when `d ≡ 7 (mod 8)`, the place above 2 splits in `K(i)` with
  residue field `𝔽₂`. Integral unit vectors reduce to the edge `(1, 1)` of the
  graph `H₂ = Cay(𝔽₂², {(1, 1)})`, a perfect matching: with `D` odd every graph
  is bipartite.
- **At 3**, when `d ≡ 2 (mod 3)`, the place above 3 has residue field `𝔽₉` and
  splits in `K(i)`. Integral graphs map to
  `H₉ = Cay(𝔽₉², {(x, x⁻¹) : x ∈ 𝔽₉*})`, which is 3-colourable. So `3 | D` is
  needed, and for `d = 83` it is not enough (below: `D = 1230`).
- **At 5**, when `d ≡ 0, 1, 4 (mod 5)`, the places above 5 have residue field
  `𝔽₅` and split in `K(i)`. Integral graphs map to `H₅`, which is
  3-colourable. So `5 | D` is needed; for `d = 35` and `D = 174`, every
  direction is integral there.

Every denominator of the published graphs meets these conditions. (Corrected on
2 October: we wrote before that with the 44 unit vectors `s/s̄` for `s` with
coefficients in `{−1, 0, 1}` the graphs over `ℚ(√47)` were bipartite because of
the gate at 2. Eight of these vectors, `±(23 ± i√47)/24` and their quarter
turns, are not integral at 2, and `12z + 12z̄ + 23 = 0` for
`z = (−23 + i√47)/24` is a closed walk of odd length 47: the graphs we built
were bipartite only because they were small.)

The graphs `H_q` are those of §12 of `notes/local_colourings.md`. Their
chromatic number is 3 for `q = 3, 5, 9`, and at least 5 for
`q = 13, 17, 25, 29, 37, 41, 49` (kissat, this note).

When `d ≡ 11 (mod 12)`, the vector `s/s̄` with `s = 1 + i√d` is integral at
neither gate. For example:
- for `d = 11`, `s = 1 + i√11` gives `(−5 + i√11)/6`;
- for `d = 47`, `s = 7 + i√47` gives `(1 + 7i√47)/48`. Its length is 1
  because `48² − 47·7² = 1`, the Pell solution of Moorhouse's Theorem 8.6.

The gates do not explain everything. On 25 September we tried the balls of
radius 5 for the 48 unit vectors of height at most 40 over `ℚ(√47)`, up to
2 248 121 points, and they were 3-colourable (`docs/research-log.md`). That set
contains `(−23 + i√47)/24`, from `s = 1 + i√47`, which is integral at neither
gate. With all 108 unit vectors of denominator dividing 240, which include
`(1 + 7i√47)/48`, the colouring-guided growth reached a graph that is not
3-colourable in 30 rounds, at 9 139 points.

**Short odd cycles.** There are no triangles, so the shortest odd cycle of a
graph built from `U_D` has length at least 5, and length 5 exactly when five
vectors of `U_D` have sum 0. For thirteen of the twenty-seven graphs `U_D` has
such a 5-cycle
(`d = 11, 23, 35, 71, 119, 131, 191, 251, 263, 455, 599, 935, 959`), and the
graph contains one; the other fourteen have no odd cycle shorter
than 7 (`data/quadratic_planes/scripts/odd_published.py`). The three smallest
graphs (71, 76 and 96 vertices) are among the thirteen; so is the graph over
`ℚ(√23)`, which a 5-cycle denominator (`D = 156`) shrank from 660 vertices to
393. `scan5.py` lists the denominators with a 5-cycle and every gate open; on
1–2 October they gave every success of the scan (`d = 251`, `455`, `935`,
`263`, `599`, `959`, in 21 to 616 seconds). It is not necessary: on 2 October `D = 1020` settled
`d = 611` and `791`, `D = 1560` settled `d = 911`, `D = 2652` settled `d = 443`, `D = 2340` settled
`d = 491`, and `D = 2460` settled `d = 851`, in 65 to 816 seconds,
although these `U_D` have no 5-cycle (the graphs' shortest odd cycles have length 7). A 5-cycle is not
enough: with such denominators `d = 299`
and `407` (and `263` with `D = 408` and `816`) stopped with no blocked
candidate. For `d = 83`, `107` and `203` there is none with `D ≤ 4 000`,
prime factors at most 61, 40 to 400 directions and every gate open. More
generally (`pentagons.py`), for `d = 47, 83, 107, 203` no unit pentagon at all
has three consecutive sides in `U_D` for the `D` we tried (510, 1020, 1530
and 2142 for `d = 83`), whatever the denominators of the other two sides.

**Spindles.** If two points `A`, `B` of a graph `G` with `|AB| = r` have the
same colour in every 3-colouring, and `4r² − 1 = t²` with `t ∈ K`, then `u = (t
+ i)/(t − i)` has modulus 1, the rotation `z ↦ A + u(z − A)` maps `K²` to
itself, and it moves `B` to a point `B'` with `|BB'| = r |1 − u| = 1`. So `G`,
its rotated copy and the edge `BB'` have no 3-colouring: this is the argument of
the Moser spindle, where `r = √3` and `t = √11`. `targets.py` lists the vectors
of the module spanned by `U_D` at such distances, `growforce.py` grows a graph
until a chosen pair is forced (the growth above, with the extra condition that
the two points have different colours, until that is impossible), and `spin.py`
builds the spindle, checks it exactly and shrinks it. A pair is forced only if
`G` itself is 3-colourable, and at first `growforce.py` did not check that. Over
`ℚ(√11)` with `D = 30` it reported the pair `0`, `m = [−15, 0, 24, −9]` (`|m|² =
(47 − 12√11)/25`, `t = (8 − 3√11)/5`) forced at once, on the 11 612 points
within two steps of `0` or of `m`; but that graph has no 3-colouring by itself,
so nothing was forced. Indeed already the 5 941 points within two steps of `0`
have none: that is where the plain growth stopped, in its first round. It still
gave the published graph. Shrinking the 11 612 points directly, over 400 random
orders, gave 94 vertices; `spin.py` shrank their union with the copy rotated by
`u` (22 481 points, denominator 750), over 300 orders, to 76 vertices and 172
edges, all of them in the rotated copy (2 000 orders found the same graph).
Rotated back by `ū` they lie in the first graph, with `D = 30` (75 of them
within two steps of `0`), and they replace the 94-vertex graph over `ℚ(√11)`.
The 82 vertices over `ℚ(√455)` and 598 over `ℚ(√119)` came the same way. Now
`growforce.py` solves the graph again without the pair, and `spin.py` refuses a
graph with no 3-colouring unless `PLAIN=1`, and writes the result rotated back
when it lies in the copy. Genuinely forced pairs came out over `ℚ(√455)` (5 690
points), `ℚ(√119)` (23 084) and `ℚ(√935)` (9 082 and 9 937); their spindles
shrank to 139, 802, 488 and 521 vertices, more than the published graphs (71,
399 and 257). Over `ℚ(√83)` no pair at such a distance was forced: there is none
two steps from `0` for any `D` tried (two steps `u₁ + u₂` are at such a distance
exactly when `(15 − s²)(1 + s²)` is a square in `K`, where `u₁ū₂ = (1 − s² +
2is)/(1 + s²)`, a curve of genus 1), `growforce.py` stopped with no blocked
candidate for the three kinds of targets three steps away, and unions of a 17
222-point 3-colourable graph with copies rotated so as to share up to 4 305
points stayed 3-colourable (`rotunion.py`, `overlap.py`).

## 4. The certificates

For each `d`, `data/quadratic_planes/q{d}.json` holds the points, the edges, the
fixed edge, the 4-colouring and the 3-colourings of the vertex-deleted graphs.
`q{d}.cnf` is the formula "the graph has a proper 3-colouring in which the fixed
edge gets the colours 0 and 1". This loses no generality: rename the colours.

**Not 3-colourable**, twice for each `d`:
- `q{d}.cnf`, with the variable of vertex `v` and colour `c` at `3v + c + 1`;
- a second formula written by separate code, with the variable at `cn + v + 1`.

In both cases kissat 4.0.4 answered UNSATISFIABLE and drat-trim VERIFIED its
DRAT proof (`q{d}.logs/`), each in a few seconds.

`python3 scripts/verify_quadratic_planes.py` checks the points, the exact unit
lengths, the 4-colourings, the vertex-critical 3-colourings, that each stored
formula is exactly its graph's formula, and the upper bounds of §5. It uses no
library and runs in about a second. With
`--kissat PATH --drat-trim PATH` it also writes the second formula, solves it
and checks the proof, in about a minute for all the fields. With
`--cake-lpr PATH` as well, drat-trim also writes each proof in LRAT form and
cake_lpr checks it, and the same is done for the stored formula `q{d}.cnf`
itself. cake_lpr (Tan, Heule and Myreen) is a proof checker verified in the
HOL4 theorem prover and compiled by the verified CakeML compiler, so this
check does not rest on drat-trim. On 2 October we ran it for every graph:
each proof `VERIFIED UNSAT`, for both formulas
(`data/quadratic_planes/cake_lpr_checks.txt`). The log records the SHA-256 hash
of each stored formula it checked, and a test checks that it is the hash of the
file in the repository.
`tests/test_quadratic_planes.py` runs the fast checks, and checks that this log
covers every field.

**Formal proofs for ten fields.** For `d = 11, 119, 131, 179, 191, 251, 431,
455, 911` and `935`, `lean/Sqrt{d}.lean` proves `χ(ℚ(√d)²) = 4` in Lean 4 with
Mathlib, and depends only on Lean's three standard axioms. For the lower bound
the kernel checks the unit distances of `q{d}.json` and, through Mathlib's
`lrat_proof`, an LRAT proof (`data/quadratic_planes/q{d}.lrat`, from kissat and
`drat-trim -L`) that `q{d}.cnf` is unsatisfiable; a 3-colouring of the plane,
with its colours renamed, would satisfy every clause. Neither the SAT solver nor
drat-trim is trusted. The upper bounds follow §5 (`lean/QuadraticPlanes.lean`).
A valuation subring of `ℚ(√d)` with 7 in its maximal ideal (Chevalley's theorem)
has residue field `𝔽₇`: by Hensel's lemma for `√d ≡ ±s` when `d ≡ s² (mod 7)`,
and by descent on the power of 7 in the denominators when `d = 7d′` (`d = 119, 455`). Unit
vectors are integral there, and the 4-colouring of `𝔽₇²` in
`finite_planes.json` colours the plane. For `d ≡ 3 (mod 8)` (`d = 131, 251`) a
valuation subring with 2 in its maximal ideal has residue field `𝔽₂`, with
`√d − 1` in the maximal ideal; in the coordinates `a = x + y/√d`, `b = 2y/√d`
the squared length is `a² − ab + ((d + 1)/4)b²`, which has no nontrivial zero
over `𝔽₂`, so the residues in `𝔽₂²` colour the plane. The ten fields cover the
three kinds of upper bound and take from half a minute to six minutes each.

**The other graphs in Lean.** For the other fields, and for the lower bound
over `ℚ(√47)`, the LRAT proofs are larger (from 0.7 MB to 70 MB), beyond what
`lrat_proof` checks in reasonable time and memory. Their files prove the
theorem from the hypothesis that the graph's formula is unsatisfiable
(`Sqrt{d}.chromaticNumber_eq_four_of_unsatisfiable`; for `d = 47`,
`not_colorable_three_of_unsatisfiable`). The formula is a Lean object,
`colourCNF E u₀ v₀` (`lean/ColouringFormula.lean`), and
`not_colorable_of_unsatisfiable` proves that a 3-colouring would satisfy it.
When a file is built, `#guard` checks that `q{d}.cnf` is exactly this formula,
and cake_lpr checked an LRAT proof that `q{d}.cnf` is unsatisfiable. So for
every graph the kernel checks the points, the unit distances, the upper bound
(except the bound 5 for `d = 47`) and the step from the unsatisfiability of the
formula to the theorem; `#guard` (an evaluation, not the kernel) checks that the
formula is the stored one; and its unsatisfiability is checked by the kernel
(ten fields) or by cake_lpr (the others).
`lean/tools/field_lean.py` writes each file from the data, and a test checks
that they are up to date.

## 5. The upper bounds

- **`d = 11, 23, 35, 71, 95, 119, 155, 179, 191, 239, 263, 359, 431, 443, 455, 491, 599, 611, 791, 851, 911, 935, 959`.** `d` is a nonzero square
  modulo 7 (`d = 11, 23, 71, 95, 155, 179, 191, 239, 263, 359, 431, 443, 491, 599, 611, 851, 911, 935`) or `d ≡ 0 (mod 7)` (`d = 35 = 5 · 7`,
  `d = 119 = 7 · 17`, `d = 455 = 5 · 7 · 13`, `d = 791 = 7 · 113`, `d = 959 = 7 · 137`), and
  `7 ≡ 3 (mod 4)`. Moorhouse's Lemma 8.2 reduces the plane modulo a prime of
  norm 7 into the unit-distance graph of `𝔽₇²`. That graph is 4-colourable, so `χ(K²) ≤ 4`.
  `data/quadratic_planes/finite_planes.json` holds the colouring, and the
  checker checks it. In Lean: `QuadraticPlanes.colorable_four` and, for `d ≡ 0 (mod 7)`,
  `colorable_four_ramified`.
- **`d = 59, 131, 251`.** `d ≡ 3 (mod 8)`, so `χ(K²) ≤ 4` by Fischer's Theorem 10 (which
  also covers `d = 11`, `35`, `155`, `179`, `443` and `611`). The
  proof in `notes/local_colourings.md` reduces modulo the place above 2, whose
  residue plane is `K₄`. In Lean: `QuadraticPlanes.colorable_four_two`.
- **`d = 47`.** `47 ≡ 3 = 5² (mod 11)` and `11 ≡ 3 (mod 4)`, so reduction modulo
  a prime of norm 11 gives `χ(K²) ≤ χ(𝔽₁₁²) = 5`. Madore's Lemma 4.5 prints a
  5-colouring of `𝔽₁₁²`; `finite_planes.json` holds one, checked by the
  checker. Reduction at 121 does not help (`notes/local_colourings.md`; the
  level-2 plane at 11 is not 4-colourable).
- **Larger fields.** Let `L ⊂ ℝ` be a number field that contains `√d` for one of
  the twenty-six `d` and has a place with residue field `𝔽₇`. Then
  `χ(L²) = 4`: the graph over `ℚ(√d)` lies in `L²`, and since `−1` is not a
  square in `𝔽₇`, unit vectors are integral at the place and reduction gives
  `χ(L²) ≤ χ(𝔽₇²) = 4`. For example `L = ℚ(√11, 7^{1/m})`: 7 splits in
  `ℚ(√11)`, and `x^m − 7` is Eisenstein at each prime above it, so `[L : ℚ] = 2m`
  and those primes are totally ramified in `L`, with residue field `𝔽₇`. So
  there are real number fields of every even degree with `χ(L²) = 4`; fields of
  odd degree have `χ(L²) = 2` (Moorhouse's Theorem 7.1).
- **Planes over `ℚ_p`.** The same two arguments work over the `p`-adic numbers, for the graph on `ℚ_p²` in
  which `(x, y)` and `(x', y')` are adjacent when `(x − x')² + (y − y')² = 1`. Bardestani and Mallahi-Karai
  ([arXiv 1507.05300](https://arxiv.org/abs/1507.05300)) show that its Borel chromatic number is finite exactly when
  `x² + y²` is anisotropic over `ℚ_p`, that is for `p = 2` and `p ≡ 3 (mod 4)`; we know of no exact value in
  the literature. For these `p`, `−1` is not a square in `ℚ_p`, so unit vectors have `p`-adic integer
  coordinates, adjacent points lie in one coset of `ℤ_p²`, and their difference reduces to a nonzero point of
  the circle of `𝔽_p²`. Colouring every coset through a proper colouring of `𝔽_p²` colours `ℚ_p²`, and the
  colouring is locally constant. In the other direction, a unit-distance graph over a number field `K ⊂ ℚ_p` is
  a subgraph of `ℚ_p²`, and `ℚ(√d) ⊂ ℚ_p` when `d` is a nonzero square modulo an odd `p` (Hensel). Hence:
  - `χ(ℚ₂²) = 2` (colour `x + y mod 2` on each coset) and `χ(ℚ₃²) = 3` (`x + y mod 3`; a 5-cycle over
    `ℚ(√7)`, and `7 ≡ 1 (mod 3)`, lies in `ℚ₃²`);
  - `χ(ℚ₇²) = 4`: `𝔽₇²` is 4-colourable, and `11 ≡ 2² (mod 7)`, so the 76-vertex graph over `ℚ(√11)` lies in
    `ℚ₇²`;
  - `4 ≤ χ(ℚ_p²) ≤ 5` for `p = 11` and `p = 19`, since `χ(𝔽₁₁²) = χ(𝔽₁₉²) = 5` and `47` is a square modulo
    both;
  - `χ(ℚ_p²) ≥ 4` for every `p ≡ 3 (mod 4)` from 7 to 79, and `χ(ℚ₈₃²) ≥ 5`: `3`, `11` and `247` are squares
    modulo 83, so the 5-chromatic graph `data/five_247_c.json` over `ℚ(√3, √11, √247)` lies in `ℚ₈₃²`.

  The Borel chromatic numbers have the same values and bounds. For `p ≡ 1 (mod 4)` the graphs over `ℚ(√d)`
  still give `χ(ℚ_p²) ≥ 4` (for example `ℚ(√11) ⊂ ℚ₅`); the Borel chromatic number is infinite there. Since
  `ℚ₁₁` contains `ℚ(√47)` and also `ℚ(√3, √5)` (Voronov's case), `χ(ℚ₁₁²) = 4` would give `χ(ℚ(√47)²) = 4` and
  `χ(ℚ(√3, √5)²) ≤ 4`, while a 5-chromatic unit-distance graph over any number field inside `ℚ₁₁` would give
  `χ(ℚ₁₁²) = 5`. `data/quadratic_planes/scripts/padic_planes.py` checks the ingredients (the stored colourings,
  the 5-cycle, the quadratic residues) and prints the table up to `p = 83`. In Lean (`lean/PadicPlanes.lean`):
  `PadicPlanes.padicSeven_chromaticNumber : (QuadraticPlanes.sumSqGraph ℚ_[7]).chromaticNumber = 4`, and the same
  with 3 for `ℚ_[3]` and 2 for `ℚ_[2]`; the lower bound at 7 holds in any field of characteristic 0 that contains
  a square root of 11 (`not_colorable_three_of_sq_eq_eleven`).

## 6. Open

- **`ℚ(√47)`.** Is `χ = 4` or `5`? It is the smallest open case of Moorhouse's
  table that this note does not settle. The plane is triangle-free, so `χ = 5`
  would give a triangle-free 5-chromatic unit-distance graph. The bound 5 is
  the reduction at a prime above 11, where `ℚ(√47)` embeds in `ℚ₁₁`: every
  finite level of the 11-adic plane bounds `χ`, but none we could test has a
  4-colouring. `(ℤ/121)²` has none (September), and on 2 October neither has
  `(ℤ/1331)²`. Obstructions lift: a 4-colouring of a level restricts to the
  preimage of any set of the level below, so it is enough to refute that
  preimage. A 69-point set of level 1 with no proper 4-colouring has a preimage
  at level 2 with none, which shrinks to 244 points (vertex-critical), and the
  29 524 points of level 3 above those have none either (kissat, drat-trim, and
  an independent second program; `data/quadratic_planes/padic11.json`,
  `scripts/liftcore.py`, `check_lift.py`). `ℚ(√3, √5)` embeds in `ℚ₁₁` too, and its bound
  5 comes from the same plane. Splitting the edges of the 816-vertex graph by
  whether their vector is integral at 2 does not help either: the edges that
  are not integral there contain an odd cycle.

  No other place can give 4 by reduction. Where `−1` is a square in `K_v` (above
  2, 3, 5, 7, 13, 17, …), no colouring of `K_v²` with finitely many colours is
  locally constant: with `z = x + iy`, `w = x − iy` the unit steps are
  `(z, w) ↦ (z + t, w + 1/t)`, so on a line `w = const` each colour class is
  bounded. The remaining places lie above 47 and above the primes
  `p ≡ 3 (mod 4)` with `(47/p) = 1`: 11, 19, 23, 31, 43, 67, …. For `p ≥ 23`
  Hoffman's bound excludes every level at once. The eigenvalues of level 1 are
  Kloosterman sums, at most `2√p`, and those new at a level `ℓ ≥ 2` are at most
  `2/(p + 1)` of the degree (sum over the rotations `≡ 1 (mod p^⌈ℓ/2⌉)` first),
  so the independence ratio of every level is at most `−μ/(1 − μ)`, where `μ` is
  the least level-1 eigenvalue over the degree: 0.2491 for `p = 23`, 0.249997
  for 31, 0.214 for 43, and below `1/4` for every `p ≥ 37`, while a 4-colouring
  has a colour class with ratio at least `1/4`. The ramified place above 47 is
  excluded the same way (residue field `𝔽₄₇`, ratio 0.211). For `p = 11` and 19
  the ratio is 0.286 and 0.276. So a 4-colouring by reduction could only come
  from the 11-adic plane at level 4 or more, or the 19-adic plane at level 2 or
  more; on 2 October that level had no 4-colouring invariant under the rotations
  `t ≡ 1 (mod 19)`, nor one invariant under the translations by `19·(1, 0)`
  (`data/quadratic_planes/scripts/hoffman_padic.py`, `spectrum.py`,
  `levelp.py`, `levelt.py`).
- **All of `d ≡ 11 (mod 12)`.** Does every such field need four colours? For
  `d = 83, 107, 203` our growth has so far found only 3-colourable graphs,
  with up to 180 directions, and on 2 October with `D = 1020` and `2040`
  (`d = 83`), `1560` (`d = 107`) and `1020` (`d = 203`), and for `d = 83` with
  `D = 1722, 3444, 5166`, which contain the factor 82 of its Pell vector
  `(1 + 9i√83)/82` and the factor 7. That may be a limit of the search: over `ℚ(√35)`,
  `ℚ(√131)`, `ℚ(√155)` and `ℚ(√179)` it failed in the same way at first, and
  succeeded with the directions of denominator dividing 390 (also 510 for
  `d = 155`, 480 for `d = 95`) and the variant of the growth.
  Some failed direction sets fail at a gate of §3: for `d = 35` and `D = 174`
  every direction is integral at the places above 5, and for `d = 83` with
  `D = 410` or `1230`, and `d = 107` with `D = 870`, at the places above 3 (a
  divisor 3 of `D` is not enough: over `ℚ(√83)` a direction with 3 in its
  denominator also has 7, 11, 17 or 31 there, and over `ℚ(√107)` 11, 13, 19 or
  23, because the primes above 3 are not principal in `ℚ(√−d)`). For the others
  we found no such reason:
  for `d = 83` with `D = 510` (108 directions) and `D = 1530` (180), no proper
  3-colouring is periodic modulo `mM`, where `M` is the lattice the directions
  span, for `m ≤ 10` or `m = 12` (`D = 510`), `m ≤ 8` (`D = 1530`) and `m ≤ 13`
  (`D = 2958`, 324 directions, where the growth stopped at 2 788 points with no
  blocked candidate) (`data/quadratic_planes/scripts/periodicq.py`).
- **Smaller witnesses.** The graph over `ℚ(√455)` has 71 vertices, and the
  one over `ℚ(√11)` 76. How small can a 4-chromatic unit-distance graph over a
  real quadratic field be?

## 7. Files

| file | contents |
|---|---|
| `data/quadratic_planes/q{d}.json` | the graph over `ℚ(√d)`: points, edges, fixed edge, 4-colouring, 3-colourings of `G − v` |
| `data/quadratic_planes/q{d}.cnf`, `q{d}.logs/` | the formula of §4 and the kissat and drat-trim logs for both encodings |
| `data/quadratic_planes/cake_lpr_checks.txt` | the checker's log with kissat, drat-trim and cake_lpr, for every field |
| `data/quadratic_planes/finite_planes.json` | proper colourings of the unit-distance graphs of `𝔽₇²` (4 colours) and `𝔽₁₁²` (5 colours) |
| `scripts/verify_quadratic_planes.py`, `tests/test_quadratic_planes.py` | the checks of §4 and §5 |
| `lean/QuadraticPlanes.lean`, `lean/ColouringFormula.lean`, `lean/Sqrt{d}.lean`, `data/quadratic_planes/q{d}.lrat` | the formal proofs of §4: for ten fields, and for the other graphs given the unsatisfiability of their formulas |
| `data/quadratic_planes/scripts/` | the code of the search of §3 (growth, shrinking, certification), as it was run, with a README |

## References

- D. Cohen, *The ℂ unit distance graph*, University of Chicago REU paper (2007),
  [math.uchicago.edu/~may/VIGRE/VIGRE2007/REUPapers/FINALFULL/Cohen.pdf](https://www.math.uchicago.edu/~may/VIGRE/VIGRE2007/REUPapers/FINALFULL/Cohen.pdf).
- Y. K. Tan, M. J. H. Heule and M. O. Myreen, *cake_lpr: verified propagation
  redundancy checking in CakeML*, TACAS 2021, LNCS 12652, 223–241.
- K. G. Fischer, *Additive K-colorable extensions of the rational plane*,
  Discrete Math. 82 (1990), 181–195.
- *A 2-adic obstruction to 5-chromatic unit-distance graphs*, public repository
  [MildlyMeticulous/hn-2adic-obstruction](https://github.com/MildlyMeticulous/hn-2adic-obstruction)
  (July 2026), not refereed.
- P. D. Johnson Jr., *Two-colorings of real quadratic extensions of ℚ² that forbid
  many distances*, Congr. Numer. 60 (1987), 51–58 (we have not seen it; its
  result as summarised by Payne).
- D. A. Madore, *The Hadwiger–Nelson problem over certain fields*,
  [arXiv:1509.07023](https://arxiv.org/abs/1509.07023) (2015).
- G. E. Moorhouse, *On the chromatic numbers of planes*, draft of 3 March 2010,
  [ericmoorhouse.org/pub/chromatic.pdf](https://www.ericmoorhouse.org/pub/chromatic.pdf).
- M. S. Payne, *Unit distance graphs with ambiguous chromatic number*, Electron.
  J. Combin. 16 (2009), Note 31; [arXiv:0707.1177](https://arxiv.org/abs/0707.1177).
