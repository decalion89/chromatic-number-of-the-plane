<!--
DRAFT until scripts/verify_g13_chi.py ends with CONFIRMED on the final archive (FINALIZE.md).
Placeholders: each name below stands in the text between double braces.
From `python3 scripts/verify_g13_chi.py --stats` (`--fill notes/g13_chi.md` writes them in):
  F36_TOP_VERIFIED, F35_TOP_VERIFIED, F34_TOP_VERIFIED   leaves of level 0 with a VERIFIED line
  F36_TOP_RESPLIT, F35_TOP_RESPLIT, F34_TOP_RESPLIT      leaves of level 0 split again
  F36_CERTIFIED, F35_CERTIFIED, F34_CERTIFIED            leaves certified in all, re-splits included
  F36_DEPTH, F35_DEPTH, F34_DEPTH                        the deepest level of re-splits
  COLOURING_CERTIFIED     the leaves of F36, F35 and F34 together
  DRAT_PROOFS             the DRAT proofs, checked by drat-trim, that the proof uses: one per certified leaf or
                          formula, the 4 826 formulas of E37 included
  RESPLITS                the re-split formulas in the certified trees
  TREE_TABLE              the rows of the table of re-splits in section 6, level by level
  KISSAT_HOURS, DRAT_TRIM_HOURS, MAX_KISSAT_S, PROOF_GB  on the log lines that certify the leaves
  ALL_LINES, TIMEOUTS, TIMEOUTS_120, TIMEOUTS_1200, TIMEOUTS_3600, ALL_KISSAT_HOURS, ALL_DRAT_TRIM_HOURS
                          over every line of every log
  UNUSED_RESPLITS         cube files of re-splits that the certified trees do not use
  REGEN_GB                the disk space that `share_driver.py regen` takes
  ARCHIVE_FILES, ARCHIVE_MB, ARCHIVE_SHA256              certificates/g13_chi_certlogs.tar.gz
By hand:
  RUN_END                 the day the last share finished, as "28 September"
  MACHINES                the number of machines the shares ran on
  DATE                    the day scripts/verify_g13_chi.py confirmed the archive
  VERIFY_MINUTES          how long that run took
-->

# The anisotropic plane over `𝔽₁₃` needs six colours

`G₁₃` is the graph on `𝔽₁₃²` in which `z ~ w` when `N(z − w) = 1`, where
`N(x, y) = x² − 2y²`. It is the finite plane `G_q` of `local_colourings.md` §4
for `q = 13`, and `notes/g13.md` §1 describes it. This note proves:

**Theorem.** `χ(G₁₃) = 6`.

A proper 6-colouring, found by a SAT solver, is stored in
`data/small_plane_colourings.json`; `tests/test_small_plane_colourings.py` and
`scripts/verify_g13_chi.py` check it on every edge. So `χ(G₁₃) ≤ 6`. The rest of
the note shows that `G₁₃` has no proper 5-colouring. Counting does not show it:
`α(G₁₃) = 36` (`notes/g13.md`), and five classes of 36 points would hold
`180 ≥ 169` vertices.

The proof is by computer: a case split on the size of the largest colour class
(§2), one formula for each case (§3), and cube and conquer, with a DRAT proof
checked by drat-trim for every leaf (§5). The note states the argument, what the
computer checked, how to check it again, and what has to be trusted. The files are
listed in §9. `docs/research-log.md` tells how the proof was found.

This is a computer proof. Nobody outside the project has refereed it.

`G₁₃` is not the finite Euclidean plane `𝔽₁₃²` of Moorhouse (*On the chromatic
numbers of planes*, draft of 3 March 2010, Table 6.1), whose unit circle is
`x² + y² = 1`. That form is isotropic modulo 13, and each vertex there has 12
neighbours, not 14. We have not found the chromatic number of `G₁₃` in the
literature. The sources that `notes/g13.md` names for its independence number do
not give it: Vinh (Electron. J. Combin. 15 (2008), R5, Theorem 12) bounds the
independence number of these graphs only asymptotically, and Vinh's note on
unit-quadrance graphs (arXiv:math/0510092) computes a chromatic number only for
the Euclidean plane over `𝔽₇`.

## 1. The graph

- `G₁₃` has 169 vertices and 14 neighbours per vertex, so 1 183 edges. The code
  writes the vertex `(x, y)` as `v = 13x + y`.
- Identify `𝔽₁₃²` with `𝔽₁₃(√2)`. The 4 732 maps `z ↦ λz + a` and
  `z ↦ λz̄ + a` with `N(λ) = 1` are automorphisms (`notes/g13.md` §1). They form
  a group.
- `lex_order()` in `scripts/g13/chi/g13cnf.py` orders the vertices: 0, then the
  circles `N = 2, 3, …, 12` and last the unit circle `N = 1`, each circle as
  `r, gr, g²r, …, g¹³r`, where `r` is its first point and `g = 3 + 2√2`
  generates the 14 units. Its first 25 vertices are 0, the 14 points of the
  circle `N = 2`, and 10 points of the circle `N = 3`.

## 2. The case split

A *5-colouring* is a map `χ` from the vertices to the colours `{0, 1, 2, 3, 4}`
with `χ(z) ≠ χ(w)` whenever `z ~ w`. Its classes `χ⁻¹(c)` are independent sets;
some may be empty. Let `M(χ)` be the size of its largest class, and, if `G₁₃` has
a 5-colouring, let `M*` be the largest value of `M(χ)` over all 5-colourings.

**Lemma (the case split).** Suppose that `G₁₃` has a 5-colouring, and let
`s = M*`. Then `34 ≤ s ≤ 36`, and there is a 5-colouring with classes
`C₀, …, C₄` (`C_c` the vertices of colour `c`) such that:

1. `C₀` is an independent dominating set of exactly `s` points;
2. read on the first 25 vertices of `lex_order()`, the 0/1 vector of `C₀` is
   lex-greater than or equal to that of `g(C₀)`, for every automorphism `g`
   (*lex-leader*);
3. each of `C₁, …, C₄` has at most `s` points, and each of them with `s` points
   is dominating;
4. if `s = 35` or `36`: the colours `1 < 2 < 3 < 4` have *value precedence*
   along `lex_order()`: colour `c + 1` is used at a vertex only if colour `c` is
   used at an earlier vertex;
5. if `s = 34`: `C₀, …, C₃` have 34 points each and `C₄` has 33, and the colours
   `1 < 2 < 3` have value precedence along `lex_order()`.

This is the lemma in the docstring of `scripts/g13/chi/plan_C.py`, items (1) to
(5).

*Proof.* **The bounds.** The five classes of a 5-colouring cover the 169
vertices, and `5 · 33 < 169`, so some class has at least 34 points: `s ≥ 34`.
Every class is independent, so `s ≤ α(G₁₃) = 36` (`notes/g13.md`).

**Largest classes dominate.** Take a 5-colouring `χ` with `M(χ) = s`, and a class
`C` of `χ` with `s` points. If a vertex `z` outside `C` had no neighbour in `C`,
give `z` the colour of `C`. The result is again a 5-colouring: `C ∪ {z}` is
independent, and the class that `z` leaves stays independent. Its class
`C ∪ {z}` has `s + 1` points, against the choice of `s`. So every class of `χ`
with `s` points is dominating. Give one of them the colour 0: it is `C₀`. Then
(1) and (3) hold, since no class has more than `s` points. If `s = 34`, the
sizes are at most 34 and add up to `169 = 4 · 34 + 33`, so they are 34, 34, 34,
34 and 33; take `C₀` among the classes of 34 points, and give the class of 33
points the colour 4.

**Lex-leader.** An automorphism `g` sends `χ` to the 5-colouring `χ ∘ g⁻¹`,
whose classes are the images `g(C_c)`: they are independent, they keep their
sizes, and dominating classes stay dominating, since `g` maps closed
neighbourhoods to closed neighbourhoods. So (1), (3) and the sizes of (5) hold
for every such image. The images `g(C₀)` are finitely many; choose `g` so that
the vector of `g(C₀)` on the first 25 vertices of `lex_order()` is the greatest
among them, and replace `χ` by `χ ∘ g⁻¹`. For every automorphism `h`, `h(g(C₀))`
is again an image of `C₀`, since the automorphisms form a group, so its vector
is not greater. This is (2).

**Value precedence.** Last, rename the colours `1, …, 4` (if `s = 34`, only
`1, 2, 3`) by the first vertex of their classes along `lex_order()`: the class
that appears first gets the colour 1, the next one 2, and empty classes come
last. Renaming moves no vertex, so (1), (2) and (3) still hold. After it, if
colour `c + 1` is used at the vertex in position `i` of `lex_order()`, its class
first appears at a position `j ≤ i`, and the class of colour `c` first appears at
a position before `j`. This is (4), or (5). ∎

*Why each step loses nothing.* Each step replaces a 5-colouring by another one
whose largest class still has `M*` points: moving one vertex is ruled out, and
an automorphism and a renaming of colours give 5-colourings of the same shape.
The order of the last two steps matters. The automorphism moves every class, so
it comes first; the renaming keeps every class in its place, so it cannot undo
(2). The two do not conflict: the lex-leader condition is on class 0 only, and
value precedence on the colours 1 to 4 only.

So if `G₁₃` had a 5-colouring, one of the three formulas of §3, `F36`, `F35` or
`F34`, would have a solution: the one for `s = M*`. All three are unsatisfiable
(§5, §6). So `G₁₃` has no proper 5-colouring, and with the stored 6-colouring,
`χ(G₁₃) = 6`. ∎

The computation also covers `s ≥ 37`, a class of 37 independent points, with the
five formulas `E37_A6`, `E37_A7`, `E37_A9`, `E37_A11` and `E37_B` of
`notes/g13.md` §3. They are the formulas of that note, with the same SHA-256, and
their logs are the ones in `certificates/g13_alpha_*`: they came from the first
share of this computation. Given `α(G₁₃) = 36`, this part adds nothing new, but
`verify_plan_D.py` checks it too.

## 3. The formulas

The variable `x(v, c) = 1 + 5v + c` says that the vertex `v` has the colour `c`;
auxiliary variables come after the 845 of these. `scripts/g13/chi/g13cnf.py`
writes the clauses:

- *a colouring*: `x(v, 0) ∨ … ∨ x(v, 4)` and `¬x(v, c) ∨ ¬x(v, d)` for each vertex
  `v` and colours `c ≠ d`; and `¬x(u, c) ∨ ¬x(v, c)` for each edge `uv` and each
  colour `c`;
- *counts*, with the totalizer of `notes/g13.md` §4 (Bailleux and Boufkhad
  2003): "class `c` has at least `t` points" is a totalizer over the variables
  `x(·, c)` with the unit clause `o_t` at its root, and "at most `t` points" is
  "at least `169 − t`" of the negations `¬x(·, c)`;
- *class `c` is dominating*: `x(v, c) ∨ x(w₁, c) ∨ … ∨ x(w₁₄, c)` for each vertex
  `v` and its neighbours `w₁, …, w₁₄`;
- *a class `c` of `s` points is dominating*: a counter over `x(·, c)` whose
  outputs are forced both ways, `o_j` if and only if at least `j` points have the
  colour `c`, and the clauses `¬o_s ∨ x(v, c) ∨ x(w₁, c) ∨ … ∨ x(w₁₄, c)`;
- *lex-leader*: the chains of `notes/g13.md` §4 (Crawford, Ginsberg, Luks and
  Roy 1996), on the variables `x(·, 0)` of class 0, one for each of the 4 731
  automorphisms other than the identity, along the first 25 vertices of
  `lex_order()`;
- *value precedence* of colours `a < b` along `lex_order() = (v₀, v₁, …)` (Law and
  Lee 2004): variables `y_i`, defined by clauses to mean that `a` is used among
  `v₀, …, v_i`, and the clauses `¬x(v₀, b)` and `¬x(v_i, b) ∨ y_{i−1}` for
  `i ≥ 1`.

| formula | `g13cnf.py OUT` with | what it says | variables | clauses | SHA-256 |
|---|---|---|---|---|---|
| `F36` | `--big0 36 --cap 36 --dom0 --lex0 25 --vp 1234 --domcap 36` | a colouring; class 0 has at least 36 points, every class at most 36; class 0 is dominating and lex-leader; a class of colours 1 to 4 with 36 points is dominating; value precedence of `1 < 2 < 3 < 4` | 125 416 | 484 970 | `6ebebab0…` |
| `F35` | the same with 35 | the same with 35 | 125 386 | 483 890 | `3341550c…` |
| `F34` | `--rigid34 --lex0 25` | a colouring; classes 0 to 3 have exactly 34 points and are dominating, class 4 at most 33; class 0 is lex-leader; value precedence of `1 < 2 < 3` | 124 192 | 451 390 | `44784afc…` |

`scripts/g13/chi/cases/SHA256SUMS` has the full SHA-256 of these three and of the
five formulas `E37_*`. A 5-colouring as in the lemma, with `s = 36`, 35 or 34,
satisfies every clause of `F36`, `F35` or `F34` once the auxiliary variables get
their intended values: each output `o_j` of a count "at least `j` of its
literals are true", each chain variable "the pairs compared so far agree", each
`y_i` "colour `a` is used up to `v_i`". In
`F34`, class 4 has at most 33 points, and so exactly 33, since every vertex has
exactly one colour.

## 4. The formulas say what they should

- **The encoders.** `g13cnf.py` here is the file of `scripts/g13/`, apart from
  one line of its docstring; its totalizer and its chains also wrote the formulas
  of `α(G₁₃) = 36`. `tests/test_g13.py` checks its encoders by brute force
  against their definitions: the totalizer both ways, the counter forced both
  ways, the chains (also with the automorphisms of `G₁₃`), and value precedence.
- **Every clause in its place.** `tests/test_g13_chi.py` reads the text of
  `F36`, `F35` and `F34`. Each clause on the colour variables alone is one of
  the kinds above, in the right number: at least one colour, at most one, an
  edge clause, the domination clause of a class that must dominate, the first
  comparison of a chain, value precedence at the first vertex. Every other clause
  has auxiliary variables of one count, of the chains, or of value precedence
  only.
- **The counts and the chains, read by the audit of `α(G₁₃)`.** The same test
  hands the counts and the chains of each formula to the functions of
  `scripts/g13/g13_audit.py`, which read clauses alone. Every count of at least
  `t` is a standard totalizer over the 169 variables of its class, with the unit
  clause `o_t` at its root; a count of at most `t` is one of at least `169 − t`
  over the negations. The chains compare the positions of one common order, the
  first 25 vertices of `lex_order()`, each position with its image under an
  automorphism that fixes the positions the chain skips. There is one chain for
  each of the 4 731 automorphisms other than the identity, and the chains are the
  same in the three formulas.
- **Intended models.** The test builds colourings from maximal independent sets
  of 36, 35 and 34 points: the set, moved to its lex-leader image, is class 0; the
  other vertices go greedily to the other classes, which are then renamed by
  first appearance. With every auxiliary variable at its intended value, the
  clauses that fail are exactly the edge clauses of the monochromatic edges and
  the domination clauses of the vertices that a class which must dominate does
  not dominate. For colourings far from this form (random ones, a class 0 that is
  not lex-leader, colours out of order, a class 0 one point short, one colour on
  every vertex), the clauses that fail are exactly those that their meaning
  predicts. For each of 23 such sets (the 15 known sets of 36 points of
  `notes/g13.md`, and four each of 35 and 34 points, in different orbits), the
  lex-leader image satisfies every lex-leader clause and the lex-least image does
  not.
- **Wrong clauses are caught.** The same test changes one clause of `F36` of each
  kind (an edge clause, a clause of a count and its root, a lex-leader clause, a
  clause of value precedence, a guarded domination clause, and a stray unit
  clause) and checks that one of the checks above fails.
- **The review.** An independent review of the plan, written before the run,
  is `scripts/g13/chi/review_plan_D.md`. It found the lemma sound, and noted
  that the lemma as then written did not justify value precedence or the normal
  form of `s = 34`; items (4) and (5) of `plan_C.py` now do. It found the
  encodings sound, with brute-force tests of the encoders, the same code run on
  the smaller planes `G₅` and `G₇`, and `F36`, `F35` and `F34` without the edge
  clauses of colours 1 to 4 found satisfiable with a class 0 in normal form. It
  found that the planned final check of the logs could pass with leaves missing,
  and proposed the tree-aware check that `verify_plan_D.py` makes. The review
  came before plan E (§5): it did not see the change to `cuber2.py`, the
  re-splits named `CASE_splitI`, or `verify_plan_D.py` as it is now.

No program reads `F36`, `F35` and `F34` alone and shows every clause sound, as
`scripts/g13/g13_audit.py` does for the formulas of `α`. The test above reads the
counts, the chains and the clauses on colour variables from the text; for value
precedence and the counters of a class of `s` points, it takes the numbering of
the auxiliary variables from its own reading of `g13cnf.py`, and checks their
clauses on the colourings it builds.

## 5. Cube and conquer

**The cube trees.** `scripts/g13/chi/cuber2.py` cut each case formula into the
leaves of a binary decision tree, with at most 14 decisions: on the variables
`x(v, 0)`, for the vertices `v` along `lex_order()`, for `F36`, `F35` and `F34`,
and on the vertex variables along `lex_order()` for `E37_B`
(`scripts/g13/chi/work/`). The trees are committed:

| cube file | leaves | closed by unit propagation |
|---|---|---|
| `cases/E37_B.icnf` | 4 822 | 564 |
| `cases/F36.icnf` | 4 823 | 565 |
| `cases/F35.icnf` | 4 823 | 564 |
| `cases/F34.icnf` | 4 823 | 564 |

`E37_B.icnf` is `certificates/g13_alpha_part_b_cubes.icnf`; `F35.icnf` and
`F34.icnf` are the same file. Every assignment satisfies exactly one leaf of each
tree: `check_cover()` in `cuber2.py` checks that the leaves are those of a binary
tree that branches on one variable both ways at each node. The closed leaves are
certified like the others.

**A leaf.** The formula of a leaf is the case formula with the literals of the
leaf as unit clauses. `scripts/g13/chi/certify.py` runs kissat on it, which writes
a DRAT proof, and drat-trim, which checks the proof against the formula; then it
deletes the proof. It writes one log line per leaf: the leaf's name, the SHA-256
of its formula, kissat's verdict and time, the size of the proof, drat-trim's
verdict and time. `certify.py` names kissat 4.0.4 and drat-trim 2e3b2dc, which
`scripts/worker_setup.sh` builds; the lines do not record the binaries.

**Re-splits.** kissat had 120 s per leaf. A leaf on which it reached the limit was
split again: its formula became a case of its own, which `cuber2.py` cut with 8
more decisions, and whose leaves were certified in the same way, with 1 200 s
each, recursively. The leaves of a re-split cover every assignment of the leaf
formula, so the leaf formula is unsatisfiable when they all are.

The first re-splits (plan D, named `CASE_leafI`) did not split. `cuber2.py`
branched on the next variable of its list that propagation had not assigned, but
the solver did not report what the unit clauses of the leaf imply. So each of the
8 decisions was on a variable whose value the leaf already implied: one branch
contradicted the leaf and was closed at once, and the other added a unit clause
that the leaf implies. Each such re-split has 8 closed leaves and one open leaf
(`verify_g13_chi.py --stats` counts them), and it only solved the leaf again,
with 1 200 s and then 3 600 s. Plan E (from 27 September, 21:20) passes the unit
clauses to the solver as assumptions, and its re-splits, named `CASE_splitI`,
have up to 256 leaves on new variables. Both kinds are valid certificates: the
leaves of each cover the leaf.

**The rule.** `scripts/g13/chi/verify_plan_D.py` accepts a case formula when its
cube file is a cover and every leaf counts. A leaf counts if a line of the log
with drat-trim VERIFIED has its name and the SHA-256 of its formula, or if it was
re-split, in either way, into a case that is accepted, recursively; the formula of
the re-split must be, byte for byte, the leaf formula. Before this,
`verify_plan_D.py` writes each of the eight case formulas with the code and
compares it with the file from which it hashes the leaves, and it checks that
the four formulas `E37_A_c` have a VERIFIED line. It ends with
`PLAN D FULLY CERTIFIED` only if all of this holds.

## 6. What was computed

The shares ran on {{MACHINES}} machines, from 27 September to {{RUN_END}} 2026;
`scripts/g13/chi/README.md` lists them.

| case | formula | leaves of level 0 | refuted there | split again | leaves certified in all | deepest re-split |
|---|---|---|---|---|---|---|
| `s ≥ 37` | `E37_A_c`, `E37_B` | 4 + 4 822 | 4 + 4 822 | 0 | 4 826 | none |
| `s = 36` | `F36` | 4 823 | {{F36_TOP_VERIFIED}} | {{F36_TOP_RESPLIT}} | {{F36_CERTIFIED}} | level {{F36_DEPTH}} |
| `s = 35` | `F35` | 4 823 | {{F35_TOP_VERIFIED}} | {{F35_TOP_RESPLIT}} | {{F35_CERTIFIED}} | level {{F35_DEPTH}} |
| `s = 34` | `F34` | 4 823 | {{F34_TOP_VERIFIED}} | {{F34_TOP_RESPLIT}} | {{F34_CERTIFIED}} | level {{F34_DEPTH}} |

Level 0 is the cube file of the case. A leaf of level `k` that was split again
became a formula of level `k + 1`, whose leaves are of level `k + 1`:

| case | level | re-split formulas | their leaves | refuted | split again |
|---|---|---|---|---|---|
{{TREE_TABLE}}

- In all, {{DRAT_PROOFS}} DRAT proofs make the proof, one for each certified
  leaf or formula: {{COLOURING_CERTIFIED}} for the leaves of `F36`, `F35` and
  `F34`, in {{RESPLITS}} re-split formulas, and 4 826 for the formulas of
  `α(G₁₃) ≤ 36`.
- On the lines that certify them, kissat took {{KISSAT_HOURS}} hours and drat-trim
  {{DRAT_TRIM_HOURS}} hours; the longest refutation took {{MAX_KISSAT_S}} s. The
  DRAT proofs took {{PROOF_GB}} GB and are not stored.
- The logs have {{ALL_LINES}} lines in all, runs repeated by more than one share
  and timeouts included: kissat reached its limit {{TIMEOUTS}} times
  ({{TIMEOUTS_120}} at 120 s, {{TIMEOUTS_1200}} at 1 200 s, {{TIMEOUTS_3600}} at
  3 600 s). All runs together took {{ALL_KISSAT_HOURS}} hours of kissat and
  {{ALL_DRAT_TRIM_HOURS}} hours of drat-trim. No run found a leaf satisfiable.
- {{UNUSED_RESPLITS}} cube files of re-splits are not used by the certified
  trees. Some leaves were refuted or re-split by more than one share, and
  `verify_plan_D.py` takes a VERIFIED line first, then a re-split of plan D, then
  one of plan E.
- The formulas are not stored: the code writes them again, and each log line
  records the SHA-256 of its formula. The logs and the cube files of the
  re-splits are in
  `certificates/g13_chi_certlogs.tar.gz` ({{ARCHIVE_FILES}} files,
  {{ARCHIVE_MB}} MB, SHA-256 `{{ARCHIVE_SHA256}}`).
- `scripts/verify_g13_chi.py` checked the archive on {{DATE}} (§7).

## 7. How to check it, and what has to be trusted

From the root of the repository, with Python, numpy and python-sat
(`requirements.txt`):

    python3 -m pytest -q tests/test_g13_chi.py         # half a minute
    python3 scripts/verify_g13_chi.py --stats          # {{VERIFY_MINUTES}} minutes
    python3 scripts/verify_g13.py --no-solve           # alpha(G_13) = 36 (notes/g13.md): seconds

The second command checks the stored 6-colouring, checks the archive against
`certificates/g13_chi_SHA256SUMS.txt`, and unpacks it into a temporary copy of
`scripts/g13/chi`. There it runs `share_driver.py regen`, which writes the eight
case formulas with the code, checks them against `cases/SHA256SUMS`, and writes
the formula of every re-split from its parent; and then `verify_plan_D.py` (§5).
It ends with `CONFIRMED` only if `verify_plan_D.py` ends with
`PLAN D FULLY CERTIFIED`, and with `--stats` it prints the numbers of §6. It runs
no solver. The formulas of the re-splits take {{REGEN_GB}} GB of disk while it
runs; `--workdir DIR` puts the copy on another disk. To refute leaves again,
keep the copy with `--keep`, and there run, for example,

    KISSAT=kissat DRAT_TRIM=drat-trim python3 certify.py --log again.log \
        --base cases/F34.cnf --cubes cases/F34.icnf --only 0,1,2

for the leaves 0, 1 and 2 of `F34`, or with `--base cases/F34_splitI.cnf` and
`--cubes cases/F34_splitI.icnf` for the leaves of a re-split. Refuting every leaf
again would take about as long as the run (§6).

What has to be trusted:

- **The argument** of §2 and §3, which is short and uses no computation beyond
  the formulas, and `α(G₁₃) = 36` (`notes/g13.md`).
- **The formulas**: that `F36`, `F35` and `F34` say what §3 says. The tests and
  the review of §4 check this; unlike the formulas of `α`, they have not been
  audited from their text alone.
- **The check of the logs**: `verify_plan_D.py` (about 100 lines) and
  `check_cover()` in `cuber2.py` (a dozen lines). The code that ran the shares
  does not have to be trusted: `verify_plan_D.py` writes every formula again and
  recomputes every SHA-256.
- **The logs.** A VERIFIED line records that drat-trim verified a proof that
  kissat wrote for the formula with that SHA-256. The proofs are not stored; to
  do without the logs, refute the leaves again.
- **drat-trim**, which checked every proof. It is not verified; unlike the
  proofs of `α(G₁₃) = 36`, these have not been checked again by cake_lpr.
- **Python and SHA-256**, which tie each line of the logs to its formula.

This proof has no formal version in Lean.

## 8. Remarks

- The fractional chromatic number of `G₁₃` is `169/36 ≈ 4.69` (`notes/g13.md`),
  so `χ(G₁₃)` exceeds it by more than one: the sizes of the independent sets
  alone would allow five colours.
- The case `s = 34` took most of the work: {{F34_TOP_RESPLIT}} of the 4 823
  leaves of `F34` had to be split again, against {{F35_TOP_RESPLIT}} of `F35` and
  {{F36_TOP_RESPLIT}} of `F36`.

## 9. Files

| file | contents |
|---|---|
| `scripts/g13/chi/` | the code of the computation, copied unchanged from the run; its `README.md` lists the files and the shares |
| `scripts/g13/chi/plan_C.py` | the lemma of §2, in its docstring |
| `scripts/g13/chi/g13.py`, `g13cnf.py`, `enum_cert.py` | the graph and the writers of the formulas: the code of `scripts/g13/` |
| `scripts/g13/chi/cuber2.py` | the cube trees, and the cover check |
| `scripts/g13/chi/certify.py` | the refutation of each leaf: kissat, then drat-trim; one log line each |
| `scripts/g13/chi/share_driver.py` | one share of the computation, and `regen` |
| `scripts/g13/chi/verify_plan_D.py` | the check of §5 |
| `scripts/g13/chi/cases/` | the four cube files, and the SHA-256 of the eight case formulas |
| `scripts/g13/chi/review_plan_D.md` | the independent review of the plan (§4) |
| `scripts/g13/chi/pack_certificates.py` | packs the merged logs into the archive below |
| `scripts/verify_g13_chi.py` | one command that checks everything again (§7) |
| `certificates/g13_chi_certlogs.tar.gz` | the logs of every leaf and the cube files of the re-splits |
| `certificates/g13_chi_SHA256SUMS.txt` | the SHA-256 of the archive and of every file in it |
| `data/small_plane_colourings.json` | the 6-colouring |
| `tests/test_g13_chi.py` | the tests of §4, and of the checks on a small made-up plan |

## References

- O. Bailleux and Y. Boufkhad, *Efficient CNF encoding of Boolean cardinality
  constraints*, CP 2003, LNCS 2833, 108–122.
- J. Crawford, M. Ginsberg, E. Luks and A. Roy, *Symmetry-breaking predicates
  for search problems*, KR 1996, 148–159.
- Y. C. Law and J. H. M. Lee, *Global constraints for integer and set value
  precedence*, CP 2004, LNCS 3258, 362–376.
- M. J. H. Heule, O. Kullmann, S. Wieringa and A. Biere, *Cube and conquer:
  guiding CDCL SAT solvers by lookaheads*, HVC 2011, LNCS 7261, 50–65.
- N. Wetzler, M. J. H. Heule and W. A. Hunt, *DRAT-trim: efficient checking and
  trimming using expressive clausal proofs*, SAT 2014, LNCS 8561, 422–429.
- A. Biere, K. Fazekas, M. Fleury and M. Heisinger, *CaDiCaL, Kissat,
  Paracooba, Plingeling and Treengeling entering the SAT Competition 2020*.
- G. E. Moorhouse, *On the chromatic numbers of planes*, draft, 3 March 2010.
- Le Anh Vinh, *On chromatic number of unit-quadrance graphs (finite Euclidean
  graphs)*, 2005, arXiv:math/0510092.
