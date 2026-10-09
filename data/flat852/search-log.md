Search log: 5-chromatic unit-distance graphs with all edges in the 126 directions U of Q(zeta21)

The search that found the graph, as recorded at the time (September 2026). File names below are those of the
working directory; the published files are in `data/flat852/` (see its README). The solvers were kissat 4.0.4 and
drat-trim; the search took about 3 hours.

HEADLINE
Two certified graphs; neither is vertex-critical yet (deletion minimisation was not run).
1. core852 (best)
   - 852 vertices, 4487 edges, every edge vector in U = mu42 ∪ omega·mu42 ∪ conj(omega)·mu42. All 126 directions occur.
   - Not 4-colourable: kissat 4.0.4 returned UNSAT in 1513 s with a 910 MiB DRAT proof. drat-trim printed "s VERIFIED" in 2037 s (16516 of 18803 clauses and 6.30 M of 10.38 M lemmas in the core).
   - A proper 5-colouring is verified, so the chromatic number is exactly 5.
2. G1023, the graph core852 was cut from
   - 1023 vertices, 5440 edges.
   - kissat UNSAT in 1577 s, 982 MiB proof; drat-trim VERIFIED in 2233 s.
   - Verified 5-colouring; chromatic number 5.
- Both graphs have coordinates in Q(zeta84)^+ and no Moser spindle. A spindle would need sqrt(-11) in Q(zeta21), and a direct search in G1023 finds none.
- Sizes used for comparison: 509 (Parts), 1299 (smallest known without Moser spindles), 2131 (Haugland). 852 is below the last two; the literature comparison is in notes/flat852.md.

1. EXACT MODEL (task 1: flat.py, task1_model.py, task1_model.log, directions.json)
- Representation:
  - An element is an integer 12-vector c meaning (c0 + c1 z + ... + c11 z^11)/7, with z = zeta21, reduced mod Phi21.
  - Conjugation z -> z^20 is an integer 12x12 matrix.
  - 7·omega = (4,-3,1,6,1,-5,1,2,-7,4,1,-2). The identity 7w^12 - 13w^6 + 7 = 0 holds exactly, and omega matches the closed formula to 6e-16.
  - D: u_{2m} = zeta42^m = (-1)^m z^{11m} and u_{2m+1} = zeta42^m·omega.
  - U = D ∪ conj(D) has 126 distinct vectors (84 + 84, with the 42 of mu42 shared).
- Exact checks on the directions:
  - u·conj(u) = 1 for all 126.
  - U is closed under negation, and conj(U) = U.
  - The smallest distance between two directions is 0.0228.
- Mod-2 claim (numerators mod 2, 1/7 read as 1), counted twice (affine solving and brute force over all 4096^2 pairs):
  - D: 42 distinct nonzero residues. 252 ordered pairs (a,b) have (a.r, b.r) != (0,0) for all residues, i.e. exactly 42 row spaces {a, b, a+b}.
  - conj(D): also 42 row spaces, none shared with D.
  - U: 63 residues and 0 pairs, so there is no linear 4-colouring. Confirmed.

2. FIRST TESTS (task 2: task2_balls.py, mkcnf_ball.py; triangle 0, 1, zeta6 fixed)
| Graph | Size | Result |
|---|---|---|
| B1(D) ∪ conj B1(D) = B1(U) | 127 vertices, 252 edges | 4-core is empty, so trivially 4-colourable |
| B2(D) ∪ conj B2(D) | 6175 points (883 shared), 31626 edges; 4-core 5923 / 30996 | Undecided: pysat CaDiCaL no answer in 11 min, kissat UNKNOWN at 1200 s, tabucol best 572 conflicts |
| B2(D) alone (colourable by Theorem A) | 3529 points, 17472 edges | kissat UNKNOWN at 600 s; tabucol stuck at 326 conflicts |
| B2(U) | 7939 points, 4-core 7687, 38052 edges | CNF written (BU2.cnf), not run |
| R5 = {x : d(0,x) + d(i√3,x) ≤ 5} in Cay(M,U), M = L + conj(L) (analogue of Haugland's T5) | 1900 points, 4-core 1886, 8653 edges | SAT in about 21 s; with col(0) = col(i√3) forced, UNKNOWN at 600 s |
- The B2(D)-alone row is a calibration: near-rigid patches are very hard for CDCL even when they are satisfiable, because the solutions are hidden F2-linear maps.

3. GROWTH (task 3: grow.py rounds with pysat CaDiCaL; grow3.py / grow4.py with tabucol, local repair and kissat)
- Growth rule: add candidates x+u that are not in the graph and whose neighbours see all 4 colours. Priority: most neighbours, then in Z[zeta21], then smallest trace norm. At most max(30, min(400, n/12)) per round.
- Seed H ∪ conj(H): Haugland's heptagon centred at the origin plus its mirror, 42 points, 84 edges.
  - Rounds 1–30 (CaDiCaL): 42 -> 806 vertices, solves ≤ 18 s. At 806, CaDiCaL stalled for more than 10 min.
  - kissat rounds:
    - 806: SAT in 47 s, +67 vertices.
    - 873: SAT in 27 s, +72.
    - 945: SAT in 169 s, +78.
    - 1023: time-out at 1500 s inside the loop.
  - At 1023, starting from the 945 colouring:
    - Local repair of balls of radius 0, 1 and 2 around the new vertices (78, 444 and 1004 of 1023 vertices free) was UNSAT each time, in under 1 s.
    - Phase-guided CaDiCaL found nothing in 6.5 min.
    - A standalone kissat run gave UNSAT; this graph is G1023.
  - The kissat colouring at 945 is far from linear: 551 of 1692 edges in the L-coset carry a minority label for their direction.
- Seed G1 ∪ conj(G1) (Haugland's G1 from record_hept, converted from the zeta42 basis): 1441 points, 39 shared, 7896 edges. No first colouring (tabucol failed, kissat time-out at 900 s), so growth could not start.
- Seed B2 ∪ conj(B2): no first colouring (see task 2), so growth could not start.
- G2 ∪ conj(G2) was not tried: record_hept contains only G1, T5 and T6.

4. MINIMISATION (task 4)
- Stage A, clause core:
  - drat-trim -c on the verified G1023 proof keeps 16719 of 22786 clauses.
  - The vertices whose at-least-one clause is in the core (850), plus the fixed triangle, give 853. The 4-core then gives 852 vertices and 4487 edges (core_vertices.py).
  - This subgraph was then certified from scratch (section 5).
- Stage B, one selector per vertex with CaDiCaL (minimize.py): the first full solve did not finish in 33 min, so I stopped it.
- Stage C, deletion-based minimisation to vertex-criticality: not run for lack of time. The code is ready in minimize_k.py (kissat plus drat-trim cores, tabucol, unique-neighbour marking).
- Sizes: 1023 -> 852. Vertex-criticality is not established.

5. CERTIFICATE (task 5: certify.py, gpcheck.gp, five_colour.py)
The same checks were run on both graphs:
- Vertices are distinct.
- Every edge difference is one of the 126 directions, and d·conj(d) = 1 exactly (integer arithmetic mod Phi21).
- An independent PARI/gp check with polmods in Q(zeta21) agrees (gpcheck_g1023.log, gpcheck_core852.log).
- The unit-distance pairs among the vertices are exactly the edges, checked in floats and at 60 digits.
- The CNF rebuilt from the JSON is identical to the certified CNF.
- The kissat and drat-trim results are as in the headline, and each proof was deleted after checking.
- The 5-colourings are verified.
- For core852, a second kissat run with a different fixed triangle was stopped once this one finished.

Shape of core852 in M/O = F7·omega ⊕ F7·conj(omega):
- 758 vertices lie in the L-coset of H, 265 of them in the one O-coset shared with the conj(L)-coset of the mirror.
- 92 more lie in that conj(L)-coset, and 2 elsewhere.
- Edges: 2252 in mu42, 1848 in omega·mu42, 387 in conj(omega)·mu42.
- Degrees 4–34 (mean 10.5); all 42 seed points are kept; the drawing fits in a disc of radius 2.53.
- This is exactly the "L-patch plus conj(L)-patch sharing a Z[zeta21]-coset" picture that the construction was designed around.

6. OBSERVATIONS
- In a flat rhombus (two unit triangles), opposite sides carry the same colour-pair label exactly when the tips, at distance √3, get different colours.
- In the graph with both unit and √3 edges on a triangular lattice, every star neighbourhood is an octahedron K_{2,2,2}. That forces c(x+s) = c(x-s), i.e. the colouring is 2-periodic.
- So Cay(M, U ∪ √3-vectors) is not 4-colourable by a purely local argument: it would 4-colour H64, whose chromatic number is at least 6.
- The expensive part in unit distance is forcing a √3 pair to be bichromatic: T6 does this with 12840 vertices.
- SAT hardness is the real bottleneck. Growth solve times jumped from 27 s to more than 25 min between 873 and 1023 vertices, and each certificate costs about 25 min of kissat plus 35 min of drat-trim on this loaded machine.

7. OPEN
- Vertex-critical minimisation (the deletion stage), and how far below 852 it can go.
- B2 ∪ conj(B2), B2(U) and G1 ∪ conj(G1) remain undecided.
- Other seeds, and symmetric (zeta7-invariant) versions.
- A literature check of the size comparison.

FILES (working directory of the search; see data/flat852/README.md for what is published)
- Best graph: core852a.json (vertices as "exact7" = [c0..c11] over 7 plus float x, y for z -> exp(2 pi i/21); edges [a, b, j] with v_b - v_a = U[j]; U list; fixed triangle; certificate fields), core852a.edges, core852a.cnf, core852a.kissat.log, core852a.drat-trim.log, core852.5colouring.json, core852_pts.npy, gpcheck_core852.log, certify_core852a.out.
- G1023: g1023.json, g1023.edges, g1023.cnf, g1023.kissat.log, g1023.drat-trim.log, g1023.core.cnf, g1023.5colouring.json, g1023_cert_pts.npy, gpcheck_g1023.log.
- Logs: task1_model.log, grow_H1.log, grow_H3.log, grow_G1.log, BD2.kissat.log, R5M*.kissat.log, g806.kissat.log, g1023_lns.log, g1023_phase.log, min_g1023m.log.
- Code: flat.py, seeds.py, cosets.py, labels.py, tabucol.c, lns.py, grow*.py, minimize*.py, core_vertices.py, certify.py, five_colour.py.
- The kissat log files contain only the s-line and exit code (this kissat build prints no statistics); the run times above come from certify.py (core852a) and from file timestamps (G1023).
