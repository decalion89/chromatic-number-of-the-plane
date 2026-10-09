# Haugland's heptagon construction: reproduction and a small test of a third vector orbit

Paper: J. K. Haugland, arXiv 2608.04542v4, Sections 2–4.
The heavy runs took at most about 20 minutes each.
Notation: z = exp(iπ/21), a primitive 42nd root of unity; K = Q(z), of degree 12; L = Z⟨u_0..u_83⟩; U = {u_j}; μ_42 = the 42nd roots of unity.

## Summary
| Item | Paper | Here |
|---|---|---|
| Unit edges of H | 42 | 42 (7 heptagon, 7 + 7 heptagram, 21 triangle); the other 168 pairs have \|d−1\| ≥ 0.0737 |
| Arcs = u_0..u_83, Table 1 | yes | yes, to 1e-60 and exactly in K |
| θ | 0.42363201413287 | 0.423632014132868551324320426773 |
| 3-ball of Cay(L, U) | 83581 | 83581 (spheres 1, 84, 3444, 80052) |
| \|V(T5)\| | 1042 | 1042 (4597 edges) |
| \|V(T6)\| | 12856 | 12856 (87619 edges) |
| G1 (7-core of G0) | 740 v, 3985 e | 740 v, 3985 e (G0: 1294 v, 6727 e) |
| T5, col(A) = col(B) | satisfiable | SAT (colouring verified) |
| T6, col(A) = col(B) | unsatisfiable | UNSAT (486 s CPU) |
| G1, col(A) = col(B) | unsatisfiable (CaDiCaL) | UNKNOWN (plain CNF: 1046 s CPU; triangle fixed: 931 s CPU); not reproduced |
| Other unit vectors in L | – | none (complete enumeration) |

## 1. H and the 84 vectors (step1_H.py, step1.log, step1.json)
The points were built from the paper's formulas with mpmath at 80 digits. There are exactly 42 unit pairs: the heptagon, the heptagrams {7/2} and {7/3}, and the triangles P_jQ_jR_j.

Of the 84 arc angles (in units of π/21), 42 are integers and 42 have fractional part θ, where θ = 21φ/π − 14 as in the paper. Each arc equals exactly one u_j, and the assignment is Table 1.

Exact form in K:
- P_0 = −1/(z³−z⁻³), Q_0 = z¹⁴/(z⁶−z⁻⁶), R_0 = z²⁸/(z¹²−z⁻¹²), and P_j = z^{6j}·P_0, likewise for Q_j and R_j. These agree with the formulas to 1e-60.
- The 42 pairs are exactly unit (x·x̄ = 1 in K) and no others are.
- u_{2j} = z^j and u_{2j+1} = z^j·u_1, with
  - u_1 = (2+z+3z²−2z⁴+z⁵+5z⁶+z⁷+z⁸−6z⁹−4z¹⁰+4z¹¹)/7 = e^{−2πi/3}·u_29,
  - u_29 = P_0 − Q_0 = i(2cos(π/7) + e^{2πi/3})/(2 sin(2π/7)).
- So U = μ_42 ∪ μ_42·u_1. The common denominator on the power basis is 7, N(u_1) = 1 and N(7u_1) = 7¹².

## 2. Point identification, T5, T6 (tgraph.py, step2_T.py, step2.log)
**Construction.**
- V(T_n) = {p : d_A(p) + d_B(p) ≤ n}, with A = 0 and B = u_14 + u_28 = (0, √3) (checked exactly).
- Since d(A,B) = 2, d_A ≤ d_B + 2 and vice versa. So only the radius-3 BFS ball R3 is needed: T5 ⊂ R3 ∩ (B+R3), and T6 also uses R2 and B+R2, where the unknown distance must be exactly 4.
- Edges are the pairs whose difference is in U.

**Identification (two independent methods, asserted to agree at every one of about 1.8 M lookups).**
1. Numeric: 60-digit fixed-point integer coordinates from the paper's angle formula; two points are identified iff both coordinates agree within 1e-30.
   - Rounding: each vector is rounded once (≤ 0.5e-60), and a compared difference involves at most 11 roundings, so the error is ≤ 6e-60.
   - Separation: if 0 ≠ x ∈ L is a sum of n of the u_j, then 7x is an algebraic integer of K. Every conjugate of a u_j has modulus 1 (K is CM), so 1 ≤ |N(7x)| ≤ |7x|²(7n)¹⁰, hence |x| ≥ 7⁻⁶n⁻⁵.
   - With n ≤ 11, distinct points are more than 5e-11 apart. The tolerance 1e-30 is about 20 orders of magnitude from both failure modes.
   - Observed minimum distance between distinct vertices: 0.0129 (T5) and 0.00348 (T6).
2. Exact: integer 12-tuples, 7 × the power-basis coordinates.

**Results.**
- T5 (d_A, d_B) profile: (0,2):1, (1,1):2, (1,2):2, (1,3):80, (2,0):1, (2,1):2, (2,2):158, (2,3):358, (3,1):80, (3,2):358.
- T6 adds (2,4):2925, (3,3):5964, (4,2):2925.
- A float all-pairs check (|d²−1| < 1e-9) finds exactly the Cayley edges in T5 and T6, with none extra and none missed. Section 4 shows no others can exist anywhere in L.
- G0 and G1 were rebuilt as in the paper's Section 3. G1 has 740 vertices and 3985 edges and contains A and B.

## 3. SAT (step3_sat.py, *.cnf, *.kissat.log, kissat 4.0.4)
**Encoding.**
- Variable 4v+c+1 means "vertex v has colour c".
- Per vertex: one at-least-one clause and 6 at-most-one clauses.
- Per edge and colour: (¬x_uc ∨ ¬x_wc).
- A = B: (¬x_Ac ∨ x_Bc) for c = 0..3.

**Sound options.**
- --core4 repeatedly deletes vertices other than A and B of degree ≤ 3; any colouring of the rest extends to them.
- --symbreak fixes colours 0, 1, 2 on a triangle A, p, q.

Every SAT model is decoded and checked on the whole graph (the 4-core colouring is extended greedily).

**Results.**
- (a) T5, A = B, plain CNF (4168 variables, 25686 clauses): SAT in 140 s wall (CPU not recorded). The colouring is in T5_AeqB.colouring.json. The pair is not forced apart in T5, as the paper says.
- (b) T6, A = B, 4-core (12840 vertices, 87577 edges), triangle fixed (51360 variables, 440195 clauses): UNSAT in 516 s wall / 486 s CPU. The pair is forced apart in T6, as the paper says. No DRAT proof was checked.
- Extra, G1 with A = B: UNKNOWN (1200 s wall / 1046 s CPU plain; 1150 s wall / 931 s CPU with the triangle fixed). The paper's G1 claim is not reproduced; the T6 result does not depend on it.

## 4. Other unit vectors, and the third-orbit idea (lat.py, step4_units.py, step4_T5.py, step4_hom.py, sat_cadical.py)
**Short sums.**
- For a vector set closed under ±60° rotation, a length-1 sum of 2 or 3 vectors is never new.
  - Two vectors: a + b = v unit forms an equilateral triangle, so v = a·e^{±iπ/3} ∈ U.
  - Three vectors: a+b+c−v = 0 with four unit vectors forces them to cancel in pairs, so v ∈ {a, b, c}.
- Computed: the 3-ball has exactly 84 points of length 1, the u_j. All sums of at most 5 of the u_j with length 1 (float filter, then exact check) are again only the 84 u_j.

**Complete: L has no other unit vectors.**
- A unit vector x ∈ L has Tr_{K/Q}(x·x̄) = 12. On power-basis coordinates this form has the integer Gram matrix G_jk = c_42(j−k) (Ramanujan sums).
- Fincke–Pohst over an LLL-reduced basis of L finds exactly 85 points with Tr ≤ 12: the origin and the 84 u_j. Also [L : Z[z]] = 7.
- So every unit-distance graph on points of L (T5, T6, G1, …) is exactly the Cayley graph of U.
- Explanation, consistent with the computation:
  - 7O_K = 𝔭⁶𝔭̄⁶, and complex conjugation swaps 𝔭 and 𝔭̄.
  - (u_1) = 𝔭𝔭̄⁻¹ (with a suitable naming of 𝔭), the exponent being pinned by the index 7, so L = 𝔭̄⁻¹.
  - A unit vector of K with 7-power denominator is ζ·u_1^m with ζ ∈ μ_42 and m ∈ Z (Kronecker's theorem). In L only m ∈ {0, 1} occur.

**Natural third orbits.**
- Rotating H by a multiple of π/21 gives the same 84 arcs. Rotating H by the angle of an odd u_j does not: it gives 42 old arcs plus W_2 = μ_42·u_1², at angles π(k+2θ)/21.
- The mirror image of H gives 42 old arcs plus W_−1 = μ_42·ū_1, at angles π(k−θ)/21.
- Both statements were checked numerically.
- By the complete enumeration, the unit vectors of each enlarged lattice are exactly its generators:

| Lattice | Index over L | Points with Tr ≤ 12 | Unit vectors |
|---|---|---|---|
| L + W_2 | 7 | 631 | exactly the 126 generators |
| L + W_−1 | 7 | 631 | exactly the 126 generators |
| L + W_2 + W_−1 | 49 | 4453 | exactly the 168 generators |

**T5 with the enlarged sets.**
- Same A and B. The float all-pairs check found no unit pairs outside the edges.
- The A = B runs used the 4-core plus a fixed triangle.

| Vectors | 3-ball | Vertices | Edges | Edges using W | col(A) = col(B) |
|---|---|---|---|---|---|
| U (84) | 83581 | 1042 | 4597 | – | SAT |
| U+W_2 (126) | 292909 | 1552 | 6631 | 1110 | SAT, 59 s CPU, verified |
| U+W_−1 (126) | 292909 | 1900 | 8689 | 2322 | UNKNOWN (see below) |
| U+W_2+W_−1 (168) | 718705 | 2458 | 10903 | 3612 | UNKNOWN (kissat 958 s CPU) |

- U+W_−1 attempts: kissat 1097 s CPU; kissat --unsat 937 s CPU; CaDiCaL 1.9.5 via pysat was killed at 1250 s wall because its interrupt did not take effect.
- Without the A = B constraint, U+W_−1 T5 is 4-colourable: SAT in 147 s CPU, verified, with A and B coloured differently.
- The U+W_−1 graph is an induced subgraph of the 168-vector graph. Its vertices lie in L+W_−1, whose only unit vectors are the 126 generators, so a SAT answer for the larger graph would have settled W_−1.

**Homomorphic colourings (step4_hom.py).**
- Take colours in F_2². A homomorphism φ: M → F_2² that is nonzero on every generator is a proper 4-colouring x ↦ φ(x) of the whole lattice graph. Its edge labels are the paper's A, B, C, and it always has φ(B) − φ(A) = φ(u_0) ≠ 0, because u_14 = u_0 + u_28.
- All six labellings of Table 2, extended by the shift A→B→C→A for u_{j+14}, are given by such homomorphisms. Each is therefore a 4-colouring of all of L. On T6, all 87619 edges are proper and carry the Table 2 label.
- Counts: L has 252 such φ (42 up to label permutation). L+W_2, L+W_−1 and L+W_2+W_−1 have none.

**Assessment.**
- Inside L the idea cannot work: there are no new unit vectors.
- W_2 does not shorten the forcing path to length 5: T5 with 1552 vertices is SAT with A = B.
- W_−1 (1900 vertices) and W_2+W_−1 (2458 vertices) are unresolved; the timeouts are not evidence of UNSAT.
- The loss of all homomorphic colourings in the enlarged lattices suggests stronger forcing there, possibly even failure of 4-colourability on larger balls. This was not tested.
- In these runs the larger T6 was refuted in 486 s while its subgraph G1 was not, so a T6-sized test for U+W_−1, or cube-and-conquer, is the natural next step.
- Side fact: every unit-distance graph with vertices in K is Moser-spindle-free. A spindle needs a rotation with cos φ = 5/6, which would put √−11 in K, but K's only quadratic subfields are Q(√−3), Q(√−7) and Q(√21).

## 5. Files (in the working directory of this reproduction, not published)
**Modules.**
- field.py: exact arithmetic in K.
- hept.py: H and the u_j, numeric and exact, plus Table 1.
- tgraph.py: balls, T_n, dual identification, all-pairs check, JSON writer.
- lat.py: HNF, LLL, Fincke–Pohst, unit vectors of a lattice.

**Scripts** (each with a .log): step1_H.py, step2_T.py, step3_sat.py, step4_units.py, step4_T5.py, step4_hom.py, sat_cadical.py.

**Graphs:** T5.json, T6.json, G1.json, T5_W2.json, T5_W-1.json, T5_W2_-1.json.
- Each vertex has x, y (32 digits), exact (7 × power-basis coordinates), dA, dB, and word (a shortest path from A as vector indices).
- Each edge is [a, b, j] with v_b − v_a = generator j. Generators 0..83 are the paper's u_j, and 84+42i+k is z^k·u_1^{m_i}.

**Also:** *.cnf (the exact instances), *.kissat.log, *.colouring.json (verified SAT models), step1.json, step4_units.json.

**Not verified here:**
- No DRAT proofs.
- Lemma 2.2's ε values were not re-derived; the separation bound replaces them.
- G1's forcing claim was not reproduced, and G2, G3 and the 2131-vertex graph were not built.
