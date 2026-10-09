# Haugland's heptagon family at the places above 2 (28 September 2026)

Subject: the 84 unit vectors of J. K. Haugland, arXiv 2608.04542v4, §§2–3 (graph H, Tables 1–2, the graphs G1, G2, G3 with 740, 1066 and 2131 vertices). The tool is the reduction principle of notes/local_colourings.md (Propositions A and C). The computations were made in a working directory that is not published. §7 separates what is proved from what is not.

## 0. Summary
- **Fields.** The coordinates of H generate K = ℚ(ζ84)⁺ = ℚ(cos 2π/7, √3, √7), of degree 12, and K(i) = ℚ(ζ84). H itself (as complex numbers) and all 84 unit vectors lie in the index-2 subfield F = ℚ(ζ21), which does not contain i.
- **Directions.** u_{2m} = ζ42^m and u_{2m+1} = ζ42^m·ω, where ω = u1 = Q0 − R0 = ζ3/(ζ7 − ζ7⁻¹) − ζ3²/(ζ7² − ζ7⁻²).
  - minpoly(ω) = 7x¹² − 13x⁶ + 7.
  - ω⁶ = (13 + 3√−3)/14, and θ = (7/2π)·arctan(3√3/13).
  - (ω) = 𝔔/𝔔̄ with 𝔔 | 7.
- **Places above 2.** 2 places of K(i), with e = 2 and f = 6, swapped by conjugation; K has one, and it splits in K(i). The residue field is 𝔽64, not 𝔽4 as in Parts' family.
- **Flatness.** All 84 directions are units at both places. They reduce to μ21 (even) and r·μ21 (odd), with r = res ω of order 9. The coset r⁻¹μ21 is missed.
- **Reduction principle.** Every unit-distance graph with edges in the 84 directions is 4-colourable. The residue Cayley graph has exactly 42 proper 4-colourings, all linear. Up to symmetry there are 6, and they are exactly Haugland's six Table 2 labellings, so those extend to all of L.
- **Consistency.** G3 escapes through ρ = (7 + i√15)/8 ∉ ℚ(ζ84), which has valuation ±2 at all four places above 2 of ℚ(ζ21, √5). This is the same flat/deep split as Parts' 509.
- **Design rule.** Forced-monochromatic pairs inside L differ by elements of 2L, so spindles sit on 2L-arms, and their rotations are automatically deep. Unit rotations are not screened out automatically here (q = 64): at least 6 missing-coset classes are needed to silence the test.

## 1. The field and ω (field.gp, field.log)

- **Points of H as complex numbers.** They agree with the paper's formulas to 1e−76, and all lie in F.
  - P_j = ζ7^j/(ζ7⁴ − ζ7⁻⁴)
  - Q_j = ζ3ζ7^j/(ζ7 − ζ7⁻¹)
  - R_j = ζ3²ζ7^j/(ζ7² − ζ7⁻²)
- **The graph and Table 1.**
  - Exactly 42 unit-distance pairs among the 21 points (7 + 7 + 7 + 21).
  - Table 1 (indices mod 7) gives exactly u_k = ζ84^k for k even, and u_k = ζ84^(k−1)·ω for k odd. Checked for all 84 arcs by exact arithmetic.
  - So the direction set is D = μ42 ∪ ωμ42.
- **ω.**
  - minpoly(ω) = 7x¹² − 13x⁶ + 7. The stabiliser of ω in (ℤ/84)^× is {1, 43}, so ℚ(ω) = F.
  - ω⁶ = (13 + 3√−3)/14. This is also the spindle rotation for arm √7, since cos = 13/14.
  - θ = (7/2π)·arctan(3√3/13) = 0.42363201413286855…, matching the paper and arg(u29)·21/π − 14.
  - In the power basis of ζ = ζ21: ω = (4 − 3ζ + ζ² + 6ζ³ + ζ⁴ − 5ζ⁵ + ζ⁶ + 2ζ⁷ − 7ζ⁸ + 4ζ⁹ + ζ¹⁰ − 2ζ¹¹)/7.
- **Coordinate field.** The stabiliser of all 42 coordinates of H is {±1}, so K = ℚ(ζ84)⁺. A random combination of the coordinates has degree 12.
- **The lattice** L = ℤ[ζ21] + ℤ[ζ21]ω (lattice.gp).
  - Rank 12, index 7 over ℤ[ζ21]. As a fractional ideal, L = 𝔔̄⁻¹.
  - A unit vector x of F has |σx| = 1 in all 12 embeddings, so Tr(x x̄) = 12.
  - qfminim finds exactly 84 vectors with Tr ≤ 12: the 84 directions, which are also the minimal vectors of L.
  - So every unit-distance graph with vertices in one translate of L has all its edges in D.

## 2. Places (places.gp; e, f and counts from idealprimedec)

"Swapped" means complex conjugation exchanges the places.

| p | F = ℚ(ζ21) | F⁺ | K(i) = ℚ(ζ84) | K = ℚ(ζ84)⁺ | residue field of K(i) |
|---|---|---|---|---|---|
| 2 | 2 places, e=1, f=6, swapped | 1, e=1, f=6 | 2 places, e=2, f=6, swapped | 1, e=2, f=6 (splits in K(i)) | 𝔽64 |
| 3 | 1 place, e=2, f=6, fixed | 1, e=2, f=3 | 2 places, e=2, f=6, each fixed | 2, e=2, f=3 (each inert) | 𝔽729 (𝔽27 for K) |
| 7 | 2 places, e=6, f=1, swapped | 1, e=6, f=1 | 2 places, e=6, f=2, swapped | 1, e=6, f=2 | 𝔽49 |

- **Names.** 𝔓1 = (2, f1(ζ21)) and 𝔓2 = 𝔓̄1 = (2, f2(ζ21)), with f1 = x⁶+x⁵+x⁴+x²+1 and f2 = x⁶+x⁴+x²+x+1.
- **Normalisation.** v(2) = 1 in F. In K(i) all 2-adic valuations double; for elements of F the residues are the same.

## 3. The 84 directions above 2 (classify2.py, crosscheck.gp, classify2.json)

- **Units.** v(u_j) = 0 at both places, for all j.
- **Residues.** With X = res ζ21 (order 21) and r = res ω: res u_{2m} = X^{11m} and res u_{2m+1} = r·X^{11m}. The same form holds at 𝔓2, and res ū = (res u)⁻¹.
  - r has order 9 at both places, so r ∉ μ21 and r³ = X¹⁴.
  - Exactly: ω³ ≡ ζ3² (mod 2) with v = 1, and ω⁶ ≡ ζ3 (mod 4) with v = 2, because ω⁶ − ζ3 = 4(2 − ζ3)/7.
- **Cosets.** 𝔽64^×/μ21 ≅ ℤ/3.
  - The even directions fill μ21.
  - The odd directions fill r·μ21.
  - The coset r⁻¹μ21 (containing res ω̄) is hit by nothing.
- **Classes.** 42 distinct residues, and each class is an antipodal pair. Independent PARI check over all 3486 pairs: v(u_j − u_k) = 1 for the 42 antipodal pairs and 0 for the other 3444.
- **Depth.** No direction is deep.
  - Even u: v(u − w) > 0 for w ∈ μ42 only when w = ±u, with v = 1. In K(i) also w = ±iu, with v(1 − i) = 1 against v(2) = 2.
  - Odd u: max over w ∈ μ42 of v(u − w) is 0.
  - Only u0 = 1 and u42 = −1 have residue 1.
- **At 3** (non-split, so Proposition A applies).
  - All 84 are units. The even ones reduce to μ14; the odd ones to the other coset of μ14 in N1 = μ28 (res ω has order 4). Together they fill all of N1.
  - λ_min(G27) = −8, so Hoffman gives χ(G27) ≥ 4.5; SAT also finds no 4-colouring. The 3-adic screen is silent.
  - Whether χ(G27) = 5 is undecided (SAT stopped at about 6 minutes).
- **At 7.** The odd directions have v = ±1. The even directions fill 𝔽7^×, so the even ones alone map to H7 with χ(H7) = 4 (SAT); by the reflection, so do the odd ones alone.

## 4. The reduction principle (cayley2.py, orbits2.py, sat2.py, certificate.py)

Let T = μ21 ∪ rμ21, and Γ_T = Cay(𝔽64², {(t, 1/t) : t ∈ T}), with 4096 vertices and degree 42.

**Theorem A.** Every unit-distance graph whose edge vectors all lie in the 84 directions is 4-colourable.

*Proof.*
1. L/2L ≅ O/2O ≅ 𝔽64², because L/ℤ[ζ21] is 7-torsion and 2 is unramified in F. The map g(z) = (res z, res z̄) induces this isomorphism, and g(u) = (t, 1/t).
2. A linear φ: 𝔽64² → 𝔽2² that vanishes on none of the 42 points gives a homomorphism c = φ∘g : L → 𝔽2² with c(u) ≠ 0 for every u ∈ D.
3. Colour each translate p + L by z ↦ c(z − p).

Exactly 42 such φ exist. They were counted twice by independent code: via the residue fields at both places, and field-free via power-basis coordinates mod 2. ∎

- **Certificate.** c(z) = M·(power-basis coordinates of z mod 2), with 1/7 read as 1 and M a 2×12 matrix over 𝔽2.
  - Haugland's labelling 1 has M = [110100111101; 001111111010], listing coefficients of ζ⁰…ζ¹¹.
  - Properness reduces to M·(u_j mod 2) ≠ 0 for the 84 vectors. All six matrices are in certificate.json.
- **No other colourings.** With a triangle pinned and the 42 affine colourings blocked, CaDiCaL returns UNSAT in 0.9 s. So Γ_T has exactly 42 4-colourings up to colour permutation.
  - Hoffman is tight: λ_min(Γ_T) = −14, λ_max = 42.
  - The full H64 has λ_min = −13, so χ(H64) ≥ 6: unrestricted Proposition C at q = 64 is useless for four colours.
- **Symmetry and Table 2.** The isometries preserving D are rotations by ζ42^k (j ↦ j+2k) and reflections z ↦ ζ42^k ω z̄ (j ↦ 1+2k−j, swapping even and odd). On residues they form a dihedral group of order 42.
  - The 42 colourings fall into 6 orbits of 7.
  - In every one, the ζ6 rotation (j ↦ j+14) acts on the labels {A, B, C} = 𝔽2²∖0 as a 3-cycle. This is Haugland's universal shift. Antipodes get equal labels (characteristic 2).
  - Each of Haugland's six Table 2 rows, extended by the shift A→B→C, is a homomorphism L → 𝔽2², and the six rows lie in the six different orbits.
  - So his colourings extend to all of L. With his "at most 6" claim (his SAT run), L has exactly 6 4-colourings up to isometry and colour permutation.
- **Corollary B (forced pairs).** Let G be a finite unit-distance graph with vertices in one translate of L.
  1. If x, y are monochromatic in every 4-colouring of G, then y − x ∈ 2L.
  2. If they are bichromatic in every 4-colouring, then y − x ∈ u_j + 2L for some j.

  *Proof.* The only class in all 42 kernels is 0, and the classes in no kernel are exactly the 42 edge classes (design2.log). ∎

  **Converse**, assuming Haugland's G2 claim:
  - Each d ∈ 2L is forced monochromatic by a finite subgraph of L: 2L is generated by the 2u_j, each 2u_j is a symmetric image of his pair 2 = (1,0) − (−1,0), and a chain of forced-monochromatic pairs is forced.
  - Each d ∈ u_j + 2L is forced bichromatic, by such a chain plus one edge.
  - Haugland's pairs fit: B − A = √−3 = 1 + 2ζ3 ≡ u0.

## 5. Consistency with the 2131-vertex graph (rho.gp)

- **G1** (in L, 740 vertices). Not 4-colourable with A = B. Every linear colouring separates A and B, since g(B − A) = (1,1) = g(u0). Consistent.
- **G2** (in L, 1066 vertices). (−1,0) and (1,0) are monochromatic in every 4-colouring. Every linear colouring agrees, since g(2) = 0. G2 itself is 4-colourable by Theorem A.
- **G3 = G2 ∪ (ρ(G2+1) − 1)**, with ρ = (7 + i√15)/8. The matrix (1/8)[[7, √15], [−√15, 7]] acting on row vectors is multiplication by ρ. G3 lies in E = ℚ(ζ21, √5), of degree 24.
  - E has four places above 2, each with e = 1 and f = 6, swapped in pairs by conjugation.
  - v(ρ) = ±2 at all four places.
  - All 6264 edges of the rotated copy have direction ρu_j, so they are deep at every place above 2.
  - The cross edge (1,0)–(3/4, √15/4) has vector 2(ρ−1) = (−1+i√15)/4, with v = ±1. The second cross edge was not identified.
  - G2's 6264 edges are units at all four places.

  This is Parts' mechanism (L374 flat, ρS136 deep) with the same ρ. Adjoining √5 also splits the place above 3; this is not needed, since G27 is not 4-colourable.
- **The requested check.** The 84 directions do reduce to a 4-colourable configuration. That is correct, not an error: Haugland's final graph is not a graph in those directions, and his Table 2 is exactly this set of 2-adic colourings.

## 6. What a 5-chromatic graph from these directions must contain

**Proposition D** (Proposition C with the residues actually used). Let G have its vertices in a CM field E ⊇ F, let W be a place of E above 2 with residue field k_W, and suppose all edge vectors are W-units with residue set T_G. If Cay(k_W², {(t,1/t) : t ∈ T_G}) is 4-colourable, then so is G.
- For k_W = 𝔽64 and T_G ⊇ T, the residue graph is 4-colourable **iff** one of the 42 kernels avoids all (t, 1/t) with t ∈ T_G∖T. This uses the UNSAT result of §4.
- T_G∖T lies in the missing coset r⁻¹μ21.
- Each missing-coset class lies in 12 of the 42 kernels.
- At least 6 classes are needed to kill all 42 (5103 such 6-sets among 54264). Three sampled 6-sets were confirmed UNSAT by SAT.
- Grouped into the 7 ζ3-orbits (the 7 triangles of the mirror image H̄), 6 of the 7 are needed.

So, at every place W above 2 of its field, a 5-chromatic graph built from the 84 directions needs one of:
- **(a)** a W-deep edge vector;
- **(b)** unit edges covering at least 6 missing-coset classes in a killing pattern;
- **(c)** residues outside 𝔽64 (residue degree above 1 over 𝔓). Not analysed.

Inside F the two places are conjugate, so a vector deep at one is deep at both. In larger fields, every place above 2 must be covered (as for Parts' p₊ and p₋).

**Design rule.**
1. **Spindles on 2L-arms are automatically deep.** By Corollary B, a base inside L can force monochromatic pairs only at d ∈ 2L.
   - If |d|·|σ − 1| = 1 and σ were a W-unit, then e = d(σ − 1) and ē would both have valuation at least 1, contradicting e ē = 1. So σ is non-integral at every place above 2.
   - σ ∈ F exactly when 1 − 4d d̄ is a square in F. For |d| = 2 this is −15, not a square in F: Haugland's ρ needs √5.
   - Among the 897 arm lengths with d = 2·(sum of at most 3 directions), exactly three have σ ∈ F (spindles.gp):

     | arm d | \|d\|² | σ | v(σ) |
     |---|---|---|---|
     | 2(u0+u8+u40) | 14−2√21 ≈ 4.835 | in F | ∓2 |
     | 4 | 16 | (31+3i√7)/32 | ±4 |
     | 2(u0+u4+u20) | 14+2√21 ≈ 23.17 | in F | ∓2 |

   - Each arm is forced by a chain of 2 or 3 copies of G2, so a Haugland-type graph inside ℚ(ζ21) is possible in principle, with a base 2 to 3 times larger. Not built.
2. **Unit rotations with residue in μ21 are useless on their own.** Examples: roots of unity, ω³, and ω⁶ = (13+3i√3)/14 (the √7 spindle). The rotated pieces keep their residues in T, so the 42 colourings extend unless cross edges bring in at least 6 missing classes. An arm of √7 is never forced inside L.
3. **Unit rotations are not automatically excluded.** In Parts' family q = 4 and χ(H4) = 4; here χ(H64) ≥ 6.
   - A unit rotation with residue outside μ21 (ω^k with 3 ∤ k, ω̄, …), applied to a piece with both even and odd directions, brings in the whole missing coset, and the level-1 test is then silent.
   - The same holds after adding the directions of 6 of the 7 mirror triangles.
   - These are necessary conditions only. The 3-adic screen is silent for everything in ℚ(ζ84), and the 7-adic screen is silent once ω is used.

## 7. Status
**Proved or exactly computed.** Exact arithmetic in PARI and Python, and linear algebra over 𝔽2:
- §§1–2.
- All valuations and residues, checked twice.
- The unit vectors of L (qfminim).
- Theorem A with its certificate, and the count of 42 (two independent programs).
- The 6 orbits and the extension of Table 2 to L.
- Corollary B, parts 1 and 2.
- The deep-spindle argument.
- The valuations of ρ and of the cross edge in E.
- χ(H64) ≥ 6 and χ(G27) ≥ 5, by Hoffman with exact character sums.
- The hitting-set numbers (12, 6, 5103, 6 of 7).

**SAT only, no DRAT certificate** (CaDiCaL via pysat):
- Every 4-colouring of Γ_T is linear.
- The sampled 6-sets and the 6 mirror orbits are UNSAT.
- G27 has no 4-colouring (Hoffman already proves this).
- χ(H7) = 4.

**Not established.**
- Haugland's SAT claims for G1, G2 and G3 were used, not re-run, and G1 was not rebuilt from Appendix A. The converse of Corollary B and the spindle chains depend on his G2.
- Whether every 4-colouring of the infinite graph on L is linear. This is Haugland's claim; we are consistent with it but do not prove it.
- Finer 2-adic levels, and residue fields larger than 𝔽64.
- No graph suggested in §6 was built or tested.
- Whether χ(G27) = 5.

## Files
- **field.gp**: H, Table 1, ω, the coordinate field.
- **places.gp**: places, factorisation of (ω), valuations.
- **lattice.gp**: L, its unit vectors, the depths of ω³ and ω⁶.
- **gf.py, classify2.py**: the residue classification.
- **cayley2.py**: the 42 linear colourings; 2912 each for the even-only and odd-only sets.
- **orbits2.py**: orbits, the ζ6 shift, the match with Table 2.
- **certificate.py**: the matrices M and the vectors u_j mod 2.
- **design2.py**: kernel statistics, Hoffman bounds, hitting sets.
- **sat2.py**: the SAT checks.
- **g27.py, h7.log**: the targets at 3 and 7.
- **rho.gp**: E = ℚ(ζ21, √5), ρ and the cross edge.
- **spindles.gp**: 2L-arm spindles inside ℚ(ζ21).
- **crosscheck.gp**: pairwise valuations and residues at 7.

Each script has a matching .log, and the JSON files where produced.
