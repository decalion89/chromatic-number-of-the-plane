# The gap (3, 7/2): χ_c(F²) ≥ 7/2 whenever χ(F²) ≥ 4

Working note, probe2 folder. References to lemma and proposition numbers are to
`papers/three-colours/three-colours.tex` (draft of 4 October 2026). The labels used there are
Lemma 1 (exact relations), Lemma B1 to B5, Lemma WA, Lemma "values" and Lemma "radii" (Section 9,
`lem:values`, `lem:radii`), Propositions 1 to 6, Theorems B, C, D and W⁺. New results here are numbered
Proposition 7, Lemma 12 and Theorem E.

**Status.** The proof below is complete, with one computer-assisted step: Proposition 7, an exact
finite computation modulo 65. It is certified by an exact branch tree with Farkas certificates
(`tree_r7_25.txt`) and exact LP dual certificates (`prop7_certificates.txt`). Two short checkers,
`verify_tree.py` and `verify_prop7.py`, accept the certificates, share no code with the generators, and
reject corrupted copies. Two further exact enumerations written separately (`probe2d.py`, `indep2d.py`)
agree. The rest of the proof is by hand. It uses, unchanged, Theorem W⁺, Lemma 1, the proofs of
Lemmas B1 to B5 and Theorem B, and Lemmas "values" and "radii" of the paper. No mathematician has
checked this note yet.

---

## Statement

> **Theorem E.** Let F be a number field. If F² has a homomorphism to a circular clique K_{p/q} with
> p/q < 7/2, then χ(F²) ≤ 3. Consequently χ_c(F²) ≥ 7/2 whenever χ(F²) ≥ 4. For every number
> field F we therefore have χ_c(F²) ∈ {2, 3} ∪ [7/2, ∞]. By Theorem C (χ_c(ℚ(√11)²) = 7/2) the
> bound 7/2 is attained, so the gap (3, 7/2) is sharp.

This supersedes Theorem D, whose bounds were 56/17 and, with a computer-assisted step, 3.3315. It also
answers the part of Question 1 of the paper that lies below 7/2. It reproves the lower bounds of
Theorem C, χ_c(ℚ(√11)²) ≥ 7/2 and χ_c(ℚ(√35)²) ≥ 7/2, from χ ≥ 4 alone, without the branch-and-bound
certificates of Section 10. The same method with the window (2,1) gives all of Question 1 below 4
(Theorem F, `THEOREM_F_local_global_below_4.md`): χ_c(F²) ∈ {2, 3, 7/2} ∪ [4, ∞], with 7/2 exactly
when (a) and (b) fail and some place above 7 has residue degree 1.

**The idea.** The probe of Section 9 uses only the rotations G_N = {iᵃρʲ}, which have 5-power
denominators. Below 3/10 its characters include "5-adic" ones: glued type-c characters, and torsion
characters of order 41, 76, … (see the section on the one-prime probe below). These have no bounded
type, so the exact-relations lemma has nothing to work with. The plane over F, however, contains every
rational rotation. Adding a single rotation with a different prime, σ = (3+2i)/(3−2i) = (5+12i)/13,
removes all of them. Already in the window |j|, |l| ≤ 1 (the 36 rotations iᵃρʲσˡ), for θ > 2/7 every
character is of type c or q at every position of the 5-probe. The only characters lost at θ = 2/7 are
the 7-adic ones of Proposition C1, which is why 7/2 is the exact limit of the method.

---

## 1. Notation

As in the paper: L = F(i) with i ∉ F, T = {z ∈ L : z z̄ = 1}, Tr = Tr_{L/ℚ(i)}, ρ = (3+4i)/5 = (2+i)/(2−i),
h = (1+i)/2, E_c = h + ℤ[i], E_q = {(α+βi)/3 : α, β ∈ {1,2}} + ℤ[i], E = E_c ∪ E_q. For
1/4 < θ ≤ 1/3 put s = 1/2 − θ, δ = 1/3 − θ, B_s = {a+bi : |a|, |b| ≤ s}. For N = 5^k,
P_k = {x : x̄ρʲ ∈ B_s for |j| ≤ k}, and for ε ∈ E_q the elements η_j(ε) ∈ {±1±i} and the sets Y_k(ε)
are as in Section 9 (Lemma "values"(iii)).

New notation: σ = (3+2i)/(3−2i) = (5+12i)/13 and, for k, m ≥ 0,

    G(k, m) = { iᵃ ρʲ σˡ : a ∈ ℤ, |j| ≤ k, |l| ≤ m },
    S^θ(k, m) = { c ∈ ℂ : Re(c̄ γ) ∈ [θ, 1−θ] + ℤ for every γ ∈ G(k, m) }.

G(k, m) consists of rational rotations. Since i ∉ F, the automorphism of L/F restricts to complex
conjugation on ℚ(i), so G(k, m) ⊂ T. S^θ(k, 0) is the probe S_N^θ of the paper. For c ∈ S^θ(k, m) and
γ ∈ G(k, m), the **strip index** n_γ(c) ∈ ℤ is the integer with
Re(c̄γ) ∈ [n_γ(c) + θ, n_γ(c) + 1 − θ]. It is unique because θ > 0. Since 65γ ∈ ℤ[i] for
γ ∈ G(1,1), the set S^θ(1,1) is invariant under translation by 65ℤ[i]. The five **type points** at
modulus 65 are T_c = 65h and T_{αβ} = (65/3)(α+βi), α, β ∈ {1, 2}.

---

## 2. The two-prime window (computer-assisted)

> **Proposition 7.** Let 2/7 < θ ≤ 1/3 and c ∈ S^θ(1,1). Then there are a type point
> T ∈ {T_c, T_11, T_12, T_21, T_22} and m ∈ ℤ[i] with n_γ(c) = n_γ(T + 65m) for every γ ∈ G(1,1).

*Proof (computer-assisted).* Put r₀ = 7/25 < 2/7. Since [n+θ, n+1−θ] ⊂ [n+r₀, n+1−r₀], we have
S^θ(1,1) ⊆ S^{r₀}(1,1), and a point of S^θ(1,1) has the same strip indices for θ and for r₀. Since
Re(c̄·iγ) = −Im(c̄γ) and [r, 1−r] + ℤ is symmetric, only the 18 functionals c ↦ Re(c̄γ), Im(c̄γ) for the
nine rotations γ = ρʲσˡ, |j|, |l| ≤ 1, need to be considered.

**(i) Completeness** (`tree_r7_25.txt`, checked by `verify_tree.py`). Let c ∈ S^{r₀}(1,1). After a
translation by 65ℤ[i], c lies in a cell [a+r₀, a+1−r₀] × [b+r₀, b+1−r₀] with 0 ≤ a, b ≤ 64; this cell
fixes the indices at γ = 1. For each of the other 16 functionals, in a fixed order, the possible strip
indices are the integers n for which [n+r₀, n+1−r₀] meets the range of the functional on the cell; the
checker computes these from the four corners. The tree lists, for every cell and every branch, each
possible index either as a branch to follow or as closed. A closed branch carries an exact Farkas
certificate: three of the valid inequalities on that branch (cell sides, and strip sides of the
functionals fixed so far) with rational multipliers λ ≥ 0 whose combination reads 0 ≤ (a negative
number). So c follows open branches only, and these end in 13 leaves. Every c ∈ S^{r₀}(1,1) therefore
has, up to translation by 65ℤ[i], one of 13 explicit index vectors I₁, …, I₁₃. The tree has 4 225 cells,
15 620 branch nodes, 4 364 closed branches and 13 leaves (288 kB). The checker runs in 2 s.

**(ii) The 13 index vectors** (`prop7_certificates.txt`, checked by `verify_prop7.py`). Five of the
I_j are the index vectors of the five type points. Each of the other eight comes with three of its 36
margin functions

    μ⁺_γ(c) = Re(c̄γ) − n_γ,    μ⁻_γ(c) = n_γ + 1 − Re(c̄γ),

and rationals λ₁, λ₂, λ₃ ≥ 0 with Σλ_i = 1 such that Σ λ_i μ_i ≡ 2/7 identically (the linear parts
cancel exactly). A point of S^θ(1,1) has every margin ≥ θ > 2/7. So no point of S^θ(1,1) has one of
these eight index vectors (up to translation), and Proposition 7 follows. ∎

*Remarks.* (1) The eight discarded index vectors contain the 7-torsion points 65(a+bi)/7 with
(a, b) ∈ {(1,2), (1,5), (2,1), (2,6), (5,1), (5,6), (6,2), (6,5)}, at which the smallest margin is
exactly 2/7. They are the characters of the 7-adic colourings of Proposition C1; for instance
λ(a+bi) = 2a+3b corresponds to the point 65(1+5i)/7. So 2/7 is optimal for this window, and for every
window, since these are characters of all of ℚ(i).
(2) Cross-checks: `probe2d.py` (lifting, in both orders 5 then 13 and 13 then 5) and `indep2d.py`
(direct over all 65² cells, rotations obtained by solving x² + y² = 65², geometry in H-representation
written separately) both find 13 components at r = 7/25 and 5 components at r = 2857143/10⁷ (just above
2/7). `kappa2d.py` computes the exact κ of the eight extra components by certified LP. A floating-point
Monte Carlo run (2·10⁷ points) finds no point of S^θ(1,1) far from a type point.
(3) The prime 13 can be replaced by 17 (σ' = (15+8i)/17, N = 85; same result just above 2/7), but the
prime 5 cannot be dropped: with the primes 13 and 17 the window (1,1) keeps 12 extra components just
above 2/7. The mixed rotations matter. The "cross" {1, ρ^{±1}, σ^{±1}} alone keeps 228 extra components,
some with κ = 55/144 > 1/3. Removing one of the nine rotation classes can still suffice
(`subwin.py`).

---

## 3. From the window to the probe at every level

> **Lemma 12.** Let 2/7 < θ ≤ 1/3, k ≥ 1, N = 5^k, and let C ∈ ℂ satisfy Re(C̄γ) ∈ [θ, 1−θ] + ℤ for
> every γ ∈ G(k, 1). Then C ∈ Nh + P_k + Nℤ[i], or C ∈ Nε + Y_k(ε) + Nℤ[i] for some ε ∈ E_q. In
> particular C = Nε' + x with ε' ∈ E and |x| ≤ R := max(s√10/3, √2·δ) (< 0.23).

*Proof.* **(a) Values and digits.** For |j| ≤ k the conditions at γ = ρʲ and γ = −iρʲ say that Re
and Im of C̄ρʲ lie in [θ, 1−θ] + ℤ (as Re(C̄·(−iρʲ)) = Im(C̄ρʲ)). So

    C̄ρʲ = h + g_j + n_j,   g_j ∈ B_s,   n_j ∈ ℤ[i],

uniquely since s < 1/2, and n_j = n_{ρʲ}(C) + i·n_{−iρʲ}(C). For −k < j ≤ k put ν_j = g_j − ρ g_{j−1}.
Since C̄ρʲ = ρ·C̄ρ^{j−1},

    ν_j = (ρ − 1)h + ρ n_{j−1} − n_j,

so ν_j depends only on the strip indices of C at ρ^{j−1}, −iρ^{j−1}, ρʲ and −iρʲ.

**(b) The type points.** Let m ∈ ℤ[i] and t ∈ {−1, 0, 1}. For T = 65h + 65m: 65ρᵗ is a Gaussian
integer of odd norm, so T̄ρᵗ ≡ h (mod ℤ[i]) (the argument of Lemma "values"(ii)), g_t(T) = 0, and
ν₀(T) = ν₁(T) = 0. For T = (65/3)(α+βi) + 65m: T̄ρᵗ ≡ (α−βi)u_t/3 (mod ℤ[i]) with u_t ≡ 65ρᵗ (mod 3).
Since 65 ≡ −1 and ρ ≡ −i (mod 3), u_t = −(−i)ᵗ ∈ μ₄. Both coordinates of (α−βi)u_t/3 lie in
{1/3, 2/3} mod 1, so g_t(T) = η_t/6 with η_t ∈ {±1±i}, and η_{t+1} = −iη_t (because −ih ≡ h). Hence
ν_t(T) = η_t/6 − ρη_{t−1}/6 = −η_{t−1}(1+3i)/10 ≠ 0, and ν₁(T) = −iν₀(T).

**(c) Windows.** Fix j with −k < j < k and put c = ρ^{−j}C, so that c̄ = ρʲC̄. For γ ∈ G(1,1) we have
ρʲγ ∈ G(k,1), so Re(c̄γ) = Re(C̄ρʲγ) ∈ [θ, 1−θ] + ℤ and c ∈ S^θ(1,1). By Proposition 7, c has the
strip indices of some T + 65m. Since c̄ρᵗ = C̄ρ^{j+t}, the strip indices of c at ρᵗ, −iρᵗ
(t ∈ {−1,0,1}) are those of C at ρ^{j+t}, −iρ^{j+t}. By (a), ν_{j+t}(C) = ν_t(c) = ν_t(T + 65m) for
t ∈ {0, 1}. By (b), (ν_j, ν_{j+1}) is either (0, 0) or (ν, −iν) with ν = −η(1+3i)/10 for some
η ∈ {±1±i}.

**(d) Pure words.** The pairs (ν_j, ν_{j+1}), −k < j < k, cover ν_{−k+1}, …, ν_k, and consecutive
pairs share a digit. So either every ν_j is 0, or every ν_j is non-zero and ν_{j+1} = −iν_j throughout.

**(e) Identification.** *All digits 0.* Then g_j = ρʲg₀. Let x be the number with x̄ = g₀. Then
x̄ρʲ = g_j ∈ B_s, so x ∈ P_k, and C̄ρʲ ≡ h + x̄ρʲ ≡ (Nh + x)‾ρʲ (mod ℤ[i]) by Lemma "values"(ii).
Multiplying by units gives Re((C − Nh − x)‾γ) ∈ ℤ for every γ ∈ G_N. The ℤ-span of G_N is (1/N)ℤ[i]:
it is a ℤ[i]-module, and it contains ρ^{±k}, with Nρ^{±k} = (2±i)^{2k} coprime in ℤ[i]. Hence
C − Nh − x ∈ Nℤ[i] (take a = 1/N and i/N in Re(w̄a) ∈ ℤ).

*The q-pattern.* By (c), every ν_j = −η'(1+3i)/10 with η' ∈ {±1±i}. Put η_{j−1} = −10ν_j/(1+3i) for
−k < j ≤ k, and η_k = −iη_{k−1}. Then η_j ∈ {±1±i} and η_j = −iη_{j−1}, and

    η_j/6 − ρη_{j−1}/6 = −η_{j−1}(i+ρ)/6 = −η_{j−1}(1+3i)/10 = ν_j.

So g_j − η_j/6 = ρ(g_{j−1} − η_{j−1}/6). Let y be the number with ȳ = g₀ − η₀/6; then
g_j = η_j/6 + ȳρʲ for |j| ≤ k. Choose ε ∈ E_q with η₀(ε) = η₀; the four classes of E_q give the four
values (Lemma "values"(iii)). Since also η_{j+1}(ε) = −iη_j(ε), we get η_j(ε) = η_j for |j| ≤ k. Then
η_j(ε)/6 + ȳρʲ = g_j ∈ B_s, so y ∈ Y_k(ε), and C̄ρʲ ≡ h + η_j(ε)/6 + ȳρʲ ≡ (Nε + y)‾ρʲ (mod ℤ[i]).
As before, C − Nε − y ∈ Nℤ[i].

The bounds |x| ≤ s√10/3 and |y| ≤ √2δ are Lemma "radii". Both proofs hold for every 1/4 < θ ≤ 1/3.
Finally, s < 3/14 gives R < 0.23. ∎

---

## 4. Proof of Theorem E

Let F² have a homomorphism to K_{p/q} with p/q < 7/2. If i ∈ F, then χ(F²) = ∞ (Davies), which is
impossible because K_{p/q} has p vertices; so i ∉ F. If p/q ≤ 3, then χ(F²) ≤ χ(K_{p/q}) ≤ 3.
Otherwise θ := q/p ∈ (2/7, 1/3) and 2q ≤ p < 4q, so Theorem W⁺ applies.

**Claim.** For every finite V ⊂ T there is β ∈ L with Tr(βv) ∈ E for every v ∈ V. This is the
conclusion of Lemma B1.

*Proof of the claim.* Enlarge V so that it contains a ℚ(i)-basis b₁, …, b_n of L inside T (exists by
the proof of Lemma B1). For v ∈ V write D_v v = Σ_j μ_{v,j} b_j with D_v ∈ ℤ∖{0} and μ_{v,j} ∈ ℤ[i].
Choose k ≥ 1 with N = 5^k > 6R·max_v(|D_v| + Σ_j |μ_{v,j}|).

Let U = G(k, 1)·V. It is a finite symmetric set of unit vectors (G(k,1) ⊂ T and T is a group), so
Cay(ℤU, U) is a subgraph of the plane Cay(L, T) and has a homomorphism to K_{p/q}. By Theorem W⁺ there
is a character ξ of ℤU with ξ(U) ⊆ [θ, 1−θ] + ℤ.

Exactly as in the proof of Lemma B1, there is a ℚ(i)-linear φ : L → ℂ with ξ(w) = Re φ(w) mod 1 for w in
a lattice Λ ⊇ U. Then ξ(γv) = Re(γφ(v)) for v ∈ V and γ ∈ G(k,1). Put C_v = φ(v)‾. Then
Re(C̄_v γ) = Re(γφ(v)) ∈ [θ, 1−θ] + ℤ for every γ ∈ G(k,1).

By Lemma 12, C_v = Nε'_v + x'_v with ε'_v ∈ E and |x'_v| ≤ R. So φ(v) = Nε_v + x_v with ε_v = ε̄'_v ∈ E
(E is stable under conjugation: h̄ = h − i and (α−βi)/3 = (α+(3−β)i)/3 − i) and |x_v| ≤ R. Now apply
Lemma 1 to D_v φ(v) − Σ_j μ_{v,j} φ(b_j) = 0, with R in place of r (the proof uses only 6E ⊆ ℤ[i] and
the bound on |x|, as noted in the proof of Theorem D). It gives D_v ε_v = Σ_j μ_{v,j} ε_{b_j}. Let β ∈ L
be the element with Tr(βb_j) = ε_{b_j}. Then D_v Tr(βv) = D_v ε_v, so Tr(βv) = ε_v ∈ E. ∎ (claim)

The proof of Lemma B2 uses the 3-colouring only through the conclusion of Lemma B1, and the rest of the
proof of Theorem B uses only (c) or (q); the paper notes this in the proof of Theorem D. So (a) or (b)
holds, and χ(F²) ≤ 3 by Theorem B.

For the second sentence: if χ(F²) ≥ 4 and χ_c(F²) < 7/2, the definition of χ_c as an infimum gives a
homomorphism to some K_{p/q} with p/q < 7/2, which we have excluded. With Corollary B6 (χ_c = χ when
χ ≤ 3), χ_c(F²) ∈ {2, 3} ∪ [7/2, ∞]. Theorem C shows that 7/2 is attained. ∎

*Remark.* With θ = 1/3, Proposition 7 and Lemma 12 also give the case of Lemma B1 used for Theorem B
(Proposition 1 is then not needed), but Proposition 1 is a short human proof there.

---

## 5. Why the one-prime probe cannot reach 7/2

This section records context; nothing in Sections 2 to 4 depends on it.

**Digit dynamics** (`gdyn.py`). A character of A = ℤ[1/5][i] whose values on G_∞ = {iᵃρʲ} stay in
[θ, 1−θ] is the same thing as a two-sided sequence g_j ∈ B_s with g_{j+1} − ρg_j ∈ Λ₊ = ((2+i)/5)ℤ[i].
For s < 3/14 the "digit" g_{j+1} − ρg_j is one of five: 0, ±(2+i)/5, ±i(2+i)/5. The components of S_N^θ
correspond one-to-one to the digit words of length 2k with non-empty polygon. This formulation
reproduces the component counts of `levels3.py` exactly at every level k ≤ 18 at θ = 0.3001611
(`test_counts.py`). Type c is the zero word; type q is the word with ν_{j+1} = −iν_j ≠ 0.

**At θ = 3/10.** The 5-cores of the words stabilise from level 13 to level 22 (`cores.py`). Besides c
and q, only words with one "glue" survive: a 2-digit q-burst between two zero-runs of radius 1/5. The
glue can sit at every position. The sequence

    g_j = ρ^{j−t}/5 (j < t),   g_t = −(1+i)/5,   g_j = iρ^{j−t}/5 (j > t)

is such a character (with κ = 3/10) for every t. Its point at level K needs a type with denominator
exactly 5^{K−t+1} (`glue_denominators.py` → `glue_denominators.txt`, exact). So no fixed finite set of
types covers S_N^θ for large N once θ ≤ 3/10. The exact-relations lemma then has no fixed lattice, and
"localisation at 5" has nothing to start from.

**Below 3/10.** At θ = 0.2999 the 14-cores stabilise from level 22 to 40 with at most one glue each:
c→c glues through a 2-digit q-burst, and q→q glues through a run of five zeros (`cores2.py`). But
one-sided words with several glues exist (`forward_glues.py`; at θ = 0.299, words with five glues of
length 60). At θ = 0.295 the number of words grows quickly (7 145 at level 14).

**Torsion types** (`torsion_exact.py` → `torsion_exact_400.txt`). Contrary to
`probe/torsion_260.txt`, there are torsion characters of A with κ > 2/7 on G_∞ other than c and q:
order 41 with κ = 12/41 ≈ 0.2927 (c = (12+12i)/41, where ρ has order 5 mod 41), and orders 76, 82,
123, 152, 164, 228, 246, 287, 304, 328, 369 (up to 400). The largest value found is 12/41. The original
program (`torsion_search_orig_copy.py`, an unmodified copy) reports M = 41 when run with arguments
`45 0.29`, so the recorded output `torsion_260.txt` ("top: []") is inconsistent with its own program and
should be corrected. For the two-prime group ⟨i, ρ, σ⟩ mod d none of these survives: no character of
order ≤ 300 prime to 65 other than c and q has κ > 2/7 (`torsion_twoprime.py` →
`torsion_twoprime_300.txt`); the order-41 character takes the value 0. By Proposition 7 this holds for
every order.

**What σ does (heuristic).** A glue forces the archimedean part to sit at a specific angle relative to
the 5-adic digit. Multiplying by σ^{±1} rotates the archimedean part by about ±67.4° and leaves the
5-adic digit unchanged (σ is a 5-adic unit), so the rows l = ±1 of the window violate the glue
condition. Torsion types with small ρ-orbits, such as order 41, get large ⟨ρ, σ⟩-orbits and reach 0.

---

## 6. Files (all in this folder)

| file | role |
|---|---|
| `probe2d.py` | two-prime probe S^θ(k, m) by exact lifting (5- and 13-lifts), labels C/Q/X |
| `kappa2d.py`, `lpexact.py` | exact κ of every extra component (certified LP: primal point plus dual multipliers) |
| `prop7_certificates.py` → `prop7_certificates.txt` | the 13 index vectors at r = 7/25, type points, dual certificates (κ = 2/7) |
| `verify_prop7.py` | independent checker of `prop7_certificates.txt` (rotations rebuilt from Gaussian integers) |
| `gen_tree.py` → `tree_r7_25.txt` | Farkas-certified branch tree over all 65² cells (completeness) |
| `verify_tree.py` | independent checker of the tree (2 s); rejects corrupted copies |
| `indep2d.py`, `indepN.py`, `hgeom_ref.py` | direct enumeration with the referee's geometry (copied from `probe/indep/hgeom.py`) |
| `subwin.py` | sub-windows (cross, 7 or 8 rotation classes) |
| `mc2d.py` | floating-point Monte Carlo (sanity check only) |
| `kappa2d_seven.py` → `seven_r2858_51.txt` | check of Lemma 12's conclusion (not used in the proof): windows (1,1) to (5,1) at θ = 0.2858 have exactly the 5 main components |
| `gdyn.py`, `test_counts.py`, `runstats.py`, `cores.py`, `cores2.py`, `forward_glues.py` | one-prime probe as digit words (Section 5) |
| `glue_denominators.py` → `glue_denominators.txt` | unbounded type denominators at θ = 3/10 |
| `torsion_exact.py` → `torsion_exact_400.txt`; `torsion_twoprime.py` → `torsion_twoprime_300.txt` | torsion types, one prime and two primes |

To re-check the computer-assisted step (about 2 s each):

    python3 verify_tree.py
    python3 verify_prop7.py

To regenerate the certificates (about 10 s each):

    python3 prop7_certificates.py
    python3 gen_tree.py
