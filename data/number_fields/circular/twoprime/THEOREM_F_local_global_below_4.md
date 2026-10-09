# A local–global principle for circular colourings below 4

Working note, probe2 folder; continuation of `THEOREM_E_seven_halves.md` (notation and Lemma 12 there).
Lemma and proposition numbers refer to `papers/three-colours/three-colours.tex`. Labels used: Lemma 1,
Lemmas B1–B5, Lemma WA, Lemmas "values" and "radii", Propositions 1–6 (Proposition C1 = `prop:C1`,
Proposition 4 = `prop:levels`, Proposition 5 = `prop:residue`), Corollary B6, Theorems B, C, W⁺.

**Status.** The proof below is complete, with one computer-assisted step: Proposition 7′, an exact
finite computation modulo 325. It is certified by `tree_K2M1.txt` (Farkas branch tree) and
`cert_K2M1.txt` (type points and LP duals), both checked by `verify_window.py`, which shares no code with
the generator. Small exact checks of finite facts about F₄₉ are in `seven_local.py` and
`seven_patterns.py`. Everything else is by hand. It reuses the paper's Lemmas B2–B5 with p = 7 in place of
3, Proposition 4 and Proposition 5. This note is newer than the one for Theorem E and has had less
re-reading. No mathematician has checked either note yet.

---

## Statement

Let F be a number field. Besides conditions (a) and (b) of Theorem B, consider

> **(7)** some place of F above 7 has residue degree 1 (residue field 𝔽₇).

> **Theorem F.** If F² has a homomorphism to a circular clique K_{p/q} with p/q < 4, then (a), (b) or (7)
> holds. Consequently, for every number field F:
> * χ_c(F²) = 2 if (a) holds;
> * χ_c(F²) = 3 if (b) holds and (a) fails;
> * χ_c(F²) = 7/2 if (7) holds and (a), (b) fail;
> * χ_c(F²) ≥ 4 if (a), (b) and (7) all fail.
>
> In particular χ_c(F²) ∈ {2, 3, 7/2} ∪ [4, ∞]. Every value below 4 comes from a locally constant
> colouring of a single completion (Corollary 5).

This answers Question 1 of the paper affirmatively, in its precise form. Examples:
* χ_c(ℚ(√23)²) = 7/2, since 7 splits.
* χ_c(ℚ(√71)²) = 7/2, since 71 ≡ 1 mod 7.
* χ_c(ℚ(√59)²) = 4, since 7 is inert and χ = 4.
* χ_c(ℚ(√47)²) ≥ 4, since (47/7) = −1; this is Question 2 of the note `circular_planes.md`.

In general, for squarefree d ≡ 11 (mod 12), χ_c(ℚ(√d)²) = 7/2 if (d/7) ≠ −1 and χ_c ≥ 4 otherwise.

*Deduction of the list from the first sentence.* Under (a), χ = χ_c = 2. Under (b) without (a), χ = 3
(Theorem B), so χ_c = 3 (Corollary B6). Under (7) without (a) and (b), χ ≥ 4 (Theorem B), so
χ_c ≥ 7/2 (Theorem E). Also χ_c ≤ 7/2 by Proposition C1, since a place above 7 of residue degree 1 has
residue field 𝔽₇. If (a), (b) and (7) all fail, no K_{p/q} with p/q < 4 receives a homomorphism, so
χ_c ≥ 4. If i ∈ F, then (a), (b) and (7) all fail: every residue field above 7 contains √−1, hence 𝔽₄₉;
and χ_c = χ = ∞.

---

## 1. Types at 7

Let 𝔽₄₉ = ℤ[i]/7. Write μ₈ for the elements of norm 1 in 𝔽₄₉; 7 ≡ 3 mod 4, so these are the image of
the rational rotations. Define

    A′₇ = { e ∈ 𝔽₄₉ : Re(ē z) ∈ {2, 3, 4, 5} for every z ∈ μ₈ }   (Re(a+bi) = a),
    E₇ = { ε ∈ (1/7)ℤ[i] : 7ε mod 7 ∈ A′₇ },     𝒜₇ = { e/7 : e ∈ ℤ₄₉, e mod 7 ∈ A′₇ } ⊂ ℚ₄₉.

**Facts** (exact enumeration of 𝔽₄₉, `seven_local.py` → `seven_local.txt`):
* (F1) A′₇ = {2+3i, 2+4i, 3+2i, 3+5i, 4+2i, 4+5i, 5+3i, 5+4i}, a single μ₈-orbit with 8 elements;
  0 ∉ A′₇.
* (F2) A′₇ contains no coset of a non-zero additive subgroup of 𝔽₄₉.
* (F3) ℓ(a+bi) = a maps A′₇ into {2, 3, 4, 5}.

For ε ∈ E₇ the 7-torsion character a ↦ Re(ε̄ ã) (ã = a mod 7) of ℤ[1/5][i] takes values in {2,…,5}/7
at every rational rotation; it is the character of the 7-adic colourings of Proposition C1. Then
6E_c, 3E_q, 7E₇ ⊂ ℤ[i], so 42E′ ⊆ ℤ[i] for E′ := E_c ∪ E_q ∪ E₇.

**Digit patterns** (`seven_patterns.py` → `seven_patterns.txt`). In the digit language of Lemma 12, the
eight 7-adic characters have g-values in {±1/14, ±3/14}² (so they lie in B_s when θ ≤ 2/7). Their digit
words are the 8 shifts of a single word S of period 8 with S_{j+2} = −iS_j: non-zero digits
u(2+i)/5 alternate with zeros. Every 3 consecutive digits of S contain a zero and a non-zero digit, and
determine the shift mod 8. Type c has digits 000…, and type q has all digits non-zero. So any 3
consecutive digits determine the family (c, q or 7) and, within q or 7, the phase.

---

## 2. The window (2, 1) (computer-assisted)

Let N₀ = 5²·13 = 325 and W = G(2, 1) = {iᵃρʲσˡ : |j| ≤ 2, |l| ≤ 1}. The type points modulo N₀ are:
* N₀h (family c);
* (N₀/3)(α+βi), α, β ∈ {1, 2} (family q);
* the 8 points N₀(a+bi)/7 with (a+bi) ≡ 5e (mod 7), e ∈ A′₇ (family 7; 5 ≡ 325⁻¹ mod 7).

> **Proposition 7′.** Let 1/4 < θ ≤ 1/3 and c ∈ S^θ(2, 1). Then c has, for every γ ∈ W, the strip
> index of T + 325m for some type point T (of family c, q or 7) and some m ∈ ℤ[i].

*Proof (computer-assisted).* Put r₀ = 249/1000 < 1/4. Then S^θ(2,1) ⊆ S^{r₀}(2,1), and strip indices
agree for θ and r₀. The proof has two certified parts.

1. `tree_K2M1.txt`, checked by `verify_window.py 2 1`. Up to translation by 325ℤ[i], every
   c ∈ S^{r₀}(2,1) has one of 29 index vectors over the 30 functionals of the 15 rotations ρʲσˡ. The
   certificate is a branch tree over all 325² cells in which every closed branch carries an exact Farkas
   certificate: 550 474 branch nodes, 199 920 closed branches, 29 leaves (12.6 MB; 50 s to check).
2. `cert_K2M1.txt`, same checker. Of the 29 vectors:
   * 5 are those of N₀h and of the four q-points;
   * 8 are those of the 7-torsion points listed above, each with all margins ≥ 2/7;
   * each of the remaining 16 has a dual certificate: λ ≥ 0, Σλ = 1, Σλ_iμ_i ≡ 1/4 identically on three
     margin functions.

   A point of S^θ(2,1) has all margins ≥ θ > 1/4, so it cannot have one of these 16 vectors.

Cross-checks:
* the lift computation (`kappa2d_seven.py`, file `seven_r249_21.txt`);
* the run at r = 1/4 + 10⁻⁷, which leaves only the 13 index vectors of families c, q and 7 (file
  `kappa2d_r25000001_21.txt`);
* the direct enumeration `indepN.py 249/1000 325` (file `indep_325_r249.txt`), which shares no code with
  the lift. It finds the 60 = 4·15 rotations of denominator dividing 325, hence exactly the window, and
  29 components: 1 c, 4 q and 24 others.

As a check of the conclusion of Lemma 12′ (not used in the proof), the windows (2,1), (3,1) and (4,1) at
θ = 0.26 have exactly 13 components, all of families c, q and 7 (`seven_r26_41.txt`). ∎

*Remark.* The window (1,1) does not suffice below 2/7. At r = 0.2501 it keeps extra components with
κ = 11/42 and 25/96 (`kappa2d_r2501_11.txt`). The window (1,2) is worse (`kappa2d_r2501_12.txt`); the
window (2,1) removes them all.

---

## 3. From the window to the 5-probe at every level

> **Lemma 12′.** Let 1/4 < θ ≤ 1/3, k ≥ 2 with k ≡ 0 (mod 6) (so N = 5^k ≡ 1 mod 7), and let C ∈ ℂ
> satisfy Re(C̄γ) ∈ [θ, 1−θ] + ℤ for every γ ∈ G(k, 1). Then C ≡ Nε + x (mod Nℤ[i]) with one of:
> * ε = h and x ∈ P_k;
> * ε ∈ E_q and x ∈ Y_k(ε);
> * ε ∈ E₇ and x ∈ X_k(ε) := {x : ζ_j(ε) + x̄ρʲ ∈ B_s for |j| ≤ k}, where h + ζ_j(ε) ≡ (Nε)‾ρʲ (mod ℤ[i]).
>
> In all cases |x| ≤ R′ := (1/4 + 3/14)√2 < 0.66.

*Proof.* As for Lemma 12, with these changes.

(a) Define g_j, n_j and the digits ν_j exactly as there.

(b) For −k+2 ≤ j ≤ k−2, the point c = ρ^{−j}C lies in S^θ(2,1), because ρʲW ⊆ G(k,1). By
Proposition 7′, c has the strip indices of a type point T + 325m at ±ρᵗ, ±iρᵗ (|t| ≤ 2). So the four
digits ν_{j−1}, ν_j, ν_{j+1}, ν_{j+2} of C equal the digits of T at t = −1, 0, 1, 2. These are:
* family c: 0000, because 325ρᵗ is a Gaussian integer of odd norm;
* family q: (ν, −iν, −ν, iν) with ν ≠ 0, as in Lemma 12(b), now with 325 ≡ 1 and ρ ≡ −i (mod 3);
* family 7: four consecutive digits of a shift of S. Indeed, T is a 7-torsion point all of whose margins
  are ≥ 2/7. Its values at the rotations iᵃρᵗ, |t| ≤ 2, which cover μ₈ mod 7 (ρ ≡ 2+5i has order 8 and
  ρ² ≡ −i), lie in {2,…,5}/7. So its class is in 325⁻¹E₇ and its pattern is a shift of S.

(c) The windows j = −k+2, …, k−2 cover ν_{−k+1}, …, ν_k, and consecutive windows share three digits.
By the digit patterns of §1, the shared digits fix the family and the phase. So the whole word is all
zeros, a q-pattern, or 2k consecutive digits of a shift of S.

(d) Identification. The first two cases are as in Lemma 12(e). In the third case, let ε ∈ E₇ be the class
whose 7-adic character has this digit pattern; the 8 classes give the 8 shifts. Because N ≡ 1 (mod 7),
the point Nε has the values of that character at the ρʲ, so g_j(Nε) = ζ_j(ε) with the same digits. Then
g_j − ζ_j = ρ(g_{j−1} − ζ_{j−1}). Let x be the number with x̄ = g₀ − ζ₀; then g_j = ζ_j + x̄ρʲ, so
x ∈ X_k(ε) and C̄ρʲ ≡ (Nε + x)‾ρʲ (mod ℤ[i]) for |j| ≤ k, whence C − Nε − x ∈ Nℤ[i] as in Lemma 12.

Bounds: |x| ≤ s√10/3 (c), √2δ (q), |g₀| + |ζ₀| ≤ (s + 3/14)√2 (7). Each is less than R′ because
s < 1/4. ∎

---

## 4. Three primes

> **Claim (Lemma B1 with E′).** Let 1/4 < θ ≤ 1/3 and suppose F² has a homomorphism to K_{p/q} with
> q/p = θ. Then for every finite V ⊂ T there is β ∈ L with Tr(βv) ∈ E′ for every v ∈ V.

*Proof.* As the claim in the proof of Theorem E. Choose k ≡ 0 (mod 6) with
N = 5^k > 42R′·max_v(|D_v| + Σ_j|μ_{v,j}|) and U = G(k,1)V. Apply Theorem W⁺ (2q ≤ p < 4q) and
Lemma 12′, then Lemma 1 with 42 in place of 6 (its proof uses only 42E′ ⊆ ℤ[i]). ∎

> **Lemma B2′.** Under the hypothesis of the claim, one of the following holds:
> * (c) as in Lemma B2;
> * (q) as in Lemma B2;
> * (7′) there is β₇ ∈ L₇ with Tr(β₇z) ∈ 𝒜₇ for every z ∈ T₇ = ∏_{u|7} T(F_u).

*Proof.* Follow the proof of Lemma B2 with a third prime. The set
K₇ = {β ∈ L₇ : Tr(βb_j) ∈ (1/7)ℤ₄₉ ∀j} is compact. On {c, q, 7}^T × K₂ × K₃ × K₇ the conditions are:
* τ(v) = c: v_{1+i}(Tr β₂v) = −1, Tr β₃v ∈ O₃, Tr β₇v ∈ ℤ₄₉;
* τ(v) = q: Tr β₂v ∈ O₂, Tr β₃v ∈ 𝒜, Tr β₇v ∈ ℤ₄₉;
* τ(v) = 7: Tr β₂v ∈ O₂, Tr β₃v ∈ O₃, Tr β₇v ∈ 𝒜₇.

They are clopen and are satisfied by the global β of the claim (E_c, E_q and E₇ satisfy them locally).
By compactness some (τ, β₂, β₃, β₇) satisfies them for all v ∈ T.

Let a, a′ (at 2), b, b′ (at 3) and c, c′ (at 7) be the local predicates. Within each pair the two are
mutually exclusive: for c, c′ because 0 ∉ A′₇. The clopen set

    Z = {(a∧b∧c) ∨ (a′∧b′∧c) ∨ (a′∧b∧c′)} ⊆ T₂ × T₃ × T₇

contains the dense diagonal image of T (Lemma WA), so Z is everything. Now three cases:
* If a(z₂⁰) holds for some z₂⁰, then b and c hold everywhere. So b′ and c′ hold nowhere, and a holds
  everywhere: case (c).
* Otherwise a′ holds everywhere. If b′(z₃⁰) holds for some z₃⁰, then c holds everywhere, so b′ holds
  everywhere: case (q).
* Otherwise b holds everywhere and c′ holds everywhere: case (7′). ∎

Cases (c) and (q) give (a) and (b) exactly as in the paper (Section 4 with [TC]; Lemmas B3–B5).

---

## 5. The places above 7

For u | 7 put g_u(z) = Tr_{L_u/ℚ₄₉}(β_u z), where β₇ = (β_u). Then Tr(β₇y) = Σ_u g_u(y_u).

*Even residue degree.* Then i ∈ F_u and T(F_u) = {(s, s⁻¹)}; exactly as in "Even f" (Section 5, with
7^m in place of 3^m), g_u = 0.

*Odd residue degree f, q = 7^f.* L_u = F_u(i) = F_uℚ₄₉ is unramified of degree 2 over F_u, and σ acts on
ℚ₄₉ as Frobenius. T(F_u) ⊂ O_{L_u}^×; it contains the Teichmüller lifts of μ_{q+1}, in particular μ₈
(8 | q+1 since f is odd), and reduces onto μ_{q+1}. Also g_u(ζz) = ζg_u(z) for ζ ∈ μ₈ ⊂ ℚ₄₉.

> **Lemma B3′ (one place).** In case (7′) there are a place u₀ | 7 of odd residue degree and
> β ∈ L_{u₀} with Tr_{L_{u₀}/ℚ₄₉}(βz) ∈ 𝒜₇ for all z ∈ T(F_{u₀}).

*Proof.* As Lemma B3. Varying one component gives g_u(z) − g_u(z′) ∈ 𝒜₇ − 𝒜₇ ⊆ (1/7)ℤ₄₉. With z′ = 1
and z = i, and since i − 1 is a unit at 7, we get g_u(T(F_u)) ⊆ (1/7)ℤ₄₉. Let Z_u ⊆ 𝔽₄₉ be the set of
residues of 7g_u; it is μ₈-stable, and Σ_u z_u ∈ A′₇ for every choice z_u ∈ Z_u.

Suppose u ≠ u′ have non-zero w ∈ Z_u and w′ ∈ Z_{u′}. Fix the other components and put
c = ζ′w′ + Σz_{u″}. The 8 distinct elements c + ζw (ζ ∈ μ₈) lie in A′₇, which has 8 elements, so they
are all of A′₇. Summing gives 8c = c = ΣA′₇ (characteristic 7; Σ_{ζ∈μ₈} ζ = 0). This holds for every
ζ′ ∈ μ₈, so w′ = 0: a contradiction. Since 0 ∉ A′₇, exactly one place u₀ has Z_{u₀} ≠ {0}. ∎

> **Lemma B4′ (level one).** Let u | 7 have odd residue degree, and let β ∈ L_u with
> Tr_{L_u/ℚ₄₉}(βz) ∈ 𝒜₇ for all z ∈ T(F_u). Then λ(x) = 7 Tr(βx) mod 7 is defined on O_{L_u} and vanishes
> on 𝔪, and the induced 𝔽₄₉-linear λ₁ : 𝔽_{q²} → 𝔽₄₉ satisfies λ₁(μ_{q+1}) ⊆ A′₇.

*Proof.* The proof of Lemma B4 holds verbatim with 3, 𝔽₉, A′ replaced by 7, 𝔽₄₉, A′₇. It uses that p
is odd; that the order of 7 mod q+1 is 2f; and fact (F2) in place of "A′ contains no coset of a non-zero
subgroup of 𝔽₉". ∎

> **Lemma B5′ (residue fields).** If f is odd and an 𝔽₄₉-linear λ₁ : 𝔽_{q²} → 𝔽₄₉ satisfies
> λ₁(μ_{q+1}) ⊆ A′₇, then f = 1.

*Proof.* By (F3), ℓ∘λ₁ : 𝔽_{q²} → 𝔽₇ is additive with values in {2,3,4,5} on μ_{q+1}, so
κ₁ ≥ 2/7 > 0 for the corresponding field. Proposition 5 gives κ₁ = 0 for p = 7 and f ≥ 3 (Weil's bound:
every additive map has a zero on μ_{q+1}). Hence f = 1. ∎

So case (7′) gives a place u₀ | 7 of residue degree 1, which is (7).

---

## 6. Proof of Theorem F

Let F² → K_{p/q} with p/q < 4. If i ∈ F this is impossible (χ = ∞). If p/q ≤ 3, then χ ≤ 3 and (a) or
(b) holds by Theorem B. Otherwise θ = q/p ∈ (1/4, 1/3), and Theorem W⁺ applies. By the claim of §4 and
Lemma B2′ we are in case (c), (q) or (7′), which give (a), (b) or (7) by §4 and §5. ∎

---

## 7. Files

| file | role |
|---|---|
| `certify_window.py` → `cert_K2M1.txt`, `tree_K2M1.txt` | certificates for Proposition 7′ (r₀ = 249/1000, N₀ = 325) |
| `verify_window.py` | independent checker (50 s); rejects corrupted copies (changed multiplier, changed 7-point, dropped branch) |
| `kappa2d_seven.py` → `seven_r249_21.txt` | lift computation with labels c/q/7 and exact κ |
| `kappa2d.py` → `kappa2d_r2501_11.txt`, `kappa2d_r2501_12.txt`, `kappa2d_r2501_21.txt`, `kappa2d_r25000001_21.txt` | windows (1,1), (1,2), (2,1) below 2/7 |
| `seven_local.py` → `seven_local.txt` | facts (F1)–(F3) |
| `seven_patterns.py` → `seven_patterns.txt` | digit patterns of the 7-adic characters |
| `indepN.py` → `indep_325_r249.txt` | direct enumeration at modulus 325 (referee geometry) |
