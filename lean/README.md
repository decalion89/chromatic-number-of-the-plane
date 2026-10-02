# Formal proofs in Lean 4

This directory holds proofs in Lean 4, with Mathlib, of:
- the two theorems of the note
  [`papers/planes-4-chromatic/planes-4-chromatic.pdf`](../papers/planes-4-chromatic/planes-4-chromatic.pdf): the planes
  over ℚ(√2, √3) and ℚ(√3, √11) have chromatic number 4. The second is K. G. Fischer's theorem (1994); the note
  gives a short proof of both, and these files check that proof;
- χ(ℚ(√d)²) = 4 for d = 11, 119, 131, 179, 191, 251, 431, 455, 911 and 935, ten of the real quadratic fields of
  [`notes/quadratic_planes.md`](../notes/quadratic_planes.md), one file per field. Each lower bound is the graph
  of `data/quadratic_planes/q{d}.json`, whose SAT certificate the kernel checks. The upper bounds are
  reductions at a place: at 7 when d is a nonzero square modulo 7 (Moorhouse's Lemma 8.2) or d = 7d′ (d = 119,
  455), and at 2 when d ≡ 3 (mod 8) (Fischer's Theorem 10; d = 131, 251);
- the same for the other fourteen fields of the note, d = 23, 35, 59, 71, 95, 155, 239, 263, 359, 443, 599, 611,
  791 and 959, and χ(ℚ(√47)²) ≥ 4, given that the graph's colouring formula `q{d}.cnf` is unsatisfiable. Their
  LRAT proofs (from 0.7 MB to 70 MB) are beyond what `lrat_proof` checks in reasonable time and memory; the
  verified checker cake_lpr checked them (`data/quadratic_planes/cake_lpr_checks.txt`). The kernel checks the
  points, the unit distances and the upper bound as for the first ten; that `q{d}.cnf` is the graph's formula is
  checked by evaluation (`#guard`) when the file is built.

## The theorems

```lean
theorem Q23.chromaticNumber_eq_four : (LocalColouring.unitDistGraph Q23.L).chromaticNumber = 4
theorem Q311.chromaticNumber_eq_four : (LocalColouring.unitDistGraph Q311.L).chromaticNumber = 4
theorem Sqrt11.chromaticNumber_eq_four :
    (LocalColouring.unitDistGraph (QuadraticPlanes.L 11)).chromaticNumber = 4
-- and the same for Sqrt119, Sqrt131, Sqrt179, Sqrt191, Sqrt251, Sqrt431, Sqrt455, Sqrt911 and Sqrt935
theorem Sqrt95.chromaticNumber_eq_four_of_unsatisfiable (h : QuadraticPlanes.Unsatisfiable Sqrt95.formula) :
    (LocalColouring.unitDistGraph (QuadraticPlanes.L 95)).chromaticNumber = 4
-- and the same for the other thirteen fields; for d = 47 the lower bound only:
theorem Sqrt47.not_colorable_three_of_unsatisfiable (h : QuadraticPlanes.Unsatisfiable Sqrt47.formula) :
    ¬(LocalColouring.unitDistGraph (QuadraticPlanes.L 47)).Colorable 3
```

Here `Q23.L`, `Q311.L` and `QuadraticPlanes.L d` are the subfields ℚ(√2, √3), ℚ(√3, √11) and ℚ(√d) of ℝ,
and `unitDistGraph K` joins two points of `K × K` at Euclidean distance 1:

```lean
noncomputable def Q23.L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√2, √3}
noncomputable def Q311.L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√3, √11}
noncomputable def QuadraticPlanes.L (d : ℕ) : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√(d : ℝ)}

def LocalColouring.unitDistGraph (K : IntermediateField ℚ ℝ) : SimpleGraph (K × K) where
  Adj p q := ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1
```

`chromaticNumber` is Mathlib's `SimpleGraph.chromaticNumber`, with values in `ℕ∞`. In the conditional theorems,
`Sqrt{d}.formula` is `QuadraticPlanes.colourCNF E u₀ v₀`, the colouring formula of the graph as a list of
clauses (lists of nonzero integers, as in the DIMACS format; variable `3v + c + 1` says that vertex `v` has
colour `c`), and

```lean
def QuadraticPlanes.Unsatisfiable (F : List (List ℤ)) : Prop := ∀ σ : ℕ → Prop, ∃ c ∈ F, ∀ l ∈ c, ¬LitHolds σ l
```

says that no assignment satisfies every clause. When a file is built, `#guard` checks that
`data/quadratic_planes/q{d}.cnf` is exactly `formula`, and `scripts/verify_quadratic_planes.py --cake-lpr` has
cake_lpr check an LRAT proof that this file is unsatisfiable.

The twenty-seven theorems depend only on Lean's three standard axioms, `propext`, `Classical.choice` and
`Quot.sound`: `axioms.expected` records the output of `#print axioms`, and CI compares them.

## Checking the proofs

With [elan](https://github.com/leanprover/elan) installed:

    cd lean
    lake exe cache get                       # Mathlib's prebuilt files, a few minutes
    lake build                               # every file: about forty minutes, one file at a time
    lake env lean PrintAxioms.lean           # compare with axioms.expected

`lean-toolchain` pins Lean v4.34.1 and `lakefile.toml` pins Mathlib v4.34.1. The GitHub Actions workflow
[`lean.yml`](../.github/workflows/lean.yml) runs these steps on every push to `main` and every pull
request. It then replays each file in Lean's kernel with `leanchecker`, independently of the tactics that
produced the proofs. The field files with an LRAT proof take from half a minute (d = 11, 191) to six minutes
(d = 119) and up to 2.5 GB of memory (d = 179) each; the time grows with the number of clauses of the formula, the
memory with the size of the LRAT proof. The conditional files take about a minute each.

## From the note to the files

| the note | Lean |
|---|---|
| the graph on L² | `LocalColouring.unitDistGraph`, `Q23.L`, `Q311.L` |
| §2, the identity x² + y² = α² − αβ + β² | `LocalColouring.toEisHom` |
| Lemma 3 | `LocalColouring.core_lemma` |
| §3, Theorem 2 (the criterion) | `LocalColouring.colorable_four_of_residue_two`, from `qGraph_colorable` |
| a place above 2 | `LocalColouring.exists_valuationSubring_two_mem_nonunits` (Chevalley's extension theorem) |
| §4, residue field 𝔽₂ for ℚ(√2, √3) | `Q23.piL_eisenstein`, `Q23.two_eq_pi_pow_four_mul`, `Q23.mem_nonunits_or_sub_one_mem` |
| §4, residue field 𝔽₂ for ℚ(√3, √11) | `Q311.two_eq`, `Q311.hensel`, `Q311.mem_nonunits_or_sub_one_mem` |
| §4, upper bounds | `Q23.colorable_four`, `Q311.colorable_four` |
| §4, lower bounds | `Q23.chain23_not_colorable`, `Q311.moser_not_colorable`, `not_colorable_three` |
| Theorem 1 | `Q23.chromaticNumber_eq_four`, `Q311.chromaticNumber_eq_four` |

For the real quadratic fields (`notes/quadratic_planes.md`; `QuadraticPlanes` holds what the fields share):

| the note | Lean |
|---|---|
| §2, the points `[a, b, c, e]` and the identities (1) | `QuadraticPlanes.pt`, `unitPairB`, `adj_pt` |
| the graph `data/quadratic_planes/q{d}.json` | `Sqrt{d}.P`, `Sqrt{d}.E`, `Sqrt{d}.checkEdges_E` (every edge has length 1) |
| §4, the formula `q{d}.cnf` and its proof | `Sqrt{d}.refuted` (Mathlib's `lrat_proof`, from `q{d}.cnf` and `q{d}.lrat`) |
| no 3-colouring | `QuadraticPlanes.clause_facts`, `Sqrt{d}.graph_not_colorable`, `Sqrt{d}.not_colorable_three` |
| §5, Moorhouse's Lemma 8.2 at 7 | `QuadraticPlanes.exists_valuationSubring_seven_mem_nonunits`, `exists_int_sub_mem` (residue field 𝔽₇ when d ≡ s² mod 7, by Hensel's lemma for √d), `exists_int_sub_mem_ramified` (when d = 7d′), `core_lemma_seven`, `sumSqGraph_colorable`, `colorable_four`, `colorable_four_ramified` |
| §5, Fischer's Theorem 10 (d ≡ 3 mod 8) | `QuadraticPlanes.mem_nonunits_or_sub_one_mem_two` (residue field 𝔽₂), `toKHom` (x² + y² = a² − ab + ((d + 1)/4)b²), `core_lemma_two`, `kGraph_colorable`, `colorable_four_two` |
| the theorem for d | `Sqrt{d}.chromaticNumber_eq_four` |
| §4, the formula as a list of clauses, for the other fields | `QuadraticPlanes.colourCNF`, `Unsatisfiable`, `not_colorable_of_unsatisfiable`, `Sqrt{d}.formula` |
| the theorem for the other fields, given the formula's unsatisfiability | `Sqrt{d}.chromaticNumber_eq_four_of_unsatisfiable`; `Sqrt47.not_colorable_three_of_unsatisfiable` |

The proofs follow the notes with one difference. The notes take a prime of the ring of integers above 2 or 7
and read its residue field from the decomposition of the prime. Here the place is a valuation subring with 2
(or 7) in its maximal ideal, given by Chevalley's extension theorem, and its residue field is shown directly:
- for ℚ(√2, √3), with the uniformiser π = (√2 + √6)/2 − 1, a root of the 2-Eisenstein polynomial
  X⁴ + 4X³ + 2X² − 4X − 2;
- for ℚ(√3, √11), where 2 has two places, with π = √3 − 1 and a Hensel-type argument for (1 + √33)/2;
- for ℚ(√d) at 7, with Hensel's lemma for the root of X² − d near ±s when d ≡ s² (mod 7), and by descent on
  the power of 7 in the denominators when d = 7d′; at 2, with π = √d − 1 when d ≡ 3 (mod 4).

The valuation subring and the coset representatives come from Zorn's lemma and `Quotient.out`, so the
colourings exist but cannot be computed.

## Files

| file | contents |
|---|---|
| `LocalColouring.lean` | the graph, the criterion at 2 for fields containing √3, and the tools the fields share |
| `Q23.lean` | ℚ(√2, √3): the place above 2, and the 10-vertex graph of [`data/chain23.json`](../data/chain23.json) |
| `Q311.lean` | ℚ(√3, √11): the place above 2, and the Moser spindle |
| `QuadraticPlanes.lean` | ℚ(√d): the three upper bounds (at 7, at 7 ramified, at 2), and the tools for the lower bounds |
| `ColouringFormula.lean` | the colouring formula as a Lean object (`colourCNF`), the lower bound from its unsatisfiability, a reader for DIMACS files, and a balanced tree of points (`PtTree`) in which the kernel finds a point in logarithmic time |
| `Sqrt{d}.lean` | ℚ(√d), for d = 11, 119, 131, 179, 191, 251, 431, 455, 911, 935: the graph of `data/quadratic_planes/q{d}.json`. It reads `data/quadratic_planes/q{d}.cnf` and the LRAT proof `data/quadratic_planes/q{d}.lrat` (kissat, then `drat-trim -L`) with `include_str`. For the other fifteen graphs (d = 23, 35, 47, 59, 71, 95, 155, 239, 263, 359, 443, 599, 611, 791, 959): the graph, `formula`, and the check that `q{d}.cnf` is `formula` |
| `PrintAxioms.lean`, `axioms.expected` | the axioms of the twenty-seven theorems |
| `tools/q23_coefficients.py`, `tools/q311_coefficients.py` | sympy scripts that produce the coefficients of the `linear_combination` steps and the edge lists |
| `tools/field_lean.py` | writes `Sqrt{d}.lean` from the data (`--check` compares, and checks the colouring of 𝔽₇² in `QuadraticPlanes.lean` against `finite_planes.json`; the tests run it) |

## What is trusted

Lean's kernel; Mathlib's definitions of `IntermediateField.adjoin`, `Real.sqrt` (written `√`) and
`SimpleGraph.chromaticNumber`; and the definitions of `Q23.L`, `Q311.L`, `QuadraticPlanes.L` and
`unitDistGraph` quoted above. For the real quadratic fields, Mathlib's `lrat_proof` turns each LRAT proof into a
proof term that the kernel checks; nothing runs as native code, so neither the SAT solver nor drat-trim is
trusted. The conditional theorems state their hypothesis; to use them one also trusts cake_lpr's check of the
LRAT proof of `q{d}.cnf` (cake_lpr is verified in HOL4 and compiled by the verified compiler CakeML), and the
`#guard` that `q{d}.cnf` is `formula` (run by Lean's interpreter when the file is built; the Python checker also
compares the file with the graph).
