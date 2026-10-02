# Formal proofs in Lean 4

This directory holds proofs in Lean 4, with Mathlib, of three theorems:
- the two theorems of the note
  [`papers/planes-4-chromatic/planes-4-chromatic.pdf`](../papers/planes-4-chromatic/planes-4-chromatic.pdf): the planes
  over ℚ(√2, √3) and ℚ(√3, √11) have chromatic number 4. The second is K. G. Fischer's theorem (1994); the note
  gives a short proof of both, and these files check that proof;
- the plane over ℚ(√11) has chromatic number 4, the first real quadratic field of
  [`notes/quadratic_planes.md`](../notes/quadratic_planes.md). The lower bound is the graph with 76 vertices of
  [`data/quadratic_planes/q11.json`](../data/quadratic_planes/q11.json), whose SAT certificate the kernel checks;
  the upper bound is Moorhouse's reduction at 7.

## The theorems

```lean
theorem Q23.chromaticNumber_eq_four : (LocalColouring.unitDistGraph Q23.L).chromaticNumber = 4
theorem Q311.chromaticNumber_eq_four : (LocalColouring.unitDistGraph Q311.L).chromaticNumber = 4
theorem Q11.chromaticNumber_eq_four : (LocalColouring.unitDistGraph Q11.L).chromaticNumber = 4
```

Here `Q23.L`, `Q311.L` and `Q11.L` are the subfields ℚ(√2, √3), ℚ(√3, √11) and ℚ(√11) of ℝ, and
`unitDistGraph K` joins two points of `K × K` at Euclidean distance 1:

```lean
noncomputable def Q23.L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√2, √3}
noncomputable def Q311.L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√3, √11}
noncomputable def Q11.L : IntermediateField ℚ ℝ := IntermediateField.adjoin ℚ {√11}

def LocalColouring.unitDistGraph (K : IntermediateField ℚ ℝ) : SimpleGraph (K × K) where
  Adj p q := ((p.1 : ℝ) - q.1) ^ 2 + ((p.2 : ℝ) - q.2) ^ 2 = 1
```

`chromaticNumber` is Mathlib's `SimpleGraph.chromaticNumber`, with values in `ℕ∞`. The three theorems depend only on
Lean's three standard axioms, `propext`, `Classical.choice` and `Quot.sound`: `axioms.expected` records the
output of `#print axioms`, and CI compares them.

## Checking the proofs

With [elan](https://github.com/leanprover/elan) installed:

    cd lean
    lake exe cache get                       # Mathlib's prebuilt files, a few minutes
    lake build                               # the four files, a few minutes
    lake env lean PrintAxioms.lean           # compare with axioms.expected

`lean-toolchain` pins Lean v4.34.1 and `lakefile.toml` pins Mathlib v4.34.1. The GitHub Actions workflow
[`lean.yml`](../.github/workflows/lean.yml) runs these steps on every push to `main` and every pull
request. It then replays each file in Lean's kernel with `leanchecker`, independently of the tactics that
produced the proofs.

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

For ℚ(√11) (`notes/quadratic_planes.md`):

| the note | Lean |
|---|---|
| §2, the points `[a, b, c, e]` and the identities (1) | `Q11.pt`, `Q11.unitPairB`, `Q11.adj_pt` |
| the graph `data/quadratic_planes/q11.json` | `Q11.P`, `Q11.E`, `Q11.checkEdges_E` (every edge has length 1) |
| §4, the formula `q11.cnf` and its proof | `Q11.q11Refuted` (Mathlib's `lrat_proof`, from `q11.cnf` and `q11.lrat`) |
| no 3-colouring | `Q11.graph_not_colorable`, `Q11.not_colorable_three` |
| §5, Moorhouse's Lemma 8.2 at 7 | `Q11.exists_valuationSubring_seven_mem_nonunits`, `Q11.exists_int_sub_mem` (residue field 𝔽₇, by Hensel's lemma for √11), `Q11.core_lemma_seven`, `Q11.sumSqGraph_colorable`, `Q11.colorable_four` |
| the theorem for d = 11 | `Q11.chromaticNumber_eq_four` |

The proofs follow the note with one difference. The note takes a prime of the ring of integers above 2 and
reads its residue field from the decomposition of 2. Here the place above 2 is a valuation subring with 2 in
its maximal ideal, given by Chevalley's extension theorem, and its residue field is shown to be 𝔽₂ directly:
- for ℚ(√2, √3), with the uniformiser π = (√2 + √6)/2 − 1, a root of the 2-Eisenstein polynomial
  X⁴ + 4X³ + 2X² − 4X − 2;
- for ℚ(√3, √11), where 2 has two places, with π = √3 − 1 and a Hensel-type argument for (1 + √33)/2.

The valuation subring and the coset representatives come from Zorn's lemma and `Quotient.out`, so the
colourings exist but cannot be computed.

## Files

| file | contents |
|---|---|
| `LocalColouring.lean` | the graph, the criterion, and the tools both fields use |
| `Q23.lean` | ℚ(√2, √3): the place above 2, and the 10-vertex graph of [`data/chain23.json`](../data/chain23.json) |
| `Q311.lean` | ℚ(√3, √11): the place above 2, and the Moser spindle |
| `Q11.lean` | ℚ(√11): the graph of `data/quadratic_planes/q11.json` and the place above 7. It reads `data/quadratic_planes/q11.cnf` and the LRAT proof `data/quadratic_planes/q11.lrat` (kissat, then `drat-trim -L`) with `include_str` |
| `PrintAxioms.lean`, `axioms.expected` | the axioms of the three theorems |
| `tools/q23_coefficients.py`, `tools/q311_coefficients.py` | sympy scripts that produce the coefficients of the `linear_combination` steps and the edge lists |
| `tools/q11_lean.py` | writes `Q11.lean` from the data (`--check` compares; the tests run it) |

## What is trusted

Lean's kernel; Mathlib's definitions of `IntermediateField.adjoin`, `Real.sqrt` (written `√`) and
`SimpleGraph.chromaticNumber`; and the definitions of `Q23.L`, `Q311.L`, `Q11.L` and `unitDistGraph` quoted above.
For ℚ(√11), Mathlib's `lrat_proof` turns the LRAT proof into a proof term that the kernel checks; nothing
runs as native code, so neither the SAT solver nor drat-trim is trusted.
