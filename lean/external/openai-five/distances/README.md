# OpenAI's transfer for any set of algebraic forbidden distances, in Lean

OpenAI's proof that the plane is not 5-colourable first turns an arbitrary proper colouring of the plane into a
*weak measurable* one (Theorem 1.3 of their paper, `proper_to_borel_weak` in their Lean library). Theorem A of
[`papers/forbidden-distances/`](../../../../papers/forbidden-distances/README.md) says that the same proof works for
every set `D` of positive real algebraic forbidden distances. This folder holds a Lean 4 proof of that statement,
built on OpenAI's formalization. Like [`../fields/`](../fields/README.md), it is not part of the build of `lean/`.

## What is proved

The hypotheses are `hpos : ∀ r ∈ D, 0 < r`, `halg : ∀ r ∈ D, IsAlgebraic ℚ r` and a colouring `c : ℂ → Fin k` with
`hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y`. `D : Set ℝ` is arbitrary otherwise; that it is countable is proved
(`countable_of_algebraic`), not assumed. [`OAI/Ext/Check.lean`](OAI/Ext/Check.lean) restates the results and prints
their axioms:

```lean
theorem borel_transfer_distances {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ d : ℂ → Fin k, Measurable d ∧ ∀ r ∈ D, ∀ᵐ p ∂(volume : Measure ℂ).prod angularMeasure,
      d p.1 ≠ d (p.1 + (r:ℂ) * direction p.2)

theorem measurable_classes_avoiding {k : ℕ} [NeZero k] {D : Set ℝ}
    (hpos : ∀ r ∈ D, 0 < r) (halg : ∀ r ∈ D, IsAlgebraic ℚ r)
    (c : ℂ → Fin k) (hc : ∀ x y : ℂ, ‖x - y‖ ∈ D → c x ≠ c y) :
    ∃ B : Fin k → Set ℂ, (∀ i, MeasurableSet (B i)) ∧ Pairwise (fun i j => Disjoint (B i) (B j)) ∧
      (∀ i, ∀ x ∈ B i, ∀ y ∈ B i, ‖x - y‖ ∉ D) ∧ volume (⋃ i, B i)ᶜ = 0
```

The first is the weak measurable colouring of Theorem A: for every `r ∈ D`, almost every pair `(x, x + r e^{it})` gets
two colours. The second is the form used in the paper's two-scale theorem: `k` pairwise disjoint measurable sets, none
containing two points at a distance in `D`, covering almost all of the plane. `Check.lean` also proves
`borel_transfer_distances_density` (no two density-one points of a colour class are at a distance in `D`),
`measurable_proper_on_conull`, `borel_transfer_all_vectors` (every displacement of length in `D` is separated almost
everywhere) and `borel_transfer_unit`, which is OpenAI's `proper_to_borel_weak` recovered as the case `D = {1}`.

The density half of the paper (Theorem B, and the bounds `χ ≥ 1/m_D`) is not formalized.

## The files

| file | lines | contents |
|---|---:|---|
| `OAI/Ext/LabelModelD.lean` | 159 | OpenAI's label model (`PlaneFiveColor.LabelModel` in `Representations.lean`) with the distance `1` replaced by the set `D` |
| `OAI/Ext/BorelTransferD.lean` | 324 | OpenAI's `BorelTransfer.lean` generalised: correlations at every radius in `D`, the label field, one sample for countably many radii, the main theorem |
| `OAI/Ext/WeakColoringD.lean` | 213 | OpenAI's `WeakColoring.lean` at radius `r`: density-one points, measurable classes, the case `D = {1}` |
| `OAI/Ext/Check.lean` | 78 | restatements and `#print axioms` |

The changes are confined to the properness predicate (`‖x - y‖ = 1` becomes `‖x - y‖ ∈ D`), the correlation lemmas,
which now hold at every radius `r ∈ D` (OpenAI's `residual_fourier` already covers every nonzero algebraic
displacement and is used unchanged), and the sampling step, which serves countably many radii at once. The spectral
and ergodic part of OpenAI's proof, including their rigidity theorem, is imported unchanged.

## What was checked (10 October 2026)

- **Sources.** OpenAI's 33 modules in the import closure of `OAI.Ext.Check`, from openai/math at commit
  `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, unmodified, and the four files here. None of the four contains `sorry`,
  `admit`, `axiom`, `native_decide`, `implemented_by` or `extern`.
- **Toolchain.** Lean 4.34.1 and Mathlib at commit d13f23b7, as for [`../README.md`](../README.md).
- **Build.** The 37 modules compiled with no error and no warning, once on top of an existing build of OpenAI's modules
  and twice from source with only Mathlib on the search path: by hand (1036 s on 2 cores), and with
  `build_distances.sh` itself on a fresh copy of OpenAI's sources (1383 s with one job; same axioms).
- **Axioms.** Every theorem listed in `Check.lean` depends only on `propext`, `Classical.choice` and `Quot.sound`:

  ```
  'borel_transfer_distances' depends on axioms: [propext, Classical.choice, Quot.sound]
  'borel_transfer_distances_density' depends on axioms: [propext, Classical.choice, Quot.sound]
  'measurable_proper_on_conull' depends on axioms: [propext, Classical.choice, Quot.sound]
  'borel_transfer_all_vectors' depends on axioms: [propext, Classical.choice, Quot.sound]
  'measurable_classes_avoiding' depends on axioms: [propext, Classical.choice, Quot.sound]
  'borel_transfer_unit' depends on axioms: [propext, Classical.choice, Quot.sound]
  ```

  and the same for the eight `OAI.Ext` theorems they restate.
- **Kernel replay.** `leanchecker` replays all 37 modules of the closure without error.

## To reproduce

With a checkout of openai/math at the commit above and the Mathlib build of `lean/`, as in
[`../README.md`](../README.md):

```sh
LP=$(cd lean && lake env printenv LEAN_PATH)
LEAN=$(cd lean && lake env which lean)
sh lean/external/openai-five/distances/build_distances.sh math/lean "$LEAN" "$LP" /tmp/oai-distances 2
```

The script copies OpenAI's `OAI/` tree to a new directory, adds `OAI/Ext/`, compiles the closure with
[`../build_closure.py`](../build_closure.py) and prints the axioms. One compilation can take 6 GB of memory.

## Licence

The four `.lean` files adapt files of [openai/math](https://github.com/openai/math), which are under the Apache
License 2.0; a copy is in [`LICENSE-openai-math`](LICENSE-openai-math), and these files are distributed under the same
licence. `build_distances.sh` and this README are part of this repository and under its licence.
