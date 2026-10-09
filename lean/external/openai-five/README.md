# A check of OpenAI's proof that χ(ℝ²) ≥ 6

In September 2026 OpenAI published a proof that the plane is not 5-colourable
([*The Euclidean plane is not five-colorable*](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf),
family 158 of [openai/math](https://github.com/openai/math)), formalized in Lean as
`OAI.EuclideanFiveColor.no_proper_five_coloring` in `lean/OAI/Geometry/PlaneColoring/Five.lean` of that repository.
The proof is not constructive and gives no finite graph. This directory records how we checked the formal proof
on 9 October 2026. It is not part of the build of `lean/`: the proof is OpenAI's, and its sources are not copied here.

## What was checked

- **Sources.** The import closure of `OAI.Geometry.PlaneColoring.Five` in openai/math at commit
  `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`: 70 modules, 30 336 lines, in `Analysis/PlaneSpectrum`,
  `Geometry/PlaneColoring` and `MeasureTheory/CompactFactors`, with no library besides Mathlib. None contains
  `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by`, `extern`, `unsafe` or new syntax, and none declares
  a `Norm` instance.
- **Toolchain.** Lean 4.34.1 and Mathlib at commit d13f23b7, the versions of `lean/` here and of openai/math.
- **Build.** `build_closure.py` compiled the 70 modules in dependency order against our Mathlib build: 70 built, 0
  failed, in 18 minutes on 2 cores, every compiler message file empty.
- **Statement.** `PlaneSix.lean` restates the theorem without OpenAI's definitions,
  `¬ ∃ c : ℂ → Fin 5, ∀ p q : ℂ, ‖p - q‖ = 1 → c p ≠ c q`, and the same for every `k ≤ 5` colours, and proves both
  from theirs. Its output:

  ```
  'plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'plane_not_k_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  ```

## χ(ℝ²) is 6 or 7, in the statement of formal-conjectures

Google DeepMind's [formal-conjectures](https://github.com/google-deepmind/formal-conjectures) states the bounds of
Erdős Problem 508 (`FormalConjectures/ErdosProblems/508.lean`) for
`χ(ℝ²) = SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ)`, the graph on `EuclideanSpace ℝ (Fin 2)` that
joins two points at distance 1; on 9 October 2026 `HadwigerNelsonAtMostSeven`, `HadwigerNelsonAtLeastFive` and
`HadwigerNelsonAtLeast4` had no proof there. [`../../PlaneSeven.lean`](../../PlaneSeven.lean), part of the build of
`lean/`, proves `χ(ℝ²) ≤ 7` by a brick colouring. `PlaneSixOrSeven.lean` here moves OpenAI's theorem to `χ(ℝ²)`
along the isometry `Complex.orthonormalBasisOneI.repr : ℂ ≃ₗᵢ[ℝ] EuclideanSpace ℝ (Fin 2)`: a 5-colouring of the graph
would give one of `ℂ`. It proves

```lean
theorem PlaneSeven.hadwigerNelson_atLeastSix : 6 ≤ SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ)
theorem PlaneSeven.hadwigerNelson_six_or_seven :
    SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) = 6 ∨
      SimpleGraph.chromaticNumber (UnitDistancePlaneGraph Set.univ) = 7
```

and the statements `hadwigerNelson_atLeastFive` and `hadwigerNelson_atLeast4`. The output of the build:

```
'PlaneSeven.hadwigerNelson_atLeastSix' depends on axioms: [propext, Classical.choice, Quot.sound]
'PlaneSeven.hadwigerNelson_six_or_seven' depends on axioms: [propext, Classical.choice, Quot.sound]
```

`leanchecker` replays `PlaneSeven` and `PlaneSixOrSeven` without error (9 October; OpenAI's 70 modules were replayed
before, see `fields/README.md`). The definition of `UnitDistancePlaneGraph` is formal-conjectures' with `dist` written
`Dist.dist` and `dist_comm` written `_root_.dist_comm`: under `import Mathlib` these names inside `namespace SimpleGraph`
mean Mathlib's graph distance. Our toolchain is Lean 4.34.1 and Mathlib d13f23b7; formal-conjectures used Lean 4.33.1
and Mathlib 0df444a3 on that date, and we have not built the files there.

## To reproduce

```sh
git clone https://github.com/openai/math && git -C math checkout fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
cd lean && lake build && cd ..          # the Mathlib build of this repository
LP=$(cd lean && lake env printenv LEAN_PATH)
LEAN=$(cd lean && lake env which lean)
python3 lean/external/openai-five/build_closure.py math/lean "$LEAN" "$LP" /tmp/oai-build 2
cp lean/external/openai-five/PlaneSix.lean math/lean/
(cd math/lean && LEAN_PATH=/tmp/oai-build:$LP "$LEAN" PlaneSix.lean)
(cd lean && lake build PlaneSeven)
cp lean/external/openai-five/PlaneSixOrSeven.lean math/lean/
(cd math/lean && LEAN_PATH=/tmp/oai-build:../../lean/.lake/build/lib/lean:$LP "$LEAN" PlaneSixOrSeven.lean)
```

## Smaller fields

[`fields/`](fields/README.md) runs the same proof over three smaller fields: the constructible numbers (ruler and
compass), the origami numbers and the numbers expressible by radicals. Already the constructible plane is not
5-colourable, and some finite set of constructible points has no proper 5-colouring. The folder holds the patches,
the statements, the axioms printed and a build script. The mathematics is in
[`notes/six_over_fields.md`](../../../notes/six_over_fields.md) and
[`papers/ruler-compass/`](../../../papers/ruler-compass/README.md).
