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

## To reproduce

```sh
git clone https://github.com/openai/math && git -C math checkout fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
cd lean && lake build && cd ..          # the Mathlib build of this repository
LP=$(cd lean && lake env printenv LEAN_PATH)
LEAN=$(cd lean && lake env which lean)
python3 lean/external/openai-five/build_closure.py math/lean "$LEAN" "$LP" /tmp/oai-build 2
cp lean/external/openai-five/PlaneSix.lean math/lean/
(cd math/lean && LEAN_PATH=/tmp/oai-build:$LP "$LEAN" PlaneSix.lean)
```

## Smaller fields

[`fields/`](fields/README.md) runs the same proof over three smaller fields: the constructible numbers (ruler and
compass), the origami numbers and the numbers expressible by radicals. Already the constructible plane is not
5-colourable, and some finite set of constructible points has no proper 5-colouring. The folder holds the patches,
the statements, the axioms printed and a build script. The mathematics is in
[`notes/six_over_fields.md`](../../../notes/six_over_fields.md) and
[`papers/ruler-compass/`](../../../papers/ruler-compass/README.md).
