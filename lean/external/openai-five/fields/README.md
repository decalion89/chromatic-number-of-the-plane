# OpenAI's proof over smaller fields

OpenAI's Lean proof that the plane is not 5-colourable (see [`../README.md`](../README.md)) works on the plane over
the algebraic numbers `ℚ̄`. These patches run the same proof on three smaller countable fields:

| patch | field `E` | statement file | theorem |
|---|---|---|---|
| `constructible.patch` | the constructible numbers: the smallest subfield of ℂ closed under square roots (ruler and compass) | `ConstructiblePlane.lean` | `sqrtClosure_plane_not_five_colourable` |
| `origami.patch` | the origami numbers: closed under square and cube roots (paper folding) | `OrigamiPlane.lean` | `sqrtCbrtClosure_plane_not_five_colourable` |
| `radicals.patch` | the numbers expressible by radicals (`solvableByRad ℚ ℂ`) | `RadPlane.lean` | `solvableByRad_plane_not_five_colourable` |

Each theorem says that every colouring of `E` with five colours gives the same colour to two points at
distance 1. For example:

```lean
theorem sqrtClosure_plane_not_five_colourable :
    ¬ ∃ c : (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ) → Fin 5,
      ∀ x y : (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ),
        ‖(x : ℂ) - (y : ℂ)‖ = 1 → c x ≠ c y
```

`ConstructibleFinite.lean` proves the finite form: some finite set of constructible points has no proper 5-colouring.
The step from the infinite statement is the de Bruijn–Erdős compactness theorem, which Mathlib has as
`SimpleGraph.nonempty_hom_of_forall_finite_subgraph_hom`:

```lean
theorem sqrtClosure_finite_set_not_five_colourable :
    ∃ S : Set (sInf {s : IntermediateField ℚ ℂ | ∀ x : ℂ, x ^ 2 ∈ s → x ∈ s} : IntermediateField ℚ ℂ),
      S.Finite ∧ ¬ ∃ c : S → Fin 5, ∀ x y : S, ‖((x : _) : ℂ) - ((y : _) : ℂ)‖ = 1 → c x ≠ c y
```

The fields are nested, constructible ⊂ origami ⊂ radicals, so the constructible case implies the other two. They were
compiled first, as steps towards it.

[`notes/six_over_fields.md`](../../../../notes/six_over_fields.md) explains the mathematics:
- what OpenAI's proof uses of the field;
- why square roots suffice;
- what the patches change.

In brief:
- a new file `Field.lean` defines `E` and proves the few facts the proof uses about it;
- for the origami and constructible fields, the finite-image lemma used in Lemma 2.4 is replaced by one that needs
  only square roots;
- for the constructible field, the triple average of the radial inequality is taken at the exponents −1, 0, 1 instead
  of 0, 1, 2. The middle factor is fixed and handled separately, which needs two new lemmas.

## What was checked (9 October 2026)

- **Sources.** The patches apply to openai/math at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, the commit
  checked in [`../README.md`](../README.md). Each one adds `Field.lean` and changes 10 files of the import closure
  (radicals) or 11 (origami, constructible).
- **Forbidden constructs.** No added line contains `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by`,
  `extern`, `unsafe`, `instance`, `notation`, `macro`, `syntax` or `set_option`.
- **Toolchain.** Lean 4.34.1 and Mathlib at commit d13f23b7, as for the unpatched proof.
- **Build.** Each field built its 72 modules (OpenAI's 70, `Field.lean` and the statement file) with no error. The
  first builds of the radical and origami cases were incremental: only the modules that import a changed file were
  recompiled (22 and 24), and the others were the unpatched ones.
- **Axioms.** The axioms printed for all seven theorems are:

  ```
  'OAI.EuclideanFiveColor.constructible_plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'OAI.EuclideanFiveColor.sqrtClosure_plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'OAI.EuclideanFiveColor.origami_plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'OAI.EuclideanFiveColor.sqrtCbrtClosure_plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'OAI.EuclideanFiveColor.radical_plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'OAI.EuclideanFiveColor.solvableByRad_plane_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  'OAI.EuclideanFiveColor.sqrtClosure_finite_set_not_five_colourable' depends on axioms: [propext, Classical.choice, Quot.sound]
  ```

- **Clean rebuild.** The constructible case was also rebuilt from a fresh copy of the pristine sources, the patch
  and `build_closure.py`, into an empty directory: 72 built, 0 failed, in 1 174 s on 2 cores. The compiler printed
  three linter warnings about style, in `Field.lean` and `Basic.lean`, and no error.
- **Kernel replay.** `leanchecker -v M` ran for each of the 72 modules of that clean build. Every module exited 0, in
  687 s in all. It also passed on `ConstructibleFinite.lean`.
- **Script test.** `build_field.sh` was run from scratch twice, as a reader would run it:
  - on the origami field with one job: 72 built, 0 failed, in 1 550 s, with the axioms above;
  - on the constructible field with two jobs: 72 built, 0 failed, in 1 196 s, with the axioms above for its three
    theorems; the compiler printed the same three style warnings and no error.
- **Search path.** In these builds the search path also held the unpatched build, after the new one. It is never read:
  `build_closure.py` compiles each module only after every module it imports is in the new directory, which comes
  first. To make sure, the replay of the constructible script test used a search path of the new build and Mathlib
  only: all 72 modules and `ConstructibleFinite.lean` exited 0, in 627 s for the 72.

## To reproduce

```sh
git clone https://github.com/openai/math && git -C math checkout fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb
cd lean && lake build && cd ..          # the Mathlib build of this repository
LP=$(cd lean && lake env printenv LEAN_PATH)
LEAN=$(cd lean && lake env which lean)
sh lean/external/openai-five/fields/build_field.sh constructible math/lean "$LEAN" "$LP" /tmp/oai-constructible 2
```

The last lines printed are the axioms of the two theorems of the field, and for the constructible field of the
finite form too. The build takes about 20 minutes on 2 cores.

## Licence

The patches modify files of [openai/math](https://github.com/openai/math), which are under the Apache License 2.0.
A copy of that licence is in [`LICENSE-openai-math`](LICENSE-openai-math). The patches, including the new
`Field.lean` files, are distributed under the same licence. The three statement files and `build_field.sh` are part
of this repository.
