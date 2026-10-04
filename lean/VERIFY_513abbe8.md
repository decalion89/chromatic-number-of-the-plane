# Lean verification of commit 513abbe8

Run in place of the `lean` GitHub Actions job (`.github/workflows/lean.yml`), whose minutes are exhausted.

- Date: 2026-10-04 (UTC)
- Commit: 513abbe8 ("Merge the Lean proof of Proposition B9 (branch claude/lean-two-roots) [skip ci]"), head of branch `claude/quadratic-planes`; record on branch `claude/lean-verify-513abbe8`
- Lean: 4.34.1 (`leanprover/lean4:v4.34.1`, from `lean-toolchain`); Lake 5.0.0-src+5045d00
- Mathlib: tag v4.34.1, rev d13f23b723b8a846827a245b89c10fc7d3f11612 (from `lake-manifest.json`)
- Machine: Linux x86_64 container, 4 cores, 15 GB RAM

## Results

| Step | Command | Wall-clock | Result |
|---|---|---|---|
| 1a. Toolchain | elan-init.sh (`--default-toolchain none`), toolchain fetched via `lean-toolchain` | 43 s | passed |
| 1b. Mathlib cache | `lake exe cache get` (8908 files) | 134 s | passed |
| 2. Build | `lake build` (all 41 default targets) | 756 s | passed: "Build completed successfully (9005 jobs).", no warnings or errors in the log |
| 3. Axioms | `lake env lean PrintAxioms.lean > axioms.txt && diff axioms.expected axioms.txt` | 10 s | passed |
| 4. Kernel replay | `lake env leanchecker -v $m` for each of the 41 modules of the workflow's loop | 1637 s | passed: every module exited 0 |
| 5. Forbidden words | grep over the 42 `.lean` files under `lean/` (excluding `.lake/`) | <1 s | passed: no match |

## Axiom diff

Output of `diff axioms.expected axioms.txt`: empty.

## Kernel replay, per module

```
LocalColouring rc=0 11s
Q23 rc=0 11s
Q311 rc=0 11s
QuadraticPlanes rc=0 17s
Sqrt11 rc=0 15s
Sqrt119 rc=0 257s
Sqrt131 rc=0 149s
Sqrt179 rc=0 129s
Sqrt191 rc=0 18s
Sqrt251 rc=0 124s
Sqrt431 rc=0 147s
Sqrt455 rc=0 15s
Sqrt911 rc=0 146s
Sqrt935 rc=0 92s
ColouringFormula rc=0 11s
Sqrt23 rc=0 17s
Sqrt35 rc=0 20s
Sqrt47 rc=0 25s
Sqrt59 rc=0 17s
Sqrt71 rc=0 21s
Sqrt95 rc=0 40s
Sqrt155 rc=0 38s
Sqrt239 rc=0 15s
Sqrt263 rc=0 16s
Sqrt359 rc=0 24s
Sqrt443 rc=0 25s
Sqrt491 rc=0 25s
Sqrt599 rc=0 23s
Sqrt611 rc=0 23s
Sqrt791 rc=0 27s
Sqrt851 rc=0 19s
Sqrt959 rc=0 19s
PadicPlanes rc=0 11s
TheoremW rc=0 7s
Recurrence rc=0 8s
TheoremWplus rc=0 7s
DistLiu rc=0 7s
TheoremWInf rc=0 8s
FourColours rc=0 20s
PadicFour rc=0 10s
TwoRoots rc=0 12s
```

The leanchecker output was one line per module (`replaying <module>`), with nothing else.

## Forbidden-word grep

```
grep -rnwE 'sorry|admit|native_decide|implemented_by|extern' --include='*.lean' --exclude-dir=.lake lean/
grep -rnE '^\s*(private\s+|protected\s+|noncomputable\s+)*axiom\b' --include='*.lean' --exclude-dir=.lake lean/
```

Both returned no matches (exit code 1). Files searched: the 41 module files of the default targets plus `PrintAxioms.lean`.

## Decisions

- The steps were run as separate commands, so a failure in one would not have hidden the others; within step 4 every module was run even if an earlier one failed (none did).
- The build log was checked with a case-insensitive grep for `error` and `warning`: 0 matching lines.
- The replay was run after the build and the axiom check, on the same machine, sequentially; the per-module times are longer than in `VERIFY_ec7b021a.md`, which reflects the machine, not the proofs.
- `axioms.txt` was generated in `lean/` as in the workflow and is not committed; no other file was changed.
