# Lean verification of commit ec7b021a

Run in place of the `lean` GitHub Actions job (`.github/workflows/lean.yml`), whose minutes are exhausted.

- Date: 2026-10-03 (UTC)
- Commit: ec7b021a ("Merge the Lean transfer to the p-adic planes (branch claude/lean-padic-four); Corollary 3 unconditional"), branch `claude/lean-verify`
- Lean: 4.34.1 (`leanprover/lean4:v4.34.1`, from `lean-toolchain`); Lake 5.0.0-src+5045d00
- Mathlib: tag v4.34.1, rev d13f23b723b8a846827a245b89c10fc7d3f11612 (from `lake-manifest.json`)
- Machine: Linux x86_64 container, 4 cores, 15 GB RAM

## Results

| Step | Command | Wall-clock | Result |
|---|---|---|---|
| 1a. Toolchain | elan-init.sh (`--default-toolchain none`), toolchain fetched via `lean-toolchain` | 24 s | passed |
| 1b. Mathlib cache | `lake exe cache get` (8908 files) | 100 s | passed |
| 2. Build | `lake build` (all 40 default targets) | 421 s | passed: "Build completed successfully (9003 jobs)", no warnings or errors in the log |
| 3. Axioms | `lake env lean PrintAxioms.lean > axioms.txt && diff axioms.expected axioms.txt` | 7 s | passed |
| 4. Kernel replay | `lake env leanchecker -v $m` for each of the 40 modules of the workflow's loop | 867 s | passed: every module exited 0 |

## Axiom diff

Output of `diff axioms.expected axioms.txt`: empty.

## Kernel replay, per module

```
LocalColouring rc=0 7s
Q23 rc=0 6s
Q311 rc=0 7s
QuadraticPlanes rc=0 10s
Sqrt11 rc=0 8s
Sqrt119 rc=0 132s
Sqrt131 rc=0 75s
Sqrt179 rc=0 66s
Sqrt191 rc=0 10s
Sqrt251 rc=0 63s
Sqrt431 rc=0 72s
Sqrt455 rc=0 9s
Sqrt911 rc=0 77s
Sqrt935 rc=0 50s
ColouringFormula rc=0 7s
Sqrt23 rc=0 11s
Sqrt35 rc=0 12s
Sqrt47 rc=0 14s
Sqrt59 rc=0 9s
Sqrt71 rc=0 13s
Sqrt95 rc=0 23s
Sqrt155 rc=0 22s
Sqrt239 rc=0 9s
Sqrt263 rc=0 10s
Sqrt359 rc=0 15s
Sqrt443 rc=0 13s
Sqrt491 rc=0 14s
Sqrt599 rc=0 12s
Sqrt611 rc=0 13s
Sqrt791 rc=0 14s
Sqrt851 rc=0 10s
Sqrt959 rc=0 10s
PadicPlanes rc=0 6s
TheoremW rc=0 4s
Recurrence rc=0 4s
TheoremWplus rc=0 5s
DistLiu rc=0 4s
TheoremWInf rc=0 4s
FourColours rc=0 11s
PadicFour rc=0 6s
```

The leanchecker output was one line per module (`replaying <module>`), with nothing else.

## Decisions

- The steps were run as separate commands, so a failure in one would not have hidden the others; within step 4 every module was run even if an earlier one failed (none did).
- `axioms.txt` was generated in `lean/` as in the workflow and is not committed; no other file was changed.
