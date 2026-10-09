# `scripts/g13/chi`: the code of `χ(G₁₃) = 6`

This directory holds the code and the fixed inputs of the computation that shows that `G₁₃`, the
anisotropic plane over `𝔽₁₃`, has no proper 5-colouring. The note is `notes/g13_chi.md`; it states
the argument, and says what was computed and how to check it.

The files come from the directory `g13-plan-d/` of commit `97047ff3` (the branch of plan D,
27 September 2026, in the private development repository; not in the public history), the version that the last shares ran. They are copied unchanged, byte for byte and
at the same relative paths, so that `share_driver.py` and `verify_plan_D.py` run here as they ran
there. No path had to change. New files: `pack_certificates.py`; the helpers `help_driver.py` and `sub_help.py`, and
the re-splits of leaf 430, `hard430.py` and `hard430c.py` with `work/vars_c12_lex.txt` and `work/vars_c012_lex.txt`
(added after the run, as they ran, except that `hard430.py` no longer names the machine's directory); and this README, which
replaces the README of the plan and keeps its account of the shares. The internal review
`review_plan_D.md` is added as it was written; `verify_plan_D.py` cites it.

## The files

| file | contents |
|---|---|
| `plan_C.py` | the lemma of the case split, in its docstring, items (1) to (5) (`notes/g13_chi.md` §2). The functions of the file belong to plan C, an enumeration of orbits of maximal independent sets that was not used; they need `work/misample`, `work/orbits.py` and `lists/`, which are not here |
| `g13.py`, `g13cnf.py`, `enum_cert.py` | the graph and the writers of the case formulas. `g13.py` and `enum_cert.py` are the files of `scripts/g13/`; `g13cnf.py` differs from its copy there in one line of its docstring |
| `cuber2.py` | the cube trees, and `check_cover()`. It differs from `scripts/g13/cuber2.py` by the fix of plan E: the unit clauses of a formula go to the solver as assumptions, so that a re-split branches on variables that its leaf has not fixed. `check_cover()` is the same |
| `certify.py` | kissat with a DRAT proof, then drat-trim, on a formula or on the leaves of a cube file; one log line per formula. The code of `scripts/g13/certify.py`, with two more lines of docstring |
| `share_driver.py` | runs one share of the computation; `share_driver.py regen` writes the formulas of the re-splits again from their parents |
| `help_driver.py`, `sub_help.py` | the helpers: they certify leaves of a share, or sub-leaves of one leaf, from the other end of its list, exactly as `share_driver.py` does, and stop where the share's own run starts |
| `hard430.py`, `hard430c.py` | leaf 430 of `F34` (`notes/g13_chi.md`): the re-split of its sub-leaf `F34_split430_split210_split90_split35_leaf1` on the variables of colours 1 and 2, and of its sub-leaves 211, 214, 215 and 223 cutting on colour 0 and then on colours 1 and 2 |
| `verify_plan_D.py` | the final check: it prints `PLAN D FULLY CERTIFIED` only when every leaf of every case is certified |
| `check_colouring.py` | checks a colouring of `G₁₃`, or a kissat model; a satisfiable leaf of `F36`, `F35` or `F34` would be a proper 5-colouring |
| `cnc_case.sh`, `resplit.sh` | the cube-and-conquer scripts of the first plan; `share_driver.py` does their work |
| `planE.json` | the leaves of shares 15 to 48 (plan E) |
| `cases/E37_B.icnf`, `F36.icnf`, `F35.icnf`, `F34.icnf` | the cube trees of the four case formulas, made once and committed so that every machine split alike. `E37_B.icnf` is `certificates/g13_alpha_part_b_cubes.icnf`; `F35.icnf` and `F34.icnf` are the same file |
| `cases/SHA256SUMS` | the SHA-256 of the eight case formulas |
| `work/vars_c0_lex.txt`, `work/vars_s_lex.txt` | the variables the cube trees branch on: `x(v, 0)`, or `s_v`, for the vertices `v` along `lex_order()` |
| `work/vars_c12_lex.txt`, `work/vars_c012_lex.txt` | the variables of colours 1 and 2, and of colours 0, 1 and 2, in the same order, for the re-splits of leaf 430 |
| `review_plan_D.md` | the internal review of the plan (written inside the project), written before the run and before plan E |
| `.gitignore` | the formulas and proofs that a run writes here |
| `pack_certificates.py` | new: packs the merged logs into `certificates/g13_chi_certlogs.tar.gz` |

Some docstrings name files that are not in this repository: `run_plan_D.sh`, the script of the first plan
on one machine, which `share_driver.py` replaced; `summarize.py`, the check that the review found unsound
and that `verify_plan_D.py` replaced; `certify_R.py` and `rest4.py` of plan C; and
`scripts/g17_part_b.py`, which `certify.py` and `cuber2.py` name as their model.

## The case formulas

`share_driver.py` and `verify_plan_D.py` write them with the same commands:

| case | command |
|---|---|
| `E37_A6`, `E37_A7`, `E37_A9`, `E37_A11` | `enum_cert.py OUT 37 --L 0 --rosette c` |
| `E37_B` | `enum_cert.py OUT 37 --norosette` |
| `F36` | `g13cnf.py OUT --big0 36 --cap 36 --dom0 --lex0 25 --vp 1234 --domcap 36` |
| `F35` | `g13cnf.py OUT --big0 35 --cap 35 --dom0 --lex0 25 --vp 1234 --domcap 35` |
| `F34` | `g13cnf.py OUT --rigid34 --lex0 25` |

The five formulas `E37_*` are those of `α(G₁₃) ≤ 36` (`notes/g13.md`), with the same SHA-256.

## How the work was split

Each cube file has about 4 820 leaves, made by `cuber2.py` with at most 14 decisions. The shares ran on
several machines from 27 September 2026:

| share | work |
|---|---|
| 1 | the four formulas `E37_A_c` and every leaf of `E37_B` (the certificates of `α(G₁₃) = 36`) |
| 2 | every leaf of `F36` (at first also every leaf of `F35`, which shares 6 and 8 took over) |
| 3, 4, 5 | the leaves `i < 1016` of `F34` with `i mod 3 = 0, 1, 2` |
| 6 | the leaves `i < 759` of `F35` |
| 7 | the 100 leaves `i ≥ 1016` of `F34` on which kissat reached 120 s in the first pass of shares 3 to 5, re-split at once |
| 8 | the leaves `i ≥ 1400` of `F35`; then the 38 leaves `759 ≤ i < 1400` of `F35` on which kissat reached 120 s in the first pass of share 6, re-split at once |
| 9, 10 | hard leaves of `F34` below 1016 (`i mod 3 = 0` and `2`) that shares 3 and 5 had not re-split, from the top down |
| 11 to 14 | every leaf of `F34` and `F35` still without a complete answer at 19:42 that no share above took from the other end |
| 15 to 48 | plan E: every leaf of `F34` (358) and `F35` (50) still without a complete answer at 21:09, 12 per share (`planE.json`) |

`python3 share_driver.py K` runs share `K`. It certifies each leaf with a kissat limit of 120 s. A leaf
without a VERIFIED line is re-split: its formula, the case formula with the leaf's literals as unit
clauses, becomes a case of its own, which `cuber2.py` cuts with 8 more decisions. Plan D names it
`CASE_leafI` and gives its leaves 1 200 s, then 3 600 s at up to three further levels; plan E names it
`CASE_splitI` and gives its leaves 1 200 s at up to seven levels. A satisfiable leaf would stop the
share (`SAT`).

**Plan E.** The re-splits of plan D did not split. `cuber2.py` branched on the next variable of its list
that propagation had not assigned; but a leaf formula carries its cube as unit clauses, and the solver
does not report what those imply. So each of the 8 decisions of a re-split was on a variable whose value
the leaf already implies by unit propagation: one branch contradicted the leaf and was closed at once,
and the other added a unit clause that the leaf implies. Every re-split of plan D has 8 closed leaves and
one open leaf, whose formula is the leaf's formula with 8 implied unit clauses. Such a re-split is still a
valid certificate, since its leaves cover the leaf; it only amounts to solving the leaf again with a
longer limit. From 27 September, 21:20, `cuber2.py` passes the unit clauses as assumptions, and a re-split
of plan E has up to 256 leaves on new variables.

## Checking the result

Merge every share's `cases/*.certlog`, `cases/*_leaf*.icnf`, `cases/*_split*.icnf` and `logs/E37_A.log`
into one copy of this directory. Then, in that copy:

    python3 share_driver.py regen      # the case formulas and the re-split formulas, from the code
    python3 verify_plan_D.py           # PLAN D FULLY CERTIFIED, or the first leaf that is not

`verify_plan_D.py` writes each case formula with the code and compares it with the one that was
certified. It checks each cube file for a cover, and recomputes the SHA-256 of each leaf formula. A leaf
counts if a drat-trim VERIFIED line of the log has its name and that SHA-256, or if it was re-split, in
either way, into a case that is fully certified, recursively.

`pack_certificates.py COPY` packs the logs and the cube files of the re-splits of such a copy into
`certificates/g13_chi_certlogs.tar.gz`, with their SHA-256 in `certificates/g13_chi_SHA256SUMS.txt`.
`python3 scripts/verify_g13_chi.py`, from the root of the repository, unpacks the archive into a
temporary copy of this directory and runs the two commands above; with `--stats` it prints the numbers
of the note.
