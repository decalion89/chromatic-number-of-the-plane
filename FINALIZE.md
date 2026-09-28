# Finalizing `χ(G₁₃) = 6`

This branch holds the kit for publishing `χ(G₁₃) = 6`, made while the last shares were still running:
the note `notes/g13_chi.md` with placeholders, the code of the run in `scripts/g13/chi/`, the checker
`scripts/verify_g13_chi.py`, the packer `scripts/g13/chi/pack_certificates.py`, and
`tests/test_g13_chi.py`. Nothing here claims the result until the steps below are done. Delete this file
in the last commit (step 8): it is not meant to stay.

Publishing (a push, a pull request, an announcement) is a separate decision, not one of these steps.

## 1. Merge the final logs

Merge every share branch (`g13-share-1` to `g13-share-48` on the remote) into one copy of `g13-plan-d/`
at commit `97047ff3`: all `cases/*.certlog` (lines concatenated, exact duplicates dropped), all
`cases/*_leaf*.icnf` and `cases/*_split*.icnf`, and `logs/E37_A.log`. The dry run of 28 September did
this with `assemble.py` in the session's scratchpad (`g13final/`). Before reusing it, check that its list
of shares covers every share branch, and that it reports no conflicts (two shares with different cube
files of the same name).

In the merged copy:

    python3 share_driver.py regen      # writes the re-split formulas: about 9 MB each
    python3 verify_plan_D.py

Go on only if the last line is `PLAN D FULLY CERTIFIED`. `verify_plan_D.py` stops at the first case that
fails, in the order `E37_B`, `F36`, `F35`, `F34`. In the dry run it stopped at `F35_leaf1001`, still
running then; but `F34` was further from done: 235 of its leaves of level 0 had no certificate yet
(`F34_leaf128` the first). `python3 scripts/verify_g13_chi.py --dir COPY --stats --no-verify` lists the
state of every case without writing the re-split formulas.

## 2. Pack the certificates

From the root of this branch:

    python3 scripts/g13/chi/pack_certificates.py COPY

It writes `certificates/g13_chi_certlogs.tar.gz` and `certificates/g13_chi_SHA256SUMS.txt`, and fails if a
log line is not a line of `certify.py` or the base cube files differ from the committed ones. The dry run
gave 753 files, 16.7 MB unpacked and 4.4 MB packed.

## 3. Check the archive, and fill the note

    python3 scripts/verify_g13_chi.py --stats --fill notes/g13_chi.md

It unpacks the archive into a temporary copy of `scripts/g13/chi` (give `--workdir DIR` on a disk with room:
about 9 MB for each cube file of a re-split, 3.5 GB in the dry run), runs `share_driver.py regen` and
`verify_plan_D.py` there, and prints the numbers. Go on only if:

- the last line begins with `CONFIRMED`;
- the statistics say `0 SAT` (the script fails otherwise) and name no case as `NOT fully certified`;
- the counts it prints agree with those of `verify_plan_D.py` (the script fails otherwise);
- every re-split of plan D has 9 leaves, 8 of them closed (the line "re-splits of plan D"), as §5 of the
  note says; if not, change that sentence.

`--fill` writes every placeholder it can compute into the note. Time the run: it gives `VERIFY_MINUTES`.

## 4. The placeholders left

Fill these by hand in `notes/g13_chi.md`:

| placeholder | source |
|---|---|
| `{{RUN_END}}` | the day the last share finished, as `28 September`: the last commits of the share branches |
| `{{MACHINES}}` | the number of machines the shares ran on |
| `{{DATE}}` | the day of step 3 |
| `{{VERIFY_MINUTES}}` | the time of step 3 |

Then:

- delete the HTML comment at the top of the note (the list of placeholders and the word DRAFT);
- check that `grep -n '{{' notes/g13_chi.md` prints nothing;
- check the sentence of §5 that names kissat 4.0.4 and drat-trim 2e3b2dc against the first line of each
  share's `progress.txt`, which names the binaries it ran; say so in the note if any share used others;
- read §6 once filled: the table of re-splits (`{{TREE_TABLE}}`) has one row per case and level, and
  "the deepest re-split" of the first table must agree with it.

The placeholders that `--stats` computes, for reference: `F36_TOP_VERIFIED`, `F35_TOP_VERIFIED`,
`F34_TOP_VERIFIED`, `F36_TOP_RESPLIT`, `F35_TOP_RESPLIT`, `F34_TOP_RESPLIT`, `F36_CERTIFIED`,
`F35_CERTIFIED`, `F34_CERTIFIED`, `F36_DEPTH`, `F35_DEPTH`, `F34_DEPTH`, `COLOURING_CERTIFIED`,
`DRAT_PROOFS`, `RESPLITS`, `TREE_TABLE`, `KISSAT_HOURS`, `DRAT_TRIM_HOURS`, `MAX_KISSAT_S`, `PROOF_GB`,
`ALL_LINES`, `TIMEOUTS`, `TIMEOUTS_120`, `TIMEOUTS_1200`, `TIMEOUTS_3600`, `ALL_KISSAT_HOURS`,
`ALL_DRAT_TRIM_HOURS`, `UNUSED_RESPLITS`, `REGEN_GB`, `ARCHIVE_FILES`, `ARCHIVE_MB`, `ARCHIVE_SHA256`. It
also prints `TOTAL_CERTIFIED`, `MAX_DEPTH`, `ALL_VERIFIED`, `SAT_LINES` and the numbers of `E37_B`, which the
note does not use.

## 5. Test

    python3 -m pytest -q tests/test_g13_chi.py tests/test_g13.py tests/test_g13_audit.py \
        tests/test_small_plane_colourings.py

With the archive in `certificates/`, `test_the_certificate_archive` runs too (it is skipped without it):
it checks the archive against its SHA256SUMS file and checks that the logs of `E37_*` in it are the
certificates of `α(G₁₃) = 36`.

## 6. Edits to files that the kit leaves alone

- `notes/g13.md`: the sentence "The theorem does not decide `χ(G₁₃)`, which is 5 or 6 (§7)" and §7 "Five
  or six colours?" now point to `notes/g13_chi.md`; for example, end §7 with "The computation is done:
  `χ(G₁₃) = 6` (`notes/g13_chi.md`)."
- `notes/local_colourings.md` §4, the row `13` of the table: `5–6` becomes `6`, citing `g13_chi.md`.
- `notes/README.md`: a row for `g13_chi.md`, after the row of `g13.md`; and in that row, "`χ(G₁₃)` stays 5
  or 6" goes.
- `scripts/README.md`: rows for `g13/chi/` and `verify_g13_chi.py`, after those of `g13/` and
  `verify_g13.py`.
- `tests/README.md`: a row for `test_g13_chi.py`; and the counts of files and tests, if it joins CI.
- `certificates/README.md`: rows for `g13_chi_certlogs.tar.gz` and `g13_chi_SHA256SUMS.txt`, and a sentence
  in the paragraph on the logs of `G₁₃`.
- `.github/workflows/tests.yml`: add `tests/test_g13_chi.py` to the list (half a minute; it needs numpy
  and python-sat, which CI installs).
- `notes/literature.md`, item 6 "Six colours for finite planes": `G₁₃` joins the planes that need six
  colours, by a SAT proof rather than the three-point bound.
- `docs/research-log.md`: the entry "`α(G₁₃) = 36` (27 September)" says that kissat reached 120 s on 2 to 4
  per cent of the leaves of `F36` and `F35`. In the dry run it did on 29 and 80 of 4 823 leaves (0.6 and
  1.7 per cent). If the final numbers agree (`F36_TOP_RESPLIT` and `F35_TOP_RESPLIT`), add a correction
  under that paragraph, as the log does for others, rather than changing it.

## 7. Drafts for the front page, the changelog and the research log

**`README.md`, the results table.** In the row of `α(G₁₃) = 36`, "Whether χ(G₁₃) is 5 or 6 stays open."
becomes "χ(G₁₃) = 6: see the next row." After it:

    | **χ(G₁₃) = 6.** The anisotropic plane over 𝔽₁₃ has no proper 5-colouring, although its fractional chromatic number is 169/36 < 5. | Computer proof | [`notes/g13_chi.md`](notes/g13_chi.md): a case split on the size of the largest colour class (34, 35 or 36 points, since α(G₁₃) = 36), one formula for each size with the symmetry broken, and cube and conquer: {{COLOURING_CERTIFIED}} leaves, each refuted by kissat with a DRAT proof checked by drat-trim (`certificates/g13_chi_*`). `scripts/verify_g13_chi.py`, `tests/test_g13_chi.py` |

**`README.md`, "Checking one result".** After the row of `α(G₁₃) = 36`:

    | χ(G₁₃) = 6, without a solver: the 6-colouring, the case formulas written again by the code, every cube tree, and a drat-trim VERIFIED line for every leaf | `python3 scripts/verify_g13_chi.py` | {{VERIFY_MINUTES}} minutes and {{REGEN_GB}} GB of disk |

In the row "Six colours for finite planes", `G₁₃` could be named as the one small case settled (by SAT).

**`CHANGELOG.md`, under `## [Unreleased]`, `### Added`, first:**

    - **`χ(G₁₃) = 6`** (`notes/g13_chi.md`): the anisotropic plane over 𝔽₁₃ has no
      proper 5-colouring. A computer proof. The largest class of a 5-colouring would
      have 34, 35 or 36 points (`α(G₁₃) = 36`); in a colouring where it is as large as
      possible it is dominating, an automorphism makes it lex-leader, and renaming
      the other colours gives them value precedence. One formula for each size says
      this; cube and conquer refutes the three, {{COLOURING_CERTIFIED}} leaves in all,
      each by kissat with a DRAT proof checked by drat-trim. `scripts/g13/chi/` is the
      code of the run, `certificates/g13_chi_certlogs.tar.gz` holds its logs, and
      `scripts/verify_g13_chi.py` checks them again. `tests/test_g13_chi.py` reads
      the three formulas (the counts and the chains with the audit of `α`) and tests
      them on intended models, and tests the checkers on a small made-up run. The
      proofs have not been checked by cake_lpr.

**`docs/research-log.md`, a new entry at the end:**

    ## `χ(G₁₃) = 6` ({{DATE}})

    The computation announced under "`α(G₁₃) = 36`" is done: `G₁₃` has no proper
    5-colouring, so `χ(G₁₃) = 6` (`notes/g13_chi.md`). The lower bound needed more than
    `α`: five classes of 36 points would hold 180 vertices.

    - **The case split.** In a 5-colouring whose largest class is as large as
      possible, every largest class is dominating, or a vertex could join it. Its size
      is 34, 35 or 36. An automorphism makes one such class lex-leader on the first 25
      vertices of `lex_order()`, and then the other colours are renamed by first
      appearance. One formula for each size (`F36`, `F35`, `F34`) says this. An
      independent review before the run found the lemma and the encodings sound, and
      the planned check of the logs not; `verify_plan_D.py` replaced that check.
    - **The run.** Each formula was cut into 4 823 leaves of at most 14 decisions, and
      each leaf refuted by kissat within 120 s, with a DRAT proof checked by
      drat-trim, or split again. The first re-splits did not split: `cuber2.py` did
      not see what the unit clauses of a leaf imply, so they only solved each leaf
      again with a longer limit. From 21:20 on 27 September, plan E passed those
      clauses as assumptions and cut each remaining leaf into up to 256 new ones. In
      all, {{COLOURING_CERTIFIED}} leaves of the three formulas were certified,
      {{F34_CERTIFIED}} of them for `F34`, in {{RESPLITS}} re-split formulas: kissat
      took {{KISSAT_HOURS}} hours and drat-trim {{DRAT_TRIM_HOURS}} on the lines that
      certify them, on {{MACHINES}} machines until {{RUN_END}}.
    - **The checks.** The logs of every share were merged, and `verify_plan_D.py`
      found every leaf of every tree certified. `scripts/verify_g13_chi.py` does the
      same from the archive `certificates/g13_chi_certlogs.tar.gz`. The tests read the
      three formulas: every clause is of a known kind; the counts and the chains pass
      the functions of the audit of `α`; and on colourings built from maximal
      independent sets, and on colourings far from any normal form, the clauses that
      fail are exactly those that their meaning predicts.

    Not done: an audit of `F34`, `F35` and `F36` from their text alone, as for the
    formulas of `α`; a second check of the proofs by cake_lpr; a proof in Lean.

## 8. Commit

    git add notes/g13_chi.md certificates/g13_chi_certlogs.tar.gz certificates/g13_chi_SHA256SUMS.txt \
        README.md CHANGELOG.md docs/research-log.md notes/g13.md notes/local_colourings.md notes/README.md \
        notes/literature.md scripts/README.md tests/README.md certificates/README.md .github/workflows/tests.yml
    git rm FINALIZE.md
    git commit -m "χ(G₁₃) = 6: the anisotropic plane over 𝔽₁₃ needs six colours"

Add only the files you changed.
