# Contributing

This is a research repository, and the most useful contribution is a
correction. That can be a statement that is wrong, a proof with a gap, a table
entry without support, or a result that does not reproduce.

## Reporting an error

Open an issue with one of the two templates:

- **Mathematical error**: a statement, proof or table entry that is wrong or
  unsupported. Quote it, say where it is, and show what fails. The template
  adds the label `erratum`.
- **Result does not reproduce**: a test, certificate or command that fails or
  gives a different answer. Give the command, the commit, the versions (Python,
  kissat, drat-trim) and the output. The template adds the label
  `reproducibility`.

A confirmed error is corrected where it stands and recorded in two places:
- [`CHANGELOG.md`](CHANGELOG.md), under *Fixed*;
- the research log,
  [`docs/research-log.md`](docs/research-log.md),
  which keeps the earlier statement, marked as corrected, with the reason.

## Reproducing a result

The section [*Reproducing*](README.md#reproducing) of the README gives one
command per result, with its running time. At the root of the repository:
- `requirements-lock.txt` pins the Python packages;
- `scripts/worker_setup.sh` builds kissat 4.0.4 and drat-trim at a fixed
  commit.

## Changing code or claims

A pull request follows the principles stated in the README:
- **Exact geometry.** Coordinates live in number fields. Floating point may
  prune a search, but never decides that two points are at distance 1.
- **Evidence with every claim.**
  - A claim that a graph has no proper k-colouring needs a DRAT proof checked
    by drat-trim, against a formula rebuilt from the coordinates. Any symmetry
    breaking in that formula is stated and justified.
  - Claims that rest on solvers alone say so.
  - A semidefinite bound needs a stored dual certificate checked in interval
    arithmetic.
  - A formal proof in `lean/` builds with no `sorry`,
    and its theorems use only the axioms listed in `lean/axioms.expected`; a
    change to that file says why.
- **Tests.** Every new statement the code supports gets a test. Run the test
  files you touched, and the fast subset that CI runs
  (`.github/workflows/tests.yml`).
- **Records.** Add an entry under *Unreleased* in `CHANGELOG.md`. Mathematical
  content also gets an entry in the research log.

## Conduct

Criticism of arguments is welcome and is the point of the repository.
Criticism of people is not.
