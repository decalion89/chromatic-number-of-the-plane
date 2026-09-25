# The Hadwiger–Nelson problem: an exact, certificate-driven attack

How many colours does the plane need so that no two points at distance exactly 1
share a colour? The answer, χ(ℝ²), has been known since 2018 to lie in
**{5, 6, 7}**.

| bound | value | who, when | how |
|---|---|---|---|
| lower | ≥ 4 | Nelson, 1950; L. and W. Moser, 1961 | the 7-vertex Moser spindle |
| lower | ≥ 5 | de Grey, 2018 | a 1581-vertex unit-distance graph with no 4-colouring |
| lower | ≥ 5 | Parts, 2020 | the same property on 509 vertices |
| upper | ≤ 7 | Isbell, 1950 | a hexagonal tiling |

This project aims at **χ(ℝ²) ≥ 6**. That needs one finite unit-distance graph
with no proper 5-colouring. A finite object can be searched for, and anyone can
check it once found.

**That graph has not been found.** Several results along the way are recorded
below, with their evidence.

## Principles

- **Exact arithmetic.** Coordinates live in number fields such as ℚ(√3, √11),
  never in floating point. "Distance exactly 1" is a decidable predicate.
- **Solvers, then certificates.** Colourability is decided by SAT solvers
  (kissat, CaDiCaL, Glucose, MiniSat). Positive answers ship as colourings,
  checked against edges rebuilt from the coordinates. Negative answers ship as
  DRAT proofs, checked by `drat-trim` against a formula rebuilt from the
  coordinates.
- **Nothing is announced unverified.** A claim of non-colourability needs:
  - an exact rebuild;
  - agreement between independent solvers;
  - a verified DRAT proof, with no hidden symmetry breaking.
- **Corrections stay visible.** Withdrawn claims are kept, with the reason, in
  the research log.

## Results

### Theorems

| statement | status | where |
|---|---|---|
| **χ(ℚ(√3, √11)²) = 4.** No 5-chromatic unit-distance graph has coordinates in ℚ(√3, √11). | Proved; short proof; unit tests. | `notes/local_colourings.md` §8, `hn/adelic.py`, `tests/test_q311.py` |
| χ(ℚ(√−3, √−11)) = 4 and χ(ℚ(√−3, √−11, √−247)) = 5, as whole complex fields | Proved. Upper bounds by reduction at the primes 2 and 11; lower bounds from the Moser spindle and Exoo–Ismailescu's graph. | `notes/rigidity.md`, `notes/local_colourings.md` |
| Necessary local conditions for a field to hold a 6-chromatic unit-distance graph | Proved | `notes/local_colourings.md` §5–§9, `scripts/fieldscreen.py` |
| Rigidity of coset colourings on the ρ₇ module | 3 840 exact Stiemke certificates | research log, "κ, the rotation every coset colouring is blind to" |

**On χ(ℚ(√3, √11)²) = 4.** The question was open in print:
- Moorhouse (2010) left the value undetermined.
- Madore (2015) proved 4 ≤ χ ≤ 5.
- Exoo and Ismailescu (2018) asked whether a 5-chromatic unit-distance graph
  embeds in this plane.
- Voronov (Polymath16, 2021) conjectured χ = 4 and noted that it was unproved.

The proof reduces z = x + iy modulo the place of ℚ(√3, √11) above 2, which is
inert in ℚ(i, √3, √11). Every unit vector becomes a nonzero element of 𝔽₄, so
the residue is a proper 4-colouring.

The 2-adic idea is David Speyer's: he used it in Polymath16 (2018) to
4-colour the Moser ring. The step here is its extension to the whole plane.

This result has not been refereed.

### Computations, all independently verified

| object | property | evidence |
|---|---|---|
| Moser spindle | no 3-colouring | `certificates/moser_spindle_no3coloring.json`, drat-trim |
| de Grey's 1581-vertex graph, rebuilt from his 39-point seed | no 4-colouring with one triangle pinned to 0, 1, 2 | `certificates/degrey_1581_no4coloring.json`, kissat plus drat-trim (13.1 M lemmas) |
| 19 vertices, 33 edges, vertex-critical | no 3-colouring; its forced pair is forced only jointly | `certificates/genuine_pair_19_no3coloring.json`, drat-trim |
| `data/five_247.json`, 1139 vertices, and `data/five_247_c.json`, 803 vertices, vertex-critical | 5-chromatic, in ℚ(√3, √11, √247). Not a record: Parts' 509 stands. | `tests/test_five_247.py`, research log |
| `data/W_moser_orbit_9_33.json`, 187 points | not 5-colourable with edges at 1 and at one Galois orbit of distances. Proves χ ≥ 6 once one gadget exists for that orbit. | four solvers plus drat-trim |
| `data/W_lattice_16_21_28_61.json`, 72 points | not 5-colourable with edges at 1, 4/√3, √7, √(28/3), √(61/3) | four solvers plus drat-trim |

### What is closed, and why

The research log records the approaches that cannot reach six, with the
reason for each:
- symmetrised growth;
- spindling at five colours;
- coset and circular colourings as obstructions;
- fields whose local planes are 5-colourable.

See "What has been ruled out so far", "The search over operations is closed, by
a theorem" and "Why every known construction stops at five".

## The route to six being searched now

The route is Exoo and Ismailescu's two-step reduction, which Polymath16 calls
"virtual edges":
1. a **witness**: a graph with edges at 1 and at some distances `d` that has
   no 5-colouring;
2. a **gadget** for each `d`: a unit-distance graph in which two points at
   distance `d` always get different colours.

Two observations are new here:
- **Repulsive distances.** Some distances, such as 2/√3, are coloured alike
  unusually rarely, so they are the natural gadget targets.
- **Galois orbits.** A Galois automorphism of the field maps gadgets to
  gadgets, so one gadget serves a whole orbit of distances. The 187-point
  witness above needs a single gadget.

The gadget searches run as parallel jobs, described in `notes/worker_jobs.md`.
None has succeeded yet.

## Layout

```
hn/             the library: exact fields, geometry, graphs, SAT colouring,
                certificates, local (adelic) colourings
tests/          pytest suite (about 590 tests; one slow de Grey solve)
scripts/        maintained tools: verification, growth, gates, field screens
scripts/experiments/
                726 one-off exploratory scripts, kept as a record
data/           graphs and witnesses with exact coordinates (JSON)
certificates/   colourings, DRAT verification logs, non-colourability claims
notes/          technical notes: local colourings, rigidity, literature, jobs
docs/research-log.md
                the full chronological log, including corrections
```

## Reproducing

```sh
cd research/hadwiger-nelson
python3 -m pip install -r requirements.txt   # python-sat, numpy, scipy, sympy, pytest
python3 -m pytest -q                         # the full suite takes hours
python3 -m pytest -q tests/test_q311.py      # the ℚ(√3, √11) theorem, seconds
sh scripts/worker_setup.sh                   # kissat and drat-trim, for the searches
```

`scripts/verify_pair.py` rebuilds a witness or gadget from its JSON file,
re-derives every edge exactly and runs the solvers. It should be the first step
in checking any claim made here.

## How this work was produced

This is AI-assisted research. The code, the experiments and most of the
writing were produced with Claude (Anthropic), through Claude Code, under the
direction of the repository owner.

Every computational claim is machine-checked as described above. The
mathematical arguments are backed by tests wherever that is possible. **None
has been refereed or independently checked by a mathematician.** Treat them as
preprint-level claims. Corrections are welcome.

## References

- A. D. N. J. de Grey, *The chromatic number of the plane is at least 5*,
  Geombinatorics 28 (2018); [arXiv:1804.02385](https://arxiv.org/abs/1804.02385)
- G. Exoo, D. Ismailescu, *The chromatic number of the plane is at least 5: a
  new proof*, DCG 64 (2020); [arXiv:1805.00157](https://arxiv.org/abs/1805.00157)
- G. Exoo, D. Ismailescu, *A 6-chromatic two-distance graph in the plane*;
  [arXiv:1909.13177](https://arxiv.org/abs/1909.13177)
- M. J. H. Heule, *Computing small unit-distance graphs with chromatic number
  5*; [arXiv:1805.12181](https://arxiv.org/abs/1805.12181)
- J. Parts, *The chromatic number of the plane is at least 5: a human-verifiable
  proof*; [arXiv:2010.12661](https://arxiv.org/abs/2010.12661)
- G. E. Moorhouse, *On the chromatic numbers of planes* (draft, 2010);
  [pdf](https://www.ericmoorhouse.org/pub/chromatic.pdf)
- D. A. Madore, *The Hadwiger–Nelson problem over certain fields*;
  [arXiv:1509.07023](https://arxiv.org/abs/1509.07023)
- Polymath16 threads,
  [3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/)
  (Speyer's 2-adic colourings of the Moser ring) and
  [17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/)
  (Voronov's conjecture)
- Á. Dúcz, *A note on geometric colorings of the Moser lattice*;
  [arXiv:2606.12325](https://arxiv.org/abs/2606.12325)

`notes/literature.md` compares the project with the literature in detail.
