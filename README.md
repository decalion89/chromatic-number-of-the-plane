# The Hadwiger–Nelson problem over number fields

**Sergi Galán** · research repository, 2026

[![tests](https://github.com/decalion89/chromatic-number-of-the-plane/actions/workflows/tests.yml/badge.svg)](https://github.com/decalion89/chromatic-number-of-the-plane/actions/workflows/tests.yml)
[![lean](https://github.com/decalion89/chromatic-number-of-the-plane/actions/workflows/lean.yml/badge.svg)](https://github.com/decalion89/chromatic-number-of-the-plane/actions/workflows/lean.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22976635.svg)](https://doi.org/10.5281/zenodo.22976635)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![peer review: not yet](https://img.shields.io/badge/peer%20review-not%20yet-lightgrey.svg)](#how-this-work-was-done)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)

The Hadwiger–Nelson problem asks for the chromatic number χ(ℝ²) of the plane: the least number of
colours such that no two points at distance exactly 1 share a colour. It has been known since 2018 that
χ(ℝ²) is 5, 6 or 7.

| bound | who, when | how |
|---|---|---|
| χ(ℝ²) ≥ 4 | Nelson, 1950; L. and W. Moser, 1961 | the 7-vertex Moser spindle |
| χ(ℝ²) ≥ 5 | de Grey, 2018 | a 1581-vertex unit-distance graph with no 4-colouring |
| χ(ℝ²) ≥ 5 | Parts, 2020 | the same property on 509 vertices |
| χ(ℝ²) ≤ 7 | Isbell, 1950 | a hexagonal tiling |

This repository studies the problem through exact arithmetic in number fields:
- every point has algebraic coordinates, so "distance 1" is decided exactly;
- claims that a graph cannot be coloured are decided by SAT solvers, and the main ones are certified by
  DRAT proofs checked by an independent program (`drat-trim`);
- the two theorems below are also proved in Lean 4.

**Status (September 2026):** AI-assisted research; nothing here has been refereed yet. Each result states
its evidence, with the statuses of [Results](#results).

## Main result

> **Theorem.** χ(ℚ(√2, √3)²) = 4. The same argument gives a short proof of
> K. G. Fischer's theorem χ(ℚ(√3, √11)²) = 4 (1994).

Both planes can be coloured with four colours, and both contain unit-distance graphs that need four.
ℚ(√3, √11) is the smallest field whose plane contains a Moser spindle (Moorhouse, 2010).

- **Background.** Fischer proved χ(ℚ(√3, √11)²) = 4 in 1994 (*Congressus Numerantium* 104). His result
  seems to have been overlooked: Moorhouse (2010), Madore (2015, who proved 4 ≤ χ ≤ 5), Cranston–Rabern
  (2017), Exoo–Ismailescu (2020) and Voronov (Polymath16, 2021) all treat the value as unknown. Voronov
  also thought χ = 4 likely for the plane over ℚ(√2, √3). Fischer's hypotheses, as the zbMATH review of
  his 1994 paper states them (we have not seen its full text), exclude that field, and we have not found
  its value in the literature.
- **Proof idea.** Change coordinates to α = x + y/√3, β = 2y/√3. The squared distance becomes
  α² − αβ + β². This form is anisotropic modulo a prime above 2 with residue field 𝔽₂, so Madore's
  reduction argument, which he stated for any quadratic form, applies. Reducing (α, β) modulo that prime,
  after subtracting fixed coset representatives, gives a proper 4-colouring of the whole plane. Speyer
  used reduction modulo 2 in 2018 to 4-colour the Moser ring; on that ring our colouring is one of his.
- **Status.** Proved, and formally verified in Lean 4 with Mathlib: both theorems depend only on Lean's
  three standard axioms. The colourings were also tested by computer on finite unit-distance graphs with
  up to 3 134 vertices. Not yet refereed.

**Read:** [the note (PDF)](papers/planes-4-chromatic/planes-4-chromatic.pdf) ·
[the Lean proofs](lean/README.md) ·
[full details](notes/local_colourings.md) (§8 and §10) ·
[comparison with the literature](notes/literature.md)

**Check:**

```sh
python3 -m pytest -q tests/test_q23.py tests/test_q311.py                        # seconds
cd lean && lake exe cache get && lake build && lake env lean PrintAxioms.lean   # minutes; needs elan
```

<details>
<summary><b>What was known, and what is new</b></summary>

According to the zbMATH review of his 1994 paper (*A planar geometric graph of chromatic number four*,
Congr. Numer. 104, 73–79), Fischer proved that ℚ(√p, √q)² has an additive 4-colouring for squarefree,
relatively prime p ≡ 3, q ≡ 11 (mod 16) with pq ≡ 1 (mod 32); for ℚ(√3, √11) the Moser spindle gives the
lower bound. In later work:
- Moorhouse (2010) left the value undetermined.
- Madore (2015) proved 4 ≤ χ ≤ 5.
- Exoo and Ismailescu (2018) asked whether a 5-chromatic unit-distance graph
  embeds in this plane.
- Cranston and Rabern (Combinatorica, 2017) asked for its fractional and its
  ordinary chromatic number.
- Voronov (Polymath16, 2021) wrote that χ = 4 "seems likely" for this plane and
  for the plane over ℚ(√2, √3), "[b]ut as far as I know, nobody has proved this
  yet".

We found Fischer's paper only after the first version of the note had been sent to two mathematicians;
the note now credits it.

The proof reduces z = x + iy modulo a place of ℚ(√3, √11) above 2 (there are two), which is inert in
ℚ(i, √3, √11). Every unit vector becomes a nonzero element of 𝔽₄, so the residue of z − rep(z), where
rep(z) is a fixed representative of the coset of z modulo the local ring, is a proper 4-colouring. In the
coordinates α = x + y/√3, β = 2y/√3 the proof is Madore's reduction argument (Prop. 3.2, which his ¶6.6
states for any quadratic form); in the coordinates (x, y) that argument fails at 2.

Fischer's colouring is additive, with values in ℤ/4. Speyer used reduction modulo 2 in Polymath16
(thread 2, April 2018) to 4-colour the Moser ring. The passage from a subgroup or a ring to the whole
field by cosets is Fischer's (1990, Thm 1), Moorhouse's and Madore's. Our contribution is the choice of
coordinates, which lets Madore's argument work at the places over 2 (they are inert in L(i), so the
reduction covers every unit vector of the plane), and the case ℚ(√2, √3). None of this has been refereed.

</details>

<p align="center">
  <img src="docs/figures/plane_q311.svg" width="560"
       alt="A piece of the plane over Q(sqrt3, sqrt11), 4-coloured by the 2-adic colouring">
</p>
<p align="center"><sub><b>Figure 1.</b> 163 points of the plane over ℚ(√3, √11) and their 594
unit-distance edges. Each point is coloured by the residue pair (ρ(α), ρ(β)) ∈ 𝔽₂²;
no edge joins two points of the same colour.</sub></p>

<p align="center">
  <img src="docs/figures/lower_bounds.svg" width="600"
       alt="The Moser spindle and the 10-vertex chain of unit rhombi, each needing four colours">
</p>
<p align="center"><sub><b>Figure 2.</b> The two lower bounds: the Moser spindle over ℚ(√3, √11) and a
chain of three unit rhombi over ℚ(√2, √3). In a 3-colouring the dashed edge would join two
points of the same colour.</sub></p>

## Results

None of these results has been refereed. Each has one or more of these statuses:
- **Proved**: a written proof, in a paper or the notes, backed by tests where possible.
- **Formally verified**: proved in Lean 4 with Mathlib; CI checks the axioms the proof uses.
- **Computer proof**: a finite computation whose answer comes with a certificate, checked in exact
  arithmetic or by an independent program: DRAT proofs checked by `drat-trim`, dual certificates of
  semidefinite programmes, exact linear-programming certificates.
- **Known**: an earlier result, credited; the repository reproduces it or proves it again.
- **Solvers only**: SAT solvers agree, and no certificate is stored. The notes and the research log mark
  such claims; none is listed here.

| result | status | evidence |
|---|---|---|
| **χ(ℚ(√2, √3)²) = 4.** Voronov's second case; not found in the literature. | Proved; formally verified | [The note](papers/planes-4-chromatic/planes-4-chromatic.pdf), [`lean/Q23.lean`](lean/Q23.lean), [`notes/local_colourings.md`](notes/local_colourings.md) §10, `tests/test_q23.py`. The lower bound was already implicit in Voronov–Neopryatnaya–Dergachev; a 10-vertex chain of unit rhombi gives a short one (`certificates/chain23_no3coloring.json`). |
| **χ(ℚ(√3, √11)²) = 4**, K. G. Fischer's theorem (1994) | Known; a new short proof, proved and formally verified | [The note](papers/planes-4-chromatic/planes-4-chromatic.pdf), [`lean/Q311.lean`](lean/Q311.lean), [`notes/local_colourings.md`](notes/local_colourings.md) §8, `hn/adelic.py`, `tests/test_q311.py` |
| χ(ℚ(√−3, √−11)) = 4 and χ(ℚ(√−3, √−11, √−247)) = 5, for the whole complex fields | Known; new local proofs | The first follows from Fischer's theorem, the second from Madore's reduction at 11 and Exoo–Ismailescu's graph. Our proofs reduce at the primes 2 and 11 ([`notes/local_colourings.md`](notes/local_colourings.md) §3, [`notes/rigidity.md`](notes/rigidity.md)); the lower bound graph `five_247_c` has a DRAT proof (`certificates/five_247_c_no4coloring.json`). |
| Necessary local conditions for a field to hold a 6-chromatic unit-distance graph | Proved | [`notes/local_colourings.md`](notes/local_colourings.md) §5–§9, `scripts/fieldscreen.py` |
| **Six colours for finite planes.** χ(G_q) ≥ 6 for every prime q ≥ 29 except 31, and χ(𝔽_q²) ≥ 6 for q = 37, 41, 43, 47, 59 and every prime q ≥ 67; χ(𝔽₄₁²) ∈ {6, 7}. Here 𝔽_q² is the plane x² + y² (Moorhouse's table stops at q = 17), and G_q the anisotropic plane, a local plane of number fields. | Computer proof; proved for large q | For 𝔽₃₇², 𝔽₄₁², 𝔽₄₃², 𝔽₄₇², G₂₉, G₃₇ and G₄₁, Schrijver's three-point bound gives α < q²/5, and each dual certificate is checked in interval and exact rational arithmetic, and again by an independent checker. The other cases follow from Proposition B and Hoffman's bound, with Weil's estimate or the exact spectrum. [`notes/local_colourings.md`](notes/local_colourings.md) §14, [`data/threepoint/`](data/threepoint/README.md), `scripts/threepoint_verify.py`, `scripts/threepoint_verify_indep.py` |
| χ(ℝ²) ≥ 4: the Moser spindle has no 3-colouring | Known; computer proof | L. and W. Moser (1961). `certificates/moser_spindle_no3coloring.json`, checked by drat-trim. |
| χ(ℝ²) ≥ 5: de Grey's 1581-vertex graph, rebuilt from his 39-point seed, has no 4-colouring | Known; computer proof | De Grey (2018). A DRAT proof of 13.1 M lemmas, checked by drat-trim, for the formula with the colours of one triangle fixed, which loses no generality: `certificates/degrey_1581_no4coloring.json`. |
| Two 5-chromatic unit-distance graphs in ℚ(√3, √11, √247): `five_247_c`, 803 vertices and vertex-critical, and `five_247`, 1 139 vertices. Not a record: Parts' 509 stands. | Computer proof | DRAT proofs checked by drat-trim: `certificates/five_247_c_no4coloring.json`, `certificates/data_no4_checks.txt`; `tests/test_five_247.py` |
| A vertex-critical unit-distance graph with 19 vertices and 33 edges and no 3-colouring. Unlike the Moser spindle, its obstruction combines two constraints, neither of which is forced on its own. | Computer proof | A DRAT proof checked by drat-trim: `certificates/genuine_pair_19_no3coloring.json` |
| Two multi-distance graphs with no 5-colouring: 187 points with edges at 1 and at one Galois orbit of two distances, and 72 points with edges at 1, 4/√3, √7, √(28/3), √(61/3) | Computer proof | Kissat and drat-trim; four solvers agree. `data/W_moser_orbit_9_33.json`, `data/W_lattice_16_21_28_61.json`, their drat-trim logs in [`certificates/`](certificates/README.md), and the research log |
| No twisted colouring of the module of `five_rho7` is proper | Computer proof | 3 840 exact linear-programming (Stiemke) certificates, one for each of the 960 × 4 pairs (ψ, t): [`notes/rigidity.md`](notes/rigidity.md) §2, `scripts/stiemke.py` |

**In progress.** G₁₇, the anisotropic plane over 𝔽₁₇, is a local plane of ℚ(√−3, √−7, √−11); χ(G₁₇) is 5
or 6, and α(G₁₇) ≤ 57 would make it 6. The first part of that bound is checked
(`scripts/g17_alpha.py`, `certificates/g17_part_a_checks.txt`); the certification of the second part is
running.

## Towards χ(ℝ²) ≥ 6

**χ(ℝ²) ≥ 6 has not been proved.** It needs one finite unit-distance graph with no proper 5-colouring. A
finite object can be searched for, and anyone can check it once found. That graph has not been found.

The route searched in September 2026 is the reduction of Exoo and Ismailescu, which Polymath16 calls
clamping onto "virtual edges":
1. a **witness**: a graph with edges at 1 and at some distances `d` that has
   no 5-colouring;
2. a **gadget** for each `d`: a unit-distance graph in which two points at
   distance `d` always get different colours.

The 187-point graph of [Results](#results) is such a witness: it would prove χ(ℝ²) ≥ 6 if a gadget
existed for its orbit of distances.

Two observations that guided the search:
- **Repulsive distances.** Some distances, such as 2/√3, are coloured alike
  unusually rarely in 5-colourings, so they are the natural gadget targets.
  We have not found this measurement in the literature. For four colours,
  Ismailescu found a 103-vertex graph with edges at 1 and 2/√3 and no
  4-colouring (Polymath16, thread 3).
- **Galois orbits.** A Galois automorphism that preserves unit distance maps
  gadgets to gadgets, so one gadget serves a whole orbit of distances. The
  principle is Tao's (Polymath16, thread 7, comment 4893). The 187-point
  witness above needs a single gadget.

The gadget searches ran as parallel jobs from 23 September 2026, described in
[`notes/worker_jobs.md`](notes/worker_jobs.md). None had succeeded by 25 September.

**What is closed, and why.** The [research log](docs/research-log.md) records the approaches that cannot
reach six, with the reason for each: symmetrised growth; spindling at five colours; coset and circular
colourings as obstructions; fields whose local planes are 5-colourable. See its sections "What has been
ruled out so far", "The search over operations is closed, by a theorem" and "Why every known
construction stops at five".

## How claims are checked

- **Exact arithmetic.** Coordinates live in number fields such as ℚ(√3, √11), never in floating point.
  "Distance exactly 1" is a decidable predicate.
- **Solvers, then certificates.** Colourability is decided by SAT solvers (kissat, CaDiCaL, Glucose,
  MiniSat), on formulas rebuilt from the coordinates. Colourings are checked against edges rebuilt from
  the coordinates. A claim that a graph has no proper colouring, when a result rests on it, comes with a
  DRAT proof checked by `drat-trim`; any symmetry breaking in its formula is stated and justified. Smaller
  claims that rest on solvers alone are marked as such ("SAT", or the solvers named).
- **Semidefinite and spectral bounds.** A lower bound from the three-point bound (the finite planes) needs
  a stored dual certificate, checked by a program independent of the solver, in interval arithmetic with
  an exact rational positive-definiteness test, and again by a second checker written from the
  definitions. Spectral (Hoffman) bounds need no certificate: `scripts/finite_hoffman.py` recomputes every
  eigenvalue in interval arithmetic.
- **Formal proofs.** The two theorems on ℚ(√2, √3) and ℚ(√3, √11) are also proved in Lean 4
  ([`lean/`](lean/README.md)). CI builds the proofs, checks that they use only Lean's standard axioms,
  and replays them in Lean's kernel.
- **Corrections stay visible.** Withdrawn claims are kept, with the reason, in the research log.

## Reproducing

```sh
git clone https://github.com/decalion89/chromatic-number-of-the-plane
cd chromatic-number-of-the-plane
python3 -m pip install -r requirements.txt   # python-sat, numpy, scipy, sympy, python-flint, pytest
python3 -m pytest -q                         # the full suite takes hours
sh scripts/worker_setup.sh                   # kissat and drat-trim, for the searches
```

**Checking one result.** From the root of the repository:

| result | command | time |
|---|---|---|
| χ(ℚ(√2, √3)²) = 4 | `python3 -m pytest -q tests/test_q23.py` | seconds |
| χ(ℚ(√3, √11)²) = 4 | `python3 -m pytest -q tests/test_q311.py` | seconds |
| both theorems, formally (needs [elan](https://github.com/leanprover/elan)) | `cd lean && lake exe cache get && lake build && lake env lean PrintAxioms.lean` | minutes |
| four and five colours suffice for the fields ℚ(√−3, √−11) and ℚ(√−3, √−11, √−247) | `python3 -m pytest -q tests/test_moser_field.py tests/test_reduce11.py` | seconds |
| six colours for a finite plane, e.g. 𝔽₄₇²: it prints the rigorous bound α ≤ 371.41…, below 47²/5 = 441.8, so χ ≥ 6 | `python3 scripts/threepoint_verify.py data/threepoint/std47.npz` | 3–5 minutes |
| all eight three-point certificates | `python3 -m pytest -q tests/test_threepoint_certificates.py` | 15 minutes |
| the same, with an independent checker | `python3 -m pytest -q tests/test_threepoint_indep.py` | 10 minutes |
| spectral bounds for large q | `python3 scripts/finite_hoffman.py 59 71` and `python3 scripts/finite_hoffman.py --inert 53 59 61` | seconds |
| de Grey's graph needs five colours | `python3 -m pytest -q tests/test_degrey.py` | up to four hours |

`scripts/verify_pair.py` rebuilds a unit-distance graph or gadget from its data file, recomputes every
edge exactly and runs the solvers; the multi-distance witnesses have their own checkers, listed in
[`scripts/README.md`](scripts/README.md).

The certificates in `data/threepoint/` were produced with Python 3.11.15, numpy 2.4.6, scipy 1.17.1,
cvxopt 1.3.3 (DSDP) and clarabel 0.11.1; checking them needs only numpy, scipy and mpmath 1.3.0
(installed with sympy), and python-flint for the independent checker. `data/threepoint/SHA256SUMS` fixes
their contents. `requirements-lock.txt` lists the exact versions of the Python packages used for the
results and of their dependencies, and `scripts/worker_setup.sh` builds the pinned kissat and drat-trim.

GitHub Actions runs the fast part of the suite, 395 tests in 35 files
([`tests.yml`](.github/workflows/tests.yml)), and builds and checks the Lean proofs
([`lean.yml`](.github/workflows/lean.yml)), on pushes to `main` and on pull requests.

## Repository layout

| path | contents |
|---|---|
| [`papers/`](papers/README.md) | the papers, in LaTeX and PDF |
| [`hn/`](hn/) | the Python library: exact number fields, geometry, unit-distance graphs, SAT colouring, certificates, local (adelic) colourings |
| [`lean/`](lean/README.md) | formal proofs in Lean 4 of the two theorems |
| [`notes/`](notes/README.md) | technical notes: local colourings, rigidity, the literature, the search jobs |
| [`data/`](data/README.md) | graphs and witnesses in exact coordinates (JSON), and the three-point certificates |
| [`certificates/`](certificates/README.md) | colourings, DRAT verification logs and non-colourability claims |
| [`scripts/`](scripts/README.md) | maintained tools: verification, search, figures; `scripts/experiments/` keeps the 733 one-off experiments behind the research log |
| [`tests/`](tests/README.md) | the test suite: 684 tests, 24 of them marked slow |
| [`docs/`](docs/README.md) | the research log, the full chronological record, and the figures |

Each of these folders has a README describing its contents. Until release 1.1.0 the project sat in
`research/hadwiger-nelson/`; version 5 of the note gives its paths in that layout.

## Citing

Cite a tagged release, so that the reader finds the version you read; the changes between releases are in
[`CHANGELOG.md`](CHANGELOG.md). Zenodo archives each release with its own DOI: version 1.0.0 is
[10.5281/zenodo.22976636](https://doi.org/10.5281/zenodo.22976636), version 1.1.0 is
[10.5281/zenodo.22985036](https://doi.org/10.5281/zenodo.22985036), and
[10.5281/zenodo.22976635](https://doi.org/10.5281/zenodo.22976635) always resolves to the latest version.
For the repository as a whole, use GitHub's "Cite this repository" button, which reads
[`CITATION.cff`](CITATION.cff):

```bibtex
@software{galan2026hn,
  author  = {Gal{\'a}n, Sergi},
  title   = {The {H}adwiger--{N}elson problem over number fields},
  version = {1.1.0},
  year    = {2026},
  doi     = {10.5281/zenodo.22985036},
  url     = {https://github.com/decalion89/chromatic-number-of-the-plane/releases/tag/v1.1.0},
  note    = {AI-assisted research; not peer reviewed}
}
```

For the note on the two 4-chromatic planes:

```bibtex
@misc{galan2026planes,
  author = {Gal{\'a}n, Sergi},
  title  = {A short proof that the planes over {$\mathbb{Q}(\sqrt{3},\sqrt{11})$}
            and {$\mathbb{Q}(\sqrt{2},\sqrt{3})$} are 4-chromatic},
  year   = {2026},
  note   = {Version 5, 27 September 2026. Preprint, not refereed. AI-assisted},
  url    = {https://github.com/decalion89/chromatic-number-of-the-plane}
}
```

Please also cite the original papers listed under [References](#references).

## How this work was done

This is AI-assisted research. The code, the computations and most of the text were produced with Claude
(Anthropic), through Claude Code, under the direction of Sergi Galán. The computational claims are
checked by machine: the main ones by certificates that an independent program verifies, the others by
solvers, as each one states. The mathematical arguments are backed by tests wherever that is possible.
No result has been refereed or independently checked by a mathematician yet: treat them as
preprint-level claims.

## Reporting an error

Corrections are welcome. Please open an issue with one of the two templates: *Mathematical error* (a
statement, proof or table entry that is wrong or unsupported) or *Result does not reproduce* (the command
you ran and what you saw). Corrections are recorded in the research log and in the changelog, not edited
away. [`CONTRIBUTING.md`](CONTRIBUTING.md) says what a pull request needs.

## License

Code, data and text are released under the [MIT License](LICENSE).

## References

- A. D. N. J. de Grey, *The chromatic number of the plane is at least 5*,
  Geombinatorics 28(1) (2018) 18–31; [arXiv:1804.02385](https://arxiv.org/abs/1804.02385)
- G. Exoo, D. Ismailescu, *The chromatic number of the plane is at least 5: a
  new proof*, Discrete Comput. Geom. 64(1) (2020) 216–226;
  [arXiv:1805.00157](https://arxiv.org/abs/1805.00157)
- G. Exoo, D. Ismailescu, *The Hadwiger–Nelson problem with two forbidden
  distances*; [arXiv:1805.06055](https://arxiv.org/abs/1805.06055)
- G. Exoo, D. Ismailescu, *A 6-chromatic two-distance graph in the plane*;
  [arXiv:1909.13177](https://arxiv.org/abs/1909.13177)
- M. J. H. Heule, *Computing small unit-distance graphs with chromatic number
  5*, Geombinatorics 28(1) (2018) 32–50; [arXiv:1805.12181](https://arxiv.org/abs/1805.12181)
- D. W. Cranston, L. Rabern, *The fractional chromatic number of the plane*,
  Combinatorica 37(5) (2017) 837–861; [arXiv:1501.01647](https://arxiv.org/abs/1501.01647)
- J. Parts, *The chromatic number of the plane is at least 5: a human-verifiable
  proof*; [arXiv:2010.12661](https://arxiv.org/abs/2010.12661)
- J. Parts, *Graph minimization, focusing on the example of 5-chromatic
  unit-distance graphs in the plane* (the 509-vertex graph);
  [arXiv:2010.12665](https://arxiv.org/abs/2010.12665)
- K. G. Fischer, *Additive K-colorable extensions of the rational plane*,
  Discrete Math. 82(2) (1990) 181–195; *The connected components of the graph
  ℚ(√N₁, …, √N_d)²*, Congr. Numer. 72 (1990) 213–221 (Zbl 0733.05048); and
  *A planar geometric graph of chromatic number four*, Congr. Numer. 104 (1994)
  73–79 (Zbl 0836.05030)
- P. D. Johnson Jr., *Two-colorings of real quadratic extensions of ℚ² that
  forbid many distances*, Congr. Numer. 60 (1987) 51–58
- P. D. Johnson Jr., *Problems posed in or arising from "Colorings of metric
  spaces": status report*, Geombinatorics 9 (2000) 170–179 (a survey we have
  not seen)
- M. S. Payne, *Unit distance graphs with ambiguous chromatic number*,
  Electron. J. Combin. 16 (2009), Note 31;
  [arXiv:0707.1177](https://arxiv.org/abs/0707.1177) (summarises Johnson's and
  Fischer's results on quadratic fields)
- G. E. Moorhouse, *On the chromatic numbers of planes* (draft, 2010);
  [pdf](https://www.ericmoorhouse.org/pub/chromatic.pdf)
- D. A. Madore, *The Hadwiger–Nelson problem over certain fields*;
  [arXiv:1509.07023](https://arxiv.org/abs/1509.07023)
- A. Medrano, P. Myers, H. M. Stark, A. Terras, *Finite analogues of Euclidean
  space*, J. Comput. Appl. Math. 68 (1996) 221–238
  ([doi](https://doi.org/10.1016/0377-0427(95)00261-8))
- A. Schrijver, *New code upper bounds from the Terwilliger algebra and
  semidefinite programming*, IEEE Trans. Inform. Theory 51 (2005) 2859–2866
- Polymath16 threads
  [2](https://dustingmixon.wordpress.com/2018/04/22/polymath16-second-thread-what-does-it-take-to-be-5-chromatic/)
  and [3](https://dustingmixon.wordpress.com/2018/05/01/polymath16-third-thread-is-6-chromatic-within-reach/)
  (Speyer's 2-adic colourings of the Moser ring),
  [7](https://dustingmixon.wordpress.com/2018/06/16/polymath16-seventh-thread-upper-bounds/)
  (Tao on Galois symmetry) and
  [17](https://dustingmixon.wordpress.com/2021/02/01/polymath16-seventeenth-thread-declaring-victory/)
  (Voronov's questions), and the Polymath16 wiki page
  [Algebraic formulation of Hadwiger–Nelson problem](https://web.archive.org/web/20210412075722/https://asone.ai/polymath/index.php?title=Algebraic_formulation_of_Hadwiger-Nelson_problem)
  (colourings of rings such as the Moser ring)
- V. A. Voronov, A. M. Neopryatnaya, E. A. Dergachev, *Constructing
  5-chromatic unit distance graphs embedded in the Euclidean plane and
  two-dimensional spheres*, Discrete Math. 345(12) (2022) 113106;
  [arXiv:2106.11824](https://arxiv.org/abs/2106.11824)
- Á. Dúcz, *A note on geometric colorings of the Moser lattice*;
  [arXiv:2606.12325](https://arxiv.org/abs/2606.12325)

[`notes/literature.md`](notes/literature.md) compares the project with the literature in detail.

## Resumen en español

Este repositorio estudia el problema de Hadwiger–Nelson: el número cromático del plano, que se sabe que
está entre 5 y 7. Se trabaja con aritmética exacta en cuerpos de números. Las afirmaciones de que un
grafo no se puede colorear se deciden con resolutores SAT, y las principales van acompañadas de una
prueba DRAT verificada por un programa independiente.

Resultado principal: el plano con coordenadas en ℚ(√2, √3) tiene número cromático exactamente 4, como
Voronov consideraba probable en 2021. El mismo argumento da una prueba corta de un teorema de
K. G. Fischer (1994): el plano sobre ℚ(√3, √11) también tiene número cromático 4. El resultado de
Fischer había pasado desapercibido; los trabajos posteriores lo daban por abierto. La prueba cambia de
coordenadas para que el argumento de reducción de Madore funcione módulo 2. Está explicada en una
[nota breve](papers/planes-4-chromatic/planes-4-chromatic.pdf) y verificada formalmente en Lean 4.
También se demuestra, con certificados verificados por dos programas independientes, que siete planos
finitos necesitan seis colores. Es un trabajo hecho con ayuda de IA y todavía no ha sido revisado por
pares.
