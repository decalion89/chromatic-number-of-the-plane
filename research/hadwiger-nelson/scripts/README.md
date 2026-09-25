# Scripts: the maintained tools

The 85 scripts here are the project's maintained tools; the 728 one-off experiments behind the
research log are in [`experiments/`](experiments/). Run every tool from `research/hadwiger-nelson/`:
each finds the `hn` package from its own location.

**Start here.**
- `verify_pair.py` rebuilds a claimed obstruction from its JSON file, rederives every edge in exact
  arithmetic and asks three pysat solvers, and kissat with drat-trim when their paths are given,
  whether the graph is 5-colourable.
- `grow_lean.py` is the colouring-guided growth used in the search for a 6-chromatic graph.
- `fieldscreen.py` and `module_gate.py` look for a colouring of a whole field or module through a
  finite residue plane; when one with at most five colours exists, no graph in that field or module
  is 6-chromatic.
- `worker_setup.sh` installs the Python packages, kissat, drat-trim and `tabu2` on a fresh machine.

The C programs `tabucol` and `tabu2` are compiled from `tabucol.c` and `tabu2.c` with `gcc -O2`;
the binaries are not in the repository. Eight scripts (`asym_grow.py`, `gadget.py`, `measure_fk.py`,
`search_disjunction.py`, `search_forced.py`, `tabuSaP.py`, `verifygate.py`, `verifyquot.py`) refer to
a temporary directory of the session that wrote them (`/tmp/claude-0/…/scratchpad`), as a default
output directory (`HN_OUT` overrides it), as an import path, or as the location of an input file.

## Verification

| script | what it does |
|---|---|
| `verify_pair.py` | Rebuilds a graph from its JSON file, rederives every edge in exact arithmetic and asks three pysat solvers (and kissat with drat-trim, when given) whether its unit-distance graph, or its graph with edges at 1 and at one rational distance d, is 5-colourable, optionally testing whether a target pair A, B is forced apart or alike; it cannot read `data/W_moser_orbit_9_33.json` (irrational extra distances) or `data/W_lattice_16_21_28_61.json` (lattice coordinates, no `field_generators`), which `orbit_witness_test.py` and `lattice_witness.py` check. |
| `verify6.py` | Rechecks a claimed 6-chromatic unit-distance graph: rederives every edge exactly, confirms that no edge has another length, solves 5-colourability with no colour pinned using more than one solver, and extracts a vertex-critical core. |
| `verify_gadget.py` | Rechecks a claimed distance-2 gadget (a unit-distance graph with vertices a, b at distance 2 that differ in every proper 5-colouring): rederives the edges and \|ab\|² = 4 exactly, merges a and b, and asks three solvers whether the merged graph is 5-colourable. |
| `orbit_witness_test.py` | Tests whether a graph with its unit edges and the edges at one or more Galois orbits of distances (`14_2`, `14_5`, `9_1`, `43`) is 5-colourable, by `tabucol` and then kissat; with orbit `9_1` it checks `data/W_moser_orbit_9_33.json`. |
| `lattice_witness.py` | Builds the graph on the points of ℤ[ω]/√−3 within a radius R whose edges are the pairs at given norms (squared distance N/3), asks three pysat solvers whether it is 5-colourable, shrinks it to a vertex-critical core and checks that core with kissat and drat-trim; `lattice_witness.py 3,16,21,28,61 5 out.json` reproduces `data/W_lattice_16_21_28_61.json`. |
| `verify_forced_same.py` | Rechecks from exact coordinates that a pair A, B at distance 5 is alike in every proper 5-colouring (two solvers), and that the union of the graph with its image under the rotation λ = (49 + 3√−11)/50 about A is not 5-colourable. |
| `verify_twisted.py` | Constructs a twisted coset colouring c(p) = ψ(p) + t·⌊(φ(p) − φ₀)/L⌋ on a saved graph and checks exactly that it is proper and colours the target pair alike. |
| `verifyfield.py` | Checks the arithmetic of ℚ(m, √−3, √−11), m a root of x³ − 10x² + 26x − 11, a field with both a Moser spindle and blocking: V(m) is totally positive and equals 33s² with s in ℚ(m), and 5 is inert in ℚ(m). |
| `verifygate.py` | Cross-checks a blocking verdict: the module has full rank modulo 5, random functionals vanish on some unit vector as predicted, and Glucose repeats CaDiCaL's refutation; it reads its input `gate.pkl` from the temporary directory, and that file is not in the repository. |
| `verifyquot.py` | Checks a periodic colouring through a map to ℤ/5 × ℤ/5 on actual points of the module in ℚ(√3, √5, √7, √11), with the unit edges found in exact arithmetic, and counts monochromatic edges. |
| `circverify.py` | Checks in exact rational arithmetic a circular colouring c(x) = ⌊5·frac(φ(x))⌋, with φ given on the basis of `circgate.py`: frac(φ(u)) lies in [1/5, 4/5] for every unit vector u, and a checkpoint graph has no monochromatic edge. |
| `verify229.py` | Rechecks with three solvers, and with the pairs rederived from exact geometry, that de Grey's graph G has no proper 5-colouring once 229 pairs at squared distances 15/16, 16, 17/2, 9 and 7 are added as edges, and that each of the five classes is needed. |
| `verify4.py` | Checks, without solver budgets, whether a sample of the 1548 pairs of de Grey's `Sa` that differed in all 24 sampled 4-colourings are forced to differ, after control runs on a unit rhombus and on K₄. |
| `verify49.py` | Reruns the tabu search that 5-coloured de Grey's graph G with every pair at squared distance 4/9 bichromatic, and checks that colouring against the exact geometry. |
| `verify5.py` | Checks with a newly built formula that a set of five vertices of `Sa` ∪ rot(`Sa`) (619 vertices) is rainbow-forcing at four colours, that is, every proper 4-colouring uses all four colours on it, that the set is minimal, and that the union is 4-colourable. |
| `verify7.py` | Checks with a direct formula that a set of seven vertices of de Grey's `Sa`, which induces only five edges, is rainbow-forcing at four colours, and that it is minimal. |
| `verifySa.py` | Checks against the exact geometry a 5-colouring of de Grey's `Sa`, found by Glucose, in which every pair at squared distance 4/9, 16/9, 4 or 16 is bichromatic. |
| `ei_rebuild.py` | Rebuilds Exoo and Ismailescu's graphs G, H and K (arXiv:1909.13177) in ℚ(√3, √11) and rechecks their sizes, the pair of H at distance 5 that is alike in every 5-colouring, and that K with edges at distances 1 and 2 has no proper 5-colouring. |

## Growth and search for obstructions

| script | what it does |
|---|---|
| `grow_lean.py` | The growth loop of `grow_kw.py` with the candidate pool held in numpy arrays (about ten times less memory) and with every unit-distance pair found, including pairs along unit vectors not yet known; modes `plain`, `apart` and `same`; the number of colours is `K` (5 by default). |
| `grow_kw.py` | The growth loop of `grow_ls2.py` extended to seeds with extra-distance pairs (`two_edges`): it colours with `tabu2`, keeping unit edges proper while minimising the alike extra pairs, and grows near the pairs that stay alike. |
| `grow_ls2.py` | Colouring-guided growth: repairs a 5-colouring of the enlarged graph with `tabucol` (then an incremental CDCL solver or kissat), and adds the candidate points whose unit neighbours carry all five colours. |
| `asym_grow.py` | Grows unit-distance graphs by attaching points one at a time at distance 1 from random vertices, without symmetrising, and measures the forced core as it grows. |
| `make_pair_seed.py` | Writes a growth seed for `MODE=apart` or `MODE=same` from a graph file, with the same points, unit vectors and colouring and the target pair A, B set (`make_pair_seed.py graph.json A B out.json`). |
| `mk_kseed.py` | Places Exoo–Ismailescu's graph H, turned by 150°, at a vertex of `five_rho7.json` and writes `data/K_rot.json` (K = H ∪ λ_A(H)) and `data/rho7_lambda_K.json`. |
| `backbone.py` | Looks for any pair alike in every 5-colouring of a grown graph, or any pair at distance 2 that differs in every one, by sampling colourings, filtering pairs with Kempe chains and settling the survivors with one incremental SAT query. |
| `search_forced.py` | Builds a ball of unit-vector walks in ℚ(√3, √11), searches it for a pair at a distance a spindle rotation can close that is monochromatic in every k-colouring (k set by `HN_K`), and spindles that pair. |
| `search_disjunction.py` | Searches for a forced disjunction, a pivot p and a set Q of vertices at one distance such that every k-colouring gives p the colour of some element of Q, and shrinks Q to a minimal forcing subset. |
| `gadget.py` | Enumerates small vertex sets S of a carrier graph (de Grey's `Sa` by default) and asks, one SAT call each, whether every k-colouring uses fewer than k colours on S, with a control run at k = 4 that must rediscover de Grey's gadget. |
| `six.py` | Asks whether de Grey's graph G has a pair, at a distance that a rotation of his field can close, that is monochromatic in every 5-colouring. |
| `measure_fk.py` | Measures f(4), the size of the smallest unit-distance configuration with a pair monochromatic in every 4-colouring, starting from balls around the centre of de Grey's graph. |
| `is63.py` | Decides whether ρ(G, 5) ≤ 63 for de Grey's graph G, ρ being the size of a smallest rainbow-forcing set, by alternating a cardinality-constrained hitting-set query with a search for a colour class that escapes the set. |
| `rhotight.py` | The decision loop of `is63.py`, with the hitting set shrunk to a minimal one before each search for an escaping colour class. |
| `bigfilter.py` | Finds the pairs alike in every sampled colouring by bucketing vertices by their colour signatures, in time linear in the number of vertices, for forced-pair searches on large graphs. |
| `pairsat.py` | Finds every pair at a given squared distance by the cell bucketing of the unit-edge finder, confirming each candidate in exact arithmetic; the experiments import it as a library. |
| `blocked.py` | Searches for a blocked point, a point of the field whose unit neighbours in G use all five colours in every proper 5-colouring of G, with one SAT call per candidate. |

## Colouring by local search

| script | what it does |
|---|---|
| `tabucol.c` | k-colouring by tabu search (Hertz–de Werra, with the incremental tables of Galinier–Hao), reading the graph on standard input. |
| `tabu2.c` | Tabu search for a proper k-colouring on hard edges that then minimises the number of alike soft pairs without breaking a hard edge. |
| `tabucol.py` | A Python TabuCol, calibrated on three instances with known answers (`Sa` at five colours with and without the class 4/9 forbidden, and `Sa` at four colours with 4/9 forbidden) before it is run on open instances. |
| `tabu.py` | An earlier tabu search that recolours a random conflicting vertex, run on `Sa`, Y and G with forbidden pairs; it stalled on a satisfiable calibration instance, and `tabucol.py` replaces it. |
| `tabuSaP.py` | Runs the calibrated TabuCol at five colours on each closable distance class of a grown seed of `Sa`, smallest first. |
| `worker_setup.sh` | Installs the Python packages, builds kissat and drat-trim under `$TOOLS` (default `~/hn-tools`) and compiles `tabu2`; it can be run more than once. |

## Gates: periodic, circular and local colourings that rule a search out

| script | what it does |
|---|---|
| `gate.py` | Tests coset colourings on the module spanned by a graph's edge vectors, in an integer echelon basis of that module: homomorphisms to ℤ/n for n = 2, …, 5 and a periodic screen for n = 6, …, 12. |
| `module_gate.py` | Reduces a graph's module at a place above an unramified odd prime where all its coordinates are integral, and colours the finite image, which bounds the chromatic number of every graph in that module. |
| `residuegate.py` | Embeds the unit vectors of a module in ℚ_{p²} at every place above an unramified prime p and, when they are integral there, colours the Cayley graph of their residues in 𝔽_{p²}, which colours every graph along these unit vectors. |
| `ramgate.py` | The residue gate at a prime p that may ramify, reducing the unit vectors modulo √p in ℚ_p(√r, √p). |
| `circgate.py` | Searches by mixed-integer programming for a circular colouring c(x) = ⌊5·frac(φ(x))⌋ with φ in Hom(M, ℝ/ℤ), maximising the margin by which every frac(φ(u)) lies inside [1/5, 4/5]. |
| `circrel.py` | A relaxation of the circular gate in the values t_u = frac(φ(u)), using short integer relations among the unit vectors; infeasibility shows that no circular colouring exists. |
| `circsat.py` | Finds by CP-SAT every circular k-colouring of a module whose character φ has denominator dividing a given Q. |
| `circextend.py` | Tests whether a circular colouring of a module M₁ extends to a larger module M′, such as its λ-closure, by enumerating the extensions of φ through the Smith form of the inclusion. |
| `cellsample.py` | Samples the polytope of characters around a known circular colouring of M₁ and tests each sample's extensions to M′ with `circextend.py`. |
| `kempegate.py` | Decides exactly, through the index of a sublattice of ker ψ, whether a Kempe swap of a coset colouring can join a pair a, a + 2e or split a pair a, a + 5e. |
| `prodgate.py` | Searches for 5-colourings periodic through non-cyclic quotients ℤ/n₁ × … × ℤ/n_k of a module that refute a distance-2 gadget or a forced pair a, a + 5e along some direction e. |
| `quotgate.py` | Looks for a 5-colouring of the Cayley graph of M/qM on the unit vectors, which gives a qM-periodic colouring of the unit-distance graph on M, and reports which pairs a, a + 2e and a, a + 5e such a colouring settles. |
| `quotient.py` | Looks for periodic 5-colourings, through the quotients (ℤ/m)^r, of the group of unit vectors generated in ℚ(√3, √5, √7, √11) by the 60° rotation and spindle rotations. |
| `idealquot.py` | Looks for a 5-colouring of the Cayley graph of M/αM on the unit vectors, for α = a + bω in ℤ[ω], which gives an αM-periodic colouring of the module that is not a coset colouring when N(α) is prime to 5. |
| `idealquot2.py` | Searches for colourings periodic under ker ψ ∩ αM and asks whether they can split a pair a, a + 5e. |
| `twistquot.py` | Searches for twisted periodic colourings ψ(x) + k(x mod αM) for an ideal α of ℤ[ω]. |
| `torusgate.py` | Searches for a 5-colouring of the plane, periodic under (Pℤ)² and constant on the cells of a grid on the torus, that is proper along the module's unit directions, with a conservative cell test so that any solution colours the whole module. |
| `torecore.py` | Extracts an unsatisfiable core of the 5-colouring problem on the Cayley graph of M/αM and lifts it greedily to M, reporting how many of its edges become unit steps. |
| `wallgate.py` | Searches for one-dimensional colourings c(x) = G(φ(x), ψ(x)), with φ an integer functional and ψ: M → ℤ/5, that refute a distance-2 gadget or a forced pair a, a + 5e. |
| `nearadm.py` | For a blocked module, completes nearly admissible maps ψ (zero exactly on the unit vectors in 5M) with a second coordinate χ: M → ℤ/m and asks whether the Cayley graph of ℤ/5 × ℤ/m on (ψ, χ)(U) is 5-colourable. |
| `screen_lambda.py` | Tests whether Exoo–Ismailescu's rotation λ = (49 + 3√−11)/50, which is not integral at 5, escapes coset colourings, on H, K = H ∪ λ_A(H), the 803-vertex graph and its union with a λ-rotated copy. |
| `stiemke.py` | Seeks exact Stiemke certificates that a module has no twisted colouring: for each admissible ψ and each t, a strictly positive integer relation among the unit vectors of the class D_t. |
| `farkas.py` | Seeks exact Farkas certificates that no twisted colouring keeps a given pair a, a + 2e alike, each expressing −e as a nonnegative combination of the class D_t. |
| `kappa_check.py` | Checks on a module that every admissible ψ satisfies ψ(κu) = ψ(u) whenever u and κu are both unit vectors, where κ = ωρ₇ = (−11 + 5√−3)/14 ≡ 1 mod 5. |
| `cosetfit.py` | Measures how close a grown graph's colouring is to a coset colouring, the best agreement over admissible ψ and relabellings of the colours, and where the mismatched vertices lie. |
| `periodfind.py` | Measures, for many small lattice vectors t, how often a grown colouring satisfies c(x + t) = π(c(x)) for the best relabelling π, to detect a hidden period. |
| `quot2d.py` | Computes the lattice N spanned by the translations that a grown colouring nearly respects, the quotient M/N with the images of the unit vectors, and whether its Cayley graph is 5-colourable on finite tori. |
| `rigid.py` | Colours a graph by tabu search with several seeds and measures the best coset fit R, which is 1 for a coset colouring and about 0.2 for an unstructured one. |

## Fields and unit vectors

| script | what it does |
|---|---|
| `fieldscreen.py` | Screens a multiquadratic CM field K = ℚ(√−d₁, …, √−d_n) by the places of K⁺ that do not split in K, each of which bounds the chromatic number of the unit-distance graph on K: by 3 at a ramified place, and by that of a finite plane otherwise. |
| `fieldtypes.py` | Decides, for each splitting type of 5 in a CM field of given degree, whether blocking is possible (every hyperplane of O/5 meets the norm-one group), which settles every field of that degree at once. |
| `field24.py` | Implements arithmetic in ℚ(ζ₂₁, √−11), the smallest field with unit triangles, blocking and a Moser spindle, and checks the spindle there. |
| `reduce11.py` | Colours unit-distance graphs in ℚ(√−3, √−11, √−247) by reduction at the places above 11 onto the 12-element norm-one group of 𝔽₁₂₁, checking every unit vector and edge exactly. |
| `moser2adic.py` | Checks in exact rational arithmetic the 2-adic 4-colouring of the plane over ℚ(√−3, √−11) on every unit vector and edge of the given files. |
| `msqrt.py` | Computes square roots in a multiquadratic field exactly, by recursion on the generators. |
| `sqrtK.py` | Computes square roots in a multiquadratic field from sign patterns that are characters of the Galois group, accepting a root only when its square is exact. |
| `allunits.py` | Finds every unit vector of an edge module, not only those the graph uses, by Fincke–Pohst enumeration under the trace form. |
| `lamclosure_units.py` | Writes a module file with the unit vectors of the λ-closure M′ = M₁ + λM₁ that need no search (U₁, λU₁ and the Exoo–Ismailescu units 5(1 − λ)e), for `circgate.py` and `circextend.py`. |
| `rotunits.py` | Writes a module file with the unit vectors R·(1, 0) for words R in the rotations ω, σ, λ, ρ₇, κ and optionally τ, for the circular gates. |
| `unitsinqm.py` | Reports, for graphs in `data/`, the rank of the module spanned by the edge vectors, the moduli q for which some unit vector lies in qM, and the content of each unit vector. |
| `tune.py` | Builds unit-distance graphs whose field is a parameter, by composing a forced pair with a rotated copy of its carrier so that it lands at a chosen squared distance, and writes `data/five_tuned_*.json` when the result has no proper 4-colouring. |

## Measurements on balls and growth curves

| script | what it does |
|---|---|
| `ball2.py` | Builds the complete combinatorial 2-ball {O + u + v} of `five_rho7.json`, colours it by tabu search and measures the best coset fit R. |
| `ballr.py` | Builds combinatorial r-balls of `five_rho7.json` in integer coordinates, colours their 5-cores by tabu search, and compares the pairs x + u, x + v with u − v in 5M against the other non-adjacent pairs. |
| `freecurve.py` | Measures free@5, the proportion of vertices with a spare colour, against mean degree, averaged over a dozen diverse colourings. |

## Figures

| script | what it does |
|---|---|
| `make_figures.py` | Draws the README figures in `docs/figures/` from exact data, after asserting that every colouring shown is proper. |
