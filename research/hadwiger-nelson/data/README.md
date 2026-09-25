# Data: graphs and witnesses in exact coordinates

Every graph here is stored with exact coordinates, so every distance can be recomputed exactly. The
table at the end lists each of the 85 JSON files tracked in this directory once. Growth runs also write
transient checkpoints and solver files here (`*_kw<n>.json`, `*_ls<n>.json`, `*_inc.json`,
`*_grow.json`, `blk_plain*.json`, `*.cnf`, `*.kissat`); `.gitignore` excludes them, and they are not
indexed.

## Formats

Four formats occur.

**Graphs with coordinates in a multiquadratic field** (65 files) have the keys:
- `field_generators`: `[a, b, …]`, meaning that the coordinates lie in ℚ(√a, √b, …);
- `points`: a list of points `[x, y]`, each coordinate a list of `[numerator, denominator]` pairs on
  the basis of square-free products of the generators, in the order `hn.field.Field(field_generators)`
  lists them (for `[2, 3]`: 1, √2, √3, √6).

Optional keys:
- `note`, `recipe`, `mechanism`, `construction`, `source`, `seed`, `centre`: how the graph was built;
  `n` and `m`: its numbers of vertices and of unit-distance edges;
- `A`, `B` (and `Bp`, an image of B): a target pair;
- `units`: the unit vectors a growth may use; `units2` and `dist2_exact`: the vectors and squared
  distances of extra edges (`L16_W0915_seed.json`);
- `two_edges`: pairs of vertices at other distances, treated as extra or soft edges;
- `dist2`: the squared distances at which the edges lie (`W_moser_orbit_9_33.json`); `edges`: the
  edge list (`chain23.json`);
- `colouring`, `colouring_prefix`: a 5-colouring found by a search, of all points or of the first
  ones;
- `status`, `mode`: the state of a growth run;
- `target_d2`, `radical`, `orbit`, `orbit_d2`, `forced_orbit`, `orbits`, `free_at_5`, `class`,
  `d2`, `cos`, `sin`, `rotations`, `base`, `spindle`, `basis`, `n_forcing`: parameters and
  measurements of the construction.

**The lattice witness** `W_lattice_16_21_28_61.json` stores each point as an integer pair `[a, b]`,
standing for (a + bω)/√−3 with ω = e^{2πi/3}, so that two points at difference a + bω are at squared
distance (a² − ab + b²)/3. `forbidden_norms` and `forbidden_d2` give the norms a² − ab + b² and the
squared distances of the edges, and `edges` lists them.

**Graphs with a stated claim**, `tight_four.json` and `witness_five.json`, give the field as
`field` (generators 3, 5, 7, 11 and the basis) and write coordinates as decimal-rational strings on
that basis: a point of `tight_four.json` is one list of 32 strings, 16 for x and then 16 for y, and
a point of `witness_five.json` is `[x, y]` with 16 strings each. `claim` states the property, and
`unit_edges` (with `forbidden_pairs` and `classes` in `witness_five.json`) gives the edges;
`tight_four.json` also records `chromatic_number`, `circular_chromatic_number` and `structure`.

**Results without coordinates** (16 files) record computations on the graph named in a `graph` key
or in the table below, by the vertex indices of that graph; `mu5_menu.json` is a list of such
records.

The tests and `scripts/verify_pair.py` read the first format. `scripts/orbit_witness_test.py` checks
`W_moser_orbit_9_33.json`, and `scripts/lattice_witness.py` rebuilds `W_lattice_16_21_28_61.json`.

## Naming

| prefix | meaning |
|---|---|
| `five_*` | unit-distance graphs with no proper 4-colouring, and results about them (`five_247_c_critical.json`) |
| `W_*` | witnesses: graphs with edges at several distances and no proper 5-colouring |
| `L16_*` | the searches for a 6-chromatic graph in the field L16 = ℚ(√−3, √−7, √−11, √−247), whose points have coordinates in ℚ(√3, √7, √11, √247) |
| `L16_g<d>_<A>_<B>_seed` | the L16 seed with the target pair (A, B) at a distance in a Galois orbit: `0737`, `0915`, `1568` and `2915` stand for d ≈ 0.7366, 0.9149, 1.5676 and 2.9149 |
| `ei_*`, `K_rot` | Exoo–Ismailescu's graphs, rebuilt |
| `rho7_*` | graphs containing `five_rho7.json` |
| `g2e_*` | searches for a distance-2 gadget around `five_247_c.json` |
| `*_seed`, `*_ckpt` | the start and a checkpoint of a growth search |
| `*_five_247_c.json`, `*_tight_hexagon_4159.json` | results about that graph |

## Key files

| file | what it is |
|---|---|
| `chain23.json` | 10 vertices with no proper 3-colouring: the lower bound in χ(ℚ(√2, √3)²) = 4 (`notes/local_colourings.md` §10). |
| `chain35.json` | 19 vertices with no proper 3-colouring: the lower bound in 4 ≤ χ(ℚ(√3, √5)²) ≤ 5 (`notes/local_colourings.md` §11). |
| `W_moser_orbit_9_33.json` | 187 points with no proper 5-colouring when edges are placed at distance 1 and at the Galois orbit d² = (9 ∓ √33)/6: a witness that needs one gadget. |
| `W_lattice_16_21_28_61.json` | 72 lattice points with no proper 5-colouring when edges are placed at 1, 4/√3, √7, √(28/3) and √(61/3). |
| `five_247.json` | 1 139 vertices, 5-chromatic, in ℚ(√3, √11, √247). |
| `five_247_c.json` | 803 vertices, 5-chromatic and vertex-critical, in ℚ(√3, √11, √247). |
| `five_tuned_16_1_3_7_11.json` | 4 081 vertices with no proper 4-colouring, in ℚ(√3, √7, √11): neither √5 nor √247 is needed. |
| `five_rho7.json` | `five_247_c.json` with its images under ρ₇ and ρ₇⁻¹: the module studied in `notes/rigidity.md`. |
| `ei_H214.json`, `ei_rho7.json`, `K_rot.json` | Exoo–Ismailescu's graphs, rebuilt from their paper and turned onto the project's modules. |
| `L16_seed.json` | 6 080 points: the starting graph of the searches in L16. |
| `L16_kw2.json` | 18 524 points: the largest growth from `L16_seed.json` kept. |
| `tight_four.json`, `witness_five.json` | Small graphs with a stated claim: circular chromatic number 4, and 161 pairs one of which is monochromatic in every proper 5-colouring. |

## All files

| file | points | field | description |
|---|---:|---|---|
| `K_rot.json` | 426 | ℚ(√3, √11, √247) | Exoo–Ismailescu's graph K = H ∪ λ_A(H), with H turned by 150° onto the unit vectors of `five_rho7.json` and A at its densest vertex; `two_edges` lists its 892 pairs at distance 2 (written by `scripts/mk_kseed.py`). |
| `L16_W0915_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with extra edges at the Galois orbit d² = (14 ∓ 2√33)/3 (d ≈ 0.9149, 2.9149), given by `units2` and `dist2_exact`, a seed for plain growth towards a witness that, by Galois conjugation, would need a single unit-distance gadget, at d ≈ 0.9149. |
| `L16_apart43_rank.json` | — | — | 200 pairs of `L16_seed.json` at distance 2/√3, each with the count 12 over 12 sampled 5-colourings, ranked as targets for `MODE=apart` growth. |
| `L16_apart43_survivors.json` | — | — | 1 600 pairs of `L16_seed.json` at distance 2/√3 kept as candidates for being forced apart, and the 12 of them chosen as targets for `MODE=apart` growth. |
| `L16_g0737_1021_1317_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 1021, B = 1317 at d² = (9 − √33)/6 (d ≈ 0.7366), in the Galois orbit of `W_moser_orbit_9_33.json`, for `MODE=apart` growth towards a gadget. |
| `L16_g0737_28_868_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | As `L16_g0737_1021_1317_seed.json`, with the target pair A = 28, B = 868. |
| `L16_g0915_32_141_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 32, B = 141 at d² = (14 − 2√33)/3 (d ≈ 0.9149), split in all 16 sampled tabu colourings, for `MODE=apart` growth towards a gadget. |
| `L16_g0915_36_628_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | As `L16_g0915_32_141_seed.json`, with the target pair A = 36, B = 628. |
| `L16_g1568_1048_1620_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 1048, B = 1620 at d² = (9 + √33)/6 (d ≈ 1.5676), in the Galois orbit of `W_moser_orbit_9_33.json`, for `MODE=apart` growth towards a gadget. |
| `L16_g1568_28_600_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | As `L16_g1568_1048_1620_seed.json`, with the target pair A = 28, B = 600. |
| `L16_g2915_2437_3218_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 2437, B = 3218 at d² = (14 + 2√33)/3 (d ≈ 2.9149), split in all 16 sampled tabu colourings, for `MODE=apart` growth towards a gadget. |
| `L16_g2915_41_715_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | As `L16_g2915_2437_3218_seed.json`, with the target pair A = 41, B = 715. |
| `L16_gadget43_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 9, B = 726 at distance 2/√3, for `MODE=apart` growth towards a gadget. |
| `L16_kw2.json` | 18 524 | ℚ(√3, √7, √11, √247) | A checkpoint of plain growth from `L16_seed.json` with 930 unit vectors and a proper 5-colouring of its first 18 484 points (`colouring_prefix`); with the edges at the orbit (14 ∓ 2√33)/3 added, its 5-colourability was undecided on 25 September 2026 (`notes/worker_jobs.md`). |
| `L16_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | The starting graph of the searches in L16: `five_tuned_16_1_3_7_11.json` and `five_rho7.json` glued along their common 402-point carrier in the Moser field (37 474 edges), with 918 unit vectors for growth and a proper 5-colouring. |
| `L16_skelG_seed.json` | 6 177 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with `W_moser_orbit_9_33.json` merged in and its 495 edges at d² = (9 ∓ √33)/6 stored as `two_edges`, for plain growth with those pairs as soft edges. |
| `L16_skeleton3_seed.json` | 6 293 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with three copies of `W_lattice_16_21_28_61.json` at vertices 2519, 4222 and 1130, turned by three unit vectors, and their 1 233 non-unit edges stored as `two_edges`. |
| `L16_skeleton_seed.json` | 6 151 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with `W_lattice_16_21_28_61.json` placed at vertex 2519 and its 411 non-unit edges (at 4/√3, √7, √(28/3) and √(61/3)) stored as `two_edges`. |
| `L16_two_apart73_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 2368, B = 2469 at distance √(7/3), for `MODE=apart` growth. |
| `L16_two_same133_seed.json` | 6 080 | ℚ(√3, √7, √11, √247) | `L16_seed.json` with the target pair A = 2135, B = 2692 at distance √(13/3), for `MODE=same` growth. |
| `W_lattice_16_21_28_61.json` | 72 | ℤ[ω]/√−3 | A vertex-critical graph on 72 points of ℤ[ω]/√−3 with 553 edges at distances 1, 4/√3, √7, √(28/3) and √(61/3) and no proper 5-colouring, rebuilt by `scripts/lattice_witness.py 3,16,21,28,61 5`. |
| `W_moser_orbit_9_33.json` | 187 | ℚ(√3, √11) | A vertex-critical graph on 187 points with 508 edges at distance 1 and 495 at the Galois orbit d² = (9 ∓ √33)/6 (d ≈ 0.7366, 1.5676) and no proper 5-colouring, checked by `scripts/orbit_witness_test.py` with orbit `9_1`. |
| `blocked_32312_coloured.json` | 32 312 | ℚ(√3, √11, √247) | Plain growth from `g2e_blocked_seed.json` to 32 312 points with a 5-colouring found by kissat, proper on the 282 909 edges along the file's 134 unit vectors but with 43 monochromatic pairs among the 223 further unit-distance pairs in other directions. |
| `blocking_requirements.json` | — | — | For a point p with 15 unit neighbours next to the 2 689-vertex carrier obtained by gluing `Sa` at the 12-point orbit of Sa[265], four pairs of those neighbours, with their squared distances, such that each 4-colouring found in three rounds (66 in all) that uses only three colours on the neighbours makes one of the pairs monochromatic (`experiments/fullcover.py`). |
| `chain23.json` | 10 | ℚ(√2, √3) | A chain of three unit rhombi (10 vertices, 16 edges) with no proper 3-colouring, the lower bound in χ(ℚ(√2, √3)²) = 4 (`tests/test_q23.py`, `certificates/chain23_no3coloring.json`). |
| `chain35.json` | 19 | ℚ(√3, √5) | A chain of six unit rhombi (19 vertices, 31 edges) joining the origin to the unit vector e^{iπ/6}, with no proper 3-colouring, the lower bound in 4 ≤ χ(ℚ(√3, √5)²) ≤ 5 (`tests/test_q35.py`, `scripts/experiments/chain35.py`). |
| `closure_five_247_c.json` | 1 851 | ℚ(√3, √11, √247) | The closure of `five_247_c.json` under adding the centre of every unit triangle, reached after four rounds. |
| `coset_unsplit_five_247_c.json` | — | — | The pairs of `five_247_c.json` at squared distance 5/9, 4/3, 3 or 25 that no coset 5-colouring splits: there are none (`experiments/cosetfilter.py`). |
| `coset_unsplit_tight_hexagon_4159.json` | — | — | The 24 pairs of `tight_hexagon_4159.json` at squared distance 5/9, 4/3, 3 or 25 that no coset 5-colouring splits, all at distance 5 (`experiments/cosetfilter.py`). |
| `deficiency_five_247_c.json` | — | — | Ten non-unit pairs of `five_247_c.json` whose addition as edges leaves no proper 5-colouring, so that at most ten such pairs are needed (δ ≤ 10; `experiments/deficiency.py`). |
| `disc_1099.json` | 1 099 | ℚ(√3, √11, √247) | The smallest disc about vertex 788 of `five_247.json` whose points still have no proper 4-colouring (`experiments/disc.py`). |
| `disc_1113.json` | 1 113 | ℚ(√3, √11, √247) | The smallest disc about the centroid of `five_247.json` whose points still have no proper 4-colouring (`experiments/disc.py`). |
| `ei_H214.json` | 214 | ℚ(√3, √11) | Exoo–Ismailescu's graph H as a point set (1 004 unit edges), with their pair A, B at distance 5, which is alike in every 5-colouring that has no monochromatic pair at distance 1 or 2. |
| `ei_rho7.json` | 638 | ℚ(√3, √11) | Exoo–Ismailescu's H with its images under ρ₇ and ρ₇⁻¹ about A (3 012 unit edges), on which the 2-adic 4-colouring of ℚ(√3, √11)² is proper (`notes/local_colourings.md`). |
| `escape5_five_247_c.json` | — | — | For a point p outside `five_247_c.json` with 12 unit neighbours in it, 32 sampled 5-colourings that use only two colours on them, each of which makes the pair (281, 655) at squared distance 3 monochromatic (`experiments/escape5.py`). |
| `five_23.json` | 7 141 | ℚ(√3, √11, √23) | `Sa` glued at the D₆-orbit of Sa[25] and spindled at squared distance 64/3 over the C₆-orbit of a pivot (sine √759/128, 759 = 3·11·23), with no proper 4-colouring (`tests/test_five_247.py`). |
| `five_23_v7.json` | 10 333 | ℚ(√3, √11, √23) | A larger graph built by the recipe of `five_23.json`, with no proper 4-colouring. |
| `five_247.json` | 1 139 | ℚ(√3, √11, √247) | A 5-chromatic graph: `Sa` glued to its 60° turn about the vertex Sa[25], then spindled, by a rotation through 2 arcsin(3/16), at a pair at distance 8/3 that is alike in every 4-colouring (`tests/test_five_247.py`). |
| `five_247_b.json` | 951 | ℚ(√3, √11, √247) | A 5-chromatic graph: `Sa` cut to its 340 highest-degree vertices, glued to its 120° image (overlap 178) and spindled at a forced pair at distance 8/3 (`tests/test_five_247.py`). |
| `five_247_c.json` | 803 | ℚ(√3, √11, √247) | A 5-chromatic, vertex-critical graph with 4 065 edges: `Sa` cut to its 327 highest-degree vertices, glued to its 120° image about a vertex, cut to the part that still forces, and spindled at distance 8/3 (`tests/test_five_247.py`). |
| `five_247_c_critical.json` | — | — | The check that `five_247_c.json` is vertex-critical, that is, deleting any vertex leaves a 4-colourable graph, with its ten slowest solver calls (`experiments/critical803.py`). |
| `five_247_c_lambda.json` | 1 605 | ℚ(√3, √11, √247) | `five_247_c.json` with its image under Exoo–Ismailescu's rotation λ = (49 + 3√−11)/50 about its densest vertex, a substrate whose edge module blocks coset colourings modulo 2, 3, 4 and 5 (`experiments/mk_lambda.py`). |
| `five_dense_10.json` | 12 469 | ℚ(√3, √11, √247) | A dense graph (mean degree 18.48) with no proper 4-colouring: ten stacked orbits of glue centres, with the composite forced pair tuned to distance 1 (`tests/test_dense_family.py`). |
| `five_dense_2.json` | 6 925 | ℚ(√3, √11, √247) | The construction of `five_dense_10.json` with two stacked orbits (mean degree 16.81; `tests/test_dense_family.py`). |
| `five_order6_free.json` | 12 011 | ℚ(√3, √11, √23, √247) | A graph with no proper 4-colouring: an order-6 carrier as in `order6_carrier.json` with its images under a rotation about each point of the forced orbit, through the free angle that brings one orbit point to distance 1 from another, which needs √23 (`experiments/freeorder6.py`). |
| `five_rho7.json` | 2 403 | ℚ(√3, √11, √247) | `five_247_c.json` with its images under ρ₇ = (1 + 4√−3)/7 and ρ₇⁻¹ about its densest vertex; `notes/rigidity.md` studies the coset colourings of its module. |
| `five_symmetric.json` | 7 141 | ℚ(√3, √11, √247) | A C₆-invariant graph with no proper 4-colouring: `Sa` glued at all six points of the C₆-orbit of Sa[25] and spindled at distance 8/3 over the whole orbit (`tests/test_five_247.py`). |
| `five_tuned_16_1_3_7_11.json` | 4 081 | ℚ(√3, √7, √11) | A graph of the tuned family (a forced pair at distance 8/3 composed with a turned copy of its carrier so that it lands at a chosen squared distance, here 16, then spindled), with no proper 4-colouring; its radical √1792 = 16√7 keeps it in ℚ(√3, √7, √11) (`tests/test_tuned_family.py`). |
| `five_tuned_1_1.json` | 2 041 | ℚ(√3, √5, √7, √11, √247) | The tuned family with the composite forced pair tuned to distance 1, so that no spindle is needed (radical √247). |
| `five_tuned_1_3.json` | 4 061 | ℚ(√3, √5, √7, √11, √23) | The tuned family with target squared distance 1/3, then spindled (radical √759, 759 = 3·11·23). |
| `five_tuned_2_1.json` | 4 081 | ℚ(√3, √5, √7, √11, √17) | The tuned family with target squared distance 2, then spindled (radical √476 = 2√(7·17)). |
| `five_tuned_4_1.json` | 4 081 | ℚ(√3, √5, √7, √11) | The tuned family with target squared distance 4, then spindled (radical √880 = 4√55), inside de Grey's field ℚ(√3, √5, √7, √11). |
| `five_twotune.json` | 12 240 | ℚ(√3, √11, √23) | The forced distance tuned to 1/√3 and then the angle to 60°: six copies whose pivot orbit has six pairs at distance 1, with no proper 4-colouring (`tests/test_two_tunings.py`). |
| `five_twotune_small.json` | 4 081 | ℚ(√3, √11, √23) | The construction of `five_twotune.json` with two copies instead of six, already without a proper 4-colouring. |
| `forced_certificates.json` | — | — | For the four unit neighbours of vertex 315 of `five_247_c.json`, the smallest ball about each pair, in the graph with vertex 315 deleted, in which the pair differs in every 4-colouring: the whole 802-vertex graph for every non-adjacent pair (`experiments/ball4.py`). |
| `forcing_set_five_247_c.json` | — | — | Two squared distances, (9 − √33)/6 and (7 − √33)/6 with 2 540 and 2 188 pairs, whose pairs added as edges leave `five_247_c.json` without a proper 5-colouring, a minimal such set of distances (`experiments/multidist.py`). |
| `free_angle.json` | 21 169 | ℚ(√3, √11, √23) | A graph with no proper 4-colouring from the free-angle spindle: in a carrier built on Sa[199] the vertices 726, 1526 and 2730 are alike in every 4-colouring (at squared distances 64/3 and 256/9 from 726), and a rotation about each point of the orbit of 726 brings 1526 to distance 1 from 2730 (`experiments/freeangle.py`). |
| `free_angle_both.json` | 39 313 | ℚ(√3, √11, √23) | The construction of `free_angle.json` with both roots of the rotation equation, with no proper 4-colouring. |
| `g2e_base.json` | 804 | ℚ(√3, √11, √247) | `five_247_c.json` with the point b = a + 2e added, where a, a + e is its only edge along e = (−3/16, √247/16): the pair A = a, B = b at distance 2. |
| `g2e_blocked_ckpt.json` | 4 425 | ℚ(√3, √11, √247) | A checkpoint of the distance-2 gadget growth from `g2e_blocked_seed.json`. |
| `g2e_blocked_seed.json` | 3 345 | ℚ(√3, √11, √247) | `five_247_c_lambda.json` with a pair A, B at distance 2 along a λ-rotated edge, the half-turn about their midpoint and both unit circles filled: a seed for a distance-2 gadget. |
| `g2e_nu_ckpt.json` | 5 927 | ℚ(√3, √11, √247) | A checkpoint of the distance-2 gadget growth from `g2e_nu_seed.json`. |
| `g2e_nu_seed.json` | 4 947 | ℚ(√3, √11, √247) | `five_247_c.json` with its images under ν = (−1 + 3√−11)/10 and ν̄ about its densest vertex, a pair A, B at distance 2, the half-turn about their midpoint and both unit circles filled. |
| `g2e_sym_ckpt.json` | 2 185 | ℚ(√3, √11, √247) | A checkpoint of the distance-2 gadget growth from `g2e_sym_seed.json`. |
| `g2e_sym_seed.json` | 1 705 | ℚ(√3, √11, √247) | `five_247_c.json` with a pair a, b = a + 2e along its √247 edge e, the half-turn about their midpoint and both unit circles filled along every edge direction. |
| `g2e_v2_ckpt.json` | 1 684 | ℚ(√3, √11, √247) | A checkpoint of the distance-2 gadget growth from `g2e_base.json`. |
| `hub_disjunction.json` | — | — | For the carrier and point of `blocking_requirements.json`, a hub disjunction: in each sampled 4-colouring with only three colours on the point's unit neighbours, the neighbour 495 shares its colour with one of three vertices W, at the listed squared distances (`experiments/hubdisj.py`). |
| `hunt_811.json` | 811 | ℚ(√3, √11, √247) | A subgraph of `five_247.json` with no proper 4-colouring, found by vertex deletion (seed 1). |
| `idealquot_ei_H214_3_1.json` | — | — | A proper 5-colouring of the Cayley graph of M/αM ≅ (ℤ/13)² on the unit vectors, for α = 3 + ω and M the module of `ei_H214.json` (`scripts/idealquot.py`). |
| `mu5_menu.json` | — | — | For 30 candidate points p outside `five_247_c.json` (images of its vertices under 60° and 30° rotations), the pairs of unit neighbours that must be forced apart to raise μ₅(p) above 2, with the size of N(p) and the squared distances of the pairs (`experiments/requirements.py`). |
| `mu5_requirement.json` | — | — | The search of `escape5_five_247_c.json` iterated: once the pair (281, 655) at squared distance 3 is required to differ, the solver, limited to 4 million conflicts per call, found no further 5-colouring with only two colours on the 12 neighbours (`experiments/escape5b.py`). |
| `narrow_disjunction.json` | — | — | A minimal set of 133 pairs of `five_247_c.json` at squared distances (9 − √33)/6 and (7 − √33)/6, one of which is monochromatic in every 5-colouring (`experiments/narrow.py`). |
| `order6_carrier.json` | 6 006 | ℚ(√3, √11, √247) | A 4-colourable carrier in which the composite isometry τ is a rotation of order 6, with a six-point orbit that is monochromatic in every 4-colouring, at squared distances 64/9, 64/3 and 256/9 from one another (`tests/test_order_six_tuning.py`). |
| `rarest_v139.json` | 24 781 | ℚ(√3, √5, √7, √11, √29) | A graph with no proper 4-colouring, spindled over six pivots at the rare forced squared distance 256/9 of the glue centre Sa[139], which needs √1015 = √(5·7·29) (`experiments/rarest.py`). |
| `rarest_v199.json` | 39 313 | ℚ(√3, √5, √7, √11, √23, √29) | A graph with no proper 4-colouring built on the glue centre Sa[199] from its rare forced squared distances 64/3 and 256/9, in a field that adds √5, √7 and √29. |
| `rho7_H.json` | 2 606 | ℚ(√3, √11, √247) | `five_rho7.json` with Exoo–Ismailescu's H, turned by 150° with A at the densest vertex; `two_edges` lists H's 446 pairs at distance 2. |
| `rho7_lambda_K.json` | 5 210 | ℚ(√3, √11, √247) | `five_rho7.json` ∪ λ_O(`five_rho7.json`) ∪ K with O = A, where K is Exoo–Ismailescu's graph turned onto the module; `two_edges` lists K's 892 pairs at distance 2 (written by `scripts/mk_kseed.py`). |
| `rho7_pair2_seed.json` | 4 981 | ℚ(√3, √11, √247) | `five_rho7.json` with a pair a, a + 2e along a unit edge e, the half-turn about their midpoint and both unit circles filled: a seed for a distance-2 gadget. |
| `rho7_pair5_ckpt.json` | 6 878 | ℚ(√3, √11, √247) | A checkpoint of the growth from `rho7_pair5_seed.json`. |
| `rho7_pair5_seed.json` | 4 998 | ℚ(√3, √11, √247) | The seed of `rho7_pair2_seed.json` with the pair a, a + 5e at distance 5 instead. |
| `shrunk_1135.json` | 1 135 | ℚ(√3, √11, √247) | A subgraph of `five_247.json` with no proper 4-colouring, extracted from an unsatisfiable core (seed 0; `experiments/shrink_core.py`). |
| `tight_four.json` | 24 | ℚ(√3, √5, √7, √11) | A unit-distance graph with 54 edges whose chromatic number and circular chromatic number are both 4 (`tests/test_tight_four.py`). |
| `tight_hexagon_4159.json` | 4 159 | ℚ(√3, √11, √247) | Six copies of `five_247_c.json` with their densest vertices on the vertices of a hexagon, each turned by a different multiple of 60°, with no proper 4-colouring. |
| `two_distance_five_247_c.json` | — | — | No squared distance among the 60 most frequent in `five_247_c.json` is forced at five colours: with the pairs at any one of them added as edges, the graph stays 5-colourable (`experiments/twodist.py`). |
| `witness_five.json` | 63 | ℚ(√3, √5, √7, √11) | 63 points with 119 unit edges and 161 listed pairs at squared distances 15/16, 16, 17/2, 9 and 7, such that no proper 5-colouring makes every listed pair bichromatic (`tests/test_homcol.py`). |
