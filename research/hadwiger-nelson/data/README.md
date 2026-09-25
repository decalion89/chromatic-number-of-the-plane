# Data: graphs and witnesses in exact coordinates

Every graph here is stored with **exact** coordinates, so every distance can be recomputed exactly.

## Format

Each file is JSON with at least:
- `field_generators`: `[a, b, …]`, meaning the coordinates lie in ℚ(√a, √b, …);
- `points`: a list of points `[x, y]`. Each coordinate is a list of `[numerator, denominator]`
  pairs on the basis of square-free products of the generators, in the order
  `hn.field.Field(field_generators)` lists them (for `[2, 3]`: 1, √2, √3, √6).

Optional keys include:
- `note`: what the graph is;
- `units`: the unit vectors a growth may use;
- `A`, `B`: a target pair;
- `two_edges`, `dist2`, `dist2_exact`: extra distances for two-distance graphs;
- `colouring`, `colouring_prefix`: a proper colouring found by a search;
- `status`, `mode`: the state of a growth run.

`scripts/verify_pair.py` and the tests rebuild graphs from these files alone.

## Naming

| prefix | meaning |
|---|---|
| `five_*` | 5-chromatic unit-distance graphs (no proper 4-colouring) |
| `W_*` | witnesses: graphs with extra forbidden distances and no proper 5-colouring |
| `L16_*` | the field ℚ(√3, √7, √11, √247): seeds for the searches for six |
| `L16_g<d>_<A>_<B>_seed` | a seed targeting the pair (A, B) at a Galois-orbit distance d |
| `ei_*`, `K_rot` | Exoo–Ismailescu's graphs, rebuilt |
| `g2e_*`, `rho7_*` | distance-2 gadget searches, in the ρ₇ module |
| `*_seed`, `*_ckpt` | the start and a checkpoint of a growth search |

## Key files

| file | what it is |
|---|---|
| `chain23.json` | 10 vertices, not 3-colourable: the lower bound in χ(ℚ(√2, √3)²) = 4 (`notes/local_colourings.md` §10). |
| `W_moser_orbit_9_33.json` | 187 points, not 5-colourable with edges at 1 and at the Galois orbit d² = (9 ∓ √33)/6: a witness that needs one gadget. |
| `W_lattice_16_21_28_61.json` | 72 lattice points, not 5-colourable with edges at 1, 4/√3, √7, √(28/3), √(61/3). |
| `five_247.json` | 1139 vertices, 5-chromatic, in ℚ(√3, √11, √247). |
| `five_247_c.json` | 803 vertices, 5-chromatic and vertex-critical, in ℚ(√3, √11, √247). |
| `five_tuned_16_1_3_7_11.json` | 4081 vertices, 5-chromatic, in ℚ(√3, √7, √11): no √5 and no √247 needed. |
| `five_rho7.json` | 803-type graph together with its ρ₇-images: the module whose coset colourings are rigid. |
| `ei_H214.json, ei_rho7.json, K_rot.json` | Exoo–Ismailescu's graphs, rebuilt from their paper and turned onto the project's modules. |
| `L16_seed.json` | 6080 points: the starting graph of the searches for six in the field L16. |
| `L16_kw2.json` | 18 524 points: the largest L16 growth kept, at the edge of 5-colourability. |
| `tight_four.json, witness_five.json` | Small certificate-style graphs (circular chromatic number 4; a forced-pair witness). |

## All files

| file | points | field | note |
|---|---:|---|---|
| `K_rot.json` | 426 | ℚ(√3, √11, √247) | E-I K = H u lambda_A(H), H turned 150 degrees onto five_rho7's units, A at five_rho7's densest vertex |
| `L16_W0915_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed as a {1, d, sigma(d)}-graph, d^2 = (14 -+ 2 sqrt33)/3 (d = 0.9149, 2.9149): grow (MODE=plain) to a witness that is not 5-coloura… |
| `L16_apart43_rank.json` |  |  |  |
| `L16_apart43_survivors.json` |  |  |  |
| `L16_g0737_1021_1317_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=1021, B=1317 at d^2 = 3/2 + -1/6*sqrt(33) (Galois orbit (9 -+ sqrt33)/6, the orbit of the verified 187-poin… |
| `L16_g0737_28_868_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=28, B=868 at d^2 = 3/2 + -1/6*sqrt(33) (Galois orbit (9 -+ sqrt33)/6, the orbit of the verified 187-point w… |
| `L16_g0915_32_141_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=32, B=141 at d^2 = 14/3 + -2/3*sqrt(33) (Galois orbit (14 -+ 2 sqrt33)/3), split in all 16 tabu colourings;… |
| `L16_g0915_36_628_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=36, B=628 at d^2 = 14/3 + -2/3*sqrt(33) (Galois orbit (14 -+ 2 sqrt33)/3), split in all 16 tabu colourings;… |
| `L16_g1568_1048_1620_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=1048, B=1620 at d^2 = 3/2 + 1/6*sqrt(33) (Galois orbit (9 -+ sqrt33)/6, the orbit of the verified 187-point… |
| `L16_g1568_28_600_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=28, B=600 at d^2 = 3/2 + 1/6*sqrt(33) (Galois orbit (9 -+ sqrt33)/6, the orbit of the verified 187-point wi… |
| `L16_g2915_2437_3218_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=2437, B=3218 at d^2 = 14/3 + 2/3*sqrt(33) (Galois orbit (14 -+ 2 sqrt33)/3), split in all 16 tabu colouring… |
| `L16_g2915_41_715_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed with a target pair A=41, B=715 at d^2 = 14/3 + 2/3*sqrt(33) (Galois orbit (14 -+ 2 sqrt33)/3), split in all 16 tabu colourings; … |
| `L16_gadget43_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed: five_tuned_16 U five_rho7 glued along their common Moser-field carrier; units = F8 lambda-closure (666) U the 247 module (134);… |
| `L16_kw2.json` | 18524 | ℚ(√3, √7, √11, √247) |  |
| `L16_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed: five_tuned_16 U five_rho7 glued along their common Moser-field carrier; units = F8 lambda-closure (666) U the 247 module (134);… |
| `L16_skelG_seed.json` | 6177 | ℚ(√3, √7, √11, √247) | L16 seed with the verified 187-point Galois witness data/W_moser_orbit_9_33.json merged in; two_edges are W's 495 edges at d^2 = (9 -+ sq… |
| `L16_skeleton3_seed.json` | 6293 | ℚ(√3, √7, √11, √247) | L16 seed: five_tuned_16 U five_rho7 glued along their common Moser-field carrier; units = F8 lambda-closure (666) U the 247 module (134);… |
| `L16_skeleton_seed.json` | 6151 | ℚ(√3, √7, √11, √247) | L16 seed: five_tuned_16 U five_rho7 glued along their common Moser-field carrier; units = F8 lambda-closure (666) U the 247 module (134);… |
| `L16_two_apart73_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed: five_tuned_16 U five_rho7 glued along their common Moser-field carrier; units = F8 lambda-closure (666) U the 247 module (134);… |
| `L16_two_same133_seed.json` | 6080 | ℚ(√3, √7, √11, √247) | L16 seed: five_tuned_16 U five_rho7 glued along their common Moser-field carrier; units = F8 lambda-closure (666) U the 247 module (134);… |
| `W_lattice_16_21_28_61.json` | 72 |  |  |
| `W_moser_orbit_9_33.json` | 187 | ℚ(√3, √11) | Not 5-colourable with edges at distance 1 and at the Galois orbit d^2 = (9 -+ sqrt33)/6 (d = 0.7366, 1.5676): a witness needing one gadge… |
| `blocked_32312_coloured.json` | 32312 | ℚ(√3, √11, √247) |  |
| `blocking_requirements.json` | 2689 |  |  |
| `chain23.json` | 10 | ℚ(√2, √3) | chi(Q(sqrt2, sqrt3)^2) >= 4: a chain of three unit rhombi O -> sqrt3 u1 -> sqrt3 (u1+u2) -> sqrt3 (u1+u2+u3) with u1 = (1, 0), u2 = (0, 1… |
| `closure_five_247_c.json` | 1851 | ℚ(√3, √11, √247) |  |
| `coset_unsplit_five_247_c.json` |  |  |  |
| `coset_unsplit_tight_hexagon_4159.json` |  |  |  |
| `deficiency_five_247_c.json` |  |  |  |
| `disc_1099.json` | 1099 | ℚ(√3, √11, √247) |  |
| `disc_1113.json` | 1113 | ℚ(√3, √11, √247) |  |
| `ei_H214.json` | 214 | ℚ(√3, √11) | Exoo-Ismailescu H (214 points), as a point set; A, B their forced pair at distance 5 |
| `ei_rho7.json` | 638 | ℚ(√3, √11) | E-I H with its rho7^{+-1} images about A; unit edges only |
| `escape5_five_247_c.json` |  |  |  |
| `five_23.json` | 7141 | ℚ(√3, √11, √23) |  |
| `five_23_v7.json` | 10333 | ℚ(√3, √11, √23) |  |
| `five_247.json` | 1139 | ℚ(√3, √11, √247) |  |
| `five_247_b.json` | 951 | ℚ(√3, √11, √247) |  |
| `five_247_c.json` | 803 | ℚ(√3, √11, √247) |  |
| `five_247_c_critical.json` | 803 |  |  |
| `five_247_c_lambda.json` | 1605 | ℚ(√3, √11, √247) | five_247_c u lambda_c(five_247_c), lambda=(49+3sqrt-11)/50 about the densest vertex |
| `five_dense_10.json` | 12469 | ℚ(√3, √11, √247) | ten stacked orbits, tuned chain to distance 1 |
| `five_dense_2.json` | 6925 | ℚ(√3, √11, √247) |  |
| `five_order6_free.json` | 12011 | ℚ(√3, √11, √23, √247) | angle tuned to order 6, then the free angle spent |
| `five_rho7.json` | 2403 | ℚ(√3, √11, √247) | five_247_c u ['rho7', 'rho7inv'] about its densest vertex |
| `five_symmetric.json` | 7141 | ℚ(√3, √11, √247) |  |
| `five_tuned_16_1_3_7_11.json` | 4081 | ℚ(√3, √7, √11) | composed forcing tuned, then spindled |
| `five_tuned_1_1.json` | 2041 | ℚ(√3, √5, √7, √11, √247) | composed forcing tuned to 1 |
| `five_tuned_1_3.json` | 4061 | ℚ(√3, √5, √7, √11, √23) | composed forcing tuned, then spindled |
| `five_tuned_2_1.json` | 4081 | ℚ(√3, √5, √7, √11, √17) | composed forcing tuned, then spindled |
| `five_tuned_4_1.json` | 4081 | ℚ(√3, √5, √7, √11) | composed forcing tuned, then spindled |
| `five_twotune.json` | 12240 | ℚ(√3, √11, √23) | distance tuned to 1/3, then angle tuned to order 6 |
| `five_twotune_small.json` | 4081 | ℚ(√3, √11, √23) | distance tuned to 1/3, angle 60, TWO copies |
| `forced_certificates.json` |  |  |  |
| `forcing_set_five_247_c.json` |  |  |  |
| `free_angle.json` | 21169 | ℚ(√3, √11, √23) |  |
| `free_angle_both.json` | 39313 | ℚ(√3, √11, √23) |  |
| `g2e_base.json` | 804 | ℚ(√3, √11, √247) | five_247_c plus b = a + 2e, e = (-3/16, sqrt247/16) the single sqrt247 edge a -- a+e |
| `g2e_blocked_ckpt.json` | 4425 | ℚ(√3, √11, √247) |  |
| `g2e_blocked_seed.json` | 3345 | ℚ(√3, √11, √247) | blocked substrate (803 u lambda 803) u its 180-degree turn about m u both unit circles; pair along a lambda-rotated edge |
| `g2e_nu_ckpt.json` | 5927 | ℚ(√3, √11, √247) |  |
| `g2e_nu_seed.json` | 4947 | ℚ(√3, √11, √247) | five_247_c.json + ['nu', 'nubar'] about its densest vertex |
| `g2e_sym_ckpt.json` | 2185 | ℚ(√3, √11, √247) |  |
| `g2e_sym_seed.json` | 1705 | ℚ(√3, √11, √247) | five_247_c u R_m(five_247_c) u both unit circles along every edge direction; a, b = a + 2e, e the sqrt247 edge |
| `g2e_v2_ckpt.json` | 1684 | ℚ(√3, √11, √247) |  |
| `hub_disjunction.json` | 2689 |  |  |
| `hunt_811.json` | 811 | ℚ(√3, √11, √247) |  |
| `idealquot_ei_H214_3_1.json` |  |  |  |
| `mu5_menu.json` |  |  |  |
| `mu5_requirement.json` |  |  |  |
| `narrow_disjunction.json` |  |  |  |
| `order6_carrier.json` | 6006 | ℚ(√3, √11, √247) | tuned to angle 60, tau of order 6 |
| `rarest_v139.json` | 24781 | ℚ(√3, √5, √7, √11, √29) |  |
| `rarest_v199.json` | 39313 | ℚ(√3, √5, √7, √11, √23, √29) |  |
| `rho7_H.json` | 2606 | ℚ(√3, √11, √247) | five_rho7 u H (E-I, turned 150 degrees, A at the densest vertex); two_edges = H's 446 distance-2 pairs |
| `rho7_lambda_K.json` | 5210 | ℚ(√3, √11, √247) | five_rho7 u lambda_O(five_rho7) u K (E-I's graph, turned onto the module), O = A |
| `rho7_pair2_seed.json` | 4981 | ℚ(√3, √11, √247) | five_rho7.json, pair a, a+2e, turn about the midpoint, circles filled |
| `rho7_pair5_ckpt.json` | 6878 | ℚ(√3, √11, √247) |  |
| `rho7_pair5_seed.json` | 4998 | ℚ(√3, √11, √247) | five_rho7.json, pair a, a+5e, turn about the midpoint, circles filled |
| `shrunk_1135.json` | 1135 | ℚ(√3, √11, √247) |  |
| `tight_four.json` | 24 | ℚ(√3, √5, √7, √11) | this unit-distance graph has circular chromatic number exactly 4: it is 4-colourable and not 3-colourable, and it admits no homomorphism … |
| `tight_hexagon_4159.json` | 4159 | ℚ(√3, √11, √247) |  |
| `two_distance_five_247_c.json` |  |  |  |
| `witness_five.json` | 63 | ℚ(√3, √5, √7, √11) | no proper 5-colouring of the unit-distance graph on these points leaves every listed pair bichromatic |
