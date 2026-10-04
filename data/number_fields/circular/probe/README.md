# The probe below `[1/3, 2/3]`: computations for Theorem D

Companion to `notes/circular_planes.md` §5 (Theorem D, Proposition D1) and to Section 9 of
`papers/three-colours/` (Theorem D, Proposition 6). All programs use exact rational arithmetic
(`fractions.Fraction`); linear-programming optima are found in floating point and then certified by exact primal
and dual solutions. In the scripts, *Theorem A* is Proposition D1 of the note (Proposition 6 of the paper), `r` is
the end of the interval `[r, 1 − r]` (`θ` in the paper), `Sq_s = [−s, s]²` with `s = 1/2 − r`, `c*_N = N(1 + i)/2`,
and `Q_N` is the set of type q points `(N/3)(a + bi)`, `a, b ∈ {1, 2}`.

The subfolder [`indep/`](indep/) holds a referee's implementation, written from scratch and sharing no code with this
folder, and its outputs; the two agree on every number quoted in the note and the paper.

| file | what it does |
|---|---|
| `snr.py` | the core: `S_N^(r)` for `N = 5^k` as a union of convex rational polygons, by lifting from level `k − 1` (cut along the strips of the conditions at `ρ^{±k}`) or directly |
| `lp3.py` | exact 3-variable linear programming with certificates |
| `check_pattern.py` → `check_pattern.txt` | Lemma 1 of the note: the values of `c*_N` and of the points of `Q_N` at `ρ^j` |
| `radii.py` → `radii.txt` | the 12-gon `P_1^(1)` and the hexagons `Z(q)` (radii `√10/3` and `√2`) |
| `base_case.py` → `base_case.txt` | the 25 classes of the base case `k = 1` at `r = 17/56` |
| `base_caseC.py` → `base_caseC.txt` | the classes with both shifts non-integral: least `s` and three-term certificates (minimal shifts) |
| `lemma_lp.py`, `lemma_more.py` → `lemma_2_14.txt` | exact thresholds of the induction lemmas C and Q with the full hypotheses `|j| ≤ k − 1`, for `2 ≤ k ≤ 14` |
| `certificates.py` → `certificates.txt` | the LP certificates behind Lemmas C and Q at `k = 2`, in readable form |
| `check_thmA.py` → `check_thmA.txt` | Proposition D1 for `r` slightly above `17/56`: only the 5 main components at every level checked, with the exact shapes |
| `check_components.py` → `check_components.txt` | every computed polygon is a whole connected component of `S_N^(r)` |
| `cross_check.py` → `cross_k1.txt`, `cross_k2.txt` | floating-point random sampling against the exact components |
| `first_new.py` → `first_new.txt` | the extra isolated points at `r = 17/56` for `N = 5, 25, 125, 625` (sharpness of `17/56`) |
| `levels2.py`, `levels3.py` → `levels_r17_56.txt`, `levels_r3001.txt`, `levels3_r3001_18.txt`, `levels3_r3001611_19.txt` | level-by-level enumeration with genealogy and the exact `κ` of every extra component; at `r = 0.3001611` only the main components remain from `k = 17` on (Theorem D(b)) |
| `far_points.py` → `far_points.txt` | the 5-adic points `c_k = c*_N + ρ^k/5 − 5^{k−1}(2 − i)` of `S_N^(3/10)`, far from all type points (the limit of the method at `10/3`) |
| `q1_tables.py`, `q1_compact.py` → `q1_k12.txt`, `q1_k3.txt`, `q1_compact.txt` | `S_N^(r)` for `N = 5, 25, 125` at several `r`, by two enumerations compared component by component |
| `chain.py` → `chain_17_56.txt`, `genealogy.py` | following the extra components through the levels |
| `torsion_search.py` → `torsion_260.txt` | torsion characters of `ℤ[1/5][i]` of order up to 260: none with `κ ≥ 0.29` |
| `chartypes_k3_002.txt` | output of `../chartypes.py` (the earlier grid search, grid `0.02`), kept for the record |

Referee's implementation ([`indep/`](indep/)):

| file | what it does |
|---|---|
| `hgeom.py`, `sn.py` | polygons in H-representation and `S_N^(r)` (lift and direct enumerations) |
| `run_structure.py` → `out_structure_k4.txt`, `out_structure_r1756.txt` | Proposition D1 for `k ≤ 4` at five values of `r`, and the extra points at `r = 17/56` |
| `lemmas.py`, `lemmaQ2.py`, `lp3x.py`, `lpcert.py` → `out_lemmas.txt`, `out_lemmaQ2.txt` | Lemmas 1, 2, C and Q, and the exact LP solver |
| `base_case_ref.py` → `out_base_case.txt` | the 25 classes of the base case with all shifts (least `s`: 0, `11/56` for 8 classes, and `1/4`, `2/7`, `1/4`, `1/6` for the four orbits), and the 12 certificates of `../base_caseC.txt` re-evaluated |
| `sample_check.py` → `out_sample.txt` | exact membership test at 350 000 random points |
| `r0check.py`, `deep_lift.py` → `out_r0.txt`, `out_lift_333.txt`, `out_deep.txt` | the thresholds `333/1106` (level 5) and `3303/10981` (level 6), and the levels up to 19 at `r = 0.3001611` and `0.3001609` |
| `thresholds.py` → `out_thresholds.txt` | the thresholds of Lemmas C and Q with the full hypotheses for `2 ≤ k ≤ 22` |
| `padic_family.py` → `out_padic.txt` | the 5-adic points `c_k` for `k ≤ 10` |

To rerun the main checks (a few minutes each):

```sh
python3 check_thmA.py 3036/10000 9     # Proposition D1 at r = 0.3036 > 17/56, levels 1 to 9
python3 levels3.py 3001611/10000000 19 > levels3_r3001611_19.txt
cd indep && python3 deep_lift.py 19 3001611/10000000 3001609/10000000 && python3 thresholds.py 22
```
