# Circular chromatic number 4: finite witnesses

`notes/circular_planes.md` §6.9 (in the paper `papers/three-colours/`, the subsection "Finite witnesses at four")
proves:

- at four colours the winding argument of Theorem W⁺ fails only on a *tight square*
  `x → x + s → x + s + t → x + t → x`; so `Cay(Γ, S)` has a 4-colouring without tight squares exactly when some
  character maps `S` into `[1/4, 3/4]` (Lemma F13);
- if `S` is finite, `Cay(Γ, S)` is 4-colourable and `κ(S) ≤ 1/4` (no character maps `S` into the open interval
  `(1/4, 3/4)`), then some finite subgraph has circular chromatic number 4 (Corollary F14);
- whenever `χ_c(F²) = 4` for a number field `F`, some finite unit-distance graph in `F²` has `χ_c = 4`
  (Theorem F16), for example over `ℚ(√59)`, `ℚ(√83)` and `ℚ(√3, √11)`; with Corollary F12, `χ_c(F²)` is attained
  by a finite unit-distance graph whenever it is at most 4.

The proof gives no bound on the size of such a graph, and we have not found one. The files here concern `ℚ(√59)`.

A unit vector with denominator `D` is `u = ((a + b√59)/D, (c + e√59)/D)` with `a² + 59b² + c² + 59e² = D²` and
`ab + ce = 0`; `U_D` is the set of all of them. For `θ ∈ ℚ⁴`, `x ↦ θ₁a + θ₂b + θ₃c + θ₄e mod 1`, on the points `x`
of the plane with these coordinates over `D`, is a character of that group `(1/D)ℤ⁴`.

| File | What it shows |
|---|---|
| `theta59_210.json` | `θ = (89/118, 1/2, 89/118, 1/2)` keeps all 108 vectors of `U_210` at distance at least `15/59` from `ℤ`. So `κ(U_210) ≥ 15/59 > 1/4`, and every graph on such points with edges in `U_210` maps to `K_{59/15}`: none has `χ_c = 4`. The growth with `D = 210` of the research log (6 October) could not succeed. |
| `theta59_1050.json` | `θ = (3/8, 5/8, 5/8, 5/8)` keeps all 300 vectors of `U_1050` at distance at least `1/4`, and exactly `1/4` for some: `κ(U_1050) ≥ 1/4`. A floating-point MILP (`kapparel.py`) finds no character with margin `0.2501`; if `κ(U_1050) = 1/4`, Corollary F14 puts a finite witness inside `Cay(ℤU_1050, U_1050)`. An exact certificate for the open interval is being computed; it is not claimed here. |
| `check_theta.py d D THETA.json` | lists `U_D` by brute force over `(b, e, a, c)` and computes the least margin of `θ` exactly. |
| `indep/` | the referee's programs: `at4check.c` and `at4check_py.py`, which share no code, check Lemma F13, Corollary F14 and the four-colour form of Theorem W⁺ on every abelian group of order at most 16 (and `ℤ/18`, `ℤ/3 × ℤ/6`; `ℤ/20` and `ℤ/2 × ℤ/10` were stopped for time), with outputs in `out_c/`, `out_py/` and the logs; `compare.py`, `aggregate.py`; `example_z8.py`, a tight 4-cycle that is not a square; `kcert210_verify.py`, another character with margin `15/59` on `U_210`. |

The finite-group runs cover 474 classes of connection sets with a 4-colourable Cayley graph (105 with `κ < 1/4`, 131
with `κ = 1/4`, 238 with `κ > 1/4`) and 163 613 685 proper 4-colourings, with no failure: every colouring in which the
two sides of a square differ has a tight square there, the averaged slopes of every colouring without tight squares
form a character, a colouring without tight squares exists exactly when `κ ≥ 1/4`, and a colouring without tight
cycles exactly when `κ > 1/4`. (In finite groups the last fact also follows from Theorem W⁺ and Guichard's lemma; the
step of Corollary F14 that needs an infinite group, Gordan's theorem, was checked by hand.)

Check:

```sh
python3 check_theta.py 59 210 theta59_210.json        # least margin 15/59
python3 check_theta.py 59 1050 theta59_1050.json      # least margin 1/4 (about 30 s)
cd indep && python3 example_z8.py && python3 at4check_py.py 2 4     # a few seconds
```
