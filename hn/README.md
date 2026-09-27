# `hn`: the library

`hn` is the Python package behind every result in this project. Coordinates are exact elements of
number fields such as ℚ(√3, √11), never floats, so "distance exactly 1" is a decidable test.
Colourability questions go to SAT solvers.

| module | what it provides |
|---|---|
| **Arithmetic and geometry** | |
| `field.py` | exact arithmetic in multiquadratic fields ℚ(√d₁, …, √dₖ) |
| `geometry.py` | points and rotations of the plane with exact coordinates |
| `cyclotomic.py`, `cyclograph.py` | points in ℤ[ζₙ] and the unit-distance graphs they span |
| `quadext.py`, `realext.py` | quadratic extensions of cyclotomic and multiquadratic fields |
| `fast.py` | a vectorised integer representation, for searching at scale |
| **Graphs, colouring and certificates** | |
| `graph.py` | unit-distance graphs: construction, reduction, input and output |
| `generate.py` | vertex-set generators |
| `coloring.py` | k-colourability by SAT, and shrinking a witness |
| `certify.py` | certificates: exactly what a third party has to check |
| `density.py` | independence ratios, as finite certificates for measurable bounds |
| **Forcing and constructions** | |
| `forced.py` | pairs forced apart or alike in every colouring, and pressure |
| `spindle.py`, `multispindle.py` | the spindle argument, automated and in general form |
| `degrey.py` | de Grey's 1581-vertex graph, rebuilt from the published recipe |
| `mixed.py`, `transfer.py`, `transversal.py` | rotations built from conflicts, gluing along interfaces, and the spindle method's ceiling |
| `blocked.py`, `homcol.py` | the μ invariant and homomorphism (coset) colourings |
| **Local colourings** | |
| `adelic.py` | colourings through one prime, including the 2-adic 4-colourings of ℚ(√3, √11)² and ℚ(√2, √3)² |
| `cli.py` | command-line entry points |

The tests in [`../tests/`](../tests/) exercise every module.
