# A Moser-spindle-free 5-chromatic unit-distance graph on 852 vertices: the data

The data of `notes/flat852.md` (29 September 2026). The graph is not known to be vertex-critical.

- `core852a.json`: 852 vertices, 4487 edges. A vertex is (c0 + c1 z + ... + c11 z^11)/7 with `exact7` = [c0..c11],
  z = zeta21, embedded by z -> exp(2 pi i/21). Every edge vector lies in U = mu42 ∪ omega*mu42 ∪ conj(omega)*mu42
  (126 unit vectors; `directions.json`), omega = zeta3/(zeta7 - zeta7^-1) - zeta3^2/(zeta7^2 - zeta7^-2), the
  second direction orbit of J. K. Haugland's heptagon lattice (arXiv 2608.04542). Coordinates lie in Q(zeta84)^+;
  there is no Moser spindle (it would need sqrt(-11) in Q(zeta21)).
- Not 4-colourable: `core852a.cnf` (one triangle fixed), kissat 4.0.4 UNSAT, drat-trim VERIFIED (`logs/`).
  A proper 5-colouring is `core852.5colouring.json`, so the chromatic number is 5.
- `g1023.json`: the 1023-vertex graph it was cut from (same certificates).
- `independent-check/verify852.py`: a second, independent check (exact unit distances with sympy in Q(zeta21),
  50-digit numerics, its own CNF encoding).
- `scripts/`: the growth, minimisation and certification code (paths written as SC = the working directory).
- `REPORT-growth-agent.md`: how the graph was found (Haugland's lattice L plus its mirror image; colouring-guided
  growth from H ∪ conj(H)).
