"""Sympy helper for `Q311.lean`.

It generates the edge list and the coefficients of the `linear_combination` proofs. A
polynomial identity that holds modulo `r3^2 = 3`, `r11^2 = 11` is written as
`E = Q0*(r3^2 - 3) + Q1*(r11^2 - 11)`; the script prints `Q0` and `Q1`.

  1. `Q311.moserEdges` and `Q311.adj..`. The certificate
     `research/hadwiger-nelson/certificates/moser_spindle_no3coloring.json` lists only the seven
     vertices (basis `1, √3, √11, √33`). The edges are all pairs at distance 1, found exactly
     here (eleven of them), with the coefficients of each unit-distance proof.
  2. Checks of the identities behind `Q311.piL_sq`, `Q311.two_eq`, `Q311.beta_sq` and
     `Q311.exists_int_combo_beta`, for both signs `t = ±√11`.

Run: `python3 lean/tools/q311_coefficients.py` from `research/hadwiger-nelson` (needs sympy).
"""
import json
import os

from sympy import Rational, expand, reduced, symbols

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
r3, r11 = symbols('r3 r11')


def quotients(E):
    """Return (Q0, Q1), or None if E does not vanish modulo r3^2 = 3, r11^2 = 11."""
    Q, R = reduced(expand(E), [r3**2 - 3, r11**2 - 11], r3, r11)
    if expand(R) != 0:
        return None
    Q = list(Q) + [0] * (2 - len(Q))
    return [expand(q) for q in Q]


def lean(e):
    return '(' + str(expand(e)).replace('**', '^').replace('r11', '√11').replace('r3', '√3') + ')'


# 1. the Moser spindle
data = json.load(open(f'{REPO}/certificates/moser_spindle_no3coloring.json'))
basis = [1, r3, r11, r3*r11]
pts = [tuple(sum(Rational(c) * b for c, b in zip(v[k], basis)) for k in ('x', 'y'))
       for v in data['vertices']]
edges = []
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        (x1, y1), (x2, y2) = pts[i], pts[j]
        q = quotients((x1 - x2)**2 + (y1 - y2)**2 - 1)
        if q is not None:
            edges.append((i, j))
            print(f'adj{i}{j}: linear_combination {lean(q[0])} * sqrt3_sq + {lean(q[1])} * sqrt11_sq')
assert len(edges) == data['m'] == 11
print('moserEdges =', edges)

# 2. the elements pi and beta, for t = sqrt(11) and t = -sqrt(11)
pi = r3 - 1
assert quotients(pi**2 - (-2*pi + 2)) is not None                   # piL_sq
assert expand(2 - pi**2*(pi + 3) + (pi + 1)*(pi**2 + 2*pi - 2)) == 0  # two_eq
a, b, c, e = symbols('a b c e')
for t in (r11, -r11):
    beta = (1 + r3*t)/2
    assert quotients(beta**2 - (beta + 8)) is not None               # beta_sq
    lhs = 3*(a + b*r3 + c*t + e*r3*t)                                # exists_int_combo_beta
    rhs = (3*a + 3*b - c - 3*e) + (3*b - c)*pi + ((2*c + 6*e) + 2*c*pi)*beta
    assert quotients(lhs - rhs) is not None
print('piL_sq, two_eq, beta_sq, exists_int_combo_beta: identities verified')
