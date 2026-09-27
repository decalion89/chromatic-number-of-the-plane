"""Sympy helper for `LocalColouring.lean` and `Q23.lean`.

It generates the coefficients of the `linear_combination` proofs. A polynomial identity that
holds modulo `u^2 = m`, `v^2 = n` is written as `E = Q0*(u^2 - m) + Q1*(v^2 - n)`; the script
prints `Q0` and `Q1`, which are the coefficients of the hypotheses `u^2 = m`, `v^2 = n`.

  1. `LocalColouring.exists_int_combo`, multiplication step: the product of
     `a1 + b1 u + c1 v + e1 uv` and `a2 + b2 u + c2 v + e2 uv`, for symbolic `m`, `n`.
  2. `Q23.piL_eisenstein`, `Q23.s2_eq`, `Q23.s3_eq`, `Q23.s6_eq`: identities in `s2 = √2`,
     `s3 = √3` for `π = (√2 + √6)/2 - 1`.
  3. `Q23.epsL` and `Q23.two_eq_pi_pow_four_mul`: `ε = 2/π⁴` in `ℤ[π]`, and the cofactor of the
     Eisenstein polynomial in `2 - π⁴ε`.
  4. `Q23.adj..`: the sixteen unit distances of `data/chain23.json`
     (coordinates on the basis `1, √2, √3, √6`).

Run: `python3 lean/tools/q23_coefficients.py` from the root of the repository (needs sympy).
"""
import json
import os

from sympy import Rational, div, expand, invert, reduced, rem, symbols

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
s2, s3, x = symbols('s2 s3 x')


def quotients(E, u, v, m, n):
    """Return (Q0, Q1) with E = Q0*(u^2 - m) + Q1*(v^2 - n); fail if the remainder is not 0."""
    Q, R = reduced(expand(E), [u**2 - m, v**2 - n], u, v)
    Q = list(Q) + [0] * (2 - len(Q))
    assert expand(R) == 0, R
    return [expand(q) for q in Q]


def lean(e, real=False):
    txt = str(expand(e)).replace('**', '^')
    if real:  # the unit-distance proofs are stated in ℝ, with `√2`, `√3`
        txt = txt.replace('s2', '√2').replace('s3', '√3')
    return '(' + txt + ')'


# 1. the multiplication step, checked symbolically in m, n
u, v, m, n = symbols('u v m n')
a1, b1, c1, e1, a2, b2, c2, e2 = symbols('a1 b1 c1 e1 a2 b2 c2 e2')
N1 = a1 + b1*u + c1*v + e1*u*v
N2 = a2 + b2*u + c2*v + e2*u*v
prod = (a1*a2 + m*b1*b2 + n*c1*c2 + m*n*e1*e2) + (a1*b2 + b1*a2 + n*(c1*e2 + e1*c2))*u \
    + (a1*c2 + c1*a2 + m*(b1*e2 + e1*b2))*v + (a1*e2 + e1*a2 + b1*c2 + c1*b2)*u*v
Q0 = b1*b2 + e1*e2*v**2 + v*(b1*e2 + b2*e1)
Q1 = c1*c2 + m*e1*e2 + u*(c1*e2 + c2*e1)
assert expand(N1*N2 - prod - Q0*(u**2 - m) - Q1*(v**2 - n)) == 0
print('exists_int_combo: product formula verified; coefficients', lean(Q0), lean(Q1))

# 2. the uniformizer pi and the basis change
pi = (s2 + s2*s3)/2 - 1
identities = {
    'piL_eisenstein': pi**4 + 4*pi**3 + 2*pi**2 - 4*pi - 2,
    's2_eq': s2 - (pi**3 + 3*pi**2 - 2),
    's3_eq': s3 - (pi**2 + 2*pi - 1),
    's6_eq': s2*s3 - (-pi**3 - 3*pi**2 + 2*pi + 4),
}
for name, E in identities.items():
    q0, q1 = quotients(E, s2, s3, 2, 3)
    print(f'{name}: linear_combination {lean(q0)} * h2 + {lean(q1)} * h3')

# 3. epsilon = 2 / pi^4 in Z[pi] = Z[x]/(g)
g = x**4 + 4*x**3 + 2*x**2 - 4*x - 2
eps = expand(2 * invert(x**4, g, x))
cof, r = div(expand(2 - x**4 * eps), g, x)
assert r == 0 and rem(expand(x**4 * eps - 2), g, x) == 0
print('epsL =', eps, '; 2 - pi^4 eps = (', cof, ') * g(pi)')

# 4. the rhombus chain
data = json.load(open(f'{REPO}/data/chain23.json'))
basis = [1, s2, s3, s2*s3]
pts = [[sum(Rational(p, q) * b for (p, q), b in zip(c, basis)) for c in P]
       for P in data['points']]
for i, j in data['edges']:
    (x1, y1), (x2, y2) = pts[i], pts[j]
    q0, q1 = quotients((x1 - x2)**2 + (y1 - y2)**2 - 1, s2, s3, 2, 3)
    print(f'adj{i}{j}: linear_combination {lean(q0, True)} * sqrt2_sq + {lean(q1, True)} * sqrt3_sq')
