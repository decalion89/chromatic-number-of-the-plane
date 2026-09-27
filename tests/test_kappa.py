"""kappa = omega * rho7 = (-11 + 5 sqrt-3)/14 is congruent to 1 mod 5.

rho7 = (1 + 4 sqrt-3)/7 is the rotation of Q(sqrt-3) at the prime 7; mod 5 it
is omega^-1, so kappa = omega * rho7 is an irrational rotation with
kappa - 1 = 5 (omega - 3)/7.  On the module of five_rho7 every unit vector u
whose kappa-image is again a unit vector has kappa u - u in 5M, so no
homomorphism M -> Z/5 can tell u from kappa u -- while omega and rho7 alone
move some unit off its class.
"""
import os
import json
from fractions import Fraction as Fr
from math import gcd

from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _module(units):
    """Integer echelon basis of the Z-module spanned by `units`; returns vec, coords."""
    raw = [tuple(u.x.c) + tuple(u.y.c) for u in units]
    den = 1
    for v in raw:
        for q in v:
            den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
    vec = lambda p: [int(Fr(q) * den) for q in tuple(p.x.c) + tuple(p.y.c)]
    basis = []
    for u in units:
        v = vec(u)
        while any(v):
            p = next(i for i, t in enumerate(v) if t)
            row = next((r for r in basis if next(i for i, t in enumerate(r) if t) == p), None)
            if row is None:
                basis.append(v if v[p] > 0 else [-t for t in v])
                break
            a, b = row[p], v[p]
            x0, x1, y0, y1, aa, bb = 1, 0, 0, 1, a, b
            while bb:
                q = aa // bb
                aa, bb = bb, aa - q * bb
                x0, x1, y0, y1 = x1, x0 - q * x1, y1, y0 - q * y1
            new = [x0 * r + y0 * w for r, w in zip(row, v)]
            v = [(a // aa) * w - (b // aa) * r for r, w in zip(row, v)]
            row[:] = new if new[p] > 0 else [-t for t in new]
    basis.sort(key=lambda r: next(i for i, t in enumerate(r) if t))

    def coords(v):
        v, c = list(v), []
        for r in basis:
            p = next(i for i, t in enumerate(r) if t)
            assert v[p] % r[p] == 0
            k = v[p] // r[p]
            c.append(k)
            v = [a - k * b for a, b in zip(v, r)]
        assert not any(v)
        return c
    return vec, coords, len(basis)


def test_kappa_arithmetic():
    F = Field((3,))
    r3 = F.sqrt(3)
    half, q = F.rational(Fr(1, 2)), lambda a, b: F.rational(Fr(a, b))
    omega = Rotation(half, r3 * half)
    rho7 = Rotation(q(1, 7), r3 * q(4, 7))
    kappa = omega * rho7
    assert kappa == Rotation(q(-11, 14), r3 * q(5, 14))
    # kappa - 1 = 5 (omega - 3) / 7, coordinate by coordinate
    assert kappa.cos - F.one() == q(5, 7) * (half - q(3, 1))
    assert kappa.sin == q(5, 7) * (r3 * half)
    # not a root of unity: cos = -11/14 is not in {0, +-1/2, +-1}
    assert kappa.cos not in {q(0, 1), q(1, 2), q(-1, 2), q(1, 1), q(-1, 1)}


def test_kappa_is_invisible_mod_5_on_five_rho7():
    d = json.load(open(f"{ROOT}/data/five_rho7.json"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
         for x, y in d["points"]]
    g = build_graph(P)
    units = list({u: 1 for a, b in g.edges() for u in (P[b] - P[a], P[a] - P[b])})
    vec, coords, rank = _module(units)
    assert (len(units), rank) == (186, 8)
    r3 = F.sqrt(3)
    q = lambda a, b: F.rational(Fr(a, b))
    rots = {"kappa": Rotation(q(-11, 14), r3 * q(5, 14)),
            "omega": Rotation(q(1, 2), r3 * q(1, 2)),
            "rho7": Rotation(q(1, 7), r3 * q(4, 7))}
    uset = set(units)
    moved = {}
    for name, R in rots.items():
        pairs = [(u, R(u)) for u in units if R(u) in uset]
        moved[name] = (len(pairs), sum(1 for u, w in pairs
                                       if any(c % 5 for c in coords([a - b for a, b in zip(vec(w), vec(u))]))))
    assert moved["kappa"] == (120, 0)          # every kappa u - u lies in 5M
    assert moved["omega"][1] > 0 and moved["rho7"][1] > 0


def test_stiemke_certificates_on_five_rho7_sample():
    """For a few admissible psi and every t, the class D_t has an exact strictly positive
    dependency (so no twisted colouring), recomputed here as in scripts/stiemke.py."""
    import itertools
    import numpy as np
    from scipy.optimize import linprog
    d = json.load(open(f"{ROOT}/data/five_rho7.json"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
         for x, y in d["points"]]
    g = build_graph(P)
    units = list({u: 1 for a, b in g.edges() for u in (P[b] - P[a], P[a] - P[b])})
    vec, coords, rank = _module(units)
    Cm5 = [[c % 5 for c in coords(vec(u))] for u in units]
    raw = [vec(u) for u in units]                          # small integer field coordinates
    adm = [psi for psi in itertools.product(range(5), repeat=rank)
           if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in Cm5)]
    assert len(adm) == 960
    for psi in adm[:: 240]:                                # four psi, all t
        val = [sum(a * b for a, b in zip(psi, c)) % 5 for c in Cm5]
        for t in range(1, 5):
            D = [raw[k] for k in range(len(units)) if val[k] == t]
            A = np.array(D, dtype=float).T
            res = linprog(np.zeros(len(D)), A_eq=A, b_eq=np.zeros(len(raw[0])),
                          bounds=[(1.0, None)] * len(D), method="highs")
            assert res.status == 0
            lam = res.x
            # exact: fix all but a basis, solve the basis with Fractions, check positivity
            import sympy
            order = sorted(range(len(D)), key=lambda q: -lam[q])
            basis = []
            for q in order:
                if sympy.Matrix([D[j] for j in basis + [q]]).rank() > len(basis):
                    basis.append(q)
                if len(basis) == rank:
                    break
            rest = [q for q in range(len(D)) if q not in basis]
            lr = {q: Fr(float(lam[q])).limit_denominator(10 ** 6) for q in rest}
            rhs = sympy.Matrix([-sum(lr[q] * D[q][c] for q in rest) for c in range(len(raw[0]))])
            Mb = sympy.Matrix([[D[q][c] for q in basis] for c in range(len(raw[0]))])
            sol, params = Mb.gauss_jordan_solve(rhs)
            assert not params
            full = {**lr, **{q: Fr(str(v)) for q, v in zip(basis, list(sol))}}
            assert all(v > 0 for v in full.values())
            assert all(sum(full[q] * D[q][c] for q in range(len(D))) == 0 for c in range(len(raw[0])))
