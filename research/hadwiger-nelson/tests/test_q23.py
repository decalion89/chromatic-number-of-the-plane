"""The plane over Q(sqrt2, sqrt3) is 4-colourable and holds a 4-chromatic graph, so its chromatic number is 4.

Voronov (Polymath16, thread 17, July 2021) wrote that chi(Q(i, sqrt3, sqrt11)) = 4 "seems likely, and the
same is true in case (2, 3)", but that nobody had proved it. The first case is tests/test_q311.py.
Here L = Q(sqrt2, sqrt3) has one place over 2, totally ramified with residue field F_2; i is not in L_v,
so the place does not split in K = L(i), and K_w = L_v(w) has residue field F_4. Every unit vector is a
w-adic unit with a nonzero residue in F_4, so the residue of z - rep(z) is a proper 4-colouring
(hn/adelic.py). A chain of three unit rhombi closing at distance 1 gives the lower bound.
"""
import json
import os
import random
from fractions import Fraction as Fr

from pysat.solvers import Solver

from hn.adelic import q23_ab, q23_colour, q23_valuation
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
F = Field((2, 3))
S2, S3, S6 = F.sqrt(2), F.sqrt(3), F.sqrt(6)
ZERO, ONE = F.zero(), F.one()


def _q(x):
    return F.rational(Fr(x))


def _square_class(n):
    """The class of a nonzero integer in Q_2* / squares, as (parity of v_2, odd part mod 8)."""
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return (k % 2, n % 8)


def _mul(a, b):
    return ((a[0] + b[0]) % 2, (a[1] * b[1]) % 8)


R15 = Rotation((S6 + S2) * _q(Fr(1, 4)), (S6 - S2) * _q(Fr(1, 4)))


def _roots():
    """The 24th roots of unity, as unit vectors."""
    out, p = [], Point(ONE, ZERO)
    for _ in range(24):
        out.append(p)
        p = R15(p)
    return out


def _unit_vectors():
    r15 = R15
    base = [Point(ONE, ZERO)]
    for x, y in ((_q(Fr(3, 5)), _q(Fr(4, 5))), (_q(Fr(1, 3)), S2 * _q(Fr(2, 3))),
                 (S3 * _q(Fr(1, 3)), S6 * _q(Fr(1, 3))), (_q(Fr(7, 9)), S2 * _q(Fr(4, 9))),
                 (_q(Fr(1, 7)), S3 * _q(Fr(4, 7))), (_q(Fr(5, 7)), S6 * _q(Fr(2, 7)))):
        base += [Point(x, y), Point(x, -y)]
    out = set()
    for p in base:
        for _ in range(24):
            out.add(p)
            p = r15(p)
    return sorted(out, key=lambda p: (p.fx, p.fy))


def _residue_nonzero(u):
    """a and b of u = a + b w are integral and not both in the maximal ideal (valuation of 0 is inf)."""
    va, vb = (float("inf") if v is None else v for v in map(q23_valuation, q23_ab(u)))
    return min(va, vb) == 0


def _proper(g):
    col = [None] * len(g.vertices)
    for s in range(len(g.vertices)):
        if col[s] is not None:
            continue
        comp, stack = [], [s]
        col[s] = -1
        while stack:
            i = stack.pop()
            comp.append(i)
            for j in g.adj[i]:
                if col[j] is None:
                    col[j] = -1
                    stack.append(j)
        for i in comp:
            col[i] = q23_colour(g.vertices[i], g.vertices[s])
    return all(col[i] != col[j] for i in range(len(col)) for j in g.adj[i])


def _colourable(g, k):
    var = lambda i, c: i * k + c + 1
    with Solver(name="cd19") as s:
        for i in range(len(g.vertices)):
            s.add_clause([var(i, c) for c in range(k)])
            for j in g.adj[i]:
                if i < j:
                    for c in range(k):
                        s.add_clause([-var(i, c), -var(j, c)])
        return s.solve()


def test_one_place_over_2_totally_ramified_and_not_split_in_K():
    two, three, minus_one, five = (_square_class(n) for n in (2, 3, -1, 5))
    group = {(0, 1), two, three, _mul(two, three)}
    # 2, 3 and 6 are independent square classes: L_v = Q_2(sqrt2, sqrt3) has degree 4, one place.
    assert len(group) == 4
    # The class of 5 (= -3) gives the unramified quadratic extension; L_v does not contain it,
    # so v is totally ramified with residue field F_2.
    assert five not in group and _square_class(-3) == five
    # -1 is not in the group: i is not in L_v, so v does not split in K = L(i).
    assert minus_one not in group
    # K_w contains sqrt(-3) = i sqrt3, so it is L_v(sqrt-3): unramified over L_v, residue field F_4.
    assert _mul(minus_one, three) == five


def test_the_valuation_is_the_2_adic_valuation_of_the_norm():
    assert q23_valuation(S2) == 2 and q23_valuation(S3 - ONE) == 2
    assert q23_valuation(S3) == 0 and q23_valuation(_q(2)) == 4


def test_every_listed_unit_vector_has_a_nonzero_residue():
    units = _unit_vectors()
    assert len(units) == 312
    for u in units:
        assert u.x * u.x + u.y * u.y == ONE
        assert _residue_nonzero(u)


def test_random_unit_vectors_have_a_nonzero_residue():
    # Hilbert 90: every norm-one element of K is t / tbar, i.e. ((a^2 - b^2), 2ab) / (a^2 + b^2).
    rng = random.Random(23)
    basis = (ONE, S2, S3, S6)
    for _ in range(300):
        a = sum((_q(Fr(rng.randint(-9, 9), rng.randint(1, 12))) * e for e in basis), ZERO)
        b = sum((_q(Fr(rng.randint(-9, 9), rng.randint(1, 12))) * e for e in basis), ZERO)
        d = a * a + b * b
        if d == ZERO:
            continue
        inv = F.one() / d
        u = Point((a * a - b * b) * inv, _q(2) * a * b * inv)
        assert u.x * u.x + u.y * u.y == ONE
        assert _residue_nonzero(u)


def test_the_residue_colouring_is_proper_on_unit_distance_graphs():
    pts, layer = {Point(ZERO, ZERO)}, {Point(ZERO, ZERO)}
    for _ in range(3):
        layer = {p + u for p in layer for u in _roots()} - pts
        pts |= layer
    g = build_graph(pts)
    assert (g.n, g.m) == (2089, 10296) and _proper(g)
    units = _unit_vectors()[::4]
    mixed = build_graph({Point(ZERO, ZERO)} | set(units) | {u + v for u in units for v in units})
    assert mixed.m > 2000 and _proper(mixed)


def test_a_chain_of_three_rhombi_needs_four_colours():
    d = json.load(open(os.path.join(ROOT, "data", "chain23.json")))
    pts = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(pts)
    assert (g.n, g.m) == (10, 16)
    # O, sqrt3 u1, sqrt3 (u1 + u2), sqrt3 (u1 + u2 + u3) are the chain's joints; the last is at distance 1.
    u3 = Point(_q(Fr(-2, 3)) - S2 * _q(Fr(1, 6)), _q(Fr(-2, 3)) + S2 * _q(Fr(1, 6)))
    assert u3.x * u3.x + u3.y * u3.y == ONE
    end = Point(S3 * (ONE + u3.x), S3 * (ONE + u3.y))
    assert end.x * end.x + end.y * end.y == ONE and end in set(pts)
    assert not _colourable(g, 3)
    assert _colourable(g, 4) and _proper(g)
