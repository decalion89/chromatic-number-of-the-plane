"""The planes over Q(sqrt3, sqrt q): chi = 3 if q = 1 (mod 3), and chi >= 4 if q = 2 (mod 3).

Lemma. If sqrt3 is in L and 1/3 is a sum of unit vectors of L(i), then L^2 has no proper
3-colouring: the sums of unit vectors form a ring C_0 containing w = e^{i pi/3}, so (1 + w)/3 is a
sum of unit vectors u_k, and the unit rhombi with long diagonals sqrt3 u_k chain 0 to
sqrt3 (1 + w)/3 = e^{i pi/6}, a unit vector, forcing it to share the colour of 0.

For squarefree q = 2 (mod 3), 3 splits in Q(sqrt -q). A generator alpha = a + b sqrt(-q) of the
h-th power of a prime above 3 (h the class number) has norm 3^h, and 3 does not divide
Tr(alpha^2); otherwise the prime would divide conj(alpha). So r = alpha/conj(alpha) + conj(alpha)/alpha
= Tr(alpha^2)/3^h is a sum of two unit vectors with exact denominator 3^h, and 1/3 lies in C_0.

Together with reduction at 3 (Madore, Cor. 3.4) for q = 1 (mod 3), and with Theorem 2 of the note
for q = 2 or q = 1, 3 (mod 8): chi(Q(sqrt3, sqrt q)^2) = 3 if q = 1 (mod 3), = 4 if q = 11, 17
(mod 24) or q even, and >= 4 if q = 5, 23 (mod 24). notes/local_colourings.md, section 11.
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pysat.solvers import Solver

from hn.certify import exact_edges
from hn.field import Field
from hn.geometry import Point

# (q, class number h of Q(sqrt -q), a, b) with alpha = a + b sqrt(-q) generating P^h, P above 3.
# Computed with PARI/GP (bnfinit, bnfisprincipal); the tests below check the arithmetic.
ALPHAS = [
    (2, 1, Fr(-1), Fr(1)), (5, 2, Fr(-2), Fr(-1)), (11, 1, Fr(1, 2), Fr(1, 2)),
    (17, 4, Fr(-8), Fr(1)), (23, 3, Fr(-2), Fr(1)), (29, 6, Fr(-2), Fr(5)),
    (41, 8, Fr(40), Fr(11)), (47, 5, Fr(14), Fr(-1)), (53, 6, Fr(-26), Fr(1)),
    (59, 3, Fr(-7, 2), Fr(-1, 2)), (71, 7, Fr(-46), Fr(-1)), (83, 3, Fr(-5, 2), Fr(1, 2)),
    (89, 12, Fr(640), Fr(37)), (101, 14, Fr(-338), Fr(215)), (107, 3, Fr(1, 2), Fr(1, 2)),
    (113, 8, Fr(-32), Fr(-7)),
]


def test_every_listed_q_puts_one_third_among_sums_of_unit_vectors():
    for q, h, a, b in ALPHAS:
        assert q % 3 == 2
        norm = a * a + q * b * b
        assert norm == 3 ** h, q
        trace = 2 * (a * a - q * b * b)          # Tr(alpha^2), an integer
        assert trace.denominator == 1 and trace.numerator % 3 != 0, q
        # alpha/conj(alpha) = ((a^2 - q b^2) + 2ab sqrt(-q))/N is a unit vector of Q(sqrt q)(i)
        x, y2 = (a * a - q * b * b) / norm, (2 * a * b) ** 2 * q / norm ** 2
        assert x * x + y2 == 1, q
        # r = 2x = trace/3^h has exact denominator 3^h, so 1/3^h and hence 1/3 lie in C_0
        r = 2 * x
        assert r.denominator == 3 ** h


def test_q_17_a_chain_of_ninety_rhombi_has_no_proper_3_colouring():
    """u = (8 + i sqrt17)/9 gives 16/9 = u + conj(u) and 1/3 = 12 (u + conj(u)) - 21."""
    F = Field((3, 17))
    r = F.rational
    s3, s17 = F.sqrt(3), F.sqrt(17)
    P = Point
    mul = lambda p, t: P(p.x * t.x - p.y * t.y, p.x * t.y + p.y * t.x)
    u = P(r(Fr(8, 9)), s17 * r(Fr(1, 9)))
    ub = P(u.x, -u.y)
    w = P(r(Fr(1, 2)), s3 * r(Fr(1, 2)))
    one = P(r(1), r(0))
    third = [u] * 12 + [ub] * 12 + [P(-one.x, -one.y)] * 21
    steps = third + [mul(t, w) for t in third]            # (1 + w)/3 as 90 unit vectors
    pts = [P(r(0), r(0))]
    cur = pts[0]
    for t in steps:
        assert t.x * t.x + t.y * t.y == F.one()
        nxt = P(cur.x + s3 * t.x, cur.y + s3 * t.y)
        half = r(Fr(1, 2))
        mid = P((cur.x + nxt.x) * half, (cur.y + nxt.y) * half)
        pts += [P(mid.x - t.y * half, mid.y + t.x * half), P(mid.x + t.y * half, mid.y - t.x * half), nxt]
        cur = nxt
    assert cur.x * cur.x + cur.y * cur.y == F.one()        # the chain ends at a unit vector
    V = list(dict.fromkeys(pts))
    E = exact_edges(V)
    k = 3
    X = lambda v, c: 1 + v * k + c
    s = Solver(name="cd19")
    for v in range(len(V)):
        s.add_clause([X(v, c) for c in range(k)])
    for a_, b_ in E:
        for c in range(k):
            s.add_clause([-X(a_, c), -X(b_, c)])
    assert not s.solve()
    s.delete()


def test_q_1_mod_3_has_a_residue_field_f3_above_3():
    for q in (7, 13, 19, 31, 37, 43, 61, 67, 73):
        assert pow(q, (3 - 1) // 2, 3) == 1                 # 3 splits in Q(sqrt q)
