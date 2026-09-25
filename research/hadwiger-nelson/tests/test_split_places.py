"""Split places colour the edges they see (notes/local_colourings.md, section 12).

Let L be totally real and P a prime of L(i) over a prime p of L that splits in L(i), with residue field
F_q.  For a unit vector u we have u * conj(u) = 1, so v_P(u) = -v_Pbar(u).  If every edge vector of a
unit-distance graph G in L^2 is a unit at P, then z -> (z mod P, conj(z) mod P) sends G into the
Cayley graph H_q on F_q x F_q with connection set {(t, 1/t) : t in F_q^*}, and chi(G) <= chi(H_q).
Points need not be integral: colour z through its class modulo the ring of P- and Pbar-integral
elements, since every edge vector lies in that ring.

Over L = Q(sqrt3, sqrt5) the primes above 2 and 3 split in L(i), with residue fields F_4 and F_9, and
chi(H_4) = 4, chi(H_9) = 3.  So a unit-distance graph over L with no proper 3-colouring has an edge
vector that is not a unit above 3 (in data/chain35.json it is tau = (2 + i sqrt5)/3), and one with no
proper 4-colouring must also have an edge vector that is not a unit above 2.

The reductions are written on the integral basis zeta^j phi^f (j < 4, f < 2) of Z[zeta12, phi], the
ring of integers of L(i) = Q(zeta12, sqrt5), zeta = e^{i pi/6}, phi = (1 + sqrt5)/2.
"""
import sys, os, json, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pysat.solvers import Solver

from hn.certify import exact_edges
from hn.fast import IntBasis, fast_edges_complete
from hn.field import Field
from hn.geometry import Point

F = Field((3, 5))
r = F.rational
S3, S5, S15 = F.sqrt(3), F.sqrt(5), F.sqrt(15)
DATA = os.path.join(os.path.dirname(__file__), "..", "data", "chain35.json")


# ---- finite fields F_q (q = p or p^2, p^3) as polynomials over F_p --------------------------------
def gf(q):
    """(elements, add, mul, zero, one) for F_q, q in {2, 3, 4, 5, 7, 8, 9, 11, 13, 16}"""
    mods = {4: (2, (1, 1, 1)), 8: (2, (1, 1, 0, 1)), 9: (3, (1, 0, 1)), 16: (2, (1, 1, 0, 0, 1))}
    # x^2+x+1, x^3+x+1, x^2+1, x^4+x+1
    if q not in mods:
        els = list(range(q))
        return els, (lambda a, b: (a + b) % q), (lambda a, b: a * b % q), 0, 1
    p, poly = mods[q]
    n = len(poly) - 1

    def mulp(a, b):
        out = [0] * (2 * n - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
        for d in range(2 * n - 2, n - 1, -1):
            c = out[d]
            if c:
                for k in range(n + 1):
                    out[d - n + k] = (out[d - n + k] - c * poly[k]) % p
        return tuple(out[:n])

    els = list(itertools.product(range(p), repeat=n))
    return els, (lambda a, b: tuple((x + y) % p for x, y in zip(a, b))), mulp, \
        tuple([0] * n), tuple([1] + [0] * (n - 1))


def hyperbola_graph(q):
    els, add, mul, zero, one = gf(q)
    nz = [e for e in els if e != zero]
    inv = {a: next(b for b in nz if mul(a, b) == one) for a in nz}
    V = [(a, b) for a in els for b in els]
    idx = {v: i for i, v in enumerate(V)}
    E = set()
    for a, b in V:
        for t in nz:
            j = idx[(add(a, t), add(b, inv[t]))]
            i = idx[(a, b)]
            E.add((min(i, j), max(i, j)))
    return V, idx, sorted(E)


def clique_to_pin(n, edges):
    """a triangle through vertex 0 if there is one, else an edge, else vertex 0: pinning its colours
    breaks the symmetry of permuting colours"""
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    for b in sorted(adj[0]):
        for c in sorted(adj[0] & adj[b]):
            return [0, b, c]
    return [0, min(adj[0])] if adj[0] else [0]


def colouring(n, edges, k):
    X = lambda v, c: 1 + v * k + c
    s = Solver(name="cd19")
    for v in range(n):
        s.add_clause([X(v, c) for c in range(k)])
    for a, b in edges:
        for c in range(k):
            s.add_clause([-X(a, c), -X(b, c)])
    pin = clique_to_pin(n, edges)
    if len(pin) > k:
        s.delete()
        return None
    for c, v in enumerate(pin):
        s.add_clause([X(v, c)])
    if not s.solve():
        s.delete()
        return None
    pos = set(l for l in s.get_model() if l > 0)
    s.delete()
    return [next(c for c in range(k) if X(v, c) in pos) for v in range(n)]


def chromatic(n, edges, kmax=6):
    return next(k for k in range(1, kmax + 1) if colouring(n, edges, k) is not None)


def test_chromatic_numbers_of_the_hyperbola_graphs():
    expect = {2: 2, 3: 3, 4: 4, 5: 3, 7: 4, 8: 4, 9: 3, 11: 4, 16: 4}
    for q, chi in expect.items():
        V, _, E = hyperbola_graph(q)
        assert chromatic(len(V), E) == chi, q


def test_h13_is_5_chromatic():
    """A proper 5-colouring of H_13 is linear: (a, b) -> col[a + 2b mod 13], since {t + 2/t} =
    {2, 3, 5, 8, 10, 11} and col is a proper colouring of that circulant.  CDCL finds no 4-colouring
    once a triangle is pinned (finding a 5-colouring by search is slow, so it is written down)."""
    V, idx, E = hyperbola_graph(13)
    col = [4, 4, 3, 3, 2, 2, 1, 4, 0, 3, 1, 2, 0]
    c = [col[(a + 2 * b) % 13] for a, b in V]
    assert all(c[u] != c[v] for u, v in E)
    assert colouring(len(V), E, 4) is None


def test_the_f17_plane_has_an_interval_6_colouring():
    """Moorhouse's Table 6.1 gives 5, 6 or 7 for F_17^2 ((x - x')^2 + (y - y')^2 = 1).  The line 3x + 6y = r
    meets the unit circle iff 45 - r^2 is a square, and 11, 10, 7 are not squares mod 17, so on the circle
    3x + 6y avoids 0, +-1 and +-2.  Blocks of three consecutive parallel lines are therefore independent:
    (x, y) -> floor(((3x + 6y) mod 17) / 3) is a proper 6-colouring."""
    q = 17
    squares = {t * t % q for t in range(q)}
    assert all((45 - r * r) % q not in squares for r in range(3))
    U = [(a, b) for a in range(q) for b in range(q) if (a * a + b * b) % q == 1]
    assert len(U) == 16 and not {(3 * a + 6 * b) % q for a, b in U} & {0, 1, 2, 15, 16}
    c = lambda x, y: (3 * x + 6 * y) % q // 3
    assert {c(x, y) for x in range(q) for y in range(q)} == set(range(6))
    assert all(c(x, y) != c((x + a) % q, (y + b) % q) for x in range(q) for y in range(q) for a, b in U)


# ---- reduction of z = x + i y, (x, y) in L^2, on the integral basis of L(i) ----------------------
# images of 1, sqrt3, sqrt5, sqrt15 and of i, i sqrt3, i sqrt5, i sqrt15 on the basis zeta^j phi^f,
# indexed j + 4 f; zeta^4 = zeta^2 - 1, phi^2 = phi + 1, i = zeta^3, sqrt3 = 2 zeta - zeta^3,
# sqrt5 = 2 phi - 1
REAL = [{0: 1}, {1: 2, 3: -1}, {4: 2, 0: -1}, {5: 4, 1: -2, 7: -2, 3: 1}]
IMAG = [{3: 1}, {2: 2, 0: -1}, {7: 2, 3: -1}, {6: 4, 2: -2, 4: -2, 0: 1}]


def obasis(p):
    """the 8 rational coordinates of x + i y on the basis zeta^j phi^f"""
    c = [Fr(0)] * 8
    for coeff, img in zip(p.x.c, REAL):
        for k, v in img.items():
            c[k] += coeff * v
    for coeff, img in zip(p.y.c, IMAG):
        for k, v in img.items():
            c[k] += coeff * v
    return c


def reducer(q, zeta_img, phi_img):
    """the ring map Z[zeta12, phi]_(p) -> F_q with zeta -> zeta_img, phi -> phi_img"""
    els, add, mul, zero, one = gf(q)
    p = q if isinstance(zero, int) else {4: 2, 8: 2, 9: 3, 16: 2}[q]
    powers = []
    for f in range(2):
        for j in range(4):
            x = one
            for _ in range(j):
                x = mul(x, zeta_img)
            for _ in range(f):
                x = mul(x, phi_img)
            powers.append(x)
    order = [j + 4 * f for f in range(2) for j in range(4)]
    base = dict(zip(order, powers))

    def scal(n, x):
        out = zero
        for _ in range(n % p):
            out = add(out, x)
        return out

    def h(point):
        out = zero
        for k, c in enumerate(obasis(point)):
            assert c.denominator % p, "not integral at the primes above p"
            out = add(out, scal(c.numerator * pow(c.denominator, -1, p), base[k]))
        return out
    return h, (els, add, mul, zero, one)


def load(path):
    d = json.load(open(path))
    fld = Field(tuple(d["field_generators"]))
    mk = lambda xy: Point(fld.element([Fr(a, b) for a, b in xy[0]]),
                          fld.element([Fr(a, b) for a, b in xy[1]]))
    return [mk(xy) for xy in d["points"]]


def is_p_unit(u, p):
    """both u and 1/u = conj(u) have coordinates prime to p on the integral basis"""
    return all(c.denominator % p for c in obasis(u)) and \
        all(c.denominator % p for c in obasis(Point(u.x, -u.y)))


def test_the_2_adic_reduction_four_colours_chain35():
    """chain35 has no proper 3-colouring and all its edge vectors are units above 2: z -> z mod P,
    a value in F_4, is a proper 4-colouring of it (and chi(H_4) = 4 is attained)."""
    pts = load(DATA)
    E = exact_edges(pts)
    assert colouring(len(pts), E, 3) is None
    w = (0, 1)                                   # a primitive cube root of unity in F_4 = F_2[x]/(x^2+x+1)
    h, _ = reducer(4, w, w)
    for a, b in E:
        u = Point(pts[b].x - pts[a].x, pts[b].y - pts[a].y)
        assert is_p_unit(u, 2)
    col = [h(p) for p in pts]
    assert all(col[a] != col[b] for a, b in E)
    assert len(set(col)) == 4


def test_chain35_has_an_edge_vector_that_is_not_a_unit_above_3():
    """chi(H_9) = 3, so the 3-chromatic bound would apply if every edge vector were a unit above 3."""
    pts = load(DATA)
    E = exact_edges(pts)
    deep = [(a, b) for a, b in E
            if not is_p_unit(Point(pts[b].x - pts[a].x, pts[b].y - pts[a].y), 3)]
    assert deep
    tau = Point(r(Fr(2, 3)), S5 * r(Fr(1, 3)))
    assert not is_p_unit(tau, 3) and is_p_unit(tau, 2)


def test_the_3_adic_reduction_three_colours_a_graph_whose_edges_are_units_above_3():
    """U = the zeta-orbit of 1, o1 and conj(o1), o1 = ((sqrt5 - sqrt3) - i (sqrt3 + sqrt5))/4 (a unit
    vector with odd valuation above 2 and a unit above 3).  G = {0} u U u (U + U) has triangles, and
    z -> (z mod P3, conj(z) mod P3) followed by a proper 3-colouring of H_9 is a proper 3-colouring."""
    zeta = Point(S3 * r(Fr(1, 2)), r(Fr(1, 2)))
    o1 = Point(S5 * r(Fr(1, 4)) - S3 * r(Fr(1, 4)), -(S3 * r(Fr(1, 4))) - S5 * r(Fr(1, 4)))
    mul = lambda a, b: Point(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x)
    base = [Point(r(1), r(0)), o1, Point(o1.x, -o1.y)]
    U = set()
    for b in base:
        p = b
        for _ in range(12):
            U.add(p)
            p = mul(p, zeta)
    U = sorted(U, key=lambda p: (float(p.x), float(p.y)))
    for u in U:
        assert u.x * u.x + u.y * u.y == F.one() and is_p_unit(u, 3)
    O = Point(F.zero(), F.zero())
    pts = [O] + U
    seen = set(pts)
    for i, u in enumerate(U):
        for v in U[i:]:
            s = Point(u.x + v.x, u.y + v.y)
            if s not in seen:
                seen.add(s)
                pts.append(s)
    basis = IntBasis.covering(pts)
    E = fast_edges_complete(basis, basis.rows(pts))                # every unit pair, checked exactly
    assert colouring(len(pts), E, 2) is None                       # triangles
    j = (0, 1)                                                     # j^2 = -1 in F_9 = F_3[x]/(x^2+1)
    h, (els, add, mulq, zero, one) = reducer(9, j, (2, 2))         # phi -> (1 + j)/2 = 2 + 2j
    V, idx, EH = hyperbola_graph(9)
    colH = colouring(len(V), EH, 3)
    col = [colH[idx[(h(p), h(Point(p.x, -p.y)))]] for p in pts]
    assert all(col[a] != col[b] for a, b in E)


def test_a_unit_that_is_not_a_unit_above_2():
    """phi2 = (1 + i sqrt15)/4: phi2 + conj(phi2) = 1/2, and phi2 is not a unit above 2."""
    phi2 = Point(r(Fr(1, 4)), S15 * r(Fr(1, 4)))
    assert phi2.x * phi2.x + phi2.y * phi2.y == F.one()
    assert phi2.x + phi2.x == r(Fr(1, 2))
    assert not is_p_unit(phi2, 2) and is_p_unit(phi2, 3)
