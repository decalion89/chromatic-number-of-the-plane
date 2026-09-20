"""Homomorphism colourings: the structural screen, and what it says."""
from fractions import Fraction
from hn.degrey import build_G
from hn.graph import build_graph
from hn.homcol import edge_vectors, has_homomorphism, screen
from hn.mixed import joint_core_union, three_hexagon_gadget


def _verify(phi, vecs, n):
    """Check phi directly rather than trusting the solver."""
    return all(sum(p * x for p, x in zip(phi, d)) % n for d in vecs)


def test_edge_vectors_include_a_real_extension_s_radical_part():
    """`.c` on a RealExtElement gives only its base component.

    Using it alone silently truncated the edge vectors of anything living over
    sqrt(v) and ran the homomorphism search on the wrong module. The joint-core
    union lives over Q(sqrt3, sqrt11)(sqrt v), so its coordinates need both
    halves: dimension 16, not 8.
    """
    g = build_graph(joint_core_union())
    vecs = edge_vectors(g)
    assert len(vecs[0]) == 16


def test_common_content_is_divided_out():
    """Clearing denominators multiplies every vector by one integer, and if
    that integer shares a factor with n the whole set looks divisible by n --
    a scaling artefact reported as an obstruction until the content was
    removed. After removing it, the vectors have no common factor."""
    from math import gcd

    for pts in (three_hexagon_gadget()[2], joint_core_union()):
        vecs = edge_vectors(build_graph(pts))
        content = 0
        for v in vecs:
            for x in v:
                content = gcd(content, abs(x))
        assert content == 1


def test_everything_here_is_five_colourable_by_cosets():
    """A homomorphism from the edge-vector module to Z/5 that avoids every edge
    vector IS a proper 5-colouring, c(p) = phi(p), constant on cosets of its
    kernel. It exists for all three graphs, which is why every rigidity
    measurement in this package came back flat -- a coset colouring composes
    with the automorphisms of Z/5 and with any translation, so colour classes
    move freely and forcing sets are enormous.

    It is also a screen: such a graph can never be 6-chromatic.
    """
    for name, g in (("gadget", build_graph(three_hexagon_gadget()[2])),
                    ("joint core", build_graph(joint_core_union())),
                    ("de Grey G", build_G())):
        vecs = edge_vectors(g)
        phi, why = has_homomorphism(vecs, 5)
        assert phi is not None, f"{name}: {why}"
        assert _verify(phi, vecs, 5), f"{name}: phi does not avoid every vector"
        assert "never 6-chromatic" in screen(g, 5)["verdict"]


def test_rank_two_lattices_can_never_block_a_coset_colouring():
    """A plane lattice's unit vectors occupy exactly three of the six
    projective points mod 5, and which three is decided by the square class of
    their norm. Three hyperplanes cover 13 of the 25 points of (Z/5)^2, so a
    homomorphism always survives -- blocking needs rank at least 4.
    """
    def proj(v, p=5):
        v = (v[0] % p, v[1] % p)
        if v == (0, 0):
            return None
        for lam in range(1, p):
            w = ((v[0] * lam) % p, (v[1] * lam) % p)
            if w[0] == 1 or (w[0] == 0 and w[1] == 1):
                return w
        return v

    def norm(a, b):
        return a * a - a * b + b * b

    squares, nonsquares = set(), set()
    for N in (1, 3, 7, 13, 21, 49, 91):
        vecs = [(a, b) for a in range(-12, 13) for b in range(-12, 13)
                if norm(a, b) == N]
        assert vecs
        cls = {proj(v) for v in vecs} - {None}
        assert len(cls) == 3, f"norm {N} gave {len(cls)} classes"
        (squares if N % 5 in (1, 4) else nonsquares).add(frozenset(cls))
    assert len(squares) == 1 and len(nonsquares) == 1
    assert squares != nonsquares
    assert set(next(iter(squares))) == {(0, 1), (1, 0), (1, 1)}
    assert set(next(iter(nonsquares))) == {(1, 2), (1, 3), (1, 4)}


def test_unit_directions_can_never_contain_a_projective_line():
    """The cheapest blocking set in PG(r-1,5) is a line, six points, and unit
    vectors cannot supply one.

    They all satisfy Q(d) = N for a single N, so their directions carry a form
    value in one square class; and no projective line over F_5 has all six of
    its points in one square class. Checked exhaustively over every binary
    form and every target: the maximum is five, reached only when the form
    degenerates to rank one, the sixth point being its radical.
    """
    p = 5
    squares = {(t * t) % p for t in range(1, p)}

    def direction(v):
        for lam in range(1, p):
            w = tuple((t * lam) % p for t in v)
            for i in range(len(w)):
                if w[i] == 1 and all(w[j] == 0 for j in range(i)):
                    return w
        return None

    best = 0
    for a in range(p):
        for b in range(p):
            for c in range(p):
                pts = {}
                for x in range(p):
                    for y in range(p):
                        if (x, y) == (0, 0):
                            continue
                        d = direction((x, y))
                        if d is None:
                            continue
                        pts.setdefault(d, set()).add(
                            (a * x * x + b * x * y + c * y * y) % p)
                for N in range(1, p):
                    cls = {(N * s) % p for s in squares}
                    best = max(best, sum(1 for vals in pts.values()
                                         if vals & cls))
    assert best == 5, f"a line should never be fully covered; got {best}"


def test_cyclotomic_step_set_admits_no_coset_five_colouring():
    """Q(zeta_7) blocks, where every multiquadratic family fails.

    Multiquadratic fields give at most 3^t directions against PG(2t-1,5)'s
    (5^2t - 1)/4 points, because 5 splits there and the norm-one group mod 5
    is a product of factors of size 2 or 3. A cyclic Galois group with 5
    inert avoids that: in Q(zeta_7), 5 is a primitive root mod 7, O/5 is
    F_5^6, and the norm-one subgroup has 126 elements.

    Modulus-one elements come free from Hilbert 90 -- u = alpha/conj(alpha) --
    and 63 of their directions admit no homomorphism to Z/5 at all. Checked by
    brute force over all 15625 maps rather than on the solver's word.
    """
    import itertools
    import random

    from hn.homcol import CYCLOTOMIC_BLOCKING

    D, N = 6, 7

    def mul(x, y):
        raw = [0] * (2 * D + 1)
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    if b:
                        raw[i + j] += a * b
        for k in range(2 * D, D - 1, -1):
            c = raw[k]
            if c:
                raw[k] = 0
                for j in range(D):
                    raw[k - D + j] -= c
        return tuple(raw[:D])

    def conj(x):
        out = [0] * D
        for i, a in enumerate(x):
            if not a:
                continue
            k = (-i) % N
            if k < D:
                out[k] += a
            else:
                for j in range(D):
                    out[j] -= a
        return tuple(out)

    def inv_mod5(x):
        e, r, b = 5 ** D - 2, tuple([1] + [0] * (D - 1)), tuple(t % 5 for t in x)
        while e:
            if e & 1:
                r = tuple(t % 5 for t in mul(r, b))
            b = tuple(t % 5 for t in mul(b, b))
            e >>= 1
        return r

    def direction(v):
        v = tuple(t % 5 for t in v)
        if not any(v):
            return None
        for lam in range(1, 5):
            w = tuple((t * lam) % 5 for t in v)
            for i in range(D):
                if w[i] == 1 and all(w[j] == 0 for j in range(i)):
                    return w
        return None

    random.seed(7)
    seen, vecs = set(), []
    for _ in range(4000):
        a = tuple(random.randint(-3, 3) for _ in range(D))
        if not any(a):
            continue
        ca = conj(a)
        if all(t % 5 == 0 for t in ca):
            continue
        u = tuple(t % 5 for t in mul(a, inv_mod5(ca)))
        d = direction(u)
        if d and d not in seen:
            seen.add(d)
            vecs.append(u)

    assert len(vecs) == CYCLOTOMIC_BLOCKING["directions"] == 63
    assert not any(all(t % 5 == 0 for t in v) for v in vecs)
    survivors = sum(1 for f in itertools.product(range(5), repeat=D)
                    if all(sum(p * t for p, t in zip(f, v)) % 5 for v in vecs))
    assert survivors == 0, f"{survivors} coset colourings survive"


def test_minimum_blocking_set_finds_the_projective_line():
    """Six slopes at rank 2 cover the dual, and six is the proven minimum."""
    from hn.homcol import has_homomorphism, minimum_blocking_set

    six = [(1, 0), (0, 1), (1, 1), (1, 2), (1, 3), (1, 4)]
    best = minimum_blocking_set(six, 5)
    assert len(best) == 6, "a projective line has no redundant point"
    assert has_homomorphism(best, 5)[0] is None
    # and dropping any one of them lets a colouring back in
    for d in six:
        rest = [e for e in six if e != d]
        assert has_homomorphism(rest, 5)[0] is not None


def test_triangular_lattice_directions_never_cover():
    """Unit vectors of a plane lattice give three of the six slopes."""
    from hn.homcol import minimum_blocking_set

    assert minimum_blocking_set([(1, 0), (0, 1), (1, -1)], 5) == []


def test_blocking_is_monotone_in_the_step_set():
    """Adding a direction can only help: the claim the construction rests on."""
    from hn.homcol import has_homomorphism

    six = [(1, 0), (0, 1), (1, 1), (1, 2), (1, 3), (1, 4)]
    assert has_homomorphism(six, 5)[0] is None
    assert has_homomorphism(six + [(2, 3), (1, 4)], 5)[0] is None


def test_blocking_denominator_threshold_is_the_first_split_prime():
    """29 = 4*7 + 1 is why the directions suddenly cover PG(5,5)."""
    from hn.homcol import BLOCKING_DENOMINATOR

    d = BLOCKING_DENOMINATOR["threshold"]
    assert d % 7 == 1
    assert all(p % 7 == 1 for p in BLOCKING_DENOMINATOR["later_denominators"])
    assert all(q % 7 != 1 for q in range(2, d) if _prime(q))


def _prime(q):
    return q > 1 and all(q % i for i in range(2, int(q ** 0.5) + 1))


def test_unit_distance_triangle_needs_zeta6():
    """Equilateral triangles are 60-degree rotations, so Q(zeta_7) has none."""
    import itertools
    from fractions import Fraction

    from hn.cyclotomic import CycloField

    F = CycloField(7)
    one = F.rational(1)

    # A primitive sixth root of unity is a root of p^2 - p + 1, which pins it
    # down exactly; testing p^6 == 1 does not, since -1 passes that too.
    zeta = F.zeta(1)
    roots, z = [], one
    for _ in range(7):
        z = F.mul(z, zeta)
        roots += [z, F.neg(z)]
    assert one in roots and len(set(roots)) == 14, "the roots of unity are mu_14"
    assert not any(F.add(F.sub(F.mul(p, p), p), one) == F.zero() for p in roots)

    # so no two unit-apart modulus-one elements exist: no triangle through 0
    us = set()
    for c in itertools.product(range(-1, 2), repeat=3):
        a = tuple(Fraction(x) for x in c) + (Fraction(0),) * 3
        if any(a) and F.norm2(a) == one:
            us.add(a)
    assert len(us) > 1
    assert not any(F.norm2(F.sub(u, v)) == one for u in us for v in us if u != v)


def test_eisenstein_lattice_does_have_triangles():
    """The converse half: zeta_6 present, triangles present."""
    from fractions import Fraction

    from hn.field import Field
    from hn.geometry import Point

    fld = Field((3,))
    h = fld.rational(Fraction(1, 2))
    a = Point(fld.zero(), fld.zero())
    b = Point(fld.one(), fld.zero())
    c = Point(h, fld.sqrt(3) * h)       # multiplication by zeta_6
    assert a.is_unit_apart(b) and b.is_unit_apart(c) and a.is_unit_apart(c)


def test_blocking_is_cheap_and_changes_no_chromatic_number():
    """A pendant edge blocks without helping: the caveat, made explicit.

    Six slopes at rank 2 block Z/5.  Hang them off a triangle as degree-one
    vertices and the graph blocks while staying 3-chromatic, because a
    degree-one vertex extends any colouring of the rest.
    """
    from hn.homcol import has_homomorphism

    triangle = [(1, 0), (0, 1), (1, -1)]           # rank 2, does not block
    assert has_homomorphism(triangle, 5)[0] is not None

    pendants = [(1, 1), (1, 2), (1, 3), (1, 4)]    # completes the line
    assert has_homomorphism(triangle + pendants, 5)[0] is None

    # the graph itself: K3 plus one pendant per added direction
    edges = [(0, 1), (1, 2), (0, 2)] + [(0, 3 + i) for i in range(len(pendants))]
    n = 3 + len(pendants)
    from pysat.solvers import Solver
    for k in (2, 3):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in edges:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            assert s.solve() is (k == 3), f"k={k}"


def test_a_plane_lattice_can_never_block():
    """At most three unit directions against PG(1,5)'s six points."""
    import itertools

    from hn.homcol import has_homomorphism

    # PG(1,5) has six points and its hyperplanes are single points, so six
    # directions are needed; a lattice supplies at most three.
    pg1 = [p for p in itertools.product(range(5), repeat=2)
           if any(p) and next(x for x in p if x) == 1]
    assert len(pg1) == 6
    for k in range(1, 6):
        for dirs in itertools.combinations(pg1, k):
            assert has_homomorphism(list(dirs), 5)[0] is not None, (
                f"{k} directions should never block at rank 2")
    assert has_homomorphism(pg1, 5)[0] is None


def test_degree_efficiency_falls_as_blocking_arrives():
    """The measured cost of blocking: forty times the steps for the same ten."""
    from hn.homcol import DEGREE_EFFICIENCY as E

    order = ["triangular patch", "de Grey Sa", "de Grey Y", "de Grey G",
             "Q(zeta_7) ball"]
    vals = [E[k] for k in order]
    assert vals == sorted(vals, reverse=True), "efficiency falls monotonically"
    for name in E["blocked_at"]:
        assert E[name] < 0.05, "blocking only appears at the bottom"
    assert E["de Grey G"] > 3 * E["Q(zeta_7) ball"]


def test_blocking_needs_the_full_rank():
    """Proper submodules hold too few unit vectors to cover anything."""
    from hn.homcol import BLOCKING_NEEDS_FULL_RANK as B

    assert B["blocks_at"] == 6
    assert all(r < B["blocks_at"] for r in B["submodule_ranks_tested"])
    # covering PG(r-1,5) needs at least six directions even in the best case
    for r, n in B["unit_directions_found"].items():
        assert n < (5 ** r - 1) // 4, f"rank {r}: {n} directions is far short"
    assert max(B["unit_directions_found"].values()) < 87


def test_no_eisenstein_spindle_over_zeta21():
    """12m - 1 is 11 mod 12, and neither 3t^2 nor 7t^2 ever is."""
    from math import isqrt

    from hn.homcol import NO_EISENSTEIN_SPINDLE_OVER_ZETA21 as N

    assert {3 * t * t % 12 for t in range(12)} == {0, 3}
    assert 11 not in {7 * t * t % 12 for t in range(12)}
    hits = [m for m in range(1, N["checked_to"] + 1)
            for d in (3, 7)
            if (12 * m - 1) % d == 0
            and isqrt((12 * m - 1) // d) ** 2 == (12 * m - 1) // d]
    assert len(hits) == N["solutions"] == 0
    assert 12 * 1 - 1 == 11, "the classical case is sqrt(-11)"


def test_the_degree_24_field_carries_both_properties():
    from hn.homcol import FIRST_BLOCKED_AND_CHROMATIC as B

    assert B["spindle"] == {"points": 7, "edges": 11, "chi": 4}
    assert B["coset_colourings"] == 0 and B["chi"] == 4
    assert "pendant" in B["honest"], "the blocking is not load bearing"
    assert B["needed_forty_vs"] == 144, "forty steps left a coset colouring"


def test_coset_colouring_is_a_colouring():
    """The gate: a coset n-colouring really is a proper n-colouring.

    This is the whole reason blocking matters -- it makes chi >= 6 imply
    blocked -- so it is walked out on a real graph rather than asserted.
    Colours are propagated along a spanning tree by phi of the edge vector;
    that the result is consistent around every cycle and proper on every
    non-tree edge is exactly the claim.
    """
    from collections import deque
    from fractions import Fraction
    from math import gcd

    from hn.degrey import build_Sa
    from hn.graph import build_graph
    from hn.homcol import _coords, edge_vectors, has_homomorphism

    g = build_graph(build_Sa())
    edges = list(g.edges())
    vecs = edge_vectors(g)
    phi, _ = has_homomorphism(vecs, 5)
    assert phi is not None, "Sa is multiquadratic, so it cannot block"

    raw = {}
    for a, b in edges:
        d = (g.vertices[b].x - g.vertices[a].x, g.vertices[b].y - g.vertices[a].y)
        raw[(a, b)] = _coords(d[0]) + _coords(d[1])
    den = 1
    for v in raw.values():
        for q in v:
            den = den * Fraction(q).denominator // gcd(den, Fraction(q).denominator)
    ints = {e: tuple(int(Fraction(q) * den) for q in v) for e, v in raw.items()}
    content = 0
    for v in ints.values():
        for x in v:
            content = gcd(content, abs(x))
    if content > 1:
        ints = {e: tuple(x // content for x in v) for e, v in ints.items()}

    def step(e):
        return sum(c * x for c, x in zip(phi, ints[e])) % 5

    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    colour = {0: 0}
    q = deque([0])
    while q:
        v = q.popleft()
        for w in adj.get(v, ()):
            if w in colour:
                continue
            colour[w] = (colour[v] + (step((v, w)) if (v, w) in ints
                                      else -step((w, v)))) % 5
            q.append(w)

    assert len(colour) == len(g.vertices), "Sa should be connected"
    for a, b in edges:
        assert colour[a] != colour[b], (a, b)


def test_denominator_29_directions_block():
    from hn.homcol import denominator_29_directions, has_homomorphism

    d = denominator_29_directions()
    assert len(d) == 300 and len(set(d)) == 300
    assert has_homomorphism(d, 5)[0] is None, "the 300 must block"


def test_blocking_is_monotone_in_the_direction_set():
    """Used to read U's blocking off a subset of its directions."""
    from hn.homcol import denominator_29_directions, has_homomorphism

    d = denominator_29_directions()
    assert has_homomorphism(d[:6], 5)[0] is not None
    assert has_homomorphism(d, 5)[0] is None


def test_every_unit_triangle_is_a_sixty_degree_pair():
    """|u| = |v| = |u-v| = 1 forces u/v to be a primitive sixth root."""
    import itertools
    from fractions import Fraction
    from hn.cyclotomic import CycloField

    K = CycloField(21)
    one = K.rational(1)
    z6 = K.neg(K.mul(K.zeta(7), K.zeta(7)))
    units = []
    for k in range(6):
        p, z = K.rational(1), one
        for _ in range(k):
            p = K.mul(p, z6)
        units.append(p)
    for u, v in itertools.permutations(units, 2):
        if K.norm2(K.sub(u, v)) != one:
            continue
        t = K.mul(u, K.conj(v))
        assert K.add(t, K.conj(t)) == one
        assert K.norm2(t) == one


def test_y_keeps_the_hexagonal_edges_at_the_origin():
    """What lets U inherit a blocking subset from G."""
    from hn.degrey import build_Y
    from hn.geometry import DEGREY_FIELD
    from hn.graph import build_graph

    Y = build_Y()
    g = build_graph(Y)
    zero = DEGREY_FIELD.zero()
    o = next(i for i, p in enumerate(Y) if p.x == zero and p.y == zero)
    nb = [j for i, j in g.edges() if i == o] + [i for i, j in g.edges() if j == o]
    half, root3 = DEGREY_FIELD.rational(Fraction(1, 2)), DEGREY_FIELD.sqrt(3)
    want = {(DEGREY_FIELD.rational(1), zero),
            (DEGREY_FIELD.rational(-1), zero)}
    for sx in (1, -1):
        for sy in (1, -1):
            want.add((DEGREY_FIELD.rational(Fraction(sx, 2)),
                      root3 * DEGREY_FIELD.rational(Fraction(sy, 2))))
    got = {(Y[j].x - Y[o].x, Y[j].y - Y[o].y) for j in nb}
    assert want <= got, "Y must keep all six hexagonal unit steps"
    assert len(nb) == 60


def test_moser_spindle_does_not_block():
    """One measurement settles every 7-vertex 4-critical graph.

    The spindle's arm rotation is forced to (5 +- sqrt-11)/6 and blocking is
    invariant under a global rotation, so the class has a single member up to
    isometry.
    """
    from hn.geometry import DEGREY_FIELD, Point, Rotation
    from hn.graph import build_graph
    from hn.homcol import edge_vectors, has_homomorphism

    F = DEGREY_FIELD
    half = F.rational(Fraction(1, 2))
    rhombus = [Point(F.zero(), F.zero()), Point(F.rational(1), F.zero()),
               Point(half, F.sqrt(3) * half),
               Point(F.rational(Fraction(3, 2)), F.sqrt(3) * half)]
    rho = Rotation(F.rational(Fraction(5, 6)),
                   F.sqrt(11) * F.rational(Fraction(1, 6)))
    g = build_graph(rhombus + [rho(q) for q in rhombus[1:]])
    vecs = edge_vectors(g)
    assert len(g.vertices) == 7 and len(list(g.edges())) == 11
    assert has_homomorphism(vecs, 5)[0] is not None
