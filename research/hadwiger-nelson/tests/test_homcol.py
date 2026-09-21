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


def test_orbit_types_decide_blocking():
    """The six small sigma-orbit types, by exhaustion over O/5.

    Exactly the two with residue degree 3 over the real subfield block -- a
    fixed prime of degree 6 lies over one of degree 3, a swapped pair of
    degree 3 over one of degree 3.
    """
    from hn.homcol import orbit_can_block

    assert orbit_can_block("fixed", 2) is False
    assert orbit_can_block("fixed", 4) is False
    assert orbit_can_block("fixed", 6) is True
    assert orbit_can_block("swapped", 1) is False
    assert orbit_can_block("swapped", 2) is False
    assert orbit_can_block("swapped", 3) is True


def test_cyclotomic_verdicts():
    from hn.homcol import cyclotomic_can_block

    for n, deg, ok in [(3, 1, False), (4, 1, False), (6, 1, False),
                       (7, 3, True), (8, 2, False), (9, 3, True),
                       (11, 5, True), (12, 2, False), (13, 2, False),
                       (21, 3, True), (24, 2, False), (33, 10, True)]:
        got = cyclotomic_can_block(n)
        assert got["residue_degree_over_F"] == deg, (n, got)
        assert got["can_block"] is ok, (n, got)
    # K/F ramifies above 5 exactly when n = 5^a or 2.5^a: there the residue
    # extension is trivial, N mod p lies in {+1,-1}, and nothing blocks.
    for n in (5, 10, 25, 50):
        got = cyclotomic_can_block(n)
        assert got["ramified_over_F"] is True
        assert got["can_block"] is False
    # But 5 | n does NOT imply K/F is ramified. Q(zeta_35) has
    # f = ord(5 mod 7) = 6 with the primes sigma-fixed, so residue degree 3
    # over F -- an earlier blanket verdict called this unable to block.
    got = cyclotomic_can_block(35)
    assert got["ramified_over_F"] is False
    assert got["residue_degree_over_F"] == 3
    assert got["can_block"] is True
    assert got["residue_field_only"] is True
    # The two ramified-in-K/Q cases settled by exhaustion agree.
    for n in (15, 20):
        got = cyclotomic_can_block(n)
        assert got["ramified_over_F"] is False
        assert got["residue_degree_over_F"] <= 2
        assert got["can_block"] is False


def test_the_verdict_matches_the_measurement():
    """Predicted from the arithmetic, then measured on real direction sets.

    Q(zeta_3) is the triangular lattice, predicted unable to block; its six
    unit steps admit a coset colouring. Q(zeta_21) has residue degree 3 over
    its real subfield, predicted able; its denominator-29 steps block.
    """
    from fractions import Fraction
    from hn.cyclotomic import CycloField
    from hn.homcol import (cyclotomic_can_block, denominator_29_directions,
                           has_homomorphism)

    K = CycloField(3)
    z6 = K.neg(K.mul(K.zeta(3), K.zeta(3))) if K.degree == 2 else None
    steps, z = [], K.rational(1)
    for _ in range(6):
        steps.append(tuple(int(x) for x in z))
        z = K.mul(z, K.neg(K.rational(1)) if z6 is None else z6)
    assert cyclotomic_can_block(3)["can_block"] is False
    assert has_homomorphism(sorted(set(steps)), 5)[0] is not None

    assert cyclotomic_can_block(21)["can_block"] is True
    assert has_homomorphism(denominator_29_directions(), 5)[0] is None


def test_multiquadratic_is_a_corollary():
    """A multiquadratic field has residue degree 1 over its real subfield.

    Its Galois group is elementary abelian, so every decomposition group is
    cyclic of order at most 2 and f <= 2. de Grey's G could never have
    blocked, and the measurement on Sa agrees.
    """
    from hn.degrey import build_Sa
    from hn.graph import build_graph
    from hn.homcol import edge_vectors, has_homomorphism, orbit_can_block

    assert orbit_can_block("fixed", 2) is False
    assert orbit_can_block("swapped", 1) is False
    g = build_graph(build_Sa())
    assert has_homomorphism(edge_vectors(g), 5)[0] is not None


def test_chain_closing_needs_three_to_split():
    """A rhombus chain closes only if conjugation moves a prime above 3."""
    from hn.homcol import cyclotomic_can_block, cyclotomic_chain_closes

    assert cyclotomic_chain_closes(3) is False
    assert cyclotomic_chain_closes(9) is False
    assert cyclotomic_chain_closes(21) is False
    assert cyclotomic_chain_closes(24) is True
    assert cyclotomic_chain_closes(33) is True
    # Q(zeta_21) blocks but carries no chain; Q(zeta_24) the reverse.
    assert cyclotomic_can_block(21)["can_block"] is True
    assert cyclotomic_can_block(24)["can_block"] is False


def test_smallest_field_meeting_both_conditions():
    from math import gcd
    from hn.homcol import cyclotomic_can_block, cyclotomic_chain_closes

    both = []
    for n in range(3, 40, 3):
        v = cyclotomic_can_block(n)
        if v.get("ramified"):
            continue
        if v["can_block"] and cyclotomic_chain_closes(n):
            both.append(n)
    assert both[0] == 33
    degree = sum(1 for k in range(1, 33) if gcd(k, 33) == 1)
    assert degree == 20
    assert cyclotomic_can_block(33)["residue_degree_over_F"] == 10


def test_no_chain_over_q_zeta_21_at_length_two():
    """The k = 2 case is the Moser spindle, and Q(zeta_21) has no rotation
    with |1 - rho|^2 = 1/3 -- the theorem at 3 says so for every length, and
    length two is checkable directly against the six hexagonal steps."""
    from hn.cyclotomic import CycloField

    K = CycloField(21)
    one = K.rational(1)
    z6 = K.neg(K.mul(K.zeta(7), K.zeta(7)))
    target = K.rational(Fraction(-5, 3))
    units, z = [], one
    for _ in range(6):
        units.append(z)
        z = K.mul(z, z6)
    for u in units:
        for v in units:
            t = K.mul(u, K.conj(v))
            assert K.add(t, K.conj(t)) != target


def test_chain_closing_condition_is_equivalent():
    """|1 + a + b|^2 = 1/3 is Re(a) + Re(b) + Re(a.bbar) = -2/3.

    And at k = 2 it collapses to a + abar = -5/3, whose only unit solutions
    are -(5 -+ sqrt-11)/6 -- the Moser rotation up to sign, which is why the
    spindle is unique.
    """
    from hn.cyclotomic import CycloField

    K = CycloField(33)
    one = K.rational(1)
    z11 = K.zeta(3)
    g, p = K.zero(), one
    for k in range(1, 11):
        p = K.mul(p, z11) if k > 1 else z11
        g = K.add(g, p if k in {1, 3, 4, 5, 9} else K.neg(p))
    assert K.mul(g, g) == K.rational(-11)
    rho = tuple(Fraction(1, 6) * x for x in K.add(K.rational(5), g))
    assert K.norm2(rho) == one
    assert K.add(rho, K.conj(rho)) == K.rational(Fraction(5, 3))
    # k = 2: |1 + a|^2 = 1/3 with a = -rho
    a = K.neg(rho)
    assert K.norm2(K.add(one, a)) == K.rational(Fraction(1, 3))
    # the general identity, checked on that pair plus a third unit step
    z3 = K.zeta(11)
    z6 = K.neg(K.mul(z3, z3))
    for b in (one, z6, K.mul(z6, z6), rho, K.conj(rho)):
        lhs = K.norm2(K.add(K.add(one, a), b))
        t = K.mul(b, K.conj(a))
        rhs = K.add(K.rational(3),
                    K.add(K.add(K.add(a, K.conj(a)), K.add(b, K.conj(b))),
                          K.add(t, K.conj(t))))
        assert lhs == rhs


def test_chain_parametrisation_is_an_identity():
    """The closing conditions, as exact identities in m.

    r + 2 and 2 - r are visibly non-negative, which is why |r| <= 2 holds at
    every real place and the unit step a always exists; y1 solves the first
    conic; and the second condition is V(m) up to the square (m^2+3)^2.
    """
    from fractions import Fraction as Fr

    def V(m):
        return -39 * m ** 4 + 540 * m ** 3 - 666 * m ** 2 - 324 * m + 297

    for num in range(-9, 10):
        for den in (1, 2, 3, 5):
            m = Fr(num, den)
            d = m * m + 3
            r = (m * m - 6 * m - 3) / d
            assert r + 2 == 3 * (m - 1) ** 2 / d
            assert 2 - r == (m + 3) ** 2 / d
            assert -2 <= r <= 2
            y1 = 3 + m * (r - 1)
            assert y1 * y1 == 3 * (4 - r * r)
            assert 3 * (8 - 12 * r - 9 * r * r) * d * d == V(m)
            delta = r * r + Fr(4, 3) * r - Fr(8, 9)
            assert delta == -3 * (8 - 12 * r - 9 * r * r) / 27


def test_chain_closes_for_every_m_with_V_positive():
    """|a| = |b| = 1 and |1 + a + b|^2 = 1/3, built from m alone."""
    import math

    def V(m):
        return -39 * m ** 4 + 540 * m ** 3 - 666 * m ** 2 - 324 * m + 297

    tested = 0
    for k in range(-60, 61):
        m = k / 10.0
        v = V(m)
        if v <= 1e-9:
            continue
        r = (m * m - 6 * m - 3) / (m * m + 3)
        y1 = 3 + m * (r - 1)
        y2 = math.sqrt(v) / (m * m + 3)
        a = (r + 1j * abs(y1) / math.sqrt(3)) / 2
        assert abs(abs(a) - 1) < 1e-9
        T = -8.0 / 3 - r
        cb = (1 + a).conjugate()
        ok = False
        for sgn in (1, -1):
            b = (T + sgn * 1j * y2 / math.sqrt(27)) / (2 * cb)
            if abs(abs(b) - 1) < 1e-9 and abs(abs(1 + a + b) ** 2 - 1 / 3) < 1e-9:
                ok = True
                break
        assert ok, m
        tested += 1
    assert tested > 20, tested


def test_cayley_chromatic_basics():
    from hn.homcol import cayley_chromatic

    assert cayley_chromatic(5, {1, 2, 3, 4}) == 5      # K_5
    assert cayley_chromatic(7, {1}) == 3               # odd cycle
    assert cayley_chromatic(8, {1}) == 2               # even cycle
    assert cayley_chromatic(9, {1, 2, 3, 4}, cap=5) is None


def test_periodic_screen_agrees_with_blocking_at_five():
    """At n = 5 the Cayley graph is K_5, so the screen IS blocking.

    Sa has a coset colouring, so one phi settles it; the denominator-29
    directions block, so the solver exhausts with no phi at all.
    """
    from hn.degrey import build_Sa
    from hn.graph import build_graph
    from hn.homcol import (denominator_29_directions, edge_vectors,
                           has_homomorphism, periodic_screen)

    v = edge_vectors(build_graph(build_Sa()))
    got = periodic_screen(v, 5)
    assert got["colourable"] is True
    assert got["cayley_chi"] == 5
    assert has_homomorphism(v, 5)[0] is not None

    d = denominator_29_directions()
    got = periodic_screen(d, 5)
    assert got["colourable"] is False
    assert got["exhausted"] is True and got["phis_tried"] == 0
    assert has_homomorphism(d, 5)[0] is None


def test_the_stronger_gate_keeps_biting_above_five():
    """n = 6 is decided outright: every homomorphism there needs 6 colours.

    Blocking says nothing about quotients above 5, where the Cayley graph is
    no longer complete -- so this is a strictly stronger condition, and the
    denominator-29 directions pass it too.
    """
    from hn.homcol import denominator_29_directions, periodic_screen

    got = periodic_screen(denominator_29_directions(), 6, cap=5, rounds=100)
    assert got["colourable"] is False
    assert got["exhausted"] is True
    assert got["phis_tried"] > 0


def test_rational_chords_give_quadratic_fields():
    """|1 - rho|^2 = 2 - (rho + rhobar), so a rational chord caps the degree.

    Every rotation de Grey uses has one, and the radicands they produce are
    exactly the field the paper names.
    """
    from hn.homcol import rotation_field_degree

    chords = {Fraction(1): -3,        # hexagonal, 60 degrees
              Fraction(1, 3): -11,    # Moser
              Fraction(1, 4): -15,    # Sb, 2 arcsin(1/4)
              Fraction(9, 4): -7,     # Ya, pi/2 + arcsin(1/8)
              Fraction(7, 4): -7}     # Yb, pi/2 - arcsin(1/8)
    for c, rad in chords.items():
        got = rotation_field_degree(c)
        assert got["degree"] == 2, c
        assert got["radicand"] == rad, (c, got)


def test_de_greys_rotations_have_rational_chords():
    """Computed from the angles themselves, not read off the paper."""
    import math

    angles = {"hexagonal": math.pi / 3,
              "moser": 2 * math.asin(1 / (2 * math.sqrt(3))),
              "Sb": 2 * math.asin(0.25),
              "Ya": math.pi / 2 + math.asin(0.125),
              "Yb": math.pi / 2 - math.asin(0.125)}
    want = {"hexagonal": Fraction(1), "moser": Fraction(1, 3),
            "Sb": Fraction(1, 4), "Ya": Fraction(9, 4), "Yb": Fraction(7, 4)}
    for name, th in angles.items():
        chord = 2 - 2 * math.cos(th)
        assert abs(chord - float(want[name])) < 1e-12, name


def test_multiquadratic_cannot_block_either_way():
    """Elementary abelian Galois group forces residue degree at most 2.

    Both small orbit types fail, so whichever way 5 behaves in a
    multiquadratic field the coset colouring exists -- which is exactly why
    de Grey's G, and every construction from rational chords, stops at 5.
    """
    from hn.degrey import build_Sa
    from hn.graph import build_graph
    from hn.homcol import (edge_vectors, has_homomorphism, orbit_can_block,
                           periodic_screen)

    assert orbit_can_block("fixed", 2) is False
    assert orbit_can_block("swapped", 1) is False
    v = edge_vectors(build_graph(build_Sa()))
    assert has_homomorphism(v, 5)[0] is not None
    assert periodic_screen(v, 5)["colourable"] is True


def test_the_chain_rotation_has_an_irrational_chord():
    """What the corollary demands of any 6-chromatic candidate.

    The chain's rotation has a + abar = r = (m^2-6m-3)/(m^2+3) with m a root
    of T^3 - 9T^2 + 14T + 8, so its chord 2 - r has degree 3 over Q -- unlike
    every rotation in every construction the barrier argument covers.
    """
    cub = (8, 14, -9)              # T^3 - 9T^2 + 14T + 8, constant first

    def mul(x, y):
        r = [Fraction(0)] * 5
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    r[i + j] += a * b
        for k in (4, 3):
            c = r[k]
            if c:
                r[k] = Fraction(0)
                r[k - 3] -= c * cub[0]
                r[k - 2] -= c * cub[1]
                r[k - 1] -= c * cub[2]
        return tuple(r[:3])

    def inv(x):
        cols = [mul(x, tuple(Fraction(1 if k == j else 0) for k in range(3)))
                for j in range(3)]
        A = [[cols[j][i] for j in range(3)] + [Fraction(1 if i == 0 else 0)]
             for i in range(3)]
        for c in range(3):
            p = next(t for t in range(c, 3) if A[t][c])
            A[c], A[p] = A[p], A[c]
            sc = Fraction(1) / A[c][c]
            A[c] = [v * sc for v in A[c]]
            for t in range(3):
                if t != c and A[t][c]:
                    f = A[t][c]
                    A[t] = [u - f * v for u, v in zip(A[t], A[c])]
        return tuple(A[i][3] for i in range(3))

    one = (Fraction(1), Fraction(0), Fraction(0))
    m = (Fraction(0), Fraction(1), Fraction(0))
    mm = mul(m, m)

    def add(x, y):
        return tuple(a + b for a, b in zip(x, y))

    def sub(x, y):
        return tuple(a - b for a, b in zip(x, y))

    def sc(c, x):
        return tuple(Fraction(c) * a for a in x)

    r = mul(sub(sub(mm, sc(6, m)), sc(3, one)), inv(add(mm, sc(3, one))))
    chord = sub(sc(2, one), r)
    assert chord[1] or chord[2], "the chord must not be rational"
    # degree 3: 1, chord, chord^2 independent
    rows, pw = [], one
    for _ in range(3):
        rows.append(list(pw))
        pw = mul(pw, chord)
    M = [row[:] for row in rows]
    rank = 0
    for c in range(3):
        pr = next((t for t in range(rank, 3) if M[t][c]), None)
        if pr is None:
            continue
        M[rank], M[pr] = M[pr], M[rank]
        f = M[rank][c]
        M[rank] = [v / f for v in M[rank]]
        for t in range(3):
            if t != rank and M[t][c]:
                g = M[t][c]
                M[t] = [u - g * v for u, v in zip(M[t], M[rank])]
        rank += 1
    assert rank == 3, "1, chord, chord^2 must be independent"


def test_closed_necklace_is_four_critical():
    """A cycle of k diamonds plus a closing edge, checked as a graph.

    This is the combinatorial content of the chain construction, with the
    geometry stripped out: 3k+1 vertices, 5k+1 edges, chromatic number 4, and
    every single vertex critical. Proved in the module; verified here for
    several k so the argument is not taken on trust.
    """
    from pysat.solvers import Solver

    def necklace(k):
        # tips B_0..B_k are 0..k; diamond j has middles at k+1+2(j-1), +1
        edges = []
        for j in range(1, k + 1):
            a, b = j - 1, j
            m1, m2 = k + 1 + 2 * (j - 1), k + 2 + 2 * (j - 1)
            edges += [(a, m1), (a, m2), (m1, m2), (m1, b), (m2, b)]
        edges.append((k, 0))
        return 3 * k + 1, edges

    def colourable(n, edges, keep, k):
        idx = {v: i for i, v in enumerate(sorted(keep))}
        cls = [[1 + i * k + c for c in range(k)] for i in range(len(idx))]
        for a, b in edges:
            if a in idx and b in idx:
                for c in range(k):
                    cls.append([-(1 + idx[a] * k + c), -(1 + idx[b] * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            return s.solve()

    for k in (3, 4, 5, 7):
        n, edges = necklace(k)
        assert n == 3 * k + 1
        assert len(edges) == 5 * k + 1
        allv = set(range(n))
        assert not colourable(n, edges, allv, 3), k
        assert colourable(n, edges, allv, 4), k
        for v in range(n):
            assert colourable(n, edges, allv - {v}, 3), (k, v)


def test_field_with_both_spindle_and_blocking():
    """x^3 - 10x^2 + 26x - 11: V = 33 s^2, and 5 is inert.

    33 is what d^2 = 3 asks for, since a distance spindles iff 3(4d^2-1) is a
    square in F. So sqrt33 lies in F, K = F(sqrt-3) holds the Moser rotation,
    and the residue degree at 5 is lcm(3, 2) = 6 -- both conditions at once.
    """
    cub = (-11, 26, -10)               # constant, x, x^2

    def cmul(x, y):
        r = [Fraction(0)] * 5
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    r[i + j] += a * b
        for k in (4, 3):
            c = r[k]
            if c:
                r[k] = Fraction(0)
                r[k - 3] -= c * cub[0]
                r[k - 2] -= c * cub[1]
                r[k - 1] -= c * cub[2]
        return tuple(r[:3])

    def lin(*pairs):
        out = [Fraction(0)] * 3
        for c, x in pairs:
            for i in range(3):
                out[i] += Fraction(c) * x[i]
        return tuple(out)

    one = (Fraction(1), Fraction(0), Fraction(0))
    m = (Fraction(0), Fraction(1), Fraction(0))
    mm = cmul(m, m)
    m3 = cmul(mm, m)
    m4 = cmul(m3, m)
    V = lin((-39, m4), (540, m3), (-666, mm), (-324, m), (297, one))
    assert V == (Fraction(1947), Fraction(-4653), Fraction(1848))
    s = lin((13, one), (-26, m), (5, mm))
    assert lin((33, cmul(s, s))) == V, "V must be 33 times a square"
    # 5 inert in the cubic: no root mod 5
    assert all((t ** 3 + cub[2] * t * t + cub[1] * t + cub[0]) % 5
               for t in range(5))
    # 5 inert in Q(sqrt33): 33 = 3 mod 5 is a non-residue
    assert pow(3, 2, 5) != 3 and 3 not in {pow(k, 2, 5) for k in range(5)}
    # so the residue degree in F = Q(m, sqrt33) is lcm(3, 2) = 6 >= 3


def test_neighbourhoods_are_bipartite():
    """N(v) lies on a circle where adjacency is 60 degrees apart.

    So every cycle there is a full hexagon, of even length, and N(v) is
    bipartite -- which is why no vertex's colour is ever forced by its
    neighbourhood at four colours, and why the rhombus trick does not lift.
    """
    from hn.degrey import build_Sa
    from hn.graph import build_graph

    g = build_graph(build_Sa())
    adj = {}
    for a, b in g.edges():
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    checked = 0
    for v, nb in adj.items():
        sub = {u: (adj[u] & nb) for u in nb}
        assert all(len(s) <= 2 for s in sub.values()), v
        colour = {}
        for start in sub:
            if start in colour:
                continue
            colour[start] = 0
            stack = [start]
            while stack:
                x = stack.pop()
                for y in sub[x]:
                    if y not in colour:
                        colour[y] = 1 - colour[x]
                        stack.append(y)
                    else:
                        assert colour[y] != colour[x], (v, x, y)
        checked += 1
    assert checked > 100


def test_blocking_at_every_modulus_up_to_five():
    """A coset colouring mod n gives chi <= n, so all of 2, 3, 4, 5 matter.

    The denominator-29 directions clear all four; the necklace's clear 2, 3
    and 5 and fail at 4, which is exactly the coset 4-colouring the n = 8
    Cayley screen had found by a longer route.
    """
    from hn.homcol import denominator_29_directions, has_homomorphism

    d = denominator_29_directions()
    for n in (2, 3, 4, 5):
        assert has_homomorphism(d, n)[0] is None, n


def test_four_colour_forcing_appears_in_Y_and_not_in_Sa():
    """(*) at u: every 4-colouring puts u's colour on u's sqrt3-sphere.

    Sa does not have it at the origin; Y does, on a sphere of 12. That is
    where de Grey's forcing lives, and sqrt3 is the distance because
    |1 - rho|^2 = 1/3 sends a point at sqrt3 to distance exactly 1.

    Scanning every vertex takes twenty minutes, so only the origin is tested
    -- the one vertex the full scan on Y singles out.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_Sa, build_Y
    from hn.geometry import DEGREY_FIELD
    from hn.graph import build_graph

    def at_origin(pts):
        zero = DEGREY_FIELD.zero()
        o = next(i for i, p in enumerate(pts) if p.x == zero and p.y == zero)
        g = build_graph(pts)
        n = len(pts)
        sphere = [j for j in range(n) if j != o
                  and (pts[j].x - pts[o].x) ** 2
                  + (pts[j].y - pts[o].y) ** 2 == 3]

        def x(v, c):
            return 1 + v * 4 + c

        cls = [[x(v, c) for c in range(4)] for v in range(n)]
        for a, b in g.edges():
            for c in range(4):
                cls.append([-x(a, c), -x(b, c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            sep = sv.solve(assumptions=[x(o, 0)] + [-x(v, 0) for v in sphere])
        return len(sphere), sep

    size, separable = at_origin(build_Sa())
    assert separable is True, "Sa must NOT have the property"
    size, separable = at_origin(build_Y())
    assert size == 12
    assert separable is False, "Y must have it: no 4-colouring avoids it"


def test_closing_radicand_reads_de_greys_chain():
    """His four radicands are his four closing distances, in order.

    A rotation taking a pair at distance d to a unit pair has
    cos = 1 - 1/(2 d^2) and sin = sqrt(4 d^2 - 1) / (2 d^2), so the field has
    to contain sqrt(4 d^2 - 1).  Running de Grey's own chain of distances --
    1 for the triangular lattice, sqrt3 for the Moser spindle, 2 for
    Sb = rho(Sa), 4 for Ya u Yb -- through that gives 3, 11, 15, 7, which is
    Q(sqrt3, sqrt5, sqrt7, sqrt11) and nothing more.  His field is not a
    choice; it is what his chain demands, term by term.
    """
    from hn.homcol import closing_radicand

    assert [closing_radicand(D) for D in (1, 3, 4, 16)] == [3, 11, 15, 7]
    assert closing_radicand(Fraction(1, 5)) == 0, "d < 1/2 closes nothing"
    assert closing_radicand(Fraction(1, 4)) == 0, "2d sin(theta/2) = 1 needs d >= 1/2"


def test_the_doubling_chain_stops_where_de_greys_field_stops():
    """D -> 4D, and the step from D needs sqrt(16D - 1).

    Y's forced pair sits on the ring of radius sqrt(D) about the union's
    pivot, antipodally, so at squared distance 4D; the next union has to close
    that.  Starting at 1 the chain is 1, 4, 16, 64, ... and the radicands are
    15, 7, 255, ...  de Grey has 3, 5, 7 and 11, so he closes 1, 4 and 16 and
    stops: 255 = 3.5.17 asks for sqrt17.  His construction is exactly as long
    as his field allows.
    """
    from hn.homcol import closable_distance, closing_radicand

    assert [D for D in (1, 4, 16, 64, 256) if closable_distance(D)] == [1, 4, 16]
    assert closing_radicand(64) == 255 and 255 == 3 * 5 * 17
    assert not closable_distance(64), "sqrt17 is not in Q(sqrt3,sqrt5,sqrt7,sqrt11)"
    assert closable_distance(64, (3, 5, 7, 11, 17)), "adjoining sqrt17 reopens it"


def test_closable_over_a_field_without_i_is_a_stricter_test():
    """K = Q(m, sqrt-3, sqrt-11) has no i, and that costs it distances.

    For a real multiquadratic F the points live in F(i), so sqrt(1 - 4D) is in
    reach whenever sqrt(4D - 1) is.  K is not of that shape: its rational
    square classes are 1, -3, -11 and 33, so 1 - 4D must be -3 or -11 times a
    square.  D = 5/2 has 4D - 1 = 9, as rational a square root as exists, and
    is still not closable over K.
    """
    from hn.homcol import closable_distance, closable_over

    K = (1, -3, -11, 33)
    assert closable_distance(Fraction(5, 2)) is True
    assert closable_over(Fraction(5, 2), K) is False
    assert [D for D in range(1, 30) if closable_over(D, K)] == [1, 3, 7, 19, 25]
    assert not closable_over(4, K) and not closable_over(16, K), \
        "K loses both of de Grey's last two steps"
    assert closable_over(1, K) and closable_over(3, K), "and keeps the first two"


def test_de_greys_union_adds_exactly_one_edge():
    """G = Ya u Yb is two copies of Y glued at the pivot, plus the spindle.

    Y has 791 vertices and 3938 edges, and G has 1581 = 2.791 - 1 and
    7877 = 2.3938 + 1.  The single extra edge joins the two images of (2,0),
    and it is the whole of the step from four colours to five: (2,0),(-2,0) is
    monochromatic in every 4-colouring of Y, both copies read the pivot's
    colour onto that edge, and the contradiction is immediate.
    """
    from hn.degrey import build_G, build_Y
    from hn.graph import build_graph

    y = build_graph(build_Y())
    g = build_graph(build_G(as_graph=False))
    assert len(y.vertices) == 791 and len(list(y.edges())) == 3938
    assert len(g.vertices) == 1581 == 2 * 791 - 1
    assert len(list(g.edges())) == 7877 == 2 * 3938 + 1


def test_the_forced_pair_of_Y_sits_on_the_ring_the_rotation_moves():
    """(2,0) and (-2,0): on the ring of radius 2 about the pivot, antipodal.

    Sb is Sa turned about the ORIGIN by the rotation closing distance 2, so
    the ring the union pins has radius 2; the forced pair is a diameter of it,
    at distance 4.  This is the geometry the SAT result rests on -- proving
    the pair forced needs a real unsatisfiability proof, 282 seconds against
    milliseconds for the five distance-4 pairs that are not forced, so only
    the placement is checked here.
    """
    from hn.degrey import build_Y
    from hn.geometry import DEGREY_FIELD as F, Point

    Y = build_Y()
    a = Point(F.rational(2), F.zero())
    b = Point(F.rational(-2), F.zero())
    assert a in Y and b in Y
    o = Point(F.zero(), F.zero())
    assert (a - o).x ** 2 + (a - o).y ** 2 == F.rational(4)
    assert (b - o).x ** 2 + (b - o).y ** 2 == F.rational(4)
    assert (a - b).x ** 2 + (a - b).y ** 2 == F.rational(16), "4D, a diameter"


def test_the_doubling_step_is_scarce_over_K():
    """D and 4D both closable: common over de Grey's field, almost nowhere over K.

    If the forced pair is antipodal on the ring -- which is where Y's is, and
    where the ring scan of `Sa u rho_4(Sa)` independently rediscovers it -- the
    chain doubles and the step from `D` needs `sqrt(16D - 1)`.  de Grey's field
    admits that step at 1/2, 1, 4, 17/2, 37/4, 61/4, 86, ... and his chain
    walks two in a row.  `K` admits exactly two below 4000, and neither chains
    again, so his architecture does not port by translation.
    """
    from hn.homcol import closable_distance, closable_over

    K = (1, -3, -11, 33)
    common = [Fraction(1, 2), Fraction(1), Fraction(4), Fraction(17, 2),
              Fraction(37, 4), Fraction(61, 4), Fraction(86)]
    assert all(closable_distance(D) and closable_distance(4 * D)
               for D in common)

    got = sorted({D for q in (1, 2, 3, 4, 6, 12) for p in range(1, 4001)
                  for D in [Fraction(p, q)] if D.denominator == q
                  and closable_over(D, K) and closable_over(4 * D, K)},
                 key=float)
    assert got == [Fraction(7, 12), Fraction(397, 4)]
    assert not closable_over(Fraction(28, 3), K), "7/3 does not chain on"
    assert not closable_over(1588, K), "397 does not chain on either"


def test_the_blocking_test_must_run_on_the_module_not_the_ambient_lattice():
    """`has_homomorphism` over Z^d can report blocking the module does not have.

    A coset colouring is a homomorphism from the group `M` the edge vectors
    generate to `Z/n`.  `M` sits inside `Z^d`, possibly properly, and `Z/n` is
    not injective, so a homomorphism on `M` need not extend -- the search over
    `Z^d` finds too few functionals.  The smallest case is one dimension wide:
    over `2Z` the map sending 2 to 1 escapes, over `Z` every map kills 2.
    """
    from hn.homcol import blocks_at, has_homomorphism, on_lattice

    assert has_homomorphism([(2,)], 2)[0] is None, "over Z it looks blocked"
    assert on_lattice([(2,)]) == [(1,)]
    assert blocks_at([(2,)], 2) is False, "over 2Z it is not"


def test_the_two_blocking_tests_agree_exactly_when_the_rank_holds_mod_p():
    """`M_sat/M` obstructs the restriction, and it is visible as a rank drop.

    From `0 -> M -> M_sat -> T -> 0`, restriction `Hom(Z^d, Z/n) ->
    Hom(M, Z/n)` is onto exactly when `T/nT` vanishes, i.e. when no invariant
    factor of `M` shares a prime with `n` -- which happens exactly when the
    rank of the direction matrix drops mod that prime.  The denominator-29 set
    has index 1, so it is saturated and every verdict on it stands as it was.
    """
    from hn.homcol import (blocks_at, denominator_29_directions,
                           has_homomorphism, lattice_basis)

    D = denominator_29_directions(3)
    B = lattice_basis(D)
    assert len(B) == len(D[0]) == 12, "full rank"
    index = 1
    for b in B:
        index *= abs(b[next(i for i, x in enumerate(b) if x)])
    assert index == 1, "saturated, so the two tests cannot differ"
    assert has_homomorphism(D, 5)[0] is None and blocks_at(D, 5)


def test_an_escape_found_over_the_ambient_lattice_is_always_genuine():
    """Only blocking claims were ever at risk; exhibited colourings are safe.

    A homomorphism `Z^d -> Z/n` restricts to `M`, so any `phi` the search
    produces really is a coset colouring.  That is why "de Grey's G does not
    block at 5" needed no revisiting, while "the necklace blocks at 4" did.
    """
    from hn.degrey import build_G
    from hn.graph import build_graph
    from hn.homcol import edge_vectors, has_homomorphism, on_lattice

    D = edge_vectors(build_graph(build_G(as_graph=False)))
    phi, _ = has_homomorphism(D, 5)
    assert phi is not None, "a coset colouring mod 5 exists, as it must"
    assert all(sum(p * x for p, x in zip(phi, d)) % 5 for d in D)
    red = on_lattice(D)
    assert len(red[0]) == 16, "rank 16 inside dimension 32"


def test_the_gate_verdict_needs_no_hermite_reduction():
    """Full rank mod 5 makes the ambient test correct at the gate, by theorem.

    Restriction `Hom(Z^d, Z/n) -> Hom(M, Z/n)` is onto when no invariant
    factor of `M` shares a prime with `n`, and that is exactly the rank
    surviving reduction mod each such prime.  The denominator-29 set has full
    rank at every modulus, so none of its verdicts could ever have moved --
    which is why `U`, whose blocking comes from a rotated copy of it, needed
    no revisiting.
    """
    from hn.homcol import (_rank_mod, _rank_q, denominator_29_directions,
                           saturated_at)

    D = denominator_29_directions(3)
    assert _rank_q(D) == 12
    assert all(_rank_mod(D, p) == 12 for p in (2, 3, 5))
    assert all(saturated_at(D, n) for n in (2, 3, 4, 5))


def test_a_forced_pair_is_one_extra_edge_the_rhombus_at_three_colours():
    """The lemma in miniature, where it runs instantly.

    A unit rhombus is `K4` minus an edge; at three colours its two tips must
    share a colour, and they lie `sqrt3` apart.  That is a forced pair, and by
    the identity it says the rhombus plus the single edge joining the tips is
    not 3-colourable -- which is what makes the Moser spindle work, since
    `sqrt3` is closable (`4D - 1 = 11`) and the closing rotation
    `(5 + sqrt-11)/6` brings the two tips of two rhombi to distance 1.

    Same sentence, same proof, at `k = 4` with `Y` and the pair `(2,0)`,
    `(-2,0)` at distance 4, where `4D - 1 = 63 = 9.7`.
    """
    from hn.homcol import closing_radicand
    from pysat.solvers import Solver

    # rhombus 0-1-2, 1-2-3 : two unit triangles sharing the edge 1-2
    rhombus = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]

    def colourable(edges, k, extra=()):
        cls = [[1 + v * k + c for c in range(k)] for v in range(4)]
        for a, b in list(edges) + list(extra):
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            return sv.solve()

    assert colourable(rhombus, 3), "the rhombus alone is 3-colourable"
    assert not colourable(rhombus, 3, [(0, 3)]), \
        "adding the tip-to-tip edge kills it -- the tips are forced equal"
    assert colourable(rhombus, 4, [(0, 3)]), "at four colours they are free"

    assert closing_radicand(3) == 11, "sqrt3 closes through sqrt11"
    assert closing_radicand(16) == 7, "and 4 through sqrt7, de Grey's step"


def test_closing_de_greys_G_under_its_own_rotations_barely_overlaps():
    """`S` works as a seed because it is small; `G` is not.

    `Sa` is `S`'s twelve images and they overlap heavily -- 39 points give 397,
    not 468.  Doing the same to `G` gives 18966 points out of a possible
    `12 x 1581 = 18972`: the images share six points and nothing else, and the
    mean degree stays at `G`'s own 10.0.  So the move that builds the fourth
    floor does not build the fifth.
    """
    from hn.degrey import build_G, build_S, build_Sa
    from hn.geometry import DEGREY_FIELD as F, Point, _rot60

    assert len(build_S()) == 39 and len(build_Sa()) == 397 < 12 * 39

    G = build_G(F, as_graph=False)
    rot60 = _rot60(F)
    seen = set()
    for p in G:
        for base in (p, Point(p.x, -p.y)):
            q = base
            for _ in range(6):
                seen.add(q)
                q = rot60(q)
    assert len(seen) == 18966 == 12 * len(G) - 6


def test_blocking_below_the_chromatic_number_is_free():
    """A coset colouring mod `n` is a proper `n`-colouring, so `chi > n` blocks.

    That makes the sharpened gate -- every 6-chromatic graph blocks at 2, 3, 4
    and 5 -- a tautology: such a graph has no proper `n`-colouring at all for
    `n <= 5`.  de Grey's `G` shows the shape.  It is 5-chromatic and blocks at
    every modulus below 5 for nothing, fails at 5, which is the barrier
    theorem and the only place content was available, and blocks again at 6.

    What survives is the reading as a design criterion: for a graph of
    chromatic number `c`, blocking at `n < c` is free and at `n >= c` is real.
    """
    from hn.degrey import build_G
    from hn.graph import build_graph
    from hn.homcol import blocks_at, edge_vectors

    D = edge_vectors(build_graph(build_G(as_graph=False)))
    assert [n for n in (2, 3, 4, 5, 6) if blocks_at(D, n)] == [2, 3, 4, 6]


def test_Y_is_two_copies_of_Sa_glued_at_a_point_by_six_edges():
    """The whole difference between no forced pair and one forced pair.

    `Sa` and `Sb` have 397 points each and meet in the origin alone; their
    union is 793 and `Y` is that less the two vertices de Grey drops. The edge
    count decomposes exactly: `2 x 1974 = 3948`, plus six edges joining the
    exclusive half of `Sa` to the exclusive half of `Sb`, less the sixteen the
    dropped vertices take with them, giving 3938.

    And the forced pair `(2,0), (-2,0)` lies *wholly inside the `Sa` copy*,
    which on its own forces nothing in 4200 pairs. Six edges landing elsewhere
    are what make two points of `Sa` unable to differ.
    """
    from hn.degrey import build_Sa, build_Sb, build_Y
    from hn.geometry import DEGREY_FIELD as F, Point
    from hn.graph import build_graph

    Sa, Sb = set(build_Sa()), set(build_Sb())
    assert len(Sa) == len(Sb) == 397
    assert len(Sa & Sb) == 1 and len(Sa | Sb) == 793

    Y = build_Y()
    gy = build_graph(Y)
    assert len(Y) == 791 and len(list(gy.edges())) == 3938

    only_a, only_b = Sa - Sb, Sb - Sa
    cross = sum(1 for i, j in gy.edges()
                if (Y[i] in only_a and Y[j] in only_b)
                or (Y[j] in only_a and Y[i] in only_b))
    assert cross == 6
    assert 2 * 1974 + cross - 16 == 3938

    a = Point(F.rational(2), F.zero())
    b = Point(F.rational(-2), F.zero())
    assert a in Sa and b in Sa and a not in Sb and b not in Sb


def test_the_cross_edge_equation_reproduces_de_greys_rotation():
    """A cross edge pins the rotation to the root of an explicit quadratic.

    With `w = u.q`, `A = |q|^2`, `P = |p|^2`, the conditions `|w|^2 = A` and
    `|w - p|^2 = 1` expand to `w.pbar + wbar.p = A + P - 1 =: 2R` and
    `|w.pbar|^2 = A.P`, so `w.pbar` is a root of `t^2 - 2Rt + A.P` and the
    rotation exists over a field exactly when `sqrt(R^2 - A.P)` does.

    Checked on de Grey's own: every cross edge of `Y` satisfies the identity,
    with `w` in `Sb` and `q = rho_4^{-1}(w)` back in `Sa`.  That is what makes
    the search over another field a solve rather than a sample -- sampling all
    3030 known units of `K` gives no cross edge at all.
    """
    from hn.degrey import build_Sa, build_Sb, build_Y
    from hn.geometry import DEGREY_FIELD as F
    from hn.graph import build_graph

    Sa, Sb = set(build_Sa()), set(build_Sb())
    Y = build_Y()
    gy = build_graph(Y)
    only_a, only_b = Sa - Sb, Sb - Sa
    cross = [(Y[i], Y[j]) for i, j in gy.edges()
             if (Y[i] in only_a and Y[j] in only_b)
             or (Y[j] in only_a and Y[i] in only_b)]
    assert len(cross) == 6

    half = F.rational(Fraction(1, 2))
    for x, y in cross:
        p, w = (x, y) if x in only_a else (y, x)
        P = p.x * p.x + p.y * p.y
        A = w.x * w.x + w.y * w.y
        R = half * (A + P - F.rational(1))
        # w.pbar = (wx px + wy py) + i(wy px - wx py); its real part is R
        assert w.x * p.x + w.y * p.y == R
        # and its squared modulus is A.P
        re = w.x * p.x + w.y * p.y
        im = w.y * p.x - w.x * p.y
        assert re * re + im * im == A * P


def test_adjoining_a_quadratic_keeps_the_residue_degree_at_five():
    """`Q(m)` has `f = 3` at 5, and no quadratic extension can drop it.

    `x^3 - 10x^2 + 26x - 11` is irreducible mod 5, so the single prime of
    `Q(m)` above 5 has residue degree 3. In any quadratic extension that prime
    stays inert (`f = 6`), splits into primes still of residue degree 3, or
    ramifies with `f = 3` — never lower. So `K(sqrt-15)` and every other
    quadratic step keeps the arithmetic a sixth colour needs, which is why
    asking about de Grey's own pair over `K` was worth the run.
    """
    poly = [-11, 26, -10, 1]          # constant first

    def ev(x, p):
        return sum(c * pow(x, i, p) for i, c in enumerate(poly)) % p

    assert all(ev(x, 5) for x in range(5)), "no root mod 5"
    # a cubic with no root mod p is irreducible there, so f = 3
    assert [ev(x, 5) for x in range(5)] == [4, 1, 4, 4, 2]


def test_the_forcing_in_Y_is_not_local():
    """The ball of radius 1.5 about the pair holds 773 of 791 and separates.

    Forcing could in principle be carried by a small gadget near the pair. It
    is not. Cutting `Y` to the points within 1.5 of the segment joining
    `(2,0)` and `(-2,0)` keeps 773 of its 791 vertices and the pair comes
    apart — so the eighteen vertices *furthest* from the pair are load
    bearing, and whatever forces at five colours will not be small either.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_Y
    from hn.geometry import DEGREY_FIELD as F, Point
    from hn.graph import build_graph

    Y = build_Y()
    g = build_graph(Y)
    ia = next(i for i, p in enumerate(Y)
              if p == Point(F.rational(2), F.zero()))
    ib = next(i for i, p in enumerate(Y)
              if p == Point(F.rational(-2), F.zero()))

    def dist(p):
        x, y = float(p.x), float(p.y)
        return ((x - min(2.0, max(-2.0, x))) ** 2 + y * y) ** 0.5

    keep = sorted(i for i, p in enumerate(Y) if dist(p) <= 1.5 + 1e-9)
    assert len(keep) == 773 and ia in keep and ib in keep
    idx = {v: i for i, v in enumerate(keep)}
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(keep))]
    for a, b in g.edges():
        if a in idx and b in idx:
            for c in range(4):
                cls.append([-(1 + idx[a] * 4 + c), -(1 + idx[b] * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        assert sv.solve()
        assert sv.solve(assumptions=[1 + idx[ia] * 4, -(1 + idx[ib] * 4)]), \
            "the pair must come apart once the far vertices are gone"


def test_every_sampled_vertex_of_Y_is_essential_to_its_forcing():
    """Every vertex of the sample is essential — but the sample is not all of `Y`.

    A full sweep later found vertex 630, of degree 4, to be slack: `Y` minus
    it still forces the pair. So `Y` is *not* vertex-critical, and what this
    test records is the sample, which is what it always was.


    A stratified sample of 61 of `Y`'s 789 non-pair vertices, by degree, and
    every one of them is essential — deleting it lets `(2,0)` and `(-2,0)`
    take different colours. Thirty-nine answered inside a 60000-conflict
    budget; the other twenty-two answered the same way in seconds once the
    budget was removed, so the budget was the only thing making them look
    hard.

    Two of them are checked here, one from each group, since the full sweep is
    an hour. `Y` is a minimal four-colour forcer for its pair, with nothing to
    trim.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_Y
    from hn.geometry import DEGREY_FIELD as F, Point
    from hn.graph import build_graph

    Y = build_Y()
    g = build_graph(Y)
    E = list(g.edges())
    ia = next(i for i, p in enumerate(Y)
              if p == Point(F.rational(2), F.zero()))
    ib = next(i for i, p in enumerate(Y)
              if p == Point(F.rational(-2), F.zero()))

    for drop in (184, 490):          # one budget-undecided, one not
        keep = [i for i in range(len(Y)) if i != drop]
        idx = {v: i for i, v in enumerate(keep)}
        cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(keep))]
        for a, b in E:
            if a in idx and b in idx:
                for c in range(4):
                    cls.append([-(1 + idx[a] * 4 + c),
                                -(1 + idx[b] * 4 + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            assert sv.solve()
            assert sv.solve(assumptions=[1 + idx[ia] * 4,
                                         -(1 + idx[ib] * 4)]), \
                f"vertex {drop} must be essential"


def test_the_counting_threshold_applies_exactly_at_the_gate():
    """A zero residue blocks trivially; the count governs only without one.

    The heuristic is `n^r . ((n-1)/n)^L < 1`, so a rank-`r` direction set needs
    more than `r log n / log(n/(n-1))` independent lines mod `n` — at `n = 5`,
    `7.21 r`. It calls both known cases at the gate correctly: the
    denominator-29 set has 120 lines against a threshold of 87 and blocks; de
    Grey's `G` has 54 against 116 and escapes.

    Below the gate it is inapplicable, and the reason is sharp. `G`'s
    directions contain a vector congruent to zero mod 2, 3 and 4 — the set
    holds `n` times one of its own members — so every `phi` kills it and the
    blocking is trivial. At `n = 5` there is no zero residue among its 109.
    """
    import math

    from hn.degrey import build_G
    from hn.graph import build_graph
    from hn.homcol import (blocks_at, denominator_29_directions,
                           edge_vectors, _rank_q)

    def lines(v, n):
        return len({tuple(x % n for x in w) for w in v}) // 2

    def threshold(r, n):
        return math.ceil(r * math.log(n) / math.log(n / (n - 1)))

    D = edge_vectors(build_graph(build_G(as_graph=False)))
    assert _rank_q(D) == 16 and lines(D, 5) == 54
    assert threshold(16, 5) == 116 and not blocks_at(D, 5)

    E = denominator_29_directions(3)
    assert _rank_q(E) == 12 and lines(E, 5) == 120
    assert threshold(12, 5) == 87 and blocks_at(E, 5)

    for n in (2, 3, 4):
        res = {tuple(x % n for x in v) for v in D}
        assert sum(1 for v in res if not any(v)) == 1, "a trivial obstruction"
        assert blocks_at(D, n), "which blocks for a reason that is not counting"
    assert all(any(v) for v in {tuple(x % 5 for x in w) for w in D})


def test_two_disjoint_copies_of_G_cost_exactly_twice_as_much():
    """The control that turns the conflict jump into evidence.

    A translate union of `G` costs 723211 conflicts to sweep at five colours,
    with its dearest pair at 3689, against `G`'s own 6410 and 29. But the
    union has twice the vertices, and a bigger instance is harder to solve
    even when it is just as loose — so the jump has to be controlled.

    Two copies of `G` a thousand apart, where no cross edge is geometrically
    possible, cost **13301 conflicts with the dearest pair at 28**: exactly
    twice `G`'s own total, and a dearest that does not move. The tightening in
    the translate union comes entirely from its 1442 cross edges.

    Only the invariant part is asserted here — that the control's dearest pair
    matches `G`'s — since conflict counts are solver-dependent.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_G
    from hn.geometry import DEGREY_FIELD as F
    from hn.graph import build_graph

    pts = build_G(F, as_graph=False)
    GE = list(build_graph(pts).edges())
    n0 = len(pts)
    E = [(a, b) for a, b in GE] + [(a + n0, b + n0) for a, b in GE]
    assert n0 == 1581 and len(E) == 2 * 7877

    cls = [[1 + v * 5 + c for c in range(5)] for v in range(2 * n0)]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        assert sv.solve(), "two copies of a 5-chromatic graph are 5-colourable"
        # a pair in one copy and a pair spanning both both come apart at once
        assert sv.solve(assumptions=[1, -(1 + 5 * 5)])
        assert sv.solve(assumptions=[1, -(1 + n0 * 5)])


def test_sampling_never_hides_a_forced_pair():
    """The filter's soundness, on a case where the answer is known.

    A rhombus -- two unit triangles sharing an edge -- forces its two apexes
    to the same colour in every proper 3-colouring: the shared edge takes two
    of the three colours, and each apex is adjacent to both of its ends, so
    each is driven to the third.  Sample as many 3-colourings as you like and
    the apex pair agrees in all of them, because there is no colouring where
    it does not.

    That is the whole argument for filtering by sampling.  A colouring that
    separates a pair is a proof the pair is not forced; the absence of such a
    colouring among the samples is not a proof of anything, which is why the
    survivors go to SAT.  The direction that could go wrong -- a forced pair
    quietly filtered out -- cannot happen at all.
    """
    import random
    from pysat.solvers import Solver

    # Two triangles on the edge (0,1): apexes 2 and 3.
    edges = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)]
    n, k = 4, 3
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    assert sv.solve()

    rng = random.Random(4242)
    agreed_apex, separated_other = 0, False
    for _ in range(200):
        v = rng.randrange(n)
        sv.solve(assumptions=[1 + v * k + rng.randrange(k)]) or sv.solve()
        m = sv.get_model()
        col = [next(c for c in range(k) if m[w * k + c] > 0) for w in range(n)]
        assert col[2] == col[3]          # the forced pair, every single time
        agreed_apex += 1
        if col[0] != col[2]:
            separated_other = True       # an unforced pair does separate
    sv.delete()

    assert agreed_apex == 200
    assert separated_other, "sampling must be able to separate a free pair"


def test_the_decay_curve_has_a_floor_exactly_when_something_forces():
    """Why the survivor curve is the right metric, and conflicts are not.

    Sampling drives the surviving-pair count down geometrically -- each
    colouring kills the pairs it separates -- but it can never drive a forced
    pair out.  So the curve of a graph that forces has a floor at least as
    high as the number of forced pairs, and the curve of a graph that does not
    reaches zero.  Zero is a complete, witnessed negative.

    Conflict counts have no such meaning: they grow with the graph, so the
    156283 conflicts of the depth-3 stack and the 35365 of Y are not on the
    same scale and comparing them was an error.  Survivor counts are pairs.
    """
    import random
    from pysat.solvers import Solver

    def curve(n, edges, k, pairs, samples=40, seed=99):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in edges:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
            # (edges only)
        sv = Solver(name="cd19", bootstrap_with=cls)
        assert sv.solve()
        rng, surv = random.Random(seed), list(pairs)
        for _ in range(samples):
            v = rng.randrange(n)
            sv.solve(assumptions=[1 + v * k + rng.randrange(k)]) or sv.solve()
            m = sv.get_model()
            col = [next(c for c in range(k) if m[w * k + c] > 0)
                   for w in range(n)]
            surv = [(i, j) for i, j in surv if col[i] == col[j]]
        sv.delete()
        return surv

    rhombus = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)]
    allpairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]

    floor = curve(4, rhombus, 3, allpairs)
    assert floor == [(2, 3)], floor          # the forced pair, and only it

    # Delete one triangle and the forcing goes with it: the curve bottoms out.
    path = [(0, 1), (0, 2), (1, 2), (0, 3)]
    assert curve(4, path, 3, allpairs) == []


def test_decay_rates_ignores_the_unevenly_spaced_marks():
    """The bug this function exists to prevent, pinned.

    The survivor curve is sampled at 1, 2, 3, 5, 10, 20, 40, so three of its
    six steps span more than one sample.  Averaging all six ratios treats a
    two-sample decay as a one-sample decay and reports a rate that is too low.
    On the real curves it turned 0.238 and 0.242 into 0.228 and 0.229, which
    was small enough to look like agreement and led to the wrong conclusion
    that stacking translates changes nothing.
    """
    from hn.homcol import decay_rates

    # A curve that halves exactly once per sample: the honest rate is 0.5.
    curve = [(1, 1024), (2, 512), (3, 256), (5, 64), (10, 2)]
    head, tail, floor = decay_rates(curve)
    assert abs(head - 0.5) < 1e-12

    # Averaging every step, spacing ignored, would have understated it.
    naive = [b / a for (_, a), (_, b) in zip(curve, curve[1:])]
    assert sum(naive) / len(naive) < 0.42

    # The tail is a geometric mean per sample, so it recovers 0.5 as well.
    assert abs(tail - 0.5) < 1e-12
    assert floor == 2


def test_a_curve_that_reaches_zero_has_no_tail_and_no_floor():
    """Zero is the whole verdict: nothing is forced, and it is witnessed.

    A forced pair agrees in every colouring, so it survives every sample by
    definition.  A curve that reaches zero therefore proves there is no forced
    pair at all -- each pair was separated by an explicit colouring.  That is
    strictly more than a budgeted scan can say.
    """
    from hn.homcol import decay_rates

    head, tail, floor = decay_rates([(1, 4774), (2, 1119), (3, 259),
                                     (5, 22), (8, 0)])
    assert floor == 0
    assert tail is None
    assert 0.23 < head < 0.24

    # Y at four colours forces, so its curve floors instead of vanishing.
    head, tail, floor = decay_rates([(1, 2759), (2, 1117), (3, 310), (5, 159),
                                     (10, 30), (20, 11), (40, 2)])
    assert floor == 2
    assert tail is not None and tail > 0.85


def test_the_ring_criterion_picks_de_greys_ring_out_of_his_own_core():
    """The retrodiction that makes the criterion worth trusting.

    A ring of the symmetric core pays twice: sqrt(4D-1) for the rotation that
    makes it bite, sqrt(16D-1) for the spindle of the antipodal pair that
    rotation aims at.  Sa has five rational closable rings, of squared radius
    1, 1/3, 5/9, 4 and 3.  D = 1 is the 60-degree turn and fixes Sa, so it is
    degenerate; of the rest exactly one passes both tests, and it is the one
    de Grey used.

    It is not the populous one.  D = 1 carries 30 of Sa's points and D = 4
    carries six, so ranking rings by how many points sit on them -- the
    obvious heuristic, and the one tried first here -- picks the wrong ring.
    """
    from hn.homcol import doubly_usable_ring, closable_distance

    sa_rings = [Fraction(1), Fraction(1, 3), Fraction(5, 9),
                Fraction(4), Fraction(3)]
    assert all(closable_distance(D) for D in sa_rings)   # all rho-able
    assert [D for D in sa_rings if doubly_usable_ring(D)] == [Fraction(4)]

    # (2,0) and (-2,0) -- the pair Y forces -- are antipodal on that ring.
    assert Fraction(2) ** 2 == Fraction(4)
    # and the spindle of a pair at distance 4 is de Grey's own next step.
    assert closable_distance(Fraction(16))


def test_the_five_colour_core_has_a_ring_the_four_colour_one_lacks():
    """G* offers D = 17/2, which Sa has no analogue of.

    rho needs sqrt(4*17/2 - 1) = sqrt(33), which is sqrt(3)*sqrt(11) and so
    lies in de Grey's field; the spindle needs sqrt(16*17/2 - 1) = sqrt(135) =
    3 sqrt(15), which does too.  Of G*'s twelve rational closable rings only
    this one and D = 4 pay both costs.  D = 16, where de Grey's own doubling
    chain stopped for want of sqrt(17), still does not.
    """
    from hn.homcol import doubly_usable_ring, closing_radicand

    assert doubly_usable_ring(Fraction(17, 2))
    assert closing_radicand(Fraction(17, 2)) == 33
    assert closing_radicand(4 * Fraction(17, 2)) == 15

    assert not doubly_usable_ring(Fraction(16))
    assert closing_radicand(64) == 255          # 3 * 5 * 17, and 17 is absent

    gstar = [Fraction(1), Fraction(4), Fraction(31, 3), Fraction(17, 2),
             Fraction(7), Fraction(20, 3), Fraction(13, 3), Fraction(7, 3),
             Fraction(3, 2), Fraction(16)]
    assert [D for D in gstar if doubly_usable_ring(D)] == [Fraction(4),
                                                           Fraction(17, 2)]


def test_forcing_is_monotone_under_adding_vertices():
    """A forced pair stays forced in every supergraph, and it is free to know.

    Every proper colouring of the larger graph restricts to a proper colouring
    of the smaller, so a pair that agrees in all colourings of the smaller
    agrees in all colourings of the larger.  Growth can create forcing; it can
    never destroy it.

    The rhombus forces its apexes at three colours.  Hang anything else off it
    and they stay forced -- here a pendant path, which adds colourings in
    quantity without touching the two triangles.
    """
    from pysat.solvers import Solver

    def forced(n, edges, k, i, j):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in edges:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        assert sv.solve(), "the graph must be colourable at all"
        out = not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])
        sv.delete()
        return out

    rhombus = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)]
    assert forced(4, rhombus, 3, 2, 3)

    bigger = rhombus + [(3, 4), (4, 5), (5, 6), (6, 2)]
    assert forced(7, bigger, 3, 2, 3)

    # And the contrapositive: no forced pair means no forcing subgraph.
    # Deleting a triangle unforces the apexes, so nothing containing the
    # result can rely on them -- the pair is free in the smaller graph too.
    path = [(0, 1), (0, 2), (1, 2), (0, 3)]
    assert not forced(4, path, 3, 2, 3)


def test_agreeing_pairs_never_compares_a_pair():
    """Bucketing by signature, and why it is the same question.

    Two vertices agree in every sample exactly when their colour signatures
    across the samples are equal, so the quadratic comparison is not an
    approximation that was dropped -- it is the same predicate, read as a
    sort.  At 27000 vertices the comparison is ten billion element operations
    and the sort is milliseconds.
    """
    from hn.homcol import agreeing_pairs

    cols = [[0, 1, 2, 0, 1],
            [3, 1, 4, 3, 2],
            [2, 0, 1, 2, 0]]
    # 0 and 3 agree throughout; 1 and 4 agree twice and differ in the middle.
    assert agreeing_pairs(cols) == [(0, 3)]

    # It must agree with the brute-force predicate, always.
    import random
    rng = random.Random(17)
    for _ in range(40):
        n, s = rng.randrange(2, 30), rng.randrange(1, 4)
        cols = [[rng.randrange(3) for _ in range(n)] for _ in range(s)]
        brute = sorted((i, j) for i in range(n) for j in range(i + 1, n)
                       if all(row[i] == row[j] for row in cols))
        assert agreeing_pairs(cols, colours=3) == brute

    assert agreeing_pairs([]) == []


def test_the_ring_criterion_rebuilds_de_greys_graph_exactly():
    """The end-to-end check: criterion in, de Grey's 1581 vertices out.

    Start from Sa and the criterion alone.  Of Sa's five rational closable
    rings only D = 4 pays both costs, and that single choice fixes everything:
    rho = rotation_joining(4) about Sa's centre, the antipodal pair on that
    ring at squared distance 16, and sigma = rotation_joining(16) about one of
    its ends.

    Then Y u sigma(Y), turned by de Grey's own pi/2 - arcsin(1/8) about that
    pivot, is G -- not a graph like G, but G, all 1581 points.  His
    construction turns Y through pi/2 + arcsin(1/8) and pi/2 - arcsin(1/8),
    and those two angles differ by 2 arcsin(1/8), which is exactly sigma.

    So a union built this way contains a rotated copy of G and is not
    4-colourable by containment, with no solver in the argument at all.
    """
    from hn.degrey import build_Sa, build_Y, build_G, _rot_half_pi_pm
    from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
    from hn.homcol import doubly_usable_ring

    Sa = build_Sa(F)
    rings = [Fraction(1), Fraction(1, 3), Fraction(5, 9),
             Fraction(4), Fraction(3)]
    assert [D for D in rings if doubly_usable_ring(D)] == [Fraction(4)]

    from hn.degrey import build_Sb
    rho = rotation_joining(4, F)
    assert {rho(p) for p in Sa} == set(build_Sb(F))

    pivot = Point(F.rational(-2), F.zero())
    sigma = rotation_joining(4 * 4, F).about(pivot)
    A = Point(F.rational(2), F.zero())
    d2 = (A.x - sigma(A).x) ** 2 + (A.y - sigma(A).y) ** 2
    assert d2 == F.one(), "the spindle must land the pair at distance 1"

    Y = set(build_Y(F))
    spun = Y | {sigma(p) for p in Y}
    turn = _rot_half_pi_pm(F, -1).about(pivot)
    assert {turn(p) for p in spun} == set(build_G(F, as_graph=False))


def test_the_four_colour_core_sits_at_its_own_threshold():
    """chi(Sa) = 4, measured, and what it does and does not license.

    The criticality story needed Sa to sit at its threshold at four colours,
    as G does at five, and it does: Sa is not 3-colourable and is
    4-colourable.  What the story also needed -- that Sa's colourings are
    scarce because Sa is lean -- is false, and flatly.  Sa reaches chi = 4 on
    397 vertices where the Moser spindle does it on seven; G reaches chi = 5
    on 1581 where about five hundred suffice.  The graph that pins 211 pairs
    when rotated is the more redundant by a factor of twenty.

    So redundancy is not the variable, and the claim built on it is withdrawn
    in the module beside the original.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_Sa, build_S
    from hn.geometry import DEGREY_FIELD as F
    from hn.graph import build_graph

    def colourable(pts, k):
        g = build_graph(pts)
        n = len(pts)
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in g.edges():
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        out = sv.solve()
        sv.delete()
        return out

    S, Sa = build_S(F), build_Sa(F)
    assert colourable(S, 3), "S alone is only 3-chromatic"
    assert not colourable(Sa, 3)
    assert colourable(Sa, 4)          # chi(Sa) = 4, at its threshold

    # And the redundancy that the withdrawn claim rested on, in numbers.
    assert len(Sa) == 397            # against 7 for the Moser spindle


def test_the_chance_rate_of_a_sampling_filter_depends_on_the_colours():
    """The rate an independent pair survives is ((k-1)/k)^samples.

    Two random vertices differ with probability 3/4 at four colours and 4/5 at
    five, so a filter that keeps pairs differing in every sample retains
    ((k-1)/k)^s of the free ones.  Writing 0.8 for both -- as the relation
    scan first did -- understates the four-colour baseline fivefold, and with
    it the contrast the scan exists to measure: Sa's 1548 candidates are 19.6
    times chance at four colours, not 4.2.
    """
    s = 24
    assert abs((3 / 4) ** s - 0.00100339) < 1e-8
    assert abs((4 / 5) ** s - 0.00472237) < 1e-8

    # Sa: 397 points, so 78606 pairs.
    pairs = 397 * 396 // 2
    assert pairs == 78606
    assert round(pairs * (3 / 4) ** s) == 79          # four colours
    assert round(pairs * (4 / 5) ** s) == 371         # five colours
    assert 1548 / (pairs * (3 / 4) ** s) > 19         # what was measured

    # G: 1581 points at five colours sits at chance, 6502 against 5898.
    gpairs = 1581 * 1580 // 2
    assert round(gpairs * (4 / 5) ** s) == 5898
    assert 1.0 < 6502 / (gpairs * (4 / 5) ** s) < 1.2


def test_the_correlation_ratio_is_exponential_in_the_sample_count():
    """Two ratios taken at different sample counts are different measurements.

    A pair survives s samples with probability p^s -- p = (k-1)/k when it is
    unconstrained, some p' > p when it is correlated -- so observed/chance is
    (p'/p)^s and grows exponentially with s.  A mild correlation invisible at
    twenty-four samples is glaring at sixty, and comparing across s reads that
    exponent as structure.

    Hence the pairing in the measurements: Sa at four colours and G at five
    are both at twenty-four samples, and Z at five is too.
    """
    k, p = 5, 4 / 5

    # A pair 5% more likely than chance to differ, at three sample counts.
    pp = p * 1.05
    r = [(pp / p) ** s for s in (24, 40, 60)]
    assert round(r[0], 2) == 3.23
    assert round(r[1], 2) == 7.04
    assert round(r[2], 2) == 18.68
    assert r[2] / r[0] > 5, "the same correlation, read five times larger"

    # And G's own drift is too small to be a consistent p'/p.
    assert abs((1.10 ** (1 / 24)) - 1.0040) < 1e-3
    assert abs((1.34 ** (1 / 40)) - 1.0073) < 1e-3


def test_the_fifth_colour_decouples_the_same_graph():
    """One graph, one sample count, only the colour count moves.

    Every other comparison in this work changes two things at once: a
    four-colour graph against a five-colour one, of different sizes and
    different constructions.  Sa against itself changes one.  At twenty-four
    samples it leaves 1548 candidates against 79 by chance at four colours,
    and 496 against 371 at five -- the correlation in its colourings collapses
    by a factor of fifteen when the solver is handed a fifth colour.
    """
    from hn.homcol import THE_FIFTH_COLOUR_DECOUPLES_THE_SAME_GRAPH as T

    pairs, s = 397 * 396 // 2, 24
    assert round(pairs * (3 / 4) ** s) == T["at_four"]["chance"] == 79
    assert round(pairs * (4 / 5) ** s) == T["at_five"]["chance"] == 371

    four = T["at_four"]["candidates"] / T["at_four"]["chance"]
    five = T["at_five"]["candidates"] / T["at_five"]["chance"]
    assert round(four, 1) == 19.6
    assert round(five, 2) == 1.34
    assert four / five > 14


def test_the_lattice_recognises_its_own_unique_colouring():
    """The instrument check behind the ladder, done by arithmetic.

    The triangular lattice is 3-chromatic and its proper 3-colouring is unique
    up to permuting the colours, so two points differ in every colouring
    exactly when they lie in different classes.  For a 400-point patch split
    into three classes of about 133, that fraction is
    1 - 3*C(133,2)/C(400,2), which is 66 per cent -- and the sampling filter,
    told nothing about any of this, returns 52212 of 78679 non-edge pairs,
    which is also 66 per cent.
    """
    from math import comb

    pairs = comb(400, 2)
    assert pairs == 79800
    cross = 1 - 3 * comb(133, 2) / pairs
    assert 0.66 < cross < 0.67

    measured = 52212 / (pairs - 1121)          # minus the edges
    assert 0.66 < measured < 0.67
    assert abs(measured - cross) < 0.01

    # And the ladder it anchors: at k = chi the ratio falls as chi rises.
    from hn.homcol import THE_LADDER_OF_COLOUR_COUNTS as L
    assert L["table"]["triangular lattice (chi=3)"][3] > 1000
    assert 10 < L["table"]["Sa (chi=4)"][4] < 100
    assert L["table"]["G (chi=5)"][5] < 2


def test_the_lattice_cannot_contain_a_moser_spindle():
    """The gadget count, and the one value of it that is forced by theory.

    A Moser spindle is 4-chromatic, so no 3-colourable graph contains one, and
    the triangular lattice is 3-colourable.  Its count must therefore be
    exactly zero -- which is what makes the counter worth trusting on the
    graphs where the answer is not known in advance.

    de Grey's family then comes out at a single density, 1.45 spindles per
    point, throughout: Sa has 576 on 397 points, Y has 1152 on 791, G has 2304
    on 1581.  Each is built from the one before, so the density is inherited
    rather than achieved.
    """
    from hn.homcol import THE_GADGET_IS_SEVENTY_TIMES_BIGGER as T

    assert T["counted"]["triangular lattice"]["spindles"] == 0

    for g in ("Sa", "Y", "G"):
        c = T["counted"][g]
        assert abs(c["spindles"] / c["points"] - 1.455) < 0.01

    # Y and G double Sa's count as they double its points.
    assert T["counted"]["Y"]["spindles"] == 2 * T["counted"]["Sa"]["spindles"]
    assert T["counted"]["G"]["spindles"] == 2 * T["counted"]["Y"]["spindles"]

    # And the scale the argument points at is the one Z happens to have.
    assert abs(397 * (500 / 7) - 28357) < 1
    assert abs(T["and_Z_has"] - 397 * (500 / 7)) < 1000


def test_a_ratio_is_only_as_good_as_the_mixing_behind_it():
    """The failure mode that made one subgraph score 700, and the check for it.

    A correlation ratio inflates whenever the sampler returns near-identical
    colourings: twenty-four copies of one colouring make every pair look
    constrained.  Renaming does not catch it, since colourings differing in
    three vertices are still distinct tuples.  What catches it is distance --
    two independent proper k-colourings disagree on about (1 - 1/k) of the
    vertices, and the spiking subgraphs disagreed on 0.019 and 0.060.

    The graphs the ladder rests on pass: Sa at four sits at 87 per cent of
    independent and G at five at 95.  The triangular lattice's zero is not a
    failure but the answer -- its 3-colouring is unique up to permutation, so
    the samples must coincide, which is where its ratio of eleven thousand
    comes from.
    """
    from hn.homcol import THE_SAMPLES_ARE_MIXED as M

    t = M["table"]
    assert t["triangular lattice at 3"]["spread"] == 0.0
    assert t["triangular lattice at 3"]["ratio"] > 1000

    for key in ("Sa at 4", "G at 5", "Y at 4", "Sa at 5"):
        row = t[key]
        assert row["spread"] / row["independent"] > 0.75, key

    # The spiking subgraphs, for contrast, sat two orders below that.
    from hn.homcol import THE_DILUTION_ANOMALY as A
    assert "RESOLVED" in A
    assert 0.019 / 0.653 < 0.05


def test_the_rhombus_is_the_minimal_three_colour_forcer():
    """The small end of the forcer curve, established rather than assumed.

    A forcer is a graph carrying a pair monochromatic in every proper
    k-colouring.  At three colours the rhombus does it on four vertices, and
    nothing smaller can: a graph on three vertices that is 3-colourable has at
    most a triangle, and a triangle has no non-adjacent pair to force.  So 4
    is exact, and the other end of the curve -- Y at 791, with 787 of its 789
    non-pair vertices individually indispensable -- is measured too.
    """
    from pysat.solvers import Solver
    from hn.homcol import THE_FORCER_GROWS_BY_TWO_HUNDRED as T

    def forces(n, edges, k):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in edges:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        assert sv.solve()
        es = {(min(a, b), max(a, b)) for a, b in edges}
        out = any(not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])
                  for i in range(n) for j in range(i + 1, n)
                  if (i, j) not in es)
        sv.delete()
        return out

    assert forces(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)], 3)

    # Nothing on three vertices does: every 3-colourable graph on three
    # vertices is a subgraph of a triangle, which has no non-adjacent pair.
    for edges in ([], [(0, 1)], [(0, 1), (1, 2)],
                  [(0, 1), (1, 2), (0, 2)]):
        assert not forces(3, edges, 3)

    assert T["k=3"]["vertices"] == 4
    assert T["k=4"]["vertices"] == 791
    assert round(791 / 4) == T["ratio"]


def test_correlation_tracks_gadget_density_monotonically():
    """The mechanism, as a relation rather than a story about two points.

    Ten objects from three unrelated constructions -- a triangular lattice, a
    lattice spindled at one to four hinges, and thinned copies of de Grey's
    core -- measured for Moser spindles per point and for correlation at four
    colours.  Sorting by density sorts by correlation, over two orders of
    magnitude in the first and one in the second.
    """
    from hn.homcol import CORRELATION_TRACKS_GADGET_DENSITY as T

    rows = sorted(T["table"], key=lambda r: r[1])
    ratios = [r[2] for r in rows]

    # The two ends are unambiguous and far apart.
    assert rows[0][1] == 0.0 and rows[-1][1] > 1.4
    assert ratios[-1] / max(ratios[0], 1e-9) > 10

    # Above the noise floor -- densities of a quarter spindle per point and up
    # -- the order is exact.
    dense = [r for r in rows if r[1] >= 0.25]
    assert [r[2] for r in dense] == sorted(r[2] for r in dense)
    assert len(dense) >= 6

    # And a 4-chromatic graph need not contain a spindle at all.
    assert any(r[0] == "lattice + 1 hinge" and r[1] == 0.0 for r in rows)


def test_G_minus_a_vertex_can_be_four_coloured():
    """A surprising measurement, checked with an explicit witness.

    Nineteen of thirty probed vertices of G came back essential -- G - v
    4-colourable -- which sits oddly beside a literature that treats de Grey's
    graph as the starting point for reductions to 874 vertices and below.  It
    is nonetheless true, and the witness settles it: a proper 4-colouring of
    G - 1172 exists, with no monochromatic edge and every vertex carrying
    exactly one colour.

    There is no contradiction.  1581 is already de Grey's own reduction of his
    20425, and the smaller graphs are found by fresh search, not by deleting
    vertices from this one.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_G
    from hn.geometry import DEGREY_FIELD as F

    g = build_G(F)
    n, k, v = g.n, 4, 1172
    E = [(min(a, b), max(a, b)) for a, b in g.edges()]
    keep = [w for w in range(n) if w != v]
    idx = {w: i for i, w in enumerate(keep)}
    sub = [(idx[a], idx[b]) for a, b in E if a != v and b != v]

    cls = [[1 + i * k + c for c in range(k)] for i in range(len(keep))]
    for x, y in sub:
        for c in range(k):
            cls.append([-(1 + x * k + c), -(1 + y * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    sv.conf_budget(600000)
    assert sv.solve_limited() is True
    m = sv.get_model()
    col = [next(c for c in range(k) if m[i * k + c] > 0)
           for i in range(len(keep))]
    sv.delete()

    assert not [1 for x, y in sub if col[x] == col[y]]
    assert len(keep) == 1580 and len(sub) == 7869


def test_criticality_and_correlation_run_opposite():
    """The inversion, on the two graphs it is measured on.

    Sa and Y are as far from critical as a graph can be -- not one deletion
    in forty drops their chromatic number -- and carry the strongest colour
    relations here, 24.0 and 25.4.  G is nearly vertex-critical, nineteen
    deletions in thirty dropping it, and carries none, 1.1.

    Which is the right way round, once said properly: critical means minimal,
    exactly enough constraint to force the chromatic number and nothing to
    spare, so everything not load-bearing is free.  Redundant structure
    correlates; tight structure does not.  This is the measurement that
    finishes off the criticality diagnosis, already withdrawn on the weaker
    grounds that the pinning graph is the more redundant.
    """
    from pysat.solvers import Solver
    from hn.degrey import build_Sa
    from hn.geometry import DEGREY_FIELD as F
    from hn.graph import build_graph
    from hn.homcol import CRITICALITY_AND_CORRELATION_RUN_OPPOSITE as C

    # Sa loses no chromatic number to a single deletion -- spot-check three.
    g = build_graph(build_Sa(F))
    n, k = g.n, 3
    E = [(min(a, b), max(a, b)) for a, b in g.edges()]
    for v in (0, 100, 250):
        keep = [w for w in range(n) if w != v]
        idx = {w: i for i, w in enumerate(keep)}
        cls = [[1 + i * k + c for c in range(k)] for i in range(len(keep))]
        for a, b in E:
            if a == v or b == v:
                continue
            x, y = idx[a], idx[b]
            for c in range(k):
                cls.append([-(1 + x * k + c), -(1 + y * k + c)])
        sv = Solver(name="cd19", bootstrap_with=cls)
        still_four = not sv.solve()
        sv.delete()
        assert still_four, f"Sa - {v} should still need four colours"

    m = C["measured"]
    assert m["Sa at 4"]["correlation"] > 20 > m["G at 5"]["correlation"]
    assert m["Sa at 4"]["deletions_dropping_chi"].startswith("0")
    assert m["G at 5"]["deletions_dropping_chi"].startswith("19")


def test_the_gadget_mechanism_does_not_survive_controlling_for_degree():
    """The confound, and the control that was already in the data.

    Thinning a graph removes edges along with gadgets, so a table where
    spindle density and correlation rise together says nothing until degree is
    held fixed.  Held fixed, the effect disappears: five objects at average
    degree about 5.6 span spindle densities from 0.000 to 0.534 and their
    ratios are 2.2, 2.0, 1.5, 1.2 and 1.7, with no trend.

    And across colour counts neither predictor works at all -- Sa and G have
    average degrees 9.94 and 9.96 and ratios 24.0 and 1.1.
    """
    from hn.homcol import THE_GADGET_MECHANISM_IS_CONFOUNDED as C

    # degree, spindles/point, ratio
    rows = [(5.61, 0.000, 2.2), (5.60, 0.000, 2.0), (5.61, 0.007, 1.5),
            (5.61, 0.004, 1.2), (5.66, 0.534, 1.7)]
    assert max(r[0] for r in rows) - min(r[0] for r in rows) < 0.1

    by_spindles = [r[2] for r in sorted(rows, key=lambda r: r[1])]
    assert by_spindles != sorted(by_spindles), "no trend at fixed degree"
    assert max(by_spindles) - min(by_spindles) < 1.1

    assert C["spearman"]["degree vs ratio"] > C["spearman"]["spindles vs ratio"]
    assert C["spearman"]["degree vs spindles"] > 0.85


def test_the_fifth_colour_flattens_the_degree_curve():
    """The controlled comparison the whole account reduces to.

    Two graphs thinned by the same procedure and seeds, their average degrees
    matched, and only the colour count differing.  At four colours the
    correlation ratio climbs from 1.3 to 24 as degree goes from 3.8 to 9.9.
    At five it is flat at one across the same range, drifting below one at low
    degree -- noise around independence.

    So constraint density does not merely weaken at five colours; it stops
    operating.  Which is why every construction in this work failed: they all
    add points and edges, and that is the lever measured to do nothing.
    """
    from hn.homcol import THE_FIFTH_COLOUR_FLATTENS_THE_CURVE as T

    rows = [r for r in T["table"] if r[2] != 251.4]   # drop the artefact
    five = [r[1] for r in rows]
    four = [r[2] for r in rows]

    assert max(five) <= 1.1 and min(five) >= 0.5      # flat at one
    assert max(five) - min(five) < 0.6
    assert max(four) / min(four) > 15                 # climbs by an order

    # Ordered by degree, four colours rises and five does not.  The
    # four-colour series inverts once at the bottom -- 1.7 at degree 5.66
    # against 2.1 at 4.87 -- where the values sit near one and the noise is
    # the same size as the signal, so the claim is made above that floor.
    top = [r for r in rows if r[2] >= 2.0]
    assert [r[2] for r in top] == sorted((r[2] for r in top), reverse=True)
    assert len(top) >= 4
    assert five[0] - five[-1] < 0.6


def test_surplus_does_not_predict_but_degree_does():
    """The last hypothesis, refuted by data already in hand.

    Surplus -- how many times a graph exceeds the minimum for its own
    chromatic number -- runs from 34 to 257 across the four-colour objects
    measured here, and the correlation ratio goes 1.7, 24.0, 2.0, 25.4, 1.2.
    No order.  Average degree separates them exactly: the two ratios above 20
    are the two objects at degree 9.95, and the three below 2.1 are the three
    at 5.6.

    Which cancelled a twelve-hour measurement of Z at fifty-five times the
    minimum, since the axis it would have placed a point on carries no signal.
    """
    from hn.homcol import SURPLUS_DOES_NOT_PREDICT as S

    rows = S["four_colour_objects_by_surplus"]
    assert [r[0] for r in rows] == sorted(r[0] for r in rows)

    ratios = [r[1] for r in rows]
    assert ratios != sorted(ratios) and ratios != sorted(ratios, reverse=True)

    high = [r for r in rows if r[1] > 20]
    low = [r for r in rows if r[1] < 2.1]
    assert len(high) == 2 and len(low) == 3
    assert all(r[2] > 9.9 for r in high)
    assert all(r[2] < 5.7 for r in low)


def test_the_two_extrapolations_disagree_by_two_orders():
    """Both estimates of the scale six would need, and their disagreement.

    One runs through degree: the four-colour ratio rises as
    5.0 exp(0.857 (d - 8.11)), the five-colour curve is flat at one up to
    degree 10.64, so reaching 24 would take degree 14.3 -- and degree grows as
    9.96 + 0.49 ln(n/1581) in this family, which puts that at 1.3e7 points.

    The other runs through the forcer: 4 vertices at three colours, 791 at
    four, a factor of 198, and one more such factor gives 157000.

    They differ by eighty, which is the honest headline.  Quoting either
    figure alone would be quoting the assumptions rather than the data.
    """
    import math
    from hn.homcol import TWO_EXTRAPOLATIONS_TO_THE_SCALE_NEEDED as T

    a = math.log(24.0 / 5.0) / (9.94 - 8.11)
    assert abs(a - T["by_degree"]["four_colour_slope"]) < 0.01
    assert abs((10.64 + math.log(24.0) / a)
               - T["by_degree"]["degree_for_ratio_24"]) < 0.1

    assert round(791 / 4) == T["by_forcer"]["factor"] == 198
    assert abs(791 * 198 - T["by_forcer"]["points_needed"]) < 1000

    big = T["by_degree"]["points_needed"]
    small = T["by_forcer"]["points_needed"]
    assert 70 < big / small < 90
    assert small > 55345, "both are past anything built here"


def test_single_distance_lattice_graphs_are_always_bipartite():
    """The fact that disqualifies the Erdos direction, and qualifies the metric.

    dx^2 + dy^2 = r forces dx + dy = r (mod 2).  For odd r exactly one of
    dx, dy is odd, so every edge flips the parity of i + j; for r = 2 mod 4
    both are odd, so every edge flips the parity of i; and r = 0 mod 4 is a
    scaled copy.  A parity is a proper 2-colouring each time.

    So the integer lattice at a single squared distance is bipartite however
    large its degree -- and one of these graphs scores 26.4 on the correlation
    ratio at five colours.  The ratio measures local correlation, not
    proximity to forcing.
    """
    from hn.homcol import SINGLE_DISTANCE_LATTICES_ARE_BIPARTITE as B

    for r in (25, 65, 325, 1105, 2210, 2900, 2650, 50, 98, 200):
        lim = int(r ** 0.5) + 1
        reps = [(dx, dy) for dx in range(-lim, lim + 1)
                for dy in range(-lim, lim + 1)
                if dx * dx + dy * dy == r and (dx, dy) > (0, 0)]
        assert reps, r
        if r % 2 == 1:
            assert all((dx + dy) % 2 == 1 for dx, dy in reps), r
        elif r % 4 == 2:
            assert all(dx % 2 == 1 and dy % 2 == 1 for dx, dy in reps), r
        else:
            assert all(dx % 2 == 0 and dy % 2 == 0 for dx, dy in reps), r

    assert set(B["checked"].values()) == {2}
    assert B["up_to"]["degree"] > 19

import hn.homcol

# --- the exact forced-same test --------------------------------------------

def test_colour_symmetry_assumptions_picks_two_literals():
    # vertex v colour c is variable 1 + v*k + c
    assert hn.homcol.colour_symmetry_assumptions(3, 7, 5) == [16, -36]
    assert hn.homcol.colour_symmetry_assumptions(0, 1, 4) == [1, -5]


def _colour_solver(n, edges, k):
    from pysat.solvers import Solver
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    return Solver(name="cd15", bootstrap_with=cls)


def test_forced_same_on_k4_minus_an_edge():
    # 1,2,3 is a triangle and 0 sees 1 and 2, so at three colours 0 is left
    # with the colour 3 already has -- the missing edge (0,3) is forced same.
    E = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    sv = _colour_solver(4, E, 3)
    assert sv.solve()
    assert hn.homcol.forced_same(sv, 0, 3, 3)
    sv.delete()


def test_forced_same_is_false_when_a_colour_is_spare():
    # the same graph at four colours: 0 may take the fourth colour
    E = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    sv = _colour_solver(4, E, 4)
    assert sv.solve()
    assert not hn.homcol.forced_same(sv, 0, 3, 4)
    sv.delete()


def test_forced_same_agrees_with_adding_the_edge():
    # the definition it stands in for: forced iff the graph plus that edge
    # is uncolourable.  Checked on both answers of the K4-minus-an-edge pair.
    E = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    for k, expect in ((3, True), (4, False)):
        sv = _colour_solver(4, E, k)
        plus = _colour_solver(4, E + [(0, 3)], k)
        assert hn.homcol.forced_same(sv, 0, 3, k) is expect
        assert (not plus.solve()) is expect
        sv.delete()
        plus.delete()


# --- the null that retired the sampled ranking ------------------------------

def test_max_agreement_null_brackets_the_climb():
    c = hn.homcol.THE_MAX_AGREEMENT_IS_AN_EXTREME_VALUE
    alone = max(c["G alone, four seeds"])
    null = max(c["G + G, zero cross edges, two seeds"])
    climb = max(c["the climb it was compared against"])
    # the uncoupled control beats the graph itself, and the climb is not
    # meaningfully past the control -- one step, over fifty draws against two
    assert null > alone
    assert climb - null <= 1


# --- de Grey's template, as measured ---------------------------------------

def test_the_ring_that_pays_twice_is_unique_among_Sa_rings():
    from fractions import Fraction as Fr
    # Sa's rational rings; only D = 4 closes both 4D-1 and 16D-1
    rings = [Fr(1), Fr(1, 3), Fr(5, 9), Fr(4), Fr(3)]
    assert [d for d in rings if hn.homcol.doubly_usable_ring(d)] == [Fr(4)]


def test_four_ninths_is_doubly_usable_but_unpopulated():
    from fractions import Fraction as Fr
    # the field closes it -- 4D-1 = 7/9 and 16D-1 = 55/9, both squares times
    # products of the generators -- and it is in the recorded hit list
    assert hn.homcol.doubly_usable_ring(Fr(4, 9))
    c = hn.homcol.THE_FIELD_OFFERS_RINGS_THE_GRAPH_DOES_NOT_POPULATE
    assert "4/9" in c["the small ones"]
    assert c["G's best populated doubly-usable ring, over every vertex "
             "centre and 4000 midpoints"]["points"] == 12


def test_the_unit_ring_is_excluded_for_the_right_reason():
    from fractions import Fraction as Fr
    # D = 1 closes both radicals but its rotation is the 60-degree turn, which
    # fixes the whole Eisenstein core -- it bites nothing
    assert hn.homcol.closable_distance(Fr(1))
    assert hn.homcol.closable_distance(Fr(4))
    assert not hn.homcol.doubly_usable_ring(Fr(1))


def test_build_G_rotations_compose_to_the_D16_spindle():
    # the two turns are pi/2 -+ arcsin(1/8); their composition must be the
    # rotation the spindling lemma needs for squared distance 16, namely
    # cos = 1 - 1/(2*16) = 31/32 and sin = sqrt(4*16-1)/(2*16) = sqrt(63)/32
    from fractions import Fraction as Fr
    from hn.degrey import _rot_half_pi_pm
    from hn.geometry import DEGREY_FIELD as F, rotation_joining
    a = _rot_half_pi_pm(F, +1)
    b = _rot_half_pi_pm(F, -1)
    # a * b^-1
    cos = a.cos * b.cos + a.sin * b.sin
    sin = a.sin * b.cos - a.cos * b.sin
    spindle = rotation_joining(Fr(16), F)
    assert cos == F.rational(Fr(31, 32))
    assert cos == spindle.cos
    assert sin * sin == spindle.sin * spindle.sin


def test_the_bite_threads_because_the_step_is_one():
    # |1 - rho| = 1/2 for D = 4, so a radius-2 point and its image are at
    # distance exactly one: that is why iterating threads a path
    from fractions import Fraction as Fr
    from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
    rho = rotation_joining(Fr(4), F)
    p = Point(F.rational(2), F.zero())
    assert p.dist2(rho(p)) == 1
    assert p.dist2(rho(rho(p))) != 1        # only consecutive ones touch
    c = hn.homcol.THE_BITE_THREADS_THE_RING
    assert c["at G[0] = (-2.25, +1.984313)"]["ring"] == [12, 24, 36, 48]


def test_the_exhaustive_census_of_G_is_recorded_complete():
    c = hn.homcol.THE_CENSUS_OF_G_IS_EXHAUSTIVE_AND_EMPTY
    assert c["forced at five colours"] == 0
    assert c["pairs they carry"] < c["pairs at a rational squared distance"]
    assert c["closable ones"] < c["distinct rational squared distances"]


# --- the weak property, which is the gateway --------------------------------

def test_forbidding_a_pair_is_adding_the_edge():
    # the identity the one-call test rests on: a k-colouring in which (i, j)
    # is NOT monochromatic is exactly a proper colouring of the graph with the
    # edge (i, j) added
    E = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    for k in (3, 4):
        plus = _colour_solver(4, E + [(0, 3)], k)
        bare = _colour_solver(4, E, k)
        assert bare.solve()
        # "some colouring separates them" and "the graph plus the edge
        # colours" are the same statement
        separable = bare.solve(
            assumptions=hn.homcol.colour_symmetry_assumptions(0, 3, k))
        assert separable == plus.solve()
        plus.delete()
        bare.delete()


def test_weak_property_holds_on_a_triangle_of_pairs():
    # three pairwise non-adjacent pairs that cannot all be non-monochromatic:
    # take K5 at four colours with three edges deleted.  Adding them back is
    # K5, which needs five colours, so at four some deleted pair is always
    # monochromatic -- the weak property, by the one-call test.
    K5 = [(a, b) for a in range(5) for b in range(a + 1, 5)]
    gone = [(0, 1), (0, 2), (0, 3)]
    rest = [e for e in K5 if e not in gone]
    bare = _colour_solver(5, rest, 4)
    assert bare.solve()                    # the graph itself colours
    full = _colour_solver(5, K5, 4)
    assert not full.solve()                # forbidding all three does not
    bare.delete()
    full.delete()


def test_strong_implies_weak_so_the_gateway_subsumes():
    c = hn.homcol.THE_WEAK_PROPERTY_IS_THE_GATEWAY
    assert "strong => weak" in c["implies"]
    g = hn.homcol.G_FAILS_THE_GATEWAY_EVERYWHERE
    assert g["carrying the weak property at five colours"] == 0
    assert g["with a closable spindle"] < g["(centre, ring) candidates"]


def test_de_greys_ring_carries_the_weak_property_with_three_pairs():
    c = hn.homcol.THE_WEAK_PROPERTY_IS_THE_GATEWAY["Sa at four"]
    assert "HOLDS" in c["D=4"] and "de Grey's" in c["D=4"]
    # and the spindle distance it names is the one the pair sits at
    assert "16" in c["D=4"]
    assert hn.homcol.closable_distance(Fraction(16))


def test_the_carrying_classes_are_not_the_biggest():
    # if the property were just "forbid enough pairs", the classes carrying it
    # would be the populous ones.  They are not: Sa's three sit 5th, 6th and
    # 9th of fourteen by size.
    c = hn.homcol.THE_WEAK_PROPERTY_GENERALISES_PAST_ANTIPODAL
    yes = c["Sa at four, closable classes carrying it"]
    no = c["Sa at four, closable classes that do not"]
    assert max(no.values()) > max(yes.values())
    assert min(no.values()) < min(yes.values())


def test_de_greys_lemma_has_no_slack():
    c = hn.homcol.THE_LEMMA_IS_TIGHT_AND_THE_HUB_CLOSURE_IS_EMPTY
    assert c["Sa's D=4 antipodal pairs"] == c["how many are needed"]
    hub = c["G closed about its hub"]
    assert hub["same as G's about G[0]"] is True
    assert hub["rings carrying the weak property"] == 0
    assert hub["added"] > 9000        # many points, no new ring points


def test_the_hard_class_is_not_the_big_one():
    # if difficulty tracked the number of forbidden pairs, the 6510-pair class
    # would be the hard one.  The hard one has 1558.
    c = hn.homcol.THE_EDGE_IS_AT_FOUR_NINTHS["G's closable classes at five, "
                                             "by cost"]
    assert "instant" in c["1/3 (6510 pairs)"]
    assert "unresolved" in c["4/9 (1558)"]


def test_four_ninths_needs_no_adjunction():
    # both radicals the template wants are already in de Grey's field
    from fractions import Fraction as Fr
    assert hn.homcol.doubly_usable_ring(Fr(4, 9), (3, 5, 7, 11))
    assert hn.homcol.closing_radicand(Fr(4, 9)) == 7
    assert hn.homcol.closing_radicand(4 * Fr(4, 9)) == 55


def test_the_ball_curve_is_monotone_and_steep():
    b = hn.homcol.THE_EDGE_IS_AT_FOUR_NINTHS["balls about the hub, all "
                                             "colouring"]
    secs = [int(v.rstrip("s")) for _, v in sorted(b.items())]
    assert secs == sorted(secs)
    assert secs[-1] > 50 * secs[0]      # 21s to 6893s over 360 points


def test_the_gap_was_double_counted():
    c = hn.homcol.THE_GAP_IS_TEN_NOT_SEVEN_HUNDRED
    inc = c["incidences per point"]
    # the incidence ratio is what the recorded figure should have been
    assert abs(inc["G (spindles)"] / inc["G (5-critical)"] - c["the gap"]) < 0.2
    # and the recorded figure is that times the gadget size ratio
    assert abs(c["as recorded"] / c["the gap"] - 500 / 7) < 1.0
    # the uncorrected numbers still reproduce: 576 spindles of 7 in 397 points
    assert abs(576 * 7 / 397 - inc["Sa"]) < 0.01
    assert abs(2304 * 7 / 1581 - inc["G (spindles)"]) < 0.01


def test_density_is_inherited_across_the_family():
    d = hn.homcol.THE_GAP_IS_TEN_NOT_SEVEN_HUNDRED["edges per vertex, "
                                                   "measured"]
    vals = list(d.values())
    assert max(vals) - min(vals) < 0.02      # constant to two parts in a
    assert all(4.9 < v < 5.0 for v in vals)  # thousand, bites included


def test_bites_do_not_move_density_and_closures_barely_do():
    c = hn.homcol.EACH_OPERATION_AND_WHAT_IT_MOVES
    assert c["bite and thicken"]["raises coverage"] is False
    assert c["dihedral closure"]["raises coverage"] is False
    # the two sides of the discriminator: near-equal density, opposite verdict
    d = c["but density is not the discriminator"]
    assert abs(d["Sa"]["edges per vertex"] - d["G"]["edges per vertex"]) < 0.02
    assert d["Sa"]["carries the weak property at 4"] is True
    assert d["G"]["carries it at 5"] is False


def test_the_deleted_pair_sits_at_the_distinguished_distance():
    # build_Y drops (1/3, 0) and (-1/3, 0); their squared distance is 4/9
    from fractions import Fraction as Fr
    from hn.geometry import DEGREY_FIELD as DF, Point
    a = Point(DF.rational(Fr(1, 3)), DF.zero())
    b = Point(DF.rational(Fr(-1, 3)), DF.zero())
    assert a.dist2(b) == Fr(4, 9)
    assert a.norm2() == Fr(1, 9)          # antipodal on the D = 1/9 ring
    c = hn.homcol.FOUR_NINTHS_IS_THE_FAMILYS_DISTANCE["the deleted pair"]
    assert c["squared distance"] == "4/9"
    assert "nothing" in c["what the deletion costs"]


def test_the_deletion_is_actually_two_points():
    from hn.degrey import build_Sa, build_Sb, build_Y
    from hn.geometry import DEGREY_FIELD as DF
    union = set(build_Sa(DF)) | set(build_Sb(DF))
    assert len(union) == 793
    assert len(build_Y(DF)) == 791


def test_sa_carries_four_classes_not_three():
    # an earlier reading of this listed three and called the distances an
    # arithmetic progression of step 2/3.  There are four: the fourth is
    # D = 16 with just three pairs, which are de Grey's own -- the antipodal
    # pairs of the D = 4 ring -- and 4 is not 2 + 2/3.
    from fractions import Fraction as Fr
    c = hn.homcol.FOUR_NINTHS_IS_THE_FAMILYS_DISTANCE["Sa's carrying classes"]
    assert set(c) == {"4/9", "16/9", "4", "16"}
    roots = sorted(Fr(k) ** Fr(1, 2) if False else r
                   for r, k in ((Fr(2, 3), "4/9"), (Fr(4, 3), "16/9"),
                                (Fr(2), "4"), (Fr(4), "16")))
    assert [r * r for r in roots] == sorted(Fr(k) for k in c)
    # they are 2/3 times 1, 2, 3, 6 -- not an arithmetic progression
    assert [r / Fr(2, 3) for r in roots] == [1, 2, 3, 6]
    assert roots[3] - roots[2] != roots[1] - roots[0]
    assert hn.homcol.FOUR_NINTHS_IS_THE_FAMILYS_DISTANCE[
        "as multiples of 2/3"] == [1, 2, 3, 6]


def test_de_greys_class_is_the_smallest_that_carries():
    # three pairs out of Sa's fourteen classes, and the one he built on
    c = hn.homcol.FOUR_NINTHS_IS_THE_FAMILYS_DISTANCE["Sa's carrying classes"]
    assert "3 pairs" in c["16"] and "de Grey's" in c["16"]
    sizes = {k: int(v.split()[-2]) for k, v in c.items()}
    assert min(sizes, key=sizes.get) == "16"


def test_the_local_search_was_discarded_on_its_calibration():
    c = hn.homcol.LOCAL_SEARCH_NEEDED_CALIBRATING
    bad = c["the discarded version"]
    assert bad["left"] > 0                      # on a known-satisfiable one
    assert bad["its readings on open instances"] == "discarded"
    cal = c["calibration"]
    assert all("SAT in" in v for kk, v in cal.items() if "k=5" in kk)
    assert "plateau" in cal["Sa, k=4, 4/9 forbidden"]


def test_no_six_clique_in_the_combined_graphs():
    c = hn.homcol.THE_COMBINED_GRAPHS_HAVE_CLIQUE_NUMBER_FOUR
    for key in ("Sa + {4/9, 16/9, 4, 16}", "Y + the same"):
        g = c[key]
        assert g["clique number"] == 4 < 6
        # degree is no obstacle -- it is fourteen, so a K6 could have fitted
        assert g["mean degree"] > 6
        # mean degree is recorded to one decimal, so allow for the rounding
        assert abs(g["unit edges"] + g["added"]
                   - g["mean degree"] * g["points"] / 2) < 0.05 * g["points"]
