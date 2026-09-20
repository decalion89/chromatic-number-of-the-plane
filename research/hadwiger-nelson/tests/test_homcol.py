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
