"""Forced colour relations, and the ladder chi >= 6 reduces to."""
import pytest

from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.forced import (ColourRelations, centre_pressure, find_rainbow,
                       min_colours_on)
from hn.graph import build_graph
from hn.mixed import joint_core_union


def test_edges_are_forced_different():
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        assert rel.colourable
        for u, v in list(g.edges())[:20]:
            assert rel.different(u, v)
            assert not rel.can_share(u, v)
    finally:
        rel.close()


def test_min_colours_recovers_the_chromatic_number():
    """On the whole vertex set the measurement is just chi."""
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        assert min_colours_on(rel, range(g.n)) == 4
        assert min_colours_on(rel, []) == 0
        assert min_colours_on(rel, [0]) == 1
    finally:
        rel.close()


def test_vacuous_queries_are_refused():
    """A graph with no k-colouring makes every query UNSAT; say so."""
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 3)
    try:
        assert not rel.colourable
    finally:
        rel.close()
    with pytest.raises(ValueError):
        find_rainbow(g, 3, 3)


def test_unit_circle_graphs_are_bipartite():
    """Two points of a unit circle are adjacent only at 60 degrees, and the
    60-degree orbit closes after six steps, so every component of a circle
    graph is a path or a 6-cycle and no circle ever needs three colours on its
    own account. Any third colour has to be forced by the ambient graph --
    which is why every neighbourhood of de Grey's G comes back at exactly two,
    however large: the two hubs of degree 60 included.
    """
    for g in (build_graph(joint_core_union()), build_G()):
        for v in range(g.n):
            circle = set(g.adj[v])
            colour = {}
            for a in sorted(circle):
                if a in colour:
                    continue
                colour[a], stack = 0, [a]
                while stack:
                    x = stack.pop()
                    for y in g.adj[x] & circle:
                        if y not in colour:
                            colour[y] = 1 - colour[x]
                            stack.append(y)
                        else:
                            assert colour[y] != colour[x], \
                                f"circle of {v} is not bipartite at {x},{y}"


@pytest.mark.slow
def test_degrey_G_has_no_local_pressure_at_five():
    """Every one of the 1581 neighbourhoods is 2-colourable at k=5.

    The ladder needs three: at three the pivot's colour is confined to two,
    which is the joint core the two-orbit block closes. This pins the measured
    fact that G supplies none of it.
    """
    g = build_G()
    rel = ColourRelations(g, 5)
    try:
        assert rel.colourable
        seen = {min_colours_on(rel, sorted(g.adj[v])) for v in range(g.n)}
        assert seen == {2}
    finally:
        rel.close()


def test_centre_pressure_reads_the_ladder():
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        hub = max(range(g.n), key=lambda v: len(g.adj[v]))
        r = centre_pressure(rel, sorted(g.adj[hub]))
        assert r["pressure"] + r["free_for_centre"] == 4
        assert r["smallest_possible_core"] == 4 - r["pressure"]
    finally:
        rel.close()


# -- the pressure theorem -------------------------------------------------

def test_pressure_bound_matches_the_construction():
    """A core of r needs pressure >= k - r; the built core of two saturates it.

    On the four points at three colours the pivot sees one neighbour, so
    pressure is 1 = k - 2, the theorem permits a core of two, and a separation
    query confirms one. The bound is tight exactly where the construction is.
    """
    from hn.forced import max_core_bound, pressure
    from hn.mixed import joint_core_configuration
    from hn.spindle import SeparationTest

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    ia, ib = g.vertices.index(pts[2]), g.vertices.index(pts[3])
    rel = ColourRelations(g, 3)
    try:
        assert pressure(rel, bp) == 1
        assert max_core_bound(rel, bp) == 2
    finally:
        rel.close()
    st = SeparationTest(g, 3, bp, [ia, ib])
    try:
        assert not st.run(subset=[ia, ib])[0]
    finally:
        st.close()


def test_pressure_forbids_small_cores_where_it_is_low():
    """The contrapositive, checked directly: squeeze the neighbourhood into
    k - r - 1 colours and every set of r vertices fails to be a core."""
    from itertools import combinations

    from hn.forced import pressure
    from hn.spindle import SeparationTest

    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        hub = max(range(g.n), key=lambda v: len(g.adj[v]))
        assert pressure(rel, hub) == 1          # so no core of size 1 or 2
    finally:
        rel.close()
    others = [v for v in range(g.n) if v != hub and v not in g.adj[hub]]
    st = SeparationTest(g, 4, hub, others)
    try:
        for pair in list(combinations(others, 2))[:40]:
            assert st.run(subset=list(pair))[0], f"{pair} should not be a core"
    finally:
        st.close()


# -- criticality makes a graph inert ---------------------------------------

def test_critical_graphs_have_no_core_at_all():
    """Colour G - p with k-1 and give p the kth: p is then alone in its
    colour, so c(p) lies outside c(T) for EVERY target set at once.

    The Moser spindle is 4-critical, so at four colours every one of its
    vertices is separable from all its non-neighbours simultaneously. This is
    the corollary that accounts for all 880 pivots of fruitless searching on
    de Grey's G, which is 5-critical, before any solver runs.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.spindle import SeparationTest

    pts, _doc = load_certificate("certificates/moser_spindle_no3coloring.json")
    g = build_graph(pts)
    assert not is_k_colorable(g, 3)[0] and is_k_colorable(g, 4)[0]
    for p in range(g.n):
        rest = build_graph([q for j, q in enumerate(pts) if j != p])
        assert is_k_colorable(rest, 3)[0], f"{p} not critical"
        others = [v for v in range(g.n) if v != p and v not in g.adj[p]]
        if not others:
            continue
        st = SeparationTest(g, 4, p, others)
        try:
            assert st.run(subset=others)[0], f"pivot {p} should be separable"
        finally:
            st.close()


def test_union_of_two_copies_is_not_critical_and_has_a_core():
    """The escape the corollary points at: in W = H union f(H), removing a
    vertex leaves a whole chi-chromatic copy, so no vertex is critical, no
    pivot can be given a colour of its own, and the full target set is a core.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.geometry import Point
    from hn.spindle import SeparationTest

    pts, _doc = load_certificate("certificates/moser_spindle_no3coloring.json")
    dx = pts[1].x - pts[0].x
    dy = pts[1].y - pts[0].y
    both = list(dict.fromkeys(pts + [Point(q.x + dx, q.y + dy) for q in pts]))
    w = build_graph(both)
    assert is_k_colorable(w, 4)[0]
    for v in range(w.n):
        rest = build_graph([q for j, q in enumerate(both) if j != v])
        if is_k_colorable(rest, 3)[0]:
            break
    else:
        # no vertex is critical, so some pivot must fail to be separable
        hub = max(range(w.n), key=lambda v: len(w.adj[v]))
        others = [v for v in range(w.n) if v != hub and v not in w.adj[hub]]
        st = SeparationTest(w, 4, hub, others)
        try:
            assert not st.run(subset=others)[0], "the full set must be a core"
        finally:
            st.close()


def test_general_block_agrees_with_the_two_sat_one():
    from hn.mixed import (blocks_targets, blocks_two_targets,
                          joint_core_configuration, joint_core_copies)

    E, pts, rho, sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    pair = [g.vertices.index(pts[2]), g.vertices.index(pts[3])]
    copies = joint_core_copies(E, rho, sigma)
    assert blocks_targets(g, bp, pair, copies) is True
    assert blocks_two_targets(g, bp, pair, copies) is True
    # one orbit alone cannot close it: three copies leave an escape
    assert not blocks_targets(g, bp, pair, copies[:3])


# -- the alignment condition ----------------------------------------------

def test_core_must_be_a_rainbow_disjoint_from_the_neighbourhood():
    """When pressure(p) = k - r exactly, a core of size r has a shape.

    In a colouring attaining the minimum the pivot has exactly r free colours
    and can take any of them, so c(T) must contain all r while holding at most
    r -- hence c(T) = free(p): the targets are a rainbow AND none of them ever
    shares a colour with a neighbour of the pivot. The built pair satisfies
    both, and its legs are the two non-pivot corners of a unit triangle, so
    they are forced apart by adjacency.
    """
    from hn.forced import core_must_be_rainbow, free_colours, pressure
    from hn.mixed import joint_core_configuration

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    legs = [g.vertices.index(pts[2]), g.vertices.index(pts[3])]
    rel = ColourRelations(g, 3)
    try:
        assert pressure(rel, bp) == 1
        assert free_colours(rel, bp) == 2 == len(legs)
        assert not rel.can_share(*legs)
        assert core_must_be_rainbow(rel, bp, legs)
        # a repeated target is one distinct point, fewer than the free
        # colours, so it cannot be a core and the test says so
        assert not core_must_be_rainbow(rel, bp, [legs[0], legs[0]])
    finally:
        rel.close()


# -- cores built forwards --------------------------------------------------

def test_core_condition_is_a_pressure_measurement():
    """T is a core of p exactly when N(p) union T uses every colour.

    If some colouring left a colour free on both, recolouring p to it stays
    proper -- properness at p asks only that its colour avoid c(N(p)) -- and
    puts c(p) outside c(T). The two characterisations must agree with the
    separation test on the built construction, where the answer is known.
    """
    from hn.forced import is_core, min_colours_on
    from hn.mixed import joint_core_configuration
    from hn.spindle import SeparationTest

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    y, z = g.vertices.index(pts[2]), g.vertices.index(pts[3])
    rel = ColourRelations(g, 3)
    try:
        assert is_core(rel, bp, [y, z])
        assert not is_core(rel, bp, [y])
        assert not is_core(rel, bp, [z])
        circle = sorted(g.adj[bp])
        assert min_colours_on(rel, circle + [y, z]) == 3
        assert min_colours_on(rel, circle + [y]) < 3
    finally:
        rel.close()
    st = SeparationTest(g, 3, bp, [y, z])
    try:
        assert not st.run(subset=[y, z])[0]      # forced, i.e. a core
        assert st.run(subset=[y])[0]             # not forced alone
        assert st.run(subset=[z])[0]
    finally:
        st.close()


def test_cegar_core_finds_the_built_pair():
    from hn.forced import cegar_core, is_core
    from hn.mixed import joint_core_configuration

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    rel = ColourRelations(g, 3)
    try:
        T, ok = cegar_core(rel, bp)
        assert ok and is_core(rel, bp, T)
        assert set(T) == {g.vertices.index(pts[2]), g.vertices.index(pts[3])}
    finally:
        rel.close()


@pytest.mark.slow
def test_cegar_core_fails_on_a_critical_graph():
    """A k-vertex-critical graph has no core, so the construction must run out.

    This is the criticality corollary as an executable statement: de Grey's G
    is 5-vertex-critical, so no amount of counterexample-killing can ever
    close a core there, and the builder has to come back empty.
    """
    from hn.forced import cegar_core

    g = build_G()
    rel = ColourRelations(g, 5)
    try:
        _T, ok = cegar_core(rel, 0, limit=40)
        assert not ok, "a critical graph must yield no core"
    finally:
        rel.close()


# -- the mechanism behind the only pressure 3 ------------------------------

def test_pressure_gadget_is_an_odd_cycle_against_k_minus_two():
    """Take the 47-vertex witness apart and check the machine causally.

    The circle splits into five hexagons, every one of the sixteen further
    points is adjacent to exactly two circle points lying in DIFFERENT
    hexagons, and those sixteen form a graph that is 3-chromatic but not
    bipartite. Squeezing the circle into two colours confines them to k-2, so
    the odd cycle bites at four colours and not at five -- and deleting three
    of them to kill the odd cycle drops the pressure from 3 to 2, which is the
    causal half rather than a coincidence of numbers.
    """
    import itertools

    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.forced import PRESSURE_GADGET, min_colours_on

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    circle = set(g.adj[piv])
    extras = [v for v in range(g.n) if v != piv and v not in circle]
    assert len(circle) == 30 and len(extras) == PRESSURE_GADGET["confined_points"]

    comps, seen = [], set()
    for s in sorted(circle):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in g.adj[x] & circle:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(set(comp))
    assert len(comps) == 5 and all(len(c) == 6 for c in comps)

    for v in extras:
        nb = [c for c in circle if c in g.adj[v]]
        assert len(nb) == 2
        assert sum(1 for c in comps if nb[0] in c) == 1
        assert not any(nb[0] in c and nb[1] in c for c in comps), \
            "both circle neighbours in one hexagon would read nothing"

    sub = build_graph([g.vertices[v] for v in extras])
    assert (sub.n, sub.m) == PRESSURE_GADGET["confined_graph"]
    assert is_k_colorable(sub, 3)[0] and not is_k_colorable(sub, 2)[0]

    def pressure_of(keep, k):
        h = build_graph([g.vertices[v] for v in sorted(keep)])
        j = h.vertices.index(g.vertices[piv])
        rel = ColourRelations(h, k)
        try:
            return min_colours_on(rel, sorted(h.adj[j])) if rel.colourable else None
        finally:
            rel.close()

    assert pressure_of(range(g.n), 4) == 3
    assert pressure_of(range(g.n), 5) == 2
    for drop in itertools.combinations(extras, 3):
        rest = [v for v in extras if v not in drop]
        if is_k_colorable(build_graph([g.vertices[v] for v in rest]), 2)[0]:
            assert pressure_of(set(range(g.n)) - set(drop), 4) == 2, \
                "killing the odd cycle must kill the pressure"
            break
    else:
        raise AssertionError("no three extras make the gadget bipartite")


def test_the_odd_cycle_covers_only_half_the_orientations():
    """The correction to the tempting one-line story.

    "An odd cycle against k-2 colours" is load-bearing but incomplete: over
    the 32 orientations of the five hexagons the confined set is bipartite in
    sixteen of them, so those are killed by something longer-range, through
    the confined points' other neighbours rather than among themselves.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    circle = set(g.adj[piv])
    extras = [v for v in range(g.n) if v != piv and v not in circle]
    seen, comps = set(), []
    for s in sorted(circle):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in g.adj[x] & circle:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(set(comp))
    par, comp_of = {}, {}
    for ci, comp in enumerate(comps):
        s = min(comp)
        par[s], front = 0, [s]
        for v in comp:
            comp_of[v] = ci
        while front:
            x = front.pop()
            for y in g.adj[x] & comp:
                if y not in par:
                    par[y] = 1 - par[x]
                    front.append(y)

    bipartite = 0
    for bits in range(1 << len(comps)):
        def col(v):
            return par[v] ^ ((bits >> comp_of[v]) & 1)
        barred = [v for v in extras
                  if len({col(c) for c in g.adj[v] & circle}) == 2]
        sub = build_graph([g.vertices[v] for v in barred])
        if sub.n == 0 or is_k_colorable(sub, 2)[0]:
            bipartite += 1
    assert bipartite == 16, f"expected half, got {bipartite}/32"


def test_the_mechanism_is_list_colouring():
    """Named completely: the auxiliary graph versus the lists the circle gives.

    Squeezing the pivot's circle into two colours hands every other vertex a
    list -- k-2 colours where it sees two differently-coloured circle points,
    k-1 where it sees two of the same. Over all 32 orientations of the five
    hexagons the auxiliary graph fails to be list-colourable in ALL of them at
    four colours, and in NONE at five, where the same lists grow to sizes 3
    and 4. That is exactly the measured pressure, 3 and then 2.
    """
    from hn.certify import load_certificate
    from hn.forced import circle_hexagons, induced_lists, list_colourable

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    comps, _par, _co = circle_hexagons(g, piv)
    assert len(comps) == 5 and all(len(c) == 6 for c in comps)
    extras = [v for v in range(g.n) if v != piv and v not in g.adj[piv]]

    for k, expected in ((4, 32), (5, 0)):
        refused = 0
        sizes = set()
        for bits in range(1 << len(comps)):
            lists, _ = induced_lists(g, piv, k, bits)
            sizes |= {len(lists[v]) for v in extras}
            if not list_colourable(g, extras, lists):
                refused += 1
        assert refused == expected, f"k={k}: {refused} refusals"
        assert sizes == {k - 2, k - 1}


def test_degrey_hexagon_offsets_are_generated_by_the_moser_angle():
    """The five hexagons sit at 0, +-theta/2 and +-(60 deg - theta), exactly,
    where cos theta = 5/6 is the Moser rotation -- the angle with
    |1 - rho|^2 = 1/3, the split-prime quotient at 3 in Q(sqrt -11).

    Checked in the field, not in degrees: the cosines between hexagon 0 and
    the rest are 1, sqrt(33)/6 twice and (5 + sqrt(33))/12 twice, and
    2*(sqrt(33)/6)^2 - 1 = 5/6 is the double-angle identity while
    (5 + sqrt 33)/12 = cos(60) cos(theta) + sin(60) sin(theta) with
    sin theta = sqrt(11)/6.
    """
    from fractions import Fraction

    from hn.certify import load_certificate
    from hn.forced import circle_hexagons

    pts, _doc = load_certificate("certificates/pressure3_witness_47.json")
    g = build_graph(pts)
    piv = 0
    comps, _par, _co = circle_hexagons(g, piv)
    pv = g.vertices[piv]
    field = pv.x.field
    reps = [(g.vertices[min(c)].x - pv.x, g.vertices[min(c)].y - pv.y)
            for c in comps]
    a = reps[0]
    cosines = [a[0] * b[0] + a[1] * b[1] for b in reps]

    half = field.sqrt(33) * field.rational(Fraction(1, 6))
    other = field.rational(Fraction(5, 12)) \
        + field.sqrt(33) * field.rational(Fraction(1, 12))
    assert sorted(map(float, cosines))[::-1] == sorted(
        map(float, [field.rational(1), half, half, other, other]))[::-1]
    assert cosines.count(half) == 2 and cosines.count(other) == 2

    cos_theta = field.rational(Fraction(5, 6))
    sin_theta = field.sqrt(11) * field.rational(Fraction(1, 6))
    assert cos_theta * cos_theta + sin_theta * sin_theta == field.rational(1)
    # theta is the Moser rotation: the chord it subtends on the unit circle
    # squares to 2 - 2 cos theta = 1/3
    two = field.rational(2)
    assert two - two * cos_theta == field.rational(Fraction(1, 3))
    # half-angle
    assert two * half * half - field.rational(1) == cos_theta
    # and 60 degrees minus theta
    assert other == field.rational(Fraction(1, 2)) * cos_theta \
        + field.sqrt(3) * field.rational(Fraction(1, 2)) * sin_theta


def test_minimise_core_shrinks_what_cegar_builds():
    """Counterexample construction is greedy, so what it returns is a core but
    rarely a small one, and the blockable sizes are exactly one, two, three."""
    from hn.forced import cegar_core, is_core, minimise_core
    from hn.mixed import joint_core_configuration

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    rel = ColourRelations(g, 3)
    try:
        T, ok = cegar_core(rel, bp)
        assert ok
        small = minimise_core(rel, bp, T)
        assert is_core(rel, bp, small)
        assert len(small) <= len(T) == 2
        for t in small:
            assert not is_core(rel, bp, [x for x in small if x != t])
    finally:
        rel.close()


@pytest.mark.slow
def test_three_hexagon_gadget_is_four_chromatic_with_no_small_core():
    """Pressure 3 permits a core of one; the gadget's smallest is seven.

    The pressure bound is necessary, not tight -- worth pinning, because it is
    the difference between "a core of one is possible here" and "there is one".
    """
    from hn.forced import cegar_core, minimise_core, pressure
    from hn.coloring import is_k_colorable
    from hn.mixed import three_hexagon_gadget

    _field, pivot, pts = three_hexagon_gadget()
    g = build_graph(pts)
    assert not is_k_colorable(g, 3)[0] and is_k_colorable(g, 4)[0]
    rel = ColourRelations(g, 4)
    try:
        best = None
        for p in range(g.n):
            if pressure(rel, p) < 3:
                continue
            T, ok = cegar_core(rel, p, limit=60)
            if ok:
                T = minimise_core(rel, p, T)
                best = len(T) if best is None else min(best, len(T))
        assert best is not None and best >= 4, \
            f"smallest core {best}; a core of 3 or less would be blockable"
    finally:
        rel.close()


def test_adding_one_point_to_a_critical_graph_gives_a_useless_core():
    """W = G + v with G k-critical: every core of v is all of V(G) \\ N(v).

    For any u outside N(v), G - u is (k-1)-colourable and so is W - u, since
    v's neighbours all lie in it; colour it with k-1 and give u the kth, and u
    is alone in that colour while none of v's neighbours carries it, so v can
    be recoloured to k as well. A target set avoiding u then misses c(v). The
    obvious resolution of the rigidity tension -- take the most rigid graph
    there is and add a single point to escape criticality -- is dead on
    arrival, since blocking caps at a core of three.
    """
    from hn.certify import load_certificate
    from hn.coloring import is_k_colorable
    from hn.forced import (cegar_core, forced_into_every_core, is_core,
                           minimise_core)
    from hn.mixed import circle_intersections

    pts, _doc = load_certificate("certificates/moser_spindle_no3coloring.json")
    g = build_graph(pts)
    one = pts[0].x.field.rational(1)
    assert all(is_k_colorable(build_graph([q for j, q in enumerate(pts)
                                           if j != i]), 3)[0]
               for i in range(g.n)), "the spindle must be 4-critical"
    v = next(q for i in range(g.n) for j in range(i + 1, g.n)
             if float(pts[i].dist2(pts[j])) <= 3.99
             for q in circle_intersections(pts[i], one, pts[j], one)
             if q not in set(pts))
    w = build_graph(pts + [v])
    piv = w.vertices.index(v)
    outside = [u for u in range(w.n) if u != piv and u not in w.adj[piv]]
    rel = ColourRelations(w, 4)
    try:
        assert rel.colourable
        T, ok = cegar_core(rel, piv, limit=40)
        assert ok
        T = minimise_core(rel, piv, T)
        assert sorted(T) == sorted(outside)
        assert sorted(forced_into_every_core(rel, piv)) == sorted(outside)
    finally:
        rel.close()


def test_unique_colourability_gives_a_core_of_k_minus_pressure():
    """The triangular lattice is uniquely 3-colourable, and that is exactly
    why sqrt(3) forces two points to agree at three colours.

    Its colouring is the Eisenstein residue modulo (1 - omega), so the colour
    classes cannot move; a pivot of pressure q therefore has a core of exactly
    k - q, one representative per class its neighbourhood misses. Measured:
    pressure 2, core of ONE, at squared distance 3 -- the classical rhombus.

    At five colours with the pressure 2 every graph here measures, the same
    statement would give a core of three, which is exactly the size the
    blocking bounds allow. That is what the remaining object has to be.
    """
    from fractions import Fraction

    from hn.field import Field
    from hn.forced import cegar_core, is_core, minimise_core, pressure
    from hn.geometry import Point

    f = Field((3,))
    r3, half = f.sqrt(3), f.rational(Fraction(1, 2))
    steps = [Point(f.rational(1), f.zero()), Point(half, r3 * half),
             Point(-half, r3 * half)]
    lat = {Point(f.zero(), f.zero())}
    for _ in range(3):
        new = set()
        for q in lat:
            for s in steps:
                new.add(Point(q.x + s.x, q.y + s.y))
                new.add(Point(q.x - s.x, q.y - s.y))
        lat |= new
    pts = sorted(lat, key=lambda q: float(q.x * q.x + q.y * q.y))
    g = build_graph(pts)
    rel = ColourRelations(g, 3)
    try:
        assert rel.colourable
        assert pressure(rel, 0) == 2
        T, ok = cegar_core(rel, 0, limit=60)
        assert ok
        T = minimise_core(rel, 0, T)
        assert len(T) == 3 - 2 == 1
        assert is_core(rel, 0, T)
        assert g.vertices[0].dist2(g.vertices[T[0]]) == f.rational(3)
    finally:
        rel.close()


def test_the_ladder_of_free_pressure_against_blockable_cores():
    """Where the method works, and where it stops, with no search involved.

    A unit circle is bipartite, so a neighbourhood containing an edge has
    pressure 2 for free and nothing more without ambient help. The pressure
    theorem makes the smallest possible core k - 2, and blocking reaches 3 and
    no further. So the free core is 1 at three colours, 2 at four, 3 at five --
    exactly saturating the blocking bound -- and 4 at six, where nothing
    blocks. Five is the last k at which the two meet.
    """
    from hn.forced import FREE_PRESSURE, LADDER, pressure
    from hn.transversal import MAX_BLOCKABLE_CORE

    assert FREE_PRESSURE == 2 and MAX_BLOCKABLE_CORE == 3
    for k in (3, 4, 5):
        assert k - FREE_PRESSURE <= MAX_BLOCKABLE_CORE, f"k={k} should work"
    assert 6 - FREE_PRESSURE > MAX_BLOCKABLE_CORE, "six must be out of reach"
    assert set(LADDER) == {3, 4, 5, 6}

    # and the free pressure is really free: any neighbourhood with an edge
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        for v in range(g.n):
            circle = set(g.adj[v])
            has_edge = any(b in g.adj[a] for a in circle for b in circle
                           if b != a)
            if has_edge:
                assert pressure(rel, v) >= FREE_PRESSURE
    finally:
        rel.close()
